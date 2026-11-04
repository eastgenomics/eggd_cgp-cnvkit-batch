"""L3: cnvkit_sex_qc.py on the 226 validation sequencings (UMI-consensus and non-UMI), in the cgp-cnvkit:1.0.0 image,
through the command lines pipeline.sh runs.

Each sequencing (keyed by run and sample: replicates share a sample ID across runs) gets its own directory with copies
of its <s>.targetcoverage.cnn and <s>.cnr (arm A0), where pipeline.sh's commands run with its file names:
  cnvkit.py sex <s>.targetcoverage.cnn --diploid-parx-genome grch38 -o <s>.sex_raw.tsv
  cnvkit.py sex <s>.cnr --diploid-parx-genome grch38 -o <s>.sex_cnr.tsv
  cnvkit_sex_qc.py resolve --sample-id <s> --bam-name <s>.bam --override ''
  cnvkit_parfix.py decide ... ; for a treated male cnvkit_parfix.py shift ... (par_fix_applied, par_span_segments)
  cnvkit_sex_qc.py qc --sample-id <s> --raw-sex <s>.sex_raw.tsv --cnr-sex <s>.sex_cnr.tsv ... --out <s>
The pinned reference `sex_check.py`, on run-prefixed tables built from those per-sample tables as in the validation,
must reproduce the validated sex QC (UMI: final_config.sex_qc.tsv; non-UMI: the raw tables equal tier1 sex_raw). The
port's <s>.sex_qc.tsv must then equal the reference in every sex_check.py column, except applied_sex /
applied_source, which follow the app's resolution (D8), and run, which the app leaves empty (pipeline.sh passes no
--run); its added columns are the decide / shift values and 0 masked probes (arm A0: no mask). Both tables must
have their full column lists (round 7, W4). sex_mismatch = 1 in 0 of the 173 LIMS M/F.
"""
import csv
import shlex
import shutil
import subprocess

import pytest

from test_equiv_parfix import image_jobs, in_image
from vdata import DATA, REFERENCE, REPO, SHIFT_KEYS, covers, index, parallel, require, script, sequencings, stats

pytestmark = [pytest.mark.validation, pytest.mark.docker]

SEX_HEADER = "sample\tsex\tX_logratio\tY_logratio\n"
# the columns sex_check.py (31a62fc) writes, and the port's, in order (round 7, W4: both are required in full)
SEX_CHECK_COLUMNS = ["run", "sample", "lims_sex", "applied_sex", "applied_source", "raw_Y_logratio", "chrY_present",
                     "chrY_sex", "raw_X_logratio", "raw_cnvkit_sex", "cnr_X_logratio", "cnr_cnvkit_sex", "sex_mismatch"]
ADDED_COLUMNS = ["par_fix_applied", "par_span_segments", "dropout_probes"]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_no_parx_raw_and_cnr_inference_are_distinct(tmp_path):
    # A no-PARX fix run on a LIMS-M sample: the raw-coverage and .cnr sex calls are not interchangeable.
    rows = read_tsv(REPO / "tests" / "harness" / "local_samples.tsv")
    row = next(r for r in rows if "25328S0023" in r["sample"])
    cnn = require(DATA / "nonumi" / "cnn" / row["run"] / f"{row['sample']}.targetcoverage.cnn")
    pon = require(DATA / "nonumi" / "pons" / f"{row['cn_reference']}.cnn")
    empty = tmp_path / "empty.cnn"
    empty.write_text("chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n")
    cnr, raw_sex, cnr_sex = (tmp_path / name for name in ("no-parx.cnr", "raw-sex.tsv", "cnr-sex.tsv"))
    q = shlex.quote
    in_image(tmp_path, f"""set -e
        cnvkit.py fix {q(str(cnn))} {q(str(empty))} {q(str(pon))} -o {q(str(cnr))}
        cnvkit.py sex {q(str(cnn))} -o {q(str(raw_sex))}
        cnvkit.py sex {q(str(cnr))} -o {q(str(cnr_sex))}
    """)
    assert read_tsv(raw_sex)[0]["sex"] == "Female"
    assert read_tsv(cnr_sex)[0]["sex"] == "Male"


def one_row(path, sample_file):
    """The one row of a per-sample `cnvkit.py sex` table, which must name sample_file (pipeline.sh's file name)."""
    rows = read_tsv(path)
    if [r["sample"] for r in rows] != [sample_file]:
        raise AssertionError(f"{path}: rows {rows}, expected one for {sample_file}")
    return rows[0]


@pytest.mark.parametrize("dataset", ["umi", "nonumi"])
def test_sex_qc_equivalence(tmp_path, dataset):
    seqs = sequencings()
    wd = {(run, s): tmp_path / run / s for run, s in seqs}
    jobs = []
    for run, s in seqs:
        d = wd[(run, s)]
        d.mkdir(parents=True)
        shutil.copyfile(require(DATA / dataset / "cnn" / run / f"{s}.targetcoverage.cnn"), d / f"{s}.targetcoverage.cnn")
        for ext in (".cnr", ".call.cns"):
            shutil.copyfile(DATA / dataset / "calls" / "A0" / run / f"{s}{ext}", d / f"{s}{ext}")
        px = ["--diploid-parx-genome", "grch38"]
        jobs.append((d, shlex.join(["cnvkit.py", "sex", f"{s}.targetcoverage.cnn", *px, "-o", f"{s}.sex_raw.tsv"])
                     + " && " + shlex.join(["cnvkit.py", "sex", f"{s}.cnr", *px, "-o", f"{s}.sex_cnr.tsv"])))
    image_jobs(tmp_path, jobs, "cnvkit.log")
    raw = {k: one_row(wd[k] / f"{k[1]}.sex_raw.tsv", f"{k[1]}.targetcoverage.cnn") for k in seqs}
    cnrs = {k: one_row(wd[k] / f"{k[1]}.sex_cnr.tsv", f"{k[1]}.cnr") for k in seqs}

    # the reference on run-prefixed tables, as in the validation
    with open(tmp_path / "ref_raw.tsv", "w") as a, open(tmp_path / "ref_cnr.tsv", "w") as b:
        a.write(SEX_HEADER); b.write(SEX_HEADER)
        for run, s in seqs:
            r, c = raw[(run, s)], cnrs[(run, s)]
            a.write(f"{run}__{s}.targetcoverage.cnn\t{r['sex']}\t{r['X_logratio']}\t{r['Y_logratio']}\n")
            b.write(f"{run}__{s}.cnr\t{c['sex']}\t{c['X_logratio']}\t{c['Y_logratio']}\n")
    subprocess.run(["python3", str(REFERENCE / "sex_check.py"), "--raw-sex", tmp_path / "ref_raw.tsv",
                    "--cnr-sex", tmp_path / "ref_cnr.tsv", "--out", tmp_path / "ref"], check=True, capture_output=True,
                   timeout=600)
    ref = index(read_tsv(tmp_path / "ref.sex_qc.tsv"), lambda r: (r["run"], r["sample"]), "reference sex_qc")
    covers(ref, seqs, "reference sex_qc")
    assert all(list(r) == SEX_CHECK_COLUMNS for r in ref.values())

    if dataset == "umi":     # the reference on our cnvkit.py sex tables reproduces the validated sex QC exactly
        assert (tmp_path / "ref.sex_qc.tsv").read_text() == require(DATA / "umi" / "sex_qc.tsv").read_text()
    else:                    # the raw-coverage sex calls equal the DI-4118 tier 1 tables
        tier1 = {}
        for run in sorted({run for run, _ in seqs}):
            tier1.update(index(read_tsv(require(DATA / "nonumi" / "sex_raw" / run / "sex.tsv")),
                               lambda r, run=run: (run, r["sample"].split("__", 1)[1][:-len(".targetcoverage.cnn")]),
                               f"sex_raw/{run}"))
        covers(tier1, seqs, "nonumi sex_raw")
        for key, r in tier1.items():
            mine = raw[key]
            assert (mine["sex"], mine["X_logratio"], mine["Y_logratio"]) == \
                   (r["sex"], r["X_logratio"], r["Y_logratio"]), key

    def per_sample(key):
        """pipeline.sh's sex resolution, PAR decision / shift and sex QC for one sequencing; (failures, row)."""
        run, s = key
        d, fails = wd[key], []
        p = script("cnvkit_sex_qc.py", "resolve", "--sample-id", s, "--bam-name", f"{s}.bam", "--override", "",
                   cwd=d)
        if p.returncode:
            return [f"{run}/{s}: resolve: {p.stderr.strip()}"], None
        sex_arg, source, lims, declared = p.stdout.rstrip("\n").split("\t")
        if sex_arg != {"M": "male", "F": "female", "U": ""}[lims] or declared != lims:
            fails.append(f"resolve {p.stdout!r}")
        p = script("cnvkit_parfix.py", "decide", "--sex-arg", sex_arg, "--cnr-sex", f"{s}.sex_cnr.tsv",
                   "--parx", "grch38", cwd=d)
        if p.returncode:
            return [f"{run}/{s}: decide: {p.stderr.strip()}"], None
        applied, span = p.stdout.split("\t", 1)[0], "0"
        if applied == "1":
            p = script("cnvkit_parfix.py", "shift", f"{s}.cnr", f"{s}.call.cns", "parfix", cwd=d)
            if p.returncode:
                return [f"{run}/{s}: shift: {p.stderr.strip()}"], None
            span = stats(p.stdout, SHIFT_KEYS)["span_segs"]
        p = script("cnvkit_sex_qc.py", "qc", "--sample-id", s, "--raw-sex", f"{s}.sex_raw.tsv",
                   "--cnr-sex", f"{s}.sex_cnr.tsv", "--lims-sex", lims, "--override-sex", "U",
                   "--applied-source", source, "--par-fix-applied", applied, "--par-span-segments", span,
                   "--dropout-probes", "0", "--out", s, cwd=d)
        if p.returncode:
            return [f"{run}/{s}: qc: {p.stderr.strip()}"], None
        rows = read_tsv(d / f"{s}.sex_qc.tsv")
        if len(rows) != 1 or not (d / f"{s}.sex_check_mqc.tsv").exists():
            return [f"{run}/{s}: qc wrote {len(rows)} rows"], None
        mine, theirs = rows[0], ref[key]
        if list(mine) != SEX_CHECK_COLUMNS + ADDED_COLUMNS:    # every column, so none is compared by omission
            return [f"{run}/{s}: qc columns {list(mine)}"], None
        for k, v in theirs.items():
            if k == "applied_source":
                expect = {"LIMS": "sample_name", "CNVkit inference (.cnr)": "cnvkit_inference"}[v]
            elif k == "applied_sex" and v == "cnvkit_inferred":
                expect = f"cnvkit_inferred({theirs['cnr_cnvkit_sex']})"
            elif k == "run":
                expect = ""
            else:
                expect = v
            if mine[k] != expect:
                fails.append(f"{k} {mine[k]!r} != {expect!r}")
        for k, expect in (("par_fix_applied", applied), ("par_span_segments", span), ("dropout_probes", "0")):
            if mine[k] != expect:
                fails.append(f"{k} {mine[k]!r} != {expect!r}")
        return [f"{run}/{s}: {f}" for f in fails], mine

    results = parallel(per_sample, seqs)
    failures = [f for fails, _ in results for f in fails]
    assert not failures, "\n".join(failures[:20])
    lims_mf = [row for _, row in results if row["lims_sex"] in ("M", "F")]
    assert (len(lims_mf), sum(row["sex_mismatch"] == "1" for row in lims_mf)) == (173, 0)
