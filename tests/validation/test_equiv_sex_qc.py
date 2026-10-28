"""L3: cnvkit_sex_qc.py on the 226 validation sequencings (UMI-consensus and non-UMI), in the cgp-cnvkit:1.0.0 image.

`cnvkit.py sex --diploid-parx-genome grch38` is run on every .targetcoverage.cnn and .cnr. The pinned reference
`sex_check.py` on those tables must reproduce the validated sex QC (UMI: final_config.sex_qc.tsv; non-UMI: the raw
tables equal tier1 sex_raw). The per-sample port must then equal the reference in every sex_check.py column, except
applied_sex / applied_source, which follow the app's resolution (D8); sex_mismatch = 1 in 0 of the 173 LIMS M/F.
"""
import argparse
import csv
import shlex
import subprocess

import pytest

from test_equiv_parfix import docker_available, in_image
from vdata import BIN, DATA, REFERENCE, REPO, load_module, require, sequencings

pytestmark = [pytest.mark.validation, pytest.mark.docker]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def batch_sex(tmp_path, files, out):
    """cnvkit.py sex on {(run, sample): path}; returns {(run, sample): row}. Replicates share a sample name across
    runs, so rows are keyed by the full path, not the file name."""
    listing = tmp_path / f"{out}.files"
    listing.write_text("\n".join(str(f) for f in files.values()) + "\n")
    in_image(tmp_path, f"xargs -a {listing} cnvkit.py sex --diploid-parx-genome grch38 -o {tmp_path / out}")
    by_path = {r["sample"]: r for r in read_tsv(tmp_path / out)}
    return {key: by_path[str(path)] for key, path in files.items()}


def one_row_table(path, sample_file, r):
    path.write_text("sample\tsex\tX_logratio\tY_logratio\n"
                    f"{sample_file}\t{r['sex']}\t{r['X_logratio']}\t{r['Y_logratio']}\n")
    return path


@pytest.mark.skipif(not docker_available(), reason="docker or cgp-cnvkit:1.0.0 not available")
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


@pytest.mark.skipif(not docker_available(), reason="docker or cgp-cnvkit:1.0.0 not available")
@pytest.mark.parametrize("dataset", ["umi", "nonumi"])
def test_sex_qc_equivalence(tmp_path, dataset):
    port = load_module(BIN / "cnvkit_sex_qc.py", "cnvkit_sex_qc_port")
    seqs = sequencings()
    cnn = {(run, s): require(DATA / dataset / "cnn" / run / f"{s}.targetcoverage.cnn") for run, s in seqs}
    cnr = {(run, s): DATA / dataset / "calls" / "A0" / run / f"{s}.cnr" for run, s in seqs}
    raw = batch_sex(tmp_path, cnn, "raw.tsv")
    cnrs = batch_sex(tmp_path, cnr, "cnr.tsv")

    # the reference on run-prefixed tables, as in the validation
    with open(tmp_path / "ref_raw.tsv", "w") as a, open(tmp_path / "ref_cnr.tsv", "w") as b:
        a.write("sample\tsex\tX_logratio\tY_logratio\n"); b.write("sample\tsex\tX_logratio\tY_logratio\n")
        for run, s in seqs:
            r, c = raw[(run, s)], cnrs[(run, s)]
            a.write(f"{run}__{s}.targetcoverage.cnn\t{r['sex']}\t{r['X_logratio']}\t{r['Y_logratio']}\n")
            b.write(f"{run}__{s}.cnr\t{c['sex']}\t{c['X_logratio']}\t{c['Y_logratio']}\n")
    subprocess.run(["python3", str(REFERENCE / "sex_check.py"), "--raw-sex", tmp_path / "ref_raw.tsv",
                    "--cnr-sex", tmp_path / "ref_cnr.tsv", "--out", tmp_path / "ref"], check=True, capture_output=True)
    ref = {(r["run"], r["sample"]): r for r in read_tsv(tmp_path / "ref.sex_qc.tsv")}
    assert set(ref) == set(seqs)

    if dataset == "umi":     # the reference on our cnvkit.py sex tables reproduces the validated sex QC exactly
        assert (tmp_path / "ref.sex_qc.tsv").read_text() == require(DATA / "umi" / "sex_qc.tsv").read_text()
    else:                    # the raw-coverage sex calls equal the DI-4118 tier 1 tables
        for run in sorted({run for run, _ in seqs}):
            for r in read_tsv(require(DATA / "nonumi" / "sex_raw" / run / "sex.tsv")):
                s = r["sample"].split("__", 1)[1][:-len(".targetcoverage.cnn")]
                mine = raw[(run, s)]
                assert (mine["sex"], mine["X_logratio"], mine["Y_logratio"]) == \
                       (r["sex"], r["X_logratio"], r["Y_logratio"]), (run, s)

    failures, mismatches, lims_mf = [], 0, 0
    for run, s in seqs:
        sex_arg, source, lims, declared = port.resolve(s, f"{s}.bam", "")
        d = tmp_path / "per_sample" / run
        d.mkdir(parents=True, exist_ok=True)
        a = argparse.Namespace(
            sample_id=s, run=run, lims_sex=lims, override_sex="U", applied_source=source,
            raw_sex=one_row_table(d / f"{s}.sex_raw.tsv", f"{s}.targetcoverage.cnn", raw[(run, s)]),
            cnr_sex=one_row_table(d / f"{s}.sex_cnr.tsv", f"{s}.cnr", cnrs[(run, s)]),
            par_fix_applied="0", par_span_segments=0, dropout_probes=0)
        mine, theirs = port.qc_row(a), ref[(run, s)]
        for k, v in theirs.items():
            if k == "applied_source":
                expect = {"LIMS": "sample_name", "CNVkit inference (.cnr)": "cnvkit_inference"}[v]
            elif k == "applied_sex" and v == "cnvkit_inferred":
                expect = f"cnvkit_inferred({theirs['cnr_cnvkit_sex']})"
            else:
                expect = v
            if str(mine[k]) != expect:
                failures.append(f"{run}/{s}: {k} {mine[k]!r} != {expect!r}")
        if lims in ("M", "F"):
            lims_mf += 1
            mismatches += mine["sex_mismatch"] == "1"
        assert sex_arg == {"M": "male", "F": "female", "U": ""}[lims]
    assert not failures, "\n".join(failures[:20])
    assert (lims_mf, mismatches) == (173, 0)
