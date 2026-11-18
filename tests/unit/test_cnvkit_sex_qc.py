"""L1: cnvkit_sex_qc.py — sex resolution (DI-4098/4099/4100, D9) and the per-sample sex QC (DI-4103, D8),
including a differential against the pinned reference `sex_check.py` (31a62fc)."""
import csv
import subprocess
import sys

import pytest

from cnvfiles import BIN, REFERENCE

pytestmark = pytest.mark.unit

SCRIPT = BIN / "cnvkit_sex_qc.py"
M = "140609591-25328S0023-26CGPV2-8475-M-96770770"
F = "138076232-25178S0009-26CGPV1-9197-F-96770770"
U = "135024088-25021S0072-26CGPVal1SHv150kb-9999-U-1111"


def run(*argv):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, argv)], capture_output=True, text=True)


def resolve(sample_id, bam_name="", override=""):
    p = run("resolve", "--sample-id", sample_id, "--bam-name", bam_name, "--override", override)
    assert p.returncode == 0, p.stderr
    return tuple(p.stdout.rstrip("\n").split("\t"))


# ── resolve ──────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("sample_id, bam, override, expected", [
    (M, f"{M}.bam", "", ("male", "sample_name", "M", "M")),
    (F, f"{F}.bam", "", ("female", "sample_name", "F", "F")),
    (U, f"{U}.bam", "", ("", "cnvkit_inference", "U", "U")),
    # field 5 present but not M/F: U, no fallback to the BAM name
    ("1-2-3-4-X-6", f"{M}.bam", "", ("", "cnvkit_inference", "U", "U")),
    ("1-2-3-4-m-6", f"{M}.bam", "", ("", "cnvkit_inference", "U", "U")),   # case-sensitive, as sex_check.py
    # sample ID without a 5th field: the BAM name stem is used
    ("25328S0023", f"{M}.bam", "", ("male", "bam_name", "M", "M")),
    ("25328S0023", f"{F}_markdup.bam", "", ("female", "bam_name", "F", "F")),
    ("25328S0023", "sample.bam", "", ("", "cnvkit_inference", "U", "U")),
    ("25328S0023", "", "", ("", "cnvkit_inference", "U", "U")),
    # override wins and is reported; LIMS is still reported
    (U, "", "male", ("male", "override", "U", "M")),
    (M, "", "F", ("female", "override", "M", "F")),
    (M, "", " Female ", ("female", "override", "M", "F")),
    (F, "", "m", ("male", "override", "F", "M")),
    (M, "", "", ("male", "sample_name", "M", "M")),
])
def test_resolve(sample_id, bam, override, expected):
    assert resolve(sample_id, bam, override) == expected


@pytest.mark.parametrize("bad", ["X", "U", "unknown", "males", "1", "Male, female"])
def test_resolve_rejects_bad_override(bad):
    p = run("resolve", "--sample-id", M, "--bam-name", "", "--override", bad)
    assert p.returncode == 2
    assert "not accepted" in p.stderr


# ── qc ───────────────────────────────────────────────────────────────────────

def sex_table(path, sample_file, sex, x, y):
    path.write_text(f"sample\tsex\tX_logratio\tY_logratio\n{sample_file}\t{sex}\t{x}\t{y}\n")
    return path


def qc(tmp_path, sample, raw_y, lims, override="U", source=None, cnr_sex="Male", raw_sex="Male", **extra):
    raw = sex_table(tmp_path / "raw.tsv", f"{sample}.targetcoverage.cnn", raw_sex, "-0.97", raw_y)
    cnr = sex_table(tmp_path / "cnr.tsv", f"{sample}.cnr", cnr_sex, "-0.94", "nan")
    source = source or ("override" if override != "U" else "sample_name" if lims != "U" else "cnvkit_inference")
    args = dict(par_fix_applied="1", par_span_segments="0", dropout_probes="2")
    args.update(extra)
    p = run("qc", "--sample-id", sample, "--raw-sex", raw, "--cnr-sex", cnr, "--lims-sex", lims,
            "--override-sex", override, "--applied-source", source, "--par-fix-applied", args["par_fix_applied"],
            "--par-span-segments", args["par_span_segments"], "--dropout-probes", args["dropout_probes"],
            "--out", tmp_path / "s")
    return p, tmp_path / "s.sex_qc.tsv", tmp_path / "s.sex_check_mqc.tsv"


def row(path):
    return next(csv.DictReader(open(path), delimiter="\t"))


@pytest.mark.parametrize("raw_y, present", [("-1.13", 1), ("-4.99", 1), ("-5", 0), ("-5.0", 0), ("-19.9", 0),
                                            ("nan", 0), ("", 0)])
def test_chry_presence(tmp_path, raw_y, present):
    p, out, _ = qc(tmp_path, U, raw_y, "U")
    assert p.returncode == 0, p.stderr
    r = row(out)
    assert r["chrY_present"] == str(present) and r["chrY_sex"] == ("Male" if present else "Female")


@pytest.mark.parametrize("lims, override, raw_y, mismatch", [
    ("M", "U", "-1.1", "0"), ("M", "U", "-20", "1"), ("F", "U", "-20", "0"), ("F", "U", "-1.1", "1"),
    ("U", "U", "-1.1", "NA"), ("U", "U", "-20", "NA"),
    ("U", "M", "-1.1", "0"), ("M", "F", "-1.1", "1"), ("F", "M", "-1.1", "0"),   # override is the declared sex
])
def test_sex_mismatch(tmp_path, lims, override, raw_y, mismatch):
    p, out, _ = qc(tmp_path, M, raw_y, lims, override)
    assert p.returncode == 0, p.stderr
    assert row(out)["sex_mismatch"] == mismatch


def test_qc_columns_and_applied_fields(tmp_path):
    p, out, mqc = qc(tmp_path, U, "-1.13", "U", cnr_sex="Male", par_span_segments="1")
    assert p.returncode == 0, p.stderr
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == ["run", "sample", "lims_sex", "applied_sex", "applied_source", "raw_Y_logratio",
                                    "chrY_present", "chrY_sex", "raw_X_logratio", "raw_cnvkit_sex", "cnr_X_logratio",
                                    "cnr_cnvkit_sex", "sex_mismatch", "par_fix_applied", "par_span_segments",
                                    "dropout_probes"]
    r = row(out)
    assert (r["run"], r["applied_sex"], r["applied_source"]) == ("", "cnvkit_inferred(Male)", "cnvkit_inference")
    assert (r["par_fix_applied"], r["par_span_segments"], r["dropout_probes"]) == ("1", "1", "2")
    assert p.stdout.startswith("[sex] chrY Male")


def test_mqc_format(tmp_path):
    _, _, mqc = qc(tmp_path, M, "-1.1", "M")
    text = mqc.read_text().splitlines()
    assert text[:4] == ["# id: 'cnvkit_sex_check'", "# section_name: 'CNVkit sex check'",
                        "# description: 'LIMS sex vs chrY presence in raw CNVkit coverage (sex_mismatch = 1 needs "
                        "review); X log-ratios from raw coverage and from the .cnr for context.'",
                        "# plot_type: 'table'"]
    assert text[4] == ("Sample\tlims_sex\tchrY_sex\tsex_mismatch\traw_Y_logratio\traw_X_logratio\tcnr_X_logratio\t"
                       "applied_sex\tdropout_probes")
    assert text[5] == f"{M}\tM\tMale\t0\t-1.1\t-0.97\t-0.94\tmale\t2"
    assert len(text) == 6


@pytest.mark.parametrize("bad", [
    dict(sample_file="other.targetcoverage.cnn"),
    dict(two_rows=True),
])
def test_qc_rejects_wrong_table(tmp_path, bad):
    raw = tmp_path / "raw.tsv"
    name = bad.get("sample_file", f"{M}.targetcoverage.cnn")
    body = f"{name}\tMale\t-1\t-1\n" * (2 if bad.get("two_rows") else 1)
    raw.write_text("sample\tsex\tX_logratio\tY_logratio\n" + body)
    cnr = sex_table(tmp_path / "cnr.tsv", f"{M}.cnr", "Male", "-1", "nan")
    p = run("qc", "--sample-id", M, "--raw-sex", raw, "--cnr-sex", cnr, "--lims-sex", "M", "--override-sex", "U",
            "--applied-source", "sample_name", "--par-fix-applied", "1", "--par-span-segments", "0",
            "--dropout-probes", "0", "--out", tmp_path / "s")
    assert p.returncode == 2


# ── differential against sex_check.py ───────────────────────────────────────

CASES = [(M, "-1.13", "-0.97", "Male", "Male"), (M, "-21.3", "-0.05", "Female", "Female"),
         (F, "-29.3", "-0.07", "Female", "Female"), (F, "-2.2", "-0.9", "Male", "Male"),
         (U, "-1.13", "-0.97", "Male", "Male"), (U, "-25", "0.01", "Female", "Female"),
         (U, "nan", "0.01", "Female", "Male"), ("1-2-3-4-m-6", "-1.0", "-1.0", "Male", "Male")]


@pytest.mark.parametrize("sample, raw_y, raw_x, raw_sex, cnr_sex", CASES)
def test_matches_reference(tmp_path, sample, raw_y, raw_x, raw_sex, cnr_sex):
    run_name = "0797_V2"
    raw = tmp_path / "raw.tsv"
    raw.write_text(f"sample\tsex\tX_logratio\tY_logratio\n{run_name}__{sample}.targetcoverage.cnn\t{raw_sex}\t{raw_x}\t{raw_y}\n")
    cnr = tmp_path / "cnr.tsv"
    cnr.write_text(f"sample\tsex\tX_logratio\tY_logratio\n{run_name}__{sample}.cnr\t{cnr_sex}\t-0.94\tnan\n")
    subprocess.run([sys.executable, str(REFERENCE / "sex_check.py"), "--raw-sex", raw, "--cnr-sex", cnr,
                    "--out", tmp_path / "ref"], check=True, capture_output=True)
    sex_arg, source, lims, declared = resolve(sample)
    p = run("qc", "--sample-id", sample, "--raw-sex", raw, "--cnr-sex", cnr, "--lims-sex", lims, "--override-sex", "U",
            "--applied-source", source, "--par-fix-applied", "0", "--par-span-segments", "0", "--dropout-probes", "0",
            "--out", tmp_path / "port", "--run", run_name)
    assert p.returncode == 0, p.stderr
    ref, mine = row(tmp_path / "ref.sex_qc.tsv"), row(tmp_path / "port.sex_qc.tsv")
    # D8: applied_sex / applied_source come from the app's resolution; every other sex_check.py column is unchanged
    assert {k: v for k, v in mine.items() if k in ref and k not in ("applied_sex", "applied_source")} == \
           {k: v for k, v in ref.items() if k not in ("applied_sex", "applied_source")}
    assert mine["applied_source"] == {"LIMS": "sample_name", "CNVkit inference (.cnr)": "cnvkit_inference"}[ref["applied_source"]]
    assert mine["applied_sex"] == (ref["applied_sex"] if ref["applied_sex"] != "cnvkit_inferred"
                                   else f"cnvkit_inferred({cnr_sex})")
    # MultiQC: same header lines and columns, plus dropout_probes
    ref_mqc = (tmp_path / "ref.sex_check_mqc.tsv").read_text().splitlines()
    port_mqc = (tmp_path / "port.sex_check_mqc.tsv").read_text().splitlines()
    assert port_mqc[:4] == ref_mqc[:4]
    assert port_mqc[4] == ref_mqc[4] + "\tdropout_probes"
    ref_vals, port_vals = ref_mqc[5].split("\t"), port_mqc[5].split("\t")
    assert port_vals[:7] == ref_vals[:7] and port_vals[-1] == "0"
