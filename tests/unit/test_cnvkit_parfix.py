"""L1: cnvkit_parfix.py (DI-4101) — PAR geometry, shift, male decision, compare; differential against the pinned
reference `genemetrics_parfix.py` (f933f7a) on synthetic inputs."""
import os
import random
import re
import subprocess
import sys

import pytest

from cnvfiles import BIN, REFERENCE
from procutil import SCRIPT_TIMEOUT

pytestmark = pytest.mark.unit

SCRIPT = BIN / "cnvkit_parfix.py"
CNR_H = "chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n"
CNS_H = "chromosome\tstart\tend\tgene\tlog2\tcn\tdepth\tprobes\tweight\n"
GM_H = "gene\tchromosome\tstart\tend\tlog2\tdepth\tweight\tcn\tprobes\tsegment_weight\tsegment_probes\n"


def load_port():
    import importlib.util
    spec = importlib.util.spec_from_file_location("cnvkit_parfix", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load_reference():
    """The reference module (top level: par_state, compare, ...) plus shift_par, which only exists inside its
    in-image WORKER source: extracted and exec'd without the cnvlib import."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("genemetrics_parfix_ref", REFERENCE / "genemetrics_parfix.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    src = re.search(r"(def shift_par\(.*?)\njob = json\.load", m.WORKER, re.S).group(1)
    ns = {"par_state": m.par_state}
    exec(src, ns)
    m.shift_par = ns["shift_par"]
    # the pinned driver's own subprocesses (docker run, cnvkit.py sex) have no timeout; the tests use only its pure
    # functions, so any call into its subprocess module fails at once instead of hanging
    m.subprocess = NoSubprocess()
    return m


class NoSubprocess:
    def __getattr__(self, name):
        raise AssertionError(f"reference genemetrics_parfix.py: subprocess.{name} is not allowed in the tests")


port, ref = load_port(), load_reference()


def python(*args, env=None):
    """sys.executable with args, bounded; env: extra environment variables."""
    return subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
                          env={**os.environ, **(env or {})})


def run(*argv):
    return python(SCRIPT, *argv)


# ── PAR geometry ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("chrom, start, end, state", [
    ("chrX", 10000, 10100, "in"), ("chrX", 9999, 10100, "span"), ("chrX", 0, 10000, "out"),
    ("chrX", 2781379, 2781479, "in"), ("chrX", 2781379, 2781480, "span"), ("chrX", 2781479, 2781600, "out"),
    ("chrX", 155701382, 155701500, "in"), ("chrX", 155701381, 155701500, "span"),
    ("chrX", 156030800, 156030895, "in"), ("chrX", 156030800, 156030896, "span"), ("chrX", 156030895, 156031000, "out"),
    ("chrX", 1190891, 1212654, "in"), ("chrX", 5000000, 6000000, "out"), ("chrX", 0, 200000000, "span"),
    ("chrY", 10000, 10100, "out"), ("chr1", 10000, 10100, "out"),
])
def test_par_state(chrom, start, end, state):
    assert port.par_state(chrom, start, end) == state == ref.par_state(chrom, start, end)


# ── shift ────────────────────────────────────────────────────────────────────

def test_shift(tmp_path):
    cnr = tmp_path / "s.cnr"
    cnr.write_text(CNR_H + "chr1\t10000\t10100\tA\t0.5\t100\t1\n"
                           "chrX\t1190891\t1191000\tCRLF2\t1.26805\t100\t1\n"
                           "chrX\t2781400\t2781600\tSPAN\t0.9\t100\t1\n"
                           "chrX\t5000000\t5000144\tNONPAR\t0.1\t100\t1\n"
                           "chrX\t155701400\t155701500\tPAR2\t-0.25\t100\t1\n")
    cns = tmp_path / "s.call.cns"
    cns.write_text(CNS_H + "chrX\t1190891\t2781479\tCRLF2\t1.26805\t3\t100\t8\t7.2\n")
    proc = run("shift", cnr, cns, tmp_path / "parfix")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "par_bins=2 span_bins=1 par_segs=1 span_segs=0"
    out = (tmp_path / "parfix" / "s.cnr").read_text().splitlines()
    assert out[0] + "\n" == CNR_H
    assert out[1] == "chr1\t10000\t10100\tA\t0.5\t100\t1"                 # untouched, byte-identical
    assert out[2] == f"chrX\t1190891\t1191000\tCRLF2\t{1.26805 - 1.0!r}\t100\t1"   # repr, as the reference
    assert repr(1.26805 - 1.0) == "0.2680499999999999"
    assert out[3] == "chrX\t2781400\t2781600\tSPAN\t0.9\t100\t1"           # spans PAR1 end: unchanged
    assert out[5] == "chrX\t155701400\t155701500\tPAR2\t-1.25\t100\t1"
    assert (tmp_path / "parfix" / "s.call.cns").read_text().splitlines()[1].split("\t")[4] == repr(1.26805 - 1.0)


# ── decide ───────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("sex_arg, inferred, parx, male", [
    ("male", "Female", "grch38", True),
    ("female", "Male", "grch38", False),
    ("", "Male", "grch38", True),
    ("", "Female", "grch38", False),
    ("male", "Male", "", False),       # no PARX: the s0 arms, post-step not applicable
    ("", "Male", "", False),
])
def test_decide(sex_arg, inferred, parx, male):
    assert port.decide(sex_arg, inferred, parx)[0] is male


def test_decide_cli(tmp_path):
    sex = tmp_path / "sex_cnr.tsv"
    sex.write_text("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\t-0.95\tnan\n")
    proc = run("decide", "--sex-arg", "", "--cnr-sex", sex, "--parx", "grch38")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == "1\tno -x; cnvkit.py sex on the .cnr: Male\n"
    assert run("decide", "--sex-arg", "female", "--cnr-sex", sex, "--parx", "grch38").stdout.startswith("0\t")
    bad = run("decide", "--sex-arg", "M", "--cnr-sex", sex, "--parx", "grch38")
    assert bad.returncode != 0 and "male, female or empty" in bad.stderr


def test_decide_cli_rejects_multi_sample_table(tmp_path):
    sex = tmp_path / "sex_cnr.tsv"
    sex.write_text("sample\tsex\tX_logratio\tY_logratio\na.cnr\tMale\t-1\tnan\nb.cnr\tFemale\t0\tnan\n")
    proc = run("decide", "--sex-arg", "", "--cnr-sex", sex, "--parx", "grch38")
    assert proc.returncode != 0 and "expected one sample row" in proc.stderr


@pytest.mark.parametrize("body, message", [
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\t-0.95\n", "the sample row has 3 fields"),   # short row
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\n", "the sample row has 1 fields"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tgarbage\t-0.95\tnan\n", "is not one of Male, Female, Unknown"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tmale\t-0.95\tnan\n", "is not one of Male, Female, Unknown"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tunknown\tNA\tNA\n", "is not one of Male, Female, Unknown"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\tNA\tNA\n", "X_logratio 'NA' is not a finite number"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\t-0.95\tNA\n", "Y_logratio 'NA' is not a finite number"),
    # R3: the same holes as in cnvkit_sex_qc.py: infinities, overflow, duplicate columns
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\t-inf\tnan\n", "X_logratio '-inf' is not a finite number"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\t-0.95\tinf\n", "Y_logratio 'inf' is not a finite number"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\t-1e999\tnan\n", "X_logratio '-1e999' is not a finite"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\tgarbage\tnan\n", "X_logratio 'garbage' is not a finite"),
    ("sample\tsex\tsex\tX_logratio\tY_logratio\ns.cnr\tgarbage\tMale\t-0.95\tnan\n", "duplicate column(s) sex"),
    ("sample\tsex\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\tgarbage\t-0.95\tnan\n", "duplicate column(s) sex"),
    ("sample\tsex\tX_logratio\ns.cnr\tMale\t-0.95\n", "missing column(s) Y_logratio"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\t\t-0.95\tnan\n", "is not one of Male, Female, Unknown"),
    ("sample\tX_logratio\tY_logratio\ns.cnr\t-0.95\tnan\n", "not a cnvkit.py sex table"),
    ("sample\tsex\tX_logratio\tY_logratio\n", "expected one sample row, found 0"),
    ("", "expected one sample row, found 0"),
    # round 6 (V1): CNVkit ends every record with a newline; an unterminated last line is a truncated table, even
    # when what is left still reads as a valid row
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\t-0.95\tnan", "not newline-terminated"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tFemale\t-0.9", "not newline-terminated"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\t-0.95\t-5.01", "not newline-terminated"),
    ("sample\tsex\tX_logratio\tY_logratio\ns.cnr\tMale\t-0.9", "not newline-terminated"),
])
@pytest.mark.parametrize("sex_arg", ["", "male", "female"])
def test_decide_cli_rejects_malformed_table(tmp_path, body, message, sex_arg):
    # a malformed .cnr sex table must fail, not read as "not Male" and silently skip the post-step
    sex = tmp_path / "sex_cnr.tsv"
    sex.write_text(body)
    proc = run("decide", "--sex-arg", sex_arg, "--cnr-sex", sex, "--parx", "grch38")
    assert proc.returncode == 1 and message in proc.stderr and proc.stdout == ""


@pytest.mark.parametrize("inferred", ["garbage", "", "male", "unknown", None])
def test_decide_rejects_invalid_inferred_sex(inferred):
    with pytest.raises(ValueError, match="is not one of Male, Female, Unknown"):
        port.decide("", inferred, "grch38")
    assert port.decide("male", inferred, "grch38")[0] is True     # -x given: the inference is not consulted


# CNVkit fc65941 writes this row when the .cnr has no usable chrX bins (no chrX targets, or only PAR ones with
# --diploid-parx-genome; checked with cnvkit.py sex in cgp-cnvkit:1.0.0): Unknown, as genemetrics' own guess (no
# shift), so the post-step is not applied. The base (8b683c8) accepted it; the expected output below is its output.
UNKNOWN = "sample\tsex\tX_logratio\tY_logratio\ns.cnr\tUnknown\tNA\tNA\n"


@pytest.mark.parametrize("sex_arg, stdout", [("", "0\tno -x; cnvkit.py sex on the .cnr: Unknown\n"),
                                             ("male", "1\t-x male\n"), ("female", "0\t-x female\n")])
def test_decide_accepts_cnvkits_unknown_sex_as_before(tmp_path, sex_arg, stdout):
    sex = tmp_path / "sex_cnr.tsv"
    sex.write_text(UNKNOWN)
    proc = run("decide", "--sex-arg", sex_arg, "--cnr-sex", sex, "--parx", "grch38")
    assert proc.returncode == 0 and proc.stdout == stdout, proc.stderr


@pytest.mark.parametrize("body", ["sample\tsex\tX_logratio\tY_logratio\ns.cnr\tFemale\tnan\tnan\n",
                                  "sample\tsex\tX_logratio\tY_logratio\ns.cnr\tUnknown\t-0.95\tnan\n",
                                  # columns in another order, or another column: not rejected (no exact header rule)
                                  "sample\tsex\tY_logratio\tX_logratio\ns.cnr\tMale\tnan\t-0.95\n",
                                  "sample\tsex\tX_logratio\tY_logratio\tother\ns.cnr\tMale\t-0.95\tnan\tx\n"])
def test_decide_accepts_well_formed_tables(tmp_path, body):
    sex = tmp_path / "sex_cnr.tsv"
    sex.write_text(body)
    proc = run("decide", "--sex-arg", "", "--cnr-sex", sex, "--parx", "grch38")
    assert proc.returncode == 0, proc.stderr


# ── compare ──────────────────────────────────────────────────────────────────

NONPAR = "TP53\tchr17\t7661779\t7687538\t-0.62\t400\t5\t1\t10\t5\t10\n"
CRLF2_ORIG = "CRLF2\tchrX\t1190891\t1212654\t1.26805\t679\t7.2\t3\t8\t7.2\t8\n"
# the post-step changes only log2: depth, weight, cn (the call.cns cn column, not shifted) and probes stay (on the
# 180 validated males, every corrected row differs from its original in log2 alone)
CRLF2_FIXED = "CRLF2\tchrX\t1190891\t1212654\t0.26805\t679\t7.2\t3\t8\t7.2\t8\n"


def gm(tmp_path, name, *rows, header=GM_H):
    p = tmp_path / name
    p.write_text(header + "".join(rows))
    return p


def call_cns(tmp_path, log2="0.26805", start=1190891, end=2781479):
    p = tmp_path / "s.call.cns"
    p.write_text(CNS_H + f"chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                         f"chrX\t{start}\t{end}\tCRLF2\t{log2}\t2\t679\t8\t7.2\n")
    return p


def test_compare_par_row_removed(tmp_path):
    # CRLF2 at 1.268 (+1) before; after the fix the segment is 0.268 < 0.3, so the row disappears
    rec = port.compare(gm(tmp_path, "o", NONPAR, CRLF2_ORIG), gm(tmp_path, "n", NONPAR), call_cns(tmp_path))
    assert rec["par_check"] == 1 and rec["rows_changed"] == 1 and rec["file_changed"] == 1


def test_compare_par_row_kept(tmp_path):
    rec = port.compare(gm(tmp_path, "o", NONPAR, CRLF2_ORIG.replace("1.26805", "1.5")),
                       gm(tmp_path, "n", NONPAR, CRLF2_FIXED.replace("0.26805", "0.5")), call_cns(tmp_path, "0.5"))
    assert rec["par_check"] == 1 and rec["par_rows_after"] == 1


def test_compare_non_par_difference_fails(tmp_path):
    proc = run("compare", gm(tmp_path, "o", NONPAR, CRLF2_ORIG), gm(tmp_path, "n", NONPAR.replace("-0.62", "-0.63")),
               call_cns(tmp_path))
    assert proc.returncode == 1 and "non-PAR rows differ" in proc.stderr


def test_compare_header_difference_fails(tmp_path):
    proc = run("compare", gm(tmp_path, "o", NONPAR), gm(tmp_path, "n", NONPAR, header=GM_H.replace("cn\t", "CN\t")),
               call_cns(tmp_path))
    assert proc.returncode == 1 and "FAILED: headers differ" in proc.stderr and proc.stdout == ""
    with pytest.raises(port.CheckFailed, match="headers differ"):
        port.compare(tmp_path / "o", tmp_path / "n", call_cns(tmp_path))


# round 4 (T2): the checks are explicit exceptions, so python -O (PYTHONOPTIMIZE) cannot remove them

def run_optimized(level, *argv):
    return python(SCRIPT, *argv, env={"PYTHONOPTIMIZE": str(level)})


@pytest.mark.parametrize("level", [1, 2])
def test_compare_checks_survive_python_optimize(tmp_path, level):
    # the reviewer's case: differing non-PAR rows used to give par_check=1 once the assert was optimized away
    o = gm(tmp_path, "o", NONPAR, CRLF2_ORIG)
    proc = run_optimized(level, "compare", o, gm(tmp_path, "n", NONPAR.replace("-0.62", "-0.63")), call_cns(tmp_path))
    assert proc.returncode == 1 and "FAILED: non-PAR rows differ" in proc.stderr and proc.stdout == ""
    proc = run_optimized(level, "compare", o, gm(tmp_path, "n2", NONPAR, CRLF2_ORIG.replace("1.26805", "1.3")),
                         call_cns(tmp_path))
    assert proc.returncode == 1 and "FAILED: par_check 0" in proc.stderr
    proc = run_optimized(level, "compare", o, gm(tmp_path, "n3", NONPAR, header=GM_H.replace("cn\t", "CN\t")),
                         call_cns(tmp_path))
    assert proc.returncode == 1 and "FAILED: headers differ" in proc.stderr
    proc = run_optimized(level, "compare", o, gm(tmp_path, "n4", NONPAR), call_cns(tmp_path))
    assert proc.returncode == 0 and proc.stdout.strip().endswith("par_check=1"), proc.stderr


@pytest.mark.parametrize("level", [1, 2])
def test_headers_differ_check_survives_python_optimize(tmp_path, level):
    o, n = gm(tmp_path, "o", NONPAR), gm(tmp_path, "n", NONPAR, header=GM_H.replace("cn\t", "CN\t"))
    code = (f"import importlib.util, sys; spec = importlib.util.spec_from_file_location('p', {str(SCRIPT)!r}); "
            "m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); "
            "m.compare(sys.argv[1], sys.argv[2], sys.argv[3])")
    proc = python("-c", code, o, n, call_cns(tmp_path), env={"PYTHONOPTIMIZE": str(level)})
    assert proc.returncode == 1 and "CheckFailed: headers differ" in proc.stderr


def test_compare_accepts_empty_table_header(tmp_path):
    rec = port.compare(gm(tmp_path, "o", CRLF2_ORIG), gm(tmp_path, "n", header=port.EMPTY_GM_HEADER),
                       call_cns(tmp_path))
    assert rec["empty_after"] == 1 and rec["par_check"] == 1


def test_compare_span(tmp_path):
    # the CRLF2 segment continues into non-PAR chrX: not shifted, genemetrics unchanged (segment + 1) -> span
    cns = call_cns(tmp_path, "0.26805", start=1190891, end=9000000)
    rec = port.compare(gm(tmp_path, "o", NONPAR, CRLF2_ORIG), gm(tmp_path, "n", NONPAR, CRLF2_ORIG), cns)
    assert rec["par_check"] == "span" and rec["file_changed"] == 0


def test_compare_span_names_the_gene_and_segment(tmp_path):
    cns = call_cns(tmp_path, "0.26805", start=1190891, end=9000000)
    rec = port.compare(gm(tmp_path, "o", NONPAR, CRLF2_ORIG), gm(tmp_path, "n", NONPAR, CRLF2_ORIG), cns)
    assert rec["span_genes"] == "CRLF2@chrX:1190891-9000000" and "unexplained" not in rec
    proc = run("compare", tmp_path / "o", tmp_path / "n", cns)
    assert proc.returncode == 0, proc.stderr
    assert "span_genes=CRLF2@chrX:1190891-9000000" in proc.stdout and proc.stdout.strip().endswith("par_check=span")


def test_compare_unexplained_par_row_fails(tmp_path):
    # no spanning segment; CRLF2 is still at +1 after the fix (1.268 != segment 0.268): par_check 0 stops the job
    o, n = gm(tmp_path, "o", NONPAR, CRLF2_ORIG), gm(tmp_path, "n", NONPAR, CRLF2_ORIG.replace("1.26805", "1.3"))
    assert port.compare(o, n, call_cns(tmp_path))["par_check"] == 0
    proc = run("compare", o, n, call_cns(tmp_path))
    assert proc.returncode == 1
    assert "FAILED: par_check 0" in proc.stderr and "par_check=0" in proc.stderr
    assert proc.stdout == ""


def test_compare_span_elsewhere_does_not_explain_a_par_row(tmp_path):
    # a segment spans the PAR2 boundary, but CRLF2's own segment lies wholly inside PAR1 and still carries the +1
    # after the fix: the reference calls this 'span' (any spanning segment + unchanged PAR rows); the port does not
    cns = tmp_path / "s.call.cns"
    cns.write_text(CNS_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                           "chrX\t1190891\t2781479\tCRLF2\t0.26805\t2\t679\t8\t7.2\n"
                           "chrX\t155000000\t155800000\tPAR2SPAN\t-0.1\t2\t500\t8\t7.2\n")
    o, n = gm(tmp_path, "o", NONPAR, CRLF2_ORIG), gm(tmp_path, "n", NONPAR, CRLF2_ORIG)
    rec = port.compare(o, n, cns)
    assert rec["span_segs"] == 1 and rec["par_check"] == 0
    assert rec["unexplained"] == "CRLF2@chrX:1190891-1212654"
    (tmp_path / "s.cnr").write_text(CNR_H)
    assert ref.compare(dict(action="dry", rundir=str(tmp_path), sample="s", orig=str(o), new=str(n)))["par_check"] \
        == "span"
    proc = run("compare", o, n, cns)
    assert proc.returncode == 1 and "CRLF2@chrX:1190891-1212654" in proc.stderr


# per-row verdicts (R1, R2): each PAR row is judged against its own segment

PAR2_ORIG = "VAMP7\tchrX\t155881000\t155943000\t1.5\t500\t6\t3\t9\t7.2\t8\n"
PAR2_FIXED = "VAMP7\tchrX\t155881000\t155943000\t0.5\t500\t6\t3\t9\t7.2\t8\n"


def cns_rows(tmp_path, *rows):
    """rows: (start, end, log2[, probes, weight]); probes 8 and weight 7.2 by default, as the rows above carry."""
    p = tmp_path / "s.call.cns"
    p.write_text(CNS_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n" + "".join(
        f"chrX\t{s}\t{e}\tG\t{l2}\t2\t500\t{n}\t{w}\n" for s, e, l2, n, w in (r + (8, 7.2)[len(r) - 3:] for r in rows)))
    return p


def verdict(tmp_path, before, after, cns):
    o, n = gm(tmp_path, "o", NONPAR, *before), gm(tmp_path, "n", NONPAR, *after)
    return port.compare(o, n, cns), run("compare", o, n, cns)


@pytest.mark.parametrize("log2", ["99", "-0.5", "0.26805", "1.3", "1.2681"])
def test_compare_span_row_must_carry_the_male_shift(tmp_path, log2):
    # R1: an unchanged row in a spanning segment is explained only by segment log2 + 1.0 (here 1.26805)
    row = CRLF2_ORIG.replace("1.26805", log2)
    rec, proc = verdict(tmp_path, [row], [row], cns_rows(tmp_path, (1190891, 9000000, "0.26805")))
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1212654"
    assert proc.returncode == 1 and "FAILED: par_check 0" in proc.stderr


def test_compare_span_row_with_large_value_is_not_excused(tmp_path):
    # R1's example: a row at log2 99 in a spanning segment at -0.5 (the +1 would give 0.5)
    row = CRLF2_ORIG.replace("1.26805", "99")
    rec, proc = verdict(tmp_path, [row], [row], cns_rows(tmp_path, (1190891, 9000000, "-0.5")))
    assert rec["par_check"] == 0 and proc.returncode == 1


@pytest.mark.parametrize("log2, check", [
    ("1.26805", "span"), ("1.268054", "span"), ("1.268046", "span"),     # within half a unit in the 6th digit
    ("1.268056", 0), ("1.268044", 0), ("1.26806", 0), ("1.26804", 0),    # beyond it (1e-5 used to pass)
])
def test_compare_span_tolerance(tmp_path, log2, check):
    # the row is CNVkit's %.6g print of segment + 1.0 (1.26805): it may differ by half a unit in the 6th digit
    row = CRLF2_ORIG.replace("1.26805", log2)
    rec, _ = verdict(tmp_path, [row], [row], cns_rows(tmp_path, (1190891, 9000000, "0.26805")))
    assert rec["par_check"] == check


# round 4 (T3): the tolerance scales with the printed precision and magnitude, as CNVkit prints with %.6g

CRLF2_ROW = "CRLF2\tchrX\t1190891\t1212654\t{}\t679\t7.2\t3\t8\t7.2\t8\n"
WHOLLY, SPANNING = (1190891, 2781479), (1190891, 9000000)


def cnvkit_rows(L, segment):
    """The genemetrics PAR rows CNVkit writes before and after the fix for a CRLF2 segment at log2 L (a %.6g string),
    reproducing its arithmetic: before, L + 1.0; after, (L - 1.0) + 1.0 when wholly in a PAR (else unchanged); a gene
    is reported when |level| >= 0.3; every value printed with %.6g."""
    x = float(L)
    before = [CRLF2_ROW.format("%.6g" % (x + 1.0))] if abs(x + 1.0) >= 0.3 else []
    if segment == SPANNING:
        return before, before
    after = (x - 1.0) + 1.0
    return before, [CRLF2_ROW.format("%.6g" % after)] if abs(after) >= 0.3 else []


def par_check_of(tmp_path, L, segment, before=None, after=None):
    b, a = cnvkit_rows(L, segment)
    before, after = b if before is None else before, a if after is None else after
    rec = port.compare(gm(tmp_path, "o", NONPAR, *before), gm(tmp_path, "n", NONPAR, *after),
                       cns_rows(tmp_path, (*segment, L)))
    return rec["par_check"]


BOUNDARY = ["9.99996", "99.9996", "999.996", "9999.96", "0.999996", "1.99999", "8.99999", "-1.99999", "-2.00001",
            "-9.99996", "-10.9999", "0.3", "-0.3", "0.300001", "-0.300001", "0.299999", "-0.699999", "-0.700001",
            "-1.29999", "-1.30001", "1e-05", "-0.999996", "4.99999", "0.0499995"]


@pytest.mark.parametrize("L", BOUNDARY)
@pytest.mark.parametrize("segment", [WHOLLY, SPANNING], ids=["wholly", "spanning"])
def test_compare_rounding_boundaries(tmp_path, L, segment):
    # the reviewer's case is L = 9.99996: the valid pre-fix value 10.99996 prints as 11 (1e-5 rejected it)
    assert par_check_of(tmp_path, L, segment) == (1 if segment == WHOLLY else "span")


def test_compare_reviewers_case(tmp_path):
    assert cnvkit_rows("9.99996", WHOLLY) == ([CRLF2_ROW.format("11")], [CRLF2_ROW.format("9.99996")])
    assert par_check_of(tmp_path, "9.99996", WHOLLY) == 1
    assert par_check_of(tmp_path, "9.99996", SPANNING) == "span"
    proc = run("compare", tmp_path / "o", tmp_path / "n", tmp_path / "s.call.cns")
    assert proc.returncode == 0 and proc.stdout.strip().endswith("par_check=span"), proc.stderr


BEYOND = [
    ("9.99996", "11.0001", "9.99996"),       # a unit off at magnitude 10 (half-unit 5e-5)
    ("9.99996", "11", "9.99986"),
    ("99.9996", "101.001", "99.9996"),
    ("0.5", "1.50001", "0.5"),               # a unit off at magnitude 1 (half-unit 5e-6)
    ("0.5", "1.5", "0.500001"),              # a unit off at magnitude 0.1 (half-unit 5e-7)
    ("-0.699999", "0.300002", "-0.699999"),
]


@pytest.mark.parametrize("L, before, after, segment", [
    *[(*c, WHOLLY) for c in BEYOND],
    # a spanning segment's row is not shifted (before = after): only the cases with a wrong before-value
    *[(*c, SPANNING) for c in BEYOND if [CRLF2_ROW.format(c[1])] != cnvkit_rows(c[0], SPANNING)[0]],
])
def test_compare_beyond_the_printed_precision_fails(tmp_path, L, before, after, segment):
    b, a = [CRLF2_ROW.format(before)], [CRLF2_ROW.format(after if segment == WHOLLY else before)]
    assert par_check_of(tmp_path, L, segment, b, a) == 0


@pytest.mark.parametrize("seed", range(5))
def test_compare_explains_every_cnvkit_print(tmp_path, seed):
    rng = random.Random(seed)   # noqa: S311 -- reproducible test values, not security
    for _ in range(80):
        L = "%.6g" % (rng.choice([-1, 1]) * rng.uniform(0.5, 9.99999) * 10.0 ** rng.randint(-3, 3))
        for segment in (WHOLLY, SPANNING):
            assert par_check_of(tmp_path, L, segment) == (1 if segment == WHOLLY else "span"), (L, segment)


def test_printed_tol():
    assert port.printed_tol(11, 10.99996) == pytest.approx(5e-5, abs=2e-12)
    assert port.printed_tol(1.26805) == pytest.approx(5e-6, abs=2e-12)
    assert port.printed_tol(0.5) == pytest.approx(5e-7, abs=2e-12)
    assert port.printed_tol(0.0) == port.FLOAT_SLACK


def test_compare_span_row_that_changed_is_not_excused(tmp_path):
    # a spanning segment is not shifted, so its gene's row must not change
    before = CRLF2_ORIG.replace("679", "680")
    rec, proc = verdict(tmp_path, [before], [CRLF2_ORIG], cns_rows(tmp_path, (1190891, 9000000, "0.26805")))
    assert rec["par_check"] == 0 and proc.returncode == 1


def test_compare_mixed_span_and_corrected(tmp_path):
    # R2: CRLF2 in a segment spanning the PAR1 boundary (unchanged, +1) and VAMP7 in a wholly-PAR2 segment (corrected)
    cns = cns_rows(tmp_path, (1190891, 9000000, "0.26805"), (155800000, 156000000, "0.5"))
    rec, proc = verdict(tmp_path, [CRLF2_ORIG, PAR2_ORIG], [CRLF2_ORIG, PAR2_FIXED], cns)
    assert rec["par_check"] == "span" and rec["span_genes"] == "CRLF2@chrX:1190891-9000000"
    assert "unexplained" not in rec
    assert proc.returncode == 0 and proc.stdout.strip().endswith("span_genes=CRLF2@chrX:1190891-9000000 par_check=span")
    # the reference requires all PAR rows unchanged for span, so it fails this sample (par_check 0)
    (tmp_path / "s.cnr").write_text(CNR_H)
    assert ref.compare(dict(action="dry", rundir=str(tmp_path), sample="s", orig=str(tmp_path / "o"),
                            new=str(tmp_path / "n")))["par_check"] == 0


def test_compare_mixed_span_and_dropped(tmp_path):
    # the corrected PAR2 gene drops below 0.3 once the +1 is gone; CRLF2 stays spanned
    cns = cns_rows(tmp_path, (1190891, 9000000, "0.26805"), (155800000, 156000000, "0.2"))
    rec, proc = verdict(tmp_path, [CRLF2_ORIG, PAR2_ORIG.replace("1.5", "1.2")], [CRLF2_ORIG], cns)
    assert rec["par_check"] == "span" and proc.returncode == 0


def test_compare_mixed_with_a_bad_corrected_row_fails(tmp_path):
    cns = cns_rows(tmp_path, (1190891, 9000000, "0.26805"), (155800000, 156000000, "0.5"))
    rec, proc = verdict(tmp_path, [CRLF2_ORIG, PAR2_ORIG], [CRLF2_ORIG, PAR2_ORIG], cns)   # VAMP7 kept its +1
    assert rec["par_check"] == 0 and rec["unexplained"] == "VAMP7@chrX:155881000-155943000"
    assert proc.returncode == 1


def test_compare_row_dropped_from_a_spanning_segment_fails(tmp_path):
    # an unshifted segment cannot lose its gene
    rec, proc = verdict(tmp_path, [CRLF2_ORIG], [], cns_rows(tmp_path, (1190891, 9000000, "0.26805")))
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1212654 (dropped)"
    assert proc.returncode == 1


# round 3 (S1): a row in a spanning segment is never 'corrected', whatever its value

@pytest.mark.parametrize("before, after", [
    ("-0.5", "-0.5"),          # the reviewer's case: unchanged, equal to the unshifted segment
    ("0.5", "-0.5"),           # changed to the unshifted segment log2, as a wholly-PAR row would be
    ("-0.5", "0.5"),
])
def test_compare_span_row_equal_to_the_segment_is_not_corrected(tmp_path, before, after):
    b, a = CRLF2_ORIG.replace("1.26805", before), CRLF2_ORIG.replace("1.26805", after)
    rec, proc = verdict(tmp_path, [b], [a], cns_rows(tmp_path, (1190891, 9000000, "-0.5")))
    assert rec["par_check"] == 0 and rec["unexplained"].startswith("CRLF2@chrX:1190891-1212654")
    assert proc.returncode == 1 and "FAILED: par_check 0" in proc.stderr


def test_compare_span_row_equal_to_the_segment_without_a_spanned_label(tmp_path):
    # R1's bypass: with the row 'corrected' and no row spanned, par_check was span with span_genes=- (now 0)
    row = CRLF2_ORIG.replace("1.26805", "-0.5")
    rec, _ = verdict(tmp_path, [row], [row], cns_rows(tmp_path, (1190891, 9000000, "-0.5")))
    assert rec["par_check"] == 0 and "span_genes" not in rec


# round 3 (S2): a before-row is replaced only by the after-row with its gene, coordinates and other fields

CRLF2_A = "CRLF2\tchrX\t1190891\t1200000\t{}\t679\t7.2\t3\t4\t7.2\t8\n"      # CRLF2 split over two segments
CRLF2_B = "CRLF2\tchrX\t1200000\t1212654\t{}\t679\t7.2\t3\t4\t7.2\t8\n"
TWO_SEGS = ((1190891, 1200000, "0.5"), (1200000, 2781479, "0.6"))


def test_compare_same_gene_in_two_segments_corrected(tmp_path):
    rec, proc = verdict(tmp_path, [CRLF2_A.format("1.5"), CRLF2_B.format("1.6")],
                        [CRLF2_A.format("0.5"), CRLF2_B.format("0.6")], cns_rows(tmp_path, *TWO_SEGS))
    assert rec["par_check"] == 1 and proc.returncode == 0


def test_compare_same_gene_one_segment_dropped_below_threshold(tmp_path):
    segs = ((1190891, 1200000, "0.2"), (1200000, 2781479, "0.6"))
    rec, _ = verdict(tmp_path, [CRLF2_A.format("1.2"), CRLF2_B.format("1.6")], [CRLF2_B.format("0.6")],
                     cns_rows(tmp_path, *segs))
    assert rec["par_check"] == 1


def test_compare_missing_row_not_excused_by_the_same_gene_in_another_segment(tmp_path):
    # the reviewer's case: CRLF2_A (|0.5| >= 0.3, so genemetrics reports it after the fix) is missing; CRLF2_B in the
    # other segment no longer excuses it
    rec, proc = verdict(tmp_path, [CRLF2_A.format("1.5"), CRLF2_B.format("1.6")], [CRLF2_B.format("0.6")],
                        cns_rows(tmp_path, *TWO_SEGS))
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1200000 (dropped)"
    assert proc.returncode == 1 and "CRLF2@chrX:1190891-1200000 (dropped)" in proc.stderr
    (tmp_path / "s.cnr").write_text(CNR_H)        # the reference accepts it: any same-name after-row
    assert ref.compare(dict(action="dry", rundir=str(tmp_path), sample="s", orig=str(tmp_path / "o"),
                            new=str(tmp_path / "n")))["par_check"] == 1


def test_compare_missing_row_not_excused_by_the_same_gene_in_par2(tmp_path):
    other = "CRLF2\tchrX\t155881000\t155943000\t0.5\t500\t6\t3\t9\t7.2\t8\n"
    cns = cns_rows(tmp_path, (1190891, 2781479, "0.5"), (155800000, 156000000, "0.5"))
    rec, _ = verdict(tmp_path, [CRLF2_ORIG.replace("1.26805", "1.5")], [other], cns)
    assert rec["par_check"] == 0
    assert rec["unexplained"] == ("CRLF2@chrX:1190891-1212654 (dropped);"
                                   "CRLF2@chrX:155881000-155943000 (not in the original)")


def test_compare_replacement_must_keep_every_field_but_log2(tmp_path):
    after = CRLF2_FIXED.replace("0.26805", "0.5").replace("\t679\t", "\t680\t")          # depth changed
    rec, _ = verdict(tmp_path, [CRLF2_ORIG.replace("1.26805", "1.5")], [after], call_cns(tmp_path, "0.5"))
    assert rec["par_check"] == 0
    assert rec["unexplained"] == "CRLF2@chrX:1190891-1212654 (dropped);CRLF2@chrX:1190891-1212654 (not in the original)"


def test_compare_replacement_multiplicity(tmp_path):
    # two identical before-rows, one corrected after-row: one before-row is unexplained
    b = CRLF2_ORIG.replace("1.26805", "1.5")
    rec, _ = verdict(tmp_path, [b, b], [CRLF2_FIXED.replace("0.26805", "0.5")], call_cns(tmp_path, "0.5"))
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1212654 (before the fix)"
    rec, _ = verdict(tmp_path, [b], [CRLF2_FIXED.replace("0.26805", "0.5")] * 2, call_cns(tmp_path, "0.5"))
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1212654 (not in the original)"


def test_compare_spanned_row_multiplicity(tmp_path):
    # one unchanged before-row explains one spanned after-row, not two
    rec, _ = verdict(tmp_path, [CRLF2_ORIG], [CRLF2_ORIG, CRLF2_ORIG],
                     cns_rows(tmp_path, (1190891, 9000000, "0.26805")))
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1212654"


def test_compare_original_row_must_carry_the_male_shift(tmp_path):
    # a wholly-PAR row equal to the segment both before and after the fix: the original was never shifted
    row = CRLF2_FIXED.replace("0.26805", "0.5")
    rec, proc = verdict(tmp_path, [row], [row], call_cns(tmp_path, "0.5"))
    assert rec["par_check"] == 0 and rec["unexplained"].startswith("CRLF2@chrX:1190891-1212654")
    assert proc.returncode == 1
    # dropped below the threshold, but the original was not at segment + 1
    rec, _ = verdict(tmp_path, [CRLF2_ORIG.replace("1.26805", "0.9")], [], call_cns(tmp_path))
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1212654 (dropped)"


# round 5 (U1): a before-row must be at a level genemetrics reports (|L + 1.0| >= 0.3), whatever the segment kind

@pytest.mark.parametrize("L, before", [("-1.2", "-0.2"), ("-1", "1e-12"), ("-0.9", "0.1"), ("-1.29", "-0.29")])
def test_compare_before_row_below_the_threshold_is_unexplained(tmp_path, L, before):
    # the reviewer's case (L = -1.2): an original row at -0.2 and a corrected row at -1.2 used to pass (par_check 1);
    # genemetrics -t 0.3 cannot have reported the original
    b, a = CRLF2_ROW.format(before), CRLF2_ROW.format(L)
    rec, proc = verdict(tmp_path, [b], [a], cns_rows(tmp_path, (*WHOLLY, L)))
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1212654 (before the fix)"
    assert proc.returncode == 1 and "FAILED: par_check 0" in proc.stderr
    # the same corrected row with no original is what CNVkit writes
    rec, _ = verdict(tmp_path, [], [a], cns_rows(tmp_path, (*WHOLLY, L)))
    assert rec["par_check"] == 1


@pytest.mark.parametrize("L, row", [("-1.2", "-0.2"), ("-1", "1e-12"), ("-1.29", "-0.29")])
def test_compare_unchanged_span_row_below_the_threshold_is_unexplained(tmp_path, L, row):
    # the reviewer's second case: an unchanged row at -0.2 in a spanning segment at -1.2 used to pass (span)
    r = CRLF2_ROW.format(row)
    rec, proc = verdict(tmp_path, [r], [r], cns_rows(tmp_path, (*SPANNING, L)))
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1212654" and "span_genes" not in rec
    assert proc.returncode == 1 and "FAILED: par_check 0" in proc.stderr
    rec, _ = verdict(tmp_path, [], [], cns_rows(tmp_path, (*SPANNING, L)))     # what CNVkit writes: no row
    assert rec["par_check"] == "span"


@pytest.mark.parametrize("segment", [WHOLLY, SPANNING], ids=["wholly", "spanning"])
@pytest.mark.parametrize("L", ["-0.7", "-1.3", "-0.699999", "-1.30001"])
def test_compare_before_row_at_the_threshold_is_accepted(tmp_path, segment, L):
    # |L + 1.0| at (or just past) 0.3: reported before the fix, as CNVkit does
    assert cnvkit_rows(L, segment)[0]
    assert par_check_of(tmp_path, L, segment) == (1 if segment == WHOLLY else "span")


def test_compare_new_row_after_the_fix(tmp_path):
    # segment -0.8: before the fix -0.8 + 1 = 0.2 (< 0.3, not reported); after it -0.8 is reported
    after = CRLF2_FIXED.replace("0.26805", "-0.8")
    rec, proc = verdict(tmp_path, [], [after], call_cns(tmp_path, "-0.8"))
    assert rec["par_check"] == 1 and proc.returncode == 0
    # segment 0.5: the original (at 1.5) should have reported it
    rec, proc = verdict(tmp_path, [], [CRLF2_FIXED.replace("0.26805", "0.5")], call_cns(tmp_path, "0.5"))
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1212654 (not in the original)"
    assert proc.returncode == 1


# Z1 (round 9): overlapping segments. CNVkit writes overlapping segments for overlapping targets, and genemetrics
# reports a gene once per segment whose range overlaps its bins (by_ranges, mode 'outer'), at that segment's log2;
# a row is explained when any segment it overlaps explains it.

BASE_STDOUT = ("par_segs=2 span_segs=0 empty_after=0 nonpar_identical=1 par_rows_before=1 par_rows_after=1 "
               "rows_changed=2 file_changed=1 par_before=CRLF2:1.5 par_after=CRLF2:0.5 par_check=1\n")


@pytest.mark.parametrize("segments", [((1100000, 1300000, "0.5"), (1190891, 2781479, "0.5")),
                                      ((1190891, 2781479, "0.5"), (1100000, 1300000, "0.5"))], ids=["fixture", "swapped"])
def test_compare_gene_in_two_overlapping_segments_matches_the_base(tmp_path, segments):
    # the committed round-9 fixture: CRLF2 lies in two overlapping wholly-PAR segments, both at 0.5, and both
    # explain the corrected row. The base (8b683c8), whose compare is the reference's, returns par_check 1; so must
    # the port, with the same statistics and the same stdout byte for byte
    cns = cns_rows(tmp_path, *segments)
    o = gm(tmp_path, "s.genemetrics.noparfix.tsv", NONPAR, CRLF2_ORIG.replace("1.26805", "1.5"))
    n = gm(tmp_path, "s.genemetrics.tsv", NONPAR, CRLF2_FIXED.replace("0.26805", "0.5"))
    (tmp_path / "s.cnr").write_text(CNR_H)
    mine = port.compare(o, n, cns)
    theirs = ref.compare(dict(action="dry", rundir=str(tmp_path), sample="s", orig=str(o), new=str(n)))
    assert mine["par_check"] == theirs["par_check"] == 1
    assert {k: theirs[k] for k in mine} == mine
    proc = run("compare", o, n, cns)
    assert proc.returncode == 0 and proc.stdout == BASE_STDOUT and proc.stderr == ""


# genemetrics output of cgp-cnvkit:1.0.0 (cnvkit.py genemetrics -t 0.3 --min-probes 3 --diploid-parx-genome grch38
# -x male), before and after `shift`, on a .cnr whose CRLF2 bins 1191000-1191300 and 1191200-1191500 overlap and a
# call.cns with two overlapping segments meeting there: the reported row 1191000-1212654 starts before its own
# segment (1191200-1212654), so no segment contains it. The base reports par_check 0 here (a warning); the port,
# as par_check 0 now fails the job, must explain the row.
CRLF2_SPLIT = "CRLF2\tchrX\t1191000\t1212654\t{}\t600\t4.5\t2\t5\t3.6\t4\n"


def split_segments(tmp_path, second):
    p = tmp_path / "s.call.cns"
    p.write_text(CNS_H + "chr17\t7000000\t7009200\tTP53\t-0.6\t1\t400\t10\t9\n"
                         "chrX\t1190891\t1191300\tCRLF2\t0.5\t2\t600\t2\t1.8\n"
                         f"chrX\t1191200\t1212654\tCRLF2\t{second}\t2\t600\t4\t3.6\n"
                         "chrX\t5000000\t5004200\tNONPAR\t0.0\t2\t500\t5\t4.5\n")
    return p


@pytest.mark.parametrize("second, before, after", [
    ("0.8", ["1.8"], ["0.8"]),      # the second segment: corrected
    ("0.5", ["1.5"], ["0.5"]),      # both segments at 0.5
    ("-0.8", [], ["-0.8"]),         # not reported before the fix (0.2), reported after it
])
def test_compare_row_extending_past_its_segment_is_explained(tmp_path, second, before, after):
    cns = split_segments(tmp_path, second)
    rec, proc = verdict(tmp_path, [CRLF2_SPLIT.format(v) for v in before], [CRLF2_SPLIT.format(v) for v in after], cns)
    assert rec["par_check"] == 1 and "unexplained" not in rec
    assert proc.returncode == 0 and proc.stdout.strip().endswith("par_check=1"), proc.stderr


@pytest.mark.parametrize("before, after", [
    (["1.8"], ["0.7"]),             # neither overlapping segment (0.5, 0.8) explains 0.7
    (["1.8"], ["1.8"]),             # still +1 after the fix
    (["1.8"], []),                  # dropped, although 0.8 is reported
    ([], ["0.8"]),                  # new, although 1.8 was reported before the fix
])
def test_compare_overlapping_segments_still_reject_an_unexplained_row(tmp_path, before, after):
    cns = split_segments(tmp_path, "0.8")
    rec, proc = verdict(tmp_path, [CRLF2_SPLIT.format(v) for v in before], [CRLF2_SPLIT.format(v) for v in after], cns)
    assert rec["par_check"] == 0 and rec["unexplained"].startswith("CRLF2@chrX:1191000-1212654")
    assert proc.returncode == 1 and "FAILED: par_check 0" in proc.stderr


def test_compare_two_rows_of_one_gene_from_two_overlapping_segments(tmp_path):
    # one row per overlapping segment, each at its own segment's log2 (all of them would reject this valid table)
    cns = cns_rows(tmp_path, (1100000, 1300000, "0.5"), (1190891, 2781479, "0.8", 9, 9))
    rows = [CRLF2_ORIG.replace("1.26805", "1.5"), CRLF2_ORIG.replace("1.26805", "1.8").replace("\t7.2\t8\n", "\t9\t9\n")]
    rec, proc = verdict(tmp_path, rows, [r.replace("\t1.5\t", "\t0.5\t").replace("\t1.8\t", "\t0.8\t") for r in rows], cns)
    assert rec["par_check"] == 1 and proc.returncode == 0, proc.stderr


# F1 (round 2): genemetrics copies a segment's probes and weight into each row it reports from it (segment_probes,
# segment_weight), and --min-probes 3 drops the rows of a segment with fewer probes: an overlapping segment whose
# fields the row does not carry cannot explain it

@pytest.mark.parametrize("probes, weight", [(1, 7.2), (8, 1), (1, 1)])
def test_compare_overlapping_segment_without_the_rows_fields_does_not_explain_it(tmp_path, probes, weight):
    # CRLF2 (segment_probes 8) comes from its wholly-PAR segment and still carries the +1 after the fix; an
    # overlapping segment spanning the PAR1 boundary at the same log2, with 1 probe (filtered out) or another weight
    cns = cns_rows(tmp_path, WHOLLY + ("0.26805",), (1200000, 9000000, "0.26805", probes, weight))
    rec, proc = verdict(tmp_path, [CRLF2_ORIG], [CRLF2_ORIG], cns)
    assert rec["par_check"] == 0 and rec["unexplained"] == "CRLF2@chrX:1190891-1212654" and "span_genes" not in rec
    assert proc.returncode == 1 and "FAILED: par_check 0" in proc.stderr
    # the same segment with the row's probes and weight explains it (span), at their %.6g print
    cns = cns_rows(tmp_path, WHOLLY + ("0.26805",), (1200000, 9000000, "0.26805", 8, "7.200001"))
    assert verdict(tmp_path, [CRLF2_ORIG], [CRLF2_ORIG], cns)[0]["par_check"] == "span"


def test_compare_corrected_row_needs_its_segments_fields(tmp_path):
    # a corrected row is judged against the segment it carries the fields of: the overlapping 1-probe segment at
    # 0.5 does not explain CRLF2 at 0.5 when CRLF2's own segment (8 probes) is at 0.8
    cns = cns_rows(tmp_path, (1100000, 1300000, "0.5", 1, 1), WHOLLY + ("0.8",))
    rec, _ = verdict(tmp_path, [CRLF2_ROW.format("1.5")], [CRLF2_ROW.format("0.5")], cns)
    assert rec["par_check"] == 0
    rec, _ = verdict(tmp_path, [CRLF2_ROW.format("1.8")], [CRLF2_ROW.format("0.8")], cns)
    assert rec["par_check"] == 1


def test_row_segments():
    segs = {"chrX": [(1100000, 1300000, 0.5), (1190891, 2781479, 0.5), (2781479, 9000000, 0.1)]}
    assert port.row_segments(segs, 1190891, 1212654) == segs["chrX"][:2]
    assert port.row_segments(segs, 2781000, 2782000) == segs["chrX"][1:]      # a row across two adjacent segments
    assert port.row_segments(segs, 1000, 2000) == []
    # half-open: a segment ending where the row starts does not overlap it
    assert port.row_segments({"chrX": [(1000, 1190891, 0.5)]}, 1190891, 1212654) == []
    # every segment the base's lookup finds (containment) is a candidate, a zero-length row at a segment start too
    zero = {"chrX": [(1190891, 1212654, 0.5), (1000, 1190891, 0.4)]}
    assert port.row_segments(zero, 1190891, 1190891) == zero["chrX"]
    assert [x for x in zero["chrX"] if x[0] <= 1190891 and 1190891 <= x[1]] == zero["chrX"]


def test_compare_zero_length_par_row_as_the_base(tmp_path):
    # a zero-length row at its segment's start: the base (containment) explains it, so must the port
    cns = cns_rows(tmp_path, (1190891, 2781479, "0.5"))
    row = "CRLF2\tchrX\t1190891\t1190891\t{}\t679\t7.2\t3\t8\t7.2\t8\n"
    o = gm(tmp_path, "s.genemetrics.noparfix.tsv", NONPAR, row.format("1.5"))
    n = gm(tmp_path, "s.genemetrics.tsv", NONPAR, row.format("0.5"))
    (tmp_path / "s.cnr").write_text(CNR_H)
    mine = port.compare(o, n, cns)
    theirs = ref.compare(dict(action="dry", rundir=str(tmp_path), sample="s", orig=str(o), new=str(n)))
    assert mine["par_check"] == theirs["par_check"] == 1 and {k: theirs[k] for k in mine} == mine


# the minimum the comparison needs to read a genemetrics table (round 3, S3): a header naming the columns it reads,
# one field per column, integer coordinates and a numeric log2; an empty or truncated table fails

GM_TRUNCATED = "CRLF2\tchrX\t1190891\t1212654\t0.2"


@pytest.mark.parametrize("orig, new, message", [
    ("", "", "empty file, no header"),
    (GM_H + NONPAR, "", "empty file, no header"),
    ("", GM_H + NONPAR, "empty file, no header"),
    ("garbage\n", "garbage\n", "not a genemetrics header"),
    ("gene\tchromosome\tstart\tend\n", "gene\tchromosome\tstart\tend\n", "not a genemetrics header"),
    ("chromosome\tgene\tstart\tend\tlog2\n", "chromosome\tgene\tstart\tend\tlog2\n", "not a genemetrics header"),
    (GM_H + NONPAR, GM_H + NONPAR + GM_TRUNCATED, "line 3 is not a 11-field genemetrics row"),     # truncated row
    (GM_H + NONPAR, GM_H + NONPAR.replace("\t5\t10\n", "\n"), "line 2 is not a 11-field genemetrics row"),
    (GM_H + NONPAR, GM_H + NONPAR + "\n", "line 3 is not a 11-field genemetrics row"),
    (GM_H + NONPAR, GM_H + NONPAR.replace("-0.62", "NA"), "line 2 is not"),
    (GM_H + NONPAR, GM_H + NONPAR.replace("-0.62", ""), "line 2 is not"),
    (GM_H + NONPAR, GM_H + NONPAR.replace("7661779", "x"), "line 2 is not"),
    (GM_H + NONPAR, GM_H + NONPAR.replace("7687538", "7.5e6"), "line 2 is not"),
    (GM_H + NONPAR + "CRLF2\tchrX\t1190891\t1212654\tgarbage\t679\t7.2\t3\t8\t7.2\t8\n", GM_H + NONPAR,
     "line 3 is not"),
])
def test_compare_rejects_an_unreadable_genemetrics(tmp_path, orig, new, message):
    (tmp_path / "o").write_text(orig)
    (tmp_path / "n").write_text(new)
    with pytest.raises(ValueError, match=re.escape(message)):
        port.compare(tmp_path / "o", tmp_path / "n", call_cns(tmp_path))
    proc = run("compare", tmp_path / "o", tmp_path / "n", call_cns(tmp_path))
    assert proc.returncode == 1 and "cnvkit_parfix.py compare: FAILED" in proc.stderr and message in proc.stderr
    assert proc.stdout == ""


def test_compare_does_not_check_the_genemetrics_schema(tmp_path):
    # identical headers with other columns after log2 are compared as before: no exact-schema rule
    other = "gene\tchromosome\tstart\tend\tlog2\tdepth\textra\n"
    rec = port.compare(gm(tmp_path, "o", "TP53\tchr17\t7661779\t7687538\t-0.62\t400\tx\n", header=other),
                       gm(tmp_path, "n", "TP53\tchr17\t7661779\t7687538\t-0.62\t400\tx\n", header=other),
                       call_cns(tmp_path))
    assert rec["par_check"] == 1


def test_compare_accepts_header_only_tables(tmp_path):
    # CNVkit's empty table (no gene at all) before and after the fix
    rec = port.compare(gm(tmp_path, "o", header=port.EMPTY_GM_HEADER), gm(tmp_path, "n", header=port.EMPTY_GM_HEADER),
                       call_cns(tmp_path))
    assert rec["par_check"] == 1 and rec["empty_after"] == 1
    rec = port.compare(gm(tmp_path, "o"), gm(tmp_path, "n"), call_cns(tmp_path))      # full header, no rows
    assert rec["par_check"] == 1


def test_compare_accepts_a_call_cns_nan_segment(tmp_path):
    cns = tmp_path / "s.call.cns"
    cns.write_text(CNS_H + "chr1\t100\t200\tG\tnan\t2\t0\t1\t0\n"
                           "chrX\t1190891\t2781479\tCRLF2\t0.26805\t2\t679\t8\t7.2\n")
    assert port.compare(gm(tmp_path, "o", NONPAR, CRLF2_ORIG), gm(tmp_path, "n", NONPAR), cns)["par_check"] == 1


@pytest.mark.parametrize("log2", ["inf", "-inf", "nan"])
def test_compare_non_finite_segment_explains_no_row(tmp_path, log2):
    rec, proc = verdict(tmp_path, [CRLF2_ORIG], [CRLF2_FIXED], call_cns(tmp_path, log2))
    assert rec["par_check"] == 0 and proc.returncode == 1 and "Traceback" not in proc.stderr


def test_compare_accepts_the_empty_header_before_the_fix(tmp_path):
    # no gene reported before the fix (5-column header), a PAR gene reported after it (segment -0.8: at 0.2 before)
    rec = port.compare(gm(tmp_path, "o", header=port.EMPTY_GM_HEADER),
                       gm(tmp_path, "n", CRLF2_FIXED.replace("0.26805", "-0.8")), call_cns(tmp_path, "-0.8"))
    assert rec["par_check"] == 1 and rec["file_changed"] == 1


def test_compare_cli_output(tmp_path):
    proc = run("compare", gm(tmp_path, "o", NONPAR, CRLF2_ORIG), gm(tmp_path, "n", NONPAR), call_cns(tmp_path))
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip().endswith("par_check=1")


# ── differential against the reference ──────────────────────────────────────

def random_rows(rng, n):
    rows, pos = [], 0
    for _ in range(n):
        chrom = rng.choice(["chrX", "chrX", "chrX", "chr7", "chrY"])
        pos = rng.choice([rng.randint(0, 3_000_000), rng.randint(155_600_000, 156_100_000), rng.randint(0, 160_000_000)])
        end = pos + rng.choice([1, 120, 144, 5000, 3_000_000])
        rows.append(f"{chrom}\t{pos}\t{end}\tG\t{rng.uniform(-2, 2)!r}\t100\t1\n")
    return rows


@pytest.mark.parametrize("seed", range(10))
def test_shift_byte_identical_to_reference(tmp_path, seed):
    src = tmp_path / "in.cnr"
    src.write_text(CNR_H + "".join(random_rows(random.Random(seed), 400)))
    a, b = tmp_path / "port.cnr", tmp_path / "ref.cnr"
    assert port.shift_par(src, a, 1, 4) == ref.shift_par(src, b, 1, 4)
    assert a.read_bytes() == b.read_bytes()


@pytest.mark.parametrize("case", ["removed", "kept", "span", "empty"])
def test_compare_matches_reference(tmp_path, case):
    if case == "removed":
        o, n, cns = [NONPAR, CRLF2_ORIG], [NONPAR], call_cns(tmp_path)
    elif case == "kept":
        o, n = [NONPAR, CRLF2_ORIG.replace("1.26805", "1.5")], [NONPAR, CRLF2_FIXED.replace("0.26805", "0.5")]
        cns = call_cns(tmp_path, "0.5")
    elif case == "span":
        o, n, cns = [NONPAR, CRLF2_ORIG], [NONPAR, CRLF2_ORIG], call_cns(tmp_path, "0.26805", end=9000000)
    else:
        o, n, cns = [CRLF2_ORIG], [], call_cns(tmp_path)
    orig = gm(tmp_path, "s.genemetrics.noparfix.tsv", *o)
    new = gm(tmp_path, "s.genemetrics.tsv", *n, header=port.EMPTY_GM_HEADER if case == "empty" else GM_H)
    (tmp_path / "s.cnr").write_text(CNR_H)
    mine = port.compare(orig, new, cns)
    theirs = ref.compare(dict(action="dry", rundir=str(tmp_path), sample="s", orig=str(orig), new=str(new)))
    assert set(mine) - set(theirs) == ({"span_genes"} if case == "span" else set())   # port addition (A1)
    for k, v in mine.items():
        if k != "span_genes":
            assert theirs[k] == v, k
