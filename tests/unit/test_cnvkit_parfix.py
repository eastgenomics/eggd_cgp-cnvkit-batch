"""L1: cnvkit_parfix.py (DI-4101) — PAR geometry, shift, male decision, compare; differential against the pinned
reference `genemetrics_parfix.py` (f933f7a) on synthetic inputs."""
import random
import re
import subprocess
import sys

import pytest

from cnvfiles import BIN, REFERENCE

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
    return m


port, ref = load_port(), load_reference()


def run(*argv):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, argv)], capture_output=True, text=True)


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


# ── compare ──────────────────────────────────────────────────────────────────

NONPAR = "TP53\tchr17\t7661779\t7687538\t-0.62\t400\t5\t1\t10\t5\t10\n"
CRLF2_ORIG = "CRLF2\tchrX\t1190891\t1212654\t1.26805\t679\t7.2\t3\t8\t7.2\t8\n"
CRLF2_FIXED = "CRLF2\tchrX\t1190891\t1212654\t0.26805\t679\t7.2\t2\t8\t7.2\t8\n"


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
    assert proc.returncode == 1 and "headers differ" in proc.stderr


def test_compare_accepts_empty_table_header(tmp_path):
    rec = port.compare(gm(tmp_path, "o", CRLF2_ORIG), gm(tmp_path, "n", header=port.EMPTY_GM_HEADER),
                       call_cns(tmp_path))
    assert rec["empty_after"] == 1 and rec["par_check"] == 1


def test_compare_span(tmp_path):
    # the CRLF2 segment continues into non-PAR chrX: not shifted, genemetrics unchanged -> par_check span
    cns = call_cns(tmp_path, "1.26805", start=1190891, end=9000000)
    rec = port.compare(gm(tmp_path, "o", NONPAR, CRLF2_ORIG), gm(tmp_path, "n", NONPAR, CRLF2_ORIG), cns)
    assert rec["par_check"] == "span" and rec["file_changed"] == 0


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
        o, n, cns = [NONPAR, CRLF2_ORIG], [NONPAR, CRLF2_ORIG], call_cns(tmp_path, "1.26805", end=9000000)
    else:
        o, n, cns = [CRLF2_ORIG], [], call_cns(tmp_path)
    orig = gm(tmp_path, "s.genemetrics.noparfix.tsv", *o)
    new = gm(tmp_path, "s.genemetrics.tsv", *n, header=port.EMPTY_GM_HEADER if case == "empty" else GM_H)
    (tmp_path / "s.cnr").write_text(CNR_H)
    mine = port.compare(orig, new, cns)
    theirs = ref.compare(dict(action="dry", rundir=str(tmp_path), sample="s", orig=str(orig), new=str(new)))
    for k, v in mine.items():
        assert theirs[k] == v, k
