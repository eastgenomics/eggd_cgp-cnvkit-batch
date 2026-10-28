"""L1: dropout_mask.py (DI-4139) — rule boundaries, output format, CLI, and byte-identity with the pinned
reference `deadprobe_mask.py mask M` (ddad79b) on synthetic inputs."""
import random
import subprocess
import sys

import pytest

from cnvfiles import BIN, CNN_HEADER, CNR_HEADER, REFERENCE, bins, cnn_rows, cnr_rows, write_table

pytestmark = pytest.mark.unit

SCRIPT = BIN / "dropout_mask.py"
QC_HEADER = "arm\tchromosome\tstart\tend\tgene\tlog2\traw_reads\treason\n"


def run_mask(tmp_path, bin_list, depth=0.0, script=SCRIPT, tag="new"):
    cnr = write_table(tmp_path / "s.cnr", CNR_HEADER, cnr_rows(bin_list))
    cnn = write_table(tmp_path / "s.targetcoverage.cnn", CNN_HEADER, cnn_rows(bin_list, depth))
    out, qc = tmp_path / f"{tag}.seg_input.cnr", tmp_path / f"{tag}.dropout_qc.tsv"
    if script == SCRIPT:
        args = ["mask", str(cnr), str(cnn), str(out), str(qc)]
    else:  # the pinned reference retains the validation CLI
        args = ["mask", "M", str(cnr), str(cnn), "spec", "-", str(out), str(qc)]
    proc = subprocess.run([sys.executable, str(script), *args], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return cnr, out, qc, proc.stdout


def masked(out):
    """Indices of rows with weight 0 in a masked .cnr."""
    rows = out.read_text().splitlines()[1:]
    return [i for i, r in enumerate(rows) if r.split("\t")[6] == "0"]


def qc_rows(qc):
    return [l.split("\t") for l in qc.read_text().splitlines()[1:]]


# ── rule ─────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("log2, is_masked", [(-3.0, False), (-3.0001, True), (-2.99, False), (-28.8, True)])
def test_deep_threshold(tmp_path, log2, is_masked):
    _, out, _, _ = run_mask(tmp_path, bins([0.1, log2, 0.1]))
    assert masked(out) == ([1] if is_masked else [])


@pytest.mark.parametrize("depth, is_masked", [
    (0.0, True),
    (9.99, True),    # 9.99 read-equivalents (length 144)
    (10.0, False),   # exactly 10 is not < 10
    (25.0, False),
])
def test_raw_read_threshold(tmp_path, depth, is_masked):
    # with length 144, raw read-equivalents = depth * 144 / 144 = depth
    _, out, _, _ = run_mask(tmp_path, bins([0.1, -10.0, 0.1]), depth=depth)
    assert masked(out) == ([1] if is_masked else [])


def test_raw_reads_scale_with_bin_length(tmp_path):
    # depth 6 on a 288 bp bin = 12 read-equivalents (kept); on a 144 bp bin = 6 (masked)
    long_bin = bins([0.1, -10.0, 0.1], length=288)
    _, out, _, _ = run_mask(tmp_path, long_bin, depth=6.0)
    assert masked(out) == []


@pytest.mark.parametrize("run, expect", [
    ([0.1, -9, 0.1], [1]),                       # run of 1
    ([0.1, -9, -9, 0.1], [1, 2]),                # run of 2
    ([0.1, -9, -9, -9, 0.1], []),                # run of 3: a callable deletion, kept
    ([-9, 0.1, -9, -9, 0.1, -9, -9, -9], [0, 2, 3]),
])
def test_run_length(tmp_path, run, expect):
    _, out, _, _ = run_mask(tmp_path, bins(run))
    assert masked(out) == expect


def test_run_broken_by_chromosome_change(tmp_path):
    # 2 deep bins at the end of chr1 and 2 at the start of chr2: two runs of 2, all masked
    b = bins([0.1, -9, -9], chrom="chr1") + bins([-9, -9, 0.1], chrom="chr2")
    _, out, qc, _ = run_mask(tmp_path, b)
    assert masked(out) == [1, 2, 3, 4]
    assert {r[7] for r in qc_rows(qc)} == {"isolated_run2"}


@pytest.mark.parametrize("chrom", ["chrX", "chrY"])
def test_sex_chromosomes_never_masked(tmp_path, chrom):
    _, out, _, _ = run_mask(tmp_path, bins([0.1, -20, 0.1], chrom=chrom))
    assert masked(out) == []


def test_sex_chromosome_bin_ends_an_autosomal_run(tmp_path):
    # the deep chrX bin is on another chromosome, so the chr22 bin is a run of 1
    b = bins([0.1, -9], chrom="chr22") + bins([-9, 0.1], chrom="chrX")
    _, out, _, _ = run_mask(tmp_path, b)
    assert masked(out) == [1]


def test_bin_absent_from_coverage_counts_as_zero_reads(tmp_path):
    b = bins([0.1, -9, 0.1])
    cnr = write_table(tmp_path / "s.cnr", CNR_HEADER, cnr_rows(b))
    cnn = write_table(tmp_path / "s.cnn", CNN_HEADER, cnn_rows([b[0], b[2]], 500.0))   # deep bin missing
    out, qc = tmp_path / "o.cnr", tmp_path / "q.tsv"
    subprocess.run([sys.executable, str(SCRIPT), "mask", str(cnr), str(cnn), str(out), str(qc)],
                   check=True, capture_output=True)
    assert masked(out) == [1]
    assert qc_rows(qc)[0][6] == "0.0"


def test_deep_bin_with_enough_reads_keeps_its_run_partner_masked(tmp_path):
    # a run of 2 where one bin has >= 10 reads: only the other one is masked
    b = bins([0.1, -9, -9, 0.1])
    _, out, _, _ = run_mask(tmp_path, b, depth={(b[2][0], b[2][1]): 50.0, (b[1][0], b[1][1]): 0.0})
    assert masked(out) == [1]


# ── output format ────────────────────────────────────────────────────────────

def test_only_weight_changes(tmp_path):
    cnr, out, _, _ = run_mask(tmp_path, bins([0.1, -9, 0.2, -12, -15, 0.3]))
    before, after = cnr.read_text().splitlines(), out.read_text().splitlines()
    assert before[0] == after[0]
    for i, (a, b) in enumerate(zip(before[1:], after[1:])):
        a, b = a.split("\t"), b.split("\t")
        assert a[:6] == b[:6]
        assert b[6] == ("0" if i in (1, 3, 4) else a[6])


def test_qc_rows(tmp_path):
    b = bins([0.1, -28.7962, -10.6188, 0.1], chrom="chr12", start=11873311, length=272, gap=0, gene="ETV6_")
    _, _, qc, stdout = run_mask(tmp_path, b, depth={("chr12", 11873583): 0.0, ("chr12", 11873855): 0.2647})
    assert qc.read_text() == QC_HEADER + (
        "M\tchr12\t11873583\t11873855\tETV6_1\t-28.7962\t0.0\tisolated_run2\n"
        "M\tchr12\t11873855\t11874127\tETV6_2\t-10.6188\t0.5\tisolated_run2\n")
    assert stdout.strip() == "[dropout] masked probes: 2"


def test_header_only_qc_when_nothing_masked(tmp_path):
    cnr, out, qc, stdout = run_mask(tmp_path, bins([0.1, -1.0, 0.2]))
    assert qc.read_text() == QC_HEADER
    assert out.read_bytes() == cnr.read_bytes()
    assert stdout.strip() == "[dropout] masked probes: 0"


def test_qc_header_subcommand_matches_mask_output(tmp_path):
    _, _, qc, _ = run_mask(tmp_path, bins([0.1, 0.2]))
    hdr = tmp_path / "header_only.tsv"
    subprocess.run([sys.executable, str(SCRIPT), "qc-header", str(hdr)], check=True)
    assert hdr.read_bytes() == qc.read_bytes()


# ── CLI ──────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("argv, message", [
    (["mask", "M", "a", "b", "s", "-", "o", "q"], "Usage"),  # old validation CLI
    (["mask", "a", "b", "o"], "Usage"),  # missing QC file
    (["mask", "a", "b", "o", "q", "extra"], "Usage"),
    (["blacklist", "out.tsv", "map.tsv", "a.cnr"], "Usage"),
    ([], "Usage"),
])
def test_cli_rejects(tmp_path, argv, message):
    proc = subprocess.run([sys.executable, str(SCRIPT), *argv], capture_output=True, text=True, cwd=tmp_path)
    assert proc.returncode != 0
    assert message in proc.stderr


# ── differential: byte-identical to the pinned reference on synthetic data ──

def random_case(rng):
    """A few chromosomes of bins with a mix of normal, deep, borderline and long-run values."""
    out = []
    for chrom in ["chr1", "chr12", "chrX", "chr22", "chrY"]:
        spec, n = [], rng.randint(5, 40)
        while len(spec) < n:
            kind = rng.random()
            if kind < 0.75:
                spec.append(round(rng.gauss(0, 0.3), 6))
            else:
                spec += [rng.choice([-3.0, -3.0001, -4.2, -9.5, -28.7962, -20.0])] * rng.choice([1, 1, 2, 2, 3, 4])
        out += bins(spec, chrom=chrom, length=rng.choice([120, 144, 272]))
    depth = {(c, s): rng.choice([0.0, 3.3, 9.9, 10.0, 12.5, 500.0]) for c, s, *_ in out}
    return out, depth


@pytest.mark.parametrize("seed", range(25))
def test_byte_identical_to_reference(tmp_path, seed):
    b, depth = random_case(random.Random(seed))
    _, new_out, new_qc, _ = run_mask(tmp_path, b, depth)
    _, ref_out, ref_qc, _ = run_mask(tmp_path, b, depth, script=REFERENCE / "deadprobe_mask.py", tag="ref")
    assert new_out.read_bytes() == ref_out.read_bytes()
    assert new_qc.read_bytes() == ref_qc.read_bytes()
