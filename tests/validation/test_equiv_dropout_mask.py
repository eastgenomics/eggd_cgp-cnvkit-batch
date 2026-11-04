"""L3: dropout_mask.py, run through its command line as pipeline.sh runs it, is byte-identical to the pinned
`deadprobe_mask.py mask M` (ddad79b) on the 226 validation sequencings, on UMI-consensus and non-UMI coverage, and
reproduces the validated mask_qc.tsv and totals.

Each sequencing (keyed by run and sample: replicates share a sample ID across runs) writes into its own directory:
  dropout_mask.py mask <s>.cnr <s>.targetcoverage.cnn <s>.seg_input.cnr <s>.dropout_qc.tsv
and its masked-probe count is the QC table's rows, as pipeline.sh counts DROPOUT_PROBES (the printed count must agree).
"""
import pytest

from vdata import DATA, DATASETS, QC_FAILED, REFERENCE, load_module, parallel, require, script, sequencings

pytestmark = pytest.mark.validation

# masked probes / sequencings with >= 1 masked probe, excluding the QC-failed library (DI-4139 comment 72496)
EXPECTED_TOTALS = {"umi": (82, 42), "nonumi": (78, 41)}
EXPECTED_QC_FAILED_BINS = {"umi": 2091, "nonumi": 2075}


@pytest.fixture(scope="module")
def reference():
    return load_module(REFERENCE / "deadprobe_mask.py", "deadprobe_mask_ref")


@pytest.mark.parametrize("dataset", DATASETS)
def test_dropout_mask_equivalence(tmp_path, reference, dataset):
    require(DATA / dataset / "calls" / "M")
    seqs = sequencings()

    def one(key):
        run, s = key
        d = tmp_path / run / s
        d.mkdir(parents=True)
        cnr = DATA / dataset / "calls" / "M" / run / f"{s}.cnr"
        cnn = DATA / dataset / "cnn" / run / f"{s}.targetcoverage.cnn"
        expected_qc = DATA / dataset / "calls" / "M" / run / f"{s}.mask_qc.tsv"
        fails = []
        # the arm M .cnr is the unmasked file: identical to A0
        if cnr.read_bytes() != (DATA / dataset / "calls" / "A0" / run / f"{s}.cnr").read_bytes():
            fails.append("M .cnr differs from A0 .cnr")
        p = script("dropout_mask.py", "mask", cnr, cnn, f"{s}.seg_input.cnr", f"{s}.dropout_qc.tsv", cwd=d)
        if p.returncode:
            return [f"{run}/{s}: dropout_mask.py mask: exit {p.returncode}: {p.stderr.strip()}"], 0
        out_p, qc_p = d / f"{s}.seg_input.cnr", d / f"{s}.dropout_qc.tsv"
        out_r, qc_r = d / "ref.cnr", d / "ref.qc.tsv"
        reference.mask("M", str(cnr), str(cnn), s, "-", str(out_r), str(qc_r))
        n = len(qc_p.read_text().splitlines()) - 1       # pipeline.sh: $(wc -l < qc) - 1
        if p.stdout != f"[dropout] masked probes: {n}\n":
            fails.append(f"printed {p.stdout!r}, the QC table has {n} rows")
        if out_p.read_bytes() != out_r.read_bytes():
            fails.append("masked .cnr differs from the reference")
        if qc_p.read_bytes() != qc_r.read_bytes():
            fails.append("QC differs from the reference")
        if qc_p.read_bytes() != expected_qc.read_bytes():
            fails.append("QC differs from the validated mask_qc.tsv")
        return [f"{run}/{s}: {f}" for f in fails], n

    results = dict(zip(seqs, parallel(one, seqs), strict=True))
    mismatches = [f for fails, _ in results.values() for f in fails]
    assert not mismatches, "\n".join(mismatches[:20])
    counted = [n for (_, s), (_, n) in results.items() if QC_FAILED not in s and n]
    assert (sum(counted), len(counted)) == EXPECTED_TOTALS[dataset]
    assert sum(n for (_, s), (_, n) in results.items() if QC_FAILED in s) == EXPECTED_QC_FAILED_BINS[dataset]
