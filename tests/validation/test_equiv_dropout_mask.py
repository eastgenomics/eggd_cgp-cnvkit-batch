"""L3: dropout_mask.py is byte-identical to the pinned `deadprobe_mask.py mask M` (ddad79b) on the 226 validation
sequencings, on UMI-consensus and non-UMI coverage, and reproduces the validated mask_qc.tsv and totals."""
import pytest

from vdata import BIN, DATA, DATASETS, QC_FAILED, REFERENCE, load_module, require, sequencings

pytestmark = pytest.mark.validation

# masked probes / sequencings with >= 1 masked probe, excluding the QC-failed library (DI-4139 comment 72496)
EXPECTED_TOTALS = {"umi": (82, 42), "nonumi": (78, 41)}
EXPECTED_QC_FAILED_BINS = {"umi": 2091, "nonumi": 2075}


@pytest.fixture(scope="module")
def ported():
    return load_module(BIN / "dropout_mask.py", "dropout_mask_port")


@pytest.fixture(scope="module")
def reference():
    return load_module(REFERENCE / "deadprobe_mask.py", "deadprobe_mask_ref")


@pytest.mark.parametrize("dataset", DATASETS)
def test_dropout_mask_equivalence(tmp_path, ported, reference, dataset):
    require(DATA / dataset / "calls" / "M")
    probes = seqs = qc_failed_bins = 0
    mismatches = []
    for run, s in sequencings():
        cnr = DATA / dataset / "calls" / "M" / run / f"{s}.cnr"
        cnn = DATA / dataset / "cnn" / run / f"{s}.targetcoverage.cnn"
        expected_qc = DATA / dataset / "calls" / "M" / run / f"{s}.mask_qc.tsv"
        # the arm M .cnr is the unmasked file: identical to A0
        if cnr.read_bytes() != (DATA / dataset / "calls" / "A0" / run / f"{s}.cnr").read_bytes():
            mismatches.append(f"{run}/{s}: M .cnr differs from A0 .cnr")
        out_p, qc_p = tmp_path / "port.cnr", tmp_path / "port.qc.tsv"
        out_r, qc_r = tmp_path / "ref.cnr", tmp_path / "ref.qc.tsv"
        n = ported.mask(str(cnr), str(cnn), str(out_p), str(qc_p))
        reference.mask("M", str(cnr), str(cnn), s, "-", str(out_r), str(qc_r))
        if out_p.read_bytes() != out_r.read_bytes():
            mismatches.append(f"{run}/{s}: masked .cnr differs from the reference")
        if qc_p.read_bytes() != qc_r.read_bytes():
            mismatches.append(f"{run}/{s}: QC differs from the reference")
        if qc_p.read_bytes() != expected_qc.read_bytes():
            mismatches.append(f"{run}/{s}: QC differs from the validated mask_qc.tsv")
        if QC_FAILED in s:
            qc_failed_bins += n
        elif n:
            probes += n
            seqs += 1
    assert not mismatches, "\n".join(mismatches[:20])
    assert (probes, seqs) == EXPECTED_TOTALS[dataset]
    assert qc_failed_bins == EXPECTED_QC_FAILED_BINS[dataset]
