"""L1 (round 6, V3): the L3 helpers key every per-sequencing result by (run, sample) and require complete, unique
coverage, so a replicate (same sample ID, another run) cannot overwrite or stand in for another."""
import csv

import pytest

import vdata

pytestmark = pytest.mark.unit

S = "135024088-25021S0072-26CGPVal1SHv150kb-9999-U-1111"
REPLICATES = [{"run": "0734_Val1_50kb", "sample": S, "sex": "Male"},
              {"run": "0811_V3", "sample": S, "sex": "Female"}]     # conflicting results for one sample ID


def test_replicates_with_conflicting_results_are_kept_apart():
    by_seq = vdata.index(REPLICATES, lambda r: (r["run"], r["sample"]), "sex calls")
    assert by_seq[("0734_Val1_50kb", S)]["sex"] == "Male" and by_seq[("0811_V3", S)]["sex"] == "Female"
    vdata.covers(by_seq, [("0734_Val1_50kb", S), ("0811_V3", S)], "sex calls")


def test_a_sample_id_key_cannot_silently_overwrite_a_replicate():
    # keyed by the sample ID alone (the round-5 hole), the second replicate would replace the first
    with pytest.raises(AssertionError, match=f"sex calls: duplicate entry for {S}"):
        vdata.index(REPLICATES, lambda r: r["sample"], "sex calls")


@pytest.mark.parametrize("keys, message", [
    ([("0734_Val1_50kb", S)], "missing \\[\\('0811_V3'"),
    ([("0734_Val1_50kb", S), ("0811_V3", S), ("0786_V1", S)], "extra \\[\\('0786_V1'"),
])
def test_incomplete_or_extra_coverage_fails(keys, message):
    with pytest.raises(AssertionError, match=message):
        vdata.covers(keys, [("0734_Val1_50kb", S), ("0811_V3", S)], "log")


def test_sequencings_requires_all_226(tmp_path, monkeypatch):
    manifest = tmp_path / "manifest" / "11_final_per_sample_files.csv"
    manifest.parent.mkdir()
    with open(manifest, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["record_id", "call_set"])
        w.writerows([[f"0734_Val1_50kb/{S}", "main"], [f"0811_V3/{S}", "main"], [f"0811_V3/{S}", "main"]])
    monkeypatch.setattr(vdata, "DATA", tmp_path)
    with pytest.raises(AssertionError, match="2 sequencings, expected 226"):
        vdata.sequencings()


# round 7 (W4): the L3 comparisons require every statistic, not only those both sides happen to have

COMPARE_OUT = ("par_segs=1 span_segs=0 empty_after=1 nonpar_identical=1 par_rows_before=1 par_rows_after=0 "
               "rows_changed=1 file_changed=1 par_before=CRLF2:1.26805 par_after= par_check=1\n")
LOG_ROW = {"rundir": "/x/d00/R1", "sample": "S", "par_segs": "1", "span_segs": "0", "empty_after": "1",
           "nonpar_identical": "1", "par_rows_before": "1", "par_rows_after": "0", "rows_changed": "1",
           "file_changed": "1", "par_before": "CRLF2:1.26805", "par_after": "", "par_check": "1"}


def test_stats_parses_the_full_set():
    got = vdata.stats(COMPARE_OUT, vdata.COMPARE_KEYS)
    assert got == vdata.logged(LOG_ROW, vdata.COMPARE_KEYS) and got["par_after"] == ""


@pytest.mark.parametrize("key", vdata.COMPARE_KEYS)
def test_a_missing_statistic_fails(key):
    # the reviewer's case: removing e.g. rows_changed or nonpar_identical from the output must not pass
    out = " ".join(kv for kv in COMPARE_OUT.split() if kv.split("=", 1)[0] != key)
    with pytest.raises(AssertionError, match="expected exactly the statistics"):
        vdata.stats(out, vdata.COMPARE_KEYS)


@pytest.mark.parametrize("out", [
    COMPARE_OUT + " span_genes=CRLF2@chrX:1-2",           # a port statistic not expected here
    COMPARE_OUT + " par_check=1",                          # repeated
    COMPARE_OUT + " stray",                                # not key=value
    COMPARE_OUT + " unexplained=CRLF2@chrX:1-2",
])
def test_an_unexpected_or_repeated_statistic_fails(out):
    with pytest.raises(AssertionError, match="expected exactly the statistics"):
        vdata.stats(out, vdata.COMPARE_KEYS)


def test_the_port_statistic_is_required_when_expected():
    out = COMPARE_OUT.replace("par_check=1", "span_genes=CRLF2@chrX:1-9 par_check=span")
    assert vdata.stats(out, vdata.COMPARE_KEYS, vdata.COMPARE_PORT_KEYS)["span_genes"] == "CRLF2@chrX:1-9"
    with pytest.raises(AssertionError):
        vdata.stats(COMPARE_OUT, vdata.COMPARE_KEYS, vdata.COMPARE_PORT_KEYS)


@pytest.mark.parametrize("key", vdata.COMPARE_KEYS)
def test_a_statistic_missing_from_the_log_fails(key):
    with pytest.raises(AssertionError, match=f"the validation log has no {key}"):
        vdata.logged({k: v for k, v in LOG_ROW.items() if k != key}, vdata.COMPARE_KEYS)


def test_the_l3_key_sets_are_what_cnvkit_parfix_prints(tmp_path):
    # the constants stay tied to the script: its shift and compare output parse with exactly these keys
    cnr_h = "chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n"
    cns_h = "chromosome\tstart\tend\tgene\tlog2\tcn\tdepth\tprobes\tweight\n"
    gm_h = "gene\tchromosome\tstart\tend\tlog2\tdepth\tweight\tcn\tprobes\tsegment_weight\tsegment_probes\n"
    (tmp_path / "s.cnr").write_text(cnr_h + "chrX\t1190891\t1191000\tCRLF2\t1.26805\t100\t1\n")
    for end, check in ((2781479, "1"), (9000000, "span")):
        (tmp_path / "s.call.cns").write_text(cns_h + f"chrX\t1190891\t{end}\tCRLF2\t0.26805\t2\t679\t8\t7.2\n")
        crlf2 = "CRLF2\tchrX\t1190891\t1212654\t1.26805\t679\t7.2\t3\t8\t7.2\t8\n"
        (tmp_path / "o.tsv").write_text(gm_h + crlf2)
        (tmp_path / "n.tsv").write_text(gm_h + (crlf2 if check == "span" else ""))
        p = vdata.script("cnvkit_parfix.py", "shift", "s.cnr", "s.call.cns", "parfix", cwd=tmp_path)
        assert p.returncode == 0, p.stderr
        assert set(vdata.stats(p.stdout, vdata.SHIFT_KEYS)) == set(vdata.SHIFT_KEYS)
        p = vdata.script("cnvkit_parfix.py", "compare", "o.tsv", "n.tsv", "s.call.cns", cwd=tmp_path)
        assert p.returncode == 0, p.stderr
        got = vdata.stats(p.stdout, vdata.COMPARE_KEYS, vdata.COMPARE_PORT_KEYS if check == "span" else ())
        assert got["par_check"] == check
