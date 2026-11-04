"""L4: the local harness (run_local.py) on the 16 local test sequencings, from the validated coverage, in the four
data x dropout-mask configurations. The BAM-start run (4.8 GB download) is manual:
    tests/harness/run_local.py --data nonumi --start-from bam --samples 25328S0023

Each configuration must compare all 16 sequencings, each once with at least one PASS and no FAIL (round 7, W2): a run
that exercised fewer, or none, cannot pass.
"""
import csv

import pytest

import run_local
from vdata import DATA, require

pytestmark = [pytest.mark.validation, pytest.mark.docker]


@pytest.mark.parametrize("data, mask", [("umi", "true"), ("umi", "false"), ("nonumi", "false"), ("nonumi", "true")])
def test_local_harness(data, mask):
    require(DATA / data / "calls")
    label = f"pytest-{data}-mask_{mask}"
    assert run_local.main(["--data", data, "--dropout-mask", mask, "--label", label]) == 0
    with open(run_local.HARNESS / label / "report.tsv") as fh:
        report = list(csv.DictReader(fh, delimiter="\t"))
    expected = {r["sample"] for r in run_local.local_samples()}
    assert len(expected) == run_local.N_LOCAL == 16
    assert {r["sample"] for r in report} == expected
    assert {r["sample"] for r in report if r["result"] == "PASS"} == expected
    assert not [r for r in report if r["result"] not in ("PASS", "n/a")]
