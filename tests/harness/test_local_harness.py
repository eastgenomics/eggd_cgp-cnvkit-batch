"""L4: the local harness (run_local.py) on the 16 local test sequencings, from the validated coverage, in the four
data x dropout-mask configurations. The BAM-start run (4.8 GB download) is manual:
    tests/harness/run_local.py --data nonumi --start-from bam --samples 25328S0023
"""
import pytest

import run_local
from test_equiv_parfix import docker_available
from vdata import DATA, require

pytestmark = [pytest.mark.validation, pytest.mark.docker]


@pytest.mark.skipif(not docker_available(), reason="docker or cgp-cnvkit:1.0.0 not available")
@pytest.mark.parametrize("data, mask", [("umi", "true"), ("umi", "false"), ("nonumi", "false"), ("nonumi", "true")])
def test_local_harness(data, mask):
    require(DATA / data / "calls")
    assert run_local.main(["--data", data, "--dropout-mask", mask, "--label", f"pytest-{data}-mask_{mask}"]) == 0
