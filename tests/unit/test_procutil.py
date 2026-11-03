"""L1: the bounded-subprocess helper used by the tests (procutil.run_group)."""
import subprocess
import sys

import pytest

from procutil import run_group

pytestmark = pytest.mark.unit


def test_completes_like_subprocess_run():
    p = run_group([sys.executable, "-c", "import sys; print('out'); print('err', file=sys.stderr); sys.exit(3)"],
                  30, capture_output=True, text=True)
    assert (p.returncode, p.stdout, p.stderr) == (3, "out\n", "err\n")
    with pytest.raises(subprocess.CalledProcessError):
        run_group([sys.executable, "-c", "raise SystemExit(1)"], 30, check=True)
    assert run_group(["cat"], 30, input="in", capture_output=True, text=True).stdout == "in"


def test_timeout_kills_the_group_and_calls_on_timeout(tmp_path):
    called = []
    with pytest.raises(subprocess.TimeoutExpired):
        run_group(["bash", "-c", "sleep 300 & wait"], 1, on_timeout=lambda: called.append(1), capture_output=True)
    assert called == [1]
