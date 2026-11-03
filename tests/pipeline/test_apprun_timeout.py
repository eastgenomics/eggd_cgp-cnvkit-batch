"""L2 harness check: a timed-out app run is killed with every descendant, not only the bash leader."""
import subprocess
import time

import pytest

from apprun import run_app

pytestmark = pytest.mark.pipeline


def alive(pid):
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] != "Z"    # a zombie has exited
    except FileNotFoundError:
        return False


def test_timeout_kills_descendants(tmp_path):
    pids = tmp_path / "pids"
    script = tmp_path / "hang.sh"
    script.write_text('main() {\n'
                      '    sleep 300 & echo "$!" >> "$PIDFILE"\n'
                      '    bash -c \'sleep 300 & echo "$!" >> "$PIDFILE"; wait\' & echo "$!" >> "$PIDFILE"\n'
                      '    wait\n'
                      '}\n')
    with pytest.raises(subprocess.TimeoutExpired):
        run_app(tmp_path, script=script, stub_env={"PIDFILE": str(pids)}, timeout=3)
    started = [int(p) for p in pids.read_text().split()]
    assert len(started) == 3                      # a child, a grandchild's parent and the grandchild
    deadline = time.monotonic() + 10
    while any(alive(p) for p in started) and time.monotonic() < deadline:
        time.sleep(0.1)
    assert not [p for p in started if alive(p)]
