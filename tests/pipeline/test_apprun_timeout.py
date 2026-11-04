"""L2 harness check: a timed-out app run is killed with every descendant, not only the bash leader."""
import io
import os
import subprocess
import time

import pytest

from apprun import run_app

pytestmark = pytest.mark.pipeline


def alive(pid, opener=open):
    """False once the process has exited (zombie, or /proc/<pid> gone). /proc/<pid>/stat can vanish between open and
    read, which raises ProcessLookupError (ESRCH): that also means exited."""
    try:
        with opener(f"/proc/{pid}/stat") as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] != "Z"    # a zombie has exited
    except (FileNotFoundError, ProcessLookupError):
        return False


READY_TIMEOUT = 60   # seconds for the three descendants to start and report their PIDs
APP_TIMEOUT = 3      # seconds the app then runs before it is timed out
EXIT_TIMEOUT = 10    # seconds for the killed descendants to exit


def proc_state(pid):
    """/proc/<pid>/stat state and command, for the failure report."""
    try:
        with open(f"/proc/{pid}/stat") as fh:
            head, rest = fh.read().rsplit(")", 1)
        return f"{pid}: state {rest.split()[0]} {head.split('(', 1)[1]}"
    except (FileNotFoundError, ProcessLookupError):
        return f"{pid}: gone"


def start_timeout_after_ready(monkeypatch, ready):
    """Readiness handshake: the app's timeout (the first Popen.communicate) starts only once ready() is true, so
    the run is never timed out before every descendant has started. Bounded by READY_TIMEOUT."""
    real, state = subprocess.Popen.communicate, {}

    def communicate(self, *a, **k):
        if "waited" not in state:
            deadline = time.monotonic() + READY_TIMEOUT
            while not ready() and time.monotonic() < deadline and self.poll() is None:
                time.sleep(0.05)
            state["waited"] = time.monotonic() - (deadline - READY_TIMEOUT)
            state["ready"] = ready()
        return real(self, *a, **k)
    monkeypatch.setattr(subprocess.Popen, "communicate", communicate)
    return state


def test_timeout_kills_descendants(tmp_path, monkeypatch):
    pids, ready = tmp_path / "pids", tmp_path / "ready"
    script = tmp_path / "hang.sh"
    script.write_text('main() {\n'
                      '    sleep 300 & echo "$!" >> "$PIDFILE"\n'
                      '    bash -c \'sleep 300 & echo "$!" >> "$PIDFILE"; wait\' & echo "$!" >> "$PIDFILE"\n'
                      '    while [ "$(wc -l < "$PIDFILE")" -lt 3 ]; do sleep 0.05; done\n'
                      '    : > "$READY"\n'
                      '    wait\n'
                      '}\n')
    handshake = start_timeout_after_ready(monkeypatch, ready.exists)
    with pytest.raises(subprocess.TimeoutExpired) as e:
        run_app(tmp_path, script=script, stub_env={"PIDFILE": str(pids), "READY": str(ready)}, timeout=APP_TIMEOUT)
    started = [int(p) for p in pids.read_text().split()] if pids.exists() else []
    deadline = time.monotonic() + EXIT_TIMEOUT
    while any(alive(p) for p in started) and time.monotonic() < deadline:
        time.sleep(0.1)

    def report():
        return (f"handshake {handshake}; pids {started}; "
                f"{'; '.join(proc_state(p) for p in started)}; notes {getattr(e.value, 'cleanup_notes', None)}; "
                f"stderr {(e.value.stderr or b'')[-2000:]!r}")
    assert handshake.get("ready"), f"descendants not ready within {READY_TIMEOUT} s: {report()}"
    assert len(started) == 3, report()            # a child, a grandchild's parent and the grandchild
    assert not [p for p in started if alive(p)], report()


class VanishingStat(io.StringIO):
    def read(self, *a):
        raise ProcessLookupError(3, "No such process")


def test_alive_treats_a_vanishing_proc_entry_as_exited():
    # R5: the race seen in CI, /proc/<pid>/stat opened and then gone before the read
    assert alive(12345, opener=lambda path: VanishingStat()) is False
    assert alive(2 ** 22 + 1) is False                     # above pid_max: no /proc entry
    assert alive(os.getpid()) is True
