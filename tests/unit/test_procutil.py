"""L1: the bounded-subprocess helper used by the tests (procutil.run_group, remove_container)."""
import os
import re
import signal
import subprocess
import sys
import time

import pytest

import procutil
from procutil import remove_container, run_group

pytestmark = pytest.mark.unit


def test_completes_like_subprocess_run():
    p = run_group([sys.executable, "-c", "import sys; print('out'); print('err', file=sys.stderr); sys.exit(3)"],
                  30, capture_output=True, text=True)
    assert (p.returncode, p.stdout, p.stderr) == (3, "out\n", "err\n")
    with pytest.raises(subprocess.CalledProcessError):
        run_group([sys.executable, "-c", "raise SystemExit(1)"], 30, check=True)
    assert run_group(["cat"], 30, input="in", capture_output=True, text=True).stdout == "in"


def test_timeout_kills_the_group_and_calls_cleanup():
    called = []
    with pytest.raises(subprocess.TimeoutExpired) as e:
        run_group(["bash", "-c", "sleep 300 & wait"], 1, cleanup=lambda: called.append(1), capture_output=True)
    assert called == [1] and e.value.cleanup_notes == []


def test_group_kill_failure_falls_back_and_still_cleans_up(monkeypatch):
    # R4: killpg raising something other than ProcessLookupError must not skip the cleanup or hang
    def denied(pid, sig):
        raise PermissionError(1, "Operation not permitted")
    monkeypatch.setattr(procutil.os, "killpg", denied)
    called = []
    start = time.monotonic()
    with pytest.raises(subprocess.TimeoutExpired) as e:
        run_group(["sleep", "300"], 1, cleanup=lambda: called.append(1), capture_output=True)
    assert time.monotonic() - start < 1 + procutil.GRACE
    assert called == [1]
    assert any("killpg" in n and "PermissionError" in n for n in e.value.cleanup_notes)
    assert any("killpg" in n for n in getattr(e.value, "__notes__", []))


def test_unkillable_process_is_bounded_and_reported(tmp_path, monkeypatch):
    # both kills fail: collecting output is bounded by GRACE, the cleanup still runs, and both failures are reported
    def denied(*a):
        raise PermissionError(1, "Operation not permitted")
    monkeypatch.setattr(procutil.os, "killpg", denied)
    monkeypatch.setattr(subprocess.Popen, "kill", lambda self: denied())
    monkeypatch.setattr(procutil, "GRACE", 1)
    pidfile = tmp_path / "pid"
    called = []
    start = time.monotonic()
    try:
        with pytest.raises(subprocess.TimeoutExpired) as e:
            run_group(["bash", "-c", f'echo $$ > "{pidfile}"; exec sleep 60'], 1,
                      cleanup=lambda: called.append(1), capture_output=True)
        assert time.monotonic() - start < 10
        assert called == [1]
        notes = e.value.cleanup_notes
        assert any("killpg" in n for n in notes) and any("kill(" in n for n in notes)
        assert any("still running" in n for n in notes)
    finally:
        os.kill(int(pidfile.read_text()), signal.SIGKILL)


def test_cleanup_failure_is_reported_without_masking_the_timeout():
    def fail():
        raise RuntimeError("docker rm --force x: exit 1: permission denied")
    with pytest.raises(subprocess.TimeoutExpired) as e:
        run_group(["sleep", "300"], 1, cleanup=fail)
    assert e.value.cleanup_notes == ["cleanup unconfirmed: docker rm --force x: exit 1: permission denied"]


@pytest.fixture
def docker_stub(tmp_path, monkeypatch):
    """A `docker` on PATH whose `rm` behaviour is set per test: write(script body)."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")

    def write(body):
        stub = bin_dir / "docker"
        stub.write_text("#!/bin/sh\n" + body + "\n")
        stub.chmod(0o755)
    return write


@pytest.mark.parametrize("body", ["exit 0", 'echo "Error response from daemon: No such container: $3" >&2; exit 0',
                                  'echo "Error response from daemon: No such container: $3" >&2; exit 1'])
def test_remove_container_confirmed(docker_stub, body):
    docker_stub(body)
    remove_container("cnvkit-l3-x")


@pytest.mark.parametrize("body, message", [
    ('echo "permission denied while trying to connect" >&2; exit 1', "exit 1: permission denied"),
    ("exec sleep 30", "timed out after 1 s"),
])
def test_remove_container_failure_raises(docker_stub, body, message):
    docker_stub(body)
    with pytest.raises(RuntimeError, match=message):
        remove_container("cnvkit-l3-x", timeout=1)


def test_remove_container_without_docker(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(RuntimeError, match="FileNotFoundError"):
        remove_container("cnvkit-l3-x")


def test_failed_container_removal_is_noted_on_the_timeout(docker_stub):
    docker_stub('echo "daemon unavailable" >&2; exit 1')
    with pytest.raises(subprocess.TimeoutExpired) as e:
        run_group(["sleep", "300"], 1, cleanup=lambda: remove_container("cnvkit-l3-x", timeout=5))
    assert e.value.cleanup_notes == ["cleanup unconfirmed: docker rm --force cnvkit-l3-x: exit 1: daemon unavailable"]


# ── interruption (round 3, S4): the group is killed, the cleanup runs, the original exception is re-raised ──────

def interrupt_once(monkeypatch, exc, ready=None):
    """Make the first Popen.communicate raise exc (after ready() is true, bounded), as a Ctrl-C while waiting would;
    later calls (e.g. the cleanup's own `docker rm`) run normally."""
    real, fired = subprocess.Popen.communicate, []

    def communicate(self, *a, **k):
        if fired:
            return real(self, *a, **k)
        fired.append(1)
        deadline = time.monotonic() + 10
        while ready is not None and not ready() and time.monotonic() < deadline:
            time.sleep(0.05)
        raise exc
    monkeypatch.setattr(subprocess.Popen, "communicate", communicate)


def exited(pid):
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] == "Z"
    except (FileNotFoundError, ProcessLookupError):
        return True


@pytest.mark.parametrize("exc", [KeyboardInterrupt(), RuntimeError("pytest internal error")])
def test_interrupt_kills_the_group_cleans_up_and_reraises(tmp_path, monkeypatch, exc):
    pidfile = tmp_path / "pid"
    interrupt_once(monkeypatch, exc, ready=lambda: pidfile.exists() and pidfile.read_text().strip())
    called = []
    with pytest.raises(type(exc)) as e:
        run_group(["bash", "-c", f'sleep 300 & echo "$!" > "{pidfile}"; wait'], 60,
                  cleanup=lambda: called.append(1), capture_output=True)
    assert e.value is exc and called == [1] and e.value.cleanup_notes == []
    grandchild = int(pidfile.read_text())
    deadline = time.monotonic() + 10
    while not exited(grandchild) and time.monotonic() < deadline:
        time.sleep(0.05)
    assert exited(grandchild)


def test_interrupt_cleanup_failure_is_noted_without_masking(monkeypatch):
    exc = KeyboardInterrupt()
    interrupt_once(monkeypatch, exc)

    def fail():
        raise RuntimeError("docker rm --force x: exit 1: permission denied")
    with pytest.raises(KeyboardInterrupt) as e:
        run_group(["sleep", "300"], 60, cleanup=fail)
    assert e.value is exc
    assert e.value.cleanup_notes == ["cleanup unconfirmed: docker rm --force x: exit 1: permission denied"]
    assert e.value.__notes__ == e.value.cleanup_notes


def test_interrupt_with_an_unkillable_leader_is_bounded(tmp_path, monkeypatch):
    def denied(*a):
        raise PermissionError(1, "Operation not permitted")
    pidfile = tmp_path / "pid"
    interrupt_once(monkeypatch, KeyboardInterrupt(), ready=lambda: pidfile.exists())
    monkeypatch.setattr(procutil.os, "killpg", denied)
    monkeypatch.setattr(subprocess.Popen, "kill", lambda self: denied())
    monkeypatch.setattr(procutil, "GRACE", 1)
    called = []
    start = time.monotonic()
    try:
        with pytest.raises(KeyboardInterrupt) as e:
            run_group(["bash", "-c", f'echo $$ > "{pidfile}"; exec sleep 60'], 60, cleanup=lambda: called.append(1))
        assert time.monotonic() - start < 20 and called == [1]
        notes = e.value.cleanup_notes
        assert any("killpg" in n for n in notes) and any("still running" in n for n in notes)
    finally:
        os.kill(int(pidfile.read_text()), signal.SIGKILL)


@pytest.mark.parametrize("rm, note", [("exit 0", None),
                                      ('echo "daemon unavailable" >&2; exit 1', "exit 1: daemon unavailable")])
def test_interrupted_l3_run_removes_its_container(tmp_path, docker_stub, monkeypatch, rm, note):
    # the L3 helper: a Ctrl-C while `docker run` is waiting kills the client and `docker rm --force`s its container
    import test_equiv_parfix as l3
    log = tmp_path / "docker.log"
    docker_stub(f'echo "$@" >> "{log}"\nif [ "$1" = run ]; then exec sleep 300; fi\n' + rm)
    exc = KeyboardInterrupt()
    interrupt_once(monkeypatch, exc, ready=lambda: log.exists())
    with pytest.raises(KeyboardInterrupt) as e:
        l3.in_image(tmp_path, "true")
    assert e.value is exc
    run_argv, rm_argv = [line.split() for line in log.read_text().splitlines()]
    name = run_argv[run_argv.index("--name") + 1]
    assert name.startswith("cnvkit-l3-") and rm_argv == ["rm", "--force", name]
    assert e.value.cleanup_notes == ([] if note is None else
                                     [f"cleanup unconfirmed: docker rm --force {name}: {note}"])


# ── round 4 (T4): a failed wait or output collection during the cleanup is a note, never the exception raised ───

def test_interrupt_with_a_failing_wait_reraises_the_interrupt(monkeypatch):
    exc = KeyboardInterrupt()
    interrupt_once(monkeypatch, exc)
    real_wait, failed = subprocess.Popen.wait, []

    def wait(self, *a, **k):
        if not failed:
            failed.append(self)
            raise OSError(10, "No child processes")
        return real_wait(self, *a, **k)
    monkeypatch.setattr(subprocess.Popen, "wait", wait)
    called = []
    with pytest.raises(KeyboardInterrupt) as e:
        run_group(["sleep", "300"], 60, cleanup=lambda: called.append(1))
    assert e.value is exc and called == [1]
    assert len(e.value.cleanup_notes) == 1 and "waiting for process" in e.value.cleanup_notes[0]
    assert "No child processes" in e.value.cleanup_notes[0] and e.value.__notes__ == e.value.cleanup_notes
    assert real_wait(failed[0], timeout=10) == -signal.SIGKILL       # the group was killed all the same


@pytest.mark.parametrize("error", [OSError(5, "Input/output error"), ValueError("I/O operation on closed file")])
def test_timeout_with_a_failing_output_collection_raises_the_timeout(monkeypatch, error):
    real, calls = subprocess.Popen.communicate, []

    def communicate(self, *a, **k):
        calls.append(self)
        if len(calls) == 2:          # the bounded collection after the kill
            raise error
        return real(self, *a, **k)
    monkeypatch.setattr(subprocess.Popen, "communicate", communicate)
    called = []
    with pytest.raises(subprocess.TimeoutExpired) as e:
        run_group(["bash", "-c", "sleep 300 & wait"], 1, cleanup=lambda: called.append(1), capture_output=True)
    assert called == [1] and len(calls) == 2
    assert len(e.value.cleanup_notes) == 1 and "collecting the output of process" in e.value.cleanup_notes[0]
    assert repr(error) in e.value.cleanup_notes[0] and e.value.__notes__ == e.value.cleanup_notes
    assert calls[0].wait(timeout=10) == -signal.SIGKILL


# ── round 5 (U2): the original TimeoutExpired is re-raised, with the output it captured ─────────────────────────

def recording_communicate(monkeypatch, fail_second=None):
    """Record the TimeoutExpired the first Popen.communicate raises; make the second (the collection after the kill)
    raise fail_second, if given."""
    real, raised, calls = subprocess.Popen.communicate, [], []

    def communicate(self, *a, **k):
        calls.append(self)
        if len(calls) == 2 and fail_second is not None:
            raise fail_second
        try:
            return real(self, *a, **k)
        except subprocess.TimeoutExpired as e:
            raised.append(e)
            raise
    monkeypatch.setattr(subprocess.Popen, "communicate", communicate)
    return raised, calls


SPEAK_THEN_HANG = "echo diagnostic; echo warning >&2; sleep 300 & wait"


def test_timeout_reraises_the_original_exception_with_the_complete_output(monkeypatch):
    raised, calls = recording_communicate(monkeypatch)
    with pytest.raises(subprocess.TimeoutExpired) as e:
        run_group(["bash", "-c", SPEAK_THEN_HANG], 1, capture_output=True, text=True)
    assert len(raised) == 1 and e.value is raised[0] and len(calls) == 2
    assert (e.value.output, e.value.stderr) == ("diagnostic\n", "warning\n")     # collected after the kill
    assert e.value.cmd == ["bash", "-c", SPEAK_THEN_HANG] and e.value.timeout == 1


@pytest.mark.parametrize("error", [OSError(5, "Input/output error"), ValueError("I/O operation on closed file")])
def test_timeout_keeps_the_partial_output_when_the_collection_fails(monkeypatch, error):
    # the reviewer's case: a command that printed `diagnostic` timed out with output=None
    raised, _ = recording_communicate(monkeypatch, fail_second=error)
    with pytest.raises(subprocess.TimeoutExpired) as e:
        run_group(["bash", "-c", SPEAK_THEN_HANG], 1, capture_output=True, text=True)
    assert e.value is raised[0]
    assert (e.value.output, e.value.stderr) == (b"diagnostic\n", b"warning\n")   # as communicate captured it
    assert "collecting the output of process" in e.value.cleanup_notes[0]
    assert "stderr: warning" in procutil.failure_text(e.value)


def test_timeout_keeps_the_partial_output_when_a_descendant_holds_the_pipes(tmp_path, monkeypatch):
    # a descendant in another session survives the group kill and holds the pipes: the bounded collection times out,
    # and the output read so far is kept
    monkeypatch.setattr(procutil, "GRACE", 1)
    pidfile = tmp_path / "pid"
    raised, _ = recording_communicate(monkeypatch)
    try:
        with pytest.raises(subprocess.TimeoutExpired) as e:
            run_group(["bash", "-c", f'echo diagnostic; setsid sleep 300 & echo $! > "{pidfile}"; wait'], 1,
                      capture_output=True, text=True)
        assert e.value is raised[0] and e.value.output == b"diagnostic\n"
        assert any("still running after the kill" in n for n in e.value.cleanup_notes)
    finally:
        if pidfile.exists() and pidfile.read_text().strip():
            try:
                os.kill(int(pidfile.read_text()), signal.SIGKILL)
            except ProcessLookupError:
                pass


def test_failure_text_keeps_stderr_and_cleanup_notes():
    e = subprocess.TimeoutExpired(["dx", "find"], 5, stderr=b"partial\n")
    e.cleanup_notes = ["killpg(1) failed: PermissionError()"]
    assert procutil.failure_text(e) == (f"{e}; {procutil.UNCONFIRMED}; stderr: partial; "
                                        "killpg(1) failed: PermissionError()")
    assert procutil.timed_out("dx find", 5, e) == (f"dx find timed out after 5 s; {procutil.UNCONFIRMED}; "
                                                  "stderr: partial; killpg(1) failed: PermissionError()")
    assert procutil.cleanup_unconfirmed(e)
    assert procutil.failure_text(OSError(2, "No such file")) == "[Errno 2] No such file"
    assert not procutil.cleanup_unconfirmed(subprocess.TimeoutExpired(["x"], 1))
    assert procutil.details(subprocess.CalledProcessError(1, ["x"], stderr="  boom \n")) == ["stderr: boom"]


# ── round 5 (U4): every named container is removed, an unconfirmed removal reported, the original failure kept ──

def test_removing_containers_removes_each_name(docker_stub, tmp_path):
    log = tmp_path / "rm.log"
    docker_stub(f'echo "$*" >> "{log}"; [ "$3" = gone ] && echo "Error: No such container: $3" >&2 && exit 1; exit 0')
    with procutil.removing_containers("a", "gone", "b", timeout=5):
        pass
    assert log.read_text().splitlines() == ["rm --force a", "rm --force gone", "rm --force b"]


def test_removing_containers_fails_a_passing_block_on_unconfirmed_removal(docker_stub, tmp_path):
    log = tmp_path / "rm.log"
    docker_stub(f'echo "$3" >> "{log}"; [ "$3" = a ] && echo "daemon unavailable" >&2 && exit 1; '
                '[ "$3" = b ] && exec sleep 30; exit 0')
    with pytest.raises(RuntimeError) as e:
        with procutil.removing_containers("a", "b", "c", timeout=1):
            pass
    assert log.read_text().splitlines() == ["a", "b", "c"]               # every name attempted
    msg = str(e.value)
    assert msg.startswith("container removal unconfirmed: docker rm --force a: exit 1: daemon unavailable")
    assert "docker rm --force b: timed out after 1 s" in msg and "--force c" not in msg


@pytest.mark.parametrize("body", ['echo "daemon unavailable" >&2; exit 1', "exec sleep 30"])
def test_removing_containers_keeps_the_original_failure(docker_stub, body):
    docker_stub(body)
    original = AssertionError("the test's own failure")
    with pytest.raises(AssertionError) as e:
        with procutil.removing_containers("a", timeout=1):
            raise original
    assert e.value is original
    assert len(e.value.__notes__) == 1 and e.value.__notes__[0].startswith("container removal unconfirmed: "
                                                                           "docker rm --force a: ")


def test_removing_containers_without_docker(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(RuntimeError, match="container removal unconfirmed: docker rm --force a: FileNotFoundError"):
        with procutil.removing_containers("a"):
            pass


# ── round 6 (V2): no exception from a removal escapes, every name is removed, the original failure is kept ─────────

INVALID_UTF8 = r"printf 'daemon \377\376 unavailable\n' >&2"


def test_text_output_that_is_not_utf8_never_raises():
    p = run_group(["sh", "-c", r"printf 'out \377\n'; printf 'err \376\n' >&2"], 30, capture_output=True, text=True)
    assert (p.stdout, p.stderr) == ("out \ufffd\n", "err \ufffd\n")
    with pytest.raises(subprocess.TimeoutExpired) as e:   # the output collected after a timeout kill, too
        run_group(["sh", "-c", r"printf 'out \377\n'; exec sleep 300"], 1, capture_output=True, text=True)
    assert e.value.output == "out \ufffd\n"


def test_an_explicit_errors_argument_is_kept():
    with pytest.raises(UnicodeDecodeError):
        run_group(["sh", "-c", r"printf '\377'"], 30, capture_output=True, text=True, errors="strict")
    assert run_group(["sh", "-c", r"printf '\377'"], 30, capture_output=True).stdout == b"\xff"   # bytes untouched


@pytest.mark.parametrize("rm_exit", [1, 0])
def test_remove_container_with_invalid_utf8_from_docker(docker_stub, rm_exit):
    docker_stub(f"{INVALID_UTF8}; exit {rm_exit}")
    if rm_exit:
        with pytest.raises(RuntimeError, match="exit 1: daemon \ufffd\ufffd unavailable"):
            remove_container("cnvkit-l3-x", timeout=5)
    else:
        remove_container("cnvkit-l3-x", timeout=5)
    docker_stub(r"printf 'Error: No such container: \377\n' >&2; exit 1")
    remove_container("cnvkit-l3-x", timeout=5)       # absent, confirmed despite the invalid byte


def test_remove_container_raises_only_runtime_error(monkeypatch):
    def broken(*a, **k):
        raise ValueError("unexpected")
    monkeypatch.setattr(procutil, "run_group", broken)
    with pytest.raises(RuntimeError, match=r"docker rm --force x: ValueError\('unexpected'\)"):
        remove_container("x")


def test_remove_containers_reports_any_exception_and_continues(monkeypatch):
    attempted = []

    def remove(name, timeout):
        attempted.append(name)
        if name == "a":
            raise UnicodeDecodeError("utf-8", b"\xff", 0, 1, "invalid start byte")
        if name == "b":
            raise RuntimeError("docker rm --force b: exit 1: unavailable")
    monkeypatch.setattr(procutil, "remove_container", remove)
    problems = procutil.remove_containers(["a", "b", "c"])
    assert attempted == ["a", "b", "c"]
    assert problems[0].startswith("docker rm --force a: UnicodeDecodeError(") and problems[1:] == [
        "docker rm --force b: exit 1: unavailable"]
    original = AssertionError("the test's own failure")
    with pytest.raises(AssertionError) as e:
        with procutil.removing_containers("a", "b", "c"):
            raise original
    assert e.value is original and len(e.value.__notes__) == 2
    with pytest.raises(RuntimeError, match="container removal unconfirmed: docker rm --force a: UnicodeDecodeError"):
        with procutil.removing_containers("a", "c"):
            pass


def test_removing_containers_with_invalid_utf8_from_docker(docker_stub, tmp_path):
    log = tmp_path / "rm.log"
    docker_stub(f'echo "$3" >> "{log}"; [ "$3" = a ] && {INVALID_UTF8} && exit 1; exit 0')
    original = AssertionError("the test's own failure")
    with pytest.raises(AssertionError) as e:
        with procutil.removing_containers("a", "b", timeout=5):
            raise original
    assert e.value is original and log.read_text().splitlines() == ["a", "b"]
    assert e.value.__notes__ == ["container removal unconfirmed: docker rm --force a: exit 1: daemon \ufffd\ufffd "
                                 "unavailable"]


@pytest.mark.parametrize("where", ["timeout", "interrupt"])
def test_an_unexpected_kill_failure_is_a_note(monkeypatch, where):
    # kill_group raising something other than OSError is noted; the original exception is still raised, and the
    # cleanup still runs
    real_kill = procutil.kill_group

    def kill(proc):
        real_kill(proc)
        raise ValueError("kill bookkeeping failed")
    monkeypatch.setattr(procutil, "kill_group", kill)
    exc = KeyboardInterrupt()
    if where == "interrupt":
        interrupt_once(monkeypatch, exc)
    called = []
    expected = subprocess.TimeoutExpired if where == "timeout" else KeyboardInterrupt
    with pytest.raises(expected) as e:
        run_group(["sleep", "300"], 1 if where == "timeout" else 60, cleanup=lambda: called.append(1))
    assert called == [1] and (where == "timeout" or e.value is exc)
    assert len(e.value.cleanup_notes) == 1
    assert re.fullmatch(r"killing process group \d+ failed: ValueError\('kill bookkeeping failed'\)",
                        e.value.cleanup_notes[0])


def test_a_non_runtime_cleanup_failure_is_a_note():
    def fail():
        raise UnicodeDecodeError("utf-8", b"\xff", 0, 1, "invalid start byte")
    with pytest.raises(subprocess.TimeoutExpired) as e:
        run_group(["sleep", "300"], 1, cleanup=fail)
    assert len(e.value.cleanup_notes) == 1 and e.value.cleanup_notes[0].startswith(
        "cleanup unconfirmed: UnicodeDecodeError(")


# round 7 (W3): a pipe that fails to close after the kill is recorded on the original exception, not discarded; it is
# not unconfirmed cleanup (the process was killed), so it cannot turn a skip into a failure

class _ClosingFails:
    """A pipe whose close() raises once `armed` (the cleanup ran): only run_group's own close after the kill fails."""

    def __init__(self, f, armed):
        self._f, self._armed = f, armed

    def close(self):
        self._f.close()
        if self._armed:
            raise OSError(5, "Input/output error")

    def __getattr__(self, name):
        return getattr(self._f, name)


@pytest.mark.parametrize("where", ["timeout", "interrupt"])
def test_a_pipe_close_failure_is_a_note_on_the_original_exception(monkeypatch, where):
    exc = KeyboardInterrupt()
    if where == "interrupt":
        interrupt_once(monkeypatch, exc)
    armed = []
    real = subprocess.Popen

    class Popen(real):
        def __init__(self, *a, **k):
            super().__init__(*a, **k)
            self.stdout = _ClosingFails(self.stdout, armed)
    monkeypatch.setattr(procutil.subprocess, "Popen", Popen)
    expected = subprocess.TimeoutExpired if where == "timeout" else KeyboardInterrupt
    with pytest.raises(expected) as e:
        run_group(["sleep", "300"], 1 if where == "timeout" else 60, cleanup=lambda: armed.append(1),
                  stdout=subprocess.PIPE)
    assert armed == [1] and (where == "timeout" or e.value is exc)
    assert e.value.cleanup_notes == [] and not procutil.cleanup_unconfirmed(e.value)
    assert len(e.value.close_notes) == 1
    assert re.fullmatch(r"closing a pipe of process \d+ failed: OSError\(5, 'Input/output error'\)",
                        e.value.close_notes[0])
    assert e.value.__notes__ == e.value.close_notes
    text = procutil.failure_text(e.value)
    assert e.value.close_notes[0] in text and procutil.UNCONFIRMED not in text


def test_pipes_that_close_leave_no_close_notes():
    with pytest.raises(subprocess.TimeoutExpired) as e:
        run_group(["sleep", "300"], 1, capture_output=True)
    assert e.value.close_notes == [] and e.value.cleanup_notes == [] and not hasattr(e.value, "__notes__")
