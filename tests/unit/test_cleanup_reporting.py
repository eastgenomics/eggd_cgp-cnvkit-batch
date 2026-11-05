"""L1 (round 5, U3): the callers of procutil.run_group keep its cleanup notes and diagnostics.

The Docker check fails a test, rather than skipping it, when its own cleanup is unconfirmed; the dx and git helpers
of the DNAnexus driver keep stderr and the notes in their messages. Nothing here reaches the real `dx` or `docker`:
both are failing stubs first on PATH, and every command is mocked.
"""
import importlib.util
import os
import subprocess
from pathlib import Path

import pytest

import procutil

pytestmark = pytest.mark.unit

REPO = Path(__file__).resolve().parents[2]
NOTE = "process 7 (or a descendant holding its output) still running after the kill and 10 s"


@pytest.fixture(autouse=True)
def no_live_tools(tmp_path_factory, monkeypatch):
    """`dx` and `docker` first on PATH that fail: every call here is mocked, and must never reach the real CLIs."""
    bin_dir = tmp_path_factory.mktemp("no_live_tools")
    for name in ("dx", "docker"):
        stub = bin_dir / name
        stub.write_text(f'#!/bin/sh\necho "unit test reached the {name} CLI: $*" >&2\nexit 99\n')
        stub.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def timeout(notes=(), stderr=None):
    def raise_timeout(cmd, t, **kwargs):
        e = subprocess.TimeoutExpired(cmd, t, stderr=stderr)
        e.cleanup_notes = list(notes)
        raise e
    return raise_timeout


@pytest.fixture
def conftest(monkeypatch):
    return load(REPO / "tests" / "conftest.py", "conftest_reporting_under_test")


class Item:
    def get_closest_marker(self, name):
        return name == "docker" or None


def test_docker_check_with_unconfirmed_cleanup_fails(conftest, monkeypatch):
    monkeypatch.setattr(conftest, "run_group", timeout([NOTE]))
    action, reason = conftest.docker_check.__wrapped__("cgp-cnvkit:1.0.0")
    assert action == "fail"
    assert reason.startswith("docker image inspect cgp-cnvkit:1.0.0 timed out after 30 s; ")
    assert procutil.UNCONFIRMED in reason and NOTE in reason
    monkeypatch.setattr(conftest, "docker_check", lambda: (action, reason))
    with pytest.raises(pytest.fail.Exception, match="cleanup unconfirmed"):
        conftest.pytest_runtest_setup(Item())


def test_docker_check_with_a_confirmed_timeout_skips(conftest, monkeypatch):
    monkeypatch.setattr(conftest, "run_group", timeout(stderr=b"slow daemon\n"))
    action, reason = conftest.docker_check.__wrapped__("cgp-cnvkit:1.0.0")
    assert action == "skip" and reason.endswith("timed out after 30 s; stderr: slow daemon")
    monkeypatch.setattr(conftest, "docker_check", lambda: (action, reason))
    with pytest.raises(pytest.skip.Exception, match="timed out"):
        conftest.pytest_runtest_setup(Item())


def test_docker_check_with_only_a_pipe_close_note_still_skips(conftest, monkeypatch):
    # round 7 (W3): a pipe that failed to close after a confirmed kill is reported, but is not unconfirmed cleanup
    note = "closing a pipe of process 7 failed: OSError(5, 'Input/output error')"

    def raise_timeout(cmd, t, **kwargs):
        e = subprocess.TimeoutExpired(cmd, t)
        e.cleanup_notes, e.close_notes = [], [note]
        raise e
    monkeypatch.setattr(conftest, "run_group", raise_timeout)
    action, reason = conftest.docker_check.__wrapped__("cgp-cnvkit:1.0.0")
    assert action == "skip" and reason == f"docker image inspect cgp-cnvkit:1.0.0 timed out after 30 s; {note}"
    assert procutil.UNCONFIRMED not in reason


def test_docker_check_with_an_oserror_and_unconfirmed_cleanup_fails(conftest, monkeypatch):
    # PR #6 (CodeRabbit): run_group notes an unconfirmed cleanup on any exception while waiting, OSError included
    def raise_oserror(cmd, t, **kwargs):
        e = OSError(5, "Input/output error")
        e.cleanup_notes = [NOTE]
        raise e
    monkeypatch.setattr(conftest, "run_group", raise_oserror)
    action, reason = conftest.docker_check.__wrapped__("cgp-cnvkit:1.0.0")
    assert action == "fail" and reason.startswith("docker not runnable: ") and procutil.UNCONFIRMED in reason


def test_docker_check_with_a_plain_oserror_skips(conftest, monkeypatch):
    def raise_oserror(cmd, t, **kwargs):
        raise PermissionError(13, "Permission denied")
    monkeypatch.setattr(conftest, "run_group", raise_oserror)
    action, reason = conftest.docker_check.__wrapped__("cgp-cnvkit:1.0.0")
    assert action == "skip" and reason.startswith("docker not runnable: ") and procutil.UNCONFIRMED not in reason


def test_docker_check_missing_image_skips_with_stderr(conftest, monkeypatch):
    monkeypatch.setattr(conftest, "run_group", lambda cmd, t, **kw: subprocess.CompletedProcess(
        cmd, 1, stdout="[]", stderr="Error: No such image: cgp-cnvkit:1.0.0\n"))
    assert conftest.docker_check.__wrapped__("cgp-cnvkit:1.0.0") == (
        "skip", "cgp-cnvkit:1.0.0 not available (Error: No such image: cgp-cnvkit:1.0.0)")


def test_docker_check_usable(conftest, monkeypatch):
    monkeypatch.setattr(conftest, "run_group", lambda cmd, t, **kw: subprocess.CompletedProcess(cmd, 0, "", ""))
    assert conftest.docker_check.__wrapped__("cgp-cnvkit:1.0.0") == (None, None)
    monkeypatch.setattr(conftest, "docker_check", lambda: (None, None))
    conftest.pytest_runtest_setup(Item())       # neither skip nor fail


def test_docker_check_reaches_only_the_stub(conftest):
    # the failing stub is what runs (exit 99): a skip with its stderr, never the real docker
    action, reason = conftest.docker_check.__wrapped__("cgp-cnvkit:1.0.0")
    assert action == "skip" and "unit test reached the docker CLI" in reason


def test_docker_check_with_invalid_utf8_skips(conftest, tmp_path, monkeypatch):
    # round 6 (V2): a failing stub first on PATH whose stderr is not valid UTF-8: a skip with the stderr replaced,
    # never UnicodeDecodeError out of the check
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "docker").write_text("#!/bin/sh\nprintf 'Error: no such image \\377\\n' >&2\nexit 1\n")
    (bin_dir / "docker").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    assert conftest.docker_check.__wrapped__("cgp-cnvkit:1.0.0") == (
        "skip", "cgp-cnvkit:1.0.0 not available (Error: no such image \ufffd)")


@pytest.fixture
def vj(monkeypatch):
    return load(REPO / "tests" / "dx" / "validation_jobs.py", "validation_jobs_reporting_under_test")


def test_dx_timeout_keeps_the_cleanup_notes(vj, monkeypatch):
    monkeypatch.setattr(vj, "run_group", timeout([NOTE], stderr="partial progress\n"))
    with pytest.raises(RuntimeError) as e:
        vj.dx("describe", "file-X")
    msg = str(e.value)
    assert msg.startswith(f"dx describe file-X timed out after {vj.DX_TIMEOUT} s; {procutil.UNCONFIRMED}; ")
    assert "stderr: partial progress" in msg and NOTE in msg


def test_dx_reaches_only_the_stub(vj):
    with pytest.raises(RuntimeError, match="unit test reached the dx CLI: describe file-X"):
        vj.dx("describe", "file-X")


@pytest.mark.parametrize("failure", ["timeout", "nonzero"])
def test_git_failure_keeps_stderr_and_notes(vj, monkeypatch, failure):
    if failure == "timeout":
        monkeypatch.setattr(vj, "run_group", timeout([NOTE]))
    else:
        def nonzero(cmd, t, **kw):
            assert kw["check"] is True and t == vj.GIT_TIMEOUT
            raise subprocess.CalledProcessError(128, cmd, stderr="fatal: not a git repository\n")
        monkeypatch.setattr(vj, "run_group", nonzero)
    with pytest.raises(RuntimeError) as e:
        vj.git("rev-parse", "--short", "HEAD")
    msg = str(e.value)
    if failure == "timeout":
        assert msg.startswith("git rev-parse --short HEAD timed out after 60 s; ") and NOTE in msg
        assert procutil.UNCONFIRMED in msg
    else:
        assert msg.startswith("git rev-parse --short HEAD: Command ") and "stderr: fatal: not a git repository" in msg
