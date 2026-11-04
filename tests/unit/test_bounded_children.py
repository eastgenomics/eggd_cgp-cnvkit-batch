"""L1 (round 4, T5): every bounded helper command kills its whole process group on timeout.

subprocess.run(..., timeout=...) kills only the process it started, so a child of a timed-out CLI (or of a wrapper
script standing in for it) survives. Each helper below runs a stub `docker` / `dx` that starts a long-running child
and waits; on timeout the child must be gone too.
"""
import importlib.util
import io
import os
import signal
import subprocess
import time
from pathlib import Path

import pytest

import procutil
import run_local

pytestmark = pytest.mark.unit

REPO = Path(__file__).resolve().parents[2]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def exited(pid):
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] == "Z"
    except (FileNotFoundError, ProcessLookupError):
        return True


@pytest.fixture
def hanging(tmp_path, monkeypatch):
    """Put stub `docker` and `dx` on PATH that start `sleep 300` (sharing their output) and wait for it. Yields a
    function returning the child's pid; the child is killed at teardown if it survived."""
    bin_dir, pidfile = tmp_path / "bin", tmp_path / "child.pid"
    bin_dir.mkdir()
    for name in ("docker", "dx"):
        stub = bin_dir / name
        stub.write_text(f'#!/bin/sh\nsleep 300 &\necho "$!" > "{pidfile}"\nwait\n')
        stub.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")

    def child():
        deadline = time.monotonic() + 10
        while not (pidfile.exists() and pidfile.read_text().strip()) and time.monotonic() < deadline:
            time.sleep(0.05)
        return int(pidfile.read_text())

    def gone(pid):
        deadline = time.monotonic() + 10
        while not exited(pid) and time.monotonic() < deadline:
            time.sleep(0.05)
        return exited(pid)
    yield child, gone
    if pidfile.exists() and pidfile.read_text().strip() and not exited(int(pidfile.read_text())):
        os.kill(int(pidfile.read_text()), signal.SIGKILL)


def test_remove_container_kills_the_child(hanging):
    child, gone = hanging
    start = time.monotonic()
    with pytest.raises(RuntimeError, match="timed out after 1 s"):
        procutil.remove_container("cnvkit-l3-x", timeout=1)
    assert time.monotonic() - start < 1 + procutil.GRACE
    assert gone(child())


def test_docker_check_kills_the_child(hanging, monkeypatch):
    child, gone = hanging
    conftest = load(REPO / "tests" / "conftest.py", "conftest_under_test")
    monkeypatch.setattr(conftest, "DOCKER_CHECK_TIMEOUT", 1)
    assert conftest.docker_check.__wrapped__("cgp-cnvkit:1.0.0") == (
        "skip", "docker image inspect cgp-cnvkit:1.0.0 timed out after 1 s")     # the kill was confirmed
    assert gone(child())


def test_harness_stop_container_kills_the_child(hanging, monkeypatch):
    child, gone = hanging
    monkeypatch.setattr(run_local, "CLEANUP_TIMEOUT", 1)
    log = io.StringIO()
    assert not run_local.stop_container("cnvkit-harness-x", log)
    assert "Container cleanup not confirmed" in log.getvalue()
    assert gone(child())


def test_harness_bed_download_kills_the_child(hanging, tmp_path, monkeypatch):
    child, gone = hanging
    monkeypatch.setattr(run_local, "HARNESS", tmp_path / "harness")
    monkeypatch.setattr(run_local, "BED_TIMEOUT", 1)
    with pytest.raises(subprocess.TimeoutExpired):
        run_local.fetch_bed()
    assert gone(child())
    assert not list((tmp_path / "harness" / "inputs").iterdir())


def test_validation_jobs_dx_kills_the_child(hanging, monkeypatch):
    child, gone = hanging
    vj = load(REPO / "tests" / "dx" / "validation_jobs.py", "validation_jobs_under_test")
    monkeypatch.setattr(vj, "DX_TIMEOUT", 1)
    with pytest.raises(RuntimeError, match="timed out after 1 s"):
        vj.dx("describe", "file-X")
    assert gone(child())


def test_fetch_validation_data_download_kills_the_child(hanging, tmp_path, monkeypatch):
    child, gone = hanging
    fv = load(REPO / "tests" / "validation" / "fetch_validation_data.py", "fetch_validation_data_under_test")
    monkeypatch.setattr(fv, "DOWNLOADS", tmp_path / "downloads")
    monkeypatch.setattr(fv, "DOWNLOAD_TIMEOUT", 1)
    with pytest.raises(subprocess.TimeoutExpired):
        fv.download("x.tar.gz", "file-X")
    assert gone(child())
