"""L1: the provenance guards of tests/dx/validation_jobs.py (the DNAnexus validation driver), with dx mocked."""
import argparse
import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO = Path(__file__).resolve().parents[2]
JOB_COLS = ["key", "run", "sample", "sample_id", "extra", "job", "applet", "git_commit"]


@pytest.fixture
def vj(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("validation_jobs", REPO / "tests" / "dx" / "validation_jobs.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    monkeypatch.setattr(m, "RUNS", tmp_path / "runs")
    monkeypatch.setattr(m, "DL", tmp_path / "dl")
    (tmp_path / "runs").mkdir()
    (tmp_path / "runs" / "applet.txt").write_text("applet-NEW\tabc1234\n")
    m.write_tsv(tmp_path / "runs" / "inputs.tsv",
                [dict(run="0797_V2", sample="S1", bam="p:file-B", bai="p:file-I", pon_id="D_GC__shared__gcclip__parx",
                      cn_reference="p:file-P")], ["run", "sample", "bam", "bai", "pon_id", "cn_reference"])
    calls = []

    def fake_dx(*args, check=True):
        calls.append(list(args))
        return "job-NEW\n" if args[0] == "run" else ""
    monkeypatch.setattr(m, "dx", fake_dx)
    m.calls = calls
    return m


def record(vj, applet):
    vj.write_tsv(vj.jobs_path("E3"), [dict(key="0797_V2/S1", run="0797_V2", sample="S1", sample_id="S1",
                                           extra="{}", job="job-OLD", applet=applet, git_commit="old")], JOB_COLS)


def test_build_refuses_untracked_files(vj, monkeypatch):
    real_run = subprocess.run

    def fake_run(cmd, **kw):
        if cmd[:2] == ["git", "status"]:
            assert "--untracked-files=no" not in cmd
            return subprocess.CompletedProcess(cmd, 0, stdout="?? resources/usr/local/bin/new.py\n", stderr="")
        return real_run(cmd, **kw)
    monkeypatch.setattr(vj.subprocess, "run", fake_run)
    with pytest.raises(SystemExit, match="commit the working tree first"):
        vj.build(argparse.Namespace())
    assert vj.calls == []                         # nothing built


def test_submit_refuses_jobs_from_another_applet(vj):
    record(vj, "applet-OLD")
    with pytest.raises(SystemExit, match="not run with applet-NEW"):
        vj.submit(argparse.Namespace(set="E3", reset=False))
    assert not [c for c in vj.calls if c[0] == "run"]


def test_submit_resumes_jobs_from_the_same_applet(vj):
    record(vj, "applet-NEW")
    vj.submit(argparse.Namespace(set="E3", reset=False))
    assert not [c for c in vj.calls if c[0] == "run"]
    assert [r["job"] for r in vj.read_tsv(vj.jobs_path("E3"))] == ["job-OLD"]


def test_submit_reset_starts_again_with_the_current_applet(vj):
    record(vj, "applet-OLD")
    vj.submit(argparse.Namespace(set="E3", reset=True))
    runs = [c for c in vj.calls if c[0] == "run"]
    assert len(runs) == 1 and runs[0][1] == "applet-NEW"
    rows = vj.read_tsv(vj.jobs_path("E3"))
    assert [(r["job"], r["applet"], r["git_commit"]) for r in rows] == [("job-NEW", "applet-NEW", "abc1234")]


def test_check_fails_a_job_run_by_another_executable(vj, monkeypatch):
    record(vj, "applet-NEW")
    desc = {"id": "job-OLD", "state": "failed", "executable": "applet-OTHER", "project": "p", "output": None}
    monkeypatch.setattr(vj, "prefetch", lambda which, rows: {"job-OLD": (desc, "")})
    monkeypatch.setattr(vj.run_local, "treated_male_log", lambda data, mask: {})   # no validation data in CI
    assert vj.check(argparse.Namespace(set="E3")) == 1
    report = vj.read_tsv(vj.RUNS / "E3_report.tsv")
    assert any(r["check"] == "applet" and r["result"] == "FAIL" and "applet-OTHER" in r["detail"] for r in report)
    assert json.loads(vj.read_tsv(vj.jobs_path("E3"))[0]["extra"]) == {}


# ── dx failure paths (CodeRabbit review) ─────────────────────────────────────

def load_fresh():
    spec = importlib.util.spec_from_file_location("validation_jobs_raw", REPO / "tests" / "dx" / "validation_jobs.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_dx_non_zero_exit_raises_with_stderr(monkeypatch):
    m = load_fresh()
    monkeypatch.setattr(m.subprocess, "run", lambda cmd, **kw: subprocess.CompletedProcess(
        cmd, 3, stdout="", stderr="dxpy.exceptions.PermissionDenied: VIEW permission required in project-x"))
    with pytest.raises(RuntimeError, match="PermissionDenied: VIEW permission required"):
        m.dx("describe", "file-x")


def test_dx_missing_binary_is_a_clear_error(monkeypatch):
    m = load_fresh()

    def missing(cmd, **kw):
        raise FileNotFoundError(2, "No such file or directory", "dx")
    monkeypatch.setattr(m.subprocess, "run", missing)
    with pytest.raises(RuntimeError, match="dx CLI not found"):
        m.dx("describe", "file-x")


def test_dx_timeouts(monkeypatch):
    m = load_fresh()
    seen = []

    def hang(cmd, **kw):
        seen.append((cmd[1], kw["timeout"]))
        raise subprocess.TimeoutExpired(cmd, kw["timeout"])
    monkeypatch.setattr(m.subprocess, "run", hang)
    for args in (("describe", "job-x"), ("download", "file-x")):
        with pytest.raises(RuntimeError, match="timed out"):
            m.dx(*args)
    assert seen == [("describe", m.DX_TIMEOUT), ("download", m.LONG_TIMEOUT)]


def test_failed_build_records_no_applet(vj, monkeypatch):
    def fake_dx(*args, check=True, timeout=None):
        if args[0] == "build":
            raise RuntimeError("dx build /repo: dxpy.app_builder.AppBuilderException: invalid dxapp.json")
        return ""
    monkeypatch.setattr(vj, "dx", fake_dx)
    monkeypatch.setattr(vj.subprocess, "run", lambda cmd, **kw: subprocess.CompletedProcess(
        cmd, 0, stdout="" if cmd[:2] == ["git", "status"] else "abc1234\n", stderr=""))
    before = (vj.RUNS / "applet.txt").read_text()
    with pytest.raises(SystemExit, match="invalid dxapp.json"):
        vj.build(argparse.Namespace())
    assert (vj.RUNS / "applet.txt").read_text() == before


@pytest.mark.parametrize("state, log, cached", [("done", "[sex] applied ...\n", True), ("done", "", False),
                                                ("running", "partial\n", False)])
def test_job_log_cached_only_when_finished_and_not_empty(vj, tmp_path, monkeypatch, state, log, cached):
    monkeypatch.setattr(vj, "dxjson", lambda *a: {"id": "job-1", "state": state, "project": "p", "output": {}})
    monkeypatch.setattr(vj, "dx", lambda *a, **kw: log)
    d, text = vj.fetch_job("job-1", tmp_path / "job")
    assert text == log
    assert (tmp_path / "job" / "job.log").exists() is cached
