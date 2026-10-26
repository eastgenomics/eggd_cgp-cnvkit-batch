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
    assert vj.check(argparse.Namespace(set="E3")) == 1
    report = vj.read_tsv(vj.RUNS / "E3_report.tsv")
    assert any(r["check"] == "applet" and r["result"] == "FAIL" and "applet-OTHER" in r["detail"] for r in report)
    assert json.loads(vj.read_tsv(vj.jobs_path("E3"))[0]["extra"]) == {}
