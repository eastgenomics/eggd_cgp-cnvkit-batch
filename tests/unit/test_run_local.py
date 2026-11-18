"""L1: shared BED and timeout cleanup; no live DNAnexus or Docker calls."""
import argparse
import concurrent.futures as cf
import csv
import io
import signal
import subprocess
from pathlib import Path

import pytest

import run_local

pytestmark = pytest.mark.unit


def test_shared_bed_is_published_once(tmp_path, monkeypatch):
    monkeypatch.setattr(run_local, "HARNESS", tmp_path)
    commands = []

    def download(cmd, **kwargs):
        commands.append(cmd)
        assert kwargs["timeout"] == run_local.BED_TIMEOUT
        assert not (tmp_path / "inputs" / "panel.bed").exists()
        Path(cmd[-1]).write_text("complete\n")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(run_local.subprocess, "run", download)
    with cf.ThreadPoolExecutor(max_workers=8) as ex:
        beds = list(ex.map(lambda _: run_local.fetch_bed(), range(8)))
    assert len(commands) == 1 and all(p == beds[0] for p in beds)
    assert beds[0].read_text() == "complete\n"
    assert not list(beds[0].parent.glob("*.part"))


@pytest.mark.parametrize("failure", ["timeout", "nonzero", "missing"])
def test_shared_bed_failure_removes_partial(tmp_path, monkeypatch, failure):
    monkeypatch.setattr(run_local, "HARNESS", tmp_path)

    def download(cmd, **kwargs):
        assert kwargs["timeout"] == run_local.BED_TIMEOUT
        Path(cmd[-1]).write_text("partial")
        if failure == "timeout":
            raise subprocess.TimeoutExpired(cmd, kwargs["timeout"])
        if failure == "missing":
            raise FileNotFoundError("dx")
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(run_local.subprocess, "run", download)
    with pytest.raises((subprocess.SubprocessError, OSError)):
        run_local.fetch_bed()
    assert not (tmp_path / "inputs" / "panel.bed").exists()
    assert not list((tmp_path / "inputs").glob("*.part"))


@pytest.mark.parametrize("failure", ["timeout", "nonzero"])
def test_shared_bed_failure_reports_all_without_workers(tmp_path, monkeypatch, failure):
    monkeypatch.setattr(run_local, "HARNESS", tmp_path)
    rows = [{"run": "R1", "sample": "S1"}, {"run": "R1", "sample": "S2"}]
    monkeypatch.setattr(run_local, "local_samples", lambda: rows)
    monkeypatch.setattr(run_local, "treated_male_log", lambda *args: {})
    wd = tmp_path / "test" / "R1" / "S1"
    wd.mkdir(parents=True)
    marker = wd / ".cleanup-unconfirmed"
    marker.write_text("cnvkit-harness-old\n")
    log = wd / "pipeline.log"
    log.write_text("previous diagnostics\n")
    calls = []

    def download(cmd, **kwargs):
        calls.append(cmd)
        Path(cmd[-1]).write_text("partial")
        if failure == "timeout":
            raise subprocess.TimeoutExpired(cmd, kwargs["timeout"])
        raise subprocess.CalledProcessError(1, cmd)

    def must_not_run(*args):
        pytest.fail("a sample worker started after shared BED failure")

    monkeypatch.setattr(run_local.subprocess, "run", download)
    monkeypatch.setattr(run_local, "run_sample", must_not_run)
    assert run_local.main(["--data", "nonumi", "--start-from", "bam", "--label", "test"]) == 1
    assert len(calls) == 1
    assert not list((tmp_path / "inputs").glob("*.part"))
    assert marker.read_text() == "cnvkit-harness-old\n"
    assert log.read_text() == "previous diagnostics\n"
    with (tmp_path / "test" / "report.tsv").open() as fh:
        report = list(csv.DictReader(fh, delimiter="\t"))
    assert [(r["sample"], r["result"]) for r in report] == [("S1", "FAIL"), ("S2", "FAIL")]
    assert all("shared BED download failed" in r["detail"] for r in report)
    assert "cleanup unconfirmed" in report[0]["detail"] and str(marker) in report[0]["detail"]


def test_bed_preflight_runs_before_workers(tmp_path, monkeypatch):
    monkeypatch.setattr(run_local, "HARNESS", tmp_path)
    rows = [{"run": "R1", "sample": "S1"}, {"run": "R1", "sample": "S2"}]
    monkeypatch.setattr(run_local, "local_samples", lambda: rows)
    monkeypatch.setattr(run_local, "treated_male_log", lambda *args: {})
    calls = []
    bed = tmp_path / "panel.bed"

    def preflight():
        calls.append("bed")
        return bed

    def sample(row, args, out):
        assert calls[0] == "bed" and args.bed == bed
        return 0

    monkeypatch.setattr(run_local, "fetch_bed", preflight)
    monkeypatch.setattr(run_local, "run_sample", sample)
    monkeypatch.setattr(run_local, "check_sample", lambda row, *args: [(row["sample"], "pipeline", "PASS", "")])
    assert run_local.main(["--data", "nonumi", "--start-from", "bam", "--label", "test"]) == 0
    assert calls == ["bed"]


def coverage_sample(tmp_path, monkeypatch):
    monkeypatch.setattr(run_local, "DATA", tmp_path / "data")
    pon = tmp_path / "data" / "nonumi" / "pons" / "pon.cnn"
    pon.parent.mkdir(parents=True)
    pon.write_text("reference\n")
    cnn = tmp_path / "data" / "nonumi" / "cnn" / "R1" / "S1.targetcoverage.cnn"
    cnn.parent.mkdir(parents=True)
    cnn.write_text("coverage\n")
    return ({"run": "R1", "sample": "S1", "cn_reference": "pon"},
            argparse.Namespace(data="nonumi", start_from="coverage", nproc=2, mask=True))


def test_normal_completion_uses_unique_container_names(tmp_path, monkeypatch):
    row, args = coverage_sample(tmp_path, monkeypatch)
    names = []

    class Done:
        def wait(self, timeout):
            assert timeout == run_local.PIPELINE_TIMEOUT
            return 0

    def popen(cmd, **kwargs):
        assert kwargs["start_new_session"] is True
        assert '--name "$HARNESS_CONTAINER_NAME"' in cmd[2]
        names.append(kwargs["env"]["HARNESS_CONTAINER_NAME"])
        return Done()

    monkeypatch.setattr(run_local.subprocess, "Popen", popen)
    assert run_local.run_sample(row, args, tmp_path / "first") == 0
    assert run_local.run_sample(row, args, tmp_path / "second") == 0
    assert len(set(names)) == 2 and all(n.startswith("cnvkit-harness-") for n in names)


def test_pipeline_timeout_stops_process_group_and_exact_container(tmp_path, monkeypatch):
    row, args = coverage_sample(tmp_path, monkeypatch)
    events, seen = [], {}

    class FakeProcess:
        pid = 12345

        def wait(self, timeout):
            events.append(("wait", timeout))
            if timeout == run_local.PIPELINE_TIMEOUT:
                raise subprocess.TimeoutExpired("bash", timeout)
            return -signal.SIGKILL

    def popen(cmd, **kwargs):
        assert kwargs["start_new_session"] is True
        seen["container"] = kwargs["env"]["HARNESS_CONTAINER_NAME"]
        return FakeProcess()

    def docker(cmd, **kwargs):
        events.append(("docker", cmd, kwargs["timeout"]))
        return subprocess.CompletedProcess(cmd, 0, stdout=cmd[-1], stderr="")

    monkeypatch.setattr(run_local.subprocess, "Popen", popen)
    monkeypatch.setattr(run_local.os, "killpg", lambda pid, sig: events.append(("killpg", pid, sig)))
    monkeypatch.setattr(run_local.subprocess, "run", docker)
    assert run_local.run_sample(row, args, tmp_path / "out") == 124
    container = seen["container"]
    assert events == [("wait", run_local.PIPELINE_TIMEOUT), ("killpg", 12345, signal.SIGKILL),
                      ("wait", run_local.CLEANUP_TIMEOUT),
                      ("docker", ["docker", "rm", "--force", container], run_local.CLEANUP_TIMEOUT)]
    wd = tmp_path / "out" / "R1" / "S1"
    assert "TIMEOUT" in (wd / "pipeline.log").read_text()
    assert "process group 12345" in (wd / "pipeline.log").read_text()
    assert not (wd / ".cleanup-unconfirmed").exists()


@pytest.mark.parametrize("failure", ["kill_error", "wait_timeout", "wait_error", "docker_error", "docker_timeout"])
def test_unconfirmed_cleanup_blocks_reuse_and_always_attempts_container_cleanup(tmp_path, monkeypatch, failure):
    row, args = coverage_sample(tmp_path, monkeypatch)
    commands = []

    class FakeProcess:
        pid = 12345
        first = True

        def wait(self, timeout):
            if self.first:
                self.first = False
                raise subprocess.TimeoutExpired("bash", timeout)
            if failure == "wait_timeout":
                raise subprocess.TimeoutExpired("bash", timeout)
            if failure == "wait_error":
                raise OSError("wait failed")
            return 0

    def killpg(*args):
        if failure == "kill_error":
            raise PermissionError("kill failed")

    def docker(cmd, **kwargs):
        commands.append(cmd)
        assert kwargs["timeout"] == run_local.CLEANUP_TIMEOUT
        if failure == "docker_timeout":
            raise subprocess.TimeoutExpired(cmd, kwargs["timeout"])
        return subprocess.CompletedProcess(cmd, int(failure == "docker_error"), stdout="", stderr="unavailable")

    monkeypatch.setattr(run_local.subprocess, "Popen", lambda *a, **kw: FakeProcess())
    monkeypatch.setattr(run_local.os, "killpg", killpg)
    monkeypatch.setattr(run_local.subprocess, "run", docker)
    out = tmp_path / "out"
    assert run_local.run_sample(row, args, out) == 124
    wd = out / "R1" / "S1"
    marker = wd / ".cleanup-unconfirmed"
    container = marker.read_text().strip()
    assert commands == [["docker", "rm", "--force", container]]
    assert "cleanup not confirmed" in (wd / "pipeline.log").read_text()
    assert "process group 12345" in (wd / "pipeline.log").read_text()
    (wd / "evidence").write_text("preserve")

    def must_not_start(*args, **kwargs):
        pytest.fail("sample with unconfirmed cleanup was restarted")

    monkeypatch.setattr(run_local.subprocess, "Popen", must_not_start)
    assert run_local.run_sample(row, args, out) == 1
    assert (wd / "evidence").read_text() == "preserve"
    assert marker.read_text().strip() == container
    assert "Refusing to reuse" in (wd / "pipeline.log").read_text()


def test_cleanup_missing_docker_is_unconfirmed(monkeypatch):
    def missing(*args, **kwargs):
        raise FileNotFoundError("docker")

    monkeypatch.setattr(run_local.subprocess, "run", missing)
    log = io.StringIO()
    assert not run_local.stop_container("cnvkit-harness-test", log)
    assert "cleanup not confirmed" in log.getvalue()
