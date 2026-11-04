"""Real-Docker check that a timed-out harness run stops only its own container."""
import argparse
import uuid

import pytest

import run_local
from procutil import removing_containers, run_group

pytestmark = [pytest.mark.validation, pytest.mark.docker]


def test_timeout_stops_only_run_container(tmp_path, monkeypatch):
    token, other_token = uuid.uuid4().hex, uuid.uuid4().hex
    name, other = f"cnvkit-harness-{token}", f"cnvkit-harness-control-{other_token}"
    monkeypatch.setattr(run_local.uuid, "uuid4", lambda: type("Token", (), {"hex": token})())
    monkeypatch.setattr(run_local, "PIPELINE_TIMEOUT", 5)
    monkeypatch.setattr(run_local, "DATA", tmp_path / "data")
    pon = tmp_path / "data" / "nonumi" / "pons" / "pon.cnn"
    pon.parent.mkdir(parents=True)
    pon.write_text("reference\n")
    cnn = tmp_path / "data" / "nonumi" / "cnn" / "R1" / "S1.targetcoverage.cnn"
    cnn.parent.mkdir(parents=True)
    cnn.write_text("coverage\n")
    lib = tmp_path / "pipeline.sh"
    lib.write_text('run_pipeline() { docker run --rm --name "$HARNESS_CONTAINER_NAME" '
                   f'{run_local.IMAGE} /bin/sleep 120; }}\n')
    monkeypatch.setattr(run_local, "LIB", lib)
    row = {"run": "R1", "sample": "S1", "cn_reference": "pon"}
    args = argparse.Namespace(data="nonumi", start_from="coverage", nproc=1, mask=True)
    # Only the two UUID-named containers created by this test are removed, each confirmed (removed or already gone);
    # an unconfirmed removal fails the test, and never replaces a failure of its own (procutil.removing_containers)
    with removing_containers(name, other):
        started = run_group(["docker", "run", "--rm", "--detach", "--name", other, run_local.IMAGE, "/bin/sleep",
                             "120"], 30, capture_output=True, text=True)
        assert started.returncode == 0, started.stderr
        assert run_local.run_sample(row, args, tmp_path / "out") == 124
        log = (tmp_path / "out" / "R1" / "S1" / "pipeline.log").read_text()
        assert f"Removed timed-out container {name}" in log
        gone = run_group(["docker", "inspect", name], 30, capture_output=True, text=True)
        control = run_group(["docker", "inspect", "--format", "{{.State.Running}}", other], 30,
                            capture_output=True, text=True)
        assert gone.returncode != 0
        assert control.returncode == 0 and control.stdout.strip() == "true", control.stderr
