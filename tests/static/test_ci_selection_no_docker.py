"""L0: collecting and running the Docker-free CI selection never runs `docker`.

A recording `docker` that always fails is put first on PATH. Docker availability is checked only when a
`docker`-marked test is set up (tests/conftest.py), so neither collection nor the L0–L2 selection may call it.
"""
import os
import sys

import pytest

from apprun import REPO
from procutil import run_group

pytestmark = pytest.mark.static

SELECTION = "static or unit or pipeline"
NESTED = "CNVKIT_TEST_NESTED"


@pytest.fixture
def recording_docker(tmp_path):
    if os.environ.get(NESTED):
        pytest.skip("nested pytest run")
    bin_dir, calls = tmp_path / "bin", tmp_path / "docker_calls"
    bin_dir.mkdir()
    stub = bin_dir / "docker"
    stub.write_text(f'#!/bin/sh\necho "$@" >> "{calls}"\nexit 99\n')
    stub.chmod(0o755)
    env = dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}", **{NESTED: "1"})
    return env, calls


def pytest_run(env, *args):
    return run_group([sys.executable, "-m", "pytest", "-p", "no:cacheprovider", "-q", *args], 300,
                     cwd=REPO, env=env, capture_output=True, text=True)


def test_collection_does_not_run_docker(recording_docker):
    env, calls = recording_docker
    proc = pytest_run(env, "--collect-only", "-m", SELECTION)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    assert not calls.exists(), calls.read_text()


def test_ci_selection_over_docker_modules_does_not_run_docker(recording_docker):
    # the modules with docker-marked tests, run with the CI marker expression (their tests are deselected)
    env, calls = recording_docker
    proc = pytest_run(env, "-m", SELECTION, "tests/validation", "tests/harness", "tests/static/test_shellcheck.py")
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    assert " deselected" in proc.stdout
    assert not calls.exists(), calls.read_text()


def test_docker_tests_check_docker_lazily(recording_docker):
    # a docker-marked test that is selected is skipped at setup with a reason, after one bounded check
    env, calls = recording_docker
    proc = pytest_run(env, "-rs", "-m", "docker", "tests/harness/test_timeout_cleanup.py")
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    assert "cgp-cnvkit:1.0.0 not available" in proc.stdout
    assert calls.read_text().splitlines() == ["image inspect cgp-cnvkit:1.0.0"]
