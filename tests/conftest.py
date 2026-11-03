"""Shared pytest hooks.

Docker availability is checked lazily, when a `docker`-marked test is set up, and never at import or collection:
the Docker-free CI selection (-m "static or unit or pipeline") deselects those tests and so never runs `docker`.
The check is bounded and cached for the session.
"""
import functools
import shutil
import subprocess

import pytest

IMAGE = "cgp-cnvkit:1.0.0"
DOCKER_CHECK_TIMEOUT = 30   # seconds for `docker image inspect`


@functools.lru_cache(maxsize=None)
def docker_unavailable(image=IMAGE):
    """None when docker and the image are usable, else the reason to skip."""
    if shutil.which("docker") is None:
        return "docker not installed"
    try:
        proc = subprocess.run(["docker", "image", "inspect", image], capture_output=True,
                              timeout=DOCKER_CHECK_TIMEOUT)
    except subprocess.TimeoutExpired:
        return f"docker image inspect {image} timed out after {DOCKER_CHECK_TIMEOUT} s"
    except OSError as e:
        return f"docker not runnable: {e}"
    return None if proc.returncode == 0 else f"{image} not available"


def pytest_runtest_setup(item):
    if item.get_closest_marker("docker"):
        reason = docker_unavailable()
        if reason:
            pytest.skip(reason)
