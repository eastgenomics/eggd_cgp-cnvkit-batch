"""Shared pytest hooks.

Docker availability is checked lazily, when a `docker`-marked test is set up, and never at import or collection:
the Docker-free CI selection (-m "static or unit or pipeline") deselects those tests and so never runs `docker`.
The check is bounded (process group killed on timeout: procutil.run_group) and cached for the session. A missing
docker or image skips the test; a check whose own cleanup is unconfirmed (a `docker` process may remain) fails it, as
a skip would hide that (round 5, U3).
"""
import functools
import shutil
import subprocess

import pytest

from procutil import cleanup_unconfirmed, failure_text, run_group, timed_out

IMAGE = "cgp-cnvkit:1.0.0"
DOCKER_CHECK_TIMEOUT = 30   # seconds for `docker image inspect`


@functools.lru_cache(maxsize=None)
def docker_check(image=IMAGE):
    """(action, reason): (None, None) when docker and the image are usable, else ('skip' | 'fail', reason)."""
    if shutil.which("docker") is None:
        return "skip", "docker not installed"
    try:
        proc = run_group(["docker", "image", "inspect", image], DOCKER_CHECK_TIMEOUT, capture_output=True, text=True)
    except subprocess.TimeoutExpired as e:
        return "fail" if cleanup_unconfirmed(e) else "skip", timed_out(f"docker image inspect {image}",
                                                                         DOCKER_CHECK_TIMEOUT, e)
    except OSError as e:
        return "skip", f"docker not runnable: {failure_text(e)}"
    if proc.returncode == 0:
        return None, None
    err = proc.stderr.strip()[-300:]
    return "skip", f"{image} not available" + (f" ({err})" if err else "")


def pytest_runtest_setup(item):
    if item.get_closest_marker("docker"):
        action, reason = docker_check()
        if action == "fail":
            pytest.fail(reason, pytrace=False)
        if action == "skip":
            pytest.skip(reason)
