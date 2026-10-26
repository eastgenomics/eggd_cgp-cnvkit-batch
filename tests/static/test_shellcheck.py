"""L0: shellcheck on the app's bash.

SC2154 (variable referenced but not assigned) is excluded: DNAnexus injects the job inputs as environment
variables, so shellcheck can't see them assigned. Everything else must be clean.
"""
import shutil
import subprocess

import pytest

from apprun import REPO

pytestmark = pytest.mark.static

SCRIPTS = [REPO / "src" / "code.sh", *sorted((REPO / "resources").rglob("*.sh")),
           *sorted((REPO / "tests").rglob("*.sh"))]

# v2.0.2 code.sh is not clean (unquoted word splitting: SC2046, SC2086; SC2002). The Phase 2 refactor fixes it;
# strict, so the suite fails as soon as it is clean and this entry must be removed.
KNOWN_UNCLEAN = {"src/code.sh": "v2.0.2: SC2046/SC2086/SC2002, fixed by the Phase 2 refactor"}


def _param(p):
    rel = str(p.relative_to(REPO))
    marks = [pytest.mark.xfail(reason=KNOWN_UNCLEAN[rel], strict=True)] if rel in KNOWN_UNCLEAN else []
    return pytest.param(p, id=rel, marks=marks)


@pytest.mark.skipif(shutil.which("shellcheck") is None, reason="shellcheck not installed")
@pytest.mark.parametrize("script", [_param(p) for p in SCRIPTS])
def test_shellcheck(script):
    proc = subprocess.run(["shellcheck", "-s", "bash", "-e", "SC2154", str(script)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout
