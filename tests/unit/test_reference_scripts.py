"""L1: the pinned reference scripts are the unmodified validation originals."""
import hashlib

import pytest

from cnvfiles import REFERENCE

pytestmark = pytest.mark.unit


def test_reference_scripts_unmodified():
    """Every pinned reference script matches its MD5 in MD5SUMS (an identity check, not a security use)."""
    for line in (REFERENCE / "MD5SUMS").read_text().splitlines():
        md5, name = line.split()
        assert hashlib.md5((REFERENCE / name).read_bytes(), usedforsecurity=False).hexdigest() == md5, name
