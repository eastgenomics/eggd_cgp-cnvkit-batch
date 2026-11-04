"""L1: the pinned reference scripts are the unmodified validation originals."""
import hashlib

import pytest

from cnvfiles import REFERENCE

pytestmark = pytest.mark.unit


# the originals the ports name in their docstrings (round 7: all three are checked, so an emptied or partial MD5SUMS
# cannot pass)
EXPECTED = {"deadprobe_mask.py": "382acdecb5be92d6074f8957362bdfd9",
            "genemetrics_parfix.py": "51ada0e597420655c89c86e1a2e13ea6",
            "sex_check.py": "a6106d50c43cd96ebde9ae3f24360e2b"}


def test_reference_scripts_unmodified():
    listed = dict(reversed(line.split()) for line in (REFERENCE / "MD5SUMS").read_text().splitlines())
    assert listed == EXPECTED
    assert {p.name for p in REFERENCE.glob("*.py")} == set(EXPECTED)
    for name, md5 in EXPECTED.items():
        assert hashlib.md5((REFERENCE / name).read_bytes()).hexdigest() == md5, name
