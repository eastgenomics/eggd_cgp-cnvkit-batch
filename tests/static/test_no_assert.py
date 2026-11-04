"""L0: the bundled scripts use no `assert` for control flow (round 4, T2).

`python -O` (or PYTHONOPTIMIZE) removes assert statements, so a check written as one would silently pass. Every check
in the resource scripts must raise an explicit exception.
"""
import ast

import pytest

from apprun import REPO

pytestmark = pytest.mark.static

SCRIPTS = sorted((REPO / "resources" / "usr" / "local" / "bin").glob("*.py"))


def test_scripts_found():
    assert {p.name for p in SCRIPTS} >= {"cnvkit_parfix.py", "cnvkit_sex_qc.py", "dropout_mask.py"}


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_no_assert_statements(script):
    tree = ast.parse(script.read_text(), filename=str(script))
    found = [n.lineno for n in ast.walk(tree) if isinstance(n, ast.Assert)]
    assert not found, f"{script.name}: assert on line(s) {found}; raise an explicit exception instead"
