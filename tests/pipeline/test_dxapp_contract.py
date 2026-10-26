"""L0/L2: dxapp.json is consistent with what the app reads and produces."""
import fnmatch
import re

import pytest

from apprun import DXAPP, REPO, run_app

APP_SOURCES = [REPO / "src" / "code.sh", *sorted((REPO / "resources").rglob("*.sh"))]
OUTPUTS = {o["name"]: o for o in DXAPP["outputSpec"]}


def app_source():
    return "\n".join(p.read_text() for p in APP_SOURCES if p.exists())


@pytest.mark.static
def test_input_names_unique():
    names = [i["name"] for i in DXAPP["inputSpec"]]
    assert len(names) == len(set(names))


@pytest.mark.static
@pytest.mark.parametrize("name", [i["name"] for i in DXAPP["inputSpec"]])
def test_every_input_is_read(name):
    assert re.search(r"\$\{?" + re.escape(name) + r"\b", app_source()), f"input {name} is never read"


@pytest.mark.static
def test_defaults_are_optional():
    for spec in DXAPP["inputSpec"]:
        if "default" in spec:
            assert spec.get("optional") is True, spec["name"]


@pytest.mark.pipeline
def test_outputs_match_spec(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert set(run.outputs) <= set(OUTPUTS), "job adds an output not in outputSpec"
    required = {n for n, o in OUTPUTS.items() if not o.get("optional")}
    assert required <= set(run.outputs), f"required outputs missing: {required - set(run.outputs)}"
    for name, out in run.outputs.items():
        spec = OUTPUTS[name]
        assert out["class"] == spec["class"], name
        uploaded = out["value"].removeprefix("file-")
        assert any(fnmatch.fnmatch(uploaded, p) for p in spec.get("patterns", ["*"])), (name, uploaded)
