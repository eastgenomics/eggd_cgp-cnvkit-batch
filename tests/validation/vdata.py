"""Shared helpers for the L3 validation tests: data location, sample lists, loading scripts as modules."""
import csv
import importlib.util
import os
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("CNVKIT_VALIDATION_DATA", REPO / "tests" / ".data" / "validation"))
BIN = REPO / "resources" / "usr" / "local" / "bin"
REFERENCE = REPO / "tests" / "validation" / "reference_scripts"
DATASETS = ("umi", "nonumi")
QC_FAILED = "25330S0104"     # the failed 0.5x library, reported separately in the masked-probe totals


def require(path):
    if not path.exists():
        pytest.skip(f"validation data absent: {path} (run tests/validation/fetch_validation_data.py)")
    return path


def sequencings():
    """The 226 validation sequencings as (run, sample), from the final-configuration manifest."""
    manifest = require(DATA / "manifest" / "11_final_per_sample_files.csv")
    with open(manifest) as fh:
        return sorted({tuple(r["record_id"].split("/", 1)) for r in csv.DictReader(fh) if r["call_set"] == "main"})


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
