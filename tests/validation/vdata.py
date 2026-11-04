"""Shared helpers for the L3 validation tests: data location, sample lists, loading scripts as modules, and running
the bundled scripts through their command lines, as pipeline.sh does.

Every per-sequencing result is keyed by (run, sample), never by the sample ID alone: the 226 sequencings have only 150
distinct sample IDs (replicates share an ID across runs), so a sample-keyed table would let one replicate's result
overwrite another's (round 6, V3). index() rejects a duplicate key, and covers() an incomplete or extra key set.
"""
import concurrent.futures as cf
import csv
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("CNVKIT_VALIDATION_DATA", REPO / "tests" / ".data" / "validation")).resolve()
BIN = REPO / "resources" / "usr" / "local" / "bin"
REFERENCE = REPO / "tests" / "validation" / "reference_scripts"
DATASETS = ("umi", "nonumi")
QC_FAILED = "25330S0104"     # the failed 0.5x library, reported separately in the masked-probe totals
N_SEQUENCINGS = 226          # in the final-configuration manifest (150 distinct sample IDs)
SCRIPT_TIMEOUT = 600         # seconds for one bundled-script command on one sequencing's validation tables


def require(path):
    if not path.exists():
        pytest.skip(f"validation data absent: {path} (fetch it: python3 tests/validation/fetch_validation_data.py all)")
    return path


def sequencings():
    """The 226 validation sequencings as (run, sample), from the final-configuration manifest."""
    manifest = require(DATA / "manifest" / "11_final_per_sample_files.csv")
    with open(manifest) as fh:
        seqs = sorted({tuple(r["record_id"].split("/", 1)) for r in csv.DictReader(fh) if r["call_set"] == "main"})
    if len(seqs) != N_SEQUENCINGS:
        raise AssertionError(f"{manifest}: {len(seqs)} sequencings, expected {N_SEQUENCINGS}")
    return seqs


def index(rows, key, where):
    """{key(row): row}; raises AssertionError on a duplicate key (a replicate overwriting another)."""
    out = {}
    for r in rows:
        k = key(r)
        if k in out:
            raise AssertionError(f"{where}: duplicate entry for {k}")
        out[k] = r
    return out


def covers(keys, expected, where):
    """Raise AssertionError unless keys is exactly the expected set (no sequencing missing, none extra)."""
    keys, expected = set(keys), set(expected)
    if keys != expected:
        raise AssertionError(f"{where}: missing {sorted(expected - keys)[:5]}, extra {sorted(keys - expected)[:5]}")


# the statistics cnvkit_parfix.py prints and the validation's parfix_log.tsv records (round 7, W4: each is required)
SHIFT_KEYS = ("par_bins", "span_bins", "par_segs", "span_segs")
COMPARE_KEYS = ("par_segs", "span_segs", "empty_after", "nonpar_identical", "par_rows_before", "par_rows_after",
                "rows_changed", "file_changed", "par_before", "par_after", "par_check")
COMPARE_PORT_KEYS = ("span_genes",)   # compare adds it with par_check span (and `unexplained` with 0, which fails)


def stats(stdout, keys, extra=()):
    """`key=value key=value ...` (shift and compare output) as a dict, which must have exactly the statistics keys
    + extra, each once (round 7, W4): a missing, repeated or unexpected one raises AssertionError, so a comparison
    cannot pass on only the keys both sides happen to have."""
    pairs = [kv.split("=", 1) for kv in stdout.split()]
    names = [kv[0] for kv in pairs]
    if any(len(kv) != 2 for kv in pairs) or len(set(names)) != len(names) or set(names) != {*keys, *extra}:
        raise AssertionError(f"expected exactly the statistics {', '.join([*keys, *extra])}, got {stdout.strip()!r}")
    return dict(pairs)


def logged(rec, keys):
    """{key: rec[key]} for a validation-log row; raises AssertionError unless the log records every key."""
    missing = [k for k in keys if rec.get(k) is None]
    if missing:
        raise AssertionError(f"the validation log has no {', '.join(missing)} for {rec.get('rundir')} "
                             f"{rec.get('sample')}")
    return {k: rec[k] for k in keys}


def script(name, *args, cwd):
    """Run a bundled script's command line (BIN/name ARGS, in cwd), as pipeline.sh does; a CompletedProcess. No
    descendants, so a plain bounded subprocess.run (procutil's policy for bundled scripts)."""
    return subprocess.run([sys.executable, str(BIN / name), *map(str, args)], cwd=cwd, capture_output=True,
                          text=True, errors="replace", timeout=SCRIPT_TIMEOUT)


def parallel(fn, items):
    """[fn(item) for item in items], run in threads (each call mostly waits on a subprocess)."""
    with cf.ThreadPoolExecutor(os.cpu_count() or 2) as ex:
        return list(ex.map(fn, items))


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
