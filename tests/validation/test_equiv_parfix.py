"""L3: cnvkit_parfix.py + the CNVkit CLI reproduce the validated PAR post-step on the 226 validation sequencings
(UMI-consensus A0, non-UMI A0 and non-UMI M), in the cgp-cnvkit:1.0.0 image.

For every sequencing: the male decision (port `decide` on `cnvkit.py sex` of the .cnr) equals the validation log, and
re-running genemetrics on the unmodified inputs reproduces the original read-out (argument check, also for -x). For
every treated male: `shift` + `cnvkit.py genemetrics` in a per-sample directory + `compare` reproduce the validated
post-step genemetrics.tsv byte for byte, with the same par_check.
"""
import csv
import os
import shutil
import subprocess

import pytest

from vdata import BIN, DATA, load_module, require, sequencings

pytestmark = [pytest.mark.validation, pytest.mark.docker]

IMAGE = "cgp-cnvkit:1.0.0"
GM_ARGS = ["-t", "0.3", "--min-probes", "3", "--diploid-parx-genome", "grch38"]
# (dataset, arm label, call directory, arm name in parfix_log, rundir filter)
SETS = [
    ("umi", "A0", "umi/final", "D_gcclip_on_s2", "/phaseB/calls/D_gcclip_on_s2/d00/"),
    ("nonumi", "A0", "nonumi/calls/A0", "A0_s2", "/d00/"),
    ("nonumi", "M", "nonumi/calls/M", "M_s2", "/d00/"),
]


def docker_available():
    return shutil.which("docker") and subprocess.run(["docker", "image", "inspect", IMAGE],
                                                     capture_output=True).returncode == 0


def lims_sex_arg(sample):
    f = sample.split("-")
    return {"M": "male", "F": "female"}.get(f[4] if len(f) > 4 else "", "")


def in_image(tmp_path, script):
    """Run a shell script in the image with the data (read-only) and tmp_path mounted at their own paths."""
    cmd = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp", "-e", "MPLCONFIGDIR=/tmp",
           "-v", f"{DATA.resolve()}:{DATA.resolve()}:ro", "-v", f"{tmp_path}:{tmp_path}", "-w", str(tmp_path),
           IMAGE, "bash", "-c", script]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr[-3000:]


def genemetrics_jobs(tmp_path, jobs):
    """jobs: [(workdir, cnr name, cns name, out name, sex_arg)], each run as `cnvkit.py genemetrics` in workdir."""
    listing = tmp_path / "jobs.tsv"
    listing.write_text("".join("\t".join(j) + "\n" for j in jobs))
    gm = " ".join(GM_ARGS)
    in_image(tmp_path, f"""set -euo pipefail
        while IFS=$'\\t' read -r wd cnr cns out sex; do
          printf '%s\\0%s\\0%s\\0%s\\0%s\\0' "$wd" "$cnr" "$cns" "$out" "$sex"
        done < {listing} | xargs -0 -n 5 -P "$(nproc)" sh -c \\
          'cd "$0" && cnvkit.py genemetrics "$1" -s "$2" {gm} ${{4:+-x $4}} -o "$3" > /dev/null 2> "$3.log"'
    """)


@pytest.fixture(scope="module")
def port():
    return load_module(BIN / "cnvkit_parfix.py", "cnvkit_parfix_port")


@pytest.mark.skipif(not docker_available(), reason=f"docker or {IMAGE} not available")
@pytest.mark.parametrize("dataset, arm, calls, log_arm, rundir_filter", SETS, ids=[f"{s[0]}-{s[1]}" for s in SETS])
def test_parfix_equivalence(tmp_path, port, dataset, arm, calls, log_arm, rundir_filter):
    calls = require(DATA / calls)
    log = {}
    with open(require(DATA / dataset / "parfix_log.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["arm"] == log_arm and rundir_filter in r["rundir"]:
                log[(r["rundir"].rsplit("/", 1)[1], r["sample"])] = r
    seqs = sequencings()
    assert set(log) == set(seqs), "parfix_log does not cover the 226 sequencings"

    # cnvkit.py sex on every .cnr, in one call
    cnrs = [str(calls / run / f"{s}.cnr") for run, s in seqs]
    (tmp_path / "cnrs.txt").write_text("\n".join(cnrs) + "\n")
    in_image(tmp_path, f"xargs -a cnrs.txt cnvkit.py sex --diploid-parx-genome grch38 -o {tmp_path}/sex_cnr.tsv")
    with open(tmp_path / "sex_cnr.tsv") as fh:
        inferred = {os.path.basename(r["sample"])[:-4]: r["sex"] for r in csv.DictReader(fh, delimiter="\t")}

    jobs, males = [], []
    for run, s in seqs:
        rec = log[(run, s)]
        sex_arg = lims_sex_arg(s)
        assert (rec["plan_sex"] if rec["plan_sex"] != "-" else "") == sex_arg, (run, s, "LIMS -x")
        assert inferred[s] == rec["inferred"], (run, s, "cnvkit.py sex")
        male, _ = port.decide(sex_arg, inferred[s], "grch38")
        assert int(male) == int(rec["treated_male"]), (run, s, "decision")
        src = calls / run
        original = src / f"{s}.genemetrics.noparfix.tsv" if male else src / f"{s}.genemetrics.tsv"
        # argument check: genemetrics on the unmodified inputs reproduces the original read-out
        rerun = tmp_path / "rerun" / run / s
        rerun.mkdir(parents=True)
        for ext in (".cnr", ".call.cns"):
            shutil.copyfile(src / f"{s}{ext}", rerun / f"{s}{ext}")
        jobs.append((str(rerun), f"{s}.cnr", f"{s}.call.cns", f"{s}.genemetrics.tsv", sex_arg))
        if male:
            fix = tmp_path / "parfix" / run / s
            fix.mkdir(parents=True)
            port.shift_par(src / f"{s}.cnr", fix / f"{s}.cnr", 1, 4)
            port.shift_par(src / f"{s}.call.cns", fix / f"{s}.call.cns", 1, 4)
            jobs.append((str(fix), f"{s}.cnr", f"{s}.call.cns", f"{s}.genemetrics.tsv", sex_arg))
            males.append((run, s, original, fix, rec))
    genemetrics_jobs(tmp_path, jobs)

    failures = []
    for run, s in seqs:
        male = int(log[(run, s)]["treated_male"])
        original = calls / run / (f"{s}.genemetrics.noparfix.tsv" if male else f"{s}.genemetrics.tsv")
        if (tmp_path / "rerun" / run / s / f"{s}.genemetrics.tsv").read_bytes() != original.read_bytes():
            failures.append(f"{run}/{s}: genemetrics re-run differs from the original read-out")
    for run, s, original, fix, rec in males:
        new = fix / f"{s}.genemetrics.tsv"
        if new.read_bytes() != (calls / run / f"{s}.genemetrics.tsv").read_bytes():
            failures.append(f"{run}/{s}: post-step genemetrics differs from the validated read-out")
        got = port.compare(original, new, calls / run / f"{s}.call.cns")
        if str(got["par_check"]) != rec["par_check"]:
            failures.append(f"{run}/{s}: par_check {got['par_check']} != {rec['par_check']}")
    assert not failures, "\n".join(failures[:20])
    assert len(males) == 60
