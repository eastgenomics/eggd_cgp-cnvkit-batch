"""L3: the PAR post-step, run through cnvkit_parfix.py's command line exactly as pipeline.sh runs it, with the CNVkit
CLI in the cgp-cnvkit:1.0.0 image, reproduces the validated post-step on the 226 validation sequencings of each set
(UMI-consensus A0, non-UMI A0 and non-UMI M).

Each sequencing (keyed by run and sample: replicates share a sample ID across runs) gets its own directory with copies
of its <s>.cnr and <s>.call.cns, where pipeline.sh's commands run with its file names:
  cnvkit.py sex <s>.cnr --diploid-parx-genome grch38 -o <s>.sex_cnr.tsv
  cnvkit.py genemetrics <s>.cnr -s <s>.call.cns -t 0.3 --min-probes 3 --diploid-parx-genome grch38 [-x SEX]
  cnvkit_parfix.py decide --sex-arg SEX --cnr-sex <s>.sex_cnr.tsv --parx grch38
and, for a treated male,
  cnvkit_parfix.py shift <s>.cnr <s>.call.cns parfix
  (cd parfix && cnvkit.py genemetrics <the same arguments>)
  cnvkit_parfix.py compare <s>.genemetrics.noparfix.tsv parfix/<s>.genemetrics.tsv <s>.call.cns
For every sequencing: the sex table names <s>.cnr and its call equals the validation log, decide accepts it (its
sex-table checks), the decision equals the log, and genemetrics on the unmodified inputs reproduces the original
read-out byte for byte (argument check, also for -x). For every treated male: shift and compare succeed (compare's
genemetrics checks and per-row PAR verdicts included), the shift counts and every compare statistic equal the log
(par_check included), and the post-step genemetrics.tsv is
byte-identical to the validated one. Each command must print exactly its full set of statistics (SHIFT_KEYS,
COMPARE_KEYS, and span_genes for a span), and the log must record each: a missing statistic fails (round 7, W4).
"""
import csv
import os
import shlex
import shutil
import uuid

import pytest

from procutil import remove_container, run_group
from vdata import (COMPARE_KEYS, COMPARE_PORT_KEYS, DATA, REFERENCE, SHIFT_KEYS, covers, index, load_module, logged,
                   parallel, require, script, sequencings, stats)

pytestmark = [pytest.mark.validation, pytest.mark.docker]

IMAGE = "cgp-cnvkit:1.0.0"
IMAGE_TIMEOUT = 3600   # seconds per `docker run` (the batched genemetrics / sex calls over the 226 sequencings)
GM_ARGS = ["-t", "0.3", "--min-probes", "3", "--diploid-parx-genome", "grch38"]   # pipeline.sh gm_args, before -x
# (dataset, arm label, call directory, arm name in parfix_log, rundir filter)
SETS = [
    ("umi", "A0", "umi/final", "D_gcclip_on_s2", "/phaseB/calls/D_gcclip_on_s2/d00/"),
    ("nonumi", "A0", "nonumi/calls/A0", "A0_s2", "/d00/"),
    ("nonumi", "M", "nonumi/calls/M", "M_s2", "/d00/"),
]
# treated males per set: (par_check 1, par_check span, span samples with no PAR row before or after the fix)
PAR_CHECKS = {"umi-A0": (51, 9, 4), "nonumi-A0": (50, 10, 4), "nonumi-M": (50, 10, 4)}


def lims_sex_arg(sample):
    f = sample.split("-")
    return {"M": "male", "F": "female"}.get(f[4] if len(f) > 4 else "", "")


def in_image(tmp_path, script_text):
    """Run a shell script in the image with the data (read-only) and tmp_path mounted at their own paths. The
    container is named, so a timed-out or interrupted run (client process group killed) also removes its container."""
    name = f"cnvkit-l3-{uuid.uuid4().hex}"
    cmd = ["docker", "run", "--rm", "--name", name, "-u", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp",
           "-e", "MPLCONFIGDIR=/tmp",
           "-v", f"{DATA.resolve()}:{DATA.resolve()}:ro", "-v", f"{tmp_path}:{tmp_path}", "-w", str(tmp_path),
           IMAGE, "bash", "-c", script_text]
    proc = run_group(cmd, IMAGE_TIMEOUT, cleanup=lambda: remove_container(name), capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr[-3000:]


def image_jobs(tmp_path, jobs, log):
    """jobs: [(workdir, shell command)], run in the image in parallel, each in its workdir with its output in
    <workdir>/<log>; a failing job fails the run, named with the tail of its log."""
    listing = tmp_path / f"{log}.jobs"
    listing.write_bytes(b"".join(f"{wd}\0{cmd}\0".encode() for wd, cmd in jobs))
    inner = (f'cd "$0" && eval "$1" > {log} 2>&1 || '
             f'{{ echo "FAILED in $0: $1" >&2; tail -n 5 {log} >&2; exit 1; }}')
    in_image(tmp_path, f'xargs -0 -n 2 -P "$(nproc)" -a {shlex.quote(str(listing))} sh -c {shlex.quote(inner)}')


def genemetrics_cmd(s, sex_arg):
    """pipeline.sh's genemetrics command line, in the directory holding <s>.cnr and <s>.call.cns."""
    x = ["-x", sex_arg] if sex_arg else []
    return shlex.join(["cnvkit.py", "genemetrics", f"{s}.cnr", "-s", f"{s}.call.cns", *GM_ARGS, *x,
                       "-o", f"{s}.genemetrics.tsv"])


@pytest.mark.parametrize("dataset, arm, calls, log_arm, rundir_filter", SETS, ids=[f"{s[0]}-{s[1]}" for s in SETS])
def test_parfix_equivalence(tmp_path, dataset, arm, calls, log_arm, rundir_filter):
    calls = require(DATA / calls)
    with open(require(DATA / dataset / "parfix_log.tsv")) as fh:
        log = index((r for r in csv.DictReader(fh, delimiter="\t")
                     if r["arm"] == log_arm and rundir_filter in r["rundir"]),
                    lambda r: (r["rundir"].rsplit("/", 1)[1], r["sample"]), "parfix_log")
    seqs = sequencings()
    covers(log, seqs, "parfix_log")

    # per sequencing: its own directory, pipeline.sh's file names; cnvkit.py sex + genemetrics (unmodified inputs)
    wd = {(run, s): tmp_path / run / s for run, s in seqs}
    jobs = []
    for run, s in seqs:
        wd[(run, s)].mkdir(parents=True)
        for ext in (".cnr", ".call.cns"):
            shutil.copyfile(calls / run / f"{s}{ext}", wd[(run, s)] / f"{s}{ext}")
        sex = shlex.join(["cnvkit.py", "sex", f"{s}.cnr", "--diploid-parx-genome", "grch38", "-o", f"{s}.sex_cnr.tsv"])
        jobs.append((wd[(run, s)], f"{sex} && {genemetrics_cmd(s, lims_sex_arg(s))}"))
    image_jobs(tmp_path, jobs, "cnvkit.log")

    def before_the_fix(key):
        run, s = key
        d, rec, sex_arg, fails = wd[key], log[key], lims_sex_arg(s), []
        if (rec["plan_sex"] if rec["plan_sex"] != "-" else "") != sex_arg:
            fails.append("LIMS -x differs from the log")
        with open(d / f"{s}.sex_cnr.tsv") as fh:
            sex_rows = list(csv.DictReader(fh, delimiter="\t"))
        if [(r["sample"], r["sex"]) for r in sex_rows] != [(f"{s}.cnr", rec["inferred"])]:
            fails.append(f"cnvkit.py sex {sex_rows} != {s}.cnr {rec['inferred']}")
        original = calls / run / (f"{s}.genemetrics.noparfix.tsv" if rec["treated_male"] == "1" else
                                  f"{s}.genemetrics.tsv")
        if (d / f"{s}.genemetrics.tsv").read_bytes() != original.read_bytes():
            fails.append("genemetrics re-run differs from the original read-out")
        shutil.copyfile(d / f"{s}.genemetrics.tsv", d / f"{s}.genemetrics.noparfix.tsv")
        p = script("cnvkit_parfix.py", "decide", "--sex-arg", sex_arg, "--cnr-sex", f"{s}.sex_cnr.tsv",
                   "--parx", "grch38", cwd=d)
        decision = p.stdout.split("\t", 1)[0]
        if p.returncode or decision != rec["treated_male"]:
            fails.append(f"decide: {p.stdout.strip()!r} (exit {p.returncode}: {p.stderr.strip()}), "
                         f"log treated_male {rec['treated_male']}")
        elif decision == "1":
            p = script("cnvkit_parfix.py", "shift", f"{s}.cnr", f"{s}.call.cns", "parfix", cwd=d)
            if p.returncode or stats(p.stdout, SHIFT_KEYS) != logged(rec, SHIFT_KEYS):
                fails.append(f"shift: {p.stdout.strip()!r} (exit {p.returncode}: {p.stderr.strip()}), log "
                             + " ".join(f"{k}={rec[k]}" for k in SHIFT_KEYS))
        return [f"{run}/{s}: {f}" for f in fails], decision == "1" and not fails

    results = dict(zip(seqs, parallel(before_the_fix, seqs), strict=True))
    failures = [f for fails, _ in results.values() for f in fails]
    males = [key for key, (_, male) in results.items() if male]
    assert not failures, "\n".join(failures[:20])

    image_jobs(tmp_path, [(wd[(run, s)] / "parfix", genemetrics_cmd(s, lims_sex_arg(s))) for run, s in males],
               "cnvkit.log")

    def after_the_fix(key):
        run, s = key
        d, rec, fails = wd[key], log[key], []
        if (d / "parfix" / f"{s}.genemetrics.tsv").read_bytes() != (calls / run / f"{s}.genemetrics.tsv").read_bytes():
            fails.append("post-step genemetrics differs from the validated read-out")
        p = script("cnvkit_parfix.py", "compare", f"{s}.genemetrics.noparfix.tsv", f"parfix/{s}.genemetrics.tsv",
                   f"{s}.call.cns", cwd=d)
        if p.returncode:
            return [f"{run}/{s}: compare (exit {p.returncode}: {p.stderr.strip()})"], {}
        # every statistic the log records, each compared (round 7, W4); span_genes, a port addition, iff span
        got = stats(p.stdout, COMPARE_KEYS, COMPARE_PORT_KEYS if "par_check=span" in p.stdout.split() else ())
        want = logged(rec, COMPARE_KEYS)
        diff = {k: (got[k], want[k]) for k in COMPARE_KEYS if got[k] != want[k]}
        if diff:
            fails.append(f"compare differs from the log (port, log): {diff}")
        return [f"{run}/{s}: {f}" for f in fails], got

    after = dict(zip(males, parallel(after_the_fix, males), strict=True))
    failures = [f for fails, _ in after.values() for f in fails]
    assert not failures, "\n".join(failures[:20])
    assert len(males) == sum(r["treated_male"] == "1" for r in log.values()) == 60
    checks = [got["par_check"] for _, got in after.values()]
    span_without_par_rows = sum(got["par_check"] == "span" and got["par_rows_before"] == got["par_rows_after"] == "0"
                                for _, got in after.values())
    assert (checks.count("1"), checks.count("span"), span_without_par_rows) == PAR_CHECKS[f"{dataset}-{arm}"]
    assert "0" not in checks


# Z1 (round 9): overlapping segments, with CNVkit's own genemetrics. CNVkit writes overlapping segments for
# overlapping targets, and reports a gene once per segment whose range overlaps its bins (by_ranges, mode 'outer').
# A synthetic .cnr with overlapping CRLF2 bins in PAR1 and a call.cns with two overlapping segments, run through
# pipeline.sh's PAR commands for a male (-x male): the post-step must complete and every PAR row be explained. Where
# both segments contain the reported row ('both'), the statistics and stdout must equal the base's (8b683c8), whose
# compare is the pinned reference's.
CNR_H = "chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n"
CNS_H = "chromosome\tstart\tend\tgene\tlog2\tcn\tdepth\tprobes\tweight\n"
OVERLAP_CNR = CNR_H + "".join(
    [f"chr17\t{7000000 + k * 1000}\t{7000200 + k * 1000}\tTP53\t-0.6\t400\t0.9\n" for k in range(10)]
    + [f"chrX\t{s}\t{e}\tCRLF2\t1.5\t600\t0.9\n" for s, e in ((1190891, 1191100), (1191000, 1191300),
                                                             (1191200, 1191500), (1191600, 1191900),
                                                             (1200000, 1200300), (1212000, 1212654))]
    + [f"chrX\t{5000000 + k * 1000}\t{5000200 + k * 1000}\tNONPAR\t0.0\t500\t0.9\n" for k in range(5)])
OVERLAP_CNS = {   # case -> the two overlapping CRLF2 segments (start, end, log2, probes, weight)
    "both": ((1100000, 1300000, "0.5", 6, 5.4), (1190891, 2781479, "0.5", 6, 5.4)),
    "split": ((1190891, 1191300, "0.5", 2, 1.8), (1191200, 1212654, "0.8", 4, 3.6)),
    "split_same": ((1190891, 1191300, "0.5", 2, 1.8), (1191200, 1212654, "0.5", 4, 3.6)),
    "split_new": ((1190891, 1191300, "0.5", 2, 1.8), (1191200, 1212654, "-0.8", 4, 3.6)),
}
BASE_BOTH = ("par_segs=2 span_segs=0 empty_after=0 nonpar_identical=1 par_rows_before=2 par_rows_after=2 "
             "rows_changed=2 file_changed=1 par_before=CRLF2:1.5 par_after=CRLF2:0.5 par_check=1\n")


def test_parfix_overlapping_segments_in_the_image(tmp_path):
    s = "S9-RUN-PANEL-0001-M-1"
    for case, segs in OVERLAP_CNS.items():
        d = tmp_path / case
        d.mkdir()
        (d / f"{s}.cnr").write_text(OVERLAP_CNR)
        (d / f"{s}.call.cns").write_text(CNS_H + "chr17\t7000000\t7009200\tTP53\t-0.6\t1\t400\t10\t9\n" + "".join(
            f"chrX\t{a}\t{b}\tCRLF2\t{lg}\t2\t600\t{n}\t{w}\n" for a, b, lg, n, w in segs)
            + "chrX\t5000000\t5004200\tNONPAR\t0.0\t2\t500\t5\t4.5\n")
        p = script("cnvkit_parfix.py", "shift", f"{s}.cnr", f"{s}.call.cns", "parfix", cwd=d)
        assert p.returncode == 0 and p.stdout.startswith("par_bins=6 span_bins=0 par_segs=2 span_segs=0"), p.stderr
    image_jobs(tmp_path, [(tmp_path / case / sub, genemetrics_cmd(s, "male")) for case in OVERLAP_CNS
                          for sub in ("", "parfix")], "cnvkit.log")
    ref = load_module(REFERENCE / "genemetrics_parfix.py", "genemetrics_parfix_ref_overlap")
    for case in OVERLAP_CNS:
        d = tmp_path / case
        os.replace(d / f"{s}.genemetrics.tsv", d / f"{s}.genemetrics.noparfix.tsv")
        p = script("cnvkit_parfix.py", "compare", f"{s}.genemetrics.noparfix.tsv", f"parfix/{s}.genemetrics.tsv",
                   f"{s}.call.cns", cwd=d)
        assert p.returncode == 0 and p.stdout.strip().endswith("par_check=1"), (case, p.stdout, p.stderr)
        assert "CRLF2" in (d / "parfix" / f"{s}.genemetrics.tsv").read_text(), case    # a PAR row was judged
        theirs = ref.compare(dict(action="dry", rundir=str(d), sample=s,
                                  orig=str(d / f"{s}.genemetrics.noparfix.tsv"),
                                  new=str(d / "parfix" / f"{s}.genemetrics.tsv")))
        if case == "both":    # the row lies inside both segments: the base's answer, byte for byte
            assert p.stdout == BASE_BOTH and theirs["par_check"] == 1
        else:                 # the row starts before its own segment: the base (containment) cannot explain it
            assert theirs["par_check"] == 0
