#!/usr/bin/env python3
"""L4 local harness (PLAN_v2.1.0_build.md Phase 5): the app's pipeline.sh with the real cgp-cnvkit:1.0.0 image,
on the local test sequencings, compared byte for byte with the validation baselines.

Each sequencing runs in its own directory with its validated PoN (leave-out for PoN members, shared otherwise) and
--diploid-parx-genome grch38, the Atlas setting, and starts from the validated coverage (or from the BAM).

Baselines (§0b): mask on -> U-M / N-M; mask off -> U-A0 (umi/final) / N-A0. U-M genemetrics.tsv is from before the
PAR post-step, so for UMI with the mask on the app's genemetrics.noparfix.tsv is compared with it, and the post-step
genemetrics.tsv with U-M (not treated as male) or U-A0 (treated male, nothing masked); a treated male with a masked
probe has no UMI post-step baseline (reported as n/a; the non-UMI N-M covers that path).

Usage
  run_local.py --data umi|nonumi [--dropout-mask true|false] [--start-from coverage|bam] [--samples S ...]
               [--jobs N] [--label NAME]
  --start-from bam downloads the non-UMI BAM/BAI and the bins BED (cached under tests/.data/harness/inputs/).
Writes tests/.data/harness/<label>/<run>/<sample>/ and <label>/report.tsv; exit 1 on any FAIL.
"""
import argparse
import concurrent.futures as cf
import csv
import os
import shutil
import subprocess
import sys
import threading
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "validation"))
from vdata import DATA  # noqa: E402

HARNESS = REPO / "tests" / ".data" / "harness"
LIB = REPO / "resources" / "usr" / "local" / "lib" / "cgp_cnvkit_batch" / "pipeline.sh"
BIN = REPO / "resources" / "usr" / "local" / "bin"
SAMPLES = REPO / "tests" / "harness" / "local_samples.tsv"
IMAGE = "cgp-cnvkit:1.0.0"
BED = "project-Fkb6Gkj433GVVvj73J7x8KbV:file-JBxvjg8433GxJf9Xypj5qyXf"   # atlas_cnvkit_bins_for_cgp_hg38_v1.1.bed
PIPELINE_TIMEOUT = 2 * 3600   # seconds per sequencing; a hung container is reported as a FAIL (exit 124)
_BED_LOCK = threading.Lock()
QC_HEADER = "arm\tchromosome\tstart\tend\tgene\tlog2\traw_reads\treason\n"
SEX_FIELDS_UMI = ["lims_sex", "raw_Y_logratio", "chrY_present", "chrY_sex", "raw_X_logratio", "raw_cnvkit_sex",
                  "cnr_X_logratio", "cnr_cnvkit_sex", "sex_mismatch"]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def local_samples():
    return read_tsv(SAMPLES)


def treated_male_log(data, mask):
    """{(run, sample): '0'|'1'} from the validation PAR post-step log of the matching arm."""
    arm = {("umi", False): "D_gcclip_on_s2", ("umi", True): "D_gcclip_on_s2",
           ("nonumi", False): "A0_s2", ("nonumi", True): "M_s2"}[(data, mask)]
    return {(r["rundir"].rsplit("/", 1)[1], r["sample"]): r["treated_male"]
            for r in read_tsv(DATA / data / "parfix_log.tsv") if r["arm"] == arm and "/d00/" in r["rundir"]}


def fetch_bam_inputs(run, s):
    """Non-UMI BAM + BAI (the eggd_chr_prefix output the coverage app used) and the bins BED, cached."""
    inv = {(r["run"], r["sample_id"]): r for r in read_tsv(DATA / "manifest" / "nonumi_inventory.tsv")}
    bam_ref = inv[(run, s)]["new_bam"]
    d = HARNESS / "inputs" / run / s
    d.mkdir(parents=True, exist_ok=True)
    bam, bai, bed = d / "tumour.bam", d / "tumour.bam.bai", HARNESS / "inputs" / "panel.bed"
    project = bam_ref.split(":")[0]
    if not bam.exists():
        subprocess.run(["dx", "download", "-f", bam_ref, "-o", f"{bam}.part"], check=True)
        os.replace(f"{bam}.part", bam)
    if not bai.exists():
        name = subprocess.run(["dx", "describe", bam_ref, "--name"], check=True, capture_output=True,
                              text=True).stdout.strip()
        found = subprocess.run(["dx", "find", "data", "--project", project, "--name", name + ".bai", "--brief"],
                               check=True, capture_output=True, text=True).stdout.split()
        if len(found) != 1:
            raise RuntimeError(f"{s}: expected one {name}.bai, found {found}")
        subprocess.run(["dx", "download", "-f", found[0], "-o", str(bai)], check=True)
    with _BED_LOCK:   # shared by all workers: download once, publish the complete file atomically
        if not bed.exists():
            subprocess.run(["dx", "download", "-f", BED, "-o", f"{bed}.part"], check=True)
            os.replace(f"{bed}.part", bed)
    return bam, bai, bed


def run_sample(row, a, out):
    run, s = row["run"], row["sample"]
    wd = out / run / s
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True)
    shutil.copyfile(DATA / a.data / "pons" / f"{row['cn_reference']}.cnn", wd / "reference.cnn")
    if a.start_from == "coverage":
        shutil.copyfile(DATA / a.data / "cnn" / run / f"{s}.targetcoverage.cnn", wd / f"{s}.targetcoverage.cnn")
    else:
        bam, bai, bed = fetch_bam_inputs(run, s)
        os.symlink(bam, wd / "tumour.bam")
        os.symlink(bai, wd / "tumour.bam.bai")
        shutil.copyfile(bed, wd / "panel.bed")
    env = {**os.environ, "PATH": f"{BIN}:{os.environ['PATH']}",
           "SAMPLE_ID": s, "BAM_NAME": f"{s}.bam", "START_FROM": a.start_from, "NPROC": str(a.nproc),
           "DIPLOID_PARX_GENOME": "grch38", "SEX_OVERRIDE": "", "DROP_LOW_COVERAGE": "true",
           "DROPOUT_MASK": "true" if a.mask else "false"}
    # BAM inputs are symlinks outside the work directory: mount their directory at the same path
    extra = f'-v "{HARNESS / "inputs"}:{HARNESS / "inputs"}:ro"' if a.start_from == "bam" else ""
    script = (f'set -eo pipefail; source "{LIB}"; '
              f'run_cnvkit() {{ docker run --rm -u {os.getuid()}:{os.getgid()} -e HOME=/tmp -e MPLCONFIGDIR=/tmp '
              f'{extra} -v "$(pwd)":/work -w /work {IMAGE} cnvkit.py "$@"; }}; run_pipeline')
    with open(wd / "pipeline.log", "w") as log:
        try:
            proc = subprocess.run(["bash", "-c", script], cwd=wd, env=env, stdout=log, stderr=subprocess.STDOUT,
                                  timeout=PIPELINE_TIMEOUT)
        except subprocess.TimeoutExpired:
            log.write(f"\nTIMEOUT after {PIPELINE_TIMEOUT} s\n")
            return 124
    return proc.returncode


def same(a, b):
    return a.exists() and b.exists() and a.read_bytes() == b.read_bytes()


def check_sample(row, a, out, rc, male_log):
    run, s = row["run"], row["sample"]
    wd = out / run / s
    if rc != 0:
        return [(s, "pipeline", "FAIL", f"exit {rc}; see {wd / 'pipeline.log'}")]
    res = []

    def chk(name, ok, detail=""):
        res.append((s, name, "PASS" if ok else "FAIL", detail))

    if a.mask:
        base = DATA / a.data / "calls" / "M" / run
    else:
        base = DATA / "umi" / "final" / run if a.data == "umi" else DATA / "nonumi" / "calls" / "A0" / run
    for ext in (".cnr", ".cns", ".call.cns"):
        chk(ext, same(wd / f"{s}{ext}", base / f"{s}{ext}"))
    if a.start_from == "bam":
        chk("targetcoverage.cnn", same(wd / f"{s}.targetcoverage.cnn",
                                       DATA / a.data / "cnn" / run / f"{s}.targetcoverage.cnn"),
            "batch-app coverage vs the coverage app (PoN input)")

    male = male_log[(run, s)]
    qc = read_tsv(wd / f"{s}.sex_qc.tsv")[0]
    chk("post-step decision", qc["par_fix_applied"] == male, f"app {qc['par_fix_applied']}, validation {male}")
    gm, nop = wd / f"{s}.genemetrics.tsv", wd / f"{s}.genemetrics.noparfix.tsv"
    if a.data == "umi" and a.mask:          # U-M has no post-step
        chk("genemetrics.noparfix.tsv", same(nop, base / f"{s}.genemetrics.tsv"), "vs U-M genemetrics (pre-post-step)")
        masked = len((wd / f"{s}.dropout_qc.tsv").read_text().splitlines()) > 1
        if male == "0":
            chk("genemetrics.tsv", same(gm, base / f"{s}.genemetrics.tsv"), "not treated as male: vs U-M")
        elif not masked:
            chk("genemetrics.tsv", same(gm, DATA / "umi" / "final" / run / f"{s}.genemetrics.tsv"),
                "treated male, nothing masked: vs U-A0 post-step")
        else:
            res.append((s, "genemetrics.tsv", "n/a", "treated male with a masked probe: no UMI post-step baseline"))
    else:
        chk("genemetrics.tsv", same(gm, base / f"{s}.genemetrics.tsv"))
        expected_nop = base / f"{s}.genemetrics.noparfix.tsv"
        chk("genemetrics.noparfix.tsv", same(nop, expected_nop if expected_nop.exists() else base / f"{s}.genemetrics.tsv"))
    dq = (wd / f"{s}.dropout_qc.tsv").read_text()
    if a.mask:
        chk("dropout_qc.tsv", dq == (base / f"{s}.mask_qc.tsv").read_text(), f"{len(dq.splitlines()) - 1} masked")
    else:
        chk("dropout_qc.tsv", dq == QC_HEADER, "header only")

    if a.data == "umi":
        truth = {(r["run"], r["sample"]): r for r in read_tsv(DATA / "umi" / "sex_qc.tsv")}[(run, s)]
        diff = [f for f in SEX_FIELDS_UMI if qc[f] != truth[f]]
    else:
        raw = {r["sample"].split("__", 1)[1][:-len(".targetcoverage.cnn")]: r
               for r in read_tsv(DATA / "nonumi" / "sex_raw" / run / "sex.tsv")}[s]
        diff = [f for f, t in (("raw_Y_logratio", "Y_logratio"), ("raw_X_logratio", "X_logratio"),
                               ("raw_cnvkit_sex", "sex")) if qc[f] != raw[t]]
        y = float(raw["Y_logratio"]) if raw["Y_logratio"] not in ("", "nan") else float("nan")
        if qc["chrY_sex"] != ("Male" if y == y and y > -5 else "Female"):
            diff.append("chrY_sex")
    chk("sex_qc", not diff, ",".join(diff) or f"chrY {qc['chrY_sex']}, mismatch {qc['sex_mismatch']}")
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, choices=("umi", "nonumi"))
    ap.add_argument("--dropout-mask", default="true", choices=("true", "false"))
    ap.add_argument("--start-from", default="coverage", choices=("coverage", "bam"))
    ap.add_argument("--samples", nargs="*", help="sample IDs or specimen IDs (default: all local test sequencings)")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--nproc", type=int, default=2)
    ap.add_argument("--label")
    a = ap.parse_args(argv)
    a.mask = a.dropout_mask == "true"
    if a.start_from == "bam" and a.data != "nonumi":
        sys.exit("--start-from bam is for the non-UMI (production) BAMs")
    rows = local_samples()
    if a.samples:
        rows = [r for r in rows if any(x in r["sample"] for x in a.samples)]
        if not rows:
            sys.exit("no local test sequencing matches --samples " + " ".join(a.samples))
    label = a.label or f"{a.data}-mask_{a.dropout_mask}-{a.start_from}"
    out = HARNESS / label
    out.mkdir(parents=True, exist_ok=True)
    male_log = treated_male_log(a.data, a.mask)
    print(f"{label}: {len(rows)} sequencings, {a.jobs} in parallel -> {out}", flush=True)
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        rcs = dict(zip([(r["run"], r["sample"]) for r in rows], ex.map(lambda r: run_sample(r, a, out), rows)))
    results = []
    for r in rows:
        results += check_sample(r, a, out, rcs[(r["run"], r["sample"])], male_log)
    with open(out / "report.tsv", "w") as fh:
        fh.write("sample\tcheck\tresult\tdetail\n")
        for x in results:
            fh.write("\t".join(x) + "\n")
    by_sample = {}
    for s, check, result, detail in results:
        by_sample.setdefault(s, []).append((check, result, detail))
    for r in rows:
        checks = by_sample[r["sample"]]
        fails = [f"{c} ({d})" if d else c for c, res, d in checks if res == "FAIL"]
        na = sum(res == "n/a" for _, res, _ in checks)
        status = "FAIL" if fails else "PASS"
        print(f"{status}  {r['run']:15s} {r['sample']:52s} {len(checks) - len(fails) - na} pass"
              f"{', ' + str(na) + ' n/a' if na else ''}{'; ' + '; '.join(fails) if fails else ''}")
    n_fail = sum(res == "FAIL" for _, _, res, _ in results)
    print(f"{label}: {len(results)} checks, {n_fail} FAIL; report {out / 'report.tsv'}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
