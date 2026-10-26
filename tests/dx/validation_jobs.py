#!/usr/bin/env python3
"""L5/L6 on DNAnexus (PLAN_v2.1.0_build.md Phase 7): build the applet, run the T matrix and the E3 equivalence test
(226 sequencings from the non-UMI BAMs, dropout mask on, validated non-UMI PoNs), and check the outputs against the
baselines with the same checks as the local harness.

Everything runs in the validation project 004_261008_cgp_cnvkit_pon_comparison under /batch_app_v2.1.0/. Job IDs
are recorded in tests/dx/runs/ (committed: they are the evidence); downloads go to tests/.data/dx/.

Usage
  validation_jobs.py prepare            BAI lookup, PoN uploads -> tests/dx/runs/inputs.tsv
  validation_jobs.py build              build the applet from the committed HEAD -> tests/dx/runs/applet.txt
  validation_jobs.py submit T|E3        submit the set (skips rows already submitted)
  validation_jobs.py status T|E3        job states
  validation_jobs.py check T|E3         download outputs and logs, check against the baselines -> report
"""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO / "tests" / "validation"), str(REPO / "tests" / "harness")]
from vdata import DATA  # noqa: E402
import run_local  # noqa: E402

PROJECT = "project-JF3p1f04Yggp0vf3PYjy7Q0F"
ROOT = "/batch_app_v2.1.0"
RUNS = REPO / "tests" / "dx" / "runs"
DL = REPO / "tests" / ".data" / "dx"
BED = "project-Fkb6Gkj433GVVvj73J7x8KbV:file-JBxvjg8433GxJf9Xypj5qyXf"         # bins BED v1.1
CONTROLLED_PON = f"{PROJECT}:file-JFP3fzQ455KjXVZ23JgX1gF2"   # = validated non-UMI shared PoN (md5 b4a37372...)
PLAN = f"{PROJECT}:file-JFPV4pQ4Yggg5f8K30VGV02F"              # real_plan_nonumi.tsv: the PoN of every call
OUTPUTS = ["copy_ratios", "segments", "call_segments", "genemetrics", "genemetrics_noparfix", "sex_qc",
           "sex_check_mqc", "dropout_qc", "target_coverage", "scatter_png", "diagram_pdf"]


def dx(*args, check=True):
    p = subprocess.run(["dx", *args], capture_output=True, text=True)
    if check and p.returncode:
        raise RuntimeError(f"dx {' '.join(args)}: {p.stderr.strip()}")
    return p.stdout


def dxjson(*args):
    return json.loads(dx(*args, "--json"))


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ── prepare ──────────────────────────────────────────────────────────────────

def bai_index(bam_refs):
    """{bam file ID: bai project:file ID}, one `dx find` per BAM folder."""
    folders, out = {}, {}
    for ref in bam_refs:
        d = dxjson("describe", ref)
        folders.setdefault((d["project"], d["folder"]), []).append((d["id"], d["name"]))
    for (project, folder), bams in folders.items():
        found = json.loads(dx("find", "data", "--project", project, "--folder", folder, "--name", "*.bai", "--json"))
        by_name = {f["describe"]["name"]: f["id"] for f in found}
        for fid, name in bams:
            out[fid] = f"{project}:{by_name[name + '.bai']}"
    return out


def upload_pons(names, local_dir, folder):
    """Upload the named PoNs (once) into 004:<folder>; returns {name: project:file ID}."""
    dx("mkdir", "-p", f"{PROJECT}:{folder}")
    existing = {f["describe"]["name"]: f["id"]
                for f in json.loads(dx("find", "data", "--project", PROJECT, "--folder", folder, "--json"))}
    out = {}
    for n in sorted(names):
        fname = f"{n}.cnn"
        if fname not in existing:
            existing[fname] = dx("upload", str(local_dir / fname), "--destination", f"{PROJECT}:{folder}/",
                                 "--wait", "--brief").strip()
        out[n] = f"{PROJECT}:{existing[fname]}"
    return out


def prepare(_):
    plan_path = DL / "real_plan_nonumi.tsv"
    DL.mkdir(parents=True, exist_ok=True)
    dx("download", "-f", PLAN, "-o", str(plan_path))
    plan = {(r["run"], r["sample_id"]): r for r in read_tsv(plan_path)}
    inv = {(r["run"], r["sample_id"]): r for r in read_tsv(DATA / "manifest" / "nonumi_inventory.tsv")}
    assert set(plan) == set(inv), "plan and inventory differ"
    bams = [inv[k]["new_bam"] for k in sorted(inv)]
    bais = bai_index(bams)
    loo = {r["pon_id"] for r in plan.values() if "shared" not in r["pon_id"]}
    pons = upload_pons(loo, DATA / "nonumi" / "pons", f"{ROOT}/inputs/nonumi_pons")
    rows = []
    for k in sorted(plan):
        run, s = k
        bam = inv[k]["new_bam"]
        pon_id = plan[k]["pon_id"]
        rows.append(dict(run=run, sample=s, bam=bam, bai=bais[bam.split(":")[1]], pon_id=pon_id,
                         cn_reference=CONTROLLED_PON if "shared" in pon_id else pons[pon_id]))
    write_tsv(RUNS / "inputs.tsv", rows, ["run", "sample", "bam", "bai", "pon_id", "cn_reference"])
    print(f"inputs.tsv: {len(rows)} sequencings; {len(pons)} leave-out PoNs in {ROOT}/inputs/nonumi_pons")


# ── build ────────────────────────────────────────────────────────────────────

def build(_):
    status = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=REPO,
                            capture_output=True, text=True).stdout
    if status.strip():
        sys.exit("commit the working tree first: the applet must be built from a recorded commit")
    sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO, capture_output=True,
                         text=True).stdout.strip()
    folder = f"{ROOT}/applet/{sha}"
    dx("mkdir", "-p", f"{PROJECT}:{folder}")
    out = subprocess.run(["dx", "build", str(REPO), "--destination", f"{PROJECT}:{folder}/", "--overwrite"],
                         capture_output=True, text=True)
    if out.returncode:
        sys.exit(out.stderr)
    applet = json.loads(out.stdout)["id"]
    dx("set_properties", f"{PROJECT}:{applet}", f"git_commit={sha}")
    (RUNS / "applet.txt").write_text(f"{applet}\t{sha}\n")
    print(f"built {applet} from {sha}")


def applet():
    applet_id, sha = (RUNS / "applet.txt").read_text().split()
    return applet_id, sha


# ── job sets ─────────────────────────────────────────────────────────────────

U = "135024088-25021S0072-26CGPVal1SHv150kb-9999-U-1111"
M = "140609591-25328S0023-26CGPV2-8475-M-96770770"


def t_matrix(inputs):
    """The T cases as dicts: case, run, sample, sample_id, extra (app inputs beyond the defaults), row.

    Default inputs: diploid_parx_genome=grch38, the validated non-UMI PoN, dropout_mask as given. N-A0 cases run with
    the mask off (the baseline without the mask); E3 covers the mask on for all 226. T15 of the plan (controlled
    PoN) is covered by E3: every non-member runs with the controlled PoN. T14 (female) is T2. T13 (UMI BAMs) is
    dropped with E2: UMI data is no longer produced, and L3/L4 cover it.
    """
    def find(spec, run=None):
        hits = [r for r in inputs if spec in r["sample"] and (run is None or r["run"] == run)]
        assert len(hits) == 1, (spec, run, len(hits))
        return hits[0]
    m = find("25328S0023", "0797_V2")
    off = {"dropout_mask": False}
    on = {"dropout_mask": True}
    cases = [
        ("T1_lims_M", m, None, off),
        ("T2_lims_F", find("25178S0009", "0786_V1"), None, off),
        ("T3_lims_U", find("25021S0072"), None, off),
        ("T4_xxy", find("26127S0031"), None, off),
        ("T5_override_M_on_U", find("25050S0022"), None, {**on, "sex_override": "male"}),
        ("T6_override_F_on_M", m, None, {**on, "sex_override": "F"}),
        ("T7_bad_override", m, None, {**on, "sex_override": "X"}),
        ("T8_letter_from_bam_name", m, "25328S0023", on),
        ("T9_no_letter", find("25021S0072"), "25021S0072", on),
        ("T10_parx_empty", m, None, {**on, "diploid_parx_genome": None}),
        ("T11_mask_off_masked_sample", find("25330S0048"), None, off),
        ("T12_deprecated_inputs", m, None, {**off, "purity": 0.5, "ploidy": 3.1, "sample_sex": "female"}),
    ]
    return [dict(case=c, run=r["run"], sample=r["sample"], sample_id=sid or r["sample"], extra=x, row=r)
            for c, r, sid, x in cases]


def run_inputs(row, sample_id, extra):
    args = ["-itumour_bam=" + row["bam"], "-itumour_bai=" + row["bai"], "-ibaits=" + BED,
            "-icn_reference=" + row["cn_reference"], "-isample_id=" + sample_id]
    extra = dict(extra)
    parx = extra.pop("diploid_parx_genome", "grch38")   # None: input not supplied (empty)
    if parx:
        args.append(f"-idiploid_parx_genome={parx}")
    for k, v in extra.items():
        if v is None:
            continue
        args.append(f"-i{k}={str(v).lower() if isinstance(v, bool) else v}")
    return args


def jobs_path(which):
    return RUNS / f"{which}_jobs.tsv"


def submit(a):
    applet_id, sha = applet()
    inputs = read_tsv(RUNS / "inputs.tsv")
    done = {r["key"]: r for r in read_tsv(jobs_path(a.set))} if jobs_path(a.set).exists() else {}
    todo = []
    if a.set == "T":
        for c in t_matrix(inputs):
            todo.append((c["case"], c["row"], c["sample_id"], c["extra"], f"T/{c['case']}"))
    else:
        for r in inputs:
            todo.append((f"{r['run']}/{r['sample']}", r, r["sample"], {"dropout_mask": True},
                         f"E3/{r['run']}/{r['sample']}"))
    rows = list(done.values())
    for key, row, sample_id, extra, folder in todo:
        if key in done:
            continue
        cmd = ["run", applet_id, *run_inputs(row, sample_id, extra), "--destination", f"{PROJECT}:{ROOT}/{folder}",
               "--name", f"cnvkit-batch v2.1.0 {a.set} {key}", "--property", f"validation_set={a.set}",
               "--property", f"git_commit={sha}", "--brief", "-y"]
        job = dx(*cmd).strip()
        rows.append(dict(key=key, run=row["run"], sample=row["sample"], sample_id=sample_id,
                         extra=json.dumps(extra), job=job))
        write_tsv(jobs_path(a.set), rows, ["key", "run", "sample", "sample_id", "extra", "job"])
        print(f"{key}\t{job}", flush=True)
    print(f"{a.set}: {len(rows)} jobs recorded in {jobs_path(a.set)}")


def states(which):
    rows = read_tsv(jobs_path(which))
    out = {}
    for i in range(0, len(rows), 100):
        for d in json.loads(dx("describe", *[r["job"] for r in rows[i:i + 100]], "--json", "--multi")):
            out[d["id"]] = d
    return rows, out


def status(a):
    rows, desc = states(a.set)
    from collections import Counter
    print(Counter(desc[r["job"]]["state"] for r in rows))


# ── check ────────────────────────────────────────────────────────────────────

def download_outputs(job_desc, dest):
    dest.mkdir(parents=True, exist_ok=True)
    for name, link in (job_desc.get("output") or {}).items():
        fid = link["$dnanexus_link"]
        fid = fid if isinstance(fid, str) else fid["id"]
        target = dest / dx("describe", fid, "--name").strip()
        if not target.exists():
            dx("download", f"{job_desc['project']}:{fid}", "-o", str(target), "-f")
    log = dest / "job.log"
    if not log.exists():
        log.write_text(dx("watch", job_desc["id"], "--get-streams", "--no-timestamps", "--no-job-info", "-q",
                          check=False))
    return log.read_text()


def check(a):
    rows, desc = states(a.set)
    results = []
    male = {**run_local.treated_male_log("nonumi", True)}
    male_a0 = run_local.treated_male_log("nonumi", False)
    for r in rows:
        d = desc[r["job"]]
        key, s, run = r["key"], r["sample"], r["run"]
        extra = json.loads(r["extra"])
        expect = "fail" if extra.get("sex_override") == "X" else None

        def add(check_name, ok, detail=""):
            results.append((key, r["job"], check_name, "PASS" if ok else "FAIL", detail))

        if expect == "fail":
            log = dx("watch", r["job"], "--get-streams", "--no-timestamps", "--no-job-info", "-q", check=False)
            add("job failed", d["state"] == "failed", d["state"])
            add("before coverage", "not accepted" in log and "[coverage]" not in log)
            continue
        if d["state"] != "done":
            add("job done", False, d["state"])
            continue
        add("job done", True)
        sid = r["sample_id"]
        dest = DL / a.set / key.replace("/", "__")
        log = download_outputs(d, dest)
        missing = [o for o in OUTPUTS if o not in (d.get("output") or {})]
        add("all outputs", not missing, ",".join(missing))
        add("edge correction (log)", "Correcting for density bias" in log)
        mask = extra.get("dropout_mask", True)
        if a.set == "E3" or key.startswith(("T1", "T2", "T3", "T4", "T11")):
            # same file names as the baseline: the harness checks, on the downloaded outputs
            ns = argparse.Namespace(data="nonumi", mask=mask, start_from="bam")
            tmp = DL / f"{a.set}_layout"
            link = tmp / run / s
            link.parent.mkdir(parents=True, exist_ok=True)
            if not link.exists():
                link.symlink_to(dest.resolve())
            for _, c, res, det in run_local.check_sample({"run": run, "sample": s}, ns, tmp, 0,
                                                          male if mask else male_a0):
                results.append((key, r["job"], c, res, det))
        elif key.startswith(("T8", "T9")):
            # sample_id without the 5th field: same content as N-M, other file names
            base = DATA / "nonumi" / "calls" / "M" / run
            for ext in (".cnr", ".cns", ".call.cns", ".genemetrics.tsv"):
                add(ext, (dest / f"{sid}{ext}").read_bytes() == (base / f"{s}{ext}").read_bytes())
            src = "bam_name" if key.startswith("T8") else "cnvkit_inference"
            add("sex source", f"(source: {src}" in log, src)
        elif key.startswith("T5"):
            qc = read_tsv(dest / f"{sid}.sex_qc.tsv")[0]
            add("override applied", qc["applied_sex"] == "male" and qc["applied_source"] == "override")
            add("post-step (male)", qc["par_fix_applied"] == "1")
            add("sex_mismatch (declared M, no chrY)", qc["sex_mismatch"] == "1", qc["chrY_sex"])
        elif key.startswith("T6"):
            qc = read_tsv(dest / f"{sid}.sex_qc.tsv")[0]
            add("override applied", qc["applied_sex"] == "female" and qc["applied_source"] == "override")
            add("no post-step", qc["par_fix_applied"] == "0")
            add("sex_mismatch (declared F, chrY present)", qc["sex_mismatch"] == "1")
        elif key.startswith("T10"):
            add("parx warning", "[parx] WARNING: diploid_parx_genome is empty" in log)
            add("no post-step", "[parfix] not applied: no --diploid-parx-genome" in log)
            add("noparfix = genemetrics", (dest / f"{sid}.genemetrics.tsv").read_bytes()
                == (dest / f"{sid}.genemetrics.noparfix.tsv").read_bytes())
        elif key.startswith("T12"):
            add("deprecation warnings", "[deprecated] purity/ploidy supplied but ignored" in log
                and "[deprecated] sample_sex supplied but ignored" in log)
            t1 = DL / "T" / "T1_lims_M"
            for ext in (".cnr", ".cns", ".call.cns", ".genemetrics.tsv", ".genemetrics.noparfix.tsv", ".sex_qc.tsv"):
                add(f"{ext} = T1", (dest / f"{sid}{ext}").read_bytes() == (t1 / f"{sid}{ext}").read_bytes())
    report = RUNS / f"{a.set}_report.tsv"
    write_tsv(report, [dict(zip(["key", "job", "check", "result", "detail"], x)) for x in results],
              ["key", "job", "check", "result", "detail"])
    fails = [x for x in results if x[3] == "FAIL"]
    if a.set == "E3":
        qcs = [read_tsv(DL / "E3" / r["key"].replace("/", "__") / f"{r['sample']}.sex_qc.tsv")[0]
               for r in rows if desc[r["job"]]["state"] == "done"]
        mf = [q for q in qcs if q["lims_sex"] in ("M", "F")]
        print(f"LIMS M/F: {len(mf)}; sex_mismatch = 1: {sum(q['sex_mismatch'] == '1' for q in mf)}")
    for x in fails[:40]:
        print("FAIL", *x)
    print(f"{a.set}: {len(rows)} jobs, {len(results)} checks, {len(fails)} FAIL -> {report}")
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("prepare")
    sp.add_parser("build")
    for c in ("submit", "status", "check"):
        p = sp.add_parser(c)
        p.add_argument("set", choices=("T", "E3"))
    a = ap.parse_args()
    return {"prepare": prepare, "build": build, "submit": submit, "status": status, "check": check}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main() or 0)
