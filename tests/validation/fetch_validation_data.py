#!/usr/bin/env python3
"""Download the validation data for the L3/L4 tests into $CNVKIT_VALIDATION_DATA (default tests/.data/validation).

Sources are pinned by DNAnexus file ID and sha256. Members of the outputs tarballs are
also checked against the sha256 lists published with them. Stdlib only; needs a logged-in `dx`.

Layout written (<run> = 0734_Val1_50kb | 0786_V1 | 0797_V2 | 0811_V3, <s> = sample ID):
  umi/cnn/<run>/<s>.targetcoverage.cnn            UMI-consensus coverage
  umi/calls/{A0,M}/<run>/<s>.*                    U-A0 / U-M raw calls, BEFORE the PAR post-step (+ M mask_qc.tsv)
  umi/final/<run>/<s>.*                           U-A0 read-out AFTER the PAR post-step (+ genemetrics.noparfix.tsv)
  umi/pons/*.cnn                                  validated UMI PoNs (shared + 124 leave-out)
  umi/sex_qc.tsv, umi/parfix_log.tsv
  nonumi/cnn/<run>/<s>.targetcoverage.cnn         non-UMI coverage (eggd_cgp-cnvkit-coverage 2.0.2)
  nonumi/calls/{A0,M}/<run>/<s>.*                 N-A0 / N-M, AFTER the PAR post-step (+ M mask_qc.tsv)
  nonumi/pons/*.cnn                               validated non-UMI PoNs (shared = the controlled PoN)
  nonumi/sex_raw/<run>/sex.tsv, nonumi/parfix_log.tsv
  pon/atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn the controlled static PoN
  manifest/{11_final_per_sample_files.csv,03_pon_samples.csv,nonumi_inventory.tsv}

Usage: fetch_validation_data.py [umi|nonumi|all]
"""
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

PROJECT = "project-JF3p1f04Yggp0vf3PYjy7Q0F"
REPO = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("CNVKIT_VALIDATION_DATA", REPO / "tests" / ".data" / "validation"))
DOWNLOADS = DATA / "_downloads"
RUNS = r"(0734_Val1_50kb|0786_V1|0797_V2|0811_V3)"

# (set, name, file ID, sha256 of the downloaded file, [(member regex, destination template)] or a destination path)
SOURCES = [
    ("umi", "cnn_shared.f7c499e87092.tar.gz", "file-JF61K7Q4YggpJkXxY6K9Jkz8", None,
     [(rf"^cnn/{RUNS}/([^/]+\.targetcoverage\.cnn)$", "umi/cnn/{0}/{1}")]),
    *[("umi", f"sim2c_real__{run}.calls.tar.gz", fid, None,
       [(rf"^calls/(A0|M)/d00/{RUNS}/([^/]+)$", "umi/calls/{0}/{1}/{2}")])
      for run, fid in [("0734_Val1_50kb", "file-JFGG4F04fzz6Fk5pqFxqPzbx"), ("0786_V1", "file-JFGG24j4By26f9j330BJV50z"),
                       ("0797_V2", "file-JFGG4BQ40z59b6k6bqxxVyVB"), ("0811_V3", "file-JFGG8Z046yXqfQkq6Z4Zvf0v")]],
    ("umi", "final_config.D_gcclip_on_s2.parfix.outputs.tar.gz", "file-JF61vJQ4Y2jj0V7FzGQ89X94", None,
     [(rf"^phaseB/calls/D_gcclip_on_s2/d00/{RUNS}/([^/]+)$", "umi/final/{0}/{1}")]),
    ("umi", "final_config.D_gcclip_on_s2.parfix.sha256.tsv", "file-JF61vJQ4Y2jy4xj3GyBp5YVy", None,
     "umi/final.sha256.tsv"),
    ("umi", "final_main__batch_000.pons.tar.gz", "file-JF61YzQ4P0z8QFVv534xJvXx", None,
     [(r"^phaseB/pons/([^/]+\.cnn)$", "umi/pons/{0}")]),
    ("umi", "final_config.sex_qc.tsv", "file-JF6X7gQ4YggvqyqFBxjq6ZXF", None, "umi/sex_qc.tsv"),
    ("umi", "parfix_log.tsv", "file-JF61vJQ4Y2jZ2g7gpXJvKPQ5", None, "umi/parfix_log.tsv"),
    ("nonumi", "cnn_shared.a96c4f793074.tar.gz", "file-JFP4JvQ4Yggg5f8K30VGPFj1", None,
     [(rf"^nonumi/cnn/{RUNS}/([^/]+\.targetcoverage\.cnn)$", "nonumi/cnn/{0}/{1}")]),
    ("nonumi", "dropout_nonumi_NM.calls.tar.gz", "file-JFPX1JQ4763FK094J6vqj1yQ", None,
     [(rf"^calls/(A0|M)/d00/{RUNS}/([^/]+)$", "nonumi/calls/{0}/{1}/{2}")]),
    ("nonumi", "dropout_nonumi_NM.sha256.tsv", "file-JFPX1J847638Qy9X54JyJxyV", None, "nonumi/calls.sha256.tsv"),
    ("nonumi", "nonumi_main__batch_000.pons.tar.gz", "file-JFP4Yvj4GJKXqfkxzk71FpYj", None,
     [(r"^nonumi/phaseN/pons/(D_GC__[^/]+\.cnn)$", "nonumi/pons/{0}")]),
    ("nonumi", "final_config_nonumi.tier1.tar.gz", "file-JFP4Q704VG07vK8B14BJgbZK", None,
     [(rf"^tier1/sex_raw/{RUNS}/sex\.tsv$", "nonumi/sex_raw/{0}/sex.tsv")]),
    ("nonumi", "parfix_log.tsv", "file-JFPX1J84763FG09Kk80y0ZJK", None, "nonumi/parfix_log.tsv"),
    ("nonumi", "cgp_cnvkit_static_pon_di4093_nonumi.cnn", "file-JFP3fzQ455KjXVZ23JgX1gF2", None,
     "pon/atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn"),
    ("all", "11_final_per_sample_files.csv", "file:/validation_manifest/11_final_per_sample_files.csv", None,
     "manifest/11_final_per_sample_files.csv"),
    ("all", "03_pon_samples.csv", "file:/validation_manifest/03_pon_samples.csv", None,
     "manifest/03_pon_samples.csv"),
    ("all", "nonumi_inventory.tsv", "file:/final_config_nonumi/repo_snapshot/e8bb628/inventory/nonumi_inventory.tsv",
     None, "manifest/nonumi_inventory.tsv"),
]

# sha256 of every downloaded source (pinned 2026-10-26 from the first download; tarball members are also checked
# against the published sha256 lists below, and the controlled PoN against its md5)
PINNED_SHA256 = {
    "umi/cnn_shared.f7c499e87092.tar.gz": "9327e1613776e6248c2fa8bdb16cc180139347baf3c8ab37866b67336dc87526",
    "umi/sim2c_real__0734_Val1_50kb.calls.tar.gz": "0a1176007c6ad7b7bc4c861d92eef805b570543c76407e74e4023a242a7c719e",
    "umi/sim2c_real__0786_V1.calls.tar.gz": "4068c55fed4809dafcacbc931ea7f27f3f21837df8c66f340bf6f54f8d94675c",
    "umi/sim2c_real__0797_V2.calls.tar.gz": "ff295083af973a8159f261e121e8f2741f57138189efcdf36f97ee425b94dc1d",
    "umi/sim2c_real__0811_V3.calls.tar.gz": "e7c949b0b19f04fa40d874880f15dcaad87bcc811edecea6709a7dd6d5a865b0",
    "umi/final_config.D_gcclip_on_s2.parfix.outputs.tar.gz": "436eb94835b3f70c1ea6fe1fe3b560ec6c28f596adc30d1107fcaf57cbde7c13",
    "umi/final_config.D_gcclip_on_s2.parfix.sha256.tsv": "1935775c005305e657dcbf4e2477aca8dbd7896b56e874f7893035d2e9e43139",
    "umi/final_main__batch_000.pons.tar.gz": "1e8b7a0b2756e32b36fcdb8602aa4b4001616870bf195d4eddd6aef2173a4a58",
    "umi/final_config.sex_qc.tsv": "6a7a1975779b0077f50e777965237a5accfe23cb034e90bfeb5cb509993abb6a",
    "umi/parfix_log.tsv": "02ff78e7f45449c81573f32eff13c577e0878348b800b6a0363072d684a9c7a1",
    "nonumi/cnn_shared.a96c4f793074.tar.gz": "bf532e205586c4e1f551f5c26dde7d4a71041235714b31dba6d4c430b9a9f7cb",
    "nonumi/dropout_nonumi_NM.calls.tar.gz": "2050b7192743069b74f7699b2da1df3293ca05a08e7a9e69018f2ab473934e55",
    "nonumi/dropout_nonumi_NM.sha256.tsv": "bbe00fafdbd2f33ae00c98fde4a962871b5962ef607c595c9e4246233f3a7ad1",
    "nonumi/nonumi_main__batch_000.pons.tar.gz": "3b124954562143c71888f84b52e6c92676df5a236d70e57c5f422d3676ba5632",
    "nonumi/final_config_nonumi.tier1.tar.gz": "f28067405c8001f2c2d2607aeb13a1a279dc8a37f8d58daf0960e6078cbea056",
    "nonumi/parfix_log.tsv": "4885c8266c4b38d4553d5d55898c1639f9caba04102d9a91f1d198f701f70105",
    "nonumi/cgp_cnvkit_static_pon_di4093_nonumi.cnn": "1fb3d8762100efbd1712bac2a1be9d5b4c05ed813f6f174e3018855788576595",
    "all/11_final_per_sample_files.csv": "1b9d5b633aab8326d624a65c985187b710d3d3b47dc8afadd65b3b9b2b1f2489",
    "all/03_pon_samples.csv": "f9fa217999958670c5bb88d90e10b1ab6a25bce971fa4bab8277273a913e768b",
    "all/nonumi_inventory.tsv": "52d7167644936f0804c9dc4f145325a7a38bc351e5c94fb9e3a99239db146c59",
}

# member sha256 lists: (list file, tarball-path prefix to strip, destination prefix)
MEMBER_CHECKS = [
    ("umi/final.sha256.tsv", rf"^phaseB/calls/D_gcclip_on_s2/d00/{RUNS}/([^/]+)$", "umi/final/{0}/{1}"),
    ("nonumi/calls.sha256.tsv", rf"^calls/(A0|M)/d00/{RUNS}/([^/]+)$", "nonumi/calls/{0}/{1}/{2}"),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download(name, fid):
    DOWNLOADS.mkdir(parents=True, exist_ok=True)
    key = fid[5:].strip("/").replace("/", "_") if fid.startswith("file:") else fid
    dest = DOWNLOADS / f"{key}__{name}"
    if not dest.exists():
        ref = f"{PROJECT}:{fid[5:]}" if fid.startswith("file:") else f"{PROJECT}:{fid}"
        print(f"  download {name}", flush=True)
        subprocess.run(["dx", "download", "-f", ref, "-o", str(dest) + ".part"], check=True,
                       stdout=subprocess.DEVNULL)
        os.replace(str(dest) + ".part", dest)
    return dest


def extract(tar_path, rules):
    n = 0
    with tarfile.open(tar_path) as tar:
        for m in tar:
            if not m.isfile():
                continue
            for rx, template in rules:
                hit = re.match(rx, m.name)
                if hit:
                    dest = DATA / template.format(*hit.groups())
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    with tar.extractfile(m) as src, open(dest, "wb") as out:
                        shutil.copyfileobj(src, out)
                    n += 1
                    break
    return n


def main(which):
    print(f"validation data -> {DATA}")
    failures = []
    for group, name, fid, _, target in SOURCES:
        if which != "all" and group not in (which, "all"):
            continue
        path = download(name, fid)
        got = sha256(path)
        key = f"{group}/{name}"
        if key in PINNED_SHA256 and PINNED_SHA256[key] != got:
            failures.append(f"{key}: sha256 {got} != pinned {PINNED_SHA256[key]}")
            continue
        if key not in PINNED_SHA256:
            failures.append(f"{key}: not pinned (sha256 {got})")
            continue
        if isinstance(target, str):
            dest = DATA / target
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, dest)
        else:
            print(f"  extract {name}: {extract(path, target)} files", flush=True)
    for listing, rx, template in MEMBER_CHECKS:
        listing = DATA / listing
        if not listing.exists():
            continue
        checked = 0
        for line in listing.read_text().splitlines():
            member, digest = line.split("\t")
            hit = re.match(rx, member)
            if not hit:
                continue
            dest = DATA / template.format(*hit.groups())
            if not dest.exists() or sha256(dest) != digest:
                failures.append(f"{dest.relative_to(DATA)}: missing or sha256 differs from {listing.name}")
            checked += 1
        print(f"  {listing.relative_to(DATA)}: {checked} members checked")
    controlled = DATA / "pon" / "atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn"
    if controlled.exists() and hashlib.md5(controlled.read_bytes()).hexdigest() != "b4a373729d9282ce0a0ed80fd29f6d0a":
        failures.append("controlled PoN md5 differs from b4a373729d9282ce0a0ed80fd29f6d0a")
    if failures:
        sys.exit("FAILED:\n  " + "\n  ".join(failures))
    print("ok")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which not in ("umi", "nonumi", "all"):
        sys.exit(__doc__)
    main(which)
