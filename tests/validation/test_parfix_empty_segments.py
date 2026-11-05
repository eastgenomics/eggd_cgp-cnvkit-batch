"""L3 (F2, round 2): the PAR post-step on an empty segmentation, in the cgp-cnvkit:1.0.0 image.

A sample whose bins all have zero depth is valid input that CNVkit handles: `segment --drop-low-coverage` drops every
bin and writes a header-only .cns, `call` a header-only .call.cns, and genemetrics, given no segments, falls back to
per-gene bin means (CNVkit fc65941 reports.do_genemetrics: `if segments:` is False, so gene_metrics_by_gene). The
commands are pipeline.sh's (-x male, --diploid-parx-genome grch38). The base (8b683c8, whose compare is the reference
genemetrics_parfix.compare) reported par_check 0 as a warning and published the corrected table; the port must
publish it too (exit 0, par_check bins), not fail the job: with no segment to judge them by, the PAR rows are not
checked (F4, round 3), the headers and non-PAR rows are. The reviewer's F4 case adds three overlapping zero-depth bins,
2781000-2781300, 2781100-2781600 (spanning the PAR1 end, so not shifted) and 2781200-2781400, of one gene: its row
lies wholly in PAR1 (first start, last end) and moves from -9 to -9.66667, not by 1.0. Synthetic input: no validation
data is needed.
"""
import os
import shlex
import uuid

import pytest

from procutil import remove_container, run_group
from vdata import REFERENCE, load_module, script

pytestmark = [pytest.mark.validation, pytest.mark.docker]

IMAGE = "cgp-cnvkit:1.0.0"
IMAGE_TIMEOUT = 600
PX = ["--diploid-parx-genome", "grch38"]
GM = ["-t", "0.3", "--min-probes", "3", *PX, "-x", "male"]     # pipeline.sh gm_args for a LIMS male


OVERLAP = [(2781000, 2781300), (2781100, 2781600), (2781200, 2781400)]   # the reviewer's F4 bins (gene OVL)


def cnr(overlap=False):
    """Six zero-depth bins (log2 -10) per gene: two autosomal, CRLF2 wholly in PAR1, one non-PAR chrX; with overlap,
    also OVL's three overlapping bins at the PAR1 end."""
    rows = ["chromosome\tstart\tend\tgene\tlog2\tdepth\tweight"]
    for chrom, start, gene in [("chr1", 1000000, "GENEA"), ("chr2", 2000000, "GENEB"),
                               ("chrX", 1190891, "CRLF2"), ("chrX", 5000000, "NONPAR")]:
        rows += [f"{chrom}\t{start + i * 200}\t{start + (i + 1) * 200}\t{gene}\t-10\t0\t0.5" for i in range(6)]
        if overlap and gene == "CRLF2":
            rows += [f"chrX\t{a}\t{b}\tOVL\t-10\t0\t0.5" for a, b in OVERLAP]
    return "\n".join(rows) + "\n"


def in_image(workdir, commands):
    name = f"cnvkit-l3-nosegs-{uuid.uuid4().hex}"
    cmd = ["docker", "run", "--rm", "--name", name, "-u", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp",
           "-e", "MPLCONFIGDIR=/tmp", "-v", f"{workdir}:{workdir}", "-w", str(workdir), IMAGE,
           "bash", "-c", "set -e; " + "; ".join(shlex.join(c) for c in commands)]
    proc = run_group(cmd, IMAGE_TIMEOUT, cleanup=lambda: remove_container(name), capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr[-3000:]


@pytest.mark.parametrize("overlap", [False, True], ids=["plain", "overlapping_bins"])
def test_parfix_on_an_empty_segmentation(tmp_path, overlap):
    s = "s"
    (tmp_path / f"{s}.cnr").write_text(cnr(overlap))
    in_image(tmp_path, [
        ["cnvkit.py", "segment", f"{s}.cnr", "--method", "cbs", "--drop-outliers", "10", "--drop-low-coverage", *PX,
         "--processes", "1", "-o", f"{s}.cns"],
        ["cnvkit.py", "sex", f"{s}.cnr", *PX, "-o", f"{s}.sex_cnr.tsv"],
        ["cnvkit.py", "call", f"{s}.cns", *PX, "-x", "male", "-o", f"{s}.call.cns"],
        ["cnvkit.py", "genemetrics", f"{s}.cnr", "-s", f"{s}.call.cns", *GM, "-o", f"{s}.genemetrics.noparfix.tsv"],
    ])
    assert (tmp_path / f"{s}.call.cns").read_text().count("\n") == 1                # no segments
    before = (tmp_path / f"{s}.genemetrics.noparfix.tsv").read_text()
    assert before.startswith("gene\tchromosome\tstart\tend\tlog2\tdepth\tweight\tprobes\n")   # bin means
    assert "CRLF2\tchrX\t1190891\t1192091\t-9\t0\t3\t6\n" in before
    assert ("OVL\tchrX\t2781000\t2781400\t-9\t0\t1.5\t3\n" in before) is overlap

    decide = script("cnvkit_parfix.py", "decide", "--sex-arg", "male", "--cnr-sex", f"{s}.sex_cnr.tsv",
                    "--parx", "grch38", cwd=tmp_path)
    assert decide.returncode == 0 and decide.stdout.startswith("1\t"), decide.stderr
    shift = script("cnvkit_parfix.py", "shift", f"{s}.cnr", f"{s}.call.cns", "parfix", cwd=tmp_path)
    assert shift.returncode == 0, shift.stderr
    assert shift.stdout.strip() == ("par_bins=8 span_bins=1" if overlap else "par_bins=6 span_bins=0") + \
        " par_segs=0 span_segs=0", shift.stderr
    in_image(tmp_path / "parfix", [["cnvkit.py", "genemetrics", f"{s}.cnr", "-s", f"{s}.call.cns", *GM,
                                    "-o", f"{s}.genemetrics.tsv"]])
    after = (tmp_path / "parfix" / f"{s}.genemetrics.tsv").read_text()
    assert after == before.replace("CRLF2\tchrX\t1190891\t1192091\t-9\t", "CRLF2\tchrX\t1190891\t1192091\t-10\t") \
        .replace("OVL\tchrX\t2781000\t2781400\t-9\t", "OVL\tchrX\t2781000\t2781400\t-9.66667\t")

    # HEAD: published, the PAR rows not checked (par_check bins; pipeline.sh warns)
    proc = script("cnvkit_parfix.py", "compare", f"{s}.genemetrics.noparfix.tsv", f"parfix/{s}.genemetrics.tsv",
                  f"{s}.call.cns", f"{s}.cnr", cwd=tmp_path)
    assert proc.returncode == 0 and proc.stdout.strip().endswith("par_check=bins"), proc.stderr
    # the base: par_check 0, which it only logged
    ref = load_module(REFERENCE / "genemetrics_parfix.py", "genemetrics_parfix_ref_nosegs")
    rec = ref.compare(dict(action="dry", rundir=str(tmp_path), sample=s, orig=str(tmp_path / f"{s}.genemetrics."
                           "noparfix.tsv"), new=str(tmp_path / "parfix" / f"{s}.genemetrics.tsv")))
    assert rec["par_check"] == 0
