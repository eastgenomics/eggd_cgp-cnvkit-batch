"""L2: the CNVkit command lines the app runs (argv contract), recorded for v2.0.2.

These tests fix the behaviour of v2.0.2 (`main` @ 7559cc2) before any change (PLAN_v2.1.0_build.md, Phase 1).
Phase 2 (refactor) must keep them passing unchanged. Phase 4 edits them, one Jira subtask at a time, into the
v2.1.0 contract. The expected argv lists are written out from reading the code, not captured from a run.
"""
import pytest

from apprun import SAMPLE, run_app

pytestmark = pytest.mark.pipeline

S = SAMPLE
IMAGE = "cgp-cnvkit:1.0.0"


def v202_calls(segment_extra=("--drop-low-coverage",), call_extra=(), diagram_extra=()):
    """The full v2.0.2 CNVkit call sequence for sample S, with --processes from the stubbed nproc (4)."""
    return [
        ["version"],
        ["coverage", "tumour.bam", "panel.bed", "--processes", "4", "--output", f"{S}.targetcoverage.cnn"],
        # DI-4096: edge correction on (no --no-edge), as the PoN was built
        ["fix", f"{S}.targetcoverage.cnn", "empty.antitargetcoverage.cnn", "reference.cnn", "-o", f"{S}.cnr"],
        ["segment", f"{S}.cnr", "--method", "cbs", "--drop-outliers", "10", *segment_extra,
         "--processes", "4", "-o", f"{S}.cns"],
        ["call", f"{S}.cns", *call_extra, "-o", f"{S}.call.cns"],
        ["scatter", f"{S}.cnr", "-s", f"{S}.cns", "--trend", "--y-min", "-3", "--y-max", "3",
         "--fig-size", "14", "4", "--title", S, "-o", f"{S}.scatter.png"],
        ["diagram", f"{S}.cnr", "-s", f"{S}.call.cns", "--threshold", "0.5", "--min-probes", "3",
         *diagram_extra, "--title", S, "-o", f"{S}.diagram.pdf"],
        ["genemetrics", f"{S}.cnr", "-s", f"{S}.call.cns", "-t", "0.3", "--min-probes", "3",
         "-o", f"{S}.genemetrics.tsv"],
    ]


ALL_OUTPUTS = {"copy_ratios": f"{S}.cnr", "segments": f"{S}.cns", "call_segments": f"{S}.call.cns",
               "genemetrics": f"{S}.genemetrics.tsv", "scatter_png": f"{S}.scatter.png",
               "diagram_pdf": f"{S}.diagram.pdf"}


@pytest.mark.parametrize("name, inputs, expected", [
    ("defaults", {}, v202_calls()),
    ("purity_above_gate", {"purity": 0.5, "ploidy": 3.1},
     v202_calls(call_extra=("--purity", "0.5", "--ploidy", "3.1"))),
    ("purity_default_ploidy", {"purity": 0.62},
     v202_calls(call_extra=("--purity", "0.62", "--ploidy", "2"))),
    ("purity_at_gate", {"purity": 0.4, "ploidy": 3.1}, v202_calls()),
    ("purity_below_gate", {"purity": 0.25}, v202_calls()),
    ("sample_sex", {"sample_sex": "male"}, v202_calls(diagram_extra=("--sample-sex", "male"))),
    ("no_drop_low_coverage", {"drop_low_coverage": False}, v202_calls(segment_extra=())),
])
def test_v202_command_lines(tmp_path, name, inputs, expected):
    run = run_app(tmp_path, inputs)
    assert run.returncode == 0, run.stderr
    assert run.cnvkit == expected


def test_v202_docker_invocation(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert {d["image"] for d in run.docker} == {IMAGE}
    for d in run.docker:
        assert d["options"] == ["--rm", "-v", f"{run.workdir}:/work", "-w", "/work"]
    # the image is loaded from its pinned DNAnexus file
    assert ["download", "project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8j7Vyj45FG1BbK26JQgQY6q",
            "-o", "cnvkit-image.tar.gz"] in run.dx


def test_v202_outputs(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert {k: v["value"] for k, v in run.outputs.items()} == {k: f"file-{v}" for k, v in ALL_OUTPUTS.items()}
    assert {v["class"] for v in run.outputs.values()} == {"file"}


def test_v202_bed_deduplicated(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert (run.workdir / "panel.bed").read_text().splitlines() == [
        "chr1\t100\t200\tGENE1", "chr2\t100\t200\tGENE2", "chrX\t1190891\t1212654\tCRLF2"]


@pytest.mark.parametrize("bai_name", [f"{S}.bam.bai", f"{S}.bai"])
def test_v202_bai_name_accepted(tmp_path, bai_name):
    run = run_app(tmp_path, names={"file-BAI": bai_name})
    assert run.returncode == 0, run.stderr


def test_v202_bai_name_mismatch_fails(tmp_path):
    run = run_app(tmp_path, names={"file-BAI": "other-sample.bam.bai"})
    assert run.returncode != 0
    assert "does not match BAM filename" in run.stdout
    assert run.cnvkit == [["version"]]       # fails before any analysis step
    assert run.outputs == {}


@pytest.mark.parametrize("failing, missing", [("scatter", "scatter_png"), ("diagram", "diagram_pdf")])
def test_v202_plot_failure_is_not_fatal(tmp_path, failing, missing):
    run = run_app(tmp_path, fail=(failing,))
    assert run.returncode == 0, run.stderr
    assert f"[{failing}] WARNING" in run.stdout
    assert set(run.outputs) == set(ALL_OUTPUTS) - {missing}


@pytest.mark.parametrize("failing", ["coverage", "fix", "segment", "call", "genemetrics"])
def test_v202_analysis_failure_is_fatal(tmp_path, failing):
    run = run_app(tmp_path, fail=(failing,))
    assert run.returncode != 0
    assert run.outputs == {}
