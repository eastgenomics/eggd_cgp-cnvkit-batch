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


def expected_calls(parx=False, segment_extra=("--drop-low-coverage",), call_extra=(), diagram_extra=()):
    """The full CNVkit call sequence for sample S, with --processes from the stubbed nproc (4).

    parx: diploid_parx_genome=grch38 (DI-4097): --diploid-parx-genome on fix, segment, call, diagram, genemetrics.
    """
    px = ("--diploid-parx-genome", "grch38") if parx else ()
    return [
        ["version"],
        ["coverage", "tumour.bam", "panel.bed", "--processes", "4", "--output", f"{S}.targetcoverage.cnn"],
        # DI-4096: edge correction on (no --no-edge), as the PoN was built
        ["fix", f"{S}.targetcoverage.cnn", "empty.antitargetcoverage.cnn", "reference.cnn", *px, "-o", f"{S}.cnr"],
        ["segment", f"{S}.cnr", "--method", "cbs", "--drop-outliers", "10", *segment_extra, *px,
         "--processes", "4", "-o", f"{S}.cns"],
        ["call", f"{S}.cns", *call_extra, *px, "-o", f"{S}.call.cns"],
        # scatter takes neither --diploid-parx-genome nor a sex argument in fc65941
        ["scatter", f"{S}.cnr", "-s", f"{S}.cns", "--trend", "--y-min", "-3", "--y-max", "3",
         "--fig-size", "14", "4", "--title", S, "-o", f"{S}.scatter.png"],
        ["diagram", f"{S}.cnr", "-s", f"{S}.call.cns", "--threshold", "0.5", "--min-probes", "3",
         *diagram_extra, *px, "--title", S, "-o", f"{S}.diagram.pdf"],
        ["genemetrics", f"{S}.cnr", "-s", f"{S}.call.cns", "-t", "0.3", "--min-probes", "3", *px,
         "-o", f"{S}.genemetrics.tsv"],
    ]


PARX_WARNING = "[parx] WARNING: diploid_parx_genome is empty"


ALL_OUTPUTS = {"copy_ratios": f"{S}.cnr", "segments": f"{S}.cns", "call_segments": f"{S}.call.cns",
               "genemetrics": f"{S}.genemetrics.tsv", "scatter_png": f"{S}.scatter.png",
               "diagram_pdf": f"{S}.diagram.pdf"}


@pytest.mark.parametrize("name, inputs, expected", [
    ("defaults", {}, expected_calls()),
    ("sample_sex", {"sample_sex": "male"}, expected_calls(diagram_extra=("--sample-sex", "male"))),
    ("no_drop_low_coverage", {"drop_low_coverage": False}, expected_calls(segment_extra=())),
])
def test_command_lines(tmp_path, name, inputs, expected):
    run = run_app(tmp_path, inputs)
    assert run.returncode == 0, run.stderr
    assert run.cnvkit == expected


# ── DI-4121: no purity/ploidy; the inputs are deprecated and ignored (D4) ──

PURITY_WARNING = "[deprecated] purity/ploidy supplied but ignored"


@pytest.mark.parametrize("inputs", [{"purity": 0.5, "ploidy": 3.1}, {"purity": 0.62}, {"purity": 0.25},
                                    {"ploidy": 2}])
def test_purity_ploidy_ignored(tmp_path, inputs):
    run = run_app(tmp_path, inputs)
    assert run.returncode == 0, run.stderr
    assert run.cnvkit == expected_calls()
    assert PURITY_WARNING in run.stdout


def test_no_purity_no_warning(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert PURITY_WARNING not in run.stdout


def test_ploidy_has_no_default():
    # a default would make DNAnexus set ploidy on every job, and the deprecation warning would always fire
    from apprun import DXAPP
    spec = {i["name"]: i for i in DXAPP["inputSpec"]}
    assert "default" not in spec["ploidy"] and spec["ploidy"]["optional"] is True
    assert spec["purity"]["optional"] is True


# ── DI-4097: diploid_parx_genome ─────────────────────────────────────────────

def test_parx_grch38(tmp_path):
    run = run_app(tmp_path, {"diploid_parx_genome": "grch38"})
    assert run.returncode == 0, run.stderr
    assert run.cnvkit == expected_calls(parx=True)
    assert PARX_WARNING not in run.stdout
    assert "[parx] --diploid-parx-genome grch38" in run.stdout


def test_parx_empty_warns_and_omits_the_flag(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert not any("--diploid-parx-genome" in c for c in run.cnvkit)
    assert PARX_WARNING in run.stdout


def test_parx_rejects_other_values(tmp_path):
    # dxapp.json restricts the value with choices; the app checks it too, before any analysis step
    run = run_app(tmp_path, {"diploid_parx_genome": "hg19"})
    assert run.returncode != 0
    assert "diploid_parx_genome" in run.stdout + run.stderr
    assert run.cnvkit == [["version"]]


def test_docker_invocation(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert {d["image"] for d in run.docker} == {IMAGE}
    for d in run.docker:
        assert d["options"] == ["--rm", "-v", f"{run.workdir}:/work", "-w", "/work"]
    # the image is loaded from its pinned DNAnexus file
    assert ["download", "project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8j7Vyj45FG1BbK26JQgQY6q",
            "-o", "cnvkit-image.tar.gz"] in run.dx


def test_outputs(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert {k: v["value"] for k, v in run.outputs.items()} == {k: f"file-{v}" for k, v in ALL_OUTPUTS.items()}
    assert {v["class"] for v in run.outputs.values()} == {"file"}


def test_bed_deduplicated(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert (run.workdir / "panel.bed").read_text().splitlines() == [
        "chr1\t100\t200\tGENE1", "chr2\t100\t200\tGENE2", "chrX\t1190891\t1212654\tCRLF2"]


@pytest.mark.parametrize("bai_name", [f"{S}.bam.bai", f"{S}.bai"])
def test_bai_name_accepted(tmp_path, bai_name):
    run = run_app(tmp_path, names={"file-BAI": bai_name})
    assert run.returncode == 0, run.stderr


def test_bai_name_mismatch_fails(tmp_path):
    run = run_app(tmp_path, names={"file-BAI": "other-sample.bam.bai"})
    assert run.returncode != 0
    assert "does not match BAM filename" in run.stdout
    assert run.cnvkit == [["version"]]       # fails before any analysis step
    assert run.outputs == {}


@pytest.mark.parametrize("failing, missing", [("scatter", "scatter_png"), ("diagram", "diagram_pdf")])
def test_plot_failure_is_not_fatal(tmp_path, failing, missing):
    run = run_app(tmp_path, fail=(failing,))
    assert run.returncode == 0, run.stderr
    assert f"[{failing}] WARNING" in run.stdout
    assert set(run.outputs) == set(ALL_OUTPUTS) - {missing}


@pytest.mark.parametrize("failing", ["coverage", "fix", "segment", "call", "genemetrics"])
def test_analysis_failure_is_fatal(tmp_path, failing):
    run = run_app(tmp_path, fail=(failing,))
    assert run.returncode != 0
    assert run.outputs == {}
