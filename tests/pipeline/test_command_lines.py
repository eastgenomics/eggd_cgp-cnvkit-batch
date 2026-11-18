"""L2: the CNVkit command lines the app runs (argv contract), recorded for v2.0.2.

These tests recorded the behaviour of v2.0.2 (`main` @ 7559cc2) before any change.
Phase 2 (refactor) must keep them passing unchanged. Phase 4 edits them, one Jira subtask at a time, into the
v2.1.0 contract. The expected argv lists are written out from reading the code, not captured from a run.
"""
import pytest

from apprun import SAMPLE, run_app

F_SAMPLE = "138076232-25178S0009-26CGPV1-9197-F-96770770"
U_SAMPLE = "135024088-25021S0072-26CGPVal1SHv150kb-9999-U-1111"

pytestmark = pytest.mark.pipeline

S = SAMPLE
IMAGE = "cgp-cnvkit:1.0.0"


def expected_calls(parx=False, sex="male", s=None, segment_extra=("--drop-low-coverage",), mask=True, parfix=False):
    """The full CNVkit call sequence for a sample (default S, LIMS M), with --processes from the stubbed nproc (4).

    parx: diploid_parx_genome=grch38 (DI-4097): --diploid-parx-genome on fix, segment, call, diagram, genemetrics.
    sex:  the applied sex (DI-4100): -x on call, diagram and genemetrics; None = no -x (CNVkit infers).
    mask: dropout_mask (DI-4139): segment reads <s>.seg_input.cnr; every other step reads the unmasked <s>.cnr.
    parfix: the PAR post-step ran (DI-4101): genemetrics re-run, same arguments, on the PAR-shifted copies in parfix/.
    """
    S = s or SAMPLE
    px = ("--diploid-parx-genome", "grch38") if parx else ()
    x = ("-x", sex) if sex else ()
    return [
        ["version"],
        ["coverage", "tumour.bam", "panel.bed", "--processes", "4", "--output", f"{S}.targetcoverage.cnn"],
        # DI-4096: edge correction on (no --no-edge), as the PoN was built
        ["fix", f"{S}.targetcoverage.cnn", "empty.antitargetcoverage.cnn", "reference.cnn", *px, "-o", f"{S}.cnr"],
        ["segment", f"{S}.seg_input.cnr" if mask else f"{S}.cnr", "--method", "cbs", "--drop-outliers", "10", *segment_extra, *px,
         "--processes", "4", "-o", f"{S}.cns"],
        # DI-4102: sex checks on raw coverage (chrY presence) and on the .cnr (chrX; post-step decision)
        ["sex", f"{S}.targetcoverage.cnn", *px, "-o", f"{S}.sex_raw.tsv"],
        ["sex", f"{S}.cnr", *px, "-o", f"{S}.sex_cnr.tsv"],
        ["call", f"{S}.cns", *px, *x, "-o", f"{S}.call.cns"],
        # scatter takes neither --diploid-parx-genome nor a sex argument in fc65941
        ["scatter", f"{S}.cnr", "-s", f"{S}.cns", "--trend", "--y-min", "-3", "--y-max", "3",
         "--fig-size", "14", "4", "--title", S, "-o", f"{S}.scatter.png"],
        ["diagram", f"{S}.cnr", "-s", f"{S}.call.cns", "--threshold", "0.5", "--min-probes", "3",
         *px, *x, "--title", S, "-o", f"{S}.diagram.pdf"],
        ["genemetrics", f"{S}.cnr", "-s", f"{S}.call.cns", "-t", "0.3", "--min-probes", "3", *px, *x,
         "-o", f"{S}.genemetrics.tsv"],
    ] + ([gm_call(S, px, x)] if parfix else [])


def gm_call(S, px, x):
    return ["genemetrics", f"{S}.cnr", "-s", f"{S}.call.cns", "-t", "0.3", "--min-probes", "3", *px, *x,
            "-o", f"{S}.genemetrics.tsv"]


PARX_WARNING = "[parx] WARNING: diploid_parx_genome is empty"


ALL_OUTPUTS = {"copy_ratios": f"{S}.cnr", "segments": f"{S}.cns", "call_segments": f"{S}.call.cns",
               "genemetrics": f"{S}.genemetrics.tsv", "scatter_png": f"{S}.scatter.png",
               "diagram_pdf": f"{S}.diagram.pdf",
               # DI-4104
               "genemetrics_noparfix": f"{S}.genemetrics.noparfix.tsv", "sex_qc": f"{S}.sex_qc.tsv",
               "sex_check_mqc": f"{S}.sex_check_mqc.tsv", "dropout_qc": f"{S}.dropout_qc.tsv",
               "target_coverage": f"{S}.targetcoverage.cnn"}


@pytest.mark.parametrize("name, inputs, expected", [
    ("defaults", {}, expected_calls()),
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


# ── DI-4098 / DI-4099 / DI-4100: sex sources, -x on call, genemetrics, diagram (D2, D9) ──

@pytest.mark.parametrize("sample, bam, override, sex, log", [
    (SAMPLE, None, None, "male", "[sex] applied: -x male (source: sample_name; LIMS M)"),
    (F_SAMPLE, None, None, "female", "[sex] applied: -x female (source: sample_name; LIMS F)"),
    (U_SAMPLE, None, None, None, "[sex] applied: no -x, CNVkit infers the sex (source: cnvkit_inference; LIMS U)"),
    (U_SAMPLE, None, "male", "male", "[sex] applied: -x male (source: override; LIMS U)"),
    (SAMPLE, None, "F", "female", "[sex] applied: -x female (source: override; LIMS M)"),
    ("25328S0023", f"{SAMPLE}.bam", None, "male", "[sex] applied: -x male (source: bam_name; LIMS M)"),
    ("25328S0023", "25328S0023.bam", None, None,
     "[sex] applied: no -x, CNVkit infers the sex (source: cnvkit_inference; LIMS U)"),
])
def test_sex_resolution(tmp_path, sample, bam, override, sex, log):
    from apprun import default_inputs
    inputs, names = default_inputs(sample_id=sample, bam_name=bam)
    if override is not None:
        inputs["sex_override"] = override
    run = run_app(tmp_path, inputs, names)
    assert run.returncode == 0, run.stderr
    assert run.cnvkit == expected_calls(sex=sex, s=sample)
    assert log in run.stdout


@pytest.mark.parametrize("bad", ["X", "unknown", "U"])
def test_bad_sex_override_fails_before_analysis(tmp_path, bad):
    run = run_app(tmp_path, {"sex_override": bad})
    assert run.returncode != 0
    assert "not accepted" in run.stdout + run.stderr
    assert run.cnvkit == [["version"]]
    assert run.outputs == {}


def test_sample_sex_deprecated_and_ignored(tmp_path):
    # sample_sex is wired from PURPLE's inferred sex in the current workflow: it must not reach CNVkit (D2)
    from apprun import default_inputs
    inputs, names = default_inputs(sample_id=U_SAMPLE)
    inputs["sample_sex"] = "male"
    run = run_app(tmp_path, inputs, names)
    assert run.returncode == 0, run.stderr
    assert run.cnvkit == expected_calls(sex=None, s=U_SAMPLE)
    assert "[deprecated] sample_sex supplied but ignored" in run.stdout


# ── DI-4139: dropout mask (D7) ───────────────────────────────────────────────

QC_HEADER = "arm\tchromosome\tstart\tend\tgene\tlog2\traw_reads\treason\n"


def masked_fixtures(tmp_path):
    """A .cnr with one isolated deep autosomal probe (0 reads) and a deep chrX probe (never masked)."""
    cnr = tmp_path / "fix.cnr"
    cnr.write_text("chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n"
                   "chr12\t11873311\t11873583\tETV6\t-0.95\t292.4\t0.886\n"
                   "chr12\t11873583\t11873855\tETV6\t-28.7962\t0\t0.743\n"
                   "chr12\t11873855\t11874127\tETV6\t0.1\t300\t0.9\n"
                   "chrX\t5000000\t5000144\tXGENE\t-20\t0\t0.5\n")
    cnn = tmp_path / "coverage.cnn"
    cnn.write_text("chromosome\tstart\tend\tgene\tdepth\tlog2\n"
                   "chr12\t11873311\t11873583\tETV6\t292.4\t8.2\n"
                   "chr12\t11873583\t11873855\tETV6\t0\t-20\n"
                   "chr12\t11873855\t11874127\tETV6\t300\t8.2\n"
                   "chrX\t5000000\t5000144\tXGENE\t0\t-20\n")
    return {"STUB_OUTPUT_FIX": str(cnr), "STUB_OUTPUT_COVERAGE": str(cnn)}


def test_dropout_mask_default_on(tmp_path):
    run = run_app(tmp_path, stub_env=masked_fixtures(tmp_path))
    assert run.returncode == 0, run.stderr
    assert run.cnvkit == expected_calls(mask=True)
    seg_input = (run.workdir / f"{S}.seg_input.cnr").read_text().splitlines()
    assert seg_input[2].split("\t")[6] == "0"                       # the isolated ETV6 probe: weight 0
    assert [r.split("\t")[6] for r in seg_input[1:]] == ["0.886", "0", "0.9", "0.5"]
    # the .cnr every other step reads (and copy_ratios) is the unmasked file
    assert (run.workdir / f"{S}.cnr").read_text() == open(masked_fixtures(tmp_path)["STUB_OUTPUT_FIX"]).read()
    assert (run.workdir / f"{S}.dropout_qc.tsv").read_text() == (
        QC_HEADER + "M\tchr12\t11873583\t11873855\tETV6\t-28.7962\t0.0\tisolated_run1\n")
    assert "[dropout] masked probes: 1" in run.stdout


def test_dropout_mask_off(tmp_path):
    run = run_app(tmp_path, {"dropout_mask": False}, stub_env=masked_fixtures(tmp_path))
    assert run.returncode == 0, run.stderr
    assert run.cnvkit == expected_calls(mask=False)
    assert not (run.workdir / f"{S}.seg_input.cnr").exists()
    assert (run.workdir / f"{S}.dropout_qc.tsv").read_text() == QC_HEADER
    assert ("[dropout] dropout_mask=false: segment reads the unmasked .cnr (validated final configuration\n"
            "[dropout] without the mask, arm A0; not v2.0.x)") in run.stdout


def test_only_segment_reads_the_masked_file(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    readers = [c[0] for c in run.cnvkit if any(a.endswith(".seg_input.cnr") for a in c)]
    assert readers == ["segment"]


# ── DI-4101: PAR post-step (D3, D11) ─────────────────────────────────────────

GM_H = "gene\tchromosome\tstart\tend\tlog2\tdepth\tweight\tcn\tprobes\tsegment_weight\tsegment_probes\n"
TP53 = "TP53\tchr17\t7661779\t7687538\t-0.62\t400\t5\t1\t10\t5\t10\n"
CRLF2 = "CRLF2\tchrX\t1190891\t1212654\t1.26805\t679\t7.2\t3\t8\t7.2\t8\n"
CALL_H = "chromosome\tstart\tend\tgene\tlog2\tcn\tdepth\tprobes\tweight\n"


def parfix_fixtures(tmp_path, parfix_rows=(TP53,)):
    """Original genemetrics with CRLF2 at +1 (male shift); the re-run in parfix/ drops it (0.268 < 0.3)."""
    (tmp_path / "gm.tsv").write_text(GM_H + TP53 + CRLF2)
    (tmp_path / "gm_parfix.tsv").write_text(GM_H + "".join(parfix_rows))
    (tmp_path / "call.cns").write_text(CALL_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                                                "chrX\t1190891\t2781479\tCRLF2\t0.26805\t2\t679\t8\t7.2\n")
    return {"STUB_OUTPUT_GENEMETRICS": str(tmp_path / "gm.tsv"),
            "STUB_OUTPUT_GENEMETRICS_PARFIX": str(tmp_path / "gm_parfix.tsv"),
            "STUB_OUTPUT_CALL": str(tmp_path / "call.cns")}


PARX = {"diploid_parx_genome": "grch38"}


@pytest.mark.parametrize("sample, cnr_sex, override, applied", [
    (SAMPLE, "Female", None, True),        # -x male: post-step, whatever the .cnr says
    (F_SAMPLE, "Male", None, False),       # -x female: never
    (U_SAMPLE, "Male", None, True),        # no -x: CNVkit's own guess on the .cnr decides
    (U_SAMPLE, "Female", None, False),
    (U_SAMPLE, "Female", "male", True),    # override male
])
def test_parfix_decision(tmp_path, sample, cnr_sex, override, applied):
    from apprun import default_inputs
    inputs, names = default_inputs(sample_id=sample)
    inputs.update(PARX)
    if override:
        inputs["sex_override"] = override
    env = dict(parfix_fixtures(tmp_path), STUB_CNR_SEX=cnr_sex)
    run = run_app(tmp_path, inputs, names, stub_env=env)
    assert run.returncode == 0, run.stderr
    sex = {SAMPLE: "male", F_SAMPLE: "female"}.get(sample) if not override else override
    assert run.cnvkit == expected_calls(parx=True, sex=sex, s=sample, parfix=applied)
    w = run.workdir
    noparfix = (w / f"{sample}.genemetrics.noparfix.tsv").read_text()
    final = (w / f"{sample}.genemetrics.tsv").read_text()
    assert noparfix == GM_H + TP53 + CRLF2                     # always the uncorrected read-out (D3)
    assert final == (GM_H + TP53 if applied else noparfix)
    assert ("[parfix] post-step applied" in run.stdout) is applied
    if applied:
        rerun = run.docker[-1]                                     # the genemetrics re-run is the last call
        assert rerun["cmd"][1] == "genemetrics"
        assert rerun["options"] == docker_options(w / "parfix")
        assert (w / "parfix" / f"{sample}.cnr").exists() and (w / "parfix" / f"{sample}.call.cns").exists()
        assert "par_check=1" in run.stdout


def test_parfix_not_applicable_without_parx(tmp_path):
    run = run_app(tmp_path, stub_env=parfix_fixtures(tmp_path))          # LIMS M, no diploid_parx_genome
    assert run.returncode == 0, run.stderr
    assert run.cnvkit == expected_calls(parx=False, parfix=False)
    assert (run.workdir / f"{S}.genemetrics.noparfix.tsv").read_bytes() == \
           (run.workdir / f"{S}.genemetrics.tsv").read_bytes()
    assert "[parfix] not applied: no --diploid-parx-genome" in run.stdout


def test_parfix_non_par_difference_fails_the_job(tmp_path):
    env = parfix_fixtures(tmp_path, parfix_rows=(TP53.replace("-0.62", "-0.70"),))
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode != 0
    assert "non-PAR rows differ" in run.stderr
    assert run.outputs == {}


# ── DI-4103: sex QC (D8) ─────────────────────────────────────────────────────

def sex_qc_row(run, sample):
    import csv
    with open(run.workdir / f"{sample}.sex_qc.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    return rows[0]


@pytest.mark.parametrize("sample, override, raw_y, cnr_sex, expect", [
    (SAMPLE, None, "-1.13", "Male",
     dict(lims_sex="M", applied_sex="male", applied_source="sample_name", chrY_sex="Male", sex_mismatch="0",
          par_fix_applied="1")),
    (SAMPLE, None, "-22.0", "Female",                 # LIMS M but no chrY: flagged
     dict(lims_sex="M", applied_sex="male", chrY_sex="Female", sex_mismatch="1", par_fix_applied="1")),
    (SAMPLE, "female", "-1.13", "Male",               # override female on a chrY-present sample: flagged
     dict(lims_sex="M", applied_sex="female", applied_source="override", sex_mismatch="1", par_fix_applied="0")),
    (U_SAMPLE, None, "-1.13", "Male",
     dict(lims_sex="U", applied_sex="cnvkit_inferred(Male)", applied_source="cnvkit_inference",
          sex_mismatch="NA", par_fix_applied="1")),
    (U_SAMPLE, None, "-25.0", "Female",
     dict(lims_sex="U", applied_sex="cnvkit_inferred(Female)", sex_mismatch="NA", par_fix_applied="0")),
])
def test_sex_qc(tmp_path, sample, override, raw_y, cnr_sex, expect):
    from apprun import default_inputs
    inputs, names = default_inputs(sample_id=sample)
    inputs.update(PARX)
    if override:
        inputs["sex_override"] = override
    env = dict(masked_fixtures(tmp_path), STUB_RAW_Y=raw_y, STUB_CNR_SEX=cnr_sex)
    run = run_app(tmp_path, inputs, names, stub_env=env)
    assert run.returncode == 0, run.stderr
    row = sex_qc_row(run, sample)
    for k, v in expect.items():
        assert row[k] == v, k
    assert row["dropout_probes"] == "1" and row["par_span_segments"] == "0" and row["run"] == ""
    mqc = (run.workdir / f"{sample}.sex_check_mqc.tsv").read_text().splitlines()
    assert mqc[0] == "# id: 'cnvkit_sex_check'"
    assert mqc[5].split("\t")[0] == sample and mqc[5].split("\t")[-1] == "1"
    assert "[sex] chrY " in run.stdout


def test_sex_qc_failure_is_fatal(tmp_path):
    # a sex table for another sample (e.g. a mis-named file) must stop the job
    other = tmp_path / "other_sex.tsv"
    other.write_text("sample\tsex\tX_logratio\tY_logratio\nother.targetcoverage.cnn\tMale\t-1\t-1\n")
    run = run_app(tmp_path, stub_env={"STUB_OUTPUT_SEX": str(other)})
    assert run.returncode != 0
    assert run.outputs == {}


def test_uploaded_files_are_the_final_versions(tmp_path):
    # genemetrics = the post-step read-out; copy_ratios = the unmasked .cnr
    env = dict(parfix_fixtures(tmp_path), **masked_fixtures(tmp_path))
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode == 0, run.stderr
    assert (run.workdir / f"{S}.genemetrics.tsv").read_text() == GM_H + TP53
    assert run.uploads.count(f"{S}.cnr") == 1 and f"{S}.seg_input.cnr" not in run.uploads
    assert "=== Done:" in run.stdout


def test_uploads_wait_for_finalisation(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    uploads = [c for c in run.dx if c[0] == "upload"]
    assert len(uploads) == len(run.outputs) == len(ALL_OUTPUTS)
    assert all(c[2:] == ["--wait", "--brief"] for c in uploads)


def test_missing_output_file_fails_the_job(tmp_path):
    # a failed dx upload must stop the job, not add an empty output
    run = run_app(tmp_path, stub_env={"STUB_SKIP_OUTPUT": "sex_qc"})
    assert run.returncode != 0
    assert "sex_qc" not in run.outputs


def test_start_from_coverage(tmp_path):
    from apprun import run_pipeline_direct
    cnn = "chromosome\tstart\tend\tgene\tdepth\tlog2\nchr1\t100\t200\tGENE1\t500.0\t9.0\n"
    run = run_pipeline_direct(tmp_path, {"SAMPLE_ID": S, "START_FROM": "coverage", "NPROC": "4",
                                         "DIPLOID_PARX_GENOME": "grch38"},
                              files={f"{S}.targetcoverage.cnn": cnn, "reference.cnn": "stub\n"})
    assert run.returncode == 0, run.stdout + run.stderr
    assert run.cnvkit == expected_calls(parx=True, parfix=True)[2:]      # no version, no coverage
    assert (run.workdir / f"{S}.targetcoverage.cnn").read_text() == cnn
    assert not (run.workdir / "panel.bed").exists()


def test_start_from_coverage_needs_the_file(tmp_path):
    from apprun import run_pipeline_direct
    run = run_pipeline_direct(tmp_path, {"SAMPLE_ID": S, "START_FROM": "coverage"}, files={"reference.cnn": "x\n"})
    assert run.returncode != 0
    assert "targetcoverage.cnn is missing" in run.stderr
    assert "ERROR:" not in run.stdout
    assert run.cnvkit == []


def test_code_sh_always_starts_from_the_bam(tmp_path):
    # a stray START_FROM in the job environment must not skip coverage
    run = run_app(tmp_path, stub_env={"START_FROM": "coverage"})
    assert run.returncode == 0, run.stderr
    assert run.subcommand("coverage")


@pytest.mark.parametrize("bad", ["a/b", "..", ".", "-starts-with-dash", "has space", ""])
def test_bad_sample_id_fails_before_analysis(tmp_path, bad):
    run = run_app(tmp_path, {"sample_id": bad})
    assert run.returncode != 0
    assert "ERROR: sample_id" in run.stderr
    assert "must not start with" in run.stderr
    assert "ERROR:" not in run.stdout
    assert run.cnvkit == []


def test_version():
    from apprun import DXAPP
    assert DXAPP["version"] == "2.1.0" and DXAPP["properties"]["githubRelease"] == "2.1.0"


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
    assert run.cnvkit == expected_calls(parx=True, parfix=True)
    assert PARX_WARNING not in run.stdout
    assert "[parx] --diploid-parx-genome grch38" in run.stdout


def test_parx_empty_warns_and_omits_the_flag(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert not any("--diploid-parx-genome" in c for c in run.cnvkit)
    assert PARX_WARNING in run.stdout
    assert "raw-coverage sex QC can call male samples female" in run.stdout
    assert "inference calls every sample female" not in run.stdout


def test_parx_rejects_other_values(tmp_path):
    # dxapp.json restricts the value with choices; the app checks it too, before any analysis step
    run = run_app(tmp_path, {"diploid_parx_genome": "hg19"})
    assert run.returncode != 0
    assert "ERROR: diploid_parx_genome 'hg19' is not supported" in run.stderr
    assert "ERROR:" not in run.stdout
    assert run.cnvkit == [["version"]]


def docker_options(workdir):
    """CNVkit runs as the executor's user, so the host can overwrite container-written files (PAR post-step);
    HOME and the Matplotlib cache point at a writable directory for a non-root user."""
    import os
    return ["--rm", "--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp", "-e", "MPLCONFIGDIR=/tmp",
            "-v", f"{workdir}:/work", "-w", "/work"]


def test_docker_invocation(tmp_path):
    run = run_app(tmp_path)
    assert run.returncode == 0, run.stderr
    assert {d["image"] for d in run.docker} == {IMAGE}
    for d in run.docker:
        assert d["options"] == docker_options(run.workdir)
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
    assert "does not match BAM filename" in run.stderr
    assert "Expected" in run.stderr
    assert "ERROR:" not in run.stdout
    assert run.cnvkit == [["version"]]       # fails before any analysis step
    assert run.outputs == {}


@pytest.mark.parametrize("failing, missing", [("scatter", "scatter_png"), ("diagram", "diagram_pdf")])
def test_plot_failure_is_not_fatal(tmp_path, failing, missing):
    run = run_app(tmp_path, fail=(failing,))
    assert run.returncode == 0, run.stderr
    assert f"[{failing}] WARNING" in run.stdout
    assert set(run.outputs) == set(ALL_OUTPUTS) - {missing}


@pytest.mark.parametrize("failing", ["coverage", "fix", "segment", "sex", "call", "genemetrics"])
def test_analysis_failure_is_fatal(tmp_path, failing):
    run = run_app(tmp_path, fail=(failing,))
    assert run.returncode != 0
    assert run.outputs == {}
