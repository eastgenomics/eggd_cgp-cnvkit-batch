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
# cn 2: its call.cns segment's cn, which genemetrics copies into the row (the .cnr has no cn column)
CRLF2 = "CRLF2\tchrX\t1190891\t1212654\t1.26805\t679\t7.2\t2\t8\t7.2\t8\n"
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


def test_parfix_unexplained_par_discrepancy_fails_the_job(tmp_path):
    # CRLF2 still carries the +1 after the fix and no segment spans a PAR boundary: par_check 0 is fatal (A1)
    env = parfix_fixtures(tmp_path, parfix_rows=(TP53, CRLF2.replace("1.26805", "1.3")))
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode != 0
    assert "FAILED: par_check 0" in run.stderr and "CRLF2@chrX:1190891-1212654" in run.stderr
    assert run.outputs == {} and run.uploads == []


def test_parfix_span_elsewhere_does_not_excuse_a_par_row(tmp_path):
    # a segment spans the PAR2 boundary, but CRLF2's segment lies wholly in PAR1 and its row is unchanged (+1 kept)
    env = parfix_fixtures(tmp_path, parfix_rows=(TP53, CRLF2))
    (tmp_path / "call.cns").write_text(CALL_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                                                "chrX\t1190891\t2781479\tCRLF2\t0.26805\t2\t679\t8\t7.2\n"
                                                "chrX\t155000000\t155800000\tP2\t-0.1\t2\t500\t8\t7.2\n")
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode != 0
    assert "FAILED: par_check 0" in run.stderr
    assert run.outputs == {}


def test_parfix_span_names_the_gene_and_segment(tmp_path):
    # CRLF2's own segment continues into non-PAR chrX: the documented span limitation, not a failure (the row keeps
    # the male +1: 1.26805 = segment 0.26805 + 1)
    env = parfix_fixtures(tmp_path, parfix_rows=(TP53, CRLF2))
    (tmp_path / "call.cns").write_text(CALL_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                                                "chrX\t1190891\t9000000\tCRLF2\t0.26805\t2\t679\t8\t7.2\n")
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode == 0, run.stderr
    assert "span_genes=CRLF2@chrX:1190891-9000000 par_check=span" in run.stdout
    assert "[parfix] WARNING: 1 call.cns segment(s) span a PAR boundary" in run.stdout
    assert sex_qc_row(run, S)["par_span_segments"] == "1"
    assert set(run.outputs) == set(ALL_OUTPUTS)


def test_parfix_span_row_without_the_male_shift_fails_the_job(tmp_path):
    # R1: an unchanged row in a spanning segment must equal segment + 1; here the segment is -0.9
    env = parfix_fixtures(tmp_path, parfix_rows=(TP53, CRLF2))
    (tmp_path / "call.cns").write_text(CALL_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                                                "chrX\t1190891\t9000000\tCRLF2\t-0.9\t2\t679\t8\t7.2\n")
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode != 0
    assert "FAILED: par_check 0" in run.stderr and "CRLF2@chrX:1190891-1212654" in run.stderr
    assert run.outputs == {}


def test_parfix_span_row_equal_to_the_unshifted_segment_fails_the_job(tmp_path):
    # round 3 (S1): an unchanged CRLF2 row at -0.5 in a spanning segment at -0.5 is not 'corrected'
    row = CRLF2.replace("1.26805", "-0.5")
    env = parfix_fixtures(tmp_path, parfix_rows=(TP53, row))
    (tmp_path / "gm.tsv").write_text(GM_H + TP53 + row)
    (tmp_path / "call.cns").write_text(CALL_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                                                "chrX\t1190891\t9000000\tCRLF2\t-0.5\t2\t679\t8\t7.2\n")
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode != 0
    assert "FAILED: par_check 0" in run.stderr and "CRLF2@chrX:1190891-1212654" in run.stderr
    assert "par_check=span" not in run.stdout
    assert run.outputs == {} and run.uploads == []


def test_parfix_missing_row_not_excused_by_the_same_gene_elsewhere(tmp_path):
    # round 3 (S2): CRLF2 (segment 0.5, so reported after the fix) is missing; a CRLF2 row in PAR2 does not excuse it
    other = "CRLF2\tchrX\t155881000\t155943000\t0.5\t500\t6\t2\t9\t7.2\t8\n"
    env = parfix_fixtures(tmp_path, parfix_rows=(TP53, other))
    (tmp_path / "gm.tsv").write_text(GM_H + TP53 + CRLF2.replace("1.26805", "1.5"))
    (tmp_path / "call.cns").write_text(CALL_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                                                "chrX\t1190891\t2781479\tCRLF2\t0.5\t2\t679\t8\t7.2\n"
                                                "chrX\t155800000\t156000000\tP2\t0.5\t2\t500\t8\t7.2\n")
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode != 0
    assert "FAILED: par_check 0" in run.stderr and "CRLF2@chrX:1190891-1212654 (dropped)" in run.stderr
    assert run.outputs == {} and run.uploads == []


@pytest.mark.parametrize("row, before, after", [
    ("CRLF2\tchrX\t1190891\t1192091\t{}\t0\t3\t6\n", "-9", "-10"),
    ("OVL\tchrX\t2781000\t2781400\t{}\t0\t1.5\t3\n", "-9", "-9.66667"),     # F4: the reviewer's overlapping bins
], ids=["crlf2", "overlapping_bins"])
def test_parfix_without_segments_warns_and_publishes(tmp_path, row, before, after):
    # F2 (round 2), F4 (round 3): a header-only call.cns (an empty segmentation of all-low-coverage bins) makes
    # genemetrics report per-gene bin means; the corrected table (cgp-cnvkit:1.0.0's output) is published, as by the
    # base, with a warning that the PAR rows were not checked (par_check bins: no segment to check them against)
    h = "gene\tchromosome\tstart\tend\tlog2\tdepth\tweight\tprobes\n"
    env = parfix_fixtures(tmp_path)
    (tmp_path / "gm.tsv").write_text(h + row.format(before))
    (tmp_path / "gm_parfix.tsv").write_text(h + row.format(after))
    (tmp_path / "call.cns").write_text("chromosome\tstart\tend\tgene\tlog2\tcn\tdepth\tweight\n")
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode == 0, run.stderr
    assert "par_check=bins" in run.stdout and "[parfix] WARNING: the call.cns has no segments" in run.stdout
    assert "rows were not checked" in run.stdout
    assert (run.workdir / f"{S}.genemetrics.tsv").read_text() == h + row.format(after)
    assert set(run.outputs) == set(ALL_OUTPUTS)


@pytest.mark.parametrize("seg_end, parfix_row, message", [
    (2781479, CRLF2.replace("1.26805", "-1.2"), "CRLF2@chrX:1190891-1212654 (before the fix)"),    # wholly in PAR1
    (9000000, CRLF2.replace("1.26805", "-0.2"), "CRLF2@chrX:1190891-1212654"),                    # spanning
], ids=["wholly", "spanning"])
def test_parfix_original_row_below_the_threshold_fails_the_job(tmp_path, seg_end, parfix_row, message):
    # round 5 (U1): segment -1.2, so genemetrics -t 0.3 cannot report CRLF2 before the fix (-1.2 + 1 = -0.2); an
    # original row at -0.2 is impossible, whether corrected to -1.2 or left unchanged in a spanning segment
    env = parfix_fixtures(tmp_path, parfix_rows=(TP53, parfix_row))
    (tmp_path / "gm.tsv").write_text(GM_H + TP53 + CRLF2.replace("1.26805", "-0.2"))
    (tmp_path / "call.cns").write_text(CALL_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                                                f"chrX\t1190891\t{seg_end}\tCRLF2\t-1.2\t2\t679\t8\t7.2\n")
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode != 0
    assert "FAILED: par_check 0" in run.stderr and f"({message})" in run.stderr
    assert "par_check=1" not in run.stdout and "par_check=span" not in run.stdout
    assert run.outputs == {} and run.uploads == []


@pytest.mark.parametrize("end, after, check", [(2781479, "9.99996", "1"), (9000000, "11", "span")])
def test_parfix_rounding_at_a_decade_completes(tmp_path, end, after, check):
    # round 4 (T3): segment 9.99996; genemetrics prints the valid pre-fix 10.99996 as 11 (%.6g)
    env = parfix_fixtures(tmp_path, parfix_rows=(TP53, CRLF2.replace("1.26805", after)))
    (tmp_path / "gm.tsv").write_text(GM_H + TP53 + CRLF2.replace("1.26805", "11"))
    (tmp_path / "call.cns").write_text(CALL_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                                                f"chrX\t1190891\t{end}\tCRLF2\t9.99996\t2\t679\t8\t7.2\n")
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode == 0, run.stderr
    assert f"par_check={check}" in run.stdout and set(run.outputs) == set(ALL_OUTPUTS)


@pytest.mark.parametrize("parfix_gm, message", [
    ("", "empty file, no header"),                                   # a producer that wrote nothing
    (GM_H + TP53[:20], "line 3 is not a 11-field genemetrics row"),  # truncated mid-row
    ("gene\tchromosome\n", "not a genemetrics header"),
])
def test_parfix_malformed_rerun_output_fails_the_job(tmp_path, parfix_gm, message):
    # round 3 (S3): a post-step genemetrics that fails the structural checks (empty, a header not starting gene
    # chromosome start end log2, a row without one field per column, integer coordinates and a numeric log2) is never
    # published; a truncation leaving well-formed fields is not caught by these checks
    env = parfix_fixtures(tmp_path)
    (tmp_path / "gm_parfix.tsv").write_text(parfix_gm.replace(TP53[:20], TP53 + TP53[:20]))
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode != 0
    assert "cnvkit_parfix.py compare: FAILED" in run.stderr and message in run.stderr
    assert run.outputs == {} and run.uploads == []


def test_parfix_empty_original_and_rerun_fails_the_job(tmp_path):
    # round 3 (S3): two empty genemetrics files used to compare equal (par_check 1)
    env = parfix_fixtures(tmp_path)
    (tmp_path / "gm.tsv").write_text("")
    (tmp_path / "gm_parfix.tsv").write_text("")
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode != 0
    assert "cnvkit_parfix.py compare: FAILED" in run.stderr and "empty file, no header" in run.stderr
    assert run.outputs == {} and run.uploads == []


@pytest.mark.parametrize("case", ["male", "no_parx"])
def test_call_cns_with_overlapping_segments_completes(tmp_path, case):
    # Z1: CNVkit writes overlapping segments for overlapping targets. CRLF2 lies in two overlapping wholly-PAR
    # segments, both at 0.5, and both explain the corrected row: the base completes with par_check 1, and so must
    # the app (the outputs are the base's)
    env = parfix_fixtures(tmp_path, parfix_rows=(TP53, CRLF2.replace("1.26805", "0.5")))
    (tmp_path / "gm.tsv").write_text(GM_H + TP53 + CRLF2.replace("1.26805", "1.5"))
    (tmp_path / "call.cns").write_text(CALL_H + "chr17\t7000000\t8000000\tTP53\t-0.62\t1\t400\t10\t5\n"
                                       "chrX\t1100000\t1300000\tCRLF2\t0.5\t2\t679\t8\t7.2\n"
                                       "chrX\t1190891\t2781479\tCRLF2\t0.5\t2\t679\t8\t7.2\n")
    run = run_app(tmp_path, PARX if case == "male" else {}, stub_env=env)
    assert run.returncode == 0, run.stderr
    assert set(run.outputs) == set(ALL_OUTPUTS)
    if case == "male":
        assert "par_check=1" in run.stdout
        assert (run.workdir / f"{S}.genemetrics.tsv").read_text() == GM_H + TP53 + CRLF2.replace("1.26805", "0.5")


def test_parfix_header_only_tables_complete(tmp_path):
    # CNVkit's legitimate empty table (no gene above the threshold) before and after the fix
    env = parfix_fixtures(tmp_path)
    empty = "gene\tchromosome\tstart\tend\tlog2\n"
    (tmp_path / "gm.tsv").write_text(empty)
    (tmp_path / "gm_parfix.tsv").write_text(empty)
    run = run_app(tmp_path, PARX, stub_env=env)
    assert run.returncode == 0, run.stderr
    assert "par_check=1" in run.stdout and (run.workdir / f"{S}.genemetrics.tsv").read_text() == empty
    assert set(run.outputs) == set(ALL_OUTPUTS)


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


@pytest.mark.parametrize("parx", [{}, PARX])
@pytest.mark.parametrize("stub_env, message", [
    ({"STUB_RAW_Y": "garbage"}, "Y_logratio 'garbage' is not a finite number"),
    ({"STUB_RAW_Y": ""}, "Y_logratio '' is not a finite number"),
    ({"STUB_RAW_SEX_ROW": "{sample}\tMale\t-0.971"}, "the sample row has 3 fields"),     # short row: no Y value
    ({"STUB_RAW_SEX": "unknown"}, "sex 'unknown' is not one of Male, Female, Unknown"),
    ({"STUB_RAW_Y": "NA"}, "Y_logratio 'NA' is not a finite number"),                   # NA only with Unknown
    ({"STUB_CNR_SEX_ROW": "{sample}\tMale\t-0.945"}, "the sample row has 3 fields"),
    ({"STUB_CNR_SEX_ROW": "{sample}\tMale\tgarbage\tnan"}, "X_logratio 'garbage' is not a finite number"),
    ({"STUB_RAW_Y": "-inf"}, "Y_logratio '-inf' is not a finite number"),               # R3: not read as no chrY
    ({"STUB_RAW_Y": "-1e999"}, "Y_logratio '-1e999' is not a finite number"),
    ({"STUB_CNR_SEX_ROW": "{sample}\tMale\t-0.945\tinf"}, "Y_logratio 'inf' is not a finite number"),
])
def test_malformed_sex_table_fails_the_job(tmp_path, parx, stub_env, message):
    # A2: malformed CNVkit sex output fails the job instead of reading as no chrY / female
    run = run_app(tmp_path, parx, stub_env=stub_env)
    assert run.returncode != 0
    assert message in run.stderr
    assert run.outputs == {} and run.uploads == []


@pytest.mark.parametrize("parx", [{}, PARX])
@pytest.mark.parametrize("stub_env", [
    {"STUB_RAW_Y": "-5.01", "STUB_RAW_SEX_END": ""},         # -5.01e-05 cut short: would read as no chrY (Female)
    {"STUB_RAW_SEX_END": ""},
    {"STUB_CNR_SEX_END": ""},
    {"STUB_CNR_SEX_ROW": "{sample}\tMale\t-0.9", "STUB_CNR_SEX_END": ""},    # -0.945 cut short
])
def test_truncated_sex_table_fails_the_job(tmp_path, parx, stub_env):
    # round 6 (V1): CNVkit ends every record with a newline; a table whose last line is unterminated is truncated,
    # and fails the job (the .cnr table in the PAR decision with grch38, else in the sex QC), not the QC alone
    run = run_app(tmp_path, parx, stub_env=stub_env)
    assert run.returncode != 0
    assert "the last line is not newline-terminated (a truncated table)" in run.stderr
    assert run.outputs == {} and run.uploads == []


@pytest.mark.parametrize("cnr_sex", ["garbage", "male", "unknown"])
def test_invalid_inferred_sex_fails_the_parfix_decision(tmp_path, cnr_sex):
    # A2: an invalid .cnr sex call must not read as "not Male" and skip the PAR post-step
    from apprun import default_inputs
    inputs, names = default_inputs(sample_id=U_SAMPLE)
    inputs.update(PARX)
    run = run_app(tmp_path, inputs, names, stub_env={"STUB_CNR_SEX": cnr_sex})
    assert run.returncode != 0
    assert "cnvkit_parfix.py decide: FAILED" in run.stderr and f"sex {cnr_sex!r} is not one of Male, Female, Unknown" in run.stderr
    assert "[parfix] not applied" not in run.stdout
    assert run.outputs == {}


@pytest.mark.parametrize("parx", [{}, PARX])
def test_cnvkit_unknown_sex_completes_as_before(tmp_path, parx):
    # CNVkit writes `Unknown NA NA` when the sample has no usable chrX bins (no chrX targets, or only PAR ones with
    # --diploid-parx-genome); the base accepted it, so the job completes as before: no post-step (genemetrics'
    # own guess is no shift), chrY absent (NA reads as NaN)
    from apprun import default_inputs
    inputs, names = default_inputs(sample_id=U_SAMPLE)
    inputs.update(parx)
    unknown = {"STUB_RAW_SEX_ROW": "{sample}\tUnknown\tNA\tNA", "STUB_CNR_SEX_ROW": "{sample}\tUnknown\tNA\tNA"}
    run = run_app(tmp_path, inputs, names, stub_env=unknown)
    assert run.returncode == 0, run.stderr
    note = "no -x; cnvkit.py sex on the .cnr: Unknown" if parx else "no --diploid-parx-genome"
    assert f"[parfix] not applied: {note}" in run.stdout
    row = sex_qc_row(run, U_SAMPLE)
    assert {k: row[k] for k in ("applied_sex", "raw_Y_logratio", "chrY_present", "chrY_sex", "raw_cnvkit_sex",
                                "cnr_X_logratio", "cnr_cnvkit_sex", "sex_mismatch", "par_fix_applied")} == {
        "applied_sex": "cnvkit_inferred(Unknown)", "raw_Y_logratio": "NA", "chrY_present": "0", "chrY_sex": "Female",
        "raw_cnvkit_sex": "Unknown", "cnr_X_logratio": "NA", "cnr_cnvkit_sex": "Unknown", "sex_mismatch": "NA",
        "par_fix_applied": "0"}


@pytest.mark.parametrize("sample", ["patient__01", "sample.cnr", "sample.cnn", "x.targetcoverage.cnn"])
def test_sample_id_with_separator_or_suffix_completes(tmp_path, sample):
    # A3: accepted sample IDs containing `__` or ending in .cnr/.cnn are matched whole in the sex tables
    from apprun import default_inputs
    inputs, names = default_inputs(sample_id=sample)
    inputs.update(PARX)
    run = run_app(tmp_path, inputs, names)
    assert run.returncode == 0, run.stderr
    assert sex_qc_row(run, sample)["sample"] == sample
    assert {k: v["value"] for k, v in run.outputs.items()} == \
        {k: f"file-{v.replace(S, sample)}" for k, v in ALL_OUTPUTS.items()}


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


@pytest.mark.parametrize("failing, missing", [("scatter", "scatter_png"), ("diagram", "diagram_pdf")])
def test_partial_plot_is_not_published(tmp_path, failing, missing):
    # A4: a plot command that writes part of its file and then fails: the partial file is removed, not uploaded
    run = run_app(tmp_path, stub_env={"STUB_FAIL_PARTIAL": failing})
    assert run.returncode == 0, run.stderr
    assert f"[{failing}] WARNING" in run.stdout
    assert not (run.workdir / ALL_OUTPUTS[missing]).exists()
    assert ALL_OUTPUTS[missing] not in run.uploads
    assert set(run.outputs) == set(ALL_OUTPUTS) - {missing}


def test_both_partial_plots_are_not_published(tmp_path):
    run = run_app(tmp_path, stub_env={"STUB_FAIL_PARTIAL": "scatter,diagram"})
    assert run.returncode == 0, run.stderr
    assert set(run.outputs) == set(ALL_OUTPUTS) - {"scatter_png", "diagram_pdf"}
    assert not any(u.endswith((".png", ".pdf")) for u in run.uploads)


@pytest.mark.parametrize("failing", ["coverage", "fix", "segment", "sex", "call", "genemetrics"])
def test_analysis_failure_is_fatal(tmp_path, failing):
    run = run_app(tmp_path, fail=(failing,))
    assert run.returncode != 0
    assert run.outputs == {}
