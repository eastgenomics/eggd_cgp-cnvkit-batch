#!/bin/bash
# cgp_cnvkit_batch/pipeline.sh — the per-sample CNVkit analysis of eggd_cgp-cnvkit-batch.
#
# Sourced by src/code.sh on the worker (installed from resources/ to /usr/local/lib/cgp_cnvkit_batch/) and by
# the tests. No DNAnexus calls: code.sh downloads the inputs and uploads the outputs.
#
# Contract
#   caller defines   run_cnvkit ARGS...   runs `cnvkit.py ARGS...` with the current directory as working directory
#   inputs (files in the current directory)
#     tumour.bam, tumour.bam.bai, panel.bed, reference.cnn
#   inputs (environment)
#     SAMPLE_ID           output file stem (required)
#     NPROC               CPUs for coverage/segment (default: nproc)
#     START_FROM          bam (default) | coverage: start from an existing <SAMPLE_ID>.targetcoverage.cnn (no BAM or
#                         BED needed; the local harness uses it to start from the validated coverage)
#     DIPLOID_PARX_GENOME grch38 | '' (DI-4097): --diploid-parx-genome on fix, segment, call, genemetrics, diagram
#     PURITY, PLOIDY      deprecated (DI-4121, D4): ignored, with a warning when supplied
#     BAM_NAME            BAM file name: LIMS-letter fallback when SAMPLE_ID has no 5th field (DI-4098)
#     SEX_OVERRIDE        male | female | M | F | '' (DI-4099, D9): wins over the LIMS letter
#     SAMPLE_SEX_DEPRECATED  the old sample_sex input (wired from PURPLE): ignored, with a warning (D2)
#     DROP_LOW_COVERAGE   true (default) | false
#     DROPOUT_MASK        true (default) | false (DI-4139, D7): segment on the dropout-masked copy of the .cnr
#   outputs (files in the current directory)
#     <SAMPLE_ID>.targetcoverage.cnn .cnr .cns .call.cns .genemetrics.tsv .genemetrics.noparfix.tsv
#     .dropout_qc.tsv .sex_raw.tsv .sex_cnr.tsv .sex_qc.tsv .sex_check_mqc.tsv, and .scatter.png / .diagram.pdf
#     when the (non-fatal) plots succeed
#   outputs (shell variables)  N_SEGS, N_GENES, SEX_ARG, SEX_SOURCE, LIMS_SEX, DECLARED_SEX, DROPOUT_PROBES,
#                              PAR_FIX_APPLIED, PAR_SPAN_SEGMENTS

dedup_bed() {
    # Deduplicate BED intervals: some panel BEDs have identical coordinates from multiple transcript
    # annotations, and CNVkit 0.9.13+ raises on duplicates.
    awk '!seen[$1 FS $2 FS $3]++' "$1" > "$1.dedup" && mv "$1.dedup" "$1"
}

run_pipeline() {
    : "${SAMPLE_ID:?SAMPLE_ID is required}"
    local s="${SAMPLE_ID}"
    local nproc="${NPROC:-$(nproc)}"

    # ── Parameters (checked before any analysis step) ─────────────────────────
    local -a px=()
    case "${DIPLOID_PARX_GENOME:-}" in
        grch38)
            px=(--diploid-parx-genome grch38)
            echo "[parx] --diploid-parx-genome grch38 on fix, segment, call, genemetrics and diagram" ;;
        "")
            echo "[parx] WARNING: diploid_parx_genome is empty: chrX PAR1/PAR2 are treated as X-specific, so"
            echo "[parx] WARNING: CNVkit's sex inference calls every sample female. Atlas requires grch38." ;;
        *)
            echo "ERROR: diploid_parx_genome '${DIPLOID_PARX_GENOME}' is not supported (grch38, or empty)"
            return 1 ;;
    esac

    # Sex (DI-4098, DI-4099, DI-4100): override, else the LIMS letter (field 5 of the sample ID, else of the BAM
    # name), else no -x (CNVkit infers). The applied sex goes to call, genemetrics and diagram as -x.
    local resolved
    resolved=$(cnvkit_sex_qc.py resolve --sample-id "${s}" --bam-name "${BAM_NAME:-}" \
        --override "${SEX_OVERRIDE:-}") || return 1
    SEX_ARG=$(cut -f1 <<< "${resolved}")
    SEX_SOURCE=$(cut -f2 <<< "${resolved}")
    LIMS_SEX=$(cut -f3 <<< "${resolved}")
    DECLARED_SEX=$(cut -f4 <<< "${resolved}")
    local -a x=()
    if [ -n "${SEX_ARG}" ]; then
        x=(-x "${SEX_ARG}")
        echo "[sex] applied: -x ${SEX_ARG} (source: ${SEX_SOURCE}; LIMS ${LIMS_SEX})"
    else
        echo "[sex] applied: no -x, CNVkit infers the sex (source: ${SEX_SOURCE}; LIMS ${LIMS_SEX})"
    fi
    if [ -n "${SAMPLE_SEX_DEPRECATED:-}" ]; then
        echo "[deprecated] sample_sex supplied but ignored (D2): it is wired from PURPLE's inferred sex;"
        echo "[deprecated] use sex_override to set the sex explicitly."
    fi

    if [ -n "${PURITY:-}${PLOIDY:-}" ]; then
        echo "[deprecated] purity/ploidy supplied but ignored (DI-4121): CNVkit calls are made without"
        echo "[deprecated] --purity/--ploidy, as validated; PURPLE is the source of absolute copy number."
    fi

    # ── Coverage (amplicon mode — panel BED directly, no autobin) ─────────────
    if [ "${START_FROM:-bam}" = "coverage" ]; then
        if [ ! -s "${s}.targetcoverage.cnn" ]; then
            echo "ERROR: START_FROM=coverage but ${s}.targetcoverage.cnn is missing"
            return 1
        fi
        echo "[coverage] START_FROM=coverage: using the existing ${s}.targetcoverage.cnn"
    else
        dedup_bed panel.bed
        echo "[coverage] Computing on-target coverage..."
        run_cnvkit coverage tumour.bam panel.bed \
            --processes "${nproc}" \
            --output "${s}.targetcoverage.cnn"
    fi

    # ── Fix: apply PoN reference ──────────────────────────────────────────────
    echo "[fix] Applying PoN reference..."
    # cnvkit.py fix requires 3 positional args even in amplicon mode: target antitarget reference.
    # Pass an empty antitarget .cnn (header only) as placeholder. Edge correction is on (DI-4096): the PoN is
    # built with CNVkit's default edge correction, so reference and fix match.
    printf "chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n" > empty.antitargetcoverage.cnn
    run_cnvkit fix \
        "${s}.targetcoverage.cnn" \
        empty.antitargetcoverage.cnn \
        reference.cnn \
        "${px[@]}" \
        -o "${s}.cnr"

    # ── Dropout mask (DI-4139) ────────────────────────────────────────────────
    # Isolated autosomal probe dropouts (log2 < -3, run of <= 2, < 10 raw read-equivalents) get weight 0 in a
    # copy of the .cnr that only segment reads; every other step reads the unmasked <s>.cnr.
    local seg_in="${s}.cnr"
    if [ "${DROPOUT_MASK:-true}" = "true" ]; then
        dropout_mask.py mask M "${s}.cnr" "${s}.targetcoverage.cnn" "${s}" - \
            "${s}.seg_input.cnr" "${s}.dropout_qc.tsv" > /dev/null
        seg_in="${s}.seg_input.cnr"
        DROPOUT_PROBES=$(( $(wc -l < "${s}.dropout_qc.tsv") - 1 ))
        echo "[dropout] masked probes: ${DROPOUT_PROBES} (listed in ${s}.dropout_qc.tsv)"
    else
        dropout_mask.py qc-header "${s}.dropout_qc.tsv"
        DROPOUT_PROBES=0
        echo "[dropout] dropout_mask=false: segment reads the unmasked .cnr (validated configuration without the mask)"
    fi

    # ── Segment: CBS ──────────────────────────────────────────────────────────
    echo "[segment] Running CBS segmentation on ${seg_in}..."
    local -a seg_args=()
    [ "${DROP_LOW_COVERAGE:-true}" = "true" ] && seg_args+=(--drop-low-coverage)
    run_cnvkit segment \
        "${seg_in}" \
        --method cbs \
        --drop-outliers 10 \
        "${seg_args[@]}" \
        "${px[@]}" \
        --processes "${nproc}" \
        -o "${s}.cns"
    N_SEGS=$(grep -vc "^chromosome" "${s}.cns" || true)
    echo "[segment] Segments: ${N_SEGS}"

    # ── Sex checks (DI-4102) ──────────────────────────────────────────────────
    # Raw coverage keeps the non-PAR chrY bins (chrY presence: the QC flag); the .cnr has chrX only (the
    # PoN masks chrY), so its call is the guess genemetrics makes without -x (the PAR post-step decision).
    echo "[sex] cnvkit.py sex on the raw coverage and on the .cnr..."
    run_cnvkit sex "${s}.targetcoverage.cnn" "${px[@]}" -o "${s}.sex_raw.tsv"
    run_cnvkit sex "${s}.cnr" "${px[@]}" -o "${s}.sex_cnr.tsv"

    # ── Call: log2 ratio thresholds (no --purity/--ploidy; DI-4121) ──────────
    echo "[call] Calling copy number from log2 ratio thresholds..."
    run_cnvkit call \
        "${s}.cns" \
        "${px[@]}" \
        "${x[@]}" \
        -o "${s}.call.cns"

    # ── Scatter plot (non-fatal) ──────────────────────────────────────────────
    # y-axis clamped to ±3: amplifications clip cleanly at the top while losses at -1 to -1.5 remain legible.
    # 14×4 inch figure matches the PURPLE plot aspect ratio for side-by-side review.
    # cnvkit.py scatter has no --sample-sex flag: chrX appears at log2R ≈ -1 for male samples.
    echo "[scatter] Generating scatter plot..."
    run_cnvkit scatter \
        "${s}.cnr" \
        -s "${s}.cns" \
        --trend \
        --y-min -3 --y-max 3 \
        --fig-size 14 4 \
        --title "${s}" \
        -o "${s}.scatter.png" 2>/tmp/scatter_err.txt || \
    echo "[scatter] WARNING: scatter plot failed (non-fatal): $(tail -1 /tmp/scatter_err.txt)"

    # ── Diagram plot (chromosome ideogram with gene labels; non-fatal) ────────
    # -x corrects the chrX/Y baseline, as in call and genemetrics, so the plot matches the calls.
    echo "[diagram] Generating chromosome diagram..."
    run_cnvkit diagram \
        "${s}.cnr" \
        -s "${s}.call.cns" \
        --threshold 0.5 \
        --min-probes 3 \
        "${px[@]}" \
        "${x[@]}" \
        --title "${s}" \
        -o "${s}.diagram.pdf" 2>/tmp/diagram_err.txt || \
    echo "[diagram] WARNING: diagram failed (non-fatal): $(tail -1 /tmp/diagram_err.txt)"

    # ── Gene metrics ──────────────────────────────────────────────────────────
    echo "[genemetrics] Computing per-gene metrics..."
    local -a gm_args=(-t 0.3 --min-probes 3 "${px[@]}" "${x[@]}")
    run_cnvkit genemetrics "${s}.cnr" -s "${s}.call.cns" "${gm_args[@]}" -o "${s}.genemetrics.tsv"

    # ── PAR post-step (DI-4101) ───────────────────────────────────────────────
    # For a sample genemetrics treats as male, CNVkit (fc65941) shifts every chrX bin +1.0, the PARs included.
    # Re-run genemetrics, with the same arguments, on copies with log2 - 1 in the bins/segments wholly inside
    # PAR1/PAR2; every non-PAR row must be byte-identical (the job fails otherwise). The uncorrected read-out is
    # kept as <s>.genemetrics.noparfix.tsv (an identical copy when the post-step does not run; D3).
    cp "${s}.genemetrics.tsv" "${s}.genemetrics.noparfix.tsv"
    local decision
    decision=$(cnvkit_parfix.py decide --sex-arg "${SEX_ARG}" --cnr-sex "${s}.sex_cnr.tsv" \
        --parx "${DIPLOID_PARX_GENOME:-}")
    PAR_FIX_APPLIED=$(cut -f1 <<< "${decision}")
    PAR_SPAN_SEGMENTS=0
    if [ "${PAR_FIX_APPLIED}" = "1" ]; then
        local shifted
        shifted=$(cnvkit_parfix.py shift "${s}.cnr" "${s}.call.cns" parfix)
        echo "[parfix] post-step applied ($(cut -f2 <<< "${decision}")): ${shifted}"
        PAR_SPAN_SEGMENTS=$(sed -n 's/.*span_segs=\([0-9]*\).*/\1/p' <<< "${shifted}")
        if [ "${PAR_SPAN_SEGMENTS}" -gt 0 ]; then
            echo "[parfix] WARNING: ${PAR_SPAN_SEGMENTS} call.cns segment(s) span a PAR boundary: PAR genes in them"
            echo "[parfix] WARNING: keep the +1 shift (known limitation): review male CRLF2 calls."
        fi
        (cd parfix && run_cnvkit genemetrics "${s}.cnr" -s "${s}.call.cns" "${gm_args[@]}" -o "${s}.genemetrics.tsv")
        # assigned first, so a failed check (non-PAR rows differ) stops the job under set -e
        local checked
        checked=$(cnvkit_parfix.py compare "${s}.genemetrics.noparfix.tsv" "parfix/${s}.genemetrics.tsv" \
            "${s}.call.cns")
        echo "[parfix] ${checked}"
        cp -f "parfix/${s}.genemetrics.tsv" "${s}.genemetrics.tsv"
    else
        echo "[parfix] not applied: $(cut -f2 <<< "${decision}")"
    fi

    N_GENES=$(( $(wc -l < "${s}.genemetrics.tsv") - 1 ))
    echo "[genemetrics] Genes with CN change (threshold 0.3): ${N_GENES}"

    # ── Sex QC (DI-4103) ──────────────────────────────────────────────────────
    # sex_mismatch: declared sex (override, else LIMS) vs chrY presence in the raw coverage (Y log-ratio > -5).
    local override_sex=U
    [ "${SEX_SOURCE}" = "override" ] && override_sex="${DECLARED_SEX}"
    cnvkit_sex_qc.py qc --sample-id "${s}" --raw-sex "${s}.sex_raw.tsv" --cnr-sex "${s}.sex_cnr.tsv" \
        --lims-sex "${LIMS_SEX}" --override-sex "${override_sex}" --applied-source "${SEX_SOURCE}" \
        --par-fix-applied "${PAR_FIX_APPLIED}" --par-span-segments "${PAR_SPAN_SEGMENTS}" \
        --dropout-probes "${DROPOUT_PROBES}" --out "${s}"
}
