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
#     DIPLOID_PARX_GENOME grch38 | '' (DI-4097): --diploid-parx-genome on fix, segment, call, genemetrics, diagram
#     PURITY, PLOIDY      optional tumour purity / ploidy
#     SAMPLE_SEX          optional chromosomal sex for the diagram
#     DROP_LOW_COVERAGE   true (default) | false
#   outputs (files in the current directory)
#     <SAMPLE_ID>.targetcoverage.cnn .cnr .cns .call.cns .genemetrics.tsv, and .scatter.png / .diagram.pdf when
#     the (non-fatal) plots succeed
#   outputs (shell variables)  N_SEGS, N_GENES

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

    dedup_bed panel.bed

    # ── Coverage (amplicon mode — panel BED directly, no autobin) ─────────────
    echo "[coverage] Computing on-target coverage..."
    run_cnvkit coverage tumour.bam panel.bed \
        --processes "${nproc}" \
        --output "${s}.targetcoverage.cnn"

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

    # ── Segment: CBS ──────────────────────────────────────────────────────────
    echo "[segment] Running CBS segmentation..."
    local -a seg_args=()
    [ "${DROP_LOW_COVERAGE:-true}" = "true" ] && seg_args+=(--drop-low-coverage)
    run_cnvkit segment \
        "${s}.cnr" \
        --method cbs \
        --drop-outliers 10 \
        "${seg_args[@]}" \
        "${px[@]}" \
        --processes "${nproc}" \
        -o "${s}.cns"
    N_SEGS=$(grep -vc "^chromosome" "${s}.cns" || true)
    echo "[segment] Segments: ${N_SEGS}"

    # ── Call: integer CN ──────────────────────────────────────────────────────
    echo "[call] Calling integer CNs..."
    local -a call_args=()
    if [ -n "${PURITY:-}" ]; then
        if [ "$(python3 -c "print('yes' if float('${PURITY}') > 0.40 else 'no')")" = "yes" ]; then
            call_args=(--purity "${PURITY}" --ploidy "${PLOIDY:-2}")
            echo "[call] Using PURPLE purity=${PURITY}, ploidy=${PLOIDY:-2}"
        else
            echo "[call] purity=${PURITY} <= 0.40 — log2R thresholds only"
        fi
    else
        echo "[call] No purity supplied — log2R thresholds only"
    fi
    run_cnvkit call \
        "${s}.cns" \
        "${call_args[@]}" \
        "${px[@]}" \
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
    # --sample-sex corrects the chrX/Y baseline; omitted when no sex is supplied (CNVkit guesses).
    echo "[diagram] Generating chromosome diagram..."
    local -a sex_args=()
    [ -n "${SAMPLE_SEX:-}" ] && sex_args=(--sample-sex "${SAMPLE_SEX}")
    run_cnvkit diagram \
        "${s}.cnr" \
        -s "${s}.call.cns" \
        --threshold 0.5 \
        --min-probes 3 \
        "${sex_args[@]}" \
        "${px[@]}" \
        --title "${s}" \
        -o "${s}.diagram.pdf" 2>/tmp/diagram_err.txt || \
    echo "[diagram] WARNING: diagram failed (non-fatal): $(tail -1 /tmp/diagram_err.txt)"

    # ── Gene metrics ──────────────────────────────────────────────────────────
    echo "[genemetrics] Computing per-gene metrics..."
    run_cnvkit genemetrics \
        "${s}.cnr" \
        -s "${s}.call.cns" \
        -t 0.3 \
        --min-probes 3 \
        "${px[@]}" \
        -o "${s}.genemetrics.tsv"
    N_GENES=$(( $(wc -l < "${s}.genemetrics.tsv") - 1 ))
    echo "[genemetrics] Genes with CN change (threshold 0.3): ${N_GENES}"
}
