#!/bin/bash
# cgp-cnvkit-batch/src/code.sh
# Per-sample CNVkit analysis: coverage → fix → segment → call → plot → genemetrics
# Accepts optional PURPLE purity/ploidy for integer CN calling.
# Only calls integers when purity is supplied AND purity > 0.40 (PURPLE reliable zone).
set -eo pipefail

main() {
    echo "=== CGP CNVkit batch: ${sample_id} ==="
    echo "Instance: $(hostname); CPUs: $(nproc)"

    # ── Install cnvkit + R DNAcopy for CBS segmentation ───────────────────────
    echo "[setup] Installing cnvkit..."
    python3 -m venv /tmp/cnvkit-env
    /tmp/cnvkit-env/bin/pip install cnvkit --quiet 2>&1 | tail -3
    export PATH="/tmp/cnvkit-env/bin:$PATH"
    cnvkit.py version

    echo "[setup] Installing R DNAcopy for CBS segmentation..."
    Rscript -e "
        options(warn=-1)
        if (!requireNamespace('BiocManager', quietly=TRUE))
            install.packages('BiocManager', repos='https://cloud.r-project.org', quiet=TRUE)
        if (!requireNamespace('DNAcopy', quietly=TRUE))
            BiocManager::install('DNAcopy', ask=FALSE, quiet=TRUE)
        cat('DNAcopy version:', as.character(packageVersion('DNAcopy')), '\n')
    "

    # ── Stage inputs ──────────────────────────────────────────────────────────
    echo "[inputs] Downloading files..."
    dx download "${tumour_bam}"    -o tumour.bam
    dx download "${tumour_bai}"    -o tumour.bam.bai
    dx download "${baits}"         -o panel.bed
    dx download "${cn_reference}"  -o reference.cnn

    # ── Coverage (amplicon mode — panel BED directly, no autobin) ─────────────
    echo "[coverage] Computing on-target coverage..."
    cnvkit.py coverage tumour.bam panel.bed \
        --processes "$(nproc)" \
        --output "${sample_id}.targetcoverage.cnn"

    # ── Fix: apply PoN reference ──────────────────────────────────────────────
    echo "[fix] Applying PoN reference..."
    # cnvkit.py fix requires 3 positional args even in amplicon mode:
    # target antitarget reference
    # Pass an empty antitarget .cnn (header only) as placeholder.
    # Note: --drop-low-coverage is a batch-mode flag, not valid for fix.
    printf "chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n" \
        > empty.antitargetcoverage.cnn

    cnvkit.py fix \
        "${sample_id}.targetcoverage.cnn" \
        empty.antitargetcoverage.cnn \
        reference.cnn \
        --no-edge \
        -o "${sample_id}.cnr"

    # ── Segment: CBS ──────────────────────────────────────────────────────────
    echo "[segment] Running CBS segmentation..."
    cnvkit.py segment \
        "${sample_id}.cnr" \
        --method cbs \
        --drop-outliers 10 \
        --processes "$(nproc)" \
        -o "${sample_id}.cns"

    N_SEGS=$(grep -vc "^chromosome" "${sample_id}.cns" || true)
    echo "[segment] Segments: ${N_SEGS}"

    # ── Call: integer CN ──────────────────────────────────────────────────────
    echo "[call] Calling integer CNs..."
    CALL_ARGS=""
    if [ -n "${purity:-}" ]; then
        PURITY_OK=$(python3 -c "print('yes' if float('${purity}') > 0.40 else 'no')")
        if [ "$PURITY_OK" = "yes" ]; then
            CALL_ARGS="--purity ${purity} --ploidy ${ploidy:-2}"
            echo "[call] Using PURPLE purity=${purity}, ploidy=${ploidy:-2}"
        else
            echo "[call] purity=${purity} <= 0.40 — log2R thresholds only"
        fi
    else
        echo "[call] No purity supplied — log2R thresholds only"
    fi

    cnvkit.py call \
        "${sample_id}.cns" \
        ${CALL_ARGS} \
        -o "${sample_id}.call.cns"

    # ── Scatter plot ──────────────────────────────────────────────────────────
    # y-axis clamped to ±3: amplifications clip cleanly at the top while losses
    # at -1 to -1.5 remain clearly legible. Amplified segment lines still show
    # at the scale edge, signalling the event continues off-screen.
    # 14×4 inch figure matches PURPLE plot aspect ratio for side-by-side review.
    # Note: cnvkit.py scatter has no --sample-sex flag; chrX will appear at
    # log2R ≈ -1 for male samples due to mixed-gender PoN (cosmetic artefact).
    echo "[scatter] Generating scatter plot..."
    cnvkit.py scatter \
        "${sample_id}.cnr" \
        -s "${sample_id}.cns" \
        --trend \
        --y-min -3 --y-max 3 \
        --fig-size 14 4 \
        --title "${sample_id}" \
        -o "${sample_id}.scatter.png" 2>/dev/null || \
    echo "[scatter] WARNING: scatter plot failed (non-fatal)"

    # ── Diagram plot (chromosome ideogram with gene labels) ───────────────────
    # --sample-sex corrects chrX/Y baseline for the given sex.
    # Omit if sample_sex not supplied — CNVkit will guess from X/Y coverage.
    echo "[diagram] Generating chromosome diagram..."
    SEX_ARG=""
    [ -n "${sample_sex:-}" ] && SEX_ARG="--sample-sex ${sample_sex}"
    cnvkit.py diagram \
        "${sample_id}.cnr" \
        -s "${sample_id}.call.cns" \
        --threshold 0.5 \
        --min-probes 3 \
        ${SEX_ARG} \
        --title "${sample_id}" \
        -o "${sample_id}.diagram.pdf" 2>/dev/null || \
    echo "[diagram] WARNING: diagram failed (non-fatal)"

    # ── Gene metrics ──────────────────────────────────────────────────────────
    echo "[genemetrics] Computing per-gene metrics..."
    cnvkit.py genemetrics \
        "${sample_id}.cnr" \
        -s "${sample_id}.call.cns" \
        -t 0.3 \
        --min-probes 3 \
        -o "${sample_id}.genemetrics.csv"

    N_GENES=$(wc -l < "${sample_id}.genemetrics.csv")
    echo "[genemetrics] Genes with CN change (threshold 0.3): $((N_GENES - 1))"

    # ── Upload outputs ────────────────────────────────────────────────────────
    copy_ratios=$(dx   upload "${sample_id}.cnr"             --brief)
    segments=$(dx      upload "${sample_id}.cns"             --brief)
    call_segments=$(dx upload "${sample_id}.call.cns"        --brief)
    genemetrics=$(dx   upload "${sample_id}.genemetrics.csv" --brief)

    dx-jobutil-add-output copy_ratios   "${copy_ratios}"   --class=file
    dx-jobutil-add-output segments      "${segments}"      --class=file
    dx-jobutil-add-output call_segments "${call_segments}" --class=file
    dx-jobutil-add-output genemetrics   "${genemetrics}"   --class=file

    if [ -f "${sample_id}.scatter.png" ]; then
        scatter_png=$(dx upload "${sample_id}.scatter.png" --brief)
        dx-jobutil-add-output scatter_png "${scatter_png}" --class=file
    fi

    if [ -f "${sample_id}.diagram.pdf" ]; then
        diagram_pdf=$(dx upload "${sample_id}.diagram.pdf" --brief)
        dx-jobutil-add-output diagram_pdf "${diagram_pdf}" --class=file
    fi

    echo "=== Done: ${sample_id} — ${N_SEGS} segments, $((N_GENES-1)) genes with CN change ==="
}
