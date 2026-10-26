#!/bin/bash
# eggd_cgp-cnvkit-batch/src/code.sh
# DNAnexus entry point: load the CNVkit image, download and check the inputs, run the per-sample analysis
# (resources/usr/local/lib/cgp_cnvkit_batch/pipeline.sh), upload the outputs.
set -eo pipefail

CGP_CNVKIT_LIB="${CGP_CNVKIT_LIB:-/usr/local/lib/cgp_cnvkit_batch}"

main() {
    echo "=== CGP CNVkit batch: ${sample_id} ==="
    echo "Instance: $(hostname); CPUs: $(nproc)"
    # shellcheck source=SCRIPTDIR/../resources/usr/local/lib/cgp_cnvkit_batch/pipeline.sh
    source "${CGP_CNVKIT_LIB}/pipeline.sh"

    # ── Load CNVkit Docker image ──────────────────────────────────────────────
    # Image stored in DNAnexus; no external internet required.
    CNVKIT_IMAGE_ID="project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8j7Vyj45FG1BbK26JQgQY6q"   # cgp-cnvkit:1.0.0
    CNVKIT_IMAGE_TAG="cgp-cnvkit:1.0.0"
    echo "[setup] Loading CNVkit image..."
    dx download "${CNVKIT_IMAGE_ID}" -o cnvkit-image.tar.gz
    docker load < cnvkit-image.tar.gz
    run_cnvkit() { docker run --rm -v "$(pwd)":/work -w /work "${CNVKIT_IMAGE_TAG}" cnvkit.py "$@"; }
    run_cnvkit version

    # ── Stage inputs ──────────────────────────────────────────────────────────
    echo "[inputs] Downloading files..."
    dx download "${tumour_bam}"    -o tumour.bam
    dx download "${tumour_bai}"    -o tumour.bam.bai
    dx download "${baits}"         -o panel.bed

    # Sanity check: verify BAM and BAI filenames are a matched pair
    BAM_NAME=$(dx describe "${tumour_bam}" --json | python3 -c "import sys,json; print(json.load(sys.stdin)['name'])")
    BAI_NAME=$(dx describe "${tumour_bai}" --json | python3 -c "import sys,json; print(json.load(sys.stdin)['name'])")
    if [[ "${BAI_NAME}" != "${BAM_NAME}.bai" && "${BAI_NAME}" != "${BAM_NAME%.bam}.bai" ]]; then
        echo "ERROR: BAI filename '${BAI_NAME}' does not match BAM filename '${BAM_NAME}'"
        echo "       Expected '${BAM_NAME}.bai' or '${BAM_NAME%.bam}.bai'"
        exit 1
    fi
    echo "[inputs] BAM/BAI pair validated: ${BAM_NAME} / ${BAI_NAME}"
    dx download "${cn_reference}"  -o reference.cnn

    # ── Analysis ──────────────────────────────────────────────────────────────
    SAMPLE_ID="${sample_id}" \
    DIPLOID_PARX_GENOME="${diploid_parx_genome:-}" \
    PURITY="${purity:-}" \
    PLOIDY="${ploidy:-}" \
    SAMPLE_SEX="${sample_sex:-}" \
    DROP_LOW_COVERAGE="${drop_low_coverage:-true}" \
        run_pipeline

    # ── Upload outputs ────────────────────────────────────────────────────────
    copy_ratios=$(dx   upload "${sample_id}.cnr"             --brief)
    segments=$(dx      upload "${sample_id}.cns"             --brief)
    call_segments=$(dx upload "${sample_id}.call.cns"        --brief)
    genemetrics=$(dx   upload "${sample_id}.genemetrics.tsv" --brief)

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

    echo "=== Done: ${sample_id} — ${N_SEGS} segments, ${N_GENES} genes with CN change ==="
}
