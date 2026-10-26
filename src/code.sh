#!/bin/bash
# eggd_cgp-cnvkit-batch/src/code.sh
# DNAnexus entry point: load the CNVkit image, download and check the inputs, run the per-sample analysis
# (resources/usr/local/lib/cgp_cnvkit_batch/pipeline.sh), upload the outputs.
set -eo pipefail

CGP_CNVKIT_LIB="${CGP_CNVKIT_LIB:-/usr/local/lib/cgp_cnvkit_batch}"

# The DNAnexus bash executor calls main(); this file must not call it.
main() {
    echo "=== CGP CNVkit batch: ${sample_id} ==="
    # sample_id is the stem of every file name: letters, digits and . _ - only, not starting with - or .
    if [[ ! "${sample_id}" =~ ^[A-Za-z0-9_][A-Za-z0-9._-]*$ ]]; then
        echo "ERROR: sample_id '${sample_id}' is not a valid file stem (letters, digits, '.', '_', '-';"
        echo "       must not start with '.' or '-')"
        exit 1
    fi
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
    START_FROM=bam \
    DIPLOID_PARX_GENOME="${diploid_parx_genome:-}" \
    BAM_NAME="${BAM_NAME}" \
    SEX_OVERRIDE="${sex_override:-}" \
    SAMPLE_SEX_DEPRECATED="${sample_sex:-}" \
    PURITY="${purity:-}" \
    PLOIDY="${ploidy:-}" \
    DROP_LOW_COVERAGE="${drop_low_coverage:-true}" \
    DROPOUT_MASK="${dropout_mask:-true}" \
        run_pipeline

    # ── Upload outputs ────────────────────────────────────────────────────────
    upload() {   # upload OUTPUT_NAME FILE: upload the file and add it as a job output
        local id
        id=$(dx upload "$2" --brief)     # assigned first, so a failed upload stops the job under set -e
        dx-jobutil-add-output "$1" "${id}" --class=file
    }
    upload copy_ratios          "${sample_id}.cnr"                       # unmasked
    upload segments             "${sample_id}.cns"                       # segmented on the dropout-masked input
    upload call_segments        "${sample_id}.call.cns"
    upload genemetrics          "${sample_id}.genemetrics.tsv"           # after the PAR post-step
    upload genemetrics_noparfix "${sample_id}.genemetrics.noparfix.tsv"
    upload sex_qc               "${sample_id}.sex_qc.tsv"
    upload sex_check_mqc        "${sample_id}.sex_check_mqc.tsv"
    upload dropout_qc           "${sample_id}.dropout_qc.tsv"
    upload target_coverage      "${sample_id}.targetcoverage.cnn"
    if [ -f "${sample_id}.scatter.png" ]; then upload scatter_png "${sample_id}.scatter.png"; fi
    if [ -f "${sample_id}.diagram.pdf" ]; then upload diagram_pdf "${sample_id}.diagram.pdf"; fi

    echo "=== Done: ${sample_id} — ${N_SEGS} segments, ${N_GENES} genes with CN change;" \
         "sex ${SEX_ARG:-inferred} (${SEX_SOURCE}); PAR post-step ${PAR_FIX_APPLIED};" \
         "masked probes ${DROPOUT_PROBES} ==="
}
