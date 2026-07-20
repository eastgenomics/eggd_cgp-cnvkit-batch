<!-- dx-header -->
# eggd_cgp-cnvkit-batch (DNAnexus Platform App)

Per-sample CNVkit analysis for tumour-only CGP panel sequencing.
Runs as **Step 3** of the three-app CGP CNVkit pipeline:

```
eggd_cgp-cnvkit-coverage (×N, parallel) → eggd_cgp-cnvkit-pon (×1) → eggd_cgp-cnvkit-batch (×N, parallel)
```

## What does this app do?

Runs the full CNVkit per-sample analysis pipeline against a pre-built PoN reference:

1. **Coverage** — `cnvkit.py coverage` in amplicon mode (panel BED, no antitarget bins)
2. **Fix** — `cnvkit.py fix` applies PoN normalisation to produce per-bin log₂ copy ratios
3. **Segment** — `cnvkit.py segment` runs CBS (circular binary segmentation) with outlier protection (`--drop-outliers 10`)
4. **Call** — `cnvkit.py call` assigns integer CN calls; uses PURPLE purity/ploidy when
   purity > 0.40, otherwise falls back to log₂R thresholds
5. **Scatter plot** — genome-wide log₂R scatter with y-axis clamped to ±3 and `--fig-size 14 4`
   for readability alongside PURPLE backbone plots
6. **Diagram** — chromosome ideogram with gene labels (sex-corrected when `sample_sex` supplied)
7. **Gene metrics** — per-gene log₂R summary filtered to |log₂R| > 0.3, ≥3 probes

R and the `DNAcopy` Bioconductor package are installed at job start for CBS segmentation.

## What are the typical use cases for this app?

- Gene-panel focal copy-number calling at exon resolution (~2–10 kb) to complement
  a backbone-based tool such as PURPLE (~50 kb resolution)
- Detecting single-gene homozygous deletions invisible to the backbone (e.g. focal
  CDKN2A deletion, ~28 kb)
- Generating per-sample scatter plots for visual QC and clinical review
- Producing cohort-level genemetrics TSVs for CN frequency analysis

**Purity routing:**

| Condition | Integer calling method |
|---|---|
| `purity` not supplied | Log₂R thresholds (amp ≥0.7, gain ≥0.2, loss <−0.25, del <−1.1) — CNVkit `fc65941d` defaults |
| `purity` ≤ 0.40 | Log₂R thresholds — PURPLE estimate unreliable for WGD at low purity |
| `purity` > 0.40 | `cnvkit.py call --purity --ploidy` (absolute CN correction) |
| NO_TUMOR sample | Log₂R thresholds — do not pass sentinel purity=1.0 |

## What are the inputs?

| Input | Class | Required | Description |
|---|---|---|---|
| `tumour_bam` | file | ✅ | Tumour BAM (chr-prefix) |
| `tumour_bai` | file | ✅ | BAM index |
| `sample_id` | string | ✅ | Output file stem |
| `cn_reference` | file | ✅ | PoN reference from `eggd_cgp-cnvkit-pon` |
| `baits` | file | ✅ | Gene panel BED (same as used for PoN build) |
| `purity` | float | ➖ | PURPLE purity (0–1); enables purity-corrected integer calling when > 0.40 |
| `ploidy` | float | ➖ | PURPLE ploidy; accepts raw float output (e.g. 3.1); used with purity (default: 2) |
| `sample_sex` | string | ➖ | `male` or `female`; corrects chrX/Y in diagram plot |
| `drop_low_coverage` | boolean | ➖ | Filter very-low-coverage bins (default: true) |

## What are the outputs?

| Output | Class | Description |
|---|---|---|
| `copy_ratios` | file | Per-bin log₂ copy ratios (`.cnr`) |
| `segments` | file | CBS segments (`.cns`) |
| `call_segments` | file | Integer CN calls (`.call.cns`) |
| `scatter_png` | file | Genome-wide scatter plot (`.scatter.png`), y = ±3 |
| `diagram_pdf` | file | Chromosome ideogram with gene labels (`.diagram.pdf`) |
| `genemetrics` | file | Per-gene gain/loss summary (`.genemetrics.tsv`) |

## How to run this app from the command line?

```bash
dx run eggd_cgp-cnvkit-batch \
  -itumour_bam=file-xxxx \
  -itumour_bai=file-xxxx \
  -isample_id="25330S0047" \
  -icn_reference=file-xxxx \
  -ibaits=file-xxxx \
  -ipurity=0.79 \
  -iploidy=3.1 \
  -isample_sex="male" \
  -idrop_low_coverage=true \
  --destination "project-xxxx:/cnvkit/results/25330S0047/" \
  --instance-type mem1_ssd1_v2_x4 \
  --priority high \
  -y
```

Typical runtime: **15–20 min** on `mem1_ssd1_v2_x4` (dominated by CBS segmentation ~10 min).
All samples can be submitted simultaneously.

## Dependencies

**Runtime environment:** pre-built Docker image (`cgp-cnvkit:1.0.0`) loaded from DNAnexus at job
start (~15 s). No internet access required at runtime.

- **CNVkit:** `0.9.13` (git `fc65941d`), Python `3.12`, installed in `/opt/cnvkit` venv
- **R / DNAcopy:** R `4.3.3` + DNAcopy `1.76.0` as pre-compiled Ubuntu 24.04 debs
- **Base image:** `ubuntu:24.04`

## Known limitations

- **chrX artefact (male samples):** `cnvkit.py scatter` does not support a
  `--sample-sex` flag. Male samples will show a broad chrX band at log₂R ≈ -1 due to
  the mixed-gender PoN median at chrX intervals. This is a cosmetic artefact and does not
  affect autosomal CN calls. The diagram plot (`diagram_pdf`) is sex-corrected via
  `--sample-sex` and does not show this artefact.
- **Tumour PoN attenuation:** When the PoN is built from tumour samples (no matched
  normals), recurrent CN events (e.g. EGFR amplification in >30% of cohort) are absorbed
  into the PoN median and appear at reduced log₂R compared to PURPLE purity-corrected
  values. Replace with a normal-sample PoN when matched-protocol normal samples are
  available.
- **No cnLOH detection:** CNVkit is a coverage-only tool. Copy-neutral LOH (biallelic
  loss at constant total CN) requires joint BAF + depth segmentation (e.g. FACETS).
