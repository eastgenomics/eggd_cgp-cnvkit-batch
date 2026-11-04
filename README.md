<!-- dx-header -->
# eggd_cgp-cnvkit-batch (DNAnexus Platform App)

Per-sample CNVkit analysis for tumour-only CGP panel sequencing. Version 2.1.0 runs the configuration validated for
Atlas CGP (`D_gcclip_on_s2` with the dropout mask; Confluence DV `4837998596` and `4850548748`, Jira DI-4094), on
production (UMI-masked, non-UMI) data with the controlled static PoN.

```
PoN build (once):  PoN specimen BAMs → eggd_cgp-cnvkit-coverage (per specimen) → eggd_cgp-cnvkit-pon → static PoN (.cnn)

Per sample:        tumour BAM + BAI, panel BED ─┐
                   static PoN (cn_reference) ───┴→ eggd_cgp-cnvkit-batch (×N, parallel)
```

The batch app computes its own coverage from the BAM; the coverage app only supplies the PoN build.

For Atlas the PoN is the controlled file `atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn` (Confluence DV `4854120471`;
md5 `b4a373729d9282ce0a0ed80fd29f6d0a`), built from 124 FFPE specimens with GC clipping and
`--diploid-parx-genome grch38`.

## What does this app do?

Runs the full CNVkit per-sample analysis pipeline against a pre-built PoN reference.

For one sample, each of the steps below is run. CNVkit steps run inside the `cgp-cnvkit:1.0.0` image (CNVkit
`fc65941`); additional processing uses Python scripts on the worker.

- `<s>` is `sample_id`.
- `PARX` is `--diploid-parx-genome grch38` when `-idiploid_parx_genome=grch38` (required for Atlas).
- `-x` is the applied sex (below).

```
coverage     tumour.bam panel.bed --processes N                        -> <s>.targetcoverage.cnn
fix          <s>.targetcoverage.cnn empty.antitarget reference.cnn PARX -> <s>.cnr     (edge correction on)
dropout mask dropout_mask.py mask <s>.cnr <s>.targetcoverage.cnn <s>.seg_input.cnr <s>.dropout_qc.tsv
segment      <s>.seg_input.cnr --method cbs --drop-outliers 10 --drop-low-coverage PARX -> <s>.cns
sex          <s>.targetcoverage.cnn PARX -> <s>.sex_raw.tsv;  <s>.cnr PARX -> <s>.sex_cnr.tsv
call         <s>.cns PARX [-x]                                         -> <s>.call.cns  (no --purity/--ploidy)
scatter      <s>.cnr -s <s>.cns --trend --y-min -3 --y-max 3 --fig-size 14 4   -> <s>.scatter.png (non-fatal)
diagram      <s>.cnr -s <s>.call.cns --threshold 0.5 --min-probes 3 PARX [-x]  -> <s>.diagram.pdf (non-fatal)
genemetrics  <s>.cnr -s <s>.call.cns -t 0.3 --min-probes 3 PARX [-x]   -> <s>.genemetrics.tsv
PAR post-step (samples treated as male)                                -> final <s>.genemetrics.tsv
sex QC       cnvkit_sex_qc.py                                          -> <s>.sex_qc.tsv, <s>.sex_check_mqc.tsv
```

Only `segment` reads the dropout-masked `<s>.seg_input.cnr`; every other step reads the unmasked `<s>.cnr`. Duplicate
BED intervals are removed before `coverage` (CNVkit ≥ 0.9.13 rejects them). The job log records each decision with a
fixed prefix: `[parx]`, `[sex]`, `[dropout]`, `[parfix]`, `[deprecated]`.

### Sex

The applied sex is passed as `-x` to `call`, `genemetrics` and `diagram`. Precedence:

1. `sex_override`: `male`, `female`, `M` or `F` (case-insensitive). Any other value fails the job before analysis.
2. The LIMS letter: field 5 of the `-`-split `sample_id` (e.g. `140609591-25328S0023-26CGPV2-8475-M-96770770` → `M`),
   or of the BAM file name when `sample_id` has fewer than 5 fields. `M` → `-x male`, `F` → `-x female`.
3. Otherwise, including a LIMS letter of `U`, no `-x` is supplied: CNVkit infers sex from chrX
   (`applied_sex = cnvkit_inferred(...)`). For Atlas, use `diploid_parx_genome=grch38`; if it is empty, the
   raw-coverage sex QC can call male samples female. An explicit override or M/F LIMS letter still takes precedence.

`sample_sex` (wired from PURPLE's inferred sex in the current workflow) is ignored, with a warning.

### PAR post-step

For a sample that `genemetrics` treats as male, CNVkit `fc65941` shifts every chrX bin and segment by +1.0 so a
normal male's single X does not look like a deletion. This includes the pseudo-autosomal regions, for which males
already have two copies and do not require correction, so PAR genes (on this panel, CRLF2) are reported 1.0 too high.
A sample is treated as male when it gets `-x male`, or no `-x` and `cnvkit.py sex` on the `.cnr` reports Male. To
correct the PAR issue, `cnvkit_parfix.py` subtracts 1.0 from the bins and segments **wholly inside** GRCh38 PAR1
`[10000, 2781479)` and PAR2 `[155701382, 156030895)` in copies under `parfix/`, `genemetrics` is re-run there with the
original arguments, and every non-PAR row must be byte-identical to the original. Each PAR gene row is then judged
against its own `call.cns` segment. In a segment wholly inside a PAR (shifted), the row must equal the segment log2
after the fix and the segment log2 + 1.0 before it, each at a level `genemetrics` reports (magnitude ≥ 0.3), so an
original row below the threshold is never accepted (`par_check=1`); a row that disappears must have dropped below the
0.3 threshold, or be replaced by the row with the same gene and coordinates that differs only in log2 (a same-name
row elsewhere does not count). In a segment that spans a PAR boundary (not shifted), the row must be unchanged and
equal the segment log2 + 1.0, the male shift, at magnitude ≥ 0.3 (`par_check=span`, also when other PAR genes were
corrected; the gene and segment are logged as `span_genes`); a row there that equals the segment log2 is not
accepted as corrected. Anything else is an unexplained PAR discrepancy (`par_check=0`). "Equal" allows for CNVkit
printing every value to 6 significant digits (`%.6g`): the difference may be up to half a unit in the 6th digit at
the value's magnitude (a segment at 9.99996 gives 10.99996 before the fix, printed as `11`). CNVkit writes
overlapping segments for overlapping targets, and `genemetrics` reports a gene once for each segment that overlaps
its bins, at that segment's log2; a PAR row is therefore explained when any segment it overlaps explains it (with
non-overlapping segments, that is the one segment containing it). `compare` reads the two genemetrics tables with
only the checks it needs to compare them: identical headers starting `gene chromosome start end log2` (CNVkit's
header-only table, when no gene is reported, may stand on either side), one field per column, integer coordinates
and a numeric log2; an empty or truncated table fails. The job fails, with no outputs, on such a table, when non-PAR
rows differ or when `par_check` is 0. The corrected read-out is `genemetrics`; the uncorrected one is always written
as `genemetrics_noparfix`. Without `diploid_parx_genome` no post-step runs.

The `cnvkit.py sex` tables are checked before use, by the sex QC and by the PAR decision: one row naming
`<s>.targetcoverage.cnn` / `<s>.cnr`, every required column present once, a value for each, a sex call CNVkit writes
(`Male`, `Female`, or `Unknown` when the sample has no usable chrX bins) and finite numeric log-ratios, and a newline
at the end of the last line, as CNVkit writes it (a raw Y log-ratio `-5.01e-05` cut to `-5.01` would still be a
number). CNVkit's literals for no value are accepted: `nan` (no usable bins on that chromosome; a `nan` raw Y
log-ratio means no chrY) and `NA`, which it writes with `Unknown`; both read as no value, as before. Any other value
fails the job instead of being read as female or as no chrY, or skipping the PAR post-step.

### Dropout mask

A probe dropout is 1–2 adjacent probes with almost no reads among normal neighbours. Probe dropouts inflate the CBS
noise estimate and hide focal events on the same chromosome (e.g. the common inherited ETV6 deletion, gnomAD-SV
DEL_chr12_9b8e999a, impairs chr12 calling). To correct this, probe dropouts are masked before being passed to
`segment`. During masking, an **autosomal** probe is assigned a weight of 0 in the segmentation input when all of the
following criteria hold: `.cnr` log2 < −3; it is in a run of ≤ 2 adjacent probes with log2 < −3 (longer runs are
callable deletions and are kept); fewer than 10 raw read-equivalents (depth × length / 144) in
`<s>.targetcoverage.cnn`. chrX/chrY are never masked. `dropout_qc` lists the masked probes. `dropout_mask=false`
segments the unmasked `.cnr`: the validated final configuration without the mask (arm A0), not v2.0.x (v2.1.0 keeps
edge correction, `--diploid-parx-genome` and the other changes).

### Copy-number calls

`call` uses CNVkit's log2 ratio thresholds to classify copy-number events (`fc65941` defaults: amplification ≥ 0.7,
gain ≥ 0.2, loss < −0.25, deletion < −1.1). No purity or ploidy correction is applied: CNVkit stays an independent
depth-based caller, and within the Atlas workflow PURPLE is the source of absolute copy number, hom/het and LOH in
the side-by-side review (DI-4121).

## What are the inputs?

| Input | Class | Required | Description |
|---|---|---|---|
| `tumour_bam` | file | ✅ | Tumour BAM (chr-prefix) |
| `tumour_bai` | file | ✅ | BAM index (`<bam>.bai` or `<stem>.bai`) |
| `sample_id` | string | ✅ | Output file stem; field 5 is the LIMS sex letter |
| `cn_reference` | file | ✅ | PoN reference (`.cnn`) built with the same BED |
| `baits` | file | ✅ | Panel BED (bins BED v1.1 for Atlas), identical to the PoN build |
| `diploid_parx_genome` | string | ➖ | `grch38` (**required for Atlas**) or empty. Empty: no `--diploid-parx-genome`, no PAR post-step, and a warning (raw-coverage sex QC can call male samples female) |
| `sex_override` | string | ➖ | `male` / `female` / `M` / `F`; wins over the LIMS letter. Leave empty in routine use |
| `drop_low_coverage` | boolean | ➖ | `--drop-low-coverage` on `segment` (default `true`) |
| `dropout_mask` | boolean | ➖ | Mask isolated probe dropouts before segmentation (default `true`) |
| `purity`, `ploidy` | float | ➖ | **Deprecated, ignored** (warning logged); kept for the current workflow wiring |
| `sample_sex` | string | ➖ | **Deprecated, ignored** (warning logged); use `sex_override` |

## What are the outputs?

| Output | Description |
|---|---|
| `copy_ratios` | `<s>.cnr`: per-bin log2 copy ratios (unmasked) |
| `segments` | `<s>.cns`: CBS segments (segmented on the dropout-masked input) |
| `call_segments` | `<s>.call.cns`: segments with copy-number calls |
| `genemetrics` | `<s>.genemetrics.tsv`: genes with \|log2\| ≥ 0.3 and ≥ 3 probes, after the PAR post-step |
| `genemetrics_noparfix` | `<s>.genemetrics.noparfix.tsv`: the read-out before the PAR post-step (identical when it did not run) |
| `sex_qc` | `<s>.sex_qc.tsv`: LIMS sex, applied sex and source, chrY presence (raw Y log-ratio > −5), X/Y log-ratios, `sex_mismatch`, `par_fix_applied`, `par_span_segments`, `dropout_probes` |
| `sex_check_mqc` | `<s>.sex_check_mqc.tsv`: MultiQC custom content (`id: cnvkit_sex_check`) |
| `dropout_qc` | `<s>.dropout_qc.tsv`: masked probes (`arm chromosome start end gene log2 raw_reads reason`); header only when none |
| `target_coverage` | `<s>.targetcoverage.cnn`: raw on-target coverage (audit, PoN rebuilds) |
| `scatter_png` | `<s>.scatter.png`: genome-wide scatter, y = ±3 (optional: plot failures are not fatal, and a partly written plot from a failed command is not uploaded) |
| `diagram_pdf` | `<s>.diagram.pdf`: chromosome ideogram with gene labels (optional, as `scatter_png`) |

Within `<s>.sex_qc.tsv`, **`sex_mismatch`** is 1 when the declared sex (override, else LIMS) disagrees with chrY
presence in the raw coverage: review the sample (possible sample swap, or a sex-chromosome variant). `NA` when the LIMS
letter is not M/F and there is no override.

## How to run this app from the command line?

```bash
dx run eggd_cgp-cnvkit-batch \
  -itumour_bam=file-xxxx -itumour_bai=file-xxxx \
  -isample_id="140609591-25328S0023-26CGPV2-8475-M-96770770" \
  -icn_reference=file-xxxx \
  -ibaits=file-xxxx \
  -idiploid_parx_genome=grch38 \
  --destination "project-xxxx:/cnvkit/" -y
```

All samples can run in parallel (one job per sample).

## Dependencies

- Docker image `cgp-cnvkit:1.0.0` (`project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8j7Vyj45FG1BbK26JQgQY6q`), loaded at job
  start; no internet access needed. CNVkit `0.9.13` (git `fc65941d`), Python 3.12, R 4.3.3 with DNAcopy 1.76.0.
- Resource scripts (standard-library Python on the worker), each a port of the validation code with its source
  pinned in the header:

  | Script | Ported from (cnv-backbone-purple-atlas) |
  |---|---|
  | `dropout_mask.py` | `deadprobe_mask.py` @ `ddad79b`, arm M |
  | `cnvkit_parfix.py` | `genemetrics_parfix.py` @ `f933f7a` |
  | `cnvkit_sex_qc.py` | `sex_check.py` @ `31a62fc` |

- `resources/usr/local/lib/cgp_cnvkit_batch/pipeline.sh`: the analysis.
- `src/code.sh`: stages inputs and uploads outputs.

## Known limitations

- **CRLF2 in males:** when CBS merges the PAR1 bins into a segment that continues into non-PAR chrX, the PAR
  post-step leaves it unchanged, and the PAR gene keeps the +1 shift (`par_span_segments` > 0, logged). Review every
  male CRLF2 call.
- **Dropout mask:** a real 1–2-probe event that meets the rule (log2 < −3, < 10 reads) is masked from segmentation; it
  is listed in `dropout_qc`. **FGFR3 probe 1** is a capture dropout that the mask does not remove at non-UMI depth
  (10 of 11 affected specimens have 11–30 raw reads), so FGFR3 loss calls from short segments need review.
- **No absolute copy number:** no purity/ploidy correction (by design). CNVkit log2 is relative to the sample's own
  average copy number: in a near-tetraploid tumour a loss at log2 −0.4 can still be ≥ 3 copies.
- **chrX in the scatter plot:** `cnvkit.py scatter` has no sex option, so male samples show chrX at log2 ≈ −1. The
  diagram and the calls are sex-corrected.
- **Tumour PoN attenuation:** the PoN is built from FFPE tumour specimens; recurrent events are partly absorbed into
  the PoN median and appear at reduced log2.
- **No copy-neutral LOH:** CNVkit uses depth only.

## Testing

Prerequisites:

- Python 3.12 and `uv`, installed into the venv below (no system `uv` needed).
- ShellCheck for L0. Without it the shellcheck tests are **skipped** locally (CI installs a pinned version and fails
  if it is missing).
- For L3/L4: Docker, with the pinned image loaded as `cgp-cnvkit:1.0.0`
  (`dx download project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8j7Vyj45FG1BbK26JQgQY6q -o cnvkit-image.tar.gz &&
  docker load < cnvkit-image.tar.gz`), and the validation data.
- To fetch the data, for the BAM start and for the DNAnexus validation runs (`tests/dx/validation_jobs.py`): `dx login`
  with access to `project-JF3p1f04Yggp0vf3PYjy7Q0F`.

```bash
python3 -m venv .venv && .venv/bin/python -m pip install --require-hashes --no-deps --only-binary :all: \
    -r requirements-bootstrap.txt                              # uv, by artifact hash
.venv/bin/uv pip install --python .venv/bin/python --require-hashes -r requirements-test.txt   # hash-pinned
.venv/bin/python -m pytest -m "static or unit or pipeline" -rs # L0–L2: no Docker, data or DNAnexus (CI)
.venv/bin/python tests/validation/fetch_validation_data.py all # ~3 GB, pinned by file ID and sha256
.venv/bin/python -m pytest -m validation -rs                   # L3 + L4: needs Docker and cgp-cnvkit:1.0.0
.venv/bin/python tests/harness/run_local.py --data nonumi --start-from bam --samples 25328S0023   # BAM start (4.8 GB)
```

Tests whose prerequisite is absent (Docker or the image, the validation data, ShellCheck) are **skipped, not
failed**, and pytest still exits 0. Read the `-rs` skip summary, and do not take a run as a pass unless the layers you
need ran with no skips (L3 + L4: 14 tests). The one exception: a Docker check that timed out and could not confirm
its own cleanup (a `docker` process may remain) fails the Docker tests instead of skipping them.

| Layer | What it proves |
|---|---|
| L0 static | shellcheck; `dxapp.json` inputs read and outputs produced; the CI selection never calls `docker`; the resource scripts use no `assert` (`python -O` would remove it) |
| L1 unit | The three resource scripts on synthetic data, including byte-identity with the pinned originals and the malformed-input failures; the test helper commands (`docker`, `dx`, `git`) are bounded, and on timeout the whole process group is killed, the original timeout raised with its output and any unconfirmed cleanup reported (a Docker check whose cleanup is unconfirmed fails rather than skips). Any exception during cleanup (invalid UTF-8 from Docker included, or from the L4 harness's group kill) is reported as unconfirmed cleanup, never in place of the original failure, and every named container is still removed; the L4 harness then keeps its `.cleanup-unconfirmed` marker and reports the timeout as a FAIL (exit 124). A pipe that fails to close after the kill is recorded on the original exception, but is not unconfirmed cleanup, so it alone does not fail a test that would pass or skip. A failing `dx` is first on `PATH` in the unit tests that mock it, so none can reach DNAnexus. One deferred exception: the L4 harness's `--start-from bam` BAM/BAI transfers and their `dx describe`/`dx find` calls have no timeout ([#5](https://github.com/eastgenomics/eggd_cgp-cnvkit-batch/issues/5)) |
| L2 pipeline | The exact CNVkit command lines for every input scenario, and the failure paths, under PATH stubs for `dx` and `docker` |
| L3 equivalence | Each port against the pinned original and the validation outputs on the 226 sequencings (UMI-consensus and non-UMI); what is compared differs per port (below) |
| L4 harness | The real image + `pipeline.sh` reproduce the validation baselines on 16 sequencings, in four data × mask configurations. The harness refuses an empty or incomplete sample list (exactly 16 sequencings, each sample ID once) and fails a sequencing with nothing compared, and each L4 test requires a PASS for all 16. The DNAnexus check (`tests/dx/validation_jobs.py check`) likewise refuses an empty or incomplete jobs file |

What L3 compares, per port. Every script runs through its command line, in a directory of its own per sequencing,
with `pipeline.sh`'s file names and arguments. Results are keyed by run and sample (the 226 sequencings have 150
distinct sample IDs: replicates share one), and each comparison must cover all 226, each once.

- `dropout_mask.py` (UMI and non-UMI, all 226): `dropout_mask.py mask`; the masked `.cnr` and the QC table are
  byte-identical to the pinned `deadprobe_mask.py mask M` and the QC table to the validated `mask_qc.tsv`, and the
  printed count equals the QC table's rows. The masked-probe totals are 82 probes in 42 sequencings (UMI) and 78 in
  41 (non-UMI), with the QC-failed library counted separately.
- `cnvkit_parfix.py` (UMI A0, non-UMI A0 and non-UMI M; all 226 in each): per sequencing, `cnvkit.py sex` on the
  `.cnr` and `genemetrics` on the unmodified inputs, then `decide`. The `cnvkit.py sex` call and the male decision
  equal the validation log, and the re-run reproduces the original read-out byte for byte. For the **60 treated
  males** in each set: `shift` (its counts equal the log), `genemetrics` in `parfix/` and `compare`; every `compare`
  statistic equals the log, and the post-step `genemetrics.tsv` is byte-identical to the validated one. `shift` and
  `compare` must print their full sets of statistics (and the log record each), so a missing statistic fails rather
  than going uncompared. `par_check` is 1 / span for 51/9 (UMI A0), 50/10 (non-UMI A0) and 50/10 (non-UMI M), never
  0; 4 span samples in each set have no PAR row. A further test, with no validation data, runs CNVkit's own
  `genemetrics` in the image on a synthetic male with overlapping PAR1 bins and two overlapping segments: the
  post-step completes with `par_check=1`, and where both segments contain the PAR row the output equals the base
  (v2.1.0 before the review fixes, whose `compare` is the pinned original's) byte for byte.
- `cnvkit_sex_qc.py` (UMI and non-UMI, all 226): per sequencing, `cnvkit.py sex` on the coverage and on the `.cnr`
  (arm A0), then `resolve`, `decide` / `shift` (for `par_fix_applied` and `par_span_segments`) and `qc`. The pinned
  `sex_check.py`, on run-prefixed tables built from those per-sample tables, reproduces the validated UMI
  `sex_qc.tsv` byte for byte. On non-UMI, the raw-coverage calls equal the DI-4118 tier 1 tables. The port is then
  compared with `sex_check.py` **field by field**, not byte for byte. Every `sex_check.py` column must match except
  three deliberate differences: `applied_source` is renamed (`LIMS` → `sample_name`, `CNVkit inference (.cnr)` →
  `cnvkit_inference`), `applied_sex` for inferred samples names the call (`cnvkit_inferred(Male)`) (D8), and `run` is
  empty, as `pipeline.sh` passes no `--run`. The port's table must have every `sex_check.py` column and its added
  columns, in order, and the added columns must equal the `decide` / `shift` values and 0 masked probes (arm A0); its
  MultiQC table is checked in L1. `sex_mismatch` is 0 for all 173 LIMS M/F sequencings.
- The pinned originals are checked by MD5: `MD5SUMS` must list exactly the three scripts, and each must match.
- L4 compares `.cnr`, `.cns`, `.call.cns`, `genemetrics(.noparfix).tsv`, `dropout_qc.tsv` and the sex QC fields with
  the baselines. In the UMI mask-on configuration, a treated male with a masked probe has no UMI post-step baseline,
  so its `genemetrics.tsv` is reported as `n/a` (one check); the non-UMI N-M run covers that path.
