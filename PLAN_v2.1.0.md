# Plan: `eggd_cgp-cnvkit-batch` v2.1.0

**Jira:**
- DI-4094, with subtasks DI-4096 to DI-4104, DI-4121 and DI-4139 (dropout mask);
- review comments of 2026-10-20 on each subtask and on DI-4094;
- decisions D1–D7 and the purity/ploidy review: DI-4094 comment 72470 (2026-10-25).

**Basis:** the Atlas CGP CNVkit validation, final configuration `D_gcclip_on_s2`.
- `FINAL_RESULTS.md` is the source of truth: `project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JF6YBbj4YggvFK194V51V7pX` (repo `woook/cnv-backbone-purple-atlas`, `analysis/pon_comparison/`).
- The technical validation page (Confluence `4837998596`) has the complete list of production changes.
- Dropout mask: Confluence `4850548748` (probe dropouts and the dropout mask); `SIM2_RESULTS.md` §5, `project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JFJQqfQ4YggqqZJFzFk0vJK7` (commit `7b12192`).
- **Non-UMI re-validation** (DI-4118 comments 72480, 72481; `analysis/pon_comparison/nonumi/README.md` @ `e8bb628`): production now masks the UMIs at demultiplexing and aligns with `sentieon-bwa`, with duplicates marked (DI-4114). The final configuration `D_gcclip_on_s2` re-run on the non-UMI data of the same 226 sequencings passes **9/9** acceptance criteria. **The dropout mask has not yet been evaluated on non-UMI data**: that analysis runs in `cnv-backbone-purple-atlas` and gates the release (§5, §8).
- **Static PoN (controlled file):** `atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn`, Confluence `4854120471` (DRAFT). Non-UMI, 124 members, `file-JFP3fzQ455KjXVZ23JgX1gF2`, md5 `b4a373729d9282ce0a0ed80fd29f6d0a`. It is currently in 004 under the build name `cgp_cnvkit_static_pon_di4093_nonumi.cnn`. It will be renamed and moved to 001 after the page is approved. It supersedes the UMI-consensus PoN (`02e1c85e…`).

**Status:** PLAN, not started; decisions D1–D7 accepted (2026-10-25). Current code: `main` at `7559cc2` (v2.0.2).

---

## 1. Goal
Make the app run exactly the validated read-out with the dropout mask (arm M), and add the sex and dropout QC outputs. The CNVkit image is unchanged: `cgp-cnvkit:1.0.0`, CNVkit `fc65941`, `project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8j7Vyj45FG1BbK26JQgQY6q`.

**Validated per-sample read-out:**
```
coverage   tumour.bam panel.bed                                   -> <s>.targetcoverage.cnn   (unchanged; = coverage app output)
fix        <s>.targetcoverage.cnn empty.antitarget reference.cnn  [edge correction ON: no --no-edge]  --diploid-parx-genome grch38 -> <s>.cnr
dropout mask  dropout_mask.py mask M <s>.cnr <s>.targetcoverage.cnn <s> - <s>.seg_input.cnr <s>.dropout_qc.tsv   (DI-4139)
segment    <s>.seg_input.cnr --method cbs --drop-outliers 10 --drop-low-coverage --diploid-parx-genome grch38 -> <s>.cns
call       <s>.cns --diploid-parx-genome grch38 [-x male|female]   (NO --purity/--ploidy)          -> <s>.call.cns
genemetrics <s>.cnr -s <s>.call.cns -t 0.3 --min-probes 3 --diploid-parx-genome grch38 [-x male|female] -> <s>.genemetrics.tsv
PAR post-step (only if genemetrics treated the sample as male)                                     -> final <s>.genemetrics.tsv
sex check: cnvkit.py sex on <s>.targetcoverage.cnn and on <s>.cnr (both --diploid-parx-genome grch38)
           + sex_check logic                                                                       -> <s>.sex_qc.tsv, <s>.sex_check_mqc.tsv
```
Every step except `segment` reads the **unmasked** `<s>.cnr` (`call` reads only the `.cns`; `genemetrics`, `scatter`, `diagram`, `cnvkit.py sex` on the `.cnr`, and the PAR post-step read `<s>.cnr`).

`-x` comes from the LIMS letter (field 5 of the `-`-split sample name: M → `-x male`, F → `-x female`, anything else → no `-x`).

## 2. Changes by subtask

| Subtask | Change in `src/code.sh` / `dxapp.json` | Notes |
|---|---|---|
| DI-4096 | Delete `--no-edge` from `fix` (code.sh L63) | The PoN is built with default edge correction, so reference and fix match |
| DI-4097 | New input `diploid_parx_genome` (string, optional). If set, add `--diploid-parx-genome <value>` to `fix`, `segment`, `call`, `genemetrics`, `diagram`, both `cnvkit.py sex` runs and the post-step genemetrics. If empty, omit everywhere | **The Atlas config must pass `grch38`**; with it empty, every sample is inferred female. Decision D1: default empty (as asked in DI-4097) or `grch38` |
| DI-4098 | `LIMS_SEX=$(parse field 5 of sample_id)`. Fall back to the BAM filename stem if `sample_id` lacks the field; record the source | Accept only `M`/`F`; anything else → U |
| DI-4099 | Sex override input. See decision D2 | **Must not be PURPLE's inferred sex.** Today `stage-cnvkit_batch.sample_sex` is wired from PURPLE |
| DI-4100 | `APPLIED_SEX` = override if set, else LIMS M/F, else none. Add `-x $APPLIED_SEX` to `call` and `genemetrics` (and `diagram`) | `fix` and `segment` have no `-x` in fc65941 |
| DI-4102 | Always run `cnvkit.py sex <s>.targetcoverage.cnn … -o <s>.sex_raw.tsv` and `cnvkit.py sex <s>.cnr … -o <s>.sex_cnr.tsv` | The raw run gives chrY presence (QC flag); the `.cnr` run gives chrX context and the U-sample post-step decision |
| DI-4101 | Post-step if and only if treated as male: `-x male`, or no `-x` and `sex_cnr.tsv` says Male. Use the **`.cnr`** result, never the raw one. Inference always resolves to M/F (indeterminate → Female, CNVkit #360) | See §3. The analysis script `genemetrics_parfix.py` is batch-oriented (plans, arm dirs, a hard-coded workdir), so the logic is ported to a per-sample resource script and proven equivalent (§5) |
| DI-4103 | Run `sex_check.py` (**not** `sex_qc.py`): `--raw-sex <s>.sex_raw.tsv --cnr-sex <s>.sex_cnr.tsv --out <s>` | Source: commit `31a62fc`, `project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JF6X7g04YggpgK1QB942gY64`. stdlib-only Python. `sex_mismatch` = LIMS vs raw chrY presence (log-ratio > −5) |
| DI-4104 | New outputs (§4) | Decision D3: always emit `genemetrics.noparfix.tsv` (identical copy for females) |
| DI-4121 | Remove purity/ploidy from `call` entirely (delete the `CALL_ARGS` block); remove or ignore the `purity`/`ploidy` inputs. See D4 | Validated = no purity. A future hom/het label would be a separate `call -m clonal` output; out of scope for 2.1.0 |
| DI-4139 | After `fix`, run `dropout_mask.py mask M <s>.cnr <s>.targetcoverage.cnn <s> - <s>.seg_input.cnr <s>.dropout_qc.tsv`. `segment` reads `<s>.seg_input.cnr`, with all other options unchanged. New input `dropout_mask` (boolean, default `true`; D7). When it is `false`, `segment` reads `<s>.cnr` (= v2.0.x / arm A0). Log the masked-probe count | See §3a. Rule: autosomal, `.cnr` log2 < −3, in a run of ≤ 2 adjacent deep probes, < 10 raw read-equivalents (depth × length / 144). chrX/Y are never masked. No site file. Only the `segment` input changes |
| DI-4094 comment | Keep `drop_low_coverage` default `true`; add the equivalence test (§5) | — |
| (consistency) | `diagram`: use `-x $APPLIED_SEX` (or `--sample-sex`) and `--diploid-parx-genome`. The scatter plot is unchanged | So the plot matches the calls |

## 3. PAR post-step (per-sample port)
New resource script `resources/usr/local/bin/cnvkit_parfix.py` (stdlib only):
1. Read `<s>.cnr` and `<s>.call.cns`. Copy them to `parfix/`, and in the copies subtract 1.0 from `log2` for every bin or segment **wholly inside** GRCh38 chrX PAR1 `[10000, 2781479)` or PAR2 `[155701382, 156030895)` (0-based half-open).
   - Leave bins or segments that cross a PAR boundary unchanged, and log them (`par_span`).
2. `code.sh` re-runs genemetrics in the image on the copies, with the original arguments (`-t 0.3 --min-probes 3 --diploid-parx-genome grch38 -x male` or no `-x`).
3. Keep the original as `<s>.genemetrics.noparfix.tsv`; the result becomes `<s>.genemetrics.tsv`.
4. Assert that every non-PAR row is byte-identical to the original. Fail the job if not.

The logic is copied from `genemetrics_parfix.py` (unchanged since `f933f7a`: `project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JF6PpjQ4Yggx4JxXFkQqZxGb`). Record the source commit in the script header.

**Known limitation:** CRLF2 can still be called in males when its segment spans the PAR1 boundary. All male CRLF2 calls need review (FINAL_RESULTS §7 item 17).

## 3a. Dropout mask (DI-4139)
New resource script `resources/usr/local/bin/dropout_mask.py` (stdlib only):
- **Source:** the `mask` function and constants of `scripts/deadprobe_mask.py` at `ddad79b`:
  - file: `project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JFGy8Vj4Yggp4Zj50GVVzF2Q`, md5 `382acdecb5be92d6074f8957362bdfd9`;
  - constants `DEEP=-3.0`, `RUN_MAX=2`, `RAW_READS=10`, `READ_LEN=144`.
- **Scope:** only arm `M` is accepted. Remove the `blacklist` mode and the `L5`/`L3`/`MB` arms, but keep the positional CLI, so the command is identical to the validation. The header cites the source commit and the evaluation page.
- **Behaviour:**
  1. Copies the `.cnr` with `weight` set to `0` on masked bins. CNVkit `segment` drops weight-0 bins before CBS, the same operation as `--drop-low-coverage`.
  2. Writes `<s>.dropout_qc.tsv` with the columns `arm chromosome start end gene log2 raw_reads reason`. It writes the header only when nothing is masked.
- **Raw depth:** read from column 5 (`depth`) of the `coverage` output. A bin absent from the `.cnn` counts as 0 reads, as in the validation.
- **Validated effect** (226 sequencings):
  - UMI-consensus (round 2c): 82 probes masked in 42 sequencings, plus 2,091 bins in the QC-failed library 25330S0104; calls change in 4 sequencings: 25330S0048 +14 genes, 25330S0109 +5, 24260S0018 +3, 25330S0044 −FGFR3;
  - non-UMI (production; DI-4139 comment 72496): 78 probes masked in 41 sequencings, plus 25330S0104; calls change in 3 sequencings: 25330S0048 +17, 25330S0109 +5, 24260S0018 +3. The 25330S0044 FGFR3 call stays (probe 1 has 11.2 raw reads). Rule and constants unchanged;
  - all other sequencings are byte-identical to the final configuration.

## 4. Outputs (`dxapp.json` `outputSpec`)
| Output | Class | Always? | Content |
|---|---|---|---|
| `copy_ratios`, `segments`, `call_segments`, `genemetrics`, `scatter_png`, `diagram_pdf` | file | as now | `genemetrics` = the post-step read-out |
| `genemetrics_noparfix` (new) | file | yes (D3) | uncorrected genemetrics; identical to `genemetrics` when no post-step ran |
| `sex_qc` (new) | file | yes | `<s>.sex_qc.tsv` from `sex_check.py`, plus added fields: `applied_sex`, `applied_sex_source` (override / sample_name / cnvkit_inference), `par_fix_applied`, `par_span_rows` |
| `sex_check_mqc` (new) | file | yes | `<s>.sex_check_mqc.tsv`, MultiQC custom content (`id: cnvkit_sex_check`) |
| `dropout_qc` (new; DI-4139) | file | yes | `<s>.dropout_qc.tsv`: one row per masked probe (arm, chromosome, start, end, gene, log2, raw_reads, reason; `arm` is always `M`, D10). Header-only if none, or if `dropout_mask=false` |
| `target_coverage` (new, optional; D5) | file | yes | `<s>.targetcoverage.cnn`, so the raw-coverage sex check and future PoN rebuilds can be audited |

`sex_check.py` writes the base columns. Add the extra fields with a small wrapper, or by extending the script in a new commit of the analysis repo. Keep the MultiQC file format unchanged, apart from one added column, `dropout_probes` (the masked-probe count, DI-4139), in `sex_check_mqc`.

`segments` (`<s>.cns`) is segmented on the masked input. `copy_ratios` (`<s>.cnr`) is always the **unmasked** file. Don't output `<s>.seg_input.cnr`.

## 5. Testing
The detail is in `PLAN_v2.1.0_build.md` §0b and §3–§4. In summary, there are four **baselines**:

| ID | Data | Configuration | Location (004) | Status |
|---|---|---|---|---|
| U-A0 | UMI-consensus | Final configuration, no mask | `/final_config/per_sample/<run>/` | Available |
| U-M | UMI-consensus | Arm M (dropout mask); `genemetrics.tsv` is **before** the PAR post-step, so compare it with the app's `genemetrics.noparfix.tsv` (build plan §0b) | `/sim2c/real/sim2c_real__<run>.calls.tar.gz` (`calls/M/d00/`) | Available |
| N-A0 | Non-UMI (production) | Final configuration, no mask | `/final_config_nonumi/outputs/final_config_nonumi.D_gcclip_on_s2.parfix.outputs.tar.gz` (`nonumi/phaseN/calls/D_gcclip_on_s2/d00/`) | Available |
| N-M | Non-UMI (production) | Arm M | `/final_config_nonumi/dropout/parfix/dropout_nonumi_NM.calls.tar.gz` (`file-JFPX1JQ4763FK094J6vqj1yQ`; `calls/M/d00/`; `calls/A0/` = N-A0) | Available (DI-4139 comment 72496) |

- **Code correctness** (byte-identity) is proven on every available baseline, each sequencing with its validated PoN (leave-out for members, shared for non-members).
- The **production read-out** is the non-UMI data with the controlled PoN.

1. **Local** (Docker, no DNAnexus), 16 sequencings (`tests/harness/local_samples.tsv`):
   - match U-A0/U-M, N-A0 and N-M: `.cnr`, `.cns`, `.call.cns` and `genemetrics.tsv` byte-identical;
   - `dropout_qc` equal to `mask_qc.tsv`;
   - the post-step applied to the same samples;
   - the `sex_qc` chrY fields equal to the validation;
   - one run from a **non-UMI BAM**, with its `.targetcoverage.cnn` byte-identical to the `eggd_cgp-cnvkit-coverage` 2.0.2 output the PoN was built from (`/coverage_nonumi/<run>/`).
   - **Dropout-mask unit test:** `dropout_mask.py` against `deadprobe_mask.py mask M` (`ddad79b`) on all 226 UMI **and** all 226 non-UMI `.cnr`/`.targetcoverage.cnn` pairs: masked `.cnr` and QC rows byte-identical. UMI totals: 82 probes / 42 sequencings, plus 25330S0104. Non-UMI totals: 78 probes / 41 sequencings, plus 25330S0104 (the N-M `mask_qc.tsv` files).
2. **Applet build** in a 003 test project: the T1–T15 matrix (build plan Phase 7), from non-UMI BAMs except the UMI dropout-mask cases.
   - Check the outputs and the job logs: `Correcting for density bias`, `--diploid-parx-genome`, `-x` and the sex source logged, the post-step decision logged.
3. **Equivalence test, full** (DI-4094 comment; part of DI-4118):
   - run the 226 through the applet or app in `project-JF3p1f04Yggp0vf3PYjy7Q0F`, with the final inputs (`diploid_parx_genome=grch38`, no purity/ploidy, `drop_low_coverage=true`);
   - **runs:**
     - N-A0 (`dropout_mask=false`, validated non-UMI PoNs);
     - U-M (`dropout_mask=true`, validated UMI PoNs);
     - N-M (**release gate**);
     - production: non-UMI, controlled PoN, mask on; differences expected for PoN members only;
   - **acceptance:**
     - all four call files byte-identical, or every difference explained;
     - `dropout_qc` rows equal the `mask_qc.tsv` rows;
     - the post-step applied to exactly the same samples;
     - `sex_mismatch` = 1 in 0/173 LIMS M/F;
     - chrY presence matching.
4. **Regression on the sex inputs:** override M/F vs sample-name M/F vs U; and a bad override value, which must fail.

## 6. Release (team lifecycle)
1. Branch `dev_2.1.0` from `main`. Bump `dxapp.json` version to `2.1.0`, and update the input/output help text.
2. README:
   - the commands as run;
   - sex handling and its sources;
   - the PAR post-step and the CRLF2 limitation;
   - the sex QC outputs;
   - the dropout mask: rule, segment-only scope, the `dropout_qc` output, the `dropout_mask` input, and the limitation that real 1–2-probe events meeting the rule are masked (listed in `dropout_qc`);
   - no purity;
   - `diploid_parx_genome` must be `grch38` for Atlas.
3. Changelog `[2.1.0]`, listing the changes under Changed / Added / Removed, with the Jira keys.
4. Resource scripts with a header citing their source commit and validation.
5. PR to `dev`, then `main` (pr-workflow skill; review-fix loop).
6. Build and publish the app (`dx build --app`, `dx publish eggd_cgp-cnvkit-batch/2.1.0`).
7. Update the Confluence dev-doc for the app (see `memory/confluence_doc.md`): the test evidence from §5 and links to the validation pages.
8. Update `memory/version_state.md`.

## 7. Decisions (D1–D7 accepted 2026-10-25, DI-4094 comment 72470; D8–D11 accepted 2026-10-25, see `PLAN_v2.1.0_build.md` §1)
| ID | Decision | Accepted option |
|---|---|---|
| D1 | `diploid_parx_genome` default: empty (DI-4097) or `grch38` | Empty in the app, with `grch38` set explicitly in the Atlas config, plus a loud log warning when empty |
| D2 | Sex override input: reuse `sample_sex` or add a new name (e.g. `sex_override`) | New name. Then `sample_sex` (fed from PURPLE) can be deprecated: ignored, with a warning. It also avoids PURPLE sex silently overriding LIMS before the workflow is updated |
| D3 | `genemetrics_noparfix` for females | Always emit (identical copy) |
| D4 | `purity`/`ploidy` inputs: remove, or keep as deprecated and ignored | Keep, ignored with a log warning, for 2.1.0. That stays compatible with the current workflow wiring. Remove in the next major version, once the workflow drops the wiring |
| D5 | Also output `<s>.targetcoverage.cnn` | Yes (small file; useful for audit and PoN rebuilds) |
| D6 | Where the resource scripts come from | Bundle in `resources/` with pinned source commits (`sex_check.py` @ `31a62fc`, post-step logic @ `f933f7a`, `dropout_mask.py` from `deadprobe_mask.py` @ `ddad79b`) |
| D7 | `dropout_mask` input: default `true`, or no input (always on) | Input, default `true`. `false` reproduces v2.0.x segmentation for the equivalence test and for troubleshooting. The Atlas config leaves it at the default |
| D8 | How to reconcile `sex_check.py` with the app's sex sources | Port it to a per-sample `cnvkit_sex_qc.py`, with the `sex_check.py` columns and values unchanged when there is no override (proved on 226). `applied_sex`/`applied_source` come from the app's resolution (`override` / `sample_name` / `bam_name` / `cnvkit_inference`). Append `par_fix_applied`, `par_span_segments` and `dropout_probes`. `sex_mismatch` = declared sex (override, else LIMS) vs chrY. The MultiQC file adds only `dropout_probes` |
| D9 | Accepted override values | `male`, `female`, `M` or `F`, case-insensitive; empty = none. Anything else fails the job before compute |
| D10 | `arm` column in `dropout_qc` | Keep it (always `M`), so the file is byte-identical to the validation `mask_qc.tsv` (8 columns) |
| D11 | Allowed `diploid_parx_genome` values | `choices: ["grch38"]`, optional, default empty (D1). When empty, skip the PAR post-step (GRCh38 coordinates) and log a warning |

## 8. Dependencies outside this app (not in DI-4094; blocking for production)
- **Static PoN** (DI-4093): **done (non-UMI).** `atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn` (Confluence `4854120471`, DRAFT; `file-JFP3fzQ455KjXVZ23JgX1gF2`, md5 `b4a37372…`), built with `eggd_cgp-cnvkit-pon` v2.2.0 (GIABv3 FASTA trio, clip 0.3/0.7, `diploid_parx_genome=grch38`) from the same 124 members (non-UMI coverage).
  - Remaining: approval of the page, then the rename and the move to 001.
  - The workflow must reference the 001 copy: check its file ID and md5 after the move.
  - The UMI-consensus PoN (`file-JFB91BQ49ZYfv92jqv6qzxpf`) is superseded and must not be used in production.
- **Non-UMI dropout-mask evaluation** (`cnv-backbone-purple-atlas`; **done 2026-10-26**, DI-4139 comment 72496): arm M meets every T3 margin on the 226 non-UMI sequencings.
  - Rule and constants unchanged (`deadprobe_mask.py` @ `ddad79b`), so D6 and D7 stand.
  - The masked-probe set changed only at the < 10-read limit: 78 vs 82 probes. FGFR3 probe 1 is mostly above the limit at non-UMI depth, so the 25330S0044 FGFR3 call is not removed.
  - The N-M baseline is published (§5).
- **`eggd_atlas_cnv_workflow` new version:**
  - the batch app 2.1.0;
  - the static PoN as `cn_reference`;
  - **remove the PURPLE → batch `purity`/`ploidy`/`sample_sex` wiring**;
  - expose `sex_qc`, `sex_check_mqc`, `genemetrics_noparfix`, `dropout_qc` (and `target_coverage`).
- **Conductor CGP config:**
  - drop the per-run PoN (analysis_7);
  - the new workflow version;
  - `diploid_parx_genome=grch38`;
  - fix the 18 malformed `$dnanexus_link` inputs.
- **MultiQC:** collect `*_mqc.tsv` from the CNV workflow outputs.
- **Assay design (panel owners):**
  - flag the two ETV6 probes inside gnomAD-SV DEL_chr12_9b8e999a in the bins BED;
  - review the FGFR3 probe 1 capture.
- **CNVkit upgrade (DI-3657):** re-check the dropout sites and the mask rule.
- **Purity/ploidy (DI-4094 comment 72470):**
  - CNVkit stays an independent depth-based caller: CBS, no `--purity`/`--ploidy`, no VCF;
  - PURPLE is the source of absolute CN, hom/het, allele-specific CN and LOH in the side-by-side review;
  - report-layer purity annotations are proposed, not agreed, and are outside this app;
  - note: a gene-level three-way concordance on the WGS cases (CNVkit vs PURPLE vs DRAGEN) before go-live.
- **Validation 2** (DI-4118, "Atlas CNV validation 2 (non-UMI)"):
  - the final configuration passes 9/9 on non-UMI data (comment 72481);
  - the sex check cut-off still holds: the lowest real male is now −3.9 (UMI: −2.5) against the −5 cut-off;
  - remaining: the dropout mask on non-UMI (above), the equivalence test above, and the prospective items (single-exon events, hom/het, high-purity homozygous deletions).
