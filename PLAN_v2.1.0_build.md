# Build plan: tests, then the v2.1.0 app updates

**Companion to:** `PLAN_v2.1.0.md` (what to build). This file says **how and in what order**.
**Jira:** DI-4094; subtasks DI-4096 to DI-4104, DI-4121 and DI-4139. D1–D7 accepted 2026-10-25 (DI-4094 comment 72470); D8–D11 (comment 72472). Non-UMI re-validation: DI-4118 (comments 72480, 72481).
**Starting point:** `main` at `7559cc2` (v2.0.2). There are no tests, test directory or CI. The app is a single `src/code.sh` (167 lines) that mixes dx I/O, Docker calls and logic.

**Progress (2026-10-26):** Phases 0–6 done on `dev_2.1.0` (local, not pushed), last commit `df4a412`.
- L0–L2: 229 tests pass.
- L3: the three ports are byte-identical to the originals on 226 × {UMI, non-UMI} sequencings.
- L4: 4 × 16 sequencings, 512 checks pass (1 n/a by design); the BAM-start run passes, including coverage = coverage app.
- The review-fix loop is clean after 3 rounds:
  - round 1: `sample_id` is now validated as a file stem; the reviewer's "code.sh never calls main" was rejected, because the DNAnexus executor calls `main`;
  - round 2: CNVkit now runs as the executor's user, so the PAR post-step can overwrite container-written files.

Next: Phase 7 (DNAnexus testing), which needs a go-ahead because it creates a 003 project and runs jobs.

---

## 0. Facts checked while writing this plan (2026-10-25)

### CNVkit flags
Checked against the local `cgp-cnvkit:1.0.0` image (CNVkit 0.9.13, `fc65941`):

| Subcommand | `--diploid-parx-genome` | `-x/--sample-sex` | Other |
|---|---|---|---|
| `fix` | yes | – | `--no-edge` (to be removed) |
| `segment` | yes | – | `--drop-low-coverage` |
| `call` | yes | yes | `--purity`/`--ploidy` (to be removed) |
| `genemetrics` | yes | yes | |
| `diagram` | yes | yes | |
| `sex` | yes | – | |
| `scatter` | **no** | **no** | Unchanged, as the plan says |

### Source scripts
All three are reachable in `project-JF3p1f04Yggp0vf3PYjy7Q0F:/repo_snapshot/<commit>/analysis/pon_comparison/scripts/`:

| Script | Commit | Size | md5 |
|---|---|---|---|
| `sex_check.py` | `31a62fc` | 4.2 kB | `a6106d50c43cd96ebde9ae3f24360e2b` |
| `genemetrics_parfix.py` | `f933f7a` | 19 kB | `51ada0e597420655c89c86e1a2e13ea6` |
| `deadprobe_mask.py` | `ddad79b` | 5.3 kB | `382acdecb5be92d6074f8957362bdfd9` (matches the plan) |

### Findings from reading the sources
They lead to the new decisions D8–D11 in §1.

- **`sex_check.py`** takes LIMS sex from field 5 of the sample name in the `cnvkit.py sex` table. It writes `applied_sex` and `applied_source` itself, but derives them from LIMS only: `male`/`female`/`cnvkit_inferred` and `LIMS`/`CNVkit inference (.cnr)`.
  - It knows nothing about the new override input (DI-4099) or the BAM-stem fallback (DI-4098).
  - Its column names and values differ from plan §4 (`applied_sex_source`: override / sample_name / cnvkit_inference).
  - Bundling it verbatim and adding fields afterwards would mean rewriting its columns.
- **`deadprobe_mask.py mask`** writes `dropout_qc` with an `arm` column first (8 columns). Jira and plan §4 list 7 columns.
- **`genemetrics_parfix.py`**:
  - The core logic is small: `par_state`, `shift_par` (log2 is column 4, written with `repr(float - 1.0)`), the male decision, and the `compare` assertions.
  - It re-runs genemetrics in-process (`cnvlib.commands.parse_args(argv).func`) with cwd set to a temp dir and relative file names.
  - It skips the post-step when `--diploid-parx-genome` is off (the "s0" arms).

### Housekeeping
`memory/version_state.md` was stale: it still said PR #3 was open, but that PR is merged at `7559cc2`. Corrected in Phase 0. v2.0.2 is also published and released.

### Phase 0 results (2026-10-25)

**Resolved IDs.** All are in the validation project `004_261008_cgp_cnvkit_pon_comparison` (`project-JF3p1f04Yggp0vf3PYjy7Q0F`; "004" in the plans refers to this project).

| Item | Location |
|---|---|
| Static PoN (DI-4093), UMI-consensus | `/pon/di4093_app_v2.2.0/cgp_cnvkit_static_pon_di4093.cnn`, `file-JFB91BQ49ZYfv92jqv6qzxpf`, md5 `02e1c85e3d1407f33b82f759e27f61c7`. Byte-identical to the validation shared PoN `D_GC__shared__gcclip__parx.cnn`. **Superseded by the non-UMI PoN (§0b); test use only** |
| Shared and leave-out PoNs | `/final_config/jobs_build_calls/final_main/final_main__batch_000.pons.tar.gz` (`file-JF61YzQ4P0z8QFVv534xJvXx`): `phaseB/pons/D_GC__shared__gcclip__parx.cnn` plus 124 `D_GC__minus__<specimen>__<patient>__gcclip__parx.cnn` |
| Final-config outputs (A0) | `/final_config/per_sample/<run>/<s>.{cnr,cns,call.cns,genemetrics.tsv}`, listed in `/validation_manifest/11_final_per_sample_files.csv` (`call_set=main`: 226 × 4 files) |
| Arm M outputs | **Tarballs**, not folders: `/sim2c/real/sim2c_real__<run>.calls.tar.gz`, members `calls/M/d00/<run>/<s>.{cnr,cns,call.cns,genemetrics.tsv,mask_qc.tsv}`. Arm M's `.cnr` is the **unmasked** file (md5 equal to A0 and to `final_config`); only `.cns` onwards differ |
| Sex QC truth | `/repo_snapshot/31a62fc/analysis/pon_comparison/results/phaseB/sex_check/final_config.sex_qc.tsv` (`file-JF6X7gQ4YggvqyqFBxjq6ZXF`) |
| PAR post-step log | `/final_config/outputs/parfix_log.tsv` (`file-JF61vJQ4Y2jZ2g7gpXJvKPQ5`); 226 main calls: 60 treated male (51 changed, 9 `span`), 166 not male |
| PoN membership | `/validation_manifest/03_pon_samples.csv` (124 members; `loo_pon_id` per member; one specimen per patient) |

**PoN membership.** Every LIMS M/F sequencing (173/173) belongs to a PoN specimen. Only 27 LIMS-U sequencings are non-members, so a non-member-only comparison can't test the M/F paths.

**Resolution.** The app doesn't depend on any particular PoN, so each sequencing can be run with its **validated** PoN as `cn_reference`: the matching leave-out PoN for members, the shared PoN for non-members. Then byte-identity can be required for all 226. A second run with the static PoN then shows the production behaviour; differences there are expected for members only.

**Local test set.** 16 sequencings in `tests/harness/local_samples.tsv`: 5 male including the XXY, 4 female, and 7 LIMS-U, including the three masked U samples and one PoN member. Each row has the PoN to use and the expected post-step and chrY results. Together they cover:
- each post-step path: changed, `span`, `-x` overriding the `.cnr` sex, and post-step via inference;
- each dropout case: 25330S0048, 24260S0018, 25330S0044 and 25330S0104.

The leave-out PoN names are the same in the UMI and non-UMI PoN tarballs (checked: 124/124, plus `D_GC__shared__gcclip__parx`), so the `cn_reference` column works for both data sets.

### 0b. Non-UMI data and the four baselines (2026-10-25)
Production now masks the UMIs at demultiplexing and aligns with `sentieon-bwa`, with duplicates marked (DI-4114).
- **Validation status:** DI-4118 re-ran the final configuration on the non-UMI data of the same 226 sequencings: 9/9 criteria pass (comment 72481).
- **Static PoN:** the controlled file is `atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn` (Confluence `4854120471`, DRAFT).
  - Today it is `file-JFP3fzQ455KjXVZ23JgX1gF2` in 004, under the build name `cgp_cnvkit_static_pon_di4093_nonumi.cnn`, md5 `b4a373729d9282ce0a0ed80fd29f6d0a`.
  - It will be renamed and moved to 001 after approval.
  - The tests use the 004 copy by file ID and md5. The app never hard-codes a PoN.
- **Dropout mask:** evaluated on non-UMI data (DI-4139 comment 72496; Confluence `4850548748` v6, test T5). The rule and constants are unchanged (`deadprobe_mask.py` @ `ddad79b`, md5 `382acdec…`), so D6 and D7 stand. Non-UMI totals: 78 probes / 41 sequencings, plus 25330S0104.

| ID | Data | Configuration | Location in 004 | PoNs | Status |
|---|---|---|---|---|---|
| **U-A0** | UMI-consensus | Final configuration, no mask | `/final_config/per_sample/<run>/<s>.{cnr,cns,call.cns,genemetrics.tsv}` | `final_main__batch_000.pons.tar.gz` (`file-JF61YzQ4P0z8QFVv534xJvXx`), `phaseB/pons/` | Available |
| **U-M** | UMI-consensus | Arm M, **before** the PAR post-step | `/sim2c/real/sim2c_real__<run>.calls.tar.gz`, `calls/M/d00/<run>/` (+ `mask_qc.tsv`) | As U-A0 | Available |
| **N-A0** | Non-UMI | Final configuration, no mask, PAR post-step applied | `/final_config_nonumi/outputs/final_config_nonumi.D_gcclip_on_s2.parfix.outputs.tar.gz` (`file-JFP4vx84k6GK39PpXZ4b439y`), `nonumi/phaseN/calls/D_gcclip_on_s2/d00/<run>/<s>.{cnr,cns,call.cns,genemetrics.tsv,genemetrics.noparfix.tsv}`; sha256 list `file-JFP4vx84k6GPZ0gxJK4vp1yg` | `nonumi_main__batch_000.pons.tar.gz` (`file-JFP4Yvj4GJKXqfkxzk71FpYj`), `nonumi/phaseN/pons/`. The shared PoN md5 `b4a37372…` = the controlled file | Available |
| **N-M** | Non-UMI | Arm M | `/final_config_nonumi/dropout/parfix/dropout_nonumi_NM.calls.tar.gz` (`file-JFPX1JQ4763FK094J6vqj1yQ`), `calls/M/d00/<run>/<s>.{cnr,cns,call.cns,genemetrics.tsv,genemetrics.noparfix.tsv,mask_qc.tsv}`; sha256 list `file-JFPX1J847638Qy9X54JyJxyV`. `calls/A0/` in the same tarball = N-A0 (286/286 genemetrics identical) | As N-A0 | Available (DI-4139 comment 72496); release gate |

**Non-UMI inputs:**
- coverage (`eggd_cgp-cnvkit-coverage` 2.0.2): `/coverage_nonumi/<run>/<s>.targetcoverage.cnn`;
- BAMs and inventory: `/final_config_nonumi/repo_snapshot/e8bb628/inventory/nonumi_inventory.tsv`, with columns `old_cnn`, `new_cnn`, `new_bam` (226/226 `ok`);
- sex truth: `final_config_nonumi.tier1.tar.gz` (`file-JFP4Q704VG07vK8B14BJgbZK`), `tier1/sex_raw/<run>/sex.tsv` and `tier1/sex_threshold/`;
- PAR post-step log: `/final_config_nonumi/outputs/parfix_log.tsv` (`file-JFP4vx84k6G2Gx5J0QP2g8yP`).

**U-M has no PAR post-step** (checked 2026-10-26).
- The `genemetrics.tsv` in the `/sim2c/real/` tarballs is the **pre-post-step** read-out. On run 0734, the 12 treated-male sequencings whose post-step changes a row differ from `/final_config/` (U-A0, post-step applied); every other `genemetrics.tsv` is identical.
- So, against U-M:
  - the app's `genemetrics.noparfix.tsv` must equal U-M `genemetrics.tsv` (all 226);
  - the app's `genemetrics.tsv` must equal U-M `genemetrics.tsv` for sequencings not treated as male, and U-A0 `genemetrics.tsv` for treated males without a masked probe;
  - for treated males with a masked probe (e.g. 24260S0018) there is no UMI post-step baseline; the `compare` invariants apply, and the non-UMI N-M covers the same path.
- N-A0 and N-M both include the post-step (60 treated males, each with a `genemetrics.noparfix.tsv`).

**N-M as delivered** (DI-4139 comment 72496):
- the tarball has `calls/{A0,M}/d00/<run>/` for all 226;
- `genemetrics.tsv` is after the PAR post-step, and `.cnr` is unmasked;
- it was built with the validated non-UMI PoNs (`file-JFP4Yvj4GJKXqfkxzk71FpYj`);
- checks: A0 = DI-4118 in 904/904 raw files and 286/286 genemetrics files; the post-step has 0 check failures.

Fetched files are verified against the sha256 list.

**Roles of the baselines:**
- U-A0, U-M and N-A0 prove the code is correct on both data types.
- N-M proves the production configuration (non-UMI, mask on) and is the release gate (E3/E4).

**Expected dropout effect on non-UMI data** (for the L4 and T13/T15 expectations):
- calls change in 3 sequencings: 25330S0048 +17, 25330S0109 +5, 24260S0018 +3;
- the 25330S0044 FGFR3 loss call **stays**: probe 1 has 11.2 raw reads, so it is not masked;
- 4 probes are masked on UMI only (FGFR3 in 25330S0044 and 25330S0054, EGFR in 25330S0057, CCNE1 in 25330S0113), and none on non-UMI only.

**Coverage consistency.** The PoN and N-A0 were built from the coverage app's `.targetcoverage.cnn`, but the batch app computes its own coverage from the BAM. The two must be byte-identical; the batch app's BED dedup and `--processes` are the candidate differences. This is tested in L4 (one non-UMI BAM) and L5 (every BAM-start run).

---

## 1. New decisions (accepted 2026-10-25)

| ID | Decision | Accepted option |
|---|---|---|
| **D8** | How to reconcile `sex_check.py` with the app's sex sources | Port it to a per-sample `cnvkit_sex_qc.py`, in the same way as the other two scripts, with explicit `--lims-sex/--lims-source/--override` arguments.<br>• Keep every `sex_check.py` column name and value for the no-override case. Prove this with an equivalence test on all 226 sequencings.<br>• Fill `applied_sex`/`applied_source` from the app's resolution: `override`, `sample_name`, `bam_name` or `cnvkit_inference`. For U samples, `applied_sex` gives the inferred value, e.g. `cnvkit_inferred(Male)`.<br>• Append `par_fix_applied`, `par_span_segments` and `dropout_probes`.<br>• `sex_mismatch` compares the declared sex (override if set, else LIMS) with chrY. It is identical to the validation whenever there is no override.<br>• The MultiQC file keeps its columns and adds `dropout_probes`. |
| **D9** | Accepted override values | `male`, `female`, `M` or `F`, case-insensitive; empty means none. Any other value fails the job before any compute. |
| **D10** | `arm` column in `dropout_qc` | Keep it. The value is always `M`, and the file stays byte-identical to the validation `mask_qc.tsv`. Change the Jira and plan §4 column lists to match. |
| **D11** | Allowed values for `diploid_parx_genome` | Set `choices: ["grch38"]` in `dxapp.json`, optional, default empty (D1).<br>• The PAR post-step coordinates are GRCh38-only.<br>• When the input is empty, skip the post-step (as for the validation s0 arms) and log a warning. |

---

## 2. Target code layout

```
src/code.sh                                    dx I/O only: download, validate inputs, call run_pipeline, upload
resources/usr/local/lib/cgp_cnvkit_batch/
    pipeline.sh                                all analysis logic as bash functions (no dx calls)
resources/usr/local/bin/
    dropout_mask.py                            port of deadprobe_mask.py mask, arm M (DI-4139)
    cnvkit_parfix.py                           port of genemetrics_parfix.py: decide / shift / compare (DI-4101)
    cnvkit_sex_qc.py                           port of sex_check.py: resolve / qc (DI-4098, DI-4100, DI-4103; D8)
tests/
    conftest.py, pytest.ini                    markers: unit, pipeline, validation, docker
    fixtures/                                  small synthetic .cnr/.cnn/.cns/sex.tsv files (a few dozen rows each)
    unit/test_dropout_mask.py
    unit/test_cnvkit_parfix.py
    unit/test_cnvkit_sex_qc.py
    pipeline/stub_cnvkit.sh                    records argv, writes fixture outputs
    pipeline/test_command_lines.py             the argv contract per scenario
    pipeline/test_dxapp_contract.py            dxapp.json ↔ code.sh inputs/outputs consistency
    validation/fetch_validation_data.sh        dx download into $CNVKIT_VALIDATION_DATA (gitignored)
    validation/test_equiv_*.py                 bundled scripts vs the analysis-repo originals (226 sequencings)
    harness/run_local.sh                       real image, real data, pipeline.sh, cmp against expected
    dx/run_matrix.sh, dx/check_outputs.py      applet tests in the 003 project
    dx/run_equivalence.sh, dx/compare_equivalence.py   the 226-sequencing test
requirements-test.txt                          pytest only (the bundled scripts stay stdlib-only)
```

### How the files fit together
- On the worker, `resources/` is unpacked to `/`. `code.sh` sources `/usr/local/lib/cgp_cnvkit_batch/pipeline.sh`, and the tests source the same file from the repo.
- `pipeline.sh` doesn't call `docker` directly; it calls `run_cnvkit`.
  - `code.sh` defines `run_cnvkit` as the Docker call, with `-w` passed so the post-step can run in `parfix/`.
  - The pipeline tests define it as `stub_cnvkit.sh`.
  - The local harness defines it as a real Docker call.
- `pipeline.sh` takes its inputs from environment variables: `SAMPLE_ID`, `BAM_NAME`, `DIPLOID_PARX_GENOME`, `SEX_OVERRIDE`, `SAMPLE_SEX_DEPRECATED`, `PURITY`, `PLOIDY`, `DROP_LOW_COVERAGE`, `DROPOUT_MASK`, `NPROC` and `START_FROM`.
  - `START_FROM` is `bam`, or `coverage` for the harness, which then starts from a given `<s>.targetcoverage.cnn`.
  - It reads and writes files in the current directory.
- Python runs on the worker host (Ubuntu 24.04, Python 3.12, stdlib only), not in the image. Local tests use the same Python 3.12.

### Step order in `run_pipeline`
Steps 1–8 are fatal on error; steps 9–10 are non-fatal, as in v2.0.x.

1. `resolve` sex (override → sample name → BAM name → U). Fail on a bad override.
   - Log the warnings: empty `diploid_parx_genome`; `purity`/`ploidy`/`sample_sex` supplied but ignored.
2. `coverage` (skipped if `START_FROM=coverage`), then `fix`, with no `--no-edge`.
3. `dropout_mask.py mask M …`, or a header-only `dropout_qc` when `dropout_mask=false`. Log the masked count.
4. `segment` on `seg_input.cnr` or `.cnr`.
5. `cnvkit.py sex` on the `.targetcoverage.cnn` and on the `.cnr`.
6. `call [-x]`, then `genemetrics [-x]`.
7. Post-step decision, then `shift`, `genemetrics` in `parfix/`, `compare` (which fails the job on any non-PAR difference), then the renames. Write `noparfix` as a copy when no post-step runs.
8. `cnvkit_sex_qc.py qc`, which writes `sex_qc.tsv` and `sex_check_mqc.tsv`.
9. `scatter` (unchanged).
10. `diagram [-x] [--diploid-parx-genome]`.

Every decision is written to the log with a fixed prefix (`[sex]`, `[dropout]`, `[parfix]`, `[deprecated]`) so the DNAnexus tests can grep for it.

---

## 3. Test layers

| Layer | What it proves | Needs | Runs |
|---|---|---|---|
| **L0 static** | `shellcheck` clean; `dxapp.json` is valid JSON; every input is used in `code.sh`; every output is uploaded; patterns match the file names | – | Every commit, CI |
| **L1 unit** | Each Python script on synthetic fixtures: rule boundaries, edge cases, error paths | pytest | Every commit, CI |
| **L2 pipeline (argv contract)** | The exact CNVkit command lines per scenario, step order, which `.cnr` each step reads, the post-step decision, fail paths. No Docker | pytest + bash | Every commit, CI |
| **L3 script equivalence** | Bundled ports byte-identical to `deadprobe_mask.py` / `genemetrics_parfix.py` / `sex_check.py` on the 226 sequencings, on **both** the UMI and the non-UMI data | Validation data downloaded locally | Before PR; on changes to `resources/` |
| **L4 local harness** | Real image + `pipeline.sh` reproduce U-A0, U-M, N-A0 and N-M byte for byte on 16 sequencings; batch-app coverage = coverage-app coverage | Docker, image, data | Before PR |
| **L5 applet matrix** | dx I/O, inputs, outputs, job logs, fail paths on DNAnexus, from non-UMI BAMs | 003 project | After build |
| **L6 full equivalence** | 226 sequencings against N-A0, U-M and N-M (the release gate); plus the production read-out with the controlled PoN | Validation project | Before release; Confluence evidence |

### L1 unit cases (minimum)

**`dropout_mask.py`:**
- Threshold boundaries:
  - log2 exactly −3 is not masked;
  - a bin with exactly 10 raw reads is not masked;
  - a bin with 9.99 raw reads is masked.
- Runs of deep probes:
  - runs of 1 and 2 are masked; a run of 3 is kept;
  - a run broken by a chromosome change counts as two runs.
- chrX and chrY are never masked.
- A bin absent from the `.cnn` counts as 0 reads.
- Only `weight` changes; every other byte is identical.
- The QC file is header-only when nothing is masked; the `qc-header` output is identical to it.
- An arm other than `M` is refused.

**`cnvkit_parfix.py`:**
- `par_state` at each PAR edge, using the 0-based half-open ends.
- A bin spanning a PAR boundary is left unchanged and counted.
- The `repr` formatting of the shifted log2.
- The decision table:
  - `-x male` → yes;
  - `-x female` → no;
  - no `-x` with `.cnr` sex Male → yes;
  - no `-x` with `.cnr` sex Female → no;
  - empty PARX → no.
- `compare`: passes on PAR-only changes; fails on a non-PAR change; accepts the 5-column empty genemetrics header.

**`cnvkit_sex_qc.py`:**
- LIMS letter parsing, including names with fewer than 5 fields, lower-case letters and extra fields.
- BAM-stem fallback only when the sample name has fewer than 5 fields.
- The override values (D9), including a bad value failing.
- The `sex_mismatch` truth table.
- chrY present at the −5 boundary and for NaN.
- MultiQC header lines unchanged.

### L2 pipeline scenarios
Each scenario asserts the full argv list of every `run_cnvkit` call.

1. **v2.0.2 characterisation** (written first, against the current code; see Phase 1).
2. LIMS M, F and U, each with PARX `grch38`.
3. An override that agrees with LIMS, and one that disagrees.
4. A bad override value: exits non-zero before `coverage`.
5. `dropout_mask` set to true and to false: `segment` input and `dropout_qc` content.
6. PARX empty: no `--diploid-parx-genome` anywhere, no post-step, a warning logged.
7. `purity`/`ploidy`/`sample_sex` supplied: no `--purity`, `--ploidy` or `--sample-sex` from them, and a warning logged.
8. The PAR post-step genemetrics runs with cwd `parfix/` and the original arguments.
9. The non-PAR assertion fails: the job fails.
10. Scatter or diagram fails: the job still succeeds.

### Validation data (L3, L4, L6)
`tests/validation/fetch_validation_data.sh [umi|nonumi|all]` downloads to `$CNVKIT_VALIDATION_DATA/{umi,nonumi}/` from `project-JF3p1f04Yggp0vf3PYjy7Q0F`. Each item is verified against its manifest md5 or the sha256 list.
- **UMI:**
  - `/coverage/<run>/` coverage;
  - U-A0 (`11_final_per_sample_files.csv`);
  - U-M (`/sim2c/real/` tarballs);
  - `final_config.sex_qc.tsv` and `parfix_log.tsv`;
  - the `final_main` PoN tarball.
- **Non-UMI:**
  - `/coverage_nonumi/<run>/` coverage;
  - N-A0 (the parfix outputs tarball);
  - the tier1 sex tables and `parfix_log.tsv`;
  - the `nonumi_main` PoN tarball;
  - the controlled PoN `file-JFP3fzQ455KjXVZ23JgX1gF2`;
  - N-M (`file-JFPX1JQ4763FK094J6vqj1yQ`, sha256 list `file-JFPX1J847638Qy9X54JyJxyV`).

In total about 1.2 GB. The data is gitignored, never committed. Tests skip with a clear message when the data is absent.

---

## 4. Phases

Each phase ends green on L0–L2. Commits are one per subtask, with the Jira key in the message.

### Phase 0: preparation (no app code) — done 2026-10-25
1. Done: branch `dev_2.1.0` from `main`; `.gitignore` updated; `memory/version_state.md` corrected.
2. Done: D8–D11 recorded on DI-4094 (comment 72472) and DI-4139 (comment 72473); plan §7 and the `dropout_qc` columns updated.
3. Done: IDs resolved (§0, Phase 0 results).
4. Done: 16 local sequencings chosen (`tests/harness/local_samples.tsv`). Members use their leave-out PoN.

### Phase 1: test scaffolding and a record of current behaviour
1. `requirements-test.txt` and `pytest.ini` (markers). Set up with `python3 -m venv .venv && .venv/bin/uv pip install -r requirements-test.txt`.
2. L0: `shellcheck src/code.sh`, and `test_dxapp_contract.py`.
3. L2 scenario 1: record the **current** v2.0.2 argv for one sample. To do this without refactoring, run `code.sh`'s `main` with `dx`, `dx-jobutil-add-output` and `docker` replaced by PATH stubs.
4. Commit: "test: add test scaffolding and v2.0.2 command-line characterisation".

### Phase 2: refactor without behaviour change
1. Move the logic into `pipeline.sh` (`run_pipeline`, using `run_cnvkit`). Keep only dx I/O in `code.sh`.
2. Scenario 1 must still pass unchanged. That shows the refactor didn't change behaviour.
3. Commit: "refactor: split dx I/O from analysis logic (no functional change)".

### Phase 3: resource scripts (test-driven, ports first)
For each script: write the L1 tests, port the code with a header (source path, commit, md5, validation page), make it pass, then run L3 equivalence.

1. **`dropout_mask.py`** (DI-4139). Copy `mask` and the constants verbatim, and drop the `blacklist`, `L5`/`L3` and `MB` paths. The positional CLI stays identical.
   - Add a `qc-header <out>` subcommand for `dropout_mask=false`, so the header text exists in one place only.
   - L3: 226/226 byte-identical (masked `.cnr` and QC) on the UMI data; totals 82 probes in 42 sequencings, plus 25330S0104.
   - L3 also runs on the 226 non-UMI pairs against `deadprobe_mask.py` itself. That proves the port on production data. Expected totals: 78 probes / 41 sequencings, plus 25330S0104; QC rows equal the N-M `mask_qc.tsv`.
2. **`cnvkit_parfix.py`** (DI-4101). Subcommands:
   - `decide --applied-sex … --cnr-sex sex_cnr.tsv --parx …` prints `1` or `0` with a reason;
   - `shift <cnr> <call.cns> <outdir>` prints the counts;
   - `compare <orig.gm> <new.gm> <call.cns>` exits 1 on a non-PAR difference and prints `par_check`.
   - L3: on the validation males, the shifted inputs and the final `genemetrics.tsv` are byte-identical to the analysis outputs, on UMI (`/final_config/outputs/`) and non-UMI (N-A0) data. Re-running on unmodified inputs reproduces the original (the `--rerun-check` equivalent).
3. **`cnvkit_sex_qc.py`** (DI-4098, DI-4100, DI-4103; D8). Subcommands:
   - `resolve --sample-id --bam-name --override` prints `applied_sex`, the source and the LIMS fields;
   - `qc …` writes `sex_qc.tsv` and `sex_check_mqc.tsv`.
   - L3: with no override, the columns shared with `sex_check.py` are byte-identical on all 226. The chrY fields equal `final_config.sex_qc.tsv` (UMI) and the tier1 `sex_raw` tables (non-UMI). `sex_mismatch=1` in 0 of the 173 LIMS M/F sequencings on both data sets.

Commits: one per script, with DI-4139, DI-4101 and DI-4103 respectively.

### Phase 4: app changes in `pipeline.sh` / `code.sh` / `dxapp.json`
In order. Each change updates the L2 scenarios first, so the test fails, and then the code.

| Step | Subtask | Change |
|---|---|---|
| 4.1 | DI-4096 | Remove `--no-edge` |
| 4.2 | DI-4097 (D1, D11) | Add the `diploid_parx_genome` input; pass it to fix/segment/call/genemetrics/diagram/sex/post-step; log a warning when empty |
| 4.3 | DI-4121 (D4) | Delete `CALL_ARGS`; `purity`/`ploidy` kept in `dxapp.json`, labelled deprecated, ignored with a warning |
| 4.4 | DI-4098, DI-4099, DI-4100 (D2, D9) | Add `sex_override`; `sample_sex` is deprecated, ignored with a warning; `-x` on call/genemetrics/diagram; source logged |
| 4.5 | DI-4102 | Both `cnvkit.py sex` runs |
| 4.6 | DI-4139 (D7) | Add the `dropout_mask` input; the mask step; `segment` reads `seg_input.cnr`; only `segment` sees the masked file |
| 4.7 | DI-4101 (D3) | The post-step; `noparfix` always written |
| 4.8 | DI-4103 | `cnvkit_sex_qc.py qc` |
| 4.9 | DI-4104 (D5) | New outputs `genemetrics_noparfix`, `sex_qc`, `sex_check_mqc`, `dropout_qc`, `target_coverage`; update the help text of every input and output; version `2.1.0` |

### Phase 5: local harness (L4)
- `tests/harness/run_local.sh <sample> --data umi|nonumi [--start-from coverage|bam] [--dropout-mask true|false]`:
  - sources `pipeline.sh` with the real Docker `run_cnvkit` and the row's `cn_reference` from the matching PoN tarball;
  - writes to `tests/.data/harness/<data>/<sample>/`, `cmp`s against the baseline and prints a per-file PASS/FAIL table.
- **Runs:**

  | Data | Mask | Sequencings | Start | Expected |
  |---|---|---|---|---|
  | UMI | on | 16 | coverage | U-M (= U-A0 where nothing is masked) |
  | UMI | off | the 4 masked | coverage | U-A0 |
  | Non-UMI | off | 16 | coverage | N-A0 |
  | Non-UMI | off | 1 (25328S0023, male) | **BAM** | N-A0, and `.targetcoverage.cnn` byte-identical to `/coverage_nonumi/` (covers `coverage`, the BED dedup and coverage consistency, §0b) |
  | Non-UMI | on | 16 | coverage | N-M; 25330S0044 keeps its FGFR3 call (not masked on non-UMI) |

- Expected in every run:
  - `.cnr`, `.cns`, `.call.cns`, `genemetrics.tsv` and `genemetrics.noparfix.tsv` byte-identical;
  - `dropout_qc` equal to `mask_qc.tsv` (mask on);
  - the post-step applied to the same samples as `parfix_log.tsv`;
  - the `sex_qc` chrY fields equal to the sex truth of that data set.
- **If the BAM-start coverage differs from the coverage app:** stop. The batch app must reproduce the coverage the PoN was built from, so find and fix the cause before Phase 7.

### Phase 6: documentation
README:
- the commands as run;
- sex sources and their precedence;
- the PAR post-step and the CRLF2 limitation;
- the dropout mask: the rule, that it only affects `segment`, its limitation (real 1–2-probe events meeting the rule are masked, listed in `dropout_qc`), and the `dropout_qc` output;
- the FGFR3 probe 1 limitation (DI-4139 comment 72496): a capture dropout that the mask does not remove at non-UMI depth (10 of 11 specimens have 11–30 raw reads), so FGFR3 loss calls from short segments need review;
- the sex QC outputs;
- no purity;
- `grch38` required for Atlas;
- the deprecated inputs;
- how to run the tests.

Also:
- Changelog `[2.1.0]` (Changed / Added / Removed / Deprecated, with Jira keys).
- Script headers.

Then run the pr-workflow review-fix loop locally.

### Phase 7: DNAnexus testing
1. **L5 matrix.** `dx build` into a 003 project (`003_<YYMMDD>_eggd_cgp-cnvkit-batch_v2.1.0`). `tests/dx/run_matrix.sh` launches the runs below, and `check_outputs.py` downloads the outputs, compares them with the expected files and greps the job logs.
   - Inputs are **non-UMI BAMs** (`nonumi_inventory.tsv` `new_bam`) with the validated non-UMI PoN of each sequencing, unless the row says otherwise.
   - Every run also checks `target_coverage` against `/coverage_nonumi/` (UMI runs: `/coverage/`).

   | # | Case | Expect |
   |---|---|---|
   | T1–T4 | LIMS M / F / U / XXY, `dropout_mask=false` | N-A0 outputs; `-x` and the source in the log |
   | T5 | Override M on a LIMS-U sample | `-x male`, source `override` |
   | T6 | Override F on a LIMS-M sample | `sex_mismatch=1` (D8), post-step not run |
   | T7 | Override `X` | Job fails before `coverage` |
   | T8 | `sample_id` without field 5, BAM name with it | Source `bam_name` |
   | T9 | Neither has field 5 | U, inference, post-step by `.cnr` sex |
   | T10 | `diploid_parx_genome` empty | Warning, no PARX flags, no post-step |
   | T11 | `dropout_mask=false` on 25330S0048 | N-A0 outputs, header-only `dropout_qc` |
   | T12 | `purity`/`ploidy`/`sample_sex` supplied | Warnings; outputs identical to T1 |
   | T13 | **UMI BAMs** and UMI PoNs: 25330S0048, 24260S0018, 25330S0044, 25330S0104, mask on | U-M outputs and `dropout_qc` rows |
   | T14 | Female-only | `noparfix` identical to `genemetrics` |
   | T15 | Production setting: the controlled PoN `file-JFP3fzQ455KjXVZ23JgX1gF2`, mask on, 25330S0048 (a PoN non-member) | N-M outputs (+17 genes vs N-A0) and the N-M `mask_qc.tsv` rows |

2. **L6 full equivalence** (the DI-4094 comment; part of DI-4118). `run_equivalence.sh` runs all 226 from BAM in `project-JF3p1f04Yggp0vf3PYjy7Q0F` with the final inputs (`diploid_parx_genome=grch38`, no purity/ploidy, `drop_low_coverage=true`). `compare_equivalence.py` writes a report TSV and Markdown.

   | Run | Data | Mask | PoN | Expected | When |
   |---|---|---|---|---|---|
   | E1 | Non-UMI | off | Validated non-UMI (leave-out / shared) | **N-A0**, byte-identical, all 226 | Now |
   | E2 | UMI | on | Validated UMI | **U-M**, byte-identical, all 226, comparing post-step files as in §0b | Now |
   | E3 | Non-UMI | on | Validated non-UMI | **N-M**, byte-identical, all 226 | Now: **release gate** |
   | E4 | Non-UMI | on | Controlled PoN | Production read-out. Non-members equal E3; differences for PoN members are listed, not failed | Now |

   - Acceptance as in `PLAN_v2.1.0.md` §5.3.
   - All four runs can go as soon as the applet passes L5.
3. The L5 and L6 reports become the Confluence test evidence: T1… in the `eggd_cgp-cnvkit-pon v2.0.0` format, as in `memory/confluence_doc.md`. Cite `atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn` (Confluence `4854120471`) as the PoN.

### Phase 8: release
**Gates:**
- E3 and E4 pass;
- the PoN page is approved.

The non-UMI dropout analysis left the rule and constants unchanged (DI-4139 comment 72496), so no change to `dropout_mask.py`, D6 or D7 is needed.

PR `dev_2.1.0` → `dev` → `main`, then `dx build --app`, `dx publish`, a GitHub release tag `v2.1.0`, and updates to Confluence and memory (`PLAN_v2.1.0.md` §6).

**Workflow handover:** the `eggd_atlas_cnv_workflow` `cn_reference` points at the **001** copy of the PoN after its rename and move. Check its file ID and md5 `b4a373729d9282ce0a0ed80fd29f6d0a` then. A `dx cp` keeps the file ID, and a rename doesn't change it.

---

## 5. CI (proposed)
Add `.github/workflows/tests.yml`. On every push and PR it runs L0–L2 on Ubuntu 24.04 with Python 3.12: `shellcheck` plus `pytest -m "unit or pipeline"`. It needs no Docker, data or dx.

L3–L6 stay manual, because they need DNAnexus data or access. Their reports are attached to the PR and to Confluence.

## 6. Risks and mitigations

| Risk | Mitigation |
|---|---|
| The refactor silently changes v2.0.2 behaviour | Phase 1 records the current argv before any change; Phase 2 must keep it |
| A port drifts from the validated script | L3 byte-identity on all 226 sequencings; source md5 recorded in each header |
| The post-step via the CLI differs from the validation's in-process call | L3 re-runs both and `cmp`s; the harness re-checks with the real image |
| PoN-member differences hide real regressions | L6 run 1 gives each sequencing its validated leave-out or shared PoN, so byte-identity holds for all 226; run 2 (static PoN) shows the production behaviour |
| The workflow still wires PURPLE `sample_sex`/`purity`/`ploidy` | They are ignored with a warning (D2, D4); T12 proves the outputs are unaffected |
| `dropout_mask=false` drifts from A0 | T11, the Phase 5 A0 reruns and E1 |
| The batch app's own coverage differs from the coverage app's, which the PoN and N-A0 were built from | L4 BAM-start `cmp`; `target_coverage` checked in every L5/L6 run; Phase 5 stops on a difference |
| The mask rule was tuned on UMI data; depth rose about 1.43× on non-UMI | **Closed:** the non-UMI evaluation (DI-4139 comment 72496) shows the rule and constants are unchanged and arm M meets every T3 margin. L3 and E3 test the port on non-UMI data |
| A UMI-only expectation is applied to non-UMI data (e.g. the 25330S0044 FGFR3 call removal) | §0b lists the non-UMI dropout effect; every L4/L5/L6 expectation names its baseline (U-* or N-*) |
| The production PoN file is renamed or moved after testing | Tests and evidence cite the file ID and md5; the workflow handover re-checks both on the 001 copy |
