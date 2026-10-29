# Changelog

## [2.1.0] - unreleased

Runs the validated Atlas CGP configuration (`D_gcclip_on_s2` with the dropout mask) on production (non-UMI) data
(DI-4094).

### Changed
- `fix`: edge correction on (`--no-edge` removed), matching the PoN build (DI-4096)
- `call`: no `--purity`/`--ploidy`; log2 ratio thresholds only, as validated (DI-4121)
- `call`, `genemetrics` and `diagram` get `-x` from the applied sex: `sex_override`, else the LIMS letter (field 5 of
  `sample_id`, else of the BAM name), else none (DI-4098, DI-4100)
- `segment` reads the dropout-masked copy of the `.cnr`; every other step reads the unmasked file (DI-4139)
- `genemetrics` output: after the PAR post-step for samples treated as male (DI-4101)
- Analysis moved to `resources/usr/local/lib/cgp_cnvkit_batch/pipeline.sh`; `src/code.sh` only stages inputs and
  uploads outputs. A failed upload now stops the job

### Added
- `diploid_parx_genome` input (`grch38` or empty): `--diploid-parx-genome` on fix, segment, call, genemetrics, diagram
  and the sex checks; warning when empty (DI-4097)
- `sex_override` input (`male`/`female`/`M`/`F`; any other value fails the job) (DI-4099)
- `dropout_mask` input (default `true`; `false` = the validated configuration without the mask, arm A0) and `dropout_mask.py`, a port of `deadprobe_mask.py mask M` @ `ddad79b`
  (DI-4139)
- `cnvkit.py sex` on the raw coverage and on the `.cnr` (DI-4102)
- PAR post-step, `cnvkit_parfix.py` (port of `genemetrics_parfix.py` @ `f933f7a`) (DI-4101)
- Sex QC, `cnvkit_sex_qc.py` (port of `sex_check.py` @ `31a62fc`) (DI-4103)
- Outputs `genemetrics_noparfix`, `sex_qc`, `sex_check_mqc`, `dropout_qc`, `target_coverage` (DI-4104)
- Tests: static checks, unit tests, command-line contract tests (CI), equivalence tests against the validation on
  226 sequencings and a local harness with the real image

### Deprecated
- `purity`, `ploidy`: ignored with a warning; the `ploidy` default of 2 is removed (DI-4121)
- `sample_sex`: ignored with a warning; it was wired from PURPLE's inferred sex (DI-4099)

### Removed
- The purity > 0.40 routing to `call --purity --ploidy` (DI-4121)

## [2.0.2] - 2026-07-17

### Changed
- `ploidy` input: `class` changed from `int` to `float`; label updated to "Tumour ploidy"; now accepts PURPLE's raw float output (e.g. 3.1) directly without manual rounding

### Removed
- `target_avg_size` input: was defined in `dxapp.json` but never used in `src/code.sh`; amplicon mode uses the panel BED directly and does not require bin-size specification

## [2.0.1] - 2026-06-16

### Fixed
- Address review comments (nadia31415)
- Correct genemetrics output file extension in README: `.csv` → `.tsv`
- Correct documented log₂R call thresholds in README to match CNVkit `fc65941d` defaults (amp ≥0.7, gain ≥0.2, loss <−0.25, del <−1.1); previous values were incorrect

## [2.0.0] - 2026-06-16

### Changed
- Replace per-job venv/pip + BiocManager install with pre-built Docker image (`cgp-cnvkit:1.0.0`, ubuntu:24.04 base); eliminates ~5 min install overhead and removes internet dependency at runtime
- CNVkit updated to git SHA `fc65941d` (post-0.9.13): pomegranate dependency removed, HMM rewritten in pure numpy/scipy, Python 3.12/3.14 CI-tested
- R / DNAcopy now installed as pre-compiled Ubuntu 24.04 debs inside Docker image

### Added
- BED interval deduplication before coverage step (`awk !seen`); guards against CNVkit ≥0.9.13 raising on duplicate genomic coordinates in panel BEDs
- `authorizedUsers` and `developers` fields in `dxapp.json` (org-emee_1)
- `allProjects` VIEW access in `dxapp.json`

### Removed
- `networkAccess` field from `dxapp.json` (no internet required at runtime)

## [1.0.0] - 2026-05-27

### Added
- Initial release
- Full CNVkit per-sample pipeline: coverage → fix → segment → call → scatter → diagram → genemetrics
- Purity/ploidy routing: PURPLE-corrected integer calls when purity > 0.40
- Scatter plot: y-axis ±3, fig-size 14×4 for readability alongside PURPLE plots
- Diagram: sex-corrected chromosome ideogram when sample_sex supplied
- Virtual environment install to avoid Ubuntu 24.04 system package conflicts
- Genemetrics filtered to |log2R| > 0.3, ≥3 probes
