# Changelog

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
