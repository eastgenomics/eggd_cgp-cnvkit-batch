# Changelog

## [1.0.0] - 2026-05-27

### Added
- Initial release
- Full CNVkit per-sample pipeline: coverage → fix → segment → call → scatter → diagram → genemetrics
- Purity/ploidy routing: PURPLE-corrected integer calls when purity > 0.40
- Scatter plot: y-axis ±3, fig-size 14×4 for readability alongside PURPLE plots
- Diagram: sex-corrected chromosome ideogram when sample_sex supplied
- Virtual environment install to avoid Ubuntu 24.04 system package conflicts
- Genemetrics filtered to |log2R| > 0.3, ≥3 probes
