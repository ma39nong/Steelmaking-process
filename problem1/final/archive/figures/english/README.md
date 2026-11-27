# English publication figures for Problem 1

This directory contains **49 figure sets**, each exported as a 300 dpi PNG, an SVG with editable text, and a PDF with embedded fonts. Figure-internal titles, axes, legends, annotations, and color-bar labels use English. Original Chinese assets remain in their original directories.

- [English figure report: captions, interpretation, and sources](../../reports/Problem1_English_Figure_Report.md)
- [Chinese manuscript with English figures and captions](../../reports/问题一_中文论文_V2_英文配图版.md)
- [Main-text figure mapping](manuscript_figure_map.csv)
- [Asset hashes](figure_manifest.csv)
- [Language audit](language_check.json)

`main/` contains nine overview, workflow, comparison, and validation figures. `appendix/` contains four figures for each heat: local repair details, full spectral matrices and differences, the final surface, and rapid-transition diagnostics. Nineteen unique figure sets are linked from the main manuscript; the remaining sets support the appendix.

SVG is the preferred editable source for labels and layout. Heat maps and surfaces may be rasterized within the SVG or PDF; text and ordinary curves remain separate elements. Numerical changes require returning to the frozen data, rather than editing the displayed curves.

The plotting script reads existing NPZ/CSV/JSON artifacts. It does not run anomaly detection, model fitting, parameter searches, or experiments. All original results and released data remain unchanged. The English convergence panels additionally display saved objective values, and the interpolation comparison separates time and wavelength groups for readability.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/final/src/plots_english.py
.venv/bin/python problem1/final/src/english_figure_report.py
```

Run from the repository root using its existing environment and frozen local data. Every visible Matplotlib text item and all SVG content are checked for CJK characters before export.

## Completed publication checks

All 147 asset hashes were verified against the manifest. All 49 SVGs retain text elements, all 49 PDFs yield readable text without CJK characters, and every PNG reports 300 dpi. Representative layouts and the main-figure contact sheet were visually inspected. The new chapter differs only in its 19 image links and 19 captions; Chinese prose, equations, and tables match the original. Twenty frozen numerical artifacts (ten interfaces and ten difference matrices) still match the existing artifact manifest.
