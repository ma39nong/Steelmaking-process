"""Write figure captions and a Chinese chapter copy using English figure assets."""
from pathlib import Path
import json
import re
import pandas as pd

P = Path(__file__).resolve().parents[1]
OUT = P / 'figures/english'
REPORTS = P / 'reports'
HEATS = pd.read_csv(P / 'results/final_summary.csv').heat_id.tolist()


def entries():
    rows = [
        (1, 'main/singular_value_energy', 'Cumulative singular-value energy of the ten raw spectral matrices.', 'Each curve uses every eighth recorded frame. The first four components account for more than 99.95% of the squared singular-value energy in every heat, supporting a rank-4 background.', 'v3/results/singular_energy.csv'),
        (2, 'main/method_workflow', 'Workflow for adaptive anomaly identification, robust point repair, and continuous spectral reconstruction.', 'The stages connect actual-coordinate quality analysis, locally calibrated combined evidence, rank-4 Huber completion, and bilinear reconstruction. Recorded top codes and conditional estimates are stored separately.', 'final/src/detectors.py; v3/src/model3.py; v3/src/stable_completion.py; v3/src/interface.py'),
        (3, 'main/interpolation_comparison', 'Held-out reconstruction errors under temporal and wavelength downsampling.', 'All methods use the same retained coordinates. Bilinear interpolation has the lowest development aggregate NRMSE and the lowest test RMSE at 1.5 s, 2.0 s, and both wavelength intervals; PCHIP is slightly more accurate at 1.0 s. Cubic spline and zero-smoothing tensor B-spline agree numerically in this configuration.', 'v2/results/interpolation_test.csv; v2/results/interpolation_independent_macro.csv'),
        (4, 'main/end_to_end_comparison', 'Anomaly detection and end-to-end repair under the paired contamination-confirmation protocol.', 'Detection scores use pooled TP, FP, and FN. RMSE is averaged within each heat and then across heats. The selected adaptive score achieves precision 89.31%, recall 70.52%, and F1 78.81%; the comparison uses the same 300-step solver budget for all strategies.', 'final/results/unified_comparison.csv: confirmation'),
        (5, 'main/solver_convergence', 'Huber objectives and relative iterate changes for four representative full-grid fits.', 'Rows correspond to Heats 6C02057, 6C02058, 6C02107, and 6C02108. Objectives and iterate changes come from the saved 1,200-step histories; the dashed line marks the fixed relative-change tolerance of 1e-4.', 'final/results/6C02057_solver.json; 6C02058_solver.json; 6C02107_solver.json; 6C02108_solver.json'),
        (6, 'main/repair_counts', 'Numbers and percentages of actual point repairs across the ten heats.', 'There are 1,043 changed observations out of 27,816,683 records (0.00375%). The plots show actual exported changes, rather than candidate counts or injected anomalies.', 'final/results/final_summary.csv'),
        (7, 'main/four_reviewed_spikes', 'Measured and repaired local spectra at four reviewed statistical spikes.', 'The four locations are 6B03205 at 14.5 s and 512.659 nm; 6B03270 at 412.5 s and 410.105 nm; 6B03271 at 208.5 s and 347.691 nm; and 6C02059 at 757.0 s and 575.499 nm. The dashed curves are model-based repair estimates.', 'v3_1/results/high_confidence_statistical_points.csv; final/results/four_reviewed_actual_repairs.csv; final/data/*_interface.npz'),
        (8, 'appendix/6C02058_repair_detail', 'Local spectral and temporal evidence for a representative actual repair in Heat 6C02058.', 'The top row shows a full spectrum, its local wavelength neighborhood, and a time trace. The lower row compares the raw and repaired neighborhoods using shared intensity limits, followed by their signed difference. The example is the actual repair with the largest absolute change in this heat.', 'final/results/actual_repairs.csv; final/data/6C02058_interface.npz'),
        (9, 'main/all_heats_continuous_surfaces', 'Recorded-scale reconstructed spectral surfaces for all ten heats.', 'Each panel uses the actual wavelength and time coordinates of one heat and the same intensity color scale from 0 to 65,535 counts. The plotting grid is thinned for readability; the continuous query interface retains the complete repaired node matrix.', 'final/data/*_interface.npz; v3/src/interface.py'),
    ]
    for i, heat in enumerate(HEATS, 10):
        rows.append((i, f'appendix/{heat}_continuous_surface', f'Recorded-scale continuous spectral surface for Heat {heat}.', 'The node values are the frozen final repair result. Axes show wavelength (nm), relative time (s), and recorded intensity (counts); the bilinear function is defined on the actual observation domain.', f'final/data/{heat}_interface.npz; v3/src/interface.py'))
    return rows


def main():
    manifest = pd.read_csv(OUT / 'figure_manifest.csv')
    audit = json.loads((OUT / 'language_check.json').read_text())
    assert audit['unique_figures'] == 49 and audit['files'] == 147
    rows = entries()
    title = '# Problem 1: English Figure Report\n\n'
    title += 'This figure report accompanies Section 5 of the manuscript. Figure-internal titles, axes, legends, annotations, and color-bar labels are in English. The Chinese manuscript and its original figures are preserved. All figures are rendered from frozen arrays, saved statistics, or saved solver histories; no detector, model fit, parameter search, or experiment is rerun.\n\n'
    title += '## Files and editing\n\n'
    title += '- `problem1/final/figures/english/main/`: main comparison, workflow, validation, and overview figures.\n- `problem1/final/figures/english/appendix/`: per-heat local repairs, full spectra and signed changes, surfaces, and rapid-transition diagnostics.\n- PNG: 300 dpi preview and raster publication output. SVG: editable text and vector elements. PDF: publication export with embedded TrueType fonts. Heat maps and surfaces may contain rasterized graphical layers; their text remains editable in SVG.\n- [Chinese chapter with English figures](问题一_中文论文_V2_英文配图版.md): the original Chinese prose, equations, and tables are retained; only figure links and captions are replaced.\n- [Asset manifest](../figures/english/figure_manifest.csv) and [language check](../figures/english/language_check.json): 49 unique figures, three formats each, no CJK characters in figure text or SVG content.\n\n'
    title += '## Main-text figure map\n\n|Figure|Subject|PNG|SVG|PDF|\n|---|---|---|---|---|\n'
    for number, stem, caption, analysis, source in rows:
        title += f'|5-{number}|{caption}|[PNG](../figures/english/{stem}.png)|[SVG](../figures/english/{stem}.svg)|[PDF](../figures/english/{stem}.pdf)|\n'
    title += '\n## Figures, captions, and interpretation\n\n'
    for number, stem, caption, analysis, source in rows:
        title += f'### Figure 5-{number}\n\n![Figure 5-{number}](../figures/english/{stem}.png)\n\n**Caption.** {caption}\n\n**Interpretation.** {analysis}\n\n**Source (relative to `problem1/`).** `{source}`\n\n'
    title += '## Supplementary per-heat evidence\n\nThe following panels use the same frozen data. The largest actual repair in each heat determines its local example. Rapid-transition spectra are selected using the saved frame diagnostics; no whole-frame deletion is introduced.\n\n|Heat|Local repair|Full spectra and differences|Rapid transition|\n|---|---|---|---|\n'
    for heat in HEATS:
        title += f'|{heat}|[PNG](../figures/english/appendix/{heat}_repair_detail.png) / [SVG](../figures/english/appendix/{heat}_repair_detail.svg) / [PDF](../figures/english/appendix/{heat}_repair_detail.pdf)|[PNG](../figures/english/appendix/{heat}_full_spectrum_difference.png) / [SVG](../figures/english/appendix/{heat}_full_spectrum_difference.svg) / [PDF](../figures/english/appendix/{heat}_full_spectrum_difference.pdf)|[PNG](../figures/english/appendix/{heat}_rapid_transition.png) / [SVG](../figures/english/appendix/{heat}_rapid_transition.svg) / [PDF](../figures/english/appendix/{heat}_rapid_transition.pdf)|\n'
    title += '\n### Wavelength-region comparison\n\n![Performance across wavelength regions](../figures/english/main/wavelength_region_comparison.png)\n\n**Caption.** Repair RMSE and mean case-wise F1 in blue, middle, and red wavelength regions under the common contamination-confirmation protocol.\n\n**Interpretation.** Regional values come from the saved per-heat, per-region results. Case-wise F1 averages in this panel have a different aggregation definition from pooled F1 in Figure 5-4.\n\n**Source.** `final/results/confirmation_per_heat_band.csv`. [SVG](../figures/english/main/wavelength_region_comparison.svg) / [PDF](../figures/english/main/wavelength_region_comparison.pdf).\n\n'
    title += '## Reproduction and checks\n\nFrom the repository root, using the existing environment and frozen local data:\n\n```bash\nOPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/final/src/plots_english.py\n.venv/bin/python problem1/final/src/english_figure_report.py\n```\n\nThe renderer checks every visible Matplotlib text item and every SVG for CJK characters. The original figures, numerical arrays, parameters, and experiment tables are preserved. Spatial/temporal repair panels retain their original example-selection rule. The English convergence figure includes the saved objectives as well as relative changes, and the interpolation comparison separates time and wavelength panels for legibility.\n'
    (REPORTS / 'Problem1_English_Figure_Report.md').write_text(title)
    chapter = (REPORTS / '问题一_中文论文_V2_正式章节版.md').read_text()
    # Change only image lines and their stand-alone captions; leave prose/math/tables intact.
    for number, stem, caption, analysis, source in rows:
        pattern = rf'^!\[图5-{number} [^\]]*\]\([^\n]*\)$'
        replacement = f'![Figure 5-{number}. {caption}](../figures/english/{stem}.png)'
        chapter, n = re.subn(pattern, replacement, chapter, flags=re.M)
        assert n == 1, (number, 'image replacement', n)
        chapter, n = re.subn(rf'^图5-{number}　[^\n]*$', f'Figure 5-{number}. {caption}', chapter, flags=re.M)
        assert n == 1, (number, 'caption replacement', n)
    (REPORTS / '问题一_中文论文_V2_英文配图版.md').write_text(chapter)
    metadata = []
    for number, stem, caption, analysis, source in rows:
        metadata.append(dict(manuscript_figure=f'5-{number}',english_stem=stem,caption=caption,source=source))
    pd.DataFrame(metadata).to_csv(OUT / 'manuscript_figure_map.csv', index=False)
    print(f'Wrote English figure report and chapter copy; {len(rows)} main-text figures, {audit["unique_figures"]} unique figure sets.')


if __name__ == '__main__':
    main()
