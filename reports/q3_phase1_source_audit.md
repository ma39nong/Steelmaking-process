# Q3 Phase 1 — repository and source-material audit

Status: complete, provisional baseline scope only
Date: 2026-11-28 (Asia/Shanghai)

## Source priority and boundary

This phase uses the GitHub repository as the authoritative source for code, first-question paper, data handoff, Release manifest, and Release asset. The local `problem2` directory is not used as an evidence source. The ten original XLSB files are read from the supplied `Data` directory and remain read-only.

All future outputs in this branch are marked `Q3 provisional baseline`. No second-question algorithm, dynamic delay result, or unfrozen second-question Release is used as a formal input.

## Repository identity

- Remote: `https://github.com/ma39nong/Steelmaking-process.git`
- Working directory: `G:\亚太杯\problem3`
- Local branch: `problem3`, created from `origin/main`
- `origin/main`: `52e648ef31b489725344ce7b9d29aff28004569b`
- `origin/problem1`: `399c4edf2ec317f23835b1e0fffc679498ff2a75`
- `origin/problem2`: `d08f9898e34c4c3c725ba35d624dd701ff664588`
- Git tag `problem1-final-data-v1`: `36c2f86a3778104425f3b99413f4b13ae3dd24fc`
- Git tag `problem2-final-data`: `cea8248ea57633b3c736fd4a52fc5d5a9e63254e`

The GitHub releases REST endpoint returned HTTP 403 rate limiting during this audit. The Release asset itself was downloaded directly and checked against its attached `SHA256SUMS`.

## Documents read from GitHub or official supplied files

- `audit/source_materials/problem_b.txt` extracted from `2026 APMCM Problem B.pdf`.
- `audit/source_materials/essay_format.txt` extracted from `2026 APMCM Essay format and submission specifications.pdf`.
- `audit/source_materials/problem1_formal_paper_github.md` read from the GitHub `origin/problem1` branch path `problem1/final/reports/问题一_中文论文_V4_机制与结果图优化版.md`.
- GitHub Problem 1 `README.md`, `TASK.md`, `final/delivery/README.md`, `data_manifest.csv`, `data_package_manifest.json`, and `SHA256SUMS` copied into `audit/source_materials/` with source-specific names.
- The official APMCM AI-tool notice was recorded in `audit/source_materials/ai_tool_use_rules.md` from the committee notice published 2026-11-04. The final contest package must disclose AI-assisted work according to that notice.
- The 2022–2025 excellent-paper PDFs were inventoried and hashed. Their extracted text was inspected for structure and relevant modeling/reporting patterns. Full extracted text was kept only in a temporary audit cache and was not committed.

## Question 3 interpretation

The problem asks for real-time prediction of delayed flue-gas information from converter-mouth flame spectra, using the complementary time-domain relationship between immediate spectral/oxygen information and delayed CO₂/CO observations. Models trained on different attachments may differ, so a combined model or grouped models must generalize across all ten attachments. The Q3 baseline therefore needs causal spectral features, explicit future horizons, cross-heat validation, and a documented distinction between predictive lead time and physical transport delay.

The supplied problem statement specifies a 0.5 s sampling step, 2047 wavelength samples, nominal spectral range wording from 340.535 nm to 1027.678 nm, and CO₂/CO as instantaneous volume-flow concentrations in m³/s. The first-question GitHub handoff records that the actual wavelength header must be used and that many Release arrays are offline reconstructed outputs. For formal causal prediction, `raw` is the source; `reconstructed` is restricted to offline comparison and diagnostics.

## Release and raw-data identity

The Release asset downloaded from `problem1-final-data-v1` is:

- ZIP bytes: `548825449`
- ZIP SHA256: `ad9d74c7e9a5be182fe2de7b998eb5762cce5585cef10bf39a38e08b4acdaa55`
- Attached Release checksum: `ad9d74c7e9a5be182fe2de7b998eb5762cce5585cef10bf39a38e08b4acdaa55`
- ZIP members: 32
- Internal manifest SHA256: `359b2937f23497b674e1ee3eda581e2dce7745ced2998f878ea6cf44e7c77bc5`
- Internal manifest: `data_source_commit=3366db522887b2ca0f5f461ea6bd52ed71def1ba`, `packaging_base_commit=36c2f86a3778104425f3b99413f4b13ae3dd24fc`, 10 heats, 1043 actual repairs, offline `true`.
- All 31 files listed by the internal manifest matched byte length and SHA256.

There is a repository-version discrepancy that must remain visible. The Git tag copy of `data_package_manifest.json` has SHA256 `e5fdbc23b5c50431e4318e80c0a5ee1ba532a0a566f19f52d826b2b6112891ca` and records a different packaging base commit and payload byte count. The downloaded Release asset is internally consistent and is the data identity used for Q3. This discrepancy will be rechecked before any final Q3 data freeze.

The ten supplied XLSB files match the GitHub `data_manifest.csv` SHA256 values for all ten heats. The detailed comparison is in `audit/source_inventory.json`. The first filename contains a duplicated `.xlsb` suffix locally, but its SHA256 matches the manifest and it maps to heat `6B03205`.

## Planned formal inputs

- From the verified Release interfaces: `raw`, `time`, `wavelength`, `reconstructed`, and quality/mask fields for audit and offline comparison.
- From the ten original XLSB files: `Sequence`, cumulative oxygen consumption, CO₂, and CO, with field names, units, time axes, heat IDs, and sampling intervals rechecked in the next phase.
- No Q2 cleaned labels, dynamic alignment outputs, or Q2-derived time-delay labels are used in the current baseline.

## Phase 1 decision

The repository and source-material audit is complete and reproducible from the recorded Git refs, hashes, and source files. The next phase is the data-interface and raw-gas label audit. It will begin only after the user explicitly approves moving to Phase 2.
