# Progress log

## 2026-11-28 Asia/Shanghai
- Created independent local clone at `G:\亚太杯\problem3`.
- Created local branch `problem3` from `origin/main` at SHA `52e648ef31b489725344ce7b9d29aff28004569b`.
- Confirmed remote refs for `main`, `problem1`, `problem2`, and tags `problem1-final-data-v1`, `problem2-final-data`.
- Started Phase 1 source-material audit.
- Corrected source priority so GitHub branch/tag/Release contents are authoritative. Local `problem2` is not used as evidence.
- Downloaded and checked the GitHub `problem1-final-data-v1` Release asset. ZIP and internal manifest are self-consistent. A tag-manifest versus Release-asset manifest discrepancy was recorded for later recheck.
- Verified all ten local original XLSB SHA256 values against GitHub `data_manifest.csv`.
- Wrote `reports/q3_phase1_source_audit.md` and `audit/source_inventory.json`.
- Fresh Phase 1 verification passed: all 10 raw hashes matched, all 31 Release manifest entries matched bytes and SHA256, ZIP checksum matched attached Release `SHA256SUMS`, and current branch was `problem3`.
