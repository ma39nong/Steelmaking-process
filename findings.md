# Findings

## Repository
- Local source workspace: `G:\亚太杯\problem3`.
- Remote: `https://github.com/ma39nong/Steelmaking-process.git`.
- `problem3` was absent remotely at audit start and was created locally from `origin/main`.
- `origin/main` initial commit SHA: `52e648ef31b489725344ce7b9d29aff28004569b`.
- Remote `problem1`: `399c4edf2ec317f23835b1e0fffc679498ff2a75`.
- Remote `problem2`: `d08f9898e34c4c3c725ba35d624dd701ff664588`.
- Tag `problem1-final-data-v1`: `36c2f86a3778104425f3b99413f4b13ae3dd24fc`.
- Tag `problem2-final-data`: `cea8248ea57633b3c736fd4a52fc5d5a9e63254e`.
- GitHub releases REST endpoint was rate-limited (HTTP 403), so release identity will be verified from tag contents and repository manifests.

## Scope
The current work must remain independent of Q2 final conclusions. Q2 reports can explain status and interfaces, but no Q2-derived labels or data are formal inputs.

## Source-priority correction
- GitHub repository contents are authoritative for first-question code, paper, handoff documentation, tag manifests, and Release assets.
- Local `G:\亚太杯\problem2` was used only as a navigation hint before this correction. It is not an evidence source for the Q3 baseline.
- The GitHub formal paper was read from `origin/problem1:problem1/final/reports/问题一_中文论文_V4_机制与结果图优化版.md` at the remote Problem 1 branch commit.
- The downloaded Release asset was independently hashed and inspected. Its internal manifest is self-consistent, but the tag copy of `data_package_manifest.json` differs from the Release asset manifest. The Release asset hash and internal manifest hash are recorded as the data-source identity for this work.
