# Q3 independent baseline plan

## Goal
Build and audit an independent APMCM Problem B Question 3 provisional baseline on branch `problem3`, using only the Problem 1 final Release and the ten original XLSB files. Complete one phase at a time, push each completed phase, report, and wait for explicit user approval before the next phase.

## Safety constraints
- Never modify or merge `main`, `problem1`, or `problem2`.
- Never overwrite existing Releases or use an unfrozen Q2 Release as a formal data source.
- Preserve raw inputs as read-only. Mark all outputs `Q3 provisional baseline`.
- Use only causal/current-history features for formal prediction inputs.

## Phases
1. Repository and source-material audit — complete, awaiting user approval to push/advance
2. Data interface and raw-gas label audit — pending user approval
3. Q3 feature/lead-time exploratory analysis — pending user approval
4. Switchable labels and classical baselines — pending user approval
5. CNN/TCN lightweight deep baselines — pending user approval
6. LOHO validation, transfer matrix, and fusion — pending user approval
7. Final audit, README, and reproducibility package — pending user approval

## Phase 1 success criteria
- `problem3` is a local branch created from origin/main with no uncommitted unsafe changes.
- Remote branch existence and tag/release identities are recorded.
- Required source documents are inventoried and relevant text extracted.
- Problem 1 commit, Release tag, and available data fingerprints are recorded.
- Q2 materials are explicitly treated as provisional context only.
- Phase 1 audit report is committed and pushed to `origin/problem3`.

## Errors encountered
| Error | Attempt | Resolution |
|---|---:|---|
| Root workspace was not a Git repository | 1 | Cloned remote repository into independent `problem3` directory |
| GitHub REST releases endpoint returned HTTP 403 rate limit | 1 | Use git refs/tags and repository contents for identity; record API limitation |
