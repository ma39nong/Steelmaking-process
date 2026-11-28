# Weak-form delay identification

Independent study based on `fd47730`, on branch `problem2`. Same first-order plus pure-delay dynamics, different identification moments. Production and frozen artifacts are protected byte for byte.

Run from repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m problem2.research.weak_integral.run
.venv/bin/python -m pytest problem2/research/weak_integral/test_weak.py -q
```

Primary sin-squared weight support: 20 s, stride 5 s. Fixed 10/40 s sensitivity is never a selection step. Evaluation seeds 9201–9210 are disjoint from prior benchmark seeds. The configuration and implementation hashes are written before evaluation. No post-evaluation tuning.

Read `results/REPORT.md` and `results/decision.json` for the pass/fail decision. Conditional parametric Bootstrap intervals and tolerance profile ranges have different meanings. Real-data exploration runs only if all controlled checks pass; no automatic formal-model replacement or flow-delay regression.
