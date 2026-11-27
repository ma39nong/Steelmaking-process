# 第一问补充研究

[补充研究报告](reports/补充研究与模型依据.md)、[物理来源](reports/physical_sources.md)、[协议](protocol.json)、[结论](conclusions.json)。现有Release与1043点基线不变；data/为本地派生数组及明确独立的敏感性矩阵。

顺序：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/analyze.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/finite_checks.py
.venv/bin/python problem1/final/supplement_v1/src/report.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/publication_figures.py --stage priority
OPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/publication_figures.py --stage remaining
.venv/bin/python problem1/final/tables_v3/build_tables.py
.venv/bin/python problem1/final/supplement_v1/src/manuscript.py
.venv/bin/python problem1/final/supplement_v1/src/validate.py
```

通道回归采用炉内时间块区间；新实验按新位置分开发/确认；历史17峰只诊断回放，不重新用于调参。输入哈希及参考册/指令哈希已登记。图与表分别在../figures/publication_v3/和../tables_v3/。
