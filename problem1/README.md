# 第一问：真实光谱异常处理与连续重构

> 当前正式模型和数据见[final研究说明](final/README.md)，中文正文见[正式章节版](final/reports/问题一_中文论文_V2_正式章节版.md)，交接见[处理后数据说明](final/reports/数据交接与后续问题使用说明.md)。下文完整保留V1的复现命令与输出定义；不作为当前正式数据入口。

本项目仅使用十炉火焰光谱，不把氧耗、CO₂、CO作为拟合输入。原始文件位于仓库根目录`Data/`，附件1以指定SHA256验证；原始记录复制到本机`data/raw/`后只读。

## 复现

从仓库根目录运行（Python3.12）：

```bash
python3 -m venv .venv
.venv/bin/pip install -r problem1/requirements.txt
.venv/bin/python references/acquire.py
.venv/bin/python references/build_index.py
.venv/bin/python problem1/src/audit.py
.venv/bin/python problem1/src/audit_extended.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/src/eda.py
.venv/bin/python problem1/src/reports_initial.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/src/experiment.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/src/final.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/src/sensitivity.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/src/benchmark.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/src/rank_audit.py
.venv/bin/python problem1/src/peak_audit.py
.venv/bin/python problem1/src/report.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q problem1/tests
```

首次文献检索须阅读根目录`references/README.md`和`reports/00_*`；`LITERATURE_REVIEW_COMPLETED`仅表示已执行检索与范围记录，不代表所有原文取得。全文缺口详见笔记。

实验固定种子20261126，开发/测试掩码分别保存。参数开发完成后锁定模型族，独立测试不参与选择。可恢复的`*_development.json`、`*_outer.json`缓存属于当前配置；改变源文件、配置或数据哈希后须移除对应缓存和COMPLETE标记后重跑，不能沿用旧结果。

## 输出

- `../data_manifest.csv`：唯一正式数据映射与SHA256。
- `results/audit`：全数据与逐帧/通道审计。
- `data/processed`、`data/masks`：逐炉独立真实轴、清洗矩阵与类型掩码。
- `results/model_comparison/model_comparison.csv`：10炉×3模型×5测试场景。
- `results/validation`：遮蔽划分、异常注入、BEADS、参数敏感性与人工截顶。
- `data/reconstructed/*_model.npz`：可在域内任意位置计算的连续模型参数。
- `*_dense.npz`：连续模型在双倍密集评价网格上的值，并非仅图片插值。
- `*_uncertainty.npz`：按遮蔽形态的经验误差半径与真实饱和不可辨识掩码。
- `figures/reconstructed`：10张独立3D图、10张热力图和总览，PNG/SVG。
- `reports/问题一_数据分析与模型建立及实验结果.md`：中文完整报告与局限。

加载任意连续模型：

```bash
.venv/bin/python problem1/src/evaluate.py problem1/data/reconstructed/6B03205_model.npz --times 100.25 100.75 --wavelengths 589.1 770.2 --output /tmp/spectrum-evaluated.npz
```

模型B仅在遮蔽处使用低秩补全，对有效观测施加保真投影，然后连续插值；保存潜在低秩层以便检查。模型A/B保存重构结点网格与真实轴，模型A使用三次张量样条，模型B先沿时间再沿波长执行形状保持PCHIP；模型C直接保存低维系数与节点。禁止域外外推。三维绘图采用较稀评价网格以便阅读，完整密集数据另存。所有图的Z值都来自对应连续模型。

## 阅读限制

R01/R02/R05/R06/R07/R11/R12本次全文未得；不得将本项目称为这些原论文的完整复现。R03对照是全局低秩+逐点稀疏的二维简化，并未实现原文非凸局部log/三维TV。真实饱和区没有上界真值，经验误差带不提供其严格覆盖保证。详见中文报告。
