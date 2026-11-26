# Problem 1 V2

从 `problem1` 分支 `d8817f2` 继续；V1目录内容和所有历史结果不变。V2全部代码/协议/结果位于本目录，原始观测直接只读加载V1 `data/processed/{heat}.npz`；不使用烟气或氧气。

本轮锁定候选并未通过全部科学性能验收，尤其整帧/整列异常、真实饱和和潜变量人工峰风险；详见报告，不能只凭平均误差宣称正式模型成立。

完整报告：`reports/问题一_V2模型深化与严格验证.md`。审计：`reports/00_V1模型审计.md`。公式：`reports/01_数学推导.md`。新增文献：`references/README.md`；原R01–R12仍在仓库根references/。

## 运行

环境沿用V1 `problem1/requirements.txt` 和根`.venv`。数值命令使用 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`。从仓库根运行：

```bash
.venv/bin/python problem1/v2/src/audit_v1.py
.venv/bin/python problem1/v2/src/develop.py
.venv/bin/python problem1/v2/src/capacity.py
.venv/bin/python problem1/v2/src/c_expand.py
.venv/bin/python problem1/v2/src/continuous.py development
# 以上开发选择锁定后才能运行测试；现有锁文件不应删除覆盖。
.venv/bin/python problem1/v2/src/test_experiments.py
.venv/bin/python problem1/v2/src/c_expanded_test.py
.venv/bin/python problem1/v2/src/continuous.py test
.venv/bin/python problem1/v2/src/v1_matched.py
.venv/bin/python problem1/v2/src/saturation.py
.venv/bin/python problem1/v2/src/real_saturation.py
.venv/bin/python problem1/v2/src/convergence.py
.venv/bin/python problem1/v2/src/full_confirmation.py
.venv/bin/python problem1/v2/src/collect_full_confirmation.py  # 并行炉次缓存完成后统一核对550行
.venv/bin/python problem1/v2/src/full_continuous.py
.venv/bin/python problem1/v2/src/final_v2.py
.venv/bin/python problem1/v2/src/structure_diagnostics.py
.venv/bin/python problem1/v2/src/time_sensitivity.py
.venv/bin/python problem1/v2/src/peak_retention.py
.venv/bin/python problem1/v2/src/solver_diagnostics.py
.venv/bin/python problem1/v2/src/summarize.py
.venv/bin/python problem1/v2/src/figures.py
.venv/bin/python problem1/v2/src/report.py
.venv/bin/pytest problem1/v2/tests -q
```

`develop.py` 和 `c_expand.py` 拒绝覆盖已锁选择。如需全新复现实验，请复制V2源代码与配置到新实验目录，指向同一只读V1输入，不要删除本轮结果。`c_develop.py`的两个极限C配置作为开发先导记录保留；最终C选择使用`c_expand.py`完整有限网格。它们均未根据测试调整。

## 结果语义

- `config/protocol.json`：独立局部片段主协议。`LOCKED.json`锁定B0–B5和检测器；`C_EXPANDED_LOCKED.json`锁定扩大C；`FINAL_SELECTION_LOCKED.json`锁定最终离散模型；`INTERPOLATION_LOCKED.json`独立锁定插值方法。
- `test_completion.csv`：B0–B5；`c_expanded_test.csv`、`v1_matched_test.csv`：相同可见性比较。按炉次先平均重复实验，再报告宏平均与炉间标准差；重复遮蔽不是独立炉次。
- `full_confirmation.csv`：全原始矩阵确认性实验，模型不再调整；完整几何和局限见`full_confirmation.json`。`full_interpolation.csv`：全矩阵降采样与组合方法确认，见独立配置。
- `artificial_censoring.csv`、`real_saturation_comparison.csv`、`focus_local_convex_saturation.csv`：人工可知真值与真实不可辨识饱和严格区分。
- `convergence_curves.csv`：四炉、两求解器、目标/变化/遮蔽误差/耗时。`capacity_development.csv`：容量/时间惩罚/Huber受控诊断。
- `data/{heat}_model.npz`：最终连续函数所需真实轴/网格/方法/潜变量/参数；`data/{heat}_masks.npz`：异常类别、训练掩码、真实饱和不可辨识标记。
- 大型模型文件只保存本地，未强行推送Git；代码、完整CSV、锁文件和报告/图片随分支推送，可用只读原始数据重新生成。

调用接口示例：

```bash
.venv/bin/python problem1/v2/src/evaluate.py 6C02058 --time 100 100.25 --wave 400 589
```

V2掩码types位：1=输入不可用，2=65535疑似饱和，4=双向MAD局部异常，32=低秩残差联合连通支持。整帧/持续通道辅助候选单独保存CSV，不直接删除。

不允许定义域外外推。有效观测保持误差不是独立预测精度。局部波长隐藏区的峰指标是该区域内峰指标，不能解释为全谱主峰；降采样时间行的谱形指标涵盖完整可用波长。
