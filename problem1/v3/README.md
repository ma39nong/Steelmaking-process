# Problem 1 V3：最终定型与预处理审核

基于 `problem1@d5b8f5e`；V1/V2的1049个已跟踪文件哈希逐项保留。正式目标是带质量标记的**仪器记录尺度离线连续曲面**，不是无噪声或超量程真实光谱。

- [中文最终研究报告](reports/问题一_最终模型与实验结果.md)
- [数据预处理专项审核](reports/数据预处理专项审核报告.md)、[异常分级](reports/异常类型与判定标准.md)
- [数学推导与来源](reports/最终数学模型推导与算法来源.md)、[文献获取边界](references/README.md)
- [科学验收](results/SCIENTIFIC_ACCEPTANCE.json)、[工程验收](results/ENGINEERING_ACCEPTANCE.json)

最终P0不删除未确认统计候选；正式V3-C使用秩4、Huber尺度0.003的缺失拟合，随后保真/非负局部支撑投影和双线性插值。当前实际无确认无效点，因此正式原网格恰为原始记录，不开展默认降噪。缺失预测由独立遮蔽检验。1,129帧候选保留，真实整帧故障确认0；31个V2负值被消除。饱和记录14,805点保留65535，条件性局部凸估计单独保存，上界未知。6C02058遮蔽拟合1200步仍不稳定，预声明回退改善稳定却损失精度，不宣称全面解决。

## 运行环境与顺序

沿用根 `.venv` 与 `problem1/requirements.txt`。所有数值命令设置 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`。原始XLSB/缓存沿用V1只读文件。以下为已执行研究流程，**不要删除锁文件或在原目录覆盖后重新宣称独立确认**；复现应复制配置到新运行目录并保留新日志。

```bash
.venv/bin/python problem1/v3/src/initialize.py
.venv/bin/python problem1/v3/src/frame_review.py
.venv/bin/python problem1/v3/src/preprocessing_experiments.py development
# 开发选择后再执行确认
.venv/bin/python problem1/v3/src/preprocessing_experiments.py confirmation
.venv/bin/python problem1/v3/src/preprocessing_report.py
.venv/bin/python problem1/v3/src/confirmation.py gate
.venv/bin/python problem1/v3/src/confirmation.py confirmation
.venv/bin/python problem1/v3/src/conditional.py
# 四炉各执行一次 stability.py HEAT；有限预算至1200
.venv/bin/python problem1/v3/src/stability.py 6C02057
.venv/bin/python problem1/v3/src/stability.py 6C02058
.venv/bin/python problem1/v3/src/stability.py 6C02107
.venv/bin/python problem1/v3/src/stability.py 6C02108
.venv/bin/python problem1/v3/src/acceptance_experiments.py
.venv/bin/python problem1/v3/src/final3.py
.venv/bin/python problem1/v3/src/figures.py
.venv/bin/python problem1/v3/src/cost_and_curves.py
.venv/bin/python problem1/v3/src/write_final_report.py
.venv/bin/python problem1/v3/src/integrity.py
.venv/bin/python -m pytest problem1/v3/tests -q
```

`stable_completion.complete_with_stability` 为后续新缺失矩阵提供600/1200步比较及预声明逐点线性回退；未改写固定100步确认结果。回退不使用隐藏真值。正式十炉没有确认缺失，不触发回退。

## 查询接口与数据对象

```bash
.venv/bin/python problem1/v3/src/interface.py 6C02058 --time 445.5 445.75 --wave 767.07 770
# 条件性估计必须显式请求
.venv/bin/python problem1/v3/src/interface.py 6C02058 --time 445.5 --wave 767.07 --target conditional
```

Python可用 `Spectrum.load('.../data/HEAT_interface.npz')`，`evaluate(times, wavelengths)` 返回两组坐标笛卡尔积，禁止域外外推。时间s、波长nm、强度ADC counts；原始波长非均匀。函数C0连续，不全局可微。`raw_until(cutoff)` 只给原始历史切片，**未执行因果预处理**，不能将离线质量诊断或离线补全作为在线无泄漏输入。

每炉本地 `data/HEAT_interface.npz` 保存：time、wavelength、raw、cleaned（仅确认无效置NaN）、reconstructed（记录尺度）、latent（低秩潜变量，不能当正式谱）、conditional（仅顶码位置替换条件估计）、valid（低秩拟合集，不含顶码）、invalid、quality、quality_grade、suspected_saturation、saturated_unidentifiable、conditional_solver_converged及模型参数。

`HEAT_dense.npz` 为双倍分辨率采样连续函数，recorded为float32显示网格；精确查询使用float64节点与公式。`HEAT_quality.npz` 补充帧/通道/点候选、等级、掩码及独立理由账本。data/沿用项目大数组不入Git约定，本地已全部生成；结果CSV、图表、源代码、哈希与参数锁入Git，可由原始附件复现。

质量位可叠加：2=顶码；4=补全；8=局部支撑投影；16=负潜变量修正；32=无轴支撑；64=非有限/负原始结构无效；128=原始顶码类型；256=点候选；512=确认无效；1024=帧候选；8192=条件性饱和局部求解未收敛。quality_grade为0保留、1疑似、2无效；不同于概率置信度。数值回退另存numerical_fallback掩码；正式数据未触发。通道候选是诊断，不自动整列删除。

## 实验记录

`config/protocol.json` 在评分前声明；`PROTOCOL_LOCKED.json`、`NUMERICAL_LOCKED.json`、`final_preprocessing_policy.json`、`FINAL_MODEL_LOCKED.json` 保存时间/哈希/选型。新确认使用十炉此前未用于V2局部实验的间隙，**不是外部炉次验证**。

`preprocessing_comparison.csv`、`frame_anomaly_candidates.csv`、`channel_anomaly_candidates.csv`、`final_preprocessing_policy.json` 在本目录根部提供专项交付副本。完整结果见results/。开发门槛实现复核及首次输出保存在results/development_gate_audit/，未覆盖历史，也未使用确认分数重选。V1/V2回归测试分别执行，避免同名历史模块路径互相污染。
