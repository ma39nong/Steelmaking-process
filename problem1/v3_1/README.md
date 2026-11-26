# Problem 1 V3.1：点级异常最终验收

从 `problem1@641ea98` 继续，V1/V2/V3文件全部保留。本轮只审查和验证点级预处理，不修改低秩、双线性插值或饱和模型。

**验收结论：保守点处理未通过新区域端到端确认，正式回退V3-P0。正式结果属于异常诊断与保真连续重构，没有实施实际点异常消除。** 这次回退依据污染恢复实验失败，不以“缺少真实故障标签”排除统计证据。仍不能因此断言4个高证据点是真实工艺信号，也不能声称所有异常已去除。

- [中文最终报告](reports/问题一_最终模型与实验结果.md)
- [数据预处理最终验收](reports/数据预处理_点级异常最终验收.md)
- [冻结策略](config/FINAL_FROZEN.json)、[预声明协议](config/protocol.json)
- [独立确认门槛结果](results/selection_gates_confirmation.csv)
- [工程/科学验收](results/ACCEPTANCE.json)

## 数据与证据

共27,816,683个记录、51,444个统计候选（以results/raw_point_audit.csv的实际合计为准）。自动强规则提出5点；逐点复核将1个持续谱线候选降级保留，剩余4点称“高可信统计异常”，不是已确认仪器故障。由于策略级确认未通过，正式仍保留这4点；其受限低秩修复仅作为敏感性附录。

独立清单：

- `results/high_confidence_statistical_points.csv`：4个逐点证据及复核依据。
- `results/general_suspect_points.csv`：一般疑似及持续谱线否决点。
- `results/retained_points_rle.csv`：其余保留点的无损游程清单；每行指定炉次、原始行号、包含两端的列区间，含时间和波长范围。行号和列号均从0开始，区间包含两端。区间可逐单元展开，避免为2700余万点生成庞大重复CSV。疑似顶码另在接口中独立标记。
- `results/audit_before_line_review/`：最初5点自动结果，完整保留。
- `config/native_review_vetoes.json`：逐点只减不增的人工复核否决账本，确认前锁定，不用于删除注入试验的假阳性。
- `results/counterfactual_repairs.csv`：仅敏感性修复，`formal_repair_applied=False`。

正式每炉 `data/HEAT_interface.npz` 保存原始矩阵、清洗矩阵（本轮未删除）、正式重构、V3潜变量、条件性饱和估计、统计候选/分级、实际修复掩码、差值矩阵、坐标和质量标记。`HEAT_difference.npz`保存完整差值矩阵（本次全零）；`HEAT_dense.npz`保存连续查询得到的密集网格。数据数组按V1–V3既有约定留在本地data/，代码、参数锁、完整CSV、哈希、报告和图表入Git。

查询复用未修改的 `problem1/v3/src/interface.py` 中 `Spectrum.load()` 与 `evaluate()`：

```python
import sys
sys.path.insert(0, 'problem1/v3/src')
from interface import Spectrum
m = Spectrum.load('problem1/v3_1/data/6C02058_interface.npz')
values = m.evaluate([445.5, 445.75], [767.07, 770.0])
```

`raw_until()`仅返回原始历史切片，未实现因果预处理。离线重构不能直接用于第三问在线性能评价。质量位沿用V3，16384预留统计点修复、32768预留数值回退；本次两位均不触发。quality_grade的2若未来采用处理表示统计掩码，不是仪器故障概率或日志确认。

## 有限验证与复现

新的32×256片段位于旧V2第3/4片段之间；列64开发，列900/1720确认。精确坐标与SHA见results/region_index.csv；此前十炉已研究，因此不是外部炉次验证。开发种子314159/271828，确认161803/141421。仅比较P0、保守、历史P1；孤立上冲/突降，各0.3/0.7幅度，污染密度0.5%，两个种子。注入前带噪原谱作为条件参考；不利用检测器剔除参考假阳性。

全部数值命令设置 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`，沿用根`.venv`。以下列出实际研究顺序，锁文件和历史结果不应覆盖；复现请另建运行目录。

```bash
.venv/bin/python problem1/v3_1/src/initialize.py
.venv/bin/python problem1/v3_1/src/audit_points.py
.venv/bin/python problem1/v3_1/src/native_review.py
.venv/bin/python problem1/v3_1/src/experiments.py development
.venv/bin/python problem1/v3_1/src/line_control.py development
.venv/bin/python problem1/v3_1/src/select.py development
.venv/bin/python problem1/v3_1/src/experiments.py confirmation
.venv/bin/python problem1/v3_1/src/line_control.py confirmation
.venv/bin/python problem1/v3_1/src/select.py confirmation
# finalize.py逐炉执行，再aggregate；native_sensitivity.py仅对4个有保留强点的炉次执行
for heat in 6B03205 6B03206 6B03270 6B03271 6C02057 6C02058 6C02059 6C02060 6C02107 6C02108; do
  .venv/bin/python problem1/v3_1/src/finalize.py "$heat"
done
.venv/bin/python problem1/v3_1/src/finalize.py aggregate
for heat in 6B03205 6B03270 6B03271 6C02059; do
  .venv/bin/python problem1/v3_1/src/native_sensitivity.py "$heat"
done
.venv/bin/python problem1/v3_1/src/figures31.py
.venv/bin/python problem1/v3_1/src/figures31.py sensitivity
.venv/bin/python problem1/v3_1/src/report31.py
.venv/bin/python -m pytest problem1/v3_1/tests -q
.venv/bin/python problem1/v3_1/src/integrity31.py
```

所有十炉的正式处理前后曲线、局部热图、全矩阵差值和独立三维图在figures/appendix；技术路线、端到端误差图和十炉总览在figures/main。正式前后重合是未实施修复的真实结果，不伪造变化。4点敏感性修复图明确标注“非正式”。

原始协议、源文件、选型及复核账本SHA分别锁定。指标澄清记录在config/METRIC_DEFINITION_NOTE.md：区分未消除的污染峰与处理新引入峰；最初开发结果保留于metric_definition_audit。开发可行结果未替代确认失败结果，确认后未调参数。
