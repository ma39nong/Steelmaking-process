# 目标激励检查位置的独立 A/B 验证

研究正式V2使用未平移烟气窗口判定激励强度的影响。[报告](results/REPORT.md)与[决策记录](results/decision.json)给出结果。正式冻结版 `d08f989`、第一问、NPZ接口和既有研究均不改动。

## 复现

在仓库根目录运行：

```bash
.venv/bin/python -m problem2.research.target_excitation.run
.venv/bin/python -m pytest -q problem2/tests problem2/research/response_broadening/test_response.py problem2/research/real_event_validation/test_review.py problem2/research/target_excitation/test_excitation.py
```

使用已有NumPy、SciPy、pandas、Matplotlib、Pillow、pyxlsb环境。严格字体配置需要Fontconfig、真实Times New Roman及pdffonts。输入沿用十炉原始XLSB和前轮冻结NPZ，逐文件核对哈希；论文冻结清单中的100个文件在研究前后全部核验。字体、原始XLSB和本地正式NPZ不提交Git。

配置见 [config.json](config.json)。不调原阈值、不调窗口或搜索范围。结果始终原位更新，不创建平行版本目录。

## 比较定义

- A先检验源和未平移目标窗口的标准差，再搜索正标准差目标lag。与生产函数逐字段比对，同时复现冻结NPZ的24个合格窗口。
- B只保留对应平移目标窗口达到原0.05全局标准差阈值的lag。接受后置信权重使用该平移窗口的强度，置信分数不参与主指标筛选。
- B剔除弱lag会改变竞争峰集合。因此同时记录保留所有原始有限竞争峰时，同一B胜出lag是否仍满足原分离度；只作机制诊断，不引入第三种正式模型。
- 不用零或负无穷填充掩蔽的相关曲线来计算峰宽。B的峰宽、显著性在完整原始相关曲线上计算，仅作诊断。

## 实验单位与口径

合成试验30个固定种子。源为8秒紧支撑、峰值1的非负瞬态；基线1000，噪声标准差为ratio乘完整瞬态序列的标准差。中心响应与晚响应分别将源放在中心和中心后20秒，5、15、30、55秒延迟及3种噪声交叉；另有极弱复制、两类无关信号与零/60秒边界参考。不重跑前轮18场景全套实验。

真实背景注入直接在V2入口的冻结流量和烟气上添加同形非负瞬态，隔离本次代码假设，不重新求导累计氧量，也不解释为实际氧枪操作。每炉气体每种子只均匀抽取一个共同中心，12种剂量延迟组合共享背景。正确恢复为接受且误差≤1秒；误接受和错误时滞分开。7200种配对条件不是7200次独立工业事件。

600次源单独注入对照每炉气体每种子只计一次，不重复计入12种条件。无注入响应接受数不能等同于工业误报率。对原始数据不计算事件Recall、Precision或F1。

## 输出

| 文件 | 内容 |
|---|---|
| `baseline_replay.csv` | A对生产函数的12字段、20条序列精确/数值重放核查 |
| `real_windows.csv`、`real_changed_windows.csv` | 原始窗口的A/B输出与接受变化；空变化文件明确保留表头 |
| `synthetic_trials.csv`、`injection_trials.csv` | 所有指定中心、接受和拒绝、已知延迟、数值误差与激励比值 |
| `synthetic_summary.csv`、`injection_summary.csv`、`scenario_summary.csv` | 各场景、剂量、延迟的恢复、拒绝及错误指标 |
| `paired_transitions.csv`、`transition_summary.csv`、`injection_by_heat.csv` | 配对救回、损失、逐炉汇总及竞争峰掩蔽诊断 |
| `original_gate_diagnosis.json` | 原始CO被未平移目标门拒绝的34窗口及重新搜索最大相关；CO₂对应数量为0 |
| `figures/`、`figure_captions.csv`、`font_validation.json` | 英文图、真实Times New Roman、PNG/SVG/PDF及独立可编辑图注 |
| `layout_preview.png` | 唯一版宽检查预览 |
| `effective_config.json`、`protocol.json`、`input_manifest.json`、`provenance.json` | 实际配置、比较定义、输入与输出SHA256、代码指纹 |

研究结果不自动进入正式论文或第三问接口。决定是否更改冻结实现，必须明确区分代码门控行为改进和原始工业时滞证据改善。
