# 十炉真实数据事件级验证

以 `dea6637` 为冻结基准，只开展真实数据的候选重放、竞争图、采样更新、表示尺度及机会配对审查。不开发新合成算法，不更改正式模型或NPZ接口。

当前结论见 [报告](results/REPORT.md)，所有候选可在 [本地检索页](results/review.html) 筛选。用户已确认中文附件1与现有英文命名文件为同一原始文件，仅作改名；本轮核对SHA256后关闭该身份事项。前轮报告作为历史证据保留，其“待核验”状态由本轮确认取代。

## 复现

仓库根目录使用现有环境执行：

```bash
.venv/bin/python -m problem2.research.real_event_validation.run
.venv/bin/python -m pytest -q problem2/tests problem2/research/response_broadening/test_response.py problem2/research/real_event_validation/test_review.py
```

配置为 [config.json](config.json)，默认全部参数显式保存。不联网下载，不依赖新增包。需要原始十炉XLSB、本地前轮基线NPZ，以及已提交的前轮结果。原始文件和NPZ均逐文件比对前轮SHA256。当前机器路径由前轮输入清单继承；移动到另一台机器时需在同一已核验数据上重新生成前轮清单，不可静默替换。

九种尺度组合为事后敏感性分析，不用于根据匹配数修改正式参数。候选阈值的参考区间严格复用前轮：第120个采样点至45%记录位置，源与目标共用同一时间索引；不能改成从记录起点开始后仍声称精确重放。

## 结果解释

- `source_candidates.csv`：共享氧枪候选分别引用到CO和CO₂。候选编号不等于独立工艺事件。时间为秒，流量为m³/s，未确认事件不能进入时滞回归。
- `candidate_edges.csv`：同极性、非负0–60秒内的全部边，保留源与目标竞争度、三种时间差和响应宽度。唯一边仍是候选，不是物理标签。
- `baseline_24_window_review.csv`：24个冻结窗口的时滞、位置及所含变化候选数。12个源窗口重叠组单独绘图，不称12次独立响应。
- `previous_unique_pairs.csv`：逐条复核之前V2适配器与直接事件规则的候选配对，包含邻近源数与尺度重现数。V2只估计局部平移量，表中其他描述符来自公共候选检测器。
- `scale_sensitivity.csv`：同极性峰一对一时间一致性，不是准确率；九种配置全重现也不证明物理独立性。
- `circular_controls.csv`：每条序列99次循环移位，固定种子。接缝波瓣丢弃，不算正式假设检验或误匹配率。
- `sampling_cadence.csv`：原始累计计数重复比例、正增量与更新间隔。不凭此认定仪器工作机制。
- `figures/`、`figure_captions.csv`、`font_validation.json`：一套当前波形页，PNG/SVG/PDF及可编辑图注、字体嵌入记录；`layout_preview.png`是唯一总览。
- `attachment1_identity.json`：用户来源确认和原始文件SHA256。
- `effective_config.json`、`input_manifest.json`、`provenance.json`：实际配置、输入与代码指纹、完整输出哈希。

重新运行原位替换结果。原始数据、第一问、正式V2、旧研究目录及第三问接口均冻结。没有独立事件标签时，不报告真实数据Precision、Recall、F1，也不按候选数假造统计样本量。
