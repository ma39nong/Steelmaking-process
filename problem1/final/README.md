# 第一问最终异常修复与中文正文

[中文论文初稿](reports/问题一_中文论文初稿.md)按问题分析、数据、异常识别、低秩修复、求解、十炉结果、检验与评价组织；[历史复核](reports/历史复核与最终决策依据.md)单独保存，不以软件版本流水账替代正文。

最终采用开发阶段锁定的**局部尺度自适应多证据评分（阈值4）＋原V3固定秩4 Huber补全＋保真/局部非负支撑投影＋双线性连续函数**。十炉实际修复1,043点，原始值、估计值、坐标和统计证据独立保存。四个历史高证据尖峰均被新评分重新识别并实际修复。整帧及持续通道只诊断、不盲目整行/整列删除。疑似顶码保持原记录并单独作条件性估计。

V1至V3.1的代码、数据、结果、锁与报告完整保留。历史数据和新比较协议不同，不直接用不同版本历史分数判优。当前目录为独立研究结果，不覆盖旧文件。

## 已执行实验

- 三波长区域、十炉、开发/确认独立时间片段与种子；每阶段240污染案例，各策略共同输入。
- A原始保留、B原V1双向MAD、C严格多条件、D自适应评分；D仅阈值4、6两候选。
- 统一300步补全预算选择检测策略；另对同一固定确认案例检查实际发布的600/1200步及数值回退配置，不称第二批独立验证。
- 不再使用全局Recall/NRMSE绝对门槛一票否决所有局部修复。预声明连续代价衡量污染恢复与误改代价，确认后不调参。
- 同一十炉历史已研究，本轮为新区域/新注入确认，不是外部炉次验证。参考原谱带噪，注入实验不等于真实故障准确率。

主要数据：`results/unified_comparison.csv`、`confirmation_per_heat_band.csv`、`exact_final_pipeline_confirmation.csv`、`actual_repairs.csv`、`native_mask_summary.csv`、`final_summary.csv`、`conditional_tiles.csv`、`ACCEPTANCE.json`。历史失败的尺度及条件交集归因见`V31_failure_attribution.csv`。

## 运行顺序

沿用根`.venv`。所有数值命令设置 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`。以下为已执行顺序；已有锁及实验结果不要覆盖，复现应复制到同级新运行目录。初始化只用于在基线提交登记协议；保留已有配置用于核对即可，不删除历史锁后宣称新独立测试。

```bash
.venv/bin/python problem1/final/src/initialize.py
.venv/bin/python problem1/final/src/history.py
.venv/bin/python problem1/final/src/experiment_final.py development
.venv/bin/python problem1/final/src/experiment_final.py confirmation
.venv/bin/python problem1/final/src/summarize.py
.venv/bin/python problem1/final/src/native_masks.py
.venv/bin/python problem1/final/src/controls_line.py
.venv/bin/python problem1/final/src/check_final_pipeline.py
.venv/bin/python problem1/final/src/run_native.py
.venv/bin/python problem1/final/src/conditional_final.py
.venv/bin/python problem1/final/src/export_results.py
.venv/bin/python problem1/final/src/plots.py
.venv/bin/python problem1/final/src/paper.py
.venv/bin/python -m pytest problem1/final/tests -q
.venv/bin/python problem1/final/src/integrity.py
```

`run_native.py`只并行独立炉次，使用原V3求解器，不创建新模型；每炉有限1200步、360s上限。原始ADC顶码、非负支撑、回退逻辑和双线性实现不变。条件性分块只对受修复掩码影响的块重算，其余相同目标复用V3结果；未达停止准则的块保留标记。

## 正式数据对象与接口

每炉本地`data/HEAT_interface.npz`保存真实time、wavelength、raw、cleaned（统计修复位置置NaN）、reconstructed、latent、conditional、valid、invalid、quality、quality_grade、statistical_repair、suspected_saturation、saturated_unidentifiable、conditional_solver_converged、numerical_fallback、score及模型参数。latent数组保存counts尺度的UL；公式中的L已除以65535。

`invalid`表示被排除后进行统计修复的位置，不表示已确认硬件故障。`valid`为低秩拟合集，不含疑似顶码；顶码记录在正式曲面中保留。`conditional`只在顶码节点替换潜在估计，其他位置等于处理矩阵。quality位：2顶码，4修复，8局部支撑投影，16负潜变量修正，32弱支撑，256一般疑似，1024帧候选，8192条件估计未达停止标准，16384统计修复，32768数值回退。grade0保留、1疑似、2统计修复，不是概率。

掩码文件为`HEAT_masks.npz`，含评分分量、自动/复核/实际修复及保留疑似。`HEAT_difference.npz`是完整处理−原始差值矩阵；`HEAT_dense.npz`为双倍分辨率的连续函数显示网格。`results/actual_repairs.csv`逐点列出原始测量与模型估计，行/列从0开始。

大数组沿用项目既有约定留在本地data/，均已生成；代码、全部指标、哈希、参数、论文、图表及参考索引入Git，可从原始附件复现。图表显示降采样不改变保存的全矩阵或精确查询。

```python
import sys
sys.path.insert(0, 'problem1/v3/src')
from interface import Spectrum
m = Spectrum.load('problem1/final/data/6C02058_interface.npz')
recorded = m.evaluate([445.5, 445.75], [767.07, 770.0])
conditional = m.evaluate([445.5], [767.07], target='conditional')
```

查询为坐标笛卡尔积，只在实际观测定义域内；时间s、波长nm、强度counts。函数C0连续、非负、无单元内超出四角值域的过冲，不全局可微。`raw_until`只返回原始历史切片，未执行因果预处理。第三问实时预测不能直接使用本轮含未来信息的离线结果。

## 文献与图表

[方法文献](references/README.md)沿用合法全文与核读笔记，区分固定秩Huber改造与原PCP/Soft-Impute。[写作组织核读](references/writing/组织方式阅读笔记.md)记录2022/2023 APMCM和2024/2025全国赛官方B题页面，后两者不混称APMCM、不转载图像正文。

figures/main供正文使用；figures/appendix包含十炉实际修复曲线、原始/处理热图、差值及三维曲面，另有帧级快变诊断。工程测试状态和科学收益/代价分别记录于验收文件，新增峰、漏检、数值回退及条件性估计局限保留在论文检验/评价部分。
