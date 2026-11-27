# Response Broadening and Event-based Delay Identification

独立研究，以 `9b8513a` 为基线。当前结论见 [研究报告](results/REPORT.md)。
正式V2、第一问、第三问接口及前轮审查成果保持冻结；本研究不拟合流量—时滞函数。

## 复现

在仓库根目录、安装现有第二问依赖以及 `threadpoolctl`、Pillow 的Python环境中运行：

```bash
.venv/bin/python -m problem2.research.response_broadening.run
.venv/bin/python -m pytest -q problem2/tests problem2/research/response_broadening/test_response.py
```

配置为 [config.json](config.json)。研究入口显式读取配置。
原始十炉文件需与前轮 `input_audit.json` 的哈希相同；基线NPZ读取自
`problem2/final/delay_audit/baseline/data/`，并与前轮输出清单逐文件核对。
若缺少该本地NPZ，先按第二问README复现前轮审查；不得用空数组或Release文件静默替代输入。
运行固定写入 `results/`，不自动创建版本或备份目录。

分阶段选项为 `--stage development`、`synthetic`、`real`、`finish`。
`finish` 要求前三阶段的代码与配置签名一致；参数或数值代码改变后需重跑相应数值阶段。
`all` 默认完整运行。开发与评价种子分开，未运行纯平移实验；前轮7种非理想场景复用原生成器和种子。
任何方法都不能从真实候选配对数直接推断事件Precision或物理独立性。

## Release证据

GitHub连接器未提供Release资产下载，使用仓库已有 `gh` 进行只读获取。
首次运行前将以下资产下载到本机缓存；ZIP和字体文件均不提交Git。

```bash
gh release download problem2-final-data --repo ma39nong/Steelmaking-process --pattern '*' --dir /home/user/.cache/apmcm-b/releases/problem2-final-data --skip-existing
gh release download problem1-final-data-v1 --repo ma39nong/Steelmaking-process --pattern data_package_manifest.json --pattern SHA256SUMS --dir /home/user/.cache/apmcm-b/releases/problem1-final-data-v1 --skip-existing
.venv/bin/python -m problem2.research.response_broadening.release_audit
```

先创建缓存目录。其他机器将 `/home/user/.cache` 换为自身缓存位置，代码使用 `Path.home()`。
`release_audit` 从GitHub刷新元数据；完整研究在结束阶段重用已保存的Release元数据快照，并重新验证缓存资产及570个NPZ字段。
已验证处理接口的一致性不等于已取得指定的新原始附件1：两个Release的清单均未列出原始XLSB，最终冻结门槛仍为待核验。

## 方法与评价口径

- 公共输入：q和气体瞬时流量，0.5 s采样。事件/核方法使用同一9点SG变化率；正式V2保留自己的差分表示。
- FIR：121个有符号系数、二阶差分与L2惩罚、训练期背景项；α只按验证段选择。测试段不参与拟合或缩放，块之间有全响应范围和SG半窗的保护区。
- 有效起点：累计绝对响应质量5%；峰值：绝对值最大；质心：绝对质量加权；宽度：5–95%。它们不等于严格反应起点、FWHM或化学收率。
- 事件法：正负极值、训练期MAD/幅值阈值；源和目标均只有一条候选边才接受。竞争候选保留歧义，孤立无响应保留未匹配。
- V2事件接口：公共源候选峰附近的可靠窗口中位时滞筛选目标；MAE使用V2时滞，不用观测峰差替代。起点和质心估计为空。
- 核筛选：按验证、正则化敏感性、单主导极性和尾部质量检查核，再限制候选目标峰。事件时刻仍来自可见事件；核描述符与事件描述符分开保存。
- 真实研究：20条炉次—气体序列分开建模，99次相位负对照并做20重BH校正；这些诊断依赖零假设，不能直接解释为因果置信度。
- 常数背景对照属于观察到线性外推失败后的事后敏感性分析，未用于根据测试集重新选择主方法。未执行流量模型、炉际/炉内效应估计。

## 结果与字段

| File or directory | Content |
|---|---|
| `results/REPORT.md` | 当前研究报告；合成与真实结果分章节报告 |
| `results/prior_review.json` | 18场景清单、11场景不重跑、7场景重放与冻结引用指纹 |
| `results/prior_nonideal_kernel.json` | 原7类连续非理想场景的核重放；无事件真值时PR/F1为空 |
| `results/development/` | 6种子开发对照，不与30种子评价混合 |
| `results/synthetic/` | 逐种子指标、事件真值、匹配表、汇总与核描述符误差 |
| `results/real/` | 十炉候选事件、核验证、SG/背景敏感性和相位负对照 |
| `results/examples/` | 每类合成场景首个评价种子的曲线与核，用于复核绘图 |
| `results/real_examples/` | 20条真实序列的输入、响应与预测，用于独立研究图，不是第三问接口 |
| `results/release_inventory.json`、`release_identity.json`、`release_npz_comparison.csv` | Release快照、资产/字段校验与处理接口证据 |
| `results/attachment1_identity.json` | 本地原始文件指纹及指定新原始附件1的待核验状态 |
| `results/figure_captions.csv`、`font_validation.json`、`layout_preview.png` | 可编辑图注、真实Times New Roman/PDF嵌入记录、唯一版宽预览 |
| `results/effective_config.json`、`stage_status.json`、`provenance.json` | 实际参数、阶段签名、代码/依赖/输入来源与全部输出SHA256 |

事件表的联合键为 `(heat_id, gas或scenario, method, event_id)`；`event_id`是候选编号，
`independence_status=not_established`不表示已确认独立事件。
时间字段单位均为秒；真实 `source_flow` 为候选源区间内流量中位数，单位m³/s。
`status`可取matched、ambiguous、unmatched；非matched行的gas与delay字段为空，不可填零。
`delay_*`是可见候选事件的时间差；`estimator_*`才是对应方法用于MAE的估计，尤其V2只提供peak。
`candidate_target_ids`保留竞争目标列表；`confidence`是启发式分数，不是概率。
`kernel_validation_screen_pass`不等于独立物理标签，`validated_physical_match`对本轮真实数据均为false。

F1按真实事件TP/FP/FN汇总；MAE条件于正确匹配，必须与召回一起读。
`test_skill`是相对训练期背景外推的MSE改善；`test_r2`相对测试段真实均值，不可互换。
真实数据没有独立响应真值，不能用上述合成准确率替代真实准确率。
所有图内英文、数字均采用共享真实Times New Roman配置；PDF嵌入检查失败时停止。
图注不烧录进图像，SVG保留可编辑文字并要求查看端安装Times New Roman。
