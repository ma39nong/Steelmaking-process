# Problem 2

当前定稿入口为 `python -m problem2.src.final_delivery`，配置见 [final_analysis.json](configs/final_analysis.json)。[中文正式章节](final/reports/问题二_中文论文_正式章节版.md)包含九张统计表、六组核心图；[最终验收](final/reports/第二问最终验收报告.md)给出十炉接口检查。时滞算法开发暂停，正式V2保留，尚未验证流量—时滞物理函数。

本目录实现 APMCM 2026 B 题第二问的工艺信号分析。输入是十炉原始 XLSB；输出是连续工艺信号、氧枪流量、累计烟气排放、可靠性筛选后的局部响应时滞、时间映射和可供后续问题调用的 NPZ 接口。原始 XLSB 只读读取，所有结果写入 `problem2/final/`。

## 目录

```text
problem2/
├── configs/                 # 固定参数
├── src/                     # 数据读取、重构、流量、关系及时滞算法
├── tests/                   # 合成信号与接口回归测试
├── final/
│   ├── data/                # 运行生成的每炉 NPZ（默认不提交大文件）
│   ├── results/             # JSON/CSV 实验结果
│   ├── figures/             # 真实数据和模拟验证图
│   ├── tables/              # 表格说明与导出位置
│   └── reports/             # 建模、接口、复现和图表说明
└── README.md
```

## 一键复现

在仓库根目录执行：

```powershell
uv run --with 'numpy<2' --with pandas --with pyxlsb --with scipy --with matplotlib --with pillow `
  python -m problem2.src.run_all --data-dir G:\亚太杯\Data
```

该命令依次完成原始字段与时间轴核验、十炉连续信号处理、关系分析、留一炉次验证、微积分一致性、原始数据异常注入、遮蔽实验和敏感性图。结果写入 `problem2/final/results/`，每炉接口写入 `problem2/final/data/`。

分模块运行示例：

```powershell
uv run --with pandas --with pyxlsb python -m problem2.src.data_validation `
  --data-dir G:\亚太杯\Data --output-dir problem2\final\results

uv run --with 'numpy<2' --with pandas --with pyxlsb --with scipy `
  python -m problem2.src.q2_pipeline_v2 `
  --data-dir G:\亚太杯\Data --output-dir problem2\final\data `
  --summary-dir problem2\final\results

uv run --with 'numpy<2' python -m problem2.src.v2_validation `
  --data-dir problem2\final\data --output problem2\final\results\cross_validation.json

uv run --with 'numpy<2' --with scipy --with pytest pytest -q problem2/tests
```

## 当前模型

1. Hampel 稳健异常识别在局部尺度为零时使用平台/邻域回退规则；PCHIP 在保留点上重构三类信号。
2. 氧枪流量采用累计耗氧量连续曲线的 Savitzky–Golay 近似导数，并截断负值；同时保存由流量梯形积分得到的累计量用于一致性检查。
3. CO 和 CO₂ 的累计排放由清洗后的瞬时体积流量积分得到。
4. 局部时滞只对动态激励、相关阈值、峰值分离和边界条件均合格的窗口返回观测值；峰宽与峰显著性仅记录为诊断量，不参与筛选。不可辨识窗口保存 NaN 和掩膜。
5. 正式接口保留基于可靠窗口的既有固定与约束反平方根探索结果；不足五个合格窗口时模型不可用。24个重叠窗口不构成已确认独立事件，不重新开展流量—时滞回归。
6. 时间映射遵循 `t_measurement = t_source + delay(t_source)`，并保存观测映射、模型映射和单调有效掩膜。

## 结果与接口

- `final/results/process_validation.csv`：十炉字段、Sequence、采样间隔和缺失值核验。
- `final/results/q2_summary.json`：每炉处理摘要和可靠窗口覆盖率。
- `final/relationships/`：当前四类关系的原值/增量相关、等耗氧进度阶段、趋势与留一炉次结果、流量复核及接口验收。
- `final/results/relationship_summary.json`：旧版四类关系诊断（包含分位数阶段），与当前等耗氧进度阶段区别保存。
- `final/results/cross_validation.json`：留一炉次时滞模型验证和阈值敏感性。
- `final/results/reconstruction_evaluation.json`：重构、反积分和全字段质量检查。
- `final/results/raw_injection_validation.json`：原始 XLSB 注入尖峰后的端到端检测结果。
- `final/results/masking_validation.json`：随机点和连续区块遮蔽结果。
- `final/results/calculus_comparison.json`：导数、积分和解析积分对照。
- `final/data/{heat_id}_q2_processed.npz`：后续问题的逐炉接口。

NPZ 字段、维度、单位、NaN 语义和读取示例见 `final/reports/数据字典与接口说明.md`。第一问光谱结果仅通过 `heat_id` 和逐炉时间轴只读连接，第二问不修改第一问文件。

## 数据解释边界

CO/CO₂ 原始字段按题面作为瞬时体积流量处理。统计响应时滞同时包含反应、传输和仪器响应成分，不能直接视为独立测得的管道传输时间。动态时滞的可靠窗口覆盖率较低，后续问题必须读取有效掩膜和置信度，不能把无效区间填零或把模型时滞当作实测标签。

## 第一轮时滞专项审查

审查基于 `cea8248`，不修改正式时滞模型和第三问接口。独立结果位于
`final/delay_audit/`；默认实验配置为 `configs/delay_audit_config.json`，
预处理和时滞参数显式读取 `configs/q2_config.json`。

完整复现（仓库根目录，已安装 README 上述依赖的 Python 环境）：

```bash
python -m problem2.src.delay_audit --data-dir Data
python -m pytest -q problem2/tests
```

本机已有环境可将 `python` 换成 `.venv/bin/python`。完整运行会重新生成审查目录中的
基线 NPZ、实验、图表与报告第6节，覆盖当前输出，不替换 `final/data/` 或原有正式结果。
正式源文件必须与 `cea8248` 一致，否则拒绝把本次运行称为原基线复现。

分阶段运行使用同一配置和输出目录：

```bash
python -m problem2.src.delay_audit --stage real
python -m problem2.src.delay_audit --stage synthetic
python -m problem2.src.delay_audit --stage figures
python -m problem2.src.delay_audit --stage report
```

`real` 生成十炉重放、独立拒绝原因、窗口重叠分组、双向/全局对照、窗口敏感性及已有流量模型验证；
`synthetic` 运行30个固定种子的18种场景；`figures` 需要前两阶段结果；`report` 更新现有审查报告。
更改配置后必须完整运行，禁止混用不同配置的分阶段结果。正式原入口的配置接入不属于本轮修改。

| 审查目录文件 | 内容 |
|---|---|
| `effective_config.json`、`provenance.json` | 有效配置、输入及源代码来源、依赖版本、输出SHA256 |
| `input_audit.json` | 十炉哈希、时间轴、工艺字段、数据质量与积分总量检查 |
| `baseline/data/`、`baseline/results/` | 重放NPZ、历史JSON比较及逐数组SHA256；NPZ按既有策略不提交Git，可重新生成 |
| `window_diagnostics.csv`、`real_summary.json` | 所有4482个炉次—气体窗口的独立检查和顺序筛选计数 |
| `accepted_windows.csv` | 24个可靠窗口、流量/阶段、三种重叠分组、跨气体数据共享；事件身份保留不确定 |
| `signed_comparison.csv`、`global_correlations.json` | 相同有效源中心的双向/非负搜索，以及全局描述性互相关 |
| `window_sensitivity.csv` | 十炉81/121/161点窗口对照；不同窗口的中心范围不同 |
| `synthetic_runs.csv`、`synthetic_windows.csv` | 每种子/场景指标与逐窗口估计、真实起点/峰值/重心 |
| `synthetic_summary.json`、`synthetic_peak_calibration.csv` | 合成误差、拒绝率、方向错误及峰宽/显著性诊断 |
| `flow_model_evidence.json`、`flow_model_cross_validation.csv`、`flow_stage_coverage.csv` | 既有模型训练误差、固定/反平方根留一炉次比较和流量/阶段覆盖 |
| `figures/`、`figure_captions.csv` | 24张逐窗口检查图和2张总览，PNG/SVG/PDF；图注独立可编辑 |
| `font_validation.json`、`layout_preview.png` | 四种真实Times New Roman字形与PDF嵌入检查；唯一版宽检查预览 |

报告见 `final/reports/独立数学审查与修复报告.md` 第6节。
窗口、重叠组均不能直接作为独立事件；合成窗口也高度重叠，30个种子才是重复实验单位。
同事历史NPZ未包含在当前检出中，故逐数组比较属于配置重放对未改动基线代码的验证，
不表示已经直接核对其历史数组。旧论文章节的峰宽筛选描述暂未改写，应以本轮审查更正为准。

新图要求真实 Times New Roman、fontconfig 的 `fc-match` 和 Poppler 的 `pdffonts`。
共享配置在 `src/figure_fonts.py`；缺失时明确失败，不能替代字体继续导出。
PDF嵌入字体，SVG保留可编辑文字并要求查看设备安装Times New Roman。

## 响应展宽与事件匹配独立研究

在 `9b8513a` 基线之后，新增 [独立研究入口](research/response_broadening/README.md) 与
[当前研究报告](research/response_broadening/results/REPORT.md)。新研究比较V2事件接口、保留歧义的一对一事件匹配及正则化FIR响应核筛选，
合成与十炉真实证据分开保存。正式V2、第三问接口及前轮审查报告均保持冻结；本轮未启动流量—时滞拟合。

## 当前真实事件验证

以 `dea6637` 为基准的[真实事件验证](research/real_event_validation/results/REPORT.md)及[复现入口](research/real_event_validation/README.md)保存十炉候选、原24窗口复核、重叠与尺度检查。用户已确认中文附件1只是现有英文原始文件的改名；本轮文件指纹核验关闭了前轮报告中的数据身份待确认项。历史报告仍保留当时状态，正式模型与NPZ接口保持冻结。

当前正式复现步骤与输入要求见 [实验配置与复现说明](final/reports/实验配置与复现说明.md)。历史研究脚本按各自冻结提交复现；其旧待办不覆盖本轮数据身份确认和时滞开发暂停决定。

## 最终论文审查与冻结

建模与实验结束，当前只维护正式论文和离线交接。九项审查见 [最终论文审查报告](final/reports/第二问最终论文审查.md)。正文在原八表基础上补入流量分层候选时滞表，现为九表、六组核心图；窗口统计不代表独立物理响应。正式生成后运行 `python -m problem2.src.final_review` 核对数值、图表、输入输出指纹并更新 `final/relationships/paper_freeze_manifest.json`。第三问可使用未作时滞修正的原时间轴作为基线。
