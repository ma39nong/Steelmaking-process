# Problem 2

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
uv run --with 'numpy<2' --with pandas --with pyxlsb --with scipy --with matplotlib `
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
4. 局部时滞只对具有动态激励、显著相关、峰值分离、峰宽和边界条件均合格的窗口返回观测值；不可辨识窗口保存 NaN 和掩膜。
5. 流量—时滞关系使用可靠窗口中心样本进行固定模型、反平方根模型和约束模型比较。样本不足时保留固定模型或不可用状态。
6. 时间映射遵循 `t_measurement = t_source + delay(t_source)`，并保存观测映射、模型映射和单调有效掩膜。

## 结果与接口

- `final/results/process_validation.csv`：十炉字段、Sequence、采样间隔和缺失值核验。
- `final/results/q2_summary.json`：每炉处理摘要和可靠窗口覆盖率。
- `final/results/relationship_summary.json`：四类累计耗氧量—烟气关系。
- `final/results/cross_validation.json`：留一炉次时滞模型验证和阈值敏感性。
- `final/results/reconstruction_evaluation.json`：重构、反积分和全字段质量检查。
- `final/results/raw_injection_validation.json`：原始 XLSB 注入尖峰后的端到端检测结果。
- `final/results/masking_validation.json`：随机点和连续区块遮蔽结果。
- `final/results/calculus_comparison.json`：导数、积分和解析积分对照。
- `final/data/{heat_id}_q2_processed.npz`：后续问题的逐炉接口。

NPZ 字段、维度、单位、NaN 语义和读取示例见 `final/reports/数据字典与接口说明.md`。第一问光谱结果仅通过 `heat_id` 和逐炉时间轴只读连接，第二问不修改第一问文件。

## 数据解释边界

CO/CO₂ 原始字段按题面作为瞬时体积流量处理。统计响应时滞同时包含反应、传输和仪器响应成分，不能直接视为独立测得的管道传输时间。动态时滞的可靠窗口覆盖率较低，后续问题必须读取有效掩膜和置信度，不能把无效区间填零或把模型时滞当作实测标签。
