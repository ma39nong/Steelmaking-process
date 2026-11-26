# Problem 2

第二问的代码和正式报告位于本目录。原始 XLSB 文件只读存放在本地资料目录，代码不会覆盖它们。

## 当前入口

```powershell
uv run --with pandas --with pyxlsb python -m problem2.src.data_validation `
  --data-dir G:\亚太杯\Data `
  --output-dir problem2\final\results
```

`data_validation` 逐炉读取唯一工作表，按表头识别累计耗氧量、CO₂ 和 CO，验证 2051 列、连续 Sequence、0.5 s 时间步、缺失值和负值，并写出 JSON 与 CSV 审计结果。

## 基线重构与评价

```powershell
uv run --with 'numpy<2' --with pandas --with pyxlsb --with scipy `
  python -m problem2.src.q2_pipeline `
  --data-dir G:\亚太杯\Data `
  --output-dir problem2\final\data `
  --summary problem2\final\results\q2_summary.json

uv run --with 'numpy<2' --with scipy python -m problem2.src.evaluation
uv run --with 'numpy<2' --with matplotlib python -m problem2.src.visualization
```

基线输出已在本地十炉数据上运行。NPZ 和 PNG 按 `.gitignore` 保留在实验工作区，JSON、CSV、代码和论文材料进入版本库。

## 数据接口

后续重构模块将从原始工艺列生成 `problem2/final/data/{heat_id}_q2_processed.npz`。接口至少保留原始 Sequence、相对时间、原始和清洗后的三类工艺变量、氧枪流量、总排放量、CO/CO₂ 时滞、时间映射和元数据。第一问光谱接口通过同一 `heat_id` 只读引用。

## 已知数据核验

十炉均为一个工作表、2051 列、Sequence 从 1 连续递增，时间步均为 0.5 s，工艺三列没有缺失值或负值。表头存在英文、中文和缩写差异，代码按实际表头识别。单位按题面采用累计氧量 m³、CO₂/CO 瞬时体积流量 m³/s。
