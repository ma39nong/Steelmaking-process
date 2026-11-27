# 第一问当前正式图集

本目录是唯一当前取图入口，配套[唯一当前正式正文](../../reports/问题一_中文论文_V4_机制与结果图优化版.md)和[统一图表附录](../../reports/问题一_V3_图表与研究附录.md)。保留既有目录名称，不再增加版本目录或正文副本。

|图号|科学问题与结构|图片|
|---|---|---|
|5-1|分组谱形、低秩及暗亮帧结构|[F01](main/F01_data_structure.png)|
|5-2|整体算法流程；质量分类、三证据、修复及截顶独立分支|[F02](main/F02_model_mechanism.png)|
|5-3|2×2双向证据：第一行修复孤立异常，第二行保留持续谱峰|[F03](main/F03_local_evidence.png)|
|5-4|1×2连续方法对照：时间、波长降采样，共享图例|[F04](main/F04_continuous_reconstruction.png)|
|5-5|共同协议下的检测与恢复检验|[F05](main/F05_joint_validation.png)|
|5-6|全部1043点的实际改动分布|[F06](main/F06_actual_repair_results.png)|
|5-7|十炉完整连续光谱总览|[F07](main/F07_ten_heat_surfaces.png)|

所有正式图片均保存PNG/SVG/PDF。main/为七幅正文图，appendix/为十九幅附录图：十炉独立曲面、七幅已有研究附图，以及本轮从正文移出的[A11定位和其他三例尖峰](appendix/A11_case_locator_and_spikes.png)、[A12双线性单元](appendix/A12_bilinear_cell.png)。已有十七幅附图从历史来源原字节集中到这里；统一附录索引只引用本目录。

## 三张核心图的简化

F02采用横向输入、三个并行证据分支、纵向修复步骤及横向连续输出；框内只有必要英文方法名称。疑似截顶分为记录尺度保留和条件性潜在估计，记录量程路径进入正式矩阵。没有单个案例数值或大段公式。

F03为两行两列：上行为6B03205的14.5 s、512.659 nm孤立异常，下行为6C02108的304.0 s、766.413 nm持续谱峰；左列为时间，右列为波长。原始、修复、邻域参考统一视觉编码，纵轴按每个面板实际数据分别设置。灰色叉号是目标节点的实际邻域参考，非整段参考曲线。

F04仅保留时间和波长降采样两面板，共用图例，所有十五个RMSE值及采样间隔沿用原独立实验。原插值单元移到A12，原公式子图删除；凸组合、非负、连续与无额外过冲性质由正文5.5推导。

## 真实版宽检查与复现

按6.8英寸（约173 mm）正文版宽制作，基本文字8.2 pt；逐张打开实际渲染检查了框内换行、公共图例、坐标、面板间距及留白。F02框内文字还通过实际渲染包围盒核对，避免溢出框体。600 dpi仅作为输出规格，不代替版面审查。

- [F02全页检查](review/F02_model_mechanism_page.png)
- [F03全页检查](review/F03_local_evidence_page.png)
- [F04全页检查](review/F04_continuous_reconstruction_page.png)
- [简化与归档说明](../../archive/reports/三张核心图_重设计说明.md)
- [图注](figure_captions.md)、[来源哈希](source_manifest.csv)、[输出哈希](asset_manifest.csv)、[数值及文件核验](validation.json)

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/final/figures/publication_v4/src/redesign.py
.venv/bin/python problem1/final/figures/publication_v4/src/package.py
```

simplify.py是F02/F03/F04及A11/A12的唯一当前绘图实现；package.py只核验当前正文、数据与图片，不从旧稿重新生成正文。当前26张图统一采用真实Times New Roman，所有子图标题居中；redesign.py复用已有绘图数据并调用当前简化实现，历史图像保持归档。历史publication_v3只作研究来源，模型、1043点修复结果、实验指标及Release不变。

## 论文图题与字体

参照仓库内2022—2025年B题优秀论文，整图题在图片下方以正文排版，图内保留英文标签。七条英文图题见[english_captions.md](english_captions.md)，实际样本核对与字体规范见[说明](../../reports/优秀论文图表样式核对与采用规范.md)。带外置图题的A4检查页在`review/*_caption_page.png`；正式正文仍使用`main/`原图。

字体安装：需本机安装真实Times New Roman（Regular/Bold/Italic/Bold Italic），脚本使用`findfont(..., fallback_to_default=False)`核实，缺失时停止，避免静默替换字体。字体文件不随仓库分发。

本次字体与标题统一已完成402项核验：26张PDF均嵌入TimesNewRoman字体，数值源文件和模型保持冻结，逐图绘图数据在样式处理前后相同。全部主图的带图题版宽预览见`review/F*_caption_page.png`。
