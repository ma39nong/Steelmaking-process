# Steelmaking-process · APMCM B题

本分支完成**第一问：火焰光谱异常修复与连续曲面重构**。新同事按下面的入口阅读即可，不需要先浏览各历史版本。

## 从这里开始

|你要做什么|直接打开|
|---|---|
|了解方法与结果、撰写论文|[当前论文正文](problem1/final/reports/问题一_中文论文_V4_机制与结果图优化版.md)|
|找论文图片、图注和数据来源|[当前图集](problem1/final/figures/publication_v4/README.md)|
|找可编辑表格和完整附录|[表格](problem1/final/tables_v3/README.md) · [图表附录](problem1/final/reports/问题一_V3_图表与研究附录.md)|
|取得处理后数据，接着做第二、三问|[数据交接说明](problem1/final/reports/数据交接与后续问题使用说明.md) · [Release下载](https://github.com/ma39nong/Steelmaking-process/releases/tag/problem1-final-data-v1)|
|查模型代码、参数和复现方法|[正式交付说明](problem1/final/README.md)|

## 最终采用的模型

**自适应多证据异常识别 → 秩4 Huber鲁棒低秩补全 → 观测保真与局部非负支撑 → 双线性连续重构。**

十炉共27,816,683个观测，实际修复1,043点；14,805个疑似65535截顶位置独立标记。论文第5.2—5.4节介绍异常检测与修复，第5.5节给出连续函数，第5.6节验证，第5.7节展示十炉结果。

## 数据怎么拿

- **原始数据**：使用题目附件，按[data_manifest.csv](data_manifest.csv)核对版本和SHA256。
- **处理后数据**：从[Release](https://github.com/ma39nong/Steelmaking-process/releases/tag/problem1-final-data-v1)下载数据包，在仓库根目录解压，恢复`problem1/final/data/`。
- **第二问**：可使用修复光谱，同时保留原始光谱与质量掩码作对照。
- **第三问**：实时预测应从原始数据构造因果输入；第一问离线重构使用前后时刻信息。

只克隆代码仓库不会获得大型数据数组。图像用于展示，数值分析使用数据包中的矩阵和坐标。

<details>
<summary>历史研究、文献与复现实验（需要时再展开）</summary>

- [研究档案与复现说明](problem1/final/reports/研究档案与复现说明.md)：V1—V3.1索引、模型实现位置与完整运行流程。
- [文献索引](references/README.md)：原文获取状态、书目信息及中文笔记。
- [补充研究](problem1/final/supplement_v1/README.md)：持续通道、物理候选、SVD与修复可靠性。

历史代码、数据和报告保持原位。当前阅读与取数统一使用上方入口。

</details>
