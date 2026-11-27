# 第一问机制与结果图重设计

本轮针对publication_v3的F02、F06、F07重新组织信息，保留原图集。沿用b837d0d的研究数据，并保留b144fe8的符号及坐标修正。算法、参数、1043点实际修复及数据Release不变。

[配套七节正文](../../reports/问题一_中文论文_V4_机制与结果图优化版.md) · [设计说明与核对](../../reports/三张核心图_重设计说明.md) · [中文图注](figure_captions.md)

|图|原图问题|新图的论证作用|直接查看|
|---|---|---|---|
|5-2|六幅实验图加公式，缺少机制关系|按有向流程连接原始观测、三类证据、评分/连通门控、修复集合、Huber补全、保真与支撑、连续函数；虚线区分截顶分支|[F02](main/F02_model_mechanism.png)|
|5-6|重复同一尖峰，一像元差值信息少|全部1043个节点的波长—炉次—相对时间图，配绝对改变量与时间分布，区分465回退与578非回退节点|[F06](main/F06_actual_repair_results.png)|
|5-7|三维坐标区内部留白及重复标注挤占曲面|扩大投影、正交视角、共享单位及色标，保留十炉实际终点、完整前缘和全局峰值|[F07](main/F07_ten_heat_surfaces.png)|

[F01](main/F01_data_structure.png)、[F03](main/F03_local_evidence.png)、[F04](main/F04_continuous_reconstruction.png)、[F05](main/F05_joint_validation.png)与当前publication_v3逐字节相同，仅复制至独立完整七图目录。F01暗亮差异、F03案例组织、F04几何表达可作为后续排版选择，本轮不重画。原[17幅附录图](../publication_v3/appendix/)及[完整图表附录](../../reports/问题一_V3_图表与研究附录.md)继续适用。

所有图片保留600 dpi PNG、英文可编辑SVG及PDF。review/提供A4纸上6.8英寸最大内容宽的检查页；PNG预览用于审图，正文可优先使用PDF/SVG。图注、正文引用及分析已随三图同步更新。

## 数据与表达边界

F02的局部记录、残差得分与修复值取自正式接口和实际评分；流程关系为模型机制表达，未生成新实验。图中顶码分支表示另一观测机制，不将该普通尖峰节点解释为饱和。

F06逐点展示实际修复集合，每个节点绘制一次；小改变量在共同线性色标中接近白色，其幅度分布另由经验CDF呈现。圆/菱形区分数值分支，回退子集不与修复总数相加。纵向带内时间t/T仅用于比较分布，连续模型仍以实际秒坐标查询。分布图并非新增检测或修复精度实验。

F07保持统一0–65535记录尺度，显示网格与正式数据分离。放大投影后关闭集合的坐标区裁剪，另核对所有投影顶点都位于整图范围内，防止截掉前缘。全局峰值及案例坐标仍保留；时间起点均为0，前四行省略重复的零时刻与400 nm刻度，末行显示完整参考刻度。

## 复现与检查

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/final/figures/publication_v4/src/redesign.py
.venv/bin/python problem1/final/figures/publication_v4/src/package.py
```

依赖已有正式接口NPZ及supplement_v1/results中的真实评分证据。来源哈希见source_manifest.csv，输出哈希见asset_manifest.csv。repair_atlas_points.csv保存1043点绘图坐标、实际差值及回退标记；surface_sampling_audit.csv保存各炉峰值及投影范围。checks.json与validation.json记录数值、文件及图表一致性检查，旧publication_v3图片逐项与Git对象核对。
