> 历史图集：当前唯一有效图集为[publication_v4](../publication_v4/README.md)，唯一正式正文见[当前章节](../../reports/问题一_中文论文_V4_机制与结果图优化版.md)。本目录保留研究来源及原始审查记录。

# 第一问出版图集

本轮以英文配图版原19图为对象，重构为7幅正文主图、十炉独立曲面及7幅研究附图；旧图完整保留。推荐配套阅读[七节新正文](../../reports/问题一_中文论文_V3_研究深化与图表重构版.md)和[完整图表附录](../../reports/问题一_V3_图表与研究附录.md)。

|正文图|研究作用|文件|
|---|---|---|
|5-1|分组谱形、秩剩余能量、暗亮通道结构|main/F01_data_structure.*|
|5-2|同一真实异常的局部参考、评分及修复机制|main/F02_model_mechanism.*|
|5-3|四例尖峰与持续谱线保留的对照证据|main/F03_local_evidence.*|
|5-4|连续方法比较、实际单元及凸组合性质|main/F04_continuous_reconstruction.*|
|5-5|有限策略、检测收益、修复代价和区域表现|main/F05_joint_validation.*|
|5-6|十炉实际修复与回退子集、局部差值|main/F06_actual_repair_results.*|
|5-7|统一记录尺度的全部十炉连续曲面|main/F07_ten_heat_surfaces.*|

附图A01—A10为十炉独立曲面；B01为求解稳定性，B02为宽峰回退案例，B03四图为两炉原始/中心化SVD，B04为持续通道时间轨迹及亮度关系。PNG为600 dpi，SVG保留英文文字，PDF嵌入字体；热图与曲面高分辨率栅格化，标注仍为矢量。中文图注在[figure_captions.md](figure_captions.md)，论文中置于对应论证位置。

[19图实际审查](figure_audit_19.md)、[旧图映射及来源](figure_mapping.csv)、[图片哈希](asset_manifest.csv)、[显示采样审计](surface_sampling_audit.csv)、[绘图参数](render_config.json)、[尺寸核验](render_dimensions.csv)。review/保存按A4纸、6.8英寸内容宽缩放的正文全页检查图，仅供排版检查，不用于论文插图。所有主图及附图均已打开渲染检查；T1表头已改为两行以避免相邻列标题拥挤，F02和F03的证据注释已避开/遮底保障可读性。

十炉展示网格均保留实际全局峰值和已核实案例，避免规则抽样漏掉窄峰；曲面仍为冻结连续函数的采样展示，原矩阵及参数没有改变。统一0–65535 counts记录尺度不用于表达超量程潜在强度。条件性估计及质量标记仍见正式接口。

复现数据、协议和脚本：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/analyze.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/finite_checks.py
.venv/bin/python problem1/final/supplement_v1/src/report.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/publication_figures.py --stage priority
OPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/publication_figures.py --stage remaining
.venv/bin/python problem1/final/tables_v3/build_tables.py
.venv/bin/python problem1/final/supplement_v1/src/manuscript.py
.venv/bin/python problem1/final/supplement_v1/src/validate.py
```

两阶段绘图须依次执行，避免并行写图片清单。输入包括原始附件、已发布正式接口、历史实验CSV及本轮真实补充计算。需先按仓库交接说明取得正式数据包，supplement_v1/data/是可复现的本地派生数组，不作为正式数据Release。外部参考册只提供图表布局参考；不分发第三方图像，也不将其数值作为本题结果。
