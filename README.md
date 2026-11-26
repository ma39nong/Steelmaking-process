# Steelmaking-process · 2026 APMCM B题

当前研究分支为 **`problem1`**。第一问已完成数据分析、实际异常点修复、连续光谱建模和中文正式章节，最终采用 **自适应多证据异常识别＋秩4 Huber鲁棒低秩补全＋观测保真及局部非负支撑＋双线性连续重构**。第二、三问尚未在本分支实现。

## 从这里开始

|你要找什么|正式入口|
|---|---|
|第一问完整中文正文，直接纳入整篇论文第5章|[问题一_中文论文_V2_正式章节版.md](problem1/final/reports/问题一_中文论文_V2_正式章节版.md)|
|正文每幅图、每张表的位置及数据来源|[问题一_V2_正文图表清单.md](problem1/final/reports/问题一_V2_正文图表清单.md)|
|摘要素材、假设、统一符号、文献和附录入口|[问题一_V2_全文整合材料.md](problem1/final/reports/问题一_V2_全文整合材料.md)|
|最终模型、代码、运行顺序、数据字段|[最终研究说明](problem1/final/README.md)|
|给第二、三问同事的数据选择与读取方法|[数据交接与后续问题使用说明](problem1/final/reports/数据交接与后续问题使用说明.md)|
|处理后数据包与下载/发布方式|[数据交付目录](problem1/final/delivery/README.md)|
|原始附件名称、炉次、范围、SHA256|[data_manifest.csv](data_manifest.csv)|

## 最后定下来的是什么

最终计算结果保存在 **`problem1/final/`**。这与“中文论文V2”的文稿编号不同：文稿V2是正式章节重写，不是将算法改回历史V2。

1. 用真实时间和非均匀波长坐标构造去当前点的双向预测，按局部波动尺度校准残差。
2. 融合时间偏离、波长偏离、秩4背景残差、谱线持续性与双向参考一致性；综合阈值为4，局部8连通区域跨度不超过3帧、7通道。
3. 将有依据的局部统计异常设为不可见，采用固定秩4、Huber尺度0.003的低秩模型估计；保留其余记录，并对修复位置施加局部非负支撑约束及既有数值稳定性处理。
4. 疑似顶码65535在正式记录曲面中保留量程信息，条件性潜在强度另用局部核范数下界约束估计。
5. 在修复矩阵上建立双线性函数，输出十炉连续三维光谱。

十炉共27,816,683个观测，**实际修复1043点（0.00375%）**；14,805个疑似顶码独立标记，整帧删除0。正式流程在固定污染确认中污染区RMSE由6510.69降至2359.80 counts，降低63.76%；全部观测误差降低56.96%。[逐炉结果](problem1/final/results/final_summary.csv)、[逐点实际修复](problem1/final/results/actual_repairs.csv)、[确认明细](problem1/final/results/exact_final_pipeline_confirmation.csv)可逐项追溯。

### 模型具体写在什么地方

|模块|正式实现|论文位置|
|---|---|---|
|自适应尺度、评分、持续性、连通范围|[detectors.py](problem1/final/src/detectors.py)|第5.2节，式(5-3)—(5-8)|
|原网格修复集合与逐点复核|[native_masks.py](problem1/final/src/native_masks.py)|第5.2.4节|
|固定秩Huber目标、梯度、回溯及支撑投影|[model3.py](problem1/v3/src/model3.py)|第5.3—5.4节，式(5-10)—(5-19)|
|600/1200步输出稳定性处理|[stable_completion.py](problem1/v3/src/stable_completion.py)|第5.4节|
|条件性截顶约束与ADMM|[convex_censor.py](problem1/v2/src/convex_censor.py)、[conditional_final.py](problem1/final/src/conditional_final.py)|第5.3.4及5.4节|
|正式逐炉计算与输出|[fit_native.py](problem1/final/src/fit_native.py)、[export_results.py](problem1/final/src/export_results.py)|第5.7节|
|连续查询接口|[interface.py](problem1/v3/src/interface.py)|第5.5节，式(5-21)—(5-23)|
|参数、开发选型及确认协议|[protocol.json](problem1/final/config/protocol.json)、[SELECTION_LOCKED.json](problem1/final/config/SELECTION_LOCKED.json)|第5.6节|

最终实现复用历史目录中的稳定模块；代码位于`v2/`或`v3/`并不表示该目录内旧的预处理策略仍是当前正式策略。

## 图片在哪，分别用于什么分析

|分析内容|图像或清单|正文位置|
|---|---|---|
|光谱低秩性质|[奇异值能量](problem1/final/figures/main/低秩能量贡献.png)|5.1|
|整体方法|[技术路线](problem1/final/figures/main/第一问技术路线.png)|5.2|
|线性、PCHIP、三次样条、张量B样条选型|[插值对照](problem1/v2/figures/interpolation_test.png)|5.5—5.6|
|同一污染协议的检测与修复收益|[端到端比较](problem1/final/figures/main/统一端到端异常处理比较.png)|5.6|
|目标函数及迭代变化|[收敛图](problem1/final/figures/main/最终求解收敛.png)|5.6|
|实际修复数量与四个典型尖峰|[逐炉统计](problem1/final/figures/main/十炉实际修复统计.png)、[四点对照](problem1/final/figures/main/四个高证据尖峰实际修复.png)|5.7|
|十炉正式曲面|[三维总览](problem1/final/figures/main/十炉最终连续光谱总览.png)、[独立曲面及附录](problem1/final/figures/appendix/)|5.7|

[正文图表清单](problem1/final/reports/问题一_V2_正文图表清单.md)逐一登记19幅图、10张表和对应数据文件。各炉`_全谱处理与差值.png`展示原始/处理热图及差值，`_异常修复曲线与局部热图.png`展示局部依据，`_帧级快变诊断.png`分析宏观变化。图中三维高度来自连续函数；数值共享请使用NPZ，而非从图片读取。

## V1、V2、V3、V3.1保留了什么

|研究阶段|主要研究作用|目录及文字报告|当前定位|
|---|---|---|---|
|V1|双向MAD、三类模型、最初掩码补全与连续光谱|[说明](problem1/README.md)、[报告](problem1/reports/问题一_数据分析与模型建立及实验结果.md)|历史基线|
|V2|低秩消融、独立插值、饱和约束及严格遮蔽验证|[说明](problem1/v2/README.md)、[报告](problem1/v2/reports/问题一_V2模型深化与严格验证.md)|实验与数学依据|
|V3|保真、非负支撑、记录/条件性饱和双输出、双线性及数据诊断|[说明](problem1/v3/README.md)、[报告](problem1/v3/reports/问题一_最终模型与实验结果.md)|最终核心模块来源|
|V3.1|点级统计证据、端到端污染恢复与保守策略复核|[说明](problem1/v3_1/README.md)、[报告](problem1/v3_1/reports/数据预处理_点级异常最终验收.md)|检测改进依据|
|最终研究|自适应多证据＋实际1043点修复＋已定型重构模块|[final说明](problem1/final/README.md)、[研究初稿](problem1/final/reports/问题一_中文论文初稿.md)|**正式算法和数据**|
|正式章节写作|按整篇竞赛论文第5章重构七节正文|[正式章节版](problem1/final/reports/问题一_中文论文_V2_正式章节版.md)|**当前论文入口**|

历史报告中的“最终”“P0保留”“修复0点”等均属于各自研究阶段。后续同事选择数据时统一以`problem1/final/data/*_interface.npz`和当前正式章节为准。历史实验协议不同，方法比较采用正文明确登记的同协议结果。

## 数据在哪，别人如何取得

|数据|本机/仓库相对位置|Git是否包含|用途|
|---|---|---|---|
|原始十炉XLSB|`Data/`；只读副本`problem1/data/raw/`|不包含|按官方附件或团队原始来源下载，再按根manifest校验|
|原始解析缓存|`problem1/data/processed/`|不包含|真实轴、原始Y及仅审计的工艺附加列|
|最终数据接口|`problem1/final/data/*_interface.npz`|不包含|十炉原始/修复光谱、真实轴、质量掩码、条件估计|
|修复差值与密集网格|`problem1/final/data/*_difference.npz`、`*_dense.npz`|不包含|修改追溯、连续曲面展示|
|论文、代码、统计CSV、图像、哈希|`problem1/final/reports/`、`src/`、`results/`、`figures/`|包含|阅读、验证与复现|

**只克隆仓库不会得到大数组。** 已发布523.40 MiB处理后数据交付包，含十炉接口、十炉差值及校验清单：[Release下载页](https://github.com/ma39nong/Steelmaking-process/releases/tag/problem1-final-data-v1)、[直接下载ZIP](https://github.com/ma39nong/Steelmaking-process/releases/download/problem1-final-data-v1/problem1-final-data.zip)、[打包与Release说明](problem1/final/delivery/README.md)。下载并在仓库根目录解压后恢复正式数据路径。包内包含十炉原始网格上的处理结果、原始光谱副本、质量标记及函数接口；原始XLSB中的工艺附加列仍由后两问同事从原附件读取。

截图中的 **Releases适合发布这个数据包**，单个附件须小于2 GiB；Packages主要用于npm、容器等软件包注册表，本项目采用Release附件更直接。[GitHub官方Release说明](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases)、[Packages说明](https://docs.github.com/en/packages/working-with-a-github-packages-registry)。正式数据附件已通过上述Release单独发布；提交代码不会自动上传被忽略的数据。

## 给第二、三问

应交付**修复矩阵＋原始矩阵＋时间/波长轴＋异常及截顶质量掩码＋读取接口＋数据说明**。第二问的离线特征分析可使用`reconstructed`，同时保留`raw`做处理前后对照；氧气、烟气等变量按原始附件、实际采样口径对齐。第三问实时预测须从`raw`建立仅使用当前及历史信息的处理流程，本问的离线异常判断、掩码、低秩估计和连续曲面均可能利用未来观测。训练/测试划分和输入可见范围由各任务单独固定。详见[交接说明](problem1/final/reports/数据交接与后续问题使用说明.md)。

## 文稿交付与核对记录

正式章节版完整保存低秩、Huber、约束、数值算法及双线性推导，并补充普通缺失和独立插值选型证据；按七个二级小节组织，复用19幅图、整理10张表。图表位于对应论证位置，主要局限集中在5.6节末尾。摘要素材、假设、符号、参考文献及附录单独保存供全文整合。

文稿提交`77df8a7`已检查数学公式与实现对应、关键数值来源、全部本地链接、公式/图表编号和十炉结果覆盖；原初稿、模型、参数、实验及1043点结果均未修改。具体来源登记在[全文整合材料](problem1/final/reports/问题一_V2_全文整合材料.md)。

## 原始数据与文献入口

[原始数据映射与SHA256](data_manifest.csv)固定十炉版本，6B03205使用已校验的替换附件。第三方全文按既有获取约定管理。[文献索引](references/README.md)、[相关综述](reports/00_相关文献综述.md)、[方法适用性](reports/00_文献方法与本题适用性分析.md)、[原始数学来源](reports/00_数学模型推导与算法来源.md)完整保留。早期[复现命令](problem1/README.md)对应V1；当前正式运行及输出以[final说明](problem1/final/README.md)为准。
