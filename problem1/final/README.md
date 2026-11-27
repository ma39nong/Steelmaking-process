# 第一问正式交付

本目录是**最终算法、结果和当前论文**的交付入口。首次阅读从正文开始；只有复现实验时才需要进入历史目录。

## 当前交付物

|内容|位置|
|---|---|
|唯一推荐正文|[七节中文论文](reports/问题一_中文论文_V4_机制与结果图优化版.md)|
|论文图片|[当前正式图集](figures/publication_v4/README.md)：7幅正文图、附图与完整研究索引；英文标注，PNG/SVG/PDF|
|论文表格|[tables_v3](tables_v3/README.md)：五组正文表与完整附表，Markdown/CSV/LaTeX|
|完整图表附录|[图表与研究附录](reports/问题一_V3_图表与研究附录.md)|
|全文整合素材|[假设、符号与参考文献](reports/问题一_V2_全文整合材料.md)|
|处理后数据|[Release下载](https://github.com/ma39nong/Steelmaking-process/releases/tag/problem1-final-data-v1) · [数据交接说明](reports/数据交接与后续问题使用说明.md)|

旧版正文与四套旧图集已集中到[archive/](archive/README.md)，正式`figures/`只保留当前图集。当前写作、翻译和排版统一取上表文件。

## 模型与结果

最终流程为**自适应多证据检测＋秩4 Huber低秩补全＋观测保真/局部非负支撑＋双线性连续函数**。疑似截顶另存质量标记及条件性估计，正式记录曲面保留量程信息。

- [逐炉修复统计](results/final_summary.csv)：十炉实际修复1,043点。
- [逐点原始值与修复值](results/actual_repairs.csv)。
- [同协议实验比较](results/unified_comparison.csv)与[发布流程确认](results/exact_final_pipeline_confirmation.csv)。

## 代码从哪里看

|作用|实现入口|
|---|---|
|异常评分与局部掩码|[detectors.py](src/detectors.py) · [native_masks.py](src/native_masks.py)|
|低秩拟合、保真与支撑约束|[model3.py](../v3/src/model3.py) · [稳定性处理](../v3/src/stable_completion.py)|
|条件性饱和估计|[conditional_final.py](src/conditional_final.py)|
|逐炉求解与结果导出|[fit_native.py](src/fit_native.py) · [export_results.py](src/export_results.py)|
|连续函数查询|[interface.py](../v3/src/interface.py)|
|参数与选择锁|[protocol.json](config/protocol.json) · [SELECTION_LOCKED.json](config/SELECTION_LOCKED.json)|

历史目录中的模块是最终实现的依赖，不能按目录版本号判断当前采用的预处理策略。

## 数据怎么用

数据包解压后，每炉`data/HEAT_interface.npz`包含真实时间、波长、原始/修复矩阵及质量掩码。完整字段与读取示例见[交接说明](reports/数据交接与后续问题使用说明.md)。第二问可使用离线修复矩阵；第三问实时预测需从原始记录建立因果处理。

需要重跑实验、追溯旧结果或核对文献时，再打开[研究档案与复现说明](docs/研究与维护说明.md#研究档案与复现说明)。请保留已有数据、参数锁及实验结果，在独立目录复现。
