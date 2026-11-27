# 第二问正式统计结果

[中文正式章节](../reports/问题二_中文论文_正式章节版.md)给出九张统计表与六组核心图，[验收报告](../reports/第二问最终验收报告.md)区分计算验收和物理结论。[复现说明](../reports/实验配置与复现说明.md)说明配置、输入和统计口径。

| 正文表 | 逐炉或完整证据 |
|---|---|
| 表1 十炉总量 | heat_summary.csv |
| 表2 原值相关 | level_correlations.csv |
| 表3 增量相关 | increment_correlations.csv |
| 表4 五阶段 | stage_statistics.csv、stage_aggregate.csv |
| 表5 趋势参数 | trend_parameters.csv、pooled_trends.csv |
| 表6 留一炉次 | trend_loho.csv |
| 表7 流量复核 | flow_review.csv |
| 表8 流量分层时滞 | delay_flow_summary.csv、delay_flow_definition.json |
| 表9 接口验收 | interface_checks.csv、interface_field_comparison.csv |

原始计数更新机制未被确证，24个合格窗口没有转化为独立事件，不拟合流量—时滞规律。当前入口不改变正式V2、第一问或第三问字段。最终论文审查见 [九项审查](../reports/第二问最终论文审查.md)，冻结文件指纹见 `paper_freeze_manifest.json`；仅描述已保存窗口的流量分层，不再拟合新时滞模型。研究代码的历史冻结检查应在对应Git提交上复现；本目录是后续正式统计收尾，不覆盖科研证据。
