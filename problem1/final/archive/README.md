# 第一问写作与图集归档

这里保存早期写作和制图成果，供研究追溯使用。**当前论文、图表和数据入口见[正式交付](../README.md)。**

|目录|保留内容|当前替代入口|
|---|---|---|
|`reports/`|初稿、章节V2/V3、英文配图报告、早期图表清单及校对记录，共9份|[当前正式正文](../reports/问题一_中文论文_V4_机制与结果图优化版.md)|
|`figures/main/`、`figures/appendix/`|最初中文图及逐炉图|[当前图集](../figures/publication_v4/README.md)|
|`figures/english/`|早期英文19图与完整逐炉扩展图|[当前图集](../figures/publication_v4/README.md)|
|`figures/publication_v3/`|上一套组合图及原始核查材料|[当前Times New Roman图集](../figures/publication_v4/README.md)|

此次归档移动345个已跟踪文件，原总大小95,924,165字节；没有删除历史实验、修复数据或选择锁。图像与数值材料按字节保留，Markdown仅修正搬迁后的相对链接。历史路径可通过[迁移清单](relocation_manifest.csv)逐文件查询；清单记录原路径、新路径、原SHA256、归档SHA256及来源提交。目录映射见[relocation_map.json](relocation_map.json)。

旧脚本继续留在原有源代码位置，输出路径已指向本归档区。正式流程仍复用`v2/`、`v3/`及`supplement_v1/`的部分实现，不能按版本号批量删除它们。归档区中的旧结论、旧字体和旧实验协议属于历史记录，不覆盖当前正式结果。

原始字节可从迁移清单的`source_commit`与`original_path`恢复，例如：

```bash
git show <source_commit>:<original_path>
```

不覆盖已有归档内容重跑历史研究。重复实验请使用独立输出目录。
