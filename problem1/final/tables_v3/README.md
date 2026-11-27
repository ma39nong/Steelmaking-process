# 第一问出版表格V3

正文为T1、T2、T3(a)/(b)、T4、T5五个论证主题；附录A1—A7保留完整数值。每表有可编辑Markdown、CSV、booktabs LaTeX片段；preview为排版检查，不作为唯一表源。

数值从源CSV/JSON生成：P/R/F1三位小数，RMSE均值及SD两位，强度统计整数，时间相邻相关四位、波长相邻相关六位，修复率ppm两位。先比较未舍入值再标最优，平局共同加粗。相对比例按总计分母计算。

表格内所有英文字母、数字、表头及统计标记统一使用真实Times New Roman。Markdown与CSV保存可编辑内容，显示字体由阅读器决定；正式输出以预览PDF和LaTeX排版字体为准。

LaTeX片段需使用XeLaTeX或LuaLaTeX，在全文模板中加载booktabs及fontspec，并安装Times New Roman；不是独立LaTeX文档。所有文件的中文标题、协议、单位和舍入说明在同名Markdown中。预览为7.2英寸单栏宽，8 pt字；数据表可作为双栏跨栏表。

复现：`.venv/bin/python problem1/final/tables_v3/build_tables.py`。来源见[source_mapping.csv](source_mapping.csv)，原十表逐项审核见[table_review_v3.md](../reports/table_review_v3.md)。
