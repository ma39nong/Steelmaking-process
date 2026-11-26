# 第一问文献索引（2026-11-26）

入口：`literature_index.csv`、`references.bib`、`notes/R01.md` 至 `notes/R12.md`。

## 获取与阅读范围

已获得并阅读 R03 作者预印本、R04 出版社OA全文、R09厂家手册、R10作者预印本；R08已阅读出版社HTML正文，保存浏览工具提取的正文。R01、R02、R05、R06、R07、R11、R12 尚未取得可用全文，笔记逐项标为待核实，绝不以摘要冒充全文。

补充来源 S01：Eilers/Marx 的 JOPS 官方CRAN手册（`fulltext/R01_cran.pdf`，ps2DNormal pp.56–57，rowtens pp.80–81）；S06：Schafer (2011), DOI 10.1109/MSP.2011.941097，大学教学公开PDF（`fulltext/S06.pdf`，pp.111–117）；S11：SciPy官方 `special.huber` 文档。补充来源支持独立数学推导和实现，并不等于已读相应原始文献。

R03 的PMC和arXiv标题及核心模型相符，但2022预印本与最终刊发版本作者名单不同，公式编号一律以本地15页预印本为准，不混用版本。R12的DOI含2025，但Crossref及出版社列出首次上线2026-07-06，不擅自写作2025发表。

## 文件管理与合法来源

`acquire.py` 仅请求作者机构、出版社、官方存档公开资源；失败及HTTP验证页面不计全文。`metadata/*_crossref.json` 保存书目原始记录；`*_acquisition.json` 和 `acquisition.log` 保存获取失败与来源证据。HTTP 200 仍可能是验证页面，须检查PDF签名及正文。个别作者站证书链无法验证，日志已记录；文件保留SHA256。

`fulltext/` 为本机阅读缓存，默认不提交第三方论文到Git，以免将合法阅读下载误当成有权再分发。索引、链接、笔记与获取程序可提交；OA文件也保留原版权。数据同样保留本机，使用哈希索引复现。

运行：项目根目录 `.venv/bin/python references/acquire.py`（部分外站仍可能拒绝访问）。`build_index.py` 重建书目索引。未取得的关键文献可由用户提供合法作者稿或机构订阅下载稿继续补充，当前研究不伪造原文细节。

## 核心算子数值复核

`reproductions/r03_log_prox.py`复现R03预印本式(22)–(23)的列log近端及奇异值log近端，不是完整L3S3TV。`verify.py`对100组随机半径/参数与20001点密网格目标比较，并检查零惩罚恒等性；同时核对Huber标准定义。结果存`metadata/method_reproduction_checks.json`。原稿证明式(25)的1/4与式(22)的1/2不一致，实现以明示目标重新推导，不照抄疑似排版误差。

补充S13：SciPy官方PchipInterpolator完整文档，https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html 。用于修正普通双三次插值在快速时间变化间的负过冲；是确定性连续求值改造，不归为R03原算法。
