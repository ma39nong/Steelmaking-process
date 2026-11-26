# 最终方法的文献边界

延续根目录R01–R12及V2的V13–V17索引，不覆盖原文、笔记或获取状态。只保留最终模型需要的论据，不为了引用数量增加算法。

|来源|全文/核对状态|最终用途及限制|
|---|---|---|
|V13 Candès, Li, Ma, Wright (2011), Robust principal component analysis? JACM58(3),article11, DOI10.1145/1970392.1970395|V2作者公开PDF与中文笔记|低秩稀疏思想；本项目非凸固定秩Huber不是凸PCP|
|V14 Mazumder, Hastie, Tibshirani (2010), JMLR11:2287–2322|V2 JMLR全文与笔记|缺失矩阵谱正则、迭代SVD；本项目没有直接使用Soft-Impute目标|
|R11 Huber (1964), Ann. Math. Statist.35(1):73–101, DOI10.1214/aoms/1177703732|原始全文未合法获取；不假称已读|Huber历史归属；实际分段损失核对官方SciPy定义|
|R09 USB4000操作手册|根目录公开手册/笔记|16位ADC与顶码解释；本题具体仪器型号仍不能由数据独立证明|
|V15/V16 Fritsch–Carlson1980/Fritsch–Butland1984|V2元数据、官方SciPy实现说明；原文全文未获取|PCHIP历史对照，不能声称阅读全文复现|

- [V2原文及笔记](../../v2/references/README.md)，[原R01–R12](../../../references/README.md)。
- [JMLR公开全文](https://jmlr.org/papers/volume11/mazumder10a/mazumder10a.pdf)。
- [RPCA作者全文](https://www.columbia.edu/~jw2966/papers/CLMW11-JACM.pdf)。
- [SciPy Huber官方定义](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.huber.html)。其页面的历史论文参考条目存在不匹配，本文使用已核实的1964书目信息，不复制该错误条目。
- [SciPy规则网格插值官方说明](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html)。

2026-11-26核对官方实现说明。完整数学目标、优化更新及文献原方法/本题改造的区别见最终数学推导。R07/R12只保留后续研究背景，本轮不进行预测网络实验。
