from common import *
R=P/'references';items=json.loads((R/'retrieval.json').read_text())
notes={
'V13':'''原文核心阅读范围：§1问题、可辨识性与定理，§4数值实验、§5求解器；附录证明未逐行复核。PDF p4式(1.1)为 min ||L||*+λ||S||1 s.t.L+S=M；p5式(1.2)–(1.3)约束奇异向量非相干性，p8式(1.5)推广到缺失观测。p5明确整列稀疏项可能同时低秩，从而不可分。§4随机矩阵与视频/人脸验证，原论文λ=1/sqrt(max(n,m))的结论有严格抽样与非相干条件。

本题不是等式约束核范数PCP复现：保留V1固定秩+带噪平方项+稀疏罚，额外引入可见数据自适应时间项和条件饱和项。矩阵本身的低秩与随机缺失可借鉴，整帧、整列异常及结构缺失不能套用其精确恢复保证。真实光谱峰也可能稀疏，S不自动等于仪器异常。''',
'V14':'''出版物JMLR 11(80):2287–2322，2010；作者Rahul Mazumder, Trevor Hastie, Robert Tibshirani。JMLR页面与原PDF核实。10.5555/1756006.1859931是文献数据库标识，Crossref404，索引DOI置空，不能伪称出版社DOI。

核心阅读§1–7与实验§8–9。印刷p2288式(3)：.5||PΩ(X−Z)||F²+λ||Z||*。p2291式(8)–(9)的近端为U diag((σ−λ)+)Vᵀ，p2292 Algorithm1每轮补上缺失再软阈值SVD。§4给凸目标的下降及收敛分析；p2297–2298 §6及Algorithm2 HARD-IMPUTE以硬谱阈值处理非凸秩罚。Netflix数据及模拟实验验证规模与计算效率，不能直接转移其时间或误差数字。

本项目B1更接近固定秩HARD-IMPUTE，不是核范数SOFT-IMPUTE。V2小矩阵使用精确截断SVD与目标回溯；全炉随机SVD仅保证每个接受步的已计算目标不增加，不继承凸全局最优保证。V2另作核范数近端的小型数学核对，不把其视为本题最优模型。''',
'V15':'''已核实SIAM Journal on Numerical Analysis17(2):238–246，1980，F.N.Fritsch与R.E.Carlson，DOI10.1137/0717021。出版社403，未取得原论文全文；不伪造原公式编号、实验结果或声称复现该算法。出版社摘要确认其研究分段三次函数单调性的充要条件。SciPy PCHIP实现主要引用1984论文，不能等同于1980原算法。可作历史背景；本轮不能提供该原文页码级复核，这是明确文献缺口。''',
'V16':'''已核实SIAM Journal on Scientific and Statistical Computing5(2):300–304，1984，F.N.Fritsch与J.Butland。出版社403，未取得作者/出版社授权全文，不把第三方镜像默认视为授权；原文公式和原实验未声称复核。

补充依据：SciPy官方PchipInterpolator Notes及本环境源码。相邻割线d_{k−1},d_k异号或有零则导数m_k=0；同号时(w1+w2)/m_k=w1/d_{k−1}+w2/d_k，w1=2h_k+h_{k−1},w2=h_k+2h_{k−1}；端点用单侧三点规则。这是官方实现依据，不冒用原论文公式编号。每段Hermite多项式由端值/端导数唯一确定。一维C1与区间保形不意味着二维全局C2；先时间再波长为本项目工程组合，顺序可能影响二维值。矩形内角点界、非负和过冲均需独立实测。''',
'V17':'''取得arXiv1209.3672v3，2014-07-01，共31页；最终期刊Information and Inference3(3):189–223，2014，DOI10.1093/imaiai/iau006（用户给定iau004对应另一论文，已保留其元数据核查记录）。核心阅读预印本§1–4；理论证明未逐行复核。p3式(1)为sign(M+Z)；p6–7 §2及§3.1式(3)最大化二元对数似然，约束||X||*≤αsqrt(r d1 d2)及||X||∞≤α。p11–13 §4.1含谱投影梯度、核范数球投影与ADMM求交。§4.2–4.3为模拟与评分量化实验。

不能把其正负二元噪声模型改名成65535上截顶。这里只借鉴受限观测的可辨识性思想；本题在条件性截顶假设下有L≥1的单侧信息，并保有其他点幅值。若噪声模型已知，可用−logPr(Z≥U−F)建立删失似然；当前没有标定噪声分布，V2使用可解释单侧二次罚和约束分裂，不能宣称统计MLE或给出真实峰高上界。'''}
rows=[];bib=[]
for item in items:
 id=item['id'];meta=R/f'{id}_metadata.json';d=json.loads(meta.read_text()) if meta.exists() else {};title=(d.get('title') or ['Spectral Regularization Algorithms for Learning Large Incomplete Matrices'])[0];authors=' and '.join(a.get('given','')+' '+a['family'] for a in d.get('author',[])) or 'Rahul Mazumder and Trevor Hastie and Robert Tibshirani';year=2010 if id=='V14' else d.get('published',{}).get('date-parts',[[2014]])[0][0]
 body=f'# {id} {title}\n\n作者：{authors}。年份：{year}。DOI：{item["doi"] if id!="V14" else "无已核实出版社DOI"}。\n获取状态：{item["status"]}。来源：{item["url"]}。\n\n'+notes[id]+'\n\n适用范围：本项目时间—真实非均匀波长矩阵；无需氧气/烟气输入。原文结果与本项目改造、实验分别陈述。\n'
 (R/'notes'/f'{id}.md').write_text(body);rows.append(dict(id=id,title=title,authors=authors,year=year,doi=item['doi'] if id!='V14' else '',status=item['status'],url=item['url'],pages=item.get('pages',''),sha256=item.get('sha256',''),notes=f'notes/{id}.md'))
 bib.append('@article{'+id+',\n title={'+title+'},\n author={'+authors+'},\n year={'+str(year)+'},\n url={'+item['url']+'}'+(',\n doi={'+item['doi']+'}' if id!='V14' else '')+'\n}\n')
import pandas as pd
pd.DataFrame(rows).to_csv(R/'literature_index_v2.csv',index=False);(R/'references_v2.bib').write_text('\n'.join(bib));(R/'README.md').write_text('''# V2补充文献

沿用根目录 references/ 的R01–R12，绝不修改V1笔记。新增V13–V17见本目录索引和中文笔记。公开全文缓存及逐页提取在本地fulltext/，不默认有再分发权，不随Git推送。三篇公开全文已获取，两篇SIAM全文未获取；详细阅读范围如实列明，不声称逐行验证所有证明。

PCHIP公式补充来源：[SciPy官方文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html)。SOFT-IMPUTE来源：[JMLR正式全文](https://jmlr.org/papers/v11/mazumder10a.html)。原文公式和本题推导分开，见 reports/01_数学推导.md。
''')
