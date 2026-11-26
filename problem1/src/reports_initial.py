from audit import ROOT,P
import pandas as pd
q=pd.read_csv(P/'results/audit/data_quality_summary.csv');w=pd.read_csv(P/'results/audit/wavelength_summary.csv');e=pd.read_csv(P/'results/eda/statistics.csv')
def table(df):
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(str(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
(P/'reports/01_问题重述与分析.md').write_text('''# 问题重述与分析

第一问以十炉次独立的离散时间—波长光谱矩阵为输入，识别不可靠观测并建立可在实际观测域内任意坐标计算的连续光强函数。输出异常类别、有效掩码、连续模型参数及十炉次曲面。

异常处理旨在减小仪器误差影响而保留真实谱峰与工艺状态变化；连续重构不等于将图片像素增密。每个炉次具有自己的实际Sequence时间轴、波长校准网格和持续时间，分别拟合。第一问仅用光谱，不将氧耗、CO2、CO作为拟合输入。

该问题有两个不同尺度：光谱相邻通道往往高度相关，而时间变化可能迅速且非平稳。必须保留孤立真实窄谱线和整体工艺突变，不能按百分位直接裁剪。65535上量程候选与一般离群值分开处理。

统计异常没有独立人工真值，故用保守掩码、注入实验和遮蔽重构验证分别评价。真实超量程光强可能不可识别，不能把平滑曲面当作已测真值。
''')
(P/'reports/02_模型假设与符号说明.md').write_text('''# 模型假设与符号说明

1. 十个文件Sequence均连续，按t=(Sequence−Sequence1)×0.5s构造时间。0.5s来自题目，表格本身未记录硬件时间戳；若将来有缺号，禁止以行号补时间。
2. λ来自逐文件表头，不要求均匀；十炉实际有两类校准网格，分别保存。
3. 信号局部相关但可非平稳并有真实窄峰；相关性由EDA核验，不能推导时间处处光滑。
4. 仪器异常与真实极端光谱可能统计重叠，局部检测仅为候选证据。
5. 65535结合16bit ADC和局部平台视为高度疑似截顶，单独遮蔽/删失敏感性处理。缺少采集软件配置，不能确定全部限幅来源。
6. 不假设所有噪声Gaussian；采用MAD/Huber，标准化尺度仅由可见训练数据估计。
7. 所有连续函数只在各自观测矩形内计算，不外推。独立遮蔽误差针对参考观测，不是干净真值误差。

|符号|含义|单位|
|---|---|---|
|h|炉次1至10|无|
|t,λ|实际时间、波长|s、nm|
|Y,F,E,S|观测、连续信号、测量误差、异常分量|原始记录单位|
|M|有效观测掩码|0/1|
|B_t,B_λ|实际坐标处B样条基矩阵|无|
|C|归一化强度系数矩阵|无|
|K_t,K_λ|基函数个数|无|
|α_t,α_λ,α_tλ|分方向及混合粗糙度系数|对应归一化目标|
|δ,ν|Huber转折点、训练噪声尺度|标准化后无量纲|
|U|疑似上量程|65535记录单位|

题目写W/cm²，整数上限又呈ADC特征，缺乏辐照度校准文件，图中谨慎标为recorded units；不擅自作绝对物理标定。
''')
(P/'reports/03_数据审计报告.md').write_text('# 十炉次数据审计报告\n\n新版附件1文件名有重复.xlsb后缀，但SHA256精确匹配指定值，因此按内容而非文件名判断版本。原文件复制到只读raw，本轮没有使用旧版附件1。各文件读2051列，其中2047为光谱，最后3列仅审计保存；未加入拟合。\n\n'+table(q[['heat_id','records','minimum','maximum','sat65535','sat_fraction','saturation_components']])+'\n\n全部Sequence连续；NaN/Inf/零/负值/重复光谱均为0；未发现整行或整列恒定。相邻观测个别相同值是量化/平稳现象，不能直接判异常。逐时间帧及逐通道统计见results/audit。\n\n'+table(w[['heat_id','start_nm','end_nm','delta_min_nm','delta_max_nm','uniform']])+'\n\n前4炉340.535–1027.963nm，后6炉339.390–1023.303nm；题目正文1027.678与自身示例表末列1027.963不一致，使用实际表头。6C02057只有775条，不能按约1400行强制补齐。8炉有65535，6B03205和6B03270没有。重复上界与16bit支持截顶判断，但长饱和区真实峰值未知。\n')
(P/'reports/04_探索性数据分析.md').write_text('# 十炉次探索性分析\n\n'+table(e[['heat_id','adjacent_time_correlation','adjacent_wavelength_correlation','median_abs_time_derivative','median_abs_wavelength_derivative']])+'\n\n相邻波长相关极强，时间相关较弱，特别6C02058约0.59；因此不能把整条快速变化的帧直接删除，时间强平滑必须经验证。两方向变化率单位不同，不能直接比较数值大小。\n\n每炉raw下保存真实轴热力图、早中晚光谱、400/589/770/950nm最近实际通道时间曲线、光强分布与最大值波长曲线；有截顶炉次另存饱和局部图。曲线显示多种宽带与局部峰结构，单一固定裁剪会损失信息。\n\n'+''.join(f'## {h}\n\n![原始谱及统计](../figures/raw/{h}_eda.png)\n\n![原始热力图](../figures/raw/{h}_heatmap.png)\n\n' for h in q.heat_id))

if (P/'results/audit/extended_quality.csv').exists():
 with (P/'reports/03_数据审计报告.md').open('a') as f:
  f.write('\n## 连续恒值、重复通道和整体跳变补充审计\n\n'+table(pd.read_csv(P/'results/audit/extended_quality.csv'))+'\n\n重复通道保留原始表头与观测，避免擅自减少2047通道；恒值区间和整体跳变是诊断候选，不直接删除。逐通道组号及最长恒值段另存。\n')

eda_report=P/'reports/04_探索性数据分析.md'
eda_report.write_text(eda_report.read_text().rstrip()+'\n')
