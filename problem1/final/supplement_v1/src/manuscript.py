"""Assemble independent seven-section manuscript, audits and source-to-publication map."""
from pathlib import Path
import re,csv,json
import pandas as pd
from PIL import Image
S=Path(__file__).resolve().parents[1];F=S.parent;R=F/'reports';O=F/'figures/publication_v3';T=F/'tables_v3'
old=(R/'问题一_中文论文_V2_英文配图版.md').read_text()
images=re.findall(r'!\[([^]]*)\]\(([^)]+)\)',old)
names=['F01_data_structure','F02_model_mechanism','F03_local_evidence','F04_continuous_reconstruction','F05_joint_validation','F06_actual_repair_results','F07_ten_heat_surfaces']
caps=[
'光谱分组结构及低秩依据。上图将逐炉时间中位谱与10–90%分位曲线映射至共同波长网格，再按6B四炉、6C六炉等权平均；阴影为炉内分位曲线的平均。下图分别给出原始抽样矩阵的未解释能量，以及两组目标通道暗、亮各15%帧的平均谱形。参考竖线指示Na I、O₂ A带和K I候选区域。宽带共同结构与局部尺度差异共同支持低秩背景及自适应检测。',
'自适应多证据异常修复机制。以6B03205的14.5 s、506.707 nm节点为例，在相同局部观测上展示时间参考、波长参考、标准化证据、综合评分及实际修复谱形；下部概括评分、秩4 Huber拟合和支撑保真的作用关系。图中数值来自冻结矩阵和最终评分掩码，尖峰节点由1767修复至523.89 counts。',
'局部修复与持续谱峰保留的证据链。第一行定位6B03205的代表尖峰并给出前后时间参考；第二行展示四个高证据尖峰的实际变化与综合评分；第三行展示6C02108在304.0 s、766.413 nm的持续特征及邻域时间变化。前者说明局部结构不一致的修复，后者说明持续窄谱形的保留。',
'连续重构选型及双线性单元。上排为独立真实坐标降采样检验的炉次宏平均RMSE，三次样条与该边界下的零平滑张量B样条合并表示；下排使用真实四角节点展示非均匀坐标单元、权重和查询值。PCHIP在1.0 s下误差较小，双线性在其余场景更小，并具有非负和凸包性质。',
'异常识别与修复的共同协议验证。左上以精确率—召回率及等F1曲线展示有限策略；右上展示共同300步配置的污染区收益与未污染区变化；左下按波长区给出逐案例F1均值；右下单独展示发布600/1200步流程在同一批案例上的恢复变化。两种迭代配置分别解释，区域F1均值也不等同于汇总计数F1。',
'十炉实际修复分布及局部变化。左侧给出1043个修复节点和其中465个回退节点的逐炉计数；右侧给出同一真实局部的修复前后曲线和有符号差值热图。回退是修复集合的子集，差值只作用于实际修复节点。',
'十炉最终连续光谱总览。各炉使用实际时间和波长定义域，统一视角、0–65535 counts纵轴及共享色标；展示网格除规则采样外保留全局最大值和已核实案例坐标。各炉完整连续函数以原生节点参数定义，图像采样只用于展示。'
]
# Replace factual coordinate/value in mechanism caption directly from frozen source.
pt=pd.read_csv(S/'results/repair_reliability_points.csv').query("heat_id=='6B03205' and row==29 and column==466").iloc[0]
caps[1]=caps[1].replace('506.707',f'{pt.wavelength_nm:.3f}' if 'wavelength_nm' in pt.index else '506.707').replace('523.89',f'{pt.repaired:.2f}')
# Source image audit and figure mapping: all 19 originals were visually inspected.
reasons=[
'原能量曲线集中在0.8–1且图例拥挤；改为未解释能量对数坐标，并与分组光谱、暗亮目标窗口组合。',
'原流程为大留白文字框；用同一真实异常的双向参考、评分与输出构成机制证据链。',
'原连续对照四组柱中三次与张量值重复；合并等价表示，用折线显示采样跨度并增加真实单元。',
'原RMSE共轴柱使未污染误差几乎不可见；改为检测权衡、恢复代价散点及分区域对照。',
'原8小图重复坐标且信息密集；收敛移附录，按各炉初始记录归一化目标并保留相对变化。',
'原计数与比例双柱重复；改点线分布并明确回退子集，配同一案例曲线与差值。',
'原四尖峰重复图例且缺时间证据；改空间定位—时间参考—四例—持续谱线保留对照。',
'原6C02058宽峰凹口不能仅作成功示例；移附录，联查潜变量、双向参考与回退输出。',
'原5×2总览留白大、各面板z轴自动缩放；统一色标/视角并强制保留全局峰值坐标。'
]
newfor=[names[0],names[1],names[3],names[4],'B01_solver_stability',names[5],names[2],'B02_wide_peak_tradeoff',names[6]]
rows=[]
for k,(alt,rel) in enumerate(images,1):
 q=(R/rel).resolve();im=Image.open(q);dest=newfor[k-1] if k<=9 else f'A{k-9:02d}_{["6B03205","6B03206","6B03270","6B03271","6C02057","6C02058","6C02059","6C02060","6C02107","6C02108"][k-10]}_surface'
 sec='main' if dest.startswith('F') else 'appendix'
 rows.append(dict(old_number=f'5-{k}',old_path=q.relative_to(F).as_posix(),width_px=im.width,height_px=im.height,new_asset=f'{sec}/{dest}',placement='正文' if sec=='main' else '附录',reason=reasons[k-1] if k<=9 else '原独立曲面规则60×90采样可能漏过窄峰；保留实际定义域，以全局最大值和案例坐标补充展示网格。',review='已逐张打开原图检查'))
for k,r in enumerate(rows,1):
 r['chapter']={1:'5.1',2:'5.2-5.4',3:'5.5',4:'5.6',5:'5.6/appendix',6:'5.7',7:'5.2/5.7',8:'5.6/appendix',9:'5.7'}.get(k,'5.7/appendix')
 r['data_source']=('problem1/final/data/*_interface.npz; problem1/final/supplement_v1/results; problem1/v3/results/singular_energy.csv' if k==1 else 'problem1/v2/results/interpolation_independent_macro.csv' if k==3 else 'problem1/final/results/unified_comparison.csv; exact_final_pipeline_confirmation.csv; confirmation_band_summary.csv' if k==4 else 'problem1/final/results/*_solver.json' if k==5 else 'problem1/final/data/*_interface.npz; supplement_v1/results/repair_reliability_points.csv')
 r['script']='problem1/final/supplement_v1/src/publication_figures.py'
 r['parameters']='render_config.json; protocol.json; frozen final config'
pd.DataFrame(rows).to_csv(O/'figure_mapping.csv',index=False)
(O/'figure_audit_19.md').write_text('# 原19图逐张审查与映射\n\n审查对象为英文配图版实际引用的19张PNG，逐张打开检查布局、坐标、峰值采样及论证作用。参考册仅用于布局与表达，所有新图数值来自项目数组/CSV。原文件保留。\n\n'+'|原图|新位置|审查及处理|\n|---|---|---|\n'+''.join(f"|{r['old_number']}|{r['new_asset']}|{r['reason']}|\n" for r in rows)+'\n\n主图统一6.8英寸设计宽、600 dpi PNG及SVG/PDF，英文图内标注；正文中文图注解释统计对象、计算方式和发现。两组光谱先映射实际共同网格再等炉汇总。全部曲面统一记录尺度，展示采样保留全局最大值。SVG保留可编辑文字。\n\n优先完成F01/F02/F03/F05并实际打开审查后，生成其余图。F02证据面板已扩大下侧空白，使P、c、S注释避开证据点线。发布曲面采样与原始矩阵分开审计，见surface_sampling_audit.csv。\n')
(O/'figure_captions.md').write_text('# 正文与附录图注\n\n'+ '\n\n'.join(f'图5-{k+1}　{c}\n\n文件：main/{names[k]}.{{png,svg,pdf}}' for k,c in enumerate(caps))+'\n\n附图A1—A10　十炉独立连续光谱。统一记录尺度及视角，实际定义域逐炉取值，完整参数和质量标记沿用正式接口。\n\n附图B1　四炉归一化目标函数及相对迭代变化，竖线表示600与1200步。\n\n附图B2　宽峰修复的权衡案例：同一节点的记录、潜变量、正式回退、双向参考及差值热图。\n\n附图B3(a–d)　两炉原始/中心化抽样SVD的谱基、时间系数、秩1贡献、秩4背景和残差；谱基单位范数，最大绝对载荷取正，系数同步翻转。\n\n附图B4　两组目标持续通道的逐炉残差轨迹和与排除目标窗口亮度的关系。\n')
# Extract old prose, retain all original equations and algorithm pseudocode.
s=old
s=re.sub(r'^表5-\d+[^\n]*\n\n(?:\|[^\n]*\n)+',lambda m:f'@@TABLE{re.search(r"表5-(\d+)",m[0]).group(1)}@@\n',s,flags=re.M)
fig_counter=iter(range(1,20))
s=re.sub(r'^!\[[^]]*\]\([^)]+\)\n',lambda m:f'@@FIG{next(fig_counter)}@@\n',s,flags=re.M)
s=re.sub(r'^(?:Figure 5-\d+\b|图5-\d+\s)[^\n]*\n','',s,flags=re.M)
tmap={1:'表5-1',2:'表5-1',3:'附表A1',4:'表5-2',5:'附表A3',6:'表5-3',7:'表5-4',8:'附表A6',9:'表5-5',10:'附表A5'}
fmap={1:'图5-1',2:'图5-2',3:'图5-4',4:'图5-5',5:'附图B1',6:'图5-6',7:'图5-3',8:'附图B2',9:'图5-7',**{k:f'附图A{k-9}' for k in range(10,20)}}
s=re.sub(r'表5-(\d+)',lambda m:tmap[int(m[1])],s);s=re.sub(r'图5-(\d+)',lambda m:fmap[int(m[1])],s)
for k in range(1,11):
 dest={1:'T1',4:'T2',6:'T3a',7:'T4',9:'T5'}.get(k)
 text=(T/f'{dest}.md').read_text() if dest else ''
 if k==6:text+='\n\n'+(T/'T3b.md').read_text()
 s=s.replace(f'@@TABLE{k}@@',text)
def img(n):return f'![图5-{n}](../figures/publication_v3/main/{names[n-1]}.png)\n\n图5-{n}　{caps[n-1]}\n'
for k in range(1,20):s=s.replace(f'@@FIG{k}@@',img({1:1,2:2,3:4,4:5,6:6,9:7}[k]) if k in [1,2,3,4,6,9] else '')
s=s.replace('## 5.3 鲁棒低秩异常修复模型建立',img(3)+'\n## 5.3 鲁棒低秩异常修复模型建立')
# Additional scientific evidence in the corresponding sections, not a version chronology.
ins1='''为比较两组不同波长网格上的记录结构，先逐炉计算时间中位谱及10%—90%分位谱，再线性映射至共同范围340.535—1023.303 nm内的0.5 nm查询网格，最后分别对四炉和六炉等权平均。图5-1的阴影表示炉内分位曲线的平均，因而既呈现各组典型幅度，也保留随时间变化的强度范围。两组的宽带形状与窄通道结构存在差异，说明跨炉分析应统一物理坐标而非列号，局部异常阈值也应相对于所在波段的波动尺度定义。

为进一步区分平均背景与共同动态，在6B03205、6C02058上对同一每8帧抽样矩阵分别作原始SVD与逐通道去时间均值的中心化SVD。两炉原始前四阶能量分别为99.952402%、99.959361%，中心化后分别为99.900322%、99.957039%。中心化减少固定平均谱的影响，剩余变化仍集中于少量模态。这里的中心化结果对应同一原始抽样矩阵，与历史掩码填充全网格的中心化结果分开解释。每个谱基取单位范数，最大绝对载荷取正，时间系数同步翻转；谱基、系数、秩1及秩4背景见附图B3。

'''
s=s.replace('## 5.2 数据预处理与异常识别',ins1+'## 5.2 数据预处理与异常识别')
ins2='''持续通道诊断使用同一实际波长参考，以残差符号持续率不低于95%、残差中位数绝对值不低于250 counts、标准化绝对残差中位数不低于4筛查全部20,470个炉次—通道组合，得到9个统计候选，其中6个为383.293 nm正残差通道，另3个为邻峰负残差响应；邻域预测共享相邻值，因此这些候选并非9个独立故障。6C组383.293 nm在六炉均通过筛查，其正残差持续率为99.48%—100%；6B组835.887 nm虽有98.02%—99.53%的正残差持续率，标准化中位偏离仅2.08—3.11，未通过该规则。由此可见，持续性须与幅度、局部尺度共同判断。

对两个目标窗口进一步取宽带亮度最低和最高各15%的帧，并以10%、20%作敏感性分析；亮度由排除目标±2 nm窗口的实际波长梯形积分除以有效带宽得到。局部稳健回归 $d(t)=\\alpha+\\beta B(t)+\\varepsilon(t)$ 将固定残差与亮度相关分量分开，区间由16帧移动块、64次重采样估计。6B四炉斜率均为正，6C斜率则正负并存，说明同一突出通道可能包含不同的背景或响应成分。其时间轨迹、峰宽描述及回归区间见附图B4及补充通道表。

物理结构的检查以NIST空气参考为索引：Na I双线588.9950/589.5924 nm、K I双线766.48991/769.89645 nm，并检查约760 nm的O₂ A带。[补充物理来源](../supplement_v1/reports/physical_sources.md)记录参考约定和原文状态。6C02108在304.0 s、766.413 nm的持续特征原始与正式值均为39674 counts，与K I第一参考位置相近；其邻帧结构支持保留。波长接近用于定位候选，实际仪器配置及空气/真空约定统一列入模型评价。

'''
s=s.replace('### 5.2.2 双向局部预测及自适应尺度估计',ins2+'### 5.2.2 双向局部预测及自适应尺度估计')
s=s.replace('### 5.3.2 固定秩Huber鲁棒矩阵补全','''### 5.3.2 固定秩Huber鲁棒矩阵补全

低秩近似可写为 $L=\\sum_{k=1}^{4}a_kv_k^{\\mathsf T}$。每个 $v_k$ 表示跨时间共享的谱形，$a_k$ 描述其时间系数；共同表示使一个缺失节点能够受到其他时刻及其他波长可见数据的共同约束。SVD用于刻画这种相关结构，实际补全则在可见集合上最小化鲁棒损失，两者的用途分别是结构分析与缺失估计。对持续谱线和细节的保留由后续观测保真完成。
''')
# New validation separated from historical masking and actual release.
ins6='''#### 持续通道扰动的有限补充验证

在历史局部片段间隙另固定32×96片段，开发与确认采用不同时间槽，两个波长起点分别为800和1600。在中心候选范围选择局部标准化残差较小的普通通道，注入两种幅度的固定偏置、增益或冻结扰动。通道修复使用已知注入列掩码，结构检测另行计分。确认集炉次宏平均中，冻结、增益、偏置的修复前RMSE分别为2296.31、6622.75、7520.75 counts；双向线性修复均为126.85 counts，既有Huber低秩为117.66 counts，未污染区变化均为0。结果说明在这一普通局部背景及正确定位条件下，两种已有方法均可恢复通道污染，其中低秩略优。该已知掩码实验与自动异常定位检验承担不同作用。自动通道筛查在确认偏置/增益扰动上的Precision为1/3、Recall为1、F1为0.5，冻结扰动Recall为0；邻域参考使突出通道两侧也出现负残差，提示自动定位和已知掩码恢复应分开评价。

'''
s=s.replace('### 5.6.3 连续插值与降采样重构验证',ins6+'### 5.6.3 连续插值与降采样重构验证')
ins_lim='''#### 修复可靠性与适用范围

实际掩码核验得到1043点由465个回退点和578个非回退点组成。高亮定义为各炉积分最高15%的帧，快变定义为相邻积分绝对变化最高15%的帧，近截顶定义为11×21邻域存在顶码。回退点位于高亮帧的比例为100%，快变比例58.28%，近截顶比例27.96%；非回退点对应82.18%、52.77%、14.19%。这些分布说明，高亮与局部参考分歧值得作为修复可靠性的重点检查对象。

6C02058在339.5 s、674.478 nm的原始值为65067 counts，潜在低秩值为63601.74 counts，而正式回退值为34559.91 counts。时间参考34288.00与波长参考65339.49明显分歧；综合评分15.746中波长得分仅0.02385，时间和背景得分较大。该位置位于宽峰内，主要凹口出现在稳定性回退阶段，附图B2以完整曲线和差值揭示这一代价。预先固定的近顶码保护规则在新开发片段减少了部分非污染区变化，在确认片段未触发，因此本轮保留既有正式输出，不据单个案例反向调整规则。

历史人工确认中17个新增峰的诊断回放逐例复现原计数，其中15处位于未注入节点、2处位于注入节点，14处使用数值回退。这些计数属于人工污染验证，不属于真实1043点的修复账本。持续通道整列替换仅作为独立敏感性输出，平均相对积分变化逐炉为0.0033%—0.0200%，平均峰高变化最高2.09%；结合亮度关系差异及谱线候选证据，正式处理继续保留这些通道。记录网格半高宽只描述采样谱形，仪器实际分辨率及波长标定配置需结合设备配置解释。模型的主要适用范围是具有局部连续性与共同背景的离线光谱重构，整段时间缺失、超量程真实峰高及弱污染与窄谱线的区分仍是主要限制。

'''
s=s.replace('## 5.7 十炉次最终求解结果',ins_lim+'## 5.7 十炉次最终求解结果')
s=re.sub(r'本模型适用于具有局部连续性和共同光谱背景的离线重构。模型评价集中于以下范围：[^\n]+\n','',s)
s=s.replace('仍是主要限制。','仍是主要限制。条件性饱和计算的104个分块中1块保留有限迭代近似标记，相关质量信息随接口输出。')
# Clarify updated display, full function remains frozen.
s=s.replace('在每炉实际定义域内分别构造含 $2n_h-1$ 个时间点及4093个波长点的等距查询坐标，生成密集网格，绘制总览及各炉独立三维图。节点模型仍采用原始非均匀波长坐标，等距坐标仅用于曲面展示。','已有密集查询网格包含 $2n_h-1$ 个时间点及4093个波长点，并完整保存。图5-7及独立附图进一步在规则展示采样中加入原网格全局最大值和已核实案例坐标，使窄峰与典型修复位置在图中得到表达；十炉展示最大值均与正式原网格最大值一致。节点模型仍采用原始非均匀波长坐标，绘图采样不改变连续函数参数或正式修复结果。')
s += '\n\n完整图表、逐炉数值及补充研究入口见[图表与研究附录](问题一_V3_图表与研究附录.md)；统一假设、符号和参考文献沿用[全文整合材料](问题一_V2_全文整合材料.md)，新增物理来源单独核实并在附录链接。\n'
# Normalize empty space and ensure no assembly markers remain.
s=re.sub(r'\n{4,}','\n\n\n',s);assert '@@' not in s
(R/'问题一_中文论文_V3_研究深化与图表重构版.md').write_text(s)
# Comprehensive appendix, editable tables and scientific sources.
a='# 第一问图表与研究附录\n\n正文图5-1—5-7服务于数据—机制—证据—连续模型—验证—修复—结果的论证；历史正文和数据接口完整保留。本附录保存完整表与独立曲面，数值来源及统计口径随表给出。\n\n## 可编辑完整表\n\n'
for k in range(1,8):a+=(T/f'A{k}.md').read_text()+f'\n\n[CSV](../tables_v3/data/A{k}.csv) · [LaTeX片段](../tables_v3/latex/A{k}.tex)\n\n'
a+='## 补充科学证据\n\n[补充研究报告](../supplement_v1/reports/补充研究与模型依据.md) · [物理来源核实](../supplement_v1/reports/physical_sources.md) · [固定协议](../supplement_v1/protocol.json)\n\n全部20,470条通道诊断、9个通道候选、暗亮回归、峰宽/邻域支持、逐点修复分级、17个新增峰定位及新片段验证CSV保存在[results](../supplement_v1/results/)。派生数组由脚本复现，正式接口仍使用已发布冻结版本。\n\n## 十炉独立连续曲面\n\n'
heats=pd.read_csv(F/'results/final_summary.csv').heat_id
for k,h in enumerate(heats,1):a+=f'![附图A{k}](../figures/publication_v3/appendix/A{k:02d}_{h}_surface.png)\n\n附图A{k}　{h}连续光谱。实际定义域、统一色标和视角，展示采样保留全局最大值及案例位置。\n\n'
for name,cap in [('B01_solver_stability','附图B1　四炉归一化目标函数及相对变化。'),('B02_wide_peak_tradeoff','附图B2　339.5 s宽峰内的修复权衡及回退机制。'),('B03_6B03205_svd_raw','附图B3(a)　6B03205原始抽样SVD。'),('B03_6B03205_svd_centered','附图B3(b)　6B03205同样本中心化SVD。'),('B03_6C02058_svd_raw','附图B3(c)　6C02058原始抽样SVD。'),('B03_6C02058_svd_centered','附图B3(d)　6C02058同样本中心化SVD。'),('B04_persistent_channel_evidence','附图B4　持续通道的时间轨迹与亮度关系。')]:a+=f'![{cap}](../figures/publication_v3/appendix/{name}.png)\n\n{cap}\n\n'
(R/'问题一_V3_图表与研究附录.md').write_text(a.rstrip()+'\n')
# Ten original tables audited separately, with precise destination and rationale.
tr=[('数据规模','T1/表5-1','2047通道、0.5 s及两组波长范围移表注；逐炉保留帧数、时长及质量统计。'),('强度及相关','T1及A1','合并中位强度和相关性，P99及完整统计移A1；相邻对展平、隔4抽样口径保留。'),('SVD四阶能量','F01及A1','避免并排炉号和99.9%挤压，正文展示1−E(r)；原始与中心化严格分开。'),('正式参数','T2及A2','数学意义与依据进正文，时限/halo等实现参数进附录。'),('连续插值','F04及A3','三次/张量等价列合并；1.0 s PCHIP最优按真实未舍入值高亮。'),('共同污染','T3(a)/(b)及A7','纳入S=6和新增峰代价；P0精确率—*；300步与发布600/1200分开。'),('局部缺失','T4及A4','正文三个代表场景；全部11场景附录，平方低秩最优项如实加粗；全矩阵独立协议不拼表。'),('求解检查点','A6及B1','归一化曲线表达趋势，同一运行目标及δ完整保存；不把跨炉绝对目标用于优劣排名。'),('实际修复','T5','修复率改ppm；797下/246上及465回退子集明确，总率按27,816,683加权。'),('四例','F03及A5','图给时空证据，表留精确坐标、变化量、S及回退；持续峰为保留对照。')]
review='# 原十表审查与论文级重构\n\n所有数字由源CSV/JSON导出，历史表不改动。正文五个主题对应六个物理表块（T3拆为a/b），三线LaTeX源和单栏预览均保留，Markdown提供编辑与检索。\n\n|原表|主题|新位置|处理依据|\n|---|---|---|---|\n'+''.join(f'|5-{i}|{x}|{d}|{why}|\n' for i,(x,d,why) in enumerate(tr,1))+'\n## 数值与公平性复核\n\n共同300步：S4污染2390.39、整体198.50、未污染68.94，新增48；S6污染3440.22、新增17。发布S4在同一确认案例上污染2359.80、整体195.80、未污染66.68，新增17，单独A7说明。两个17不合并解释。P/R/F1汇总计数，误差炉次宏平均，区域图逐案例F1均值，口径分别标注。\n\n未污染0不能独立代表修复能力；零预测阳性的P0精确率未定义。各表在舍入前判定最优，平局共同高亮。T5修复1043、向下797、向上246、回退465；总率37.50 ppm而非逐炉百分比均值；14805顶码和1129帧候选不与点修复相加。\n\n源与计算入口：[source_mapping.csv](../tables_v3/source_mapping.csv)、[build_tables.py](../tables_v3/build_tables.py)。列名统一量纲、每指标统一精度。LaTeX为booktabs片段，三横线无竖线；预览按7.2英寸、8pt跨双栏宽设计，正式模板中应作跨栏表，不能直接压缩至半栏。预览已渲染检查，完整表保留附录，不靠栅格截图交付。\n'
(R/'table_review_v3.md').write_text(review)
print('MANUSCRIPT_AUDITS_COMPLETE',len(images))
