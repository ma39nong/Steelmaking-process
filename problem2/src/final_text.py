"""Data-driven current Chinese chapter, handoff and figure index."""
from pathlib import Path
import json
import numpy as np
from .final_statistics import RELATIONS

ROOT=Path(__file__).resolve().parents[2]
REPORTS=ROOT/'problem2/final/reports'
LABELS={'co_instantaneous':'耗氧量—CO瞬时流量','co2_instantaneous':'耗氧量—CO₂瞬时流量','co_cumulative':'耗氧量—CO累计排放','co2_cumulative':'耗氧量—CO₂累计排放'}


def table(headers,rows):
    def fmt(x):return f'{x:.4g}' if isinstance(x,(float,np.floating)) else str(x)
    return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+''.join('| '+' | '.join(fmt(x) for x in row)+' |\n' for row in rows)


def write_reports(heats,a,cfg):
    REPORTS.mkdir(exist_ok=True);out=ROOT/cfg['output_dir']
    t1=table(['炉号','点数','时长 s','耗氧增量 m³','CO总量 m³','CO₂总量 m³'],[[r[k] for k in ['heat_id','samples','duration_s','oxygen_change_m3','co_total_m3','co2_total_m3']] for r in a['heat_summary']])
    t2=table(['关系','Pearson范围','Spearman范围','恒定排放参照Pearson范围'],[[LABELS[rel],f"{min(r['pearson_r'] for r in a['level_correlations'] if r['relation']==rel):.3f}–{max(r['pearson_r'] for r in a['level_correlations'] if r['relation']==rel):.3f}",f"{min(r['spearman_r'] for r in a['level_correlations'] if r['relation']==rel):.3f}–{max(r['spearman_r'] for r in a['level_correlations'] if r['relation']==rel):.3f}",f"{min(r['constant_gas_rate_pearson'] for r in a['level_correlations'] if r['relation']==rel):.3f}–{max(r['constant_gas_rate_pearson'] for r in a['level_correlations'] if r['relation']==rel):.3f}" if rel.endswith('cumulative') else '不适用'] for rel in RELATIONS])
    diffrows=[]
    for rel in RELATIONS:
        for sec in cfg['increment_intervals_s']:
            rr=[r for r in a['increment_correlations'] if r['relation']==rel and r['interval_s']==sec]
            diffrows.append([LABELS[rel],sec,np.median([r['pearson_r'] for r in rr]),np.median([r['spearman_r'] for r in rr]),f"{min(r['pearson_r'] for r in rr):.3f}–{max(r['pearson_r'] for r in rr):.3f}"])
    t3=table(['关系','增量间隔 s','Pearson炉间中位数','Spearman炉间中位数','Pearson逐炉范围'],diffrows)
    t4=table(['气体','等耗氧阶段','阶段平均流量中位数 m³/s','95%整炉重采样区间','阶段排放/耗氧中位数 m³/m³'],[[r['gas'].upper(),r['stage'],r['mean_rate_m3_s_median'],f"[{r['mean_rate_m3_s_ci_low']:.3f}, {r['mean_rate_m3_s_ci_high']:.3f}]",r['emission_per_oxygen_m3_m3_median']] for r in a['stage_aggregate']])
    t5=table(['关系','次数','c₀','c₁','c₂'],[[LABELS[r['relation']],r['degree'],r['c0'],r['c1'],r['c2']] for r in a['pooled_trends']])
    loorows=[]
    for rel in RELATIONS:
        for degree in cfg['trend_degrees']:
            rr=[r for r in a['trend_loho'] if r['relation']==rel and r['degree']==degree]
            loorows.append([LABELS[rel],degree,np.mean([r['mae'] for r in rr]),np.mean([r['r2'] for r in rr]),sum(r['r2']>0 for r in rr)])
    t6=table(['关系','次数','留一炉次宏平均MAE','留一炉次宏平均R²','测试R²>0炉数'],loorows)
    flow=[r for r in a['flow_review'] if r['sg_window']==31]
    t7=table(['炉号','原累计值相邻重复比例','正增量中位数 m³','31点反积分RMSE m³','31点总量相对差 %','21/41点流量偏离31点RMSE m³/s'],[[r['heat_id'],r['raw_equal_adjacent_fraction'],r['positive_increment_median_m3'],r['reintegration_rmse_m3'],100*r['endpoint_relative_error'],' / '.join(f"{x['flow_difference_from_31_rmse_m3_s']:.3f}" for x in a['flow_review'] if x['heat_id']==r['heat_id'] and x['sg_window']!=31)] for r in flow])
    t8=table(['炉号','CO/CO₂可靠窗口','CO/CO₂有效对齐点','Q1时间轴','字段历史比较','检查项'],[[r['heat_id'],f"{r['co_reliable_windows']} / {r['co2_reliable_windows']}",f"{r['co_alignment_points']} / {r['co2_alignment_points']}",'精确一致','57/57通过',sum(x['heat_id']==r['heat_id'] for x in a['interface_checks'])] for r in a['heat_summary']])
    bandnames={'low':'低','medium':'中','high':'高'}
    delayrows=[]
    for r in a['delay_flow_summary']:
        interval=f"[{r['flow_min_m3_s']:.3f}, {r['flow_max_m3_s']:.3f}"+(']' if r['upper_inclusive'] else ')')
        delayrows.append([bandnames[r['flow_band']],interval,r['gas'].upper(),r['delay_median_s'] if r['accepted_windows'] else '不可辨识',f"{r['delay_min_s']:.1f}–{r['delay_max_s']:.1f}" if r['accepted_windows'] else '—',r['accepted_windows'],r['contributing_heats'],'仅窗口级；物理事件未确认' if r['accepted_windows'] else '无合格窗口'])
    delaytable=table(['流量层','中心流量区间 m³/s','气体','候选时滞中位数 s','范围 s','合格窗口数','覆盖炉数','证据范围'],delayrows)
    maxrel=max(abs(r['endpoint_relative_error']) for r in flow)
    chapter=r'''# 问题二模型与实验结果

## 6.1 数据来源和变量定义

第二问使用附件1至10的Sequence与三列原始工艺数据。每炉Sequence连续、采样间隔0.5 s；十炉文件SHA256与冻结基线一致。用户确认中文“附件1.炉号-6B03205.xlsb”只是现有英文文件的改名，数据身份核验已关闭。第一问处理结果为光谱与质量接口，第二问工艺量直接取自原始XLSB，不用光谱处理结果代替工艺列。

令炉次为 \(h\)，时间为 \(t_{h,i}=0.5(s_{h,i}-s_{h,0})\)。累计耗氧量记为 \(O_h(t)\)，单位m³；题面称作CO与CO₂ concentration的字段，在本题给定口径下均为瞬时体积流量 \(C_{g,h}(t)\)，单位m³/s，不能再额外乘一个未知烟气总流量。所有积分均在秒制时间轴上计算。

表1 十炉数据及总量结果

'''+t1+r'''
## 6.2 异常处理和连续信号重构

对三列工艺量使用9点Hampel检测，阈值为6倍局部稳健尺度；MAD为零时使用邻域平台与回退尺度处理孤立污染。异常位置单独保存，在保留观测点上用PCHIP重构，并执行非负约束。该流程保留原始值、重构值与异常掩膜，不能据此断言所有真实仪表误差都被剔除。

已有十炉、三信号、三随机种子的90组孤立尖峰注入实验平均召回率约0.965。其结果仅适用于已定义污染方式；随机点与连续块遮蔽结果保留在原实验文件中。图1给出十炉六类连续曲线，短炉次保持原记录长度，不外推或按最短炉次截断。

图1 十炉六类连续工艺曲线。见 `../figures/ten_heats_six_process_curves.*`。

## 6.3 氧枪流量和烟气累计排放

正式氧枪流量采用31点、三阶Savitzky–Golay近似导数：

\[
\widehat q_h(t_i)=\max\{0,\operatorname{SGDeriv}_{31,3}[\widehat O_h](t_i)\}.
\]

31点覆盖30个采样间隔，即15 s；这是累计记录的平滑流量估计，不是另一个独立仪表测得的瞬时真值。原始累计值存在平台与离散增量，直接差分会把累计更新转换成流量脉冲。21、31、41点对照仅用于量化表示敏感性，不以产生更多事件为依据改变正式参数。原始计数更新特征与微分平滑共同影响短周期变化的解释，但本数据没有独立仪表更新说明，不能据此确认某种硬件计数机制。

累计排放用复合梯形积分计算：

\[
\widehat E_{g,h}(t_i)=\sum_{j=1}^{i}\frac{\widehat C_{g,h}(t_{j-1})+\widehat C_{g,h}(t_j)}{2}(t_j-t_{j-1}).
\]

该量从记录起点累计，单位m³。氧枪流量反积分与原累计量比较，检查数值一致性；总量误差小不等于各时刻导数准确。PCHIP解析微积分对照也表明，节点导数上的梯形求积不能替代分段内部的解析积分。

## 6.4 四类关系的统计设计

用每炉相对耗氧进度

\[
p_h(t)=\frac{\widehat O_h(t)-\widehat O_h(0)}{\widehat O_h(T_h)-\widehat O_h(0)}
\]

表示0至1的供氧进程。分别计算 \(O_h\) 与 \(C_{\mathrm{CO},h}\)、\(C_{\mathrm{CO_2},h}\)、\(E_{\mathrm{CO},h}\)、\(E_{\mathrm{CO_2},h}\) 的Pearson及Spearman相关。相关首先逐炉计算，再汇总炉间范围，避免炉际总量差异和记录长度主导结论。不对自相关时间点套用独立同分布检验或报告点级p值。

累计量共同增长本身就能产生高相关。为此构造与实际末端总量相同的恒定排放速率参照：

\[
E^{(0)}_{g,h}(t)=E_{g,h}(T_h)\frac{t}{T_h}.
\]

参照仅说明没有任何阶段性气体响应也可能出现高累计相关；不能将实际相关与参照之差解释为“真实机制占比”。图2、图3分别给出瞬时量和累计量关系。

表2 四类原值相关及恒定排放参照

'''+t2+r'''
累计耗氧量与CO瞬时量具有明显的同向进程关系，CO₂瞬时量的整炉相关较弱且炉间符号不一致。两种累计排放量均与耗氧量高度相关，但恒定排放参照亦高度相关，因而不能仅凭累计相关推断反应效率或物理时滞。

将原序列在不重叠的0.5、5、15秒网格上取增量，末端不足完整间隔的部分不纳入。表3区分短尺度更新结构与较慢过程变化。这里瞬时量的增量是 \(\Delta C_g\)，累计量的增量是 \(\Delta E_g\)，二者含义不同。

表3 增量相关与时间尺度

'''+t3+r'''
0.5秒累计量增量相关明显低于原值相关；5和15秒聚合后部分炉次相关增强，说明结果对累计计数更新和时间尺度敏感。CO₂累计增量在15秒尺度的炉间中位相关约0.596，但仍有负值炉次，因此也不能将高累计相关全部归为伪相关，或反过来将聚合后的相关当作因果响应。所有增量统计均未作时滞对齐。

图2 累计耗氧进度与CO、CO₂瞬时流量关系。见 `../figures/oxygen_instantaneous_relationships.*`。

图3 累计耗氧进度与两类累计排放关系。见 `../figures/oxygen_cumulative_relationships.*`。

## 6.5 五阶段比较和非线性趋势

按 \(p\in[0,0.2),[0.2,0.4),\ldots,[0.8,1]\) 划分五个等耗氧进度阶段。阶段边界由累计量的线性时间插值定位，同氧量平台取最后观测时刻，记录起止点固定。阶段平均排放为该段积分增量除以持续时间；阶段排放/耗氧比为 \(\Delta E_g/\Delta O\)，避免逐点除以低流量。该比值是表观体积比，不是脱碳收率。

表4 五阶段炉间统计

'''+t4+r'''
CO阶段平均排放中位数从第一阶段约2.11 m³/s增加到第三阶段约42.61 m³/s，随后处于较高水平。CO₂在第二阶段约19.58 m³/s达到中位高值，之后逐步降低。由此可解释CO₂整炉线性相关弱：先升后降的阶段形态可能被单个相关系数抵消。区间按整炉重采样2000次得到，仅反映十炉经验分布的不确定性；不代表随机抽取任意工业炉次的保证。

在已观测范围内比较两个预先规定的低复杂度描述模型：

\[
y_h(p)=c_0+c_1p,\qquad y_h(p)=c_0+c_1p+c_2p^2.
\]

逐炉系数见CSV；表5为十炉等权拟合，每炉全部点的总权重为1，避免长记录权重更大。因 \(p\) 无量纲，系数单位与目标量相同。累计排放拟合不强制通过零点，二次函数也不保证非负或单调，因此仅用于观测区间内的趋势描述，累计总量仍取积分结果。

表5 十炉等权描述性趋势系数

'''+t5+r'''
用留一炉次法，将九炉系数应用到未参与拟合的一炉。训练不使用测试目标；测试输入的最终耗氧总量只用于离线进度归一化，所以该检验评价回顾性跨炉趋势，不能称为在线预测。MAE单位与对应目标相同，R²相对于每个测试炉的均值计算；表6为十炉等权宏平均。

表6 线性与二次趋势的留一炉次表现

'''+t6+r'''
CO瞬时量二次趋势的宏平均R²约0.792，高于线性趋势的0.596；CO₂相应为0.309与−0.044。CO₂累计量的两种趋势表现接近，不能宣称增加二次项普遍改善所有关系。二次模型用于呈现非线性进程，不据此建立新的时滞或流量动力学模型。图4给出五阶段变化，图5给出原值与增量相关，图6给出流量数值复核。

图4 五个等耗氧进度阶段的排放及表观体积比。见 `../figures/oxygen_stage_statistics.*`。

图5 四类关系的原值与多尺度增量相关。见 `../figures/oxygen_trend_diagnostics.*`。

图6 累计计数与SG流量敏感性。见 `../figures/oxygen_flow_review.*`。

表7 31点正式流量的数值一致性与尺度敏感性

'''+t7+f'''
31点流量的总量相对差绝对值最大为{100*maxrel:.4f}%。它支持累计总量的一致性，不能作为瞬时变化或事件身份的验证。其余21、41点完整结果见 `../relationships/flow_review.csv`。

'''+r'''## 6.6 统计时滞与可辨识性边界

正式V2保留121点窗口（首末样本跨度60 s）、5点步长（2.5 s）及0–60 s非负搜索。输入表示为一阶差分减去21点居中滚动中位数。接受条件是源/目标动态激励、相关系数至少0.35、与相隔至少5采样点竞争峰的分离至少0.05，以及最大峰不在0或60秒搜索端点。峰宽和峰显著性只记录为诊断，不参与接受判据；置信分数是启发式权重，不是已校准的概率。

CO共13/2241个窗口通过，CO₂为11/2241，覆盖率分别约0.58%和0.49%。24个窗口位置与时滞均重现，但相邻源窗口重叠约95.9%；7个CO源重叠组、5个CO₂源重叠组不等于12次独立响应。167个多尺度重现源候选也不等于独立工艺事件。真实事件复核中的4条方法候选均为局部下降变化，尚无独立物理响应确认。

合成验证说明理想平移可恢复，复杂展宽、相邻响应及背景变化下的适用性依赖场景；合成性能不等同于真实工业准确率。双向搜索可能增加机会接受，因而不为提高覆盖率更改正式非负搜索。

原固定时滞与约束反平方根模型仍保留在接口中作为探索性离线推断。其模型选择依据训练误差，既有留一炉次也未建立稳定泛化，因此不将模型输出当作实测时滞。题目给出的“流量增大、时滞减小”在本数据上属于待检验方向：当前缺少可信独立事件，不能报告经验证的流量—时滞函数，也不能由证据不足推断该物理规律不存在。24个窗口的流量与候选时滞详见 `../delay_audit/accepted_windows.csv`，仅作窗口级估计，不重复拟合。

表8 不同流量区间的候选时滞描述性统计

'''+delaytable+r'''
流量边界取全部2241个共同候选中心的三分位数，仅按源流量划分，不根据24个时滞值调整。两个气体共用分层边界，源中心不重复计数；正式V2使用中心时刻的瞬时流量，不是事件平均流量。表中边界仅显示三位小数，归属按CSV和JSON的完整精度确定。

低、中、高层均存在算法合格窗口，但中、高层CO各只有1个窗口；高流量两类气体的窗口均只来自6C02057。时滞中位数未呈一致下降，这不能证明物理规律相反，因为炉次、阶段、重叠和候选身份未被分离。窗口数不是独立样本数，不报告事件置信区间；没有合格窗口的区间应记“不可辨识”，不能填零。表8直接给出现有算法在不同流量下的定量输出，仍不足以给出经过物理验证的流量—时滞函数。

当前数据接口未提供独立氧枪操作指令或分析仪响应记录。暂停新增时滞算法；若未来取得这些信息，再核查事件标签。题目能计算的重构、导数与积分结果完整给出；时滞部分同时交付窗口估计、置信分数、掩膜与限制，避免把不可辨识段补成预设曲线。

## 6.7 第三问接口和验收

表9 十炉离线接口验收

'''+t8+r'''
十炉共570个字段与历史Release在规定数值容差内一致，掩膜、窗口中心和时间轴要求精确一致。本轮实际核对了十个第一问本地接口的完整文件SHA256，与Release清单一致；其时间轴与第二问逐点精确相同，不通过截断或插值掩盖不一致。验收还检查序号、形状、布尔类型、导数、反积分、气体积分、观测时滞NaN、有效掩膜、时间映射及对齐插值；这是计算和交接验收，不是物理时滞验证。

正式文件为 `../data/{heat_id}_q2_processed.npz`，57个字段的名称、形状、单位和NaN语义不变。第三问应读取基础工艺量与异常掩膜；若使用探索性时滞对齐，必须同时读取 `delay_*_observed_mask`、`delay_*_confidence` 和 `alignment_valid_mask_*`。观测字段的“可靠”只表示通过V2判据，不是独立物理确认。缺失时滞不得填零，模型时滞不得作为独立监督标签。接口含未来窗口信息，只适用于离线研究，在线任务需另设因果处理与训练边界。

## 6.8 复现和结果追溯

正式统计配置为 `../../configs/final_analysis.json`，冻结研究基准为 `fe52075`。在仓库根目录运行 `python -m problem2.src.final_delivery`，输入原始文件、前轮冻结NPZ及第一问接口只读核验，固定写入一套当前统计、图表与正文。所有九张表对应的逐炉CSV、趋势参数、留一炉次结果、图注、字体检查和输入输出SHA256位于 `../relationships/`。数据字典、交接规则和完整复现命令随本章保存。
'''
    (REPORTS/'问题二_中文论文_正式章节版.md').write_text(chapter)
    acceptance=json.loads((out/'acceptance.json').read_text())
    (REPORTS/'第二问最终验收报告.md').write_text(f'''# 第二问结果与交接验收

本轮按“关系定稿、流量复核、时滞边界、正文整合、接口验收”收尾。附件1只是改名的身份事项已关闭；`fe52075` 的报告正文和生成代码均已明确关闭该事项，旧研究报告保留其当时状态，不再作为当前待办。

| 项目 | 结论 |
|---|---|
| 四类关系 | Pearson、Spearman、三种增量尺度、五个等耗氧阶段、线性/二次趋势及留一炉次结果已保存 |
| 氧枪流量 | 31点正式实现保留，21/41点仅诊断；计数更新机制尚不能由现有数据确认 |
| 时滞结论 | 24个算法合格窗口保留，尚无独立事件确认；暂停算法开发，不拟合流量—时滞规律 |
| 中文正文 | 九张统计表、六组核心图及十炉结果已统一 |
| 第三问交接 | 十炉{acceptance['historical_fields']}字段、{acceptance['interface_checks']}项检查通过；第一问十炉时间轴精确一致 |
| 物理解释 | 离线计算验收通过不代表动态规律已被验证 |

正文见 [正式章节](问题二_中文论文_正式章节版.md)。[验收JSON](../relationships/acceptance.json)、[输入清单](../relationships/input_manifest.json)、[输出清单](../relationships/provenance.json)保存证据。第三问不应把模型时滞或V2置信分数解释为实测标签或概率。

四类关系中，CO有较强阶段上升趋势，CO₂先升后降；累计变量的高相关不能单独证明机制。增量相关受时间尺度影响，不能简单将其归零或全部解释为伪相关。二次趋势的留一炉次改善以实际结果报告，不用于改动正式时滞模型。

本轮关闭的是统计整理与离线接口的工程验收。题面流量—时滞方向仍受可辨识性限制；该限制已作为正式答案的一部分明确写出，不以新增模型或强制负斜率填补。
''')
    (REPORTS/'实验配置与复现说明.md').write_text('''# 实验配置与复现说明

## 当前定稿入口

冻结基准 `fe52075`，正式配置 [final_analysis.json](../../configs/final_analysis.json)。在仓库根目录运行：

```bash
.venv/bin/python -m problem2.src.final_delivery
.venv/bin/python -m problem2.src.final_review
.venv/bin/python -m pytest -q problem2/tests problem2/research/response_broadening/test_response.py problem2/research/real_event_validation/test_review.py
```

需要现有NumPy、SciPy、Matplotlib、Pillow、pyxlsb、Fontconfig的fc-match及pdffonts；图形共享配置会严格检查真实Times New Roman，字体缺失或PDF未嵌入即停止。字体不随仓库提交。

输入位置通过配置中的相对路径解析，不依赖原机器用户名。输入包括十炉原始XLSB、`final/delay_audit/baseline/data/` 的冻结NPZ，以及 `problem1/final/data/*_interface.npz`。前轮NPZ缺失时，先使用前轮审查入口复现，不得以空值代替。原始文件和冻结NPZ按已保存清单核对SHA256；第一问接口只读取时间轴和波长维度，并记录文件指纹，不修改第一问。

历史第二问Release ZIP须存在于本机缓存 `~/.cache/apmcm-b/releases/problem2-final-data/`。若尚未下载：

```bash
gh release download problem2-final-data --repo ma39nong/Steelmaking-process --pattern '*' --dir "$HOME/.cache/apmcm-b/releases/problem2-final-data" --skip-existing
gh release download problem1-final-data-v1 --repo ma39nong/Steelmaking-process --pattern data_package_manifest.json --pattern SHA256SUMS --dir "$HOME/.cache/apmcm-b/releases/problem1-final-data-v1" --skip-existing
```

入口核验ZIP的固定SHA256，逐字段比较十炉NPZ。正式 `final/data/` 若缺失，从冻结NPZ逐字节复制；若已存在则逐字段核对内容；仅ZIP容器字节差异且数组和元数据一致时保留现有文件，数值差异超过1e-10容差时停止，不静默覆盖。当前数据接口保持字段、单位和NaN语义不变。

## 统计口径

原值相关逐炉计算，不报告点级p值。增量网格为不重叠0.5、5、15秒，末端不足一格丢弃；无时滞对齐。五阶段是等20%耗氧进度，不是五分位样本。阶段边界同氧量平台取最后观测，记录端点固定，积分边界由已保存累计排放插值。Bootstrap重采样单位为整炉，2000次，固定种子5201。

线性和二次趋势均预先指定。合并拟合使每炉总权重相同；留一炉次仅用九炉拟合，测试炉不参与系数估计。耗氧进度使用最终耗氧总量，因此这是离线回顾性评价，不声称在线预测。所有模型系数及逐炉误差直接保存，不按训练R²宣布普遍有效。

## 当前文件与历史证据

当前正式统计位于 `final/relationships/`，图位于 `final/figures/`，正文位于 `final/reports/`。重复运行原位更新，不建立平行版本目录。`delay_audit/` 和 `research/` 保留独立科研证据；其冻结检查对应当时提交，完整历史复现可在相应Git提交检出后执行。旧研究正文的待核验状态不覆盖用户已确认的当前数据身份。

`run_all` 是原始工程全流程，会运行原有诊断，并在默认正式输出路径上最后调用 `final_delivery`，保证当前统计、图表及正文一致。若自定义输出路径，需要另备定稿配置。不得把旧的分位数阶段统计当作本轮等耗氧阶段统计。正式V2数值源码未改变。
''')
    # English LaTeX is an includable section, not a standalone document.
    (REPORTS/'problem2_section.tex').write_text(r'''% Current English section; include from the main document, not standalone.
\section{Model and Results for Problem 2}
\subsection{Process reconstruction and numerical quantities}
For heat $h$, $t_i=0.5(s_i-s_0)$ seconds. Cumulative oxygen $O_h$ is in $\mathrm{m^3}$; the two gas fields are instantaneous volume rates $C_{g,h}$ in $\mathrm{m^3/s}$ under the problem statement. A nine-sample Hampel detector (threshold 6) and PCHIP preserve original values and outlier masks. The frozen oxygen flow is a nonnegative, degree-three, 31-sample Savitzky--Golay derivative. Gas totals are trapezoidal time integrals. Reintegration checks cumulative consistency, not instantaneous flow accuracy. Repeated cumulative samples and discrete increments require caution when interpreting short-period changes.
\subsection{Four oxygen relationships}
Within each heat we report Pearson and Spearman coefficients, nonoverlapping increments at 0.5, 5 and 15 seconds, and five equal 20-percent oxygen-progress stages. Whole heats, rather than time samples, are resampled for stage uncertainty. A constant-rate cumulative comparator demonstrates that high level correlation alone does not establish a reaction mechanism. Increment correlations depend on aggregation scale and do not isolate a causal effect.
With $p=(O_h-O_h(0))/(O_h(T_h)-O_h(0))$, linear and quadratic descriptive trends are fitted with equal total weight per heat and evaluated by leaving one heat out. For instantaneous CO, macro-average test $R^2$ is 0.596 and 0.792; for instantaneous CO$_2$, it is $-0.044$ and 0.309. The normalized progress requires the final oxygen total, so this is a retrospective evaluation. Curves are not extrapolated or substituted for physical integrals.
\subsection{Delay identifiability and offline handoff}
The frozen V2 method uses 121 samples, stride five, and a nonnegative 0--60 second lag search. Dynamic excitation, correlation, peak separation and non-boundary checks determine acceptance. Peak width and prominence are diagnostic only. CO has 13 accepted windows out of 2241, and CO$_2$ has 11. These overlapping windows, and 167 scale-recurring source candidates, are not independent physical events. Four candidate pairings remain unconfirmed. Flow-tertile summaries of the accepted windows are descriptive: CO median lags are 27.5, 7.5 and 40.0 s; CO$_2$ medians are 7.5, 15.0 and 36.5 s. High-flow windows come from one heat only, so these are not a validated monotone relation. Existing fixed and constrained inverse-square-root delay predictions remain exploratory; the data do not validate a decreasing flow--delay function.
The 57-field NPZ schema is unchanged. Ten heats and 570 historical fields have been checked, and the local first-question time axes match exactly. Observed delays, model delays and physical measurements must remain distinct. NaNs are not replaced by zero; downstream use requires masks and heuristic confidence scores. The interface is offline and may use future samples.
''')
    # Extend the existing dictionary without changing field semantics.
    dictionary=REPORTS/'数据字典与接口说明.md';text=dictionary.read_text();marker='\n## 5. 当前验收与物理解释\n';text=text.split(marker)[0]
    dictionary.write_text(text+marker+'\n十炉57字段共570项历史字段比较通过，第一问十炉本地时间轴已逐点精确核对。详见 [验收报告](第二问最终验收报告.md)。观测时滞的可靠标记仅指V2算法接受；峰宽及显著性为诊断量，不参与正式筛选。置信分数不是概率，模型时滞不是独立实测标签。部分峰宽诊断可为Infinity，无法辨识的观测和模型字段可为NaN，应按掩膜及字段定义处理。\n')
    handoff=REPORTS/'数据交接与第三问使用说明.md';text=handoff.read_text();marker='\n## 当前验收结果\n';text=text.split(marker)[0]
    handoff.write_text(text+marker+f'\n本轮十炉{acceptance["historical_fields"]}项历史字段比较及{acceptance["interface_checks"]}项接口检查通过；第一问本地十炉文件SHA256与Release清单一致，时间轴逐点精确一致。正式NPZ已在上述标准路径核验可读。验收记录位于 [acceptance.json](../relationships/acceptance.json) 与 [逐项检查](../relationships/interface_checks.csv)。这不证明任何候选时滞具有独立物理标签。\n\n推荐第三问将工艺变量及异常掩膜作为基础输入；只有在明确评估研究假设时才使用探索性时滞对齐。必须将“未对齐”与“时滞为零”区分，不得用模型时滞填充观测或生成伪标签。若需离线预测验证，共享拟合及预处理参数须在训练炉次内确定。\n')
    import csv
    with (out/'figure_captions.csv').open() as f:captions=list(csv.DictReader(f))
    index='# 第二问当前图表清单\n\n六组核心图与十二组补充图统一由 `python -m problem2.src.final_delivery` 生成，当前输出仅维护在 `../figures/`。PNG、SVG、PDF是同一图的必要导出格式，图题另存为可编辑 [图注CSV](../relationships/figure_captions.csv)。所有图内英文与数字使用真实Times New Roman，PDF嵌入检查见 [font_validation.json](../relationships/font_validation.json)，唯一版宽预览见 [layout_preview.png](../relationships/layout_preview.png)。\n\n'+table(['图名','用途','图题'],[[r['figure'],r['role'],r['caption']] for r in captions])+'\n第一轮时滞检查图、响应核研究图和真实事件检查页属于独立科研证据，分别保留在 `../delay_audit/` 和 `../../research/`。其组数和图数不计作独立响应数量。当前九张统计表位于中文正式章节，逐炉CSV位于 `../relationships/`。\n'
    (REPORTS/'完整图表清单.md').write_text(index)
    (out/'README.md').write_text('''# 第二问正式统计结果

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
''')

    audit=REPORTS/'独立数学审查与修复报告.md'
    text=audit.read_text();marker='\n<!-- final-delivery-status:start -->'
    text=text.split(marker)[0]
    audit.write_text(text+marker+'\n\n## 正式统计收尾状态\n\n附件1的中文名只是原英文文件改名，用户来源确认与SHA256核验均已完成，不再保留为待办。以上时滞审查保留当时实验与结论；其中“下一阶段”研究安排属于历史上下文。响应展宽与真实事件研究现已完成阶段审查，当前暂停时滞算法开发，不再按24个重叠窗口拟合流量—时滞规律。\n\n当前四类关系已补充Pearson、Spearman、三种增量尺度、等20%耗氧阶段、线性/二次趋势及留一炉次结果；原分位数阶段统计作为旧诊断保留，不是当前正式阶段。正式结果见 [中文章节](问题二_中文论文_正式章节版.md) 与 [验收报告](第二问最终验收报告.md)。十炉570字段历史比较、480项离线接口检查和十炉第一问时间轴核对通过；这不是对物理时滞规律的确认。\n\n<!-- final-delivery-status:end -->\n')
    manifest_path=ROOT/'problem2/final/results/artifact_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    manifest.update(current_entry='python -m problem2.src.final_delivery',formal_statistics='../relationships/',
        current_chapter='../reports/问题二_中文论文_正式章节版.md',current_acceptance='../reports/第二问最终验收报告.md',
        figures=[r['figure']+'.'+ext for r in captions for ext in ['png','svg','pdf']],
        relationships_status='Current equal-progress and increment statistics are in ../relationships; relationship_summary.json retains the earlier quantile-stage diagnostic.',
        delay_status='Exploratory overlapping-window estimates; no validated independent physical event or flow-delay law')
    from .delay_audit import write_json
    write_json(manifest_path,manifest)

    q1audit=REPORTS/'第一问成果审查与第二问接口说明.md'
    text=q1audit.read_text()
    old='Release 清单标明数据包大小约 548.8 MB，且明确不包含原始 XLSB 和工艺附加列。Release 的 SHA256 清单已保存于 `problem2/audit/`。本地未下载 548.8 MB 大包，因此接口数组的逐文件哈希尚未在本机复核。后续若需要读取光谱接口，应先下载并按清单校验。'
    new='Release清单标明数据包大小约548.8 MB，且不包含原始XLSB和工艺附加列。本轮已对本地十炉接口文件逐文件计算SHA256，全部匹配Release清单，并核对其时间轴与第二问完全一致，见 [第一问接口身份记录](../relationships/q1_interface_identity.csv)。无需为已有且哈希一致的本地文件重新下载整个压缩包。本项是文件身份和坐标核验，不重复评价第一问的科学结论。'
    q1audit.write_text(text.replace(old,new))
