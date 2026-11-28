"""Competition manuscript: model, evidence, answer. No numerical pipeline changes."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import pandas as pd
from .delay_audit import write_csv

ROOT=Path(__file__).resolve().parents[2]
REPORTS=ROOT/'problem2/final/reports'


def research_evidence():
    folder=ROOT/'problem2/research/signal_fidelity/results'
    gas=pd.read_csv(folder/'gas_fidelity.csv');source=pd.read_csv(folder/'source_fidelity.csv');control=pd.read_csv(folder/'negative_controls.csv')
    evidence=[]
    def value(name,frame,field,condition,file):
        v=float(frame[field].median());path=folder/file
        evidence.append(dict(metric=name,aggregation='median',value=v,conditions=len(frame),condition=condition,
                             source=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        return v
    result={}
    for g in ['co','co2']:
        for method in ['acquisition_only','formal']:
            q=gas[(gas.domain=='raw_background')&(gas.gas==g)&(gas.duration_s==2)&(gas.amplitude_multiplier==.5)&(gas.acquisition=='grid_0.01')&(gas.method==method)]
            result[g+'_'+method]=value(g+'_'+method,q,'peak_ratio','2 s triangle; 0.5 whole-record gas SD amplitude; raw heat backgrounds; stored grid 0.01','gas_fidelity.csv')
    for method in ['formal','SG11']:
        q=source[(source.domain=='raw_background')&(source.duration_s==4)&(source.acquisition=='grid_0.01')&(source.method==method)]
        result[method+'_peak']=value(method+'_peak',q,'peak_ratio','4 s triangle; median background flow amplitude; integrate flow before cumulative acquisition; stored grid 0.01','source_fidelity.csv')
        q=control[(control.signal=='oxygen')&(control.acquisition=='empirical_hold')&(control.method==method)]
        result[method+'_ripple']=value(method+'_ripple',q,'constant_flow_ripple_cv','constant true flow 5 m3/s; ten recorded nonzero-update schedules as competing hold hypotheses','negative_controls.csv')
    spike=control[(control.signal.isin(['co','co2']))&(control.acquisition=='grid_0.01')]
    for method in ['formal','hampel_window5']:
        q=spike[spike.method==method];result[method+'_spike_detected']=int(q.corrupted_point_flagged.sum());result[method+'_spike_conditions']=len(q)
        evidence.append(dict(metric=method+'_spike_detected',aggregation='count',value=result[method+'_spike_detected'],conditions=len(q),
                             condition='one 10-SD gas logging spike per heat and gas; stored grid 0.01',source=str((folder/'negative_controls.csv').relative_to(ROOT)),sha256=hashlib.sha256((folder/'negative_controls.csv').read_bytes()).hexdigest()))
    path=ROOT/'problem2/final/results/raw_injection_validation.json';r=json.loads(path.read_text());rows=[v for h in r for v in h['rows']]
    result['spike_recall']=float(np.mean([v['recall'] for v in rows]));assert len(rows)==90
    evidence.append(dict(metric='historical_spike_recall',aggregation='mean',value=result['spike_recall'],conditions=len(rows),condition='10 heats x 3 signals x 3 seeds; 20 spikes per trial',source=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    path=ROOT/'problem2/research/measurement_chain/results/background_trials.csv'
    chain=pd.read_csv(path);selected=chain[(chain.scenario=='known_addition')&(chain.duration_s==4)]
    for method in ['direct_flow','raw_counter_chain']:
        q=selected[selected.method==method]
        assert len(q)==1200 and np.array_equal(q.correct.to_numpy(),(q.accepted & (q.error_s.abs()<=1)).to_numpy())
        result[method+'_correct']=int(q.correct.sum());result[method+'_conditions']=len(q)
        evidence.append(dict(metric=method+'_4s_delay_correct',aggregation='count',value=result[method+'_correct'],conditions=len(q),
                             condition='4 s triangular paired addition; same gas injection and cleaning; accepted and absolute delay error <=1 s; conditions are not independent physical events',
                             source=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return result,evidence


def render(a,cfg,t):
    from .final_text import table
    e,evidence=research_evidence()
    compare=table(['导数估计','4 s供氧脉冲峰值保留中位数 / %','恒定流量伪波动中位数 / %'],[
        ['31点SG',f"{100*e['formal_peak']:.1f}",f"{100*e['formal_ripple']:.1f}"],['11点SG',f"{100*e['SG11_peak']:.1f}",f"{100*e['SG11_ripple']:.1f}"]])
    chapter=r'''# 基于稳健重构的转炉工艺量计算与烟气响应分析

针对累计耗氧量与烟气观测中的异常和离散更新，采用稳健重构、平滑微分及数值积分建立六类工艺曲线，并从耗氧进程和局部响应两个尺度分析供氧与排放的关系。

## 6.1 工艺变量与连续信号重构

采用附件1至附件10提供的十炉炼钢工艺数据，以0.5 s为时间采样间隔，提取累计耗氧量、CO和CO₂瞬时体积流量三个变量。对第 \(h\) 炉，令 \(t_i=0.5(s_i-s_0)\)，以 \(O_h(t)\) 表示累计耗氧量，\(C_{g,h}(t)\) 表示气体 \(g\in\{\mathrm{CO},\mathrm{CO_2}\}\) 的瞬时体积流量，单位分别为m³和m³/s。十炉记录时长为387–777 s，分别在各自观测区间内重构。

表1 十炉数据规模及耗氧、排放总量

'''+t['t1']+r'''
对每个工艺序列 \(x_i\)，在9点邻域 \(\mathcal N_i\) 中计算中位数与稳健尺度：

\[
m_i=\operatorname{median}_{j\in\mathcal N_i}x_j,\qquad
\sigma_i=1.4826\operatorname{median}_{j\in\mathcal N_i}|x_j-m_i|.
\]

当 \(|x_i-m_i|>6\sigma_i\) 时标记异常；局部尺度为零时，结合两侧平台和邻域尺度识别孤立尖峰，非有限值及负值同时剔除。利用其余观测构造PCHIP连续函数 \(\widehat O_h(t)\) 和 \(\widehat C_{g,h}(t)\)。该插值保持节点取值及局部单调形态，避免普通三次样条在平台和突变附近产生额外过冲，并对重构值施加非负约束。

'''+f'''90组孤立尖峰注入试验的平均召回率为{100*e['spike_recall']:.1f}%。在2 s低幅烟气响应试验中，CO与CO₂的峰值保留中位数由采集后的{100*e['co_acquisition_only']:.1f}%和{100*e['co2_acquisition_only']:.1f}%降至清洗后的{100*e['co_formal']:.1f}%和{100*e['co2_formal']:.1f}%。缩小至5点窗口虽能改善短峰保留，但另一组20个已知强孤立尖峰仅检出{e['hampel_window5_spike_detected']}个，9点窗口全部检出。故保留9点方案用于连续曲线与累计排放，并在时滞分析中考察短响应损失；具体注入条件见补充材料。

'''+r'''## 6.2 氧枪流量与烟气排放量

氧枪流量由累计耗氧量的平滑导数估计：

\[
\widehat q_h(t_i)=\max\left\{0,\operatorname{SGDeriv}_{31,3}[\widehat O_h](t_i)\right\}.
\]

采用31点、三阶Savitzky–Golay局部多项式，窗口首末跨度为15 s。原始累计耗氧记录约一半相邻值重复，正增量中位数约4.73 m³；直接差分会将离散更新放大为流量脉冲。平滑窗口因而需要兼顾动态保留与更新波动抑制。表2比较31点与11点方案的两项受控结果。

表2 SG导数估计的短峰保留与伪波动权衡

'''+compare+r'''
表注：两项指标来自不同受控条件，分别衡量4 s脉冲峰值保留与恒定流量下的伪波动；后者为估计流量标准差与真流量之比。采集保持机制为实验假设，详细设置见补充材料。

'''+f'''31点SG的恒定流量伪波动中位数为{100*e['formal_ripple']:.1f}%，11点SG虽提高短峰保留率，却使伪波动增至{100*e['SG11_ripple']:.1f}%。综合流量稳定性与累计量一致性，选取31点SG进行流量重构，十炉反积分相对误差绝对值均小于0.08%。该估计用于连续流量与总量分析，不保证短瞬态时序准确。

'''+r'''烟气累计排放量定义为

\[
\widehat E_{g,h}(t)=\int_0^t\widehat C_{g,h}(u)\,\mathrm du,
\]

在采样节点采用复合梯形公式：

\[
\widehat E_{g,h}(t_i)=\sum_{j=1}^{i}
\frac{\widehat C_{g,h}(t_{j-1})+\widehat C_{g,h}(t_j)}{2}(t_j-t_{j-1}).
\]

累计量从记录起点计，单位为m³，十炉末端总量见表1。图1给出氧枪流量、累计耗氧量、CO与CO₂瞬时流量及累计排放量六类连续曲线；图2展示累计记录的更新结构和平滑导数的尺度差异。

![六类连续工艺曲线](../figures/ten_heats_six_process_curves.png)

图1 十炉氧枪流量、累计耗氧量与两种烟气的瞬时流量及累计排放量

![累计耗氧量与流量估计](../figures/oxygen_flow_review.png)

图2 累计耗氧记录的离散更新特征与不同SG窗口下的流量估计

## 6.3 累计耗氧量与四类烟气变量的关系

以累计耗氧量为自变量，分别研究CO、CO₂瞬时流量及其累计排放量。先逐炉计算Pearson与Spearman系数，再汇总炉间范围，避免记录长度和炉际总量差异主导统计。为识别累计变量共同增长的影响，构造末端总排放量相同的恒定速率参照

\[
E^{(0)}_{g,h}(t)=E_{g,h}(T_h)\frac{t}{T_h}.
\]

表3 四类耗氧关系的相关系数及恒定排放参照

'''+t['t2']+r'''
CO瞬时流量与累计耗氧量的Pearson系数为0.537–0.930，表现出明显的同向进程关系；CO₂相应为−0.199–0.249，单一线性关系较弱。两类累计排放与耗氧量高度相关，但恒定排放参照也达到0.997–1.000，表明高累计相关本身不足以区分共同进程与动态响应，也不能据此确定响应时滞。图3和图4分别展示瞬时流量与累计排放的关系。

![耗氧进度与瞬时烟气流量](../figures/oxygen_instantaneous_relationships.png)

图3 耗氧进度与CO、CO₂瞬时流量

![耗氧进度与累计排放](../figures/oxygen_cumulative_relationships.png)

图4 耗氧进度与CO、CO₂累计排放

进一步在不重叠的0.5、5和15 s区间计算 \(\Delta O\) 与 \(\Delta C_g\)、\(\Delta E_g\) 的相关。累计量在0.5 s尺度上的增量相关明显低于原值相关；聚合至15 s后，CO与CO₂累计增量的Pearson炉间中位数分别为0.356和0.596，但逐炉仍有明显差异（图6、附录表S1）。因此，排放与供氧的关联包含时间尺度差异，原值相关不足以描述局部动态。

为比较不同炉次的排放进程，定义无量纲耗氧进度

\[
p_h(t)=\frac{\widehat O_h(t)-\widehat O_h(0)}{\widehat O_h(T_h)-\widehat O_h(0)},
\]

将其划分为五个等耗氧阶段。每段以排放积分增量除以持续时间得到平均流量，以 \(\Delta E_g/\Delta O\) 得到排放与耗氧的表观体积比。阶段边界由累计耗氧量插值定位，边界处同氧量平台取最后时刻。

表4 五个等耗氧阶段的烟气排放统计

'''+t['t4']+r'''
CO阶段平均流量中位数由第一阶段2.11 m³/s升至第三阶段42.61 m³/s，随后保持高位；CO₂在第二阶段达到19.58 m³/s，之后下降。这一先升后降形态解释了CO₂整炉线性相关较弱的原因。表中区间以整炉为重采样单位进行2000次Bootstrap，避免把相邻时间点作为独立样本。

在 \(p\in[0,1]\) 内，用线性与二次多项式描述上述进程：

\[
f_1(p)=c_0+c_1p,\qquad f_2(p)=c_0+c_1p+c_2p^2.
\]

合并拟合时令每炉总权重相等，并以留一炉次法评价跨炉趋势。由于耗氧进度使用末端耗氧总量，该检验针对离线进程描述。系数见附录表S2，验证结果见表5。

表5 线性与二次耗氧趋势的留一炉次验证

'''+t['t6']+r'''
CO瞬时流量的二次趋势将宏平均MAE由8.923降至6.204 m³/s，R²由0.596提高至0.792；CO₂的R²由−0.044提高至0.309，正R²炉数由2增至9。二次项有效刻画了CO₂阶段转折。累计CO₂的两种趋势MAE接近，因此仍以积分结果确定累计排放，多项式只用于解释观测区间内的趋势。

![耗氧阶段统计](../figures/oxygen_stage_statistics.png)

图5 五个等耗氧阶段的烟气流量与排放耗氧体积比

![原值与增量相关](../figures/oxygen_trend_diagnostics.png)

图6 四类耗氧关系的原值与多尺度增量相关

## 6.4 不同氧枪流量下的烟气响应时滞

以供氧变化与随后烟气变化的局部对应关系估计时滞。对 \(\widehat q_h\) 和 \(\widehat C_{g,h}\) 分别作一阶差分，并减去21点居中滚动中位数，得到动态分量 \(u_h\) 与 \(v_{g,h}\)。在121点窗口 \(W_k\) 内计算

\[
\widehat\tau_{g,h,k}=\underset{\tau\in\{0,0.5,\ldots,60\}}{\arg\max}
\operatorname{corr}\!\left(u_h(t),v_{g,h}(t+\tau)\right)_{t\in W_k}.
\]

正迟延表示烟气变化晚于供氧变化。窗口首末跨度60 s，中心每2.5 s移动一次。要求源窗口与未平移烟气窗口的动态标准差均不低于各自全炉标准差的5%，并保留最大相关系数不低于0.35、与相隔至少2.5 s竞争峰的相关系数差不低于0.05、且最优值不在搜索端点的估计。峰宽与显著性用于诊断，不参与筛选。

十炉共有2241个候选中心，CO通过13个，CO₂通过11个，比例分别为0.58%和0.49%。以全部候选中心的氧枪流量三分位数划分低、中、高流量，分层边界约为4.652和5.199 m³/s，两种气体使用相同边界。表6给出各层时滞的中位数、范围与覆盖情况。

表6 不同氧枪流量下的候选烟气时滞

'''+t['delaytable']+r'''
表注：各项均为窗口级候选估计，尚未确认独立物理响应事件。

CO在低、中、高流量下的时滞中位数分别为27.5、7.5和40.0 s，CO₂分别为7.5、15.0和36.5 s。中、高流量的CO各仅有1个窗口；高流量两种气体均只覆盖同一炉。相邻窗口共享约95.9%的源数据，24个合格窗口对应的7个CO与5个CO₂源重叠组也不能作为12次独立响应。

'''+f'''在相同烟气注入与清洗条件下，4 s已知供氧脉冲直接叠加至流量序列时，时滞正确恢复数为{e['direct_flow_correct']}/{e['direct_flow_conditions']}；先积分成累计耗氧量并经过重构链后，降至{e['raw_counter_chain_correct']}/{e['raw_counter_chain_conditions']}，表明短时动态信息的衰减显著影响时滞辨识。结合表6，CO与CO₂候选时滞均未呈现一致的单调递减关系。因此保留各流量区间的统计时滞估计，尚不建立统一的流量—时滞函数；第三问采用原始时间轴建立预测基线，避免不可靠的时滞修正引入额外误差。
'''
    supplement=r'''# 问题二补充统计与验证

以下结果补充正文中的方法选择和关系分析。正文仅保留直接回答题目所需的模型与主要结果。

异常检测召回率来自十炉、三信号、三随机种子的90组孤立尖峰注入。短烟气响应试验向十炉背景注入宽2 s、峰高为原序列全程标准差0.5倍的三角响应，施加0.01 m³/s存储网格后比较清洗前后的峰值；强尖峰对照则向每炉、每种气体加入一个幅度为原序列标准差10倍的孤立记录尖峰。峰值保留统计为中位数，召回率统计为试验平均值。

正文表2的峰值保留试验在十炉原始背景上加入4 s三角流量脉冲，幅值取各炉流量中位数，先积分成累计量，再施加0.01 m³存储网格并完成重构。恒定流量对照以5 m³/s真流量和各炉记录的非零更新时刻构造保持输入，伪波动为估计流量标准差与真流量之比。记录更新时刻仅用于构造保持假设，不据此认定仪表更新机制。

正文6.4节的1200组配对条件覆盖十炉、两种气体、30个固定种子及15/30 s已知迟延。两入口使用相同供氧脉冲与烟气注入、清洗结果；正确恢复要求通过局部互相关筛选且迟延误差不超过1 s。累计量入口将已知流量先积分，再加到原始累计耗氧量并经过Hampel、PCHIP和31点SG。条件共享抽样背景，不代表1200个独立物理事件。完整数值见[测量链注入结果](../../research/measurement_chain/results/background_trials.csv)。

表S1 四类关系的多尺度增量相关

'''+t['t3']+r'''
表S2 十炉等权线性与二次趋势系数

'''+t['t5']+r'''
表S3 逐炉流量反积分一致性及窗口敏感性

'''+t['t7']+r'''
表S4 积分辨识对照的关键验证结果

| 验证条件 | 普通积分 | 弱形式积分 | 解释 |
|---|---|---|---|
| 混合宽度、独立种子加噪，正确恢复 | 44/60 | 60/60 | 同一一阶加纯迟延模型，误差≤1 s且通过筛选 |
| 4 s响应加噪，正确恢复 | 30/60 | 60/60 | 短响应点估计改善 |
| 30 s响应无噪声，正确恢复 | 6/6 | 5/6 | 弱形式1例剖面过宽被拒绝 |
| 4 s响应、迟延15 s、T=8 s，名义95%区间覆盖 | 0/10 | 5/10 | 短响应区间未得到可靠校准 |
| 两级响应无噪声，真纯迟延15 s | 16.5 s | 17.0 s | 模型失配下点估计仍偏离纯迟延 |
| 十炉真实数据，弱形式通过炉数 | — | CO：0/10；CO₂：0/10 | 未通过可辨识性和连续留出验收 |

弱形式改善受控噪声下的点估计，但未通过完整采用验收，故不列为正文采用的方法，也不用于流量—时滞回归。其置信区间以一阶模型与噪声假设为条件，不能解释为仪表延迟和工艺迟延的物理分离。该对照与正文累计排放的数值积分用途不同，不替换后者。

工程数据字典、接口检查和复现命令见[接口说明](数据字典与接口说明.md)、[交接说明](数据交接与第三问使用说明.md)及[复现说明](实验配置与复现说明.md)。完整受控条件见[信号保真验证](../../research/signal_fidelity/results/REPORT.md)与[弱形式验证](../../research/weak_integral/results/REPORT.md)。
'''
    return chapter,supplement,evidence


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--review',action='store_true');args=parser.parse_args()
    if args.review:
        from .final_manuscript_review import run
        run();return
    from .final_text import write_reports
    cfg=json.loads((ROOT/'problem2/configs/final_analysis.json').read_text());a={}
    for path in (ROOT/cfg['output_dir']).glob('*.csv'):
        def convert(v):
            if v=='True':return True
            if v=='False':return False
            if v=='':return None
            try:return float(v) if any(c in v for c in '.eE') else int(v)
            except ValueError:return v
        with path.open() as f:a[path.stem]=[{k:convert(v) for k,v in row.items()} for row in csv.DictReader(f)]
    write_reports([],a,cfg,manuscript_only=True)




def english_section():
    return r'''% Includable section, not a standalone document.
\section{Process Reconstruction and Flue-Gas Response}
\subsection{Process variables and continuous reconstruction}
We use ten heats sampled every 0.5 s. Cumulative oxygen $O_h$ is measured in $\mathrm{m^3}$ and the CO and CO$_2$ fields represent instantaneous volume rates $C_{g,h}$ in $\mathrm{m^3/s}$. A nine-point Hampel detector with threshold six removes isolated outliers, followed by nonnegative PCHIP reconstruction. Mean recall in 90 specified spike-injection trials is 96.5\%. For 2 s triangular gas responses of amplitude one-half of the record standard deviation, median peak retention after the stored-value grid is 86.9\% for CO and 87.5\% for CO$_2$; after cleaning it is 28.6\% and 37.7\%. A five-point detector preserves more short responses but detects only 15 of 20 strong logging spikes, versus 20 for the nine-point method. We retain the latter for process curves and cumulative quantities.
\subsection{Flow and emission calculations}
Oxygen flow is estimated using a nonnegative 31-point cubic Savitzky--Golay derivative, spanning 15 s. In separate controlled conditions, shortening the window to 11 points raises median 4 s pulse retention from 35.0\% to 77.3\%, but increases spurious variation for constant flow under hypothesized recorded holds from 1.6\% to 11.8\%. Reintegration errors for all ten heats are below 0.08\%. This supports total-volume consistency, while short-transient timing remains affected by smoothing. Gas emissions are evaluated by trapezoidal integration of $C_{g,h}$ from the record start.
\subsection{Four oxygen--emission relationships}
Within-heat Pearson and Spearman coefficients, nonoverlapping increments, and five equal oxygen-progress stages distinguish overall progress from local dynamics. The normalized progress is $p_h=(O_h-O_h(0))/(O_h(T_h)-O_h(0))$. CO increases from a first-stage median of 2.11 to a third-stage median of 42.61 $\mathrm{m^3/s}$, then stays high. CO$_2$ peaks in stage two at 19.58 $\mathrm{m^3/s}$. Equal-heat-weight quadratic trends improve leave-one-heat-out $R^2$ from 0.596 to 0.792 for CO and from $-0.044$ to 0.309 for CO$_2$. A constant-emission comparator also has cumulative correlation above 0.997, so cumulative correlation alone does not identify a dynamic response. Progress uses the final oxygen total; validation therefore describes retrospective trends.
\subsection{Delay estimates at different oxygen flows}
Local correlations of detrended first differences use 121-point windows, 2.5 s center spacing, and a 0--60 s lag search. Dynamic excitation, correlation at least 0.35, peak separation at least 0.05, and an interior optimum define acceptance. CO has 13 accepted windows and CO$_2$ has 11 out of 2241 centers. Shared flow tertiles give CO median lags of 27.5, 7.5 and 40.0 s, and CO$_2$ medians of 7.5, 15.0 and 36.5 s. These are overlapping-window estimates: high-flow results come from one heat and the windows are not independent physical events. In paired 4 s pulse injections with the same gas input and cleaning, accepted delay recovery within 1 s falls from 1199/1200 for direct flow addition to 3/1200 after cumulative-oxygen reconstruction. Short-response attenuation therefore materially limits delay identification. The ten heats do not establish a reliable decreasing flow--delay function; subsequent analyses retain the original time axis where alignment is unsupported.
'''


if __name__=='__main__':main()
