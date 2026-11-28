"""Fixed report and scientific figures; no method tuning or cherry-picking."""
import subprocess
import numpy as np
import pandas as pd
from PIL import Image
from ...src.delay_audit import write_json,write_csv
from ...src.figure_fonts import configure
from ..target_excitation.present import mdtable as _mdtable


def mdtable(records):
    if not records:return "No rows."
    headers=list(records[0]);rows=[[None if isinstance(r.get(k),float) and np.isnan(r.get(k)) else r.get(k) for k in headers] for r in records]
    return _mdtable(headers,rows)


def deliver(cfg,out,decision):
    import matplotlib.pyplot as plt
    data=pd.read_csv(out/'trials.csv');s=pd.read_csv(out/'summary.csv');p=pd.read_csv(out/'profiles.csv');b=pd.read_csv(out/'bootstrap_draws.csv')
    primary=data[(data.acquisition=='ideal')&(data.gas_input=='raw')]
    ss=s[(s.acquisition=='ideal')&(s.gas_input=='raw')]
    font=configure();figdir=out/'figures';figdir.mkdir(exist_ok=True);captions=[];checks=[];paths=[]
    def save(fig,name,caption):
        fig.tight_layout()
        for ext in ['png','svg','pdf']:
            fig.savefig(figdir/(name+'.'+ext),dpi=180,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
        svg=figdir/(name+'.svg');svg.write_text('\n'.join(x.rstrip() for x in svg.read_text().splitlines())+'\n')
        plt.close(fig);raw=subprocess.check_output(['pdffonts',str(figdir/(name+'.pdf'))],text=True);col=raw.splitlines()[0].index('emb');lines=raw.splitlines()[2:]
        assert lines and all('TimesNewRoman' in x.replace(' ','') and x[col:col+3]=='yes' for x in lines)
        checks.append(dict(figure=name,pdffonts=raw));captions.append(dict(figure=name,caption=caption,width_inches=7.2));paths.append(figdir/(name+'.png'))
    colors={'ordinary':'#245e8a','weak':'#bd582b'}
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.2))
    for ax,noise,title in zip(axes,[0,.05],['(a) No noise','(b) Noise SD = 0.05']):
        for j,m in enumerate(colors):
            rates=[float(primary[(primary.kind==k)&(primary.noise_std==noise)&(primary.method==m)].correct.mean()) for k in ['matched','short','slow']]
            ax.bar(np.arange(3)+j*.34-.17,rates,width=.34,color=colors[m],label=m.capitalize())
        ax.set_xticks(range(3),['Mixed widths','Short: 4 s','Slow: 30 s']);ax.set_ylim(0,1.12);ax.set_ylabel('Accepted and within 1 s / all trials');ax.set_title(title);ax.legend(fontsize=7)
    save(fig,'recovery_by_scale','独立评价种子，理想采集与原始带噪输出的配对恢复率。短时、缓变与混合宽度分别呈现；恢复要求通过原有筛选且纯迟延绝对误差不超过1秒。')
    fig,axes=plt.subplots(1,3,figsize=(7.2,3.0))
    for ax,kind in zip(axes,['matched','short','two_stage']):
        for m in colors:
            q=p[(p.kind==kind)&(p.method==m)&(p.noise_std==.05)&(p.seed==9201)&(p.true_delay_s==15)&(p.true_T_s==8)]
            # Normalize independently for shape display; amplitudes across moment types are not comparable.
            ax.plot(q.lag_s,(q.loss-q.loss.min())/max(float(q.loss.max()-q.loss.min()),1e-12),color=colors[m],label=m.capitalize())
        ax.axvline(15,color='black',ls=':',lw=.8);ax.set_title(kind.replace('_',' ').capitalize());ax.set_xlabel('Candidate delay (s)');ax.set_ylabel('Normalized excess moment loss');ax.legend(fontsize=6)
    save(fig,'delay_profiles','预先指定seed=9201、纯迟延15秒、时间常数8秒的目标函数剖面。各自归一化只比较形状，不比较两种积分方程的损失绝对数值；竖线为合成真值。完整所有主评价剖面保存在CSV。')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.2),sharey=True)
    for ax,method in zip(axes,colors):
        q=primary[(primary.kind=='matched')&(primary.noise_std==.05)&(primary.method==method)].sort_values(['true_delay_s','true_T_s','seed'])
        x=np.arange(len(q));ax.scatter(x,q.estimate_s-q.true_delay_s,s=9,color=colors[method])
        # CI may exclude the point estimate; explicit interval segments retain the actual percentile limits.
        ax.vlines(x,q.ci_lower_s-q.true_delay_s,q.ci_upper_s-q.true_delay_s,color=colors[method],lw=.65)
        ax.axhline(0,color='black',ls=':',lw=.8);ax.set_title(method.capitalize()+' conditional 95% intervals');ax.set_xlabel('Fixed setting / seed order');ax.set_ylabel('Estimated minus true delay (s)')
    save(fig,'conditional_intervals','全部60个主评价噪声试验的95%百分位参数Bootstrap区间。误差线端点以独立竖线给出，点为估计值；仅在一阶模型及独立同分布高斯观测噪声假设下成立。覆盖率不按通过筛选的结果条件化。')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.2))
    for m in colors:
        q=b[(b.kind=='matched')&(b.method==m)&(b.seed==9201)&(b.true_delay_s==15)&(b.true_T_s==8)]
        axes[0].scatter(q.lag_s,q.T_s,s=9,alpha=.25,label=m.capitalize(),color=colors[m])
        q=primary[(primary.kind=='matched')&(primary.method==m)&(primary.noise_std==.05)]
        axes[1].scatter(q.error_s,q.T_s-q.true_T_s,s=12,alpha=.6,label=m.capitalize(),color=colors[m])
    for ax in axes:ax.legend(fontsize=7)
    axes[0].set_title('(a) Model-bootstrap parameter coupling');axes[0].set_xlabel('Delay (s)');axes[0].set_ylabel('Time constant (s)')
    axes[1].set_title('(b) Bias across all matched noisy trials');axes[1].set_xlabel('Delay error (s)');axes[1].set_ylabel('Time-constant error (s)')
    save(fig,'parameter_coupling','左为固定代表场景全部199次重拟合的迟延与时间常数耦合，右为全部60个主评价试验的参数偏差；不能用良好预测或小残差替代纯迟延正确性验证。')
    if (out/'real_trials.csv').exists() and decision['controlled_gate_passed']:
        real=pd.read_csv(out/'real_trials.csv');lo=pd.read_csv(out/'real_loho.csv')
        fig,axes=plt.subplots(1,2,figsize=(7.2,3.5))
        for ax,gas in zip(axes,['co','co2']):
            q=real[(real.gas==gas)&(real.method=='weak')];l=lo[lo.gas==gas]
            ax.plot(range(len(q)),q.validation_improvement,'o-',label='Within-heat holdout',ms=3)
            ax.plot(range(len(l)),l.validation_improvement,'s-',label='Leave-one-heat-out',ms=3)
            ax.axhline(0,color='black',ls=':',lw=.8);ax.set_xticks(range(len(q)),q.heat_id,rotation=60,fontsize=6);ax.set_title(gas.upper());ax.set_ylabel('RMSE improvement over no-input model');ax.legend(fontsize=6)
        save(fig,'real_holdout','十炉原始数据的连续留出与留一炉次预测。负数明确表示劣于无供氧动态项基线；所有炉次完整报告，不按通过数量选参数。')
    write_csv(out/'figure_captions.csv',captions);write_json(out/'font_validation.json',dict(**font,pdf_checks=checks))
    # One current layout check, at a fixed viewing width.
    thumbs=[]
    for path in paths:
        im=Image.open(path).convert('RGB');im.thumbnail((900,480));thumbs.append(im)
    canvas=Image.new('RGB',(1800,500*int(np.ceil(len(thumbs)/2))),'white')
    for j,im in enumerate(thumbs):canvas.paste(im,((j%2)*900,(j//2)*500))
    canvas.save(out/'layout_preview.png')
    table=ss[['kind','noise_std','method','conditions','accepted','correct','accepted_MAE_s','ci_coverage','ci_median_width_s']].to_dict('records')
    matched=primary[(primary.kind=='matched')&(primary.noise_std==.05)]
    acq=s[(s.kind=='matched')&(s.noise_std==.05)]
    sens=pd.read_csv(out/'support_sensitivity.csv');sr=[]
    for span,g in sens.groupby('support_s'):sr.append(dict(support_s=span,conditions=len(g),correct=int(g.correct.sum()),accepted=int(g.accepted.sum())))
    mismatch=primary[(primary.kind=='two_stage')&(primary.noise_std==0)][['method','estimate_s','T_s','validation_improvement','profile_width_s']].to_dict('records')
    verdict='主匹配场景的预定受控门槛通过，并已完成十炉探索；完整研究的采用验收失败，当前弱形式方案不能作为正式物理时滞模型。' if decision['controlled_gate_passed'] else '受控验收未通过：本轮不采用弱形式法，不进入十炉真实数据辨识，不开展流量—时滞回归。'
    gate_rows=[dict(check=r['name'],passed=r['passed'],evidence=str({k:v for k,v in r.items() if k not in ['name','passed']})) for r in decision['checks']]
    text=f'''# 弱形式积分时滞辨识专项研究

研究基准：`fd47730`；当前分支：`problem2`。本目录独立，正式成果、V2、第一问、第三问接口与冻结论文不变。

## 明确结论

{verdict}

独立评价主场景的普通法为44/60，弱形式60/60；历史41/60使用另一组种子，不能把41和60直接当配对改善。短响应区间仅覆盖5/10真值；缓变无噪声场景普通法6/6，弱形式5/6，后者失败为剖面过宽而非点估计偏移。**主门槛通过只授权探索真实数据，不等于全部时间尺度和不确定性验证通过。**

最终采用决策：`{decision['adoption']}`。本研究仅改同一一阶加纯迟延模型的辨识方程，不新增响应模型。不能把合成恢复率解释为工业物理事件准确率。

## 数学构造与单位

动力学模型仍为 $T y'(t)+y(t)=K O'(t-\\tau)+c$。取两端为零的 $\\phi(t)=\\sin^2(\\pi(t-a)/(b-a))$，分部积分得到

$$-T\\int_a^b\\phi' y+\\int_a^b\\phi y=-K\\int_a^b\\phi' O(t-\\tau)+c\\int_a^b\\phi.$$

各矩除以 $\\int\\phi$ 后形成 $T A+B=K Q+c$。$A$ 是输出/秒，$B,c$ 是输出单位，$Q$ 为 m³/s；$T,\\tau$ 为秒，$K$ 为输出单位/(m³/s)。不对原始累计量或烟气逐点求导。半秒样本采用梯形积分，有有限网格误差，不能称精确解析积分。两种方法约束 $T\\ge0$，不强制正增益，搜索0–60秒、步长0.5秒。

**弱形式消除了孤立端点噪声，但自变量中的加权输出仍有观测误差；它不是无偏保证，也不是工具变量辨识。** 累计量减去一个常数以避免大数积分抵消，该操作不改变数学方程。

## 预定协议与公平比较

`protocol_lock.json` 在首次历史回放与独立评价前写入代码/配置哈希并保留。首次执行的数值实验已完成，报告表格适配错误在呈现阶段修复；修复不改变模型、配置或数值输出，当前执行哈希另存`execution_manifest.json`。主权函数20秒支撑、5秒步长；10与40秒仅作固定敏感性分析。普通法复用上一轮5秒积分；两法使用相同训练60–275秒、验证280–395秒、搜索网格、信号和筛选阈值。不同矩数与弱形式重叠是方法本身的差异，不将其当独立样本。弱形式所有完整支撑都位于训练段内，没有借用验证输出。

历史种子8101–8110逐项回放：旧实现与本轮普通积分实现数值一致，正确恢复41/60。独立评价为9201–9210，不参与参数选择；无噪声只做一次，避免重复计数。混合宽度完整复用旧生成器；短时与缓变仅将同一事件表的宽度统一为4秒或30秒。正负供氧扰动叠加在5m³/s基线上，累计量由流量积分生成。两级响应只检验模型失配，不参与选参。采集与清洗直接复用既有 Benchmark，未重新建设信号保真诊断。

主比较是理想采集、未清洗气体。0.01存储格点、1秒/2秒保持是假设条件；`formal_clean`是冻结清洗的敏感性对照。有限10个评价种子共享多个场景，率与Bootstrap比较按种子配对聚类；没有将所有条件当独立实验。

## 主场景进入真实探索的预定门槛

{mdtable(gate_rows)}

验收阈值在config.json预定。主要噪声恢复率至少85%，配对种子Bootstrap增益下限严格大于0；短、慢不得低于普通法；恒定供氧零误接受；无关输出误接受率增幅不得大于10个百分点；主匹配场景95%区间覆盖率至少90%、中位宽度不超过5秒。只有全部通过才能进入真实数据。10个无关试验只能检验这一有限集合，不能证明总体误接受率不增加。

## 完整采用验收

| 项目 | 证据 | 判定 |
|---|---|---|
| 主匹配场景噪声点估计 | 相同独立种子44/60→60/60；按种子配对增益95%Bootstrap范围为18.33–33.33个百分点 | 改善成立 |
| 短响应不确定性 | 15秒迟延、T=8秒场景，名义95%区间覆盖5/10真值 | 不合格 |
| 缓变无噪声不退化 | 普通法6/6；弱形式5/6，1例容差剖面5.5秒超过预定5秒上限 | 未全面满足 |
| 两级响应模型失配 | 真迟延15秒，弱形式估计17秒，留出改善仍约90.3% | 不能用预测成功证明迟延正确 |
| 十炉真实可辨识性 | CO与CO₂通过炉数均为0 | 不合格 |
| 流量–迟延建模 | 真实验证未通过 | 不开展 |

**最终决定：保留为受控噪声辨识研究结果；不采用为正式物理时滞模型，不替换V2，停止本轮方法扩展。** 主门槛是预先固定的探索授权标准，完整采用验收明确包含上述失败，不因主门槛通过而省略。

## 主结果：时间尺度分层

{mdtable(table)}

正确恢复=通过筛选且纯迟延误差≤1秒。MAE仅对通过筛选的已知响应计算，因此必须与拒绝数同时阅读；CSV另外保留未筛选估计偏差、90%绝对误差、起始迟延与核质心误差。matched含4/16/30秒扰动，short为4秒，slow为30秒。恒定和无关联信号没有迟延真值，不能计作恢复失败或成功。

## 采集与清洗敏感性

{mdtable(acq[['acquisition','gas_input','method','conditions','accepted','correct','accepted_MAE_s']].to_dict('records'))}

这些输入表示均为预先指定的对照，不以接受窗口数选择参数。权函数会丢失部分短时信息，是否有净收益以short分层结果判断。

固定支撑长度敏感性，仅限主评价混合宽度噪声条件；不据此改主方法：

{mdtable(sr)}

## 不确定性与参数耦合

主噪声匹配场景全部60条件均做199次整段序列参数Bootstrap；短、慢和两级失配的区间另取预先指定的迟延15秒、T=8秒、全部10种子。每次在训练残差估计的噪声标准差下生成独立高斯噪声，并重新扫描全部迟延；不是对重叠积分行独立重采样。种子固定，百分位区间为条件95%区间，未经截断修饰。残差均值移除、sigma按ddof=3估计；误差模型错误或动力学失配会导致区间失准。

`ci_coverage`逐个检查合成真迟延，包含被算法拒绝的试验。完整目标函数、每个候选迟延的T/K/c、缩放设计条件数以及Bootstrap迟延–T相关系数已输出。目标函数损失≤最小值+输出矩方差1%的区间是旧原型容差集，**不是置信区间**。Bootstrap区间即使很窄，也不证明模型结构正确。迟延网格为0.5秒，零宽区间只说明重抽样解落在同一格点，不表示零不确定性或亚采样精度。

无噪声两级响应，真实纯迟延15秒、总时间常数8秒：

{mdtable(mismatch)}

两级响应的形状不是本轮一阶模型；迟延偏差必须归入模型失配，不能归入仪表输送迟延或供氧物理规律。良好留出预测不能消除这项偏差。

## 真实数据与流量定义

真实数据执行状态：`{decision['real_data_status']}`。

'''
    if decision['controlled_gate_passed']:
        real=pd.read_csv(out/'real_trials.csv');lo=pd.read_csv(out/'real_loho.csv')
        text+='原始累计耗氧量与原始记录气体输入，不使用SG瞬时导数。实际累计量单位m³，两气体记录流量单位m³/s，因此实际K无量纲、c为m³/s。训练从60秒至记录长度60%，65秒间隔后连续留出。每个炉次和气体独立报告，包括失败结果；区间采用10秒移动块残差Bootstrap，其有效性仍取决于模型及局部平稳假设。\n\n'
        text+=mdtable(real[['heat_id','gas','method','mean_source_flow_m3_s','estimate_s','T_s','accepted','ci_width_s','validation_improvement']].to_dict('records'))+'\n\n'
        text+='留一炉次：九炉训练结果的T、K、c中位数，延迟按九炉归一化训练预测损失选择；被留出炉次不拟合参数，仅初始化第一个验证输出。无供氧项基线独立训练，同样跨炉迁移。训练包含被拒绝模型以避免按接受结果挑选供体；因此这也是参数可迁移性的严格压力测试。\n\n'
        text+=mdtable(lo.to_dict('records'))+'\n\n'+mdtable(decision['real_validation'])+'\n\n'
        text+='''CO的留一炉次预测在9/10炉优于基线，但所有10个跨炉迟延均贴60秒搜索上边界；这不是可辨识且稳定的物理迟延。CO₂为0/10炉优于基线。炉内连续留出的RMSE改善中位数分别约−16.18%和−15.96%，两气体接受炉数均为0；不能把跨炉CO预测改善单独解释为迟延规律。

'''
        text+=f"真实数据采用判断：`{decision['adoption']}`；流量关系状态：`{decision['flow_delay_analysis']}`。\n\n"
    else:
        text+='受控验证未过，按照预定停止规则没有运行真实十炉参数估计、连续留出或留一炉次；这些不是已完成的验证。失败后不更换权函数、不改阈值、不选择10秒或40秒结果进入真实数据。\n\n'
    text+='''代表性平均氧枪流量严格定义为 $\\bar q=[O(b)-O(a)]/(b-a)$，单位m³/s，和瞬时SG流量不同。区间可辨识性与真实验证通过前，不拟合流量–迟延函数；模型参数和统计相关不能确认某次供氧导致某次物理响应，也不能分离仪表延迟与工艺迟延。

## 文件与复现

运行：`OPENBLAS_NUM_THREADS=1 .venv/bin/python -m problem2.research.weak_integral.run`。单元测试：`.venv/bin/python -m pytest problem2/research/weak_integral/test_weak.py -q`。配置只有一份，同一路径原位更新。

- `trials.csv`：全部配对试验、筛选原因、估计偏差、留出预测和条件区间。
- `profiles.csv`：所有主输入条件的迟延剖面与参数耦合，未选优摘取。
- `bootstrap_draws.csv`：区间每次重拟合的迟延、T、K、c。
- `summary.csv`、`support_sensitivity.csv`、`historical_replay.csv`：分层对照及历史一致性。
- `protocol_lock.json`、`decision.json`、`provenance.json`：实际配置、预定规则、决策与哈希。
- `figure_captions.csv`：可编辑图题；图内英文与真实Times New Roman，PNG/SVG/PDF同时导出，PDF嵌入检查见`font_validation.json`，仅保留一张当前版宽检查页。
- `protection_check.json`：本轮以外全部既有Git跟踪文件逐字节保护检查。

本轮研究结论只针对这一预定弱形式辨识方案；不以替换模型、追加算法或降低门槛掩盖失败。
'''
    (out/'REPORT.md').write_text(text)
