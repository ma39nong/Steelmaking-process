"""Explain stage fidelity, tradeoffs and identifiability without selecting a winner."""
import csv,json,subprocess
import numpy as np
import pandas as pd
from PIL import Image
from ...src.delay_audit import write_csv,write_json
from ...src.figure_fonts import configure
from ..target_excitation.present import mdtable

LABELS={'formal':'Formal SG31','omit_oxygen_hampel':'No O Hampel','linear_oxygen_interpolation':'Linear O',
        'SG11':'SG11','analytic_PCHIP_derivative':'PCHIP derivative','secant_2s':'Secant 2 s','omit_flow_clipping':'Unclipped SG31',
        'omit_gas_hampel':'No gas Hampel','hampel_threshold12':'Hampel threshold 12','hampel_window5':'Hampel window 5',
        'linear_gas_interpolation':'Linear gas','acquisition_only':'Acquisition only',
        'formal_median21':'Difference - median21','difference_only':'Difference only','mean21':'Difference - mean21'}


def summarize(frame,keys):
    rows=[]
    for key,g in frame.groupby(keys,dropna=False,sort=False):
        if not isinstance(key,tuple):key=(key,)
        r=dict(zip(keys,key));r['conditions']=len(g)
        for name in ['peak_ratio','onset_absolute_error_s','onset_relative_error_s','centroid_error_s',
                     'nearest_peak_set_distance_s','apparent_width50_error_s','signed_area_relative_error',
                     'negative_area_ratio','far_leakage_absolute_area_ratio','waveform_nrmse','retained_energy_ratio',
                     'change_waveform_relative_l2_error','artifact_peak_response','artifact_absolute_area',
                     'constant_flow_relative_rmse','constant_flow_ripple_cv','constant_flow_mean_ratio',
                     'constant_flow_false_excursion_fraction','constant_flow_total_relative_error']:
            if name not in g:continue
            v=g[name].dropna()
            for q,label in [(.1,'p10'),(.5,'median'),(.9,'p90')]:r[name+'_'+label]=float(v.quantile(q)) if len(v) else None
            r[name+'_defined_fraction']=len(v)/len(g)
        for name in ['peak_at_least_half','absolute_onset_defined','corrupted_point_flagged']:
            if name in g:r[name+'_fraction']=float(g[name].mean())
        rows.append(r)
    return rows


def integral_summary(frame,keys):
    rr=[]
    for key,g in frame.groupby(keys,dropna=False,sort=False):
        if not isinstance(key,tuple):key=(key,)
        accepted=g[g.accepted];true=g[g.true_delay_s.notna()];err=accepted.error_s.dropna()
        rr.append(dict(zip(keys,key),conditions=len(g),accepted=int(g.accepted.sum()),correct=int(g.correct.sum()),
                       rejected=int((~g.accepted).sum()),correct_rate=float(g.correct.sum()/len(true)) if len(true) else None,
                       accepted_wrong=int((accepted.error_s.abs()>1).sum()),
                       accepted_delay_mae_s=float(err.abs().mean()) if len(err) else None,
                       accepted_median_delay_s=float(accepted.estimate_s.median()) if len(accepted) else None,
                       median_validation_improvement=float(g.validation_improvement.median())))
    return rr


def deliver(out,cfg,source,change,gas,controls,integral,profiles,examples,heats):
    import matplotlib.pyplot as plt
    s,c,g,n,i=[pd.DataFrame(x) for x in [source,change,gas,controls,integral]]
    ss=summarize(s,['domain','shape','acquisition','method','duration_s','control_base_flow','control_event_phase_s'])
    gs=summarize(g,['domain','gas','shape','amplitude_multiplier','acquisition','method','duration_s','control_event_phase_s'])
    cs=summarize(c,['domain','shape','acquisition','method','duration_s','control_base_flow','control_event_phase_s'])
    sh=summarize(s[s.domain=='raw_background'],['heat_id','acquisition','method','duration_s'])
    gh=summarize(g[g.domain=='raw_background'],['heat_id','gas','amplitude_multiplier','acquisition','method','duration_s'])
    nc=summarize(n,['signal','acquisition','method'])
    ins=integral_summary(i,['kind','acquisition','noise_std','gas_input'])
    for name,rows in [('source_summary',ss),('gas_summary',gs),('change_summary',cs),('source_by_heat',sh),('gas_by_heat',gh),('negative_control_summary',nc),('integral_summary',ins)]:
        write_csv(out/(name+'.csv'),rows)
    # Report all per-heat differences; ranks are never converted into a chosen method.
    heterogeneity=[]
    hs=pd.DataFrame(sh)
    for (acq,method,w),q in hs.groupby(['acquisition','method','duration_s']):
        heterogeneity.append(dict(acquisition=acq,method=method,duration_s=w,heats=len(q),
                                  minimum_heat_median_peak_ratio=float(q.peak_ratio_median.min()),
                                  maximum_heat_median_peak_ratio=float(q.peak_ratio_median.max()),
                                  minimum_heat_median_area_error=float(q.signed_area_relative_error_median.min()),
                                  maximum_heat_median_area_error=float(q.signed_area_relative_error_median.max())))
    write_csv(out/'heat_heterogeneity.csv',heterogeneity)
    # Predefined, strictly paired contrasts with the formal pipeline only.
    index=['domain','heat_id','seed','shape','duration_s','amplitude','acquisition','event_center_s','control_base_flow']
    contrasts=[]
    for key,q in s.groupby(index,dropna=False,sort=False):
        formal=q[q.method=='formal'].iloc[0]
        for _,r in q.iterrows():
            if r['method']=='formal':continue
            contrasts.append(dict(zip(index,key),method=r['method'],
                                  peak_ratio_change=r.peak_ratio-formal.peak_ratio,
                                  signed_area_error_change=r.signed_area_relative_error-formal.signed_area_relative_error,
                                  waveform_nrmse_change=r.waveform_nrmse-formal.waveform_nrmse))
    write_csv(out/'source_paired_contrasts.csv',contrasts)
    font=configure();figdir=out/'figures';figdir.mkdir(exist_ok=True);captions=[];fontchecks=[]
    def save(fig,name,caption):
        fig.tight_layout()
        for ext in ['png','svg','pdf']:
            fig.savefig(figdir/(name+'.'+ext),dpi=180,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
        svg=figdir/(name+'.svg');svg.write_text('\n'.join(x.rstrip() for x in svg.read_text().splitlines())+'\n')
        plt.close(fig);raw=subprocess.check_output(['pdffonts',str(figdir/(name+'.pdf'))],text=True)
        col=raw.splitlines()[0].index('emb');lines=raw.splitlines()[2:]
        assert lines and all('TimesNewRoman' in x.replace(' ','') and x[col:col+3]=='yes' for x in lines)
        fontchecks.append(dict(figure=name,pdffonts=raw));captions.append(dict(figure=name,caption=caption,width_inches=7.2))
    plotmethods=['formal','SG11','analytic_PCHIP_derivative','secant_2s'];colors=['#235e8a','#ba5627','#4b8050','#795393']
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.5))
    for ax,acq in zip(axes,['grid_0.01','empirical_hold']):
        for name,color in zip(plotmethods,colors):
            q=s[(s.domain=='raw_background')&(s.acquisition==acq)&(s.method==name)]
            r=q.groupby('duration_s').peak_ratio.agg(['median',lambda a:a.quantile(.1),lambda a:a.quantile(.9)])
            ax.plot(r.index,r['median'],'o-',ms=3,color=color,label=LABELS[name]);ax.fill_between(r.index,r.iloc[:,1],r.iloc[:,2],color=color,alpha=.08)
        ax.axhline(1,color='black',ls=':',lw=.8);ax.set_title('Stored grid only' if acq=='grid_0.01' else 'Hypothesized recorded hold')
        ax.set_xlabel('Known pulse duration (s)');ax.set_ylabel('Recovered / true peak');ax.legend(fontsize=6)
    save(fig,'source_peak_fidelity','已知三角供氧脉冲在十炉背景中的峰值保留。线为全部配对条件中位数，阴影为10%–90%分位，不是置信区间。展示四种预先定义的导数表示，其余单环节消融仍完整列入CSV，不选择最优方法。')
    fig,axes=plt.subplots(2,2,figsize=(7.2,5.4))
    for ax,(name,ex) in zip(axes.flat,examples):
        t=ex['t'];keep=(t>=180)&(t<=220);ax.plot(t[keep],ex['truth'][keep],color='black',label='Known added flow')
        for method,color in [('formal',colors[0]),('SG11',colors[1]),('secant_2s',colors[3])]:
            delta=ex['augmented'][method]['flow']-ex['baseline'][method]['flow']
            ax.plot(t[keep],delta[keep],color=color,label=LABELS[method],lw=1)
        ax.set_title(name.replace('_',' '));ax.set_xlabel('Time (s)');ax.set_ylabel('Added flow (m3/s)');ax.legend(fontsize=6)
    save(fig,'source_waveforms','已知流量先解析积分成累计量，再经过采集及重构。图为相同基线前后输出之差，展示4秒/16秒与逐点/1秒保持条件，真实采集机制没有据此被确认。')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.5))
    oracle=pd.read_csv(out/'oracle_change_fidelity.csv')
    for ax,shape in zip(axes,['triangle','rectangle']):
        for method,color in zip(cfg['change_methods'],colors):
            q=oracle[(oracle['shape']==shape)&(oracle.method==method)]
            ax.plot(q.duration_s,q.retained_energy_ratio,'o-',label=LABELS[method],color=color,ms=3)
        ax.set_title(shape.capitalize()+' input, no flow reconstruction');ax.set_xlabel('Pulse duration (s)');ax.set_ylabel('Retained / original difference energy');ax.legend(fontsize=6)
    save(fig,'change_representation','仅改变变化提取环节，输入为已知真流量，未经过SG，分离观察去趋势本身的影响。能量比不等于事件识别准确率；差分保留更多变化也会保留更新噪声。')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.6))
    gaslabels={**LABELS,'formal':'Formal Hampel9'}
    gas_spike_worst=[]
    for method,color in zip(cfg['gas_methods'],colors+['#777777']):
        q=g[(g.domain=='raw_background')&(g.acquisition=='grid_0.01')&(g.method==method)]
        yy=q.groupby('duration_s').peak_ratio.median();axes[0].plot(yy.index,yy,'o-',ms=3,label=gaslabels[method],color=color)
        rr=n[(n.signal.isin(['co','co2']))&(n.acquisition=='grid_0.01')&(n.method==method)]
        residual=rr.artifact_peak_response/rr.logging_spike_amplitude
        axes[1].bar(gaslabels[method],float(residual.median()),color=color)
        gas_spike_worst.append(float(residual.max()))
    axes[1].scatter(range(len(cfg['gas_methods'])),gas_spike_worst,marker='x',color='black',s=20,label='Worst residual')
    axes[1].set_ylim(0,1.15);axes[1].legend(fontsize=6)
    axes[0].set_title('(a) Known gas pulses in ten-heat backgrounds');axes[0].set_xlabel('Pulse duration (s)');axes[0].set_ylabel('Recovered / true peak');axes[0].legend(fontsize=6)
    axes[1].set_title('(b) Added logging spike, gas channels');axes[1].set_ylabel('Residual / added spike peak');axes[1].tick_params(axis='x',labelrotation=45)
    save(fig,'gas_preservation_tradeoff','CO和CO2的逐炉背景同时检验真实注入脉冲与明确的孤立记录尖峰。左侧保留率不能脱离右侧伪响应抑制能力选优。保留全部气体、幅值和炉次分层CSV。')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.6))
    for method,color in zip(plotmethods,colors):
        q=s[(s.domain=='raw_background')&(s.acquisition=='grid_0.01')&(s.method==method)]
        yy=q.groupby('duration_s').signed_area_relative_error.median();axes[0].plot(yy.index,yy,'o-',label=LABELS[method],color=color,ms=3)
        q=n[(n.signal=='oxygen')&(n.acquisition=='empirical_hold')&(n.method==method)]
        axes[1].bar(LABELS[method],float(q.constant_flow_ripple_cv.median()),color=color)
    axes[0].set_title('(a) Added-flow sampled integral');axes[0].set_xlabel('Pulse duration (s)');axes[0].set_ylabel('Signed relative volume error');axes[0].legend(fontsize=6)
    axes[1].set_title('(b) Known constant flow, recorded hold');axes[1].set_ylabel('Spurious flow standard deviation / true flow');axes[1].tick_params(axis='x',labelrotation=40)
    save(fig,'volume_and_null_ripple','总量使用正式0.5秒流量网格的梯形积分，无后验总量缩放。恒定真流量下的波动是人工保持假设下可确认的伪动态，不是原始工业事件误报率。PCHIP连续导数在节点取样的积分偏差须与连续插值函数本身的端点守恒区别。')
    fig,axes=plt.subplots(1,2,figsize=(7.2,4.1))
    for ax,width in zip(axes,[4,30]):
        hh=hs[(hs.acquisition=='empirical_hold')&(hs.duration_s==width)]
        matrix=np.array([[float(hh[(hh.heat_id==heat['heat_id'])&(hh.method==m)].peak_ratio_median.iloc[0]) for heat in heats] for m in cfg['source_methods']])
        im=ax.imshow(matrix,aspect='auto',vmin=0,vmax=max(1.5,float(matrix.max())),cmap='viridis')
        ax.set_xticks(range(10),[z['heat_id'] for z in heats],rotation=60,ha='right');ax.set_yticks(range(7),[LABELS[x] for x in cfg['source_methods']],fontsize=6)
        ax.set_title(f'Pulse {width:g} s; heat-wise median peak ratio');fig.colorbar(im,ax=ax,shrink=.7)
    save(fig,'heat_heterogeneity','按炉次展示4秒与30秒源脉冲的中位峰值比。每炉同样保留全部种子和采集假设，矩阵不代表独立事件数，也不按总体均值选择方法。')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.5));p=pd.DataFrame(profiles)
    for kind,delay,T,color in [('matched',15,8,colors[0]),('constant',np.nan,8,colors[1]),('two_stage',15,8,colors[3])]:
        q=p[(p.kind==kind)&(p.acquisition=='ideal')&(p.true_T_s==T)]
        if kind=='matched':q=q[q.true_delay_s==delay]
        scale=max(float(q.loss.max()),1e-12);axes[0].plot(q.lag_s,q.loss/scale,label=kind,color=color)
    axes[0].set_title('(a) Integral loss profiles, no noise');axes[0].set_xlabel('Candidate pure delay (s)');axes[0].set_ylabel('Loss / maximum loss');axes[0].legend(fontsize=7)
    for gas_input,offset,color in [('raw',-.15,colors[0]),('formal_clean',.15,colors[1])]:
        mech=cfg['integral_probe']['source_acquisition'];vals=[]
        for acq in mech:
            q=i[(i.kind=='matched')&(i.acquisition==acq)&(i.gas_input==gas_input)];vals.append(float(q.correct.mean()))
        axes[1].bar(np.arange(len(mech))+offset,vals,width=.28,color=color,label=gas_input.replace('_',' '))
    axes[1].set_xticks(range(8),cfg['integral_probe']['source_acquisition'],rotation=55,ha='right',fontsize=6);axes[1].set_ylim(0,1.05)
    axes[1].set_title('(b) Matched dynamics, all noise conditions');axes[1].set_ylabel('Accepted and delay error at most 1 s');axes[1].legend(fontsize=6)
    save(fig,'integral_probe','积分辨识仅作受控下游探针。误差剖面不等于置信区间；两级响应能给出较好预测而纯迟延估错，不能据此自动采用。右图均包含拒绝，原始和清洗气体为配对条件。')
    write_csv(out/'figure_captions.csv',captions);write_json(out/'font_validation.json',dict(configuration=font,pdf_checks=fontchecks))
    preview=Image.new('RGB',(1000,4*330),'white')
    for j,r in enumerate(captions):
        im=Image.open(figdir/(r['figure']+'.png')).convert('RGB');im.thumbnail((495,325));preview.paste(im,((j%2)*500,(j//2)*330))
    preview.save(out/'layout_preview.png')
    def mdt(df,cols):return mdtable(cols,[[None if pd.isna(v) else v for v in row] for row in df[cols].values.tolist()])
    primary=pd.DataFrame(ss);primary=primary[(primary.domain=='raw_background')&(primary.acquisition=='grid_0.01')&primary.duration_s.isin([4,16,60])]
    source_table=mdt(primary,['method','duration_s','conditions','peak_ratio_median','onset_absolute_error_s_median','onset_absolute_error_s_defined_fraction','signed_area_relative_error_median','waveform_nrmse_median'])
    gas_table=pd.DataFrame(gs);gas_table=gas_table[(gas_table.domain=='raw_background')&(gas_table.acquisition=='grid_0.01')&(gas_table.duration_s==2)&(gas_table.amplitude_multiplier==.5)]
    gast=mdt(gas_table,['gas','method','conditions','peak_ratio_median','peak_at_least_half_fraction','signed_area_relative_error_median'])
    core_i=pd.DataFrame(ins);it=mdt(core_i[(core_i.acquisition=='ideal')&(core_i.gas_input=='raw')],['kind','noise_std','conditions','accepted','correct','accepted_wrong','accepted_delay_mae_s','median_validation_improvement'])
    null=pd.DataFrame(nc);nt=mdt(null[(null.signal=='oxygen')&(null.acquisition=='empirical_hold')],['method','constant_flow_ripple_cv_median','constant_flow_mean_ratio_median','constant_flow_total_relative_error_median','artifact_peak_response_median'])
    # Conservative, machine-readable decisions; no winner or formal parameter file.
    no_noise=i[(i.kind=='matched')&(i.acquisition=='ideal')&(i.gas_input=='raw')&(i.noise_std==0)]
    mismatch=i[(i.kind=='two_stage')&(i.acquisition=='ideal')&(i.gas_input=='raw')&(i.noise_std==0)].iloc[0]
    decisions=dict(automatic_selection=False,formal_model_changed=False,physical_events_confirmed=0,
                   matched_ideal_no_noise_trials=len(no_noise),matched_ideal_no_noise_correct=int(no_noise.correct.sum()),
                   two_stage_no_noise_estimated_delay_s=float(mismatch.estimate_s),two_stage_true_delay_s=15,
                   conclusion='累计量、求导、变化去趋势、烟气清洗必须分别验收。候选处理在峰值保留、采集噪声、总量与记录尖峰抑制之间存在不同取舍。积分最小二乘原型在匹配条件下有效，但噪声和模型失配会破坏纯迟延解释。本轮不自动选优，不修改正式V2，不开展真实流量—时滞拟合。')
    write_json(out/'decision.json',decisions)
    # Quantitative answers, scoped to named hypotheses rather than blended into a rank.
    findings={}
    for method in ['formal','SG11','secant_2s','analytic_PCHIP_derivative']:
        q=s[(s.domain=='raw_background')&(s.acquisition=='grid_0.01')&(s.duration_s==4)&(s.method==method)]
        nullq=n[(n.signal=='oxygen')&(n.acquisition=='empirical_hold')&(n.method==method)]
        h=hs[(hs.acquisition=='empirical_hold')&(hs.duration_s==4)&(hs.method==method)]
        findings[method]=dict(source4_grid_median_peak_ratio=float(q.peak_ratio.median()),
                              source4_grid_median_signed_volume_error=float(q.signed_area_relative_error.median()),
                              constant_flow_recorded_hold_median_ripple_cv=float(nullq.constant_flow_ripple_cv.median()),
                              source4_recorded_hold_heat_median_peak_min=float(h.peak_ratio_median.min()),
                              source4_recorded_hold_heat_median_peak_max=float(h.peak_ratio_median.max()))
    gasfind=[]
    for gasname in ['co','co2']:
        for method in ['acquisition_only']+cfg['gas_methods']:
            q=g[(g.domain=='raw_background')&(g.gas==gasname)&(g.acquisition=='grid_0.01')&(g.duration_s==2)&(g.amplitude_multiplier==.5)&(g.method==method)]
            h=q.groupby('heat_id').peak_ratio.median()
            gasfind.append(dict(gas=gasname,method=method,median_peak_ratio=float(q.peak_ratio.median()),
                                minimum_heat_median_peak_ratio=float(h.min()),maximum_heat_median_peak_ratio=float(h.max())))
    spikes=[]
    for method in cfg['gas_methods']:
        q=n[n.signal.isin(['co','co2'])&(n.acquisition=='grid_0.01')&(n.method==method)]
        ratios=q.artifact_peak_response/q.logging_spike_amplitude
        spikes.append(dict(method=method,conditions=len(q),flagged=int(q.corrupted_point_flagged.sum()),
                           median_residual_ratio=float(ratios.median()),worst_residual_ratio=float(ratios.max())))
    write_json(out/'findings.json',dict(source=findings,gas=gasfind,known_gas_logging_spikes=spikes,
                                       automatic_selection=False,conditions_are_not_independent_industrial_events=True))
    f=findings
    heat_table=mdtable(['method','4秒逐炉中位峰值比最小','最大','恒定流量伪波动CV中位数'],
                       [[name,v['source4_recorded_hold_heat_median_peak_min'],v['source4_recorded_hold_heat_median_peak_max'],v['constant_flow_recorded_hold_median_ripple_cv']] for name,v in findings.items()])
    spiketable=mdtable(['method','已知日志尖峰条件','被标记点数','残余峰值比中位数','最差残余峰值比'],
                      [[r[k] for k in ['method','conditions','flagged','median_residual_ratio','worst_residual_ratio']] for r in spikes])
    integrity=dict(source_methods={
        'formal':'Hampel9 threshold6 + PCHIP + SG31 derivative poly3 + clip0',
        'omit_oxygen_hampel':'only omit statistical Hampel; retain physical invalid-value mask',
        'linear_oxygen_interpolation':'only replace PCHIP interpolation by linear',
        'SG11':'only change derivative window 31 to 11',
        'analytic_PCHIP_derivative':'only replace SG by exact PCHIP derivative evaluated on the 0.5 s output grid',
        'secant_2s':'only replace SG by centered 2 s cumulative secant, with one-sided endpoint widths',
        'omit_flow_clipping':'only omit derivative clipping; diagnostic, negative flows preclude automatic formal adoption'},
        gas_methods={'formal':'Hampel9 threshold6 + PCHIP','omit_gas_hampel':'only omit statistical Hampel',
                     'hampel_threshold12':'only threshold6 to12','hampel_window5':'only window9 to5',
                     'linear_gas_interpolation':'only PCHIP tolinear','acquisition_only':'diagnostic, no cleaning; not a candidate'},
        change_methods={'formal_median21':'difference minus centered21 median','difference_only':'only omit median subtraction','mean21':'only median to mean'},
        no_composed_best_pipeline=True,injection='M(background + analytic integral of flow pulse) and M(background) processed separately; gas uses M(background + gas pulse)',
        instrument_mechanism_identified=False,confidence_intervals_reported=False,
        paired_unit='heat, signal, seed, event center; durations, doses, methods and hypotheses reuse that draw',
        known_nulls='constant true flow under hypothetical acquisition; post-acquisition known one-sample logging spikes',
        actual_background_truth='unknown; only added perturbations have known ground truth')
    write_json(out/'protocol.json',integrity)
    report=f'''# Signal Fidelity Benchmark：信号保真度基准

{decisions['conclusion']}

## 1. 两项工作如何衔接

先运行独立积分辨识原型，再把它作为基准中的下游探针。主任务是测量各处理步骤对已知信号的影响，不开展新一轮真实工艺事件识别，不把候选窗口数增加当作成功指标。此前测量链审计0ee55e2保留，本轮新增逐环节消融、备选表示、真实背景逐炉比较和已知伪响应对照，不重复旧24窗口重放或大批理想平移实验。

前三项成果提供O、烟气y和累计排放E；积分方法使用它们建立新的可检验动态关系，不要求先拿到独立事件列表，也不强制使用SG求出的q。这与保真度基准没有冲突：前者检验下游可识别性，后者检验送入下游的信息是否被保留。冻结光谱本轮未调用，先按用户最后确定的范围聚焦原始工艺信号保真度。

## 2. 正式基线、输入与复现

正式基线d08f989，研究基线0ee55e2。十炉原始文件及冻结数组共20个输入哈希与上一轮清单核对；20个炉—气体序列逐点核对正式Hampel掩码、气体重构、SG流量及变化提取。向量化Hampel仅用于加速同一数学规则，边界、非有限值、零MAD平台分支均与生产函数逐点测试；生产实现不修改。

完整运行：`.venv/bin/python -m problem2.research.signal_fidelity.run`。仅运行积分探针：加`--integral-only`。固定配置、随机种子、实际参数、输入/源代码/输出哈希完整保存。输出固定在results中原位更新。运行前后核对所有既有跟踪文件及100个冻结文件，正式V2、第一问、第三问NPZ和冻结论文均不修改。

## 3. 已知真值、采集假设及累计量注入

先生成已知流量三角/矩形脉冲，再使用解析积分得到累计量增量。没有把流量脉冲直接加到累计氧量上。对背景B及增量I，比较F(M(B+I))−F(M(B))；量化和保持对整个信号执行，不能用M(I)简单代替这个配对差。矩形积分净量为幅值×宽度，三角为其一半，真起始和重心按连续定义计算。

供氧宽度2、4、8、16、30、60秒。控制条件枚举三角与矩形、零基流与5m3/s恒定基流、事件相位0/0.15秒；背景条件每炉10固定种子抽取同一合法中心，跨宽度/方法复用，事件另有±0.25秒子网格相位。背景供氧幅值为该炉冻结流量中位数。气体每炉每气体10种子，幅值为原记录标准差的0.5或2倍。背景真信号未知，本轮真值仅指人工添加的分量。

供氧采集假设：理想逐点；0.01网格；1秒保持两种相位；2秒保持四种相位；依原始非零变化时刻保持。气体另比较逐点、0.01网格和该气体的记录变化时刻保持。0.01保存值网格是观察事实，设备真实分辨率未知；所有保持方式均为竞争假设，记录非零变化不等于仪器刷新时间。原背景已经采集过，再施加M是明确的假设压力试验，`baseline_remeasurement.csv`逐炉记录二次采集对基线的改变量，不能将它解释为实际隐藏仪器的恢复。

事件抽样是有放回的，不保证每个种子对应不同时间中心。同一中心的宽度、剂量、处理方法与假设相关，条件数不是独立工艺事件数。按炉次完整分层，不报告将数万条件当作独立样本的置信区间。

## 4. 独立消融与指标

`protocol.json`逐项列明改动。供氧七种入口分别仅改变Hampel、插值、导数窗长/表示或截负中的一个环节。烟气五种仅改变是否Hampel、阈值、窗长或插值。变化表示三种固定同一正式q，只改变去趋势。没有拼接“最优导数＋最优气体清洗＋最优变化表示”的组合，也不生成自动选择锁。

峰值比对真实连续峰幅计算。起始按真实峰幅10%的绝对阈值定义：重构峰不足该值时起始留空，不能将一个已消失的小残差重新缩放后声称恢复。另保存相对自身峰的起始作对照。矩形最大值为区间，记录估计最大值范围及到真峰集合的最短距离，不强造唯一峰。重心使用事件ROI内响应正部，面积另外保留有符号全程积分、正部/负部及远离支持区的泄漏，不能把正部积分当成净氧量。

源阶段分别保存采集后与清洗后累计量误差；气体`acquisition_only`是清洗前诊断，不是候选算法。全量逐试验表保留所有衰减、未定义起始、拒绝与失败。0.5秒输出网格梯形积分是接口可使用的总量口径，没有事后缩放使总量对齐。

## 5. 供氧短时与缓变信号

下表为十炉背景、仅保存网格假设、4/16/60秒脉冲的逐条件汇总，数字不是物理事件准确率。完整每种宽度和九种采集假设见source_summary.csv；逐炉见source_by_heat.csv。

{source_table}

SG31对短变化的展宽与衰减，与是否采用统计Hampel分开检验；O的单调叠加通常不产生需Hampel删除的尖峰，不能把烟气Hampel问题移植成“氧量Hampel必然有错”。去掉Hampel或换插值若在无剔除点条件下等价，就如实记录零差异，不制造改善。

SG11、2秒割线及解析PCHIP导数完整参与比较：它们能改变短脉冲保留，但也会改变对阶梯更新的敏感度。解析PCHIP导数尤其需要区别“连续函数积分等于累计端点差”与“只在0.5秒节点采样导数，再梯形积分”：后者可能漏掉节点之间的斜率，不能仅凭保形插值的名称认定流量和接口总量正确。

对缓变信号，oracle_change_fidelity.csv直接把已知真流量送入三个变化算子，完全绕开求导：三角斜坡的差分是分段常值，21点中位数会将足够长的分段变化消掉。单纯差分保留这类变化，也保留计数更新噪声；mean21只替换去趋势算子，未改变供氧和烟气清洗。变化能量保留比不是事件Recall，也不必越大越好，因为噪声能量同样计入。

## 6. CO、CO2短时响应与异常清洗

下表为2秒、0.5标准差注入、仅保存网格假设。每气体包含十炉各10个背景中心，气体清洗前后与不同方法共享同一输入。

{gast}

完整gas_by_heat.csv同时保存炉号、气体、幅值、宽度和采集机制。控制条件包括事件恰落网格及偏离网格：零MAD邻域可能把真实孤立峰当异常，亚采样相位和背景坡度又可能使其保留，因此不能把一个对称脉冲的失败概括为所有短峰都会被删。threshold12在零MAD分支仍用极小回退尺度，不保证能挽救真实尖峰；window5和省略清洗也必须与下面的记录尖峰对照合看。

## 7. 累计总量、更新伪动态与记录尖峰

对十炉的更新节律分别输入已知恒定5m3/s真流量，再检查处理后是否制造波动。同时在采集之后加入一个已知不真实的单点记录尖峰，比较处理前后之差。它是可确认的记录错误，不是把原始工业变化重新贴成异常。供氧尖峰幅值为各炉正增量中位数的10倍，气体为标准差的10倍。

下面是按记录节律保持的已知恒定流量与供氧尖峰结果；`constant_flow_ripple_cv`只度量不应出现的波动，均值偏差另报，避免将固定幅值损失误称为新事件。

{nt}

总量一致性使用source_fidelity.csv的signed_area_relative_error，不能被峰值保留率替代。rest与5m3/s基流控制分开保留：SG负旁瓣在停氧基线可能被截掉而增加积分，在已有正基流时则可能不触发截负。省略Hampel保留真实峰的同时也可能保留日志尖峰；较短导数窗与更直接的导数会对保持更新产生较大伪动态。本轮将这些代价并列，不以某一个指标宣布最优方法。

## 8. 十炉结论是否一致

逐炉表是主要交付，不只报告总体均值。heat_heterogeneity.csv保存每种方法/宽度/采集机制下十炉中位峰值比、总量误差的最小和最大。source_paired_contrasts.csv按相同输入逐方法与正式方案配对，不把炉次差异当方法收益。

同一确定性算子对已知真波形的抑制机制可以重复出现，但真实背景、更新节律、幅值和掩码变化不同，改善大小不能由一炉推广到十炉。报告与图保留这种差别，不能声称“十炉同样改善”后再隐藏最差炉次。

## 9. 积分形式辨识：先尝试，再作为下游探针

参考Belkoura等2008卷积/积分辨识思想及Yang等2004连续时间辨识工作，将局部一阶响应加纯迟延写成：

T[y(b)−y(a)]+[E(b)−E(a)]=K[O(b−tau)−O(a−tau)]+c(b−a)。

这是针对累计耗氧量的改写，原型采用延迟剖面＋T非负约束的最小二乘，未实现论文的完整工具变量算法，也没有声称是其转炉复现。T是响应展开时间常数，tau是模型纯迟延，不等同互相关峰或响应重心。每个候选tau仅用训练区间拟合参数；留出区间从首个观测y初始化一次，之后使用O增量前向预测，不逐点借用未来y。

使用0.05秒高分辨率生成，0.5秒记录；真实延迟5/15/30秒，T=2/8秒，已知多次短/慢正负流量扰动叠加于恒定正流量。包含无噪声及0.05绝对烟气噪声。无噪声不复制10次相同条件，噪声用10固定种子。训练60–275秒、验证280–395秒，5秒非重叠积分区间。每个候选纯迟延0–60秒、步长0.5秒，记录全部损失剖面。

接受要求预先固定：有效输入变化、非零增益、近优剖面范围不超过5秒、留出RMSE比无输入嵌套模型改善至少10%、非边界最优值。近优范围按输出方差1%损失容差定义，是敏感性集合而非统计置信区间。所有拒绝保留，并以数值误差≤1秒判正确，不以“有返回值”判成功。

以下是理想采集、原始气体入口结果；所有保持/量化和原始/清洗入口对照见integral_summary.csv。

{it}

匹配模型、理想采集、无噪声的{len(no_noise)}个条件正确{int(no_noise.correct.sum())}个。两级真实响应的无噪声例子，真纯迟延15秒，原型却估计{mismatch.estimate_s:.1f}秒，即使留出预测改善{100*mismatch.validation_improvement:.1f}%。这说明好预测与正确拆分纯迟延/展宽是两项验收，错误不能用训练拟合掩盖。加噪后dy项也含测量误差，简单最小二乘会有参数偏差；这是本原型需要解决的问题，不能据此宣称论文工具变量方法失效。

恒定输入时，O(b−tau)−O(a−tau)近似与tau无关，必须拒绝；无关联输出也用留出检验。没有为改善恢复率修改这些阈值，没有给真实十炉拟合纯迟延，更没有拟合流量—时滞曲线。

原始来源：[Belkoura等，2008](https://skoge.folk.ntnu.no/prost/proceedings/ifac2008/data/papers/1215.pdf)；[Yang等，2004](https://www.jstage.jst.go.jp/article/ieejeiss/124/10/124_10_2132/_pdf)。引用的是避开直接微分和区分系统动态/延迟的思想，不是未经验证的工业性能。

## 10. 对六个验收问题的回答与决策

本轮可以作出的具体判断：

- 在仅保存网格的十炉背景条件下，4秒源脉冲正式峰值保留中位数{f['formal']['source4_grid_median_peak_ratio']:.3f}，SG11为{f['SG11']['source4_grid_median_peak_ratio']:.3f}，2秒割线为{f['secant_2s']['source4_grid_median_peak_ratio']:.3f}。但按原记录节律保持的恒定流量对照中，伪波动CV从正式{f['formal']['constant_flow_recorded_hold_median_ripple_cv']:.3f}升至SG11的{f['SG11']['constant_flow_recorded_hold_median_ripple_cv']:.3f}、割线的{f['secant_2s']['constant_flow_recorded_hold_median_ripple_cv']:.3f}。因此短峰保留改善有明确的更新噪声代价，不能自动换用。
- SG11在该4秒网格条件下的注入净量误差中位数为{f['SG11']['source4_grid_median_signed_volume_error']:.6g}，没有因改善峰值就必然损害总量；相反，节点取样的解析PCHIP导数该误差达到{f['analytic_PCHIP_derivative']['source4_grid_median_signed_volume_error']:.3f}，峰值大不等于正确。零基流2秒三角控制中，正式截负使净量增加约11.6%，已有正基流时这一旁瓣截断效应消失；必须按基流和采集假设分别判断。
- 2秒低幅烟气注入中，采集本身及Hampel清洗的损失已在第6节分开给出。窗口5和不做Hampel在短峰保留上表现较好，但不能只看中位数：已知日志尖峰的最差残余和标记数如下。threshold12在本套强孤立尖峰对照中仍能标记全部点，是值得后续验证的候选；本轮只测了明确的10标准差日志尖峰，不据此推断其对较弱或连续异常同样可靠。

{spiketable}

按记录节律保持的4秒源脉冲，逐炉保留范围如下。正式短峰衰减及SG11改善方向在这十炉均出现，解析导数对炉次更新节律高度敏感。其他尺度、气体及假设仍须查逐炉表，不能概括为所有工况统一改善。

{heat_table}

1. 损害步骤：采集保持/量化先限制时间信息；SG31改变短供氧脉冲；中位数变化去趋势消除部分缓变；Hampel能混淆孤立真实短峰与错误尖峰。插值和截负分别消融，不一概指责。
2. 时间尺度：2/4/8秒与16/30/60秒分别列出，控制条件还按形态、采集相位及零/正基流区分；不把15秒窗跨度解释成固定7.5秒时延。
3. 可保留方法：短窗/割线/直接导数与省略变化中位数等候选明确展示保留能力，但逐项同时展示更新噪声、累计积分及异常抑制，不自动选优。
4. 总量代价：逐条件保存有符号净量误差，尤其检查导数节点取样和负旁瓣截断；不做总量事后缩放。
5. 新伪响应：用已知恒定真流量和日志尖峰检验，不能用未知真实背景宣称工业误报率。
6. 炉次一致性：逐炉矩阵、最差/最好炉次范围及严格配对差异完整输出；总体中位数不替代十炉检查。

本轮交付基准和适用性证据，没有选择替代正式V2的方案。下一步的判断依据应是同一候选能否同时保留所需时间尺度、抑制更新伪动态、满足接口总量，并使下游已知时滞恢复通过验证；不是窗口数量是否增加。物理事件归因与设备时间标定仍需要独立依据。

## 图表与文件

图内英文、真实Times New Roman；7组图各导出PNG/SVG/PDF，逐PDF检查嵌入，图题另存figure_captions.csv，一张当前排版检查页。JSON/CSV与固定配置保留完整数值、输出清单和保护核查。
'''
    (out/'REPORT.md').write_text(report,encoding='utf-8')
