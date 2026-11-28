"""Report separately what is observed, simulated, and unidentified."""
import csv,json,subprocess
from pathlib import Path
import numpy as np
from PIL import Image
from ...src.delay_audit import write_csv,write_json
from ...src.figure_fonts import configure
from ..target_excitation.present import mdtable


def summary(rows,keys):
    groups={}
    for r in rows:groups.setdefault(tuple(r[k] for k in keys),[]).append(r)
    result=[]
    for key,rr in groups.items():
        true=[r for r in rr if r['truth_s'] is not None]
        ee=[abs(r['error_s']) for r in true if r['error_s'] is not None]
        result.append(dict(zip(keys,key),conditions=len(rr),accepted=sum(r['accepted'] for r in rr),
                           correct=sum(r['correct'] for r in true),wrong=sum(r['wrong'] for r in true),
                           correct_rate=sum(r['correct'] for r in true)/len(true) if true else None,
                           reject_rate=sum(not r['accepted'] for r in rr)/len(rr),
                           conditional_mae_s=float(np.mean(ee)) if ee else None,
                           conditional_max_error_s=float(max(ee)) if ee else None))
    return result


def deliver(out,cfg,heats,stats,flows,model,lags,examples,background,timing,ar,loo,boot):
    from .run import ROOT,HERE,sha
    from .core import pulse,errors,held,features
    from ...src.q2_pipeline import hampel_mask,pchip_clean,_standardized_change
    import matplotlib.pyplot as plt
    font=configure();figdir=out/'figures';figdir.mkdir(exist_ok=True);captions=[];fontchecks=[]
    def save(fig,name,caption):
        fig.tight_layout()
        for ext in ['png','svg','pdf']:
            fig.savefig(figdir/(name+'.'+ext),dpi=180,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
        svg=figdir/(name+'.svg')
        svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
        plt.close(fig)
        text=subprocess.check_output(['pdffonts',str(figdir/(name+'.pdf'))],text=True)
        ix=text.splitlines()[0].index('emb');ff=text.splitlines()[2:]
        assert ff and all('TimesNewRoman' in x.replace(' ','') and x[ix:ix+3]=='yes' for x in ff)
        fontchecks.append(dict(figure=name,pdffonts=text));captions.append(dict(figure=name,caption=caption,width_inches=7.2))
    signals=['oxygen_cumulative','co','co2'];labels=['Cumulative oxygen','CO','CO2'];colors=['#236394','#ba592b','#478055']
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.5))
    for j,k in enumerate(signals):
        yy=[r for r in stats if r['signal']==k and r['scope']=='whole']
        axes[0].plot(range(10),[r['repeat_fraction'] for r in yy],'o-',color=colors[j],label=labels[j],ms=3)
        with (out/'update_intervals.csv').open() as f:
            import csv
            ir=[r for r in csv.DictReader(f) if r['signal']==k and r['scope']=='whole']
        vals=[sum(int(r['count']) for r in ir if float(r['interval_s'])==gap) for gap in [.5,1,1.5,2,2.5]]
        other=sum(int(r['count']) for r in ir)-sum(vals);vals.append(other)
        axes[1].bar(np.arange(6)+(j-1)*.25,np.array(vals)/sum(vals),width=.24,color=colors[j],label=labels[j])
    axes[0].set_xticks(range(10),[z['heat_id'] for z in heats],rotation=55,ha='right');axes[0].set_ylabel('Equal adjacent values / pairs');axes[0].set_title('(a) Exact repeats in raw records');axes[0].legend(fontsize=7)
    axes[1].set_xticks(range(6),['0.5','1','1.5','2','2.5','Other']);axes[1].set_xlabel('Interval between value changes (s)');axes[1].set_ylabel('Fraction of change intervals');axes[1].set_title('(b) Recorded changes, pooled');axes[1].legend(fontsize=7)
    save(fig,'raw_updates','十炉原始记录的精确重复及变化间隔。变化间隔不等于设备采样周期，重复值也不单独证明输出保持机制。')
    fig,axes=plt.subplots(2,2,figsize=(7.2,5.2),sharex=True)
    for ax,e in zip(axes.flat,examples):
        keep=(e['t']>=185)&(e['t']<=215)
        ax.plot(e['t'][keep],e['truth'][keep],label='Known input',color='black')
        ax.plot(e['t'][keep],e['flow'][keep],label='Hampel + PCHIP + SG31',color=colors[0])
        ax.set_title(f"Pulse {e['duration']:g} s; counter update {e['period']:g} s");ax.set_xlabel('Time (s)');ax.set_ylabel('Flow (m3/s)');ax.legend(fontsize=7)
    save(fig,'controlled_flow','已知三角脉冲经过假设计数保持和完整正式重构链后的形态。显示的更新周期是实验条件，不是对真实设备的识别结果；无量化、相位为零。')
    controlled_summary=summary(lags,['shape','duration_s','counter_period_s','quantum_m3','gas_period_s'])
    write_csv(out/'controlled_delay_summary.csv',controlled_summary)
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.4))
    for ax,shape in zip(axes,['triangle','rectangle']):
        for p,color in zip(cfg['simulation']['counter_update_periods_s'],['#236394','#478055','#ba592b','#775397']):
            yy=[]
            for w in cfg['simulation']['durations_s']:
                rr=[r for r in lags if r['shape']==shape and r['duration_s']==w and r['counter_period_s']==p]
                yy.append(np.mean([r['correct'] for r in rr]))
            ax.plot(cfg['simulation']['durations_s'],yy,'o-',label=f'Counter {p:g} s',color=color,ms=3)
        ax.set_xlabel('Known pulse duration (s)');ax.set_ylabel('Accepted and within 1 s');ax.set_ylim(-.03,1.03);ax.set_title(shape.capitalize()+' input');ax.legend(fontsize=7)
    save(fig,'controlled_delay','假设计数与烟气保持条件下，正式A的指定中心正确恢复率。周期内每个离散相位枚举一次；烟气周期、相位、延迟和量化交叉，不将这些条件当作独立工艺事件。')
    bs=summary(background,['scenario','method','duration_s']);bh=summary(background,['heat_id','gas','scenario','method','duration_s'])
    write_csv(out/'background_summary.csv',bs);write_csv(out/'background_by_heat.csv',bh)
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.4))
    for method,offset,color in [('direct_flow',-.15,colors[0]),('raw_counter_chain',.15,colors[1])]:
        for j,scenario in enumerate(['known_addition','source_only']):
            yy=[next(r for r in bs if r['scenario']==scenario and r['method']==method and r['duration_s']==w) for w in cfg['background']['durations_s']]
            val=[r['correct_rate'] if j==0 else r['accepted']/r['conditions'] for r in yy]
            axes[j].bar(np.arange(2)+offset,val,width=.28,color=color,label=method.replace('_',' '))
            for x,v,r in zip(np.arange(2)+offset,val,yy):
                count=r['correct'] if j==0 else r['accepted']
                axes[j].text(x,v+(.012 if j==0 else .00035),f"{count}/{r['conditions']}",ha='center',fontsize=6,rotation=90)
    for ax in axes:ax.set_xticks(range(2),['4','16']);ax.set_xlabel('Added pulse duration (s)');ax.set_ylim(0,1);ax.legend(fontsize=7)
    axes[0].set_title('(a) Raw-background known additions');axes[0].set_ylabel('Accepted and within 1 s')
    axes[1].set_title('(b) Source-only controls');axes[1].set_ylabel('Acceptance with no added gas response')
    axes[0].set_ylim(0,1.2);axes[1].set_ylim(0,.02)
    save(fig,'background_chain','直接叠加流量与从原始累计量重新重构的配对对照。两者目标信号均从同一原始烟气叠加后清洗。仅源注入接受率不是工业误报率。')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.6))
    axes[0].bar(range(10),[100*r['gain_rate'] for r in ar],color=colors[0]);axes[0].set_xticks(range(10),[r['heat_id'] for r in ar],rotation=55,ha='right');axes[0].set_ylabel('Correct recovery gain (percentage points)');axes[0].set_title('(a) Previous A/B, per heat')
    axes[1].plot(range(10),[100*r['gain_rate'] for r in loo],'o',color=colors[0]);axes[1].axhline(100*boot['gain_rate'],ls='--',color='black');axes[1].set_xticks(range(10),[r['omitted_heat'] for r in loo],rotation=55,ha='right');axes[1].set_ylabel('Gain after omitting one heat (pp)');axes[1].set_title('(b) Leave-one-heat-out sensitivity')
    save(fig,'prior_cluster_sensitivity','只重新汇总3790d5f已有注入结果，不再研发B。每次留出整炉；Bootstrap同样以整炉为重采样单元，保留炉内配对和重复条件。')
    write_csv(out/'figure_captions.csv',captions);write_json(out/'font_validation.json',dict(configuration=font,pdf_checks=fontchecks))
    preview=Image.new('RGB',(900,3*350),'white')
    for i,r in enumerate(captions):
        im=Image.open(figdir/(r['figure']+'.png')).convert('RGB');im.thumbnail((450,345));preview.paste(im,((i%2)*450,(i//2)*350))
    preview.save(out/'layout_preview.png')
    # Recorded facts are independent of any assumed instrument model.
    pooled=[]
    for k,label in zip(signals,labels):
        rr=[r for r in stats if r['signal']==k and r['scope']=='whole'];den=sum(r['adjacent_pairs'] for r in rr)
        pooled.append(dict(signal=k,repeated_fraction=sum(r['repeated_pairs'] for r in rr)/den,
                           minimum_repeat_fraction=min(r['repeat_fraction'] for r in rr),maximum_repeat_fraction=max(r['repeat_fraction'] for r in rr)))
    write_json(out/'recorded_facts.json',pooled)
    with (out/'raw_increments.csv').open() as f:
        positive_increments=np.array([float(r['increment']) for r in csv.DictReader(f)
                                      if r['signal']=='oxygen_cumulative' and float(r['increment'])>0])
    oxygen_increment_stats={str(q):float(np.quantile(positive_increments,q)) for q in [0,.05,.25,.5,.75,.95,1.]}
    write_json(out/'oxygen_increment_quantiles.json',dict(quantiles=oxygen_increment_stats,unit='m3',meaning='recorded positive differences, not inferred count quanta'))
    model_features=[]
    for shape in cfg['simulation']['shapes']:
        for width in cfg['simulation']['durations_s']:
            for period in cfg['simulation']['counter_update_periods_s']:
                rr=[r for r in model if r['shape']==shape and r['duration_s']==width and r['counter_period_s']==period]
                result=dict(shape=shape,duration_s=width,counter_period_s=period,conditions=len(rr))
                for key in ['onset10_s_error','rise50_s_error','fall50_s_error','centroid_s_error','peak_set_distance_s','peak_ratio','area_relative_error','final_integral_error_relative']:
                    v=[r[key] for r in rr if r[key] is not None]
                    for name,fun in [('min',np.min),('median',np.median),('max',np.max)]:result[key+'_'+name]=float(fun(v)) if v else None
                model_features.append(result)
    write_csv(out/'controlled_timing_summary.csv',model_features)
    delay_summary=[]
    for shape in cfg['simulation']['shapes']:
        for width in cfg['simulation']['durations_s']:
            for period in cfg['simulation']['counter_update_periods_s']:
                rr=[r for r in lags if r['shape']==shape and r['duration_s']==width and r['counter_period_s']==period]
                result=dict(shape=shape,duration_s=width,counter_period_s=period,conditions=len(rr))
                for name in ['onset','peak','centroid']:
                    values=[r['feature_delay_'+name+'_error_s'] for r in rr if r['feature_delay_'+name+'_error_s'] is not None]
                    result[name+'_defined_conditions']=len(values)
                    result[name+'_error_min_s']=float(min(values)) if values else None
                    result[name+'_error_max_s']=float(max(values)) if values else None
                    result[name+'_absolute_error_p95_s']=float(np.quantile(abs(np.asarray(values)),.95)) if values else None
                delay_summary.append(result)
    write_csv(out/'controlled_feature_delay_summary.csv',delay_summary)
    # A complete chain must also audit whether a valid short gas peak is erased.
    sc=cfg['simulation'];dt=cfg['dt_s'];tt=np.arange(sc['samples'])*dt;gas_rows=[]
    for shape in sc['shapes']:
        for width in sc['durations_s']:
            for period in sc['gas_update_periods_s']:
                for phase in range(round(period/dt)):
                    for delay in sc['delays_s']:
                        y=sc['amplitude_m3_s']*pulse(tt,sc['source_center_s']+delay,width,shape)
                        obs=held(y,round(period/dt),phase);mask=hampel_mask(obs);clean=pchip_clean(tt,obs,mask)
                        maxima=tt[np.isclose(clean,clean.max(),rtol=1e-10,atol=1e-12)] if clean.max()>0 else []
                        gas_rows.append(dict(shape=shape,duration_s=width,gas_period_s=period,gas_phase_s=phase*dt,truth_delay_s=delay,
                                             flagged_points=int(mask.sum()),raw_peak=float(obs.max()),clean_peak=float(clean.max()),
                                             clean_peak_first_s=float(maxima[0]) if len(maxima) else None,
                                             clean_peak_last_s=float(maxima[-1]) if len(maxima) else None,
                                             clean_peak_sample_count=len(maxima),
                                             clean_to_measured_peak_ratio=float(clean.max()/obs.max()) if obs.max()>0 else None,
                                             **errors(tt,y,clean)))
    write_csv(out/'controlled_gas_cleaning.csv',gas_rows)
    sharp=next(r for r in gas_rows if r['shape']=='triangle' and r['duration_s']==2 and r['gas_period_s']==.5 and r['truth_delay_s']==15)
    representation=[]
    for shape in ['triangle','rectangle']:
        for width in [2.,4.,5.,15.,16.,30.]:
            known=pulse(tt,sc['source_center_s'],width,shape)
            difference=np.diff(known,prepend=known[0]);change=_standardized_change(known)
            representation.append(dict(shape=shape,duration_s=width,difference_energy=float(np.sum(difference**2)),
                                       retained_change_energy=float(np.sum(change**2)),
                                       retained_energy_ratio=float(np.sum(change**2)/np.sum(difference**2)),
                                       max_absolute_retained_change=float(np.max(abs(change)))))
    write_csv(out/'change_representation.csv',representation)
    # Metadata ledger states exactly which missing observations prevent identification.
    task=ROOT/'2026 APMCM Problem B.pdf'
    tracked=subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
    ledger=dict(problem_statement=dict(path=str(task.relative_to(ROOT)),sha256=sha(task),
                 established=['nominal adjacent sample spacing 0.5 s','gas travels through sampling path before measurement','Raman gas analyzer described'],
                 not_established=['pipe length / flow and calibrated transport delay','analyzer internal response time','independent acquisition timestamps / clock offset','PLC setpoints and commands']),
                local_original_inputs=[str(p.relative_to(ROOT)) for p in sorted((ROOT/'Data').glob('*')) if p.is_file()],
                tracked_operation_log_filename_matches=[p for p in tracked if any(s in Path(p).name.lower() for s in ['plc','clock','操作记录','启停记录','同步记录'])],
                workbook_schema='one sequence index, 2047 wavelength fields, cumulative oxygen, CO2, CO; no device timestamps or operator command fields',
                acquisition_mechanism='unidentified: repeated values admit held acquisition, quantized smooth input, and genuine plateaus',
                physical_delay_identifiability='transport, analyzer, reaction, and clock effects cannot be separately identified from these columns',
                precision='0.5 s is the lag grid, not a validated 0.5 s physical timing accuracy',
                known_processing_support=dict(SG31_past_s=7.5,SG31_future_s=7.5,Hampel9_past_s=2.,Hampel9_future_s=2.,
                                              change_median21_past_s=5.,change_median21_future_s=5.),
                support_interpretation='centered SG has no automatic 7.5 s causal delay; symmetric smoothing can shift threshold crossings and broaden short pulses',
                independent_events_confirmed_by_this_audit=0,event_protocol=cfg['event_definition'])
    write_json(out/'measurement_evidence.json',ledger)
    gap=max(abs(r['final_relative_error']) for r in flows)
    peak2=next(r for r in model_features if r['shape']=='triangle' and r['duration_s']==2 and r['counter_period_s']==.5)
    short={r['method']:r for r in bs if r['scenario']=='known_addition' and r['duration_s']==4}
    incremental=[r for r in timing if r['duration_s']==4]
    decision=dict(adopt_new_formal_model=False,physical_events_confirmed=0,physical_timing_accuracy_identified=False,
                  SG31_window_span_s=15.,SG31_future_support_s=7.5,
                  maximum_actual_final_integral_relative_error=gap,
                  short_triangle_peak_ratio=peak2['peak_ratio_median'],
                  four_second_background_direct_correct=short['direct_flow']['correct'],
                  four_second_background_full_chain_correct=short['raw_counter_chain']['correct'],
                  four_second_background_conditions=short['direct_flow']['conditions'],
                  conclusion=f"测量链审计找到了明确的软件失真：同样的4秒已知供氧脉冲，在真实背景中直接叠加流量的正确恢复为{short['direct_flow']['correct']}/{short['direct_flow']['conditions']}，从原始累计量重新重构则为{short['raw_counter_chain']['correct']}/{short['raw_counter_chain']['conditions']}。31点SG显著衰减和展宽短变化，21点变化去趋势能完全消除某些缓变三角供氧激励，Hampel也能误删构造的真实烟气尖峰。原始十炉还缺少设备时间戳、PLC事件及分析仪标定，因此不能给出物理时滞的校准精度，也不能拆分工艺、采样输送和时钟偏移。反积分一致性及0.5秒网格都不能证明时序准确。保留正式V2，不用候选窗口拟合物理流量—时滞函数。")
    write_json(out/'decision.json',decision)
    observed_table=mdtable(['信号','重复相邻值比例（按记录对数汇总）','逐炉最小','逐炉最大'],[[r['signal'],r['repeated_fraction'],r['minimum_repeat_fraction'],r['maximum_repeat_fraction']] for r in pooled])
    flowtable=mdtable(['炉号','氧量变化 m3','反积分相对误差','最大路径残差 m3','截负增加 m3'],[[r['heat_id'],r['oxygen_span_m3'],r['final_relative_error'],r['max_path_residual_m3'],r['zero_clipping_added_m3']] for r in flows])
    bt=mdtable(['场景','入口','宽度 s','条件数','接受','±1s正确','错误时滞','拒绝比例','接受后MAE s'],[[r[k] for k in ['scenario','method','duration_s','conditions','accepted','correct','wrong','reject_rate','conditional_mae_s']] for r in bs])
    mt=mdtable(['形态','宽度 s','计数周期 s','起始10%误差范围 s','50%上升误差范围 s','峰值集合距离最大 s','重心误差范围 s','峰值比中位'],[[r['shape'],r['duration_s'],r['counter_period_s'],f"{r['onset10_s_error_min']:.3f} to {r['onset10_s_error_max']:.3f}",f"{r['rise50_s_error_min']:.3f} to {r['rise50_s_error_max']:.3f}",r['peak_set_distance_s_max'],f"{r['centroid_s_error_min']:.3f} to {r['centroid_s_error_max']:.3f}",r['peak_ratio_median']] for r in model_features])
    abt=mdtable(['炉号','A正确率','B正确率','净增加','改善百分点'],[[r['heat_id'],r['A_rate'],r['B_rate'],r['net_correct'],100*r['gain_rate']] for r in ar])
    report=f'''# 原始测量链与时滞可辨识性审计

{decision['conclusion']}

## 1. 输入与保护边界

正式基线d08f989，上一轮研究3790d5f。十个XLSB SHA256与原清单逐一一致；读取全部实际表头和原始工艺列，与冻结NPZ逐点精确核对。30组Hampel掩码和PCHIP结果、十炉SG流量重放一致。所有原有受Git跟踪文件以及100个冻结文件在执行前后核对哈希；本轮只在独立measurement_chain目录写入。第一问、正式V2、第三问接口和冻结论文均不修改。

本轮没有开发新的事件算法。实验比较的时滞器始终是正式A，生产参数和阈值不变。合成采集模型、原始背景叠加和实际记录统计分别报告。

## 2. 原始记录揭示的更新特征

{observed_table}

`raw_signal_summary.csv`逐炉报告全记录及首末正氧增量之间的活动区间，避免停吹平台单独支配结论。`update_intervals.csv`统计相邻非零变化之间的间隔；`raw_increments.csv`保留每一次正、负变化；`constant_runs.csv`明确区分相等转换持续时间(n−1)×0.5秒与网格占用n×0.5秒，并标注两端截断。重复用精确相等判定，无人为容差合并。

多数炉次耗氧量的变化间隔众数为1秒，同时存在0.5、1.5秒等间隔；CO/CO₂也不是每个0.5秒点都变化。`update_autocorrelation.csv`给出变化指示序列0.5–20秒的相关，`representation_grids.csv`检验保存值的小数网格。它们分别描述记录节律和数值表示，不能据此唯一识别设备周期或计数分辨率。变化大小随过程输入变化，正增量的中位数也不能解释为固定计数量子。

正耗氧增量的中位数为{oxygen_increment_stats['0.5']:.2f}m3，5%–95%范围为{oxygen_increment_stats['0.05']:.2f}–{oxygen_increment_stats['0.95']:.2f}m3，最小/最大为{oxygen_increment_stats['0']:.2f}/{oxygen_increment_stats['1.0']:.2f}m3。三个工艺列保存值均落在0.01的数值网格上；这是存储表示证据，不能等同于设备真实分辨率。

`update_coincidence.csv`记录跨列同时变化的数量及Jaccard；同一行索引不保证设备同时测量，同步变化也可能来自共享导出节律。仅靠这三列，固定保持、平滑量经量化和真实平台都能产生重复值。本轮结论是“存在记录更新限制”，没有把一种仪表内部机制当成事实。

## 3. 时间轴、测量延迟和证据清单

实际2051列全部核验：序列、2047个数值波长及三个工艺变量。没有设备绝对时间戳、时间质量标志、氧枪开关、流量设定或PLC指令列。当前数据入口以t=0.5(Sequence−1)构造炉内时间。题面支持0.5秒记录间距，也明确提到烟气输送采样延迟；没有提供管长、输送流量、分析仪响应标定或共钟证书。检索本地Data文件及受跟踪路径，没有取得独立PLC/时钟标定记录。详见`raw_headers.json`和`measurement_evidence.json`。没有凭通用仪器资料给本项目虚构固定输送延迟。

必须区分名义记录间隔、数值变化间隔、仪器采集周期和时滞估计精度。前两者可从题面和文件计算，后两者不能由重复值唯一反演。若记录形如y(t)=s(t−τ−δ)，工艺/输送延迟τ与设备时钟偏移δ只以和出现；缺少独立校准时，增加互相关算法不能将它们拆开。即使找到较高相关峰，也只能报告复合的统计响应时差。

31点SG的首尾跨度是15秒，中心两侧各用7.5秒；这不是固定7.5秒输出延迟。9点Hampel两侧各用2秒，变化信号的21点居中中位数两侧各用5秒；PCHIP在剔除区间依赖两侧有效点，实际修复跨度见`repair_spans.csv`。这些是离线处理的数据依赖范围，不能直接相加当作物理时滞误差条，也不能据此宣称真实测量链有10秒误差。

## 4. 流量与累计量的一致性

{flowtable}

反积分最大末端相对误差绝对值为{gap:.6g}。`flow_integral.csv`另报全程残差和截负贡献。这只验证总量与曲线计算的一致性，不验证供氧变化时刻。首末正增量也不是氧枪启停记录：记录可能从已供氧状态开始，首末样本之间的短间隔变化无法单独定义开关事件。

## 5. 已知测量链条件下的时序传递

不再复跑理想平移算法基准，而专门将已知流量积分为O，然后经过假设计数采集、原Hampel、PCHIP、31点SG；气体先加入已知延迟，再经过假设保持与原清洗。三角与矩形输入宽度2、5、15、30秒，幅值5m3/s。计数更新周期0.5、1、2、5秒，枚举每个离散相位；量化为0或0.01m3。烟气周期0.5、1、2秒，枚举相位；延迟5、15、30秒。所有周期是预先设定的压力条件，不声称是真实设备参数；没有根据真实时滞结果择优选参数。

共有{len(model)}种源链重构条件和{len(lags)}种源—气体采集时滞条件，无随机重复伪增样本。起始定义为各自峰幅10%的上升穿越，另报50%上升/下降、半高宽、正信号重心；这些是明确的操作定义，不等于仪器检出限。矩形输入有一段相同最大值，峰值误差以估计峰到真实最大值集合的最短距离计算，不能拿任选的矩形峰时刻作唯一真值。所有拒绝保留。完整互相关时滞数值误差、拒绝率见`controlled_delay.csv`及分层汇总。同表的`feature_delay_onset/peak/centroid_s`直接计算重构源与气体的三种时间差及其对真值的误差，独立于互相关是否接受；矩形的单一峰时滞留空。`controlled_feature_delay_summary.csv`汇总三者误差，不能拿重心准确代替起始准确。

{mt}

最短2秒三角响应即使每点更新、无额外保持，其重构峰值比中位数仅{peak2['peak_ratio_median']:.4f}。这直接证明“累计量积分正确”与“短时流量形态正确”不是同一验收指标。居中平滑能在真实输入前产生正响应，起始偏早并不说明实际先发生反应；重心、峰和阈值起始也不能混用。时间指标相对各自幅值计算，幅值衰减另列，避免隐藏短信号损失。

完整链也不能忽略烟气异常处理。`controlled_gas_cleaning.csv`逐条记录假设采集后、原Hampel/PCHIP前后的峰值：2秒三角脉冲在0.5秒逐点更新、15秒延迟的控制条件下，真实有效脉冲有{sharp['flagged_points']}点被标异常，清洗峰值只剩观测峰的{sharp['clean_to_measured_peak_ratio']:.3f}。这是异常检测与短真实瞬态发生混淆的构造证据，不是对原始工业异常标签的重新认定，也不能由此自动放宽正式Hampel。清洗还会将唯一三角峰变成平台，该例清洗最大值位于{sharp['clean_peak_first_s']:.1f}–{sharp['clean_peak_last_s']:.1f}秒，真峰为215秒；时间特征中的估计峰使用最早最大值这一明确约定，不能将其当成唯一物理峰，CSV另报近似并列最大值范围。源链另存`SG_only_vs_full_max_difference_m3_s`，用于区分SG自身与氧量清洗引起的差异。

## 6. 从真实原始O开始的完整软件链注入

每炉—气体—30固定种子抽取一个有效中心，复用4/16秒三角脉冲与15/30秒延迟。源幅值为该炉冻结流量中位数，气体幅值为原始气体全程标准差。源脉冲积分后加到原始O，再跑原Hampel/PCHIP/SG；气体脉冲加到原始气体，再跑原清洗。配对入口对照直接将同一已知流量脉冲加到冻结q，目标端使用同一清洗序列。气体延迟真值只指添加分量，不为未知原背景贴标签。另做仅源注入对照，每个中心和宽度一次，不随延迟重复计算。

这覆盖当前软件处理链，仍不包含未知实际设备重新测量、时钟、烟气运输或分析仪传递函数。没有人为把5秒保持强加到真实仪器。600个炉—气体—种子抽样背景重复用于多条件，抽样有放回，不保证600个不同时间中心，也不作为600个独立物理事件。

{bt}

4秒配对条件中，从直接流量入口的{short['direct_flow']['correct']}/{short['direct_flow']['conditions']}正确恢复降至完整源链的{short['raw_counter_chain']['correct']}/{short['raw_counter_chain']['conditions']}；气体输入和清洗相同，差异来自源表示经过累计量积分和正式重构。源增量的峰值比中位数{np.median([r['peak_ratio'] for r in incremental]):.4f}、起始误差中位数{np.median([r['onset10_s_error'] for r in incremental]):.3f}秒、半高宽误差中位数{np.median([r['width50_s_error'] for r in incremental]):.3f}秒，源Hampel掩码改变为{sum(r['changed_oxygen_mask_points'] for r in incremental)}。这给出了源重构使短输入难以与未同样平滑的烟气匹配的受控证据；它没有证明所有原始工艺变化都属于4秒脉冲。16秒脉冲两种入口本来都恢复很少，说明恢复也受去趋势变化表示和背景条件制约，不能将全部失败归给SG。

为核对16秒失败是否只是噪声，另对无背景、未经过SG的已知输入直接审计正式变化表示，不新增识别方法。`change_representation.csv`记录差分能量以及减去21点滚动中位数后的能量：4秒三角输入保留比为1，16秒三角输入为0，15/30秒也降到数值舍入量级；矩形输入的边缘差分则保留。这不是激励不存在，而是三角信号的分段常值差分被滚动中位数当作趋势消掉。它证明正式表示偏向尖锐边缘，不能覆盖所有真实供氧变化形态。真实背景中还存在更新噪声，因此不能据此认定任何一个具体原始候选就是该机制产生。

这里“正确”严格指接受且距添加延迟≤1秒；接受后MAE不能代替正确率。`background_source_timing.csv`比较加脉冲前后流量之差与真脉冲，记录掩码变化。该差可以有负值，时间特征使用正部，面积/重心的定义与控制条件一致；不能将负部截掉后的面积当作净耗氧守恒量。原始记录与合成控制的结果不得互相替代。

## 7. 上轮A/B的炉次异质性与聚类统计

本节只读取已保存结果，不再研究B。7200是配对条件，600是有放回背景抽样单元，12种条件共享背景；Bootstrap将整炉所有气体、种子和条件一起重采样，A/B不拆开。以10炉为聚类单位，固定种子、10000次百分位Bootstrap，改善率为{100*boot['gain_rate']:.3f}个百分点，描述性95%区间[{100*boot['percentile_95_lower']:.3f}, {100*boot['percentile_95_upper']:.3f}]个百分点。这个区间仅针对注入设计中的炉次差异，不是物理事件准确率的置信区间，10炉也不能代表所有工业工况。

{abt}

留一炉改善率范围{100*min(r['gain_rate'] for r in loo):.3f}–{100*max(r['gain_rate'] for r in loo):.3f}个百分点，见`prior_AB_leave_one_heat_out.csv`。逐炉最差A、B恢复率分别为{min(r['A_rate'] for r in ar):.4f}、{min(r['B_rate'] for r in ar):.4f}。保留全部炉次，不用总体均值掩盖改善集中现象。正式V2仍不替换。

## 8. 独立事件与下一步决策

本轮确认的独立物理事件为0。累计量数值更新、SG导数峰和滑动窗口都不直接授予事件身份。`config.json`给出后续事件约定：需要独立PLC命令/实际流量记录、可校准时间戳，再记录输入启停/设定变化、气体起始/峰/重心、各自定义、时间区间、输入隔离情况及匹配依据；允许unmatched、ambiguous、candidate及confirmed_by_independent_record。PLC命令时刻与实际流量变化时刻也必须区分，仅有命令时估计的是包含执行机构响应在内的命令—烟气时差。没有记录时维持statistical_candidate，不强制一对一指派，不开始新事件算法。

下一步的必要输入是同炉PLC实际流量或启停/设定时间、采集时钟对齐资料、采样输送和分析仪动态标定。获得这些输入后才能给物理时滞设置校准误差，并决定1秒级估计是否有意义；仍需结合事件隔离和工况信息检验供氧变化的归因，补齐时钟本身不能证明因果。当前十炉可以继续支持窗口级统计时差，不能证明流量增加必然使独立物理时滞减小。本轮停止在测量链审计，不降低0.35、不扩大事件集、不拟合新流量函数。

## 9. 复现与输出

在项目根目录运行 `.venv/bin/python -m problem2.research.measurement_chain.run`，固定目录原位更新。原工程环境用于准确复用生产scipy算法及pyxlsb读取（文档运行时不含这两个依赖）。完整参数、输入哈希、源代码哈希、执行提交和输出哈希分别保存。图内英文与真实Times New Roman，PNG/SVG/PDF齐全，PDF嵌入逐图检查，图题另存`figure_captions.csv`，一张当前检查页。所有原始数据和此前实验目录只读。
'''
    (out/'REPORT.md').write_text(report,encoding='utf-8')
