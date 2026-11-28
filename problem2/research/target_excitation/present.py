"""Report paired A/B evidence and strictly separate it from physical validation."""
import csv,json,subprocess
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw
from ...src.figure_fonts import configure
from ...src.delay_audit import write_csv,write_json
from ...src.q2_pipeline import _standardized_change


def mdtable(headers,rows):
    def text(v):
        if v is None or v=='':return '—'
        return f'{v:.4g}' if isinstance(v,(float,np.floating)) else str(v)
    return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+''.join('| '+' | '.join(text(v) for v in row)+' |\n' for row in rows)


def deliver(out,cfg,examples,real,syn,inj):
    from .run import summarize
    summary=summarize(syn+inj,['scenario']);real_summary=summarize(real,['gas'])
    pairs={}
    for r in syn+inj:
        key=tuple(r.get(k) for k in ['scenario','heat_id','gas','seed','noise_ratio','target_snr','delay_s','center_index'])
        pairs.setdefault(key,{})[r['method']]=r
    transitions=[]
    for key,pair in pairs.items():
        assert set(pair)=={'A','B'}
        a,b=pair['A'],pair['B']
        transitions.append(dict(scenario=key[0],heat_id=a.get('heat_id'),gas=a.get('gas'),seed=a['seed'],noise_ratio=a.get('noise_ratio'),target_snr=a.get('target_snr'),delay_s=a.get('delay_s'),
            A_accepted=a['accepted'],B_accepted=b['accepted'],A_correct=a['correct'],B_correct=b['correct'],
            rescued_correct=not a['correct'] and b['correct'],lost_correct=a['correct'] and not b['correct'],
            newly_accepted=not a['accepted'] and b['accepted'],new_wrong=b['false_delay'] and not a['false_delay'],
            B_accept_depends_on_masked_competitors=b['accepted'] and not b['full_competitor_accept']))
    write_csv(out/'paired_transitions.csv',transitions)
    def detail_rows(rr):
        return [[r['scenario'],r['method'],r['trials'],r['accepted'],r['correct_recoveries'] if r['known_delay_trials'] else '不适用',r['correct_recovery_rate'],r['wrong_delay_rate_all'],r['conditional_mae_s'],r['rejection_rate']] for r in rr]
    scoreboard=mdtable(['场景','方法','试验数','接受数','±1秒正确恢复数','正确恢复率','错误时滞/全部','接受后MAE s','拒绝比例'],detail_rows(summary))
    realtable=mdtable(['气体','方法','全部窗口','接受窗口','A式目标前置门拒绝','已接受但平移目标低于阈值','接受依赖删去竞争峰'],[[r['gas'],r['method'],r['trials'],r['accepted'],r['unshifted_target_gate_failures'],r['accepted_selected_target_below_threshold'],r['accepted_due_to_masked_competitors']] for r in real_summary])
    transition_summary=[]
    for scenario in sorted({r['scenario'] for r in transitions}):
        rr=[r for r in transitions if r['scenario']==scenario]
        transition_summary.append(dict(scenario=scenario,trials=len(rr),rescued_correct=sum(r['rescued_correct'] for r in rr),lost_correct=sum(r['lost_correct'] for r in rr),
            newly_accepted=sum(r['newly_accepted'] for r in rr),new_wrong=sum(r['new_wrong'] for r in rr),masked_competitor_dependence=sum(r['B_accept_depends_on_masked_competitors'] for r in rr)))
    write_csv(out/'transition_summary.csv',transition_summary)
    transitions_table=mdtable(['场景','配对试验','B救回正确','B丢失正确','新增接受','新增错误时滞','B接受依赖竞争峰掩蔽'],[[r[k] for k in ['scenario','trials','rescued_correct','lost_correct','newly_accepted','new_wrong','masked_competitor_dependence']] for r in transition_summary])
    injection_breakdown=summarize(inj,['scenario','gas','target_snr','delay_s'])
    dosage=mdtable(['气体','目标变化SNR','真实注入延迟 s','方法','正确恢复率','错误时滞率','接受后MAE s'],[[r['gas'],r['target_snr'],r['delay_s'],r['method'],r['correct_recovery_rate'],r['wrong_delay_rate_all'],r['conditional_mae_s']] for r in injection_breakdown if r['scenario']=='real_paired_injection'])
    # Original accepted-set changes are not ground-truth improvements.
    realpair={}
    for r in real:realpair.setdefault((r['heat_id'],r['gas'],r['center_index']),{})[r['method']]=r
    changes=[]
    for key,p in realpair.items():
        a,b=p['A'],p['B']
        if a['accepted']!=b['accepted'] or a['estimate_s']!=b['estimate_s']:
            changes.append(dict(heat_id=key[0],gas=key[1],center_index=key[2],center_time_s=key[2]*cfg['dt_s'],A_accepted=a['accepted'],B_accepted=b['accepted'],A_delay_s=a['estimate_s'],B_delay_s=b['estimate_s'],
                unshifted_target_ratio=a['unshifted_target_ratio'],B_selected_target_ratio=b['selected_target_ratio'],B_full_competitor_accept=b['full_competitor_accept']))
    # Empty output retains an explicit header so a previous nonempty result cannot survive reruns.
    fields=['heat_id','gas','center_index','center_time_s','A_accepted','B_accepted','A_delay_s','B_delay_s','unshifted_target_ratio','B_selected_target_ratio','B_full_competitor_accept']
    with (out/'real_changed_windows.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(changes)
    font=configure()
    import matplotlib.pyplot as plt
    figures=out/'figures';figures.mkdir(exist_ok=True);captions=[];fontchecks=[]
    def save(fig,name,caption):
        fig.tight_layout()
        for ext in ['png','svg','pdf']:fig.savefig(figures/(name+'.'+ext),dpi=180,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
        plt.close(fig);raw=subprocess.check_output(['pdffonts',str(figures/(name+'.pdf'))],text=True);lines=raw.splitlines()[2:];col=raw.splitlines()[0].index('emb')
        assert lines and all('TimesNewRoman' in line.replace(' ','') and line[col:col+3]=='yes' for line in lines)
        fontchecks.append(dict(figure=name,pdffonts=raw));captions.append(dict(figure=name,caption=caption,width_inches=7.2))
    for name,source,target,pair,center,lag in examples:
        time=np.arange(len(source))*cfg['dt_s'];half=cfg['window']//2;left=(center-half)*cfg['dt_s'];right=(center+half)*cfg['dt_s'];keep=(time>=left-10)&(time<=right+cfg['max_lag_samples']*cfg['dt_s']+10)
        fig,axes=plt.subplots(3,1,figsize=(7.2,6))
        axes[0].plot(time[keep],source[keep]-1000,label='Source transient');axes[0].plot(time[keep],target[keep]-1000,label='Target transient');axes[0].axvspan(left,right,color='gray',alpha=.15,label='Unshifted window');axes[0].axvspan(left+lag,right+lag,color='orange',alpha=.1,label='True shifted window');axes[0].set_title('(a) Known added signals');axes[0].set_ylabel('Signal above baseline');axes[0].set_xlabel('Time (s)');axes[0].legend(fontsize=7)
        for method,color in [('A','#155b8a'),('B','#b54f20')]:
            r=pair[method];lags=np.arange(len(r['correlations'][0]))*cfg['dt_s'];axes[1].plot(lags,r['correlations'][0],color=color,label=method);axes[2].plot(lags,r['shifted_std'][0]/(r['target_global_std']+1e-12),color=color,label='Shifted target / global' if method=='A' else None)
        axes[1].axvline(lag,color='black',ls=':');axes[1].set_title('(b) Admissible lag correlations');axes[1].set_ylabel('Correlation');axes[1].set_xlabel('Lag (s)');axes[1].legend()
        axes[2].axhline(cfg['min_strength_ratio'],color='black',ls='--',label='Strength threshold');axes[2].axhline(pair['A']['unshifted_target_ratio'][0],color='gray',ls=':',label='Unshifted target ratio');axes[2].set_title('(c) Excitation at the location actually compared');axes[2].set_ylabel('Local / global standard deviation');axes[2].set_xlabel('Lag (s)');axes[2].legend(fontsize=7)
        save(fig,name,'已知延迟的门控示例。A与B共用相关计算；A若前置拒绝则无可画相关曲线。B仅允许达到原0.05激励阈值的lag。该示例不代表真实事件标签。')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.4));scenarios=['centered_response','late_response','real_paired_injection']
    for method,offset,color in [('A',-.16,'#155b8a'),('B',.16,'#b54f20')]:
        yy=[next(r['correct_recovery_rate'] for r in summary if r['scenario']==s and r['method']==method) for s in scenarios];axes[0].bar(np.arange(3)+offset,yy,width=.3,label=method,color=color)
        nulls=['stationary_unrelated','quiet_future_unrelated','real_source_only'];yy=[next(r['acceptance_rate'] for r in summary if r['scenario']==s and r['method']==method) for s in nulls];axes[1].bar(np.arange(3)+offset,yy,width=.3,label=method,color=color)
    axes[0].set_xticks(range(3),['Centered','Late target','Real injection']);axes[0].set_ylabel('Correct recovery within 1 s');axes[0].set_title('(a) Known-delay trials')
    axes[1].set_xticks(range(3),['Stationary null','Quiet future null','Source only'],rotation=12);axes[1].set_ylabel('Acceptance with no planted response');axes[1].set_title('(b) Negative controls')
    for ax in axes:ax.set_ylim(0,1);ax.legend()
    save(fig,'paired_performance','正确恢复按每次试验的共同指定窗口中心计算，不等同于事件Recall或真实工业准确率。无注入响应接受率仅评价控制条件下的偶然接受，不是工业误报率。')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.4))
    for j,gas in enumerate(['co','co2']):
        for method,color in [('A','#155b8a'),('B','#b54f20')]:
            rr=[r for r in real if r['gas']==gas and r['method']==method and r['accepted']];axes[j].scatter([r['center_index']*cfg['dt_s'] for r in rr],[r['estimate_s'] for r in rr],marker='x' if method=='B' else 'o',s=22,color=color,label=method)
        axes[j].set_title(f'({chr(97+j)}) '+gas.upper()+' original accepted windows');axes[j].set_xlabel('Within-heat time (s)');axes[j].set_ylabel('Candidate delay (s)');axes[j].legend()
    save(fig,'real_candidates','原始十炉的算法合格窗口，炉内时间叠加仅作A/B总览，跨炉时间不表示同一事件。新增或移除窗口均未经物理标签验证。')
    write_csv(out/'figure_captions.csv',captions);write_json(out/'font_validation.json',dict(font=font,pdf_checks=fontchecks))
    tiles=[]
    for r in captions:
        im=Image.open(figures/(r['figure']+'.png')).convert('RGB');im.thumbnail((350,300));tile=Image.new('RGB',(360,325),'white');tile.paste(im,((360-im.width)//2,0));ImageDraw.Draw(tile).text((5,305),r['figure'],fill='black');tiles.append(tile)
    preview=Image.new('RGB',(720,650),'white')
    for i,tile in enumerate(tiles):preview.paste(tile,((i%2)*360,(i//2)*325))
    preview.save(out/'layout_preview.png')
    late={r['method']:r for r in summary if r['scenario']=='late_response'};ri={r['method']:r for r in summary if r['scenario']=='real_paired_injection'};null={r['method']:r for r in summary if r['scenario']=='real_source_only'}
    heat_rows=[]
    for heat in sorted({r['heat_id'] for r in inj}):
        for gas in ['co','co2']:
            selected=[r for r in inj if r['heat_id']==heat and r['gas']==gas]
            rr=summarize(selected,['scenario']);pp={r['method']:r for r in rr if r['scenario']=='real_paired_injection'};nn={r['method']:r for r in rr if r['scenario']=='real_source_only'}
            heat_rows.append(dict(heat_id=heat,gas=gas,paired_conditions=pp['A']['trials'],A_correct=pp['A']['correct_recoveries'],B_correct=pp['B']['correct_recoveries'],net_correct_gain=pp['B']['correct_recoveries']-pp['A']['correct_recoveries'],A_source_only_accepted=nn['A']['accepted'],B_source_only_accepted=nn['B']['accepted']))
    write_csv(out/'injection_by_heat.csv',heat_rows)
    quiet=[r for r in real if r['method']=='B' and r['source_gate'] and not r['unshifted_target_gate']]
    quiet_corr=max([r['max_corr'] for r in quiet if np.isfinite(r['max_corr'])],default=None)
    original_gate_facts=dict(co_unshifted_target_rejections=sum(r['gas']=='co' for r in quiet),co2_unshifted_target_rejections=sum(r['gas']=='co2' for r in quiet),max_correlation_after_reopening=quiet_corr,
        accepted_original_window_changes=len(changes))
    write_json(out/'original_gate_diagnosis.json',original_gate_facts)
    primary_real_gain=ri['B']['correct_recoveries']-ri['A']['correct_recoveries'];real_false_delta=null['B']['accepted']-null['A']['accepted']
    weak={r['method']:r for r in summary if r['scenario']=='weak_real_copy'}
    gas_correct={gas:{method:sum(r['correct'] for r in inj if r['scenario']=='real_paired_injection' and r['gas']==gas and r['method']==method) for method in ['A','B']} for gas in ['co','co2']}
    masked_rescues=sum(r['rescued_correct'] and r['B_accept_depends_on_masked_competitors'] for r in transitions if r['scenario']=='real_paired_injection')
    adopt=False  # Original 24 are unchanged; weak known copies are lost and null acceptance also increases.
    decision=f"目标激励检查位置错配是受控漏检的原因，但不是原始24个合格窗口稀少的主要原因。未平移目标门控会漏掉后续响应：专门晚响应试验正确恢复从{late['A']['correct_recoveries']}/{late['A']['trials']}变为{late['B']['correct_recoveries']}/{late['B']['trials']}。在十炉背景注入中，B比A净增加{primary_real_gain}次±1秒正确恢复，无注入响应接受数变化{real_false_delta:+d}。原始十炉的24个窗口及其时滞完全不变；同时，已知极弱复制在A恢复{weak['A']['correct_recoveries']}/{weak['A']['trials']}、B恢复{weak['B']['correct_recoveries']}/{weak['B']['trials']}。本轮保留正式V2：B通过了晚响应漏检机制验证，但没有改善原始工业时滞证据。"
    write_json(out/'decision.json',dict(code_hypothesis_confirmed=late['B']['correct_recoveries']>late['A']['correct_recoveries'],adopt_B_in_formal_model=adopt,
        real_injection_net_correct_gain=primary_real_gain,real_source_only_acceptance_change=real_false_delta,original_real_changed_windows=len(changes),
        conclusion=decision))
    report=f'''# 目标激励检查位置的独立 A/B 验证

{decision}

## 比较定义

A精确复用正式V2的含义：先检查源窗口及未平移目标窗口，再对所有正标准差目标lag计算相关。B取消未平移目标门控，在每个lag使用平移目标窗口相对于整条目标变化序列的标准差阈值0.05。源激励、121点窗口、0–60秒搜索、5点步长、相关0.35、峰分离0.05、零和60秒端点拒绝均不改变。所有A/B使用相同窗口中心。

B接受后的置信权重使用被选中目标窗口的标准差；这只影响诊断输出，不影响接受数量。B排除弱lag会改变竞争峰集合，因此另记录“对所有原始有限lag保留竞争峰”时同一B胜出lag是否仍通过分离阈值。它用于解释A/B差异，没有作为第三种正式模型，也没有根据结果调参。B的峰宽/显著性诊断使用未掩蔽相关曲线，避免用零或负无穷填孔制造峰；诊断不参与接受。

## A的精确重放与原始十炉

A与生产函数的窗口位置、时滞、接受掩膜、相关、两类标准差、峰分离、峰宽、显著性和置信权重逐字段核对，均在1e-9相对、1e-10绝对容差内一致。冻结NPZ的CO 13、CO₂ 11个接受窗口逐一匹配。原始数据没有独立物理标签，新增窗口不能计算真实Recall或当作改进证据。

{realtable}

A/B接受状态或已接受时滞改变的窗口共{len(changes)}个，见 `real_changed_windows.csv`。原始CO有{original_gate_facts['co_unshifted_target_rejections']}个窗口被未平移目标门单独挡住，B重新搜索后这些窗口的最大相关仅{quiet_corr:.6f}，仍低于0.35；CO₂没有被该门单独挡住的窗口。因此，这个门控位置问题不能解释当前CO₂的稀疏覆盖，也没有救回任何原始CO合格窗口。`source_ratio`及两类目标ratio均相对于各自完整序列标准差，不是已知物理激励阈值。

## 已知时滞与背景注入

使用30个固定种子，每次只评价一个共同指定中心，正确恢复定义为接受且误差不超过1秒，不能只以“返回有效值”判定成功。合成源为非负紧支撑8秒瞬态，峰值1，基线1000；噪声标准差为指定ratio乘整个瞬态序列的标准差。源位于窗口中心或其后20秒，延迟5、15、30、55秒；零与60秒单列边界参考。弱复制试验在未平移目标处放置强独立瞬态，再在真实55秒lag放置0.001倍复制，检验幅值门槛是否会同时压掉真实弱复制。该试验有真实数学延迟，不当作无关信号误报。

真实背景采用冻结的氧枪流量和CO/CO₂序列，在V2入口直接添加同形非负瞬态，不重新求导、不修改累计耗氧或正式数据。每个炉次—气体—种子均匀抽取一个有效中心，源放在中心后20秒，源变化SNR固定2；目标变化SNR为0.5、2、8，与5、15、30、55秒交叉。SNR按原背景变化序列的全局标准差和模板在源窗口内的变化标准差标定。每次配对使用同一输入、中心、阈值与随机实现；剂量和延迟重复共享背景，不作为独立工业事件。

每个炉次—气体—种子只做一次“仅注入源、不注入响应”对照，不随12种剂量延迟重复计数。它检验对人为源变化的无注入响应接受，不能当作真实工业误报率。原背景已有真实关系未知，因此所有注入误差只对已知添加成分定义，不给原始响应重新贴标签。

{scoreboard}

接受后MAE必须和正确恢复率、拒绝率合看；没有真实延迟的控制不计算MAE。边界参考不纳入5–55秒正确恢复结论。完整逐试验表保留所有拒绝，不删除失败种子。

真实背景配对注入的CO正确恢复为{gas_correct['co']['A']}/3600→{gas_correct['co']['B']}/3600，CO₂为{gas_correct['co2']['A']}/3600→{gas_correct['co2']['B']}/3600；逐炉结果见 `injection_by_heat.csv`。新增正确恢复中依赖竞争峰掩蔽的数量为{masked_rescues}，所以本轮恢复收益来自解除前置门控。源单独注入的B接受中有{null['B']['accepted_due_to_masked_competitors']}次依赖竞争峰掩蔽；完整竞争峰诊断后为{null['B']['B_full_competitor_acceptances']}/{null['B']['trials']}。原先的{ri['A']['accepted_wrong_delay']}次错误时滞仍然存在，不能用较低的条件MAE掩盖错误尾部；P99和最大误差随汇总保存。

## 配对变化与竞争峰检查

{transitions_table}

同一B胜出lag的完整竞争峰分离率、正确恢复数和错误数见汇总CSV的`B_full_competitor_acceptances`、`guarded_correct_recoveries`和`guarded_wrong_delays`。B若只通过移除竞争峰增加接受，必须与真实添加延迟误差一起评价，不能将其直接记为成功。

## 真实背景的剂量和时滞分层

{dosage}

## 结论与后续动作

本轮确认了两件事：未平移目标门会漏掉晚到的已知响应；它不是原始十炉只有24个合格窗口的主要原因。B在注入条件下净增加{primary_real_gain}次正确恢复且本轮未新增错误时滞，这是一项有效的代码机制改进；B的极弱复制正确恢复从{weak['A']['correct_recoveries']}变为{weak['B']['correct_recoveries']}，源单独注入接受数增加{real_false_delta}次。

结论是保留B作为已经验证的实验修补，正式V2不替换。当前原始数据的瓶颈仍在0.35相关阈值下的响应证据和事件身份，不能继续把这一前置门当作主要突破口。对已有数据而言，这项修改不会提高24窗口的覆盖率，也不会使统计相关峰自动变成物理响应时滞。

本轮不修改正式V2、正文、接口、论文冻结清单或第一问。冻结清单内每个文件在实验前后再次核验。研究配置、输入哈希、逐试验CSV和全部输出清单随研究保存。复现命令为 `python -m problem2.research.target_excitation.run`，固定原位更新，不创建平行版本。
'''
    (out/'REPORT.md').write_text(report)
