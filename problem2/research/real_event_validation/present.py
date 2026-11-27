"""Fixed review pages, editable captions, report and font verification."""
import csv,html,json,subprocess
from pathlib import Path
import numpy as np
from PIL import Image,ImageOps,ImageDraw
from ...src.figure_fonts import configure
from ...src.delay_audit import write_json,write_csv
from .run import ROOT,PRIOR,RESEARCH,readcsv,events_at


def render(out,cfg):
    font=configure()
    import matplotlib.pyplot as plt
    cases=json.loads((out/'review_cases.json').read_text());summary=json.loads((out/'summary.json').read_text())
    rows=readcsv(out/'heat_summary.csv');scales=readcsv(out/'scale_sensitivity.csv');cadence=readcsv(out/'sampling_cadence.csv')
    pairs=readcsv(out/'previous_unique_pairs.csv');controls=readcsv(out/'circular_controls.csv');sources=readcsv(out/'source_candidates.csv')
    figures=out/'figures';figures.mkdir(exist_ok=True);captions=[];pdfs=[]
    base=json.loads((ROOT/cfg['event_configuration']).read_text())
    def save(fig,name,caption):
        fig.tight_layout()
        for ext in ['png','svg','pdf']:
            fig.savefig(figures/(name+'.'+ext),dpi=180,metadata={'Creator':'APMCM real event audit'} if ext=='pdf' else None)
        plt.close(fig);captions.append(dict(figure=name,caption=caption));pdfs.append(figures/(name+'.pdf'))
    for case in cases:
        heat,gas=case['heat_id'],case['gas']
        with np.load(PRIOR/'baseline/data'/f'{heat}_q2_processed.npz',allow_pickle=False) as z:
            t=z['time'];q=z['oxygen_flow'];y=z[gas+'_clean'];src,tgt=events_at(z,gas,base,31,9)
            idx=(t>=case['left'])&(t<=case['right']);u=np.load(RESEARCH/'results/real_examples'/f'{heat}_{gas}.npz')
            fig,axes=plt.subplots(3,1,figsize=(7.2,6.0),sharex=True)
            axes[0].plot(t[idx],q[idx],color='#1f557d',lw=1,label='Oxygen flow')
            axes[1].plot(t[idx],y[idx],color='#a8522e',lw=1,label=gas.upper()+' signal')
            for ax,ev,values,label in [(axes[0],src,q,'Source candidate'),(axes[1],tgt,y,'Gas candidate')]:
                pp=[s['peak_index'] for s in ev if case['left']<=s['peak']<=case['right']]
                ax.scatter(t[pp],values[pp],marker='x',s=17,label=label,zorder=3)
            for key,color,label in [('source','#1f557d','Oxygen change'),('response','#a8522e',gas.upper()+' change')]:
                val=u[key];den=max(np.max(abs(val[idx])),1e-12)
                axes[2].plot(t[idx],val[idx]/den,color=color,lw=.9,label=label+' (local scaled)')
            u.close()
            if 'source_peak' in case:
                for ax in axes:
                    ax.axvline(case['source_peak'],ls='--',color='#1f557d',lw=.8)
                    ax.axvline(case['gas_peak'],ls=':',color='#a8522e',lw=.8)
            else:
                for ax in axes:
                    for center,lag in zip(case['centers'],case['lags']):
                        ax.axvline(center,color='gray',alpha=.4,lw=.8)
                        ax.axvline(center+lag,color='#6d4291',alpha=.4,lw=.8,ls=':')
            for ax,title in zip(axes,['(a) Oxygen flow and source candidates','(b) Gas signal and response candidates','(c) Change signals and competing fluctuations']):
                ax.set_title(title,loc='center');ax.grid(alpha=.15);ax.legend(loc='upper right',fontsize=7)
            axes[0].set_ylabel('Flow (m3/s)');axes[1].set_ylabel('Gas signal (m3/s)');axes[2].set_ylabel('Scaled derivative');axes[-1].set_xlabel('Time (s)')
            axes[-1].set_xlim(case['left'],case['right'])
            save(fig,case['case_id'],f"{heat} {gas.upper()}：{case['kind']}。叉号表示检测候选，不代表已确认工艺事件；虚线仅定位已有方法的源峰或窗口中心，点线定位目标峰或中心加时滞。底图各自按显示区间绝对最大值缩放，不比较幅值。")
    fig,ax=plt.subplots(1,2,figsize=(7.2,3.2));x=np.arange(len(cadence))
    ax[0].bar(x,[float(r['equal_adjacent_fraction']) for r in cadence],color='#1f557d');ax[0].set_ylabel('Equal adjacent fraction');ax[0].set_title('(a) Raw cumulative oxygen updates')
    ax[1].bar(x,[float(r['positive_increment_band_fraction']) for r in cadence],color='#a8522e');ax[1].set_ylabel('Fraction of positive increments');ax[1].set_title('(b) Increments within 4.5 to 5.0 m3')
    for a in ax:a.set_xticks(x,[r['heat_id'] for r in cadence],rotation=60);a.set_ylim(0,1);a.grid(axis='y',alpha=.15)
    save(fig,'sampling_cadence','原始累计耗氧量的重复采样与近固定增量统计。该现象提示计数更新结构，但不足以证明仪器机制或认定全部候选为伪事件。')
    fig,ax=plt.subplots(1,2,figsize=(7.2,3.2))
    for a,gas in zip(ax,['co','co2']):
        for qw in cfg['flow_sg_windows']:
            a.plot(cfg['change_sg_windows'],[sum(int(r['source_candidates']) for r in scales if r['gas']==gas and int(r['flow_sg_samples'])==qw and int(r['change_sg_samples'])==cw) for cw in cfg['change_sg_windows']],marker='o',label=f'Flow SG {qw}')
        a.set_title('('+('a' if gas=='co' else 'b')+') '+gas.upper()+' source representation');a.set_xlabel('Change SG window (samples)');a.set_ylabel('Source candidates');a.legend();a.grid(alpha=.15)
    save(fig,'scale_sensitivity','九种预设表示尺度下的源候选数量。两个气体共用氧枪输入，左右相同是预期结果，不构成两份独立样本；本轮不按配对数量选择新参数。')
    pdfrecords=[]
    for pdf in pdfs:
        output=subprocess.check_output(['pdffonts',str(pdf)],text=True)
        lines=output.splitlines()[2:]
        if not lines or any('TimesNewRoman' not in line.replace(' ','') or 'yes' not in line for line in lines):raise AssertionError('PDF Times New Roman embedding failed: '+str(pdf))
        # Poppler columns: name type [encoding] emb sub uni object ID; embedding is parsed from the right.
        if any(line.split()[-5]!='yes' for line in lines):raise AssertionError('PDF font not embedded')
        pdfrecords.append(dict(path=str(pdf.relative_to(out)),pdffonts=output))
    write_json(out/'font_validation.json',dict(configuration=font,pdf_checks=pdfrecords))
    write_csv(out/'figure_captions.csv',captions)
    thumbs=[]
    for record in captions:
        im=Image.open(figures/(record['figure']+'.png')).convert('RGB');im.thumbnail((355,295))
        tile=Image.new('RGB',(365,320),'white');tile.paste(im,((365-im.width)//2,0));ImageDraw.Draw(tile).text((5,300),record['figure'],fill='black');thumbs.append(tile)
    preview=Image.new('RGB',(365*4,320*((len(thumbs)+3)//4)),'#dddddd')
    for i,tile in enumerate(thumbs):preview.paste(tile,((i%4)*365,(i//4)*320))
    preview.save(out/'layout_preview.png')
    # All source rows are searchable locally; no records or arrays leave the workstation.
    headers=['heat_id','gas','event_id','source_peak','source_flow_m3_s','oxygen_stage','status','other_sources_within_60s','scale_agreements','target_ids']
    table=''.join('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in headers)+'</tr>' for r in sources)
    links=''.join('<li><a href="figures/'+r['figure']+'.png">'+html.escape(r['figure'])+'</a> '+html.escape(r['caption'])+'</li>' for r in captions)
    page='''<!doctype html><meta charset="utf-8"><title>十炉真实事件审查</title><style>body{font:15px sans-serif;margin:24px;max-width:1500px}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:5px}input{padding:8px;width:450px}th{position:sticky;top:0;background:white}</style><h1>十炉真实事件审查</h1><p>所有行都是候选。未确认物理匹配或独立性；CO 和 CO₂ 使用同一氧枪源候选，不得相加。</p><input id="filter" placeholder="输入炉号、气体、source 编号或状态筛选"><table><thead><tr>'''+''.join('<th>'+k+'</th>' for k in headers)+'</tr></thead><tbody>'+table+'</tbody></table><h2>波形检查页</h2><ul>'+links+r'''</ul><script>document.querySelector('#filter').oninput=e=>{let q=e.target.value.toLowerCase().split(/\s+/);document.querySelectorAll('tbody tr').forEach(r=>r.hidden=!q.every(v=>r.textContent.toLowerCase().includes(v)))}</script>'''
    (out/'review.html').write_text(page)
    totals={gas:{key:sum(int(r[key]) for r in rows if r['gas']==gas) for key in ['source_candidates','gas_candidates','same_polarity_edges','unique_pairs','ambiguous_sources','unmatched_sources','all_scales_recurring_sources']} for gas in ['co','co2']}
    def table_md(headers,data):return '| '+' | '.join(headers)+' |\n|'+ '|'.join(['---']*len(headers))+'|\n'+''.join('| '+' | '.join(map(str,r))+' |\n' for r in data)
    pairtab=table_md(['炉号','气体','方法','源峰 s','目标峰 s','起始差 s','峰值差 s','邻近其他源数','尺度重现数'],[[r['heat_id'],r['gas'],r['method'],r['source_peak'],r['gas_peak'],r['delay_onset'],r['delay_peak'],r['other_sources_within_60s'],r['scale_agreements']+'/9'] for r in pairs])
    heattab=table_md(['炉号','气体','源候选','目标候选','唯一边','歧义源','无匹配','隔离源'],[[r[k] for k in ['heat_id','gas','source_candidates','gas_candidates','unique_pairs','ambiguous_sources','unmatched_sources','isolated_source_candidates']] for r in rows])
    controltotal={gas:[sum(int(r['unique_pairs']) for r in controls if r['gas']==gas and int(r['replicate'])==k) for k in range(cfg['circular_controls'])] for gas in ['co','co2']}
    controltab=table_md(['气体','观测唯一边','循环移位总数中位数','移位总数范围'],[[gas,totals[gas]['unique_pairs'],float(np.median(controltotal[gas])),f"{min(controltotal[gas])}–{max(controltotal[gas])}"] for gas in ['co','co2']])
    notes=json.loads((Path(__file__).parent/'review_notes.json').read_text())
    write_json(out/'case_review_notes.json',notes)
    notetab=table_md(['炉号与气体','波形复核观察','当前判定'],[[r['heat_id']+' '+r['gas'],r['observation'],r['assessment']] for r in notes['cases']])
    text=f'''# 十炉真实数据事件级验证

以 `{cfg['baseline_commit']}` 为基准，本轮重放十炉真实候选和24个可靠窗口，检查源事件竞争、表示尺度、采样更新与机会配对。当前证据仍不足以构造具有统计依据的流量—时滞函数。没有独立响应标签，不能给真实数据计算Precision、Recall或F1，也不能把“唯一候选边”当作已确认物理响应。

## 数据身份和核验范围

十炉原始XLSB及基线NPZ逐文件SHA256与已冻结清单一致，见 [input_manifest.json](input_manifest.json)。原始工艺列没有被第一问光谱接口替代。Release历史处理结果570字段的一致性已在上轮核验；用户已明确确认中文附件1只是本地英文文件的改名，不存在另一个待取得的数据版本。文件哈希与既有清单一致，附件1身份核验项关闭，见 [attachment1_identity.json](attachment1_identity.json)。此前报告中的待核验为当时状态，由本轮来源确认取代。6B03205现有工艺三列缺失计数均为零；没有证据表明存在需要补回的另一版缺失数据。

## 事件和窗口重放

共重现959个共享氧枪源候选，CO和CO₂各引用这959个源，不得合计为1918次工艺事件。原24个窗口全部保留位置和时滞，逐窗候选数见 [baseline_24_window_review.csv](baseline_24_window_review.csv)。CO源窗口重叠组7个、CO₂为5个；对应12个检查页仍称窗口组，不称独立响应。另有{len(pairs)}个已有方法的配对检查页，共{len(cases)}个波形检查页。

{heattab}
已有V2事件适配器在全记录中有3条候选配对（CO 1条、CO₂ 2条），直接事件法1条；均为下降方向的局部变化配对，不能当成主要气体上升响应。之前的事件汇总按CSV直接复核如下，不能据此宣称真实准确率改善。

{pairtab}
起始差、峰值差和重心差是可见变化信号的描述性时间差；起始仍为绝对变化质量5%定义。负起始差不能用正峰值差掩盖，也不能直接解释为气体真正先于氧枪反应。V2自身只估计局部平移量，不估计物理响应起点。上述配对的描述符来自公共候选检测器，而非V2的新输出。

{notetab}
这些是对指定波形页的定性复核，不是独立工艺标注。6C02057的7.5秒唯一边匹配的是小幅下降变化，不能将其移用于图中主要CO上升段。当前2/9尺度重现及邻近源竞争也不足以确认这条局部下降的物理时滞。

## 重叠和歧义

在预设60秒响应范围内，每个源峰的前后60秒邻域内都存在其他检测候选；满足前后各60秒无其他候选的数量为0。这个条件只检验在该时间范围下能否隔离，并不证明所有候选都是同一工艺事件。源支持区间及其后60秒连成的组见 [source_candidates.csv](source_candidates.csv)，组编号不是独立样本编号。

候选图逐边保留0–60秒内同极性的所有目标。源和目标度数均为1才称“唯一候选边”；多源竞争、多个目标均保留歧义。CO的歧义源为{totals['co']['ambiguous_sources']}个，CO₂为{totals['co2']['ambiguous_sources']}个。另存不限制极性的候选数，提醒“同极性”是匹配假设，不能视为化学机制。全部边的起始差、峰值差、重心差和响应宽度在 [candidate_edges.csv](candidate_edges.csv)；歧义边不能进入回归。

## 计数更新和信号表示

十炉原始累计耗氧量相邻相等比例为{min(float(r['equal_adjacent_fraction']) for r in cadence):.1%}–{max(float(r['equal_adjacent_fraction']) for r in cadence):.1%}，正增量中位数为{min(float(r['positive_increment_median_m3']) for r in cadence):.2f}–{max(float(r['positive_increment_median_m3']) for r in cadence):.2f} m³。原始更新间隔和4.5–5.0 m³增量比例逐炉见 [sampling_cadence.csv](sampling_cadence.csv)。这些现象与离散累计计数更新相符，但采集机制尚未独立确认；不能直接认定某个事件是仪器伪影。

采用21、31、41点流量SG与9、21、41点变化SG的九种组合，只做事后敏感性审查。原31/9配置精确重放既有候选；其他尺度的同极性峰在±5秒内用一对一分配比较时刻重现，不是检测准确率。CO共享源中{totals['co']['all_scales_recurring_sources']}个在全部九配置中重现。稳定的时刻也不证明是独立物理事件。不依据匹配条数选择平滑参数，不替换正式31点流量。

## 循环移位机会配对对照

固定源候选，将目标候选整体循环平移，每条炉次气体序列使用99个固定随机移位。位移避开±60秒；跨接缝的完整目标波瓣丢弃，实际保留数逐次保存。复用同一同极性0–60秒候选图，不重算V2或核模型。因此这里只评价直接事件规则的机会配对数量。

{controltab}
这项对照保留目标候选间隔结构，不保留非平稳工艺阶段对齐，不能解释为误匹配率、显著性或因果检验。唯一边少本身既不能证明准确，也不能证明无效。结果与全部候选竞争一起表明：目前不能凭唯一配对数量确认真实响应。

## 对流量与时滞问题的回答

目前没有已确认且可作为独立样本的源响应配对，四类流量时滞模型的有效事件样本量尚不能建立。现有24个窗口既重叠又未确认响应身份；直接用它们拟合固定、线性、反平方根或单调曲线，会将窗口重复和匹配不确定性误作统计证据。

因此，本轮对“氧枪流量增加，CO和CO₂时滞是否减小”的回答是：现有证据不足以判定；不能据此肯定，也不能声称已否定该物理关系。固定模型与反平方根的既有训练和留一炉次结果仍由前轮报告提供，保持探索性地位。正式接口继续使用冻结实现，不能把研究中的空值替换成零或导入新时滞。

## 定稿路径和停止条件

附件1身份已确认。下一步结合采集说明判定累计计数更新是否会造成短周期流量波动，再利用这些波形页选择有外部依据的供氧启停或调整段，由工艺记录或独立复核提供源与气体响应标注；不要求每炉都有可辨识时滞。若仍无法形成可信事件样本，第二问应如实报告时滞可辨识性和证据边界，并完成四类耗氧关系的趋势、差分与阶段统计，无需继续增加合成算法。

四类关系统计、最终中文章节和最终接口冻结属于后续定稿工作，本轮未将它们宣称为完成。流量—时滞函数只有在事件证据足够时才比较；炉内和炉际效应需分开，不能用总体相关替代。

## 复现和检查入口

运行命令见 [复现说明](../README.md)。[本地检索页](review.html)覆盖所有源候选和全部波形页；[图注](figure_captions.csv)可编辑，图内为英文与真实Times New Roman，PNG、SVG、PDF逐图导出，PDF嵌入记录见 [font_validation.json](font_validation.json)。所有合成实验保持冻结，第一问、正式V2和第三问NPZ字段均未修改。配置、原始输入与基线NPZ指纹、代码指纹和输出SHA256随结果保存。
'''
    (out/'REPORT.md').write_text(text)
