"""Current Q2 figures with strict Times New Roman and external captions."""
from pathlib import Path
import json,subprocess
import numpy as np
from PIL import Image,ImageDraw
from .figure_fonts import configure
from .delay_audit import write_csv,write_json
from .final_statistics import RELATIONS,target,progress
from .q2_pipeline import oxygen_flow

ROOT=Path(__file__).resolve().parents[2]


def render(heats,tables,cfg):
    font=configure()
    import matplotlib.pyplot as plt
    folder=ROOT/cfg['figure_dir'];folder.mkdir(parents=True,exist_ok=True);out=ROOT/cfg['output_dir'];captions=[];checks=[]
    colors=plt.get_cmap('tab10').colors
    def finish(fig,name,caption,role='core',global_legend=False):
        if global_legend:
            handles,labels=fig.axes[0].get_legend_handles_labels();fig.legend(handles,labels,ncol=5,loc='upper center',fontsize=7,frameon=False)
        for ax in fig.axes:ax.grid(alpha=.16)
        fig.tight_layout(rect=(0,0,1,.94) if global_legend else (0,0,1,1))
        for ext in ['png','svg','pdf']:fig.savefig(folder/(name+'.'+ext),dpi=220,metadata={'CreationDate':None,'ModDate':None,'Creator':'APMCM Q2 final analysis'} if ext=='pdf' else None)
        plt.close(fig)
        pdf=subprocess.check_output(['pdffonts',str(folder/(name+'.pdf'))],text=True);lines=pdf.splitlines()[2:];emb=pdf.splitlines()[0].index('emb')
        if not lines or any('TimesNewRoman' not in line.replace(' ','') or line[emb:emb+3]!='yes' for line in lines):raise RuntimeError('Times New Roman PDF embedding failed '+name)
        checks.append(dict(figure=name,pdffonts=pdf));captions.append(dict(figure=name,role=role,width_inches=7.2,caption=caption))
    labels={'co_instantaneous':'CO rate (m3/s)','co2_instantaneous':'CO2 rate (m3/s)','co_cumulative':'Cumulative CO (m3)','co2_cumulative':'Cumulative CO2 (m3)'}
    for cumulative,name in [(False,'oxygen_instantaneous_relationships'),(True,'oxygen_cumulative_relationships')]:
        rr=RELATIONS[2:] if cumulative else RELATIONS[:2]
        fig,axes=plt.subplots(1,2,figsize=(7.2,3.6))
        for j,(ax,relation) in enumerate(zip(axes,rr)):
            for i,d in enumerate(heats):ax.plot(progress(d),target(d,relation),color=colors[i],lw=.8,alpha=.8,label=d['heat_id'])
            coef=next(r for r in tables['pooled_trends'] if r['relation']==relation and r['degree']==2)
            x=np.linspace(0,1,201);pred=coef['c0']+coef['c1']*x+coef['c2']*x*x
            ax.plot(x,pred,color='black',ls='--',lw=1.3)
            ax.set_title(f'({chr(97+j)}) '+('CO2' if relation.startswith('co2') else 'CO'));ax.set_xlabel('Oxygen progress');ax.set_ylabel(labels[relation])
        finish(fig,name,'累计耗氧进度与'+('累计排放' if cumulative else '瞬时排放')+'的炉次曲线。颜色区分十炉，黑色虚线为每炉等权的二次描述性趋势，不作反应动力学或因果解释。','core',True)
    fig,axes=plt.subplots(2,2,figsize=(7.2,6.0));xs=np.arange(1,6)
    for col,gas in enumerate(['co','co2']):
        rows=[r for r in tables['stage_aggregate'] if r['gas']==gas]
        for row,key,label in [(0,'mean_rate_m3_s','Mean rate (m3/s)'),(1,'emission_per_oxygen_m3_m3','Stage gas / oxygen (m3/m3)')]:
            ax=axes[row,col]
            for d in heats:
                vals=[r[key] for r in tables['stage_statistics'] if r['heat_id']==d['heat_id'] and r['gas']==gas]
                ax.plot(xs,vals,color='gray',lw=.55,alpha=.35)
            ax.fill_between(xs,[r[key+'_ci_low'] for r in rows],[r[key+'_ci_high'] for r in rows],color=colors[col],alpha=.18,label='95% heat bootstrap interval')
            ax.plot(xs,[r[key+'_median'] for r in rows],color=colors[col],marker='o',label='Across-heat median')
            ax.set_title(f'({chr(97+row*2+col)}) '+gas.upper()+(' rate' if row==0 else ' stage ratio'));ax.set_ylabel(label);ax.set_xticks(xs);ax.set_xlabel('Oxygen stage (equal 20% progress)')
        axes[0,col].legend(fontsize=6.5)
    finish(fig,'oxygen_stage_statistics','五个等耗氧进度阶段的排放变化。细灰线为逐炉结果，粗线为十炉中位数，区间按整炉重采样2000次得到；不把时间采样点视为独立样本。阶段排放除以阶段耗氧量为表观体积比，不是化学收率。')
    fig,axes=plt.subplots(2,2,figsize=(7.2,5.8))
    for j,(ax,relation) in enumerate(zip(axes.flat,RELATIONS)):
        rr=[r for r in tables['level_correlations'] if r['relation']==relation]
        for metric,offset,color,marker in [('pearson_r',-.13,colors[0],'o'),('spearman_r',.13,colors[1],'x')]:
            groups=[rr]+[[r for r in tables['increment_correlations'] if r['relation']==relation and r['interval_s']==seconds] for seconds in cfg['increment_intervals_s']]
            for k,gg in enumerate(groups):
                vals=[r[metric] for r in gg if r[metric] is not None];ax.scatter(np.full(len(vals),k+offset)+np.linspace(-.04,.04,len(vals)),vals,color=color,s=13,marker=marker,label=('Pearson' if metric=='pearson_r' else 'Spearman') if k==0 else None)
                if vals:ax.plot([k+offset-.06,k+offset+.06],[np.median(vals)]*2,color='black',lw=1.4)
        ax.axhline(0,color='gray',lw=.7);ax.set_ylim(-1.05,1.05);ax.set_xticks(range(4),['Level','Delta 0.5 s','Delta 5 s','Delta 15 s'],rotation=15)
        ax.set_title(f'({chr(97+j)}) '+labels[relation]);ax.set_ylabel('Within-heat correlation')
    axes[0,0].legend(fontsize=7)
    finish(fig,'oxygen_trend_diagnostics','四类关系的原值相关与不重叠0.5、5、15秒增量相关。每个标记是一炉的描述性相关，黑线为炉间中位数。差分用于检查共同趋势与时间尺度，不是因果检验，不报告基于采样点独立性的p值。')
    d=heats[0];t=d['time'];idx=t<=100
    fig,axes=plt.subplots(2,2,figsize=(7.2,5.8))
    axes[0,0].step(t[idx],d['oxygen_cumulative_raw'][idx]-d['oxygen_cumulative_raw'][0],where='post',lw=1);axes[0,0].set_title('(a) Raw cumulative oxygen: '+d['heat_id']);axes[0,0].set_ylabel('Oxygen change (m3)');axes[0,0].set_xlabel('Time (s)')
    for window in cfg['flow_sg_windows']:axes[0,1].plot(t[idx],oxygen_flow(t,d['oxygen_cumulative_clean'],window)[idx],lw=.8,label=f'SG {window}')
    axes[0,1].set_title('(b) Derivative sensitivity');axes[0,1].set_xlabel('Time (s)');axes[0,1].set_ylabel('Oxygen flow (m3/s)');axes[0,1].legend(fontsize=7)
    for w,marker in [(21,'o'),(31,'x'),(41,'s')]:
        rr=[r for r in tables['flow_review'] if r['sg_window']==w]
        axes[1,0].plot(range(10),[r['reintegration_rmse_m3'] for r in rr],marker=marker,lw=.8,label=f'SG {w}')
        axes[1,1].plot(range(10),[100*r['endpoint_relative_error'] for r in rr],marker=marker,lw=.8,label=f'SG {w}')
    for ax in axes[1]:ax.set_xticks(range(10),[d['heat_id'] for d in heats],rotation=60);ax.legend(fontsize=7)
    axes[1,0].set_title('(c) Reintegration consistency');axes[1,0].set_ylabel('Cumulative RMSE (m3)');axes[1,1].set_title('(d) End-point volume consistency');axes[1,1].set_ylabel('Relative error (%)')
    finish(fig,'oxygen_flow_review','原始累计计数及21、31、41点SG导数的数值敏感性。31点为冻结正式流量，其他尺度仅诊断；反积分接近累计量不代表瞬时流量已被独立验证。典型区间固定为6B03205前100秒。')
    fig,axes=plt.subplots(3,2,figsize=(7.2,8.0))
    keys=['oxygen_cumulative_clean','oxygen_flow','co_clean','co2_clean','total_co_emission','total_co2_emission'];names=['Cumulative oxygen (m3)','Oxygen flow (m3/s)','CO rate (m3/s)','CO2 rate (m3/s)','Cumulative CO (m3)','Cumulative CO2 (m3)']
    for j,(ax,key,label) in enumerate(zip(axes.flat,keys,names)):
        for i,d in enumerate(heats):ax.plot(d['time'],d[key],color=colors[i],lw=.65,label=d['heat_id'])
        ax.set_title(f'({chr(97+j)}) '+label);ax.set_ylabel(label);ax.set_xlabel('Time (s)')
    finish(fig,'ten_heats_six_process_curves','十炉六类工艺信号随相对时间的变化。烟气字段按题面单位解释为体积流量，累计排放由其时间积分得到。各炉记录长度不同，不截断或外推补齐。','core',True)
    for d in heats:
        fig,axes=plt.subplots(2,2,figsize=(7.2,5.6))
        for j,(ax,relation) in enumerate(zip(axes.flat,RELATIONS)):
            p=progress(d);ax.plot(d['oxygen_cumulative_clean'],target(d,relation),lw=.8,color=colors[0],label='Clean signal')
            for degree,style in [(1,':'),(2,'--')]:
                r=next(r for r in tables['trend_parameters'] if r['heat_id']==d['heat_id'] and r['relation']==relation and r['degree']==degree)
                ax.plot(d['oxygen_cumulative_clean'],r['c0']+r['c1']*p+r['c2']*p*p,ls=style,lw=.9,label='Linear' if degree==1 else 'Quadratic')
            ax.set_title(f'({chr(97+j)}) '+labels[relation]);ax.set_xlabel('Cumulative oxygen (m3)');ax.set_ylabel(labels[relation]);ax.legend(fontsize=6.5)
        finish(fig,d['heat_id']+'_oxygen_relationships',d['heat_id']+'的四类关系及线性、二次描述性趋势。仅在本炉已观测耗氧区间内解释，不将回归外推或拟合优度当作物理规律。','supplement')
    d=heats[0];fig,axes=plt.subplots(3,1,figsize=(7.2,6.0),sharex=True)
    axes[0].plot(d['time'],d['oxygen_cumulative_raw'],color='gray',lw=.8,label='Raw oxygen');axes[0].plot(d['time'],d['oxygen_cumulative_clean'],lw=.9,label='Clean oxygen');axes[0].set_ylabel('Oxygen (m3)');axes[0].set_title('(a) Raw and reconstructed cumulative oxygen')
    for gas in ['co','co2']:axes[1].plot(d['time'],d[gas+'_clean'],lw=.9,label=gas.upper());axes[2].scatter(d[f'delay_{gas}_window_time'],d[f'delay_{gas}_observed'],s=12,label=gas.upper())
    axes[1].set_title('(b) Reconstructed gas signals');axes[1].set_ylabel('Gas rate (m3/s)');axes[2].set_title('(c) Accepted statistical delay windows');axes[2].set_ylabel('Delay (s)');axes[2].set_xlabel('Time (s)')
    if not sum(d[f'delay_{gas}_valid_mask'].sum() for gas in ['co','co2']):axes[2].text(.5,.5,'No accepted delay windows',ha='center',va='center',transform=axes[2].transAxes)
    for ax in axes:ax.legend(fontsize=7)
    finish(fig,d['heat_id']+'_q2_processed_signals_and_delay',d['heat_id']+'重构及统计时滞诊断。本炉没有通过正式V2判据的窗口；空白不填零，不绘制未经支持的时滞曲线。','supplement')
    old=json.loads((ROOT/'problem2/final/results/cross_validation.json').read_text());fig,ax=plt.subplots(figsize=(7.2,3.5))
    for gas in ['co','co2']:
        ts=sorted({r['min_corr'] for r in old['threshold_sensitivity'] if r['gas']==gas});values=[]
        for threshold in ts:
            rr=[r for r in old['threshold_sensitivity'] if r['gas']==gas and r['min_corr']==threshold];values.append(sum(r['accepted_windows'] for r in rr)/sum(r['candidate_windows'] for r in rr))
        ax.plot(ts,values,marker='o',label=gas.upper())
    ax.set_title('Accepted windows under diagnostic thresholds');ax.set_xlabel('Minimum correlation threshold');ax.set_ylabel('Accepted-window fraction');ax.legend()
    finish(fig,'delay_threshold_sensitivity','既有阈值敏感性结果的统一排版。正式阈值仍为0.35，不依据覆盖率调整判定或增加物理事件数。','supplement')
    write_csv(out/'figure_captions.csv',captions);write_json(out/'font_validation.json',dict(configuration=font,pdf_checks=checks))
    tiles=[]
    for row in captions:
        image=Image.open(folder/(row['figure']+'.png')).convert('RGB');image.thumbnail((345,300));tile=Image.new('RGB',(355,325),'white');tile.paste(image,((355-image.width)//2,0));ImageDraw.Draw(tile).text((5,305),row['figure'],fill='black');tiles.append(tile)
    preview=Image.new('RGB',(355*4,325*((len(tiles)+3)//4)),'#ddd')
    for i,tile in enumerate(tiles):preview.paste(tile,((i%4)*355,(i//4)*325))
    preview.save(out/'layout_preview.png')
