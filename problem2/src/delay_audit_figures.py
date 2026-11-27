"""English audit figures with editable external captions and embedded TNR PDFs."""
from __future__ import annotations
import csv
import json
import subprocess
import numpy as np
from .figure_fonts import configure
from .delay_audit_core import scan
from .delay_audit import settings,write_json,write_csv
from .q2_pipeline import _standardized_change


def render(out,cfg):
    font_info=configure()
    import matplotlib.pyplot as plt
    from PIL import Image,ImageDraw,ImageFont
    figures=out/'figures';figures.mkdir(exist_ok=True)
    accepted=list(csv.DictReader((out/'accepted_windows.csv').open()))
    expected={f'window_{r["heat_id"]}_{r["gas"]}_{int(r["center_index"]):04d}' for r in accepted}|{'accepted_window_distribution','synthetic_comparison'}
    for path in figures.iterdir():
        if path.suffix in ['.png','.svg','.pdf'] and (path.stem.startswith('window_') or path.stem in expected) and path.stem not in expected:
            path.unlink()
    captions=[];font_checks=[];pngs=[]
    def save(fig,name,caption):
        fig.set_size_inches(7.2,fig.get_figheight())
        for ext in ['png','svg','pdf']:
            meta={'CreationDate':None,'ModDate':None} if ext=='pdf' else ({'Date':None} if ext=='svg' else {})
            fig.savefig(figures/(name+'.'+ext),dpi=200,metadata=meta)
        fonts=subprocess.check_output(['pdffonts',str(figures/(name+'.pdf'))],text=True)
        lines=fonts.splitlines()[2:]
        if not lines or any('TimesNewRoman' not in l.replace('-','') or 'yes' not in l for l in lines):
            raise RuntimeError('unexpected PDF fonts: '+fonts)
        # pdffonts fixed columns: font name, type, encoding, emb/sub/uni, object ID.
        emb_column=fonts.splitlines()[0].index('emb')
        if any(l[emb_column:emb_column+3]!='yes' for l in lines): raise RuntimeError('unembedded font: '+fonts)
        font_checks.append(dict(figure=name,pdffonts=fonts))
        captions.append(dict(figure=name,caption=caption,width_inches=7.2))
        pngs.append(figures/(name+'.png'));plt.close(fig)
    for row in accepted:
        heat,gas=row['heat_id'],row['gas'];c=int(row['center_index'])
        with np.load(out/'baseline/data'/f'{heat}_q2_processed.npz',allow_pickle=False) as z:
            t=z['time'];q=z['oxygen_flow'];y=z[gas+'_clean'];raw=z[gas+'_raw']
        s=scan(q,y,**settings(cfg));i=int(np.flatnonzero(s['center_index']==c)[0]);lag=float(row['best_delay_s'])
        left=int(row['search_left_index']);right=int(row['search_right_index']);half=cfg['delay']['window']//2
        ix=np.arange(left,right+1); si=np.arange(c-half,c+half+1); ti=si+int(round(lag/cfg['time_step_s']))
        fig,axes=plt.subplots(2,2,figsize=(7.2,5.5),layout='constrained')
        ax=axes[0,0];ax.plot(t[ix],q[ix],color='#245a81',lw=1);ax.axvspan(t[si[0]],t[si[-1]],alpha=.13,color='#245a81');ax.axvline(t[c],color='black',ls=':',lw=.8)
        ax.set(title='(a) Oxygen-lance flow',xlabel='Source time (s)',ylabel='Flow (m3/s)')
        ax=axes[0,1];ax.plot(t[ix],raw[ix],color='.7',lw=.7,label='Raw');ax.plot(t[ix],y[ix],color='#b34c32',lw=1,label='Clean');ax.axvspan(t[ti[0]],t[ti[-1]],alpha=.13,color='#b34c32');ax.legend(loc='best');ax.set(title='(b) Candidate gas response',xlabel='Measurement time (s)',ylabel=gas.upper()+' flow (m3/s)')
        ax=axes[1,0];a=_standardized_change(q)[si];b=_standardized_change(y)[ti]
        a=(a-a.mean())/(a.std()+1e-12);b=(b-b.mean())/(b.std()+1e-12)
        ax.plot(t[si],a,lw=1,label='Source change');ax.plot(t[si],b,lw=1,label='Shifted gas change');ax.legend(loc='best');ax.set(title='(c) Changes after candidate alignment',xlabel='Source time (s)',ylabel='Standardized change')
        ax=axes[1,1];ax.plot(s['lags_s'],s['correlations'][i],lw=1,color='#245a81');ax.axhline(cfg['delay']['min_corr'],ls=':',color='.5',lw=.8);ax.axvline(lag,color='#b34c32',ls='--',lw=.8)
        ax.set(title=f'(d) Correlation peak at {lag:g} s',xlabel='Target delay (s)',ylabel='Correlation',ylim=(-.65,1.05));ax.text(.03,.96,'Event identity: uncertain',transform=ax.transAxes,va='top',fontsize=8)
        for ax in axes.flat: ax.grid(alpha=.15);ax.title.set_horizontalalignment('center')
        name=f'window_{heat}_{gas}_{c:04d}'
        save(fig,name,f'{heat}, {gas.upper()}, source center {t[c]:g} s; candidate delay {lag:g} s. Shading marks the source window and shifted target window. This correlation observation does not identify a unique physical event. Group IDs and cross-gas overlap are in accepted_windows.csv.')
    summary=json.loads((out/'real_summary.json').read_text());heats=sorted({r['heat_id'] for r in summary})
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.4),layout='constrained',sharex=True)
    for ax,gas in zip(axes,['co','co2']):
        for row in accepted:
            if row['gas']!=gas:continue
            idx=heats.index(row['heat_id']);ax.scatter(float(row['center_time_s']),idx,c=float(row['best_delay_s']),vmin=0,vmax=60,cmap='viridis',s=28)
        ax.set_yticks(range(len(heats)),heats);ax.set(title=gas.upper()+' accepted windows',xlabel='Source time (s)',ylabel='Heat');ax.grid(alpha=.15);ax.set_xlim(0,800)
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,60),cmap='viridis'),ax=list(axes),label='Delay (s)',shrink=.8)
    save(fig,'accepted_window_distribution','Locations of all 24 baseline accepted windows. Color encodes delay from 0 to 60 s. Points are overlapping window observations, not independent events.')
    aggregate=json.loads((out/'synthetic_summary.json').read_text());names={r['scenario'] for r in aggregate}
    cases=sorted([n for n in names if n.startswith('fixed_')],key=lambda n:float(n[6:-1]))+sorted(n for n in names if not n.startswith('fixed_'))
    def label(name):
        if name.startswith('fixed_'): return 'Delay '+name[6:-1]+' s'
        if name.startswith('broadened_theta_'): return 'Broadening '+name[16:-1]+' s'
        return {'adjacent_events':'Adjacent events','drifting_5_to_25s':'Drift 5-25 s','noise_0.5':'Noise 0.5','noise_1':'Noise 1.0','unrelated':'Unrelated'}[name]
    fig,axes=plt.subplots(2,1,figsize=(7.2,6.0),layout='constrained',sharex=True)
    for mode,color,offset in [('nonnegative','#245a81',-.18),('signed','#b34c32',.18)]:
        rs={r['scenario']:r for r in aggregate if r['mode']==mode and r['truth_definition']=='peak'}
        axes[0].bar(np.arange(len(cases))+offset,[rs[k]['acceptance_rate'] for k in cases],width=.35,label=mode,color=color)
        axes[1].scatter(np.arange(len(cases))+offset,[rs[k].get('mae_s') if rs[k].get('mae_s') is not None else np.nan for k in cases],label=mode,color=color,s=20)
    axes[0].set(title='(a) Acceptance on shared source centers',ylabel='Accepted fraction',ylim=(0,1.1));axes[0].legend()
    axes[1].set(title='(b) Error conditional on acceptance',ylabel='MAE to response peak (s)');axes[1].set_xticks(range(len(cases)),[label(n) for n in cases],rotation=55,ha='right',fontsize=7)
    for ax in axes:ax.grid(axis='y',alpha=.15)
    save(fig,'synthetic_comparison','Thirty fixed seeds per scenario. Nonnegative and signed searches use identical source centers. Error is conditional on acceptance; the null pair has no true delay. Broadened-response errors refer to the known impulse-response peak, not its onset or centroid. Adjacent-event metrics score centers within 5 s of a known source event.')
    # A single contact sheet for layout inspection, kept separate from formal figures.
    thumbs=[]
    font=ImageFont.truetype(font_info['fonts'][0]['file'],18)
    for path in pngs:
        im=Image.open(path).convert('RGB');im.thumbnail((550,450))
        tile=Image.new('RGB',(570,490),'white');tile.paste(im,((570-im.width)//2,30));ImageDraw.Draw(tile).text((10,5),path.stem,font=font,fill='black');thumbs.append(tile)
    preview=Image.new('RGB',(570*3,490*((len(thumbs)+2)//3)),'#dddddd')
    for i,tile in enumerate(thumbs):preview.paste(tile,((i%3)*570,(i//3)*490))
    preview.save(out/'layout_preview.png')
    write_csv(out/'figure_captions.csv',captions)
    write_json(out/'font_validation.json',dict(configuration=font_info,pdf_checks=font_checks))
