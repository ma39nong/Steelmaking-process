"""Research plots, strict Times New Roman and external editable captions."""
from pathlib import Path
import csv,json,subprocess
import numpy as np
from ...src.figure_fonts import configure
from ...src.delay_audit import write_csv,write_json


def render(out,cfg):
    fonts=configure()
    import matplotlib.pyplot as plt
    from PIL import Image,ImageDraw,ImageFont
    folder=out/'figures';folder.mkdir(exist_ok=True);captions=[];checks=[];pngs=[]
    def save(fig,name,caption):
        fig.set_size_inches(7.2,fig.get_figheight())
        for ax in fig.axes:
            ax.grid(alpha=.15);ax.title.set_horizontalalignment('center')
        for ext in ['png','svg','pdf']:
            meta={'CreationDate':None,'ModDate':None} if ext=='pdf' else ({'Date':None} if ext=='svg' else {})
            fig.savefig(folder/(name+'.'+ext),dpi=200,metadata=meta)
        text=subprocess.check_output(['pdffonts',str(folder/(name+'.pdf'))],text=True);header=text.splitlines()[0];lines=text.splitlines()[2:];emb=header.index('emb')
        if not lines or any('TimesNewRoman' not in line.replace('-','') or line[emb:emb+3]!='yes' for line in lines):raise RuntimeError('PDF font embedding check failed: '+name)
        checks.append(dict(figure=name,pdffonts=text));captions.append(dict(figure=name,width_inches=7.2,caption=caption));pngs.append(folder/(name+'.png'));plt.close(fig)
    model_list=json.loads((out/'real/kernel_validation.json').read_text())
    expected={'synthetic_'+n for n in cfg['synthetic']['scenarios']}|{'real_'+r['heat_id'] for r in model_list}|{'event_comparison','kernel_diagnostics','background_sensitivity'}
    for path in folder.iterdir():
        if path.suffix in ['.png','.svg','.pdf'] and path.stem not in expected:path.unlink()
    eventrows=list(csv.DictReader((out/'synthetic/event_matches.csv').open()));seed=cfg['synthetic']['evaluation_seeds'][0]
    colors={'v2':'#245a81','events':'#c78d24','kernel':'#3d805e'}
    for name in cfg['synthetic']['scenarios']:
        with np.load(out/'examples'/(name+'.npz')) as z:d={k:z[k] for k in z.files}
        t=d['time'];start=int(len(t)*cfg['kernel']['validation_end_fraction'])+int(cfg['max_response_s']/cfg['dt_s'])+cfg['change_smoothing_samples']//2
        fig,axes=plt.subplots(2,2,figsize=(7.2,5.5),layout='constrained')
        axes[0,0].plot(t,d['source'],lw=.9,color='#245a81');axes[0,0].set(title='(a) Input change signal',xlabel='Time (s)',ylabel='Change rate')
        axes[0,1].plot(t[start:],d['response'][start:],lw=.9,label='Response');axes[0,1].plot(t[start:],d['prediction'][start:],lw=.9,label='FIR prediction');axes[0,1].legend();axes[0,1].set(title='(b) Held-out response',xlabel='Time (s)',ylabel='Change rate')
        tt=np.arange(len(d['estimated_kernel']))*cfg['dt_s'];axes[1,0].plot(tt,d['estimated_kernel'],label='Estimated',lw=1)
        if len(d['true_kernel']):axes[1,0].plot(tt,d['true_kernel'],ls='--',label='Known kernel',lw=1)
        else:axes[1,0].text(.05,.95,'No stationary kernel truth',transform=axes[1,0].transAxes,va='top',fontsize=8)
        axes[1,0].legend(loc='best');axes[1,0].set(title='(c) Response kernel',xlabel='Lag (s)',ylabel='Discrete response weight')
        ax=axes[1,1]
        for j,method in enumerate(['v2','events','kernel']):
            rs=[r for r in eventrows if int(r['seed'])==seed and r['scenario']==name and r['method']==method and r['evaluation_partition']=='test']
            for r in rs:
                marker={'matched':'o','ambiguous':'s','unmatched':'x'}[r['status']]
                ax.scatter(float(r['source_peak']),j,marker=marker,color=colors[method],s=18)
                if r['gas_peak']:ax.plot([float(r['source_peak']),float(r['gas_peak'])],[j,j],color=colors[method],lw=.8)
        ax.set_yticks(range(3),['V2 adapter','Events','Kernel screen']);ax.set(title='(d) Source-response candidates',xlabel='Source / response time (s)',xlim=(t[start],t[-1]),ylim=(-.5,2.5))
        save(fig,'synthetic_'+name,f'{name}, evaluation seed {seed}. Parameters are fixed across scenarios. The FIR is trained on the first block and selected using validation only. Timeline markers: circle = matched, square = ambiguous, cross = unmatched; a match decision is not necessarily correct. Kernel onset uses 5% absolute cumulative response mass; width uses 5-95%.')
    models=json.loads((out/'real/kernel_validation.json').read_text());heats=sorted({r['heat_id'] for r in models})
    for heat in heats:
        fig,axes=plt.subplots(2,2,figsize=(7.2,5.5),layout='constrained')
        for j,gas in enumerate(['co','co2']):
            with np.load(out/'real_examples'/f'{heat}_{gas}.npz') as z:d={k:z[k] for k in z.files}
            row=next(r for r in models if r['heat_id']==heat and r['gas']==gas);mask=d['time']>=row['test_start_s']
            if j==0:axes[0,0].plot(d['time'],d['source'],lw=.8,color='#245a81')
            ax=axes[0,1] if j==0 else axes[1,0]
            ax.plot(d['time'][mask],d['response'][mask],lw=.9,label='Observed change');ax.plot(d['time'][mask],d['prediction'][mask],lw=.9,label='FIR prediction');ax.legend(fontsize=7)
            ax.set(title=f'({"b" if j==0 else "c"}) {gas.upper()} held-out response',xlabel='Time (s)',ylabel='Change rate')
            h=d['kernel'];den=max(abs(h).sum(),1e-12);axes[1,1].plot(np.arange(len(h))*cfg['dt_s'],h/den,label=gas.upper(),lw=1)
        axes[0,0].set(title='(a) Oxygen-flow changes',xlabel='Time (s)',ylabel='Change rate')
        axes[1,1].axhline(0,color='.5',lw=.5);axes[1,1].legend();axes[1,1].set(title='(d) Candidate response kernels',xlabel='Lag (s)',ylabel='Weight / absolute mass')
        save(fig,'real_'+heat,f'{heat}: CO and CO2 are modelled separately. Panels b and c are held-out samples after a full-horizon guard. Panel d displays signed, absolute-mass-normalized fitted FIRs, not validated physical transport kernels. No independent event labels are available for this heat.')
    summary=json.loads((out/'synthetic/summary.json').read_text());names=cfg['synthetic']['scenarios'];labels=['Gaussian','Asymmetric','Slow','Adjacent','Noise / trend','Mixed strength','Null','Biphasic']
    fig,axes=plt.subplots(2,1,figsize=(7.2,5.5),layout='constrained')
    for j,method in enumerate(['v2','events','kernel']):
        rs=[next(r for r in summary if r['scenario']==n and r['method']==method) for n in names]
        values=[r['f1'] if r['truth_events'] else np.nan for r in rs]
        axes[0].bar(np.arange(len(names))+(j-1)*.23,values,width=.23,label={'v2':'V2 adapter','events':'Events','kernel':'Kernel screen'}[method],color=colors[method])
        values=[r['false_match_rate'] if r['false_match_rate'] is not None else np.nan for r in rs]
        axes[1].scatter(np.arange(len(names))+(j-1)*.18,values,s=20,color=colors[method])
    axes[0].text(names.index('null_response'),.04,'N/A',ha='center',fontsize=8,fontstyle='italic')
    axes[0].legend(ncol=3);axes[0].set(title='(a) Event F1 on held-out sources',ylabel='Event F1',ylim=(0,1.15));axes[0].set_xticks(range(len(names)),labels,rotation=25,ha='right')
    axes[1].set(title='(b) False matches among accepted pairs',ylabel='False match fraction',ylim=(-.05,1.05));axes[1].set_xticks(range(len(names)),labels,rotation=25,ha='right')
    save(fig,'event_comparison','Thirty evaluation seeds per scenario. F1 is pooled from event TP/FP/FN, not overlapping windows. All methods share source/response event detection; V2 uses an explicit event adapter, and Kernel screen uses an estimated FIR to screen candidate edges. Null cases have no true response events and therefore no recall/F1 comparison. Missing false-match fractions indicate no accepted pairs, not perfect precision.')
    km=json.loads((out/'synthetic/kernel_metrics.json').read_text());known=[n for n in names if any(r['scenario']==n and r['truth'] is not None for r in km)]
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.8),layout='constrained')
    for j,key in enumerate(['onset','peak','centroid']):
        vals=[np.mean([r['errors'][key] for r in km if r['scenario']==n]) for n in known]
        axes[0].bar(np.arange(len(known))+(j-1)*.23,vals,width=.23,label=key)
    axes[0].set_xticks(range(len(known)),[labels[names.index(n)] for n in known],rotation=40,ha='right',fontsize=7);axes[0].legend();axes[0].set(title='(a) Kernel timing error',ylabel='MAE (s)')
    values=[np.mean([r['errors']['width'] for r in km if r['scenario']==n]) for n in known]
    axes[1].bar(range(len(known)),values,color='#3d805e');axes[1].set_xticks(range(len(known)),[labels[names.index(n)] for n in known],rotation=40,ha='right',fontsize=7)
    axes[1].set(title='(b) Kernel width error',ylabel='MAE of 5-95% width (s)')
    save(fig,'kernel_diagnostics','Kernel timing errors use a single known stationary kernel only; adjacent/null cases have none. Widths and effective onsets use cumulative absolute kernel mass, not strict support or FWHM. Errors include every fitted kernel, without selecting successful runs.')
    fig,axs=plt.subplots(2,1,figsize=(7.2,5.5),layout='constrained');ax=axs[0]
    nuisance=list(csv.DictReader((out/'real/nuisance_sensitivity.csv').open()))
    for order,color in [(0,'#245a81'),(1,'#b34c32')]:
        vals=[float(r['test_r2']) for r in nuisance if int(r['nuisance_order'])==order]
        ax.plot(range(len(vals)),vals,'o-',lw=.7,ms=3,label='Constant background' if order==0 else 'Linear background',color=color)
    ax.axhline(0,color='.5',lw=.7);ax.legend();ax.set(title='(a) Held-out sensitivity to background specification',xlabel='Heat-gas series index',ylabel='Test R2')
    axs[1].bar(range(len(models)),[r['sg_peak_range_s'] for r in models],color='#3d805e');axs[1].set(title='(b) Kernel peak sensitivity to SG window',xlabel='Heat-gas series index',ylabel='Peak-delay range (s)')
    save(fig,'background_sensitivity','Post hoc diagnostic of constant versus linear nuisance background within the same FIR model. This comparison was added after observing poor real-data extrapolation, and is not used to choose the reported main estimator or claim held-out validation for a revised model.')
    font=ImageFont.truetype(fonts['fonts'][0]['file'],17);tiles=[]
    for p in pngs:
        im=Image.open(p).convert('RGB');im.thumbnail((540,430));tile=Image.new('RGB',(560,460),'white');tile.paste(im,((560-im.width)//2,28));ImageDraw.Draw(tile).text((8,4),p.stem,font=font,fill='black');tiles.append(tile)
    preview=Image.new('RGB',(1680,460*((len(tiles)+2)//3)),'#dddddd')
    for i,tile in enumerate(tiles):preview.paste(tile,((i%3)*560,(i//3)*460))
    preview.save(out/'layout_preview.png')
    write_csv(out/'figure_captions.csv',captions);write_json(out/'font_validation.json',dict(configuration=fonts,pdf_checks=checks))
