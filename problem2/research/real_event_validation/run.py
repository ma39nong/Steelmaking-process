"""Run the real-only event review from the repository root."""
from __future__ import annotations
import argparse,copy,csv,hashlib,importlib.metadata,json,platform,subprocess
from pathlib import Path
import numpy as np
from .core import graph,intervals,temporal_recurrence,shifted_targets,unique_count
from ..response_broadening.core import changes,detect_events
from ...src.q2_pipeline import oxygen_flow
from ...src.delay_audit import write_json,write_csv

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
PRIOR=ROOT/'problem2/final/delay_audit'
RESEARCH=ROOT/'problem2/research/response_broadening'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def readcsv(path):
    with path.open() as f:return list(csv.DictReader(f))
def events_at(z,gas,base,qwindow,cwindow):
    cfg=copy.deepcopy(base);cfg['change_smoothing_samples']=cwindow
    q=z['oxygen_flow'] if qwindow==31 else oxygen_flow(z['time'],z['oxygen_cumulative_clean'],window=qwindow)
    train=np.arange(int(round(base['max_response_s']/base['dt_s'])), int(len(q)*base['kernel']['train_fraction']))
    return detect_events(changes(q,cfg),cfg,train),detect_events(changes(z[gas+'_clean'],cfg),cfg,train)

def execute(cfg,out):
    base=json.loads((ROOT/cfg['event_configuration']).read_text())
    prior=json.loads((PRIOR/'provenance.json').read_text());npzhashes={r['path']:r['sha256'] for r in prior['outputs']}
    inventory=json.loads((PRIOR/'input_audit.json').read_text())
    olddetections=readcsv(RESEARCH/'results/real/event_detections.csv')
    matches=readcsv(RESEARCH/'results/real/event_matches.csv')
    accepted=readcsv(PRIOR/'accepted_windows.csv')
    source_rows=[];edge_rows=[];summaries=[];cadence=[];scales=[];controls=[];window_rows=[];reproduced=[];cases=[];inputs=[]
    for heatinfo in sorted(inventory,key=lambda x:x['heat_id']):
        heat=heatinfo['heat_id'];raw=Path(heatinfo['path']);npz=PRIOR/'baseline/data'/f'{heat}_q2_processed.npz'
        if sha(raw)!=heatinfo['sha256'] or sha(npz)!=npzhashes[str(npz.relative_to(PRIOR))]:raise ValueError('input identity changed: '+heat)
        inputs.extend([dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in [raw,npz]])
        with np.load(npz,allow_pickle=False) as z:
            t=z['time'];duration=float(t[-1]+cfg['dt_s']);raw_o=z['oxygen_cumulative_raw'];diff=np.diff(raw_o);finite=np.isfinite(diff);positive=diff[finite & (diff>0)]
            nz=np.flatnonzero(finite & (diff>0));gap=np.diff(nz)*cfg['dt_s'];band=cfg['counter_increment_band_m3']
            cadence.append(dict(heat_id=heat,samples=len(t),time_end_s=t[-1],raw_missing=int(np.sum(~np.isfinite(raw_o))),
                equal_adjacent_fraction=float(np.mean(diff[finite]==0)),negative_increment_count=int(np.sum(diff[finite]<0)),
                positive_increment_median_m3=float(np.median(positive)),positive_increment_iqr_m3=float(np.subtract(*np.percentile(positive,[75,25]))),
                positive_increment_band_fraction=float(np.mean((positive>=band[0])&(positive<=band[1]))),update_gap_median_s=float(np.median(gap)),
                co_equal_fraction=float(np.mean(np.diff(z['co_raw'])==0)),co2_equal_fraction=float(np.mean(np.diff(z['co2_raw'])==0))))
            for gas in ['co','co2']:
                src,tgt=events_at(z,gas,base,31,9)
                old=[r for r in olddetections if r['heat_id']==heat and r['gas']==gas]
                for role,es in [('source',src),('gas',tgt)]:
                    oo=[r for r in old if r['role']==role]
                    if len(oo)!=len(es) or any(abs(float(a['peak'])-b['peak'])>1e-10 or int(a['polarity'])!=b['polarity'] for a,b in zip(oo,es)):
                        raise ValueError('candidate replay differs: '+heat+gas+role)
                edges,sd,td=graph(src,tgt,cfg['response_horizon_s']);all_edges,asd,atd=graph(src,tgt,cfg['response_horizon_s'],False)
                groups=intervals(src,cfg['response_horizon_s']);recurrence=np.zeros(len(src),int)
                for qw in cfg['flow_sg_windows']:
                    for cw in cfg['change_sg_windows']:
                        ss,gg=events_at(z,gas,base,qw,cw)
                        rr=temporal_recurrence(src,ss,cfg['event_time_tolerance_s'])
                        recurrence[list(rr)]+=1
                        scales.append(dict(heat_id=heat,gas=gas,flow_sg_samples=qw,change_sg_samples=cw,source_candidates=len(ss),gas_candidates=len(gg),
                            baseline_source_agreements=len(rr),baseline_source_count=len(src),unique_same_polarity_pairs=unique_count(ss,gg,cfg['response_horizon_s']),
                            role='posthoc sensitivity; no configuration selected by pair count'))
                source_lookup={}
                for i,s in enumerate(src):
                    future=[j for j,a in enumerate(src) if j!=i and abs(a['peak']-s['peak'])<=cfg['response_horizon_s']]
                    complete=s['start']>=0 and s['end']+cfg['response_horizon_s']<=t[-1]
                    candidates=[j for ii,j in edges if ii==i]
                    unique=len(candidates)==1 and td[candidates[0]]==1
                    state='unique_candidate' if unique else ('ambiguous' if candidates else 'unmatched')
                    lo=int(s['start']/cfg['dt_s']);hi=int(s['end']/cfg['dt_s'])+1
                    p=(z['oxygen_cumulative_clean'][s['peak_index']]-z['oxygen_cumulative_clean'][0])/(z['oxygen_cumulative_clean'][-1]-z['oxygen_cumulative_clean'][0])
                    row=dict(heat_id=heat,gas=gas,event_id=f'source_{i:04d}',source_index=i,source_start=s['start'],source_peak=s['peak'],source_end=s['end'],source_polarity=s['polarity'],
                        source_flow_m3_s=float(np.median(z['oxygen_flow'][lo:hi])),oxygen_progress=float(p),oxygen_stage=min(5,max(1,int(p*5)+1)),
                        source_horizon_overlap_group=groups[i],other_sources_within_60s=len(future),complete_60s_followup=complete,
                        same_polarity_targets=int(sd[i]),all_polarity_targets=int(asd[i]),unique_candidate=unique,status=state,
                        target_ids=';'.join(map(str,candidates)),scale_agreements=int(recurrence[i]),scale_configurations=len(cfg['flow_sg_windows'])*len(cfg['change_sg_windows']),
                        independent_event_confirmed=False,physical_match_confirmed=False,flow_regression_eligible=False,
                        evidence_status='overlapping source context' if future else 'isolated candidate without independent label')
                    source_rows.append(row);source_lookup[i]=row
                for i,j in edges:
                    s,g=src[i],tgt[j]
                    edge_rows.append(dict(heat_id=heat,gas=gas,source_id=f'source_{i:04d}',target_id=j,source_peak=s['peak'],gas_peak=g['peak'],
                        delay_onset_s=g['start']-s['start'],delay_peak_s=g['peak']-s['peak'],delay_centroid_s=g['centroid']-s['centroid'],response_width_s=g['width'],
                        source_degree=int(sd[i]),target_degree=int(td[j]),unique_edge=bool(sd[i]==1 and td[j]==1),physical_match_confirmed=False))
                rng=np.random.default_rng(cfg['random_seed']+2*len(cadence)+(gas=='co2'))
                unique=unique_count(src,tgt,cfg['response_horizon_s'])
                for k in range(cfg['circular_controls']):
                    shift=float(rng.uniform(cfg['response_horizon_s'],duration-cfg['response_horizon_s']))
                    tt=shifted_targets(tgt,shift,duration)
                    controls.append(dict(heat_id=heat,gas=gas,replicate=k,shift_s=shift,target_count_after_seam_drop=len(tt),
                        unique_pairs=unique_count(src,tt,cfg['response_horizon_s']),observed_unique_pairs=unique))
                summary=dict(heat_id=heat,gas=gas,source_candidates=len(src),gas_candidates=len(tgt),same_polarity_edges=len(edges),all_polarity_edges=len(all_edges),
                    unique_pairs=unique,ambiguous_sources=int(sum(sd[i]>0 and not(sd[i]==1 and td[[j for ii,j in edges if ii==i][0]]==1) for i in range(len(src)))),
                    unmatched_sources=int(sum(sd==0)),source_horizon_groups=len(set(groups)),isolated_source_candidates=sum(r['other_sources_within_60s']==0 for r in source_rows if r['heat_id']==heat and r['gas']==gas),
                    median_source_gap_s=float(np.median(np.diff([s['peak'] for s in src]))),complete_followup_sources=sum(s['end']+60<=t[-1] for s in src),
                    all_scales_recurring_sources=int(sum(recurrence==9)),confirmed_physical_pairs=0,independent_confirmed_events=0)
                summaries.append(summary)
                for old in [r for r in matches if r['heat_id']==heat and r['gas']==gas and r['status']=='matched']:
                    i=int(old['source_index']);row=source_lookup[i]
                    reproduced.append({**old,**{key:row[key] for key in ['other_sources_within_60s','scale_agreements','complete_60s_followup']},
                        'review_status':'candidate_only; no independent physical confirmation'})
                    cases.append(dict(case_id=f"pair_{heat}_{gas}_{old['method']}_{i:04d}",heat_id=heat,gas=gas,left=max(0,float(old['source_start'])-15),
                        right=min(t[-1],float(old['source_end'])+75),kind='previous unique pair',source_peak=float(old['source_peak']),gas_peak=float(old['gas_peak'])))
                ww=[r for r in accepted if r['heat_id']==heat and r['gas']==gas]
                for w in ww:
                    pos=np.flatnonzero(z[f'delay_{gas}_center_index']==int(w['center_index']))
                    if len(pos)!=1 or not z[f'delay_{gas}_valid_mask'][pos[0]] or abs(z[f'delay_{gas}_observed'][pos[0]]-float(w['best_delay_s']))>1e-10:
                        raise AssertionError('accepted window does not match frozen NPZ: '+heat+gas)
                    left=float(t[int(w['source_left_index'])]);right=float(t[int(w['source_right_index'])]);lag=float(w['best_delay_s'])
                    ss=[s for s in src if left<=s['peak']<=right];gg=[g for g in tgt if left+lag<=g['peak']<=right+lag]
                    window_rows.append(dict(heat_id=heat,gas=gas,center_index=int(w['center_index']),source_overlap_group=int(w['source_overlap_group']),
                        window_time_s=float(w['center_time_s']),delay_s=lag,source_candidates_in_window=len(ss),gas_candidates_in_shifted_window=len(gg),
                        physical_event_identity='uncertain',reason='correlation window is not a labelled source-response pair'))
                for group in sorted({int(w['source_overlap_group']) for w in ww}):
                    groupww=[w for w in ww if int(w['source_overlap_group'])==group]
                    cases.append(dict(case_id=f'window_{heat}_{gas}_{group:02d}',heat_id=heat,gas=gas,
                        left=max(0,float(t[min(int(w['source_left_index']) for w in groupww)])-15),
                        right=min(t[-1],float(t[max(int(w['search_right_index']) for w in groupww)])+15),kind='baseline window overlap group',
                        centers=[float(w['center_time_s']) for w in groupww],lags=[float(w['best_delay_s']) for w in groupww]))
                print(heat,gas,'sources',len(src),'unique',unique,'scale9',summary['all_scales_recurring_sources'],flush=True)
    for name,rows in [('source_candidates',source_rows),('candidate_edges',edge_rows),('heat_summary',summaries),('sampling_cadence',cadence),('scale_sensitivity',scales),('circular_controls',controls),('baseline_24_window_review',window_rows),('previous_unique_pairs',reproduced)]:
        write_csv(out/(name+'.csv'),rows)
    write_json(out/'review_cases.json',cases);write_json(out/'input_manifest.json',inputs)
    if len(window_rows)!=24 or len([r for r in source_rows if r['gas']=='co'])!=959:raise AssertionError('baseline replay counts changed')
    write_json(out/'summary.json',dict(baseline_commit=cfg['baseline_commit'],heats=len(inventory),shared_source_candidates=sum(r['source_candidates'] for r in summaries if r['gas']=='co'),
        replayed_windows=len(window_rows),previous_unique_pairs=len(reproduced),review_pages=len(cases),physical_confirmed_events=0,flow_regression_eligible_events=0,
        note='No independent ground-truth labels. Zero confirmed events is an evidence status, not proof of absent physical responses.'))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,default=HERE/'config.json');args=ap.parse_args()
    cfg=json.loads(args.config.read_text());out=HERE/'results';out.mkdir(exist_ok=True)
    protected=['problem2/src','problem2/final','problem2/research/response_broadening']
    for folder in protected:
        result=subprocess.check_output(['git','diff',cfg['baseline_commit'],'--',folder],cwd=ROOT)
        if result:raise ValueError('protected baseline modified: '+folder)
    write_json(out/'effective_config.json',cfg)
    write_json(out/'attachment1_identity.json',dict(status='resolved_same_file_renamed',user_confirmation='是原始的数据，就是这个，我只是换了一个名字罢了；就是这个',
        confirmed_original_path='Data/Attachment 1. Heat No. 6B03205.xlsb',
        sha256=sha(ROOT/'Data/Attachment 1. Heat No. 6B03205.xlsb'),
        designated_attachment_is_existing_original=True,distinct_new_data_version=False,data_identity_gate_pass=True,
        scope='User confirmed filename alias; byte hash rechecked against frozen input manifest. This does not certify all scientific conclusions.'))
    execute(cfg,out)
    from .present import render
    render(out,cfg)
    write_json(out/'provenance.json',dict(baseline_commit=cfg['baseline_commit'],execution_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        config_sha256=sha(args.config),python=platform.python_version(),source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted(list(HERE.glob('*.py'))+list(HERE.glob('*.json')))},
        dependency_versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','matplotlib','Pillow']},
        outputs=[dict(path=str(p.relative_to(out)),sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='provenance.json']))

if __name__=='__main__':main()
