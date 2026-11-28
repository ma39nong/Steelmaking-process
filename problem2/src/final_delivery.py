"""Freeze final descriptive statistics, figures, text and offline handoff checks."""
from __future__ import annotations
import argparse,csv,hashlib,importlib.metadata,io,json,platform,shutil,subprocess,zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid,trapezoid
from scipy.signal import savgol_filter
from .delay_audit import write_csv,write_json,compare_tree
from .q2_pipeline import oxygen_flow
from .final_statistics import analyze

ROOT=Path(__file__).resolve().parents[2]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def audit_interface(d,q1path):
    n=len(d['time']);t=d['time'];checks=[];heat=d['heat_id']
    def check(name,value):
        checks.append(dict(heat_id=heat,check=name,pass_check=bool(value)))
    check('sequence_time_identity',np.array_equal(t,.5*(d['sequence']-d['sequence'][0])))
    check('continuous_sequence',np.all(np.diff(d['sequence'])==1))
    check('time_step',np.allclose(np.diff(t),.5,atol=1e-12,rtol=0))
    check('57_field_schema',len(d)==57)
    for key in ['oxygen_cumulative_raw','oxygen_cumulative_clean','oxygen_flow','oxygen_from_flow','co_raw','co_clean','co2_raw','co2_clean','total_co_emission','total_co2_emission']:
        check(key+'_finite_shape',d[key].shape==(n,) and np.isfinite(d[key]).all())
    for key in ['outlier_oxygen','outlier_co','outlier_co2']:
        check(key+'_bool_shape',d[key].dtype.kind=='b' and d[key].shape==(n,))
    check('nonnegative_flow',np.all(d['oxygen_flow']>=0))
    check('formal_flow_unchanged',np.allclose(d['oxygen_flow'],oxygen_flow(t,d['oxygen_cumulative_clean']),atol=1e-10,rtol=1e-10))
    check('oxygen_reintegration_array',np.allclose(d['oxygen_from_flow'],d['oxygen_cumulative_clean'][0]+np.r_[0,cumulative_trapezoid(d['oxygen_flow'],t)],atol=1e-9,rtol=1e-10))
    metadata=json.loads(str(d['metadata']))
    for gas in ['co','co2']:
        centers=d[f'delay_{gas}_center_index'];valid=d[f'delay_{gas}_valid_mask'];w=len(centers)
        windowkeys=['observed','max_corr','source_std','target_std','peak_separation','peak_prominence','peak_width_samples','boundary_hit','valid_mask','confidence','support_left_index','support_right_index']
        check(gas+'_window_shapes',all(d[f'delay_{gas}_{k}'].shape==(w,) for k in windowkeys))
        check(gas+'_centers',np.all((centers>=0)&(centers<n)) and np.all(np.diff(centers)>0))
        check(gas+'_window_time',np.array_equal(d[f'delay_{gas}_window_time'],t[centers]))
        obs=d[f'delay_{gas}_observed'];mask=d[f'delay_{gas}_observed_mask'];align=d[f'alignment_valid_mask_{gas}']
        check(gas+'_nan_semantics',np.array_equal(np.isfinite(obs),valid) and np.isnan(obs[~valid]).all())
        expected=np.zeros(n,bool);expected[centers[valid]]=True
        check(gas+'_observed_masks',np.array_equal(mask,expected) and valid.dtype.kind=='b' and mask.dtype.kind=='b' and align.dtype.kind=='b')
        mapping=d[f'time_mapping_{gas}'];expected_mapping=np.full(n,np.nan);expected_mapping[centers[valid]]=t[centers[valid]]+obs[valid]
        check(gas+'_mapping',np.allclose(mapping,expected_mapping,equal_nan=True,rtol=1e-12,atol=1e-10))
        check(gas+'_alignment_domain',np.all(~align|mask) and np.all((mapping[align]>=t[0])&(mapping[align]<=t[-1])) and np.all(np.diff(mapping[align])>0))
        aligned=d[gas+'_aligned'];check(gas+'_aligned_nan',np.array_equal(np.isfinite(aligned),align))
        check(gas+'_aligned_values',np.allclose(aligned[align],np.interp(mapping[align],t,d[gas+'_clean']),rtol=1e-10,atol=1e-10))
        model=d[f'delay_{gas}_model'];check(gas+'_model_bounds',np.all((model[np.isfinite(model)]>=0)&(model[np.isfinite(model)]<=60)))
        check(gas+'_model_mapping',np.allclose(d[f'time_mapping_{gas}_model'],t+model,equal_nan=True,atol=1e-10,rtol=1e-10))
        emission=np.r_[0,cumulative_trapezoid(d[gas+'_clean'],t)]
        check(gas+'_emission_integral',np.allclose(d['total_'+gas+'_emission'],emission,atol=1e-8,rtol=1e-10))
        check(gas+'_metadata_counts',metadata[gas+'_coverage']['reliable_window_count']==int(valid.sum()) and metadata[gas+'_coverage']['alignment_valid_count']==int(align.sum()))
    with np.load(q1path,allow_pickle=False) as q1:
        check('q1_time_exact',np.array_equal(q1['time'],t));check('q1_wavelength_2047',q1['wavelength'].shape==(2047,))
    if not all(r['pass_check'] for r in checks):raise AssertionError([r for r in checks if not r['pass_check']])
    return checks


def flow_audit(d,cfg):
    rows=[];t=d['time'];o=d['oxygen_cumulative_clean'];raw=d['oxygen_cumulative_raw'];dd=np.diff(raw);positive=dd[dd>0]
    for window in cfg['flow_sg_windows']:
        unclip=savgol_filter(o,window,3,deriv=1,delta=cfg['time_step_s']);q=np.maximum(unclip,0);rebuilt=o[0]+np.r_[0,cumulative_trapezoid(q,t)]
        rows.append(dict(heat_id=d['heat_id'],sg_window=window,window_span_s=(window-1)*cfg['time_step_s'],
            raw_equal_adjacent_fraction=float(np.mean(dd==0)),positive_increment_median_m3=float(np.median(positive)),
            positive_increment_4p5_to_5_fraction=float(np.mean((positive>=4.5)&(positive<=5))),negative_unclipped_fraction=float(np.mean(unclip<0)),
            oxygen_change_m3=float(o[-1]-o[0]),flow_integral_m3=float(trapezoid(q,t)),
            endpoint_relative_error=float((rebuilt[-1]-o[-1])/(o[-1]-o[0])),reintegration_rmse_m3=float(np.sqrt(np.mean((rebuilt-o)**2))),
            flow_difference_from_31_rmse_m3_s=float(np.sqrt(np.mean((q-d['oxygen_flow'])**2))),flow_variation_std_m3_s=float(np.std(q)),
            role='numerical sensitivity, no instantaneous flow ground truth'))
    return rows


def execute(cfg,config_path=None):
    out=ROOT/cfg['output_dir'];out.mkdir(parents=True,exist_ok=True)
    prior=ROOT/'problem2/final/delay_audit';pr=json.loads((prior/'provenance.json').read_text());known={r['path']:r['sha256'] for r in pr['outputs']}
    rawinfo=json.loads((prior/'input_audit.json').read_text());inputs=[];heats=[];checks=[];flow=[];heatrows=[];release=[]
    archive=Path.home()/'.cache/apmcm-b/releases/problem2-final-data/Steelmaking-process-q2-final-data.zip'
    expected='c1c444c56ea5827d5ea91465e86b7b98818883d6487e839f67c23524fe68c8a1'
    if not archive.exists() or sha(archive)!=expected:raise ValueError('Download verified problem2-final-data Release ZIP to local cache first')
    inputs.append(dict(role='historical_q2_release',path=str(archive),sha256=sha(archive)))
    q1_manifest=Path.home()/'.cache/apmcm-b/releases/problem1-final-data-v1/data_package_manifest.json'
    q1_sums=q1_manifest.parent/'SHA256SUMS'
    if not q1_manifest.exists() or not q1_sums.exists():raise FileNotFoundError('Q1 Release manifest and SHA256SUMS required in local cache')
    q1sha=next(line.split()[0] for line in q1_sums.read_text().splitlines() if line.split() and line.split()[-1].lstrip('*')=='data_package_manifest.json')
    if sha(q1_manifest)!=q1sha:raise ValueError('Q1 Release manifest checksum mismatch')
    q1known={r['path']:r['sha256'] for r in json.loads(q1_manifest.read_text())['files']}
    q1identity=[]
    inputs.extend(dict(role=role,path=str(p),sha256=sha(p)) for role,p in [('q1_release_manifest',q1_manifest),('q1_release_checksums',q1_sums)])
    canonical=ROOT/'problem2/final/data';canonical.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive) as package:
        for info in sorted(rawinfo,key=lambda r:r['heat_id']):
            heat=info['heat_id'];raw=ROOT/cfg['raw_data_dir']/Path(info['path']).name;path=ROOT/cfg['data_dir']/f'{heat}_q2_processed.npz';q1=ROOT/cfg['q1_data_dir']/f'{heat}_interface.npz'
            if sha(raw)!=info['sha256'] or sha(path)!=known['baseline/data/'+path.name]:raise ValueError('frozen input identity changed '+heat)
            if not q1.is_file():raise FileNotFoundError('Q1 local interface required for actual cross-question validation: '+str(q1))
            q1hash=sha(q1)
            if q1hash!=q1known[str(q1.relative_to(ROOT))]:raise ValueError('Q1 local interface does not match Release '+heat)
            q1identity.append(dict(heat_id=heat,path=str(q1.relative_to(ROOT)),sha256=q1hash,release_sha256_match=True))
            inputs.extend(dict(role=role,path=str(p.relative_to(ROOT)),sha256=sha(p)) for role,p in [('raw_xlsb',raw),('baseline_npz',path),('q1_read_only',q1)])
            with np.load(path,allow_pickle=False) as z:d={k:z[k].copy() for k in z.files}
            d['heat_id']=str(d['heat_id']);heats.append(d);checks.extend(audit_interface(d,q1));flow.extend(flow_audit(d,cfg))
            member=next(name for name in package.namelist() if Path(name).name==path.name)
            with np.load(io.BytesIO(package.read(member)),allow_pickle=False) as historical:
                for key in d:
                    a=historical[key];b=np.asarray(d[key]);numeric=a.dtype.kind in 'biufc'
                    equal=bool(np.array_equal(a,b,equal_nan=True)) if numeric else bool(np.array_equal(a,b))
                    near=bool(np.allclose(a,b,atol=1e-8,rtol=1e-8,equal_nan=True)) if numeric else equal
                    if key=='metadata':near=not compare_tree(json.loads(str(a)),json.loads(str(b)))
                    strict='mask' in key or key in ['time','sequence'] or key.endswith(('center_index','window_time'))
                    if a.shape!=b.shape or not near or(strict and not equal):raise AssertionError('Release field changed '+heat+key)
                    release.append(dict(heat_id=heat,field=key,shape=str(b.shape),dtype=str(b.dtype),exact_historical_match=equal,numerical_historical_match=near,strict_identity_required=strict))
            destination=canonical/path.name
            if destination.exists() and sha(destination)!=sha(path):
                with np.load(destination,allow_pickle=False) as existing:
                    if set(existing.files)!=set(d):raise ValueError('canonical NPZ fields differ: '+str(destination))
                    for field in existing.files:
                        aa=existing[field];bb=np.asarray(d[field])
                        same=aa.shape==bb.shape and aa.dtype==bb.dtype
                        if field=='metadata':same=same and not compare_tree(json.loads(str(aa)),json.loads(str(bb)))
                        elif 'mask' in field or field in ['time','sequence'] or field.endswith(('center_index','window_time')):same=same and np.array_equal(aa,bb)
                        elif aa.dtype.kind in 'biufc':same=same and np.allclose(aa,bb,equal_nan=True,atol=1e-10,rtol=1e-10)
                        else:same=same and np.array_equal(aa,bb)
                        if not same:raise ValueError('canonical NPZ differs; investigate before replacing '+str(destination)+' '+field)
            if not destination.exists():shutil.copyfile(path,destination)
            if not destination.is_file():raise AssertionError('handoff file absent')
            metadata=json.loads(str(d['metadata']));row=dict(heat_id=heat,samples=len(d['time']),duration_s=float(d['time'][-1]),oxygen_start_m3=float(d['oxygen_cumulative_clean'][0]),
                oxygen_change_m3=float(d['oxygen_cumulative_clean'][-1]-d['oxygen_cumulative_clean'][0]),co_total_m3=float(d['total_co_emission'][-1]),co2_total_m3=float(d['total_co2_emission'][-1]),
                co_reliable_windows=int(d['delay_co_valid_mask'].sum()),co2_reliable_windows=int(d['delay_co2_valid_mask'].sum()),co_alignment_points=int(d['alignment_valid_mask_co'].sum()),co2_alignment_points=int(d['alignment_valid_mask_co2'].sum()),
                mask_oxygen=int(d['outlier_oxygen'].sum()),mask_co=int(d['outlier_co'].sum()),mask_co2=int(d['outlier_co2'].sum()),npz_sha256=sha(destination),baseline_npz_sha256=sha(path))
            heatrows.append(row)
    tables=analyze(heats,cfg)
    from .final_delay_summary import summarize
    delay_bands,delay_definition=summarize(heats,cfg['delay_flow_summary_quantiles'])
    tables['delay_flow_summary']=delay_bands
    write_json(out/'delay_flow_definition.json',delay_definition)
    tables.update(heat_summary=heatrows,flow_review=flow,interface_checks=checks,interface_field_comparison=release,q1_interface_identity=q1identity)
    for name,rows in tables.items():write_csv(out/(name+'.csv'),rows)
    write_json(out/'effective_config.json',cfg);write_json(out/'input_manifest.json',inputs)
    write_json(out/'acceptance.json',dict(data_identity_resolved=True,heats=len(heats),q1_time_axes_exact=True,q1_release_hashes_match=True,historical_fields=len(release),interface_checks=len(checks),all_interface_checks_pass=True,
        observed_co_windows=sum(r['co_reliable_windows'] for r in heatrows),observed_co2_windows=sum(r['co2_reliable_windows'] for r in heatrows),
        scientific_flow_delay_law_validated=False,delay_method_development_paused=True,interface_scope='offline computation and schema validation, not physical validation of delays'))
    from .final_figures import render
    from .final_text import write_reports
    render(heats,tables,cfg);write_reports(heats,tables,cfg)
    protected=['problem2/src/q2_pipeline.py','problem2/src/q2_pipeline_v2.py','problem2/src/delay_v2.py']
    if subprocess.check_output(['git','diff',cfg['baseline_commit'],'--',*protected],cwd=ROOT):raise AssertionError('formal numeric model changed')
    outputs=[p for p in out.rglob('*') if p.is_file() and p.name!='provenance.json' and not p.name.startswith('paper_') and p.name!='statement_answer_mapping.csv']
    outputs+=list((ROOT/cfg['figure_dir']).glob('*'))
    outputs.append(ROOT/'problem2/final/results/artifact_manifest.json')
    reports=['第一问成果审查与第二问接口说明.md','独立数学审查与修复报告.md','问题二_中文论文_正式章节版.md','问题二_补充统计与验证.md','第二问最终验收报告.md','数据交接与第三问使用说明.md','数据字典与接口说明.md','完整图表清单.md','实验配置与复现说明.md','problem2_section.tex']
    outputs+=[ROOT/'problem2/final/reports'/name for name in reports]
    write_json(out/'provenance.json',dict(baseline_commit=cfg['baseline_commit'],execution_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        execution_source_state='Working-tree source hashes record the executed code; execution_commit is the pre-commit HEAD.',
        git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        python=platform.python_version(),dependencies={n:importlib.metadata.version(n) for n in ['numpy','scipy','matplotlib','Pillow']},
        sources={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'problem2/src').glob('final_*.py'))},
        protected_sources={name:sha(ROOT/name) for name in protected},config_sha256=sha(config_path or ROOT/'problem2/configs/final_analysis.json'),
        outputs=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in sorted(set(outputs))]))
    print('Final analysis complete:',len(heats),'heats;',len(release),'historical fields;',len(checks),'interface checks',flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,default=ROOT/'problem2/configs/final_analysis.json');a=p.parse_args();execute(json.loads(a.config.read_text()),a.config)

if __name__=='__main__':main()
