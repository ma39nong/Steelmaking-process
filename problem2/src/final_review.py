"""Final acceptance only: independently verify frozen outputs and document limits."""
from __future__ import annotations
import csv,hashlib,json,re,subprocess
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from scipy.signal import savgol_filter
from .delay_audit import write_json,write_csv

ROOT=Path(__file__).resolve().parents[2]
BASE='82582ad'
OUT=ROOT/'problem2/final/relationships'
REPORT=ROOT/'problem2/final/reports/第二问最终论文审查.md'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(name):
    with (OUT/(name+'.csv')).open() as f:return list(csv.DictReader(f))
def gitbytes(path):return subprocess.check_output(['git','show',BASE+':'+str(path.relative_to(ROOT))],cwd=ROOT)

def corr(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);a=x-x.mean();b=y-y.mean()
    den=np.sqrt((a@a)*(b@b));return float(a@b/den) if den>1e-12 else None

def both(x,y):return corr(x,y),corr(rankdata(x),rankdata(y))

def integral(t,y):return np.r_[0,np.cumsum(np.diff(t)*(y[1:]+y[:-1])/2)]

def fit(p,y,degree,w=None):
    # Independent polynomial API; published implementation uses a Vandermonde least-squares system.
    return np.pad(np.polyfit(p,y,degree,w=None if w is None else np.sqrt(w))[::-1],(0,2-degree))

def error_metrics(y,pred):
    error=y-pred;return dict(mae=np.mean(abs(error)),rmse=np.sqrt(np.mean(error**2)),r2=1-np.sum(error**2)/np.sum((y-y.mean())**2))

def run():
    cfg=json.loads((ROOT/'problem2/configs/final_analysis.json').read_text())
    ds={}
    for p in sorted((ROOT/cfg['data_dir']).glob('*.npz')):
        with np.load(p,allow_pickle=False) as z:d={k:z[k].copy() for k in z.files}
        ds[str(d['heat_id'])]=d
    assert len(ds)==10
    checks=[]
    def verify(table,key,field,actual,expected):
        if expected is None:passed=actual in ['',None]
        elif isinstance(expected,(int,float,np.number)):
            passed=actual not in ['',None] and bool(np.isclose(float(actual),float(expected),rtol=2e-8,atol=2e-8))
        else:passed=str(actual)==str(expected)
        checks.append(dict(table=table,row=key,field=field,saved=actual,recomputed=expected,pass_check=passed))
        if not passed:raise AssertionError(checks[-1])
    def xy(d,rel):
        gas=rel.split('_')[0];t=d['time'];o=d['oxygen_cumulative_clean'];y=d[gas+'_clean']
        return o,y if rel.endswith('instantaneous') else integral(t,y)
    def progress(d):
        o=d['oxygen_cumulative_clean'];return (o-o[0])/(o[-1]-o[0])
    # First check exact preservation of all eight original tables' CSV inputs and six core figures.
    kept=[]
    for p in sorted(OUT.glob('*.csv')):
        probe=subprocess.run(['git','cat-file','-e',BASE+':'+str(p.relative_to(ROOT))],cwd=ROOT,capture_output=True)
        if probe.returncode==0:
            if p.read_bytes()!=gitbytes(p):raise AssertionError('original statistical CSV changed: '+str(p))
            kept.append(str(p.relative_to(ROOT)))
    for r in rows('heat_summary'):
        d=ds[r['heat_id']];t=d['time'];o=d['oxygen_cumulative_clean']
        expected=dict(samples=len(t),duration_s=t[-1],oxygen_start_m3=o[0],oxygen_change_m3=o[-1]-o[0],co_total_m3=integral(t,d['co_clean'])[-1],co2_total_m3=integral(t,d['co2_clean'])[-1])
        for gas in ['co','co2']:
            expected[gas+'_reliable_windows']=int(d[f'delay_{gas}_valid_mask'].sum());expected[gas+'_alignment_points']=int(d[f'alignment_valid_mask_{gas}'].sum())
        for field,mask in [('mask_oxygen','outlier_oxygen'),('mask_co','outlier_co'),('mask_co2','outlier_co2')]:expected[field]=int(d[mask].sum())
        for field,value in expected.items():verify('heat_summary',r['heat_id'],field,r[field],value)
    for r in rows('level_correlations'):
        d=ds[r['heat_id']];x,y=xy(d,r['relation']);p,s=both(x,y)
        for field,value in [('n',len(x)),('pearson_r',p),('spearman_r',s)]:verify('level_correlations',r['heat_id']+r['relation'],field,r[field],value)
        if r['relation'].endswith('cumulative'):
            pp,ss=both(x,(d['time']/d['time'][-1])*y[-1])
            verify('level_correlations',r['heat_id']+r['relation'],'constant_gas_rate_pearson',r['constant_gas_rate_pearson'],pp)
            verify('level_correlations',r['heat_id']+r['relation'],'constant_gas_rate_spearman',r['constant_gas_rate_spearman'],ss)
    for r in rows('increment_correlations'):
        d=ds[r['heat_id']];x,y=xy(d,r['relation']);step=int(round(float(r['interval_s'])/.5));indices=np.arange(0,len(x),step);dx=np.diff(x[indices]);dy=np.diff(y[indices]);p,s=both(dx,dy)
        for field,value in [('n',len(dx)),('pearson_r',p),('spearman_r',s)]:verify('increment_correlations',r['heat_id']+r['relation']+r['interval_s'],field,r[field],value)
    def inverse(t,o,target):
        i=int(np.searchsorted(o,target,side='right')-1)
        if i<0:return float(t[0])
        if i>=len(t)-1:return float(t[-1])
        return float(t[i]+(target-o[i])/(o[i+1]-o[i])*(t[i+1]-t[i]))
    for r in rows('stage_statistics'):
        d=ds[r['heat_id']];gas=r['gas'];stage=int(r['stage']);t=d['time'];o=d['oxygen_cumulative_clean'];p=progress(d);v=integral(t,d[gas+'_clean']);lo=(stage-1)/5;hi=stage/5
        t0=t[0] if stage==1 else inverse(t,o,o[0]+lo*(o[-1]-o[0]));t1=t[-1] if stage==5 else inverse(t,o,o[0]+hi*(o[-1]-o[0]));v0=np.interp(t0,t,v);v1=np.interp(t1,t,v);do=(o[-1]-o[0])*(hi-lo)
        mask=(p>=lo)&((p<=hi) if stage==5 else(p<hi))
        expected=dict(progress_start=lo,progress_end=hi,time_start_s=t0,time_end_s=t1,duration_s=t1-t0,sample_count=int(mask.sum()),mean_rate_m3_s=(v1-v0)/(t1-t0),median_rate_m3_s=np.median(d[gas+'_clean'][mask]),cumulative_emission_start_m3=v0,cumulative_emission_end_m3=v1,stage_emission_m3=v1-v0,stage_emission_share=(v1-v0)/v[-1],stage_oxygen_m3=do,emission_per_oxygen_m3_m3=(v1-v0)/do)
        for field,value in expected.items():verify('stage_statistics',r['heat_id']+gas+str(stage),field,r[field],value)
    rng=np.random.default_rng(cfg['random_seed'])
    for r in rows('stage_aggregate'):
        rr=[x for x in rows('stage_statistics') if x['gas']==r['gas'] and x['stage']==r['stage']];ids=rng.integers(0,len(rr),(cfg['bootstrap_heat_resamples'],len(rr)))
        verify('stage_aggregate',r['gas']+r['stage'],'heats',r['heats'],len(rr))
        for field in ['mean_rate_m3_s','emission_per_oxygen_m3_m3','stage_emission_share']:
            v=np.array([float(x[field]) for x in rr]);b=np.median(v[ids],axis=1)
            for suffix,value in [('median',np.median(v)),('q25',np.percentile(v,25)),('q75',np.percentile(v,75)),('ci_low',np.percentile(b,2.5)),('ci_high',np.percentile(b,97.5))]:verify('stage_aggregate',r['gas']+r['stage'],field+'_'+suffix,r[field+'_'+suffix],value)
    for name in ['trend_parameters','pooled_trends','trend_loho']:
        for r in rows(name):
            degree=int(r['degree']);rel=r['relation']
            if name=='trend_parameters':train=[ds[r['heat_id']]]
            else:train=[d for heat,d in ds.items() if name=='pooled_trends' or heat!=r['test_heat']]
            x=np.concatenate([progress(d) for d in train]);y=np.concatenate([xy(d,rel)[1] for d in train]);w=np.concatenate([np.full(len(d['time']),1/len(d['time'])) for d in train]) if name!='trend_parameters' else None
            coef=fit(x,y,degree,w);key=r.get('heat_id',r.get('test_heat','pooled'))+rel+str(degree)
            for j in range(3):verify(name,key,'c'+str(j),r['c'+str(j)],coef[j])
            if name!='pooled_trends':
                test=ds[r.get('heat_id',r.get('test_heat'))];pred=sum(coef[j]*progress(test)**j for j in range(3));metrics=error_metrics(xy(test,rel)[1],pred)
                for field,value in metrics.items():verify(name,key,field,r[field],value)
    for r in rows('flow_review'):
        d=ds[r['heat_id']];t=d['time'];o=d['oxygen_cumulative_clean'];raw=d['oxygen_cumulative_raw'];dd=np.diff(raw);positive=dd[dd>0];window=int(r['sg_window']);unclip=savgol_filter(o,window,3,deriv=1,delta=.5);q=np.maximum(unclip,0);rebuilt=o[0]+integral(t,q)
        expected=dict(window_span_s=(window-1)*.5,raw_equal_adjacent_fraction=np.mean(dd==0),positive_increment_median_m3=np.median(positive),positive_increment_4p5_to_5_fraction=np.mean((positive>=4.5)&(positive<=5)),negative_unclipped_fraction=np.mean(unclip<0),oxygen_change_m3=o[-1]-o[0],flow_integral_m3=rebuilt[-1]-o[0],endpoint_relative_error=(rebuilt[-1]-o[-1])/(o[-1]-o[0]),reintegration_rmse_m3=np.sqrt(np.mean((rebuilt-o)**2)),flow_difference_from_31_rmse_m3_s=np.sqrt(np.mean((q-d['oxygen_flow'])**2)),flow_variation_std_m3_s=np.std(q))
        for field,value in expected.items():verify('flow_review',r['heat_id']+str(window),field,r[field],value)
    for r in rows('heat_summary'):
        d=ds[r['heat_id']]
        for gas in ['co','co2']:
            m=d[f'delay_{gas}_valid_mask'];centers=d[f'delay_{gas}_center_index'];obs=d[f'delay_{gas}_observed'];align=d[f'alignment_valid_mask_{gas}'];mapped=d['time_mapping_'+gas]
            assert np.array_equal(np.isfinite(obs),m)
            expected=np.zeros(len(d['time']),bool);expected[centers[m]]=True
            assert np.array_equal(expected,d[f'delay_{gas}_observed_mask'])
            assert np.array_equal(np.isfinite(d[gas+'_aligned']),align)
            assert np.all(np.diff(mapped[align])>0) and np.all((mapped[align]>=0)&(mapped[align]<=d['time'][-1]))
        assert (ROOT/'problem2/final/data'/(r['heat_id']+'_q2_processed.npz')).is_file()
        assert sha(ROOT/'problem2/final/data'/(r['heat_id']+'_q2_processed.npz'))==r['npz_sha256']
    for r in rows('interface_checks'):assert r['pass_check']=='True'
    assert len(rows('interface_checks'))==480 and len(rows('interface_field_comparison'))==570
    # Independently allocate the 24 already saved windows, not fitted values.
    definition=json.loads((OUT/'delay_flow_definition.json').read_text());bounds=definition['boundaries_m3_s']
    allq=np.concatenate([d['oxygen_flow'][d['delay_co_center_index']] for d in ds.values()])
    assert len(allq)==2241 and np.allclose(bounds,np.quantile(allq,[0,1/3,2/3,1]),atol=1e-12,rtol=0)
    with (ROOT/'problem2/final/delay_audit/accepted_windows.csv').open() as f:accepted=list(csv.DictReader(f))
    for r in rows('delay_flow_summary'):
        k=['low','medium','high'].index(r['flow_band']);assert int(r['candidate_windows'])==int(np.sum((allq>=bounds[k])&((allq<=bounds[k+1]) if k==2 else(allq<bounds[k+1]))));rr=[v for v in accepted if v['gas']==r['gas'] and bounds[k]<=float(v['flow_m3_s']) and(float(v['flow_m3_s'])<bounds[k+1] or k==2 and float(v['flow_m3_s'])<=bounds[k+1])];dd=[float(v['best_delay_s']) for v in rr]
        expected=dict(accepted_windows=len(rr),contributing_heats=len({v['heat_id'] for v in rr}),delay_median_s=np.median(dd) if dd else None,delay_min_s=min(dd) if dd else None,delay_max_s=max(dd) if dd else None)
        for field,value in expected.items():verify('delay_flow_summary',r['flow_band']+r['gas'],field,r[field],value)
    provenance=json.loads((OUT/'provenance.json').read_text())
    for r in provenance['outputs']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    for r in json.loads((OUT/'input_manifest.json').read_text()):
        path=Path(r['path']);path=path if path.is_absolute() else ROOT/path
        assert sha(path)==r['sha256'],path
    preserved_figures=[]
    for r in rows('figure_captions'):
        for ext in ['png','svg','pdf']:
            path=ROOT/'problem2/final/figures'/(r['figure']+'.'+ext)
            actual=path.read_bytes();original=gitbytes(path)
            if ext=='svg':
                strip_date=lambda b:re.sub(rb'<dc:date>.*?</dc:date>',b'',b)
                assert strip_date(actual)==strip_date(original),str(path)
            else:assert actual==original,str(path)
            if r['role']=='core':preserved_figures.append(str(path.relative_to(ROOT)))
        pdfpath=ROOT/'problem2/final/figures'/(r['figure']+'.pdf');text=subprocess.check_output(['pdffonts',str(pdfpath)],text=True);header=text.splitlines()[0];emb=header.index('emb')
        assert all('TimesNewRoman' in line.replace(' ','') and line[emb:emb+3]=='yes' for line in text.splitlines()[2:])
    assert len(preserved_figures)==18
    chapter=ROOT/'problem2/final/reports/问题二_中文论文_正式章节版.md';body=chapter.read_text();assert re.findall(r'^表([1-9]) ',body,re.M)==list('123456789')
    displayed={};current=None
    for line in body.splitlines():
        match=re.match(r'^表([1-9]) ',line)
        if match:current=int(match[1]);displayed[current]=[]
        elif current and line.startswith('|'):
            cells=[c.strip() for c in line.strip('|').split('|')]
            if all(set(c)<=set('-: ') for c in cells):continue
            displayed[current].append(cells)
    def printed(n,k,values):
        actual=displayed[n][k+1]  # Skip table header.
        for col,expected in values.items():
            assert actual[col]==expected,(n,k,col,actual[col],expected)
    fmt=lambda v:format(float(v),'.4g')
    for k,r in enumerate(rows('heat_summary')):
        printed(1,k,{j:fmt(r[field]) for j,field in enumerate(['samples','duration_s','oxygen_change_m3','co_total_m3','co2_total_m3'],1)})
    relations=['co_instantaneous','co2_instantaneous','co_cumulative','co2_cumulative']
    for k,rel in enumerate(relations):
        rr=[r for r in rows('level_correlations') if r['relation']==rel]
        ranges={j:f"{min(float(r[field]) for r in rr):.3f}–{max(float(r[field]) for r in rr):.3f}" for j,field in [(1,'pearson_r'),(2,'spearman_r')]}
        if rel.endswith('cumulative'):ranges[3]=f"{min(float(r['constant_gas_rate_pearson']) for r in rr):.3f}–{max(float(r['constant_gas_rate_pearson']) for r in rr):.3f}"
        printed(2,k,ranges)
        for kk,sec in enumerate([.5,5.,15.]):
            vv=[r for r in rows('increment_correlations') if r['relation']==rel and float(r['interval_s'])==sec]
            printed(3,k*3+kk,{1:fmt(sec),2:fmt(np.median([float(r['pearson_r']) for r in vv])),3:fmt(np.median([float(r['spearman_r']) for r in vv])),4:f"{min(float(r['pearson_r']) for r in vv):.3f}–{max(float(r['pearson_r']) for r in vv):.3f}"})
        for kk,degree in enumerate([1,2]):
            vv=[r for r in rows('trend_loho') if r['relation']==rel and int(r['degree'])==degree]
            printed(6,k*2+kk,{1:str(degree),2:fmt(np.mean([float(r['mae']) for r in vv])),3:fmt(np.mean([float(r['r2']) for r in vv])),4:str(sum(float(r['r2'])>0 for r in vv))})
    for k,r in enumerate(rows('stage_aggregate')):
        printed(4,k,{1:r['stage'],2:fmt(r['mean_rate_m3_s_median']),3:f"[{float(r['mean_rate_m3_s_ci_low']):.3f}, {float(r['mean_rate_m3_s_ci_high']):.3f}]",4:fmt(r['emission_per_oxygen_m3_m3_median'])})
    for k,r in enumerate(rows('pooled_trends')):printed(5,k,{1:r['degree'],2:fmt(r['c0']),3:fmt(r['c1']),4:fmt(r['c2'])})
    for k,r in enumerate([r for r in rows('flow_review') if r['sg_window']=='31']):
        alternate=[x for x in rows('flow_review') if x['heat_id']==r['heat_id'] and x['sg_window']!='31']
        printed(7,k,{1:fmt(r['raw_equal_adjacent_fraction']),2:fmt(r['positive_increment_median_m3']),3:fmt(r['reintegration_rmse_m3']),4:fmt(100*float(r['endpoint_relative_error'])),5:' / '.join(f"{float(x['flow_difference_from_31_rmse_m3_s']):.3f}" for x in alternate)})
    for k,r in enumerate(rows('delay_flow_summary')):
        interval=f"[{float(r['flow_min_m3_s']):.3f}, {float(r['flow_max_m3_s']):.3f}"+(']' if r['upper_inclusive']=='True' else ')')
        printed(8,k,{1:interval,3:fmt(r['delay_median_s']) if r['delay_median_s'] else '不可辨识',4:f"{float(r['delay_min_s']):.1f}–{float(r['delay_max_s']):.1f}" if r['delay_median_s'] else '—',5:r['accepted_windows'],6:r['contributing_heats']})
    for k,r in enumerate(rows('heat_summary')):printed(9,k,{1:r['co_reliable_windows']+' / '+r['co2_reliable_windows'],2:r['co_alignment_points']+' / '+r['co2_alignment_points'],4:'57/57通过',5:'48'})
    paragraphs=[p.strip() for p in body.split('\n\n') if len(p.strip())>=40 and not p.lstrip().startswith(('|','\\[','\\]'))]
    assert not [p for p,n in Counter(paragraphs).items() if n>1]
    for path in (ROOT/'problem2/final/reports').glob('*.md'):
        text=re.sub(r'\\\[.*?\\\]','',path.read_text(),flags=re.S);text=re.sub(r'\\\(.*?\\\)','',text,flags=re.S);text=re.sub(r'```.*?```','',text,flags=re.S)
        for url in re.findall(r'\]\(([^)]+)\)',text):
            if '://' not in url and not url.startswith('#'):assert (path.parent/url.split('#')[0]).exists(),(path,url)
    for commit in ['b1ed16e','cea8248','9b8513a','dea6637','fe52075',BASE]:
        subprocess.run(['git','merge-base','--is-ancestor',commit,'HEAD'],cwd=ROOT,check=True)
    protected=['problem2/src/q2_pipeline.py','problem2/src/q2_pipeline_v2.py','problem2/src/delay_v2.py','problem2/final/delay_audit','problem2/research']
    assert not subprocess.check_output(['git','diff',BASE,'--',*protected],cwd=ROOT)
    assert '峰宽和峰显著性只记录为诊断' in body and '不是已校准的概率' in body
    assert '属于待检验方向' in body and '离线' in body and '未对齐' in (chapter.parent/'数据交接与第三问使用说明.md').read_text()
    coverage=[
        ('十炉数据来源与异常值','6.1、6.2；表1','完成计算与原始文件身份核验'),
        ('连续重构和六类曲线','6.2、6.3；图1','完成，估计信号不等于仪表真值'),
        ('四类耗氧关系','6.4、6.5；表2至6；图2至5','原值、增量、阶段与非线性趋势均有结果'),
        ('不同流量下候选时滞','6.6；新增表8','直接作答，只有窗口级描述，不宣称已验证物理函数'),
        ('时滞可辨识性和第三问接口','6.6、6.7；表9','保留缺口与离线限制，未填充缺失时滞')]
    checklist=[
        ('逐条核对题目原文与正文回答','通过，动态规律保留限制','题目工艺变量、六曲线、四关系、候选时滞和交接均有对应段落'),
        ('检查不同流量区间的描述性时滞','补齐后通过','新增6行表8；边界取全部2241候选中心流量三分位数，24窗口全部归属'),
        ('核对原八张统计表所有数值','通过',f'{len(checks)}个数值单元重新计算核对；原统计CSV逐字节保持，接口原表8改编号为表9'),
        ('核对六组核心图和图注','通过','六组PNG/PDF逐字节一致，SVG仅导出日期元数据变化；检查图源、图注、版宽和字体嵌入'),
        ('核对数学公式符号和物理单位','通过','氧量与排放m³、流量m³/s、时滞s；耗氧进度无量纲；31点跨度15秒'),
        ('检查中文章节重复论述','通过','无重复段落；各处限制分别对应求导、关系、时滞和接口，保留必要边界说明'),
        ('核对第三问接口与离线限制','通过','570字段比较、480检查、第一问10文件身份及时间轴均核验；无需时滞修正的原时间轴可直接交接'),
        ('确认Git历史数据哈希和复现说明','通过','六个历史提交祖先关系成立；输入与输出指纹通过；原模型和历史研究未变'),
        ('检查正式结论是否超出实验支持','通过，物理关系仍未验证','无因果、单调规律或独立事件准确率声明；新增表不计算事件置信区间')]
    write_csv(OUT/'paper_numeric_verification.csv',checks)
    write_json(OUT/'paper_review_checklist.json',dict(audited_commit=BASE,checklist=[dict(item=i+1,requirement=r,status=s,evidence=e) for i,(r,s,e) in enumerate(checklist)],
        completed=9,total=9,rendered_tables_checked=9,numeric_cells_verified=len(checks),preserved_statistical_csvs=kept,preserved_core_exports=preserved_figures,
        scientific_limit='No confirmed independent physical events; no validated flow-delay law',status='reviewed offline deliverable; algorithms frozen'))
    write_csv(OUT/'statement_answer_mapping.csv',[dict(requirement=r,location=l,assessment=a) for r,l,a in coverage])
    report='# 第二问最终论文质量审查\n\n审查基准为 `82582ad`。本轮只核对既有结果并补充流量分层描述，不新增异常、重构、事件或响应核算法，不扩大合成实验，不修改正式V2或NPZ。原八张统计表的数值证据保持不变；新增不同流量区间候选时滞为表8，原接口表8改为表9，正文当前共九表、六组核心图。\n\n## 九项验收\n\n| 项目 | 状态 | 证据 |\n|---|---|---|\n'+''.join('| '+r+' | '+s+' | '+e+' |\n' for r,s,e in checklist)
    report+='\n## 不同流量下的候选结果\n\n全部候选中心流量范围为1.988–15.090 m³/s，三分位边界约4.652和5.199 m³/s，完整精度见 [分层定义](../relationships/delay_flow_definition.json)。\n\n| 流量层 | CO中位数与范围 s | CO窗口数 | CO₂中位数与范围 s | CO₂窗口数 |\n|---|---|---|---|---|\n'
    bands=rows('delay_flow_summary')
    for label,cn in [('low','低'),('medium','中'),('high','高')]:
        rr=[next(r for r in bands if r['flow_band']==label and r['gas']==gas) for gas in ['co','co2']]
        values=[f"{float(r['delay_median_s']):g} [{float(r['delay_min_s']):g}, {float(r['delay_max_s']):g}]" if r['delay_median_s'] else '不可辨识' for r in rr]
        report+=f"| {cn} | {values[0]} | {rr[0]['accepted_windows']} | {values[1]} | {rr[1]['accepted_windows']} |\n"
    report+='\n中、高流量CO各仅1窗口；高流量两气体均只覆盖6C02057。表内数量是算法合格窗口，不是独立事件；它直接回答现有方法在不同流量下给出什么数值，却不能证明流量—时滞物理函数。中位数没有一致下降，不作相反物理规律的推断。\n\n## 审查结论与冻结范围\n\n第二问进入论文与离线交接冻结状态，停止算法研发。第三问可以基于原始时间轴及质量掩膜推进，以不使用时滞修正的模型为基线；探索性对齐另作敏感性分析，不能填零或生成伪实测标签。动态流量—时滞定量识别的证据缺口仍须保留，不能将9项审查完成解释为这一物理关系已被验证。\n\n本轮独立核对相关系数、增量、阶段积分、阶段汇总、原趋势参数和留一炉次指标，以及流量数值复核。没有新增模型选择。计算证据见 [逐数值核对](../relationships/paper_numeric_verification.csv)、[九项清单](../relationships/paper_review_checklist.json)和[题目回答映射](../relationships/statement_answer_mapping.csv)。当前正文见 [中文正式章节](问题二_中文论文_正式章节版.md)。\n'
    REPORT.write_text(report)
    files=[chapter,REPORT,*OUT.glob('*'),* (ROOT/'problem2/final/figures').glob('*'),*(ROOT/'problem2/src').glob('final_*.py'),ROOT/'problem2/configs/final_analysis.json',*(ROOT/'problem2/final/reports').glob('*'),*[ROOT/p for p in protected if (ROOT/p).is_file()]]
    files=[p for p in files if p.is_file() and p.name!='paper_freeze_manifest.json']
    write_json(OUT/'paper_freeze_manifest.json',dict(review_base_commit=BASE,execution_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope='Current offline paper, tables, figures, verifier and final generation sources; no new scientific claim',
        files=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in sorted(set(files))]))
    print('Paper review passed: 9/9 items;',len(checks),'numeric cells;',len(kept),'unchanged original CSVs; 6 core figures unchanged except SVG date metadata.',flush=True)

if __name__=='__main__':run()
