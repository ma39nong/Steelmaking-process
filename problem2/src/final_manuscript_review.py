"""Editorial acceptance: independent old-table checks and protected-result hashes."""
from pathlib import Path
import csv,hashlib,json,re,subprocess
import numpy as np
import pandas as pd
from .delay_audit import write_csv,write_json

ROOT=Path(__file__).resolve().parents[2];BASE='8c8cea5';EDITORIAL_BASE='94cb15f';OUT=ROOT/'problem2/final/relationships';REPORTS=ROOT/'problem2/final/reports'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def original(path):return subprocess.check_output(['git','show',BASE+':'+str(path.relative_to(ROOT))],cwd=ROOT)


def tables(body,prefix=''):
    lines=body.splitlines();found={}
    for i,line in enumerate(lines):
        match=re.match(r'^表('+prefix+r'\d+) ',line)
        if not match:continue
        j=i+1
        while j<len(lines) and not lines[j].startswith('|'):j+=1
        cells=[]
        while j<len(lines) and lines[j].startswith('|'):
            r=[v.strip() for v in lines[j].strip('|').split('|')]
            if not all(set(v)<=set('-: ') for v in r):cells.append(r)
            j+=1
        found[match[1]]=cells
    return found


def run():
    chapter=REPORTS/'问题二_中文论文_正式章节版.md';supp=REPORTS/'问题二_补充统计与验证.md'
    body=chapter.read_text();extra=supp.read_text();prior=tables(original(chapter).decode());current=tables(body);appendix=tables(extra,'S')
    assert list(current)==list('123456') and list(appendix)==['S1','S2','S3','S4']
    assert re.findall(r'^## (6\.\d+) ',body,re.M)==['6.1','6.2','6.3','6.4']
    assert not any(word in body for word in ['SHA256','用户确认','Release','NPZ','冻结','Git','V2','复现命令','暂停新增','提交'])
    before=subprocess.check_output(['git','show',EDITORIAL_BASE+':'+str(chapter.relative_to(ROOT))],cwd=ROOT,text=True)
    assert current==tables(before), 'academic expression changed a table'
    formulas=re.findall(r'\\\[(.*?)\\\]',body,re.S)
    assert len(formulas)==8 and formulas==re.findall(r'\\\[(.*?)\\\]',before,re.S), 'academic expression changed a model formula'
    assert re.findall(r'\\\((.*?)\\\)',body,re.S)==re.findall(r'\\\((.*?)\\\)',before,re.S), 'academic expression changed inline mathematics'
    checked=[dict(item='expression_only_revision',pass_check=True,evidence='All six main tables and displayed mathematical formulas match 94cb15f exactly')]
    def compare(label,a,b):
        assert a==b,label
        checked.append(dict(item=label,pass_check=True,evidence=f'{len(a)-1} rows; every displayed cell matches pre-edit source'))
    for new,old in [('1','1'),('3','2'),('4','4'),('5','6')]:compare('main_table_'+new,current[new],prior[old])
    compare('main_table_6',current['6'],[r[:-1] for r in prior['8']])
    for new,old in [('S1','3'),('S2','5'),('S3','7')]:compare('supplement_table_'+new,appendix[new],prior[old])
    source=pd.read_csv(ROOT/'problem2/research/signal_fidelity/results/source_fidelity.csv')
    controls=pd.read_csv(ROOT/'problem2/research/signal_fidelity/results/negative_controls.csv')
    for i,method in enumerate(['formal','SG11'],1):
        peak=source.query('domain=="raw_background" and duration_s==4 and acquisition=="grid_0.01" and method==@method').peak_ratio.median()
        ripple=controls.query('signal=="oxygen" and acquisition=="empirical_hold" and method==@method').constant_flow_ripple_cv.median()
        assert current['2'][i][1:]==[f'{100*peak:.1f}',f'{100*ripple:.1f}']
    checked.append(dict(item='main_table_2',pass_check=True,evidence='Two metrics recomputed from separate explicitly defined controlled conditions'))
    chain=pd.read_csv(ROOT/'problem2/research/measurement_chain/results/background_trials.csv')
    chain=chain.query('scenario=="known_addition" and duration_s==4')
    keys=['heat_id','gas','seed','center_index','duration_s','truth_s']
    paired=chain.pivot(index=keys,columns='method',values='correct')
    assert len(paired)==1200 and paired.notna().all().all()
    assert int(paired.direct_flow.sum())==1199 and int(paired.raw_counter_chain.sum())==3
    assert '1200组配对受控条件中有1199组正确恢复时滞' in body and '正确恢复数降至3组' in body
    assert '各项均为窗口级候选估计' in body and all(len(r)==7 for r in current['6'])
    checked.append(dict(item='measurement_chain_delay_recovery',pass_check=True,evidence='1200 matched controlled conditions per entry; accepted and error <=1 s: direct flow 1199, cumulative reconstruction 3'))
    bands=pd.read_csv(OUT/'delay_flow_summary.csv')
    for gas,expected in [('co',[27.5,7.5,40.0]),('co2',[7.5,15.0,36.5])]:
        values=[float(bands[(bands.gas==gas)&(bands.flow_band==level)].delay_median_s.iloc[0]) for level in ['low','medium','high']]
        assert values==expected
    accepted=pd.read_csv(ROOT/'problem2/final/delay_audit/accepted_windows.csv')
    threshold=float(bands[bands.flow_band=='high'].flow_min_m3_s.iloc[0])
    high=accepted[(accepted.gas=='co2')&(accepted.flow_m3_s>=threshold)]
    assert len(high)==5 and set(high.heat_id)=={'6C02057'} and int((high.best_delay_s==36.5).sum())==3
    assert '呈先减小后增大的非单调特征' in body and '表现为逐层增大' in body
    checked.append(dict(item='gas_specific_delay_interpretation',pass_check=True,evidence='CO nonmonotone medians; CO2 increasing medians; high-flow CO2 has 5 windows in heat 6C02057, 3 at 36.5 s; checked against existing records'))


    trials=pd.read_csv(ROOT/'problem2/research/weak_integral/results/trials.csv')
    q=trials.query('acquisition=="ideal" and gas_input=="raw"')
    counts={k:[int(q.query('kind==@k and method==@m and noise_std==.05').correct.sum()) for m in ['ordinary','weak']] for k in ['matched','short','slow']}
    assert counts=={'matched':[44,60],'short':[30,60],'slow':[48,53]}
    assert appendix['S4'][1][1:3]==['44/60','60/60'] and appendix['S4'][2][1:3]==['30/60','60/60']
    coverage=q.query('kind=="short" and method=="weak" and noise_std==.05 and true_delay_s==15 and true_T_s==8')
    assert int(coverage.ci_covers_true_delay.sum())==5
    real=pd.read_csv(ROOT/'problem2/research/weak_integral/results/real_trials.csv')
    assert not real.query('method=="weak"').accepted_with_interval.any()
    checked.append(dict(item='supplement_table_S4',pass_check=True,evidence='Paired recovery counts, short CI coverage and real-data rejection checked against saved experiment outputs'))
    # Input fingerprints make references to new research evidence verifiable.
    for r in csv.DictReader((OUT/'manuscript_evidence.csv').open()):assert sha(ROOT/r['source'])==r['sha256']
    preserved=[]
    for path in OUT.glob('*.csv'):
        if path.name in ['manuscript_evidence.csv','manuscript_review.csv','statement_answer_mapping.csv']:continue
        probe=subprocess.run(['git','cat-file','-e',BASE+':'+str(path.relative_to(ROOT))],cwd=ROOT,capture_output=True)
        if probe.returncode==0:
            assert path.read_bytes()==original(path),path
            preserved.append(str(path.relative_to(ROOT)))
    protected=['problem1','problem2/src/q2_pipeline.py','problem2/src/q2_pipeline_v2.py','problem2/src/delay_v2.py',
               'problem2/research','problem2/final/delay_audit','problem2/configs/final_analysis.json']
    assert not subprocess.check_output(['git','diff',BASE,'--',*protected],cwd=ROOT)
    summary=pd.read_csv(OUT/'heat_summary.csv')
    assert (summary.co_total_m3>summary.co2_total_m3).all(), 'CO/CO2 total ordering changed'
    for r in summary.to_dict('records'):
        assert sha(ROOT/'problem2/final/data'/(r['heat_id']+'_q2_processed.npz'))==r['npz_sha256']
    figures=re.findall(r'!\[[^\]]*\]\(([^)]+)\)',body);assert len(figures)==6
    assert re.findall(r'^图(\d+) ',body,re.M)==list('123456')
    for name in figures:
        path=(chapter.parent/name).resolve();assert path.is_file()
        for ext in ['.png','.pdf','.svg']:
            actual=path.with_suffix(ext).read_bytes();prior_bytes=original(path.with_suffix(ext))
            if ext=='.svg':
                strip=lambda b:re.sub(rb'<dc:date>.*?</dc:date>',b'',b)
                assert strip(actual)==strip(prior_bytes)
            else:assert actual==prior_bytes
        fonts=subprocess.check_output(['pdffonts',str(path.with_suffix('.pdf'))],text=True);col=fonts.splitlines()[0].index('emb');rows=fonts.splitlines()[2:]
        assert rows and all('TimesNewRoman' in r.replace(' ','') and r[col:col+3]=='yes' for r in rows)
    for path in [chapter,supp,REPORTS/'实验配置与复现说明.md',REPORTS/'完整图表清单.md']:
        link_text=re.sub(r'\\\[.*?\\\]','',path.read_text(),flags=re.S)
        link_text=re.sub(r'\\\(.*?\\\)','',link_text,flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',link_text):
            if '://' not in link and not link.startswith('#'):assert (path.parent/link.split('#')[0]).exists(),(path,link)
    write_csv(OUT/'manuscript_review.csv',checked)
    mapping=[dict(requirement=r,location=l,assessment=s) for r,l,s in [
        ('异常处理与连续重构','6.1；表1','给出稳健尺度、PCHIP和短响应保留代价'),
        ('六类连续曲线与总量','6.2；表1、2；图1、2','保留SG31与梯形积分；量化平滑权衡'),
        ('四类耗氧关系','6.3；表3至5；图3至6','相关、增量、阶段和跨炉趋势有对应结果'),
        ('不同流量下烟气时滞','6.4；表6','给出分层候选数值；未建立可靠物理函数'),
        ('补充统计与复现','附录表S1至S4；独立交接和复现文档','移出主文的证据仍可追溯')]]
    write_csv(OUT/'statement_answer_mapping.csv',mapping)
    metadata=dict(review_base_commit=BASE,editorial_base_commit=EDITORIAL_BASE,main_sections=4,main_tables=6,supplement_tables=4,embedded_core_figures=6,
                  table_checks=checked,preserved_statistical_csvs=preserved,protected_results_unchanged=True,
                  scientific_status='No validated physical flow-delay function',editorial_status='accepted',
                  note='Old numeric verification remains valid for unchanged statistical CSVs. New controlled summaries independently checked; no new model selected.')
    write_json(OUT/'paper_review_checklist.json',metadata)
    report='''# 第二问论文审查

正文保持“连续重构—流量与排放—四类关系—时滞估计”四节结构，本轮仅优化学术表达。六张主表、六组原有核心图直接回答题目；多尺度增量、趋势系数、逐炉反积分及未采用方法对照列入四张补充表。文件改名、哈希、提交和接口验收过程均移出正文。

| 核查项 | 结果 |
|---|---|
| 原统计表 | 保留的五张主表与三张补充表数值匹配修改前正文；表6删除重复说明列并合并为表注，统计CSV未变 |
| 新增方法依据 | 核对2 s烟气响应、SG窗口权衡及1200组测量链配对恢复数；区分中位数与平均召回率 |
| 实验条件 | 脉冲保留与恒定流量伪波动分别标明条件；存储网格不解释为仪表分辨率 |
| 时滞结论 | 逐气体解释分层特征：CO先降后升，CO₂中位数递增；统计范围集中说明 |
| 未采用方法 | 弱形式只在补充材料说明采用失败，不列为正文方法 |
| 六组核心图 | 图像文件未变；正文按出现顺序嵌入，PDF保持Times New Roman嵌入 |
| 模型与数据 | 第一问、正式V2、研究实验、第三问NPZ和原统计配置均未变 |
| 生成与追溯 | 正文生成模板同步更新；重跑正文入口不会恢复工程记录 |

正文采用Hampel–PCHIP、31点SG导数、梯形积分、分阶段统计和局部互相关。已有验证支持连续工艺量与总量计算。正文分别给出CO的非单调分层形态和CO₂的递增中位数，并保留1200组配对条件中正确恢复1199组与3组的受控对照；结论采用分层候选时滞，不拟合统一递减函数。本轮六张正文表及全部显示公式与94cb15f逐项一致，未重跑实验。

[正式正文](问题二_中文论文_正式章节版.md) · [补充统计与验证](问题二_补充统计与验证.md) · [逐项核对](../relationships/manuscript_review.csv) · [新增数值来源](../relationships/manuscript_evidence.csv)
'''
    (REPORTS/'第二问最终论文审查.md').write_text(report)
    # Refresh provenance only for editorially changed outputs; preserve old execution identities.
    path=OUT/'provenance.json';prov=json.loads(path.read_text())
    prov['editorial_revision']=dict(base_commit=EDITORIAL_BASE,source='python -m problem2.src.final_manuscript',review='python -m problem2.src.final_review',scope='Presentation only; numeric arrays, formal statistics and figures unchanged')
    for r in prov['outputs']:
        p=ROOT/r['path'];r['sha256']=sha(p)
    prov['sources'].update({str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'problem2/src').glob('final_*.py')})
    included={r['path'] for r in prov['outputs']}
    for p in [supp,OUT/'manuscript_evidence.csv']:
        if str(p.relative_to(ROOT)) not in included:prov['outputs'].append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
    write_json(path,prov)
    manifest=OUT/'paper_freeze_manifest.json';old=json.loads(manifest.read_text())
    files={r['path'] for r in old['files']}
    files.update(str(p.relative_to(ROOT)) for p in [supp,OUT/'manuscript_evidence.csv',OUT/'manuscript_review.csv',ROOT/'problem2/src/final_manuscript.py',ROOT/'problem2/src/final_manuscript_review.py'])
    write_json(manifest,dict(review_base_commit=BASE,editorial_base_commit=EDITORIAL_BASE,execution_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        previous_manifest_sha256=sha(manifest),scope='Editorial revision; current manuscript and supplement; existing scientific results unchanged',
        files=[dict(path=p,sha256=sha(ROOT/p)) for p in sorted(files) if (ROOT/p).is_file()]))
    print(f'Manuscript review passed: 4 sections; 6 main tables; 4 supplement tables; {len(preserved)} statistical CSVs unchanged; 6 figures unchanged.',flush=True)


if __name__=='__main__':run()
