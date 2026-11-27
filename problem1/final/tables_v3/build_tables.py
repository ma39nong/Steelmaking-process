"""Editable publication tables derived from frozen source CSVs, with explicit rounding."""
from pathlib import Path
import json,re,textwrap
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
T=Path(__file__).resolve().parent;F=T.parent;P=F.parent
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
SOURCES=[]

def mdtable(headers,rows):return '|'+ '|'.join(headers)+'|\n|'+'|'.join(['---']+['---:']*(len(headers)-1))+'|\n'+'\n'.join('|'+ '|'.join(map(str,r))+'|' for r in rows)+'\n'

def latexescape(x):
 x=str(x).replace('%',r'\%').replace('_',r'\_').replace('±',r'$\pm$').replace('↑',r'$\uparrow$').replace('↓',r'$\downarrow$').replace('—',r'---')
 return re.sub(r'\*\*(.*?)\*\*',lambda m:r'\textbf{'+m.group(1)+'}',x)

def emit(name,caption,headers,rows,note,sources,widths=None):
 (T/f'{name}.md').write_text(f'**{caption}**\n\n'+mdtable(headers,rows)+'\n'+note+'\n')
 clean=[[re.sub(r'\*\*','',str(c)) for c in r] for r in rows]
 pd.DataFrame(clean,columns=headers).to_csv(T/f'data/{name}.csv',index=False)
 tex=r'\begin{tabular}{l'+('r'*(len(headers)-1))+'}\n'+r'\toprule'+'\n'+' & '.join(map(latexescape,headers))+r' \\'+'\n'+r'\midrule'+'\n'+'\n'.join(' & '.join(map(latexescape,r))+r' \\' for r in rows)+'\n'+r'\bottomrule'+'\n'+r'\end{tabular}'+'\n'
 (T/f'latex/{name}.tex').write_text('% booktabs fragment; captions/notes remain editable in paired Markdown.\n'+tex)
 if name.startswith('T'):
  n=len(rows);height=max(1.5,.29*n+1.0)
  fig,ax=plt.subplots(figsize=(7.2,height));ax.axis('off');table=ax.table(cellText=clean,colLabels=[{'Time (s)':'Time\n(s)','Top codes':'Top\ncodes','Median (counts)':'Median\n(counts)','Time corr.':'Time\ncorr.','Wave corr.':'Wave\ncorr.'}.get(h,h) if name=='T1' else h for h in headers],loc='center',cellLoc='right',colLoc='center',colWidths=widths)
  table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,1.5)
  for (r,c),cell in table.get_celld().items():
   cell.visible_edges='';cell.set_linewidth(.45);cell.PAD=.12
   if c==0:cell.set_text_props(ha='left')
   if r==0:
    cell.visible_edges='TB';cell.set_linewidth(.8);cell.set_text_props(weight='bold')
    if name=='T1':cell.set_height(cell.get_height()*1.6)
   if r==n:cell.visible_edges='B';cell.set_linewidth(.8)
   if r and '**' in str(rows[r-1][c]):cell.set_text_props(weight='bold')
  for ext in ['png','svg','pdf']:
   q=T/f'preview/{name}.{ext}';fig.savefig(q,dpi=300,bbox_inches='tight')
   if ext=='svg':q.write_text('\n'.join(x.rstrip() for x in q.read_text().splitlines())+'\n')
  plt.close(fig)
 SOURCES.append(dict(table=name,caption=caption,sources='; '.join(str(Path(x).relative_to(P.parent)) for x in sources),note=note,rows=len(rows),columns=len(headers)))


def main():
 native=pd.read_csv(F/'results/final_summary.csv');eda=pd.read_csv(P/'results/eda/statistics.csv');coords=pd.read_csv(P.parent/'data_manifest.csv');energy=pd.read_csv(P/'v3/results/singular_energy.csv');e4=energy[energy['rank']==4][['heat_id','energy']];all_data=native.merge(eda,on='heat_id').merge(coords[['heat_id','records','time_end_s','wavelength_start_nm','wavelength_end_nm']],on='heat_id').merge(e4,on='heat_id')
 rows=[[r.heat_id,f'{r.records:,}',f'{r.time_end_s:g}',f'{r.top_codes:,}',f'{r.median:.0f}',f'{r.adjacent_time_correlation:.4f}',f'{r.adjacent_wavelength_correlation:.6f}'] for r in all_data.itertuples()]
 emit('T1','表5-1 数据结构与质量特征',['Heat','Frames','Time (s)','Top codes','Median (counts)','Time corr.','Wave corr.'],rows,'每帧2047通道，间隔0.5 s。6B组波长340.535–1027.963 nm；6C组339.390–1023.303 nm，均为非均匀网格。中位数为原始矩阵的全体记录中位数。相邻相关性使用实际相邻对展平计算，时间方向每4通道抽样、波长方向每4帧抽样；不是逐帧相关系数的平均。',[str(P/'results/eda/statistics.csv'),str(P.parent/'data_manifest.csv'),str(F/'results/final_summary.csv')],widths=[.145,.085,.095,.10,.20,.13,.145])
 rows=[['Detection','Score threshold','4','Development loss'],['Detection','Scale floors','25 / 50 counts','Local residual scale'],['Detection','Connected span','3 frames x 7 channels','Local point repair'],['Background','Rank','4','Energy and development'],['Completion','Huber scale','0.003','196.605 counts'],['Completion','Rank','4','Locked development choice'],['Solver','Relative change','1e-4','Stopping criterion'],['Continuous','Interpolation','Bilinear','Held-out / convex hull']]
 emit('T2','表5-2 核心模型参数与选取依据',['Module','Parameter','Value','Basis'],rows,'Huber尺度针对Y/65535的归一化残差。检测25/50 counts分别是双向和低秩残差尺度下限。连通跨度是包围盒上限；停止阈值衡量迭代变化。工程预算、偏移窗、ADMM及halo等完整配置转附录A2。',[str(F/'config/protocol.json'),str(P/'v3/config/protocol.json')],widths=[.16,.24,.24,.36])
 u=pd.read_csv(F/'results/unified_comparison.csv');u=u[u.split=='confirmation'].set_index('strategy');keys=['A_P0','B_V1_MAD','C_V31','D_adaptive4','D_adaptive6'];names=['Retain all','Bidirectional MAD','Conservative','Adaptive S=4','Adaptive S=6']
 def bestcol(vals,direction):
  vals=np.asarray(vals,float);v=np.nanmax(vals) if direction=='max' else np.nanmin(vals);return [np.isfinite(x) and abs(x-v)<1e-10 for x in vals]
 metrics=['precision','recall','f1','introduced_unsupported_peak_count'];directions=['max','max','max','min'];flags=[]
 for col,d in zip(metrics,directions):
  vals=u.loc[keys,col].to_numpy(float,copy=True)
  if col=='precision':vals[0]=np.nan
  flags.append(bestcol(vals,d))
 rows=[]
 for i,(key,name) in enumerate(zip(keys,names)):
  r=[name]
  for k,col in enumerate(metrics):
   x='—*' if i==0 and col=='precision' else (str(int(u.loc[key,col])) if k==3 else f'{u.loc[key,col]:.3f}')
   r.append('**'+x+'**' if flags[k][i] else x)
  rows.append(r)
 emit('T3a','表5-3(a) 共同污染确认的检测性能与新增峰',['Strategy','Precision ↑','Recall ↑','F1 ↑','New peaks ↓'],rows,'* 无预测阳性，精确率未定义；源CSV按约定置0。P/R/F1由汇总TP/FP/FN计算；新增峰是全部确认案例计数。加粗为同指标最优，平局共同高亮。共同补全上限300步。',[str(F/'results/unified_comparison.csv')],widths=[.32,.18,.17,.15,.18])
 metrics=['pollution_rmse','overall_rmse','unpolluted_rmse'];flags=[bestcol(u.loc[keys,col],'min') for col in metrics];rows=[]
 for i,(key,name) in enumerate(zip(keys,names)):
  vals=[f'{u.loc[key,col]:.2f}' for col in metrics];rows.append([name]+['**'+v+'**' if flags[k][i] else v for k,v in enumerate(vals)])
 emit('T3b','表5-3(b) 同一协议的异常恢复误差',['Strategy','Polluted RMSE ↓','Overall RMSE ↓','Unpolluted RMSE ↓'],rows,'单位counts；先炉内平均、再十炉等权平均。P0和保守策略未污染误差为0，须结合污染恢复解释。发布配置使用同一确认案例但600/1200步及回退，另报告2359.80、195.80、66.68 counts及17个新增峰，不混入本300步表。',[str(F/'results/unified_comparison.csv'),str(F/'results/exact_final_pipeline_confirmation.csv')],widths=[.28,.24,.24,.24])
 c=pd.read_csv(P/'v2/results/all_models_test_macro.csv',header=[0,1],index_col=[0,1]);scenarios=['random_0.01','random_0.05','random_0.1','time_2','time_4','time_8','time_16','wave_4','wave_8','wave_16','wave_32'];labels=['Random 1%','Random 5%','Random 10%','Time 2 frames','Time 4 frames','Time 8 frames','Time 16 frames','Wave 4 channels','Wave 8 channels','Wave 16 channels','Wave 32 channels']
 def completion_rows(indices):
  result=[]
  for s,label in zip(scenarios,labels):
   if s not in indices:continue
   vals=[c.loc[(model,s),('rmse','mean')] for model in ['B0','B1','B3']];best=min(vals);row=[label]
   for model,v in zip(['B0','B1','B3'],vals):
    cell=f'{v:.2f} ± {c.loc[(model,s),("rmse","std")]:.2f}';row.append('**'+cell+'**' if abs(v-best)<1e-10 else cell)
   result.append(row)
  return result
 note='历史独立局部96×192片段协议，RMSE/counts，均值±炉间SD；随机/时间/波长为不同遮蔽场景。加粗按未舍入均值判最优。Huber在部分场景接近平方损失，正式选择同时考虑污染恢复。全矩阵配对试验另列，不与本表混合。'
 emit('T4','表5-4 代表性局部缺失补全检验',['Mask scenario','Linear','Rank-4 squared','Rank-4 Huber'],completion_rows(['random_0.05','time_4','wave_8']),note,[str(P/'v2/results/all_models_test_macro.csv')],widths=[.24,.25,.25,.26])
 rows=[[r.heat_id,str(r.repair_points),f'{r.repair_fraction*1e6:.2f}',str(r.downward_repairs),str(r.upward_repairs),str(r.numerical_fallback_points)] for r in native.itertuples()];total=native.sum(numeric_only=True);rows.append(['Total',f'{int(total.repair_points):,}',f'{total.repair_points/total.observations*1e6:.2f}',str(int(total.downward_repairs)),str(int(total.upward_repairs)),str(int(total.numerical_fallback_points))])
 emit('T5','表5-5 十炉实际修复与数值回退子集',['Heat','Repaired','Rate (ppm)','Downward','Upward','Fallback subset'],rows,'1 ppm=每百万原始记录1点；总比例由1043/27816683加权计算，为37.4955 ppm=0.00374955%。465回退点包含在1043修复点内，不额外相加。14,805个疑似顶码独立标记，1129个帧候选保留，整帧删除0。方向是估计值相对原值的变化，不等同于检测类型。',[str(F/'results/final_summary.csv')],widths=[.18,.13,.18,.15,.14,.22])
 all_data.to_csv(T/'data/A1_full_data_statistics.csv',index=False)
 emit('A1','附表A1 全部原始统计与秩4剩余能量',['Heat','P99 (counts)','E4 (%)','1-E4 (%)'],[[r.heat_id,f'{r.p99:.0f}',f'{r.energy*100:.6f}',f'{(1-r.energy)*100:.6f}'] for r in all_data.itertuples()],'原始矩阵每8帧抽样；未中心化。中心化口径另见supplement_v1/results/svd_conventions.csv。原始完整统计保留于data/A1_full_data_statistics.csv。',[str(P/'v3/results/singular_energy.csv'),str(P/'results/eda/statistics.csv')])
 core=json.loads((P/'v3/config/protocol.json').read_text());final=json.loads((F/'config/protocol.json').read_text());params=[['Iteration budget','1200'],['Per-heat wall budget','360 s'],['Stability checkpoints','600 / 1200'],['Point stability threshold','65.535 counts'],['Backtracking attempts','12'],['Step decay / growth','0.5 / 1.1'],['Time / wave scale offsets','+-2, +-3, +-4'],['Background center window','9 x 9'],['Spectral prominence window','5 channels'],['ADMM tau / rho','0.01 / 1'],['ADMM budget / tolerance','300 / 1e-4'],['Core / halo','96 x 192 / 8 x 16']]
 emit('A2','附表A2 求解与检测完整实现参数',['Parameter','Value'],params,'所有数值来自既有配置与实现；JSON源配置完整保留。随机种子、数据协议及图像参数见各源协议，不能因省略正文工程列而删除配置。',[str(P/'v3/config/protocol.json'),str(F/'config/protocol.json'),str(P/'v3/src/model3.py')])
 i=pd.read_csv(P/'v2/results/interpolation_independent_macro.csv',header=[0,1],index_col=[0,1]);rows=[]
 for scenario,label in [('time_stride_2','Time 1.0 s'),('time_stride_3','Time 1.5 s'),('time_stride_4','Time 2.0 s'),('wave_stride_2','Wave stride 2'),('wave_stride_3','Wave stride 3')]:
  vals=[i.loc[(m,scenario),('rmse','mean')] for m in ['linear','pchip','cubic']];rows.append([label]+['**'+f'{v:.2f}'+'**' if v==min(vals) else f'{v:.2f}' for v in vals])
 emit('A3','附表A3 独立连续降采样检验',['Sampling','Bilinear','PCHIP','Cubic / tensor'],rows,'RMSE/counts。三次样条与零平滑张量B样条在该边界设置下数值等价；合为一列。1.0 s的PCHIP最优，其余按源CSV呈现。',[str(P/'v2/results/interpolation_independent_macro.csv')])
 emit('A4','附表A4 全部11种局部缺失场景',['Mask scenario','Linear','Rank-4 squared','Rank-4 Huber'],completion_rows(scenarios),note,[str(P/'v2/results/all_models_test_macro.csv')])
 r=pd.read_csv(F/'results/four_reviewed_actual_repairs.csv');emit('A5','附表A5 四个高证据尖峰的坐标与修复证据',['Heat','Time (s)','Wave (nm)','Raw -> estimate','Delta','Score','Fallback'],[[x.heat_id,f'{x.time_s:.1f}',f'{x.wavelength_nm:.3f}',f'{x.original:.0f} -> {x.repaired:.2f}',f'{x.difference:+.2f}',f'{x.score:.2f}','No'] for x in r.itertuples()],'强度及差值单位counts；修复值是模型估计。实际双向参考和持续谱线保留反例见supplement_v1/results/repair_reliability_points.csv。',[str(F/'results/four_reviewed_actual_repairs.csv')])
 rows=[]
 for h in ['6C02057','6C02058','6C02107','6C02108']:
  j=json.loads((F/f'results/{h}_solver.json').read_text());hist={x['iteration']:x for x in j['history']};f=native.set_index('heat_id').loc[h];rows.append([h]+[f'{hist[k]["objective"]:.6f}' for k in [300,600,1200]]+[f'{hist[1200]["relative_change"]:.5e}',f'{100*(hist[1200]["objective"]-hist[600]["objective"])/hist[600]["objective"]:+.3f}',str(int(f.numerical_fallback_points))])
 emit('A6','附表A6 原网格保存迭代与输出稳定性',['Heat','Phi300','Phi600','Phi1200','delta1200','Phi change (%)','Fallback'],rows,'同一正式原网格配置的检查点；目标相对变化为(Φ1200−Φ600)/Φ600。不同炉目标值不用于横向比较；满足变化率标准与局部输出稳定性分开记录。',[str(F/'results/*_solver.json'),str(F/'results/final_summary.csv')])
 common=pd.read_csv(F/'results/contamination_confirmation.csv').query("strategy=='D_adaptive4'");exact=pd.read_csv(F/'results/exact_final_pipeline_confirmation.csv');means=exact.groupby('heat_id').mean(numeric_only=True).mean();emit('A7','附表A7 同一确认案例的发布配置复核',['Configuration','Polluted RMSE','Overall RMSE','Unpolluted RMSE','New peaks','Stopped'],[['300-step S=4',f'{u.loc["D_adaptive4","pollution_rmse"]:.2f}',f'{u.loc["D_adaptive4","overall_rmse"]:.2f}',f'{u.loc["D_adaptive4","unpolluted_rmse"]:.2f}',str(int(u.loc['D_adaptive4','introduced_unsupported_peak_count'])),f'{int(common.converged.sum())}/{len(common)}'],['Released 600/1200',f'{means.pollution_rmse:.2f}',f'{means.overall_rmse:.2f}',f'{means.unpolluted_rmse:.2f}',str(int(exact.introduced_unsupported_peak_count.sum())),f'{int(exact.converged.sum())}/{len(exact)}']],'同一批固定确认案例，不是第二次独立验证；RMSE为炉次宏平均。两个配置新增峰分别48和17；后者与300步阈值6的17属于不同统计结果。',[str(F/'results/exact_final_pipeline_confirmation.csv'),str(F/'results/unified_comparison.csv')])
 pd.DataFrame(SOURCES).to_csv(T/'source_mapping.csv',index=False)
 (T/'README.md').write_text('# 第一问出版表格V3\n\n正文为T1、T2、T3(a)/(b)、T4、T5五个论证主题；附录A1—A7保留完整数值。每表有可编辑Markdown、CSV、booktabs LaTeX片段；preview为排版检查，不作为唯一表源。\n\n数值从源CSV/JSON生成：P/R/F1三位小数，RMSE均值及SD两位，强度统计整数，时间相邻相关四位、波长相邻相关六位，修复率ppm两位。先比较未舍入值再标最优，平局共同加粗。相对比例按总计分母计算。\n\nLaTeX片段需在全文模板中加载booktabs；不是独立LaTeX文档。所有文件的中文标题、协议、单位和舍入说明在同名Markdown中。预览为7.2英寸单栏宽，8 pt字；数据表可作为双栏跨栏表。\n\n复现：`.venv/bin/python problem1/final/tables_v3/build_tables.py`。来源见[source_mapping.csv](source_mapping.csv)，原十表逐项审核见[table_review_v3.md](../reports/table_review_v3.md)。\n')
 print('TABLES_COMPLETE',len(SOURCES))
if __name__=='__main__':main()
