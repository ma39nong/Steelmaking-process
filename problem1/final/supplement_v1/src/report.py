"""Summarize executed supplementary analyses; no historical result is overwritten."""
from pathlib import Path
import json
import numpy as np,pandas as pd
S=Path(__file__).resolve().parents[1];F=S.parent

def md(df):
 headers=list(df.columns);return '|'+ '|'.join(headers)+'|\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('|'+ '|'.join(map(str,r))+'|' for r in df.itertuples(index=False,name=None))+'\n'

def main():
 focus=pd.read_csv(S/'results/focus_channels.csv');diag=pd.read_csv(S/'results/channel_diagnostics.csv');follow=pd.read_csv(S/'results/channel_followup.csv');rel=pd.read_csv(S/'results/repair_reliability_points.csv');sv=pd.read_csv(S/'results/svd_conventions.csv');candidate=pd.read_csv(S/'results/channel_sensitivity.csv');extra=[]
 for h in focus.heat_id:
  with np.load(F/f'data/{h}_interface.npz') as z:y=z['raw'];w=z['wavelength']
  median=np.median(y,axis=0)
  for j in follow[follow.heat_id==h].channel.unique():
   j=int(j);a=max(0,j-5);b=min(len(w)-1,j+5);base=np.interp(w[a:b+1],[w[a],w[b]],[median[a],median[b]]);prom=median[a:b+1]-base;sign=1 if prom[j-a]>=0 else -1;prom*=sign;half=.5*prom[j-a];l=r=j
   while l>a and prom[l-a-1]>=half:l-=1
   while r<b and prom[r-a+1]>=half:r+=1
   width=w[r]-w[l]+np.median(np.diff(w[a:b+1]));reference=(y[:,j-1]+y[:,j+1])*.5
   extra.append(dict(heat_id=h,channel=j,wavelength_nm=w[j],median_feature_sign=sign,half_height_channels=r-l+1,half_height_span_nm=width,neighbor_intensity_correlation=float(np.corrcoef(y[:,j],reference)[0,1])))
 pd.DataFrame(extra).to_csv(S/'results/channel_width_and_support.csv',index=False)
 recurrence=diag.groupby(['group','channel','wavelength_nm']).agg(candidate_heats=('candidate','sum'),available_heats=('heat_id','count'),mean_positive_fraction=('positive_fraction','mean'),median_z=('median_abs_z','median')).reset_index();recurrence['candidate_heat_fraction']=recurrence.candidate_heats/recurrence.available_heats;recurrence[recurrence.candidate_heats>0].to_csv(S/'results/cross_heat_recurrence.csv',index=False)
 a=rel[rel.statistical_repair];bad=rel[(rel.heat_id=='6C02058')&(abs(rel.time_s-339.5)<1e-8)&(abs(rel.wavelength_nm-674.478)<.001)].iloc[0];retain=rel[~rel.statistical_repair].iloc[0]
 fh=focus[['heat_id','wavelength_nm','positive_fraction','median_abs_z','dark_residual','bright_residual','beta_ci_low','beta_ci_high']].copy()
 for col in fh.columns[1:]:fh[col]=fh[col].map(lambda v:f'{v:.5f}')
 rs=pd.read_csv(S/'results/repair_reliability_groups.csv');rs.iloc[:,1:]=rs.iloc[:,1:]*100
 for col in rs.columns[1:]:rs[col]=rs[col].map(lambda v:f'{v:.2f}%')
 selected=sv[sv['rank']==4][['heat_id','mode','energy','unexplained']].copy();selected.energy=selected.energy.map(lambda x:f'{x*100:.6f}%');selected.unexplained=selected.unexplained.map(lambda x:f'{x*100:.6f}%')
 cv=pd.read_csv(S/'results/channel_validation.csv');detection=cv[cv.method=='rank4_huber'].groupby(['split','kind'])[['tp','fp','fn']].sum().reset_index();detection['precision']=detection.tp/(detection.tp+detection.fp).replace(0,np.nan);detection['recall']=detection.tp/(detection.tp+detection.fn);detection['f1']=2*detection.tp/(2*detection.tp+detection.fp+detection.fn);detection.to_csv(S/'results/channel_detection_summary.csv',index=False)
 validation=pd.read_csv(S/'results/channel_validation_macro.csv').query('split=="confirmation"')[['kind','method','before_rmse','repair_rmse','unpolluted_rmse']].copy()
 for col in validation.columns[2:]:validation[col]=validation[col].map(lambda v:f'{v:.2f}')
 guard=pd.read_csv(S/'results/guard_validation.csv');gm=guard.groupby(['split','method'])[['protected_count','pollution_rmse','overall_rmse','unpolluted_rmse']].mean().reset_index()
 for col in gm.columns[2:]:gm[col]=gm[col].map(lambda x:f'{x:.3f}')
 report=f'''# 第一问补充研究 光谱结构与修复可靠性

## 研究对象与结论

本轮以4d28373及已发布1043点实际修复为基线，在独立目录完成十炉持续通道扫描、暗亮帧稳健回归、原始与中心化SVD、回退子集复核、有限通道污染验证和历史新增峰定位。原始10份附件及20个正式接口/差值文件的哈希均已验证，原始附件1沿用替换版本。正式算法、参数、发布数组及Release保持不变。

持续结构诊断增加了模型的数据依据：6B组835.887 nm表现为相邻参考上方的长期小幅结构，6C组383.293 nm在6/6炉通过固定的窄通道筛选，但各炉亮度斜率并不一致。原始抽样的秩4剩余能量很小，中心化后仍有清晰的共同变化结构。实际修复则须同时检查检测证据、潜变量与最终输出；宽峰案例说明数值稳定性回退和局部光谱保真可以发生取舍。

## 1 坐标与组间光谱

每炉在原生非均匀波长网格计算时间中位谱和10–90%分位带。6B四炉、6C六炉分别等权汇总：每炉曲线先线性映射到共同340.535–1023.303 nm范围的0.5 nm网格，再对同组炉次取均值。分位带是每炉分位曲线的等权平均，不是混合所有帧的总体分位带；末端共同网格点为1023.035 nm，不作越界外推。各炉完整时间积分、亮度、变化率保存在data/*_summary.npz，组内/逐炉共同网格CSV保存在results/。

两个组的光谱幅度和局部结构不同，提供分区尺度校准的依据。本轮不根据组间曲线推断工艺原因，亦不按相同列号混合两组不同波长坐标。

## 2 持续通道证据

对全2047通道使用去当前点的相邻线性预测Rλ，按实际距离加权；以相邻偏移±2/±3/±4的绝对残差中位数及25 counts下限校准。候选固定要求残差符号持续率≥95%、|残差中位数|≥250 counts、标准化绝对残差中位数≥4，边界4通道只诊断。总计20,470个通道记录全部保存；筛选结果为{len(diag[diag.candidate])}个“炉次—通道”候选。

{md(fh)}

表中符号持续率是d=Y−Rλ>0的帧比例；暗/亮使用排除目标±2 nm后的宽带平均强度，分别选最低/最高15%帧。10%、20%另作固定敏感性检查。6B目标的标准化中位偏离为2.08–3.11，未通过候选阈值；6C目标为6.89–21.77，全部通过。因此不能把先前截图的持续率当成本题已证实的坏像元比例。

回归模型d(t)=α+βB(t)+ε采用Huber IRLS，B为实际波长梯形积分去除目标窗口再除以有效带宽。β的区间来自64次、16帧移动块重采样，是各炉内部的条件性区间。6B四炉β均为正；6C的斜率正负并存，多炉区间覆盖0。这支持将固定偏置分量、亮度相关分量和局部谱形分开讨论，不能以同组出现同一突出通道就直接判为硬件故障。

9个候选中6个为383.293 nm正残差通道，另3个是6C02058/6C02060中紧邻该峰的负残差通道。相邻插值参考会使强窄峰的两侧出现负残差，因此不把9个候选解释为9个独立故障。分级清单channel_evidence_grades.csv将统计证据、物理解释及建议处理分别保存。

采样中位谱的半高连续通道数、跨度和相邻强度相关见channel_width_and_support.csv。这里的半高跨度用于描述记录网格上的局部结构，不替代未知的仪器光学FWHM。完整逐炉候选、暗亮谱、时间轨迹和跨炉复现率均提供CSV。

### 物理候选与仪器口径

[物理来源记录](physical_sources.md)核对了NIST官方Na I空气参考588.9950/589.5924 nm、K I空气参考766.48991/769.89645 nm，以及760 nm氧气A带背景。附件未说明空气/真空约定和光栅狭缝配置。每炉最近通道、参考距离、局部峰位置和间距见physical_feature_candidates.csv；距离不当作标定误差或物种确认。

6C02108保留案例304.0 s、766.413 nm的原始与正式值均为{retain.raw:.0f} counts，仍未进入修复集合；与K I空气参考相差−0.07691 nm。持续谱形及位置一致性支持保留这一物理候选特征。

### 独立通道敏感性与有限验证

仅作为敏感性输出，将6B目标835.887 nm及6C目标383.293 nm的整列替换为实际波长邻域线性参考，基于正式D保留其他节点和顶码，得到data/*_channel_sensitivity.npz。共改变{candidate.changed.sum():,}个节点，其中仅{candidate.overlap_with_baseline_repair.sum()}点与正式1043点集合重合。平均相对谱积分变化逐炉约{candidate.integral_relative_mae.min()*100:.4f}%–{candidate.integral_relative_mae.max()*100:.4f}%，平均峰高变化最高{candidate.peak_relative_mean.max()*100:.2f}%。这些是对基线的描述性变化，不是恢复误差。

新开发/确认片段位于历史局部时间槽的间隙3/4，列起点800/1600，32×96；位置清单new_panel_index.csv。固定偏置、增益和冻结扰动各取两幅度。通道恢复采用已知注入列掩码，分别用双向线性与既有秩4 Huber模型；结构检测结果单独报告，不把已知掩码恢复能力等同自动故障定位准确率。确认炉次宏平均为：

{md(validation)}

自动通道筛查在确认偏置/增益扰动上的Precision为1/3、Recall为1、F1为0.5，冻结扰动Recall为0（无阳性输出，Precision未定义）。邻通道参考会使单列突出同时在左右相邻列产生负残差，因此高召回不等同精确定位；完整检测计数见channel_detection_summary.csv。确认通道恢复两种方法的平均峰高变化均为0，新增峰均为0，积分误差与逐案例峰值指标保存在channel_validation.csv。

参考列按局部偏离较小的条件选取，故该验证衡量普通局部背景上的通道扰动恢复。它支持既有模型具备通道缺失估计能力，而真实持续谱线的归属和保留仍需独立证据。整列候选没有写入正式发布。

## 3 低秩含义与SVD口径

在6B03205、6C02058按每8帧抽样计算精确SVD；另对同一抽样矩阵按通道去时间均值，再单独分解。四个谱基的最大绝对载荷统一取正，并同步翻转时间系数。谱基单位范数，时间系数为σₖuₖ；模态不作化学过程归属。

{md(selected)}

现有历史中心化99.30%–99.86%来自V1经掩码填充的全原网格和随机SVD；本轮是原始记录的固定时间抽样及精确SVD。二者的数据对象与抽样口径不同。正文能量依据继续使用原始未中心化抽样；新增中心化结果用于说明共同动态结构，分别标注。谱基、时间系数、秩1/秩4背景及残差的数值保存于data/*_svd_*.npz，附录提供联合图。

## 4 实际修复与数值回退

掩码交集复核：1043点=465回退+578非回退；465点均为检测分类中的尖峰型，非回退含523尖峰型和55局部突降。该分类按原始值相对时间参考定义，与最终估计的向上/向下方向不同。

{md(rs)}

上表各比例以对应子集为分母。高亮为全谱积分最高15%帧，快变为相邻积分绝对变化最高15%，近顶码为11×21邻域存在65535。回退点全部位于高亮帧，近顶码邻域比例高于非回退子集；这是相关分布，不作因果归因。

### 339.5 s宽峰案例

6C02058第{int(bad.row)}行、第{int(bad.column)}通道（0基），339.5 s、674.478 nm：原始{bad.raw:.0f}，正式{bad.repaired:.5f}，差值{bad.difference:.5f} counts；潜变量{bad.latent:.5f}，时间参考{bad.time_reference:.5f}，波长参考{bad.wavelength_reference:.5f}，支撑区间[{bad.support_low:.0f},{bad.support_high:.0f}]。综合S={bad.score:.5f}，z_t={bad.time_z:.5f}，z_λ={bad.wave_z:.5f}，z_L={bad.lowrank_z:.5f}，P={int(bad.persistence)}。

该位置本身不是顶码，附近存在顶码；波长参考与原始宽峰接近，时间与检测背景证据却较大。实际低秩拟合的潜变量仍约63,602 counts，明显凹口主要出现在后续数值回退输出：既有双向填充同时使用时间和波长信息，快速亮度变化时两方向不一致，平均后降低局部强度。因此该图作为修复代价案例放附录，不作为单纯成功例。定位记录、真实参考及回退标记均可复核。

预先固定的一项保护敏感性只作用于回退点：原始与波长参考≥0.9U，且双向分歧>0.5|Rλ|时取波长参考并作既有支撑投影。新片段结果为：

{md(gm)}

开发片段的非污染区域变化有所减少，但确认片段没有触发该规则，因而没有建立独立确认收益；本轮不据此升级正式模型，也不反向修改保护阈值。原宽峰案例的分歧约47.5%，不满足这项50%规则，同样如实记录。

## 5 历史17个新增峰定位

在原正式确认表中筛出13个存在新增峰的案例，仅按原种子、位置、检测和600/1200步稳定流程诊断回放。逐案例重现数量与原表完全一致，共17处；位置与参考/污染/修复强度、是否注入/掩码/回退保存于introduced_peaks_17_locations.csv。这是历史案例解释，不是新独立实验，不混入真实1043点账本。

## 6 正文采用与版本结论

推荐进入正文的内容是：等权共同网格组间谱形、持续通道与物理候选分级、原始/中心化SVD口径解释、465/578子集分布、原波长参考与回退凹口的机制，以及实际17峰定位的可追溯性。新通道恢复实验与保护实验放在5.6及附录，并分别说明已知列掩码和新确认的条件。

保留正式自适应多证据＋固定秩4 Huber＋保真/非负支撑＋双线性模型和发布接口。独立整列线性候选只作敏感性结果：其峰高变化和物理可辨识性尚未形成升级依据。可选低秩可分离样条没有开展，避免在完成本轮优先任务后再次展开模型搜索。

文件入口：protocol.json固定协议；results/全部CSV；src/复现脚本；data/本地派生数组。后续问题继续使用已发布冻结基线，研究候选明确标记formal_release=False；第三问仍须从原始数据构造因果处理。
'''
 grade=pd.concat([diag[diag.candidate],diag[diag.target_channel]]).drop_duplicates(['heat_id','channel']).copy()
 grade['statistical_evidence']=np.where(grade.candidate,'persistent_residual_thresholds_pass','persistent_positive_structure_below_z_threshold')
 grade['interpretation']=np.where(grade.median_residual<0,'adjacent_reference_response_near_narrow_peak','persistent_narrow_structure_unassigned')
 grade['hardware_fault_confirmed']=False;grade['treatment']='retain_formal_observations; independent_sensitivity_only'
 grade.to_csv(S/'results/channel_evidence_grades.csv',index=False)
 (S/'reports/补充研究与模型依据.md').write_text(report)
 conclusions={'frozen_model_changed':False,'formal_repairs':1043,'fallback':465,'nonfallback':578,'persistent_channel_candidates':int(diag.candidate.sum()),'target_6C383_candidates':6,'target_6B835_candidates':0,'introduced_peak_replay':17,'historical_replay_cases':13,'guard_confirmation_activations':int(guard[(guard.split=='confirmation')&(guard.method=='guard_candidate')].protected_count.sum()),'recommendation':'retain frozen release; new channel candidates and protection are diagnostic/sensitivity only','optional_separable_spline':'not undertaken'}
 (S/'conclusions.json').write_text(json.dumps(conclusions,ensure_ascii=False,indent=2)+'\n')
 (S/'README.md').write_text('# 第一问补充研究\n\n[补充研究报告](reports/补充研究与模型依据.md)、[物理来源](reports/physical_sources.md)、[协议](protocol.json)、[结论](conclusions.json)。现有Release与1043点基线不变；data/为本地派生数组及明确独立的敏感性矩阵。\n\n顺序：\n\n```bash\nOPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/analyze.py\nOPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/finite_checks.py\n.venv/bin/python problem1/final/supplement_v1/src/report.py\nOPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/publication_figures.py --stage priority\nOPENBLAS_NUM_THREADS=1 .venv/bin/python problem1/final/supplement_v1/src/publication_figures.py --stage remaining\n.venv/bin/python problem1/final/tables_v3/build_tables.py\n.venv/bin/python problem1/final/supplement_v1/src/manuscript.py\n.venv/bin/python problem1/final/supplement_v1/src/validate.py\n```\n\n通道回归采用炉内时间块区间；新实验按新位置分开发/确认；历史17峰只诊断回放，不重新用于调参。输入哈希及参考册/指令哈希已登记。图与表分别在../figures/publication_v3/和../tables_v3/。\n')
 print('REPORT_COMPLETE',conclusions)
if __name__=='__main__':main()
