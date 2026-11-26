import json,shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ctx import *

def table(df):
 df=df.reset_index();return '| '+' | '.join(map(str,df.columns))+' |\n|'+'|'.join(['---']*len(df.columns))+'|\n'+'\n'.join('| '+' | '.join(f'{v:.5g}' if isinstance(v,(float,np.floating)) else str(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
def main():
 policy=json.loads((P/'config/final_preprocessing_policy.json').read_text());dev=pd.read_csv(P/'results/preprocessing_development.csv');test=pd.read_csv(P/'results/preprocessing_confirmation.csv');frames=pd.read_csv(P/'results/frame_diagnostics_all.csv');inj=pd.read_csv(P/'results/synthetic_diagnosis_confirmation.csv');signal=pd.read_csv(P/'results/preprocessing_signal_fidelity_development.csv');risk=pd.read_csv(P/'results/fast_change_risk_confirmation.csv')
 pd.concat([dev,test]).to_csv(P/'results/preprocessing_comparison.csv',index=False)
 for f in ['frame_anomaly_candidates.csv','channel_anomaly_candidates.csv','preprocessing_comparison.csv']:shutil.copyfile(P/'results'/f,P/f)
 shutil.copyfile(P/'config/final_preprocessing_policy.json',P/'final_preprocessing_policy.json')
 (P/'reports/异常类型与判定标准.md').write_text('''# 异常类型与分级标准

审计不等于降噪，统计候选不等于已确认故障。等级0=保留有效，1=疑似、保留质量标记，2=有明确结构性/外部证据的无效观测。饱和采用独立位标记，记录尺度保留65535。

| 类型 | 判断依据 | 默认等级及操作 |
|---|---|---|
|孤立点|双向局部中位数/MAD及低秩残差支持；可能为瞬时真峰|疑似；P1/P2掩码仅为历史统计策略对照，正式策略由预测及保真门槛选择|
|连续局部块|局部时空支持、持续长度、与前后谱对照|疑似，须外部确认范围才进入明确无效账本|
|整帧|积分、亮度、去当前帧参考光谱角、形状残差、波长覆盖、持续/恢复、饱和比例、交叉拟合低秩残差|统计候选不整行删除；确认账本的行号才全行掩码|
|持续通道|长期归一化对比与窄线持续性；真实发射线也可突出|保留标记；无外部证据不整列删除|
|疑似截顶|16bit量程及65535重复顶码|独立censored标志，保留记录值；潜在幅度条件估计，上界不可辨识|
|真实工艺突变|亮度共同改变、形状变化、前后阶段连续关系；没有故障日志|原始值保留；诊断不等同工艺因果确认|

外部确认账本config/confirmed_invalid.json支持整行和指定单元，当前为空。非有限或物理不可能的负值可作结构性无效，实际附件没有这类值。异常平坦或恒定光谱仅为候选，不能因其人工模拟可检出就认定真实仪器故障。若只确认部分波长无效，账本仅记单元，不扩大为整帧。

P0只处理明确无效，顶码从低秩拟合中排除但在记录曲面保留；P1为V2点级MAD+残差连通规则；P2增加逐帧诊断与确认账本；P3对疑似帧以0.25权重进入 ΣW hη(X−L) 的实际拟合；P4将候选整行屏蔽，仅作敏感性，不作正式故障标签。
''')
 stats=frames.groupby('heat_id').agg(frames=('candidate','size'),candidates=('candidate','sum'),structural_invalid=('structural_invalid','sum'),mean_angle=('spectral_angle_deg','mean'),mean_sat_fraction=('saturation_fraction','mean'))
 old=pd.read_csv(V2/'results/detector_test_with_raw_FPR.csv').groupby(['method','kind'])[['recall','f1','raw_false_positive_rate']].mean()
 tested=test[test.scenario!='suspect_frame_mask_sensitivity'].groupby(['heat_id','scenario','policy'])[['rmse','nrmse','masked_fraction']].mean().groupby(['scenario','policy']).mean()
 synth=inj.groupby(['heat_id','kind','level'])[['precision','recall','f1','false_positive_rate','region_iou']].mean().groupby(['kind','level']).mean()
 text='''# 数据预处理专项审核报告

基于d5b8f5e，V1/V2不改写。数据流为：原始XLSB只读审计→真实坐标与版本验证→多指标候选及证据分级→依据锁定策略建立训练掩码/权重→低秩缺失补全→观测保真与质量约束→连续插值。降噪是改变有效观测的独立操作，本轮不默认开展；候选标记并不代表故障已经确认或异常全部消除。

## 历史预处理复核

V2的孤立点召回约48%，整帧/整列/连续块更低。双向局部规则要求两轴同时偏离，天然漏掉整行增益和持久列；低秩残差也可能吸收低秩型故障。低召回不构成大范围删帧的依据，真实冶炼变化亦有相同统计特征。历史数值：
'''+table(old)+'''
## 十炉去当前帧参考诊断

参考用实际时间距离加权的前后两帧，明确排除当前帧；低秩谱基由相反奇偶帧交叉拟合，亦不含当前帧。完整清单覆盖13589帧，包括积分和实际时间差分、亮度、光谱角、形状残差、时间异常波长比例、持续长度与前后恢复、饱和占比、交叉拟合残差；sequence保留原0.5s连续序号。局部参考变化率不等于真实导数或工艺因果判断。
'''+table(stats)+'''
20%亮度和2°/5°光谱角来自历史描述性诊断，没有成为删除标准的验证依据。frame_threshold_sensitivity.csv展示10%/20%/30%与1°/2°/5°/10°的候选数量；最终多指标候选结合每炉MAD尺度，不将固定阈值当成故障标签。没有仪器日志或独立故障证据；当前没有确认真实整帧采集错误。逐帧分类和置信程度仅表达诊断等级。

## 新区域人工验证及真实快变误删风险

预声明时间间隙与波长区域避开V2局部片段，并使用7307/7907与9307/9907两组新种子；开发与确认是不同原始片段。模拟增益、偏置、整体衰减、异常平坦、多帧持续和局部波长异常，幅度0.3/0.7、持续1/3帧，并模拟ADC上限。增益/衰减可以和真实工艺变化不可区分，因此下面是候选识别性能，不是自动故障识别准确率。点级定位将帧候选展开为整行，局部波长事件可能因此出现较多误报；IoU和异或单元数同时保存。事件按连续帧匹配，事件误报率定义为未匹配事件数/片段帧数，不冒充点级TN口径。
'''+table(synth)+'''
真实快速变化片段以原积分变化上10%作风险集合，这仍不是人工真值；候选与实际自动删帧分别列出：
'''+table(risk.set_index('heat_id'))+'''
## 相同模型下P0–P4处理对照

各策略使用同一秩4、η=0.003、缺失几何和双线性插值。P3的权重真实进入损失，η不因权重改变转折点；不是事后缩放预测。P4候选帧整行重构另见suspect_frame_mask_sensitivity行，不与随机点误差混合。原始有效观测保真和积分/时间变化的开发对照：
'''+table(signal.groupby('policy')[['raw_peak_p95','integral_relative_mae','temporal_integral_difference_rmse','masked_fraction']].mean())+'''
确认片段分场景结果（先每炉平均重复种子）：
'''+table(tested)+'''
## 锁定策略与专项验收

预声明要求额外处理相对P0至少降低5%开发宏NRMSE、各场景恶化不超过10%、掩码不超过1%、原始峰变化p95不超过1%；P4不进入最终选择。开发只出现约1%的误差改善，未达到额外删除/降权的收益门槛，因此锁定 **'''+policy['policy']+'''**：保留所有没有明确无效证据的统计候选，同时输出完整质量标记。P1/P2在没有确认行时数值相同；P3亦无充分新增收益。确认实验不改选策略。

正式数据处理不是“已清除全部异常”：保留测量噪声、可疑点和真实快变，缺失/确认无效位置才补全，65535以截顶记录保留。疑似点的敏感性对照已完整保存；以后有采集日志可更新确认账本，但需要新版本再评估。工程与专项科学验收通过的含义是分级依据、可见性、风险与对照齐全，绝不声称现有数据可准确识别真实全部故障。
'''
 (P/'reports/数据预处理专项审核报告.md').write_text(text)
 # Representative candidate views, all from real data and selected by fixed metric maxima.
 for metric in ['brightness_relative','spectral_angle_deg','crossfit_lowrank_residual','saturation_fraction']:
  ix=frames[metric].abs().idxmax();row=frames.loc[ix];h=row.heat_id;y,t,w=load(h);i=int(np.argmin(abs(t-row.time_s)));lo=max(0,i-5);hi=min(len(t),i+6);fig,axs=plt.subplots(1,3,figsize=(15,4))
  for k in range(max(0,i-1),min(len(t),i+2)):axs[0].plot(w,y[k],lw=.8,label=f'{t[k]:g}s')
  axs[0].legend();axs[0].set(xlabel='Wavelength (nm)',ylabel='Intensity (recorded counts)');axs[1].plot(t[lo:hi],np.trapezoid(y[lo:hi],w,axis=1));axs[1].axvline(t[i],color='r');axs[1].set(xlabel='Time (s)',ylabel='Integrated recorded intensity');axs[2].imshow(y[lo:hi].T,origin='lower',aspect='auto',extent=[t[lo],t[hi-1],w[0],w[-1]]);axs[2].set(xlabel='Time (s)',ylabel='Wavelength (nm)');fig.suptitle(h+' '+metric+' candidate, unconfirmed');fig.tight_layout();fig.savefig(P/f'figures/appendix/frame_candidate_{metric}.png',dpi=170);plt.close(fig)
 dump(P/'results/PREPROCESSING_ACCEPTED.json',dict(policy=policy['policy'],actual_confirmed_fault_frames=0,automatic_real_frame_deletions=0,all_frames=len(frames),diagnostic_candidates=int(frames.candidate.sum()),scope='graded conservative handling, not proof of all true fault removal',protocol_sha256=sha(P/'config/protocol.json')))
if __name__=='__main__':main()
