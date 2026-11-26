"""Build Chinese report solely from stored measured results; no invented values."""
import json,hashlib
import numpy as np,pandas as pd
from audit import ROOT,P,HEATS
from reports_initial import table

def fmt(df):
 d=df.copy()
 for c in d:
  if pd.api.types.is_float_dtype(d[c]):d[c]=d[c].map(lambda v:f'{v:.5g}' if np.isfinite(v) else '不适用/无有效评分点')
 return table(d)
def main():
 required=['results/validation/COMPLETE','results/final/COMPLETE','results/validation/spline_sensitivity.csv','results/validation/artificial_clipping.csv','results/model_comparison/resources_isolated.csv']
 for r in required:
  if not (P/r).exists():raise RuntimeError('Missing '+r)
 select=json.loads((P/'results/final/model_selection.json').read_text());f=select['family'];q=pd.read_csv(P/'results/audit/data_quality_summary.csv');a=pd.read_csv(P/'results/anomaly_detection/summary.csv');r=pd.read_csv(P/'results/model_comparison/model_comparison.csv');inj=pd.read_csv(P/'results/validation/anomaly_injection.csv');beads=pd.read_csv(P/'results/validation/beads_experiment.csv');fin=pd.read_csv(P/'results/final/reconstruction_summary.csv');sat=pd.read_csv(P/'results/final/saturation_sensitivity.csv');clip=pd.read_csv(P/'results/validation/artificial_clipping.csv');mem=pd.read_csv(P/'results/model_comparison/resources_isolated.csv');dev=pd.read_csv(P/'results/validation/development_scores.csv')
 if r.n.min()<=0 or not np.isfinite(r.rmse).all():raise RuntimeError('Empty or nonfinite independent validation scores')
 if 'isolated_peak_rss_mb' not in r:r=r.merge(mem.drop(columns=['solver_converged']),on=['heat_id','model'],validate='many_to_one');r.to_csv(P/'results/model_comparison/model_comparison.csv',index=False)
 agg=r.groupby(['model','scenario'],as_index=False)[['mae','rmse','nrmse','sam_deg','seconds','empirical_band_coverage']].mean()
 perf=r.groupby(['heat_id','model'],as_index=False)[['mae','rmse','nrmse','anomaly_f1']].mean()
 chosen=r[r.model==f];selected_params=pd.DataFrame([dict(heat_id=h,**select['per_heat'][h]['models'][f]['params']) for h in HEATS]);conv=r.groupby('model').solver_converged.agg(['sum','count']).reset_index()
 doc=fr'''# 问题一：数据分析、模型建立与实验结果

研究日期：2026-11-26；分支`problem1`。本报告使用真实十炉数据，模型族在开发遮蔽上选定，再报告独立测试。当前有限候选集选中**模型{f}**。这表示本配置、本验证形态下的选择，不是全部模型空间的数学最优结论。

## 第一章 问题重述

输入十炉次离散时间—波长—光强矩阵，输出不可靠观测类型、清洗掩码、连续函数及十张3D谱面。只用光谱作为拟合输入；氧耗和烟气列仅审计。原始数据保持只读，十炉独立真实坐标，不拼接不同波长网格。

## 第二章 问题分析

EDA显示相邻λ的相关性显著高于相邻时间，宽带辐射形状与局部窄峰共存，部分炉次有剧烈整体亮度变化。强时间平滑与小时间基函数维数容易吞掉真实变化。异常需区分孤立尖峰、整帧状态变化、持续通道特征与上量程截顶。

6C02058代表性早中晚谱可见宽峰形状大体相似、幅度明显变化，770nm附近强时间跳变；6B03205热力图亦呈细密横向亮度条带，不能将全部快速时间变化当噪声。图见[EDA](04_探索性数据分析.md)。

## 第三章 模型假设与符号

t=(Sequence−Sequence1)×0.5s，λ从各文件表头读取；允许非平稳和真实极端值，不默认噪声Gaussian。只在观测域重构。Y为观测，F为连续信号，S异常分量，M有效掩码，C为样条系数，α_t/α_λ分方向粗糙度，δ为Huber阈值。完整假设、单位见[符号说明](02_模型假设与符号说明.md)。题目光强标签为W/cm²，但缺少原始辐照度校准，图中标recorded units，未自行换算。

## 第四章 数据来源与预处理

新版6B03205哈希与指定`8c73aaefadb9c8d8c3f83b193d9b4583dcfdfb239446c69b7eacd4120820783c`一致，虽然文件名带重复.xlsb后缀，按内容确认版本，未使用旧附件1。映射见[manifest](../../data_manifest.csv)。

{fmt(q[['heat_id','records','minimum','maximum','sat65535','sat_fraction']])}

十炉共{int(q.records.sum()):,}条时间记录，每条2047通道。零值、负值、NaN、Inf、Sequence断点、重复光谱均为0。前4炉实际波长340.535–1027.963nm，后6炉339.390–1023.303nm，均非均匀；与题目正文终点不一致，以表头为准。6C02057仅775条，不强行补到1400条。

## 第五章 异常识别

双向局部中位数残差分别作MAD标准化，只有时间与波长两方向同时显著偏离的点才作较强孤立候选。持续谱线和共同亮度跃变只在一向突出时保留。窗口/阈值由独立种子的人工异常开发集选择，在另一种子注入集测试。

掩码类型位：1=不可见/非有限，2=65535截顶候选，4=双向局部异常，8=低维PCA+MCD整帧候选，16=持续通道候选。**8/16只作诊断，不自动删除，也不排除真实高强度帧的验证资格。** PCA取6维后MCD，避免在2047维直接估协方差。原始MCD显著性不具有非平稳时间序列的严格p值含义。

{fmt(a)}

不按固定百分位裁剪。mask外存NaN，重构值另存，不把65535置零。局部检测是保守工程规则，不能保证完全区分真实瞬时窄峰与仪器尖峰。

### 注入实验（十炉全部报告）

下表是独立测试种子下每种异常的十炉宏平均；人工异常不是实测故障标签。

{fmt(inj[inj.split=='test'].groupby('kind',as_index=False)[['precision','recall','f1']].mean())}

整行/整列异常可能表现为真实整体变化，保守双向局部检测不保证高召回。完整逐炉逐类型数据见`results/validation/anomaly_injection.csv`；低F1必须保留，不能用易检尖峰的结果概括所有异常。

## 第六章 连续曲面模型

### A：局部多项式基线

训练掩码填补只见可见数据；均匀时间轴使用SG，波长使用真实坐标局部二次多项式：V_jk=(λ_j−λ0)^k，a=V†y。随后在真实轴上构建三次张量样条，形成任意域内坐标可计算的F_A(t,λ)。这不是对非均匀λ直接套均匀SG。

三点二次时间SG在完整点上等于恒等映射；五点窗口才有实际时间平滑。将这一退化候选保留，是为了检验强时间平滑是否必要，不能将三点情形称为已进行了有效时间去噪。

### B：低秩—稀疏补全、观测保真投影与连续插值

归一化Z=Y/65535，求解

$$\min_{{L,S}}\ \frac12\|M\odot(Z-L-S)\|_F^2+\eta\|M\odot S\|_1,\quad rank(L)\le r.$$

S由软阈值更新，L经可见残差步和截断SVD迭代投影。固定随机种子，保存目标及变化率。它是本题二维简化，**没有实现R03原始局部logdet/按列log稀疏/三维TV，不称完整L3S3TV复现**。候选设置r=4或12、η=.003或.01，最多120步，达到相对变化停止准则才记收敛。

为避免低秩近似在本来有效的通道造出窄峰，加入观测保真投影：D=M⊙Y+(1−M)⊙(65535L)，再沿真实轴执行顺序PCHIP：先逐实际波长对时间插值，再对任意目标时刻沿波长插值。标准双三次插值曾在采样点之间产生明显负值，因此改用具有区间形状保持性质的确定性连续求值，不对最终曲面作视觉裁剪。最终曲面本身不要求低秩，低秩只约束潜在补全层。有效观测在原网格上保持不变，同时也保留其中测量噪声；只有已遮蔽点由低秩层填补。数学表达为F_B(t,λ)=H_λ(λ; {{(λ_j,H_t(t; {{(t_i,D_ij)}}_i))}}_j)，H表示实际坐标上的分段三次Hermite形状保持插值。该顺序构造连续，不宣称固定二维B样条系数或全局C²。补充依据S13为SciPy官方PchipInterpolator算法文档。保存原网格D、潜在L与两轴可唯一复建该插值函数，禁止外推。

### C：鲁棒各向异性张量B样条

$$F_C(t,λ)=65535\sum_{{p=1}}^{{K_t}}\sum_{{q=1}}^{{K_λ}}c_{{pq}}B_p(u)B_q(v).$$

u/v为逐炉归一化真实坐标；在等距节点上建三次样条，但观测网格无需均匀。目标为固定训练噪声尺度ν下的Huber误差，加两个二阶系数差分惩罚和可选一阶混合惩罚。IRLS权min(1,δ/|r|)仅乘一次。每步CG用算子B_tᵀ[W⊙(B_tCB_λᵀ)]B_λ和惩罚，不建数百万行的稠密Kronecker设计矩阵。

完整推导、权重和预条件器见[算法来源与推导](../../reports/00_数学模型推导与算法来源.md)。候选K_t=64/192、K_λ=256/512，有限容量并不代表该方法全部潜力。CG与IRLS分别记录，达到最大步数不当作完全收敛。

### BEADS独立实验

以R10原式(49)/(63)–(68)为依据，调用pybaselines已实现的BEADS；明确重采样为均匀λ，再归一化。每炉早/中/晚各比较R04参数(.08,.1,1,.8)和更弱(.008,.001,.01,.008)，分别保存峰x、基线f和残差；保留x+f用于原信号估计。

共运行{len(beads)}次，满足停止准则{int(beads.converged.sum())}次。只保留峰的平均面积相对变化为{beads.peak_only_area_relative_error.mean():.4g}，中位基线面积占比为{beads.baseline_area_fraction.median():.4g}。这些结果说明去基线会显著改变原始辐射总量；不能因R04用于预测特征就强制用于本问完整光强。详见`beads_experiment.csv`与30组对照图。

## 第七章 模型求解与结果

### 参数选择

开发集同时包含随机遮蔽、4行时间带与12通道波长带；测试集单独含随机1%/5%/10%、6行时间带与16通道波长带。本问为离线重构，可以使用缺失段前后可见谱，并非未来预测。参数开发阶段不见所有外层测试遮蔽并集，外层独立评估按各场景重新拟合。整段时间/λ验证隐藏整个带，评分仅用非强局部异常且Y<64000的可参考观测，近量程保护阈值只用于评分资格，不用于裁剪有效数据或删去拟合输入。此高可信筛选使误差结论有条件，不能代表饱和/近饱和区的真实恢复误差。

开发宏平均NRMSE：{json.dumps(select['development_macro_nrmse'],ensure_ascii=False)}。因此锁定模型{f}，后续测试仅作评价。不能根据测试表倒改选型。修改过一次协议：不再将MCD整帧候选排除评分，以解决空时间测试集合；同时增加求解预算，修正依据是可评分性和收敛状态，而非追逐测试误差。旧试运行缓存标为superseded。

### 十炉×三模型独立测试（五场景宏平均）

{fmt(perf)}

此表各场景等权宏平均，随机遮蔽场景可相互重叠，不能把它们当独立样本累计标准误。完整150行和耗时、内存、SAM、峰误差位于[模型对比CSV](../results/model_comparison/model_comparison.csv)。

### 分遮蔽形态比较

{fmt(agg)}

不能只凭随机缺失误差判断完整时间段恢复。真实光谱快速变化没有同一时刻可见信息时，本质上更难。λ带为局部范围，其SAM/峰位/积分误差仅针对该评分带；时间带则对应隐藏谱帧的高可信通道集合。随机遮蔽不计算不完整峰形的峰位误差，CSV记NaN，不能解释为0。

### 低秩假设的数据证据

{fmt(pd.read_csv(P/'results/eda/low_rank_energy.csv'))}

上述能量谱是全数据描述性审计，不用于从测试误差倒选参数；未中心化能量会受整体亮度支配，因此同时报告按通道减均值的结果。高能量解释比例仍不保证窄峰和稀有状态无损。

### 最终参数及十炉重构

{fmt(selected_params)}

{fmt(fin[['heat_id','model','solver_converged','negative_fraction','dense_negative_fraction','saturated_count','dense_time','dense_wavelength']])}

最终拟合使用验证选定的参数与求解预算。密集评价网格各轴为2n−1、2m−1；它来自明确连续函数，不是仅增加图片网格。3D图用160×360评价点以保持可读性，完整密集数值另存。密集网格负值比例据实列出；模型B顺序PCHIP的区间范围性质避免普通双三次插值的负过冲，并未对图形做后处理裁剪。

含饱和/异常遮蔽区域的峰高包含模型推断；图形不是对真实超量程强度的直接测量，具体区域另有不可辨识掩码。

![十炉次连续曲面总览](../figures/reconstructed/all_heats_surfaces.png)

''' + ''.join(f'''### {h}

![连续3D曲面](../figures/reconstructed/{h}_surface.png)

![重构热力图](../figures/reconstructed/{h}_heatmap.png)

![真实曲线处理前后](../figures/comparison/{h}_before_after.png)

''' for h in HEATS)+f'''
## 第八章 验证、敏感性与不确定性

### 观测保真与峰形复核

视觉核查曾发现纯低秩曲面在6C02058低亮度t=445.5s的有效通道造出额外窄峰，因此采用上述有明确硬约束的观测保真投影。潜在优化和缺失点估计没有改变；50项模型B独立遮蔽测试已逐项重新拟合、重新评分，并检验有效观测差为0、遮蔽预测差为0，记录见`observation_projection_equivalence.csv`，未据测试误差倒选参数。

{fmt(pd.read_csv(P/'results/validation/peak_retention_descriptive.csv'))}

原始峰可能含异常，表中差异不是干净真峰误差；原网格有效点保真不能保证遮蔽的瞬时真实峰被正确识别，形状保持插值仅保证其数学区间性质，不证明未观测时刻的真实峰高；真实饱和区仍不可辨识。保留每炉最大原始峰差异图，避免只展示有利时刻。

![修正后最大原始峰差异](../figures/comparison/6C02058_largest_raw_peak_difference.png)

### 指标定义

MAE和RMSE仅在独立遮蔽参考点上统计；NRMSE按该评分点的q95−q05归一化。SAM按每帧可共同评分的通道子集计算弧度后转度，至少3通道。峰位单位nm，峰高和面积误差为相对值；λ带局部峰不可误称全谱峰。拟合点残差另存为非验证列，不用于宣称重构精度。

### 参数敏感性

在预先指定6B03205、6C02058、6C02108三炉完整网格上单因素改变α_t、α_λ、K_t、K_λ、δ、混合惩罚，固定其余参数；使用独立探索遮蔽，不反向调整已锁定最终模型。见`spline_sensitivity.csv`及六组敏感性图。异常阈值—F1图来自十炉开发注入；不能把3炉敏感性说成全部炉次置信结论。

![时间平滑敏感性](../figures/comparison/sensitivity_alpha_t.png)

![波长基函数数量敏感性](../figures/comparison/sensitivity_kw.png)

![异常阈值敏感性](../figures/comparison/anomaly_threshold_sensitivity.png)

### 饱和与删失

R09的16bit ADC与65535重复上限支持截顶解释，但缺少采集软件配置，不能确定全部限幅机制。常规路线将其掩码。另对8炉比较C掩码与单侧平方下界惩罚η/2Σ[1−F/65535]_+²，**不当作等式真值，也不称Gaussian删失似然**。

{fmt(sat)}

选定模型的掩码预测若低于65535必须标出；这种区域不能据未截顶邻域唯一推定真实超量程峰高。单侧约束减少下界违背，但不能给出有限上界或证明恢复了真峰。

人工截顶仅在已知未饱和参考的连续局部块上降低阈值，十炉比较四个处理方式，真实结果为

{fmt(clip.groupby('model',as_index=False)[['rmse','mae','nrmse','below_cap_fraction']].mean())}

这衡量已知参考场景的恢复能力，不能代表真实长饱和区的真实误差。保存原始参考、人工截顶观测、掩码与四模型预测以便复核。

### 经验预测带

分别按开发集随机/时间/λ带残差的95%绝对分位数形成经验预测半径，独立测试覆盖率已列在场景比较表。它是观测误差带，没有对非独立/异方差数据提供严格95%保证。真实饱和区额外输出不可辨识掩码，不套用这些带声称真实超量程值已被覆盖。

### 求解与资源

{fmt(conv)}

上述sum/count为达到停止准则的外层拟合数。若未收敛，结果属于记录预算下的近似解，不宣称已达到相应目标的最优点。每个CG线性残差、IRLS目标和低秩变化率见`*_solver.json`及最终模型diagnostics。

{fmt(mem.groupby('model',as_index=False)[['fit_and_grid_seconds','isolated_peak_rss_mb','pre_fit_peak_rss_mb']].mean())}

内存为每个炉次/模型启动独立进程的实际peak RSS，包括解释器、输入与拟合临时数组；不把主进程累积高水位当单模型内存。三模型共享同一检测器，异常F1因此不能被解读为模型族自身检出能力的比较。

### 工程复现测试

最终pytest为19项通过、0失败、0跳过，覆盖十炉实际XLSB读取、新版哈希、非均匀坐标、异常掩码、矩阵自由正规方程与小型稠密解一致性、参数回读、保真投影、形状保持求值及十炉图/模型/密集网格完整性。R03核心标量log近端另有100组密网格目标核对。工程测试通过不等于所有非凸拟合已达严格停止准则，二者分别报告。

## 第九章 结论与未解决问题

已实际读取十炉、保护新版附件1、独立真实网格审计/EDA、运行三类模型、150项独立遮蔽评估、六种人工异常和BEADS/截顶实验，并输出十个可评价连续模型与十张3D图。当前有限候选选择模型{f}；C的较大误差提示节点容量、噪声尺度和时间快速变化之间存在失配，并不证明张量P样条理论上不适合本题。

以下限制不能省略：

1. R01/R02/R05/R06/R07/R11/R12全文仍未取得，不能宣称全部正式原文阅读或原算法复现完成；只就已获得正文与补充官方来源作确证。
2. R03是二维简化对照，没有复现原文非凸局部/按列log稀疏/三维TV；若后续加入，需单独标为新实验。
3. 没有真实异常人工标签；注入F1不能证明所有实测异常被正确识别。整帧/列候选需结合物理/仪器记录复核。
4. 真饱和块缺少上界真值，掩码曲面不是已证实真实峰；删失敏感性不能消除不可辨识性。
5. 低秩补全保真投影已解决有效原网格点的附加窄峰，但普通双三次插值的负过冲通过顺序PCHIP避免；未观测位置的真实过程和已遮蔽真实峰仍需物理证据。部分迭代可能未满足停止准则，按表如实标记；有限候选、单组时间/λ带宽与三炉敏感性未覆盖全部参数/缺失机制。
6. 遮蔽参考含测量噪声，不代表干净真值；潜在低秩估计未显式强制非负；本次实际原网格估计为正，形状保持求值继承其非负范围，逐炉密集网格统计另有验证。不得将这一经验事实当成对任意新数据的理论保证。
7. 本轮未开展烟气预测、延迟估计或双模态网络，不据第一问结果推断后续问题性能。

因此工程产物已经形成可复现研究基线，但“12篇原文全部获得并完整复现、真实超量程精确恢复”的严格要求尚未达成，不能冒称全部验收完成。

文献来源与逐篇页码/式号见[文献索引](../../references/README.md)、[中文笔记](../../references/notes/)、[文献综述](../../reports/00_相关文献综述.md)。
'''
 (P/'reports/问题一_数据分析与模型建立及实验结果.md').write_text(doc)
 # Compact validation-facing extracts.
 (P/'reports/05_异常识别与评价.md').write_text('# 异常识别与评价\n\n'+fmt(a)+'\n\n'+fmt(inj[inj.split=='test'].groupby('kind',as_index=False)[['precision','recall','f1']].mean())+'\n\n整帧/整列仅候选，不自动删除。详细算法、局限和图见完整报告第五章。\n')
 (P/'reports/06_模型比较与敏感性.md').write_text('# 模型比较与敏感性\n\n'+fmt(agg)+'\n\n模型族由开发集选择，外层测试不参与。敏感性使用三个预定炉次完整网格，见完整报告第八章。\n')
 checks=dict(heat_count=10,outer_score_rows=len(r),finite_outer_rmse=bool(np.isfinite(r.rmse).all()),empty_outer_scores=int((r.n==0).sum()),final_model=f,surfaces=sum((P/f'figures/reconstructed/{h}_surface.png').exists() for h in HEATS),literature_fulltexts_missing=['R01','R02','R05','R06','R07','R11','R12'],strict_all_requirements_complete=False)
 (P/'results/final/acceptance_status.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
 print(json.dumps(checks,ensure_ascii=False))
if __name__=='__main__':main()
