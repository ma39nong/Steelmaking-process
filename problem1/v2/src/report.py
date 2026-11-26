"""Generate eleven-chapter Chinese report from immutable experimental CSVs."""
import json
import numpy as np,pandas as pd
from common import *

def table(df):
 if isinstance(df,pd.Series):df=df.to_frame()
 df=df.reset_index();cols=list(df.columns)
 def fmt(v):
  if isinstance(v,(float,np.floating)):return f'{v:.5g}' if np.isfinite(v) else '—'
  return str(v).replace('|','/')
 return '| '+' | '.join(map(str,cols))+' |\n|'+'|'.join(['---']*len(cols))+'|\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
def read(name):return pd.read_csv(P/'results'/name)
def macro(df,group,col='rmse'):return df.groupby(['heat_id']+group)[col].mean().groupby(group).agg(['mean','std'])
def main():
 for marker in ['TEST_COMPLETE.json','FINAL_COMPLETE.json','FULL_CONFIRMATION_COMPLETE.json','FULL_INTERPOLATION_COMPLETE.json']:assert (P/'results'/marker).exists(),marker
 lock=json.loads((P/'config/FINAL_SELECTION_LOCKED.json').read_text());inter=json.loads((P/'config/INTERPOLATION_LOCKED.json').read_text());B=json.loads((P/'config/LOCKED.json').read_text());test=read('test_completion.csv');c=read('c_expanded_test.csv');v1=read('v1_matched_test.csv');allmodels=pd.concat([test,c,v1],ignore_index=True);inj=read('detector_test_with_raw_FPR.csv');clip=read('artificial_censoring.csv');final=read('final_summary.csv');cont=read('interpolation_test.csv');full=read('full_confirmation.csv');fullcont=read('full_interpolation.csv');capacity=read('capacity_development.csv');conv=read('convergence_budgets.csv');rank=pd.read_csv(V1/'results/eda/low_rank_energy.csv');slow=read('solver_slow_change_diagnostics.csv');sat=read('real_saturation_comparison.csv');focus=read('focus_local_convex_saturation.csv');structure=read('structure_diagnostic_summary.csv');bootstrap=json.loads((P/'results/paired_heat_bootstrap.json').read_text());source=pd.read_csv(ROOT/'data_manifest.csv')
 result=[]
 result.append('''# 问题一 V2：模型深化、独立验证与严格验收

基准分支problem1，V1提交d8817f2。全部V1代码、原始数据、结果与报告保留不变；本报告仅对应problem1/v2。数字来自本目录实际CSV，单位为原始记录强度。没有开展第二、三问，没有使用烟气/氧气预测输入，没有撰写英文论文。

## 1 问题重述与问题分析

目标是识别火焰光谱中的仪器异常，利用时间和波长双向相关性补全缺失/异常区域，并建立可在实际坐标查询的连续曲面。必须分别回答“原始离散网格上能补回多少”和“留出的真实采样坐标能预测多少”，不能把观察点的硬复制当成预测准确率。

V1审计发现：检测器开发/测试共用原始片段；检测参数选择未先隔离外层区域；参考质量与随机遮蔽几何有关；150个V1分数是网格补全而不是连续插值。详见[完整审计](00_V1模型审计.md)。V2先按几何位置隐藏观测，再计算检测、噪声尺度与时间权重。清洗、补全、插值全程不可访问隐藏值。

本轮不以算法堆叠制造贡献。实际新增内容是可审核的可见性协议、带梯度核对的条件性饱和优化、有效观测保持与独立预测分离、具有失败证据的模型对比。尚未完成近邻文献全景比较，不能据此声称这些组合具有绝对原创性。

## 2 模型假设与符号说明

Y∈R^{n×m}是单炉次原始强度，t_i、λ_j均取附件真实数值。M是训练掩码，Q是只在评分时生成的参考位置，Ω_sat是可见65535编码位置；U=65535，X=Y/U。L为潜在低秩矩阵，S为稀疏偏离，r为秩上界，η为稀疏罚，γ为单侧下界罚，β为二阶时间罚，a_i为可见观测得到并冻结的时间权重。

假设包括局部时空相关和近似低秩；不假设所有窄峰都是异常，不假设整帧亮度跳变就是仪器错误。16位ADC及重复顶码支持条件性上截顶解释，但没有仪器标定/原始曝光记录可证明全部顶码的确切机制。真实饱和区域仅提供条件性下界，峰高上界不可辨识。

## 3 光谱数据审计及统计特征分析

沿用已核实的十炉13589条时间记录、2047个光谱通道。新附件1 SHA256为8c73aaefadb9c8d8c3f83b193d9b4583dcfdfb239446c69b7eacd4120820783c。前四炉波长340.535–1027.963nm，后六炉339.390–1023.303nm；不统一成虚构等间距坐标。时间间隔0.5s；无NaN/Inf、零值、负值或重复整谱。8炉共14805个65535疑似饱和点。

主协议每炉六段互不重叠的96×192原始片段，交替分配开发/测试，分别覆盖早、中、晚和三个波长区域。开发种子137、271；测试种子811、1297。共30个开发片段与30个独立原始测试片段。时间与波长区域在此设计中部分混杂，不能据此单独解释“阶段效应”；全矩阵确认提供补充但不是新炉次外部验证。

训练随机遮蔽独立于强度。测试1%、5%、10%随机缺失，连续时间2/4/8/16点，连续波长4/8/16/32通道。每次完整隐藏几何带，包括不适合评分的值。参考筛选0<Y<64000只影响评分，同时保存0<Y<65535的更宽条件误差及保留比例。它仍偏向未饱和参考区域，原观测也含噪声和自然异常，并非人工无噪声真值。

全部局部比较先在每炉内平均阶段和种子，再报告十炉宏平均与炉间标准差；重复遮蔽不是独立炉次。任何置信解释均受十炉固定样本限制。

## 4 光谱低秩性质的数学验证

以下为V1保存的真实完整矩阵奇异值能量审计，本轮只读复用。能量集中说明存在低维结构；未中心化的首分量大量反映亮度，不能直接等同于光谱形状的低维性，故同时列中心化结果。高能量覆盖不保证窄峰、异常和整段缺失的恢复。
'''+table(rank.set_index('heat_id')))
 result.append('''V13在PDF p5明确指出“整列稀疏”可能同时低秩，分离不可辨识。整行完全缺失时，行潜因子没有观测约束；时间先验只能补充假设，无法知道被完全遮蔽的突变。低秩证据因此支持候选模型，不能证明所有缺失模式都可唯一恢复。

## 5 分类异常识别模型

保留双向MAD基线，新增低秩残差、连通支持、归一化光谱角、亮度与形状变化分解，以及持续通道候选。低秩残差联合双向局部证据的连通区域才可增加自动掩码；整帧和持续通道主要保留为辅助诊断，不因单一异常分数删除。持续窄线也可能是真实发射峰。两种诊断脚本的角色区分：detector.py是锁定的训练检测器；structure_diagnostics.py仅解释全炉可见结构，不修改任何验证掩码。

人工异常六类覆盖不同原始阶段/波长、两种强度、开发持续长度1/4与独立测试长度2/6，种子与原始时间片段均分离。开发选择受新增误报率≤0.005约束。下表为测试十炉宏平均，Precision/F1为相对未注入基线的增量评价；raw_FPR把原始基线检测也作为非注入点误报计入，是更保守的另一口径。不能仅报告扣除基线后的高Precision。
'''+table(inj.groupby(['heat_id','method','kind'])[['precision','recall','f1','false_positive_rate','raw_false_positive_rate','total_deleted_fraction']].mean().groupby(['method','kind']).mean()))
 result.append('''辅助检测的提升很小，整帧、整列和短块的低召回是明确失败结果。这里不通过删除大量整帧/整列追求召回，不能声称已建立可靠的全类型异常识别器。亮度变化与形状变化的全炉辅助统计如下；候选不等同于仪器错误标签。
'''+table(structure.set_index('heat_id'))+'''
![六类异常](../figures/detector_test.png)

## 6 饱和感知低秩光谱补全模型

完整目标、消元、梯度、回溯、ADMM子问题与参数单位见[数学推导](01_数学推导.md)。V13原始PCP为核范数+ℓ1等式模型，V14 SOFT-IMPUTE为核范数正则；本轮B系列是V1固定秩平方/鲁棒损失的任务改造，不能继承凸PCP的全局恢复结论。

单侧软罚为γ/2 Σ_{Ωsat}[1−L_ij]_+²。B5额外β||A^{1/2}D²_t L||²，a_i仅由可见相邻行变化率确定。未看到的快速变化不能参与估权。

显式固定秩分裂采用不精确非凸ADMM，返回的L必须检查下界违反；“相对残差小”不等于每个点严格满足。另实现带核范数的凸可行替代，X子问题投影下界、Z子问题SVT；每轮X都是优化变量，不是在绘图前改数值；最终输出另在有效点作观测保真投影，饱和点分数不含这些复制位置，输出网格也不再保证小核范数。其有限迭代仍需原始/对偶残差验证。

人工截顶使用已知未饱和真实数据，三个阈值分位数0.65/0.8/0.9，两种区域4×16、12×48，两种测试种子，三阶段十炉，共1800次方法比较；其中n=0的条件保留为空、不伪造零误差。下表先炉内平均再列炉间离散程度。
'''+table(clip.groupby(['heat_id','method'])[['rmse','lower_bound_violation','shortfall','converged']].mean().groupby('method').agg(['mean','std']))+'''
人工截顶各条件以cap归一化为1，因此固定η对应的强度阈值为η×cap；跨阈值比较还包含这一尺度变化，不能把所有差异归因于截顶阈值本身。凸对照也将鲁棒固定秩损失改成平方项与核范数，τ=0.01、ρ=1固定，因而其优于其他方法不能单独归因于“硬约束”。

凸核范数下界方案的优势是在本轮人工截顶条件下得到的，不能直接证明实际65535以上的峰高已知。它是预设对照，没有根据测试结果反向替换锁定最终模型。

真实饱和8炉的下界违反统计如下，数值不是对未知真实峰高的误差。三个重点炉次还对饱和最密集的局部窗口运行凸约束，明确报告局部覆盖率，不能冒充全炉凸重构。
'''+table(sat[sat.heat_id.isin(['6C02057','6C02058','6C02107'])].set_index(['heat_id','method'])[['count','lower_violation','shortfall','converged']])+table(focus.set_index('heat_id'))+'''
![人工截顶](../figures/artificial_censoring.png)

## 7 连续光谱曲面重构方法

连续验证只保留时间间隔1.0/1.5/2.0s或每2/3个真实波长采样，再在被留出的真实坐标evaluate。隐藏值不参与清洗、补全或参数选择。方法比较共享同一清洗网格，分别为线性、三次样条、PCHIP和零平滑张量积插值B样条。三次样条与后者数学上是同一函数的不同基表示，数值相等合理，不能当作两个独立方法的成功证据。

PCHIP基于邻接割线的加权调和斜率，有一维保形/C1依据；本题二维先时间后波长的组合需要独立评价，不提升为全局C2保证。线性插值在每个原始矩形单元内为双线性组合，连续但跨单元导数可不连续。

插值方法的开发宏NRMSE如下，最终锁定 '''+inter['method']+'。\n'+table(pd.Series(inter['development_macro_nrmse'],name='development_macro_nrmse'))+'''
下表为独立局部片段测试RMSE，报告十炉宏平均和炉间标准差；MAE、SAM、峰位、峰高、积分、负值、角点界过冲完整见interpolation_independent_macro.csv。时间留行时是全可用波长谱形，波长留带时仅为该带内谱形；随机稀疏Q上的峰指标也仅对应评分子集，不能作为全谱主峰误差。
'''+table(macro(cont,['scenario','method']))+'''
![插值独立测试](../figures/interpolation_test.png)

## 8 参数估计与数值求解

B候选秩4/8，η=0.003/0.01，γ=0.5/4，β=0.03/0.3；先按开发数据锁定。C先做受控诊断，再在3种容量×3种Huber尺度的九配置中按相同开发片段、场景、种子和掩码选择。容量对照先保持alpha×nm/(KtKw)这个系数差分罚的权重，避免V1显式乘数混杂；但节点间距改变还会改变其对应的连续导数强度，因此不构成严格的连续粗糙度不变对照。另补做按(kt−3)^3/(kw−3)及对称空间因子的近似积分导数尺度修正，只用于开发解释，非精确导数Gram积分，也不修改锁定模型。C的PCG/IRLS是否停止与预测误差分别记录。

开发选择表：
'''+table(pd.Series(lock['development_macro_nrmse'],name='development_macro_nrmse'))+'\n开发评分锁定的离散候选为 '+lock['model']+'，参数 '+json.dumps(B['models'][lock['model']]['params'])+'。B2无饱和/时间项，故JSON中的gamma/beta为共用候选字典的未启用字段，不代表最终模型使用了它们。B2/B3缺失网格分数相同，按预设较低复杂度打破并列，选择B2，不因后来确认性结果改选。\n\n'+'''
C受控开发诊断（只对应随机5%，各对照共享MAD阈值9的训练掩码；不同于主比较选出的阈值6辅助掩码，不能直接与主表比较，也不能概括全部场景）：
'''+table(macro(capacity,['setting']))+'''
''' + table(read('capacity_physical_scaling_development.csv').groupby(['heat_id','kt','kw']).rmse.mean().groupby(['kt','kw']).agg(['mean','std'])) + '''

增大容量在系数罚固定对照下改善误差；导数尺度修正后的结果需与上表一起看，不能把全部改善都归因为节点数量。减小时间罚和增大Huber尺度也有影响，三者不能用一个原因概括。全容量、极小惩罚、二次损失的组合反而恶化，说明无约束扩容可能不稳定。C整体开发和测试表现不足，不等于B样条方法本身无效。

四炉的同可见性求解器预算对照：
'''+table(conv.set_index(['heat_id','solver','iteration'])[['relative_change','held_rmse','converged_1e4','converged_1e5']])+'''
V1的50条历史外层曲线没有观测到目标上升；全炉保存迭代的SVD种子扰动远小于实际每步变化，不能把主要未收敛原因归咎于随机SVD抖动。Huber尾部会限制梯度，部分炉明显，但并非四炉都处于同样尾部比例；结构性缺失与低曲率慢移动仍是恢复误差的关键。以下只诊断V1完整数据保存的潜变量，与前述带遮蔽预算审计的输入不同。
'''+table(slow.set_index('heat_id').drop(columns=['warning']))+'''
从120至300步，三炉达到1e−4，6C02108仍未达到；四炉均未达到1e−5。隐藏区RMSE下降有限，不能通过放宽阈值宣布问题解决。V2回溯提供实际目标不增检查，并未产生显著更快的恢复。全部目标、相对变化、误差与耗时曲线保存于convergence_curves.csv。

## 9 模型对比与消融实验

B0双轴线性；B1低秩平方项；B2低秩+稀疏；B3加观测保真；B4加饱和软罚；B5加自适应时间罚。B2与B3相同超参数下隐藏点预测完全相同，硬复制只改变有效观测保持，绝不算独立预测提高。

局部测试共B0–B5 3960行，C 660行，同掩码V1 660行。下表每个场景先炉内平均三阶段和两种子，再取十炉平均RMSE；完整炉次/场景/种子记录在CSV。
'''+table(allmodels.groupby(['heat_id','scenario','model']).rmse.mean().groupby(['scenario','model']).mean().unstack('model'))+'''
全场景宏平均及炉间标准差：
'''+table(macro(allmodels,['model']))+'''
V1/V2对比使用相同隐藏几何、相同输入/检测器，仅V1使用历史锁定参数与原求解器。其历史开发与当前数据并非全新炉次，不能作为外部泛化验证。B2−V1的平均RMSE差为 '''+f"{bootstrap['mean_rmse_difference']:.3f}"+'，以炉为簇的10000次重采样百分位区间为 '+str(bootstrap['cluster_bootstrap_percentile95'])+'''。这是十炉条件性不确定性，不把多次遮蔽当成数百个独立样本。

各炉总体RMSE（场景等权）：
'''+table(allmodels.groupby(['heat_id','model']).rmse.mean().unstack('model'))+'''
B5没有整体改善，并在部分波长缺失场景明显恶化；开发后追加的β=0、1e−4、1e−3、0.01、0.03、0.3与均匀/自适应权重对照仅用于开发敏感性解释，没有改变最终模型，也未反向使用测试调参。不能据此声称所有时间正则都无效。

![离散补全比较](../figures/completion_test.png)

## 10 连续重构结果及模型验证

最终全炉函数为固定秩鲁棒潜变量的网格及锁定的线性插值，而不是提前指定的“鲁棒B样条”。它改变有效观测：这既可能去噪也可能损伤真实窄峰，下面的有效点RMSE是描述性保真代价而不是独立预测误差。函数使用每炉真实坐标，不做外推；三维图没有绘图阶段截顶、补饱和下界或隐藏负值。
'''+table(final.set_index('heat_id')[['model','converged','valid_observation_rmse_not_prediction','masked_fraction','saturated_count','saturated_lower_violation','grid_negative_fraction']])+'''
''' + table(read('peak_retention_descriptive.csv').set_index('heat_id')) + '''

6C02058在t=445.5s相对于原始谱的峰值变化达到228.9%，仍出现与V1纯潜变量类似的人工峰风险。原始观测不是无噪声真值，但这种大幅新结构不足以支持无条件用于正式论文的重构。B2/B3在离散隐藏点并列，却在有效观测保持和后续插值上不同；以单一开发NRMSE和低复杂度打破并列的协议没有保护这一物理要求，这是V2选择准则的明确局限。不能在看过本轮测试后静默改选B3并继续声称独立验证；应在后续版本预声明物理保真/峰值约束和新的独立检验。

![6C02058最大原始峰变化](../figures/6C02058_raw_peak_change.png)

真实饱和不可辨识掩码、异常掩码、参数、潜变量和可调用函数均在data/{heat}_*.npz；大文件仅本地保留，可从只读输入和锁定源代码重建。

为检验局部结论，全原始矩阵上进行确认性测试：随机遮蔽仅在与开发不重叠的测试时间段，时间缺失完整覆盖所有波长，波长带覆盖测试时段的三个预设区域。B0与B2均跑两种子，V1同掩码对照跑811种子；配对比较必须仅取共同种子。确认结果不重新选择模型。下表为对应模型的逐炉场景宏平均：
'''+table(full.groupby(['heat_id','scenario','model']).rmse.mean().groupby(['scenario','model']).mean().unstack('model'))+'''
''' + table(full[full.seed==811].groupby(['heat_id','scenario','model']).rmse.mean().groupby(['scenario','model']).mean().unstack('model')) + '''

上一个补充表限定共同种子811，才是V1与V2的严格配对比较；前表B0/B2包含两种子，不能逐格与仅811的V1解释成配对差值。

全矩阵共同种子811、场景等权的宏平均RMSE，V1为2750.13、B2为2364.22，下降约14.0%，主要来自波长带补全。随机1%却从223.04恶化到314.78，完整时间段的改进极小。这些相反结果与人工峰风险一起说明：整体平均误差改善不能支持全场景或物理保真的优越性。

完整矩阵连续降采样另提供250行确认数据：四种共同清洗网格方法，加“锁定B2补全+锁定线性插值”组合。该确认覆盖原始历史开发时段，不能取代严格分离的局部主测试，也不用于调参。组合方法和仅清洗网格插值的差异反映潜变量修改有效点带来的代价。
'''+table(macro(fullcont,['scenario','method']))+'''
![十炉三维图总览](../figures/all_heats_surfaces.png)

十炉三维图：
''')
 for h in HEATS:result.append(f'![{h}](../figures/{h}_surface.png)\n')
 result.append('''## 11 敏感性分析与模型局限

本轮完成工程和实验交付，但不能宣称已满足所有科学性能要求。主要未解决问题：

1. 整帧、持续整列、短块异常的识别可靠性不足；没有真实仪器异常标签，人工实验的可迁移性有限。
2. B5不稳定改善，整段时间缺失仍难；缺失内部的快速变化不可凭正则唯一识别。
3. B2在少量低亮度原始网格点输出负值，未做事后截断；有效观测修改与窄峰保持存在取舍，饱和软/硬方案虽改善人工截顶，却没有被开发整体评分选为最终模型。真实上量程部分仍可能违反条件下界，不虚构恢复峰高或上置信界。
4. 四炉最终B2拟合未达锁定停止准则；预算审计显示更多迭代改善有限。目标回溯不意味着非凸全局最优。
5. 原始参考含噪声，0<Y<64000筛选造成条件性；NRMSE按评分值q95−q05归一化，在低亮度小动态范围片段可能放大相对误差。故同时报告原始RMSE、MAE、保留比例与更宽未饱和参考。
6. 主实验是分层局部原生片段；全炉确认和最终重拟合不能制造新的独立炉次。时间与波长区域部分混杂，十炉全部来自现有数据，无法保证跨炉型/仪器泛化。
7. 线性函数连续但不全局可微；三次/张量B样条表达等价。最终选择仅覆盖预设有限候选，不能证明方法类别的全局最优。
8. 两篇SIAM PCHIP原文仍未取得合法作者/出版社全文；V13/V14/V17已获取公开全文并核读核心方法，未声称逐行复核全部证明。原R01–R12未获取全文的缺口保留。V17正确DOI为iau006，原给定iau004的错配记录保存。

新增文献和原始方法边界见[文献索引](../references/README.md)。V1与V2联合测试37项通过（68.70s）；测试覆盖源文件不变、真实十炉读取、梯度差分、掩码独立、B2/B3关系、连续坐标查询与边界、显式凸下界及结果数量。首次失败是空截顶参考条件的测试断言错误，已修复并保留首次日志；科学性能未达标的事实没有修改。代码测试通过不代表科学验收通过。

代码、CSV、锁文件、模型公式及本报告共同构成可审查结果；未来改进应启用新的实验版本及新的独立测试，不能在已打开的V2测试上继续调参再宣称独立。
''')
 text='\n\n'.join(result);(P/'reports/问题一_V2模型深化与严格验证.md').write_text(text)
 dump(P/'results/acceptance_status.json',dict(engineering_deliverables_complete=True,scientific_all_requirements_pass=False,selected_model=lock['model'],selected_interpolation=inter['method'],issues=['Spurious peak risk remains: 6C02058 at 445.5s has 228.9% raw peak change','Frame/column/block anomaly recall remains low','Four final fits fail locked change criterion','Saturated latent heights unidentifiable; final selected masked model violates conditional lower bounds','Small number of unconstrained low-light latent predictions are negative','Two original PCHIP fulltexts unavailable','No independent external heats; local panel and full-native confirmation scopes differ'],report='reports/问题一_V2模型深化与严格验证.md'))
if __name__=='__main__':main()
