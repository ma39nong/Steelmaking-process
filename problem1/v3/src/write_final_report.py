from ctx import *
import pandas as pd

def table(d):
 d=d.copy()
 for c in d:
  if pd.api.types.is_float_dtype(d[c]):d[c]=d[c].map(lambda x:f'{x:.6g}' if pd.notna(x) else '—')
 return '| '+' | '.join(map(str,d.columns))+' |\n|'+'|'.join(['---']*len(d.columns))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in d.itertuples(index=False,name=None))+'\n'

def main():
 q=pd.read_csv(P/'results/model_confirmation.csv');q=q.groupby(['heat_id','scenario','model'],as_index=False).mean(numeric_only=True)
 compare=q.groupby(['scenario','model']).agg(RMSE=('rmse','mean'),between_heat_SD=('rmse','std'),NRMSE=('nrmse','mean'),MAE=('mae','mean'),SAM_deg=('sam_deg','mean'),peak_nm=('peak_position_nm','mean'),peak_relative=('peak_height_relative','mean'),area_relative=('area_relative','mean')).reset_index();compare.to_csv(P/'results/paired_model_summary.csv',index=False)
 cont=pd.read_csv(P/'results/continuous_confirmation.csv');cs=cont.groupby(['scenario','model']).agg(RMSE=('rmse','mean'),between_heat_SD=('rmse','std'),MAE=('mae','mean'),SAM_deg=('sam_deg','mean'),peak_nm=('peak_position_nm','mean'),peak_relative=('peak_height_relative','mean'),area_relative=('area_relative','mean'),negative_fraction=('negative_prediction_fraction','mean')).reset_index();cs.to_csv(P/'results/continuous_summary.csv',index=False)
 fs=pd.read_csv(P/'results/final_summary.csv');stats=pd.read_csv(P/'results/data_statistics.csv');sv=pd.read_csv(P/'results/singular_energy.csv');sta=pd.concat([pd.read_csv(P/f'results/{h}_stability.csv') for h in CFG['stability_heats']]);sta.to_csv(P/'results/solver_budget_comparison.csv',index=False)
 cost=pd.read_csv(P/'results/paired_compute_cost.csv').groupby('model').seconds.agg(['mean','std']).reset_index();fb=pd.read_csv(P/'results/numerical_fallback_audit.csv');cap=pd.read_csv(P/'results/fixed_artificial_cap_check.csv');tiles=pd.read_csv(P/'results/conditional_saturation_tiles.csv');peaks=pd.read_csv(P/'results/peak_origin_audit.csv');neg=pd.read_csv(P/'results/negative_origin_audit.csv')
 ratio=q[q.model=='V3-C'].nrmse.mean()/q[q.model=='linear_baseline'].nrmse.mean();delta=(q[q.model=='V1'].nrmse.mean()-q[q.model=='V3-C'].nrmse.mean())/q[q.model=='V1'].nrmse.mean()
 modules=pd.DataFrame([
 ('双向MAD','保留候选/历史基线','人工孤立点召回约48%；不确认真实故障'),('低秩补全','保留缺失模型','新随机缺失有收益，整段时间缺失近线性'),('稀疏分解','仅固定Huber等价表达','无独立精度收益证据，不将S作故障标签'),('观测保真投影','正式保留','消除有效位置潜变量替换造成的峰失真'),('饱和约束','条件估计单独保留','104局部块收敛，非全炉核范数解'),('自适应时间正则B5','研究对照','历史未显示足够增益，不扩大搜索'),('双线性插值','正式保留','C0连续、非负、单元凸包，无额外过冲'),('PCHIP','V1对照','保形但不解决错误输入网格'),('张量B样条/模型C','研究对照','V2容量与尺度研究已保留，不继续搜索'),('BEADS','独立历史对照','不将光谱基线等同异常，不进入正式流程')],columns=['模块','最终定位','依据/限制']);modules.to_csv(P/'results/module_decisions.csv',index=False)
 text='''# 问题一：最终模型与实验结果

## 一、问题重述

针对十炉转炉口火焰光谱，建立时间—波长定义域内可任意查询的连续光谱函数，识别和表达数据质量问题，对有依据的无效/缺失记录补全，保持有效记录及疑似上量程截顶信息。本研究不推断烟气、不建立实时预测网络，也不将ADC计数解释为已定标的绝对光强。

正式结果为 **P0分级预处理＋固定秩4有界影响缺失补全＋有效观测保真/局部非负支撑投影＋双线性连续重构（锁定候选V3-C）**。饱和潜在谱使用独立条件性局部凸估计。当前没有确认无效单元，因此十炉正式记录尺度原网格等于原始记录；低秩预测能力由遮蔽实验单独检验。这不构成对原始数据无噪声的断言。

## 二、问题分析

光谱有较强共同变化与窄线结构，时间突变、通道突出和整体亮度变化可能来自真实工艺，也可能来自采集污染。缺少仪器故障日志使统计候选无法唯一判真。直接以低秩潜变量替代观测会产生人工窄峰；对顶码作普通缺失又丢失量程信息。因此先按可辨识性分级，再对缺失作预测、对有效点作保真、对顶码分开表达。

主张的新增贡献是本题数据条件下的处理边界、分级验证和物理约束组合，不声称提出全新的低秩或插值算法。将候选识别、真实故障确认、补全和连续插值分别验收，可避免用复制观测的零误差掩盖缺失预测风险。

## 三、模型假设

时间与波长坐标正确且严格递增；相邻可见光谱给出局部支撑，但不保证隐藏过程无更高真峰。局部低秩是经验近似，不保证所有阶段统一低秩。65535为疑似顶码，在上量程截顶条件成立时才有潜在强度下界。原始记录有噪声，未确认统计异常保留。离线模型可以使用前后时刻，不能作为实时预测的因果预处理。

## 四、符号说明

|符号|含义|
|---|---|
|$t_i,\\lambda_j$|实际时间s及波长nm|
|$Y,U$|原始ADC记录及65535顶码|
|$\\Omega,\\Omega_o$|低秩拟合可见集、记录保留集|
|$L,S$|归一化低秩潜变量、等价稀疏残差|
|$w,r,\\eta$|拟合权重、秩4、Huber尺度0.003|
|$D,\\widehat F$|正式节点矩阵、双线性连续函数|
|$a,b$|缺失点最近可见轴向支撑范围|

## 五、数据分析与预处理

原始版本与历史代码由SHA256清单验证。十炉共13,589帧，每帧2047通道，时间间隔0.5s；波长非均匀，前四炉与后六炉仪器坐标范围略异。实际缓存无非有限、负值或全零帧。顶码共14,805点，分布八炉。下表同时记录候选帧数量；候选不是已确认采集错误。

'''+table(stats)+'''
去当前帧参考计算积分、亮度、光谱角、归一化形状、异常波长覆盖、持续长度/恢复、饱和比例和交叉拟合残差，得到1,129帧候选。通道清单亦保留实际波长、证据及建议。20%亮度及2°/5°光谱角未被当作已验证删除界限。

P0–P4采用同一低秩与插值流程对照。开发宏NRMSE为P0 0.185529、P1/P2 0.183708、P3 0.183569、P4 0.230261。额外统计处理仅约1%改善，未达到预先声明的5%收益门槛；P4还显著改变真实峰与积分。因此锁定P0，不自动删除未确认点、帧或列。P3的0.25帧权重实际进入损失，P4整帧屏蔽误差单独保存。六类人工污染的点/帧/事件指标及真实快变误删风险见[专项审核报告](数据预处理专项审核报告.md)，不能把人工异常候选F1当真实故障准确率。

## 六、模型建立

完整目标与推导见[数学推导](最终数学模型推导与算法来源.md)。核心为

$$\\min_{\\operatorname{rank}(L)\\le4}\\sum_{(i,j)\\in\\Omega}w_{ij}h_{0.003}(L_{ij}-Y_{ij}/65535).$$

消元稀疏残差得到Huber损失，采用截断SVD投影梯度与目标下降回溯。它与凸RPCA或Soft-Impute原算法存在区别，不混称复现。最终对有效观测精确保留，对待补单元约束到最近可见轴向支撑的非负区间；之后建立双线性函数。

低秩性使用原始记录每8帧子采样的SVD作描述性检查，四阶累计能量如下。它包含共同亮度、顶码和真实尖峰，能量集中不能证明噪声与工艺信号分离正确，也不单独决定秩。

'''+table(sv[sv['rank']==4][['heat_id','energy']])+'''
方法取舍如下；保留研究对照不意味着纳入正式流程。

'''+table(modules)+'''
稀疏项没有独立预测收益证据：无稀疏B对照随机恢复更好，有限候选中只有C同时满足预声明的观测、非负、局部支撑及记录顶码表达要求，因此采用C；不声称这种候选集证明所有无稀疏受限变体均更差。未在确认结果打开后新增“B加全部约束”搜索。

双线性函数在实际坐标上四角加权，连续且非负，无单元内超出四角范围的新峰；不全局可微，不保证隐藏窄峰真实高度。疑似顶码在正式记录曲面保留65535及不可辨识标记；条件性潜在估计另用带下界的局部凸核范数ADMM，绝不称超量程估计为实测值。

## 七、模型求解与参数锁定

研究协议先于分数锁定，见config/PROTOCOL_LOCKED.json；模型、预处理与核心数值源文件另有锁。复用V2秩4及尺度0.003，未重新展开参数搜索。开发与确认使用旧V2局部96行片段之间的32×256新片段；预处理与模型选型片段分离，确认种子9307/9907预先固定，精确索引见results/new_region_index.csv。原十炉此前已研究，因此这是新遮蔽确认，不是新炉次外部验证。

所有候选共享可见输入和P0策略，隐藏真值不参与检测、补全、支撑范围或选择。评分才筛选0<Y<64000，并同时报告全部未饱和参考RMSE；这对高强度和低强度带来条件性选择偏差。NRMSE分母是评分参考的95%–5%分位差。SAM只对每帧至少3个留出参考通道计算，峰位/峰高/积分只对至少5个留出参考通道的子谱计算；随机稀疏通道积分跨越间隙，是条件子谱指标，不等同完整光谱积分。正式全谱保真另以原始节点逐点、逐峰及预处理积分指标检查。两种子先在炉内平均，十炉宏平均及炉间SD报告；重复遮蔽不作为独立炉次。V1表示历史架构和锁定参数在共同新协议下重拟合，不是复用历史测试分数。新确认固定100步候选求解，与原网格扩预算研究分开报告。

### 计算成本

共同新确认32×256片段下，每候选分别计时拟合、投影及原网格查询；不含共同预处理。每炉一次、十炉平均，时间秒，受运行环境影响，不作普遍速度保证。原始确认CSV中共享潜变量的seconds仅是潜变量拟合时间，完整可比开销以下表为准。

'''+table(cost)+'''

## 八、重构结果与分析

以下为新确认缺失点配对结果，单位counts；V3-A与V2共用潜变量，缺失预测最大差严格为0。V3-C局部边界保护会略改缺失预测；复制观测不计入这些误差。

'''+table(compare)+f'''
V3-C综合宏NRMSE为{q[q.model=='V3-C'].nrmse.mean():.6f}，相对V1下降{100*delta:.2f}%，为可见线性基线的{ratio:.4f}倍，满足预声明≤1.25倍。相对V2并无精度提升（V2为{q[q.model=='V2'].nrmse.mean():.6f}）；主要改善是有效记录保真、非负和截顶语义。连续4帧缺失约7431 counts，几乎等于线性基线，不能宣称时间缺失问题已解决。

6C02058的445.5s、767.070nm位置原始记录741，V2原网格已为2861.247，该单元在旧掩码中为有效。故人工窄峰来自低秩潜变量替换有效记录，既非该点被检测器删除，也非离网格插值过冲。V3该点恢复为记录741；十炉原始有效峰相对变化均为0，仅描述保真。旧V2负值31点全部位于其补全区域，V3正式原网格与密集网格负值均为0。缺失点的支撑投影仅保证不超可见范围，无法证明所有隐藏伪峰或真峰被正确辨认。

'''+table(peaks)+'''
十炉连续函数、质量掩码和密集网格全部保存，接口允许原观测定义域内任意查询。

'''+table(fs[['heat_id','confirmed_invalid_count','recorded_saturation_count','valid_max_error','negative_count','dense_negative_count','raw_peak_max_relative_change','iterations','lowrank_converged']])+'''
八炉条件性饱和估计的104局部分块全部满足原ADMM停止标准，14,805顶码全部覆盖，下界违反0。正式记录仍保留顶码，真实潜在强度上界未知。局部块拼接不是全炉统一核范数解，估计不提供真实峰高保证。固定人工截顶检查如下（仅对原本未饱和参考，80%分位阈值，不调参）。

'''+table(cap)+'''
## 九、模型检验

真实坐标降采样独立于原网格补全评分，时间隔点对应1.0s，波长隔点使用实际nm坐标；在从未输入模型的真实坐标处调用evaluate。V2历史1.0/1.5/2.0s、更多遮蔽比例/块长保留作历史证据，不重新声明独立。新确认结果如下，峰高/积分误差相对于带噪记录。

'''+table(cs)+'''
原网格四炉预算研究如下。即使relative_change达标，未观测饱和潜变量仍可能漂移；记录尺度保真不依赖这些漂移，另用凸下界估计处理饱和。

'''+table(sta)+'''
6C02058遮蔽片段1200步仍未严格收敛且预测不稳定，不能因原网格四炉都达阈值就宣称整个求解器已完全稳定。执行预声明逐点线性回退后的检查如下。其精度损失直接报告，未据此调参或改写已锁定100步确认结果。该回退只作数值保障，不证明非凸优化全局最优。

'''+table(fb)+'''
工程测试共47项通过（V1 19、V2 18、V3 10）；V2保留1条全NaN参考变化率的历史警告，未改写历史代码。工程测试与科学验收分开：数据/协议/历史文件哈希、接口范围、非负、保真、顶码、凸下界、十炉产物完整性通过自动检查；真实故障检出率、连续时间缺失精度及6C02058缺失求解稳定性没有全部解决。正式验收对象限定为带质量标记的记录尺度连续曲面及条件性估计，不宣称潜在真实光谱完整恢复。

## 十、敏感性分析

预处理P1–P4与亮度/角度阈值敏感性均已保存；全帧统计屏蔽会显著损伤信号，未选择。非负/局部支撑投影使V2与C随机RMSE约834.9→838.0，代价较小；整时间段恢复未明显改善。无稀疏B随机RMSE约720.9，说明稀疏项不应宣传为新增贡献。四炉300/600/1200步揭示相对变化小不等于缺失稳定。新人工截顶只作固定检验，饱和参数敏感性及大规模历史实验继续引用V2，未重复搜索。

## 十一、模型评价及局限

本模型的可证明性质是有效记录保持、正式输出非负、双线性连续与单元凸包性质、条件估计在所解局部问题上的饱和下界。它不证明原始记录无噪声、统计候选是真故障、隐藏尖峰真实高度可恢复、全局可微或非凸全局最优。没有真实仪器日志和外部炉次验证；整帧/整列异常仍以诊断为主，人工异常表现不能外推真实故障准确率。长时间缺失恢复仍弱，6C02058有限迭代回退稳定但误差代价大。饱和真实上界不可辨识，局部估计和全炉潜变量均不作真实峰高宣称。

后续接口保存真实坐标、原始/重构矩阵、有效/异常/饱和与不可辨识标记；raw_until仅返回截止时刻原始记录，不附离线重构。第三问实时评价不得使用含未来观测的本轮离线结果；本轮未实现在线预测。

## 图表与复现入口

- [整体技术路线](../figures/main/第一问整体技术路线.png)
- [十炉数据特征](../figures/main/十炉数据特征统计.png)、[奇异值能量](../figures/main/奇异值能量贡献率.png)
- [配对比较](../figures/main/V1_V2_V3配对比较.png)、[谱峰保真](../figures/main/典型谱峰保真对比.png)
- [连续插值](../figures/main/连续插值精度比较.png)、[求解曲线](../figures/main/代表性求解器收敛.png)
- [十炉三维总览](../figures/main/十炉连续三维光谱总览.png)；独立十炉三维图、质量图与帧候选前后谱/时间序列/热图见figures/appendix。
- [预处理审核](数据预处理专项审核报告.md)、[分类标准](异常类型与判定标准.md)、[数学推导](最终数学模型推导与算法来源.md)、[文献边界](../references/README.md)。
- 数据数组位于本地data/，哈希、CSV、锁文件、图表与复现源代码版本化；data/沿用项目大数组不入Git约定。完整命令及接口说明见[README](../README.md)。
'''
 (P/'reports/问题一_最终模型与实验结果.md').write_text(text)
 dump(P/'results/SCIENTIFIC_ACCEPTANCE.json',dict(formal_target='recorded-scale quality-marked offline continuous spectrum',physical_integrity_pass=True,accuracy_ratio_to_linear=float(ratio),accuracy_gate_pass=bool(ratio<=1.25),confirmed_real_faults=0,all_faults_removed_claim=False,raw_grid_equal_to_records=True,conditional_local_tiles=int(len(tiles)),conditional_tiles_converged=bool(tiles.converged.all()),conditional_lower_violations=int(tiles.lower_violations.sum()),unresolved=['no external heat/fault-log validation','continuous temporal gaps weak','6C02058 held fit nonconverged at 1200; conservative fallback loses accuracy','latent true saturated peak and upper bound unidentified'],overall='restricted recorded-target acceptance; scientific limitations retained, not universal latent-recovery pass'))
if __name__=='__main__':main()
