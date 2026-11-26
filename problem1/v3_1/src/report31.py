from common31 import *

def table(d):
 d=d.copy()
 for c in d:
  if pd.api.types.is_float_dtype(d[c]):d[c]=d[c].map(lambda x:f'{x:.6g}' if pd.notna(x) else '—')
 return '| '+' | '.join(map(str,d.columns))+' |\n|'+'|'.join(['---']*len(d.columns))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in d.itertuples(index=False,name=None))+'\n'

def main():
 a=pd.read_csv(P/'results/raw_point_audit.csv');final=pd.read_csv(P/'results/final_summary.csv');high=pd.read_csv(P/'results/high_confidence_statistical_points.csv');dev=pd.read_csv(P/'results/selection_gates_development.csv');confirm=pd.read_csv(P/'results/selection_gates_confirmation.csv');counter=pd.read_csv(P/'results/counterfactual_repairs.csv');all_summary=[];perheat=[];strata=[]
 for split in ['development','confirmation']:
  d=pd.read_csv(P/f'results/contamination_{split}.csv');byheat=d.groupby(['heat_id','strategy']).mean(numeric_only=True).reset_index();byheat.insert(0,'split',split);perheat.append(byheat)
  summary=byheat.groupby('strategy').agg(pollution_rmse=('pollution_rmse','mean'),pollution_between_heat_SD=('pollution_rmse','std'),before_pollution_rmse=('before_pollution_rmse','mean'),unpolluted_rmse=('unpolluted_rmse','mean'),unpolluted_changed_fraction=('unpolluted_changed_fraction','mean'),overall_rmse=('overall_rmse','mean'),before_overall_rmse=('before_overall_rmse','mean'),integral_relative_mae=('integral_relative_mae','mean'),peak_p95=('peak_p95','mean'),peak_position_mae_nm=('peak_position_mae_nm','mean')).reset_index();summary.insert(0,'split',split);all_summary.append(summary)
  s=d.groupby(['kind','strength','strategy']).agg(tp=('tp','sum'),fp=('fp','sum'),fn=('fn','sum'),pollution_rmse=('pollution_rmse','mean'),overall_rmse=('overall_rmse','mean'),before_pollution_rmse=('before_pollution_rmse','mean'),before_overall_rmse=('before_overall_rmse','mean'),introduced_unsupported_peak_count=('introduced_unsupported_peak_count','sum')).reset_index();s.insert(0,'split',split);strata.append(s)
 summary=pd.concat(all_summary);summary.to_csv(P/'results/end_to_end_summary.csv',index=False);pd.concat(perheat).to_csv(P/'results/end_to_end_per_heat.csv',index=False);strata=pd.concat(strata);strata.to_csv(P/'results/end_to_end_by_kind_strength.csv',index=False)
 control=pd.read_csv(P/'results/control_confirmation.csv');controls=control.groupby('strategy').agg(peak_p95_mean=('peak_p95','mean'),peak_p95_max=('peak_p95','max'),masked_fraction=('masked_fraction','mean'),real_prominent_cells_masked=('real_prominent_cells_masked','sum'),fast_frame_cells_masked=('fast_frame_cells_masked','sum')).reset_index();line=pd.read_csv(P/'results/persistent_line_control_confirmation.csv').groupby('strategy').line_cells_masked.sum().reset_index();controls=controls.merge(line,on='strategy');controls.to_csv(P/'results/control_summary.csv',index=False)
 text='''# 数据预处理：点级异常最终验收

## 结论与正式任务边界

**本轮未通过实际点异常消除的科学验收。正式策略冻结为V3-P0，十炉实际修复点0，最终矩阵仍等于原始记录。正式结果必须称“异常诊断与保真连续重构”，不能称“已完成实际异常值消除”。**

这不是因为缺少仪器故障标签就拒绝统计异常处理。相反，建立了可实施的统计掩码、实际低秩修复和端到端污染恢复评价，并筛出有多重证据的原始孤立尖峰；但预声明规则在新区域确认中未达到部署门槛。保留未通过状态，不在确认后调低门槛或把开发收益当作确认收益。若题目硬性要求已经消除真实异常，本结果仍有不足，不能以曲面完成或工程测试通过代替这一要求。

## 1. 现有P0输出复核

十炉合计13,589帧，每帧2047通道，各炉有不同时间长度与实际波长坐标。原始记录总数27,816,683，顶码14,805。V3正式原网格保持原记录、确认无效0、全部统计候选保留，只完成了诊断、缺失模型研究与连续化。保真有助于避免潜变量制造新峰，但保真本身不等于去异常。旧P0采用普通缺失NRMSE改善5%门槛，不能充分回答污染恢复是否有价值的问题。

本轮将普通缺失补全与污染修复分开：使用注入前参考→注入污染→检测→掩码→未修改的秩4 Huber低秩拟合→非负支撑投影/观测保真→与参考比较。全谱积分和峰保持使用完整参考谱，不再用随机留出子谱指标替代整谱指标。

## 2. 固定统计证据与原始点分级

仅使用已有双向MAD、低秩残差、邻域一致性与连通区域。V2 MAD阈值6、辅助阈值10的并集为统计候选。保守强规则固定为：双轴z均≥8；对候选掩码后秩4背景计算原始点残差z≥6；去当前点的时间与波长参考同号偏离且两参考一致；超出八个轴向邻居范围至少6倍噪声下限；8连通区域≤3单元、时间/波长跨度均≤2；非边缘且距顶码至少2单元；不具备持续谱线的邻帧显著性支持；单帧修复比例不得超过2%。完整公式及参数见config/protocol.json、src/points.py。

局部MAD尺度与已有实现相同；低秩残差按通道稳健尺度计算，尺度下限50 counts。原始点与已掩码背景比较，避免先插值后只评价被抹平的残差。局部规则与阈值是本题保守启发式，不是已证明的故障概率。大块、整帧、持久通道与顶码不由本轮自动修复。

'''+table(a)+'''
合计51,444个统计候选。自动强规则提出5点；逐点核对全部五点后，将6C02108、304.0s、766.413nm的39674 counts保留：该局部峰在四个邻帧持续存在、相邻通道同时增强，不能以孤立统计核心冒充故障。自动中位数显著性比例未充分识别宽谱线，这是规则局限的直接例子。

其余4点有局部尖峰统计支持，独立清单如下。它们仍不是已确认仪器故障，强候选也不自动等于正式修复；正式行动还需通过端到端确认。

'''+table(high[['heat_id','time_s','wavelength_nm','raw','time_z','wave_z','lowrank_residual_z','temporal_reference','wave_reference','persistent_neighbor_frames','component_cells']])+'''
完整5点初始结果保留于results/audit_before_line_review；复核账本只否决、不新增，确认前锁定。不从人工试验假阳性中事后删去类似点。原始点复核只缩小掩码、不能提高确认实验的低召回，故不以复核替代规则确认。高可信、一般疑似和保留点清单互相独立，保留点用可无损展开的原始行/列区间游程清单表达；含顶码的区间另由质量位区分。未确认真实整帧错误，1,129个帧候选继续保留，整帧删除0。

## 3. 预声明协议与可见性

只比较V3-P0、保守点级、V3历史P1。补全固定秩4、Huber尺度0.003，验证最大300步、停止阈值1e-4，全部策略使用相同非负支撑投影、观测保真和双线性函数；P0无掩码时无需实际拟合。数值模型、插值和饱和代码完全不改。所有实际启动的注入验证拟合均达到原V3相对变化停止准则；该准则不保证隐藏真值恢复正确。

新的32×256原始片段位于旧V2第3/4片段之间，开发列64、确认列900/1720；与V2/V3局部片段不重叠。开发种子314159/271828、确认161803/141421。开发10片段×2种子×2类型×2幅度=80个污染案例，确认20片段×2种子×2类型×2幅度=160个案例，每例3策略。十炉均已用于历史研究，所以属于新注入确认，不是独立外部炉次验证；种子/片段重复不当作独立炉次。

上冲增加片段95%–5%范围的0.3/0.7倍（最小范围5000），突降乘以0.7/0.3，密度0.5%，注入位置远离片段边缘且相互留间隔。只选100<原始值<60000位置，不依据检测器筛掉参考异常；上冲限64000，避免混入饱和研究。参考为带噪原始记录，不是假定无噪声真值。

评价为绝对Precision/Recall/F1、污染区RMSE、未污染区RMSE/改动比例、整体RMSE前后、整谱积分/峰高/峰位及新引入无支持峰。原始已存在统计候选同样计入未污染FP，未人为排除以抬高Precision。另检查原始高显著谱线、积分最快变化帧及固定持续谱线反例。

选型规则在评分前声明：Precision≥0.85、Recall≥0.25、污染RMSE降≥20%、整体RMSE降≥10%、未污染改动≤0.1%、未污染RMSE/污染前整体RMSE≤0.1、未注入控制片段峰变化p95≤1%且掩码≤0.2%、每类每幅度整体误差不恶化、处理不新引入无支持峰。控制门槛取逐片段最坏值。合格策略优先少改未污染点；开发锁定后确认同门槛，不通过则预先约定回退P0。这套门槛不是普适异常处理标准，而是本轮在有限验证下防止误删的明确取舍。

指标澄清区分“未消除污染峰”与“修复新引入峰”，初次开发结果完整保留，未改变检测阈值或参考集合。原始4点复核仅作为原始数据人工保留否决，不用于改变任何注入评分。

## 4. 开发与确认端到端结果

开发门槛及结果：

'''+table(dev)+'''
开发选择保守策略：其Precision=1、Recall约0.590，污染RMSE从1416.72降至333.52，未污染点变化0。历史P1修复更多污染，但控制片段峰变化及部分类型的误差不满足条件。选型在确认打开前锁定。

独立确认结果：

'''+table(confirm)+'''
确认保守Precision=1并不意味着规则整体有效：Recall仅0.07094、污染RMSE仅下降3.93%，未达到0.25和20%两项门槛；因此未通过。历史P1 Precision约0.273，误删未污染点约0.277%，新引入无支持峰12个，控制片段峰p95最大变化约20.07%，亦未通过。确认各类/幅度Recall为：突降0.3为0.0281、突降0.7为0.0719、上冲0.3为0.0575、上冲0.7为0.1263；失败并非单一类型引起。完整分类型Precision/Recall/F1另存confirmation_kind_detail.csv。正式依预声明回退P0，未再搜索阈值。低召回和总体收益不足不等于已证明那4个原始强点不该修复；P0是本轮预声明策略门槛下的保守部署选择，不能扩大解释为所有统计异常处理都无价值。

以下宏误差先在炉内平均种子、片段与幅度，再十炉平均并报告炉间标准差。Precision/Recall门槛采用所有案例TP/FP/FN汇总的微平均，其样本依赖性不能消失；“总体F1”不等于真实故障识别率。

'''+table(summary)+'''
未注入真实片段与持续谱线控制：

'''+table(controls)+'''
完整逐炉及类型/幅度结果见end_to_end_per_heat.csv、end_to_end_by_kind_strength.csv，保留所有失败结果。两个确认波长区域覆盖更复杂光谱与强度范围，蓝端开发表现未外推成功，不对确认区域重新调阈值。

## 5. 四个原始强点的修复敏感性

为回答“若补全是否符合邻域结构”，用同一V3全炉掩码低秩模型对4个强点作一次固定600/1200步数值稳定性计算，使用已有逐点可见线性回退规则，不修改目标、阈值、插值或饱和。结果只为敏感性，未进入正式曲面：

'''+table(counter)+'''
投影后处于可见轴向邻居支撑范围内是数学约束结果，不是真实值正确性的独立验证。这些局部统计证据支持它们值得修复研究，但在本次策略级确认失败后，不用原始点修复后的“看起来平滑”替代独立验收。没有隐藏真实光强，不能为这4点报告真实修复RMSE。

## 6. 十炉正式处理前后与最终验收

'''+table(final[['heat_id','observations','candidate_points','high_statistical_after_review','repaired_points','repair_fraction','suspect_frames_retained','whole_frames_deleted','unrepaired_max_change','negative_nodes','negative_dense']])+'''
正式原始→最终差值全部为0；十炉曲线与热图前后重合是实际结果，未为了图表效果伪造去噪差异。完整差值矩阵、节点、密集网格及质量掩码留存；四点反事实差值另存、图中明确“非正式”。记录尺度顶码保留，潜在饱和估计完全复用V3并维持下界标记，不解释为测量真值。

工程验收与科学验收必须分开：新增8项测试与V3回归10项通过，历史文件/数据哈希核验通过；实际异常消除科学验收未通过。第一问正式函数冻结为记录尺度双线性连续曲面，能任意域内查询，但题目异常去除部分仍属于诊断而非已验证实际消除。未解决的问题明示于正式研究报告，不再增加算法规模或用工程完成掩盖科学不足。

## 7. 图表与追溯

- [持续谱线复核保留](../figures/main/持续谱线候选复核保留.png)
- [最终技术路线](../figures/main/最终预处理技术路线.png)、[污染端到端恢复](../figures/main/异常污染端到端恢复比较.png)
- [开发固定注入实例](../figures/main/固定注入案例_development.png)、[确认固定注入实例](../figures/main/固定注入案例_confirmation.png)
- [十炉最终连续光谱](../figures/main/十炉最终连续光谱总览.png)
- appendix目录包含每炉正式曲线、局部前后热图、全矩阵差值、独立三维曲面及四点非正式敏感性图。
- [运行与数据接口](../README.md)、[锁定策略](../config/FINAL_FROZEN.json)、[原V3完整研究](../../v3/reports/问题一_最终模型与实验结果.md)。
'''
 (P/'reports/数据预处理_点级异常最终验收.md').write_text(text)
 report='''# 问题一：最终模型与实验结果（V3.1冻结定稿）

## 一、问题重述

对十炉时间—波长火焰光谱进行异常审查，并建立在每炉实际观测定义域内可任意查询的连续函数。正式输出为ADC记录尺度，不等于无噪声或超量程真实光谱。

**正式模型冻结为P0分级诊断＋原V3缺失补全接口＋观测保真/非负支撑投影＋双线性连续重构。新增点级处理没有通过端到端确认，正式实际修复0点。不能将本结果称为“已完成实际异常值消除”。** 题目若要求真实异常已去除，本研究对此仍有明确不足。

## 二、问题分析

仅复制有效记录能避免潜变量伪峰，但不能代替异常去除。另一方面，瞬时真实谱峰与故障尖峰可能统计相似。本轮基于已有方法区分“候选证据”“可实施修复”“端到端确认”和“正式部署”，明确承认实际消除验收失败，不通过继续叠加算法掩盖。

## 三、模型假设

真实时间/波长坐标正确且严格递增；局部低秩仅为近似；原始记录有噪声；统计强证据不是已确认硬件故障。65535的潜在下界只在上量程截顶解释成立时有效。相邻可见支撑不是隐藏真实峰的物理上界。全部重构为离线过程，可用前后信息；不能直接用于第三问在线无泄漏评价。

## 四、符号说明

$t_i$为时间s，$\\lambda_j$为波长nm；$Y$为原始记录、$U=65535$；$L$为归一化秩≤4潜变量；$h_{0.003}$为Huber损失；$\\Omega$为可见拟合集，$D$为正式节点，$\\widehat F$为双线性函数。点级强证据集合并非实际故障真值；本轮正式修复集合为空。

## 五、数据分析与预处理

十炉13,589帧、每帧2047通道，总27,816,683记录；时间步长0.5s，真实波长非均匀；14,805顶码分布八炉。无非有限、负值或全零原始帧，并不意味着无统计异常。

本轮识别51,444个统计点候选，5个自动强候选中1个因持续谱线及邻通道共同增强被复核保留，4个孤立尖峰保留为高证据统计点。十炉1,129个帧候选继续保留，不自动删整帧/持续列。强点清单、一般疑似清单、无损保留游程清单独立保存，所有质量标记可追溯。

旧普通遮蔽NRMSE收益5%不足以决定异常处理。本轮已改为污染注入→检测→实际掩码→固定低秩修复→与注入前参考比较，含检测、污染区/未污染区误差、整谱积分、原始峰及新引入峰。开发与确认独立片段/种子、规则先声明，确认后未调参。

## 六、模型建立

复用V3目标与算法，不重新设计：

$$\\min_{\\operatorname{rank}(L)\\le4}\\sum_{(i,j)\\in\\Omega}w_{ij}h_{0.003}(L_{ij}-Y_{ij}/65535).$$

稀疏变量消元仅作为Huber等价解释，不充当故障标签或新精度收益。有效记录与顶码观测精确保留，待补点投影到最近可见轴向支撑非负范围。正式P0没有确认缺失，原网格 $D=Y$；这只是记录保真，不证明原始无噪声。

任意原始单元内 $a=(t-t_i)/(t_{i+1}-t_i)$，$b=(\\lambda-\\lambda_j)/(\\lambda_{j+1}-\\lambda_j)$：

$$\\widehat F=(1-a)(1-b)D_{ij}+a(1-b)D_{i+1,j}+(1-a)bD_{i,j+1}+abD_{i+1,j+1}.$$

四权非负且和为1，函数全域C0连续、非负、单元内不额外过冲；不全局可微，不保证精确恢复未采样过程。疑似顶码在正式记录曲面保留65535与不可辨识标记；局部凸核范数条件估计单独沿用V3，不宣称真实峰高或上界。

## 七、模型求解

Huber秩投影梯度、回溯、600/1200预算及可见线性数值回退完全复用V3。点级验证三策略固定300步，不重新搜索秩、损失尺度、插值或饱和参数。原V3四炉稳定性及6C02058缺失拟合未收敛限制继续有效，P0不作新拟合并不表示这些限制消失。数学推导与已核实文献仍见原V3文档，启发式点规则与原始论文方法明确分开。

## 八、重构结果与分析

本轮确认结果：

'''+table(confirm[['strategy','precision_micro','recall_micro','f1_micro','pollution_rmse','before_pollution_rmse','pollution_reduction','overall_reduction','unpolluted_changed_fraction','introduced_unsupported_peak_count','passed_added_processing']])+'''
保守处理开发污染RMSE降低76.46%，但确认仅降低3.93%、Recall 0.071，不通过预声明门槛；历史P1确认误删和新增峰风险亦不通过。正式P0对注入污染的修复收益为0，这一点不能回避。4个原始强点的反事实修复有邻域支持，但不作为通过确认的部署结果。

正式十炉原网格和密集曲面均无负值；原始点变化0，顶码误差0，实际删帧0。V3已解决的有效观测人工峰问题保持解决；未移除的原始统计候选仍然存在，不能与潜变量伪峰混为一谈。十炉前后矩阵、全零差值、曲线、局部热图、三维图及接口完整留存。

## 九、模型检验

新区域端到端污染恢复与普通缺失验证分开。真实原始谱是带噪条件参考，未用检测器排除假阳性；重复种子不作独立炉次；十炉历史已研究，未称外部验证。完整逐炉/类型/幅度及控制风险结果见[专项验收](数据预处理_点级异常最终验收.md)。原V3缺失预测配对、真实坐标降采样、饱和实验与求解器结果继续保留为历史证据，不改写为新确认收益。

工程测试及历史数据哈希核验通过；**点级实际异常消除的科学验收未通过**。最终函数只在记录尺度目标上定稿，异常去除要求的不足明确保留。

## 十、敏感性分析

有限三策略已经给出收益/误删取舍；未追加阈值搜索。5个自动强点逐点核对揭示持续宽谱线误判风险。4个保留强点按原固定低秩模型反事实补全后检查支撑、迭代与数值回退，单独绘图，不混入正式结果或以平滑外观宣称正确。旧V3参数和模型选择锁不改变。

## 十一、评价及局限

可证明的是记录保真、非负、双线性连续及单元凸包性质；条件性饱和估计下界在原局部问题中成立。不能证明统计候选是真故障、已清除所有异常、潜在超量程峰可识别、非凸全局最优或实时因果有效。保守策略区域迁移失败，宽谱线有统计误判，原始4点虽有证据仍无真实修复误差。正式P0未实施去除，须如实评价第一问异常处理部分的不足，不以工程图表完成掩盖。

## 定稿与入口

[预处理完整证据](数据预处理_点级异常最终验收.md) · [冻结策略](../config/FINAL_FROZEN.json) · [十炉曲面总览](../figures/main/十炉最终连续光谱总览.png) · [最终流程图](../figures/main/最终预处理技术路线.png) · [运行与接口](../README.md)

[保留的V3数学推导](../../v3/reports/最终数学模型推导与算法来源.md) · [保留的V3完整实验与文献入口](../../v3/reports/问题一_最终模型与实验结果.md)。本轮不增加与最终模型无关文献，不复现新算法，不开展第二/三问。
'''
 (P/'reports/问题一_最终模型与实验结果.md').write_text(report)
if __name__=='__main__':main()
