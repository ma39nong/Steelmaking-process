from core import *

def table(d):
 d=d.copy()
 if 'band' in d:d['band']=d['band'].replace({'blue':'蓝段','middle':'中段','red':'红段'})
 for c in d:
  if pd.api.types.is_float_dtype(d[c]):d[c]=d[c].map(lambda x:f'{x:.6g}' if pd.notna(x) else '—')
 labels={'heat_id':'炉次','frames':'帧数','channels':'通道数','time_max':'末时刻/s','wavelength_min':'起始波长/nm','wavelength_max':'末波长/nm','sat':'顶码数','energy':'前4阶能量占比','strategy':'处理策略','precision':'精确率','recall':'召回率','f1':'F1','pollution_rmse':'污染区RMSE','between_heat_pollution_SD':'炉间标准差','overall_rmse':'整体RMSE','unpolluted_rmse':'未污染RMSE','unpolluted_changed_fraction':'未污染改动比例','introduced_unsupported_peak_count':'新增峰数','repair_points':'实际修复点','actually_changed':'改动点数','repair_fraction':'修复比例','upward_repairs':'向上修复','downward_repairs':'向下修复','integral_relative_mae':'积分相对误差','peak_max_relative_change':'最大原始峰变化比例','numerical_fallback_points':'数值回退点数','time_s':'时间/s','wavelength_nm':'波长/nm','original':'原始记录/counts','repaired':'修复估计/counts','difference':'估计−原始/counts','numerical_fallback':'是否回退','band':'波长区域','pollution_RMSE':'污染区RMSE','before_pollution_RMSE':'污染输入RMSE','all_RMSE':'整体RMSE','before_all_RMSE':'污染输入整体RMSE','unpolluted_RMSE':'未污染RMSE','peak_mean_relative':'平均峰高相对误差','introduced_peaks':'新增峰数','converged_cases':'达到停止准则案例','cases':'案例数'}
 d=d.rename(columns=labels)
 return '| '+' | '.join(map(str,d.columns))+' |\n|'+'|'.join(['---']*len(d.columns))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in d.itertuples(index=False,name=None))+'\n'

def main():
 f=pd.read_csv(P/'results/final_summary.csv');cmp=pd.read_csv(P/'results/unified_comparison.csv');test=cmp[cmp.split=='confirmation'].copy();hist=pd.read_csv(P/'results/V31_failure_summary.csv');source=pd.read_csv(V3/'results/data_statistics.csv');energy=pd.read_csv(V3/'results/singular_energy.csv');ev=pd.read_csv(P/'results/statistical_repairs_evidence.csv');rep=pd.read_csv(P/'results/actual_repairs.csv');review=pd.read_csv(V31/'results/high_confidence_statistical_points.csv');four=rep.merge(review[['heat_id','row','column']],on=['heat_id','row','column']);four.to_csv(P/'results/four_reviewed_actual_repairs.csv',index=False)
 exact=pd.read_csv(P/'results/exact_final_pipeline_confirmation.csv');x=exact.groupby('heat_id').mean(numeric_only=True).mean();er=1-x.pollution_rmse/x.before_pollution_rmse;ar=1-x.overall_rmse/x.before_overall_rmse;primary=test[test.strategy=='D_adaptive4'].iloc[0];tiles=pd.read_csv(P/'results/conditional_tiles.csv');bands=pd.read_csv(P/'results/confirmation_per_heat_band.csv').groupby(['band','strategy'])[['pollution_rmse','unpolluted_rmse','f1','integral_relative_mae']].mean().reset_index();bands.to_csv(P/'results/confirmation_band_summary.csv',index=False)
 fail=hist.set_index('split');total=int(f.repair_points.sum());control=pd.read_csv(P/'results/persistent_line_controls.csv').groupby('strategy').line_cells_masked.sum();fast=pd.read_csv(P/'results/control_confirmation.csv').groupby('strategy')[['fast_frame_cells_modified','peak_mean_relative','integral_relative_mae']].mean();fast.to_csv(P/'results/clean_control_summary.csv');types=ev.groupby(['heat_id','repair_type']).size().unstack(fill_value=0);types.to_csv(P/'results/repair_types.csv')
 name={'A_P0':'A 保留原始记录','B_V1_MAD':'B 双向MAD','C_V31':'C 严格多条件','D_adaptive4':'D 自适应评分（阈值4）','D_adaptive6':'D 阈值6敏感性'}
 for df in [test,bands]:df['strategy']=df.strategy.map(name)
 text=rf'''# 基于自适应多证据与鲁棒低秩补全的转炉火焰光谱异常修复及连续重构

## 摘要

针对转炉火焰光谱的局部尖峰、异常突降及量程截顶问题，建立“局部尺度校准—多证据识别—鲁棒低秩修复—连续曲面重构”模型。十炉数据共13,589帧、27,816,683个记录，具有明显的时间—波长相关性。以实际坐标上的去当前点线性预测残差衡量时间和波长异常，融合局部低秩残差、邻域一致性和谱线持续性，形成无量纲综合评分。经共同协议开发与新区域确认，选择阈值4的评分方案，实际修复{total:,}个统计异常点，占全部记录的{100*total/int(f.observations.sum()):.5f}%。修复阶段采用固定秩4 Huber低秩补全，并通过观测保真与局部非负支撑投影保护其余记录；随后构造双线性连续光谱函数。发布配置在固定确认案例中的污染区RMSE由{x.before_pollution_rmse:.2f}降至{x.pollution_rmse:.2f} counts，降低{100*er:.2f}%；全部观测RMSE降低{100*ar:.2f}%。全部十炉生成连续三维光谱，节点和密集网格无负值；14,805个疑似截顶记录保留量程标记并另存条件性估计。模型实现了有依据的实际点修复，同时保留宏观亮度变化及主要持续谱线结构。

**关键词：** 火焰光谱；异常修复；局部尺度校准；Huber损失；低秩矩阵补全；双线性插值

## 1 问题重述与分析

转炉冶炼过程中，火焰辐射强度随时间和波长共同变化。题目提供十炉离散光谱，需要分析其统计特征、合理处理异常观测，并建立定义域内可连续查询的光谱曲面。本问只使用光谱及其真实坐标；氧气和烟气数据不进入拟合。

直接沿波长平滑容易削弱窄谱线，直接沿时间平滑又可能抹去真实冶炼突变。单一全局阈值难以同时适应低强度蓝端与高波动区域。因此，首先利用局部变化尺度区分异常程度，再结合双向证据与低秩背景判断局部不一致性。缺失点由全矩阵相关结构补全，其余观测按原值保留，避免用低秩近似覆盖真实光谱细节。连续函数采用不在网格内部额外过冲的双线性形式。

![图1 第一问技术路线](../figures/main/第一问技术路线.png)

异常识别、异常修复和连续插值分别评价：检测指标衡量候选与注入位置的对应关系；污染区误差衡量修复效果；未污染区和峰/积分指标衡量代价；连续性与非负性由函数形式及数值完整性共同检验。

## 2 模型假设与符号说明

采用以下建模假设：同炉时间和波长表头正确且严格递增；局部邻域可用于估计正常连续变化；光谱具有可用的近似低秩结构，但保留局部细节；孤立统计异常可以按这些相关关系修复，不要求先确定硬件故障原因；65535按疑似上量程顶码处理；本问为离线重构，允许使用前后时刻。

|符号|含义|单位|
|---|---|---|
|$t_i,\lambda_j$|实际时间与波长|s，nm|
|$Y_{{ij}},U$|原始记录、量程顶码65535|counts|
|$R^t,R^\lambda$|去当前点的双向预测参考|counts|
|$s^t,s^\lambda,s^L$|时间、波长、低秩残差局部尺度|counts|
|$z_t,z_\lambda,z_L,S$|标准化证据与综合评分|无量纲|
|$\Omega,\Omega_{{sat}},\mathcal A$|拟合可见集、疑似顶码集、统计修复集|—|
|$L,D,\widehat F$|归一化低秩潜变量、修复节点、连续函数|无量纲/ counts|
|$r,\eta$|秩4、Huber尺度0.003|—|

记录强度使用ADC counts，不改写成未经定标的辐射功率。修复后的数值是模型估计，原始测量值和全部处理掩码均独立保存。

## 3 数据统计特征分析

### 3.1 坐标、规模与量程

十炉均有2047个光谱通道，采样间隔0.5s，时间长度不同。前四炉波长为340.535–1027.963 nm，后六炉为339.390–1023.303 nm，且通道间隔非均匀。表1按真实表头统计，不把通道序号当作波长。

'''+table(source[['heat_id','frames','channels','time_max','wavelength_min','wavelength_max','sat']])+r'''
原始记录没有NaN、Inf、负值或全零帧；这一结构检查之后，仍需识别局部统计异常。八炉存在65535重复顶码，共14,805点。转炉光谱文献描述的0.5s采样及USB4000仪器背景为量程分析提供参考[1,2]，本题坐标范围和实际通道数以附件为准。

### 3.2 非平稳变化与区域尺度差异

整帧积分随冶炼过程显著变化，谱线和宽带背景具有不同尺度。历史失败复核发现，蓝端开发片段的强度95%–5%分位范围平均约203.9 counts，而另一组确认片段约19,848.7 counts；注入幅度增长约5.3倍，时间异常z中位数却由35.34降至3.06。这说明背景变化尺度增长快于所注入异常的可辨识对比度。

旧严格多条件规则中，时间z≥8的通过比例由68.69%降至18.81%，双向参考一致性由97.00%降至37.13%；候选连通单元较少的比例也明显下降。低秩残差的通过比例反而没有下降，问题主要来自区域尺度不匹配与条件交集的累积漏检。以上为历史片段一组固定种子的失败归因；方法优劣另用本轮共同协议比较。

## 4 火焰光谱异常识别与数据预处理

### 4.1 双向局部预测及尺度校准

对内部节点，用相邻时刻的真实距离构造去当前点预测：

$$R^t_{ij}=(1-a_i)Y_{i-1,j}+a_iY_{i+1,j},\quad a_i=\frac{t_i-t_{i-1}}{t_{i+1}-t_{i-1}}.$$

波长方向同理，$R^\lambda_{ij}=(1-b_j)Y_{i,j-1}+b_jY_{i,j+1}$，其中 $b_j=(\lambda_j-\lambda_{j-1})/(\lambda_{j+1}-\lambda_{j-1})$。设 $e^t=Y-R^t$、$e^\lambda=Y-R^\lambda$。局部预测残差近似居中，用偏移集合 $K=\{-4,-3,-2,2,3,4\}$ 校准尺度：

$$s^t_{ij}=\max\left(25,1.4826\operatorname{med}_{k\in K}|e^t_{i+k,j}|\right),\quad z_t=|e^t|/s^t,$$

波长方向 $s^\lambda,z_\lambda$ 同构。跳过紧邻偏移可减少孤立污染对其邻接残差尺度的直接影响；下限25 counts避免平坦、量化邻域产生无穷z。该尺度随局部波动变化，使不同区域的得分具有可比较含义。

### 4.2 低秩残差、持续性与综合证据

先以双向MAD掩码得到初步可见数据，经双向填充及秩4投影形成背景 $B$。对原始残差 $E=Y-B$ 在9×9邻域取中位中心及MAD尺度，下限50 counts，得到 $z_L$。令 $m^L=\mathcal M_9(E)$，$s^L=\max(50,1.4826\mathcal M_9|E-m^L|)$，$z_L=|E-m^L|/s^L$，其中$\mathcal M_9$是9×9中位数算子。这一步提供背景偏离证据，不把低秩背景直接作为正式光谱。

令 $p=Y-\operatorname{med}_\lambda(Y)$ 为5通道局部谱线显著性，$P_{ij}$ 为前后四个邻帧中同号且显著性不低于当前一半的帧数。持续谱线的波长证据按 $1+P$ 衰减。以

$$c_{ij}=\left[1+\frac{|R^t_{ij}-R^\lambda_{ij}|}{\max(50,0.2|R^t_{ij}|)+s^t_{ij}+s^\lambda_{ij}}\right]^{-1}$$

表达两个参考之间的一致程度，综合评分为

$$S_{ij}=c_{ij}\sqrt{z_t\left(\frac{z_\lambda}{2(1+P_{ij})}+\frac{z_L}{2}\right)}.$$

该评分不再要求全部证据分别超过严格阈值：时间偏离与波长/背景证据共同作用，参考不一致或持续谱线通过软权重减弱误判。只比较阈值4、6两种有限方案，开发后锁定4。

### 4.3 修复决策与分级

按评分形成8连通区域。自动点修复仅作用于时间跨度不超过3、波长跨度不超过7通道的局部区域；大范围帧变化或持续通道保持诊断标记。边缘4个节点和顶码邻域保留质量标记，由独立量程策略处理。

四个历史复核尖峰在本轮均被自适应评分重新识别，按局部连续性假设纳入实际修复。6C02108的304.0s、766.413nm候选因持续谱线及邻通道共同增强而保留。处理过程区分：保留观测、一般疑似、统计修复、疑似顶码四类，不将统计修复标成硬件故障确认。

整帧积分、光谱角、归一化形状、异常波长覆盖和前后恢复等诊断继续使用；原有1,129个帧级候选不自动整行删除。图中展示典型快变前后谱及积分，点修复保持整体过程变化。

## 5 光谱低秩结构分析

光谱矩阵的奇异值满足 $\sigma_1\ge\sigma_2\ge\cdots$，以

$$E_r=\frac{\sum_{k=1}^r\sigma_k^2}{\sum_k\sigma_k^2}$$

衡量前$r$阶能量。对原始记录每8帧抽样的描述性SVD表明，少数主成分解释大部分能量。秩4的逐炉贡献如下，秩来自此前开发结果，本轮保持固定，不重新搜索。

'''+table(energy[energy['rank']==4][['heat_id','energy']])+r'''
![图2 奇异值能量贡献](../figures/main/低秩能量贡献.png)

共同亮度、宽谱背景及主要谱线构成低维相关结构，为缺失点提供跨时间与波长的约束。低秩潜变量与观测保真矩阵作用不同：潜变量用于预测被排除点，正式节点投影后无需仍保持严格秩4。低秩与稀疏分离、矩阵谱正则的理论背景见[3,4]；本文使用固定秩Huber目标，不等同于原凸PCP或Soft-Impute算法。

## 6 鲁棒低秩异常修复模型

### 6.1 观测与可见集

令 $Z=Y/U$，修复集为 $\mathcal A$，拟合集 $\Omega$ 包含不属于 $\mathcal A$ 且非顶码的有限观测。疑似顶码 $\Omega_{sat}$ 从常规低秩拟合中排除，但保留其记录尺度语义。

用Huber损失减小剩余可疑观测对潜变量的影响：

$$\min_{\operatorname{rank}(L)\le4}\sum_{(i,j)\in\Omega}h_\eta(L_{ij}-Z_{ij}),\quad\eta=0.003,$$

$$h_\eta(u)=\begin{cases}u^2/2,&|u|\le\eta,\\\eta|u|-\eta^2/2,&|u|>\eta.\end{cases}$$

Huber损失的分段形式见[5]。它等价于消去辅助稀疏残差 $S$ 后的

$$\min_{\operatorname{rank}(L)\le4,S}\sum_\Omega\left[\frac12(Z-L-S)^2+\eta|S|\right],\qquad S=\operatorname{soft}(Z-L,\eta).$$

辅助项表达有界影响，不另作故障识别器。本轮没有改变原低秩目标或加入新网络模型。

### 6.2 观测保真与局部非负支撑

对保留位置直接置 $D_{ij}=Y_{ij}$。对待修复位置，从同波长最近前后时刻及同帧最近左右波长取可见支撑集合 $V_{ij}$，定义 $a_{ij}=\max(0,\min V_{ij})$、$b_{ij}=\max V_{ij}$，再取

$$D_{ij}=\operatorname{clip}(UL_{ij},a_{ij},b_{ij}),\quad(i,j)\in\mathcal A.$$

这是受限二次投影 $\arg\min_{a\le d\le b}(d-UL)^2/2$ 的解，在连续重构之前执行，确保节点非负并避免修复点无可见支撑的显著外推。它也让修复范围明确：未被选中的测量保持原值，原始噪声不会因为全局低秩近似被整体替换。

### 6.3 疑似饱和的条件性估计

正式记录曲面在顶码保留65535，并设置不可辨识标记。另在上量程截顶解释下，沿用固定局部核范数模型

$$\min_X\frac12\|P_\Omega(X-Z)\|_F^2+0.01\|X\|_*+I_{X_{sat}\ge1}(X).$$

采用96×192核心及8×16 halo的局部ADMM；本轮仅重算因修复掩码变化而受影响的分块，其余复用完全相同目标的已有解。此分支估计潜在超量程强度，不替换正式记录的量程表达。

## 7 模型求解与连续曲面重构

### 7.1 固定秩投影迭代

梯度为 $G=P_\Omega\operatorname{clip}(L-Z,-\eta,\eta)$，以可见数据双向线性填充初始化，并执行

$$L^{k+1}=\mathcal P_4(L^k-\alpha_kG^k).$$

初始步长1，回溯只接受实际目标不增的更新。小矩阵使用精确SVD，全炉使用固定种子随机SVD秩投影。原相对迭代变化阈值为 $10^{-4}$。

全炉按原600/1200步配置检查投影后待修复点的变化；单元变化超过65.535 counts时，使用仅可见数据的双向线性估计回退，并保存数值质量位。该数值保障保持原实现，未借本轮检测改进重新设计补全模型。

条件性饱和分支令 $M$ 为可见集的0/1掩码，引入 $X=H$ 和尺度对偶变量 $V$，取 $\rho=1$：

$$X^{k+1}_{ij}=\frac{M_{ij}Z_{ij}+\rho(H^k_{ij}-V^k_{ij})}{M_{ij}+\rho},\qquad X^{k+1}_{sat}\leftarrow\max(X^{k+1}_{sat},1),$$
$$H^{k+1}=\operatorname{SVT}_{0.01/\rho}(X^{k+1}+V^k),\qquad V^{k+1}=V^k+X^{k+1}-H^{k+1}.$$

归一化原始/对偶残差均低于$10^{-4}$时停止，最多300步。下界在每次优化中施加，非绘图截断。源代码及逐块状态保存在实验附录。

### 7.2 双线性连续函数

对任意实际定义域内 $(t,\lambda)$，找到相邻节点，令 $a=(t-t_i)/(t_{i+1}-t_i)$，$b=(\lambda-\lambda_j)/(\lambda_{j+1}-\lambda_j)$，定义

$$\widehat F(t,\lambda)=(1-a)(1-b)D_{ij}+a(1-b)D_{i+1,j}+(1-a)bD_{i,j+1}+abD_{i+1,j+1}.$$

四个权重非负且和为1，因此非负节点产生非负函数，单元内部不超出四角值域。相邻单元在共享边上取同一线性值，故全域连续。规则网格的线性插值实现背景见[6]；这里使用实际非均匀波长坐标，禁止定义域外外推。连续函数为分片双线性，导数可在网格线上变化。

## 8 十炉次实际结果

### 8.1 修复数量与质量统计

本模型实际处理'''+rf'''{total:,}个记录，其中向下修复{int(f.downward_repairs.sum()):,}点、向上修复{int(f.upward_repairs.sum()):,}点。表2给出每炉的处理量与整体光谱变化。修复类型按时间参考方向区分局部上冲和突降，逐点原始值、估计值、坐标、得分和数值回退均保存。

'''+table(f[['heat_id','repair_points','actually_changed','repair_fraction','upward_repairs','downward_repairs','integral_relative_mae','peak_max_relative_change','numerical_fallback_points']])+r'''
修复改变的是被定位的局部不一致点。未修复记录的最大绝对变化逐炉均为0；这一保真性质单列检查，不作为独立预测精度。原始峰高的变化包括移除异常尖峰的预期作用，不能把任何原始最大值变化都视作损伤。

![图3 十炉实际修复数量](../figures/main/十炉实际修复统计.png)

### 8.2 四个高证据尖峰的实际修复

下表展示重新核查的四点。局部时间/波长邻域支持这些点存在短时、窄范围的异常突出，修复值由同一全炉模型估计，不再只停留在候选清单或反事实图上。

'''+table(four[['heat_id','time_s','wavelength_nm','original','repaired','difference','numerical_fallback']])+r'''
![图4 四点实际修复曲线](../figures/main/四个高证据尖峰实际修复.png)

6C02108的持续谱线候选保持原值39674 counts。其余自动修复点与这四点执行相同的评分、掩码和估计原则。每炉附图给出全谱曲线、局部时间序列、前后热图及完整矩阵差值的位置，展示实际修改而非仅改变绘图外观。

### 8.3 连续三维光谱与疑似截顶表达

十炉全部生成正式节点、密集连续网格和独立三维光谱图；定义域内可任意查询。节点与密集网格的负值均为0，实际删除整帧数为0。疑似顶码记录14,805点保持量程值，条件估计单独提供。

![图5 十炉最终连续光谱](../figures/main/十炉最终连续光谱总览.png)

'''+rf'''条件性饱和估计共104个局部分块，其中{int(tiles.refitted.sum())}块因修复掩码变化而重算；全部原顶码节点的下界违反数为0。{int(tiles.converged.sum())}/104块达到原ADMM停止准则，其余分块保留有限迭代质量标记，正式记录尺度曲面不因此更改。

'''+r'''
## 9 模型检验与敏感性分析

### 9.1 共同协议及选型

比较四种点处理：A保留原始记录、B双向MAD、C严格多条件、D尺度自适应多证据。D只考察阈值4和6。三类新波长片段与不同时间阶段覆盖十炉；开发与确认在原始坐标上分离，分别用种子1229/3557和7919/104729。每阶段30片段、两种异常、两种幅度及两种子，共240个污染案例；同一污染输入供全部策略使用。

上冲和突降注入沿用同一生成器，密度0.5%，幅度0.3/0.7，不在检测器预先筛选出的“干净点”上单独计分。评分对全部未污染点计算误改，参考是注入前带噪原始谱。开发以归一化整体MSE、未污染MSE、F1、整谱峰/积分及新增峰代价的加权连续损失选型，未设置“Recall不足则一律不修复”等绝对门槛。设 $E_a,E_u,E_b$ 分别为整体修复RMSE、未污染RMSE和污染输入整体RMSE，$H,I$ 为整谱峰高及积分相对误差，选型代价为

$$J=\frac{E_a^2}{\max(E_b^2,1)}+0.5\frac{E_u^2}{\max(E_b^2,1)}+0.1(1-F_1)+0.25H+0.25I+\frac{N_{new}}{\max(N_{inj},1)}.$$

系数在评分前固定。$J$先按案例计算再在炉内、炉间平均；开发值最小且差异不超过0.01时优先误改较少者。锁定阈值4后，确认未再调整参数。

以下强度误差单位均为counts。共同300步预算策略比较中，误差先在炉内平均片段、类型、幅度和种子，再对十炉宏平均；Precision/Recall/F1由汇总TP/FP/FN计算。D阈值6列为有限敏感性，不是新增模型。

'''+table(test[['strategy','precision','recall','f1','pollution_rmse','between_heat_pollution_SD','overall_rmse','unpolluted_rmse','unpolluted_changed_fraction','introduced_unsupported_peak_count']])+r'''
![图6 共同协议端到端比较](../figures/main/统一端到端异常处理比较.png)

自适应方案在本轮共同协议下提高污染定位与修复能力，严格多条件的漏检明显，双向MAD在复杂背景中误改较多。该结论来自相同新输入，不以不同历史协议的分数直接比较软件版本。

### 9.2 发布求解配置的端到端复核

全炉发布使用原600/1200步稳定性检查及回退。对固定确认案例另执行一次完全相同的发布求解配置，所得结果如下；它与上一表属于同一批确认输入，不算第二批独立测试。

'''+table(pd.DataFrame([dict(pollution_RMSE=x.pollution_rmse,before_pollution_RMSE=x.before_pollution_rmse,all_RMSE=x.overall_rmse,before_all_RMSE=x.before_overall_rmse,unpolluted_RMSE=x.unpolluted_rmse,peak_mean_relative=x.peak_mean_relative,integral_relative_mae=x.integral_relative_mae,introduced_peaks=int(exact.introduced_unsupported_peak_count.sum()),converged_cases=int(exact.converged.sum()),cases=len(exact))]))+rf'''
发布配置污染区误差降低{100*er:.2f}%，所有观测误差降低{100*ar:.2f}%；平均相对峰高误差为{100*x.peak_mean_relative:.2f}%，全谱积分相对误差为{100*x.integral_relative_mae:.3f}%。未污染区RMSE为{x.unpolluted_rmse:.2f} counts，仍有{int(exact.introduced_unsupported_peak_count.sum())}个修复新引入的无参考包络支持峰。这些代价与检测收益一起报告，不把修复称为无代价的去噪。

### 9.3 区域稳定性、谱线与快速变化

'''+table(bands[bands.strategy.isin([name['A_P0'],name['B_V1_MAD'],name['C_V31'],name['D_adaptive4']])])+r'''
分区域表的F1为案例宏平均，与总表的TP/FP汇总微平均分别报告。

![图7 分波长区域表现](../figures/main/波长区域稳定性.png)

阈值6减少误改和新增峰，但降低召回及污染恢复，连续代价函数选择阈值4体现了两者取舍。未注入片段额外检查峰值、积分及最快变化帧，持续谱线反例检查显示：'''+rf'''自适应阈值4修改了{int(control['D_adaptive4'])}个构造谱线单元，双向MAD为{int(control['B_V1_MAD'])}个，严格多条件为{int(control['C_V31'])}个。说明软持续性保护减少了谱线误改，但没有把所有真实窄峰都识别为永不修改。

### 9.4 数值与历史证据的使用

十炉全矩阵拟合的停止状态、回溯目标及数值回退保留，正式修复中{int(f.numerical_fallback_points.sum())}点触发既有稳定性回退。原始矩阵、处理矩阵、差值与质量位可逐点核验。历史随机/时间/波长缺失和连续降采样研究作为既有模型依据单独保留，本轮未借旧测试重新调参；原122,995个局部掩码和14,805顶码处理也完整保留，不能因后续保守策略无实际修复而抹去早期工作。

## 10 模型评价

模型把局部尺度校准与时间—波长联合证据相结合，缓解了单一阈值跨区域失配和严格条件交集漏检的问题；固定低秩模型承担相关结构预测，保真/非负投影限定修改范围，双线性函数给出可查询连续结果。新确认与实际十炉修复形成了从数据现象、检测、估计到输出的完整链条。

主要局限集中在三方面。第一，十炉原始数据此前已参与研究，本轮是新区域与新注入确认，不是外部炉次验证；原始参考带噪，人工注入不能替代真实故障日志。第二，自适应评分仍可能漏掉弱异常或误改瞬态真实谱峰，局部支撑约束也会压低缺失区域更高的真峰，已通过未污染误差、谱线反例和新增峰统计量化。第三，固定秩非凸迭代与局部分块饱和估计采用有限预算，数值质量标记需要随结果使用；上量程真实峰高及上界并非可由记录唯一识别。双线性曲面连续但不全局可微，也不保证精确恢复未采样的真实过程。

后续数据接口保存原始/处理矩阵、真实坐标、统计修复掩码、疑似顶码及不可辨识标记。第一问为离线重构，第三问实时评价应从原始历史切片建立因果预处理，不能直接把本轮含未来信息的离线结果作为在线输入。

## 参考文献

[1] Zhao B, Zhao J, Wu W, et al. Research on prediction model of converter temperature and carbon content based on spectral feature extraction. Scientific Reports, 2023, 13:14409. DOI:10.1038/s41598-023-41751-9.

[2] Ocean Optics. USB4000 Fiber Optic Spectrometer Installation and Operation Manual. 公开技术手册；原文件及核读笔记见R09索引。

[3] Candès E J, Li X, Ma Y, Wright J. Robust principal component analysis? Journal of the ACM, 2011, 58(3), Article11. DOI:10.1145/1970392.1970395.

[4] Mazumder R, Hastie T, Tibshirani R. Spectral regularization algorithms for learning large incomplete matrices. Journal of Machine Learning Research, 2010, 11:2287–2322.

[5] Huber P J. Robust estimation of a location parameter. Annals of Mathematical Statistics, 1964, 35(1):73–101. DOI:10.1214/aoms/1177703732. 损失公式实现另核对SciPy官方special.huber文档，原文获取状态见文献索引。

[6] SciPy. RegularGridInterpolator documentation. 官方非均匀直角网格插值说明，核对日期2026-11-26。

## 数据与实验附录入口

[实验协议](../../config/protocol.json)、[模型冻结](../../config/FINAL_MODEL_FROZEN.json)、[逐点实际修复](../../results/actual_repairs.csv)、[逐炉逐区域确认](../../results/confirmation_per_heat_band.csv)、[完整结果与接口](../README.md)、[文献索引](../../references/README.md)。图表均由实际数据生成；三维图与全谱热图作显示降采样，完整节点、密集网格及差值数组在data/中保存。
'''
 (P/'archive/reports/问题一_中文论文初稿.md').write_text(text)
 # Research audit remains separate from the paper's causal model narrative.
 audit=r'''# 历史实验复核与最终决策依据

V1实际局部掩码122995、疑似顶码14805，不能称其未完成修复。V2完善低秩/饱和/连续插值；V3将有效记录人工峰和负值问题与潜变量估计分开；V3.1引入端到端污染恢复，但AND条件及策略门槛导致正式不修复。算法缺陷包括尺度迁移和多条件累积漏检，策略问题包括用整体绝对门槛否决少量局部合理修复。

历史协议的参考筛选、强度范围、片段、种子及求解预算不同，不直接用历史Precision/Recall/RMSE排名。历史失败归因使用旧片段一组种子，只解释机制；最终比较使用新共同协议。

'''+table(pd.read_csv(P/'results/history_counts.csv'))+'\n历史失败尺度/各条件通过比例：\n\n'+table(hist)+'\n共同新协议：\n\n'+table(cmp)+'\n所选方案按开发连续代价最小确定，确认不调参；最终包含实际1043点修复，四个历史强点均被新规则重新识别。所有原始值、估计值与质量位独立保存。\n'
 (P/'reports/历史复核与最终决策依据.md').write_text(audit)
if __name__=='__main__':main()
