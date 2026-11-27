"""Refresh one audit section in the existing research report, from saved evidence."""
import csv
import json
from collections import defaultdict
from pathlib import Path
import numpy as np


def update(out):
    from .delay_audit import ROOT
    read=lambda name:list(csv.DictReader((out/name).open()))
    real=json.loads((out/'real_summary.json').read_text())
    syn=json.loads((out/'synthetic_summary.json').read_text())
    inp=json.loads((out/'input_audit.json').read_text())
    models=json.loads((out/'flow_model_evidence.json').read_text())
    accepted=read('accepted_windows.csv'); signed=read('signed_comparison.csv');sens=read('window_sensitivity.csv');cv=read('flow_model_cross_validation.csv')
    def pct(x):return f'{100*x:.2f}%'
    def fmt(x):return '—' if x is None else f'{x:.3f}'
    lines=['<!-- delay-audit:start -->','## 6. 第一轮时滞专项验证（2026-11-27）','',
    '### 6.1 复现口径与证据等级','',
    '以 `cea8248` 为代码基线；十炉原始 XLSB 的 SHA256 全部匹配交付清单。第二问实际读取 Sequence、累计耗氧量、CO、CO₂ 三类工艺列，不读取第一问修复后的光谱数组。第一问光谱接口只作为后续跨模态连接依据，不应与工艺列的数据来源混淆。CO、CO₂ 按题面口径作为 m³/s 处理，不能从英文 concentration 标签推断为无量纲浓度。',
    '',
    '本轮保存的历史材料只有 JSON/CSV 摘要，没有同事运行时的历史 NPZ。因此可以直接确认的是逐炉摘要、模型参数、可靠窗口数与历史记录一致；逐点掩码、流量、窗口位置和时滞，是配置驱动的审查重放与未改动的基线代码输出核对一致。不能将其表述为已直接核对历史 NPZ。逐数组 SHA256 和比较记录保存在 `../delay_audit/baseline/results/`。',
    '',
    '审查入口显式读取 `q2_config.json`，并记录有效配置、输入指纹、代码版本、源代码指纹、依赖版本和全部输出指纹。原正式入口仍使用原有默认参数；本轮没有声称已将配置接入正式管线。实验接口在独立目录生成，原有正式 NPZ 字段和语义不变。',
    '',
    f'十炉流量积分与累计耗氧增量的终点相对误差绝对值最大为 {pct(max(abs(r["integral_relative_error"]) for r in inp))}；这支持总量一致性，但不证明局部流量及事件时间准确。',
    '', '### 6.2 为什么只有24个可靠窗口','',
    '| 气体 | 候选窗口 | 动态激励通过 | 相关阈值通过 | 最终通过 | 连续段 | 源窗口重叠组 | 完整搜索范围重叠组 |',
    '|---|---:|---:|---:|---:|---:|---:|---:|']
    for gas in ['co','co2']:
        rs=[r for r in real if r['gas']==gas]
        vals=[sum(r[k] for r in rs) for k in ['candidate_windows','accepted_windows','consecutive_groups','source_overlap_groups','search_overlap_groups']]
        strength=sum(r['check_counts']['strength']['passed'] for r in rs); corr=sum(r['check_counts']['correlation']['passed'] for r in rs)
        lines.append(f'| {gas.upper()} | {vals[0]} | {strength} | {corr} | {vals[1]} | {vals[2]} | {vals[3]} | {vals[4]} |')
    lines += ['',
    '达到相关阈值0.35的窗口全部进入最终可靠集合，峰分离和边界检查没有进一步减少这24个窗口。主要瓶颈是当前变化量表示下，局部最大相关很少达到阈值；不能据此断言真实过程没有时滞，也不能通过直接降低阈值宣称算法改善。逐项拒绝原因允许重叠，详见 `window_diagnostics.csv`；各检查独立通过数不能相加作为总样本量。',
    '',
    '峰显著性和峰宽目前只是诊断量，没有参与筛选。本报告第2、3节据此纠正；旧中文论文章节和英文片段中“峰宽参与可靠筛选”的表述尚待论文更新阶段同步，不可继续作为已实现功能引用。当前置信度是启发式评分，不是经过概率校准的置信概率。',
    '', '### 6.3 重复窗口、响应身份与气体差异','',
    '三种分组分别依据：可靠中心间距等于滑动步长；源窗口区间相交后传递合并；完整滞后搜索所用数据区间相交后传递合并。分组是数据共享诊断，不是独立性检验，也不是物理事件识别。所有24个相关观测均保守标为“事件身份不确定”，不能把7、5或4组改称独立响应次数。',
    '',
    '每个可靠窗口均有独立波形检查图（原始/清洗烟气、氧枪流量、对齐后的变化量和完整相关峰）。例如：6C02059 的 CO 与 CO₂ 都在680、682.5、685 s通过，时滞均为7.5 s，明显共享同一段数据；这也使两种气体的局部结果不能视作彼此独立的重复验证。6C02058 的 CO 在30、42.5、65 s给出35、27.5、6 s，源窗口高度重叠而估计明显改变，尚不能解释为三个唯一响应。',
    '',
    'CO可靠窗口分布在6炉，CO₂只在3炉。6C02057 的CO₂五个窗口跨越较长时间，源窗口形成两组，但完整搜索范围仍连通；相邻后段窗口流量升高而时滞保持36.5 s，并未构成普遍下降关系。上述现象反映当前估计的异质性，不能仅凭稀疏窗口认定两种气体具有确定不同的物理传输模式。',
    '', '### 6.4 双向搜索与全局对照','',
    '为避免左右边界不同造成样本变化，双向与非负搜索均在各炉相同的有效源中心上比较。该共同区域比原基线更窄，所以此处非负窗口数不同于13/11。',
    '', '| 气体 | 每种搜索的共同中心数 | 非负搜索通过 | 双向搜索通过 | 双向通过中的负时滞 |', '|---|---:|---:|---:|---:|']
    for gas in ['co','co2']:
        rs=[r for r in signed if r['gas']==gas]; pos=[r for r in rs if r['mode']=='nonnegative']; both=[r for r in rs if r['mode']=='signed'];passed=[r for r in both if r['accepted']=='True']
        lines.append(f'| {gas.upper()} | {len(pos)} | {sum(r["accepted"]=="True" for r in pos)} | {len(passed)} | {sum(float(r["best_delay_s"])<0 for r in passed)} |')
    globals_=json.loads((out/'global_correlations.json').read_text())
    lines += ['',f'全局方法使用相同去趋势变化量，在±60 s内计算逐滞后归一化相关。二十条炉次—气体序列的最大相关系数约为 {min(r["max_corr"] for r in globals_):.3f}–{max(r["max_corr"] for r in globals_):.3f}，总体均较弱，且峰方向并不一致。全局峰是描述性结果，不沿用局部阈值作显著性判定，也不能解释为整炉物理传输时间。',
    '', '| 气体/搜索 | 81点窗口通过数 | 121点窗口通过数 | 161点窗口通过数 |','|---|---:|---:|---:|']
    for gas in ['co','co2']:
        for mode in ['nonnegative','signed']:
            counts=[sum(int(r['accepted_windows']) for r in sens if r['gas']==gas and r['mode']==mode and int(r['window'])==w) for w in [81,121,161]]
            lines.append(f'| {gas.upper()}/{mode} | {counts[0]} | {counts[1]} | {counts[2]} |')
    lines += ['',
    '每种窗口长度内部，两种搜索共用有效中心；不同窗口长度之间的有效中心和样本量不同，故不能只按通过数排名。短窗口提高样本通过数也可能提高偶然匹配率，本轮不据此更改正式窗口长度。',
    '', '### 6.5 合成校准：有效条件与失效边界','',
    '固定30个种子101–130，18种情形、两种搜索，共1080次扫描。常规情形每种搜索/情形有3840个高度重叠窗口；它们不是3840个独立重复实验。逐种子结果与逐窗口结果均保存，报告误差仅条件于通过窗口，同时保留拒绝率和全体可评价窗口中误差不超过1 s的比例。',
    '',
    '正负固定时滞为±5、±15、±30、±55 s，另测0和±60 s。噪声幅度按源信号标准差缩放。缓变时滞定义在源时间轴上，并通过源到测量时间的逆映射生成响应。相邻事件场景只对距已知源事件中心5 s以内的窗口评价，防止给空白时间强加真实事件时滞。',
    '', '| 情形 | 搜索 | 通过率 | 通过窗口MAE/秒 | 全体可评价窗口中误差≤1秒 |','|---|---|---:|---:|---:|']
    wanted=['fixed_+0s','fixed_+15s','fixed_-15s','fixed_+60s','drifting_5_to_25s','broadened_theta_2s','broadened_theta_5s','adjacent_events','noise_1','unrelated']
    for name in wanted:
        for mode in ['nonnegative','signed']:
            r=next(r for r in syn if r['scenario']==name and r['mode']==mode and r['truth_definition']=='peak')
            lines.append(f'| {name} | {mode} | {pct(r["acceptance_rate"])} | {fmt(r.get("mae_s"))} | {pct(r["within_1s_rate_eligible"]) if r.get("within_1s_rate_eligible") is not None else "—"} |')
    lines += ['',
    '低噪声纯平移的±5、±15、±30、±55 s在双向搜索下均准确恢复；非负搜索准确恢复正时滞，却把负时滞中的少量偶然峰当作正值。零时滞在非负搜索中被端点规则全部排除，双向搜索则能恢复；±60 s仍处于双向搜索端点，按规则拒绝。这是边界策略的表现，不能宣称边界时滞不存在。',
    '',
    '展宽响应使用因果核 k(u) ∝ u exp(-u/θ)，截断至6θ，起点附加5 s；起点、峰值和离散响应重心分别记录。表中展宽误差对照已知核峰，不能将其当作起点恢复误差。宽核θ=5 s时两种搜索通过率均很低；缓变与相邻响应场景也明显弱于纯平移。说明现有方法适合局部近似平移且具有足够激励的信号，不能从理想平移测试推广到任意工业响应。',
    '',
    '无关信号的误接受率约为非负3.31%、双向6.77%。这是指定合成过程下的窗口级经验值，既不是普适误报率，也不能直接套用于真实数据；但足以说明扩大搜索范围会引入更多竞争峰，双向搜索不能只按覆盖率评价。峰宽与显著性的校准数据已保存，本轮未据此反复调阈值或宣称找到最优阈值。',
    '', '### 6.6 当前流量—时滞证据是否充分','',
    '| 气体 | 有可靠窗口的炉次 | 流量范围/(m³/s) | 固定模型训练RMSE/秒 | 约束反平方根训练RMSE/秒 | 无约束反平方根b |', '|---|---:|---|---:|---:|---:|']
    for r in models['pooled_training']:
        m=r['models'];lines.append(f'| {r["gas"].upper()} | {r["heats_with_windows"]} | {r["flow_min_m3_s"]:.3f}–{r["flow_max_m3_s"]:.3f} | {m["fixed"]["rmse_s"]:.3f} | {m["constrained_inverse_sqrt"]["rmse_s"]:.3f} | {m["unconstrained_inverse_sqrt"]["b_s_sqrt_m3"]:.3f} |')
    lines += ['',
    '以上只是重用现有模型的描述性拟合，不增加候选模型。CO训练误差有所降低，但CO₂无约束b为负，对应时滞随流量增大的方向；施加b≥0后退化为固定时滞。这不能证明真实物理方向相反，却说明现有观测不一致支持预设下降曲线。无约束模型使用原代码的不加权最小二乘，约束模型采用置信度权重，因此两者不是单纯“只改变约束”的公平对照。',
    '', '| 气体 | 有测试窗口的留一炉次折数 | 固定模型宏平均MAE/秒 | 约束反平方根宏平均MAE/秒 |', '|---|---:|---:|---:|']
    for gas in ['co','co2']:
        rs=[r for r in cv if r['gas']==gas and r['fixed_test_mae_s']]
        lines.append(f'| {gas.upper()} | {len(rs)} | {np.mean([float(r["fixed_test_mae_s"]) for r in rs]):.3f} | {np.mean([float(r["constrained_inverse_sqrt_test_mae_s"]) for r in rs]):.3f} |')
    lines += ['',
    '已有“按所选模型预测”的留一炉次结果已重放并与历史文件一致。上表另并列展示相同训练折中的既有固定模型与约束反平方根模型；未新增模型或独立事件拟合。宏平均只包含有可靠测试窗口的炉次，没有窗口的炉次保留缺失。训练RMSE降低不能替代跨炉泛化验证；留一炉次也不能消除同一炉内窗口重复和标签不确定性。原模型名 fixed_median 实际对应置信度加权均值，本轮保留接口名称并澄清含义。',
    '',
    '目前不能给出经过验证的独立事件数量，也不能据此进行可靠的事件级显著性或置信区间推断。阶段与流量覆盖、跨气体共享窗口见 `flow_stage_coverage.csv` 和 `accepted_windows.csv`。结论是证据不足以建立普遍的流量依赖时滞规律，而不是否定题目给出的物理背景。',
    '', '### 6.7 下一阶段决策','',
    '保留当前Hampel、PCHIP和SG基线及稀疏时滞诊断；不降低阈值、不直接换成双向搜索、不因短窗口通过更多就替换正式参数。下一阶段首先比较对事件形状更稳健的信号表示和响应展宽模型，并把SG导数时间尺度纳入敏感性分析；同步建立能表达唯一匹配、未匹配及多解的事件驱动验证。只有事件识别在合成校准和真实波形检查中成立后，才开展事件级流量模型、炉次分组验证和不确定性分析。',
    '', '### 6.8 结果入口与复现','',
    '- [当前实验目录](../delay_audit/)：固定路径重复运行，只有一套当前输出。',
    '- [运行与文件索引](../../README.md)：完整命令、分阶段命令、输出解释。',
    '- [全部窗口诊断](../delay_audit/window_diagnostics.csv)、[24个可靠窗口](../delay_audit/accepted_windows.csv)、[合成校准](../delay_audit/synthetic_summary.json)。',
    '- [图注清单](../delay_audit/figure_captions.csv)：图题与图片分离；图内使用真实Times New Roman。PNG/SVG/PDF保持一致，PDF字体嵌入检查见 `font_validation.json`；整套检查图只有一个版宽预览 `layout_preview.png`。',
    '', '<!-- delay-audit:end -->']
    report=ROOT/'problem2/final/reports/独立数学审查与修复报告.md'
    text=report.read_text()
    if '<!-- delay-audit:start -->' in text:text=text.split('<!-- delay-audit:start -->')[0].rstrip()+'\n'
    text=text.replace('增加源/目标标准差、相关阈值、竞争峰分离度、峰宽、边界和置信度；不合格窗口返回 NaN','增加源/目标标准差、相关阈值、竞争峰分离度、边界和置信度；峰宽与显著性仅记录诊断，不合格窗口返回 NaN')
    text=text.replace('只有动态标准差、最大相关、峰分离、峰宽和边界条件同时合格的窗口进入','只有动态标准差、最大相关、峰分离和边界条件同时合格的窗口进入（峰宽和显著性仅记录，不参与筛选）')
    text=text.replace('累计量之间的高秩相关主要来自共同时间趋势，不能直接解释为化学效率。','累计量之间的高秩相关可能受共同时间趋势影响，不能直接解释为化学效率；本轮未实施去趋势关系分析，尚未定量分离其来源。')
    report.write_text(text.rstrip()+'\n\n'+'\n'.join(lines)+'\n')
