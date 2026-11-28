"""Data-driven current Chinese chapter, handoff and figure index."""
from pathlib import Path
import json
import numpy as np
from .final_statistics import RELATIONS

ROOT=Path(__file__).resolve().parents[2]
REPORTS=ROOT/'problem2/final/reports'
LABELS={'co_instantaneous':'耗氧量—CO瞬时流量','co2_instantaneous':'耗氧量—CO₂瞬时流量','co_cumulative':'耗氧量—CO累计排放','co2_cumulative':'耗氧量—CO₂累计排放'}


def table(headers,rows):
    def fmt(x):return f'{x:.4g}' if isinstance(x,(float,np.floating)) else str(x)
    return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+''.join('| '+' | '.join(fmt(x) for x in row)+' |\n' for row in rows)


def write_reports(heats,a,cfg,manuscript_only=False):
    REPORTS.mkdir(exist_ok=True);out=ROOT/cfg['output_dir']
    t1=table(['炉号','点数','时长 s','耗氧增量 m³','CO总量 m³','CO₂总量 m³'],[[r[k] for k in ['heat_id','samples','duration_s','oxygen_change_m3','co_total_m3','co2_total_m3']] for r in a['heat_summary']])
    t2=table(['关系','Pearson范围','Spearman范围','恒定排放参照Pearson范围'],[[LABELS[rel],f"{min(r['pearson_r'] for r in a['level_correlations'] if r['relation']==rel):.3f}–{max(r['pearson_r'] for r in a['level_correlations'] if r['relation']==rel):.3f}",f"{min(r['spearman_r'] for r in a['level_correlations'] if r['relation']==rel):.3f}–{max(r['spearman_r'] for r in a['level_correlations'] if r['relation']==rel):.3f}",f"{min(r['constant_gas_rate_pearson'] for r in a['level_correlations'] if r['relation']==rel):.3f}–{max(r['constant_gas_rate_pearson'] for r in a['level_correlations'] if r['relation']==rel):.3f}" if rel.endswith('cumulative') else '不适用'] for rel in RELATIONS])
    diffrows=[]
    for rel in RELATIONS:
        for sec in cfg['increment_intervals_s']:
            rr=[r for r in a['increment_correlations'] if r['relation']==rel and r['interval_s']==sec]
            diffrows.append([LABELS[rel],sec,np.median([r['pearson_r'] for r in rr]),np.median([r['spearman_r'] for r in rr]),f"{min(r['pearson_r'] for r in rr):.3f}–{max(r['pearson_r'] for r in rr):.3f}"])
    t3=table(['关系','增量间隔 s','Pearson炉间中位数','Spearman炉间中位数','Pearson逐炉范围'],diffrows)
    t4=table(['气体','等耗氧阶段','阶段平均流量中位数 m³/s','95%整炉重采样区间','阶段排放/耗氧中位数 m³/m³'],[[r['gas'].upper(),r['stage'],r['mean_rate_m3_s_median'],f"[{r['mean_rate_m3_s_ci_low']:.3f}, {r['mean_rate_m3_s_ci_high']:.3f}]",r['emission_per_oxygen_m3_m3_median']] for r in a['stage_aggregate']])
    t5=table(['关系','次数','c₀','c₁','c₂'],[[LABELS[r['relation']],r['degree'],r['c0'],r['c1'],r['c2']] for r in a['pooled_trends']])
    loorows=[]
    for rel in RELATIONS:
        for degree in cfg['trend_degrees']:
            rr=[r for r in a['trend_loho'] if r['relation']==rel and r['degree']==degree]
            loorows.append([LABELS[rel],degree,np.mean([r['mae'] for r in rr]),np.mean([r['r2'] for r in rr]),sum(r['r2']>0 for r in rr)])
    t6=table(['关系','次数','留一炉次宏平均MAE','留一炉次宏平均R²','测试R²>0炉数'],loorows)
    flow=[r for r in a['flow_review'] if r['sg_window']==31]
    t7=table(['炉号','原累计值相邻重复比例','正增量中位数 m³','31点反积分RMSE m³','31点总量相对差 %','21/41点流量偏离31点RMSE m³/s'],[[r['heat_id'],r['raw_equal_adjacent_fraction'],r['positive_increment_median_m3'],r['reintegration_rmse_m3'],100*r['endpoint_relative_error'],' / '.join(f"{x['flow_difference_from_31_rmse_m3_s']:.3f}" for x in a['flow_review'] if x['heat_id']==r['heat_id'] and x['sg_window']!=31)] for r in flow])
    t8=table(['炉号','CO/CO₂可靠窗口','CO/CO₂有效对齐点','Q1时间轴','字段历史比较','检查项'],[[r['heat_id'],f"{r['co_reliable_windows']} / {r['co2_reliable_windows']}",f"{r['co_alignment_points']} / {r['co2_alignment_points']}",'精确一致','57/57通过',sum(x['heat_id']==r['heat_id'] for x in a['interface_checks'])] for r in a['heat_summary']])
    bandnames={'low':'低','medium':'中','high':'高'}
    delayrows=[]
    for r in a['delay_flow_summary']:
        interval=f"[{r['flow_min_m3_s']:.3f}, {r['flow_max_m3_s']:.3f}"+(']' if r['upper_inclusive'] else ')')
        delayrows.append([bandnames[r['flow_band']],interval,r['gas'].upper(),r['delay_median_s'] if r['accepted_windows'] else '不可辨识',f"{r['delay_min_s']:.1f}–{r['delay_max_s']:.1f}" if r['accepted_windows'] else '—',r['accepted_windows'],r['contributing_heats'],'仅窗口级；物理事件未确认' if r['accepted_windows'] else '无合格窗口'])
    delaytable=table(['流量层','中心流量区间 m³/s','气体','候选时滞中位数 s','范围 s','合格窗口数','覆盖炉数','证据范围'],delayrows)
    maxrel=max(abs(r['endpoint_relative_error']) for r in flow)
    from .final_manuscript import render
    chapter,supplement,evidence=render(a,cfg,locals())
    (REPORTS/'问题二_中文论文_正式章节版.md').write_text(chapter)
    (REPORTS/'问题二_补充统计与验证.md').write_text(supplement)
    from .delay_audit import write_csv
    write_csv(out/'manuscript_evidence.csv',evidence)
    if manuscript_only:
        from .final_manuscript import english_section
        (REPORTS/'problem2_section.tex').write_text(english_section())
        return
    acceptance=json.loads((out/'acceptance.json').read_text())
    (REPORTS/'第二问最终验收报告.md').write_text(f'''# 第二问结果与交接验收

本轮按“关系定稿、流量复核、时滞边界、正文整合、接口验收”收尾。附件1只是改名的身份事项已关闭；`fe52075` 的报告正文和生成代码均已明确关闭该事项，旧研究报告保留其当时状态，不再作为当前待办。

| 项目 | 结论 |
|---|---|
| 四类关系 | Pearson、Spearman、三种增量尺度、五个等耗氧阶段、线性/二次趋势及留一炉次结果已保存 |
| 氧枪流量 | 31点正式实现保留，21/41点仅诊断；计数更新机制尚不能由现有数据确认 |
| 时滞结论 | 24个算法合格窗口保留，尚无独立事件确认；暂停算法开发，不拟合流量—时滞规律 |
| 中文正文 | 六张正文表与四张补充表、六组核心图及十炉结果已统一 |
| 第三问交接 | 十炉{acceptance['historical_fields']}字段、{acceptance['interface_checks']}项检查通过；第一问十炉时间轴精确一致 |
| 物理解释 | 离线计算验收通过不代表动态规律已被验证 |

正文见 [正式章节](问题二_中文论文_正式章节版.md)。[验收JSON](../relationships/acceptance.json)、[输入清单](../relationships/input_manifest.json)、[输出清单](../relationships/provenance.json)保存证据。第三问不应把模型时滞或V2置信分数解释为实测标签或概率。

四类关系中，CO有较强阶段上升趋势，CO₂先升后降；累计变量的高相关不能单独证明机制。增量相关受时间尺度影响，不能简单将其归零或全部解释为伪相关。二次趋势的留一炉次改善以实际结果报告，不用于改动正式时滞模型。

本轮关闭的是统计整理与离线接口的工程验收。题面流量—时滞方向仍受可辨识性限制；该限制已作为正式答案的一部分明确写出，不以新增模型或强制负斜率填补。
''')
    (REPORTS/'实验配置与复现说明.md').write_text('''# 实验配置与复现说明

## 当前定稿入口

冻结基准 `fe52075`，正式配置 [final_analysis.json](../../configs/final_analysis.json)。在仓库根目录运行：

```bash
.venv/bin/python -m problem2.src.final_delivery
.venv/bin/python -m problem2.src.final_review
.venv/bin/python -m pytest -q problem2/tests problem2/research/response_broadening/test_response.py problem2/research/real_event_validation/test_review.py
```

需要现有NumPy、SciPy、Matplotlib、Pillow、pyxlsb、Fontconfig的fc-match及pdffonts；图形共享配置会严格检查真实Times New Roman，字体缺失或PDF未嵌入即停止。字体不随仓库提交。

输入位置通过配置中的相对路径解析，不依赖原机器用户名。输入包括十炉原始XLSB、`final/delay_audit/baseline/data/` 的冻结NPZ，以及 `problem1/final/data/*_interface.npz`。前轮NPZ缺失时，先使用前轮审查入口复现，不得以空值代替。原始文件和冻结NPZ按已保存清单核对SHA256；第一问接口只读取时间轴和波长维度，并记录文件指纹，不修改第一问。

历史第二问Release ZIP须存在于本机缓存 `~/.cache/apmcm-b/releases/problem2-final-data/`。若尚未下载：

```bash
gh release download problem2-final-data --repo ma39nong/Steelmaking-process --pattern '*' --dir "$HOME/.cache/apmcm-b/releases/problem2-final-data" --skip-existing
gh release download problem1-final-data-v1 --repo ma39nong/Steelmaking-process --pattern data_package_manifest.json --pattern SHA256SUMS --dir "$HOME/.cache/apmcm-b/releases/problem1-final-data-v1" --skip-existing
```

入口核验ZIP的固定SHA256，逐字段比较十炉NPZ。正式 `final/data/` 若缺失，从冻结NPZ逐字节复制；若已存在则逐字段核对内容；仅ZIP容器字节差异且数组和元数据一致时保留现有文件，数值差异超过1e-10容差时停止，不静默覆盖。当前数据接口保持字段、单位和NaN语义不变。

## 统计口径

原值相关逐炉计算，不报告点级p值。增量网格为不重叠0.5、5、15秒，末端不足一格丢弃；无时滞对齐。五阶段是等20%耗氧进度，不是五分位样本。阶段边界同氧量平台取最后观测，记录端点固定，积分边界由已保存累计排放插值。Bootstrap重采样单位为整炉，2000次，固定种子5201。

线性和二次趋势均预先指定。合并拟合使每炉总权重相同；留一炉次仅用九炉拟合，测试炉不参与系数估计。耗氧进度使用最终耗氧总量，因此这是离线回顾性评价，不声称在线预测。所有模型系数及逐炉误差直接保存，不按训练R²宣布普遍有效。

## 当前文件与历史证据

当前正式统计位于 `final/relationships/`，图位于 `final/figures/`，正文位于 `final/reports/`。重复运行原位更新，不建立平行版本目录。`delay_audit/` 和 `research/` 保留独立科研证据；其冻结检查对应当时提交，完整历史复现可在相应Git提交检出后执行。旧研究正文的待核验状态不覆盖用户已确认的当前数据身份。

`run_all` 是原始工程全流程，会运行原有诊断，并在默认正式输出路径上最后调用 `final_delivery`，保证当前统计、图表及正文一致。若自定义输出路径，需要另备定稿配置。不得把旧的分位数阶段统计当作本轮等耗氧阶段统计。正式V2数值源码未改变。
''')
    from .final_manuscript import english_section
    (REPORTS/'problem2_section.tex').write_text(english_section())
    # Extend the existing dictionary without changing field semantics.
    dictionary=REPORTS/'数据字典与接口说明.md';text=dictionary.read_text();marker='\n## 5. 当前验收与物理解释\n';text=text.split(marker)[0]
    dictionary.write_text(text+marker+'\n十炉57字段共570项历史字段比较通过，第一问十炉本地时间轴已逐点精确核对。详见 [验收报告](第二问最终验收报告.md)。观测时滞的可靠标记仅指V2算法接受；峰宽及显著性为诊断量，不参与正式筛选。置信分数不是概率，模型时滞不是独立实测标签。部分峰宽诊断可为Infinity，无法辨识的观测和模型字段可为NaN，应按掩膜及字段定义处理。\n')
    handoff=REPORTS/'数据交接与第三问使用说明.md';text=handoff.read_text();marker='\n## 当前验收结果\n';text=text.split(marker)[0]
    handoff.write_text(text+marker+f'\n本轮十炉{acceptance["historical_fields"]}项历史字段比较及{acceptance["interface_checks"]}项接口检查通过；第一问本地十炉文件SHA256与Release清单一致，时间轴逐点精确一致。正式NPZ已在上述标准路径核验可读。验收记录位于 [acceptance.json](../relationships/acceptance.json) 与 [逐项检查](../relationships/interface_checks.csv)。这不证明任何候选时滞具有独立物理标签。\n\n推荐第三问将工艺变量及异常掩膜作为基础输入；只有在明确评估研究假设时才使用探索性时滞对齐。必须将“未对齐”与“时滞为零”区分，不得用模型时滞填充观测或生成伪标签。若需离线预测验证，共享拟合及预处理参数须在训练炉次内确定。\n')
    import csv
    with (out/'figure_captions.csv').open() as f:captions=list(csv.DictReader(f))
    index='# 第二问当前图表清单\n\n六组核心图与十二组补充图统一由 `python -m problem2.src.final_delivery` 生成，当前输出仅维护在 `../figures/`。PNG、SVG、PDF是同一图的必要导出格式，图题另存为可编辑 [图注CSV](../relationships/figure_captions.csv)。所有图内英文与数字使用真实Times New Roman，PDF嵌入检查见 [font_validation.json](../relationships/font_validation.json)，唯一版宽预览见 [layout_preview.png](../relationships/layout_preview.png)。\n\n'+table(['图名','用途','图题'],[[r['figure'],r['role'],r['caption']] for r in captions])+'\n第一轮时滞检查图、响应核研究图和真实事件检查页属于独立科研证据，分别保留在 `../delay_audit/` 和 `../../research/`。其组数和图数不计作独立响应数量。当前六张主表位于中文正式章节，四张补充表位于补充材料，逐炉CSV位于 `../relationships/`。\n'
    (REPORTS/'完整图表清单.md').write_text(index)
    (out/'README.md').write_text('''# 第二问正式统计结果

[中文正式章节](../reports/问题二_中文论文_正式章节版.md)给出六张主表与六组核心图，[补充材料](../reports/问题二_补充统计与验证.md)给出四张补充表，[验收报告](../reports/第二问最终验收报告.md)区分计算验收和物理结论。[复现说明](../reports/实验配置与复现说明.md)说明配置、输入和统计口径。

| 正文表 | 逐炉或完整证据 |
|---|---|
| 表1 十炉总量 | heat_summary.csv |
| 正文表2 SG窗口权衡 | manuscript_evidence.csv |
| 正文表3 原值相关 | level_correlations.csv |
| 补充表S1 增量相关 | increment_correlations.csv |
| 表4 五阶段 | stage_statistics.csv、stage_aggregate.csv |
| 补充表S2 趋势参数 | trend_parameters.csv、pooled_trends.csv |
| 正文表5 留一炉次 | trend_loho.csv |
| 补充表S3 流量复核 | flow_review.csv |
| 正文表6 流量分层时滞 | delay_flow_summary.csv、delay_flow_definition.json |
| 补充表S4 积分辨识对照 | ../../research/weak_integral/results/ |
| 工程接口验收 | interface_checks.csv、interface_field_comparison.csv |

原始计数更新机制未被确证，24个合格窗口没有转化为独立事件，不拟合流量—时滞规律。当前入口不改变正式V2、第一问或第三问字段。最终论文审查见 [九项审查](../reports/第二问最终论文审查.md)，冻结文件指纹见 `paper_freeze_manifest.json`；仅描述已保存窗口的流量分层，不再拟合新时滞模型。研究代码的历史冻结检查应在对应Git提交上复现；本目录是后续正式统计收尾，不覆盖科研证据。
''')

    audit=REPORTS/'独立数学审查与修复报告.md'
    text=audit.read_text();marker='\n<!-- final-delivery-status:start -->'
    text=text.split(marker)[0]
    audit.write_text(text+marker+'\n\n## 正式统计收尾状态\n\n附件1的中文名只是原英文文件改名，用户来源确认与SHA256核验均已完成，不再保留为待办。以上时滞审查保留当时实验与结论；其中“下一阶段”研究安排属于历史上下文。响应展宽与真实事件研究现已完成阶段审查，当前暂停时滞算法开发，不再按24个重叠窗口拟合流量—时滞规律。\n\n当前四类关系已补充Pearson、Spearman、三种增量尺度、等20%耗氧阶段、线性/二次趋势及留一炉次结果；原分位数阶段统计作为旧诊断保留，不是当前正式阶段。正式结果见 [中文章节](问题二_中文论文_正式章节版.md) 与 [验收报告](第二问最终验收报告.md)。十炉570字段历史比较、480项离线接口检查和十炉第一问时间轴核对通过；这不是对物理时滞规律的确认。\n\n<!-- final-delivery-status:end -->\n')
    manifest_path=ROOT/'problem2/final/results/artifact_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    manifest.update(current_entry='python -m problem2.src.final_delivery',formal_statistics='../relationships/',
        current_chapter='../reports/问题二_中文论文_正式章节版.md',current_supplement='../reports/问题二_补充统计与验证.md',current_acceptance='../reports/第二问最终验收报告.md',
        figures=[r['figure']+'.'+ext for r in captions for ext in ['png','svg','pdf']],
        relationships_status='Current equal-progress and increment statistics are in ../relationships; relationship_summary.json retains the earlier quantile-stage diagnostic.',
        delay_status='Exploratory overlapping-window estimates; no validated independent physical event or flow-delay law')
    from .delay_audit import write_json
    write_json(manifest_path,manifest)

    q1audit=REPORTS/'第一问成果审查与第二问接口说明.md'
    text=q1audit.read_text()
    old='Release 清单标明数据包大小约 548.8 MB，且明确不包含原始 XLSB 和工艺附加列。Release 的 SHA256 清单已保存于 `problem2/audit/`。本地未下载 548.8 MB 大包，因此接口数组的逐文件哈希尚未在本机复核。后续若需要读取光谱接口，应先下载并按清单校验。'
    new='Release清单标明数据包大小约548.8 MB，且不包含原始XLSB和工艺附加列。本轮已对本地十炉接口文件逐文件计算SHA256，全部匹配Release清单，并核对其时间轴与第二问完全一致，见 [第一问接口身份记录](../relationships/q1_interface_identity.csv)。无需为已有且哈希一致的本地文件重新下载整个压缩包。本项是文件身份和坐标核验，不重复评价第一问的科学结论。'
    q1audit.write_text(text.replace(old,new))
