# 第一问处理后数据交接与后续问题使用说明

## 正式数据选哪个

正式数据统一使用`problem1/final/data/{heat}_interface.npz`。最终方案为自适应多证据评分（阈值4）、秩4 Huber补全（尺度0.003）、观测保真、局部非负支撑、条件性饱和估计和双线性连续函数。十炉实际修复1043点；历史V3/V3.1的保留原值数据不作为当前正式修复数据。

中文[正式章节](../archive/reports/问题一_中文论文_V2_正式章节版.md)解释方法与实验；[图表清单](../archive/reports/问题一_V2_正文图表清单.md)关联图片及分析；[全文整合材料](问题一_V2_全文整合材料.md)登记公式、来源和实验协议。文稿V2表示章节重构，与算法历史V2区分。

## 文件在哪里

仓库根目录为`/home/user/APMCM-B`（其他同事以自己的克隆路径为准）。原始附件在`Data/`，正式处理结果在`problem1/final/data/`。这两个数据目录均不随普通Git克隆下载。

|文件|内容|是否为后续建模首选|
|---|---|---|
|`{heat}_interface.npz`|完整原网格、坐标、质量及查询所需参数|是，每炉一个|
|`{heat}_difference.npz`|`difference=处理矩阵-原始矩阵`及修复掩码|用于追溯|
|`{heat}_dense.npz`|更密坐标上的连续曲面展示值|仅展示，建模优先原网格接口|
|`{heat}_masks.npz`|检测分量、自动/复核掩码及分级|深入检查检测依据时使用|
|`{heat}_fit.npz`|潜变量、正式矩阵和回退标记|研究求解时使用|
|`{heat}_conditional.npz`|八个有顶码炉次的条件估计|专门研究截顶时使用|

十炉接口文件合计 522.86 MiB；整个`final/data/`还包含中间量、差值和密集网格，约1.6 GiB。交接包优先包含十炉接口和差值，无须为特征建模下载重复中间矩阵。

## 接口字段与正确用法

|字段|内容|使用说明|
|---|---|---|
|`time`|长度n的相对时间/s|由附件Sequence差乘0.5构成，各炉从0开始；跨模态对齐需核对原附件时间口径|
|`wavelength`|长度2047的实际波长/nm|非均匀，各炉坐标需分别读取|
|`raw`|n×2047原始光谱/counts|测量记录，用于复核及后续因果处理|
|`cleaned`|统计修复位置置NaN的原始矩阵|检测后的掩码矩阵，尚不是完整补全结果|
|`reconstructed`|n×2047正式异常修复矩阵/counts|离线特征和原网格分析的默认矩阵|
|`latent`|低秩潜变量UL/counts|背景诊断；不替代正式矩阵|
|`statistical_repair`、`invalid`|实际排除并修复的位置|统计修复掩码，不等同于硬件故障标签|
|`valid`|低秩拟合可见集合|排除修复和顶码；顶码在正式记录矩阵中仍保留|
|`quality_grade`|0保留、1疑似、2修复|配合独立顶码标记使用|
|`suspected_saturation`|原始65535位置|不要把它当作普通缺失或真实峰高|
|`saturated_unidentifiable`|顶码处潜在强度不可完全辨识标记|随数据传递|
|`conditional`|仅顶码位置替换为条件性估计，其余与正式矩阵相同|可作敏感性分析，超量程值是估计|
|`conditional_solver_converged`|条件估计在顶码位置的求解状态|只在顶码处有意义（见下方说明）|
|`numerical_fallback`|使用既有可见数据双向线性回退的位置|保留质量信息|
|`quality`|uint16组合位标记|见下方位定义|
|`score`|无量纲综合异常评分|诊断用途|
|`offline`、`rank`、`eta`、`interpolation`|True、4、0.003、bilinear|处理元数据|

`conditional_solver_converged`只在`suspected_saturation=True`的节点判断；无截顶炉次或非截顶位置的False不表示普通修复失败。质量位：2顶码、4补全、8局部支撑投影、16负潜变量修正、32弱支撑、256一般疑似、1024帧候选、8192条件估计有限迭代、16384统计修复、32768数值回退。可以同时出现多个标记。

## 最短读取示例

从仓库根目录运行，读取只需要NumPy：

```python
import numpy as np
from pathlib import Path
p = Path('problem1/final/data/6C02058_interface.npz')
with np.load(p, allow_pickle=False) as d:
    t = d['time']
    wavelength = d['wavelength']
    spectrum = d['reconstructed']     # 后续离线特征输入
    raw = d['raw']                    # 原始对照/因果流程来源
    repaired = d['statistical_repair']
    sat = d['suspected_saturation']
    print(spectrum.shape, int(repaired.sum()), int(sat.sum()))
```

查询连续曲面：

```python
import sys
sys.path.insert(0, 'problem1/v3/src')
from interface import Spectrum
m = Spectrum.load('problem1/final/data/6C02058_interface.npz')
values = m.evaluate([445.5, 445.75], [767.07, 770.0])
# 结果为2×2坐标笛卡尔积；时间/s、波长/nm，光强/counts
historical_raw = m.raw_until(100.0)  # 仅切片原始记录，未执行因果清洗
```

[可直接运行的读取示例](../delivery/read_example.py)提供原网格统计及连续查询；接口的域外查询会被拒绝。正式矩阵保留未修复点，含疑似顶码记录；密集网格等距坐标仅供展示，原节点使用实际非均匀波长。

## 是否要交给后两问同事

**要交接，而且应连同原始值和质量标记一起交接。** 第一问提供统一的光谱数据基准，避免各同事混用不同阶段的矩阵。第二问、第三问仍需按照各自研究目标确定输入的可见范围。

|任务|建议数据|应另行完成的工作|
|---|---|---|
|第二问离线特征提取、相关分析、回顾性拟合|`reconstructed`为主，`raw`作对照，携带修复/顶码掩码|从原始XLSB读取氧气、烟气等附加量，核对名称单位和时间对齐；按炉次或明确协议划分训练/测试|
|第二问若评价某时刻即时预测|原始历史片段及因果处理结果|预处理仅用当时可见信息，独立固定模型和特征参数|
|第三问实时预测|`raw`和真实坐标；`raw_until`可提供历史切片|重新实现因果异常识别、尺度估计和特征生成；不能直接把离线补全或离线质量掩码当作在线输入|
|论文复现、曲面展示、谱峰积分分析|原网格接口＋连续函数|按实际坐标查询和积分，记录使用记录尺度还是条件性估计|

本问离线检测同时查看前后帧，低秩估计利用整炉信息。即使仅修改1043点，离线修复掩码本身也可能包含未来信息。实时预测需重新确定异常判断及补全的因果可见范围。训练/测试划分后，标准化、降维、特征和任何由数据学习的参数均按该任务协议拟合；测试炉次仅用于规定的评价。

接口包只覆盖光谱相关内容，**不包含氧气、CO、CO₂等工艺列**。原XLSB的2051列包含Sequence、2047个波长列及其余3列；已核对6C02058表头，后三列为`Cumulative oxygen consumption`（累计耗氧）、`CO₂ concentration`和`CO concentration`。后两问读取各炉原表头并核对单位、时间与变量含义。历史`problem1/data/processed/{heat}.npz`中的`process_audit_only`仅为审计保留，未经本轮作为预测输入处理。

## 原始附件与版本核对

十炉原始XLSB按赛题官方或团队已有渠道下载，文件名与SHA256见根[data_manifest.csv](../../../data_manifest.csv)。当前仓库没有登记可验证的外部原始附件直链，因此不提供猜测的下载地址。6B03205使用已核对的替换文件，SHA256为`8c73aaefadb9c8d8c3f83b193d9b4583dcfdfb239446c69b7eacd4120820783c`。应保留原文件和原表头。

最终大数组按[artifact_manifest.csv](../results/artifact_manifest.csv)校验。数据交接包将其接口、差值与校验清单一起提供；对外下载后核对包SHA256，并逐项核对清单。处理数据是从原始记录派生的估计，原始文件不被覆盖。

## Releases和Packages如何选

使用截图中的 **Releases** 上传`problem1-final-data.zip`和`SHA256SUMS`，可供同事直接下载。单个Release附件须小于2 GiB，见[GitHub官方说明](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases)。Packages是软件包注册表，见[支持类型](https://docs.github.com/en/packages/working-with-a-github-packages-registry)，本数据无需封装成npm或容器包。

正式数据已发布：[Release下载页](https://github.com/ma39nong/Steelmaking-process/releases/tag/problem1-final-data-v1)、[直接下载ZIP](https://github.com/ma39nong/Steelmaking-process/releases/download/problem1-final-data-v1/problem1-final-data.zip)。普通`git push`只同步代码、文稿和清单，不上传忽略目录中的数组。发布操作及下载流程见[交付说明](../delivery/README.md)。


当前附件1规范来源名为`Data/Attachment 1. Heat No. 6B03205.xlsb`，重复扩展名已规范化，文件SHA256不变；原接收名称和清单快照见[来源名称更正记录](../supplement_v1/results/source_name_corrections.csv)。该变更仅涉及文件名称，不涉及数据内容或已发布接口。
