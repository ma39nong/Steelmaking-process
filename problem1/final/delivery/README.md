# 第一问正式处理数据交付

## 数据包状态

已在本地生成`problem1-final-data.zip`、`data_package_manifest.json`和`SHA256SUMS`。ZIP按仓库相对路径组织，在克隆仓库根目录解压即可恢复正式数据位置。**本轮没有创建或上传GitHub Release**；下载仓库源码不会自动包含这个ZIP或大数组。

- 冻结数据来源：`3366db522887b2ca0f5f461ea6bd52ed71def1ba`对应的最终研究结果。
- 参数与实际1043点修复保持不变；本轮只编写交接说明、校验并打包现有文件。
- [清单](data_package_manifest.json)登记每个包内文件的相对路径、SHA256、字节数及来源；[SHA256SUMS](SHA256SUMS)核对ZIP和清单。
- [数据字段、后两问使用范围](../reports/数据交接与后续问题使用说明.md)是交接时应一起阅读的说明。

## 本轮交付核对

数据包大小523.40 MiB；31个载荷文件加嵌入清单，其中20个数组为十炉接口及差值。全部数组与既有artifact_manifest逐项一致，ZIP载荷CRC与SHA256检查通过。读取示例已在6C02058验证：1343×2047，477个实际修复位置、5805个疑似顶码，连续查询接口可用。本轮没有重新拟合模型、修改处理结果或重跑研究实验。

## 包内包含什么

|内容|包内位置|用途|
|---|---|---|
|十炉完整光谱接口|`problem1/final/data/*_interface.npz`|真实时间/波长、raw、reconstructed、掩码、质量及条件性估计|
|十炉差值|`problem1/final/data/*_difference.npz`|实际修改追溯|
|稳定查询接口|`problem1/v3/src/interface.py`|既有双线性函数，未修改模型|
|读取示例|`problem1/final/delivery/read_example.py`|NumPy读取及连续查询|
|交接说明及本说明|`problem1/final/reports/数据交接与后续问题使用说明.md`、本目录README|字段及离线/实时使用边界|
|逐炉、逐点统计与现有文件哈希|`problem1/final/results/`内选定CSV|1043点结果及文件验证|
|原始版本映射|`data_manifest.csv`|原始附件名称、坐标和哈希|
|数据包清单|`problem1/final/delivery/data_package_manifest.json`|包内所有实际载荷的哈希|

不打包历史V1/V2/V3/V3.1矩阵、重复中间拟合量及密集展示网格；它们在本机原目录完整保留。原始XLSB、第三方文献全文、氧气/烟气附加列也不包含在该光谱交接包。论文、图片和代码从`problem1`分支获取；使用相同原网格可直接查询生成任意展示密度。

## 数据包如何生成与检查

从仓库根目录运行；打包脚本只用Python标准库，不调用检测、补全、插值拟合或参数搜索：

```bash
.venv/bin/python problem1/final/delivery/build_data_package.py
```

脚本先将20个NPZ的SHA256与既有`artifact_manifest.csv`逐项比对，再登记载荷并写ZIP，最后检查ZIP的CRC及包内SHA256。既有NPZ已压缩，ZIP采用存储方式，避免无收益的重复压缩；输出小于Release单文件2 GiB限制。交付文件`SHA256SUMS`含ZIP和清单哈希。

单独验证下载包：

```bash
# 将ZIP、清单和SHA256SUMS下载到同一目录
sha256sum -c SHA256SUMS
```

Windows可用PowerShell `Get-FileHash .\problem1-final-data.zip -Algorithm SHA256`，与`SHA256SUMS`对应项比较。

## 同事收到后怎么用

```bash
git clone --branch problem1 https://github.com/ma39nong/Steelmaking-process.git
cd Steelmaking-process
# 替换为已下载ZIP的实际路径；不覆盖克隆目录中的已有代码/说明
unzip -n /path/to/problem1-final-data.zip
python3 -m venv .venv
.venv/bin/pip install numpy
.venv/bin/python problem1/final/delivery/read_example.py --heat 6C02058
```

读取原网格和函数只依赖NumPy；完整科研复现的依赖见`problem1/requirements.txt`。输入/输出光强均为counts。默认离线矩阵为`reconstructed`，原始矩阵为`raw`，`cleaned`在修复位置包含NaN，不能误作已经补全的矩阵。完整字段说明见[交接文档](../reports/数据交接与后续问题使用说明.md)。

## 截图中的Releases可以放数据吗

可以。**建议把处理后的ZIP作为Release附件发布**，并上传外部清单和校验文件。GitHub规定单个Release附件小于2 GiB；源代码ZIP与手工上传的数据附件是不同对象。[GitHub官方说明](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases)。Packages用于软件包注册表，本项目交接NPZ直接使用Releases。[Packages支持类型](https://docs.github.com/en/packages/working-with-a-github-packages-registry)。

发布时可采用标签`problem1-final-data-v1`，目标选择`problem1`分支，标题“问题一正式处理数据：十炉光谱与质量掩码”。标签表示数据交付批次，不重新命名已有研究版本。附件选择：

1. `problem1-final-data.zip`（本地准备的大数据包）；
2. `data_package_manifest.json`（文件级来源和SHA256）；
3. `SHA256SUMS`（完整ZIP和清单校验）。

Release说明可写：“冻结第一问自适应多证据＋秩4 Huber低秩修复结果，十炉1043点实际修复，含原始/修复光谱、坐标、质量、差值及连续查询接口。数据为离线重构；实时预测需从raw另建因果流程。原始工艺附加列从比赛附件获取。”

通过仓库右侧“Create a new release”上传上述附件即可；发布后把实际Release链接加入根README的数据入口。在没有真实发布地址前，文档不写虚假的下载链接。也可通过已登录GitHub CLI创建草稿：

```bash
gh release create problem1-final-data-v1 \
  problem1/final/delivery/problem1-final-data.zip \
  problem1/final/delivery/data_package_manifest.json \
  problem1/final/delivery/SHA256SUMS \
  --repo ma39nong/Steelmaking-process --target problem1 --draft \
  --title '问题一正式处理数据：十炉光谱与质量掩码' \
  --notes '十炉正式离线光谱接口与差值、质量掩码及SHA256；详细使用方法见分支总README。'
```

以上是发布方法，尚未执行。Git提交仅包含文档、读取/打包工具及小型清单；ZIP沿用项目忽略规则留在本地。
