# 原始测量链审计

独立研究目录，正式基线 `d08f989`、研究基线 `3790d5f`。不修改正式V2、第一问、第三问接口或冻结论文。只审计记录更新、完整软件重构链和时滞可辨识性，不新增事件算法或流量回归。

从仓库根目录运行：

```bash
.venv/bin/python -m problem2.research.measurement_chain.run
.venv/bin/python -m problem2.research.measurement_chain.run --render-only
.venv/bin/python -m pytest -q problem2/research/measurement_chain/test_chain.py
```

参数在 `config.json`，结果在固定的 `results/`。主报告 `results/REPORT.md`，其余CSV/JSON记录实际更新节律、逐条件误差、背景注入及此前A/B的整炉Bootstrap。`--render-only`先检查已有配置、输入及缓存输出哈希，再使用已保存数值更新展示；不替代完整复现命令。所有假设采集周期都明确标为实验条件，不解释成已识别的仪器参数。运行前后核对既有受跟踪文件和100个冻结文件的哈希。

实际设备时间戳、PLC指令、采样管路和分析仪标定不在当前数据字段内。0.5秒记录间距不能当作0.5秒物理时滞精度。仅源注入接受也不是工业误报率。
