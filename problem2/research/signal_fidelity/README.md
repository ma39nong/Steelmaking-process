# Signal Fidelity Benchmark

固定研究目录，正式基线d08f989，研究基线0ee55e2。先验证积分形式辨识，再在同一基准内检查测量假设与各处理环节。正式V2、第一问、第三问接口、冻结论文和此前研究输出只读。本目录不生成最优方案选择锁，不更新正式参数。

```bash
.venv/bin/python -m problem2.research.signal_fidelity.run --integral-only
.venv/bin/python -m problem2.research.signal_fidelity.run
.venv/bin/python -m pytest -q problem2/research/signal_fidelity/test_fidelity.py
```

完整配置在config.json。流量脉冲使用解析积分形成累计量，再对整个基线和整个叠加信号分别应用同一采集算子。基准同时测试短时/缓变信号、不同记录保持假设、逐炉背景、恒定真流量和孤立记录尖峰。基线重测差异单独保存，不能将二次假设采集当作已知真实仪器机制。

源码只用已有工程依赖。向量化Hampel是原规则的计算加速，经过逐点掩码等价测试及十炉重放；没有更改生产实现。数据、代码、参数及输出哈希记录在results中。重复执行原位更新。

积分探针是参考积分/连续时间辨识思想后针对累计耗氧量改写的一阶加纯迟延模型，不是原论文在转炉上的已验证实现。只做受控辨识及失配/无关联/恒定输入检验，不开展真实工业时滞拟合。具体来源及适用条件见报告。
