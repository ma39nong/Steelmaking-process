# 修复引入峰与未消除污染峰

预声明的“no new unsupported peak”指修复操作新引入的无参考支持峰，不包括检测器漏检后仍存在的注入尖峰。计算参考谱±2通道包络；显著局部峰超过包络 max(500 counts,参考全帧峰10%) 才计无支持峰。分别输出总残留unsupported_peak_count，以及相对污染输入新出现的introduced_unsupported_peak_count；选型无新增峰门槛使用后者。该指标定义澄清不改变检测阈值、注入、补全或选型门槛。在打开确认前完成，首次开发结果保存在metric_definition_audit，不覆盖原始记录。首次开发尚未据此执行模型选择。
