# C1 · NEB 与输运计算

以 LiSi 为重点对象，从结构弛豫、空位构建和路径端点出发，经 NEB 获得迁移能垒，再用 KMC 与后处理计算扩散系数和活化能。先读[阶段总览](02-Methods-and-Parameters/C1-overview.md)。

| 顺序 | 方法与参数 | 执行与证据入口 |
|---|---|---|
| C1.2.1 · 方法比较 | [计算方案及适用范围](02-Methods-and-Parameters/C1-2.1-method-comparison.md) | [阶段总览](02-Methods-and-Parameters/C1-overview.md) |
| C1.2.2 · 结构弛豫 | [弛豫原理与参数](02-Methods-and-Parameters/C1-2.2-relaxation.md) | [16 → 64 → 63 原子结构与输入](04-Execution/early-practice/C1-2.2-practice-README.md) |
| C1.2.3 · 路径与端点 | [候选路径](02-Methods-and-Parameters/C1-2.3-diffusion-paths.md) | [几何与连通性记录](04-Execution/early-practice/C1-2.3-practice-README.md) · [端点映射](04-Execution/63-atom-NEB-endpoint-mapping.md) |
| C1.2.4 · NEB | [NEB 方法](02-Methods-and-Parameters/C1-2.4-NEB.md) | [Path00 参数与能垒读取](04-Execution/early-practice/C1-2.4-practice-README.md) · [NEB 审查数据](06-Validation/neb_analysis.json) |
| C1.2.5 · KMC | [KMC 方法](02-Methods-and-Parameters/C1-2.5-KMC.md) | [数据流与结果说明](04-Execution/KMC-practice-start-here.md) · [代码](03-Code/KMC/) · [输出](05-Outputs/KMC/) |
| C1.2.6 · 处理与对照 | [数据分析](02-Methods-and-Parameters/C1-2.6-analysis-and-comparison.md) | [阶段检验](06-Validation/) · [补充检验](../Supplementary-Validation/README.md) |

`early-practice` 保存用于理解计算链的早期材料，正式模型状态见对应的运行与检验记录。输运结果需同时读取温度、空位浓度、路径集合及尝试频率假设，区分空位扩散与 Li 示踪扩散。
