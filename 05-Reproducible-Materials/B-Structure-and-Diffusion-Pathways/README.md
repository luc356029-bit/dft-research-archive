# B · 结构与扩散路径

以 [A 阶段 CIF](../A-High-Throughput-Screening/01-Inputs/CIF/) 为输入，将晶体结构转为带周期边界的 Li 候选迁移网络，再归并对称等价路径。

## 分步分析

先读[整体方法与问题记录](02-Methods-and-Parameters/B-analysis.md)，再按下表查看分步笔记与代码。

| 步骤 | 方法记录 | 代码 |
|---|---|---|
| B2.1 · 空间群 | [结构与空间群](02-Methods-and-Parameters/B2.1-space-group.md) | [脚本](03-Code/B2.1_structure_space_group.py) |
| B2.2 · 不等价 Li 位点 | [位点归类](02-Methods-and-Parameters/B2.2-inequivalent-Li-sites.md) | [脚本](03-Code/B2.2_inequivalent_Li_sites.py) |
| B2.3 · Li–Li 距离壳 | [距离壳与容差](02-Methods-and-Parameters/B2.3-distance-shells.md) | [脚本](03-Code/B2.3_Li_distance_shells.py) |
| B2.4 · 截断与扩散维度 | [周期连通性](02-Methods-and-Parameters/B2.4-cutoff-and-dimension.md) | [脚本](03-Code/B2.4_cutoff_and_dimension.py) |
| B2.5 · 物理跳跃 | [周期像与反向去重](02-Methods-and-Parameters/B2.5-physical-hops.md) | [脚本](03-Code/B2.5_physical_hops.py) |
| B2.6 · Npath 与 Ni | [对称等价路径归并](02-Methods-and-Parameters/B2.6-Npath-and-Ni.md) | [脚本](03-Code/B2.6_inequivalent_paths_Npath_Ni.py) |

## Npath 六轮对照

[六轮数据流报告](05-Outputs/Npath-six-cycle-dataflow-report.md)记录了候选规则的调整、逐轮结果与保留问题；[论文对照表](05-Outputs/Npath-paper-comparison.csv)提供数值入口。[归档索引](05-Outputs/Npath-reproduction-README.md)说明完整历史材料的结构。

上表是分步学习与分析脚本，六轮研究另有最终实现。候选计数接近论文并不等于恢复了全部路径身份；几何连通性也需要结合后续能垒判断实际输运。

下一阶段：[C1 · NEB 与输运](../C1-NEB-Calculations/README.md)。
