# A · 高通量筛选

从 Materials Project 获取二元锂合金，通过元素排除和能量高于凸包阈值筛选候选材料，为后续结构与扩散路径分析提供输入。

| 顺序 | 内容 | 文件 |
|---|---|---|
| 1 | 筛选目标、条件、学习方法与论文差异 | [筛选方法与方法论示范](02-Methods-and-Parameters/A-screening-methodology.md) |
| 2 | 代码理解与执行逻辑 | [筛选笔记](02-Methods-and-Parameters/A2.1-screening-notes.md) · [筛选脚本](03-Code/A_screen.py) |
| 3 | 候选列表与结构 | [alloy_list.csv](05-Outputs/alloy_list.csv) · [77 个 CIF](01-Inputs/CIF/) |

方法论示范按「方案验收 → 核心知识预习 → 分步实践与反馈」推进；阶段性学习记录由「3.代码执行、数据流」串联，报错、合理性检验和疑问的记录方式见该页对应章节。

历史记录中的筛选数量与论文数量存在差异，具体口径和排查过程见筛选方法。CIF 保留 Materials Project 编号，便于追溯结构来源。

下一阶段：[B · 结构与扩散路径](../B-Structure-and-Diffusion-Pathways/README.md)。
