# 补充检验与后续分析

本部分接续 A → B → C1 主线，记录内部核查、独立数据对照、差异归因及动力学模型分析。

| 阅读顺序 | 主要问题 | 入口 |
|---|---|---|
| 1 · 内部一致性 | 路径截断是否合理，正反端点是否等价，输出异常是否影响结果？ | [截断依据](01-Internal-Consistency/LiSi-3.3A-cutoff-selection.md) · [端点核查](01-Internal-Consistency/forward-reverse-NEB-endpoint-equivalence.md) · [NaN 排查](01-Internal-Consistency/NEB-POSCAR-NaN-diagnosis.md) |
| 2 · 独立数据对照 | 论文与本项目是否比较了相同路径、条件和物理量？ | [Kumar 首轮比较及更正](02-Independent-Data-Comparison/Kumar2020-first-round-comparison.md) · [统计口径](02-Independent-Data-Comparison/Kumar-statistical-definition-and-figure-differences.md) · [同口径空位 MSD](02-Independent-Data-Comparison/same-definition-vacancy-MSD-comparison.md) · [公开数据综合对照](02-Independent-Data-Comparison/LiSi-public-data-comparison.md) |
| 3 · 差异归因 | Em、D、Ea 的差异能解释到哪一步？ | [最终归因报告](03-Error-Attribution/LiSi-final-Em-D-Ea-attribution.md) · [数值总表](03-Error-Attribution/Em_D_Ea_attribution.csv) · [同模型 Ea 对照](02-Independent-Data-Comparison/LiSi-Ea-independent-data-same-model.md) |
| 4 · 动力学扩展 | KMC 时钟与真实迁移过程的关系如何，现成 ML 势能否用于检验？ | [时钟与迁移时间](04-Extended-Analysis/LiSi-KMC-clock-and-transition-time.md) · [ML 势检验](04-Extended-Analysis/LiSi-migration-time-final-validation.md) · [TST 与有效速率](04-Extended-Analysis/LiSi-TST-and-effective-hop-rate.md) |

[核查代码](03-Code/)与[来源、状态及校验记录](05-Provenance-and-Manifests/)用于追溯具体分析。[2026-11-12 阶段记录](00-Project-Record/phase-record-20261112.md)概括关键归因工作；[早期外部验证方案](01-Internal-Consistency/final-external-validation-plan.md)保留历史规划及未执行状态。
