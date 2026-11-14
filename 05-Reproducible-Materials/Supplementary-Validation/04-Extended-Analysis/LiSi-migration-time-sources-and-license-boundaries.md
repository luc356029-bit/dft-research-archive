# 2026-11-13实际验证：来源、配置与复现边界

本轮依照用户新指令实际执行两套公开势推断，不再使用前轮“没有新计算授权”的状态。继续使用Google Drive文件管理技能，在原codex/DFT复刻目录交付；不更改共享权限。

公开输入：download_receipt.json给出两种权重、ænet固定commit源码、release和许可的来源、字节与SHA；ænet权重还与上轮作者仓库Git blob身份交叉核对。MACE官方release未给独立摘要，来源/大小与本地SHA不冒称发布方独立hash。ASL链接最初路径错误404，已按作者仓库实际ASL.md解决，正确回执在license_resolution.json。MACE权重为学术非商业ASL，不能统称无限制开源权重。

软件：原版ænet Fortran/Cython/ASE接口，不自写NN内核；官方Lbfgsb PIC对象补编译后重链接，未改源文件数学公式。编译二进制及SHA保留。MACE官方0.3.14/CPU/Torch2.9.1、float64、无附加色散，完整层及元素配置在metadata。初始Windows挂载环境I/O拖慢，已仅停止本任务安装，换WSL ext4隔离依赖环境；输入/输出仍写项目独立目录，线程≤6，未用超算。

数据：reference_dataset.json以60个已有源文件/生产脚本SHA绑定246个构型及三能量/物理力，source_snapshots保留原NEB审查报告。没有重新做结构映射或DFT。逐原子力分量来自DFT物理力，非弹簧/投影力；相同组成每路径端点均值对齐，不比较不同模型绝对总能量。代码独立存在，不修改原tracer_kmc.py、transport.py、network_arrays.py或生产输运模型。

校验：本轮额外对Li/Si坐标的能量导数及两种PBC变换做直接模型调用，输出interface_checks.json；这是接口数值核验，不是MD。两模型各246个输出保留全部能量/力；ænet第一次已完成全部构型但metadata NumPy bool序列化失败，只恢复metadata并重做少量接口检查，原246输出不重算、不改动，旧脚本和失败日志保存。计时口径/RSS边界在metadata及报告中说明。

训练重叠：下载并完整流式扫描ænet公开XSF归档及MACE公开回放文件，provider checksum/本地SHA和计数保留。MACE全部MATPES微调/OMat预训练没有穷尽，不能宣称全链无泄漏；公共训练标签不当独立热态验收。

科学停止：两候选均未通过四路径静态入口；独立800K目标单空位热态DFT亦缺。没有ML MD、实际迁移时间分布、κ、返回率或记忆数据，不填0，不用1/ν。KMC时钟代数和DFT成本是已有结果/显式条件场景，不是新AIMD计时或迁移时间实测。没有新DFT、训练、HPC任务或速率更改。

复现：在新副本保存reference_dataset和脚本，执行reproduce.sh（会拒绝覆盖旧预测）；依赖按requirements-lock.txt安装。权重/固定源码和许可另包，两个大公共训练数据文件不重复上传，来源/精确SHA/扫描脚本完整提供。全部新数值原始输出、误差表及曲线随结果包，不把文件完整性通过当成动力学物理通过。
