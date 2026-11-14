# C1.2.5 · KMC 数据流与结果读取

以 LiSi 首轮 300 K 结果为例，计算链为：**NEB 能垒 → 候选事件与跳率 → 位置/时间更新 → 分块统计 → Li 示踪扩散系数**。模型含 32 个 Li 位点、31 个 Li 与 1 个空位。

## 输入与事件

[transport_model.json](../05-Outputs/KMC/transport_model.json)保存网络与能垒；[KMC_ext_00.json](../05-Outputs/KMC/KMC_ext_00.json)保存单次模拟条件、候选事件及分块结果。

| 字段 | 含义 |
|---|---|
| `T_K`、`attempt_frequency_s` | 300 K，共同尝试频率假设为 1e13 s⁻¹ |
| `nsites`、`Li_count`、`vacancy_fraction` | 32 位点、31 Li、空位比例 1/32 |
| `barriers_eV` | Path00—03 的对称化 Efr 能垒 |
| `supplied_neighbor`、`supplied_rates_s` | 每个空位位置的候选邻居与逐事件速率 |
| `supplied_displacement_A` | 对应边的未折回空位位移，单位 Å；Li 位移方向相反 |

0 号空位的邻居为 `[18,22,26,14,7,8]`，程序按速率权重选择六个具体事件；四类路径与六个出通道不是同一计数。若以两个 0.5 随机数作教学演示，则空位 0→26，原位点 26 的 Li→0；这不是历史第一跳的回放。

更新逻辑见 [tracer_kmc.py](../03-Code/KMC/tracer_kmc.py)，速率与输运处理见 [transport.py](../03-Code/KMC/transport.py)。

## 分块汇总

主结果使用 `KMC_ext_00.json`（seed 137）和 [KMC_ext_02.json](../05-Outputs/KMC/KMC_ext_02.json)：各统计 800 万次交换、160 块，每块 5 万次，另有 2 万次预运行不计入统计。

每块按 `D_b = MSD_b/(6*t_b)` 求值，再按 `sum(D_b*t_b)/sum(t_b)` 汇总；Å²/s 转为 cm²/s 乘 1e−16。得到 `D_Li = 2.012311292523e−09 cm²/s`，KMC 1SE 约 `1.967569586e−11 cm²/s`，不再额外乘空位比例。

`KMC_ext_01.json` 与 `KMC_ext_03.json` 是同种子轨迹的 320 块对照，不能作为新增独立轨迹合并。现有 JSON 保留分块 D 与时间，不含每次跳跃的完整坐标。

## 结果状态

[统计记录](../06-Validation/KMC_statistical_admission.json)、[相关误差补充](../06-Validation/correlation_error_supplement.json)和[有限模型最终记录](../06-Validation/finite_model_final_admission.json)给出检验范围。旧文件的 `PENDING` 是生成时状态，应结合后续记录与源 SHA 阅读。

结果对应固定有限空位浓度、四类路径与共同前因子假设下的 Li 示踪扩散；统计误差不覆盖 DFT 能垒与动力学模型误差。后续解释见[数据处理与对照](../02-Methods-and-Parameters/C1-2.6-analysis-and-comparison.md)。
