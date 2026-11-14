# 02 · 独立 KMC 输运计算与验证

得到四类路径的 NEB 能垒后，接下来要回答的是：这些局部跳跃，最终会让 Li 扩散得多快？

文章先计算各类路径的扩散系数，再按路径数量加权。我选择另外建立一条计算路线：**把能垒放进扩散网络，让 Li 按照相应速率迁移，再从实际产生的位移得到 D(T) 和 Ea**。

## 1. 从能垒到随机迁移

最初的 KMC 模型有 32 个 Li 位点，其中放置 31 个 Li 和 1 个空位。每次由一个相邻 Li 与空位交换，允许发生的跳跃来自已经完成 NEB 的四类路径。

第 i 类路径的跳跃速率为：

$$
\Gamma_i(T)=\nu\exp\left(-\frac{E_{m,i}}{k_BT}\right)
$$

其中，Em 决定跳跃的能量门槛，尝试频率 ν 决定时间尺度。本次计算统一取 ν = 10¹³ s⁻¹。

知道速率后，还要解决两个问题：下一次跳哪条路径，以及这次跳跃前等待多久。

$$
P_j=\frac{\Gamma_j}{\sum_j\Gamma_j},
\qquad
\Delta t=-\frac{\ln u}{\sum_j\Gamma_j},\quad 0<u<1
$$

按概率选中事件后，更新迁移 Li 的身份与位置、空位位置、累计位移和模拟时间，再继续下一次。

Li 跨出计算盒时，会从另一侧回到盒内，但累计位移保留跨越周期边界的部分。这样才能记录它逐渐走向远处的运动。

## 2. 从位移得到 D，再从温度变化得到 Ea

这里需要区分路程和位移。一个 Li 来回跳了很多次，走过的路程很长，但最后可能仍在起点附近。要描述长程扩散，就需要统计每个 Li 相对起点的净位移：

$$
\mathrm{MSD}_{\mathrm{Li}}(t)
=\frac{1}{N_{\mathrm{Li}}}
\sum_{\alpha=1}^{N_{\mathrm{Li}}}
|\mathbf r_\alpha(t)-\mathbf r_\alpha(0)|^2,
\qquad
D_{\mathrm{Li}}=\lim_{t\to\infty}\frac{\mathrm{MSD}_{\mathrm{Li}}(t)}{6t}
$$

在 250、300、350、400 K 分别计算 D(T)，再对 ln D 与 1/T 拟合：

$$
\ln D(T)=\ln D_{\mathrm{fit},0}-\frac{E_a}{k_BT}
$$

得到的主要结果为：

| 项目 | KMC 结果 |
|---|---:|
| 300 K 的 Li 示踪扩散系数 | 2.0123 × 10⁻⁹ cm²/s |
| 250–400 K 拟合的 Ea | 0.25161 eV |

Em 描述某一条路径的跳跃门槛，Ea 描述整体 D 对温度变化的响应。多条路径共同参与、回跳和连通方式随之产生影响后，Ea 要从整体输运中求出。

## 3. 计算结果是否可靠

KMC 有随机性，得到一次结果后，还需要检查采样是否足够。

初始每个随机种子运行 100 万次事件时，部分低频通道还没有被采样。后来扩展到每个种子 800 万次，并比较不同种子、不同统计块长度下的 D，同时检查是否仍有未采样的正速率边。

对不同统计块汇总时，也按各块实际经历的模拟时间加权，避免把事件数相同、时间却不同的块直接平均。

随后又用同一跳跃模型的主方程直接求长期输运，交叉检查 KMC 的随机采样结果。四个温度点的 Li 扩散系数都落在原 KMC 统计误差内，精确求值拟合得到 Ea = 0.25133 eV，与 KMC 的 0.25161 eV 相近。

两种求法得到一致结果后，这条从能垒到整体扩散的计算链就有了可以相互核查的依据。

## 4. 与文章的加权方法相比，多看到了什么

我把同一套四路径数据代入文章公式，取 Z = 6、g = 1，得到的 D300 约为 4.8621 × 10⁻⁷ cm²/s，比 KMC 大很多。

于是回到同一套 KMC 事件，分别统计每次跳跃距离的平方之和，以及每个 Li 最终净位移的平方：

| 300 K 的计算方式 | D，cm²/s |
|---|---:|
| 直接代入文章加权公式，Z = 6、g = 1 | 4.8621 × 10⁻⁷ |
| 按实际跳跃次数和跳距统计，平均到 31 个 Li | 1.5684 × 10⁻⁸ |
| 按 Li 净位移统计，即本项目 KMC | 2.0123 × 10⁻⁹ |

第一步相差约 31 倍，对应单空位模型中总跳跃量向 31 个 Li 的平均；第二步相差约 7.8 倍，体现了回跳等运动相关性对净位移的削弱。

这让我把原先笼统的“方法算出来不一样”，拆成了两个具体问题：统计量平均给了谁，以及是否考虑了前后跳跃之间的相关性。

**KMC 的价值就在这里：它把路径组合、空位交换和回跳放进了同一段运动过程，让整体 D 和 Ea 从这段过程里产生**。后续的低频路径分析与文献对照，也都以这套模型为基础。

## 原始记录与计算材料

- 学习记录：[C1 2.5 · KMC 模拟](../04-Learning-Notes/02-Stepwise-Learning/C1/C1%202.5KMC模拟.md)、[C1 2.6 · 处理与对比](../04-Learning-Notes/02-Stepwise-Learning/C1/C1%202.6处理、对比.md)。
- 计算入口：[KMC 代码与数据流](../05-Reproducible-Materials/C1-NEB-Calculations/04-Execution/KMC-practice-start-here.md)、[KMC 代码](../05-Reproducible-Materials/C1-NEB-Calculations/03-Code/KMC)、[模型与结果](../05-Reproducible-Materials/C1-NEB-Calculations/06-Validation/finite_model_final_admission.json)。
- 交叉验证：[统计检查](../05-Reproducible-Materials/C1-NEB-Calculations/06-Validation/KMC_statistical_admission_final.json)、[主方程求值与 KMC 对照](../05-Reproducible-Materials/Supplementary-Validation/02-Independent-Data-Comparison/LiSi-Ea-independent-data-same-model.md)、[精确求值代码与数据包](../05-Reproducible-Materials/Supplementary-Validation/04-Extended-Analysis/LiSi-Ea-same-model-evidence.zip)。
