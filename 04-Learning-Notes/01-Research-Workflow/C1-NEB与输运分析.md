# C1 · NEB 与输运分析

[04 · 学习记录](../README.md) · [05 · 本阶段计算材料](../../05-Reproducible-Materials/C1-NEB-Calculations/README.md)

## 0.位置
```
A.从数据库得到候选锂合金
B.分析结构、等价位点和扩散路径
C.NEB / AIMD ← 当前
D.数据分析
```
## 1.核心目的
（Li合金经过NEB计算和数据处理）得到材料总体Ea、室温D
## 2.方案
*本方案与文章的方案有一些差异，Ea最终由KMC求得*

| 分块           | 对象 → 操作 → 结果                         | 时间       |
| ------------ | ------------------------------------ | -------- |
| C1 2.1 方法核对  | 论文与现有方案 → 对照方法 → 明确区别、适用范围           | 1—2h     |
| C1 2.2 结构与弛豫 | LiSi晶体 → 完整结构/缺陷结构优化 → 合格初末态         | 2.5—3.5h |
| C1 2.3 扩散路径  | Li位点 → 候选路径分类、连通性检查 → 4类路径与网络        | 2—3h     |
| C1 2.4 NEB计算 | 初末态、路径 → NEB、收敛检查 → 各路径迁移能垒          | 3.5—5.5h |
| C1 2.5 KMC模拟 | 能垒、网络 → 连续跳跃模拟 → 各温度 D(T)            | 3.5—5h   |
| C1 2.6 数据与对照 | D(T)、NEB结果 → Ea拟合、作者方法有限对照 → 最终结果、结论 | 2.5—4h   |
| 主体合计         |                                      | 15—23h   |
**大致流程：**
前后两次弛豫 → 准备几种单步跳法 → NEB算各步Ea → KMC模拟多次跳跃得到 D(T) → 从温度变化规律求总体 Ea

## 3.代码执行、数据流
| 分步学习 | 对应计算材料（05） |
|---|---|
| [C1 2.1 · 方法对比](../02-Stepwise-Learning/C1/C1%202.1%E6%96%B9%E6%B3%95%E5%AF%B9%E6%AF%94.md) | [计算链总览](../../05-Reproducible-Materials/C1-NEB-Calculations/README.md) |
| [C1 2.2 · 结构与弛豫](../02-Stepwise-Learning/C1/C1%202.2%20%E5%BC%9B%E8%B1%AB.md) | [结构、参数与弛豫数据](../../05-Reproducible-Materials/C1-NEB-Calculations/04-Execution/early-practice/C1-2.2-practice-README.md) |
| [C1 2.3 · 扩散路径](../02-Stepwise-Learning/C1/C1%202.3%20%E6%89%A9%E6%95%A3%E8%B7%AF%E5%BE%84.md) | [路径与连通性数据](../../05-Reproducible-Materials/C1-NEB-Calculations/04-Execution/early-practice/C1-2.3-practice-README.md) |
| [C1 2.4 · NEB 计算](../02-Stepwise-Learning/C1/C1%202.4NEB%E8%AE%A1%E7%AE%97.md) | [NEB 参数与能垒读取](../../05-Reproducible-Materials/C1-NEB-Calculations/04-Execution/early-practice/C1-2.4-practice-README.md) · [NEB 检验](../../05-Reproducible-Materials/C1-NEB-Calculations/06-Validation/neb_analysis.json) |
| [C1 2.5 · KMC 模拟](../02-Stepwise-Learning/C1/C1%202.5KMC%E6%A8%A1%E6%8B%9F.md) | [KMC 代码与数据流](../../05-Reproducible-Materials/C1-NEB-Calculations/04-Execution/KMC-practice-start-here.md) |
| [C1 2.6 · 处理与对比](../02-Stepwise-Learning/C1/C1%202.6%E5%A4%84%E7%90%86%E3%80%81%E5%AF%B9%E6%AF%94.md) | [输运模型](../../05-Reproducible-Materials/C1-NEB-Calculations/06-Validation/transport_model_final.json) · [差异归因](../../05-Reproducible-Materials/Supplementary-Validation/03-Error-Attribution/LiSi-final-Em-D-Ea-attribution.md) · [同模型对照](../../05-Reproducible-Materials/Supplementary-Validation/02-Independent-Data-Comparison/LiSi-Ea-independent-data-same-model.md) |

## 4.报错、合理性检验
1. 超胞设置为4倍体积的非对称超胞，最短周期平移达到 9.3374 Å，达到了codex预设的9Å最小距离，但未进行更大超胞NEB能垒计算，不能保证这种情况下计算得到的迁移能垒一定准确（C1 2.2）
2. $\nu$ 取$10^{-13}$，可能导致D与实际值产生偏差，但不会对 Ea 造成影响
## 5.疑问
| 疑问                                                                                                                           | 根/叶 | 归属（根问题去哪深挖）                                                |
| ---------------------------------------------------------------------------------------------------------------------------- | --- | ---------------------------------------------------------- |
| 作者是否使用Li空位？<br>为什么要建立缺陷结构？<br>为什么完整结构和缺陷结构都要弛豫？只弛豫缺陷结构不行吗？<br>为什么要挖掉一个Li？真实材料中也存在这种空位吗？<br>（2-C1 2.2）                        | 根   | 《材科基》<br>5.1.1 金属中的点缺陷；<br>9.4.1 扩散机制；<br>2.3.1 原子间作用力和结合能 |
| 倒易晶格、布里渊区、波矢的数学定义、波矢和电子状态的关系<br>（2-C1 2.4）                                                                                   | 根   | 理论计算—DFT与固体电子结构                                            |
| 跳跃相关性是说Li会从A跳到B，也会重新跳回A，这不是正常情况吗？Li在B点跳跃，A就是其中一个候选点。这会有什么误差？而且MSD是在时间足够长的情况下计算的，这时候对应的Li已经静止，其位移就应该是离它初始跳跃点的距离<br>（2-C1 2.6） | 根   | 《材科基》<br>9.4.1扩散机制、9.4.2扩散通道                               |
| T越高，原子振幅越大，但频率不一定变大。但尝试频率应该和振幅有关吧，频率不变，振幅越大尝试频率不应该越大吗？<br>（2-C1 2.6）                                                         | 根   | 《材科基》<br>5.1.4原子振动、9.5.1温度                                 |



## 6.与论文相比的差异、疑点

1. 2-C1 2.2
文章未明确说明采用的什么扩散模拟机制，虽然提到了LiF的空位形成能，但不能说明其他合金也是采用这种形式
当前方案采用的是空位交换迁移扩散的模拟方式，可能会影响到和文章的可比性