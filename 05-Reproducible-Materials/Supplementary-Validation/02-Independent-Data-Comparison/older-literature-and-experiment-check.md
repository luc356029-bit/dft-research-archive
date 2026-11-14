# 旧文献与实验核验：晶态 LiSi（I41/a）

## 2026-11-09｜阶段状态：完成检索与第一手证据核验

**结论：找到独立的同相、单空位 NEB 证据；没有找到可无附加假设精确转成当前口径的旧文献 D(300 K) 或实验基准。** 最有用的是 Yao 2014 的 shortest monovacancy path 0.243 eV，以及 Moon 2016 的四类近邻空位路径 0.2348–0.4658 eV。Moon 的 0.306 eV 为 KMC 有效能垒；其具体拟合区间、逐路径振动前因子、空位浓度和 MSD 归一未全部取得，不能把它与本项目 300 K 局部网络斜率 0.250704 eV 当作相同定义。

本报告严格区分：(1)已看全文；(2)已从出版社官方图核验；(3)作者论文可补足的细节；(4)待取得全文的候选。数值接近只用于描述差值，不作为端点路径映射或验证通过的证据。未修改原代码、未运行 DFT、未上传或联系作者。

### 当前比较口径（由主任务提供，非本报告重算）

- 晶态 LiSi，mp-795 / I41/a；64 宿主位点，单中性空位 Li31Si32；PBE/PAW，固定宿主晶胞。
- 距离 [2.69496, 2.92847, 2.94600, 3.05385] Å；对应 Em [0.2534438, 0.2358302, 0.2459095, 0.4276611] eV。
- 共同假设 ν=1e13 s^-1；D_vac(300 K)=1.697203e-7 cm²/s；c_v=1/32 的 Li tracer D(300 K)=2.01231e-9 cm²/s；网络 300 K 局部斜率 Ea=0.250704 eV。
- Kumar 2020 本轮不作为最终 D 准确标杆：缺陷/粒子归一不清、Fig.23 与正文的经典/量子 D0 标签存在对调问题（沿用主任务审计结论）。

## 1. Yao et al. 2014：独立的同相单空位 NEB，适合 Em 层面核对

原文：H. R. Yao et al., *Lithium Ion Migration in Li-Si alloys: From First Principles Studies*, Int. J. Electrochem. Sci. 9, 1854–1866 (2014)。[出版社全文 PDF](https://www.electrochemsci.org/papers/vol9/90401854.pdf)，[现出版社 DOI 页面](https://doi.org/10.1016/S1452-3981(23)07897-5)。本地 `older/Yao2014.pdf`（13 页），已用 py -3.14 / fitz 全文提取为 `Yao2014.txt`，并视检 p1857 的 Fig.1。

### 方法与准确可提取的数据

| 项目 | 原文证据与位置 | 与本项目关系 |
|---|---|---|
| 目标晶相 | p1856 Table 1：Li1Si1，I41/a；a=b=9.353 Å，c=5.743 Å，直角 | 相种直接匹配，非泛称 Li–Si alloy |
| 电子结构 | p1855 Computational details：VASP，PAW，**PW91-GGA**，450 eV；relaxed atom force <0.005 eV/Å；k-mesh denser than 0.03 Å^-1 | 与本项目 PBE 不同；不可把差异只归因于收敛误差 |
| 缺陷 | p1856：从完美晶格移除一个 Li，只考虑 monovacancy | 单空位机制相符；但 bulk 超胞原子数/空位分数、NELECT/电荷设置未写明 |
| 路径 | p1857：选择 **shortest migration distance** 的一条路径，编号 1-Li1Si1；Fig.1 含结构示意与 NEB profile | 作者只采样一条最短路径，不能代表四类网络 |
| Em | p1857 正文与 Fig.1：**0.243 eV** | 可作为独立最短路径 Em 候选；不能作为网络 Ea |
| 其他机制 | p1857–1858：bulk interstitial 的最低报道障碍仍 1.172 eV；p1859 的 0.198/0.224/0.347 eV 属于(100) surface/subsurface | 表面/间隙值不混入本项目 bulk 单空位比较 |

**路径映射程度：tentative shortest-topology。** 本项目最短边 2.69496 Å / 0.2534438 eV，是合理的首选对应候选；差值为 +0.0104438 eV。依据是作者明确写“最短跳跃”，不是能垒相近。Fig.1 未给 Li 原子编号、端点分数坐标或投影方向；图的 Reaction Path(Å) 不应直接等同于端点 Li–Li 距离。故尚不能给出已验证的一一对应，也不能证明四个 Em 全部准确。

### 空位形成能与 D 不可转换的原因

p1856 **明写 bcc 金属 Li 参照**：

`E_for = E[Li_(m-1)Si_n] + E[Li_bcc] - E[Li_mSi_n]`

其中 E[Li_bcc] 被作者定义为 bcc Li metal 中每一个 Li 原子的能量；**不是文中明示的孤立 Li 原子参照**。p1856–1857 给 bulk 空位形成能约 3.0 eV。该值与 Moon 学位论文的 1.2 eV 不一致，但现有证据不足以说明原因；不能擅自改成孤立原子参考，也不能断言哪篇数值错误。

p1864 §3.3 只作数量级估算：Γ=ν0 exp(-Ea/kBT)，ν0 取 **1 THz**，D=Γl²，Li1Si1 bulk Em=0.243 eV，D 约 **1e-7 cm²/s**。这段没有给出估算温度、跳距数值、空位浓度、三维几何/相关因子，且不是多路径 KMC。即使把 ν 从 1e12 线性放大到 1e13，也仍不能获得本项目同口径 D_vac 或 Li tracer D；本报告不采用其 D 作准确标杆。

## 2. Moon et al. 2016：官方图确认四近邻障碍；KMC D 仍缺严格归一

原文：J. Moon, B. Lee, M. Cho, K. Cho, *Ab initio and kinetic Monte Carlo study of lithium diffusion in LiSi, Li12Si7, Li13Si5 and Li15Si4*, J. Power Sources 328, 558–566 (2016)，[DOI](https://doi.org/10.1016/j.jpowsour.2016.07.092)，[ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0378775316309673)，[作者机构摘要](https://khu.elsevierpure.com/en/publications/ab-initio-and-kinetic-monte-carlo-study-of-lithium-diffusion-in-l/)。标题的 Li13Si5 与摘要/正文片段的 Li13Si4 不一致；本报告只用 LiSi 部分。

完整期刊 PDF 本轮未取得；ScienceDirect 403；Elsevier API 返回元数据且标明非 OA（已留存 `Moon2016_elsevier.xml`）。常见 SI 地址 `...-mmc1.pdf` 返回 404，不能据此断言不存在 SI。

**成功取得并视检的出版社官方图（这不是二手摘录）：**

- [Fig.2，LiSi 距离/结构/五类路径](https://ars.els-cdn.com/content/image/1-s2.0-S0378775316309673-gr2.jpg)，本地 `Moon2016_gr2.jpg`。
- [Fig.4，LiSi 空位 NEB profile](https://ars.els-cdn.com/content/image/1-s2.0-S0378775316309673-gr4.jpg)，本地 `Moon2016_gr4.jpg`。
- [Fig.8，300–700 K 的 KMC D](https://ars.els-cdn.com/content/image/1-s2.0-S0378775316309673-gr8_lrg.jpg)，本地 `Moon2016_gr8_lrg.jpg`；未把图读数伪装成作者原始数值表。

另取得第一作者 Janghyuk Moon 2015 博士论文 *Multiscale mechanical modeling of high capacity anode in lithium ion batteries*：[SNU 原文](https://s-space.snu.ac.kr/bitstream/10371/118464/1/000000053290.pdf)，本地 `Moon_thesis.pdf` 149 页，全文 `Moon_thesis.txt`。其 Fig.2.4/2.6/2.9 与期刊官方 Fig.2/4/8 一致，故可补充解读；**方法细节的 thesis 页码不能冒称是 2016 期刊页码**。

### 路径证据与候选映射

官方 Fig.2(a) 的颜色/距离与 Fig.4 的颜色/Em 对应如下；已从图片直接核验：

| 颜色/拓扑 | Moon d (Å) | Moon Em (eV) | 本项目候选 d / Em | 映射状态 |
|---|---:|---:|---|---|
| Green，Li 四面体簇之间 | 2.644 | 0.2725 | 2.69496 / 0.2534438 | 最短边唯一，属几何排序候选；未验证端点分数坐标 |
| Red，簇内 | 2.900 | 0.2519 | 2.92847 / 0.2358302 | 红/棕簇内两边很接近，不能仅用长度排序断言颜色对应 |
| Brown，簇内 | 2.911 | 0.2348 | 2.94600 / 0.2459095 | 与上一行视为无序候选对；需本项目端点/位移方向才能解歧 |
| Blue，簇之间 | 3.003 | 0.4658 | 3.05385 / 0.4276611 | 较长近邻桥接边候选；未验证完整几何 |
| Purple，簇之间 | 4.128 | 0.6182 | 无本项目≤3.1 Å对应 | 超出 thesis KMC 的 <3.5 Å 跳跃候选阈值 |

这张表支持“相同晶型下存在三条约0.23–0.27 eV的低障碍近邻与一条约0.47 eV较高近邻”的独立结构化证据。它不支持“各条误差只有几meV”或“一一复现”。若按距离排序暂配，四条本项目-minus-Moon 的差值依次为 -0.0190562、-0.0160698、+0.0111095、-0.0381389 eV；其中中间两条不应当作已映射误差。不能凭“LDA通常收缩”把这些差值定量修正为零。

### 方法、KMC 与 Ea 定义

- 期刊出版社正文片段确认 **LDA/PAW/VASP，Li_sv（1s纳入价电子），500 eV**，与本项目PBE不同。完整 k-grid/缺陷 NEB 超胞细节未取得。
- thesis §2.4.1.1 p18（PDF31）：I41/a，a=9.332 Å、c=5.704 Å；但 thesis Table2.1 p30（PDF43）又写 a=9.000、c=5.319。两处不一致，不能选一个偷偷当作确定的期刊使用晶格。
- thesis §2.5.1 p26（PDF39）：从空位位置选择 **<3.5 Å** 的 Li 邻居；ν通过初态/鞍点振动频率比求得；事件驱动 KMC，D=lim <r²>/(6t)。**未给逐路径ν数值，未明确 <r²> 是所有Li、移动Li还是空位的平均、也未给Li空位个数。**
- thesis p25–27（PDF38–40）：10个样本平均，周期边界、up to 10×10×10 supercell、50 ns。无法从“up to”还原单一缺陷浓度。尤其不能擅自当作cv=1/32，或当成逐空位归一D。
- 期刊 **abstract** 明确 LiSi KMC 有效能垒 **0.306 eV**。这是由KMC宏观迁移再得到的有效量，不是Fig.4任一单步NEB障碍。已取得的官方Fig.8给300/400/500/600/700 K各点，但未取到明确的拟合公式、温区、权重、D0/误差；本轮只能说“从KMC得到”，不能具体断言0.306是300–700K等权Arrhenius回归或300K局部导数。
- Fig.8 的 LiSi 300 K 点目视是几×1e-11 cm²/s（图像估读），不是可精确引用的原始D。thesis正文仅给LixSi约1e-10的总数量级。该数量级与本项目并非完全相同，但**因ν、cv、归一和可能拟合温区不清，不得用差值直接宣告验证成败**。

## 3. Wang G. et al. 2016 应变论文：零应变0.25 eV，已获官方图，方法待全文

G. Wang, J. Shi, M. Wu, C. Ouyang, B. Xu, *The strain effect on lithium ion migration in Li-Si alloys: A first-principles study*, Solid State Commun. 247, 47–52 (2016)，[DOI](https://doi.org/10.1016/j.ssc.2016.08.008)，[出版社页面](https://www.sciencedirect.com/science/article/pii/S0038109816302022)。

全文未取得，出版社页面本轮403，API仅元数据（`WangG2016_elsevier.xml`），FULL视图401。未用二手论文中的重绘数值填补。

出版社官方图已下载并核验：

- [Fig.1晶胞](https://ars.els-cdn.com/content/image/1-s2.0-S0038109816302022-gr1.jpg)，本地 `WangG2016_gr1.jpg`。
- [Fig.2零应变NEB](https://ars.els-cdn.com/content/image/1-s2.0-S0038109816302022-gr2_lrg.jpg)：图内明确 **0%-Li1Si1，0.25 eV**，本地 `WangG2016_gr2_lrg.jpg`。
- [Fig.3应变两端](https://ars.els-cdn.com/content/image/1-s2.0-S0038109816302022-gr3.jpg)：**-10%-Li1Si1：0.56 eV；+10%-Li1Si1：0.07 eV**。

零应变0.25 eV可作为额外较低精度Em候选，表明应变条件会明显改变该研究中的障碍。尚缺全文 methods、确切晶体空间群、NEB缺陷数/电荷/超胞、应变张量定义、端点坐标和prefactor；因此不能证明它与当前Li31Si32同口径，也不能用±10%曲线为当前PBE与MoonLDA差异作数值校正。与Yao2014作者组重叠，不能把两篇简单当作完全独立计算体系计权。

## 4. Wen & Huggins 1981：相种错误和扩散定义错误均排除

C. John Wen, R. A. Huggins, *Chemical diffusion in intermediate phases in the lithium-silicon system*, J. Solid State Chem. 37, 271–278 (1981)，[DOI](https://doi.org/10.1016/0022-4596(81)90487-4)，[全文PDF镜像（原始论文内容）](https://bdt.semi.ac.cn/download/0.3553448840748338.pdf)。本地 `Wen_Huggins1981.pdf` 8页已fitz全文读取，p277 Table I已渲染视检，避免OCR把10^-7读成10^-5。

- p271 abstract、p272–273、p277 Table I：研究四相是 **Li12Si7、Li7Si3、Li13Si4、Li22Si5**。**没有晶态LiSi / Li1Si1。** 文题“lithium-silicon system”不能解释为LiSi化合物。
- 温度 **415°C = 688.15 K**，不是700K；p272–273是三电极Al/LiAl | LiCl-KCl eutectic | LiySi、potentiostatic step的长时间解。
- p277 Table I的化学扩散D~为四相分别(6.60–8.13)、(4.16–5.46)、(7.00–9.33)、(5.13–7.24)×1e-5 cm²/s；相应TEF为122–160、87–111、120–325、147–232。
- p275 Eq.(3)用TEF把化学扩散折算为Li self-diffusion，假设Li比Si扩散快；p277表的Li self D分别(4.85–5.44)、(4.16–4.96)、(2.88–5.01)、(2.99–3.45)×1e-7 cm²/s。**这是热力学折算的self D，不是同位素示踪直接测得D，也不是单位空位D。**

结论：即使进行TEF校正仍然相种不符；缺温变Ea也不能换到300K。不能把“~6e-5化学D”拿来校验LiSi单空位，更不能把由其折算的其他相~1e-7 self D当作LiSi实验证据。Moon图中EXP700K线也是跨相背景，不能解释为本项目目标相的700K实验标杆。

## 5. Stearns et al. 2003：目标相正确，NMR没有跳跃率/扩散Ea

L. A. Stearns et al., *Lithium monosilicide (LiSi), a low-dimensional silicon-based material prepared by high pressure synthesis: NMR and vibrational spectroscopy and electrical properties characterization*, J. Solid State Chem. 173, 251–258 (2003)，[DOI/出版社](https://doi.org/10.1016/S0022-4596(03)00045-8)。[原始全文镜像转录](https://www.academia.edu/22522704/Lithium_monosilicide_LiSi_a_low_dimensional_silicon_based_material_prepared_by_high_pressure_synthesis_NMR_and_vibrational_spectroscopy_and_electrical_properties_characterization)。electronicsandbooks PDF curl多次TLS失败；Academia curl403，本轮通过web工具查看整篇原文转录，**没有已下载PDF**。

p255 Table1确认I41/a、a=9.354(3)、c=5.746(3) Å。p253的7Li MAS为9.4T、155.5MHz；p256–257 Fig.9用于化学位移及四极线形，作者据室温四极线形认为Li在晶位上静止。没有变温弛豫/跳跃率或Li示踪D。p257 Fig.10在250–650K测DC电导：低于450K的**0.057eV是电子窄带隙**；450–525K上升仅被讨论为Li+贡献的可能起点。因此0.057eV不能当Li迁移Ea，电导也不能直接转D。

## 6. Chevrier/Zwanziger 与 Cui 等的追溯结论

- Chevrier, Zwanziger, Dahn 2010, J. Alloys Compd.496,25–36，[出版社原文摘要/节选](https://www.sciencedirect.com/science/article/pii/S0925838810001969)，[DOI](https://doi.org/10.1016/j.jallcom.2010.01.142)：覆盖LiSi的电荷转移、电子结构、声子与谐近似热力学，没有可供本轮使用的Li跳跃率/NEB/KMC D或Ea；不能拿总声子DOS代替跃迁态/初态振动比。
- Cui, Gao, Cui, Qu 2012, JPS207,150–159，[DOI](https://doi.org/10.1016/j.jpowsour.2012.01.145)，[原始全文镜像转录](https://www.academia.edu/117042632/A_second_nearest_neighbor_embedded_atom_method_interatomic_potential_for_Li_Si_alloys)：pp153–157的2NN-MEAM/LAMMPS，NPT，1fs；平衡1000ps、采样2000ps。D由所有扩散原子的MSD求；300K由高温Arrhenius外推，且比较crystalline/amorphous。未发现规定单空位/cv=1/32的证据；本轮没有将其图8的晶态D当精确数值取读，也不把经验势+高温外推当当前固定缺陷NEB网络的可换算基准。
- 检索“LiSi+NMR/jump/tracer”和“lithium monosilicide+diffusion”时，常见的0.18/0.32/0.55eV及NMR/PFG数据来自**Li12Si7**；Kuhn JACS2011、Dupke等的这些记录不能换到LiSi。含“LiSi”的lithium disilicate玻璃记录实为含氧硅酸盐，同样排除。未找到目标I41/a LiSi的直接示踪或弛豫跳频定量实验；这是本次检索结果，不是证明不存在。

## 最终分级及后续可执行项

| 文献 | 当前能支持的判断 | 不可支持的判断 |
|---|---|---|
| Yao2014 | 同相、最短bulk monovacancy Em=0.243eV；最短拓扑候选 | 准确D300、全部四边网络、严格端点映射、实验验证 |
| Moon2016 Fig2/4 | 同相五类跳跃，其中四近邻Em与距离已从官方图核实 | 当前PBE边与颜色已一一匹配；D/ν/cv已可精确换算 |
| Moon2016 abstract | KMC有效障碍0.306eV | 300K局部network斜率或确定区间/权重的Arrhenius回归 |
| WangG2016 Fig2 | 0%-Li1Si1障碍0.25eV；应变图显示敏感性 | 未取得methods时声称固定单空位Li31Si32准确复现 |
| Wen1981 | 其他四种晶态硅化锂的化学D与TEF | LiSi目标相的实验D/Ea |
| Stearns2003 | LiSi结构/静态MAS/电导性质 | 0.057eV为Li迁移Ea；NMR有跳频或tracerD |

关键决定：仅将旧文献保留为“Em机制/数量级的独立支持”，不设一个可转换后的准确D300金标准，也不因为D恰好接近就通过验证。影响：本项目若要声明外部定量一致，仍须补齐端点坐标映射及频率/缺陷归一，或另找明确同口径研究。

未决项/风险：Moon完整正文/SI仍缺，无法核其NEB原胞大小、prefactors和Ea拟合定义；Yao bulk尺寸/NELECT未给；WangG2016 methods未取得；实验无目标相定量跳频。建议下一步由主任务结合自身端点结构验证Moon四类几何及multiplicity，优先补齐Moon正文原始表或作者公开输入；不以缺失信息自动构造换算。
