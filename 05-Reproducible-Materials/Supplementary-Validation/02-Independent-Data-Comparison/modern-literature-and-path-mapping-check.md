# 现代晶态 LiSi 公开扩散证据核查

## 2026-11-09 阶段记录

- 状态：2025 主文、2022 SI、2025/2026 LiSi 原始 NEB 与结构数据已取得；已完成标准晶胞与路径轨道映射。2025 SI PDF 尚未取得，MSD 的严格公式、拟合区间和不确定度未完全复核。
- 范围：只讨论晶态 Li1Si1（空间群 88，I41/a）及明确缺陷模型；没有用非晶、液态、稀 Li/c-Si 或其他化学计量晶相的 D 来验证本项目 Li31Si32。
- 决策：优先将 2025 的同相、单中性 Li 空位 NEB 作为迁移势垒外部对照；800 K AIMD D 作为有口径条件的量级对照。不能声称已验证绝对 D 或热平衡本征 D。
- 分析基于公开资料与既有项目数据，未新增 DFT/MD。

## 核心结论

1. 最接近本项目的公开数据是 Kirsch–Dreßler–Sebastiani 2025 的 **Li127Si128、一个中性 Li 空位、800 K、100 ps AIMD**。主文 Table 3 的 Li 原子扩散数值为 **1.0 × 10^-7 cm²/s**，不是空位位置 MSD 的扩散系数。它的 Li-site 空位浓度为 1/128，本项目为 1/32；不能将文献数值直接同本项目 300 K 单空位 D=1.6972×10^-7 cm²/s 相比。
2. 原始 NEB 能量提供三条单原子空位交换路径。严格几何映射的结果是 **文献 (6)_def → 本项目 path00；(7)_def → path02；(8)_def → path01**。本项目第四条 path03 在该文献中没有找到孤立 Li 空位对应路径。
3. 对应迁移势垒的本项目减文献差为 **+0.033158、−0.006701、−0.049041 eV**。这支持三条几何跳跃类别和 0.2–0.3 eV 迁移势垒量级，但数十 meV 差异仍能显著改变 300 K 速率。
4. 2022 年无缺陷晶态 LiSi 在 500/800 K 的有限轨迹中没有观察到 Li 扩散；这是有限观测结论，不是“真实 D 严格为零”。2026 年的非零 LiSi D 来自 Na 加入体系，不能作为孤立 Li 空位的直接验证。
5. **800 K D 的严格归一尚缺一个关键证据**：2025 主文把定义放在 SI section S1.2 / Eq. S5。还不能确认 N_Li 分母、时间原点平均方式、中心漂移扣除、拟合 lag 区间和统计误差。2022 SI 图表强烈交叉支持通常的总 MSD/(6t)，但不能代替 2025 的原式。

## 一手来源与已保存证据

| 来源 | 原始公开入口 | 本地证据与状态 |
|---|---|---|
| 2022, Atomistic Diffusion Pathways of Lithium Ions in Crystalline Lithium Silicides from ab Initio Molecular Dynamics Simulations, JPCC 126, 12136–12149 | https://doi.org/10.1021/acs.jpcc.2c01555 | 主文全文本轮 ACS HTTP403；不把未取到的参数当成已核。SI 已取到。 |
| 2022 SI, DOI .s001 | https://acs.figshare.com/articles/journal_contribution/20317266 ; 原文件 https://ndownloader.figshare.com/files/36301458 | `modern/ACS2022_SI.pdf`，21,607,263 bytes。Figshare API 的 MD5=370b4310a1c22e2d57b0cb9304820ff2，与下载吻合。 |
| 2025, Li+ diffusion in crystalline lithium silicides: influence of intrinsic point defects, J. Phys. Energy 7 025003 | https://doi.org/10.1088/2515-7655/ada639 ; https://iopscience.iop.org/article/10.1088/2515-7655/ada639/pdf | `modern/IOP2025.pdf`，真实 PDF，15 页含期刊封面；已抽文、视觉复核 PDF 第4、6页。 |
| 2025 原始数据（作者公开） | https://doi.org/10.5281/zenodo.14065754 ; API https://zenodo.org/api/records/14065754 | `modern/2025_paths_li1si1.zip` 与 `modern/2025_systems.zip`，保留原 zip；解压本在相应同名目录。数据只有松弛结构和 NEB 路径/能量，没有 AIMD XYZ、MSD 数值或 CP2K 输入。 |
| 2026, influence of extrinsic point defects and structural complexity, J. Phys. Energy 8 015004（在线发表 2025-11-06） | https://doi.org/10.1088/2515-7655/ae175e ; https://opendata.uni-halle.de/bitstream/1981185920/123640/1/Kirsch_2026_J._Phys._Energy_8_015004.pdf | 一手仓储的搜索缓存中核对 Table 3、摘要、数据可用性；实际下载遭 Anubis/Radware，尚无可用主文 PDF。不能冒充完整本地复核。 |
| 2026 原始数据（作者公开） | https://doi.org/10.5281/zenodo.15670547 （concept）；解析到版本 https://doi.org/10.5281/zenodo.15670548 | `modern/2026_paths_li1si1.zip`、`modern/2026_systems.zip`；API 元数据保存在 `zenodo15670547.json`。 |

检索时间为 2026-11-09。完整下载 SHA256、MD5、字节数、来源 URL 在 `modern/evidence_manifest.json`。2025/2026 数据集均标 CC-BY-4.0。仓储 PDF 被返回为 HTML 的失败下载应视为访问失败，不能视为全文证据。

## 2025：晶相、缺陷和计算条件

以下均由 2025 主文直接核对：

| 项目 | 原文可核实内容 | 位置 |
|---|---|---|
| 晶相 | tetragonal Li1Si1，No.88 I41/a；常规胞16Li+16Si；实验 a=b=9.354 Å、c=5.746 Å | Table1，期刊p3 / PDF第4页 |
| 超胞与缺陷 | 常规胞2×2×2；一个Li1空位；因此 Li127Si128，255原子、128个Li宿主位点、cv=1/128。raw XYZ逐帧元素计数独立确认127Li+128Si | §2.2，Table2；Zenodo systems/li1si1_li_vac_li1.xyz；三条raw NEB各16帧255原子 |
| 电荷 | 全部体系周期边界、总电荷中性。文中 Li+ 是离子性描述，不是给空位加电荷或背景电荷的设置 | §2.3，期刊p3 |
| 软件与基组 | CP2K Quickstep/GPW；DZVP-MOLOPT-SR-GTH；GTH-PBE赝势 | §2.3 |
| 泛函与取样 | PBE + DFT-D3；电子密度平面波cutoff350Ry，relative cutoff40Ry；Γ-only | §2.3 |
| SCF | Broyden mixing α=.2；Fermi–Dirac电子温度600K；SCF accuracy10^-6 | §2.3 |
| 晶胞处理 | 首先0K、1atm全优化含缺陷/无缺陷晶胞；按实验/计算无缺陷晶格参数比值缩放；最后在各自固定校正晶胞中松弛坐标。因此**缺陷胞不是严格冻结为同一无缺陷宿主晶胞** | §2.3 |
| MD | NVT；LiSi800K；1fs步长；Nosé–Hoover chain；1ps massive thermostat τ10fs，再1ps τ100fs，然后100ps global thermostat τ100fs用于分析 | §2.3 |
| NEB | CI-NEB，16 replicas，含全松弛且固定的初末端；路径来自AIMD观测，也存在理论构造路径 | §2.3与§3.1 |
| Em | 孤立Li空位single-atom vacancy机制；三条晶体学不同Li1→邻近Li1 vacancy路径 | §3.2，期刊p9 / PDF第10页；SI TableS10与FigS9(b)被主文引用 |
| Eform | LiSi vacancy约1.0eV。缺陷形成能另用无D3、TZVP基组+k网格重新优化；不能把这个不同协议的Eform直接并入本项目固定缺陷D验证 | Table5、§2.3、§3.2 |

精确含空位校正胞常数在 SI Table S2，尚未取得，不能用2×实验胞强行替代。路径映射使用**无缺陷**理想宿主，故使用 Table1 的2×常规胞是正文明确支持的；没有将该名义胞冒充含空位胞。

本项目与以上协议共同点为PBE、晶态LiSi、中性孤立Li空位；不同点有63对255原子、cv、PAW平面波对GPW/GTH、520eV对350Ry电子密度cutoff（两个cutoff不能直接比数值）、k4对Γ、D3、晶胞约束、电子smearing、静态速率模型对有限温MD。现有证据不能将某一个差异单独指定为全部势垒差或D差的原因。

## 可用的 D 数值和可比性

单位在各文表中为10^-3 Å²/ps = 10^-7 cm²/s。

| 年份/模型 | 温度 | LiSi数值 | 观测对象与位置 | 是否可用于本项目绝对D验证 |
|---|---:|---:|---|---|
| 2022 无非化学计量缺陷 | 500、800 K | 各方向未观测到扩散 | SI FigS1、TableS1脚注a，SI pS2 | 否；有限观测，不是测得非零D或严格上限 |
| 2025 无缺陷 | 800K | 0.0 | Table3；取自2022 | 否；同上 |
| 2025 一个Li间隙 | 800K | 0.0 | Table3 | 否；缺陷机制不同 |
| **2025 一个Li空位** | **800K** | **1.0 → 1.0×10^-7 cm²/s** | **Li原子MSD**；Table3；MSD图为SI FigS5，公式EqS5；Li127Si128 | **最接近的对照；仍需cv/示踪相关/温度/归一匹配** |
| 2025 Li+Si Schottky | 800K | 0.7 | Table3；另一Si空位不动、Li空位扩散 | 否；多了Si空位及局部捕获机制 |
| 2026 substitutional Na | 800K | 0.0 | Table3，Na_Li体系 | 否；Li空位已被Na替代 |
| 2026 添加Na间隙 | 800K | 0.5 | Table3；Na_i体系；作者发现Na间隙不稳定、可置换Li | 否；有Na及Li间隙/置换机制 |
| 2026 substitutional Ge | 800K | 0.0 | Table3，Ge_Si体系 | 否；有Ge掺杂 |

2026 raw structures独立计数：Na_i 为128Li+128Si+1Na=257；Na_Li为127Li+128Si+1Na=256；Ge_Si为128Li+127Si+1Ge=256。均不是Li127Si128孤立空位模型。

### 必须保留的 MSD 归一缺口

2025正文明确说从Li MSD根据EqS5求DLi。由于SI未取到，尚未逐字核验：

- 是否严格用 (1/(6N_Li)) Σ_i <|r_i(t+τ)−r_i(t)|²>/τ；N_Li是127还是另有归一。
- 时间原点平均、unwrapping、质心漂移、适用扩散维数、拟合lag区间、拟合方法。
- 是否给独立重复、block误差、置信区间或稀少跳跃误差。
- TableS7方向项定义与三方向相加规则。

2022 SI存在一个不能误读的细节：TableS1的Li12Si7(500K)方向项4.4+2.0+5.1≈11.5，对应FigS2总项≈11.6；Li13Si4只有Dx39.5，总项也是39.5。因此方向项可能已经各用1/6，是贡献项，而不能想当然认作惯用D_xx=MSD_x/(2t)。视觉交叉核对FigS1：Li22Si5(500K)在50ps的MSD约14Å²，斜率/6≈0.047Å²/ps，与FigS2总项≈0.048相符；斜率/2会大约三倍。因此**图形支持总DLi采用标准1/6，但还不是2025 EqS5的直接证据，不建议擅加一个×3修正**。

若最终EqS5采用标准总MSD/(6t)，与模型同口径比较主值1.0×10^-7；若它定义的是trace而非trace/3，标准标量对应3.33×10^-8。该第二分支目前只是未排除的定义风险，不能写成“作者有三倍错误”。Table3只保留一位小数且没有附误差；不能给任意接近程度设成绝对D“通过线”。

在通常的平均原子MSD假设下，一空位体系应注意 D_Li* ≈ f·(N_v/N_Li)·D_v，精确浓度比是1/127，而site fraction是1/128；稀释极限近似才可互换。f随跳跃网络、速率、温度及有限盒子而变，不应仅以D_Li/cv就宣称得到观测空位D。100ps是否足够也取决于实际独立跳跃次数；公开数据集未提供MD轨迹以重做统计。

## 原始 NEB 与严格轨道映射

`modern/map_literature_orbits.py` 可复现分析；结果为 `modern/literature_orbit_mapping.json`。全部源文件仅读取，代码只写本补充目录。未做逐路径参数拟合、没有调前因子来匹配结果。

映射步骤：

1. 用2025 `li1si1_no_defect_path_1.xyz` 首帧恢复无缺陷256原子宿主，以原文明确的18.708×18.708×11.492Å超胞描述。
2. 文献结构在spglib symprec=10^-4Å识别为88(I41/a)，10^-4、10^-3、10^-2Å结果稳定；10^-5Å因数值微扰降至13，不将该过严阈值误认为不同相。本项目投影宿主在10^-5Å为88。
3. 各自转为standard conventional32原子胞：本项目a=9.33740189、c=5.73414364Å；文献a=9.354、c=5.746Å。一个全局四方setting变化+原点平移，对齐16Li+16Si，双物种一对一映射RMS0.0178158Å、最大0.0239192Å。晶格相对差约a0.178%、c0.207%；没有逐route重缩放。
4. 从每条vacancy raw XYZ找最大移动Li，逐个端点映回无缺陷宿主最近Li位点，再将宿主位移转为standard分数坐标。比较本项目四轨道在实际4/m点群操作下的全部等价向量；也核对端点的全局宿主映射。进一步枚举空间群的旋转和平移，将每条边的**一对端点**整体匹配（只允许共同整数晶格平移和正反方向），不是独立给每个端点选镜像。三条正确轨道的端点对RMS分别0.0239433/0.0239407/0.0239404Å，错误轨道均不低于1.91146Å。脚本现在直接从原始XYZ和energy计算，不依赖手写的中间映射。

| 文献route | raw Li端点位移Å | 映回文献理想宿主位点距Å | 本项目route（几何确定） | 本项目宿主距Å | 归一分数向量轨道残差Å | 次近错误轨道残差Å |
|---|---:|---:|---|---:|---:|---:|
| (6)_def | 2.254919 | 2.716313 | **path00** | 2.694958 | .047887 | 1.392700 |
| (7)_def | 2.667652 | 2.907833 | **path02** | 2.946005 | .043764 | 1.334843 |
| (8)_def | 2.788540 | 2.926164 | **path01** | 2.928470 | .036537 | 1.368235 |

注意path01/path02必须交换；仅按距离数值从小到大排序会得错配。文献7几乎纯basal方向(z≈0)，文献8有显著c方向分量，足以区分这两条几何轨道。

同一个全局setting下文献6/7/8 standard fractional vectors为(-.23292750,-.15908895,.11233966)、(.15908915,.26707264,~0)、(-.21308077,-.05399180,-.36233781)。相应本项目等价向量为(-.22928098,-.16203392,.10895073)、(.16203392,.27071902,~0)、(-.21637647,-.05434255,-.35895073)。这不是只比较一个标量距离。

raw能量文件每条16行，单位按CP2K Hartree解释并以27.211386245988eV/Ha换算，与文章表述0.22/0.25/0.28eV一致。正反端点能量小差约0.00006–0.00009eV；对称化Em定义为Emax−(Einitial+Efinal)/2。

| 文献route → 本项目 | 文献正向Em eV | 文献反向Em eV | 文献对称化Em eV | 本项目精确对称化Em eV | 本项目−文献 eV |
|---|---:|---:|---:|---:|---:|
| 6 → path00 | .2202561915 | .2203155555 | .2202858735 | .2534438050 | +.0331579315 |
| 7 → path02 | .2526526367 | .2525678211 | .2526102289 | .2459094900 | −.0067007389 |
| 8 → path01 | .2849121448 | .2848308013 | .2848714731 | .2358302050 | −.0490412681 |
| 无对应 → path03 | — | — | — | .4276610950 | — |

几何对应置信度高：全宿主双物种bijective、端点残差约0.024Å、正确轨道与其他轨道分离超过1.3Å；仍未对所有鞍点图像作共同坐标系内的机制/力场比较。因此这里的“对应”是同一宿主跳跃轨道类别，不能扩大成所有鞍点构型相同或数值势垒已收敛到同一答案。

文献没有给这三条路径可直接移植的harmonic-TST ν或D0；单个800K AIMD点也不能求独立Ea和D0。本项目统一ν=10^13s^-1仍是假设。不能拟合ν以迫使800K数据吻合，再把吻合当验证。Em也不等于任意机制或热平衡本征扩散的Ea；后者可含形成能、相关因子和机制切换。

## 缺失、风险与下一步

- 最关键未完成项是2025 SI EqS5、TableS2、TableS7、FigS5和TableS10原始页面。已尝试DOI HTML/PDF links、IOP data/suppdata/mmedia、stacks重定向、Crossref、常见媒体候选；网页常遭robots/Radware/Anubis，候选SI404/403。PDF全部get_links()没有SI外链。浏览器工具没有可用browser，未尝试破解验证码。
- 后续公开OAI-PMH进一步给出明确仓储清单：`https://opendata.uni-halle.de/oai/request?verb=GetRecord&identifier=oai:opendata.uni-halle.de:1981185920/120481&metadataPrefix=xoai`（保存`modern/halle2025_xoai.xml`）中，ORIGINAL bundle只有主文PDF，另有TEXT和缩略图，没有SI；2026 item123640亦同。不能再把该仓储作为“已包含SI”的来源。
- 2025 Zenodo记录只有systems+四种化学计量paths zip；LiSi zip含53条目录/文件项：XYZ+energy。没有可重算AIMD MSD的轨迹、输入、fit或误差数据，不能许诺公开原始数据已允许绝对D复算。
- 本项目63原子有限模型收敛不等于255原子、稀释极限或宏观材料D收敛；晶胞/基组/kmesh/D3/动力学差异未被单因素隔离，归因应写“候选差异来源”，不能归结为某一项。
- 主线程可做800K、cv=1/128的同网络计算和有限盒子Li示踪，再给通常归一和trace定义分支；保留统计和方法缺口。也可只换三条对应Em做单因素敏感性，明确这是模型敏感性而非独立实验验证。
- 若后续取得SI/MD，先核公式和单位，再查LiSi一空位轨迹是否稳定线性、跳跃数/回跳、分块置信区间，以及D_v与D_Li*关联；再判定是否值得运行更大晶胞DFT/MD。此处没有提出或启动计算任务。

### 核心原始文件 SHA256

| 文件 | SHA256 |
|---|---|
| IOP2025.pdf | 083f47da4555c3779735fe49beb4a962ada3fa80478e81afcf33676a8ce9dce4 |
| ACS2022_SI.pdf | 9f1701595cae4ad19f85cfdd9066521452584a0dbff411665886ee86b86f7252 |
| 2025_paths_li1si1.zip | 3cdb80e55658ac145a893835427871071827b6c0e919ca20c37542ad6dd60b8e |
| 2025_systems.zip | 602ab80b19520a7cb147353588009b2a3f20738a90f570682c4051b3929334ba |
| 2026_paths_li1si1.zip | 8ef65bd83dafa712196f2fde039813716e3fffcd69d754e5e0b89156a7ff987c |
| 2026_systems.zip | 477fbcb68710e478de3ff5eadd898307f1f5a3654435153802c854d4b2ce274e |

源项目POSCAR/轨道报告/精确势垒证据及文献原始无缺陷首帧的SHA均保存在`modern/literature_orbit_mapping.json`，便于后续复核。

## 2026-11-09 收尾更新

- 已完成：真实PDF页数修正为15；追加逐边端点对的完整空间群匹配、精确项目势垒读取、10个映射源文件SHA；下载zip的MD5与Zenodo/API记录吻合。
- 结论保持：NEB几何轨道对应高置信；绝对D仅条件式比较；S5统计口径与误差缺口未闭合。不将800K单点倒推D0/Ea，不声称尺寸收敛。
- 原归档文件：`modern_findings.md`、`modern/literature_orbit_mapping.json`、`modern/map_literature_orbits.py`、`modern/evidence_manifest.json`。项目另行完成的 KMC 浓度、温度及 Em 替换敏感性属于模型控制，不作为文献数据。
