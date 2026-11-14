# prepare_supercell.py：16原子晶胞如何扩展为64原子超胞

整理日期：2026-10-31。复用首轮已有结构与脚本；本次只做几何复核，没有新跑DFT。

输入是 analysis/host_symmetry_projection/POSCAR，即已作对称性投影的16原子母胞（8Li+8Si），不是直接把任意原始下载晶胞重复。脚本从命令行读取输入路径，输出固定到 analysis/host_supercell64。

ASE将晶格矢量按行存储。采用整数矩阵

```text
P = [[ 0, 1,-1],
     [-1, 1, 1],
     [ 2, 0, 0]]
A_super = P @ A_parent
```

因此新基矢分别为 b−c、−a+b+c、2a；det(P)=4，体积扩大4倍，16×4=64原子（32Li+32Si）。这是非对角整数扩胞，不是2×2×2重复；后者对16原子输入会得到128原子。

脚本先枚举母胞空间群操作，检查 inv(P).T @ R @ P.T 是否为整数矩阵，以确保所需点操作与超胞周期边界兼容。随后用 make_supercell(..., order='atom-major') 扩胞，使每个母胞原子的4个周期像相邻排列；POSCAR写出时 sort=False，保留身份顺序。

supercell_mapping.json 保存每个超胞槽位对应的 basis_indices（母胞原子索引）及 primitive_images（母胞整数周期像）。可用 r_super=(f_basis+n) @ A_parent 逐原子追溯，而不是只靠坐标近邻猜测身份。脚本检验整数像残差小于1e-8。

已有输出最短周期平移为9.337401893572688 Å（原脚本在整数平移[-3,3]^3内搜索），满足脚本采用的>9 Å几何门槛。这个门槛和扩胞矩阵是本项目选择，不是论文明确给出的LiSi参数，也不能证明缺陷相互作用已达到尺寸收敛。删除一个Li后是Li31Si32，共63原子，cv=1/32；64是完整宿主原子数，不是缺陷体系原子数。

原生成脚本输出状态仍写 REQUIRES_SUPERCELL_DFT_GRID_CHECK，属于生成当时的状态；后续DFT验收见首轮主报告。本次在内存重新执行同一ASE扩胞并逐坐标比对已有POSCAR，误差<1e-10 Å。没有覆盖旧输出目录。

附件：原脚本 prepare_supercell.py；共同证据ZIP中的母胞/超胞POSCAR、supercell_mapping.json、依赖脚本和verification.json。原脚本会因既有输出目录而拒绝重建，应在独立副本运行，避免改动首轮证据。
