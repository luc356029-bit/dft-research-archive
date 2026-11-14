from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# 1. 读取晶体结构
structure = Structure.from_file("cif_files/mp-1067_LiAl.cif")

# 2. 看看 Structure 里实际有什么
print("=== Structure 输入信息 ===")
print("化学式:", structure.composition.reduced_formula)
print("晶格常数 a,b,c:", structure.lattice.abc)
print("晶格角 α,β,γ:", structure.lattice.angles)
print("原子总数:", len(structure))

print("\n原子与分数坐标:")
for i, site in enumerate(structure):
    print(i, site.species_string, site.frac_coords)

# 3. 对称性分析
sga = SpacegroupAnalyzer(
    structure,
    symprec=0.1,
    angle_tolerance=5.0
)

# 4. 输出空间群
print("\n=== 空间群分析 ===")
print("空间群符号:", sga.get_space_group_symbol())
print("空间群编号:", sga.get_space_group_number())