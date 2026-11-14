from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def analyze_li_equivalence(cif_path):
    # 1. 读取结构
    structure = Structure.from_file(cif_path)

    # 2. 空间群分析
    sga = SpacegroupAnalyzer(
        structure,
        symprec=0.10,
        angle_tolerance=5.0
    )

    # 3. 得到按对称性分组后的结构
    symm_structure = sga.get_symmetrized_structure()

    print("=" * 60)
    print("文件:", cif_path)
    print("化学式:", structure.composition.reduced_formula)
    print("空间群:", sga.get_space_group_symbol(),
          sga.get_space_group_number())

    print("\n=== 当前 Structure 中全部位点 ===")
    for i, site in enumerate(structure):
        print(i, site.species_string, site.frac_coords)

    print("\n=== equivalent_indices 原始结果 ===")
    print(symm_structure.equivalent_indices)

    # 4. 从所有等价组中筛选 Li
    li_groups = []

    for group in symm_structure.equivalent_indices:
        group = [int(i) for i in group]

        if structure[group[0]].species_string == "Li":
            # 防御性检查：同一等价组必须全部都是 Li
            assert all(
                structure[i].species_string == "Li"
                for i in group
            )
            li_groups.append(group)

    print("\n=== Li 等价组 ===")
    for n, group in enumerate(li_groups):
        print(f"Li{n}: {group}，代表位点 = {group[0]}")

    # 5. NLi
    NLi = len(li_groups)

    print("\nNLi =", NLi)

    # 6. 覆盖检查
    all_li_indices = [
        i for i, site in enumerate(structure)
        if site.species_string == "Li"
    ]

    grouped_li_indices = sorted(
        i for group in li_groups for i in group
    )

    print("全部 Li indices:", all_li_indices)
    print("分组覆盖 indices:", grouped_li_indices)
    print("完整覆盖:", sorted(all_li_indices) == grouped_li_indices)


# 第一轮：LiAl
analyze_li_equivalence("cif_files/mp-1210753_Li2Al.cif")