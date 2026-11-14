from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

SYMPREC = 0.10
ANGLE_TOLERANCE = 5.0

SEARCH_RADIUS = 9.0   # Å
SHELL_TOL = 0.01      # Å
EPS = 1e-7            # Å


def analyze_distance_shells(cif_path):
    # 1. 读取结构
    structure = Structure.from_file(cif_path)

    # 2. 得到对称等价位点
    sga = SpacegroupAnalyzer(
        structure,
        symprec=SYMPREC,
        angle_tolerance=ANGLE_TOLERANCE
    )

    symm_structure = sga.get_symmetrized_structure()

    # 3. 提取 Li 等价组
    li_groups = []

    for group in symm_structure.equivalent_indices:
        group = [int(i) for i in group]

        if structure[group[0]].species_string == "Li":
            assert all(
                structure[i].species_string == "Li"
                for i in group
            )
            li_groups.append(group)

    print("=" * 70)
    print("文件:", cif_path)
    print("化学式:", structure.composition.reduced_formula)
    print("空间群:",
          sga.get_space_group_symbol(),
          sga.get_space_group_number())

    print("\n=== Li 等价组 ===")
    for n, group in enumerate(li_groups):
        print(f"Li{n}: {group}，代表位点 = {group[0]}")

    # 4. 从每一个不等价 Li 代表位点出发，
    #    搜索 SEARCH_RADIUS 内全部 Li 周期邻居
    all_distances = []

    print("\n=== 各代表 Li 的周期 Li 邻居 ===")

    for li_type, group in enumerate(li_groups):
        representative_index = group[0]
        center_site = structure[representative_index]

        neighbors = structure.get_neighbors(
            center_site,
            SEARCH_RADIUS
        )

        print(
            f"\nLi{li_type} 代表位点 "
            f"index={representative_index}"
        )

        for neighbor in neighbors:
            # 只保留 Li 邻居
            if neighbor.species_string != "Li":
                continue

            distance = float(neighbor.nn_distance)

            # 去掉自己本身：
            # index相同且image=(0,0,0)
            image = tuple(int(x) for x in neighbor.image)

            if (
                neighbor.index == representative_index
                and image == (0, 0, 0)
            ):
                continue

            all_distances.append(distance)

            print(
                f"邻居index={neighbor.index}, "
                f"image={image}, "
                f"distance={distance:.10f} Å"
            )

    # 5. 所有距离从小到大排序
    all_distances.sort()

    print("\n=== 全部 Li-Li 距离（排序后） ===")
    for d in all_distances:
        print(f"{d:.10f}")

    # 6. 按 dmax-dmin <= 0.01 Å 聚成距离壳
    shells = []

    remaining = all_distances.copy()

    while remaining:
        dmin = remaining[0]

        current_shell = [
            d for d in remaining
            if d - dmin <= SHELL_TOL
        ]

        shells.append(current_shell)

        # 删除已经归入当前壳的距离
        remaining = remaining[len(current_shell):]

    # 7. 输出距离壳
    print("\n=== Li-Li 距离壳 ===")

    for i, shell in enumerate(shells, start=1):
        dmin = min(shell)
        dmax = max(shell)
        center = sum(shell) / len(shell)
        cutoff = dmax + EPS

        print(f"\nShell {i}")
        print(f"成员数: {len(shell)}")
        print(f"dmin   = {dmin:.10f} Å")
        print(f"dmax   = {dmax:.10f} Å")
        print(f"center = {center:.10f} Å")
        print(f"cutoff = {cutoff:.10f} Å")


# 第一轮：LiAl
analyze_distance_shells(
    "cif_files/mp-1210753_Li2Al.cif"
)