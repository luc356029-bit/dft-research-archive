from pymatgen.core import Structure
from pymatgen.core.graphs import StructureGraph
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.dimensionality import get_dimensionality_larsen


SYMPREC = 0.10
ANGLE_TOLERANCE = 5.0

SEARCH_RADIUS = 6.0
SHELL_TOL = 0.01
EPS = 1e-7


def get_li_groups(structure):
    """
    得到对称不等价 Li 分组。
    继承 B2.2。
    """
    sga = SpacegroupAnalyzer(
        structure,
        symprec=SYMPREC,
        angle_tolerance=ANGLE_TOLERANCE
    )

    symm_structure = sga.get_symmetrized_structure()

    li_groups = []

    for group in symm_structure.equivalent_indices:
        group = [int(i) for i in group]

        if structure[group[0]].species_string == "Li":
            assert all(
                structure[i].species_string == "Li"
                for i in group
            )
            li_groups.append(group)

    return li_groups


def get_distance_shells(structure, li_groups):
    """
    B2.3：
    从每组不等价 Li 的代表位点收集 Li-Li 距离，
    再按 0.01 Å 分壳。
    """
    distances = []

    for group in li_groups:
        rep_index = group[0]
        center = structure[rep_index]

        neighbors = structure.get_neighbors(
            center,
            SEARCH_RADIUS
        )

        for neighbor in neighbors:
            if neighbor.species_string != "Li":
                continue

            image = tuple(int(x) for x in neighbor.image)

            # 排除中心位点自己
            if (
                int(neighbor.index) == rep_index
                and image == (0, 0, 0)
            ):
                continue

            distances.append(float(neighbor.nn_distance))

    distances.sort()

    shells = []
    remaining = distances.copy()

    while remaining:
        dmin = remaining[0]

        current_shell = [
            d for d in remaining
            if d - dmin <= SHELL_TOL
        ]

        shells.append(current_shell)
        remaining = remaining[len(current_shell):]

    return shells


def build_li_graph(structure, cutoff):
    """
    使用当前 cutoff 建完整 Li 周期图。

    注意：
    这里不能只从“不等价 Li 代表”建边，
    因为维度判断需要完整周期 Li 网络。
    """
    graph = StructureGraph.from_empty_graph(
        structure,
        name="Li candidate diffusion network",
        edge_weight_name="distance",
        edge_weight_units="Å"
    )

    li_indices = [
        i for i, site in enumerate(structure)
        if site.species_string == "Li"
    ]

    for i in li_indices:
        center = structure[i]

        neighbors = structure.get_neighbors(
            center,
            cutoff
        )

        for neighbor in neighbors:
            if neighbor.species_string != "Li":
                continue

            j = int(neighbor.index)
            image = tuple(int(x) for x in neighbor.image)
            distance = float(neighbor.nn_distance)

            # 排除零距离自身
            if (
                j == i
                and image == (0, 0, 0)
            ):
                continue

            graph.add_edge(
                from_index=i,
                to_index=j,
                from_jimage=(0, 0, 0),
                to_jimage=image,
                weight=distance,
                warn_duplicates=False
            )

    return graph


def analyze_d0(cif_path):
    structure = Structure.from_file(cif_path)

    print("=" * 70)
    print("文件:", cif_path)
    print("化学式:", structure.composition.reduced_formula)

    # 1. B2.2：不等价 Li
    li_groups = get_li_groups(structure)

    print("\n=== Li 等价组 ===")
    for n, group in enumerate(li_groups):
        print(
            f"Li{n}: {group}, "
            f"代表位点={group[0]}"
        )

    # 2. B2.3：距离壳
    shells = get_distance_shells(
        structure,
        li_groups
    )

    print("\n=== 逐距离壳检查长程连通 ===")

    d0 = None
    d0_dimension = None
    d0_shell = None

    # 3. B2.4：逐壳测试
    for shell_number, shell in enumerate(shells, start=1):

        dmin = min(shell)
        dmax = max(shell)

        # d0 本身定义为 dmax
        # EPS 仅用于实际数值查询
        query_cutoff = dmax + EPS

        graph = build_li_graph(
            structure,
            query_cutoff
        )

        dimension = get_dimensionality_larsen(graph)

        print(
            f"Shell {shell_number}: "
            f"{dmin:.10f} ~ {dmax:.10f} Å | "
            f"cutoff={dmax:.10f} Å | "
            f"Li边数={graph.graph.number_of_edges()} | "
            f"dimension={dimension}"
        )

        # 第一次出现长程周期连通
        if dimension >= 1:
            d0 = dmax
            d0_dimension = dimension
            d0_shell = shell_number
            

    print("\n=== B2.4 结果 ===")

    if d0 is None:
        print(
            f"在 {SEARCH_RADIUS} Å 搜索范围内 "
            "尚未找到长程 Li 网络。"
        )
        print("下一步应扩大 SEARCH_RADIUS。")

    else:
        print("首次长程连通壳:", d0_shell)
        print(f"d0 = {d0:.10f} Å")
        print("dimension(d0) =", d0_dimension)


# 第二轮：Li2Al
analyze_d0(
    "cif_files/mp-16506_Li3Al2.cif"
)