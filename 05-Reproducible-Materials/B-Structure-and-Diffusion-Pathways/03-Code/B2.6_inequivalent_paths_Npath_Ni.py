import numpy as np

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


# ============================================================
# 参数
# ============================================================

SYMPREC = 0.10
ANGLE_TOLERANCE = 5.0

D0 = 2.7614055281      # LiAl，来自 B2.4，Å
EPS = 1e-7

# 空间群操作后，判断两个端点是否重合的笛卡尔距离容差
# 与本项目空间群识别 symprec 保持一致
MATCH_TOL = 0.10       # Å


# ============================================================
# B2.5：无方向物理 hop canonical 化
# ============================================================

def canonical_hop(i, j, image):

    rx, ry, rz = (int(x) for x in image)

    forward = (
        int(i), int(j),
        rx, ry, rz
    )

    reverse = (
        int(j), int(i),
        -rx, -ry, -rz
    )

    return min(forward, reverse)


def get_physical_hops(structure, d0):

    query_cutoff = d0 + EPS

    li_indices = [
        i for i, site in enumerate(structure)
        if site.species_string == "Li"
    ]

    unique_hops = {}

    for i in li_indices:

        center = structure[i]

        neighbors = structure.get_neighbors(
            center,
            query_cutoff
        )

        for neighbor in neighbors:

            if neighbor.species_string != "Li":
                continue

            j = int(neighbor.index)
            image = tuple(
                int(x) for x in neighbor.image
            )

            # 排除真正的自己
            if (
                i == j
                and image == (0, 0, 0)
            ):
                continue

            key = canonical_hop(
                i,
                j,
                image
            )

            if key not in unique_hops:
                unique_hops[key] = float(
                    neighbor.nn_distance
                )

    return [
        {
            "key": key,
            "distance": distance
        }
        for key, distance
        in sorted(unique_hops.items())
    ]


# ============================================================
# B2.6：把压缩 hop 恢复成两个周期端点
# ============================================================

def hop_endpoints(structure, hop):

    i, j, rx, ry, rz = hop["key"]

    # 起点：参考晶胞中的 Li_i
    p1 = np.array(
        structure[i].frac_coords,
        dtype=float
    )

    # 终点：Li_j 的指定周期像
    p2 = (
        np.array(
            structure[j].frac_coords,
            dtype=float
        )
        +
        np.array(
            [rx, ry, rz],
            dtype=float
        )
    )

    return p1, p2


# ============================================================
# 判断两对端点是否只相差一个共同整数晶格平移
# ============================================================

def endpoints_match(
    structure,
    p1,
    p2,
    q1,
    q2,
):

    lattice = structure.lattice

    # 如果两条 hop 是同一条周期连接，
    # 两个端点之间应该相差同一个整数晶格平移 T。
    T = np.rint(p1 - q1)

    err1_frac = p1 - (q1 + T)
    err2_frac = p2 - (q2 + T)

    err1_cart = lattice.get_cartesian_coords(
        err1_frac
    )

    err2_cart = lattice.get_cartesian_coords(
        err2_frac
    )

    err1 = np.linalg.norm(err1_cart)
    err2 = np.linalg.norm(err2_cart)

    return (
        err1 <= MATCH_TOL
        and
        err2 <= MATCH_TOL
    )


# ============================================================
# 判断两条物理 hop 是否空间群对称等价
# ============================================================

def hops_are_symmetry_equivalent(
    structure,
    hop1,
    hop2,
    symmetry_ops,
):

    a, b = hop_endpoints(
        structure,
        hop1
    )

    c, d = hop_endpoints(
        structure,
        hop2
    )

    for op_index, op in enumerate(symmetry_ops):

        # 同一个空间群操作同时作用在 hop1 两个端点
        a_new = np.array(
            op.operate(a),
            dtype=float
        )

        b_new = np.array(
            op.operate(b),
            dtype=float
        )

        # 情况1：
        # A -> C, B -> D
        if endpoints_match(
            structure,
            a_new,
            b_new,
            c,
            d,
        ):
            return True, op_index, False

        # 情况2：
        # 因为 hop 无方向，允许
        # A -> D, B -> C
        if endpoints_match(
            structure,
            a_new,
            b_new,
            d,
            c,
        ):
            return True, op_index, True

    return False, None, None


# ============================================================
# 空间群分类
# ============================================================

def classify_hops(
    structure,
    hops,
    symmetry_ops,
):

    path_classes = []

    for hop_id, hop in enumerate(hops):

        assigned = False

        for path_class in path_classes:

            representative_id = path_class[0]
            representative = hops[
                representative_id
            ]

            equivalent, _, _ = (
                hops_are_symmetry_equivalent(
                    structure,
                    hop,
                    representative,
                    symmetry_ops,
                )
            )

            if equivalent:
                path_class.append(hop_id)
                assigned = True
                break

        # 与已有类别都不等价
        # → 新建一个路径类别
        if not assigned:
            path_classes.append([hop_id])

    return path_classes


# ============================================================
# 主程序
# ============================================================

def analyze_path_classes(cif_path):

    structure = Structure.from_file(
        cif_path
    )

    sga = SpacegroupAnalyzer(
        structure,
        symprec=SYMPREC,
        angle_tolerance=ANGLE_TOLERANCE
    )

    symmetry_ops = (
        sga.get_symmetry_operations(
            cartesian=False
        )
    )

    print("=" * 70)
    print("文件:", cif_path)
    print(
        "化学式:",
        structure.composition.reduced_formula
    )
    print(
        "空间群:",
        sga.get_space_group_symbol(),
        sga.get_space_group_number()
    )
    print(
        "空间群操作数:",
        len(symmetry_ops)
    )
    print(f"d0 = {D0:.10f} Å")

    # --------------------------------------------------------
    # 1. B2.5：重新得到 N 条物理 hop
    # --------------------------------------------------------

    hops = get_physical_hops(
        structure,
        D0
    )

    print("\n=== B2.5 物理 hops ===")

    for hop_id, hop in enumerate(hops):

        print(
            f"hop{hop_id}: "
            f"{hop['key']} | "
            f"d={hop['distance']:.10f} Å"
        )

    N = len(hops)

    # --------------------------------------------------------
    # 2. B2.6：空间群分类
    # --------------------------------------------------------

    path_classes = classify_hops(
        structure,
        hops,
        symmetry_ops
    )

    print("\n=== B2.6 对称等价分类 ===")

    Ni_list = []

    for path_id, hop_ids in enumerate(
        path_classes,
        start=1
    ):

        Ni = len(hop_ids)
        Ni_list.append(Ni)

        distances = [
            hops[i]["distance"]
            for i in hop_ids
        ]

        print(
            f"Path {path_id}: "
            f"hop_ids={hop_ids} | "
            f"Ni={Ni}"
        )

        print(
            "    distances =",
            [
                f"{d:.10f}"
                for d in distances
            ]
        )

        print(
            f"    distance spread = "
            f"{max(distances) - min(distances):.10e} Å"
        )

    Npath = len(path_classes)

    # --------------------------------------------------------
    # 3. 最终结果
    # --------------------------------------------------------

    print("\n=== B2.6 结果 ===")
    print("N =", N)
    print("Npath =", Npath)
    print("Ni =", Ni_list)

    print(
        "sum(Ni) =",
        sum(Ni_list)
    )

    print(
        "守恒检查 sum(Ni) == N:",
        sum(Ni_list) == N
    )


# 第一轮：LiAl
analyze_path_classes(
    "cif_files/mp-16506_Li3Al2.cif"
)