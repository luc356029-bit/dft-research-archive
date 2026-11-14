from pymatgen.core import Structure

# =========================
# 参数
# =========================
D0 = 2.7203208768   # LiAl，来自 B2.4，Å
EPS = 1e-7


def canonical_hop(i, j, image):
    """
    将同一无方向物理 hop：
        (i, j, R)
        (j, i, -R)
    统一成唯一 canonical key。
    """

    rx, ry, rz = (int(x) for x in image)

    forward = (
        int(i), int(j),
        rx, ry, rz
    )

    reverse = (
        int(j), int(i),
        -rx, -ry, -rz
    )

    # 固定保留字典序较小的一种表示
    return min(forward, reverse)


def analyze_hops(cif_path):
    structure = Structure.from_file(cif_path)

    query_cutoff = D0 + EPS

    print("=" * 70)
    print("文件:", cif_path)
    print("化学式:", structure.composition.reduced_formula)
    print(f"d0 = {D0:.10f} Å")
    print(f"程序搜索 cutoff = {query_cutoff:.10f} Å")

    # 当前晶胞内全部 Li
    li_indices = [
        i for i, site in enumerate(structure)
        if site.species_string == "Li"
    ]

    print("\n晶胞内全部 Li indices:", li_indices)

    raw_hops = []
    unique_hops = {}

    # ==================================================
    # 1. 从当前晶胞内每一个 Li 出发，枚举 d0 内 Li 邻居
    # ==================================================
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
            image = tuple(int(x) for x in neighbor.image)
            distance = float(neighbor.nn_distance)

            # 排除真正的自己 → 自己
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

            raw_hops.append(
                {
                    "i": i,
                    "j": j,
                    "image": image,
                    "distance": distance,
                    "key": key
                }
            )

            # canonical key 相同 → 同一物理 hop
            if key not in unique_hops:
                unique_hops[key] = distance

    # ==================================================
    # 2. 打印原始枚举
    # ==================================================
    print("\n=== 原始有方向 hop ===")

    for n, hop in enumerate(raw_hops, start=1):
        print(
            f"{n}: "
            f"Li{hop['i']} (0,0,0) -> "
            f"Li{hop['j']} {hop['image']} | "
            f"d={hop['distance']:.10f} Å | "
            f"canonical={hop['key']}"
        )

    # ==================================================
    # 3. 打印 canonical 去重后的物理 hop
    # ==================================================
    print("\n=== 去重后的无方向物理 hop ===")

    for n, key in enumerate(
        sorted(unique_hops),
        start=1
    ):
        i, j, rx, ry, rz = key

        print(
            f"{n}: "
            f"(Li{i}, Li{j}, "
            f"R=({rx},{ry},{rz})) | "
            f"d={unique_hops[key]:.10f} Å"
        )

    raw_N = len(raw_hops)
    N = len(unique_hops)

    # ==================================================
    # 4. 基本合理性检查
    # ==================================================
    print("\n=== B2.5 结果 ===")
    print("原始有方向 hop 数 =", raw_N)
    print("去重后物理 hop 数 N =", N)

    print(
        "正反向成对检查:",
        raw_N == 2 * N
    )


analyze_hops(
    "cif_files/mp-1210753_Li2Al.cif"
)