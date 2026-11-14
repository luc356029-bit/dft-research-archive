# -*- coding: utf-8 -*-
"""C1.2.2 快速读取原始 VASP 输入/结果；只读，不运行 DFT。"""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
FOLDERS = [
    "01_完整晶胞16原子_全弛豫",
    "02_完整超胞64原子_单点检查",
    "03_空位63原子_固定晶格弛豫",
    "04_NEB末态B00_单点验收",
]

def poscar(folder, filename="POSCAR"):
    lines = (folder / filename).read_text(encoding="utf-8").splitlines()
    scale = float(lines[1].strip())
    lattice = [[float(x) * scale for x in line.split()[:3]] for line in lines[2:5]]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    return dict(species=dict(zip(species, counts)), n=sum(counts), lattice=lattice)

def incar(folder):
    params = {}
    for line in (folder / "INCAR").read_text(encoding="utf-8").splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            params[k.strip().upper()] = v.split("!")[0].strip()
    return params

def max_cell_diff(a, b):
    return max(abs(a[i][j] - b[i][j]) for i in range(3) for j in range(3))

parsed = []
for name in FOLDERS:
    folder = ROOT / name
    p, q = poscar(folder), poscar(folder, "CONTCAR")
    inc = incar(folder)
    res = json.loads((folder / "analysis_result.json").read_text(encoding="utf-8"))
    parsed.append(p)
    fmax = res.get("last_frame", {}).get("Fmax_eV_A")
    print("\n", name)
    print("组成（输入）", p['species'], "; 原子总数", p['n'])
    print("组成（最终）", q['species'], "; 原子总数", q['n'])
    print("ISIF", inc.get('ISIF'), "NSW", inc.get('NSW'))
    print("输入到输出：晶格矢量最大变化", f"{max_cell_diff(p['lattice'], q['lattice']):.9f}", 'Å')
    print("SCF收敛", res.get('SCF_converged'), "离子评价帧数", res.get('n_ionic_frames'), "最终最大原子受力", round(fmax, 8) if fmax is not None else None, 'eV/Å')
print("\n完整64与缺陷63：输入晶格矢量最大差", f"{max_cell_diff(parsed[1]['lattice'], parsed[2]['lattice']):.9f}", 'Å')
print("缺陷63与末态B00：输入晶格矢量最大差", f"{max_cell_diff(parsed[2]['lattice'], parsed[3]['lattice']):.9f}", 'Å')
