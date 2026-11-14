"""只读复算本教学包path 0的能量；不运行DFT，不修改任何文件。"""
from pathlib import Path
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
PATTERN = re.compile(r"free\s+energy\s+TOTEN\s*=\s*(" + FLOAT + r")")

def read_ionic_energies(path: Path) -> list[float]:
    """只取离子能量汇总区，避免把SCF内部的电子迭代当成路径优化帧。"""
    text = path.read_text(encoding="utf-8", errors="replace")
    blocks = text.split("FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)")[1:]
    energies = []
    for block in blocks:
        match = PATTERN.search(block)
        if match is None:
            raise ValueError(f"{path}: 离子能量区不完整")
        energies.append(float(match.group(1).replace("D", "E").replace("d", "e")))
    if not energies:
        raise ValueError(f"{path}: 未找到完整离子能量汇总")
    return energies

def main() -> None:
    root = Path(__file__).resolve().parent
    index = json.loads((root / "源文件索引.json").read_text(encoding="utf-8"))
    for entry in index["files"]:
        p = root / entry["practice_path"]
        if hashlib.sha256(p.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"文件已变化或不完整：{p}；请保留原始副本")
    print(f"原始复制文件校验通过：{len(index['files'])}个。")
    print("口径：每个文件最后一个完整离子评价帧的Efr/TOTEN；单位eV。")
    energies: list[float] = []
    frames: list[int] = []
    for i in range(6):
        folder = root / "path0_CI" / f"{i:02d}"
        sequence = read_ionic_energies(folder / "OUTCAR")
        energies.append(sequence[-1]); frames.append(len(sequence))
        if 1 <= i <= 4:
            tree = ET.parse(folder / "vasprun.xml").getroot()
            calculations = tree.findall("calculation")
            if len(calculations) != len(sequence):
                raise ValueError(f"Image {i}: OUTCAR与XML帧数不一致")
            element = calculations[-1].find("energy/i[@name='e_fr_energy']")
            if element is None or element.text is None:
                raise ValueError(f"Image {i}: XML末帧缺少能量")
            if abs(float(element.text) - sequence[-1]) > 1e-7:
                raise ValueError(f"Image {i}: OUTCAR与XML末帧能量不一致")
    if len(set(frames[1:5])) != 1:
        raise ValueError("中间图像优化帧数不一致；不能直接拼成共同末帧")
    print("目录  已评价帧数    最后Efr/eV        相对初态/eV")
    for i, (e, n) in enumerate(zip(energies, frames)):
        print(f"{i:02d}    {n:4d}       {e:15.8f}     {e-energies[0]:.8f}")
    print("注：00、05的帧数属于独立端态计算，不应与01—04的NEB帧数比较。")
    peak = max(range(1, 5), key=lambda i: energies[i])
    ep, ea, eb = energies[peak], energies[0], energies[-1]
    print(f"\n最高中间图像：{peak:02d}")
    print(f"正向 E峰-EA = {ep-ea:.8f} eV")
    print(f"反向 E峰-EB = {ep-eb:.8f} eV")
    print(f"端态能差 EB-EA = {eb-ea:.8f} eV")
    print(f"本轮等价端态对称化 E峰-(EA+EB)/2 = {ep-(ea+eb)/2:.9f} eV")
    print("\n本脚本没有独立重做路径寻优、SCF收敛审计或NEB投影力重构。")
    print("最高Image是否可作为鞍点，仍须结合LCLIMB、原始收敛记录和路径核验。")

if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, IndexError, ET.ParseError) as exc:
        print(f"检查停止：{exc}", file=sys.stderr)
        sys.exit(1)
