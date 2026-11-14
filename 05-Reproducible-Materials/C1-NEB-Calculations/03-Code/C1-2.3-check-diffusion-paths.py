# -*- coding: utf-8 -*-
"""C1.2.3 实践：从真实宿主POSCAR独立复算几何边，再核对原始分类和连通报告。

仅Python标准库；不进行DFT，也不将几何候选视为已验证物理跳跃。

Windows: py .\核查路径.py
连通检查: py .\核查路径.py --network
"""
import sys, json, math, hashlib
from pathlib import Path
from itertools import product

ROOT=Path(__file__).resolve().parent
POSCAR=ROOT/'01_16原子宿主结构'/'POSCAR'
REPORT=ROOT/'02_候选路径分类'/'projection_report.json'
TOPOLOGY=ROOT/'03_周期网络连通性'/'complete_orbit_connectivity_INTEGER.json'

def load_poscar(path):
    lines=path.read_text(encoding='utf-8').splitlines()
    scale=float(lines[1].split()[0])
    A=[[float(x)*scale for x in lines[i].split()[:3]] for i in range(2,5)]
    names=lines[5].split(); nums=[int(v) for v in lines[6].split()]
    assert names[:2]==['Li','Si'] and nums[:2]==[8,8] and lines[7].lower().startswith('direct'), 'POSCAR结构格式与预期不符'
    coords=[[float(x) for x in lines[i].split()[:3]] for i in range(8,8+nums[0]+nums[1])]
    return A,coords[:nums[0]],names,nums

def canonical(i,j,R):
    fwd=(i,j,*R); rev=(j,i,*(-x for x in R))
    return min(fwd,rev)

def distance(i,j,R,A,frac):
    delta=[frac[j][k]+R[k]-frac[i][k] for k in range(3)]
    cart=[sum(delta[k]*A[k][c] for k in range(3)) for c in range(3)]
    return math.sqrt(sum(v*v for v in cart))

def enumerate_edges(A,frac,rmax=4.31):
    # 仅为本16原子LiSi胞的示例：在[-2,2]^3显式枚举周期像，并核对原档的3.3A边集。
    all_edges={}
    for i in range(len(frac)):
        for j in range(len(frac)):
            for R in product(range(-2,3),repeat=3):
                if i==j and R==(0,0,0):
                    continue
                d=distance(i,j,R,A,frac)
                if d<rmax:
                    all_edges[canonical(i,j,R)]=d
    return all_edges

def main():
    A,frac,names,nums=load_poscar(POSCAR)
    report=json.loads(REPORT.read_text(encoding='utf-8'))
    orbit=report['candidate_orbits']
    all_edges=enumerate_edges(A,frac)
    geom={e:d for e,d in all_edges.items() if d<3.3}
    saved_edges={tuple(e) for o in orbit['orbits'] for e in o['members']}
    print('== C1.2.3 第一小步：几何路径和等价归类 ==')
    print('晶胞原子数量：',dict(zip(names,nums)), '；Li位点索引为0—7（从0开始）')
    print('实际参考结构空间群：',report['projected_spacegroup_at_1e5_A'])
    print('截断距离：',orbit['cutoff_A'],'Å')
    print('独立周期像枚举的无向边数：',len(geom),'；原始报告：',orbit['n_undirected_edges'])
    print('独立重建边集合与原始报告一致：',set(geom)==saved_edges)
    print('源结构SHA-256与周期连通证据一致：',hashlib.sha256(POSCAR.read_bytes()).hexdigest()==json.loads(TOPOLOGY.read_text(encoding='utf-8'))['host_poscar_sha256'])
    print('\n路径类别 | 代表Li位点对+周期像 | 直线距离(Å) | 同类无向边Ni | 跨胞边数')
    total=0
    for o in orbit['orbits']:
        i,j,*R=o['representative']
        d=distance(i,j,R,A,frac)
        Ni=o['Ni'];total+=Ni
        n_cross=sum(any(e[2:]) for e in o['members'])
        assert abs(d-o['distance_min_A'])<1e-6
        print(f"  {o['path_id']:02}       | {i} -> {j}, R={tuple(R)} | {d:.6f} | {Ni} | {n_cross}")
    print('4类合计：',total,'条无向边；对应正反向有向边：',2*total)
    print('\ncutoff(Å) -> 独立枚举的无向候选边数：')
    for rc in [2.70,2.93,3.00,3.10,3.30,4.20,4.30]:
        count=sum(d<rc for d in all_edges.values())
        print(f'  {rc:.2f} -> {count}')
    print('\n请自己判断：3.3和4.2得到的几何候选集合一样吗？4.3新增了什么？')

    if '--network' in sys.argv:
        topology=json.loads(TOPOLOGY.read_text(encoding='utf-8'))
        print('\n== C1.2.3 第二小步：周期连通性（原始整数拓扑报告） ==')
        desired={(0,),(1,),(3,),(0,1),(0,2),(0,1,2,3)}
        for row in topology['rows']:
            subset=tuple(row['paths'])
            if subset in desired:
                ranks=[component['translation_rank'] for component in row['components']]
                print(f'路径={list(subset)}；连通分量数={row["component_count"]}；每分量周期延伸维数={ranks}；整张图单连通三维={row["single_connected_3d_periodic_graph"]}')
        print('注意：[0,2]存在2个互不相通的1D连通分量，不是2D；不能把不同分量的维数相加。')
        print('拓扑连通性≠跳跃速率≠实际D(T)：本报告没有使用NEB能垒。')

if __name__=='__main__':
    main()
