#!/usr/bin/env python3
"""C1 2.5：原始归档的只读教学核查，只用 Python 标准库。

默认 --step 1：核对已有能垒、复算速率、查看当前空位的候选事件。
--step 2：使用实际网络和人为指定随机数演示一次交换（不是历史轨迹）。
--step 3：由归档的分块 D 与时间复算 300 K 汇总和部分统计检查。
不提交任务、不运行正式 KMC、不修改任何输入或原始结果。
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parent
DATA = ROOT / '原始数据'
KB_EV = 8.617333262145e-5  # 与原始 transport.py 相同，eV/K
A2_TO_CM2 = 1e-16


def load(path: Path) -> dict[str, Any]:
    with path.open(encoding='utf-8-sig') as stream:
        obj = json.load(stream)
    if not isinstance(obj, dict):
        raise ValueError(f'需要 JSON 对象：{path.name}')
    return obj


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def close(a: float, b: float, message: str, rel: float = 1e-10) -> None:
    require(math.isfinite(a) and math.isfinite(b) and
            math.isclose(a, b, rel_tol=rel, abs_tol=1e-25), message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_sources() -> None:
    source = load(ROOT / '文件来源.json')
    for entry in source['files']:
        f = ROOT / entry['local_path']
        require(f.is_file(), f'缺少文件：{entry["local_path"]}')
        require(f.stat().st_size == entry['bytes'] and digest(f) == entry['sha256'],
                f'文件已变更：{entry["local_path"]}；请使用未编辑的原始副本。')
    model_hash = digest(DATA / 'transport_model.json')
    stats = load(DATA / 'KMC_statistical_admission.json')
    final = load(DATA / 'finite_model_final_admission.json')
    require(stats['status'] == 'FINITE_MODEL_KMC_STATISTICAL_GATES_PASSED',
            '最终统计报告状态不符。')
    require(final['status'] == 'FIRST_ROUND_FINITE_NEUTRAL_VACANCY_NEB_NETWORK_CONDITIONAL_TRACER_LOOP_PASSED',
            '有限模型最终证书状态不符。')
    require(stats['model_sha256'] == model_hash, '统计报告与网络模型 SHA 不一致。')
    for name in ['KMC_ext_00.json', 'KMC_ext_01.json', 'KMC_ext_02.json', 'KMC_ext_03.json']:
        row = load(DATA / name)
        require(row['source_model_sha256'] == model_hash, f'{name} 与模型不一致。')
        matches = [v for k, v in stats['source_raw_KMC_sha256'].items()
                   if k.replace('\\', '/').split('/')[-1] == name]
        require(matches == [digest(DATA / name)], f'{name} 未正确绑定最终统计报告。')
    for old_path, expected in final['admitted_evidence'].items():
        f = DATA / old_path.replace('\\', '/').split('/')[-1]
        if f.exists():
            require(digest(f) == expected, f'最终证书绑定不一致：{f.name}')
    print(f'[通过] {len(source["files"])}份原始文件 SHA-256；模型、KMC输出、最终证书关联一致。')


def rates_and_row() -> tuple[list[float], dict[str, Any]]:
    model = load(DATA / 'transport_model.json')
    row = load(DATA / 'KMC_ext_00.json')
    require(row['sampler'] == 'direct exchanges', '本教学核查仅适用于逐次交换采样。')
    require(row['energy_family'] == 'Efr', '能量口径不符。')
    barriers = model['families']['Efr']['barriers_sym_eV']
    require(barriers == row['barriers_eV'], '模型和KMC能垒不一致。')
    rates = []
    for i, evidence in enumerate(model['barrier_evidence']):
        require(evidence['path_id'] == i, 'path次序不一致。')
        energies = evidence['energies_eV']['Efr']
        barrier = energies['S_at_CI'] - (energies['A'] + energies['B'])/2
        close(barrier, barriers[i], '对称化能垒复算不一致。')
        rates.append(row['attempt_frequency_s'] * math.exp(-barrier/(KB_EV*row['T_K'])))
    return rates, row


def family_of(rate: float, family_rates: list[float]) -> int:
    # 原输出未直接存逐边 path_id。以唯一的4类速率匹配，并核对结果唯一。
    matches = [i for i, candidate in enumerate(family_rates)
               if math.isclose(rate, candidate, rel_tol=1e-12, abs_tol=0)]
    require(len(matches) == 1, '逐边速率不能唯一对应一类path。')
    return matches[0]


def validate_arrays(row: dict[str, Any], rates: list[float]) -> None:
    neighbors = row['supplied_neighbor']
    rr = row['supplied_rates_s']
    disp = row['supplied_displacement_A']
    require(len(neighbors) == len(rr) == len(disp) == row['nsites'], '网络维度不一致。')
    for i, neighbor_list in enumerate(neighbors):
        require(len(neighbor_list) == len(rr[i]) == len(disp[i]), '事件数不一致。')
        for q, target in enumerate(neighbor_list):
            require(0 <= target < row['nsites'] and target != i, '无效迁移位点。')
            family_of(rr[i][q], rates)
            reverse = [r for r, old in enumerate(neighbors[target])
                       if old == i and math.isclose(rr[target][r], rr[i][q], rel_tol=1e-12)
                       and all(math.isclose(disp[target][r][k], -disp[i][q][k],
                                            abs_tol=1e-12) for k in range(3))]
            require(bool(reverse), f'缺少反向边：{i}->{target}')


def step1() -> None:
    rates, row = rates_and_row()
    validate_arrays(row, rates)
    print('\n=== 第一步：输入 -> 跳跃速率 -> 候选事件 ===')
    print(f'T = {row["T_K"]:g} K；统一尝试频率 = {row["attempt_frequency_s"]:.3e} s^-1（假设）')
    print(f'Li位点 = {row["nsites"]}；Li原子 = {row["Li_count"]}；空位 = 1；空位比例 = 1/32')
    print('能量口径：Efr；四类能垒 = E_TS - (E_A+E_B)/2；沿用已核验的等价端态模型。')
    print('\npath      Em(eV)             Gamma(s^-1)')
    for i, (barrier, rate) in enumerate(zip(row['barriers_eV'], rates)):
        print(f'{i:4d}      {barrier:.9f}       {rate:.9e}')
    print('\n以空位在0号位点为例（这是状态示例，不是归档的第一条采样记录）：')
    print('候选Li所在位点      path     本次被选中的概率')
    total = math.fsum(row['supplied_rates_s'][0])
    for target, rate in zip(row['supplied_neighbor'][0], row['supplied_rates_s'][0]):
        print(f'{target:8d}              {family_of(rate, rates):d}         {rate/total:.6%}')
    print('注意：邻居编号是当前位点号，不是始终固定不变的Li身份编号。')
    print(f'候选事件合计 = {len(row["supplied_neighbor"][0])}；类别合计 = {len(rates)}')
    print(f'总速率 R = {total:.9e} s^-1；平均等待时间 1/R = {1/total:.9e} s')
    print('[通过] 速率与原始KMC输入一致；每条边均存在等速率、反位移的反向边。')
    print('本步未运行随机模拟，也未检查统计收敛。')


def step2() -> None:
    rates, row = rates_and_row()
    validate_arrays(row, rates)
    print('\n=== 第二步：一次交换的教学演示（不是历史轨迹回放） ===')
    print('使用实际网络；人为指定选事件随机数 u_event=0.5、等待随机数 u_time=0.5。')
    print('仅为说明算法：初始空位在0号位点，其余Li身份号设为初始位点号。')
    total = math.fsum(row['supplied_rates_s'][0])
    cumulative = 0.0
    choice = None
    for j, rate in enumerate(row['supplied_rates_s'][0]):
        cumulative += rate / total
        if 0.5 < cumulative:
            choice = j
            break
    require(choice is not None, '选事件失败。')
    target = row['supplied_neighbor'][0][choice]
    displacement = row['supplied_displacement_A'][0][choice]
    waiting = -math.log(0.5) / total
    print(f'选中 path {family_of(row["supplied_rates_s"][0][choice], rates)}：位点{target}的Li填入位点0。')
    print(f'空位：0 -> {target}；等待时间 = {waiting:.9e} s。')
    print('空位位移(Å)：', ', '.join(f'{v:.9f}' for v in displacement))
    print('迁移Li位移(Å)：', ', '.join(f'{-v:.9f}' for v in displacement))
    print('其他Li的本步位移为0；空位和此Li位移方向相反。')
    print('下一步应读取新空位所在行的候选事件；位点号与Li身份号不能再混用。')
    print('位移向量已包含跨周期连接；不能用折回后的位点坐标直接相减替代。')
    print('这一跳不足以估计长期扩散系数。')


def pooled(values: list[float], times: list[float]) -> tuple[float, float, float]:
    require(len(values) == len(times) and len(times) >= 2, '统计块不足或维度不一致。')
    require(all(math.isfinite(v) for v in values) and
            all(math.isfinite(t) and t > 0 for t in times), '非法分块数据。')
    total_t = math.fsum(times)
    numerators = [v*t for v,t in zip(values,times)]
    total_n = math.fsum(numerators)
    mean = total_n / total_t
    # 与原始 pooled_ratio 相同的删一块 jackknife；本阶段无需背推导。
    leave = [(total_n-n)/(total_t-t) for n,t in zip(numerators,times)]
    center = math.fsum(leave)/len(leave)
    se = math.sqrt((len(leave)-1)/len(leave)*math.fsum((x-center)**2 for x in leave))
    return mean*A2_TO_CM2, se*A2_TO_CM2, total_t


def group(names: list[str]) -> tuple[float, float, float]:
    values, times = [], []
    for name in names:
        r = load(DATA/name)['raw_KMC']
        values.extend(r['block_Dtracer_A2_s'])
        times.extend(r['block_time_s'])
    return pooled(values, times)


def step3() -> None:
    family_rates, reference = rates_and_row()
    stats = load(DATA/'KMC_statistical_admission.json')
    supplement = load(DATA/'correlation_error_supplement.json')
    final = load(DATA/'finite_model_final_admission.json')
    print('\n=== 第三步：复算32位点盒的300 K分块汇总 ===')
    print('输入保留的是每个块的D和时间；本步不从未保存的逐跳坐标重算MSD。')
    print('原算法每块：D_b=MSD_b/(6*t_b)；本步汇总：sum(D_b*t_b)/sum(t_b)。')
    print('Å^2/s -> cm^2/s：乘1e-16；不再乘空位浓度。')
    rows = []
    for name in ['KMC_ext_00.json','KMC_ext_02.json']:
        r=load(DATA/name); raw=r['raw_KMC']
        value,se,elapsed=pooled(raw['block_Dtracer_A2_s'],raw['block_time_s'])
        close(value,r['DLi_scalar_cm2_s'],'单种子D复算不一致。')
        close(se,r['DLi_standard_error_cm2_s'],'单种子SE复算不一致。',rel=1e-8)
        require(raw['seed'] in [137,991] and r['nsites']==32 and r['Li_count']==31 and
                r['T_K']==300 and raw['nblocks']==160, '不是指定主结果。')
        counts=raw['event_counts_by_directed_edge']
        flat=[x for a in counts for x in a]
        require(sum(flat)==raw['nevents']==8000000,'事件计数不一致。')
        inferred=[0,0,0,0]
        for a,rate_row in zip(counts,r['supplied_rates_s']):
            for count,rate in zip(a,rate_row):
                inferred[family_of(rate,family_rates)]+=count
        require(inferred==r['actual_observed_counts_by_family'],'路径类别计数不一致。')
        missing=sum(c==0 and rate>0 for a,b in zip(counts,r['supplied_rates_s']) for c,rate in zip(a,b))
        print(f'\n{name}：seed={raw["seed"]}，800万次计入统计的交换，160块，每块5万次。')
        print(f'采样前另有{raw["burn_events"]}次预运行，不计入本统计时间和事件数。')
        print(f'D_Li={value:.12e} cm^2/s；1SE={se:.9e}；模型累计时间={elapsed:.12e} s')
        print(f'各path交换次数={inferred}；未采到的正速率有向边={missing}/{len(flat)}')
        rows.append((value,se))
    mean,se,elapsed=group(['KMC_ext_00.json','KMC_ext_02.json'])
    primary=stats['primary_300K']
    close(mean,primary['DLi_cm2_s'],'汇总D不一致。')
    close(se,primary['DLi_SE_cm2_s'],'汇总SE不一致。',rel=1e-8)
    close(elapsed,primary['total_elapsed_time_s'],'汇总时间不一致。')
    close(mean,final['actual_results']['conditional_LiD300_cm2_s'],'最终证书D不一致。')
    seed_z=abs(rows[0][0]-rows[1][0])/math.hypot(rows[0][1],rows[1][1])
    close(seed_z,primary['seed_agreement_z'],'种子对照不一致。',rel=1e-8)
    short,short_se,_=group(['KMC_ext_01.json','KMC_ext_03.json'])
    paired_bound=abs(mean-short)/(se+short_se)
    expected=next(v for v in supplement['paired_block_checks'] if v['nsites']==32)
    close(paired_bound,expected['difference_over_SE_upper_bound'],'分块对照不一致。',rel=1e-8)
    print('\n只汇总两个不同种子的160块数据，不把同轨迹的320块版本重复当新样本。')
    print(f'D_Li = {mean:.12e} cm^2/s；KMC 1SE = {se:.9e}；相对1SE={se/mean:.4%}')
    print(f'两种子模型时间之和 = {elapsed:.12e} s（不是电脑耗时，也不是一条更长轨迹）。')
    print(f'种子差异/组合SE = {seed_z:.6f}（原门槛<=3）。')
    print(f'分块对照差异/SE保守上界 = {paired_bound:.6f}（原门槛<=3，不假定两种分块独立）。')
    require(seed_z<=3 and se/mean<=0.05 and paired_bound<=3, '当前子集未通过原门槛。')
    print('[通过] 32位点、300 K的汇总及上述统计复算与原最终记录一致。')
    print('范围：没有重新审计128位点、其他温度、NEB或模型物理真实性；统计SE不含这些误差。')
    print('旧输入与单次输出的PENDING字段是历史阶段标记；最终结论以最终证书及统计补充为准。')


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step', choices=['1','2','3','all'], default='1')
    args=parser.parse_args()
    try:
        verify_sources()
        for number, func in [('1',step1),('2',step2),('3',step3)]:
            if args.step in [number,'all']:
                func()
    except (OSError, ValueError, KeyError, TypeError, IndexError, OverflowError) as error:
        print(f'\n核查停止：{error}',file=sys.stderr)
        return 1
    return 0

if __name__=='__main__':
    raise SystemExit(main())
