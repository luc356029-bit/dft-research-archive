# NEB POSCAR尾部NaN：原因、影响与证据边界

记录日期：2026-11-01。针对LiSi首轮NEB的01/POSCAR中Cartesian块全为nan的疑问。本次仅检查已有文件和复现读写，不运行新DFT、不覆盖历史证据。

## 结论

检查到的NaN位于可选离子速度块，不在原子位置中。代表文件 cases/lisi-r1-neb00-ci-k4-v1/01/POSCAR 的前63行Direct数据是正常的原子分数坐标；随后Cartesian标记对应速度块，不是第二套原子坐标。

VASP官方格式说明明确区分位置与可选速度部分：[POSCAR格式](https://vasp.at/wiki/POSCAR)。速度块也可用Cartesian标记，不能仅凭该词认定其后是位置。

## 已确认的传递原因

1. 前一轮普通NEB的源文件 hpc/results/20261014-185852-db265980660a/20261014-203439-c82a525b9b0e/01/CONTCAR 已带NaN速度，位置仍有限。
2. neb_restart.py 的 restore_band 用ASE读取该CONTCAR。读取后的Atoms不仅保存位置，还保留速度对应的momenta数组。
3. prepare_neb.py将该Atoms写成新POSCAR时，ASE也写出了保留的速度数组，并显式加上Cartesian标记。将位置展开为最小像不会自动删除速度。
4. 本次用同一源CONTCAR直接执行ASE读写，复现了“位置正常、尾部Cartesian下全nan”的现象。

因此，已确认的问题是续算输入整理未移除无用速度数组，原有几何门控只检查位置、晶胞和路径，没有覆盖该附加数据。这是本项目输入整理遗漏。

## 检查范围与影响

扫描 cases、hpc/staging、hpc/results 中的POSCAR，共找到112份含NaN的文件；逐一按原子数解析位置块，全部位置有限，NaN全部位于尾部。这112份包含工作副本、上传暂存和结果快照，不表示112次独立计算。

这些证据排除了此次发现属于“原子坐标全NaN”的情形。它本身不能否定已验收NEB的能量与力，也不能单凭尾部问题重新证明势垒正确；科学结果仍应依据原有SCF、实际评估坐标、能量、力和输入关联证据。没有因此重算或修改Ea、D300。

## 尚未确认的部分

为何VASP最初将速度写成NaN，尚未做其内部源码级追踪。VTST输入中POTIM=0是观察到的条件；“内部除以零造成NaN”只能作为待验证猜测，不能写成确定原因。静态NEB优化中的该速度尾部也不能解释为AIMD物理速度。

## 后续输入处理建议

新建静态NEB续算输入应在独立Atoms副本中移除继承的momenta（例如 atoms.arrays.pop('momenta', None)），再写POSCAR，并检查位置、晶胞的所有数值有限以及没有意外速度块。对已经绑定SHA的历史输入输出不直接覆盖，应保留原件并另存清理版本及其来源记录。AIMD重启需要真实速度，不能套用删除规则。

此次只完成诊断、记录和发布，尚未修改生成器；不把建议表述成已实施的修复。

## 附件

- [审计脚本](https://drive.google.com/file/d/1RUxkP8HCW0vJ4BXlIddn9gn0rY2cDJ76/view)：此次扫描和ASE传递复现脚本（本地名audit.py）。
- [审计JSON](https://drive.google.com/file/d/1exzwAGYXrOKqfJgt3_orMe77f8kgqwaH/view)：112份文件的路径、SHA256、位置/尾部分项结论及源CONTCAR的SHA256（本地名audit.json）。
- 复核环境：Python3.14、ASE3.29.0、NumPy2.4.6；诊断复用已有结构文件。
