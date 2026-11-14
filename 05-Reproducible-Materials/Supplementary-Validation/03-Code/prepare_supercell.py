from pathlib import Path
import json,itertools,sys
import numpy as np
from ase.io import read,write
from ase.build import make_supercell
from hop_geometry import operation_maps
ROOT=Path(__file__).resolve().parents[1]
source=Path(sys.argv[1]).resolve()
dest=ROOT/'analysis/host_supercell64'
dest.mkdir(exist_ok=False)
atoms=read(source,format='vasp')
P=np.array([[0,1,-1],[-1,1,1],[2,0,0]])
_,ops,_=operation_maps(atoms,1e-5)
for R,t,p,m in ops:
    superR=np.linalg.inv(P).T@R@P.T
    assert np.max(abs(superR-np.rint(superR)))<1e-8,'Supercell breaks required point operation'
atoms.new_array('basis_index',np.arange(len(atoms)))
supercell=make_supercell(atoms,P,order='atom-major')
write(dest/'POSCAR',supercell,format='vasp',direct=True,sort=False,vasp5=True)
period=min(np.linalg.norm(np.array(n)@supercell.cell) for n in itertools.product(range(-3,4),repeat=3) if n!=(0,0,0))
images=supercell.positions@np.linalg.inv(atoms.cell)-atoms.get_scaled_positions()[supercell.arrays['basis_index']]
assert np.max(abs(images-np.rint(images)))<1e-8
report=dict(source=str(source),supercell_matrix=P.tolist(),nat=len(supercell),shortest_periodic_translation_A=float(period),cell_A=supercell.cell.tolist(),basis_indices=supercell.arrays['basis_index'].tolist(),primitive_images=np.rint(images).astype(int).tolist(),Li_vacancy_fraction=1/32,status='PRISTINE_GEOMETRY_REQUIRES_SUPERCELL_DFT_GRID_CHECK')
(dest/'supercell_mapping.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
assert len(supercell)==64 and sum(supercell.numbers==3)==32
assert period>9.
print(json.dumps({k:v for k,v in report.items() if k not in ['basis_indices','primitive_images']},indent=2))
