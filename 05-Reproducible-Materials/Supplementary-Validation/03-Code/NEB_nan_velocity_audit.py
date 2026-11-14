from pathlib import Path
import hashlib,json,re,sys,platform
import numpy as np
from ase.io import read,write
from io import StringIO
ROOT=Path(__file__).resolve().parents[4]
out=Path(__file__).resolve().parent
records=[]
for directory in ['cases','hpc/staging','hpc/results']:
 for p in (ROOT/directory).rglob('POSCAR'):
  text=p.read_text(errors='replace')
  if not re.search(r'\bnan\b',text,re.I):continue
  lines=text.splitlines(); n=sum(map(int,lines[6].split())); start=8 if not lines[7].lower().startswith('s') else 9
  positions=np.array([list(map(float,x.split()[:3])) for x in lines[start:start+n]])
  records.append(dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),nat=n,positions_finite=bool(np.isfinite(positions).all()),nan_in_position_block=bool(re.search(r'\bnan\b','\n'.join(lines[start:start+n]),re.I)),nan_in_tail=bool(re.search(r'\bnan\b','\n'.join(lines[start+n:]),re.I))))
assert records and all(x['positions_finite'] and not x['nan_in_position_block'] and x['nan_in_tail'] for x in records)
source=ROOT/'hpc/results/20261014-185852-db265980660a/20261014-203439-c82a525b9b0e/01/CONTCAR'
a=read(source);s=StringIO();write(s,a,format='vasp',direct=True,sort=False)
assert np.isfinite(a.positions).all() and not np.isfinite(a.get_velocities()).all()
assert 'Cartesian\n' in s.getvalue() and 'nan' in s.getvalue()
report=dict(date='2026-11-01',verdict='NAN_IN_OPTIONAL_VELOCITY_TAIL_NOT_ION_POSITIONS',files=records,source=str(source.relative_to(ROOT)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),propagation_reproduced=True,mechanism='VASP CONTCAR already contains NaN velocities; ASE read stores momenta and write emits Cartesian velocity block; restore_band preserves this array',origin_limit='POTIM=0 is present in VTST inputs; precise internal VASP division causing NaN is not proven without source-level tracing',no_original_files_changed=True,python=platform.python_version())
(out/'audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(files_with_nan=len(records),all_position_blocks_finite=True,propagation_reproduced=True)))
