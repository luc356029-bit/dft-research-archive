"""Read-only geometry audit of the historic 3.3 A Li-site candidate radius."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
import hashlib
import itertools
import json
import sys
import ase
import numpy as np
import spglib
from ase.io import read
from ase.neighborlist import neighbor_list

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from hop_geometry import classify_hops
from network_arrays import primitive_edges
from periodic_topology import topology

SHA = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
host = ROOT / 'analysis/host_symmetry_projection/POSCAR'
old_report = ROOT / 'analysis/complete_orbit_connectivity_INTEGER.json'
projection_report = ROOT / 'analysis/host_symmetry_projection/projection_report.json'
old = json.loads(old_report.read_text(encoding='utf-8-sig'))
assert SHA(host) == old['host_poscar_sha256']
assert old['additional_cutoff_check_A'] == 4.2 and old['additional_cutoff_added_orbits'] == 0
atoms = read(host)
f, cell = atoms.get_scaled_positions(), np.asarray(atoms.cell)
li = np.flatnonzero(atoms.numbers == 3).tolist()
assert len(atoms) == 16 and len(li) == 8
maximum = 5.0
# For v=(delta_f+S)A, |S_k| <= R||column_k(A^-1)||+|delta_f,k|.
# Fractional coordinates are wrapped, so |delta_f,k|<1; this finite box covers all images within R.
bounds = np.ceil(maximum*np.linalg.norm(np.linalg.inv(cell), axis=0)+1).astype(int)
independent = {}
for i, j in itertools.product(li, repeat=2):
    for s in itertools.product(*(range(-int(b), int(b)+1) for b in bounds)):
        if i == j and s == (0, 0, 0):
            continue
        distance = float(np.linalg.norm((f[j]+np.array(s)-f[i])@cell))
        if distance < maximum:
            edge = min((i,j,*s), (j,i,*(-x for x in s)))
            if edge in independent:
                assert abs(independent[edge]-distance) < 1e-10
            independent[edge] = distance
i, j, images, distances = neighbor_list('ijSd', atoms, maximum)
ase_edges = {min((int(a),int(b),*map(int,s)), (int(b),int(a),*(-int(x) for x in s)))
             for a,b,s in zip(i,j,images) if atoms.numbers[a] == atoms.numbers[b] == 3}
assert ase_edges == set(independent)
shells = []
for edge, distance in sorted(independent.items(), key=lambda p:p[1]):
    if not shells or distance-shells[-1]['distance_min_A'] > 1e-6:
        shells.append(dict(distance_min_A=distance, distance_max_A=distance, undirected_edges=1))
    else:
        shells[-1]['distance_max_A'] = distance
        shells[-1]['undirected_edges'] += 1
assert len(shells) == 10
baseline = classify_hops(atoms, 3.3, 1e-5)
baseline_edges = {tuple(e) for o in baseline['orbits'] for e in o['members']}
assert len(baseline['orbits']) == 4 and len(baseline_edges) == 24
assert [o['Ni'] for o in baseline['orbits']] == [4,8,4,8]
rows = []
for radius in [2.70,2.93,3.00,3.10,3.20,3.30,4.20,4.30,5.00]:
    result = classify_hops(atoms, radius, 1e-5)
    edges = {tuple(e) for o in result['orbits'] for e in o['members']}
    assert edges == {e for e,d in independent.items() if d < radius}
    graph = topology(len(li), primitive_edges(atoms,result['orbits'],[1.]*len(result['orbits'])))
    rows.append(dict(radius_A=radius, orbit_count=len(result['orbits']), undirected_edges=len(edges),
        same_edge_set_as_3p3_A=edges==baseline_edges, geometric_3d_connected=graph['single_connected_3d_periodic_graph'],
        component_count=graph['component_count'], component_translation_ranks=[c['translation_rank'] for c in graph['components']]))
assert all(row['same_edge_set_as_3p3_A'] for row in rows if row['radius_A'] in [3.1,3.2,3.3,4.2])
assert next(row for row in rows if row['radius_A']==4.3)['orbit_count'] == 5
assert next(row for row in rows if row['radius_A']==2.93)['geometric_3d_connected']
paper = next((ROOT/'sources/repository').rglob('论文原文提取.txt'))
b_record = ROOT/'sources/repository/高通量筛选、输运研究复刻/工作流/B分析.md'
dependencies = [ROOT/'scripts'/name for name in ['hop_geometry.py','network_arrays.py','periodic_topology.py',
                    'precheck_orbit_connectivity.py','project_host.py']]
evidence = dict(timestamp=datetime.now(timezone(timedelta(hours=8))).isoformat(),
    scope='Geometry-only explanatory supplement. No new DFT, barriers, KMC or transport result.',
    source_host=str(host), source_host_sha256=SHA(host), host_atoms=16, Li_sites=8,
    geometric_cutoff_definition='Include every periodic Li-Li pair with ||(f_j+S-f_i)A|| < r_c; retain full periodic images, count reverse pairs once.',
    independent_enumeration=dict(maximum_radius_A=maximum, integer_image_bounds_per_axis=bounds.tolist(),
        complete_by_reciprocal_lattice_bound=True, agrees_with_ASE_neighbor_list=True,
        agrees_with_symmetry_enumeration_at_all_tested_radii=True),
    first_ten_distance_shells=shells,
    four_baseline_orbits=[{k:o[k] for k in ['path_id','Ni','distance_min_A','distance_max_A']} for o in baseline['orbits']],
    cutoff_sweep=rows,
    last_included_shell_max_A=shells[3]['distance_max_A'],
    first_excluded_shell_min_A=shells[4]['distance_min_A'],
    gap_between_fourth_and_fifth_shell_A=shells[4]['distance_min_A']-shells[3]['distance_max_A'],
    lower_safety_margin_for_3p3_A=3.3-shells[3]['distance_max_A'],
    upper_safety_margin_for_3p3_A=shells[4]['distance_min_A']-3.3,
    first_geometrically_3d_connected_shell_distance_A=shells[1]['distance_max_A'],
    historic_choice_status='3.3 A is an explicitly adopted finite candidate radius. Historic scripts hard-code this value; no unique numerical derivation or paper-specified LiSi 3.3 A rule was recorded. This later audit supplies a geometric justification, not a reconstructed contemporaneous derivation.',
    old_4p2_control_interpretation='Both 3.3 and 4.2 A stop before the fifth shell at about4.2251 A. The existing no-added-orbits check is correct, but does not test the first new shell or prove energetic mechanism completeness.',
    limits=['Short Li-site endpoint distance is not a proof of a low barrier or of an elementary hop.',
            'Graph connectivity is not energetic pathway completeness.',
            'The new4.3/5.0 A candidates have no new NEB validation and are not added to the admitted Ea/D model.'],
    source_evidence_sha256={str(p):SHA(p) for p in [host,old_report,projection_report,paper,b_record,*dependencies]},
    analyzer_path=str(Path(__file__).resolve()),analyzer_sha256=SHA(__file__),
    versions=dict(python=sys.version,ASE=ase.__version__,numpy=np.__version__,spglib=spglib.__version__),
    independent_edges_to_5A=[dict(edge=list(e),distance_A=d) for e,d in sorted(independent.items())])
out = HERE/'LiSi_3p3A_cutoff_audit_20261027.json'
assert not out.exists(), 'Never overwrite an audit'
out.write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:evidence[k] for k in ['source_host_sha256','last_included_shell_max_A','first_excluded_shell_min_A',
    'gap_between_fourth_and_fifth_shell_A','lower_safety_margin_for_3p3_A','first_geometrically_3d_connected_shell_distance_A',
    'cutoff_sweep']},ensure_ascii=False,indent=2))
