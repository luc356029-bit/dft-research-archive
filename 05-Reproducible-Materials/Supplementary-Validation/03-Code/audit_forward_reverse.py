"""Audit endpoint equivalence and both directional barriers without new DFT."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
import contextlib
import hashlib
import io
import json
import sys
import xml.etree.ElementTree as ET
from unittest.mock import patch
import numpy as np
from ase.io import read
from ase.geometry import find_mic
from ase.dft.kpoints import monkhorst_pack

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from analyze_neb import analyze_band
from analyze_vasp import analyze
from network_arrays import primitive_edges, supercell_arrays

J = lambda p: json.loads(Path(p).read_text(encoding='utf-8-sig'))
SHA = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
dest = ROOT/'analysis/transport_phase34_20261015'
model_path, cert_path = dest/'transport_model.json', dest/'network_admission_certificate.json'
model, cert = J(model_path), J(cert_path)
assert SHA(model_path) == cert['model_sha256']
prep_path = ROOT/'analysis/vacancy63_endpoints_k4/endpoint_preparation.json'
prep = J(prep_path)
host_path = ROOT/'analysis/host_supercell64/POSCAR'
host = read(host_path)
cell, hf = np.asarray(host.cell), host.get_scaled_positions()
G = cell@cell.T
assert len(host) == 64 and sum(host.numbers==3) == 32

def no_write_analyze(case):
    with contextlib.redirect_stdout(io.StringIO()), patch.object(Path,'write_text',return_value=0):
        result = analyze(case)
    assert result['SCF_converged'] and result['association']['input_manifest_verified']
    assert result['association']['cross_format_association_passed']
    return result

def kp_set(points):
    return {tuple(p) for p in np.round(np.mod(points,1),12)}

def kgrid(case):
    x = ET.parse(Path(case)/'vasprun.xml').getroot()
    generation = x.find('kpoints/generation')
    mesh = [int(v) for v in generation.find("v[@name='divisions']").text.split()]
    shift = [float(v) for v in generation.find("v[@name='usershift']").text.split()]
    actual = np.array([[float(v) for v in node.text.split()] for node in x.findall("kpoints/varray[@name='kpointlist']/v")])
    params = {n.attrib['name']:n.text.strip() for n in x.findall('.//parameters//i') if n.text}
    assert generation.attrib['param']=='Monkhorst-Pack' and mesh==[4,4,4] and shift==[0.,0.,0.]
    assert params['ISYM']=='0' and params['ISPIN']=='1' and params['LSORBIT']=='F' and params['LNONCOLLINEAR']=='F'
    full = monkhorst_pack(mesh)
    assert kp_set(np.vstack([actual,-actual])) == kp_set(full)
    return full

rows, permutations, raw = [], [], {}
for q, source, stored in zip(prep['paths'],cert['actual_sources'],model['barrier_evidence']):
    pid = q['path_id']
    assert pid == source['path_id'] == stored['path_id']
    report_path = Path(source['report'])
    assert SHA(report_path) == source['report_sha256'] == stored['source_sha256']
    for n, h in source['raw_sha256'].items():
        p = Path(source['output_directory'])/n
        assert SHA(p)==h
        raw[str(p)]=h
    with contextlib.redirect_stdout(io.StringIO()), patch.object(Path,'write_text',return_value=0):
        band = analyze_band(report_path.parent)
    assert band['status']=='CI_CANDIDATE_FORCE_GATE_PASSED' and band['final_CI_raw_force_gate_passed']
    A, B = [no_write_analyze(Path(e['path'])) for e in source['actual_endpoint_sources']]
    for actual, expected in zip([A,B],source['actual_endpoint_sources']):
        assert actual['raw_sha256']==expected['sha256'] and actual['stationary_at_0p01']
        for name,h in actual['raw_sha256'].items():raw[str(Path(actual['case'])/name)]=h
    assert Path(prep['source']).resolve()==Path(A['case']).resolve()
    # Inputs carry the identical full-precision fixed cell; XML prints only 8 decimals.
    # Half a printed last unit is 5e-9 A, not a physical cell mismatch.
    for actual in [A,B]:
        assert np.max(abs(read(Path(actual['case'])/'POSCAR').cell-cell))<1e-12
        assert np.max(abs(np.array(actual['last_frame']['cell_A'])-cell))<=5.01e-9
    R=np.array(q['rotation_fractional']);t=np.array(q['translation_fractional'])
    perm=np.array(q['host_permutation']);src=np.array(q['source_atom_permutation'])
    assert sorted(perm.tolist())==list(range(64)) and sorted(src.tolist())==list(range(63))
    assert np.array_equal(host.numbers,host.numbers[perm])
    assert perm[0]==q['final_vacancy_host_slot']
    metric_error=float(np.linalg.norm(R.T@G@R-G)/np.linalg.norm(G))
    host_error=float(max(find_mic((hf@R.T+t-hf[perm])@cell,cell,pbc=True)[1]))
    assert metric_error<1e-12 and host_error<1e-10
    Af=np.array(A['last_frame']['positions_fractional']);Bf=np.array(B['last_frame']['positions_fractional'])
    mapped=(Af[src]@R.T+t)@cell
    difference=Bf@cell-mapped
    shift=find_mic(difference,cell,pbc=True)[0][0]
    endpoint_error=float(max(find_mic(difference-shift,cell,pbc=True)[1]))
    assert endpoint_error<2e-6
    assert np.array_equal(read(Path(A['case'])/'POSCAR').numbers[src],read(Path(B['case'])/'POSCAR').numbers)
    for actual in [A,B]:
        full=kgrid(actual['case'])
        assert kp_set(full@np.linalg.inv(R))==kp_set(full)
    permutations.append(perm)
    h=band['history'][-1];saddle=h['climbing_image']
    energies={}
    for label,key,ekey in [('Efr','Efr_eV','e_fr_energy'),('E0','E0_eV','e_0_energy')]:
        E=h[key];a,b,s=map(float,[E[0],E[-1],E[saddle]])
        assert abs(a-A['last_frame']['energies_eV'][ekey])<1e-10
        assert abs(b-B['last_frame']['energies_eV'][ekey])<1e-10
        assert stored['energies_eV'][label]==dict(A=a,B=b,S_at_CI=s)
        forward,reverse=s-a,s-b
        symmetric=s-(a+b)/2
        assert abs((forward-reverse)-(b-a))<1e-13
        assert abs(symmetric-stored['barriers_sym_eV'][label])<1e-13
        energies[label]=dict(E_A_eV=a,E_B_eV=b,E_saddle_eV=s,delta_B_minus_A_eV=b-a,
            forward_barrier_eV=forward,reverse_barrier_eV=reverse,symmetrized_barrier_eV=symmetric)
    rows.append(dict(path_id=pid,vacancy_transition_host_slots=[0,q['final_vacancy_host_slot']],
        operation_R=R.tolist(),operation_t=t.tolist(),host_metric_relative_error=metric_error,
        host_mapping_max_error_A=host_error,actual_evaluated_defect_mapping_max_error_A=endpoint_error,
        full_actual_MP4_grid_invariant=True,energies=energies,
        Ewo_endpoint_split_eV=B['last_frame']['energies_eV']['e_wo_entrp']-A['last_frame']['energies_eV']['e_wo_entrp']))

# The actual boundary-compatible operations generate the complete 32-Li host-site orbit.
seen={0}
while True:
    new=seen|{int(p[i]) for p in permutations for i in seen}
    if new==seen:break
    seen=new
assert seen==set(np.flatnonzero(host.numbers==3))
primitive_path=ROOT/'analysis/host_symmetry_projection/POSCAR'
primitive=read(primitive_path)
rate_row=next(r for r in model['families']['Efr']['temperature_rows'] if r['T_K']==300.)
rates=rate_row['rates_s']
edges=primitive_edges(primitive,model['orbits'],rates)
edge_rates={(i,j,*map(int,S)):k for i,j,S,k in edges}
assert len(edges)==len(edge_rates)==48
for (i,j,*S),k in edge_rates.items():assert edge_rates[(j,i,*(-x for x in S))]==k
mapping_path=ROOT/'analysis/host_supercell64/supercell_mapping.json'
arrays=supercell_arrays(primitive,host,J(mapping_path),model['orbits'],rates)
assert arrays['neighbor'].shape==(32,6)
directed={}
for i in range(32):
    for n,j in enumerate(arrays['neighbor'][i]):
        directed[(i,int(j),int(arrays['path_id'][i,n]),*map(int,arrays['integer_image'][i,n]))]=(float(arrays['rates'][i,n]),arrays['displacement'][i,n])
assert len(directed)==192
for (i,j,pid,*S),(k,v) in directed.items():
    reverse=directed[(j,i,pid,*(-x for x in S))]
    assert k==reverse[0] and np.array_equal(-v,reverse[1])
delta=max(abs(r['energies']['Efr']['delta_B_minus_A_eV']) for r in rows)
kb=8.617333262145e-5
deps=[ROOT/'scripts'/n for n in ['endpoint_mapping.py','hop_geometry.py','analyze_vasp.py','analyze_neb.py',
    'compute_round1_transport.py','network_arrays.py']]
output=dict(timestamp=datetime.now(timezone(timedelta(hours=8))).isoformat(),
    verdict='EQUIVALENT_ENDPOINT_MODEL_AND_BOTH_DIRECTIONS_REVALIDATED',
    scope='Existing four-pair finite neutral-vacancy model; not proof of all mechanisms or general inequivalent-site kinetics.',
    no_new_DFT_or_KMC=True,old_reports_not_modified=True,paths=rows,
    XML_cell_precision_note='The initial audit used a1e-10 A comparison and stopped. All actual input cells are byte-equivalent in numerical values to the host; XML prints8decimal basis components, max difference4.758011e-9 A, within half-last-digit rounding. Audit now checks identical input cells and the explicit5e-9 A XML rounding bound. No DFT parameter or scientific tolerance changed.',
    symmetry_generated_Li_host_orbit_size=len(seen),primitive_directed_edges=48,supercell_directed_edges=192,
    exact_reverse_edges_and_equal_model_rates_present=True,
    maximum_observed_Efr_endpoint_split_eV=delta,
    maximum_each_direction_symmetrization_shift_eV=delta/2,
    maximum_relative_rate_shift_from_that_numerical_split_at_300K=float(np.expm1(delta/(2*kb*300))),
    model_symmetry_convention=model['symmetry_convention'],
    directional_count_interpretation='Four undirected representative channels yield eight directional barrier entries, from four bands. Reversal of the same path on one energy surface shares its saddle; inequivalent endpoints need unequal directional barriers, not an automatically independent reverse NEB.',
    conditional_warning='If complete defect states are truly inequivalent or reconstruct differently, equal-site symmetrization is invalid: retain state energies, direction-specific rates and detailed balance; re-enumerate any newly independent environments/channels.',
    primary_reference='https://vasp.at/wiki/index.php/Nudged_elastic_bands',
    raw_output_sha256=raw,source_evidence_sha256={str(p):SHA(p) for p in [model_path,cert_path,prep_path,host_path,primitive_path,mapping_path,*deps]},
    analyzer_path=str(Path(__file__).resolve()),analyzer_sha256=SHA(__file__))
out=HERE/'forward_reverse_endpoint_audit_20261027.json'
assert not out.exists()
out.write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:output[k] for k in ['verdict','symmetry_generated_Li_host_orbit_size','primitive_directed_edges',
    'supercell_directed_edges','maximum_observed_Efr_endpoint_split_eV','maximum_each_direction_symmetrization_shift_eV',
    'maximum_relative_rate_shift_from_that_numerical_split_at_300K']},ensure_ascii=False,indent=2))
print(json.dumps([{k:p[k] for k in ['path_id','actual_evaluated_defect_mapping_max_error_A','energies']} for p in rows],ensure_ascii=False,indent=2))
