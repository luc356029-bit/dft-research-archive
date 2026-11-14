"""Construct an equivalent vacancy endpoint by an exact host space operation.

Only symmetry-compatible supercells and reciprocal grids may use the associated
energy equivalence. A real DFT single point must still validate the endpoint.
"""
import numpy as np
from ase.geometry import find_mic
from scipy.optimize import linear_sum_assignment
from hop_geometry import operation_maps

def audit_assignment(host,defect,host_slots):
    report={};assigned=np.full(len(defect),-1,dtype=int)
    for species in sorted(set(defect.numbers)):
        ai=np.flatnonzero(defect.numbers==species);hi=np.flatnonzero(host.numbers==species)
        delta=defect.positions[ai,None,:]-host.positions[None,hi,:]
        _,dist=find_mic(delta.reshape(-1,3),host.cell,pbc=True)
        dist=dist.reshape(len(ai),len(hi))
        rows,cols=linear_sum_assignment(dist**2)
        assigned[ai[rows]]=hi[cols]
        missing=sorted(set(map(int,hi))-set(map(int,hi[cols])))
        report[str(species)]=dict(missing_host_slots=missing,max_assignment_distance_A=float(dist[rows,cols].max()),rms_assignment_distance_A=float(np.sqrt(np.mean(dist[rows,cols]**2))))
    report['assigned_host_slots']=assigned.tolist()
    report['identity_mapping_preserved']=assigned.tolist()==list(host_slots)
    return report

def equivalent_endpoint(host,initial,host_slots,target,vacancy_slot=0,operations=None):
    host_slots=np.asarray(host_slots,dtype=int)
    if len(host_slots)!=len(initial) or len(set(host_slots))!=len(initial) or target not in host_slots or vacancy_slot in host_slots:
        raise ValueError('invalid vacancy identity mapping')
    cell=np.asarray(host.cell)
    _,ops,error=operation_maps(host,1e-5) if operations is None else (None,operations,0.)
    delta=find_mic(initial.positions-host.positions[host_slots],cell,pbc=True)[0]@np.linalg.inv(cell)
    desired=np.where(host_slots==target,vacancy_slot,host_slots)
    atom_of_slot={int(h):i for i,h in enumerate(host_slots)}
    mover=int(np.flatnonzero(host_slots==target)[0]);candidates=[]
    for R,t,perm,images in ops:
        if perm[vacancy_slot]!=target: continue
        inverse=np.argsort(perm)
        src=np.array([atom_of_slot[int(inverse[h])] for h in desired])
        if not np.array_equal(initial.numbers[src],initial.numbers):
            raise ValueError('chemical species permutation mismatch')
        positions=host.positions[desired]+delta[src]@R.T@cell
        positions=initial.positions+find_mic(positions-initial.positions,cell,pbc=True)[0]
        si=initial.numbers==14
        shift=np.mean(positions[si]-initial.positions[si],axis=0)
        positions-=shift
        displacement=positions-initial.positions
        score=float(np.sum(np.delete(displacement,mover,axis=0)**2))
        candidates.append((score,positions,R,t,perm,src,shift))
    if not candidates: raise ValueError('no exact supercell operation maps requested vacancy sites')
    score,positions,R,t,perm,src,shift=min(candidates,key=lambda c:c[0])
    final=initial.copy();final.positions=positions
    report=dict(status='SYMMETRY_GENERATED_REQUIRES_DFT_CHECK',initial_vacancy_host_slot=vacancy_slot,final_vacancy_host_slot=target,
        mover_atom_index=mover,destination_host_slots=desired.tolist(),source_atom_permutation=src.tolist(),
        rotation_fractional=R.tolist(),translation_fractional=t.tolist(),host_permutation=perm.tolist(),
        rigid_Si_alignment_shift_A=shift.tolist(),atomic_displacements_A=(positions-initial.positions).tolist(),
        mover_distance_A=float(np.linalg.norm(positions[mover]-initial.positions[mover])),
        max_nonmover_displacement_A=float(np.linalg.norm(np.delete(positions-initial.positions,mover,axis=0),axis=1).max()),
        minimum_periodic_pair_distance_A=float(np.min(final.get_all_distances(mic=True)+np.eye(len(final))*1e6)),
        max_host_operation_error_A=error,n_equivalent_operation_candidates=len(candidates))
    return final,report
