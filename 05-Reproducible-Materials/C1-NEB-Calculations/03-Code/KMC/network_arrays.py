"""Bind crystallographic integer images to transport and exclusion KMC.

For corrector diffusion use integer cell translations instead of floating
internal site displacements. This is an exact bounded-position gauge change,
followed by Dcart=cell.T @ Dfractional @ cell. It removes roundoff loops that
could spuriously diffuse in an extremely fast but geometrically bounded cage.
KMC uses actual jump vectors, with reverse vectors made by exact negation.
"""
import numpy as np
from ase.geometry import find_mic

def primitive_edges(atoms,orbits,rates):
    li=np.flatnonzero(atoms.numbers==3);local={int(h):i for i,h in enumerate(li)}
    edges=[]
    for orbit in orbits:
        rate=float(rates[orbit['path_id']])
        for i,j,*s in orbit['members']:
            v=np.asarray(s,dtype=int)
            edges.extend([(local[i],local[j],v,rate),(local[j],local[i],-v,rate)])
    return edges

def supercell_arrays(primitive,supercell,mapping,orbits,rates):
    indices=np.asarray(mapping['basis_indices']);images=np.asarray(mapping['primitive_images'])
    P=np.asarray(mapping['supercell_matrix']);inverse=np.linalg.inv(P)
    li=np.flatnonzero(supercell.numbers==3);local={int(h):i for i,h in enumerate(li)}
    pf=primitive.get_scaled_positions();sf=supercell.get_scaled_positions()
    rows=[[] for _ in li]
    for orbit in orbits:
        pid=orbit['path_id'];rate=float(rates[pid])
        for i,j,*S in orbit['members']:
            S=np.array(S,dtype=int);vec=(pf[j]+S-pf[i])@primitive.cell
            for h in li[indices[li]==i]:
                destination=(pf[j]+S+images[h])@inverse
                candidates=li[indices[li]==j]
                _,distance=find_mic((destination-sf[candidates])@supercell.cell,supercell.cell,pbc=True)
                h2=int(candidates[np.argmin(distance)])
                if min(distance)>1e-6: raise ValueError('supercell jump mapping not exact')
                image=np.rint(destination-sf[h2]).astype(int)
                if np.max(abs(destination-sf[h2]-image))>1e-7: raise ValueError('noninteger supercell image')
                a,b=local[int(h)],local[h2]
                rows[a].append((b,vec,rate,pid,image))
                rows[b].append((a,-vec,rate,pid,-image))
    if len(set(map(len,rows)))!=1: raise ValueError('nonuniform coordination requires a ragged-array model')
    rows=[sorted(r,key=lambda e:(e[3],e[0],*e[4])) for r in rows]
    return dict(host_slots=li,neighbor=np.array([[e[0] for e in r] for r in rows],dtype=int),
        displacement=np.array([[e[1] for e in r] for r in rows],dtype=float),
        rates=np.array([[e[2] for e in r] for r in rows],dtype=float),path_id=np.array([[e[3] for e in r] for r in rows],dtype=int),
        integer_image=np.array([[e[4] for e in r] for r in rows],dtype=int))
