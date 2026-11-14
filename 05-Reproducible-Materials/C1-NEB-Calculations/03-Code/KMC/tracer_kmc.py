"""Single-vacancy exchange KMC with particle identities and unwrapped motion.

This is a kinetic model, not AIMD. Rates come from explicitly supplied barriers.
One vacancy per finite periodic box; no implicit equilibrium vacancy concentration.
"""
import numpy as np

def simulate(neighbor,displacement,rates,nevents,nblocks,seed):
    neighbor=np.asarray(neighbor,dtype=int)
    displacement=np.asarray(displacement,dtype=float)
    rates=np.asarray(rates,dtype=float)
    nsites,degree=neighbor.shape
    if displacement.shape!=(nsites,degree,3) or rates.shape!=neighbor.shape:
        raise ValueError('inconsistent kinetic arrays')
    if nevents<nblocks or nblocks<1 or nsites<2 or not np.all(np.isfinite(rates)) or not np.all(np.isfinite(displacement)) or np.any(neighbor<0) or np.any(neighbor>=nsites):
        raise ValueError('invalid simulation controls or kinetic arrays')
    if np.any(rates<0) or np.any(rates.sum(axis=1)<=0):
        raise ValueError('positive outgoing rates required')
    total=rates.sum(axis=1)
    cdf=np.cumsum(rates,axis=1)/total[:,None]
    cdf[:,-1]=1.
    rng=np.random.default_rng(seed)
    occupancy=np.arange(nsites)
    occupancy[0]=-1
    vacancy=0
    li_disp=np.zeros((nsites,3))
    vac_disp=np.zeros(3)
    event_counts=np.zeros_like(neighbor)

    def advance(nsteps,track):
        nonlocal vacancy
        random=rng.random((nsteps,2))
        time=0.
        for u,v in random:
            old=vacancy
            choice=int(np.searchsorted(cdf[old],u,side='right'))
            target=neighbor[old,choice]
            ion=occupancy[target]
            if ion<0:
                raise ValueError('periodic box aliases jump onto vacancy itself')
            if track:
                event_counts[old,choice]+=1
                vec=displacement[old,choice]
                li_disp[ion]-=vec
                vac_disp[:]+=vec
                time+=-np.log(max(v,np.finfo(float).tiny))/total[old]
            occupancy[old]=ion
            occupancy[target]=-1
            vacancy=int(target)
        return time

    burn=min(20000,nevents//10)
    advance(burn,False)
    tensors=[];scalar=[];vscalar=[];times=[]
    for block in range(nblocks):
        li_disp.fill(0.);vac_disp.fill(0.)
        dt=advance(nevents//nblocks,True)
        D=li_disp.T@li_disp/(2*(nsites-1)*dt)
        tensors.append(D)
        scalar.append(float(np.trace(D)/3))
        vscalar.append(float(vac_disp@vac_disp/(6*dt)))
        times.append(dt)
    # Pool MSD and elapsed time. Equal averaging of MSD/random time has a
    # stopping-time bias, particularly for short fixed-event blocks.
    times=np.asarray(times);numerators=np.asarray(scalar)*times
    mean=numerators.sum()/times.sum()
    if nblocks>1:
        leave=(numerators.sum()-numerators)/(times.sum()-times)
        error=np.sqrt((nblocks-1)/nblocks*np.sum((leave-leave.mean())**2))
    else: error=None
    flags=['FINITE_TIME_ESTIMATE_REQUIRES_LAG_AND_SEED_CONVERGENCE']
    if np.any((rates>0)&(event_counts==0)): flags.append('UNSAMPLED_POSITIVE_RATE_EDGE')
    if nevents//nblocks<100: flags.append('SHORT_EVENT_BLOCKS')
    tensor=np.einsum('b,bij->ij',times,np.asarray(tensors))/times.sum()
    return dict(model='single vacancy, exclusion exchanges, supplied rates; not AIMD',estimator='pooled block MSD / total elapsed time; ratio jackknife uncertainty',long_time_convergence_claimed=False,diagnostic_flags=flags,event_counts_by_directed_edge=event_counts.tolist(),seed=seed,nsites=nsites,vacancy_fraction=1/nsites,nevents=(nevents//nblocks)*nblocks,nblocks=nblocks,burn_events=burn,Dtracer_tensor_A2_s=tensor.tolist(),Dtracer_mean_A2_s=float(mean),block_standard_error_A2_s=float(error) if error is not None else None,block_Dtracer_A2_s=scalar,block_Dvacancy_A2_s=vscalar,block_time_s=times.tolist())
