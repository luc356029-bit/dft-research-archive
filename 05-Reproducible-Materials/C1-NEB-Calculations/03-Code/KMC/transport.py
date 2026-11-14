"""Reversible periodic continuous-time jump process, including connectivity.

Lengths in A, rates in s^-1, output in A^2/s. This is vacancy diffusion;
it does not silently supply a Li tracer correlation factor or concentration.
Corrector formulation includes the cancellation of bounded internal motion.
"""
import numpy as np
import mpmath as mp

KB_EV = 8.617333262145e-5

def rates_from_barriers(barriers_eV, temperature_K, frequency_s):
    if not np.isfinite(temperature_K) or not np.isfinite(frequency_s) or temperature_K <= 0 or frequency_s <= 0:
        raise ValueError('positive T and attempt frequency required')
    barriers=np.asarray(barriers_eV,dtype=float)
    if not np.all(np.isfinite(barriers)): raise ValueError('finite barriers required')
    with np.errstate(over='ignore',under='ignore'):
        rates=frequency_s*np.exp(-barriers/(KB_EV*temperature_K))
    if not np.all(np.isfinite(rates)) or np.any(rates<=0):
        raise ValueError('rate overflow/underflow; do not silently remove physical edges')
    return rates

def periodic_diffusion(nsites, directed_edges, equilibrium=None):
    if not isinstance(nsites,int) or nsites<1:
        raise ValueError('positive integer site count required')
    edges=[]; adjacency=[set() for _ in range(nsites)]
    for i,j,vec,rate in directed_edges:
        vec=np.asarray(vec,dtype=float);rate=float(rate)
        if not (0<=i<nsites and 0<=j<nsites) or vec.shape!=(3,) or not np.all(np.isfinite(vec)) or not np.isfinite(rate) or rate<0:
            raise ValueError('invalid edge')
        if rate:
            edges.append((i,j,vec,rate));adjacency[i].add(j);adjacency[j].add(i)
    supplied=None
    if equilibrium is not None:
        supplied=np.asarray(equilibrium,dtype=float)
        if supplied.shape!=(nsites,) or not np.all(np.isfinite(supplied)) or np.any(supplied<0) or supplied.sum()<=0:
            raise ValueError('finite nonnegative probability with positive sum required')
        supplied=supplied/supplied.sum()
    if not edges:
        return np.zeros((3,3)),supplied if supplied is not None else np.full(nsites,1/nsites)
    positive=[e[3] for e in edges]
    decades=np.log10(max(positive))-np.log10(min(positive))
    if decades>100:
        raise ValueError('rate dynamic range exceeds supported 100-decade input limit')
    # Double precision SVD can erase physically real slow bridges. Construct the
    # generator from individual rates at high precision before any cancellation.
    # Rare uphill probabilities may multiply along an entire component, so
    # single-edge dynamic range alone does not bound cancellation in LU.
    digits=max(50,int(np.ceil(decades))*(nsites+1)+35)
    if digits>2000:
        raise ValueError('required working precision exceeds 2000-digit safety limit')
    with mp.workdps(digits):
        rate_scale=mp.mpf(str(max(positive)))
        Q=mp.matrix(nsites);drift=mp.matrix(nsites,3)
        converted=[]
        for i,j,vec,rate in edges:
            r=mp.mpf(str(rate))/rate_scale
            v=[mp.mpf(str(float(x))) for x in vec]
            Q[i,j]+=r;Q[i,i]-=r
            for a in range(3): drift[i,a]+=r*v[a]
            converted.append((i,j,v,r))
        components=[];unseen=set(range(nsites))
        while unseen:
            stack=[min(unseen)];component=set()
            while stack:
                i=stack.pop()
                if i in component: continue
                component.add(i);stack.extend(adjacency[i]-component)
            unseen-=component;components.append(sorted(component))
        pi=mp.matrix(nsites,1);chi=mp.matrix(nsites,3)
        tolerance=mp.mpf('1e-25')
        for ids in components:
            m=len(ids)
            if supplied is None:
                # Unspecified component populations: equal probability per site.
                weight=mp.mpf(m)/nsites
                A=mp.matrix([[Q[j,i] for j in ids] for i in ids]);b=mp.matrix(m,1)
                for j in range(m): A[m-1,j]=1
                b[m-1]=weight
                try: pc=mp.lu_solve(A,b)
                except ZeroDivisionError as error:
                    raise ValueError('nonergodic directed component; provide a valid reversible graph') from error
                for k,i in enumerate(ids): pi[i]=pc[k]
            else:
                for i in ids: pi[i]=mp.mpf(str(float(supplied[i])))
            for i in ids:
                if pi[i]<0: raise ValueError('invalid stationary distribution')
                incoming=sum(pi[j]*Q[j,i] for j in ids)
                local_scale=sum(abs(pi[j]*Q[j,i]) for j in ids)
                if local_scale and abs(incoming)>mp.mpf('1e-12')*local_scale:
                    raise ValueError('invalid stationary distribution')
            for a in range(3):
                bias=sum(pi[i]*drift[i,a] for i in ids)
                bias_scale=sum(abs(pi[i]*r*v[a]) for i,j,v,r in converted if i in ids)
                if bias_scale and abs(bias)>mp.mpf('1e-12')*bias_scale:
                    raise ValueError('nonzero ballistic drift')
            if m>1:
                # Gauge chi[first]=0. Solve a principal minor without SVD cutoff.
                A=mp.matrix([[Q[i,j] for j in ids[1:]] for i in ids[1:]])
                B=mp.matrix([[-drift[i,a] for a in range(3)] for i in ids[1:]])
                for a in range(3):
                    try: x=mp.lu_solve(A,B[:,a])
                    except ZeroDivisionError as error:
                        raise ValueError('corrector precision insufficient; no trustworthy D returned') from error
                    for k,i in enumerate(ids[1:]): chi[i,a]=x[k]
            for i in ids:
                for a in range(3):
                    residual=sum(Q[i,j]*chi[j,a] for j in ids)+drift[i,a]
                    local=sum(abs(Q[i,j]*chi[j,a]) for j in ids)+abs(drift[i,a])
                    if local and abs(residual)>tolerance*local:
                        raise ValueError('corrector residual failed')
        D=mp.matrix(3)
        for i,j,vec,r in converted:
            v=[vec[a]+chi[j,a]-chi[i,a] for a in range(3)]
            for a in range(3):
                for b in range(3): D[a,b]+=pi[i]*r*v[a]*v[b]*rate_scale/2
        return np.array(D.tolist(),dtype=float),np.array(list(pi),dtype=float)
