"""Bind raw KMC blocks to reproducibility checks and conditional transport results."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from run_tracer_validation import pooled_ratio
J=lambda p:json.loads(Path(p).read_text(encoding='utf-8-sig'))
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
dest=ROOT/'analysis/transport_phase34_20261015';certificate=J(dest/'network_admission_certificate.json');model=J(dest/'transport_model.json')
assert SHA(dest/'transport_model.json')==certificate['model_sha256']
for n,d in certificate['dependency_sha256'].items():assert SHA(ROOT/'scripts'/n)==d
files=[dest/f'KMC_ext_{i:02}.json' for i in range(12)]+[dest/f'KMC_{i:02}.json' for i in [10,11,12,13]]
assert all(p.exists() for p in files),'Do not admit an unfinished experiment set'
groups={};source_sha={};checks=[]
for p in files:
    r=J(p);raw=r['raw_KMC'];assert r['sampler']=='direct exchanges'
    assert SHA(r['source_model'])==r['source_model_sha256']==certificate['model_sha256']
    assert SHA(r['experiment_design'])==r['experiment_design_sha256'] and SHA(r['worker_analyzer_path'])==r['worker_analyzer_sha256']
    t=raw['block_time_s'];li=pooled_ratio(raw['block_Dtracer_A2_s'],t);vac=pooled_ratio(raw['block_Dvacancy_A2_s'],t)
    assert np.isclose(li['mean']*1e-16,r['DLi_scalar_cm2_s'],rtol=1e-12,atol=0)
    assert np.isclose(li['standard_error']*1e-16,r['DLi_standard_error_cm2_s'],rtol=1e-10,atol=0)
    assert np.isclose(vac['mean']*1e-16,r['Dvac_KMC_time_pooled_cm2_s'],rtol=1e-12,atol=0)
    assert np.isclose(np.trace(raw['Dtracer_tensor_A2_s'])/3*1e-16,r['DLi_scalar_cm2_s'],rtol=1e-12,atol=0)
    counts=np.asarray(raw['event_counts_by_directed_edge']);rates=np.asarray(r['supplied_rates_s'])
    assert counts.shape==rates.shape and counts.sum()==raw['nevents'] and raw['burn_events']==20000
    assert len(t)==raw['nblocks'] and r['Li_count']==r['nsites']-1 and r['vacancy_fraction']==1/r['nsites']
    key=(r['energy_family'],r['T_K'],r['nsites'],raw['nblocks']);groups.setdefault(key,[]).append((p,r))
    source_sha[str(p)]=SHA(p)
def pool(key):
    pairs=groups[key];assert len(pairs)==2 and {r['raw_KMC']['seed'] for p,r in pairs}=={137,991}
    rows=[r for p,r in pairs];ref=rows[0]
    times=np.concatenate([r['raw_KMC']['block_time_s'] for r in rows]);lis=np.concatenate([r['raw_KMC']['block_Dtracer_A2_s'] for r in rows]);vacs=np.concatenate([r['raw_KMC']['block_Dvacancy_A2_s'] for r in rows])
    li=pooled_ratio(lis,times);vac=pooled_ratio(vacs,times);D=li['mean']*1e-16;se=li['standard_error']*1e-16
    V=vac['mean']*1e-16;vse=vac['standard_error']*1e-16
    seed_z=abs(rows[0]['DLi_scalar_cm2_s']-rows[1]['DLi_scalar_cm2_s'])/np.hypot(rows[0]['DLi_standard_error_cm2_s'],rows[1]['DLi_standard_error_cm2_s'])
    rates=np.asarray(ref['supplied_rates_s']);counts=sum(np.asarray(r['raw_KMC']['event_counts_by_directed_edge']) for r in rows)
    per_seed_zero=[int(np.sum((np.asarray(r['raw_KMC']['event_counts_by_directed_edge'])==0)&(rates>0))) for r in rows]
    pop=ref['population_ratio'];dv=ref['Dvac_exact_cm2_s'];du=ref['Duncorr_cm2_s']
    return dict(family=key[0],T_K=key[1],nsites=key[2],blocks_per_seed=key[3],seeds=[137,991],events_per_seed=ref['raw_KMC']['nevents'],total_elapsed_time_s=float(times.sum()),DLi_cm2_s=D,DLi_SE_cm2_s=se,relative_SE=se/D,Dvac_KMC_cm2_s=V,Dvac_KMC_SE_cm2_s=vse,Dvac_exact_cm2_s=dv,Duncorr_cm2_s=du,vacancy_exact_z=abs(V-dv)/vse,seed_agreement_z=float(seed_z),unsampled_positive_edges_pooled=int(np.sum((counts==0)&(rates>0))),unsampled_positive_edges_per_seed=per_seed_zero,vacancy_fraction=ref['vacancy_fraction'],population_ratio=pop,f_v=dv/du,f_Li=D/(pop*du),f_eff=D/(pop*dv),f_eff_SE=se/(pop*dv),raw_counts_by_family=sum(np.asarray(r['actual_observed_counts_by_family']) for r in rows).tolist(),files=[str(p) for p,r in pairs],block_pooling='Independent seeds only; nested block-length variants never pooled together',statistical_error_scope='Ratio jackknife across long event blocks; does not include DFT, defect concentration, model completeness, or prefactor uncertainty')
pooled={str(key):pool(key) for key in groups};primary=pooled[str(('Efr',300.,32,160))]
comparisons=[];admission=[]
for sites in [32,128]:
    long=pooled[str(('Efr',300.,sites,160))];short=pooled[str(('Efr',300.,sites,320))]
    z=abs(long['DLi_cm2_s']-short['DLi_cm2_s'])/np.hypot(long['DLi_SE_cm2_s'],short['DLi_SE_cm2_s'])
    comparisons.append(dict(nsites=sites,relative_change_short_to_long=long['DLi_cm2_s']/short['DLi_cm2_s']-1,conservative_combined_SE_z=float(z),paired_trajectory=True,scope='Paired blocking diagnostic, not an independent hypothesis test'))
    admission.extend([dict(check=f'{sites}_block_agreement',passed=bool(z<=3)),dict(check=f'{sites}_primary_relative_SE',passed=bool(long['relative_SE']<=.05))])
    for label,row in [('long',long),('short',short)]:
        admission.extend([dict(check=f'{sites}_{label}_seed_agreement',passed=bool(row['seed_agreement_z']<=3)),dict(check=f'{sites}_{label}_vacancy_exact',passed=bool(row['vacancy_exact_z']<=3)),dict(check=f'{sites}_{label}_edge_coverage',passed=row['unsampled_positive_edges_pooled']==0),dict(check=f'{sites}_{label}_reversible_upper_bound',passed=bool(row['DLi_cm2_s']<=row['population_ratio']*row['Duncorr_cm2_s']+3*row['DLi_SE_cm2_s']))])
box=pooled[str(('Efr',300.,128,160))];box_comparison=dict(D128_over_D32=box['DLi_cm2_s']/primary['DLi_cm2_s'],population_ratio_128_over_32=box['population_ratio']/primary['population_ratio'],fLi_128_over_32=box['f_Li']/primary['f_Li'],normalized_difference_z=abs(box['f_Li']-primary['f_Li'])/np.hypot(box['DLi_SE_cm2_s']/(box['population_ratio']*box['Duncorr_cm2_s']),primary['DLi_SE_cm2_s']/(primary['population_ratio']*primary['Duncorr_cm2_s'])),scope='One vacancy in each box with the same64-host DFT barriers: kinetic finite-concentration/box comparison, not fixed-cv or DFT-size convergence. Resolved difference is retained, not forced to disappear.')
temps=[pooled[str(('Efr',T,32,160 if T in [250.,300.] else 40))] for T in [250.,300.,350.,400.]]
T=np.array([r['T_K'] for r in temps]);y=np.log([r['DLi_cm2_s'] for r in temps]);sy=np.array([r['DLi_SE_cm2_s']/r['DLi_cm2_s'] for r in temps]);X=np.c_[np.ones(4),1/T]
W=np.diag(1/sy**2);M=np.linalg.inv(X.T@W@X)@X.T@W;beta=M@y;residual=y-X@beta;kb=8.617333262145e-5
fit=dict(T_range_K=[250,400],apparent_Li_Ea_eV=float(-kb*beta[1]),diagonal_fit_SE_eV=float(kb*np.sqrt(np.linalg.inv(X.T@W@X)[1,1])),conservative_SE_upper_bound_for_unknown_cross_temperature_covariance_eV=float(kb*np.sum(abs(M[1])*sy)),weighted_residual_chi2=float(np.sum((residual/sy)**2)),degrees_of_freedom=2,logD_residuals=residual.tolist(),scope='Descriptive250-400K fit, not a300K local LiEa proof. Same numeric seeds acrossT may correlate estimates; diagonal SE is not claimed as an independent-temperature uncertainty. The absolute-coefficient SE bound is given conservatively; non-Arrhenius curvature remains possible.',temperature_rows=temps)
E0=pooled[str(('E0',300.,32,160))]
result=dict(status='FINITE_MODEL_KMC_STATISTICAL_GATES_PASSED' if all(c['passed'] for c in admission) else 'FINITE_MODEL_KMC_GATES_NOT_ALL_PASSED',criteria=J(dest/'KMC_design.json')['criteria'],checks=admission,primary_300K=primary,box_300K=box,blocking_comparisons=comparisons,box_comparison=box_comparison,temperature_fit=fit,E0_sensitivity_300K=E0,E0_over_Efr_DLi=E0['DLi_cm2_s']/primary['DLi_cm2_s'],pooled_estimates=pooled,source_raw_KMC_sha256=source_sha,model_sha256=SHA(dest/'transport_model.json'),network_certificate_sha256=SHA(dest/'network_admission_certificate.json'),analyzer_path=str(Path(__file__).resolve()),analyzer_sha256=SHA(__file__),pooling_analyzer_path=str(ROOT/'scripts/run_tracer_validation.py'),pooling_analyzer_sha256=SHA(ROOT/'scripts/run_tracer_validation.py'),intrinsic_D_reported=False,local_300K_Li_Ea_proven=False,frequency_sensitivity={str(v):dict(DLi300_cm2_s=primary['DLi_cm2_s']*v/1e13,DLi300_SE_cm2_s=primary['DLi_SE_cm2_s']*v/1e13,scope='Exact time-rescaling of the same assumed constant common frequency; no new physical prefactor evidence') for v in [1e12,1e13,1e14]})
out=dest/'KMC_statistical_admission.json';assert not out.exists();out.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(dict(status=result['status'],failed_checks=[c for c in admission if not c['passed']],primary_300K=primary,box_comparison=box_comparison,apparent_Li_Ea_fit=fit),indent=2))
