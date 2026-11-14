"""Finite single-vacancy tracer experiment; statistical gates stay explicit."""
from pathlib import Path
import argparse,json,hashlib,itertools
import numpy as np
from ase.io import read
from ase.build import make_supercell
from hop_geometry import classify_hops
from network_arrays import supercell_arrays
from transport import rates_from_barriers
from compute_round1_transport import network_fields,read_barriers
from tracer_kmc import simulate
from dimer_kmc import simulate_accelerated
ROOT=Path(__file__).resolve().parents[1]

def validate_model(model):
 if model['status']!='CI_CANDIDATE_FINITE_MODEL_PENDING_FINAL_SENSITIVITY':
  raise ValueError('use the registered four-CI finite-model report')
 frequency=float(model['attempt_frequency_s'])
 if not np.isfinite(frequency) or frequency<=0:raise ValueError('positive finite assumed frequency required')
 evidence=model['barrier_evidence']
 rows,cell=read_barriers([dict(path_id=r['path_id'],source=r['source']) for r in evidence])
 for old,new in zip(evidence,rows):
  if old!=new:raise ValueError('CI source evidence or provenance changed; regenerate transport model')
 for label in ['Efr','E0']:
  actual=[r['barriers_sym_eV'][label] for r in rows]
  if model['families'][label]['barriers_sym_eV']!=actual:
   raise ValueError('derived barrier arrays differ from revalidated CI evidence')
 return rows

def pooled_ratio(block_values,block_times):
 values=np.asarray(block_values,dtype=float);times=np.asarray(block_times,dtype=float)
 if len(values)!=len(times) or len(values)<2 or not np.all(np.isfinite(values)) or not np.all(np.isfinite(times)) or np.any(times<=0):
  raise ValueError('at least two finite blocks with positive times required')
 num=values*times;mean=float(num.sum()/times.sum())
 leave=(num.sum()-num)/(times.sum()-times)
 error=float(np.sqrt((len(values)-1)/len(values)*np.sum((leave-leave.mean())**2)))
 return dict(mean=mean,standard_error=error,elapsed_time_s=float(times.sum()))

def build_box(replication,rates):
 replication=np.asarray(replication,dtype=int)
 if replication.shape!=(3,) or np.any(replication<1):raise ValueError('positive integer replication required')
 a=read(ROOT/'analysis/host_symmetry_projection/POSCAR');orbits=classify_hops(a,3.3,1e-5)['orbits']
 base=json.loads((ROOT/'analysis/host_supercell64/supercell_mapping.json').read_text())
 P=np.diag(replication)@np.array(base['supercell_matrix'])
 a.new_array('basis_index',np.arange(len(a)))
 s=make_supercell(a,P,order='atom-major');indices=s.arrays['basis_index']
 images=s.positions@np.linalg.inv(a.cell)-a.get_scaled_positions()[indices]
 if np.max(abs(images-np.rint(images)))>1e-7:raise ValueError('noninteger primitive image mapping')
 mapping=dict(supercell_matrix=P.tolist(),basis_indices=indices.tolist(),primitive_images=np.rint(images).astype(int).tolist())
 arrays=supercell_arrays(a,s,mapping,orbits,rates)
 return a,s,mapping,arrays

def select_sampler(path_ids,rates,mode):
 if mode not in ['auto','direct','dimer']:raise ValueError('unknown sampler mode')
 pid=np.asarray(path_ids);rates=np.asarray(rates)
 fastest=int(np.argmax(rates));fast=None
 if np.all(np.sum(pid==fastest,axis=1)==1):fast=np.argmax(pid==fastest,axis=1)
 if mode=='dimer' and fast is None:raise ValueError('fastest family is not a one-neighbor dimer family')
 escape=float(np.sum(rates[pid[0]])-rates[fastest])
 accelerated=mode=='dimer' or (mode=='auto' and fast is not None and rates[fastest]>10*escape)
 return accelerated,fastest,fast

def experiment(model_file,temperature,replication,nevents,nblocks,seed,mode='auto',energy_family='Efr'):
 model_file=Path(model_file).resolve();model=json.loads(model_file.read_text())
 if nblocks<2:raise ValueError('tracer validation requires at least two blocks')
 validate_model(model)
 for evidence in model['barrier_evidence']:
  file=Path(evidence['source'])
  if hashlib.sha256(file.read_bytes()).hexdigest()!=evidence['source_sha256']:
   raise ValueError('CI source audit changed; regenerate the transport model')
 if hashlib.sha256((ROOT/'analysis/host_symmetry_projection/POSCAR').read_bytes()).hexdigest()!=model['primitive_structure_sha256']:
  raise ValueError('transport structure source changed')
 if hashlib.sha256((ROOT/'analysis/host_supercell64/supercell_mapping.json').read_bytes()).hexdigest()!=model['supercell_mapping_sha256']:
  raise ValueError('transport cell mapping source changed')
 barriers=model['families'][energy_family]['barriers_sym_eV'];frequency=model['attempt_frequency_s']
 rates=rates_from_barriers(barriers,temperature,frequency)
 a,s,mapping,arrays=build_box(replication,rates);orbits=classify_hops(a,3.3,1e-5)['orbits']
 exact,upper=network_fields(a,orbits,rates);Dv=float(np.trace(exact)/3);Du=float(np.trace(upper)/3)
 neighbor,disp,rr=arrays['neighbor'],arrays['displacement'],arrays['rates'];pid=arrays['path_id']
 accelerated,fastest,fast=select_sampler(pid,rates,mode)
 if accelerated:
  if fast is None:raise ValueError('fastest family does not have one reciprocal fast neighbor per site')
  result=simulate_accelerated(neighbor,disp,rr,fast,nevents,nblocks,seed)
  counts=np.array(result['escape_event_counts_by_edge']);sampler='exact dimer escape skeleton'
 else:
  result=simulate(neighbor,disp,rr,nevents,nblocks,seed)
  counts=np.array(result['event_counts_by_directed_edge']);sampler='direct exchanges'
 vac=pooled_ratio(result['block_Dvacancy_A2_s'],result['block_time_s'])
 nsites=neighbor.shape[0];population=1/(nsites-1)
 Dli=result['Dtracer_mean_A2_s'];SE=result['block_standard_error_A2_s']
 period=min(np.linalg.norm(np.array(v)@s.cell) for v in itertools.product(range(-3,4),repeat=3) if v!=(0,0,0))
 return dict(status='FINITE_TIME_ESTIMATE_PENDING_SEED_LAG_BOX_COMPARISON',
  source_model=str(model_file),source_model_sha256=hashlib.sha256(model_file.read_bytes()).hexdigest(),
  source_model_status=model['status'],
  energy_family=energy_family,T_K=temperature,attempt_frequency_s=frequency,barriers_eV=barriers,
  sampler=sampler,fast_family=fastest if accelerated else None,nsites=nsites,Li_count=nsites-1,
  vacancy_fraction=1/nsites,population_ratio=population,replication=list(replication),cell_A=s.cell.tolist(),
  shortest_periodic_translation_A=float(period),network_mapping=mapping,
  supplied_neighbor=neighbor.tolist(),supplied_displacement_A=disp.tolist(),supplied_rates_s=rr.tolist(),
  actual_observed_counts_by_family=[int(counts[pid==i].sum()) for i in range(4)],
  observed_count_scope='Escape counts; fast events are marginalized, parity/expected counts are separate' if accelerated else 'Physical exchange events',
  DLi_scalar_cm2_s=Dli*1e-16,DLi_standard_error_cm2_s=SE*1e-16,
  Dvac_exact_cm2_s=Dv*1e-16,Duncorr_cm2_s=Du*1e-16,
  Dvac_KMC_time_pooled_cm2_s=vac['mean']*1e-16,Dvac_KMC_SE_cm2_s=vac['standard_error']*1e-16,
  Dvac_KMC_relative_difference=vac['mean']/Dv-1,
  f_v=Dv/Du,f_Li=Dli/(population*Du),f_eff=Dli/(population*Dv),
  f_eff_standard_error=SE/(population*Dv),f_eff_bound_assumed=False,
  correlation_scope='trace/3 ratios; f_eff includes normalization to already correlated Dvac and need not be <=1',
  raw_KMC=result,intrinsic_D_reported=False,
  box_scope='One vacancy in every box; changing replication changes concentration, not a fixed-cv DFT size check')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('model');p.add_argument('destination')
 p.add_argument('--temperature',type=float,default=300.);p.add_argument('--replication',type=int,nargs=3,default=[1,1,1])
 p.add_argument('--events',type=int,default=1000000);p.add_argument('--blocks',type=int,default=40);p.add_argument('--seed',type=int,required=True)
 p.add_argument('--mode',choices=['auto','direct','dimer'],default='auto');p.add_argument('--energy-family',choices=['Efr','E0'],default='Efr')
 args=p.parse_args();dest=Path(args.destination)
 if dest.exists():raise ValueError('use a fresh output file; previous experiments are immutable')
 result=experiment(args.model,args.temperature,args.replication,args.events,args.blocks,args.seed,args.mode,args.energy_family)
 dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(result,indent=2),encoding='utf-8')
 print(json.dumps({k:result[k] for k in ['status','sampler','T_K','nsites','DLi_scalar_cm2_s','DLi_standard_error_cm2_s','Dvac_KMC_relative_difference','f_eff','actual_observed_counts_by_family']},indent=2))
