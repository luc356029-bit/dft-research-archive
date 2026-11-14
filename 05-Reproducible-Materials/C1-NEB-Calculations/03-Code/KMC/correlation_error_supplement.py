"""Use covariance-free SE bounds for paired box/block comparisons."""
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];dest=ROOT/'analysis/transport_phase34_20261015'
J=lambda p:json.loads(Path(p).read_text(encoding='utf-8-sig'))
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
source=dest/'KMC_statistical_admission.json';a=J(source)
assert SHA(a['analyzer_path'])==a['analyzer_sha256']
for p,d in a['source_raw_KMC_sha256'].items():assert SHA(p)==d
rows=[]
for sites in [32,128]:
    l=a['pooled_estimates'][str(('Efr',300.,sites,160))];s=a['pooled_estimates'][str(('Efr',300.,sites,320))]
    bound=l['DLi_SE_cm2_s']+s['DLi_SE_cm2_s']
    assert bound>=np.hypot(l['DLi_SE_cm2_s'],s['DLi_SE_cm2_s'])
    z=abs(l['DLi_cm2_s']-s['DLi_cm2_s'])/bound
    rows.append(dict(nsites=sites,SE_of_difference_upper_bound_cm2_s=bound,difference_over_SE_upper_bound=z,passed_3SE=bool(z<=3),covariance_assumption='No cross-blocking independence assumed; triangle/Cauchy bound SE(delta)<=SE(a)+SE(b)'))
l=a['primary_300K'];s=a['box_300K'];se_l=l['DLi_SE_cm2_s']/(l['population_ratio']*l['Duncorr_cm2_s']);se_s=s['DLi_SE_cm2_s']/(s['population_ratio']*s['Duncorr_cm2_s'])
box=dict(fLi_difference_32_minus_128=l['f_Li']-s['f_Li'],SE_difference_upper_bound=se_l+se_s,difference_over_SE_upper_bound=abs(l['f_Li']-s['f_Li'])/(se_l+se_s),normalized_ratio_128_over_32=s['f_Li']/l['f_Li'],scope='Common seeds can correlate boxes. The original9.2447 combined-SE z is descriptive only; this bound does not assume independent box estimates. It still resolves a finite-concentration/box difference, not DFT size convergence.')
assert all(r['passed_3SE'] for r in rows) and box['difference_over_SE_upper_bound']>3
result=dict(status='PAIRED_COMPARISONS_VALIDATED_WITH_COVARIANCE_FREE_ERROR_BOUND',source_statistical_audit=str(source),source_statistical_audit_sha256=SHA(source),paired_block_checks=rows,paired_box_comparison=box,temperature_fit_scope=a['temperature_fit']['scope'],temperature_SE_bound_eV=a['temperature_fit']['conservative_SE_upper_bound_for_unknown_cross_temperature_covariance_eV'],original_audit_preserved=True,analyzer_path=str(Path(__file__).resolve()),analyzer_sha256=SHA(__file__))
file=dest/'correlation_error_supplement.json';assert not file.exists();file.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
