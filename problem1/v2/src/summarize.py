import numpy as np,pandas as pd,json
from common import *

def aggregate(df,group,name):
 # First average scenarios within each heat only when explicitly forming overall summaries.
 metrics=[c for c in ['rmse','mae','nrmse','sam_deg','peak_position_nm','peak_height_relative','area_relative','negative_prediction_fraction','overshoot_fraction','valid_observation_rmse_not_prediction','seconds','converged'] if c in df]
 heat=df.groupby(['heat_id']+group)[metrics].mean();heat.to_csv(P/f'results/{name}_by_heat.csv');heat.groupby(group).agg(['mean','std']).to_csv(P/f'results/{name}_macro.csv');return heat

def main():
 d=pd.concat([pd.read_csv(P/'results/test_completion.csv'),pd.read_csv(P/'results/c_expanded_test.csv'),pd.read_csv(P/'results/v1_matched_test.csv')],ignore_index=True);aggregate(d,['model','scenario'],'all_models_test');aggregate(d,['model'],'all_models_overall')
 # Paired differences use exactly same heat, stage, seed and scene.
 a=d[d.model=='B2'].set_index(['heat_id','stage','seed','scenario']);b=d[d.model=='V1_B3_matched'].set_index(['heat_id','stage','seed','scenario']);pairs=pd.DataFrame(dict(B2_rmse=a.rmse,V1_rmse=b.rmse));pairs['rmse_difference']=pairs.B2_rmse-pairs.V1_rmse;pairs['relative_improvement']=1-pairs.B2_rmse/pairs.V1_rmse;pairs.to_csv(P/'results/v1_v2_paired_panel.csv');pairs.groupby(['heat_id','scenario']).mean().groupby('scenario').agg(['mean','std']).to_csv(P/'results/v1_v2_paired_panel_macro.csv')
 for f,group,name in [('interpolation_test.csv',['method','scenario'],'interpolation_independent'),('artificial_censoring.csv',['method'],'censor_overall'),('real_saturation_comparison.csv',['method'],'real_censor_overall'),('full_confirmation.csv',['model','scenario'],'full_confirmation'),('full_interpolation.csv',['method','scenario'],'full_interpolation')]:
  path=P/'results'/f
  if path.exists():aggregate(pd.read_csv(path),group,name)
 d=pd.read_csv(P/'results/detector_test.csv');N=int(np.prod(CFG['panel_shape']));truth=d.tp+d.fn;raw_fp=np.rint(d.total_deleted_fraction*N).astype(int)-d.tp;d['raw_false_positive_rate']=raw_fp/np.maximum(N-truth,1);d['raw_precision']=d.tp/np.maximum(d.tp+raw_fp,1);d['raw_f1']=2*d.raw_precision*d.recall/np.maximum(d.raw_precision+d.recall,1e-12);d.to_csv(P/'results/detector_test_with_raw_FPR.csv',index=False)
 d.groupby(['heat_id','method','kind'])[['precision','recall','f1','false_positive_rate','raw_false_positive_rate','raw_precision','raw_f1','total_deleted_fraction','incremental_reference_exclusion']].mean().groupby(['method','kind']).agg(['mean','std']).to_csv(P/'results/detector_test_detailed_macro.csv')
 # Heat-level bootstrap of RMSE differences; repeated masks stay within each heat cluster.
 heat=pairs.groupby('heat_id').rmse_difference.mean().to_numpy();rng=np.random.default_rng(407);means=np.mean(rng.choice(heat,size=(10000,len(heat)),replace=True),axis=1);dump(P/'results/paired_heat_bootstrap.json',dict(n_heats=len(heat),mean_rmse_difference=float(heat.mean()),cluster_bootstrap_percentile95=[float(x) for x in np.quantile(means,[.025,.975])],interpretation='Ten fixed heats only, no independent external-furnace generalization; repeated masks clustered by heat'))
if __name__=='__main__':main()
