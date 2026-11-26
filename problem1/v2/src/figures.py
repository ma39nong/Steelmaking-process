import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *

def save(fig,name):fig.tight_layout();fig.savefig(P/f'figures/{name}.png',dpi=170);plt.close(fig)
def main():
 d=pd.read_csv(P/'results/convergence_curves.csv')
 for h,g in d.groupby('heat_id'):
  fig,axs=plt.subplots(2,2,figsize=(11,7))
  for name,b in g.groupby('solver'):
   for ax,col in zip(axs.ravel(),['objective','relative_change','held_rmse','seconds']):ax.plot(b.iteration,b[col],label=name);ax.set(xlabel='Iteration',ylabel=col);ax.legend()
  axs[0,1].set_yscale('log');fig.suptitle(h+' same hidden observations / same rank and eta');save(fig,h+'_convergence')
 d=pd.read_csv(P/'results/capacity_development.csv');g=d.groupby(['heat_id','setting']).rmse.mean().groupby('setting').agg(['mean','std']);fig,ax=plt.subplots(figsize=(11,5));ax.bar(g.index,g['mean'],yerr=g['std']);ax.tick_params(axis='x',rotation=30);ax.set(ylabel='RMSE: heat macro mean +/- heat SD',title='C development controlled diagnostics');save(fig,'C_capacity')
 d=pd.read_csv(P/'results/artificial_censoring.csv');g=d.groupby(['heat_id','method']).rmse.mean().groupby('method').agg(['mean','std']);fig,ax=plt.subplots(figsize=(10,5));ax.bar(g.index,g['mean'],yerr=g['std']);ax.tick_params(axis='x',rotation=20);ax.set(ylabel='Artificial clipping RMSE',title='Matched reference, three thresholds / two areas / two seeds');save(fig,'artificial_censoring')
 d=pd.read_csv(P/'results/interpolation_test.csv');g=d.groupby(['heat_id','scenario','method']).rmse.mean().groupby(['scenario','method']).mean().unstack();fig,ax=plt.subplots(figsize=(11,5));g.plot.bar(ax=ax);ax.set(ylabel='Held coordinate RMSE: heat macro mean',title='Independent continuous interpolation test');ax.tick_params(axis='x',rotation=20);save(fig,'interpolation_test')
 if (P/'results/test_completion.csv').exists():
  ds=[pd.read_csv(P/'results/test_completion.csv')]
  if (P/'results/c_expanded_test.csv').exists():ds.append(pd.read_csv(P/'results/c_expanded_test.csv'))
  if (P/'results/v1_matched_test.csv').exists():ds.append(pd.read_csv(P/'results/v1_matched_test.csv'))
  d=pd.concat(ds);g=d.groupby(['heat_id','scenario','model']).rmse.mean().groupby(['scenario','model']).mean().unstack();fig,ax=plt.subplots(figsize=(13,6));g.plot.bar(ax=ax);ax.set(ylabel='Held native-grid RMSE: heat macro mean',title='Same masks, local test panels');ax.tick_params(axis='x',rotation=25);save(fig,'completion_test')
 if (P/'results/full_confirmation.csv').exists():
  d=pd.read_csv(P/'results/full_confirmation.csv');g=d.groupby(['heat_id','scenario','model']).rmse.mean().groupby(['scenario','model']).mean().unstack();fig,ax=plt.subplots(figsize=(13,6));g.plot.bar(ax=ax);ax.set(ylabel='RMSE',title='Full native matrices confirmation; V1 matched on seed 811');ax.tick_params(axis='x',rotation=25);save(fig,'full_confirmation')
 if (P/'results/detector_test.csv').exists():
  d=pd.read_csv(P/'results/detector_test.csv');g=d.groupby(['heat_id','kind','method']).f1.mean().groupby(['kind','method']).mean().unstack();fig,ax=plt.subplots(figsize=(10,5));g.plot.bar(ax=ax);ax.set(ylabel='F1',ylim=(0,1),title='Six injected anomaly types; disjoint source panels');ax.tick_params(axis='x',rotation=20);save(fig,'detector_test')
if __name__=='__main__':main()
