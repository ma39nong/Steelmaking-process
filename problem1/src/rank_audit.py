"""Descriptive low-rank evidence; never used to tune final selection from outer scores."""
import numpy as np,pandas as pd
from sklearn.utils.extmath import randomized_svd
import matplotlib.pyplot as plt
from audit import P,HEATS
from anomalies import fill_training
from eda import save
rows=[];fig,ax=plt.subplots(figsize=(8,5))
for h in HEATS:
 d=np.load(P/f'data/processed/{h}.npz');y,t,w=d['Y'],d['time'],d['wavelength'];mask=np.load(P/f'data/masks/{h}.npz')['valid'];z=fill_training(np.where(mask,y,np.nan),t,w)
 _,s,_=randomized_svd(z,n_components=32,n_iter=4,random_state=20261126);energy=np.cumsum(s*s)/np.sum(z*z)
 centered=z-z.mean(0);_,sc,_=randomized_svd(centered,n_components=32,n_iter=4,random_state=20261126);ec=np.cumsum(sc*sc)/np.sum(centered*centered)
 rows.append(dict(heat_id=h,raw_rank1_energy=energy[0],raw_rank4_energy=energy[3],raw_rank12_energy=energy[11],centered_rank4_energy=ec[3],centered_rank12_energy=ec[11]));ax.plot(np.arange(1,33),energy,label=h,lw=.9)
ax.set(xlabel='Rank',ylabel='Cumulative uncentered energy fraction',ylim=(.8,1.001));ax.legend(fontsize=7,ncol=2);save(fig,P/'figures/comparison/singular_energy');pd.DataFrame(rows).to_csv(P/'results/eda/low_rank_energy.csv',index=False)
