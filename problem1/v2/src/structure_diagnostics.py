"""Auxiliary interpretations, never used to alter locked training masks."""
import numpy as np,pandas as pd,json
from scipy.ndimage import median_filter
from common import *

def main():
 summary=[]
 for h in HEATS:
  y,t,w=load(h);base=(y>0)&(y<65535)&np.isfinite(y);z=fill_training(np.where(base,y,np.nan),t,w)
  norm=np.maximum(np.linalg.norm(z,axis=1),1);shape=z/norm[:,None];neighbor=np.empty_like(shape);neighbor[1:-1]=(shape[:-2]+shape[2:])*.5;neighbor[0]=shape[1];neighbor[-1]=shape[-2];neighbor/=np.maximum(np.linalg.norm(neighbor,axis=1)[:,None],1e-12)
  angle=np.degrees(np.arccos(np.clip(np.sum(shape*neighbor,axis=1),-1,1)));brightness=norm/np.maximum(median_filter(norm,size=5,mode='nearest'),1)-1
  labels=np.where((abs(brightness)>.2)&(angle<2),'brightness_dominant',np.where(angle>5,'shape_change','unresolved_or_smooth'))
  pd.DataFrame(dict(time_s=t,brightness_relative=brightness,normalized_neighbor_angle_deg=angle,diagnostic=labels,automatically_deleted=False)).to_csv(P/f'results/{h}_brightness_shape.csv',index=False)
  # Contrast to neighboring wavelengths after row brightness normalization. Real narrow emission lines also qualify.
  spectral=shape-median_filter(shape,size=(1,11),mode='nearest');scale=np.maximum(1e-5,1.4826*np.median(abs(spectral),axis=1));support=abs(spectral)>8*scale[:,None];persistence=support.mean(0);signed=np.median(spectral,axis=0);candidate=(persistence>.5)&(abs(signed)>8*np.median(scale))
  pd.DataFrame(dict(wavelength_nm=w,persistence_fraction=persistence,normalized_signed_contrast=signed,candidate=candidate,automatically_deleted=False)).to_csv(P/f'results/{h}_persistent_channel_candidates.csv',index=False)
  summary.append(dict(heat_id=h,brightness_dominant_frames=int((labels=='brightness_dominant').sum()),shape_change_frames=int((labels=='shape_change').sum()),persistent_channel_candidates=int(candidate.sum()),candidate_auto_deleted=0))
 pd.DataFrame(summary).to_csv(P/'results/structure_diagnostic_summary.csv',index=False)
if __name__=='__main__':main()
