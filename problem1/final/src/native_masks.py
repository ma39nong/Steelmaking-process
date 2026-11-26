from core import *
from detectors import adaptive,local_mask,all_masks

def main():
 chosen=json.loads((P/'config/SELECTION_LOCKED.json').read_text())['strategy'];review=pd.read_csv(V31/'results/high_confidence_statistical_points.csv');veto=json.loads((V31/'config/native_review_vetoes.json').read_text())['vetoes'];rows=[];points=[]
 for h in HEATS:
  y,t,w=load(h);e=adaptive(y,t,w)
  if chosen.startswith('D_'):mask=local_mask(e,float(chosen[-1]))
  elif chosen=='B_V1_MAD':mask=~detect(y,t,w,**V1PAR[h]['detector'])[0]&(y!=65535)
  elif chosen=='C_V31':mask=conservative_evidence(y,t,w)['high']
  else:mask=np.zeros(y.shape,bool)
  auto=mask.copy();review_mask=np.zeros(y.shape,bool)
  for r in review[review.heat_id==h].itertuples():review_mask[r.row,r.column]=True
  mask|=review_mask
  for r in veto:
   if r['heat_id']==h:mask[r['row'],r['column']]=False
  valid=np.isfinite(y)&(y>=0)&(y!=65535)&~mask;suspect=(e['score']>=2)&~mask;grade=np.zeros(y.shape,np.uint8);grade[suspect]=1;grade[mask]=2;grade[y==65535]=1
  np.savez_compressed(P/f'data/{h}_masks.npz',repair=mask,automatic=auto,reviewed=review_mask,suspect=suspect,valid=valid,grade=grade,score=e['score'].astype(np.float32),time_z=e['time_z'].astype(np.float32),wave_z=e['wave_z'].astype(np.float32),lowrank_z=e['lowrank_z'].astype(np.float32),persistence=e['persistence'])
  ii,jj=np.where(mask)
  for i,j in zip(ii,jj):points.append(dict(heat_id=h,row=int(i),column=int(j),time_s=t[i],wavelength_nm=w[j],raw=y[i,j],repair_type='isolated_spike' if y[i,j]>e['time_reference'][i,j] else 'local_drop',score=e['score'][i,j],time_z=e['time_z'][i,j],wave_z=e['wave_z'][i,j],lowrank_z=e['lowrank_z'][i,j],persistence=int(e['persistence'][i,j]),consistency=e['consistency'][i,j],reviewed_local_continuity=bool(review_mask[i,j]),hardware_fault_confirmed=False))
  rows.append(dict(heat_id=h,observations=y.size,automatic_before_veto=int(auto.sum()),reviewed_points=int(review_mask.sum()),repair_points=int(mask.sum()),repair_fraction=float(mask.mean()),suspect_retained=int(suspect.sum()),top_codes=int((y==65535).sum()),whole_frames_deleted=int(mask.all(1).sum())));pd.DataFrame(rows).to_csv(P/'results/native_mask_summary.csv',index=False);print('mask',h,int(mask.sum()),flush=True)
 pd.DataFrame(points).to_csv(P/'results/statistical_repairs_evidence.csv',index=False)
if __name__=='__main__':main()
