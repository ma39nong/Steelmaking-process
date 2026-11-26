from common31 import *
from points import evidence

def runs(mask):
 for i,row in enumerate(mask):
  changes=np.diff(np.r_[False,row,False].astype(int));starts=np.flatnonzero(changes==1);ends=np.flatnonzero(changes==-1)-1
  for a,b in zip(starts,ends):yield i,int(a),int(b)
def main():
 summary=[];high=[];suspect=[];retained=[]
 for h in HEATS:
  y,t,w=load(h);e=evidence(y,t,w)
  np.savez_compressed(P/f'data/{h}_point_evidence.npz',candidate=e['candidate'],high=e['high'],suspect=e['suspect'],grade=e['grade'])
  for tag,mask,out in [('high_statistical',e['high'],high),('suspected',e['suspect'],suspect)]:
   ix,jx=np.where(mask)
   for i,j in zip(ix,jx):
    out.append(dict(heat_id=h,row=int(i),column=int(j),time_s=t[i],wavelength_nm=w[j],raw=y[i,j],classification=tag,time_z=e['time_z'][i,j],wave_z=e['wave_z'][i,j],lowrank_residual_z=e['residual_z'][i,j],temporal_reference=e['temporal_reference'][i,j],wave_reference=e['wave_reference'][i,j],neighbor_agreement=bool(e['neighbor_agreement'][i,j]),outside_neighbor_range=bool(e['exterior'][i,j]),persistent_neighbor_frames=int(e['persistence'][i,j]),component_cells=int(e['component_size'][i,j]),near_saturation=bool(e['near_saturation'][i,j]),hardware_fault_confirmed=False))
  for i,a,b in runs(~e['candidate']):retained.append(dict(heat_id=h,row=i,time_s=t[i],column_start=a,column_end=b,wave_start_nm=w[a],wave_end_nm=w[b],count=b-a+1,classification='retained_including_separately_marked_top_codes'))
  summary.append(dict(heat_id=h,observations=y.size,statistical_candidates=int(e['candidate'].sum()),high_statistical_points=int(e['high'].sum()),general_suspect_points=int(e['suspect'].sum()),retained_points=int((~e['candidate']).sum()),sat=int((y==65535).sum()),hardware_faults_confirmed=0))
  pd.DataFrame(summary).to_csv(P/'results/raw_point_audit.csv',index=False);print('audit',h,summary[-1],flush=True)
 pd.DataFrame(high).to_csv(P/'results/high_confidence_statistical_points.csv',index=False);pd.DataFrame(suspect).to_csv(P/'results/general_suspect_points.csv',index=False);pd.DataFrame(retained).to_csv(P/'results/retained_points_rle.csv',index=False)
if __name__=='__main__':main()
