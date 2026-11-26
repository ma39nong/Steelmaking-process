from core import *
from detectors import all_masks,reconstruct
from experiment_final import metrics

def main():
 rows=[]
 for split in ['development','confirmation']:
  for h,b,y,t,w,i,j in panels(split):
   ref=y.copy();ref[:,64]=np.minimum(64000,ref[:,64]+5000);masks,e=all_masks(ref,t,w,h)
   for name,mask in masks.items():
    D,sol=reconstruct(ref,t,w,mask);rows.append(dict(split=split,heat_id=h,band=b,strategy=name,line_cells_masked=int(mask[:,64].sum()),**metrics(ref,ref,D,np.zeros(y.shape,bool),mask,w),**sol))
 pd.DataFrame(rows).to_csv(P/'results/persistent_line_controls.csv',index=False)
if __name__=='__main__':main()
