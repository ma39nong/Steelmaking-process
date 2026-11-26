"""Prespecified persistent-line false-repair counterexample, no fault labels inferred."""
from common31 import *
from points import evidence,repair,mask_for
from experiments import stats

def main(split):
 rows=[]
 for h,y,t,w,i,j in panels(split):
  reference=y.copy();reference[:,128]=np.minimum(64000,reference[:,128]+5000)
  e=evidence(reference,t,w)
  for name in CFG['strategies']:
   D,s,_=repair(reference,t,w,e,name);mask=mask_for(e,name);rows.append(dict(heat_id=h,column_start=j,strategy=name,**stats(reference,reference,D,np.zeros(y.shape,bool),mask,w),line_cells_masked=int(mask[:,128].sum()),**s))
 pd.DataFrame(rows).to_csv(P/f'results/persistent_line_control_{split}.csv',index=False)
if __name__=='__main__':main(sys.argv[1])
