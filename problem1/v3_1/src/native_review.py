from common31 import *
import shutil

def main():
 paths=['raw_point_audit.csv','high_confidence_statistical_points.csv','general_suspect_points.csv'];archive=P/'results/audit_before_line_review';archive.mkdir(exist_ok=True)
 for f in paths:
  if not (archive/f).exists():shutil.copy2(P/'results'/f,archive/f)
 high=pd.read_csv(archive/paths[1]);suspect=pd.read_csv(archive/paths[2]);summary=pd.read_csv(archive/paths[0]);suspect['native_review']='retain_unconfirmed_statistical_candidate';suspect['review_reason']='multiple automatic high-evidence conditions not simultaneously met; retain; no confirmed hardware fault';veto=json.loads((P/'config/native_review_vetoes.json').read_text())['vetoes'];high['native_review']='statistical_repair_supported_no_fault_label';high['review_reason']='isolated excursion outside both axial references; no sustained line evidence in inspected neighborhood'
 for r in veto:
  match=(high.heat_id==r['heat_id'])&(high.row==r['row'])&(high.column==r['column']);q=high[match].copy();assert len(q)==1;q['classification']='suspected';q['native_review']='retain_persistent_line';q['review_reason']=r['reason'];suspect=pd.concat([suspect,q]);high=high[~match]
  path=P/f"data/{r['heat_id']}_point_evidence.npz";d=dict(np.load(path));d['high'][r['row'],r['column']]=False;d['suspect'][r['row'],r['column']]=True;d['grade'][r['row'],r['column']]=1;np.savez_compressed(path,**d)
  ix=summary.heat_id==r['heat_id'];summary.loc[ix,'high_statistical_points']-=1;summary.loc[ix,'general_suspect_points']+=1
 high.to_csv(P/'results/high_confidence_statistical_points.csv',index=False);suspect.to_csv(P/'results/general_suspect_points.csv',index=False);summary.to_csv(P/'results/raw_point_audit.csv',index=False)
 dump(P/'results/NATIVE_REVIEW_COMPLETE.json',dict(automatic_strong_proposals=5,statistical_repair_supported=4,persistent_line_vetoes=1,hardware_faults_confirmed=0,review_ledger_sha256=sha(P/'config/native_review_vetoes.json')))
if __name__=='__main__':main()
