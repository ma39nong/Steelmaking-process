from pathlib import Path
import sys,json
import numpy as np
from r03_log_prox import log_radius_prox,column_log_prox,singular_log_prox
from scipy.special import huber
R=Path(__file__).resolve().parents[1];rng=np.random.default_rng(20261126);rows=[]
for k in range(100):
 a=float(rng.uniform(0,12));alpha=float(rng.uniform(0,15));x=float(log_radius_prox(a,alpha));grid=np.linspace(0,max(a,1e-9),20001)
 obj=lambda z:.5*(z-a)**2+alpha*np.log1p(z)
 gap=float(obj(x)-np.min(obj(grid)));assert gap<1e-6
 rows.append(dict(a=a,alpha=alpha,prox=x,objective=obj(x),grid_minimum_gap=gap))
y=rng.normal(size=(13,7));z=column_log_prox(y,2);assert np.isfinite(z).all()
assert np.linalg.norm(column_log_prox(y,0)-y)<1e-12
assert np.linalg.norm(singular_log_prox(y,0)-y)<1e-12
r=np.linspace(-6,6,101);delta=1.345;ours=np.where(abs(r)<=delta,.5*r*r,delta*abs(r)-.5*delta*delta);assert np.max(abs(ours-huber(delta,r)))<1e-12
(R/"metadata/method_reproduction_checks.json").write_text(json.dumps({"r03_scope":"Core proximal operators (22)-(23); not full local HSI/TV algorithm","r03_trials":rows,"huber_supplement":"SciPy official definition","huber_max_error":float(np.max(abs(ours-huber(delta,r))))},indent=2))
print('100 R03 scalar proximal checks; column and singular zero-penalty identities; Huber definition passed')
