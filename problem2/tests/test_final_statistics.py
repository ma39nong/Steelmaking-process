"""Numerical checks for final stage accounting and honest heat validation."""
import copy
import numpy as np
from scipy.integrate import cumulative_trapezoid
from problem2.src.final_statistics import correlations,interval_increments,stage_statistics,fit_trend,analyze


def heat(name,n=101,offset=0):
    t=np.linspace(0,10,n);o=2*t;y=3+t+offset
    return dict(heat_id=name,time=t,oxygen_cumulative_clean=o,co_clean=y,co2_clean=y*2,
        total_co_emission=np.r_[0,cumulative_trapezoid(y,t)],total_co2_emission=np.r_[0,cumulative_trapezoid(2*y,t)])


def config():
    return dict(oxygen_progress_edges=[0,.2,.4,.6,.8,1],increment_intervals_s=[.5,5,15],trend_degrees=[1,2],random_seed=5201,bootstrap_heat_resamples=30)


def test_high_cumulative_correlation_does_not_force_increment_correlation():
    t=np.arange(20.)
    assert correlations(t,5*t)['pearson_r']>.999999
    assert correlations(np.diff(t),np.diff(5*t))['pearson_r'] is None


def test_nonoverlapping_increments_drop_incomplete_end():
    t=np.arange(0,10.5,.5);x=t**2
    dx,dy=interval_increments(t,x,2*x,3)
    np.testing.assert_allclose(dx,[9,27,45]);np.testing.assert_allclose(dy,2*dx)


def test_equal_oxygen_stages_conserve_both_totals():
    d=heat('a');rows=stage_statistics(d,np.linspace(0,1,6))
    for gas in ['co','co2']:
        rr=[r for r in rows if r['gas']==gas]
        np.testing.assert_allclose([r['stage_oxygen_m3'] for r in rr],4)
        assert np.isclose(sum(r['stage_emission_m3'] for r in rr),d['total_'+gas+'_emission'][-1])
        assert np.isclose(sum(r['duration_s'] for r in rr),10)
        assert sum(r['sample_count'] for r in rr)==len(d['time'])


def test_stage_boundary_plateau_tie_is_explicit():
    d=heat('a',6);d['oxygen_cumulative_clean']=np.array([0,4,4,12,16,20.])
    rows=stage_statistics(d,np.linspace(0,1,6))
    assert rows[0]['time_end_s']==4
    assert rows[0]['time_start_s']==0


def test_equal_heat_weights_invariant_to_duplicated_samples():
    x=np.array([0.,1.,0.,1.]);y=np.array([0.,1.,10.,11.]);w=np.full(4,.5)
    base=fit_trend(x,y,1,w)
    xx=np.r_[x[:2],np.repeat(x[2:],10)];yy=np.r_[y[:2],np.repeat(y[2:],10)];ww=np.r_[np.full(2,.5),np.full(20,.05)]
    np.testing.assert_allclose(fit_trend(xx,yy,1,ww),base,atol=1e-12)


def test_heldout_target_never_enters_its_fitted_coefficients():
    ds=[heat('a'),heat('b',offset=2),heat('c',offset=-1)]
    first=analyze(ds,config());changed=copy.deepcopy(ds);changed[0]['co_clean']+=100
    second=analyze(changed,config())
    one=next(r for r in first['trend_loho'] if r['test_heat']=='a' and r['relation']=='co_instantaneous' and r['degree']==2)
    two=next(r for r in second['trend_loho'] if r['test_heat']=='a' and r['relation']=='co_instantaneous' and r['degree']==2)
    np.testing.assert_allclose([one[k] for k in ['c0','c1','c2']],[two[k] for k in ['c0','c1','c2']])
    assert two['mae']>one['mae']+90


def test_decreasing_cumulative_oxygen_rejected_not_silently_fixed():
    import pytest
    d=heat('a',6);d['oxygen_cumulative_clean'][3]=-1
    with pytest.raises(ValueError):stage_statistics(d,np.linspace(0,1,6))


def test_default_full_workflow_finishes_with_current_delivery(monkeypatch):
    import sys
    from problem2.src import run_all
    calls=[]
    monkeypatch.setattr(sys,'argv',['run_all'])
    monkeypatch.setattr(run_all,'run',lambda module,args:calls.append(module))
    run_all.main()
    assert calls[-1]=='problem2.src.final_delivery'
    assert 'problem2.src.visualization' not in calls
