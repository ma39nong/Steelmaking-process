"""Numerical acceptance and honest rejection, not existence-only tests."""
import numpy as np
import pytest
from .run import configuration,signals
from .model import Design,bootstrap
from ..signal_fidelity.integral import fit

@pytest.mark.parametrize('tau,T',[(5,2),(15,8),(30,2),(30,8)])
def test_weak_numerical_recovery(tau,T):
    cfg=configuration();t,o,y=signals(cfg,tau,T,'matched',9001)
    r,*_=Design(t,o,cfg).fit(y)
    assert r['accepted'] and abs(r['estimate_s']-tau)<=.5
    assert abs(r['T_s']-T)<.3 and abs(r['gain']-.8)<.03


def test_ordinary_replay_and_no_training_validation_overlap():
    cfg=configuration();t,o,y=signals(cfg,15,8,'matched',9001)
    a,*_=fit(t,o,y,cfg);d=Design(t,o,cfg,'ordinary');b,*_=d.fit(y)
    for k in ['estimate_s','T_s','gain','offset','validation_rmse','profile_width_s']:assert b[k]==pytest.approx(a[k],abs=1e-8)
    weak=Design(t,o,cfg);support=np.flatnonzero(np.any(weak.W!=0,axis=0))
    assert t[support].max()<cfg['integral_probe']['validation_bounds_s'][0]


def test_constant_rejected_and_cumulative_offset_invariant():
    cfg=configuration();t,o,y=signals(cfg,15,8,'constant',9001)
    r,*_=Design(t,o,cfg).fit(y)
    assert not r['accepted'] and 'insufficient_input_variation' in r['rejection_reasons']
    t,o,y=signals(cfg,15,8,'matched',9001);a,*_=Design(t,o,cfg).fit(y);b,*_=Design(t,o+1e8,cfg).fit(y)
    for k in ['estimate_s','T_s','gain']:assert a[k]==pytest.approx(b[k],abs=1e-5)


def test_weak_weight_annihilates_constant_and_integrates_linear_cumulative():
    cfg=configuration();t,o,y=signals(cfg,15,8,'matched',9001);d=Design(t,25+5*t,cfg)
    assert np.max(abs(d.D@np.ones(len(t))))<1e-12
    assert np.allclose(d.W@np.ones(len(t)),1)
    # Trapezoid quadrature error is explicit and vanishes on refinement.
    assert np.max(abs(d.Q-5))<.011


def test_bootstrap_deterministic_and_not_overlapping_row_resampling():
    cfg=configuration();cfg['bootstrap_replicates']=9
    t,o,y=signals(cfg,15,8,'matched',9001);y+=np.random.default_rng(9001).normal(0,.05,len(t))
    d=Design(t,o,cfg);r,*_=d.fit(y);a,x=bootstrap(d,y,r,cfg,42);b,z=bootstrap(d,y,r,cfg,42)
    assert a==b and np.array_equal(x,z) and x.shape==(9,4)
    assert a['ci_lower_s']<=a['ci_upper_s']


def test_unrelated_output_rejected_without_forced_event_pairing():
    cfg=configuration();t,o,y=signals(cfg,15,8,'unrelated',9001)
    r,*_=Design(t,o,cfg).fit(y)
    assert not r['accepted'] and 'no_heldout_improvement' in r['rejection_reasons']


def test_model_mismatch_good_prediction_does_not_validate_pure_delay():
    cfg=configuration();t,o,y=signals(cfg,15,8,'two_stage',9001)
    r,*_=Design(t,o,cfg).fit(y)
    assert r['validation_improvement']>.8 and abs(r['estimate_s']-15)>1


def test_correct_slow_point_can_be_rejected_for_nonidentifiable_profile():
    cfg=configuration();t,o,y=signals(cfg,15,2,'slow',9001)
    r,*_=Design(t,o,cfg).fit(y)
    assert r['estimate_s']==15 and not r['accepted']
    assert 'broad_profile' in r['rejection_reasons']


def test_evaluation_seeds_are_disjoint_and_temporal_history_is_required():
    cfg=configuration();assert set(cfg['evaluation_seeds']).isdisjoint(range(8101,8111))
    t,o,y=signals(cfg,15,8,'matched',9001)
    cfg['integral_probe']['training_bounds_s']=[30,275]
    with pytest.raises(ValueError,match='history'):Design(t,o,cfg)
