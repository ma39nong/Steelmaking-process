import numpy as np
import pytest
from problem2.src.delay_audit_core import scan, group_windows, error_metrics
from problem2.src.delay_v2 import local_delay_reliable
from problem2.src.delay_audit import synthetic_cases, compare_tree


def signals(delay=18):
    rng=np.random.default_rng(71)
    x=np.convolve(rng.normal(size=1000),np.ones(11)/11,mode='same')
    y=np.interp(np.arange(1000)-delay,np.arange(1000),x)+rng.normal(0,.01,1000)
    return x,y


def test_scan_matches_production_diagnostics_and_masks():
    x,y=signals(); old=local_delay_reliable(x,y,.5); new=scan(x,y)
    for key in old:
        np.testing.assert_allclose(new[key],old[key],rtol=1e-10,atol=1e-12,equal_nan=True)


@pytest.mark.parametrize('lag',[-110,-60,-30,-10,0,10,30,60,110])
def test_signed_known_delays(lag):
    x,y=signals(lag); centers=np.arange(180,820,5)
    a=scan(x,y,signed=True,centers=centers)
    assert a['reliable'].mean()>.95
    assert np.nanmax(abs(a['delay_observed']-lag*.5))<=.5
    b=scan(x,y,centers=centers)
    np.testing.assert_array_equal(a['center_index'],b['center_index'])
    if lag==0: assert not b['reliable'].any()


@pytest.mark.parametrize('lag',[-120,120])
def test_signed_endpoints_rejected(lag):
    x,y=signals(lag); s=scan(x,y,signed=True)
    assert s['boundary_hit'].mean()>.95
    assert s['reliable'].mean()<.05


def test_groups_are_descriptive_not_independent_events():
    centers=np.array([60,65,90,300]); left=centers-60; right=centers+180
    ids=group_windows(centers,left,right,np.ones(4,bool))
    assert [len(set(x)) for x in ids]==[3,2,1]


def test_metrics_include_rejections_and_wrong_direction():
    m=error_metrics(np.array([5.,np.nan,-4.]),np.array([5.,5.,5.]))
    assert m['accepted_windows']==2
    assert m['acceptance_rate']==2/3
    assert m['mae_s']==4.5
    assert m['wrong_direction_rate_accepted_nonzero']==.5
    assert m['within_1s_rate_eligible']==1/3


def test_broadened_response_has_distinct_truths():
    cfg=dict(fixed_delays_s=[0],noise_ratios=[.05,.5])
    cases={name:(x,y,truth,eligible) for name,x,y,truth,eligible,_ in synthetic_cases(1000,.5,101,cfg)}
    truth=cases['broadened_theta_5s'][2]
    assert np.all(truth['onset']<truth['peak'])
    assert np.all(truth['peak']<truth['centroid'])
    assert not cases['adjacent_events'][3].all()
    assert np.isnan(cases['adjacent_events'][2]['onset']).all()


def test_invalid_support_and_shape_fail():
    x,y=signals()
    with pytest.raises(ValueError): scan(x,y,signed=True,centers=[60])
    with pytest.raises(ValueError): scan(x,y[:-1])
    assert compare_tree({'a':1.},{'a':1.000000001})==[]


def test_undefined_truth_is_not_counted_as_a_rejected_estimate():
    m=error_metrics(np.array([5.,5.]),np.array([np.nan,5.]))
    assert m["eligible_windows"]==1
    assert m["accepted_windows"]==1
    assert m["acceptance_rate"]==1
