import numpy as np
import pytest
from .core import held,pulse,features,errors,integrate,reconstruct,runs,delay_features


def test_hold_phase_has_no_future_access():
    x=np.arange(15.)
    assert np.array_equal(held(x,4,1),[0,1,1,1,1,5,5,5,5,9,9,9,9,13,13])
    assert np.all(held(x,4,1)<=x)
    with pytest.raises(ValueError):held(x,2,2)


def test_known_shift_moves_all_three_definitions():
    t=np.arange(501)*.5;x=pulse(t,80.,15.);y=pulse(t,97.,15.)
    a,b=features(t,x),features(t,y)
    for k in ['onset10_s','peak_s','centroid_s']:
        assert b[k]-a[k]==pytest.approx(17.)
    result=delay_features(t,x,y,17.)
    for k in ['onset','peak','centroid']:
        assert result['feature_delay_'+k+'_error_s']==pytest.approx(0.)


def test_rectangle_peak_is_a_set_not_arbitrary_first_maximum():
    t=np.arange(501)*.5;x=pulse(t,80.,20.,'rectangle');y=x.copy();y[t==83.]=1.01
    assert errors(t,x,y)['peak_set_distance_s']==0.


def test_runs_do_not_treat_missing_values_as_held_data():
    assert runs(np.array([1.,1.,np.nan,np.nan,2.,2.]))==[(0,2),(2,3),(3,4),(4,6)]


def test_sg_short_pulse_loses_peak_despite_small_integral_error():
    t=np.arange(1001)*.5;q=5*pulse(t,200.,2.);o=integrate(t,q)
    recovered,_,_=reconstruct(t,o);r=errors(t,q,recovered)
    assert r['peak_ratio']<.4
    assert abs(integrate(t,recovered)[-1]/o[-1]-1)<.3
    assert r['onset10_s_error']<0


def test_constant_flow_reconstruction_remains_correct():
    t=np.arange(301)*.5;q=np.full_like(t,5.);r,oc,mask=reconstruct(t,100+integrate(t,q))
    assert not mask.any()
    assert np.allclose(r,q,atol=1e-10)


def test_valid_sharp_gas_peak_can_be_flagged_by_frozen_cleaner():
    from ...src.q2_pipeline import hampel_mask,pchip_clean
    t=np.arange(1001)*.5;y=5*pulse(t,215.,2.)
    mask=hampel_mask(y);clean=pchip_clean(t,y,mask)
    assert np.flatnonzero(mask).tolist()==[430]
    assert clean.max()/y.max()==pytest.approx(.5)
    assert features(t,clean)['centroid_s']==pytest.approx(215.)


def test_change_representation_erases_slow_but_real_triangular_input():
    from ...src.q2_pipeline import _standardized_change
    t=np.arange(1001)*.5
    slow=pulse(t,200.,16.);fast=pulse(t,200.,4.)
    assert np.max(slow)>0 and np.std(np.diff(slow))>0
    assert np.all(_standardized_change(slow)==0)
    assert np.allclose(_standardized_change(fast),np.diff(fast,prepend=fast[0]))
