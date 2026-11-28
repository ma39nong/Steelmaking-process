import json
from pathlib import Path
import numpy as np
import pytest
from ...src.q2_pipeline import hampel_mask
from .operators import fast_hampel,pulse,acquire,source_methods,gas_methods
from .metrics import metrics
from .integral import generate,fit

CFG=json.loads(Path(__file__).with_name('config.json').read_text())


@pytest.mark.parametrize('window,threshold',[(9,6),(9,12),(5,6)])
def test_vectorized_baseline_is_exact_on_noise_plateaus_edges_and_invalids(window,threshold):
    rng=np.random.default_rng(7)
    signals=[rng.normal(10,1,151),np.repeat([1.,2.,3.,2.,0.],31),np.zeros(151)]
    signals[0][[0,15,40,100,150]]=[30,np.nan,np.inf,-1,100]
    signals[1][75]=99.;signals[2][75]=1.
    for x in signals:assert np.array_equal(fast_hampel(x,window,threshold),hampel_mask(x,window,threshold))


def test_flow_pulse_is_integrated_into_monotone_cumulative_volume():
    t=np.arange(1001)*.5
    for shape,volume in [('triangle',25.),('rectangle',50.)]:
        q,o=pulse(t,200.15,10.,5.,shape)
        assert np.all(np.diff(o)>=0) and o[-1]==pytest.approx(volume)
        assert o[t<195.15].max()==0 and q[t<195.15].max()==0


def test_quantized_pair_applies_measurement_to_whole_augmented_signal():
    t=np.arange(10)*.5;base=np.full(10,.006);add=np.full(10,.006)
    correct=acquire(base+add,t,'grid_0.01')-acquire(base,t,'grid_0.01')
    wrong=acquire(add,t,'grid_0.01')
    assert np.all(correct==0) and np.all(wrong==.01)


def test_empirical_hold_is_causal_and_preserves_its_observed_background():
    t=np.arange(8)*.5;raw=np.array([0,0,2,2,2,5,5,7.])
    assert np.array_equal(acquire(raw,t,'empirical_hold',raw),raw)
    ramp=acquire(t,t,'empirical_hold',raw)
    assert np.all(ramp<=t)


def test_single_component_source_ablations_share_mask_and_clean_signal():
    t=np.arange(401)*.5;x=100+5*t;x[200]+=100
    r=source_methods(t,x,CFG['source_methods'])
    for name in ['SG11','analytic_PCHIP_derivative','secant_2s','omit_flow_clipping']:
        assert np.array_equal(r[name]['mask'],r['formal']['mask'])
        assert np.array_equal(r[name]['clean'],r['formal']['clean'])
    assert np.array_equal(r['formal']['flow'],np.maximum(r['omit_flow_clipping']['flow'],0))


def test_small_reconstructed_peak_cannot_define_absolute_onset_by_rescaling_itself():
    t=np.arange(801)*.5;q,_=pulse(t,200.,4.,1.)
    r=metrics(t,.05*q,200.,4.,1.,'triangle',CFG)
    assert r['peak_ratio']==pytest.approx(.05)
    assert not r['absolute_onset_defined'] and r['onset_absolute_error_s'] is None
    assert r['onset_relative_error_s'] is not None


@pytest.mark.parametrize('tau,T',[(5.,2.),(15.,8.),(30.,2.)])
def test_integral_probe_recovers_delay_and_time_constant_numerically(tau,T):
    t,o,y=generate(CFG,tau,T,'matched',8101);r,_,_=fit(t,o,y,CFG)
    assert r['accepted'] and abs(r['estimate_s']-tau)<=.5
    assert abs(r['T_s']-T)<.25 and abs(r['gain']-.8)<.02


def test_constant_input_has_no_identified_delay():
    t,o,y=generate(CFG,15.,8.,'constant',8101);r,_,_=fit(t,o,y,CFG)
    assert not r['accepted'] and 'insufficient_input_variation' in r['rejection_reasons']
    assert r['profile_width_s']==60.


def test_good_prediction_does_not_guarantee_correct_pure_delay_under_model_mismatch():
    t,o,y=generate(CFG,15.,8.,'two_stage',8101);r,_,_=fit(t,o,y,CFG)
    assert r['validation_improvement']>.8
    assert abs(r['estimate_s']-15.)>1.
