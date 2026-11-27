import json
from pathlib import Path
import numpy as np
import pytest
from .core import KernelDesign,changes,descriptors,candidate_match,evaluate_matches,detect_events,adjust_bh
from .benchmark import event_case,response_kernel


@pytest.fixture
def cfg():return json.loads(Path(__file__).with_name('config.json').read_text())


def test_kernel_descriptors_distinguish_effective_onset_peak_centroid(cfg):
    d=descriptors(response_kernel('asymmetric_gamma',.5,60),.5)
    assert d['onset']<d['peak']<d['centroid']
    assert d['width']>10
    assert descriptors(np.zeros(121),.5)['peak'] is None


def test_recovery_of_broadened_kernel_and_no_test_label_leakage(cfg):
    rng=np.random.default_rng(81);u=rng.normal(size=1801);h=response_kernel('asymmetric_gamma',.5,60)
    y=np.convolve(u,h)[:len(u)]+rng.normal(0,.001,len(u))
    design=KernelDesign(u,cfg);a=design.fit(y);changed=y.copy();changed[design.test]+=1000;b=design.fit(changed)
    np.testing.assert_array_equal(a['kernel_weights'],b['kernel_weights'])
    assert a['alpha']==b['alpha']
    assert abs(a['descriptors']['peak']-descriptors(h,.5)['peak'])<=1
    assert abs(a['descriptors']['centroid']-descriptors(h,.5)['centroid'])<=1
    assert a['test_skill']>.95
    assert b['test_skill']<a['test_skill']


def ev(peak):return dict(peak=peak,start=peak-2,centroid=peak+1,width=5,end=peak+4,polarity=1,confidence=.9)


def test_match_preserves_shared_response_ambiguity_and_unmatched(cfg):
    rows=candidate_match([ev(10),ev(20),ev(200)],[ev(40)],cfg)
    assert [r['status'] for r in rows]==['ambiguous','ambiguous','unmatched']
    assert all(r['target_index'] is None for r in rows)
    rows=candidate_match([ev(10)],[ev(40)],cfg)
    assert rows[0]['status']=='matched' and rows[0]['delay_peak']==30


def test_reverse_response_is_not_forced_into_positive_delay(cfg):
    assert candidate_match([ev(50)],[ev(40)],cfg)[0]['status']=='unmatched'


def test_zero_signal_does_not_invent_events(cfg):
    assert detect_events(np.zeros(500),cfg)==[]


def test_detect_noise_threshold_is_training_only(cfg):
    t=np.arange(500);x=3*np.exp(-.5*((t-100)/3)**2)+2*np.exp(-.5*((t-400)/3)**2)
    a=detect_events(x,cfg,np.arange(200));x[450]=1000
    b=detect_events(x,cfg,np.arange(200))
    assert a and b
    assert all(e['threshold']==a[0]['threshold'] for e in b)


def test_bh_and_null_metrics(cfg):
    np.testing.assert_allclose(adjust_bh([.01,.04,.9]),[.03,.06,.9])
    rows=candidate_match([ev(10)],[ev(40)],cfg)
    m=evaluate_matches(rows,[ev(10)],[ev(40)],[dict(source_peak=10,gas_peak=None)],0,cfg)
    assert m['tp']==0 and m['fp']==1 and m['recall'] is None
    assert m['null_match_fraction']==1


def test_new_benchmark_has_no_pure_translation(cfg):
    for name in cfg['synthetic']['scenarios']:
        case=event_case(name,1001,cfg)
        assert len(case['truth'])>=9
        if case['kernel_truth'] is not None:assert descriptors(case['kernel_truth'],.5)['width']>5


def test_short_records_rejected(cfg):
    with pytest.raises(ValueError):KernelDesign(np.ones(150),cfg)


def test_v2_adapter_does_not_invent_onset_and_uses_v2_delay(cfg):
    v2=dict(center_time=np.array([10.]),reliable=np.array([True]),delay_observed=np.array([27.]))
    r=candidate_match([ev(10)],[ev(40)],cfg,method="v2",v2=v2)[0]
    assert r["status"]=="matched"
    assert r["delay_peak"]==30 and r["estimator_peak"]==27
    assert r["estimator_onset"] is None and r["estimator_centroid"] is None


def test_surrogates_preserve_power_spectrum():
    from .run import phase_surrogate
    rng=np.random.default_rng(7);y=rng.normal(size=1000)
    other=phase_surrogate(y,rng)
    np.testing.assert_allclose(abs(np.fft.rfft(other)),abs(np.fft.rfft(y)),atol=1e-10)
