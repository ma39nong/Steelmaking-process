"""Meaningful checks for exact A replay, shifted gating and no-data semantics."""
import json
from pathlib import Path
import numpy as np
from .core import compare,packet,shift
from ...src.delay_v2 import local_delay_reliable


def config():return json.loads((Path(__file__).parent/'config.json').read_text())


def test_A_matches_production_on_noisy_known_translation():
    c=config();rng=np.random.default_rng(7);x=np.convolve(rng.normal(size=601),np.ones(5)/5,mode='same');y=shift(x,15,.5)+rng.normal(size=601)*.03
    a=compare(x,y,c)['A'];p=local_delay_reliable(x,y,.5)
    assert np.array_equal(a['center_index'],p['center_index'])
    assert np.array_equal(a['accepted'],p['reliable'])
    np.testing.assert_allclose(a['observed'],p['delay_observed'],equal_nan=True)
    np.testing.assert_allclose(a['score'],p['max_corr'],equal_nan=True,atol=1e-11)
    np.testing.assert_allclose(a['confidence'],p['confidence'],equal_nan=True,atol=1e-11)


def test_late_true_response_rejected_by_unshifted_gate_is_recovered():
    c=config();rng=np.random.default_rng(2);x=packet(1001,440,rng);y=shift(x,30,.5)
    result=compare(1000+x,1000+y,c,[400])
    assert not result['A']['accepted'][0]
    assert not result['A']['unshifted_target_gate'][0]
    assert result['B']['accepted'][0]
    assert result['B']['observed'][0]==30
    assert result['B']['best_shifted_target_ratio'][0]>=.05


def test_B_cannot_accept_target_lag_below_its_excitation_threshold():
    c=config();rng=np.random.default_rng(18);x=packet(1001,360,rng);y=40*packet(1001,350,rng)+.001*shift(x,55,.5)
    result=compare(1000+x,1000+y,c,[400])
    assert result['A']['accepted'][0] and result['A']['observed'][0]==55
    assert result['A']['best_shifted_target_ratio'][0]<.05
    b=result['B']
    assert not b['accepted'][0] or b['best_shifted_target_ratio'][0]>=.05


def test_constants_produce_rejections_and_zero_confidence():
    c=config();x=np.ones(601)
    result=compare(x,x,c,[200])
    for r in result.values():
        assert not r['accepted'].any()
        assert np.isnan(r['observed']).all()
        assert np.array_equal(r['confidence'],[0.])


def test_shared_centers_and_endpoint_rejection():
    c=config();rng=np.random.default_rng(8);x=packet(1001,400,rng)
    for lag in [0,60]:
        result=compare(x,shift(x,lag,.5),c,[400])
        assert np.array_equal(result['A']['center_index'],result['B']['center_index'])
        assert not result['A']['accepted'][0] and not result['B']['accepted'][0]


def test_masking_cannot_improve_full_competitor_separation():
    c=config();rng=np.random.default_rng(21);x=np.convolve(rng.normal(size=1001),np.ones(5)/5,mode='same');y=rng.normal(size=1001)
    b=compare(x,y,c)['B'];valid=np.isfinite(b['separation'])&np.isfinite(b['full_competitor_separation'])
    assert np.all(b['full_competitor_separation'][valid]<=b['separation'][valid]+1e-12)
    assert np.all(~b['B_full_competitor_accept']|b['accepted'])


def test_invalid_center_rejected_before_comparison():
    import pytest
    c=config()
    with pytest.raises(ValueError):compare(np.ones(601),np.ones(601),c,[0])
