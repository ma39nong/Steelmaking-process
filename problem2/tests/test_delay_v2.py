import numpy as np

from problem2.src.delay_v2 import local_delay_reliable
from problem2.src.q2_pipeline import hampel_mask


def test_hampel_zero_mad_flags_isolated_spike_and_keeps_step():
    spike = np.array([10, 10, 10, 100, 10, 10, 10], dtype=float)
    step = np.array([10, 10, 10, 100, 100, 100, 100], dtype=float)
    assert hampel_mask(spike, window=5, threshold=6)[3]
    assert not hampel_mask(step, window=5, threshold=6)[3]


def test_constant_source_is_rejected():
    source = np.ones(500)
    target = np.sin(np.linspace(0, 25, 500))
    result = local_delay_reliable(source, target, 0.5)
    assert result["center_index"].size == 0


def test_known_fixed_delay_is_recovered():
    rng = np.random.default_rng(7)
    source = np.convolve(rng.normal(size=800), np.ones(11) / 11, mode="same")
    lag = 18
    target = np.roll(source, lag) + rng.normal(0, 0.03, len(source))
    result = local_delay_reliable(source, target, 0.5, max_lag=60, min_corr=0.25)
    observed = result["delay_observed"][result["reliable"]]
    assert observed.size > 20
    assert np.median(np.abs(observed - lag * 0.5)) <= 1.0


def test_unrelated_signals_have_low_acceptance():
    rng = np.random.default_rng(9)
    source = rng.normal(size=800)
    target = rng.normal(size=800)
    result = local_delay_reliable(source, target, 0.5, max_lag=60)
    assert result["reliable"].mean() < 0.25


def test_slowly_varying_delay_has_valid_centers():
    rng = np.random.default_rng(13)
    kernel_size = 9
    source = np.convolve(rng.normal(size=1000), np.ones(kernel_size) / kernel_size, mode="same")
    time = np.arange(len(source), dtype=float)
    lag = 12 + 0.01 * time
    target = np.interp(time - lag, time, source, left=source[0], right=source[-1])
    result = local_delay_reliable(source, target, 0.5, max_lag=60, min_corr=0.25)
    assert result["reliable"].sum() > 0


def test_boundary_peaks_are_flagged():
    source = np.sin(np.linspace(0, 80, 800))
    target = np.roll(source, 30)
    result = local_delay_reliable(source, target, 0.5, max_lag=30, min_corr=0.25)
    assert result["boundary_hit"].any()
