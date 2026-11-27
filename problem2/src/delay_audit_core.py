"""Experimental delay diagnostics; production delay models remain unchanged."""
from __future__ import annotations

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.signal import peak_prominences, peak_widths
from .q2_pipeline import _standardized_change


def scan(source, target, dt=0.5, window=121, max_lag=120, stride=5,
         min_corr=0.35, min_separation=0.05, min_strength_ratio=0.05,
         signed=False, centers=None):
    """Same correlation/acceptance as V2; signed search is experimental only.

    Explicit centers must have full support in both series. Signed mode treats
    zero as an interior lag, whereas V2 rejects zero as a search endpoint.
    Target excitation uses the *unshifted* target window, matching V2 exactly.
    """
    source, target = np.asarray(source, float), np.asarray(target, float)
    if source.ndim != 1 or source.shape != target.shape or not np.isfinite([source, target]).all():
        raise ValueError('expected equal, finite one-dimensional signals')
    if dt <= 0 or window < 3 or window % 2 == 0 or max_lag < 1 or stride < 1:
        raise ValueError('invalid sampling or window parameters')
    half = window // 2
    lags = np.arange(-max_lag if signed else 0, max_lag + 1)
    lo, hi = half - lags[0], len(source) - half - lags[-1]
    if centers is None:
        centers = np.arange(lo, hi, stride, dtype=int)
    else:
        centers = np.asarray(centers, int)
    if len(centers) == 0 or np.any(centers < lo) or np.any(centers >= hi):
        raise ValueError('centers have insufficient support')
    src, dst = _standardized_change(source), _standardized_change(target)
    sw, tw = sliding_window_view(src, window), sliding_window_view(dst, window)
    left = centers - half
    a, b0 = sw[left], tw[left]
    astd, bstd = a.std(axis=1), b0.std(axis=1)
    sg, tg = src.std(), dst.std()
    strength = (sg > 0) & (tg > 0) & (astd >= min_strength_ratio*sg) & (bstd >= min_strength_ratio*tg)
    ac = (a - a.mean(axis=1)[:, None]) / (astd[:, None] + 1e-12)
    corr = np.full((len(centers), len(lags)), np.nan)
    for j, lag in enumerate(lags):
        b = tw[left + lag]
        std = b.std(axis=1)
        valid = strength & (std > 0)
        if valid.any():
            bc = (b[valid] - b[valid].mean(axis=1)[:, None]) / std[valid, None]
            corr[valid, j] = np.mean(ac[valid] * bc, axis=1)
    finite = np.isfinite(corr).any(axis=1)
    best = np.argmax(np.where(np.isfinite(corr), corr, -np.inf), axis=1)
    score = corr[np.arange(len(centers)), best]
    separation = np.full(len(centers), np.nan)
    prominence, width = separation.copy(), separation.copy()
    boundary = (best == 0) | (best == len(lags)-1)
    for i in np.flatnonzero(finite):
        other = np.isfinite(corr[i]) & (np.abs(lags - lags[best[i]]) >= 5)
        separation[i] = score[i] - (np.max(corr[i, other]) if other.any() else score[i])
        if not boundary[i]:
            prominence[i] = peak_prominences(corr[i], [best[i]])[0][0]
            width[i] = peak_widths(corr[i], [best[i]], rel_height=0.5)[0][0]
        else:
            prominence[i], width[i] = 0., np.inf
    checks = {'strength': strength, 'finite_correlation': finite,
              'correlation': finite & (score >= min_corr),
              'separation': np.isfinite(separation) & (separation >= min_separation),
              'interior_peak': ~boundary & finite}
    accepted = np.logical_and.reduce(list(checks.values()))
    best_lag = lags[best]
    confidence = np.clip((score-min_corr)/(1-min_corr), 0, 1)*np.clip(separation/.25, 0, 1)*np.minimum(np.minimum(astd/(sg+1e-12), bstd/(tg+1e-12)), 1)
    confidence = np.where(accepted, confidence, 0.)
    return dict(center_index=centers, center_time=centers*dt, lags_s=lags*dt,
                correlations=corr, best_delay_s=best_lag*dt,
                delay_observed=np.where(accepted, best_lag*dt, np.nan),
                max_corr=score, source_std=astd, target_std=bstd,
                peak_separation=separation, peak_prominence=prominence,
                peak_width_samples=width, boundary_hit=boundary,
                reliable=accepted, confidence=confidence, checks=checks,
                support_left_index=left+lags[0],
                support_right_index=left+window-1+lags[-1])


def group_windows(centers, left, right, reliable, stride=5):
    """Assign three descriptive groups, never an independence claim."""
    ids = np.full((3, len(centers)), -1, dtype=int)
    half = (centers[0]-left[0]) if len(centers) else 0
    for kind in range(3):
        group, last_i, end = -1, None, -1
        for i in np.flatnonzero(reliable):
            l = centers[i]-half if kind == 1 else left[i]
            r = centers[i]+half if kind == 1 else right[i]
            same = last_i is not None and ((centers[i]-centers[last_i] == stride) if kind == 0 else l <= end)
            if not same:
                group += 1
            ids[kind, i] = group
            end = max(end, r) if same else r
            last_i = i
    return ids


def error_metrics(observed, truth, eligible=None):
    observed, truth = np.asarray(observed), np.asarray(truth)
    eligible = np.ones(observed.shape, bool) if eligible is None else np.asarray(eligible, bool)
    # Undefined estimands (e.g. Gaussian onset) are not evaluable truths.
    eligible = eligible & np.isfinite(truth)
    valid = eligible & np.isfinite(observed)
    error = observed[valid]-truth[valid]
    absolute = np.abs(error)
    denom = int(eligible.sum())
    nonzero = valid & (truth != 0)
    wrong = int(np.sum(np.sign(observed[nonzero]) != np.sign(truth[nonzero])))
    return dict(eligible_windows=denom, accepted_windows=int(valid.sum()),
                acceptance_rate=float(valid.sum()/denom) if denom else None,
                rejection_rate=float(1-valid.sum()/denom) if denom else None,
                wrong_direction_count=wrong,
                wrong_direction_rate_accepted_nonzero=float(wrong/nonzero.sum()) if nonzero.any() else None,
                bias_s=float(error.mean()) if len(error) else None,
                mae_s=float(absolute.mean()) if len(error) else None,
                median_absolute_error_s=float(np.median(absolute)) if len(error) else None,
                p90_absolute_error_s=float(np.quantile(absolute,.9)) if len(error) else None,
                within_1s_rate_eligible=float(np.sum(absolute<=1)/denom) if denom else None)
