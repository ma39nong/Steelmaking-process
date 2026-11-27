"""Reliability-aware local delay estimation and constrained model fitting."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import peak_prominences, peak_widths

from .q2_pipeline import _standardized_change


def local_delay_reliable(
    source: np.ndarray,
    target: np.ndarray,
    time_step: float,
    window: int = 121,
    max_lag: int = 120,
    stride: int = 5,
    min_corr: float = 0.35,
    min_separation: float = 0.05,
    min_strength_ratio: float = 0.05,
    peak_distance: int = 5,
) -> dict[str, np.ndarray]:
    """Estimate positive lags and reject weak, ambiguous, or boundary peaks.

    A center is valid when the source window is fully observed and the delayed
    target window is fully observed. Consequently centers start at ``half``
    and end before ``n-half-max_lag``. No boundary value is extrapolated.
    """
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    if source.ndim != 1 or target.ndim != 1 or source.size != target.size:
        raise ValueError("source and target must be one-dimensional arrays of equal length")
    n = source.size
    half = window // 2
    if window % 2 == 0 or n <= window + max_lag:
        raise ValueError("window must be odd and the series must support the requested lag range")
    src = _standardized_change(source)
    dst = _standardized_change(target)
    src_global = float(np.nanstd(src))
    dst_global = float(np.nanstd(dst))
    fields = {
        "center_index": np.array([], dtype=int),
        "center_time": np.array([], dtype=float),
        "delay_observed": np.array([], dtype=float),
        "max_corr": np.array([], dtype=float),
        "source_std": np.array([], dtype=float),
        "target_std": np.array([], dtype=float),
        "peak_separation": np.array([], dtype=float),
        "peak_prominence": np.array([], dtype=float),
        "peak_width_samples": np.array([], dtype=float),
        "boundary_hit": np.array([], dtype=bool),
        "reliable": np.array([], dtype=bool),
        "confidence": np.array([], dtype=float),
        "support_left_index": np.array([], dtype=int),
        "support_right_index": np.array([], dtype=int),
    }
    records: list[dict[str, Any]] = []
    if not np.isfinite(src_global) or not np.isfinite(dst_global) or src_global == 0 or dst_global == 0:
        return fields
    for center in range(half, n - half - max_lag, stride):
        left = center - half
        right = center + half + 1
        a = src[left:right]
        b0 = dst[left:right]
        a_std = float(np.nanstd(a))
        b0_std = float(np.nanstd(b0))
        valid_strength = a_std >= min_strength_ratio * src_global and b0_std >= min_strength_ratio * dst_global
        correlations = np.full(max_lag + 1, np.nan)
        if valid_strength and np.isfinite(a).all():
            a_centered = (a - np.mean(a)) / (a_std + 1e-12)
            for lag in range(max_lag + 1):
                b = dst[left + lag : right + lag]
                b_std = float(np.nanstd(b))
                if b_std > 0 and np.isfinite(b).all():
                    correlations[lag] = float(np.mean(a_centered * ((b - np.mean(b)) / b_std)))
        if not np.isfinite(correlations).any():
            best_lag, best_corr = 0, np.nan
            separation, prominence, width = np.nan, np.nan, np.nan
        else:
            best_lag = int(np.nanargmax(correlations))
            best_corr = float(correlations[best_lag])
            candidates = np.flatnonzero(np.isfinite(correlations))
            competitor = candidates[np.abs(candidates - best_lag) >= peak_distance]
            second = float(np.max(correlations[competitor])) if competitor.size else best_corr
            separation = best_corr - second
            if 0 < best_lag < max_lag:
                prominence = float(peak_prominences(correlations, [best_lag])[0][0])
                width = float(peak_widths(correlations, [best_lag], rel_height=0.5)[0][0])
            else:
                prominence, width = 0.0, float("inf")
        boundary = best_lag in (0, max_lag)
        reliable = bool(
            valid_strength
            and np.isfinite(best_corr)
            and best_corr >= min_corr
            and np.isfinite(separation)
            and separation >= min_separation
            and not boundary
        )
        strength = min(a_std / (src_global + 1e-12), b0_std / (dst_global + 1e-12), 1.0)
        confidence = float(np.clip((best_corr - min_corr) / (1 - min_corr), 0, 1) * np.clip(separation / 0.25, 0, 1) * strength) if np.isfinite(best_corr) else 0.0
        records.append({
            "center_index": center,
            "center_time": center * time_step,
            "delay_observed": best_lag * time_step if reliable else np.nan,
            "max_corr": best_corr,
            "source_std": a_std,
            "target_std": b0_std,
            "peak_separation": separation,
            "peak_prominence": prominence,
            "peak_width_samples": width,
            "boundary_hit": boundary,
            "reliable": reliable,
            "confidence": confidence if reliable else 0.0,
            "support_left_index": left,
            "support_right_index": right + max_lag - 1,
        })
    for key in fields:
        fields[key] = np.asarray([row[key] for row in records], dtype=fields[key].dtype)
    return fields


def _pava(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted pool-adjacent-violators for a nondecreasing fit."""
    values: list[float] = []
    weights: list[float] = []
    counts: list[int] = []
    for value, weight in zip(x, w):
        values.append(float(value))
        weights.append(float(weight))
        counts.append(1)
        while len(values) >= 2 and values[-2] > values[-1]:
            total_weight = weights[-2] + weights[-1]
            values[-2] = (values[-2] * weights[-2] + values[-1] * weights[-1]) / total_weight
            weights[-2] = total_weight
            counts[-2] += counts[-1]
            values.pop(); weights.pop(); counts.pop()
    return np.repeat(values, counts)


def fit_delay_models(flow: np.ndarray, delay: np.ndarray, confidence: np.ndarray, max_delay_s: float) -> dict[str, Any]:
    valid = np.isfinite(flow) & np.isfinite(delay) & np.isfinite(confidence) & (flow > 0) & (confidence > 0)
    if valid.sum() < 5:
        return {"n": int(valid.sum()), "selected": "insufficient_data"}
    q = flow[valid]
    y = delay[valid]
    w = np.maximum(confidence[valid], 1e-3)
    x = 1.0 / np.sqrt(q)
    fixed = float(np.average(y, weights=w))
    unconstrained = np.linalg.lstsq(np.column_stack([np.ones_like(x), x]), y, rcond=None)[0]

    def residual(theta: np.ndarray) -> np.ndarray:
        return np.sqrt(w) * (theta[0] + theta[1] * x - y)

    constrained = least_squares(residual, x0=np.array([max(0.0, fixed / 2), max(0.0, fixed)]), bounds=([0.0, 0.0], [max_delay_s, max_delay_s * np.sqrt(max(q))]))
    pred = np.clip(constrained.x[0] + constrained.x[1] * x, 0, max_delay_s)
    fixed_rmse = float(np.sqrt(np.average((y - fixed) ** 2, weights=w)))
    dynamic_rmse = float(np.sqrt(np.average((y - pred) ** 2, weights=w)))
    order = np.argsort(q)
    iso_sorted = _pava(y[order][::-1], w[order][::-1])
    iso_pred = np.empty_like(iso_sorted)
    iso_pred[np.argsort(q)[::-1]] = iso_sorted
    iso_rmse = float(np.sqrt(np.average((y - iso_pred) ** 2, weights=w)))
    selected = "constrained_inverse_sqrt" if dynamic_rmse + 1e-9 < fixed_rmse else "fixed_median"
    return {
        "n": int(valid.sum()),
        "selected": selected,
        "fixed": {"c_s": fixed, "rmse_s": fixed_rmse},
        "unconstrained_inverse_sqrt": {"a_s": float(unconstrained[0]), "b_s_sqrt_m3": float(unconstrained[1])},
        "constrained_inverse_sqrt": {"a_s": float(constrained.x[0]), "b_s_sqrt_m3": float(constrained.x[1]), "rmse_s": dynamic_rmse},
        "monotone_flow_fit": {"rmse_s": iso_rmse},
        "prediction_min_s": float(np.min(pred)),
        "prediction_max_s": float(np.max(pred)),
    }
