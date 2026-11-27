"""Offline reconstruction, oxygen-flow and delay analysis for Problem 2."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.signal import savgol_filter

from .data_validation import ALIASES, _find_column, _norm


def _heat_id(path: Path) -> str:
    return re.search(r"Heat No\.\s*([A-Za-z0-9]+)", path.name, flags=re.I).group(1)


def read_process_file(path: Path) -> dict[str, np.ndarray | str]:
    frame = pd.read_excel(path, sheet_name=0, engine="pyxlsb")
    sequence = pd.to_numeric(frame.iloc[:, 0], errors="coerce").to_numpy(float)
    columns = {name: _find_column(list(frame.columns), aliases) for name, aliases in ALIASES.items()}
    values = {name: pd.to_numeric(frame[col], errors="coerce").to_numpy(float) for name, col in columns.items()}
    if len(sequence) == 0 or not np.allclose(sequence, np.arange(1, len(sequence) + 1)):
        raise ValueError(f"{path.name}: Sequence is not continuous")
    return {"heat_id": _heat_id(path), "sequence": sequence, "time": (sequence - 1) * 0.5, **values}


def hampel_mask(values: np.ndarray, window: int = 9, threshold: float = 6.0) -> np.ndarray:
    """Flag nonfinite, negative and robust local deviations."""
    x = np.asarray(values, float)
    mask = ~np.isfinite(x) | (x < 0)
    half = window // 2
    for i in range(len(x)):
        lo, hi = max(0, i - half), min(len(x), i + half + 1)
        local = x[lo:hi]
        local = local[np.isfinite(local)]
        if len(local) < 5:
            continue
        med = np.median(local)
        scale = 1.4826 * np.median(np.abs(local - med))
        if scale > 0:
            if abs(x[i] - med) > threshold * scale:
                mask[i] = True
            continue
        # A zero MAD is common on flat instrument plateaus. Use the two
        # one-sided neighborhoods to distinguish an isolated spike from a
        # genuine step, while keeping the exact flat-plateau test sensitive.
        left = x[lo:i]
        right = x[i + 1 : hi]
        left = left[np.isfinite(left)]
        right = right[np.isfinite(right)]
        if len(left) < 2 or len(right) < 2:
            continue
        left_med, right_med = np.median(left), np.median(right)
        neighbors = np.concatenate([left, right])
        neighbor_scale = 1.4826 * np.median(np.abs(neighbors - np.median(neighbors)))
        step_size = abs(left_med - right_med)
        if step_size <= max(1e-12, 3 * neighbor_scale):
            fallback_scale = max(neighbor_scale, 1e-12)
            if abs(x[i] - np.median(neighbors)) > threshold * fallback_scale:
                mask[i] = True
    return mask


def pchip_clean(time: np.ndarray, values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    x = np.asarray(values, float).copy()
    visible = (~mask) & np.isfinite(x)
    if visible.sum() < 2:
        raise ValueError("fewer than two visible points for reconstruction")
    result = PchipInterpolator(time[visible], x[visible], extrapolate=False)(time)
    result[~np.isfinite(result)] = np.interp(time[~np.isfinite(result)], time[visible], x[visible])
    return np.maximum(result, 0.0)


def oxygen_flow(time: np.ndarray, oxygen: np.ndarray, window: int = 31, polyorder: int = 3) -> np.ndarray:
    """Estimate dO/dt with a Savitzky-Golay derivative and enforce q >= 0."""
    n = len(oxygen)
    window = min(window, n if n % 2 else n - 1)
    window = max(polyorder + 2 + ((polyorder + 2) % 2), window)
    window = min(window, n if n % 2 else n - 1)
    derivative = savgol_filter(oxygen, window_length=window, polyorder=polyorder, deriv=1, delta=float(np.median(np.diff(time)),), mode="interp")
    return np.maximum(derivative, 0.0)


def _standardized_change(values: np.ndarray) -> np.ndarray:
    difference = np.diff(values, prepend=values[0])
    return difference - pd.Series(difference).rolling(21, center=True, min_periods=1).median().to_numpy()


def local_delay(source: np.ndarray, target: np.ndarray, time_step: float, window: int = 121, max_lag: int = 120, stride: int = 5) -> tuple[np.ndarray, np.ndarray]:
    """Estimate positive target-after-source delay from local change correlations."""
    n = len(source)
    src = _standardized_change(source)
    dst = _standardized_change(target)
    delay = np.full(n, np.nan)
    score = np.full(n, np.nan)
    half = window // 2
    for center in range(half + max_lag, n - half - max_lag, stride):
        left = center - half
        right = center + half + 1
        best_lag, best_score = 0, -np.inf
        a = src[left:right]
        a = (a - a.mean()) / (a.std() + 1e-12)
        for lag in range(max_lag + 1):
            b = dst[left + lag:right + lag]
            b = (b - b.mean()) / (b.std() + 1e-12)
            corr = float(np.mean(a * b))
            if corr > best_score:
                best_lag, best_score = lag, corr
        delay[center] = best_lag * time_step
        score[center] = best_score
    valid = np.isfinite(delay)
    if valid.any():
        delay = np.interp(np.arange(n), np.flatnonzero(valid), delay[valid])
        score = np.interp(np.arange(n), np.flatnonzero(valid), score[valid])
    return delay, score


def fit_flow_delay(flow: np.ndarray, delay: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    valid = np.isfinite(flow) & np.isfinite(delay) & (flow > 1e-6)
    if valid.sum() < 3:
        return np.full(2, np.nan), {"n": int(valid.sum()), "rmse_s": float("nan")}
    x = 1.0 / np.sqrt(flow[valid])
    design = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(design, delay[valid], rcond=None)
    residual = delay[valid] - design @ coef
    return coef, {"n": int(valid.sum()), "rmse_s": float(np.sqrt(np.mean(residual**2))), "a_s": float(coef[0]), "b_s_sqrt_m3": float(coef[1])}


def process_one(path: Path, output_dir: Path) -> dict:
    raw = read_process_file(path)
    time = raw["time"]
    masks = {name: hampel_mask(raw[name]) for name in ("oxygen_cumulative", "co2", "co")}
    clean = {name: pchip_clean(time, raw[name], masks[name]) for name in masks}
    flow = oxygen_flow(time, clean["oxygen_cumulative"])
    # The derivative is the flow estimate. Its trapezoidal integral is retained for audit.
    oxygen_from_flow = clean["oxygen_cumulative"][0] + np.concatenate([[0.0], cumulative_trapezoid(flow, time)])
    total_co2 = np.concatenate([[0.0], cumulative_trapezoid(clean["co2"], time)])
    total_co = np.concatenate([[0.0], cumulative_trapezoid(clean["co"], time)])
    delay_co, score_co = local_delay(flow, clean["co"], 0.5)
    delay_co2, score_co2 = local_delay(flow, clean["co2"], 0.5)
    coef_co, fit_co = fit_flow_delay(flow, delay_co)
    coef_co2, fit_co2 = fit_flow_delay(flow, delay_co2)
    metadata = {
        "time_step_s": 0.5,
        "oxygen_unit": "m3",
        "flow_unit": "m3/s",
        "flue_gas_unit": "m3/s",
        "emission_unit": "m3",
        "outlier_method": "Hampel local median, window=9, threshold=6",
        "reconstruction": "PCHIP on visible samples",
        "flow_method": "Savitzky-Golay derivative, window=31, polyorder=3, clipped at zero",
        "delay_method": "local correlation of detrended changes, window=121 samples, max lag=60 s",
        "delay_fit": {"CO": fit_co, "CO2": fit_co2},
        "delay_fit_coefficients": {"CO": coef_co.tolist(), "CO2": coef_co2.tolist()},
        "mask_counts": {name: int(mask.sum()) for name, mask in masks.items()},
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_dir / f"{raw['heat_id']}_q2_processed.npz", heat_id=raw["heat_id"], sequence=raw["sequence"], time=time, oxygen_cumulative_raw=raw["oxygen_cumulative"], oxygen_cumulative_clean=clean["oxygen_cumulative"], oxygen_flow=flow, oxygen_from_flow=oxygen_from_flow, co_raw=raw["co"], co_clean=clean["co"], co2_raw=raw["co2"], co2_clean=clean["co2"], total_co_emission=total_co, total_co2_emission=total_co2, outlier_oxygen=masks["oxygen_cumulative"], outlier_co=masks["co"], outlier_co2=masks["co2"], delay_co=delay_co, delay_co2=delay_co2, delay_score_co=score_co, delay_score_co2=score_co2, metadata=json.dumps(metadata, ensure_ascii=False))
    return {"heat_id": raw["heat_id"], "records": len(time), **metadata["mask_counts"], "flow_mean": float(np.mean(flow)), "flow_max": float(np.max(flow)), "delay_co_median_s": float(np.median(delay_co)), "delay_co2_median_s": float(np.median(delay_co2)), "delay_co_fit_rmse_s": fit_co["rmse_s"], "delay_co2_fit_rmse_s": fit_co2["rmse_s"]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("Data"))
    parser.add_argument("--output-dir", type=Path, default=Path("problem2/final/data"))
    parser.add_argument("--summary", type=Path, default=Path("problem2/final/results/q2_summary.json"))
    args = parser.parse_args()
    rows = [process_one(path, args.output_dir) for path in sorted(args.data_dir.glob("*.xlsb*"))]
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    for row in rows:
        print(json.dumps(row, ensure_ascii=False))


if __name__ == "__main__":
    main()
