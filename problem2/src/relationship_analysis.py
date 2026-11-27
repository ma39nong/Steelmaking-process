"""Analyze the four cumulative-oxygen / flue-gas relationships required by the statement."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.stats import pearsonr, spearmanr


def _fit_metrics(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return {"n": int(len(x)), "pearson_r": np.nan, "spearman_r": np.nan, "linear_r2": np.nan, "quadratic_r2": np.nan}
    linear = np.polyfit(x, y, 1)
    quadratic = np.polyfit(x, y, 2)
    pred_linear = np.polyval(linear, x)
    pred_quadratic = np.polyval(quadratic, x)
    total = np.sum((y - np.mean(y)) ** 2)
    return {
        "n": int(len(x)),
        "pearson_r": float(pearsonr(x, y).statistic),
        "spearman_r": float(spearmanr(x, y).statistic),
        "linear_r2": float(1 - np.sum((y - pred_linear) ** 2) / total),
        "quadratic_r2": float(1 - np.sum((y - pred_quadratic) ** 2) / total),
    }


def analyze(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        time = data["time"]
        oxygen = data["oxygen_cumulative_clean"]
        flow = data["oxygen_flow"]
        co = data["co_clean"]
        co2 = data["co2_clean"]
    emissions = {
        "co": np.concatenate([[0.0], cumulative_trapezoid(co, time)]),
        "co2": np.concatenate([[0.0], cumulative_trapezoid(co2, time)]),
    }
    relationships = {
        "oxygen_vs_co_instantaneous": _fit_metrics(oxygen, co),
        "oxygen_vs_co2_instantaneous": _fit_metrics(oxygen, co2),
        "oxygen_vs_co_cumulative_emission": _fit_metrics(oxygen, emissions["co"]),
        "oxygen_vs_co2_cumulative_emission": _fit_metrics(oxygen, emissions["co2"]),
    }
    bins = np.quantile(oxygen, np.linspace(0, 1, 6))
    stages = []
    for i in range(5):
        mask = (oxygen >= bins[i]) & ((oxygen <= bins[i + 1]) if i == 4 else (oxygen < bins[i + 1]))
        stages.append({
            "stage": i + 1,
            "oxygen_min_m3": float(bins[i]),
            "oxygen_max_m3": float(bins[i + 1]),
            "n": int(mask.sum()),
            "co_m3_s_median": float(np.median(co[mask])) if mask.any() else np.nan,
            "co2_m3_s_median": float(np.median(co2[mask])) if mask.any() else np.nan,
            "co_m3_median": float(np.median(emissions["co"][mask])) if mask.any() else np.nan,
            "co2_m3_median": float(np.median(emissions["co2"][mask])) if mask.any() else np.nan,
        })
    ratio_mask = flow > max(1e-6, np.quantile(flow, 0.05))
    ratios = {
        "co_per_oxygen_m3_per_m3": _fit_metrics(oxygen[ratio_mask], (co[ratio_mask] / flow[ratio_mask])),
        "co2_per_oxygen_m3_per_m3": _fit_metrics(oxygen[ratio_mask], (co2[ratio_mask] / flow[ratio_mask])),
    }
    return {"heat_id": path.stem.split("_q2")[0], "relationships": relationships, "oxygen_stage_statistics": stages, "unit_oxygen_rate_diagnostics": ratios}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("problem2/final/data"))
    parser.add_argument("--output", type=Path, default=Path("problem2/final/results/relationship_summary.json"))
    args = parser.parse_args()
    rows = [analyze(path) for path in sorted(args.data_dir.glob("*_q2_processed.npz"))]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, ensure_ascii=False, indent=2, allow_nan=True), encoding="utf-8")
    for row in rows:
        print(row["heat_id"], row["relationships"])


if __name__ == "__main__":
    main()
