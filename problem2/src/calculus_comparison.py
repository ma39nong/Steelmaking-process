"""Compare SG and PCHIP analytical derivatives and emission integrals."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.signal import savgol_filter


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def evaluate(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        time = data["time"]
        oxygen = data["oxygen_cumulative_clean"]
        co = data["co_clean"]
        co2 = data["co2_clean"]
        clipped_flow = data["oxygen_flow"]
    oxygen_fn = PchipInterpolator(time, oxygen)
    co_fn = PchipInterpolator(time, co)
    co2_fn = PchipInterpolator(time, co2)
    pchip_flow = oxygen_fn.derivative()(time)
    flow_antiderivative = oxygen_fn.derivative().antiderivative()
    sg_flow = savgol_filter(oxygen, 31, 3, deriv=1, delta=0.5, mode="interp")
    pchip_co_total = co_fn.antiderivative()(time) - co_fn.antiderivative()(time[0])
    pchip_co2_total = co2_fn.antiderivative()(time) - co2_fn.antiderivative()(time[0])
    trap_co_total = np.concatenate([[0.0], cumulative_trapezoid(co, time)])
    trap_co2_total = np.concatenate([[0.0], cumulative_trapezoid(co2, time)])
    return {
        "heat_id": path.stem.split("_q2")[0],
        "sg_vs_pchip_flow_rmse": rmse(sg_flow, pchip_flow),
        "clipped_flow_vs_pchip_derivative_rmse": rmse(clipped_flow, pchip_flow),
        "sg_negative_fraction": float(np.mean(sg_flow < 0)),
        "pchip_derivative_negative_fraction": float(np.mean(pchip_flow < 0)),
        "co_trapezoid_vs_pchip_integral_rmse": rmse(trap_co_total, pchip_co_total),
        "co2_trapezoid_vs_pchip_integral_rmse": rmse(trap_co2_total, pchip_co2_total),
        "pchip_exact_integral_error": float(abs((flow_antiderivative(time[-1]) - flow_antiderivative(time[0])) - (oxygen[-1] - oxygen[0]))),
        "node_sampled_pchip_trapezoid_error": float(abs(np.trapz(pchip_flow, time) - (oxygen[-1] - oxygen[0]))),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("problem2/final/data"))
    parser.add_argument("--output", type=Path, default=Path("problem2/final/results/calculus_comparison.json"))
    args = parser.parse_args()
    rows = [evaluate(path) for path in sorted(args.data_dir.glob("*_q2_processed.npz"))]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
