"""Cross-heat validation and threshold sensitivity for delay observations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .delay_v2 import fit_delay_models


def _load(path: Path, gas: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        centers = data[f"delay_{gas}_center_index"]
        flow = data["oxygen_flow"][centers]
        delay = data[f"delay_{gas}_observed"]
        confidence = data[f"delay_{gas}_confidence"]
    return flow, delay, confidence


def _predict(model: dict, flow: np.ndarray) -> np.ndarray:
    if model.get("selected") == "constrained_inverse_sqrt":
        p = model["constrained_inverse_sqrt"]
        return np.clip(p["a_s"] + p["b_s_sqrt_m3"] / np.sqrt(np.maximum(flow, 1e-9)), 0, 60)
    if model.get("selected") == "fixed_median":
        return np.full_like(flow, model["fixed"]["c_s"], dtype=float)
    return np.full_like(flow, np.nan, dtype=float)


def leave_one_heat_out(data_dir: Path) -> list[dict]:
    paths = sorted(data_dir.glob("*_q2_processed.npz"))
    rows = []
    for gas in ("co", "co2"):
        samples = {p.stem.split("_q2")[0]: _load(p, gas) for p in paths}
        for test_heat, (test_flow, test_delay, test_conf) in samples.items():
            train_flow = np.concatenate([v[0][v[2] > 0] for heat, v in samples.items() if heat != test_heat])
            train_delay = np.concatenate([v[1][v[2] > 0] for heat, v in samples.items() if heat != test_heat])
            train_conf = np.concatenate([v[2][v[2] > 0] for heat, v in samples.items() if heat != test_heat])
            model = fit_delay_models(train_flow, train_delay, train_conf, 60.0)
            valid = (test_conf > 0) & np.isfinite(test_delay)
            prediction = _predict(model, test_flow[valid]) if valid.any() else np.array([])
            mae = float(np.mean(np.abs(prediction - test_delay[valid]))) if prediction.size and np.isfinite(prediction).all() else np.nan
            rows.append({"gas": gas, "test_heat": test_heat, "train_reliable_samples": int(len(train_delay)), "test_reliable_samples": int(valid.sum()), "selected": model.get("selected"), "test_mae_s": mae})
    return rows


def threshold_sensitivity(data_dir: Path) -> list[dict]:
    rows = []
    for path in sorted(data_dir.glob("*_q2_processed.npz")):
        heat = path.stem.split("_q2")[0]
        with np.load(path, allow_pickle=False) as data:
            for gas in ("co", "co2"):
                corr = data[f"delay_{gas}_max_corr"]
                separation = data[f"delay_{gas}_peak_separation"]
                boundary = data[f"delay_{gas}_boundary_hit"]
                for threshold in (0.25, 0.30, 0.35, 0.40):
                    valid = np.isfinite(corr) & np.isfinite(separation) & (corr >= threshold) & (separation >= 0.05) & ~boundary
                    rows.append({"heat_id": heat, "gas": gas, "min_corr": threshold, "candidate_windows": int(len(corr)), "accepted_windows": int(valid.sum()), "coverage": float(valid.mean())})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("problem2/final/data"))
    parser.add_argument("--output", type=Path, default=Path("problem2/final/results/cross_validation.json"))
    args = parser.parse_args()
    result = {"leave_one_heat_out": leave_one_heat_out(args.data_dir), "threshold_sensitivity": threshold_sensitivity(args.data_dir)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=True), encoding="utf-8")
    print(json.dumps({"loho_rows": len(result["leave_one_heat_out"]), "sensitivity_rows": len(result["threshold_sensitivity"])}, ensure_ascii=False))


if __name__ == "__main__":
    main()
