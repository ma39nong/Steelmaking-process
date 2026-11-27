"""Deterministic checks for the saved Problem 2 interfaces."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def evaluate(path: Path, fraction: float = 0.05, seed: int = 20261127) -> dict:
    with np.load(path, allow_pickle=False) as data:
        numeric_fields = {}
        shape_fields = {}
        for key in data.files:
            if key == "metadata":
                continue
            array = data[key]
            shape_fields[key] = list(array.shape)
            if np.issubdtype(array.dtype, np.number):
                numeric_fields[key] = bool(np.isfinite(array).all())
        time = data["time"]
        oxygen = data["oxygen_cumulative_clean"]
        co = data["co_clean"]
        co2 = data["co2_clean"]
        flow = data["oxygen_flow"]
        integrated = data["oxygen_from_flow"]
    rng = np.random.default_rng(seed)
    mask = rng.random(len(time)) < fraction
    mask[[0, -1]] = False
    visible = ~mask
    reconstructed = {}
    for name, values in (("oxygen", oxygen), ("co", co), ("co2", co2)):
        reconstructed[name] = PchipInterpolator(time[visible], values[visible])(time[mask])
    return {
        "heat_id": path.stem.split("_q2")[0],
        "masked_points": int(mask.sum()),
        "mask_fraction": float(mask.mean()),
        "oxygen_masked_rmse": rmse(reconstructed["oxygen"], oxygen[mask]),
        "co_masked_rmse": rmse(reconstructed["co"], co[mask]),
        "co2_masked_rmse": rmse(reconstructed["co2"], co2[mask]),
        "oxygen_flow_min": float(np.min(flow)),
        "oxygen_reintegration_rmse": rmse(integrated, oxygen),
        "all_numeric_fields_finite": bool(all(numeric_fields.values())),
        "numeric_field_finite": numeric_fields,
        "field_shapes": shape_fields,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("problem2/final/data"))
    parser.add_argument("--output", type=Path, default=Path("problem2/final/results/q2_evaluation.json"))
    args = parser.parse_args()
    rows = [evaluate(path) for path in sorted(args.data_dir.glob("*_q2_processed.npz"))]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
