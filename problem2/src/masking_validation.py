"""Multi-seed random and contiguous-block masking evaluation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator


def _rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def evaluate(path: Path, seeds: tuple[int, ...] = (11, 23, 37), block_lengths: tuple[int, ...] = (2, 4, 8, 16)) -> list[dict]:
    with np.load(path, allow_pickle=False) as data:
        time = data["time"]
        signals = {name: data[f"{name}_clean"] for name in ("oxygen_cumulative", "co", "co2")}
    rows = []
    for name, values in signals.items():
        for seed in seeds:
            rng = np.random.default_rng(seed)
            random_mask = rng.random(len(values)) < 0.05
            random_mask[[0, -1]] = False
            pred = PchipInterpolator(time[~random_mask], values[~random_mask])(time[random_mask])
            rows.append({"heat_id": path.stem.split("_q2")[0], "signal": name, "scenario": "random_5pct", "seed": seed, "masked": int(random_mask.sum()), "rmse": _rmse(pred, values[random_mask])})
            for length in block_lengths:
                start = int(rng.integers(2, max(3, len(values) - length - 2)))
                mask = np.zeros(len(values), dtype=bool)
                mask[start : start + length] = True
                pred = PchipInterpolator(time[~mask], values[~mask])(time[mask])
                rows.append({"heat_id": path.stem.split("_q2")[0], "signal": name, "scenario": f"block_{length}", "seed": seed, "masked": length, "rmse": _rmse(pred, values[mask])})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("problem2/final/data"))
    parser.add_argument("--output", type=Path, default=Path("problem2/final/results/masking_validation.json"))
    args = parser.parse_args()
    rows = [row for path in sorted(args.data_dir.glob("*_q2_processed.npz")) for row in evaluate(path)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    print({"experiments": len(rows), "mean_rmse": float(np.mean([row["rmse"] for row in rows]))})


if __name__ == "__main__":
    main()
