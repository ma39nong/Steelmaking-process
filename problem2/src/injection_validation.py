"""Synthetic outlier injection checks using the real cleaned heat signals."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .q2_pipeline import hampel_mask, pchip_clean


def evaluate(path: Path, seeds: tuple[int, ...] = (11, 23, 37), count: int = 20) -> dict:
    with np.load(path, allow_pickle=False) as data:
        time = data["time"]
        signals = {name: data[f"{name}_clean"] for name in ("oxygen_cumulative", "co", "co2")}
    rows = []
    for name, reference in signals.items():
        for seed in seeds:
            rng = np.random.default_rng(seed)
            candidates = np.arange(4, len(reference) - 4)
            indices = np.sort(rng.choice(candidates, size=min(count, len(candidates)), replace=False))
            contaminated = reference.copy()
            amplitude = max(np.std(reference) * 8, 1.0)
            contaminated[indices] += amplitude
            mask = hampel_mask(contaminated)
            reconstructed = pchip_clean(time, contaminated, mask)
            rows.append({
                "heat_id": path.stem.split("_q2")[0],
                "signal": name,
                "seed": seed,
                "injected": int(len(indices)),
                "detected": int(mask[indices].sum()),
                "recall": float(mask[indices].mean()),
                "injected_rmse": float(np.sqrt(np.mean((reconstructed[indices] - reference[indices]) ** 2))),
            })
    return {"heat_id": path.stem.split("_q2")[0], "rows": rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("problem2/final/data"))
    parser.add_argument("--output", type=Path, default=Path("problem2/final/results/injection_validation.json"))
    args = parser.parse_args()
    result = [evaluate(path) for path in sorted(args.data_dir.glob("*_q2_processed.npz"))]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    flat = [row for heat in result for row in heat["rows"]]
    print({"experiments": len(flat), "mean_recall": float(np.mean([row["recall"] for row in flat]))})


if __name__ == "__main__":
    main()
