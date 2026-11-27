"""Window and lag-range sensitivity on representative heats."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .delay_v2 import local_delay_reliable


def evaluate(path: Path, windows: tuple[int, ...], max_lags: tuple[int, ...]) -> list[dict]:
    with np.load(path, allow_pickle=False) as data:
        flow = data["oxygen_flow"]
        signals = {"co": data["co_clean"], "co2": data["co2_clean"]}
    rows = []
    for window in windows:
        for max_lag in max_lags:
            for gas, signal in signals.items():
                info = local_delay_reliable(flow, signal, 0.5, window=window, max_lag=max_lag)
                rows.append({
                    "heat_id": path.stem.split("_q2")[0],
                    "gas": gas,
                    "window": window,
                    "max_lag_s": max_lag * 0.5,
                    "candidate_windows": int(len(info["center_index"])),
                    "reliable_windows": int(info["reliable"].sum()),
                    "coverage": float(info["reliable"].mean()) if len(info["reliable"]) else 0.0,
                    "median_delay_s": float(np.nanmedian(info["delay_observed"])) if np.isfinite(info["delay_observed"]).any() else np.nan,
                })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("problem2/final/data"))
    parser.add_argument("--output", type=Path, default=Path("problem2/final/results/delay_parameter_sensitivity.json"))
    args = parser.parse_args()
    wanted = {"6B03205", "6C02057"}
    paths = [path for path in sorted(args.data_dir.glob("*_q2_processed.npz")) if path.stem.split("_q2")[0] in wanted]
    rows = [row for path in paths for row in evaluate(path, (81, 121, 161), (40, 60, 80))]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, ensure_ascii=False, indent=2, allow_nan=True), encoding="utf-8")
    print({"experiments": len(rows)})


if __name__ == "__main__":
    main()
