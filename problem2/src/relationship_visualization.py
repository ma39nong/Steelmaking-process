"""Plot the four required cumulative-oxygen relationships."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid


def make_figure(path: Path, output_dir: Path) -> None:
    with np.load(path, allow_pickle=False) as data:
        oxygen = data["oxygen_cumulative_clean"]
        time = data["time"]
        co = data["co_clean"]
        co2 = data["co2_clean"]
    co_total = np.concatenate([[0.0], cumulative_trapezoid(co, time)])
    co2_total = np.concatenate([[0.0], cumulative_trapezoid(co2, time)])
    targets = [(co, "CO / m$^3$ s$^{-1}$"), (co2, "CO$_2$ / m$^3$ s$^{-1}$"), (co_total, "Total CO / m$^3$"), (co2_total, "Total CO$_2$ / m$^3$")]
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    for axis, (target, label) in zip(axes.flat, targets):
        axis.scatter(oxygen, target, s=3, alpha=0.35, rasterized=True)
        axis.set_xlabel("Cumulative oxygen / m$^3$")
        axis.set_ylabel(label)
    heat_id = path.stem.split("_q2")[0]
    fig.suptitle(f"{heat_id}: cumulative oxygen relationships")
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{heat_id}_oxygen_relationships.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("problem2/final/data"))
    parser.add_argument("--output-dir", type=Path, default=Path("problem2/final/figures"))
    args = parser.parse_args()
    for path in sorted(args.data_dir.glob("*_q2_processed.npz")):
        make_figure(path, args.output_dir)


if __name__ == "__main__":
    main()
