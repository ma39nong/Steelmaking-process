"""Generate figures from the saved Problem 2 NPZ interfaces."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files if key != "metadata"}


def make_figures(data_dir: Path, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = sorted(data_dir.glob("*_q2_processed.npz"))
    if not paths:
        raise FileNotFoundError(f"no q2 NPZ files in {data_dir}")
    fig, axes = plt.subplots(6, 1, figsize=(12, 15), sharex=False, constrained_layout=True)
    for path in paths:
        d = load(path)
        label = path.name.split("_q2")[0]
        axes[0].plot(d["time"], d["oxygen_cumulative_clean"], lw=0.8, label=label)
        axes[1].plot(d["time"], d["oxygen_flow"], lw=0.8, label=label)
        axes[2].plot(d["time"], d["co_clean"], lw=0.8, label=label)
        axes[3].plot(d["time"], d["co2_clean"], lw=0.8, label=label)
        axes[4].plot(d["time"], d["total_co_emission"], lw=0.8, label=label)
        axes[5].plot(d["time"], d["total_co2_emission"], lw=0.8, label=label)
    axes[0].set_ylabel("Cumulative O$_2$ / m$^3$")
    axes[1].set_ylabel("Lance flow / m$^3$ s$^{-1}$")
    axes[2].set_ylabel("CO / m$^3$ s$^{-1}$")
    axes[3].set_ylabel("CO$_2$ / m$^3$ s$^{-1}$")
    axes[4].set_ylabel("Total CO / m$^3$")
    axes[5].set_ylabel("Total CO$_2$ / m$^3$")
    axes[5].set_xlabel("Relative time / s")
    axes[0].legend(ncol=2, fontsize=8)
    fig.savefig(figure_dir / "ten_heats_six_process_curves.png", dpi=180)
    fig.savefig(figure_dir / "ten_heats_six_process_curves.svg")
    plt.close(fig)

    d = load(paths[0])
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True, constrained_layout=True)
    axes[0].plot(d["time"], d["oxygen_cumulative_raw"], color="0.65", lw=0.8, label="raw")
    axes[0].plot(d["time"], d["oxygen_cumulative_clean"], color="#1261A0", lw=1.2, label="PCHIP")
    axes[1].plot(d["time"], d["co2_clean"], color="#C23B22", lw=1.0, label="CO$_2$")
    axes[1].plot(d["time"], d["co_clean"], color="#2E7D32", lw=1.0, label="CO")
    co_centers = d["delay_co_center_index"]
    co2_centers = d["delay_co2_center_index"]
    axes[2].plot(d["time"][co_centers], d["delay_co_observed"], color="#2E7D32", lw=1.0, label="CO observed delay")
    axes[2].plot(d["time"][co2_centers], d["delay_co2_observed"], color="#C23B22", lw=1.0, label="CO$_2$ observed delay")
    axes[0].set_ylabel("O$_2$ / m$^3$")
    axes[1].set_ylabel("Flow / m$^3$ s$^{-1}$")
    axes[2].set_ylabel("Delay / s")
    axes[2].set_xlabel("Relative time / s")
    for axis in axes:
        axis.legend()
    fig.savefig(figure_dir / f"{paths[0].stem}_signals_and_delay.png", dpi=180)
    fig.savefig(figure_dir / f"{paths[0].stem}_signals_and_delay.svg")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("problem2/final/data"))
    parser.add_argument("--figure-dir", type=Path, default=Path("problem2/final/figures"))
    args = parser.parse_args()
    make_figures(args.data_dir, args.figure_dir)


if __name__ == "__main__":
    main()
