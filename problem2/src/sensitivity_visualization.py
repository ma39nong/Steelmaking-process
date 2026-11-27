"""Plot reliability-threshold sensitivity from saved diagnostics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("problem2/final/results/cross_validation.json"))
    parser.add_argument("--output", type=Path, default=Path("problem2/final/figures/delay_threshold_sensitivity.png"))
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    for gas, color in (("co", "#2E7D32"), ("co2", "#C23B22")):
        thresholds = sorted({row["min_corr"] for row in data["threshold_sensitivity"] if row["gas"] == gas})
        coverage = []
        for threshold in thresholds:
            rows = [row for row in data["threshold_sensitivity"] if row["gas"] == gas and row["min_corr"] == threshold]
            coverage.append(sum(row["accepted_windows"] for row in rows) / sum(row["candidate_windows"] for row in rows))
        ax.plot(thresholds, coverage, marker="o", color=color, label=gas.upper())
    ax.set_xlabel("Minimum correlation threshold")
    ax.set_ylabel("Accepted-window fraction")
    ax.grid(alpha=0.25)
    ax.legend()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
