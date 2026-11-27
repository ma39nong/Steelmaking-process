"""Run the complete Problem 2 offline workflow in documented order."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run(module: str, args: list[str]) -> None:
    command = [sys.executable, "-m", module, *args]
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("Data"))
    parser.add_argument("--output-root", type=Path, default=Path("problem2/final"))
    args = parser.parse_args()
    results = args.output_root / "results"
    run("problem2.src.data_validation", ["--data-dir", str(args.data_dir), "--output-dir", str(results)])
    run("problem2.src.q2_pipeline_v2", ["--data-dir", str(args.data_dir), "--output-dir", str(args.output_root / "data"), "--summary-dir", str(results)])
    run("problem2.src.relationship_analysis", ["--data-dir", str(args.output_root / "data"), "--output", str(results / "relationship_summary.json")])
    run("problem2.src.v2_validation", ["--data-dir", str(args.output_root / "data"), "--output", str(results / "cross_validation.json")])
    run("problem2.src.calculus_comparison", ["--data-dir", str(args.output_root / "data"), "--output", str(results / "calculus_comparison.json")])
    run("problem2.src.evaluation", ["--data-dir", str(args.output_root / "data"), "--output", str(results / "reconstruction_evaluation.json")])
    run("problem2.src.injection_validation", ["--raw-data-dir", str(args.data_dir), "--output", str(results / "raw_injection_validation.json")])
    run("problem2.src.masking_validation", ["--data-dir", str(args.output_root / "data"), "--output", str(results / "masking_validation.json")])
    run("problem2.src.delay_parameter_sensitivity", ["--data-dir", str(args.output_root / "data"), "--output", str(results / "delay_parameter_sensitivity.json")])
    run("problem2.src.visualization", ["--data-dir", str(args.output_root / "data"), "--figure-dir", str(args.output_root / "figures")])
    run("problem2.src.relationship_visualization", ["--data-dir", str(args.output_root / "data"), "--output-dir", str(args.output_root / "figures")])
    run("problem2.src.sensitivity_visualization", ["--input", str(results / "cross_validation.json"), "--output", str(args.output_root / "figures" / "delay_threshold_sensitivity.png")])


if __name__ == "__main__":
    main()
