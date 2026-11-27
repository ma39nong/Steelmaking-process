"""V2 pipeline with reliability-aware delay observations and time mappings."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .delay_v2 import fit_delay_models, local_delay_reliable
from .q2_pipeline import (
    hampel_mask,
    oxygen_flow,
    pchip_clean,
    process_one,
    read_process_file,
)


def process_one_v2(path: Path, output_dir: Path, summary_dir: Path) -> dict:
    raw = read_process_file(path)
    time = raw["time"]
    masks = {name: hampel_mask(raw[name]) for name in ("oxygen_cumulative", "co2", "co")}
    clean = {name: pchip_clean(time, raw[name], masks[name]) for name in masks}
    flow = oxygen_flow(time, clean["oxygen_cumulative"])
    # V1 remains available through q2_pipeline.process_one. V2 uses observed
    # reliable centers only and never fills delay observations at boundaries.
    delay_info = {
        gas: local_delay_reliable(flow, clean[gas], 0.5)
        for gas in ("co", "co2")
    }
    models = {
        gas: fit_delay_models(flow[info["center_index"]], info["delay_observed"], info["confidence"], 60.0)
        for gas, info in delay_info.items()
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {
        "heat_id": np.asarray(raw["heat_id"]),
        "sequence": raw["sequence"],
        "time": time,
        "oxygen_cumulative_raw": raw["oxygen_cumulative"],
        "oxygen_cumulative_clean": clean["oxygen_cumulative"],
        "oxygen_flow": flow,
        "co_raw": raw["co"],
        "co_clean": clean["co"],
        "co2_raw": raw["co2"],
        "co2_clean": clean["co2"],
        "outlier_oxygen": masks["oxygen_cumulative"],
        "outlier_co": masks["co"],
        "outlier_co2": masks["co2"],
    }
    metadata = {
        "version": "V2-audit-fix",
        "delay_definition": "target measurement time minus source event time",
        "delay_parameters": {"window": 121, "max_lag_s": 60.0, "stride": 5, "min_corr": 0.35, "min_separation": 0.05},
        "models": models,
        "mask_counts": {name: int(mask.sum()) for name, mask in masks.items()},
    }
    for gas, info in delay_info.items():
        prefix = f"delay_{gas}"
        for key, value in info.items():
            if key == "center_time":
                arrays[f"{prefix}_window_time"] = value
            elif key == "delay_observed":
                arrays[f"{prefix}_observed"] = value
            elif key == "reliable":
                arrays[f"{prefix}_valid_mask"] = value
            elif key == "confidence":
                arrays[f"{prefix}_confidence"] = value
            elif key == "center_index":
                arrays[f"{prefix}_center_index"] = value
            else:
                arrays[f"{prefix}_{key}"] = value
        valid_centers = info["reliable"]
        centers = info["center_index"]
        observed = np.full(len(time), np.nan)
        confidence = np.zeros(len(time))
        observed_mask = np.zeros(len(time), dtype=bool)
        observed[centers[valid_centers]] = info["delay_observed"][valid_centers]
        confidence[centers[valid_centers]] = info["confidence"][valid_centers]
        observed_mask[centers[valid_centers]] = True
        model = models[gas]
        if model.get("selected") == "constrained_inverse_sqrt":
            a = model["constrained_inverse_sqrt"]["a_s"]
            b = model["constrained_inverse_sqrt"]["b_s_sqrt_m3"]
            delay_model = np.clip(a + b / np.sqrt(np.maximum(flow, 1e-9)), 0, 60.0)
        elif model.get("selected") == "fixed_median":
            delay_model = np.full(len(time), model["fixed"]["c_s"])
        else:
            delay_model = np.full(len(time), np.nan)
        mapping_model = time + delay_model
        mapping_observed = time + observed
        mapping_valid = observed_mask & np.isfinite(mapping_observed) & (mapping_observed >= time[0]) & (mapping_observed <= time[-1])
        if mapping_valid.sum() >= 2:
            ordered = np.diff(mapping_observed[mapping_valid]) > 0
            monotonic = np.r_[True, ordered]
            valid_positions = np.flatnonzero(mapping_valid)
            mapping_valid[valid_positions[~monotonic]] = False
        arrays[f"delay_{gas}_model"] = delay_model
        arrays[f"delay_{gas}_observed_mask"] = observed_mask
        arrays[f"time_mapping_{gas}"] = mapping_observed
        arrays[f"time_mapping_{gas}_model"] = mapping_model
        arrays[f"alignment_valid_mask_{gas}"] = mapping_valid
        aligned = np.full(len(time), np.nan)
        aligned[mapping_valid] = np.interp(mapping_observed[mapping_valid], time, clean[gas])
        arrays[f"{gas}_aligned"] = aligned
        metadata[f"{gas}_coverage"] = {
            "window_count": int(len(centers)),
            "reliable_window_count": int(valid_centers.sum()),
            "reliable_fraction": float(valid_centers.mean()) if len(valid_centers) else 0.0,
            "observed_point_count": int(observed_mask.sum()),
            "alignment_valid_count": int(mapping_valid.sum()),
            "boundary_unavailable_count": int(len(time) - (centers[-1] - centers[0] + 1)) if len(centers) else int(len(time)),
            "between_window_gap_count": int((centers[-1] - centers[0] + 1) - len(centers)) if len(centers) else 0,
        }
    arrays["metadata"] = np.asarray(json.dumps(metadata, ensure_ascii=False))
    np.savez_compressed(output_dir / f"{raw['heat_id']}_q2_v2_processed.npz", **arrays)
    result = {"heat_id": raw["heat_id"], "records": len(time), **metadata["mask_counts"], "models": models, "coverage": {gas: metadata[f"{gas}_coverage"] for gas in ("co", "co2")}}
    (summary_dir / f"{raw['heat_id']}_q2_v2_summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("Data"))
    parser.add_argument("--output-dir", type=Path, default=Path("problem2/final/v2_data"))
    parser.add_argument("--summary-dir", type=Path, default=Path("problem2/final/v2_results"))
    args = parser.parse_args()
    rows = [process_one_v2(path, args.output_dir, args.summary_dir) for path in sorted(args.data_dir.glob("*.xlsb*"))]
    (args.summary_dir / "q2_v2_summary.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    for row in rows:
        print(json.dumps(row, ensure_ascii=False))


if __name__ == "__main__":
    main()
