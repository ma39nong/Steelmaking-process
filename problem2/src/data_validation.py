"""Validate the ten raw XLSB process files without modifying them."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")


ALIASES = {
    "oxygen_cumulative": ["cumulative oxygen consumption", "累积耗氧量", "累计耗氧量", "yqx h", "yqxh"],
    "co2": ["co₂ concentration", "co2 concentration", "co₂浓度", "co2浓度"],
    "co": ["co concentration", "co浓度"],
}


def _norm(value: object) -> str:
    text = str(value).strip().lower().replace(" ", "")
    return text.replace("（", "(").replace("）", ")")


def _find_column(columns: list[object], aliases: list[str]) -> object:
    normalized = {_norm(c): c for c in columns}
    for alias in aliases:
        key = _norm(alias)
        if key in normalized:
            return normalized[key]
    for key, original in normalized.items():
        if any(_norm(alias) in key for alias in aliases):
            return original
    raise KeyError(f"cannot identify process column from aliases {aliases}")


def heat_id_from_file(path: Path) -> str:
    match = re.search(r"Heat No\.\s*([A-Za-z0-9]+)", path.name, flags=re.I)
    if not match:
        raise ValueError(f"cannot infer heat id from {path.name}")
    return match.group(1)


def validate_file(path: Path) -> dict:
    heat_id = heat_id_from_file(path)
    excel = pd.ExcelFile(path, engine="pyxlsb")
    if len(excel.sheet_names) != 1:
        raise ValueError(f"{heat_id}: expected one sheet, found {excel.sheet_names}")
    sheet = excel.sheet_names[0]
    frame = pd.read_excel(path, sheet_name=sheet, engine="pyxlsb")
    if frame.shape[1] != 2051:
        raise ValueError(f"{heat_id}: expected 2051 columns, found {frame.shape[1]}")
    sequence_col = frame.columns[0]
    selected = {
        name: _find_column(list(frame.columns), aliases)
        for name, aliases in ALIASES.items()
    }
    sequence = pd.to_numeric(frame[sequence_col], errors="coerce").to_numpy(dtype=float)
    seq_finite = sequence[np.isfinite(sequence)]
    if len(seq_finite) != len(frame) or not np.allclose(seq_finite, np.arange(1, len(frame) + 1)):
        raise ValueError(f"{heat_id}: Sequence is not the continuous 1..n index")
    values = {}
    for name, column in selected.items():
        values[name] = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
    time = (seq_finite - 1.0) * 0.5
    summary = {
        "heat_id": heat_id,
        "sheet_name": sheet,
        "records": int(len(frame)),
        "columns": int(frame.shape[1]),
        "sequence_column": str(sequence_col),
        "process_columns": {k: str(v) for k, v in selected.items()},
        "sequence_first": int(seq_finite[0]),
        "sequence_last": int(seq_finite[-1]),
        "duplicate_sequence_count": int(len(seq_finite) - len(np.unique(seq_finite))),
        "gap_count": int(np.sum(np.diff(seq_finite) != 1)),
        "time_start_s": float(time[0]),
        "time_end_s": float(time[-1]),
        "time_step_min_s": float(np.min(np.diff(time))) if len(time) > 1 else 0.0,
        "time_step_max_s": float(np.max(np.diff(time))) if len(time) > 1 else 0.0,
        "missing": {name: int(np.sum(~np.isfinite(array))) for name, array in values.items()},
        "min": {name: float(np.nanmin(array)) for name, array in values.items()},
        "max": {name: float(np.nanmax(array)) for name, array in values.items()},
        "negative_count": {name: int(np.sum(array < 0)) for name, array in values.items()},
    }
    return summary


def validate_directory(data_dir: Path, output_dir: Path) -> list[dict]:
    files = sorted(data_dir.glob("*.xlsb*"))
    if len(files) != 10:
        raise ValueError(f"expected ten XLSB files in {data_dir}, found {len(files)}")
    summaries = [validate_file(path) for path in files]
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "process_validation.json").write_text(
        json.dumps(summaries, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    flat = []
    for item in summaries:
        row = {k: v for k, v in item.items() if k not in {"process_columns", "missing", "min", "max", "negative_count"}}
        row["process_columns"] = json.dumps(item["process_columns"], ensure_ascii=False)
        for group in ("missing", "min", "max", "negative_count"):
            for key, value in item[group].items():
                row[f"{group}_{key}"] = value
        flat.append(row)
    pd.DataFrame(flat).to_csv(output_dir / "process_validation.csv", index=False, encoding="utf-8-sig")
    return summaries


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("Data"))
    parser.add_argument("--output-dir", type=Path, default=Path("problem2/final/results"))
    args = parser.parse_args()
    summaries = validate_directory(args.data_dir, args.output_dir)
    for item in summaries:
        print(item["heat_id"], item["records"], item["process_columns"], item["time_end_s"])


if __name__ == "__main__":
    main()
