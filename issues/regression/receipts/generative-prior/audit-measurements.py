"""Audit the shared-trial measurement contrast with Python's standard library."""
import argparse
import hashlib
import itertools
import json
import math
import pathlib
import runpy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=pathlib.Path)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    args = parser.parse_args()
    base = pathlib.Path(__file__).parent
    common = runpy.run_path(str(base / "audit-linear.py"))
    read = common["read_rows"]
    controls = read(args.results / "data-checks.jsonl")
    assert [r["Repetition"] for r in controls] == list(range(1, 31))
    settings_bytes = (base / "settings/g1-measurements.json").read_bytes()
    config = json.loads(settings_bytes)
    settings_sha = hashlib.sha256(settings_bytes).hexdigest()
    counts = config["MeasurementCounts"]
    rows_by_m, cells = {}, []
    for m in counts:
        rows = read(args.results / f"m{m}/trials.jsonl")
        summary = json.loads((args.results / f"m{m}/summary.json").read_text())
        assert len(rows) == 30 and summary["GatePassed"]
        if m != 150:
            assert summary["SettingsSha256"] == settings_sha
            assert summary["Measurements"] == m
        for i, (row, check) in enumerate(zip(rows, controls), 1):
            assert row["Repetition"] == check["Repetition"] == i
            assert row["Seed"] == check["Seed"] and row["ChildSeeds"] == check["ChildSeeds"]
            assert check["BaselineScalarAndVectorBitsAgree"]
            assert check["MaxReferenceDiscrepancy"] <= 1e-10
            assert [r["Measurements"] for r in check["Matrices"]] == counts
            assert all(r["EntryBitMismatches"] == 0 for r in check["Matrices"])
            common["finite"](row)
            common["validate"](row, True)
        successes = sum(r["Success"] for r in rows)
        assert summary["Successes"] == successes
        rows_by_m[m] = rows
        cells.append({
            "m": m, "m_over_k": m / 80, "successes": successes, "n": 30,
            "wilson95": common["wilson"](successes),
            "classes": {c: sum(r["EndpointClass"] == c for r in rows)
                        for c in ("recovered", "negative-neighbourhood", "elsewhere")},
            "gradient_stops": sum(r["Stop"] == "gradient" for r in rows),
            "update_limit_stops": sum(r["Stop"] == "update-limit" for r in rows),
            "evaluated_updates": sum(r["EvaluatedUpdates"] for r in rows),
            "skipped_updates": sum(r["SkippedUpdates"] for r in rows),
        })
    baseline = read(base / "results/g1-linear-curve/k80-on/trials.jsonl")
    common["compare"](baseline, rows_by_m[150])
    assert baseline == rows_by_m[150]
    for m, rows in rows_by_m.items():
        for old, row in zip(baseline, rows):
            for key in ("Seed", "ChildSeeds", "Truth", "Start"):
                assert old[key] == row[key], (m, row["Repetition"], key)
    pairs, reversals = [], []
    for lower, higher in itertools.combinations(counts, 2):
        a, b = rows_by_m[lower], rows_by_m[higher]
        both = [i for i, (x, y) in enumerate(zip(a, b), 1) if x["Success"] and y["Success"]]
        lost = [i for i, (x, y) in enumerate(zip(a, b), 1) if x["Success"] and not y["Success"]]
        gained = [i for i, (x, y) in enumerate(zip(a, b), 1) if not x["Success"] and y["Success"]]
        neither = [i for i, (x, y) in enumerate(zip(a, b), 1) if not x["Success"] and not y["Success"]]
        n = len(lost) + len(gained)
        p = min(1, 2 * sum(math.comb(n, i) for i in range(min(len(lost), len(gained)) + 1)) / 2**n)
        pairs.append({
            "lower_m": lower, "higher_m": higher, "both": len(both), "lower_only": len(lost),
            "higher_only": len(gained), "neither": len(neither),
            "lower_only_repetitions": lost, "higher_only_repetitions": gained,
            "difference": (len(gained) - len(lost)) / 30, "mcnemar_p": p,
            "primary": (lower, higher) == (80, 240),
        })
        for i in lost:
            diagnostics = ("RelativeCodeError", "EndpointClass", "FinalObjective", "FinalGradientNorm",
                           "Stop", "Updates", "NegativeAngleDegrees", "RadiusRatio", "RayOptimalRadius")
            reversals.append({
                "lower_m": lower, "higher_m": higher, "repetition": i, "seed": a[i - 1]["Seed"],
                "lower": {key: a[i - 1][key] for key in diagnostics},
                "higher": {key: b[i - 1][key] for key in diagnostics},
            })
    primary = next(p for p in pairs if p["primary"])
    patterns = [
        {"repetition": i, "seed": baseline[i - 1]["Seed"],
         "success_pattern": [rows_by_m[m][i - 1]["Success"] for m in counts]}
        for i in range(1, 31)
    ]
    audit = {
        "settings_sha256": settings_sha, "measurement_order": counts,
        "data_and_numerical_gates_passed": True, "baseline_fields_and_bits_preserved": True,
        "primary_increase_supported": primary["difference"] > 0 and primary["mcnemar_p"] < 0.05,
        "cells": cells, "pairs": pairs, "success_patterns": patterns, "reversals": reversals,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "analysis.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"cells": cells, "primary": primary, "reversals": reversals}, indent=2))


if __name__ == "__main__":
    main()
