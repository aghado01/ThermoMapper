"""Independent endpoint scoring and exact count comparisons. Python standard library only."""
import argparse
import csv
import hashlib
import json
import math
import pathlib
import struct
from collections import Counter


def read_rows(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def fisher(first, second):
    total = first + second
    weights = {
        i: math.comb(30, i) * math.comb(30, total - i)
        for i in range(max(0, total - 30), min(30, total) + 1)
    }
    observed = weights[first]
    return sum(v for v in weights.values() if v <= observed) / sum(weights.values())


def wilson(successes, n=30):
    z = 1.959963984540054
    p = successes / n
    center = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [center - half, center + half]


def close(actual, expected):
    assert abs(actual - expected) <= 2e-12 * max(1, abs(expected)), (actual, expected)


def validate(row, sign, rho=1 / math.pi):
    truth, estimate = row["Truth"], row["Estimate"]
    truth_norm, estimate_norm = math.hypot(*truth), math.hypot(*estimate)
    error = math.dist(truth, estimate) / truth_norm
    radius = estimate_norm / truth_norm
    product = math.fsum(a * b for a, b in zip(truth, estimate))
    cosine = max(-1, min(1, product / (truth_norm * estimate_norm))) if estimate_norm else 0
    angle = math.degrees(math.acos(-cosine)) if estimate_norm else None
    orientation = angle is not None and angle <= 20
    radius_match = abs(radius - rho) <= rho / 4
    classification = "recovered" if error < 0.001 else (
        "negative-neighbourhood" if orientation and radius_match else "elsewhere"
    )
    for key, value in {
        "RelativeCodeError": error, "RadiusRatio": radius, "CosineToTruth": cosine,
        "AxialProjection": product / (truth_norm * truth_norm),
        "PopulationRadiusDeviation": radius - rho,
    }.items():
        close(row[key], value)
    # acos near its endpoints is sensitive to a final-bit cosine difference.
    if angle is not None:
        assert abs(row["NegativeAngleDegrees"] - angle) < 3e-6
    assert row["Success"] == (error < 0.001)
    assert row["NegativeOrientation"] == orientation
    assert row["PopulationRadiusMatch"] == radius_match
    assert row["EndpointClass"] == classification
    assert row["FlatNegativeRay"] == (row["RaySquaredNorm"] == 0)
    if not row["FlatNegativeRay"]:
        optimum = max(0, row["RayInnerProduct"] / row["RaySquaredNorm"])
        close(row["RayOptimalRadius"], optimum)
        close(row["FiniteRayRadiusDeviation"], radius - optimum)
    assert row["SelectedPreactivationZeros"] == 0
    assert row["MaxReferenceDiscrepancy"] <= 1e-10
    assert row["Updates"] == row["EvaluatedUpdates"] + row["SkippedUpdates"]
    assert row["SignChecks"] == row["Updates"] + (row["Stop"] == "gradient")
    if row["Stop"] == "gradient":
        assert row["FinalGradientNorm"] < math.ulp(1.0)
    else:
        assert row["Stop"] == "update-limit" and row["Updates"] == 50000
    if row["SkippedUpdates"]:
        assert row["CycleLength"] > 0 and row["SkippedUpdates"] % row["CycleLength"] == 0
        assert row["Stop"] == "update-limit"
    if not sign:
        assert row["FlipIterations"] == []


def finite(value):
    if isinstance(value, float):
        assert math.isfinite(value)
    elif isinstance(value, dict):
        for item in value.values():
            finite(item)
    elif isinstance(value, list):
        for item in value:
            finite(item)


def compare(first, second):
    ignored = {"EvaluatedUpdates", "SkippedUpdates", "CycleLength"}
    assert len(first) == len(second)
    for a, b in zip(first, second):
        assert a["Repetition"] == b["Repetition"]
        for key in a.keys() & b.keys() - ignored:
            assert a[key] == b[key], (a["Repetition"], key, a[key], b[key])
        for key in ("Truth", "Start", "Estimate"):
            assert b"".join(struct.pack("d", v) for v in a[key]) == b"".join(
                struct.pack("d", v) for v in b[key]
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=pathlib.Path)
    parser.add_argument("--baseline", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--compare", nargs=2, type=pathlib.Path)
    args = parser.parse_args()
    if args.compare:
        compare(*(read_rows(p) for p in args.compare))
        print("Common numerical fields and endpoint bits agree.")
        return
    settings_path = pathlib.Path(__file__).parent / "settings/g1-linear-curve.json"
    settings = json.loads(settings_path.read_text())
    sha = hashlib.sha256(settings_path.read_bytes()).hexdigest()
    curve, pairs, summaries, all_rows, geometry = [], [], {}, {}, {}
    component_seeds = []
    for cell in settings["Cells"]:
        k = cell["K"]
        for sign in (True, False) if cell["Paired"] else (True,):
            name = f"k{k}-{'on' if sign else 'off'}"
            rows = read_rows(args.results / name / "trials.jsonl")
            summary = json.loads((args.results / name / "summary.json").read_text())
            assert summary["SettingsSha256"] == sha and summary["GatePassed"]
            assert len(rows) == 30 and [r["Repetition"] for r in rows] == list(range(1, 31))
            for row in rows:
                finite(row)
                validate(row, sign)
            assert summary["Successes"] == sum(r["Success"] for r in rows)
            assert summary["SkippedUpdates"] == sum(r["SkippedUpdates"] for r in rows)
            failed = [r for r in rows if not r["Success"]]
            ranges = {}
            for key in ("NegativeAngleDegrees", "RadiusRatio", "RayOptimalRadius", "FiniteRayRadiusDeviation"):
                values = [r[key] for r in failed if r[key] is not None]
                ranges[key] = [min(values), max(values)] if values else None
            geometry[name] = {
                "negative_orientation_count": sum(r["NegativeOrientation"] for r in rows),
                "population_radius_count": sum(r["PopulationRadiusMatch"] for r in rows),
                "negative_neighbourhood_count": sum(r["EndpointClass"] == "negative-neighbourhood" for r in rows),
                "flat_ray_count": sum(r["FlatNegativeRay"] for r in rows),
                "failed_endpoint_ranges": ranges,
            }
            if sign:
                component_seeds.extend(n for r in rows for n in r["ChildSeeds"])
            summaries[name], all_rows[name] = summary, rows
        on = all_rows[f"k{k}-on"]
        observed = sum(r["Success"] for r in on)
        p = fisher(observed, cell["PaperSuccesses"])
        curve.append({
            "k": k, "m_over_k": 150 / k, "successes": observed, "n": 30,
            "paper_successes": cell["PaperSuccesses"], "wilson95": wilson(observed),
            "paper_wilson95": wilson(cell["PaperSuccesses"]), "fisher_p": p,
            "discrepant": p < 0.05 / 15,
        })
        if cell["Paired"]:
            off = all_rows[f"k{k}-off"]
            for a, b in zip(on, off):
                for key in ("Seed", "ChildSeeds", "Truth", "Start", "RayInnerProduct", "RaySquaredNorm", "RayOptimalRadius"):
                    assert a[key] == b[key], (k, key)
            both = sum(a["Success"] and b["Success"] for a, b in zip(on, off))
            on_only = sum(a["Success"] and not b["Success"] for a, b in zip(on, off))
            off_only = sum(not a["Success"] and b["Success"] for a, b in zip(on, off))
            discordant = on_only + off_only
            paired_p = min(1, 2 * sum(math.comb(discordant, i) for i in range(min(on_only, off_only) + 1)) / 2**discordant)
            pairs.append({
                "k": k, "both": both, "on_only": on_only, "off_only": off_only,
                "neither": 30 - both - discordant, "difference": (on_only - off_only) / 30,
                "mcnemar_p": paired_p, "significant": paired_p < 0.05 / 3,
                "on_classes": dict(Counter(r["EndpointClass"] for r in on)),
                "off_classes": dict(Counter(r["EndpointClass"] for r in off)),
            })
    assert len(set(component_seeds)) == len(component_seeds)
    compare(read_rows(args.baseline), all_rows["k10-on"])
    audit = {
        "curve_gate_passed": not any(r["discrepant"] for r in curve), "curve": curve, "pairs": pairs,
        "numerical_gates_passed": True, "k10_baseline_fields_match": True,
        "distinct_component_seeds": len(component_seeds), "cells": summaries, "geometry": geometry,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "analysis.json").write_text(json.dumps(audit, indent=2) + "\n")
    with (args.output / "curve.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=curve[0].keys())
        writer.writeheader()
        writer.writerows(curve)
    print(json.dumps({"curve_gate": audit["curve_gate_passed"], "pairs": pairs}, indent=2))


if __name__ == "__main__":
    main()
