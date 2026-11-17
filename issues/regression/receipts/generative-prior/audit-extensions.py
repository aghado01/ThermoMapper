"""Independent geometry, paired outcomes and analytic radius controls for G1."""
import argparse
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
    config = json.loads((base / "settings/g1-extensions.json").read_text())
    readings, comparisons, cycles = [], [], {}
    for cell in config["ExtendedCells"]:
        if cell.get("Bridge"):
            continue
        angle = math.pi
        for _ in range(cell["Depth"]):
            angle = math.acos(max(-1, min(1, ((math.pi - angle) * math.cos(angle) + math.sin(angle)) / math.pi)))
        rho = math.cos(angle)
        if cell["Magnitude"]:
            rho = 2 * math.sin(angle) / math.pi + (math.pi - 2 * angle) / math.pi * rho
        arms = {}
        for sign in (True, False):
            name = cell["Name"] + ("-on" if sign else "-off")
            rows = common["read_rows"](args.results / name / "trials.jsonl")
            summary = json.loads((args.results / name / "summary.json").read_text())
            assert len(rows) == 30 and summary["GatePassed"]
            for row in rows:
                common["finite"](row)
                common["validate"](row, sign, rho)
                common["close"](row["PredictedNegativeRadius"], rho)
                assert row["SelectedMeasurementZeros"] == 0
                if row["CycleLength"]:
                    cycles.setdefault(str(row["CycleLength"]), {"cell": cell["Name"], "sign": sign,
                                                               "repetition": row["Repetition"]})
            assert summary["Successes"] == sum(r["Success"] for r in rows)
            arms[sign] = rows
            failed = [r for r in rows if not r["Success"]]
            readings.append({
                "cell": name, "k": cell["K"], "depth": cell["Depth"], "magnitude": cell["Magnitude"],
                "rho": rho, "successes": summary["Successes"], "wilson95": common["wilson"](summary["Successes"]),
                "gradient_stops": summary["GradientStops"], "limit_stops": summary["UpdateLimitStops"],
                "classes": summary["EndpointClasses"], "max_reference_discrepancy": summary["MaxReferenceDiscrepancy"],
                "negative_orientation_count": sum(r["NegativeOrientation"] for r in rows),
                "radius_match_count": sum(r["PopulationRadiusMatch"] for r in rows),
                "failed_angles": [min(r["NegativeAngleDegrees"] for r in failed), max(r["NegativeAngleDegrees"] for r in failed)] if failed else None,
            })
        on, off = arms[True], arms[False]
        for a, b in zip(on, off):
            for key in ("Seed", "ChildSeeds", "Truth", "Start", "RayOptimalRadius"):
                assert a[key] == b[key]
        both = sum(a["Success"] and b["Success"] for a, b in zip(on, off))
        only_on = [a["Repetition"] for a, b in zip(on, off) if a["Success"] and not b["Success"]]
        only_off = [a["Repetition"] for a, b in zip(on, off) if b["Success"] and not a["Success"]]
        n = len(only_on) + len(only_off)
        p = min(1, 2 * sum(math.comb(n, i) for i in range(min(len(only_on), len(only_off)) + 1)) / 2**n)
        comparisons.append({"cell": cell["Name"], "both": both, "on_only": only_on, "off_only": only_off,
                            "neither": 30 - both - n, "mcnemar_p": p, "family_significant": p < 0.05 / 7})
    args.output.mkdir(parents=True, exist_ok=True)
    audit = {"numerical_gates_passed": True, "readings": readings, "paired": comparisons,
             "cycle_representatives": cycles, "controls": {"rho_linear_2": 1 / math.pi,
                                                         "rho_magnitude_2": 0.6691522942169524}}
    (args.output / "analysis.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"readings": readings, "paired": comparisons, "cycles": cycles}, indent=2))


if __name__ == "__main__":
    main()
