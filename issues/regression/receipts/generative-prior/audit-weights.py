"""Independent reading of stored G3 scalar witnesses and matrix controls."""
import argparse
import json
import math
import pathlib
import statistics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=pathlib.Path)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    args = parser.parse_args()
    rows = [json.loads(line) for line in (args.results / "readings.jsonl").read_text().splitlines()]
    summary = json.loads((args.results / "summary.json").read_text())
    assert len(rows) == 720 and summary["NumericalControlsPassed"] and summary["WitnessGatesPassed"]
    for row in rows:
        assert row["MaxGap"] + 1e-12 >= row["WitnessGap"] >= abs(row["WitnessDiagonal2"] - 0.5) - 1e-12
        assert row["MaxEigenvalueCheck"] <= 1e-9 * max(1, row["MaxGap"])
        assert row["MaxKnotCheck"] <= 1e-10
    for item in summary["Witnesses"]:
        selected = [r["WitnessDiagonal2"] for r in rows if r["K"] == item["K"] and
                    r["Law"] == item["Law"] and r["N"] == 512 * item["K"]]
        mean = statistics.fmean(selected)
        se = statistics.stdev(selected) / math.sqrt(len(selected))
        assert abs(mean - item["Mean"]) <= 1e-12 and abs(se - item["StandardError"]) <= 1e-12
        assert abs(mean - item["Expected"]) <= 4 * se + 1e-10
    controls = [json.loads(line) for line in (args.results / "control-matrices.jsonl").read_text().splitlines()]
    assert len(controls) == 180
    for control in controls:
        matrix = control["Residual"]
        assert max(abs(matrix[i][j] - matrix[j][i]) for i in range(8) for j in range(8)) < 1e-12
        assert control["Norm"] + 1e-12 >= max(abs(matrix[i][i]) for i in range(8))
        if control["Pair"] == "e1-minus-e1" and control["Law"] != "sparse-threshold":
            assert control["Norm"] == 0 and all(v == 0 for row in matrix for v in row)
    decisions = []
    keys = sorted({(r["K"], r["Law"], r["Radius"]) for r in rows})
    for k, law, radius in keys:
        groups = [g for g in summary["Groups"] if g["K"] == k and g["Law"] == law and g["Radius"] == radius]
        passing = [g["Ratio"] for g in groups if g["SampledEpsilonReached"]]
        decisions.append({"k": k, "law": law, "radius": radius, "first_sampled_ratio": min(passing) if passing else None,
                          "widest_median_gap": max(groups, key=lambda g: g["N"])["MedianGap"]})
    args.output.mkdir(parents=True, exist_ok=True)
    audit = {"all_stored_controls_passed": True, "decisions": decisions, "witnesses": summary["Witnesses"]}
    (args.output / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
