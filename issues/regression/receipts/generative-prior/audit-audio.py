"""Audit saved normalized secant tables and their scalar minimax scale fits."""
import argparse
import json
import math
import pathlib
import statistics


def interval(a, c, error):
    lower, upper = 0.0, math.inf
    for x, y in zip(a, c):
        if x == 0:
            if abs(y) > error:
                return False
        else:
            bounds = sorted(((y - error) / x, (y + error) / x))
            lower, upper = max(lower, bounds[0]), min(upper, bounds[1])
            if upper < lower:
                return False
    return True


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", type=pathlib.Path, required=True)
    p.add_argument("--output", type=pathlib.Path, required=True)
    a = p.parse_args()
    rows = [json.loads(x) for x in (a.results / "readings.jsonl").read_text().splitlines()]
    secants = [json.loads(x) for x in (a.results / "secants.jsonl").read_text().splitlines()]
    summary = json.loads((a.results / "summary.json").read_text())
    assert len(rows) == 470 and len(secants) == 640 and summary["NumericalControlsPassed"]
    assert all(s["Norm"] > 0 and math.isfinite(s["Norm"]) for s in secants)
    original = {}
    for r in rows:
        measured, reference = r["MeasuredInnerProducts"], r["OriginalInnerProducts"]
        assert len(measured) == len(reference) == 96 and reference[:64] == [1] * 64
        assert all(math.isfinite(x) for x in measured + reference)
        assert r["Scale"] >= 0 and 0 <= r["Rank"] <= r["Rows"]
        assert r["RankIsJacobian"] == (r["Map"] == "kept-log-power")
        assert r["WhiteningError"] <= 1e-5 and r["ScaleCheckError"] <= 1e-8
        assert r["Condition"] <= 1e8 * (1 + 1e-10)
        objective = max(abs(r["Scale"] * x - y) for x, y in zip(measured, reference))
        assert abs(objective - r["Epsilon"]) < 1e-12
        assert interval(measured, reference, objective + 1e-12)
        assert not interval(measured, reference, max(0, objective - 1e-9))
        key = r["Manifold"], r["Repetition"]
        assert original.setdefault(key, reference) == reference
    decisions = []
    for manifold, name in sorted({(r["Manifold"], r["Map"]) for r in rows}):
        groups = [g for g in summary["Groups"] if g["Manifold"] == manifold and g["Map"] == name]
        passing = [g["Rows"] for g in groups if g["SampledNearIsometry"]]
        decisions.append({"manifold": manifold, "map": name, "first_sampled_rows": min(passing) if passing else None,
                          "widest": max(groups, key=lambda g: g["Rows"])})
    result = {"passed": True, "readings": len(rows), "secants": len(secants), "decisions": decisions,
              "max_whitening_error": max(r["WhiteningError"] for r in rows),
              "max_eigen_residual": max(r["EigenResidual"] for r in rows),
              "max_scale_control_error": max(r["ScaleCheckError"] for r in rows)}
    a.output.mkdir(parents=True, exist_ok=True)
    (a.output / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
