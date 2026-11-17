"""Reconstruct G2 summaries from retained states; compare the fixed replay case."""
import argparse
import collections
import itertools
import json
import math
import pathlib
import statistics


def load(path):
    with path.open() as stream:
        for line in stream:
            yield json.loads(line)


def key(r):
    return r.get("K", r.get("k")), r.get("Ratio", r.get("ratio")), r["Repetition"], r.get("Method", "problem")


def norm(v):
    return math.sqrt(math.fsum(x * x for x in v))


def without_timing(value):
    if isinstance(value, dict):
        return {k: without_timing(v) for k, v in value.items() if k not in ("Seconds", "DescentSeconds", "MedianSeconds")}
    if isinstance(value, list):
        return [without_timing(v) for v in value]
    return value


def selected(r):
    k, ratio, repetition, _ = key(r)
    return (k, ratio, repetition) == (2, 4, 1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", type=pathlib.Path, required=True)
    p.add_argument("--output", type=pathlib.Path, required=True)
    p.add_argument("--replay", type=pathlib.Path)
    a = p.parse_args()
    readings = list(load(a.results / "calibrations.jsonl"))
    problems = list(load(a.results / "problems.jsonl"))
    lookup = {key(r): r for r in readings}
    summary = json.loads((a.results / "summary.json").read_text())
    assert len(readings) == (216 if summary["ExpandedTo16"] else 162)
    assert len(problems) * 3 == len(readings) == len(lookup)
    raw_counts, checked = [], set()
    for identity, states in itertools.groupby(load(a.results / "posterior-draws.jsonl"), key):
        states = list(states)
        r = lookup[identity]
        checked.add(identity)
        k = identity[0]
        weights = [d["Weight"] for d in states]
        total = math.fsum(weights)
        assert abs(total - 1) < 1e-10
        assert all(len(d["Theta"]) == k and all(abs(x) <= 3 for x in d["Theta"]) for d in states)
        assert all(d["Variances"][0] > 0 and math.isfinite(d["LogLikelihood"]) and d["Weight"] >= 0 for d in states)
        means = [math.fsum(w * d["Theta"][j] for w, d in zip(weights, states)) / total for j in range(k)]
        assert max(abs(x - y) for x, y in zip(means, r["Mean"])) < 1e-10
        relative = norm([x - y for x, y in zip(means, r["Truth"])]) / norm(r["Truth"])
        assert abs(relative - r["RelativeError"]) < 1e-10
        for j in range(k):
            for probability, q in ((0.05, r["Q05"][j]), (0.95, r["Q95"][j])):
                before = math.fsum(w for w, d in zip(weights, states) if d["Theta"][j] < q) / total
                at = math.fsum(w for w, d in zip(weights, states) if d["Theta"][j] <= q) / total
                assert before <= probability + 1e-10 and at >= probability - 1e-10
        sigma = math.sqrt(math.fsum(w * d["Variances"][0] for w, d in zip(weights, states)) / total)
        assert abs(sigma - r["NoiseDeviation"]) < 1e-10
        assert r["Covered"] == sum(lo <= t <= hi for lo, t, hi in zip(r["Q05"], r["Truth"], r["Q95"]))
        negative = math.fsum(w for w, d in zip(weights, states) if norm(d["Theta"]) > 0 and
            math.fsum(x * y for x, y in zip(d["Theta"], r["Truth"])) / (norm(d["Theta"]) * norm(r["Truth"])) <= -math.cos(math.pi / 9) and
            abs(norm(d["Theta"]) / norm(r["Truth"]) - 1 / math.pi) <= 1 / (4 * math.pi)) / total
        assert abs(negative - r["NegativeMass"]) < 1e-10
        if r["Method"] == "smc":
            populations = [[d for d in states if d["Population"] == i] for i in range(3)]
            assert all(len(pop) == 256 for pop in populations) and 1 <= r["MinimumFamilies"] <= 256
            population_weights = [[d["Weight"] / math.fsum(x["Weight"] for x in pop) for d in pop] for pop in populations]
            minimum_ess = min(1 / math.fsum(w * w for w in weights) for weights in population_weights)
            spread = 0
            for j in range(k):
                pm = [math.fsum(w * d["Theta"][j] for w, d in zip(weights, pop)) for weights, pop in zip(population_weights, populations)]
                grand = statistics.fmean(pm)
                variance = statistics.fmean(math.fsum(w * (d["Theta"][j] - grand) ** 2 for w, d in zip(weights, pop))
                                            for weights, pop in zip(population_weights, populations))
                spread = max(spread, (max(pm) - min(pm)) / math.sqrt(variance) if variance > 0 else 0)
            assert abs(minimum_ess - r["MinimumParticleEss"]) < 1e-8 and abs(spread - r["PopulationSpread"]) < 1e-9
            assert r["DiagnosticGate"] == (minimum_ess >= 64 and spread <= 0.5)
        else:
            assert set(d["Chain"] for d in states) == {0, 1, 2}
            assert "unconverged" in r["Pooling"] if not r["DiagnosticGate"] else r["Pooling"] == "converged inner ensemble"
        raw_counts.append({"key": identity, "states": len(states)})
    assert checked == set(lookup)
    for g in summary["Groups"]:
        selected_rows = [r for r in readings if (r["K"], r["Ratio"], r["Method"]) == (g["K"], g["Ratio"], g["Method"])]
        assert len(selected_rows) == 6
        median = statistics.median(r["RelativeError"] for r in selected_rows)
        passes = sum(r["DiagnosticGate"] for r in selected_rows)
        assert abs(median - g["MedianError"]) < 1e-12 and passes == g["DiagnosticPasses"]
        assert g["Passed"] == (median <= 0.1 and passes >= 5)
    grid = json.loads((a.results / "grid-control.json").read_text())
    mean_gap = max(abs(x - y) for x, y in zip(grid["Coarse"]["Mean"], grid["Fine"]["Mean"]))
    quantile_gap = max(abs(x - y) for name in ("Lower", "Upper") for x, y in zip(grid["Coarse"][name], grid["Fine"][name]))
    assert grid["Passed"] and mean_gap <= 0.005 and quantile_gap <= 0.02
    alignments = list(load(a.results / "alignment-states.jsonl"))
    alignment_summary = []
    for identity, group in itertools.groupby(sorted(alignments, key=lambda r: (r["k"], r["Method"], r["RelativeError"] > 0.1)),
                                               lambda r: (r["k"], r["Method"], r["RelativeError"] > 0.1)):
        group = list(group)
        assert all(-1 - 1e-12 <= r[n] <= 1 + 1e-12 for r in group for n in ("CleanCosine", "NoisyCosine"))
        alignment_summary.append({"k": identity[0], "method": identity[1], "above_error_0_1": identity[2], "states": len(group),
            "median_clean_cosine": statistics.median(r["CleanCosine"] for r in group),
            "median_noisy_cosine": statistics.median(r["NoisyCosine"] for r in group),
            "positive_clean_fraction": sum(r["CleanCosine"] > 0 for r in group) / len(group),
            "positive_noisy_fraction": sum(r["NoisyCosine"] > 0 for r in group) / len(group)})
    descent = []
    for k, ratio in sorted({(r["k"], r["ratio"]) for r in problems}):
        group = [r for r in problems if (r["k"], r["ratio"]) == (k, ratio)]
        errors = [norm([x - y for x, y in zip(r["Descent"]["X"], r["Truth"])]) / norm(r["Truth"]) for r in group]
        descent.append({"k": k, "ratio": ratio, "median_error": statistics.median(errors), "errors": errors,
                        "below_0_1": sum(e <= 0.1 for e in errors), "median_seconds": statistics.median(r["DescentSeconds"] for r in group),
                        "outside_prior_box": sum(any(abs(x) > 3 for x in r["Descent"]["X"]) for r in group),
                        "stops": dict(collections.Counter(r["Descent"]["Stop"] for r in group))})
    replay_checks = []
    if a.replay:
        for name in ("calibrations.jsonl", "problems.jsonl", "posterior-draws.jsonl", "alignment-states.jsonl"):
            primary = (without_timing(r) for r in load(a.results / name) if selected(r))
            replay = (without_timing(r) for r in load(a.replay / name))
            count = 0
            for left, right in itertools.zip_longest(primary, replay):
                assert left == right, (name, count)
                count += 1
            replay_checks.append({"file": name, "rows_equal_excluding_timing": count})
        assert json.loads((a.replay / "grid-control.json").read_text()) == grid
    noise = []
    for g in summary["Groups"]:
        group = [r for r in readings if (r["K"], r["Ratio"], r["Method"]) == (g["K"], g["Ratio"], g["Method"])]
        actual = 0.1 / math.sqrt(g["K"] * g["Ratio"])
        noise.append({"k": g["K"], "ratio": g["Ratio"], "method": g["Method"], "actual_deviation": actual,
                      "median_posterior_rms_deviation": statistics.median(r["NoiseDeviation"] for r in group)})
    failures = [{name: r[name] for name in ("K", "Ratio", "Repetition", "Method", "MaximumRHat", "MinimumEssPerChain",
                "MinimumParticleEss", "MinimumFamilies", "PopulationSpread")} for r in readings if not r["DiagnosticGate"]]
    result = {"passed": True, "calibrations": len(readings), "problems": len(problems), "states": sum(r["states"] for r in raw_counts),
        "raw_state_counts": raw_counts, "grid_mean_gap": mean_gap, "grid_quantile_gap": quantile_gap,
        "replay": replay_checks, "descent": descent, "alignment_diagnostic_subsets": alignment_summary,
        "noise": noise, "diagnostic_failures": failures,
        "chain_metric_scope": "reported max R-hat/min ESS cover theta and log variance; the original engine verdict also gates log likelihood"}
    a.output.mkdir(parents=True, exist_ok=True)
    (a.output / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("passed", "calibrations", "states", "grid_mean_gap", "grid_quantile_gap", "replay")}))


if __name__ == "__main__":
    main()
