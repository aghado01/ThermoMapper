"""Independent readings of saved conformal-audio predictions and deployment controls."""
import argparse
import collections
import cmath
import hashlib
import itertools
import json
import math
import pathlib
import statistics
from fractions import Fraction


def rows(path):
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            yield json.loads(line)


def mean_interval(values):
    mean = statistics.fmean(values)
    se = statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else 0
    return {"mean": mean, "se": se, "lower": mean - 1.959963984540054 * se, "upper": mean + 1.959963984540054 * se, "splits": len(values)}


def controls(path):
    saved = json.loads((path / "controls.json").read_text())
    signals = [(1, 0, 1), (1, 0, -1)]
    power = lambda x, k: abs(sum(v * cmath.exp(-2j * math.pi * k * n / len(x)) for n, v in enumerate(x))) ** 2
    assert [abs(x) for x in signals[0]] == [abs(x) for x in signals[1]]
    for start in (0, 1):
        for k in (0, 1):
            assert abs(power(signals[0][start:start + 2], k) - power(signals[1][start:start + 2], k)) < 1e-10
    oracle = []
    for signal in signals:
        visible = power(signal, 1)
        predicted = 4 if abs(visible - 1) < 1e-10 else 0
        target = power(signal, 0)
        assert abs(predicted - target) < 1e-10
        oracle.append({"visible_frequency_one_power": visible, "predicted_dc_power": predicted, "target_dc_power": target})
    coverage = sum(Fraction(1, 2) * (Fraction(9, 10) + Fraction(1, 10) * (target == 0)) for target in (0, 4))
    cardinality = Fraction(9, 10) * 2 + Fraction(1, 10)
    assert coverage == Fraction(19, 20) and cardinality == Fraction(19, 10)
    assert saved["FrameMissProbability"] == 0.25 and saved["GroupMissProbability"] == 0 and saved["LevelAt18"] == 0 and saved["LevelAt19"] == 1
    return {"passed": True, "independent_dft_oracle": oracle, "exact_randomized_coverage": str(coverage), "exact_cardinality": str(cardinality)}


def audio(path, replay):
    split_rows = list(rows(path / "audio-splits.jsonl"))
    split_key = lambda r: (r["Gap"], r["Split"], r["Predictor"], r["GroupEvent"])
    splits = {split_key(r): r for r in split_rows}
    assert len(splits) == len(split_rows) == 480
    tests = {(split_key(r), r["Conversation"]): r for r in rows(path / "audio-tests.jsonl")}
    assert len(tests) == 48000
    sequence = collections.Counter()
    matched = {}
    squared = collections.defaultdict(lambda: [0.0, 0])
    checked = 0
    for r in rows(path / "audio-predictions.jsonl"):
        identity = r["Gap"], r["Split"], r["Predictor"]
        index = sequence[identity]
        sequence[identity] += 1
        signature = (r["Mark"], r["NoiseSd"], hashlib.sha256(json.dumps(r["Targets"]).encode()).hexdigest())
        pair = r["Gap"], r["Split"], r["Seed"]
        assert matched.setdefault(pair, signature) == signature
        assert len(r["Targets"]) == len(r["Centers"]) == len(r["Scores"]) in (1, 4, 12)
        errors = [[abs(y - m) for y, m in zip(yrow, mrow)] for yrow, mrow in zip(r["Targets"], r["Centers"])]
        scores = [max(e) / r["Scale"] for e in errors]
        assert max(abs(x - y) for x, y in zip(scores, r["Scores"])) < 1e-10
        for event in ("OneCase", "WorstCase"):
            key = identity + (event,)
            s = splits[key]
            t = tests[key, index]
            assert t["MarkedCovered"] == (scores[r["Mark"]] <= s["Threshold"])
            assert t["AllCovered"] == (max(scores) <= s["Threshold"])
            assert abs(t["Width"] - 2 * s["Threshold"] * r["Scale"]) < 1e-10
            assert t["Length"] == len(scores) and t["NoiseSd"] == r["NoiseSd"]
            for yrow, mrow in zip(r["Targets"], r["Centers"]):
                for y, m in zip(yrow, mrow):
                    squared[key][0] += (y - m) ** 2
                    squared[key][1] += 1
            band = all(abs(y - m) <= s["Threshold"] * r["Scale"] + 1e-12 for yrow, mrow in zip(r["Targets"], r["Centers"]) for y, m in zip(yrow, mrow))
            assert band == t["AllCovered"]
        checked += 1
    assert checked == 24000 and len(matched) == 8000 and set(sequence.values()) == {100}
    for key, s in splits.items():
        group = [tests[key, i] for i in range(100)]
        assert s["CalibrationGroups"] == 99 and s["TestGroups"] == 100 and not s["FullSpace"]
        assert abs(sum(r["MarkedCovered"] for r in group) / 100 - s["MarkedCoverage"]) < 1e-12
        assert abs(sum(r["AllCovered"] for r in group) / 100 - s["AllCoverage"]) < 1e-12
        assert abs(math.sqrt(squared[key][0] / squared[key][1]) - s["Rmse"]) < 1e-10
    subgroup = collections.defaultdict(lambda: collections.defaultdict(list))
    for (key, _), r in tests.items():
        gap, split, predictor, event = key
        subgroup[gap, predictor, event, r["Length"], r["NoiseSd"]][split].append(r)
    diagnostics = []
    for (gap, predictor, event, length, noise), split_groups in sorted(subgroup.items()):
        name = "MarkedCovered" if event == "OneCase" else "AllCovered"
        diagnostics.append({"gap": gap, "predictor": predictor, "event": event, "length": length, "noise_sd": noise,
            "groups": sum(len(g) for g in split_groups.values()), "coverage": mean_interval([statistics.fmean(r[name] for r in g) for g in split_groups.values()])})
    replay_checks = []
    if replay:
        for name in ("audio-splits.jsonl", "audio-tests.jsonl", "audio-predictions.jsonl"):
            primary = (r for r in rows(path / name) if r["Gap"] == 64 and r["Split"] == 1)
            repeat = rows(replay / name)
            count = 0
            for left, right in itertools.zip_longest(primary, repeat):
                assert left == right, (name, count)
                count += 1
            replay_checks.append({"file": name, "exact_rows": count})
        for name in ("fit-gap64-short.json", "fit-gap64-plus.json"):
            left, right = [json.loads((p / name).read_text(encoding="utf-8")) for p in (path, replay)]
            left.pop("Seconds"); right.pop("Seconds")
            assert left == right
        assert (path / "training-gap64.json").read_bytes() == (replay / "training-gap64.json").read_bytes()
    return {"passed": True, "conversations_matched": len(matched), "prediction_rows_checked": checked,
            "test_events_checked": len(tests), "split_rows_checked": len(splits), "subgroup_diagnostics": diagnostics, "replay": replay_checks}


def transport(path):
    readings = list(rows(path / "transport-counts.jsonl"))
    summary = json.loads((path / "transport-summary.json").read_text())
    assert len(readings) == 5050
    probability = miss = size = 0
    for r in readings:
        a, b, c = r["Low"], r["PositiveHigh"], r["NegativeHigh"]
        p = math.comb(99, a) * math.comb(99 - a, b) * 0.95 ** a * 0.025 ** (b + c)
        assert abs(p - r["Probability"]) < 1e-12
        probability += p
        accepted = [10 * (190 * c + w) > 10 * a + 190 * (b + c) + w for w in (10, 190)]
        assert accepted == r["WeightedNegativeAccepted"]
        miss += p * (0 if accepted[1] else 0.25)
        size += p * (1 + sum(accepted) / 2)
    assert abs(probability - 1) < 1e-12 and abs(miss - summary["WeightedMarginalMiss"]) < 1e-12
    assert abs(size - summary["WeightedCardinality"]) < 1e-12
    return {"passed": True, "count_states": len(readings), "probability": probability, "weighted_miss": miss, "weighted_size": size}


def risk(path):
    split = list(rows(path / "risk-splits.jsonl"))
    total = collections.defaultdict(lambda: [0.0, 0.0, 0, 0.0, 0])
    policy_totals = collections.defaultdict(lambda: [0.0, 0.0, 0, 0, 0.0])
    delays = [0, 0.1, 0.2, 0.4, 0.8, None]
    for r in rows(path / "risk-conversations.jsonl"):
        delay = delays[r["Policy"]]
        starts = early = 0
        for p, b, e in zip(r["Prefix"], r["Budget"], r["Completion"]):
            if delay is not None and 0.1 + delay <= b:
                starts += 1
                early += p + 0.1 + delay < e
        n = len(r["Prefix"])
        loss, rate = (early / n, starts / n) if n else (0, 0)
        assert loss == r["Loss"] and rate == r["StartRate"] and starts == r["Starts"]
        total[r["Split"]][0] += loss / 500
        total[r["Split"]][1] += rate / 500
        total[r["Split"]][2] += starts
        total[r["Split"]][3] += starts * (0.1 + delay) if delay is not None else 0
        total[r["Split"]][4] += 1
        for policy, wait in enumerate(delays):
            count = early_count = 0
            for prefix, budget, completion in zip(r["Prefix"], r["Budget"], r["Completion"]):
                if wait is not None and 0.1 + wait <= budget:
                    count += 1
                    early_count += prefix + 0.1 + wait < completion
            value = policy_totals[r["Split"], policy]
            value[0] += early_count / n / 500 if n else 0
            value[1] += count / n / 500 if n else 0
            value[2] += count
            value[3] += early_count
            value[4] += count * (0.1 + wait) if wait is not None else 0
    for r in split:
        values = total[r["Split"]]
        assert values[4] == 500 and abs(values[0] - r["TestRisk"]) < 1e-12 and abs(values[1] - r["TestStartRate"]) < 1e-12
        assert values[2] == r["Starts"] and abs(values[3] - r["DelaySum"]) < 1e-9
    assert len(split) == 200
    curves = []
    for policy, delay in enumerate(delays):
        values = [policy_totals[rep, policy] for rep in total]
        count = sum(v[2] for v in values)
        curves.append({"policy": policy, "delay": delay, "eligible_turn_risk": mean_interval([v[0] for v in values]),
            "start_fraction": mean_interval([v[1] for v in values]), "starts": count, "early_starts": sum(v[3] for v in values),
            "pooled_delay_among_starts": sum(v[4] for v in values) / count if count else None})
    return {"passed": True, "splits": len(split), "conversations": sum(r[4] for r in total.values()),
            "risk": mean_interval([r["TestRisk"] for r in split]), "fixed_policy_curves": curves}


def alarm(path):
    counts = collections.Counter()
    alarms = collections.Counter()
    for r in rows(path / "alarm-streams.jsonl"):
        counts[r["Law"]] += 1
        alarms[r["Law"]] += r["FirstAlarm"] is not None
        assert r["FirstAlarm"] is None or 1 <= r["FirstAlarm"] <= 200
    summary = json.loads((path / "alarm-summary.json").read_text())
    for r in summary["Streams"]:
        assert counts[r["Law"]] == r["EverAlarm"]["Total"] and alarms[r["Law"]] == r["EverAlarm"]["Count"]
    return {"passed": True, "streams": sum(counts.values()), "counts": dict(counts), "alarms": dict(alarms)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("block", choices=["controls", "audio", "transport", "risk", "alarm"])
    p.add_argument("--results", type=pathlib.Path, required=True)
    p.add_argument("--replay", type=pathlib.Path)
    p.add_argument("--output", type=pathlib.Path, required=True)
    a = p.parse_args()
    result = audio(a.results, a.replay) if a.block == "audio" else globals()[a.block](a.results)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "subgroup_diagnostics"}))


if __name__ == "__main__":
    main()
