"""Saved-receipt noise analysis. No simulator or inversion calls."""
import argparse
import csv
import hashlib
import itertools
import json
import math
import pathlib
import statistics

MASK = (1 << 64) - 1


def splitmix(state):
    state = (state + 0x9E3779B97F4A7C15) & MASK
    z = state
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
    return state, z ^ (z >> 31)


def derive(master, count):
    state = master & MASK
    values = []
    for _ in range(count):
        state, value = splitmix(state)
        low = value & 0xFFFFFFFF
        values.append(low if low < 1 << 31 else low - (1 << 32))
    return values


class Xoshiro:
    def __init__(self, seed):
        state = seed & MASK
        self.words = []
        for _ in range(4):
            state, value = splitmix(state)
            self.words.append(value)

    def next(self):
        def rotate(value, bits):
            return ((value << bits) | (value >> (64 - bits))) & MASK
        a, b, c, d = self.words
        result = (rotate((a + d) & MASK, 23) + a) & MASK
        t = (b << 17) & MASK
        c ^= a
        d ^= b
        b ^= c
        a ^= d
        c ^= t
        d = rotate(d, 45)
        self.words = [a, b, c, d]
        return result

    def below(self, count):
        return ((self.next() >> 32) * count) >> 32


def percentile(values, probability):
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lo = math.floor(position)
    hi = math.ceil(position)
    return ordered[lo] + (position - lo) * (ordered[hi] - ordered[lo])


def slopes(rows, offset=0.0):
    groups = {}
    for row in rows:
        error = max(0, float(row["rmse"]) + offset)
        if error == 0:
            continue
        groups.setdefault(int(row["k"]), []).append(
            (math.log(int(row["k"]) / int(row["m"])), math.log(float(row["noise"])), 2 * math.log(error))
        )
    xx = zz = xz = xy = zy = 0.0
    count = 0
    for values in groups.values():
        means = [statistics.fmean(v[i] for v in values) for i in range(3)]
        for x, z, y in values:
            x -= means[0]
            z -= means[1]
            y -= means[2]
            xx += x * x
            zz += z * z
            xz += x * z
            xy += x * y
            zy += z * y
            count += 1
    if xx == 0:
        return {"ratio_slope": None, "noise_slope": None, "rows": count}
    if zz < 1e-20:
        return {"ratio_slope": xy / xx, "noise_slope": None, "rows": count}
    determinant = xx * zz - xz * xz
    return {"ratio_slope": (xy * zz - zy * xz) / determinant,
            "noise_slope": (zy * xx - xy * xz) / determinant, "rows": count}


def read_originals(config, root):
    rows, manifest = [], []
    for stage in config["stages"]:
        for name in stage["files"]:
            path = root / name
            raw = path.read_bytes()
            selected = 0
            with path.open(newline="") as stream:
                reader = csv.DictReader(stream)
                for row in reader:
                    if row["emulator"] != "oracle" or row["mode"] != "nested" or row["disc"] != "0":
                        continue
                    rows.append({key: row[key] for key in ("stage", "k", "rep", "noise", "ratio", "m", "converged", "runs", "seed")} |
                                {"rmse": row["recon_truth_held"], "source_file": name, "source_line": reader.line_num})
                    selected += 1
            manifest.append({"file": name, "sha256": hashlib.sha256(raw).hexdigest(), "selected_rows": selected})
    return rows, manifest


def analyze(rows, config):
    results = []
    for stage_index, stage in enumerate(config["stages"]):
        selected = [r for r in rows if r["stage"] == stage["name"]]
        cells = []
        key = lambda r: (int(r["k"]), float(r["noise"]), int(r["ratio"]))
        for (k, noise, ratio), group in itertools.groupby(sorted(selected, key=key), key=key):
            group = list(group)
            squared = [float(r["rmse"]) ** 2 for r in group]
            half = config["rmse_quantum"] / 2
            cells.append({
                "k": k, "noise": noise, "ratio": ratio, "k_over_m": 1 / ratio, "rows": len(group),
                "median_mse": statistics.median(squared), "q25_mse": percentile(squared, 0.25),
                "q75_mse": percentile(squared, 0.75),
                "rounding_lower_median": statistics.median(max(0, float(r["rmse"]) - half) ** 2 for r in group),
                "rounding_upper_median": statistics.median((float(r["rmse"]) + half) ** 2 for r in group),
                "converged_rows": sum(int(r["converged"]) == int(r["runs"]) for r in group),
            })
        by_k = {}
        for row in selected:
            by_k.setdefault(int(row["k"]), {}).setdefault(int(row["rep"]), []).append(row)
        clusters = [list(group.values()) for _, group in sorted(by_k.items())]
        stage_seed = derive(config["master_seed"], len(config["stages"]))[stage_index]
        samples = []
        for seed in derive(stage_seed, config["bootstrap_repetitions"]):
            rng = Xoshiro(seed)
            bootstrap = [r for group in clusters for _ in range(len(group)) for r in group[rng.below(len(group))]]
            samples.append(slopes(bootstrap))
        fit = slopes(selected)
        intervals = {}
        for name in ("ratio_slope", "noise_slope"):
            values = [s[name] for s in samples if s[name] is not None]
            intervals[name] = [percentile(values, 0.025), percentile(values, 0.975)] if values else None
        ci = intervals["ratio_slope"]
        converged = [r for r in selected if int(r["converged"]) == int(r["runs"])]
        results.append({
            "stage": stage["name"], "rows": len(selected), "converged_rows": len(converged),
            "zero_rounded_rmse_rows": sum(float(r["rmse"]) == 0 for r in selected),
            "fit": fit, "bootstrap95": intervals, "positive_ratio_trend": ci is not None and ci[0] > 0,
            "compatible_with_unit_ratio_slope": ci is not None and ci[0] <= 1 <= ci[1],
            "converged_only_fit": slopes(converged),
            "rounding_lower_fit": slopes(selected, -config["rmse_quantum"] / 2),
            "rounding_upper_fit": slopes(selected, config["rmse_quantum"] / 2), "cells": cells,
        })
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--settings", type=pathlib.Path, default=pathlib.Path(__file__).parent / "settings/g5-noise.json")
    parser.add_argument("--input-root", type=pathlib.Path)
    parser.add_argument("--rows", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.settings.read_text())
    if args.rows:
        with args.rows.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        manifest = [{"portable_rows": str(args.rows), "sha256": hashlib.sha256(args.rows.read_bytes()).hexdigest()}]
    else:
        rows, manifest = read_originals(config, args.input_root)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "oracle-rows.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    result = {"settings_sha256": hashlib.sha256(args.settings.read_bytes()).hexdigest(),
              "source_manifest": manifest, "stages": analyze(rows, config)}
    (args.output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps([{key: value for key, value in s.items() if key != "cells"} for s in result["stages"]], indent=2))


if __name__ == "__main__":
    main()
