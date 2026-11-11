"""Independent recomputation of F-R1's gate and rule from the per-fit CSV (the review's hand check)."""
import csv
import math
import sys

GATE, CLEARS, CEILING, MARGIN = 1.05, 0.9, 1.05, 0.05
ARMS = [f"A{i}" for i in range(8)]


def rule(rows, benchmarks, designs):
    at = {(r["arm"], r["benchmark"], int(r["k"])): r for r in rows}
    excluded = set()
    for r in rows:
        if r["benchmark"] in benchmarks and int(r["k"]) < designs:
            rhat = float(r["predictive_rhat"]) if r["predictive_rhat"] else math.inf
            if not rhat <= GATE:
                excluded.add((r["benchmark"], int(r["k"])))
    scored = {}
    for f in benchmarks:
        keep = [k for k in range(designs) if (f, k) not in excluded]
        if designs - len(keep) < 2 and keep:
            scored[f] = keep
    out = {}
    for arm in ARMS[1:]:
        lr, ok = {}, True
        for f, keep in scored.items():
            lr[f] = (sum(math.log(float(at[(arm, f, k)]["nrmse"])) for k in keep) / len(keep)
                     - sum(math.log(float(at[("A0", f, k)]["nrmse"])) for k in keep) / len(keep))
            cov = sum(float(at[(arm, f, k)]["coverage"]) for k in keep) / len(keep)
            cov0 = sum(float(at[("A0", f, k)]["coverage"]) for k in keep) / len(keep)
            ok &= math.exp(lr[f]) <= CEILING and cov0 - cov <= MARGIN
        score = math.exp(sum(lr.values()) / len(lr)) if lr else float("nan")
        out[arm] = (score, {f: math.exp(v) for f, v in lr.items()}, bool(lr) and ok and score <= CLEARS)
    clearing = [(s, a) for a, (s, _, c) in out.items() if c]
    return sorted(excluded), scored, out, (min(clearing)[1] if clearing else None)


rows = list(csv.DictReader(open(sys.argv[1], newline="", encoding="utf-8")))
designs = int(sys.argv[2]) if len(sys.argv) > 2 else 3
for name, benches in (("noise-free", ["friedman", "ishigami", "borehole", "functional"]),
                      ("noisy controls", ["friedman-noisy", "ishigami-noisy"])):
    excluded, scored, out, winner = rule(rows, benches, designs)
    print(f"== {name}: excluded pairs {excluded}; scored {scored}; winner {winner}")
    for arm, (score, rf, clears) in out.items():
        print(f"  {arm}: score {score:.6f} r_f {{{', '.join(f'{f}: {v:.6f}' for f, v in rf.items())}}} clears {clears}")
