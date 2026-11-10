"""Searches behind slice 9's open items, beside conformal-checks.py's exact assertions.

1. CV+'s tournament bound (CV-b) on a fixed bag of n + 1 units, without the further test units.
   For each left-out set L (a fold with the test, m + 1 units), a symmetric fit ranks the residuals;
   d_L(x) counts the units of L ranked strictly above x, a bijection L -> {0..m} for strict orders
   (ties only help the test). A uniform ordering makes the test x uniform and the folds a uniform
   partition of the rest, and CV+ refuses iff sum_k d_{F_k + x}(x) <= l - 1. The search maximizes
   the refused share over the rankings by annealing; each figure is an exact count for the rankings
   found, so a lower bound on the worst case. Standard library only.
2. With scipy, the same worst case exactly by a MILP at K = 2, m = 2 (every level) and K = 3, m = 2,
   l = 1, and the jackknife+'s worst case for fits that read their units' order (model A, positions
   closed up) or know their units' positions (model B), m <= 5.

Results of 2026-11-09 (recorded in conformal-lemmas.md, slice 9):
- the fixed bag stays within CV-b's bound everywhere searched: K = 3, 4, 5 with m = 2 and K = 3 with
  m = 3, every level; exactly at K = 2, m = 2 (1/3, 2/3, 1 at l = 1, 2, 3) and K = 3, m = 2, l = 1
  (1/7 against 1/4);
- model A's worst case is the symmetric (2k - 1)/m for every m <= 6 solved, now a theorem
  (jackknifePlus_orderReading); model B reaches 1/2 at k = 1 for m = 3..6, the parity construction
  that conformal-checks.py verifies.
Run: python -I -X utf8 conformal-searches.py [--milp]
"""
import itertools
import math
import random
import sys
from fractions import Fraction as Q


def unordered_partitions(items, m):
    if not items:
        yield ()
        return
    first, rest = items[0], items[1:]
    for others in itertools.combinations(rest, m - 1):
        fold = (first,) + others
        remaining = [u for u in rest if u not in others]
        for tail in unordered_partitions(remaining, m):
            yield (fold,) + tail


def anneal_bag(K, m, l, steps, seed):
    rng = random.Random(seed)
    n = K * m
    pts = list(range(n + 1))
    sets = list(itertools.combinations(pts, m + 1))
    sidx = {L: j for j, L in enumerate(sets)}
    configs = []
    for x in pts:
        others = [u for u in pts if u != x]
        for part in unordered_partitions(others, m):
            configs.append([(sidx[tuple(sorted(f + (x,)))], x) for f in part])
    member = [[] for _ in sets]
    for c, parts in enumerate(configs):
        for (j, _x) in parts:
            member[j].append(c)
    d = []
    for L in sets:
        vals = list(range(m + 1))
        rng.shuffle(vals)
        d.append(dict(zip(L, vals)))
    tot = [sum(d[j][x] for j, x in parts) for parts in configs]
    good = sum(1 for v in tot if v <= l - 1)
    best = good
    for s in range(steps):
        T = 2.0 * (1 - s / steps) + 0.05
        j = rng.randrange(len(sets))
        a, b = rng.sample(sets[j], 2)
        old = d[j]
        new = dict(old)
        new[a], new[b] = old[b], old[a]
        delta, changes = 0, []
        for c in member[j]:
            after = tot[c] + sum(new[x] - old[x] for jj, x in configs[c] if jj == j)
            changes.append((c, after))
            delta += (after <= l - 1) - (tot[c] <= l - 1)
        if delta >= 0 or rng.random() < math.exp(delta / T):
            d[j] = new
            for c, after in changes:
                tot[c] = after
            good += delta
            best = max(best, good)
    return Q(best, len(configs))


def cv_b(l, m, n):
    return Q(max(2 * l + m - 2, 0), n + m)


def main():
    for K, m in ((3, 2), (4, 2), (3, 3), (5, 2)):
        n = K * m
        for l in range(1, n + 1):
            if cv_b(l, m, n) >= 1:
                continue
            steps = 40000 if n <= 8 else 15000
            frac = max(anneal_bag(K, m, l, steps, seed) for seed in range(3))
            print(f"fixed bag, K = {K}, m = {m}, l = {l}: worst found {frac} = {float(frac):.4f}; "
                  f"CV-b {float(cv_b(l, m, n)):.4f}; one fold {float(Q(l, m + 1)):.4f}; "
                  f"{'EXCEEDS CV-b' if frac > cv_b(l, m, n) else 'within'}", flush=True)
    if "--milp" in sys.argv:
        milp_checks()


def milp_checks():
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds

    def jack(m, k, model):
        t = m - 1
        perms = list(itertools.permutations(range(m)))
        ctx, rows = {}, []
        for s in perms:
            r = []
            for i in range(t):
                tup = tuple(s[j] for j in range(t) if j != i)
                key = tup if model == "A" else (i, tup)
                ctx.setdefault(key, len(ctx))
                r.append((ctx[key], s[t] < s[i]))
            rows.append(r)
        nd, ng = len(ctx), len(perms)
        A = np.zeros((ng, nd + ng))
        ub = np.zeros(ng)
        for p, r in enumerate(rows):
            A[p, nd + p] = m - k
            for c, tmin in r:
                A[p, c] += -1 if tmin else 1
                ub[p] += 0 if tmin else 1
        res = milp(np.concatenate([np.zeros(nd), -np.ones(ng)]), constraints=LinearConstraint(A, -np.inf, ub),
                   integrality=np.ones(nd + ng), bounds=Bounds(0, 1))
        return Q(int(round(-res.fun)), ng)

    for model in ("A", "B"):
        for m, k in ((3, 1), (4, 1), (4, 2), (5, 1)):
            print(f"jackknife+, model {model}, m = {m}, k = {k}: worst {jack(m, k, model)}; "
                  f"symmetric {Q(2 * k - 1, m)}", flush=True)


if __name__ == "__main__":
    main()
