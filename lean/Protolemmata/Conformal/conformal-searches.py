"""Searches behind slice 9's open items, beside conformal-checks.py's exact assertions.

1. CV+'s tournament bound (CV-b) on a fixed bag of n + 1 units, without the further test units.
   For each left-out set L (a fold with the test, m + 1 units), a symmetric fit ranks the residuals;
   d_L(x) counts the units of L ranked strictly above x, a bijection L -> {0..m} for strict orders
   (ties only help the test). A uniform ordering makes the test x uniform and the folds a uniform
   partition of the rest, and CV+ refuses iff sum_k d_{F_k + x}(x) <= l - 1. The search maximizes
   the refused share over the rankings by annealing; each figure is an exact count for the rankings
   found, so a lower bound on the worst case. Standard library only.
2. Exact fixed-bag checks at K = 2, m = 2 and at K = 3, m = 2, l = 1 use the standard library.
   The optional scipy MILP checks the jackknife+ for fits that read their units' order (model A,
   positions closed up) or know their units' positions (model B), m <= 5.

Results of 2026-11-09 (recorded in conformal-lemmas.md, slice 9):
- the fixed bag stays within CV-b's bound everywhere searched: K = 3, 4, 5 with m = 2 and K = 3 with
  m = 3, every level; exactly at K = 2, m = 2 (1/3, 2/3, 1 at l = 1, 2, 3) and K = 3, m = 2, l = 1
  (1/7 against 1/4);
- model A's worst case is the symmetric (2k - 1)/m for every m <= 6 solved, now a theorem
  (jackknifePlus_orderReading); model B reaches 1/2 at k = 1 for m = 3..6, the parity construction
  that conformal-checks.py verifies.
3. Slice 14's follow-ups (2026-11-13): weighted K-fold CV+ under likelihood-ratio weights, by merging the
   folds. Each fold's symmetric fit ranks the fold with the test as any function of the other folds; the
   test is u with chance w_u / W given the bag. Annealing over the rankings maximizes the exact miss rate.
   Result of 2026-11-13: nothing past slice 9's unweighted merged bound or past E[2c] in 24 configurations
   (K <= 4, n + 1 <= 7, weights 1..20). The largest shares found were 0.326 and 0.292. Later that day the
   merged bound was proved with the test's weight replaced by a cap read from the bag
   (cvPlus_merged_likelihoodRatio and its refined form, follow-up 2e in conformal-lemmas.md). The rerun
   asserts it in every configuration: the largest share found is 0.322, and it is below E[2c] in 17 of the 24.
   Whether E[2c] itself bounds the miss rate stays open.
Run: python -I -X utf8 conformal-searches.py [--milp]
Quick exact fixed-bag checks only: python -I -X utf8 conformal-searches.py --exact-bags
Likelihood-ratio merged search only: python -I -X utf8 conformal-searches.py --weights
"""
import itertools
import math
import random
import sys
from collections import Counter
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


def exact_bags():
    # Each string ranks the members of one triple, in lexicographic triple/member order.
    # The one-fold upper bound is l/3. A witness attaining it certifies the optimum.
    witnesses = (
        "012 102 120 120 102 012 102 012 102 120",
        "012 021 210 210 102 120 012 201 012 201",
        "012 021 120 102 210 201 012 021 102 012",
    )
    units = range(5)
    for level, witness in enumerate(witnesses, start=1):
        ranks = {}
        for triple, order in zip(itertools.combinations(units, 3), witness.split(), strict=True):
            assert sorted(order) == list("012")
            ranks[triple] = dict(zip(triple, map(int, order)))
        refused = total = 0
        for test in units:
            for folds in unordered_partitions([i for i in units if i != test], 2):
                rank_sum = sum(ranks[tuple(sorted((*fold, test)))][test] for fold in folds)
                refused += rank_sum < level
                total += 1
        assert total == 15 and Q(refused, total) == Q(level, 3)
        print(f"two folds of two, level {level}: exact worst {refused}/{total}")

    # The translates of {0, 1, 3} modulo 7 are the lines of the Fano plane.
    # For a candidate test unit, its three lines pair up the remaining six units.
    units = range(7)
    lines = [(i, (i + 1) % 7, (i + 3) % 7) for i in units]
    planes = {tuple(sorted(tuple(sorted(labels[i] for i in line)) for line in lines))
              for labels in itertools.permutations(units)}
    counts = Counter()
    for plane in planes:
        # At level one a refused unit must win on each of its lines. Any two units
        # share a line, so at most one can do so. Ranking the least unit highest
        # on each line attains one refusal per plane.
        for x, y in itertools.combinations(units, 2):
            assert sum(x in line and y in line for line in plane) == 1
        winners = 0
        for test in units:
            incident = [line for line in plane if test in line]
            folds = tuple(sorted(tuple(i for i in line if i != test) for line in incident))
            assert len(folds) == 3
            assert sorted(itertools.chain.from_iterable(folds)) == [i for i in units if i != test]
            counts[test, folds] += 1
            winners += all(min(line) == test for line in incident)
        assert winners == 1
    # Equal multiplicities give the uniform law over tests and partitions.
    assert len(planes) == 30 and len(counts) == 105 and set(counts.values()) == {2}
    print("three folds of two, level 1: exact worst 1/7 (30 planes, 105 configurations)")


def labelled_partitions(rest, K, m):
    """Labelled partitions of rest into K folds of m."""
    if K == 0:
        yield []
        return
    for comb in itertools.combinations(rest, m):
        left = [x for x in rest if x not in comb]
        for tail in labelled_partitions(left, K - 1, m):
            yield [list(comb)] + tail


def lr_merged_anneal(K, m, w, alpha, steps, seed):
    """Weighted K-fold CV+ under likelihood-ratio weights, the abstract model: fold k's symmetric fit ranks
    F_k with the test as any function of O_k, the units outside them; given the bag the test is u with chance
    w_u / W and the folds are uniform. Returns the largest exact miss rate found and the candidate bounds."""
    rr = random.Random(seed)
    N = K * m + 1
    W = sum(w)
    arr = [(u, parts) for u in range(N)
           for parts in labelled_partitions([x for x in range(N) if x != u], K, m)]
    Os = sorted({frozenset(x for x in range(N) if x != u and x not in F) for u, parts in arr for F in parts},
                key=sorted)
    rank = {}
    for O in Os:
        comp = [x for x in range(N) if x not in O]
        rr.shuffle(comp)
        rank[O] = {x: i for i, x in enumerate(comp)}

    def miss():
        per = {}
        for u, parts in arr:
            beaten = sum(w[i] for F in parts
                         for i in F if rank[frozenset(x for x in range(N) if x != u and x not in F)][i]
                         < rank[frozenset(x for x in range(N) if x != u and x not in F)][u])
            a, b = per.get(u, (0, 0))
            per[u] = (a + (beaten >= (1 - alpha) * W), b + 1)
        return sum(Q(w[u], W) * Q(a, b) for u, (a, b) in per.items())

    cur = best = miss()
    for s in range(steps):
        T = 0.05 * (1 - s / steps) + 1e-4
        O = rr.choice(Os)
        old = dict(rank[O])
        a, b = rr.sample(list(old), 2)
        rank[O][a], rank[O][b] = old[b], old[a]
        new = miss()
        if new >= cur or rr.random() < math.exp(float(new - cur) / T):
            cur = new
            best = max(best, cur)
        else:
            rank[O] = old
    n = N - 1
    b0 = Q(2) * (alpha * (n + 1) + K - 1) / (n + K)
    e2c = sum(Q(w[u], W) * 2 * (alpha * W + (K - 1) * w[u]) / (W + (K - 1) * w[u]) for u in range(N))
    return best, b0, e2c, lr_merged_proven(K, w, alpha)


def lr_merged_proven(K, w, alpha):
    """The proven merged bound (cvPlus_merged_likelihoodRatio and its refined form): the smaller of
    min(1, 2c(w_max)) and min(1, 2L(W + (K - 1)e)/N^2), gamma the largest weight at most alpha W."""
    W = sum(w)
    gm = max(w)
    simple = min(Q(1), 2 * (alpha * W + (K - 1) * gm) / (W + (K - 1) * gm))
    small = [x for x in w if x <= alpha * W]
    if not small:
        return Q(0)
    g = max(small)
    L, N = alpha * W + (K - 1) * g, W + (K - 1) * g
    e = Q(sum(x * x for x in w), W)
    refined = min(Q(1), 2 * L * (W + (K - 1) * e) / N ** 2) if 2 * L <= N else Q(1)
    return min(simple, refined)


def weighted_search():
    """Slice 14's follow-ups: merged K-fold CV+ under likelihood-ratio weights. Searches for a case past
    slice 9's unweighted merged bound B0, past E[2c] (c from the test's weight, unproven), and checks the
    proven bound, reporting how much of it the worst case found uses."""
    rr = random.Random(11)
    worst_b0 = worst_e2c = worst_pr = 0.0
    for trial in range(24):
        K, m = rr.choice([(2, 1), (2, 2), (3, 1), (2, 3), (4, 1), (3, 2)])
        w = [rr.choice([1, 1, 2, 3, 5, 8, 20]) for _ in range(K * m + 1)]
        alpha = Q(rr.randint(1, 8), 20)
        best, b0, e2c, pr = lr_merged_anneal(K, m, w, alpha, 1500 if K * m + 1 <= 5 else 600,
                                             rr.randint(0, 10 ** 6))
        assert best <= pr, (K, m, w, alpha, best, pr)
        worst_b0 = max(worst_b0, float(best / b0))
        worst_e2c = max(worst_e2c, float(best / e2c))
        worst_pr = max(worst_pr, float(best / pr)) if pr else worst_pr
        print(f"K = {K}, m = {m}, alpha = {alpha}, w = {w}: worst found {float(best):.4f}; B0 {float(b0):.4f}; "
              f"E[2c] {float(e2c):.4f}; proven {float(pr):.4f}{'  EXCEEDS' if best > min(b0, e2c) else ''}",
              flush=True)
    print(f"largest share found: of B0 {worst_b0:.3f}, of E[2c] {worst_e2c:.3f}, of the proven bound {worst_pr:.3f}")


def main():
    if "--weights" in sys.argv:
        weighted_search()
        return
    exact_bags()
    if "--exact-bags" in sys.argv:
        return
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
