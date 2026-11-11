"""Numeric checks of conformal-lemmas.md (CB-1 to CB-4, M-1 and M-2, slice 9), labelled by item, before any Lean.

Exhaustive enumeration over every ordering at small n, in exact rational arithmetic, so each check
is an assertion and a printed figure, not an estimate. Standard library only.

Conventions, fixed in CB-1 and used throughout:
- a score is a conformity score: higher conforms better (a log or plain predictive density);
- rank_le(s, u) counts the units scoring at most s[u], u itself included;
- the accept rule at level k accepts u when k < rank_le(s, u), Fong and Holmes's pi(y) > k/(n+1);
- the strict rule accepts u when k <= rank_lt(s, u), the units scoring strictly below; it is
  Hoff's "c_{n+1} above the k-th order statistic", which he states without ties.
"""
import itertools
import math
import random
from fractions import Fraction as Q

rng = random.Random(20261106)


def rank_le(s, u, C):
    return sum(1 for j in C if s[j] <= s[u])


def rank_lt(s, u, C):
    return sum(1 for j in C if s[j] < s[u])


def rank_shift(s, u, C, d):
    """Units of C, u itself always counted, scoring at least d below s[u]."""
    return sum(1 for j in C if j == u or s[j] + d <= s[u])


def accepts(k, s, u, C):
    return k < rank_le(s, u, C)


def accepts_strict(k, s, u, C):
    return k <= rank_lt(s, u, C)


def perms_on(N, C):
    """The permutations of range(N) that move only units of C, as tuples p with p[i] = sigma(i)."""
    C = list(C)
    out = []
    for img in itertools.permutations(C):
        p = list(range(N))
        for a, b in zip(C, img):
            p[a] = b
        out.append(tuple(p))
    return out


def compose(z, p):
    """z o sigma: the unit at position i of the reordered data is z[sigma(i)]."""
    return tuple(z[p[i]] for i in range(len(z)))


def coverage_orderings(z, score, k, t, C, rule=accepts):
    G = perms_on(len(z), C)
    hits = sum(1 for p in G if rule(k, score(compose(z, p)), t, C))
    return Q(hits, len(G))


def is_equivariant(z, score, C):
    s = score(z)
    return all(score(compose(z, p))[i] == s[p[i]] for p in perms_on(len(z), C) for i in range(len(z)))


def injective_on(s, C):
    return len({s[j] for j in C}) == len(list(C))


def multisets(alphabet, N):
    return list(itertools.combinations_with_replacement(alphabet, N))


# Scores used across the checks. Each is a function of the whole data tuple, symmetric in it.

def score_sq(z):
    """Minus the squared distance to the mean of all units, the test's included."""
    m = sum(Q(x) for x in z) / len(z)
    return tuple(-(Q(x) - m) ** 2 for x in z)


def score_postpred(z, K=3, a=Q(1)):
    """In-sample Dirichlet-categorical posterior predictive (Fong and Holmes's sigma_i)."""
    return tuple((z.count(x) + a) / (len(z) + a * K) for x in z)


def score_loo(z, K=3, a=Q(1)):
    """Hoff's c_i: unit i's predictive given every other unit, the candidate included."""
    return tuple((z.count(x) - 1 + a) / (len(z) - 1 + a * K) for x in z)


def score_indexed(z):
    """Not equivariant: the score reads the position, as an asymmetric fit would, and marks the
    last unit (the test) down."""
    return tuple(Q(x + 1, i + 1) for i, x in enumerate(z))


print("=== CB-1: the rank guarantee ===")

# CB-1, the counting core: whatever the scores and ties, at most k units of C have rank_le <= k;
# at least k have rank_lt < k; with no ties both are exactly k; and liberal acceptance contains
# strict acceptance, the difference being the ties straddling the threshold.
cases = 0
for N in range(1, 6):
    C = range(N)
    for s in itertools.product(range(3), repeat=N):
        for k in range(N + 1):
            rej = [u for u in C if rank_le(s, u, C) <= k]
            rej_strict = [u for u in C if rank_lt(s, u, C) < k]
            acc = [u for u in C if accepts(k, s, u, C)]
            acc_strict = [u for u in C if accepts_strict(k, s, u, C)]
            straddle = [u for u in C if rank_lt(s, u, C) < k < rank_le(s, u, C)]
            assert len(rej) <= k and len(rej_strict) >= k
            assert set(acc_strict) <= set(acc)
            assert len(acc) == len(acc_strict) + len(straddle) and len(acc_strict) <= N - k
            cases += 1
    for s in itertools.permutations(range(N)):
        for k in range(N + 1):
            assert sum(1 for u in C if rank_le(s, u, C) <= k) == k
            assert all(accepts(k, s, u, C) == accepts_strict(k, s, u, C) for u in C)
        assert sorted(rank_le(s, u, C) for u in C) == list(range(1, N + 1))  # the rank is a bijection
print(f"CB-1 counting core: {cases} (scores, k) cases with ties up to N = 5, all bounds hold; "
      "with no ties the rank is a bijection onto 1..N")

# CB-1, the finite form: over every ordering of a fixed multiset, an equivariant score's test unit
# is accepted at rate >= 1 - k/N (liberal) and <= 1 - k/N (strict), with equality when the scores
# are distinct. Random tie-breaking, as distinct labels carried with the units, makes it exact.
worst_lib, worst_strict, n_ms = Q(1), Q(0), 0
for score in (score_sq, score_postpred, score_loo):
    for N in range(2, 6):
        C = range(N)
        t = N - 1
        for z in multisets(range(3), N):
            assert is_equivariant(z, score, C)
            inj = injective_on(score(z), C)
            labels = tuple(range(N))
            zl = tuple(zip(z, labels))

            def tie_broken(w, score=score):
                base = score(tuple(x for x, _ in w))
                return tuple((b, lab) for b, (_, lab) in zip(base, w))

            for k in range(N + 1):
                lib = coverage_orderings(z, score, k, t, C)
                strict = coverage_orderings(z, score, k, t, C, accepts_strict)
                bound = 1 - Q(k, N)
                assert strict <= bound <= lib
                if inj:
                    assert lib == strict == bound
                assert coverage_orderings(zl, tie_broken, k, t, C) == bound
                worst_lib = min(worst_lib, lib - bound)
                worst_strict = max(worst_strict, strict - bound)
            n_ms += 1
print(f"CB-1 finite form: {n_ms} (score, multiset) pairs, every k: liberal - bound >= {worst_lib}, "
      f"strict - bound <= {worst_strict}; distinct scores and tie-broken scores give the bound exactly")

# The hypothesis does the work: a score that reads positions breaks the guarantee.
viol = None
for N in range(2, 5):
    C = range(N)
    for z in multisets(range(3), N):
        for k in range(1, N + 1):
            cov = coverage_orderings(z, score_indexed, k, N - 1, C)
            if cov < 1 - Q(k, N) and (viol is None or cov - (1 - Q(k, N)) < viol[0]):
                viol = (cov - (1 - Q(k, N)), z, k)
print(f"CB-1 without equivariance: coverage falls {-viol[0]} below the bound at multiset {viol[1]}, k = {viol[2]}")


def iid_law(p, N):
    out = {}
    for z in itertools.product(range(len(p)), repeat=N):
        w = Q(1)
        for x in z:
            w *= p[x]
        out[z] = w
    return out


def polya_law(K, N, a=Q(1)):
    """A Polya urn: exchangeable, not independent."""
    out = {}
    for z in itertools.product(range(K), repeat=N):
        w, counts = Q(1), [0] * K
        for i, x in enumerate(z):
            w *= (counts[x] + a) / (i + a * K)
            counts[x] += 1
        out[z] = w
    return out


def coverage_law(law, score, k, t, C, rule=accepts):
    return sum(w for z, w in law.items() if rule(k, score(z), t, C))


# CB-1 under an exchangeable law: every exchangeable law on a finite space is a mixture of the
# finite form, so the bound holds for independent draws and for a Polya urn alike.
for name, law in (("iid (1/2, 1/3, 1/6)", iid_law([Q(1, 2), Q(1, 3), Q(1, 6)], 4)),
                  ("Polya urn", polya_law(3, 4))):
    assert sum(law.values()) == 1
    row = []
    for score in (score_sq, score_postpred, score_loo):
        for k in range(5):
            cov = coverage_law(law, score, k, 3, range(4))
            assert cov >= 1 - Q(k, 4)
        row.append(str(coverage_law(law, score, 1, 3, range(4))))
    cov_bad = coverage_law(law, score_indexed, 1, 3, range(4))
    print(f"CB-1 law, {name}, N = 4: k = 1 coverage {', '.join(row)} (bound 3/4); indexed score {cov_bad}")


# CB-1 per stratum: units exchangeable within strata, and a score that sees every stratum (partial
# pooling) but is equivariant under permutations inside the test's stratum. Ranking inside the
# stratum covers there; ranking across strata does not.
def score_pooled(z, strata, K=3, c=Q(1, 2), a=Q(1)):
    out = []
    for i, x in enumerate(z):
        own = sum(1 for j, y in enumerate(z) if strata[j] == strata[i] and y == x)
        other = sum(1 for j, y in enumerate(z) if strata[j] != strata[i] and y == x)
        n_own = sum(1 for j in range(len(z)) if strata[j] == strata[i])
        out.append((own + c * other + a) / (n_own + c * (len(z) - n_own) + a * K))
    return tuple(out)


strata = (0, 0, 0, 1, 1)
pmfs = ([Q(3, 4), Q(1, 8), Q(1, 8)], [Q(1, 8), Q(1, 8), Q(3, 4)])
law = {}
for z in itertools.product(range(3), repeat=5):
    w = Q(1)
    for i, x in enumerate(z):
        w *= pmfs[strata[i]][x]
    law[z] = w
for t in (2, 4):
    Cs = [j for j in range(5) if strata[j] == strata[t]]
    for sname, sc in (("partially pooled", lambda z: score_pooled(z, strata)), ("fully pooled", score_postpred)):
        for z in itertools.product(range(3), repeat=5):
            assert is_equivariant(z, sc, Cs)
        for k in range(len(Cs) + 1):
            assert coverage_law(law, sc, k, t, Cs) >= 1 - Q(k, len(Cs))
        k, k5 = 1, 2
        within, across = coverage_law(law, sc, k, t, Cs), coverage_law(law, sc, k5, t, range(5))
        print(f"CB-1 per stratum, {sname} score, test in a stratum of {len(Cs)}: ranked within it, k = {k}, "
              f"coverage {float(within):.3f} (bound {float(1 - Q(k, len(Cs))):.3f}); ranked across all 5 "
              f"units, k = {k5}, {float(across):.3f} (bound {float(1 - Q(k5, 5)):.3f})")
# Ranked across strata, a stratum the score finds unusual undercovers.
across = coverage_law(law, lambda z: score_pooled(z, strata), 2, 4, range(5))
assert across < 1 - Q(2, 5)

# Small strata: a stratum of m units (m - 1 replicates and the test) attains only the levels
# 1 - k/m. At alpha = 0.1, k = floor(alpha m) is 0, the set is everything, until m >= 10.
print("CB-1 small strata, alpha = 0.1: m -> (k, exact coverage with distinct scores): " + ", ".join(
    f"{m} -> ({int(Q(1, 10) * m)}, {1 - Q(int(Q(1, 10) * m), m)})" for m in (2, 5, 9, 10, 11, 19, 20, 21)))

print()
print("=== CB-2: the split form holds for any fit ===")


# The calibration and test units are scored by a function fitted on the training split. Whatever
# the fit does with the training data, and with its own randomness, coverage holds given the
# training split, under any exchangeable law of all the units.
def split_coverage(law, Ttr, fit, k, seeds=(0,), conditional=False):
    """law over tuples of Ttr training units then n + 1 calibration/test units; the test is last."""
    N = len(next(iter(law))) - Ttr
    C = range(N)
    by_train = {}
    for z, w in law.items():
        tr, ct = z[:Ttr], z[Ttr:]
        for r in seeds:
            fitted = fit(tr, r)
            s = tuple(fitted(x) for x in ct)
            hit = w / len(seeds) if accepts(k, s, N - 1, C) else Q(0)
            acc = by_train.setdefault(tr, [Q(0), Q(0)])
            acc[0] += hit
            acc[1] += w / len(seeds)
    total = sum(v[0] for v in by_train.values())
    worst_cond = min(v[0] / v[1] for v in by_train.values() if v[1] > 0)
    return total, worst_cond


fits = {
    "centre at the training mean": lambda tr, r: (lambda x, m=sum(Q(v) for v in tr) / len(tr): -abs(Q(x) - m)),
    "an unmixed sampler (first training unit)": lambda tr, r: (lambda x, c=tr[0]: -abs(x - c)),
    "a fixed wrong centre": lambda tr, r: (lambda x: -abs(x - 7)),
    "a randomized fit (seeded lookup)": lambda tr, r: (
        lambda x, tbl=random.Random(hash((tr, r))).sample(range(10), 3): Q(tbl[x])),
}
for lname, law in (("iid", iid_law([Q(1, 2), Q(1, 3), Q(1, 6)], 5)), ("Polya urn", polya_law(3, 5))):
    for fname, fit in fits.items():
        for k in range(4):
            tot, cond = split_coverage(law, 2, fit, k, seeds=(0, 1, 2))
            assert tot >= 1 - Q(k, 3) and cond >= 1 - Q(k, 3)
        tot, cond = split_coverage(law, 2, fit, 1, seeds=(0, 1, 2))
        print(f"CB-2 {lname}, 2 training, 2 calibration + test, {fname}: k = 1 coverage {tot}, "
              f"worst given the training split {cond} (bound 2/3)")


# CQR: q_lo, q_hi fitted on training (here deliberately poor), conformity min(y - q_lo, q_hi - y),
# the negative of Romano et al.'s max(q_lo - y, y - q_hi). The accepted y form one interval,
# [q_lo + c_(k), q_hi - c_(k)] with c_(k) the k-th smallest calibration conformity.
def cqr_score(qlo, qhi):
    return lambda p: min(p[1] - qlo[p[0]], qhi[p[0]] - p[1])


xy = [(x, y) for x in (0, 1) for y in range(4)]
pxy = {(0, 0): Q(1, 8), (0, 1): Q(1, 4), (0, 2): Q(1, 16), (0, 3): Q(1, 16),
       (1, 0): Q(1, 16), (1, 1): Q(1, 16), (1, 2): Q(1, 8), (1, 3): Q(1, 4)}
qlo, qhi = {0: Q(1, 2), 1: Q(2)}, {0: Q(1), 1: Q(5, 2)}
s_cqr = cqr_score(qlo, qhi)
N = 4
for k in range(N + 1):
    cov = Q(0)
    for z in itertools.product(xy, repeat=N):
        w = Q(1)
        for p in z:
            w *= pxy[p]
        if accepts(k, tuple(s_cqr(p) for p in z), N - 1, range(N)):
            cov += w
    assert cov >= 1 - Q(k, N)
    if k == 1:
        print(f"CB-2 CQR, N = 4, k = 1: coverage {cov} (bound 3/4), with quantile fits that miss")
grid = [Q(i, 4) for i in range(-8, 21)]
n_sets = 0
for cal in itertools.product(xy, repeat=3):
    cal_s = sorted(s_cqr(p) for p in cal)
    for x in (0, 1):
        for k in range(1, 4):
            accepted = [y for y in grid
                        if accepts(k, tuple(s_cqr(p) for p in cal) + (s_cqr((x, y)),), 3, range(4))]
            c = cal_s[k - 1]
            assert accepted == [y for y in grid if qlo[x] + c <= y <= qhi[x] - c]
            n_sets += 1
print(f"CB-2 CQR: {n_sets} accept sets on a grid in y, each the interval [q_lo + c_(k), q_hi - c_(k)]")

print()
print("=== CB-4: latents, and the score's symmetry ===")

# A finite model: theta in {0, 1}, a latent v in {0, 1} per unit, x in {0, 1}, y in {0, 1, 2}.
TH, V, X, Y = range(2), range(2), range(2), range(3)
prior = [Q(1, 3), Q(2, 3)]
qv = [[Q(1, 2), Q(1, 2)], [Q(3, 4), Q(1, 4)]]
# f[theta][v][x][y]: the latent moves y strongly, so integrating it matters.
f = [[[[Q(6, 10), Q(3, 10), Q(1, 10)], [Q(1, 10), Q(3, 10), Q(6, 10)]],
      [[Q(1, 10), Q(1, 10), Q(8, 10)], [Q(5, 10), Q(4, 10), Q(1, 10)]]],
     [[[Q(8, 10), Q(1, 10), Q(1, 10)], [Q(2, 10), Q(6, 10), Q(2, 10)]],
      [[Q(3, 10), Q(4, 10), Q(3, 10)], [Q(1, 10), Q(1, 10), Q(8, 10)]]]]
for th in TH:
    for v in V:
        for x in X:
            assert sum(f[th][v][x]) == 1
Z = [(x, y) for x in X for y in Y]


def g(th, z):
    return sum(qv[th][v] * f[th][v][z[0]][z[1]] for v in V)


def post(zs):
    u = [prior[th] * prod(g(th, z) for z in zs) for th in TH]
    return [w / sum(u) for w in u]


def prod(it):
    out = Q(1)
    for a in it:
        out *= a
    return out


def joint(zs):
    """The joint posterior of theta and every unit's latent, by enumeration."""
    out = {}
    for th in TH:
        for vs in itertools.product(V, repeat=len(zs)):
            out[(th, vs)] = prior[th] * prod(qv[th][v] * f[th][v][z[0]][z[1]] for v, z in zip(vs, zs))
    tot = sum(out.values())
    return {key: w / tot for key, w in out.items()}


def score_int(zs):
    p = post(zs)
    return tuple(sum(p[th] * g(th, z) for th in TH) for z in zs)


def score_lat(zs):
    J = joint(zs)
    return tuple(sum(w * f[th][vs[i]][zs[i][0]][zs[i][1]] for (th, vs), w in J.items())
                 for i in range(len(zs)))


def score_lat_closed(zs):
    p = post(zs)
    return tuple(sum(p[th] * sum(qv[th][v] * f[th][v][z[0]][z[1]] ** 2 for v in V) / g(th, z) for th in TH)
                 for z in zs)


def score_lat_from_draws(zs):
    """From the n-unit joint posterior (the retained draws, exactly): the add-one-in weight
    g(z_new | theta), and the candidate's latent from p(v | theta, z_new)."""
    zn, znew = zs[:-1], zs[-1]
    Jn = joint(zn)
    wts = {key: w * g(key[0], znew) for key, w in Jn.items()}
    tot = sum(wts.values())
    out = []
    for i in range(len(zn)):
        out.append(sum(w / tot * f[th][vs[i]][zn[i][0]][zn[i][1]] for (th, vs), w in wts.items()))
    cand = sum(w / tot * sum(qv[th][v] * f[th][v][znew[0]][znew[1]] / g(th, znew) * f[th][v][znew[0]][znew[1]]
                             for v in V) for (th, vs), w in wts.items())
    return tuple(out) + (cand,)


def score_mixed(zs):
    """The broken reading: data units by their latents (the n-unit draws, reweighted), the
    candidate by the integrated density."""
    return score_lat_from_draws(zs)[:-1] + (score_int(zs)[-1],)


def score_mixed_unweighted(zs):
    """Worse still: the data units' latents from the n-unit posterior, not reweighted."""
    zn = zs[:-1]
    Jn = joint(zn)
    data = tuple(sum(w * f[th][vs[i]][zn[i][0]][zn[i][1]] for (th, vs), w in Jn.items()) for i in range(len(zn)))
    return data + (score_int(zs)[-1],)


def score_mixed_reverse(zs):
    """The other way round: data units integrated, the candidate by its latent."""
    return score_int(zs)[:-1] + (score_lat_from_draws(zs)[-1],)


# (a) add-one-in: the posterior given n + 1 units is the n-unit posterior times g(z_new | theta),
# and it is the theta-marginal of the joint posterior over every latent.
for zs in itertools.product(Z, repeat=3):
    J = joint(zs)
    marg = [sum(w for (th, _), w in J.items() if th == t0) for t0 in TH]
    p_n = post(zs[:-1])
    aoi = [p_n[th] * g(th, zs[-1]) for th in TH]
    aoi = [w / sum(aoi) for w in aoi]
    assert marg == post(zs) == aoi
print("CB-4(a): over all 216 data sets of 3 units, the joint posterior's theta-marginal, the "
      "integrated posterior and the add-one-in reweighting agree exactly")

# (b), (c): both scores equivariant; (c)'s closed form and its computation from n-unit draws agree;
# (c) >= (b) unit by unit (Jensen), equal only where f does not vary with v.
gap = Q(0)
for zs in itertools.product(Z, repeat=3):
    assert is_equivariant(zs, score_int, range(3))
    sl = score_lat(zs)
    assert sl == score_lat_closed(zs) == score_lat_from_draws(zs)
    si = score_int(zs)
    assert all(a >= b for a, b in zip(sl, si))
    gap = max(gap, max(a - b for a, b in zip(sl, si)))
for zs in multisets(Z, 3):
    assert is_equivariant(zs, score_lat, range(3))
print(f"CB-4(b, c): both scores equivariant; (c) from the joint, its closed form and the n-unit draws "
      f"with the candidate's latent agree exactly; (c) - (b) up to {float(gap):.4f}")


# Coverage under two exchangeable laws: the model's own prior predictive, and an iid law outside it.
def model_law(N):
    out = {}
    for zs in itertools.product(Z, repeat=N):
        out[zs] = sum(prior[th] * prod(Q(1, 2) * g(th, z) for z in zs) for th in TH)
    return out


p_out = {z: Q(1, 6) for z in Z}
p_out[(0, 2)], p_out[(1, 0)] = Q(1, 4), Q(1, 12)


def out_law(N):
    return {zs: prod(p_out[z] for z in zs) for zs in itertools.product(Z, repeat=N)}


N = 4
for lname, law in (("the model's own law", model_law(N)), ("an iid law outside the model", out_law(N))):
    assert sum(law.values()) == 1
    for k in (1, 2):
        bound = 1 - Q(k, N)
        row = {}
        for sname, sc in (("(b)", score_int), ("(c)", score_lat_from_draws), ("(d)", score_mixed),
                          ("(d) unweighted", score_mixed_unweighted), ("reverse", score_mixed_reverse)):
            row[sname] = coverage_law(law, sc, k, N - 1, range(N))
        assert row["(b)"] >= bound and row["(c)"] >= bound
        assert row["(d)"] < bound and row["(d)"] <= row["(c)"]
        print(f"CB-4 {lname}, N = 4, k = {k}, bound {bound} = {float(bound):.3f}: " +
              ", ".join(f"{a} {float(c):.3f}" for a, c in row.items()))

# Hoff's score is each unit's predictive given the rest, c_i = p(z_i | z_-i); Fong and Holmes's is
# the in-sample p(z_i | z). Both are equivariant, so both cover; they need not rank alike.
def score_hoff(zs):
    out = []
    for i, z in enumerate(zs):
        p = post(zs[:i] + zs[i + 1:])
        out.append(sum(p[th] * g(th, z) for th in TH))
    return tuple(out)


disagree, total = 0, 0
for zs in itertools.product(Z, repeat=3):
    a, b = score_int(zs), score_hoff(zs)
    total += 1
    if any((a[i] < a[j]) != (b[i] < b[j]) for i in range(3) for j in range(3)):
        disagree += 1
zs = ((0, 2), (0, 2), (0, 2), (0, 0))  # Hoff's own setting: no covariate, y alone
a, b = score_int(zs), score_hoff(zs)
assert a[3] < a[0] and b[3] > b[0]
law = model_law(4)
cov_h = coverage_law(law, score_hoff, 1, 3, range(4))
cov_i = coverage_law(law, score_int, 1, 3, range(4))
assert cov_h >= Q(3, 4) and cov_i >= Q(3, 4)
print(f"CB-4 Hoff's c_i against the in-sample score: the unit orderings differ on {disagree} of {total} "
      f"data sets of 3 units; with no covariate, y = (2, 2, 2, 0) puts the lone 0 lowest in-sample "
      f"({float(a[3]):.4f} < {float(a[0]):.4f}) and highest by Hoff's ({float(b[3]):.4f} > {float(b[0]):.4f}); "
      f"at N = 4, k = 1 both cover, {float(cov_h):.3f} and {float(cov_i):.3f}")

# (d) the Lean counterexample: one theta, a fair latent bit, y the bit with probability 3/4.
# Integrated: g = 1/2 for either y. Latent: (1/2)(9/16 + 1/16)/(1/2) = 5/8. Two units, k = 1:
# the candidate always scores strictly lowest, so it is never accepted, against a bound of 1/2.
qq, ff = Q(1, 2), lambda v, y: Q(3, 4) if v == y else Q(1, 4)
g0 = sum(qq * ff(v, 0) for v in (0, 1))
h0 = sum(qq * ff(v, 0) ** 2 for v in (0, 1)) / g0
assert g0 == Q(1, 2) and h0 == Q(5, 8)
for y0, y1 in itertools.product((0, 1), repeat=2):
    s = (h0, g0)  # data unit, then candidate; both values are free of y by symmetry
    assert not accepts(1, s, 1, range(2))
print(f"CB-4(d) Lean counterexample: data unit scores {h0}, candidate {g0}; coverage 0 at k = 1 against 1/2")

print()
print("=== CB-3: what an imperfect score costs ===")


# Each unit j carries its own error bound eps[j] (its Monte Carlo error, say). With every eps[j]
# equal to e, the shifts below are 2e (eps[j] + eps[u]) and 4e (2 eps[j] + 2 eps[u]).
def kth(s, C, r):
    """The r-th smallest score over C (1-based)."""
    return sorted(s[j] for j in C)[r - 1]


def rank_shift_at(s, u, C, d):
    """Units of C, u itself always counted, scoring at least d[j] below s[u]."""
    return sum(1 for j in C if j == u or s[j] + d[j] <= s[u])


def fragile(k, s, u, C, eps):
    return k < rank_le(s, u, C) and rank_shift_at(s, u, C, [e + eps[u] for e in eps]) <= k


def accepts_inflated(k, sh, u, C, eps):
    return k < sum(1 for j in C if sh[j] - eps[j] <= sh[u] + eps[u])


def fragile_est(k, sh, u, C, eps):
    return accepts_inflated(k, sh, u, C, eps) and \
        rank_shift_at(sh, u, C, [2 * e + 2 * eps[u] for e in eps]) <= k


def check_deterministic(s, sh, eps, C):
    N = len(C)
    for k in range(N + 1):
        for u in C:
            if not accepts(k, sh, u, C):
                assert (not accepts(k, s, u, C)) or fragile(k, s, u, C, eps)  # loss only through M
            if accepts(k, s, u, C):
                assert accepts_inflated(k, sh, u, C, eps)  # the inflated rule keeps it
            if fragile(k, s, u, C, eps):
                assert fragile_est(k, sh, u, C, eps)  # M is bounded from the estimated scores
            if accepts(k, sh, u, C) and not accepts(k, s, u, C):
                assert k < sum(1 for j in C if s[j] <= s[u] + eps[j] + eps[u])  # gains, the other side
    if len(set(eps)) == 1:
        for r in range(1, N + 1):
            assert abs(kth(sh, C, r) - kth(s, C, r)) <= eps[0]


n_det = 0
grid = [Q(i, 2) for i in range(5)]
for N in range(1, 5):
    C = list(range(N))
    e0 = Q(1, 2)
    for s in itertools.product(grid, repeat=N):
        for e in itertools.product((-e0, Q(0), e0), repeat=N):
            check_deterministic(s, tuple(a + b for a, b in zip(s, e)), (e0,) * N, C)
            n_det += 1
for N in range(1, 4):
    C = list(range(N))
    for s in itertools.product(grid, repeat=N):
        for eps in itertools.product((Q(0), Q(1, 2), Q(1)), repeat=N):
            for sign in itertools.product((-1, 0, 1), repeat=N):
                check_deterministic(s, tuple(a + g * b for a, g, b in zip(s, sign, eps)), eps, C)
                n_det += 1
for _ in range(3000):
    N = rng.randint(2, 7)
    C = list(range(N))
    uniform = rng.random() < 0.5
    e0 = Q(rng.randint(1, 40), rng.randint(1, 40))
    eps = tuple(e0 if uniform else Q(rng.randint(0, 40), rng.randint(1, 40)) for _ in C)
    s = tuple(Q(rng.randint(-60, 60), rng.randint(1, 12)) for _ in C)
    sh = tuple(a + Q(rng.randint(-1000, 1000), 1000) * e for a, e in zip(s, eps))
    check_deterministic(s, sh, eps, C)
    n_det += 1
print(f"CB-3 deterministic part: {n_det} (exact scores, estimates within per-unit eps) cases, "
      "exhaustive on grids (uniform eps to N = 4, per-unit eps to N = 3) and a seeded random sweep: "
      "a decision is lost only in M, the inflated rule keeps every exact acceptance, M sits inside "
      "the estimated band, and with uniform eps order statistics move at most eps")

# Tightness, and the order-statistic picture: an adversary within the bounds (the test unit down
# by its own eps, the rest up by theirs) refuses exactly the units the bound counts; with uniform
# eps and distinct scores those beyond the k refused are the accepted units below s_(k) + 2 eps.
for _ in range(2000):
    N = rng.randint(2, 7)
    C = list(range(N))
    uniform = rng.random() < 0.5
    e0 = Q(rng.randint(1, 20), rng.randint(1, 20))
    eps = tuple(e0 if uniform else Q(rng.randint(0, 20), rng.randint(1, 20)) for _ in C)
    s = tuple(Q(v, 7) for v in rng.sample(range(-80, 80), N))
    for k in range(1, N):
        for u in C:
            adv = tuple(s[j] - eps[j] if j == u else s[j] + eps[j] for j in C)
            assert (not accepts(k, adv, u, C)) == (rank_shift_at(s, u, C, [e + eps[u] for e in eps]) <= k)
            if uniform:
                in_band = kth(s, C, k + 1) <= s[u] < kth(s, C, k) + 2 * e0
                assert fragile(k, s, u, C, eps) == in_band
print("CB-3 tightness: the adversary attains the bound for every unit, uniform or per-unit eps; with "
      "uniform eps and distinct scores M = {u : s_(k+1) <= s_u < s_(k) + 2 eps}")


# Over orderings, with an asymmetric estimate: the exact score is equivariant, and the estimate for
# each ordering p is anything within that ordering's bounds eps_of(p).
# - The inflated rule keeps >= 1 - k/N whatever the bounds: by ordering and by position.
# - The margin needs bounds that travel with the units (eps0 indexed by the original unit, so
#   position i of ordering p has eps0[p[i]]): then the plain rule covers >= 1 - (k + |M|)/N, the
#   adversary attains it with distinct scores, and every ordering's estimated band counts >= |M|.
def coverage_estimated(z, est, eps_of, k, t, C, inflated=False):
    G = perms_on(len(z), C)
    hits = 0
    for p in G:
        sh = est(compose(z, p), p)
        hits += accepts_inflated(k, sh, t, C, eps_of(p)) if inflated else accepts(k, sh, t, C)
    return Q(hits, len(G))


for N in (3, 4, 5):
    C = list(range(N))
    t = N - 1
    for eps0 in ((Q(1, 50),) * N, tuple(Q(i + 1, 100) for i in range(N))):
        by_unit = lambda p, eps0=eps0: tuple(eps0[p[i]] for i in range(len(p)))
        by_position = lambda p, eps0=eps0: eps0
        for z in multisets(range(3), N):
            s = score_sq(z)
            for eps_of in (by_unit, by_position):

                def adv(w, p, eps_of=eps_of):
                    e = eps_of(p)
                    return tuple(v - e[i] if i == t else v + e[i] for i, v in enumerate(score_sq(w)))

                def noisy(w, p, eps_of=eps_of):
                    e = eps_of(p)
                    return tuple(v + e[i] * Q(random.Random(hash((p, i))).randint(-10, 10), 10)
                                 for i, v in enumerate(score_sq(w)))

                for k in range(N + 1):
                    for est in (adv, noisy):
                        assert coverage_estimated(z, est, eps_of, k, t, C, inflated=True) >= 1 - Q(k, N)
                    if eps_of is by_position:
                        continue
                    M = sum(1 for u in C if fragile(k, s, u, C, eps0))
                    rej = sum(1 for u in C if not accepts(k, s, u, C))
                    for est in (adv, noisy):
                        assert coverage_estimated(z, est, eps_of, k, t, C) >= 1 - Q(rej + M, N) >= 1 - Q(k + M, N)
                        for p in perms_on(N, C):
                            assert M <= sum(1 for u in C if fragile_est(k, est(compose(z, p), p), u, C, eps_of(p)))
                    if injective_on(s, C):
                        assert coverage_estimated(z, adv, eps_of, k, t, C) == 1 - Q(k + M, N)
print("CB-3 over orderings: the inflated rule >= 1 - k/N with bounds by ordering and position; with "
      "bounds carried by the units, the plain rule >= 1 - (k + |M|)/N, attained by the adversary with "
      "distinct scores, and every ordering's estimated band counts at least |M|")


# A real Monte Carlo case: conformal Bayes by add-one-in importance weights over T draws from the
# n-unit posterior, on a Beta-Bernoulli model with theta on a grid. The draws depend on which units
# are the n, so the estimate is not equivariant. eps is its worst error against the exact score.
TH5 = [Q(i, 6) for i in range(1, 6)]


def bern(th, y):
    return th if y == 1 else 1 - th


def cb_exact(zs):
    u = [prod(bern(th, y) for y in zs) for th in TH5]
    tot = sum(u)
    return tuple(sum(w / tot * bern(th, y) for w, th in zip(u, TH5)) for y in zs)


def cb_draws(zs, T, seed):
    zn = zs[:-1]
    u = [prod(bern(th, y) for y in zn) for th in TH5]
    r = random.Random(hash((zn, seed)))
    draws = r.choices(TH5, weights=[float(w) for w in u], k=T)
    wts = [bern(th, zs[-1]) for th in draws]
    tot = sum(wts)
    return tuple(sum(w / tot * bern(th, y) for w, th in zip(wts, draws)) for y in zs)


N, k = 5, 1
C, t = list(range(N)), N - 1
for T in (20, 200):
    rows = []
    for z in multisets((0, 1), N):
        G = perms_on(N, C)
        ests = {p: cb_draws(compose(z, p), T, 7) for p in G}
        eps = max(abs(a - b) for p in G for a, b in zip(ests[p], cb_exact(compose(z, p))))
        s = cb_exact(z)
        M = sum(1 for u in C if fragile(k, s, u, C, (eps,) * N))
        cov = Q(sum(1 for p in G if accepts(k, ests[p], t, C)), len(G))
        cov_inf = Q(sum(1 for p in G if k < sum(1 for j in C if ests[p][j] <= ests[p][t] + 2 * eps)), len(G))
        exact = coverage_orderings(z, cb_exact, k, t, C)
        assert cov >= 1 - Q(k + M, N) and cov_inf >= 1 - Q(k, N)
        # Per-unit bounds: each unit's worst error over the orderings, carried with the unit.
        eps_u = [max(abs(ests[p][i] - cb_exact(compose(z, p))[i]) for p in G for i in C if p[i] == u) for u in C]
        M_u = sum(1 for u in C if fragile(k, s, u, C, eps_u))
        cov_inf_u = Q(sum(1 for p in G if accepts_inflated(k, ests[p], t, C, [eps_u[p[i]] for i in C])), len(G))
        # Bounds by position: each ordering's own errors, the narrowest valid inflation.
        cov_inf_p = Q(sum(1 for p in G if accepts_inflated(
            k, ests[p], t, C, [abs(a - b) for a, b in zip(ests[p], cb_exact(compose(z, p)))])), len(G))
        assert M_u <= M and cov_inf_u >= 1 - Q(k, N) and cov_inf_p >= 1 - Q(k, N)
        rows.append(f"{''.join(map(str, z))}: eps {float(eps):.3f} (per unit {float(min(eps_u)):.3f} to "
                    f"{float(max(eps_u)):.3f}), |M| {M} (per unit {M_u}), exact {exact}, draws {cov}, "
                    f"inflated {cov_inf} (per unit {cov_inf_u}, by position {cov_inf_p})")
    print(f"CB-3 conformal Bayes by draws, T = {T}, N = 5, k = 1 (bound 4/5): " + "; ".join(rows))

print()
print("=== E35: conformal test martingales (M-1, M-2) ===")

# M-1. The smoothed p-value of the nth unit, against every score so far (conformity scores):
# p_n = (#{i <= n : s_i < s_n} + theta_n #{i <= n : s_i = s_n}) / n. Its finite part is the
# randomized rank R_n = L_n + J_n, with J_n = ceil(theta_n T_n) uniform on {1..T_n}; then
# p_n = (R_n - 1 + V_n)/n with V_n uniform on (0, 1], independent of J_n.


def prefix_scores(seq, n, mode):
    """The scores of the first n units at time n: a fixed A (inductive), or recomputed from the bag."""
    w = tuple(seq[:n])
    return tuple(Q(x) for x in w) if mode == 'inductive' else score_postpred(w)


def rank_law(z, mode, tie='random', stale=False):
    """The exact joint law of (R_1, ..., R_N) over every ordering of z and every tie-break draw."""
    N = len(z)
    law = {}
    for order in itertools.permutations(range(N)):
        seq = [z[o] for o in order]
        kept = []
        dists = []
        for n in range(1, N + 1):
            s = list(prefix_scores(seq, n, mode))
            if stale:  # the full form keeping each score from its own time
                kept.append(s[-1])
                s = kept[:-1] + [s[-1]]
            L = sum(1 for v in s if v < s[-1])
            T = sum(1 for v in s if v == s[-1])
            if tie == 'random':
                dists.append({L + j: Q(1, T) for j in range(1, T + 1)})
            else:  # ties counted for the unit, as CB-1's sets do
                dists.append({L + T: Q(1)})
        for combo in itertools.product(*[list(d.items()) for d in dists]):
            key = tuple(r for r, _ in combo)
            w = Q(1, math.factorial(N))
            for _, pr in combo:
                w *= pr
            law[key] = law.get(key, Q(0)) + w
    return law


def is_uniform_ranks(law, N):
    target = Q(1, math.factorial(N))
    return all(law.get(r, Q(0)) == target for r in itertools.product(*[range(1, n + 1) for n in range(1, N + 1)]))



n_m1 = 0
for N in range(1, 5):
    for z in multisets(range(3), N):
        for mode in ('inductive', 'full'):
            assert is_uniform_ranks(rank_law(z, mode), N)
            n_m1 += 1
print(f"M-1 finite core: {n_m1} (data multiset, score) cases up to N = 4, ties included, inductive and "
      "full (recomputed) scores: the randomized ranks are exactly uniform on {1}x{1,2}x...x{1..N}")

# The premises, broken one at a time.
bad_tie = next(z for z in multisets(range(3), 3) if not is_uniform_ranks(rank_law(z, 'inductive', tie='le'), 3))
bad_stale = next((z for N in (3, 4) for z in multisets(range(3), N)
                  if not is_uniform_ranks(rank_law(z, 'full', stale=True), len(z))), None)
print(f"M-1 premises: ties counted for the unit (no theta) break uniformity at {bad_tie}; the full form "
      f"keeping each score from its own time breaks it at {bad_stale}")
assert bad_stale is not None


# Ranking against a fixed calibration set, as the split form's p-values for many test points do:
# two test points sharing one calibration set have dependent p-values.
def split_pair_law(z_cal_test):
    """Units 0..m-1 calibration, m and m+1 tests; each test ranked against calibration plus itself."""
    N = len(z_cal_test)
    m = N - 2
    law = {}
    for order in itertools.permutations(range(N)):
        seq = [Q(z_cal_test[o]) for o in order]
        cal = seq[:m]
        dists = []
        for t in (seq[m], seq[m + 1]):
            L = sum(1 for v in cal if v < t)
            T = sum(1 for v in cal if v == t) + 1
            dists.append({L + j: Q(1, T) for j in range(1, T + 1)})
        for (r1, a), (r2, b) in itertools.product(dists[0].items(), dists[1].items()):
            law[(r1, r2)] = law.get((r1, r2), Q(0)) + Q(1, math.factorial(N)) * a * b
    return law


law2 = split_pair_law((0, 1, 2, 3, 4))
m1 = {r: sum(v for (a, b), v in law2.items() if a == r) for r in range(1, 5)}
m2 = {r: sum(v for (a, b), v in law2.items() if b == r) for r in range(1, 5)}
assert all(v == Q(1, 4) for v in m1.values()) and all(v == Q(1, 4) for v in m2.values())
assert law2[(1, 1)] != m1[1] * m2[1]
print(f"M-1 against a fixed calibration set (3 calibration, 2 tests): each test's rank is uniform, but "
      f"P(both lowest) = {law2[(1, 1)]} against {m1[1] * m2[1]} if independent")


# The continuous lift, exactly: P(p_1 <= x_1, ..., p_N <= x_N) = x_1 ... x_N on a grid of x, from the
# theta-smoothed p-values directly (each theta's part is a clipped linear function).
def joint_cdf(z, mode, xs):
    N = len(z)
    tot = Q(0)
    for order in itertools.permutations(range(N)):
        seq = [z[o] for o in order]
        pr = Q(1)
        for n in range(1, N + 1):
            s = prefix_scores(seq, n, mode)
            L = sum(1 for v in s if v < s[-1])
            T = sum(1 for v in s if v == s[-1])
            pr *= min(Q(1), max(Q(0), (n * xs[n - 1] - L) / T))
        tot += pr
    return tot / math.factorial(N)


grid_x = [Q(1, 7), Q(2, 5), Q(1, 2), Q(3, 4), Q(1)]
n_lift = 0
for N in (2, 3, 4):
    for z in multisets(range(3), N)[::2]:
        for mode in ('inductive', 'full'):
            for xs in itertools.product(grid_x, repeat=N):
                assert joint_cdf(z, mode, xs) == prod(xs)
                n_lift += 1
print(f"M-1 lifted: {n_lift} exact joint CDF values of the smoothed p-values equal x_1...x_N "
      "(i.i.d. uniform), ties included")

# M-2. The capital S_n = prod_{i<=n} f_i(p_i), each f_i >= 0 with integral 1 on [0, 1], chosen from
# p_1..p_{i-1} only. Exactly on a K-point midpoint grid standing in for [0, 1] (the betting
# functions averaging 1 over it), every path enumerated: E[S_n | history] = S_{n-1}, E[S_N] = 1,
# and P(max_{n <= N} S_n >= c) <= 1/c.
K, N = 4, 4
G = [Q(2 * k - 1, 2 * K) for k in range(1, K + 1)]


def jumper_fn(eps):
    return lambda p: 1 + eps * (p - Q(1, 2))


def simple_jumper(path, J=Q(1, 100)):
    """Algorithm 1 of Vovk et al. 2021 on a path of p-values: capital after each step."""
    C = {e: Q(1, 3) for e in (-1, 0, 1)}
    out = []
    for p in path:
        tot = sum(C.values())
        C = {e: (1 - J) * C[e] + J / 3 * tot for e in C}
        C = {e: C[e] * jumper_fn(e)(p) for e in C}
        out.append(sum(C.values()))
    return out


def power_grid(eps):
    """The power betting function eps p^(eps-1), renormalized to average 1 on the grid."""
    raw = [eps * p ** (eps - 1) for p in G]
    m = sum(raw) / K
    return lambda p: eps * p ** (eps - 1) / m


def product_capital(path, fs):
    out, S = [], Q(1)
    for i, p in enumerate(path):
        S *= fs(path[:i])(p)
        out.append(S)
    return out


strategies = {
    "Simple Jumper": simple_jumper,
    "power, eps = 2 (grid-normalized)": lambda path: product_capital(path, lambda h: power_grid(2)),
    "power, eps = 3, switching on the past": lambda path: product_capital(
        path, lambda h: power_grid(3) if (h and h[-1] > Q(1, 2)) else jumper_fn(-1)),
}
for name, strat in strategies.items():
    paths = list(itertools.product(G, repeat=N))
    caps = {path: strat(list(path)) for path in paths}
    for n in range(1, N + 1):
        for hist in itertools.product(G, repeat=n - 1):
            prev = Q(1) if n == 1 else caps[hist + (G[0],) * (N - n + 1)][n - 2]
            nxt = sum(caps[hist + (p,) + (G[0],) * (N - n)][n - 1] for p in G) / K
            assert nxt == prev
    assert sum(caps[pth][N - 1] for pth in paths) / len(paths) == 1
    worst = max(sum(1 for pth in paths if max(caps[pth]) >= c) / Q(len(paths)) * c for c in (Q(3, 2), Q(2), Q(3), Q(5)))
    assert worst <= 1
    print(f"M-2 exact on the grid, {name}: E[S_n | history] = S_(n-1) on every history, E[S_4] = 1, "
          f"max over c of c P(max S >= c) = {worst}")

# The premises, broken: an f integrating to 11/10, and an f that reads the current p.
over = lambda path: product_capital(path, lambda h: (lambda p: Q(11, 10)))
look = lambda path: [Q(2) ** (i + 1) for i in range(len(path))]  # bets after seeing p: doubles surely
for name, strat, c in (("integral 11/10", over, Q(5, 4)), ("reads the current p", look, Q(4))):
    paths = list(itertools.product(G, repeat=N))
    freq = sum(1 for pth in paths if max(strat(list(pth))) >= c) / Q(len(paths))
    assert freq > 1 / c
    print(f"M-2 broken, {name}: P(max S >= {c}) = {freq} > 1/c = {1 / c}")

# Continuous: the integrals, and the one-step identity for the Simple Jumper, exactly (it is
# polynomial in p): integral_0^1 f_eps(p) dp = 1 + eps (1/2 - 1/2) = 1.
for eps in (-1, 0, 1):
    integral = 1 + eps * (Q(1, 2) - Q(1, 2))
    assert integral == 1
print("M-2 continuous: each Simple Jumper bet integrates to 1 on [0, 1]; the power bet eps p^(eps-1) "
      "integrates to [p^eps]_0^1 = 1 for every eps in (0, 1]")

# Seeded Monte Carlo: Ville's bound attained, approached, and kept.
mc = random.Random(20261108)
R = 200_000
for c in (2, 10, 100):
    hits = sum(1 for _ in range(R) if mc.random() <= 1 / c)  # all in on p_1 <= 1/c: S_1 = c there
    se = math.sqrt((1 / c) * (1 - 1 / c) / R)
    assert abs(hits / R - 1 / c) < 4 * se
    print(f"M-2 Ville attained, c = {c}: all in on p_1 <= 1/c reaches c with frequency {hits / R:.5f} "
          f"(1/c = {1 / c:.5f}, se {se:.5f})")
for lam in (Q(1, 2), Q(1, 5)):
    c, R2, hits = 10, 4000, 0
    for _ in range(R2):
        S = 1.0
        while 1e-4 < S < c:
            S *= 1 + float(lam) * (1 - 2 * mc.random())
        hits += S >= c
    print(f"M-2 Ville approached, c = 10, bet 1 + {float(lam)}(1 - 2p) until absorbed: frequency "
          f"{hits / R2:.4f} (bound 0.1; overshoot shrinks with the bet)")
    assert hits / R2 <= 0.1 + 4 * math.sqrt(0.1 * 0.9 / R2)
R3, Nmax = 2000, 500
for c in (2, 10):
    hj = hp = 0
    for _ in range(R3):
        path = [mc.random() for _ in range(Nmax)]
        C = {e: 1 / 3 for e in (-1, 0, 1)}
        Sp, peakj, peakp = 1.0, 0.0, 0.0
        for p in path:
            tot = sum(C.values())
            C = {e: 0.99 * C[e] + 0.01 / 3 * tot for e in C}
            C = {e: C[e] * (1 + e * (p - 0.5)) for e in C}
            peakj = max(peakj, sum(C.values()))
            Sp *= 0.5 * max(p, 1e-12) ** -0.5
            peakp = max(peakp, Sp)
        hj += peakj >= c
        hp += peakp >= c
    se = math.sqrt((1 / c) * (1 - 1 / c) / R3)
    assert hj / R3 <= 1 / c + 4 * se and hp / R3 <= 1 / c + 4 * se
    print(f"M-2 Ville kept over {Nmax} uniform p's, c = {c}: Simple Jumper {hj / R3:.4f}, power eps = 1/2 "
          f"{hp / R3:.4f} (bound {1 / c})")


print()
print("=== Slice 9: the jackknife+ and the cross-fitted rule ===")

# Conventions. m = n + 1 units: n training units and the test. pair[i][j] is unit i's conformity
# under the fit that leaves out units i and j. For the jackknife+ it is minus i's absolute residual
# |y_i - mu_{-ij}(x_i)|, and the fit leaving out i and the test is the leave-one-out fit mu_{-i}.
# Unit i is stranger than j when pair[i][j] < pair[j][i], strictly, so a tie makes neither
# stranger. The cross-fitted rule at level k refuses u when u is stranger than at least m - k units.
# The jackknife+ takes k = floor(alpha (n + 1)).
rng9 = random.Random(20261109)


def stranger_count(P, u, C):
    return sum(1 for j in C if P[u][j] < P[j][u])


def refuses_pairwise(k, P, u, C):
    return stranger_count(P, u, C) + k >= len(C)


def pair_bound(k):
    """At most 2k - 1 units are refused at level k, and none at k = 0."""
    return max(2 * k - 1, 0)


def comparisons(m, ties=True):
    """Every outcome of the m(m - 1)/2 comparisons: i stranger, j stranger, or (with ties) neither."""
    pairs = list(itertools.combinations(range(m), 2))
    for choice in itertools.product(range(3 if ties else 2), repeat=len(pairs)):
        P = [[0] * m for _ in range(m)]
        for (i, j), c in zip(pairs, choice):
            if c == 0:
                P[j][i] = 1
            elif c == 1:
                P[i][j] = 1
        yield P


# The counting core, Landau's argument as Barber et al. use it: whatever the comparisons, ties
# included, at most 2k - 1 units are refused at level k. The count is attained by 2k - 1 units in a
# regular tournament, each stranger than every unit outside it, so it cannot be sharpened.
n_cases = 0
for m, ties in ((1, True), (2, True), (3, True), (4, True), (5, True), (6, False)):
    worst = [0] * (m + 1)
    for P in comparisons(m, ties):
        n_cases += 1
        sc = [stranger_count(P, u, range(m)) for u in range(m)]
        for k in range(m + 1):
            r = sum(1 for c in sc if c + k >= m)
            assert r <= pair_bound(k)
            worst[k] = max(worst[k], r)
    assert worst == [min(pair_bound(k), m) for k in range(m + 1)]
print(f"Slice 9, the counting core: {n_cases} comparison outcomes (all, ties included, to m = 5; every "
      f"tournament at m = 6); at most 2k - 1 units refused at level k, and 2k - 1 attained")

INF = math.inf


def jack_ends(vlo, vhi, k):
    """The jackknife+ interval at level k: [k-th smallest lower end, k-th largest upper end]. It is
    the whole line at k = 0 and empty past k = n."""
    if k == 0:
        return -INF, INF
    if k > len(vhi):
        return INF, -INF
    return sorted(vlo)[k - 1], sorted(vhi)[len(vhi) - k]


def smallest_quantile(atoms, level):
    """The smallest q with F(q) >= level, each atom of mass 1/len(atoms)."""
    a = sorted(atoms)
    return min(q for q in a if Q(sum(1 for v in a if v <= q), len(a)) >= level)


def largest_quantile(atoms, level):
    """The largest level-quantile: the largest atom q with P(X < q) <= level."""
    a = sorted(atoms)
    return max(q for q in a if Q(sum(1 for v in a if v < q), len(a)) <= level)


# The rule's three forms agree, ties included:
# - Barber et al.'s, the largest alpha-quantile of the lower ends with an atom at -inf and the
#   smallest (1 - alpha)-quantile of the upper ends with an atom at +inf;
# - the k-th smallest lower end and the k-th largest upper end, k = floor(alpha (n + 1));
# - the counts the Lean states: refused when at least n + 1 - k upper ends lie strictly below y, or
#   at least n + 1 - k lower ends strictly above it.
# The interval is closed: a y at an end is covered.
n_cases = 0
for n in range(1, 8):
    for _ in range(150):
        vlo = [Q(rng9.randint(-6, 6), 2) for _ in range(n)]
        vhi = [Q(rng9.randint(-6, 6), 2) for _ in range(n)]
        ys = sorted(set([Q(j, 4) for j in range(-14, 15)] + vlo + vhi))
        for a in range(1, 11):
            alpha = Q(a, 20)
            k = math.floor(alpha * (n + 1))
            lo, hi = jack_ends(vlo, vhi, k)
            assert hi == smallest_quantile(vhi + [INF], 1 - alpha)
            assert lo == largest_quantile(vlo + [-INF], alpha)
            for y in ys:
                assert (y > hi) == (sum(1 for v in vhi if v < y) + k >= n + 1)
                assert (y < lo) == (sum(1 for v in vlo if v > y) + k >= n + 1)
                n_cases += 1
print(f"Slice 9, the jackknife+'s ends: {n_cases} cases (n to 7, alpha from 1/20 to 1/2, ends tied and "
      f"y at the ends); Barber et al.'s quantiles, the k-th order statistics and the counts agree")


def fit_mean(train):
    ybar = sum(p[1] for p in train) / len(train)
    return lambda x: ybar


def fit_median(train):
    ys, h = sorted(p[1] for p in train), len(train)
    med = ys[h // 2] if h % 2 else (ys[h // 2 - 1] + ys[h // 2]) / 2
    return lambda x: med


def fit_line(train):
    h = len(train)
    xb, yb = sum(p[0] for p in train) / h, sum(p[1] for p in train) / h
    sxx = sum((p[0] - xb) ** 2 for p in train)
    b = sum((p[0] - xb) * (p[1] - yb) for p in train) / sxx if sxx else Q(0)
    return lambda x: yb + b * (x - xb)


def fit_nearest(train):
    def f(x):
        d = min(abs(p[0] - x) for p in train)
        near = [p[1] for p in train if abs(p[0] - x) == d]
        return sum(near) / len(near)
    return f


def fit_size_unstable(n):
    """Symmetric, but its fit on all n training units sits 100 above its fits on n - 1."""
    def fit(train):
        ybar = sum(p[1] for p in train) / len(train)
        return lambda x: ybar + (100 if len(train) >= n else 0)
    return fit


def fit_without(z, fit, out):
    return fit(tuple(z[j] for j in range(len(z)) if j not in out))


def pair_matrix(z, fit):
    m = len(z)
    return [[-abs(z[i][1] - fit_without(z, fit, (i, j))(z[i][0])) if i != j else Q(0) for j in range(m)]
            for i in range(m)]


def jack_refuses(k, z, fit, t):
    """The jackknife+ at level k for the test unit t, by its interval; the count form must agree."""
    m, y = len(z), z[t][1]
    vlo, vhi = [], []
    for i in range(m):
        if i != t:
            f = fit_without(z, fit, (i, t))
            r = abs(z[i][1] - f(z[i][0]))
            vlo.append(f(z[t][0]) - r)
            vhi.append(f(z[t][0]) + r)
    lo, hi = jack_ends(vlo, vhi, k)
    refused = y < lo or y > hi
    assert refused == (sum(1 for v in vhi if v < y) + k >= m or sum(1 for v in vlo if v > y) + k >= m)
    return refused


def jackknife_refuses(k, z, fit, t):
    """The plain jackknife: the fit on all n training units, widened by the k-th largest
    leave-one-out residual."""
    m = len(z)
    e = abs(z[t][1] - fit_without(z, fit, (t,))(z[t][0]))
    rs = [abs(z[i][1] - fit_without(z, fit, (i, t))(z[i][0])) for i in range(m) if i != t]
    return sum(1 for r in rs if r < e) + k >= m


def as_q(pts):
    return tuple((Q(x), Q(y)) for x, y in pts)


# Over every ordering of a fixed data set, with symmetric fits: the leave-two-out matrix is
# equivariant, the interval refuses only a strange test, so the jackknife+ covers at least as often
# as the cross-fitted set, which covers at least 1 - (2k - 1)/m. The plain jackknife, which centres
# on the fit to all n training units, has no such bound: with a fit that changes with its sample size
# it covers nothing at any k >= 1, while the jackknife+ on the same fit holds.
jk_data = [as_q([(0, 0), (1, 1), (2, 1), (3, 4)]),
           as_q([(0, 0), (1, 2), (2, 1), (3, 3), (4, 9)]),
           as_q([(0, 1), (0, 1), (1, 1), (2, 1), (3, 2)]),
           as_q([(0, 0), (1, 3), (2, 0), (3, 3), (4, 0), (5, 3)])]
jk_fits = {"mean": lambda n: fit_mean, "median": lambda n: fit_median, "line": lambda n: fit_line,
           "nearest": lambda n: fit_nearest, "size-unstable": fit_size_unstable}
n_orderings, slack, plus_over_set, tight_cases = 0, None, 0, []
for z in jk_data:
    m, t = len(z), len(z) - 1
    perms = list(itertools.permutations(range(m)))
    for name, make in jk_fits.items():
        fit = make(m - 1)
        P0 = pair_matrix(z, fit)
        hits = {k: [0, 0, 0] for k in range(m + 1)}
        for p in perms:
            zp = compose(z, p)
            Pp = pair_matrix(zp, fit)
            assert all(Pp[i][j] == P0[p[i]][p[j]] for i in range(m) for j in range(m))
            for k in range(m + 1):
                ref_j = jack_refuses(k, zp, fit, t)
                ref_p = refuses_pairwise(k, Pp, t, range(m))
                assert ref_p or not ref_j
                hits[k][0] += not ref_j
                hits[k][1] += not ref_p
                hits[k][2] += not jackknife_refuses(k, zp, fit, t)
            n_orderings += 1
        for k in range(m + 1):
            cj, cp, cjk = (Q(h, len(perms)) for h in hits[k])
            assert cj >= cp >= 1 - Q(pair_bound(k), m)
            plus_over_set += cj > cp
            if 1 <= k and 2 * k - 1 < m:
                d = cp - (1 - Q(pair_bound(k), m))
                slack = d if slack is None else min(slack, d)
                if d == 0:
                    tight_cases.append(f"data set {jk_data.index(z) + 1}, {name}, k = {k}: {cj}")
            if name == "size-unstable" and k >= 1:
                assert cjk == 0
print(f"Slice 9, jackknife+ over orderings: {len(jk_data)} data sets (m = 4 to 6, ties included) x 5 "
      f"symmetric fits, {n_orderings} orderings; every bound holds, least slack {slack} where the bound "
      f"is not 0; the interval covers more than the cross-fitted set in {plus_over_set} (data, fit, k) cases")
print("Slice 9, the sharp bound met by an ordinary fit: " + "; ".join(tight_cases))
print("Slice 9, the plain jackknife with a fit that changes with its sample size: coverage 0 at every "
      "k >= 1; the jackknife+ on the same fit meets its bound")

# Under an exchangeable law, exactly: i.i.d. units on four points, and a Polya urn.
pts4 = as_q([(0, 0), (0, 2), (1, 1), (2, 3)])
for law_name, law in (("i.i.d.", iid_law([Q(1, 2), Q(1, 4), Q(1, 8), Q(1, 8)], 4)), ("Polya", polya_law(4, 4))):
    covs = []
    for name in ("mean", "line", "nearest"):
        fit = jk_fits[name](3)
        for k in range(5):
            cov = sum(w for zi, w in law.items() if not jack_refuses(k, tuple(pts4[a] for a in zi), fit, 3))
            assert cov >= 1 - Q(pair_bound(k), 4)
            if k == 1:
                covs.append(f"{name} {float(cov):.4f}")
    print(f"Slice 9, jackknife+ under the {law_name} law, m = 4, k = 1 (bound 3/4): " + ", ".join(covs))


def lookup_fit(data, R):
    """A symmetric fit made for one data set of distinct points. Given the training points in any
    order, it names the two missing units a, b and predicts y_a - R[a][b] at x_a and y_b - R[b][a] at
    x_b, so the leave-two-out residuals are R's entries, each on the upper side."""
    index = {p: i for i, p in enumerate(data)}

    def fit(train):
        a, b = sorted(set(range(len(data))) - {index[p] for p in train})
        pred = {data[a][0]: data[a][1] - R[a][b], data[b][0]: data[b][1] - R[b][a]}
        return lambda x: pred.get(x, Q(0))
    return fit


def tournament_residuals(m, k):
    """2k - 1 units in a regular tournament, each with the larger residual against every unit outside."""
    s = 2 * k - 1
    R = [[Q(1)] * m for _ in range(m)]
    for i in range(s):
        for j in range(m):
            if j >= s or 1 <= (j - i) % s <= (s - 1) // 2:
                R[i][j], R[j][i] = Q(2), Q(1)
    return R


# The bound is attained by the jackknife+ itself: a symmetric fit whose leave-two-out residuals are
# a regular tournament's covers exactly 1 - (2k - 1)/m over the orderings, below Barber et al.'s
# 1 - 2alpha only by 1/m at alpha = k/m.
attained = []
for m, k in ((3, 1), (3, 2), (4, 1), (4, 2), (5, 1), (5, 2), (5, 3), (6, 2), (6, 3), (7, 3)):
    data = as_q([(i, 0) for i in range(m)])
    fit = lookup_fit(data, tournament_residuals(m, k))
    perms = list(itertools.permutations(range(m)))
    cov = Q(sum(1 for p in perms if not jack_refuses(k, compose(data, p), fit, m - 1)), len(perms))
    assert cov == 1 - Q(pair_bound(k), m)
    attained.append(f"m = {m}, k = {k}: {cov}")
print("Slice 9, the bound attained by a symmetric fit: " + "; ".join(attained))
m, k = 20, 2
data = as_q([(i, 0) for i in range(m)])
fit = lookup_fit(data, tournament_residuals(m, k))
swap_last = lambda u: tuple(m - 1 if i == u else (u if i == m - 1 else i) for i in range(m))
cov20 = Q(sum(1 for u in range(m) if not jack_refuses(k, compose(data, swap_last(u)), fit, m - 1)), m)
assert cov20 == Q(17, 20)
print(f"Slice 9, at alpha = 0.1 and n = 19 (k = 2): the worst case covers {cov20} = {float(cov20)} "
      f"(each unit as the test), against Barber et al.'s 1 - 2alpha = 0.8")

# Ties count for neither unit. Counting a tie as strange breaks the counting core: with every
# residual tied, every unit is refused at k = 1. The closed interval covers y at its ends, here all y.
z_tied = as_q([(0, 1), (1, 1), (2, 1), (3, 1)])
Pt = pair_matrix(z_tied, fit_mean)
tie_strange = sum(1 for u in range(4) if sum(1 for j in range(4) if j != u and Pt[u][j] <= Pt[j][u]) + 1 >= 4)
assert tie_strange == 4 > pair_bound(1)
assert sum(1 for u in range(4) if refuses_pairwise(1, Pt, u, range(4))) == 0
assert not jack_refuses(1, z_tied, fit_mean, 3)
print("Slice 9, ties: all residuals tied, a tie counted as strange refuses all 4 units at k = 1 "
      "(the core allows 1); with ties counting for neither, none is refused, and the closed "
      "interval covers y at its ends")


# A fit that reads its training units' order is outside the hypothesis. Model one exactly: given the
# ordered training tuple it picks which of the two missing units is stranger. Exhaustively at m = 4
# (every such choice, 2^12 fits), none refuses the test on more than (2k - 1)/m of the orderings.
# The Lean statement needs the symmetry; whether the jackknife+'s row alone needs it is left open.
def order_reading_worst(m, k):
    t = m - 1
    perms = list(itertools.permutations(range(m)))
    tuples = sorted({tuple(s[j] for j in range(m) if j not in (i, t)) for s in perms for i in range(t)})
    idx = {tp: n for n, tp in enumerate(tuples)}
    rows = [[(idx[tuple(s[j] for j in range(m) if j not in (i, t))], s[t] < s[i]) for i in range(t)]
            for s in perms]
    best = 0
    for bits in range(2 ** len(tuples)):
        hits = sum(1 for r in rows
                   if sum(1 for (ti, tmin) in r if (((bits >> ti) & 1) == 0) == tmin) + k >= m)
        best = max(best, hits)
    return Q(best, len(perms))


for k in (1, 2):
    w = order_reading_worst(4, k)
    assert w == Q(pair_bound(k), 4)
    print(f"Slice 9, order-reading fits, m = 4, k = {k}: the worst refuses {w} of the orderings, "
          f"the symmetric bound {Q(pair_bound(k), 4)} (proved since: jackknifePlus_orderReading)")


# A fit that knows its units' positions breaks the jackknife+. Swapping the test unit with training
# unit i flips the ordering's parity. Over the orderings of a fixed bag, the fit leaving out position i
# and the test knows the two missing points and the missing positions, so it can make stranger the
# point that sits at the test position in the even completion. Then in every even ordering the test
# unit is stranger than every training unit, and in every odd one than none: coverage 1/2 at every
# level 1 <= k <= n, against 1 - (2k - 1)/(n + 1). The fit reads only its training units and where
# they sit; the residuals are realized by predictions at the two missing x, on the upper side.
def parity(seq):
    return sum(1 for a, b in itertools.combinations(range(len(seq)), 2) if seq[a] > seq[b]) % 2


def jack_refuses_positional(k, order, data):
    """The jackknife+ for the test at the last position with the parity fit; order[j] is the index
    of the point at position j."""
    m = len(order)
    t = m - 1
    y = data[order[t]][1]
    vhi = []
    for i in range(t):
        seen = {j: order[j] for j in range(m) if j not in (i, t)}
        a, b = sorted(set(range(m)) - set(seen.values()))
        comp = [None] * m
        for j, u in seen.items():
            comp[j] = u
        comp[i], comp[t] = b, a
        stranger = a if parity(comp) == 0 else b
        other = b if stranger == a else a
        pred = {data[stranger][0]: data[stranger][1] - 2, data[other][0]: data[other][1] - 1}
        mu = lambda x: pred[x]
        r_i = abs(data[order[i]][1] - mu(data[order[i]][0]))
        vhi.append(mu(data[order[t]][0]) + r_i)
    return sum(1 for v in vhi if v < y) + k >= m


rows = []
for m in range(3, 8):
    data = as_q([(j, 0) for j in range(m)])
    perms = list(itertools.permutations(range(m)))
    past_sharp, past_two_alpha = [], []
    for k in range(1, m):
        refused = sum(1 for p in perms if jack_refuses_positional(k, p, data))
        assert Q(refused, len(perms)) == Q(1, 2)
        if Q(1, 2) > Q(pair_bound(k), m):
            past_sharp.append(k)
        if Q(1, 2) > Q(2 * k, m):
            past_two_alpha.append(k)
    rows.append(f"m = {m}: past (2k - 1)/m at k in {past_sharp}, past 2k/m at k in {past_two_alpha}")
print("Slice 9, a fit that knows its units' positions refuses exactly 1/2 at every level: " + "; ".join(rows))

print()
print("=== Slice 9: CV+ with K folds ===")


# Conventions. n = K m training units in K folds of m units (by position), the test unit last.
# Fold k's fit leaves out fold k (and the test unit); unit i's CV residual is under its fold's fit,
# and so is the test unit's comparison with i. The level is l = floor(alpha (n + 1)). CV+'s interval
# is the jackknife+'s with mu_{-S_k(i)} for mu_{-i}; as counts, y is refused when at least n + 1 - l
# upper ends lie strictly below it or n + 1 - l lower ends strictly above it.
# CV-a: r_k is the test unit's rank (ties counted for it) among fold k and itself, in 1..m + 1.
#   A refusal has sum_k r_k < l + K, and P(sum_k r_k < l + K) <= 2(l + K - 1)/(n + K) whatever
#   the fits do, provided each fold's fit does not read the fold or the test unit.
# CV-b: with m - 1 further test units forming a test fold, at most 2l + m - 2 units are strange in
#   the cross-fold tournament, so P <= (2l + m - 2)/(n + m) for a symmetric fit.
def folds_of(n, K):
    m = n // K
    return [list(range(k * m, (k + 1) * m)) for k in range(K)]


def cv_parts(z, fit, K, t, train):
    """For each training unit i: its fold's fit at x_t and i's CV residual."""
    n = len(train)
    out = []
    for fold in folds_of(n, K):
        units = [train[j] for j in fold]
        f = fit(tuple(z[u] for u in train if u not in units))
        for u in units:
            out.append((u, f(z[t][0]), abs(z[u][1] - f(z[u][0])), units))
    return out


def cv_refuses(l, z, fit, K, t, train):
    n, y = len(train), z[t][1]
    parts = cv_parts(z, fit, K, t, train)
    vlo = [mu - r for _, mu, r, _ in parts]
    vhi = [mu + r for _, mu, r, _ in parts]
    lo, hi = jack_ends(vlo, vhi, l)
    refused = y < lo or y > hi
    assert refused == (sum(1 for v in vhi if v < y) + l >= n + 1 or sum(1 for v in vlo if v > y) + l >= n + 1)
    return refused


def cv_ranks(z, fit, K, t, train):
    """r_k: 1 + the units of fold k whose conformity is at most the test unit's, under fold k's fit."""
    n = len(train)
    ranks = []
    for fold in folds_of(n, K):
        units = [train[j] for j in fold]
        f = fit(tuple(z[u] for u in train if u not in units))
        st = -abs(z[t][1] - f(z[t][0]))
        ranks.append(1 + sum(1 for u in units if -abs(z[u][1] - f(z[u][0])) <= st))
    return ranks


def fit_first(train):
    """Reads its training units' order: predicts the first unit's y. Not symmetric."""
    return lambda x: train[0][1]


def fit_weighted(train):
    """Weights its training units by position. Not symmetric."""
    w = [Q(j + 1) for j in range(len(train))]
    yb = sum(wj * p[1] for wj, p in zip(w, train)) / sum(w)
    return lambda x: yb


def bound_merged(l, K, n):
    return Q(2 * (l + K - 1), n + K)


def bound_tournament(l, m, n):
    return Q(max(2 * l + m - 2, 0), n + m)


# CV-a, the merging inequality, deterministically: if the K ranks sum to at most L, then
# L <= sum_k (2L/K - r_k)+. And the layer cake: for a rank super-uniform on 1..M,
# E (c - r)+ <= (1/M) sum_j (c - j)+ <= c^2/(2M). Checked on every rank vector and on super-uniform
# laws (uniform, and mass moved up).
n_cases = 0
for K in (1, 2, 3, 4):
    for M in range(1, 6):
        for r in itertools.product(range(1, M + 1), repeat=K):
            for L in range(1, K * M + 1):
                if sum(r) <= L:
                    assert L <= sum(max(Q(2 * L, K) - rk, 0) for rk in r)
                n_cases += 1
for M in range(1, 9):
    for cnum in range(0, 4 * M + 1):
        c = Q(cnum, 2)
        riemann = sum(max(c - j, 0) for j in range(1, M + 1))
        assert riemann <= c * c / 2
        for shift in range(M):
            law = [Q(0)] * (M + 1)
            for j in range(1, M + 1):
                law[min(M, j + (shift if j > M - shift else 0))] += Q(1, M)
            assert all(sum(law[1:j + 1]) <= Q(j, M) for j in range(1, M + 1))
            assert sum(law[j] * max(c - j, 0) for j in range(1, M + 1)) <= riemann / M
print(f"CV-a, merging: the deterministic inequality on {n_cases} (K, M, rank vector, L) cases; the layer "
      f"cake E (c - r)+ <= sum_j (c - j)+/M <= c^2/(2M) for super-uniform ranks")

# K = 2 with uniform ranks, the worst coupling (a permutation, by Birkhoff) against the bound
# 2L/(K M): the merged bound is within 1/M of the worst case here.
for M in range(2, 7):
    for L in range(2, M + 2):
        worst = max(sum(1 for i in range(1, M + 1) if i + pi[i - 1] <= L)
                    for pi in itertools.permutations(range(1, M + 1)))
        assert Q(2 * L, 2 * M) - Q(1, M) <= Q(worst, M) <= Q(2 * L, 2 * M)
        if M == 6 and L == 4:
            print(f"CV-a, K = 2, M = 6, L = 4: the worst coupling of two uniform ranks gives P(r1 + r2 <= 4) = "
                  f"{Q(worst, M)}, against the merged bound {Q(2 * L, 2 * M)}")

# CV-a and CV-b over orderings. CV-a over every ordering of n + 1 units, with symmetric fits and with
# two that read their units' order; CV-b over every ordering of n + m units, the test unit's fold
# made up by m - 1 further units the interval never reads, with symmetric fits.
cv_cases = [((2, 2), as_q([(0, 0), (1, 2), (2, 1), (3, 3), (4, 8)])),
            ((3, 2), as_q([(0, 0), (1, 2), (2, 1), (3, 3), (4, 0), (5, 1), (6, 7)])),
            ((2, 3), as_q([(0, 1), (1, 1), (2, 2), (3, 0), (4, 1), (5, 2), (6, 5)])),
            ((2, 2), as_q([(0, 1), (0, 1), (1, 1), (2, 1), (3, 2)]))]
cv_fits = {"mean": fit_mean, "line": fit_line, "nearest": fit_nearest,
           "first (order-reading)": fit_first, "weighted (order-reading)": fit_weighted}
n_orderings = 0
rows_a = []
for (K, m), z in cv_cases:
    n = K * m
    t = n
    train = list(range(n))
    perms = list(itertools.permutations(range(n + 1)))
    for name, fit in cv_fits.items():
        cov = {l: 0 for l in range(n + 2)}
        for p in perms:
            zp = compose(z, p)
            ranks = cv_ranks(zp, fit, K, t, train)
            for l in range(n + 2):
                ref = cv_refuses(l, zp, fit, K, t, train)
                assert (not ref) or sum(ranks) < l + K
                cov[l] += not ref
            n_orderings += 1
        for l in range(n + 2):
            c = Q(cov[l], len(perms))
            assert c >= 1 - bound_merged(l, K, n)
            assert c >= 1 - Q(l, m + 1)  # one fold (cvPlus_oneFold)
        rows_a.append(f"{name} (K = {K}, m = {m}), l = 1: {Q(cov[1], len(perms))} >= {1 - bound_merged(1, K, n)}")
print(f"CV-a over orderings: {n_orderings} orderings, 4 data sets x 5 fits (2 of them read their units' order); "
      f"every level meets 1 - 2(l + K - 1)/(n + K), and every refusal has sum_k r_k < l + K")
for row in rows_a[:5]:
    print("  " + row)

n_orderings = 0
for (K, m), z in [((2, 2), as_q([(0, 0), (1, 2), (2, 1), (3, 3), (4, 8), (5, 1)])),
                  ((3, 1), as_q([(0, 0), (1, 2), (2, 1), (3, 9)])),
                  ((2, 2), as_q([(0, 1), (1, 1), (2, 1), (3, 2), (4, 1), (5, 1)]))]:
    n = K * m
    train = list(range(n))
    perms = list(itertools.permutations(range(n + m)))
    for name in ("mean", "line", "nearest"):
        fit = cv_fits[name]
        cov = {l: 0 for l in range(n + 2)}
        for p in perms:
            zp = compose(z, p)
            for l in range(n + 2):
                cov[l] += not cv_refuses(l, zp, fit, K, n, train)
            n_orderings += 1
        for l in range(n + 2):
            assert Q(cov[l], len(perms)) >= 1 - bound_tournament(l, m, n)
print(f"CV-b over orderings of n + m units: {n_orderings} orderings, 3 data sets x 3 symmetric fits; every level "
      f"meets 1 - (2l + m - 2)/(n + m)")


# CV-b's counting core: K + 1 folds of m units, games only across folds, ties allowed. At most
# 2l + m - 2 units are stranger than n + 1 - l others (n = K m), and the count is attained.
def cross_comparisons(F, m):
    N = F * m
    pairs = [(i, j) for i, j in itertools.combinations(range(N), 2) if i // m != j // m]
    for choice in itertools.product(range(3), repeat=len(pairs)):
        P = [[0] * N for _ in range(N)]
        for (i, j), c in zip(pairs, choice):
            if c == 0:
                P[j][i] = 1
            elif c == 1:
                P[i][j] = 1
        yield P


n_cases = 0
for F, m in ((2, 2), (3, 2), (2, 3), (4, 1)):
    N, n = F * m, (F - 1) * m
    worst = [0] * (n + 2)
    for P in cross_comparisons(F, m):
        n_cases += 1
        sc = [sum(1 for j in range(N) if j // m != u // m and P[u][j] < P[j][u]) for u in range(N)]
        for l in range(n + 2):
            r = sum(1 for c in sc if c + l >= n + 1)
            assert r <= max(2 * l + m - 2, 0)
            worst[l] = max(worst[l], r)
    print(f"CV-b counting core, {F} folds of {m}: the most refused at l = 0..{n + 1} is {worst}, "
          f"bound {[min(max(2 * l + m - 2, 0), N) for l in range(n + 2)]}")
print(f"CV-b counting core: {n_cases} outcomes of the cross-fold comparisons, ties included")

# Which bound serves, at alpha = 0.1: the merged one for few folds, the tournament for many.
rows = []
for n, K in ((20, 2), (20, 5), (20, 20), (100, 5), (100, 10), (100, 100), (1000, 5), (1000, 10), (1000, 1000)):
    l = math.floor(Q(1, 10) * (n + 1))
    a, b = bound_merged(l, K, n), bound_tournament(l, n // K, n)
    rows.append(f"n = {n}, K = {K}: merged {float(a):.3f}, tournament {float(b):.3f}")
print("CV+ miscoverage bounds at alpha = 0.1: " + "; ".join(rows))

# The counts sit inside Barber et al.'s alpha forms: merged, 2alpha + 2(1 - alpha)(K - 1)/(n + K)
# (Vovk's term); tournament, 2alpha + (1 - K/n)/(K + 1).
n_cases = 0
for K in range(1, 21):
    for m in range(1, 21):
        n = K * m
        for a in range(0, 21):
            alpha = Q(a, 40)
            l = math.floor(alpha * (n + 1))
            assert bound_merged(l, K, n) <= 2 * alpha + 2 * (1 - alpha) * Q(K - 1, n + K)
            assert bound_tournament(l, m, n) <= 2 * alpha + (1 - Q(K, n)) / (K + 1)
            n_cases += 1
print(f"CV+ bounds inside Barber et al.'s alpha forms on {n_cases} (K, m, alpha) cases")

# CV-b on a fixed bag of n + 1 units, without the further test units. A refusal leaves the test unit
# within the top l of fold 1 with itself, so P <= l/(m + 1) on any bag (cvPlus_oneFold). With K = 2 and
# m >= 2 that is within CV-b's bound at every level, and m = 1 is the jackknife+: K = 2 is settled.
# For K = 3, m = 2, l = 2 a refusal needs the test on top in two of its three triples, which it tops
# once on average: P <= 1/2, CV-b's bound. The rest is searched in conformal-searches.py.
for m in range(2, 40):
    n = 2 * m
    for l in range(0, n + 2):
        if bound_tournament(l, m, n) < 1:
            assert Q(l, m + 1) <= bound_tournament(l, m, n)
for m in range(2, 40):
    for alpha in (Q(1, 20), Q(1, 10), Q(1, 5)):
        l = math.floor(alpha * (2 * m + 1))
        assert Q(l, m + 1) < 2 * alpha
assert Q(3, 3 * 2) == bound_tournament(2, 2, 6)
print("CV-b on a fixed bag: K = 2 settled by one fold (and two-fold CV+ covers above 1 - 2alpha); "
      "K = 3, m = 2, l = 2 settled by Markov at 1/2")

print()
print("=== Slice 9: asymmetric CQR, a threshold per tail ===")


# Each tail is CB-1's rule on its own conformity, fitted on the training split: y - q_lo(x) below,
# q_hi(x) - y above, at levels k_lo and k_hi, each floor((alpha/2) m) at alpha/2. The set is
# [q_lo(x) + c_lo(k_lo), q_hi(x) - c_hi(k_hi)], with c(k) the k-th smallest calibration conformity
# of that tail, and an end is infinite at level 0. Ties count for the test, as in CB-1.
def lo_conf(qlo):
    return lambda p: p[1] - qlo[p[0]]


def hi_conf(qhi):
    return lambda p: qhi[p[0]] - p[1]


xy9 = [(x, y) for x in (0, 1) for y in range(4)]
pxy9 = {(0, 0): Q(1, 8), (0, 1): Q(1, 4), (0, 2): Q(1, 16), (0, 3): Q(1, 16),
        (1, 0): Q(1, 16), (1, 1): Q(1, 16), (1, 2): Q(1, 8), (1, 3): Q(1, 4)}
qlo9, qhi9 = {0: Q(1, 2), 1: Q(2)}, {0: Q(1), 1: Q(5, 2)}
slo9, shi9 = lo_conf(qlo9), hi_conf(qhi9)

# Under an i.i.d. law, N = 4 units, every pair of levels: each tail refuses at most k/N, and the set
# covers at least 1 - (k_lo + k_hi)/N.
N = 4
law9 = {}
for z in itertools.product(xy9, repeat=N):
    w = Q(1)
    for p in z:
        w *= pxy9[p]
    law9[z] = w
for klo in range(N + 1):
    for khi in range(N + 1):
        lo_ok = hi_ok = both = Q(0)
        for z, w in law9.items():
            a = accepts(klo, tuple(slo9(p) for p in z), N - 1, range(N))
            b = accepts(khi, tuple(shi9(p) for p in z), N - 1, range(N))
            lo_ok += w * a
            hi_ok += w * b
            both += w * (a and b)
        assert lo_ok >= 1 - Q(klo, N) and hi_ok >= 1 - Q(khi, N) and both >= 1 - Q(klo + khi, N)
        if klo == khi == 1:
            print(f"Asymmetric CQR, i.i.d. law, N = 4, k_lo = k_hi = 1: below the set {1 - lo_ok}, above it "
                  f"{1 - hi_ok} (each at most 1/4), covered {both} (at least 1/2), ties included")

# Over orderings with distinct scores in each tail: each tail refuses exactly k/m, and two-sided
# coverage lies between 1 - (k_lo + k_hi)/m and 1 - max(k_lo, k_hi)/m.
z_cont = as_q([(0, Q(3, 10)), (1, Q(17, 10)), (0, Q(22, 10)), (1, Q(4, 10)), (0, Q(11, 10))])
m = len(z_cont)
assert injective_on([slo9(p) for p in z_cont], range(m)) and injective_on([shi9(p) for p in z_cont], range(m))
G = perms_on(m, range(m))
for klo in range(m + 1):
    for khi in range(m + 1):
        a = coverage_orderings(z_cont, lambda zz: tuple(slo9(p) for p in zz), klo, m - 1, range(m))
        b = coverage_orderings(z_cont, lambda zz: tuple(shi9(p) for p in zz), khi, m - 1, range(m))
        assert a == 1 - Q(klo, m) and b == 1 - Q(khi, m)
        both = Q(sum(1 for p in G if accepts(klo, tuple(slo9(q) for q in compose(z_cont, p)), m - 1, range(m))
                     and accepts(khi, tuple(shi9(q) for q in compose(z_cont, p)), m - 1, range(m))), len(G))
        assert 1 - Q(klo + khi, m) <= both <= 1 - Q(max(klo, khi), m)
print(f"Asymmetric CQR over orderings, m = {m}, distinct scores: each tail refuses exactly k/m at every "
      f"level, and two-sided coverage lies in [1 - (k_lo + k_hi)/m, 1 - max(k_lo, k_hi)/m]")

# The ends: on a grid in y, the accepted set is [q_lo + c_lo(k_lo), q_hi - c_hi(k_hi)], ties to the
# test and an infinite end at level 0; Romano et al.'s Q_{1 - alpha_lo}, the
# ceil((1 - alpha_lo)(n2 + 1))-th smallest of q_lo(x_i) - y_i (infinite past n2), gives the same lower
# end with k_lo = floor(alpha_lo (n2 + 1)).
grid9 = [Q(i, 4) for i in range(-12, 25)]
n_sets = 0
for cal in itertools.product(xy9, repeat=3):
    clo = sorted(slo9(p) for p in cal)
    chi = sorted(shi9(p) for p in cal)
    elo = sorted(qlo9[p[0]] - p[1] for p in cal)
    for x in (0, 1):
        for klo in range(4):
            for khi in range(4):
                lo_end = -INF if klo == 0 else qlo9[x] + clo[klo - 1]
                hi_end = INF if khi == 0 else qhi9[x] - chi[khi - 1]
                acc = [y for y in grid9
                       if accepts(klo, tuple(slo9(p) for p in cal) + (slo9((x, y)),), 3, range(4))
                       and accepts(khi, tuple(shi9(p) for p in cal) + (shi9((x, y)),), 3, range(4))]
                assert acc == [y for y in grid9 if lo_end <= y <= hi_end]
                n_sets += 1
        for a in range(1, 20):
            alpha_lo = Q(a, 20)
            K = math.ceil((1 - alpha_lo) * 4)
            k = math.floor(alpha_lo * 4)
            assert K == 4 - k
            q_end = -INF if K > 3 else qlo9[x] - elo[K - 1]
            assert q_end == (-INF if k == 0 else qlo9[x] + clo[k - 1])
print(f"Asymmetric CQR, the ends: {n_sets} accept sets on a grid, each [q_lo + c_lo(k_lo), q_hi - c_hi(k_hi)]; "
      f"Romano et al.'s per-tail quantile gives the same ends")

# Symmetric CQR controls only the total. With both quantile fits below the data, its conformity
# min(y - q_lo, q_hi - y) is always the upper one, and all its miscoverage falls above the set.
qlo_low, qhi_low = {0: Q(-2), 1: Q(-2)}, {0: Q(-1), 1: Q(-1)}
s_sym = cqr_score(qlo_low, qhi_low)
s_up = hi_conf(qhi_low)
s_dn = lo_conf(qlo_low)
sym_above = sym_below = asym_above = Q(0)
for z, w in law9.items():
    sym = tuple(s_sym(p) for p in z)
    if not accepts(2, sym, N - 1, range(N)):
        if s_up(z[-1]) <= s_dn(z[-1]):
            sym_above += w
        else:
            sym_below += w
    if not accepts(1, tuple(s_up(p) for p in z), N - 1, range(N)):
        asym_above += w
assert sym_below == 0 and sym_above > Q(1, 4) >= asym_above
print(f"Symmetric CQR at k = 2 (alpha = 1/2), quantile fits below the data: above the set {sym_above}, "
      f"below 0; asymmetric at k = 1 per tail: above {asym_above} (at most 1/4)")

# Small calibration sets and rounding. A tail's end is finite only once floor((alpha/2) m) >= 1,
# m >= 2/alpha, twice the symmetric form's m >= 1/alpha. And 2 floor(alpha m / 2) is floor(alpha m)
# or one less. Equal-tail rounding adds 0 or 1/m to the symmetric finite-rank bound.
# Its excess over nominal 1 - alpha lies in [0, 2/m).
rows9 = []
for alpha in (Q(1, 5), Q(1, 10), Q(1, 20)):
    m_sym = min(mm for mm in range(1, 400) if math.floor(alpha * mm) >= 1)
    m_asym = min(mm for mm in range(1, 400) if math.floor(alpha / 2 * mm) >= 1)
    assert m_asym == math.ceil(2 / alpha) and m_sym == math.ceil(1 / alpha)
    rows9.append(f"alpha = {float(alpha)}: m >= {m_sym} symmetric, m >= {m_asym} per tail")
for mm in range(1, 200):
    for alpha in (Q(1, 5), Q(1, 10), Q(1, 20)):
        assert 2 * math.floor(alpha * mm / 2) in (math.floor(alpha * mm), math.floor(alpha * mm) - 1)
print("Asymmetric CQR, finite ends need " + "; ".join(rows9) + ". At alpha = 0.1, m = 30: "
      f"two-sided {1 - Q(2 * math.floor(Q(3, 2)), 30)} per-tail against {1 - Q(3, 30)} symmetric")


print()
print("=== Slice 10: the full step's error, completed ===")


# The CDF score's equi-tailed set accepts when the candidate is in neither tail: CB-1's rule on F
# (the lower tail) and on -F (the upper), k per tail. Estimates within per-unit bounds eps.
def neg(s):
    return tuple(-v for v in s)


def two_sided(k, s, u, C):
    return accepts(k, s, u, C) and accepts(k, neg(s), u, C)


def two_sided_inflated(k, s, u, C, eps):
    return accepts_inflated(k, s, u, C, eps) and accepts_inflated(k, neg(s), u, C, eps)


def score_ecdf(z):
    """An in-sample CDF score, symmetric in the units: the share at or below each unit."""
    return tuple(Q(sum(1 for y in z if y <= x), len(z)) for x in z)


n_ts = 0
worst_overlap = None
for N in (3, 4, 5):
    C = list(range(N))
    t = N - 1
    for eps0 in ((Q(1, 50),) * N, tuple(Q(i + 1, 60) for i in range(N)), (Q(2, 5),) * N):
        by_unit = lambda p, eps0=eps0: tuple(eps0[p[i]] for i in range(len(p)))
        for score in (score_sq, score_ecdf):
            for z in multisets(range(3), N):
                s = score(z)

                def adv(w, p, score=score, eps0=eps0):
                    e = by_unit(p)
                    base = score(w)
                    # the test unit pushed toward its nearer tail, the rest away
                    lo = rank_le(base, t, C) <= rank_le(neg(base), t, C)
                    sign = -1 if lo else 1
                    return tuple(v + sign * e[i] if i == t else v - sign * e[i] for i, v in enumerate(base))

                def noisy(w, p, score=score, eps0=eps0):
                    e = by_unit(p)
                    return tuple(v + e[i] * Q(random.Random(hash((p, i, 7))).randint(-10, 10), 10)
                                 for i, v in enumerate(score(w)))

                for k in range(0, N // 2 + 1):
                    R = {u for u in C if not two_sided(k, s, u, C)}
                    Mlo = {u for u in C if fragile(k, s, u, C, eps0)}
                    Mhi = {u for u in C if fragile(k, neg(s), u, C, eps0)}
                    assert len(R) <= 2 * k
                    for u in Mlo & Mhi:
                        near = sum(1 for j in C if j == u or abs(s[j] - s[u]) < eps0[j] + eps0[u])
                        assert near >= N - 2 * k + 2
                        worst_overlap = (N, k, near) if worst_overlap is None else worst_overlap
                    for est in (adv, noisy):
                        G = perms_on(N, C)
                        hits = sum(1 for p in G if two_sided(k, est(compose(z, p), p), t, C))
                        assert Q(hits, len(G)) >= 1 - Q(len(R | Mlo | Mhi), N)
                        inf = sum(1 for p in G if two_sided_inflated(k, est(compose(z, p), p), t, C, by_unit(p)))
                        assert Q(inf, len(G)) >= 1 - Q(2 * k, N)
                        for p in G:
                            sh = est(compose(z, p), p)
                            e = by_unit(p)
                            Mhat = {u for u in C if fragile_est(k, sh, u, C, e) or fragile_est(k, neg(sh), u, C, e)}
                            assert len(Mlo | Mhi) <= len(Mhat)
                    n_ts += 1
print(f"Slice 10, the equi-tailed set: {n_ts} (data, score, bounds, k) cases over every ordering, an "
      f"adversarial and a noisy estimate: coverage >= 1 - |R u M_lo u M_hi|/m with |R| <= 2k, the "
      f"inflated set >= 1 - 2k/m, every ordering's estimated union band at least the exact one; a unit "
      f"fragile in both tails has at least m - 2k + 2 units within its bounds (first such case {worst_overlap})")


# The PIT. Its rank R = B + J, J uniform on 1..T+1 (theta's ceil), so P(R <= k | ordering) is
# sum_{r <= k} tieprob. The exact score's shifted counts bound it from both sides:
# P(R_hat <= k) <= #{u : L(u) < k}/m and >= #{u : U(u) <= k}/m, with L(u) the units at least
# eps_j + eps_u below u and U(u) the units at most eps_j + eps_u above it (u among them).
def tieprob(s, u, C, r):
    lt, le = rank_lt(s, u, C), rank_le(s, u, C)
    return Q(1, le - lt) if lt < r <= le else Q(0)


def L_count(s, u, C, eps):
    return sum(1 for j in C if s[j] + eps[j] + eps[u] < s[u])


def U_count(s, u, C, eps):
    return sum(1 for j in C if s[j] <= s[u] + eps[j] + eps[u])


def Lhat(s, u, C, eps):
    return sum(1 for j in C if s[j] + 2 * (eps[j] + eps[u]) < s[u])


def Uhat(s, u, C, eps):
    return sum(1 for j in C if s[j] <= s[u] + 2 * (eps[j] + eps[u]))


n_pit, worst_dev = 0, Q(0)
for N in (3, 4, 5):
    C = list(range(N))
    t = N - 1
    for eps0 in ((Q(1, 50),) * N, tuple(Q(i + 1, 60) for i in range(N)), (Q(1, 4),) * N):
        by_unit = lambda p, eps0=eps0: tuple(eps0[p[i]] for i in range(len(p)))
        for z in multisets(range(3), N):
            s = score_ecdf(z)
            G = perms_on(N, C)
            for est in ("exact", "noisy"):
                def sh_of(p):
                    base = score_ecdf(compose(z, p))
                    if est == "exact":
                        return base
                    e = by_unit(p)
                    return tuple(v + e[i] * Q(random.Random(hash((p, i, 3))).randint(-10, 10), 10)
                                 for i, v in enumerate(base))
                cdf = {k: sum(sum(tieprob(sh_of(p), t, C, r) for r in range(1, k + 1)) for p in G) / len(G)
                       for k in range(N + 1)}
                for k in range(N + 1):
                    up = Q(sum(1 for u in C if L_count(s, u, C, eps0) < k), N)
                    lo = Q(sum(1 for u in C if U_count(s, u, C, eps0) <= k), N)
                    assert lo <= Q(k, N) <= up
                    if est == "exact":
                        assert cdf[k] == Q(k, N)
                    else:
                        assert lo <= cdf[k] <= up
                        worst_dev = max(worst_dev, abs(cdf[k] - Q(k, N)))
                        for p in G:
                            sh, e = sh_of(p), by_unit(p)
                            assert sum(1 for u in C if L_count(s, u, C, eps0) < k) <= \
                                sum(1 for u in C if Lhat(sh, u, C, e) < k)
                            assert sum(1 for u in C if U_count(s, u, C, eps0) <= k) >= \
                                sum(1 for u in C if Uhat(sh, u, C, e) <= k)
                # Between grid values the PIT's CDF interpolates linearly, so its largest departure
                # from u is at a grid value.
                for uu in [Q(i, 7 * N) for i in range(7 * N + 1)]:
                    kk = min(int(uu * N), N - 1)
                    th = uu * N - kk
                    val = (1 - th) * cdf[kk] + th * cdf[kk + 1]
                    assert abs(val - uu) <= max(abs(cdf[kk] - Q(kk, N)), abs(cdf[kk + 1] - Q(kk + 1, N)))
                n_pit += 1
print(f"Slice 10, the PIT: {n_pit} (data, bounds, estimate) cases: the exact score's PIT is uniform at "
      f"every grid value; the estimate's P(R <= k) lies between #{{U <= k}}/m and #{{L < k}}/m, which "
      f"bracket k/m; each ordering's estimated counts dominate; the departure between grid values is "
      f"at most the larger at its ends (largest departure seen {worst_dev})")


# The carried bound that leaves out the degenerate turns D. Bounds hold, with chance 1 - nu, on the
# orderings whose candidate is outside D; on D's turns and on failure the estimate is anything
# (here the adversary that refuses). Coverage >= 1 - |D u R u M_D|/m - nu, where M_D is the margin
# under bounds carried over the kept turns only, and the inflated rule >= 1 - |D u R|/m - nu.
nu = Q(1, 10)
n_d = 0
for N in (4, 5):
    C = list(range(N))
    t = N - 1
    for z in multisets(range(3), N)[::2]:
        s = score_sq(z)
        for Dset in ({0}, {N - 1}, {0, 1}):
            kept_eps = tuple(Q(1, 40) for _ in C)  # carried over the kept turns
            G = perms_on(N, C)
            for k in range(0, N):
                R = {u for u in C if not accepts(k, s, u, C)}
                MD = {u for u in C if fragile(k, s, u, C, kept_eps)}
                hits = Q(0)
                hits_inf = Q(0)
                for p in G:
                    base = score_sq(compose(z, p))
                    e = tuple(kept_eps[p[i]] for i in range(N))
                    good = tuple(v - e[i] if i == t else v + e[i] for i, v in enumerate(base))
                    if p[t] in Dset:
                        continue  # D's turn: refused, whatever w; on failure (chance nu) refused too
                    hits += (1 - nu) * accepts(k, good, t, C)
                    hits_inf += (1 - nu) * accepts_inflated(k, good, t, C, e)
                cov, cov_inf = hits / len(G), hits_inf / len(G)
                bound = 1 - Q(len(Dset | R | MD), N) - nu
                assert cov >= bound >= 1 - Q(k + len(MD), N) - Q(len(Dset), N) - nu
                assert cov_inf >= 1 - Q(len(Dset | R), N) - nu
                n_d += 1
print(f"Slice 10, leaving out the degenerate turns: {n_d} (data, D, k) cases: coverage >= "
      f"1 - |D u R u M_D|/m - nu (and so >= 1 - (k + |M_D|)/m - |D|/m - nu); the inflated rule >= "
      f"1 - |D u R|/m - nu")
# The same for the equi-tailed set and the PIT. D's turns and failed orderings read the worst
# case: refused for coverage, 1 for the PIT's upper bound and 0 for its lower. The equi-tailed set
# covers >= 1 - |D u R u M_lo u M_hi|/m - nu (inflated >= 1 - |D u R|/m - nu), and
# #{u not in D : U(u) <= k}/m - nu <= P(PIT <= k/m) <= #{u in D or L(u) < k}/m + nu.
n_d2 = 0
for N in (4, 5):
    C = list(range(N))
    t = N - 1
    kept_eps = tuple(Q(1, 40) for _ in C)
    G = perms_on(N, C)
    for z in multisets(range(3), N)[::2]:
        s = score_ecdf(z)
        ests = {}
        for p in G:
            base = score_ecdf(compose(z, p))
            e = tuple(kept_eps[p[i]] for i in range(N))
            sign = -1 if rank_le(base, t, C) <= rank_le(neg(base), t, C) else 1
            adv = tuple(v + sign * e[i] if i == t else v - sign * e[i] for i, v in enumerate(base))
            noisy = tuple(v + e[i] * Q(random.Random(hash((p, i, 11))).randint(-10, 10), 10)
                          for i, v in enumerate(base))
            ests[p] = (e, adv, noisy)
        for Dset in ({0}, {N - 1}, {0, 1}):
            for k in range(0, N // 2 + 1):
                R2 = {u for u in C if not two_sided(k, s, u, C)}
                Mlo = {u for u in C if fragile(k, s, u, C, kept_eps)}
                Mhi = {u for u in C if fragile(k, neg(s), u, C, kept_eps)}
                kept = [p for p in G if p[t] not in Dset]
                cov = sum((1 - nu) * two_sided(k, ests[p][1], t, C) for p in kept) / len(G)
                cov_inf = sum((1 - nu) * two_sided_inflated(k, ests[p][1], t, C, ests[p][0])
                              for p in kept) / len(G)
                assert cov >= 1 - Q(len(Dset | R2 | Mlo | Mhi), N) - nu
                assert cov_inf >= 1 - Q(len(Dset | R2), N) - nu
                # What the engine reads from one ordering whose bounds hold: k (2k) plus the
                # estimated band plus the units of D the inflated rule accepts.
                R1 = {u for u in C if not accepts(k, s, u, C)}
                for p in G:
                    e = ests[p][0]
                    for est in ests[p][1:]:
                        Mh1 = {i for i in C if fragile_est(k, est, i, C, e)}
                        Mh2 = Mh1 | {i for i in C if fragile_est(k, neg(est), i, C, e)}
                        Di1 = {i for i in C if p[i] in Dset and accepts_inflated(k, est, i, C, e)}
                        Di2 = {i for i in C if p[i] in Dset and two_sided_inflated(k, est, i, C, e)}
                        assert len(Dset | R1 | Mlo) <= k + len(Mh1) + len(Di1)
                        assert len(Dset | R2 | Mlo | Mhi) <= 2 * k + len(Mh2) + len(Di2)
                n_d2 += 1
            for k in range(N + 1):
                cdf = {p: sum(tieprob(ests[p][2], t, C, r) for r in range(1, k + 1)) for p in G}
                up = sum(1 if p[t] in Dset else (1 - nu) * cdf[p] + nu for p in G) / len(G)
                lo = sum(0 if p[t] in Dset else (1 - nu) * cdf[p] for p in G) / len(G)
                assert up <= Q(sum(1 for u in C if u in Dset or L_count(s, u, C, kept_eps) < k), N) + nu
                assert lo >= Q(sum(1 for u in C if u not in Dset and U_count(s, u, C, kept_eps) <= k), N) - nu
                for p in G:
                    e, est = ests[p][0], ests[p][2]
                    assert sum(1 for u in C if u in Dset or L_count(s, u, C, kept_eps) < k) <=                         sum(1 for i in C if p[i] in Dset or Lhat(est, i, C, e) < k)
                    assert sum(1 for i in C if p[i] not in Dset and Uhat(est, i, C, e) <= k) <=                         sum(1 for u in C if u not in Dset and U_count(s, u, C, kept_eps) <= k)
                n_d2 += 1
print(f"Slice 10, leaving out the degenerate turns, two-sided and the PIT: {n_d2} (data, D, k) cases: "
      f"the equi-tailed set >= 1 - |D u R u M_lo u M_hi|/m - nu, inflated >= 1 - |D u R|/m - nu; "
      f"#{{u not in D : U <= k}}/m - nu <= P(PIT <= k/m) <= #{{u in D or L < k}}/m + nu; from any "
      f"ordering's estimates, |D u R u M| <= k (2k) + |M_hat| + |D's units the inflated rule accepts|, "
      f"and the PIT's counts with D read from the estimates bound the exact ones")
# What the engine reads, leaving out the degenerate turns (follow-up, 2026-11-10). The engine holds
# one ordering's estimates, the realized one, and when its own test unit's turn is in D it reads
# nothing there. Per data set: a count c (at most B) that every kept ordering whose bounds hold
# reads from above, read as 0 on the left-out orderings, is at most the average reading plus B
# times D's share and the failure share; a count read from below needs no D charge. So coverage >=
# 1 - avg_sigma[1[kept] c_hat(sigma)]/m - |D|/m - 2 nu, and the PIT's counts likewise.
rr = random.Random(1213)
n_read = 0
for trial in range(400):
    N = rr.choice((3, 4))
    C = list(range(N))
    t = N - 1
    G = perms_on(N, C)
    W = rr.choice((1, 2, 3))
    pw = [Q(rr.randint(1, 5)) for _ in range(W)]
    pw = [x / sum(pw) for x in pw]
    B = Q(rr.randint(1, 6))
    lhs_up, rhs_up, lhs_lo, rhs_lo, nu_w, shareD = Q(0), Q(0), Q(0), Q(0), Q(0), Q(0)
    for w in range(W):
        D = {u for u in C if rr.random() < 0.3}
        hold = {p: rr.random() < 0.8 for p in G}
        c = Q(rr.randint(0, 12), 2) * B / 6
        up = {p: (c + Q(rr.randint(0, 4), 3) if (p[t] not in D and hold[p]) else Q(rr.randint(0, 9), 3)) for p in G}
        lo = {p: (max(Q(0), c - Q(rr.randint(0, 4), 3)) if (p[t] not in D and hold[p]) else Q(rr.randint(0, 9), 3) * B / 3)
              for p in G}
        lo = {p: min(v, B) for p, v in lo.items()}
        fail = Q(sum(1 for p in G if p[t] not in D and not hold[p]), len(G))
        lhs_up += pw[w] * c
        rhs_up += pw[w] * sum(Q(0) if p[t] in D else up[p] for p in G) / len(G)
        rhs_lo += pw[w] * sum(Q(0) if p[t] in D else lo[p] for p in G) / len(G)
        nu_w += pw[w] * fail
        shareD += pw[w] * Q(len(D), N)
    assert lhs_up <= rhs_up + B * (shareD + nu_w)
    assert rhs_lo - B * nu_w <= lhs_up
    n_read += 1
# The engine's reading itself: in the slice 10 setting, each kept ordering reads
# k + |M_hat| + |D's units the inflated rule accepts| from its own estimates; coverage is at least
# 1 - (the kept orderings' average reading)/m - |D|/m - 2 nu.
n_eng = 0
for N in (4, 5):
    C = list(range(N))
    t = N - 1
    G = perms_on(N, C)
    for z in multisets(range(3), N)[::2]:
        s = score_sq(z)
        for Dset in ({0}, {N - 1}, {0, 1}):
            kept_eps = tuple(Q(1, 40) for _ in C)
            for k in range(0, N):
                hits, read = Q(0), Q(0)
                for p in G:
                    if p[t] in Dset:
                        continue
                    base = score_sq(compose(z, p))
                    e = tuple(kept_eps[p[i]] for i in range(N))
                    good = tuple(v - e[i] if i == t else v + e[i] for i, v in enumerate(base))
                    hits += (1 - nu) * accepts(k, good, t, C)
                    Mh = sum(1 for i in C if fragile_est(k, good, i, C, e))
                    Di = sum(1 for i in C if p[i] in Dset and accepts_inflated(k, good, i, C, e))
                    read += k + Mh + Di
                cov = hits / len(G)
                assert cov >= 1 - read / len(G) / N - Q(len(Dset), N) - 2 * nu
                n_eng += 1
print(f"Slice 10, reading from the realized ordering: {n_read} random data sets with readings above (and "
      f"below) the count on kept orderings whose bounds hold: count <= average kept reading + B (D share + nu), "
      f"and >= average kept reading - B nu; on {n_eng} (data, D, k) cases the engine's own reading gives "
      f"coverage >= 1 - average kept reading/m - |D|/m - 2 nu")
# Slice 10's grid reading (2026-11-10). (a) The PIT between grid values: with the randomized rank
# R = B + J and V uniform, P(PIT <= x) = sum_r P(R = r) clamp(m x - r + 1, 0, 1); at
# x = (k + theta)/m it is (1 - theta) P(R <= k) + theta P(R <= k + 1), so its departure from x is at
# most the larger at the two grid values. (b) The margin over the grid: the engine's report at the
# realized data is the largest kept reading and the largest flagged count over every value the test
# unit can take; both are free of the test unit's own value, and they bound the per-data-set
# reading form. A grid that misses a value the test unit can take can miss the reading there.
def clamp01(a):
    return max(Q(0), min(Q(1), a))


def pit_cdf(s, u, C, x):
    m = len(C)
    return sum(tieprob(s, u, C, r) * clamp01(m * x - r + 1) for r in range(1, m + 1))


def pit_rank_le(s, u, C, k):
    return sum((tieprob(s, u, C, r) for r in range(1, k + 1)), Q(0))


n_lin = 0
for N in (3, 4, 5):
    C = list(range(N))
    for z in multisets(range(3), N):
        s = score_ecdf(z)
        for u in C:
            for k in range(N):
                for th in (Q(0), Q(1, 3), Q(1, 2), Q(5, 7), Q(1)):
                    x = (k + th) / N
                    val = pit_cdf(s, u, C, x)
                    lo, hi = pit_rank_le(s, u, C, k), pit_rank_le(s, u, C, k + 1)
                    assert val == (1 - th) * lo + th * hi
                    assert abs(val - x) <= max(abs(lo - Q(k, N)), abs(hi - Q(k + 1, N)))
                    n_lin += 1
# averaged over orderings, with a noisy estimate in each: the same interpolation
for N in (3, 4):
    C, t = list(range(N)), N - 1
    G = perms_on(N, C)
    for z in multisets(range(3), N)[::3]:
        ests = {}
        for p in G:
            base = score_ecdf(compose(z, p))
            ests[p] = tuple(v + Q(random.Random(hash((p, i, 13))).randint(-10, 10), 400) for i, v in enumerate(base))
        for k in range(N):
            for th in (Q(1, 4), Q(2, 3)):
                x = (k + th) / N
                avg = sum(pit_cdf(ests[p], t, C, x) for p in G) / len(G)
                lo = sum(pit_rank_le(ests[p], t, C, k) for p in G) / len(G)
                hi = sum(pit_rank_le(ests[p], t, C, k + 1) for p in G) / len(G)
                assert avg == (1 - th) * lo + th * hi
                n_lin += 1

nu = Q(1, 10)
n_grid, n_free, misses = 0, 0, 0
for N in (3, 4):
    C, t = list(range(N)), N - 1
    G = perms_on(N, C)
    eps = tuple(Q(1, 40) for _ in C)

    def flagged(x):
        return {i for i, v in enumerate(x) if v == 2}  # a symmetric rule: the far value's turns

    def good(x):
        return tuple(v - eps[i] if i == t else v + eps[i] for i, v in enumerate(score_sq(x)))

    def reading(x, k):
        e, Dx = good(x), flagged(x)
        return (k + sum(1 for i in C if fragile_est(k, e, i, C, eps))
                + sum(1 for i in C if i in Dx and accepts_inflated(k, e, i, C, eps)))

    def upd(x, y):
        return x[:t] + (y,) + x[t + 1:]

    def grid_read(x, k, grid):
        return max(reading(upd(x, y), k) if t not in flagged(upd(x, y)) else 0 for y in grid)

    def grid_flags(x, grid):
        return max(len(flagged(upd(x, y))) for y in grid)

    for z in multisets(range(3), N):
        for k in range(N):
            hits = sum((1 - nu) * accepts(k, good(compose(z, p)), t, C)
                       for p in G if t not in flagged(compose(z, p))) / len(G)
            rep = sum(grid_read(compose(z, p), k, range(3)) + grid_flags(compose(z, p), range(3))
                      for p in G) / len(G)
            assert hits >= 1 - rep / N - 2 * nu
            for p in G:
                x = compose(z, p)
                if t not in flagged(x):
                    assert reading(x, k) <= grid_read(x, k, range(3))
                    if reading(x, k) > grid_read(x, k, (0, 2)):
                        misses += 1
                assert len(flagged(x)) <= grid_flags(x, range(3))
                assert all(grid_read(upd(x, y), k, range(3)) == grid_read(x, k, range(3)) for y in range(3))
                n_free += 1
            n_grid += 1
assert misses > 0
print(f"Slice 10's grid reading: the PIT's CDF interpolates linearly between grid values on {n_lin} "
      f"(data, unit, level) cases, so its departure is at most the larger at the ends; on {n_grid} "
      f"(data, k) cases coverage >= 1 - avg(largest kept reading + largest flagged count over every "
      f"value)/m - 2 nu, and the report is free of the test unit's own value ({n_free} arrangements); a grid "
      f"missing a value the test unit can take misses its reading in {misses} arrangements")
# What the margin can then read at m = 100: one unit far out (in D, and in the exact low tail),
# whose turn widens every carried bound to 3 against 0.002 over the kept turns.
mr = random.Random(20261110)
m, k, nu = 100, 5, Q(1, 20)
vals = sorted(mr.uniform(0, 1) for _ in range(m - 1)) + [-50.0]
sc = tuple(Q(v).limit_denominator(10**6) for v in vals)
Cm = list(range(m))
Rm = {u for u in Cm if not accepts(k, sc, u, Cm)}
M_all = {u for u in Cm if fragile(k, sc, u, Cm, (Q(3),) * m)}
M_kept = {u for u in Cm if fragile(k, sc, u, Cm, (Q(1, 500),) * m)}
Dm = {m - 1}
assert Dm <= Rm
print(f"Slice 10 at m = 100, k = 5, nu = 1/20: over every turn the margin reads "
      f"{float(1 - Q(len(Rm | M_all), m)):.2f} ({len(M_all)} fragile units); leaving out the far unit's "
      f"turn it reads {float(1 - Q(len(Dm | Rm | M_kept), m) - nu):.2f} (D inside R, so D costs nothing; "
      f"{len(M_kept)} fragile)")


print()
print("=== Slice 12: Hoff's score and the latent posterior ===")


# CB-4's model again (theta, a latent bit per unit, x and y). marg(zs) is the units' prior
# predictive with every latent integrated, the x law left out (it is the same factor for every
# theta, so it cancels in every ratio below).
def marg(zs):
    return sum(prior[th] * prod(g(th, z) for z in zs) for th in TH)


# Hoff's score c_i = p(z_i | every other unit) = p(z) / p(z_-i), with every latent integrated
# ((b) of CB-4): equivariant, so CB-1's rule (ties counted for the candidate) covers with it.
n_h = 0
for N, sets in ((3, list(itertools.product(Z, repeat=3))), (4, multisets(Z, 4))):
    C = list(range(N))
    for zs in sets:
        c = score_hoff(zs)
        assert all(c[i] == marg(zs) / marg(zs[:i] + zs[i + 1:]) for i in range(N))
        assert is_equivariant(zs, score_hoff, C)
        for k in range(N + 1):
            assert coverage_orderings(zs, score_hoff, k, N - 1, C) >= 1 - Q(k, N)
        n_h += 1
zs = ((0, 1),) * 4
assert coverage_orderings(zs, score_hoff, 1, 3, range(4), rule=accepts_strict) == 0
print(f"Slice 12, Hoff's score: c_i = p(z)/p(z_-i) on {n_h} data sets (every ordering of 3 units, "
      f"every multiset of 4); equivariant, and CB-1's rule covers >= 1 - k/m at every k; the strict "
      f"rule covers 0 when every unit ties")

# The draw identity. From the n-unit posterior, c_{n+1} is the plain n-unit predictive
# A = E_n[g(z_new | theta)], and each data unit's c_i = A / B_i with B_i = E_n[g(z_new|th)/g(z_i|th)]:
# the weight g(z_new)/g(z_i) is add-one-in and leave-one-out at once. The ranks need only the B_i:
# the candidate is accepted at level k iff at least k data units have B_i >= 1.
n_d = 0
for zs in itertools.product(Z, repeat=4):
    zn, znew = zs[:-1], zs[-1]
    pn = post(zn)
    A = sum(pn[th] * g(th, znew) for th in TH)
    B = [sum(pn[th] * g(th, znew) / g(th, z) for th in TH) for z in zn]
    c = score_hoff(zs)
    assert c[-1] == A and all(c[i] == A / B[i] for i in range(3))
    for k in range(5):
        assert accepts(k, c, 3, range(4)) == (sum(1 for b in B if b >= 1) >= k)
    n_d += 1
# From S draws the estimator is c_hat_i = A_hat / B_hat_i, c_hat_{n+1} = A_hat, with A_hat common to
# every unit: shifting every log score of an ordering by one constant changes no rule, so the
# bounds slice 10 carries are on log B_hat_i alone, and the candidate's bound is 0.
mr = random.Random(1210)
n_mc = 0
for trial in range(300):
    zs = tuple(mr.choice(Z) for _ in range(4))
    zn, znew = zs[:-1], zs[-1]
    pn = post(zn)
    draws = mr.choices(list(TH), weights=[float(w) for w in pn], k=mr.choice((1, 3, 10)))
    Ah = sum(g(th, znew) for th in draws) / len(draws)
    Bh = [sum(g(th, znew) / g(th, z) for th in draws) / len(draws) for z in zn]
    ch = tuple(Ah / b for b in Bh) + (Ah,)
    num = [sum(g(th, znew) / g(th, z) * g(th, z) for th in draws) / len(draws) for z in zn]
    assert all(ch[i] == num[i] / Bh[i] for i in range(3))  # the self-normalized weights' reading
    shifted = tuple(1 / b for b in Bh) + (Q(1),)  # every score divided by A_hat
    for k in range(5):
        assert accepts(k, ch, 3, range(4)) == accepts(k, shifted, 3, range(4))
        for cc in (Q(-3), Q(1, 7), Q(5)):
            s_add = tuple(v + cc for v in ch)
            assert accepts(k, s_add, 3, range(4)) == accepts(k, ch, 3, range(4))
            e = (Q(1, 50),) * 4
            assert accepts_inflated(k, s_add, 3, range(4), e) == accepts_inflated(k, ch, 3, range(4), e)
    n_mc += 1
print(f"Slice 12, the draw identity on all {n_d} data sets of 3 units and a candidate: c_(n+1) = A, "
      f"c_i = A/B_i, accepted iff #{{B_i >= 1}} >= k; on {n_mc} Monte Carlo draws the estimator "
      f"A_hat/B_hat_i ranks as 1/B_hat_i against 1, and a common shift changes neither rule")


# Hoff's score under latents. (c1): E[f(z_i | theta, v_i)] under the joint posterior of theta and
# every latent given the other units; unit i's latent is then at its prior, so (c1) is (b).
def score_hoff_joint(zs):
    out = []
    for i, z in enumerate(zs):
        tot, num = Q(0), Q(0)
        for th in TH:
            for vs in itertools.product(V, repeat=len(zs)):
                w = prior[th] * prod(qv[th][v] * (f[th][v][zj[0]][zj[1]] if j != i else 1)
                                     for j, (v, zj) in enumerate(zip(vs, zs)))
                tot += w
                num += w * f[th][vs[i]][z[0]][z[1]]
        out.append(num / tot)
    return tuple(out)


# (c2): each unit's own latent conditioned on its own value, theta given the other units:
# sum_theta p(theta | z_-i) sum_v q f(z_i)^2 / g(z_i).
def score_hoff_cond(zs):
    out = []
    for i, z in enumerate(zs):
        p = post(zs[:i] + zs[i + 1:])
        out.append(sum(p[th] * sum(qv[th][v] * f[th][v][z[0]][z[1]] ** 2 for v in V) / g(th, z)
                       for th in TH))
    return tuple(out)


# (c2) from the n-unit draws: each data unit reads its retained latent with Hoff's weight
# g(z_new|theta)/g(z_i|theta); the candidate draws its latent from p(v | theta, z_new), unweighted.
def score_hoff_cond_draws(zs):
    zn, znew = zs[:-1], zs[-1]
    Jn = joint(zn)
    out = []
    for i in range(len(zn)):
        wts = {key: w * g(key[0], znew) / g(key[0], zn[i]) for key, w in Jn.items()}
        tot = sum(wts.values())
        out.append(sum(w / tot * f[th][vs[i]][zn[i][0]][zn[i][1]] for (th, vs), w in wts.items()))
    cand = sum(w * sum(qv[th][v] * f[th][v][znew[0]][znew[1]] ** 2 for v in V) / g(th, znew)
               for (th, vs), w in Jn.items())
    return tuple(out) + (cand,)


def score_hoff_mixed(zs):
    """The broken reading: data units by (c2), the candidate by (b)."""
    return score_hoff_cond(zs)[:-1] + (score_hoff(zs)[-1],)


def score_hoff_mixed_reverse(zs):
    return score_hoff(zs)[:-1] + (score_hoff_cond(zs)[-1],)


gap = Q(0)
for zs in itertools.product(Z, repeat=3):
    hb = score_hoff(zs)
    assert score_hoff_joint(zs) == hb
    hc = score_hoff_cond(zs)
    assert hc == score_hoff_cond_draws(zs)
    assert all(a >= b for a, b in zip(hc, hb))
    gap = max(gap, max(a - b for a, b in zip(hc, hb)))
for zs in multisets(Z, 3):
    assert is_equivariant(zs, score_hoff_cond, range(3))
print(f"Slice 12, latents: (c1), the joint posterior given the others, equals (b) on all 216 data sets; "
      f"(c2), each unit's latent conditioned on its own value, is equivariant, equals its reading from "
      f"the n-unit draws, and sits above (b) unit by unit (by up to {float(gap):.4f})")
N = 4
for lname, law in (("the model's own law", model_law(N)), ("an iid law outside the model", out_law(N))):
    for k in (1, 2):
        bound = 1 - Q(k, N)
        row = {}
        for sname, sc in (("(b)", score_hoff), ("(c2)", score_hoff_cond), ("mixed", score_hoff_mixed),
                          ("reverse", score_hoff_mixed_reverse)):
            row[sname] = coverage_law(law, sc, k, N - 1, range(N))
        assert row["(b)"] >= bound and row["(c2)"] >= bound
        assert row["mixed"] < bound
        print(f"Slice 12 Hoff under {lname}, N = 4, k = {k}, bound {float(bound):.3f}: " +
              ", ".join(f"{a} {float(c):.3f}" for a, c in row.items()))

# CB-4(c) from S draws, the candidate's latent drawn per draw from p(v | theta_s, z_new). Over the
# drawn latents the estimate's mean is the Rao-Blackwell reading sum_s w_s E[f | theta_s, z_new]/W,
# and its variance is sum_s (w_s/W)^2 Var(f | theta_s, z_new): what the drawn latent adds.
mr = random.Random(1211)
n_rb = 0
for trial in range(300):
    znew = mr.choice(Z)
    S = mr.choice((1, 2, 3, 4))
    ths = [mr.choice(TH) for _ in range(S)]
    w = [g(th, znew) for th in ths]
    W = sum(w)
    cond = [[qv[th][v] * f[th][v][znew[0]][znew[1]] / g(th, znew) for v in V] for th in ths]
    fv = [[f[th][v][znew[0]][znew[1]] for v in V] for th in ths]
    mean, second = Q(0), Q(0)
    for vs in itertools.product(V, repeat=S):
        pr = prod(cond[s][vs[s]] for s in range(S))
        est = sum(w[s] * fv[s][vs[s]] for s in range(S)) / W
        mean += pr * est
        second += pr * est ** 2
    rb = [sum(cond[s][v] * fv[s][v] for v in V) for s in range(S)]
    var_s = [sum(cond[s][v] * (fv[s][v] - rb[s]) ** 2 for v in V) for s in range(S)]
    assert mean == sum(w[s] * rb[s] for s in range(S)) / W
    assert second - mean ** 2 == sum((w[s] / W) ** 2 * var_s[s] for s in range(S))
    n_rb += 1
for zs in itertools.product(Z, repeat=3):
    assert score_lat_from_draws(zs) == score_lat(zs)
print(f"Slice 12, CB-4(c) from draws: on {n_rb} draw sets the mean over the candidate's drawn latents is "
      f"the Rao-Blackwell reading and the added variance is sum_s (w_s/W)^2 Var(f | theta_s, z_new); "
      f"the weighted reading of the n-unit joint is (c) for every unit, the candidate included")


# Hoff's optimality, finite. P is the model's law of the units (exchangeable, positive). A procedure
# phi gives phi(z), the chance that z's last unit is accepted given the others; its expected size is
# sum_z P_-t(z) phi(z), with P_-t(z) = sum_y P(z with t -> y) the others' marginal. Among procedures
# whose coverage over the orderings of every data set is at least that of Hoff's set, Hoff's set
# has the smallest expected size; with ties broken at random to exactly 1 - k/m, it is smallest
# among all procedures valid at 1 - k/m.
def opt_check(N, k, trials, seed):
    C, t = list(range(N)), N - 1
    P = model_law(N)
    G = perms_on(N, C)

    def Pm(zs, i):
        return sum(P[zs[:i] + (y,) + zs[i + 1:]] for y in Z)

    hoff_cache = {zs: score_hoff(zs) for zs in P}
    for zs in itertools.islice(P, 0, None, 37):
        assert all(P[zs] / Pm(zs, i) == hoff_cache[zs][i] / 2 for i in range(N))

    def size(phi):
        return sum(Pm(zs, t) * phi[zs] for zs in P)

    def conformal(score):
        return {zs: (1 if accepts(k, score(zs), t, C) else 0) for zs in P}

    def randomized(score):
        """CB-1's rule with ties broken at random: exactly 1 - k/m on every data set."""
        out = {}
        for zs in P:
            c = score(zs)
            L = sum(1 for u in C if c[u] < c[t])
            T = sum(1 for u in C if c[u] == c[t])
            out[zs] = Q(1) if L >= k else (Q(0) if L + T <= k else Q(L + T - k, T))
        return out

    orbits = {}
    for zs in P:
        orbits.setdefault(tuple(sorted(zs)), []).append(zs)

    def orbit_cov(phi, key):
        # every ordering of the data set, each member weighted by how many orderings land on it
        mem = orbits[key]
        return sum(phi[zs] for zs in mem) * Q(len(G), len(mem))

    hoff, hoff_r = conformal(hoff_cache.get), randomized(hoff_cache.get)
    for key in orbits:
        assert orbit_cov(hoff_r, key) == len(G) * Q(N - k, N)
    det, rnd, short = {"Hoff": size(hoff)}, {"Hoff": size(hoff_r)}, {}
    for sname, sc in (("in-sample", score_int), ("(c)", score_lat), ("(c2)", score_hoff_cond)):
        phi, phr = conformal(sc), randomized(sc)
        det[sname], rnd[sname] = size(phi), size(phr)
        short[sname] = sum(1 for key in orbits if orbit_cov(phi, key) < orbit_cov(hoff, key))
        if short[sname] == 0:
            assert size(hoff) <= size(phi)
        assert size(hoff_r) <= size(phr) and size(hoff_r) <= size(phi)
    rr = random.Random(seed)
    for trial in range(trials):
        for target, ref in ((hoff, "cb1"), (hoff_r, "exact")):
            phi = {zs: rr.choice((Q(0), Q(1, 2), Q(1))) for zs in P}
            for key, mem in orbits.items():
                need = orbit_cov(target, key)
                order = mem[:]
                rr.shuffle(order)
                for zs in order:
                    if orbit_cov(phi, key) >= need:
                        break
                    phi[zs] = Q(1)
                assert orbit_cov(phi, key) >= need
            assert size(target) <= size(phi)
    return det, rnd, short, len(orbits)


for N, k, trials in ((3, 1, 40), (4, 1, 8), (4, 2, 8)):
    det, rnd, short, n_orb = opt_check(N, k, trials, 1212 + N + k)
    print(f"Slice 12, Hoff's optimality at m = {N}, k = {k}, expected sizes counting (x, y) values. "
          f"Ties broken at random (every set exactly 1 - k/m): " +
          ", ".join(f"{a} {float(b):.4f}" for a, b in rnd.items()) +
          ". CB-1's tie rule: " + ", ".join(f"{a} {float(b):.4f}" for a, b in det.items()) +
          f", the others covering less than Hoff's on " +
          ", ".join(f"{b} ({a})" for a, b in short.items()) + f" of {n_orb} data sets. Hoff's set no "
          f"larger than any of {2 * trials} random procedures covering as much on every data set")


print()
print("=== Slice 13: the tournament correction ===")
import math

# (a) The index: their rule refuses when the test wins at least (1 - alpha)(n + 1) games, a win
# being strictly higher nonconformity; that is the cross-fitted rule at k = floor(alpha m), m = n + 1.
rr = random.Random(1301)
n_idx = 0
for trial in range(5000):
    m = rr.randint(1, 40)
    alpha = Q(rr.randint(0, 200), 200)
    c = rr.randint(0, m + 2)
    assert (c < (1 - alpha) * m) == (c + math.floor(alpha * m) < m)
    n_idx += 1


# (b) Coverage over every ordering, with the pair fit a fit to the rest plus the pair as an
# unordered multiset. The fit reads its training units in index order (a weighted mean with
# weights by index rank: order-reading, not symmetric), as a sampler sweeping in order does.
def fit_ordered(x, T):
    T = sorted(T)
    w = [Q(j + 1) for j in range(len(T))]
    return sum(wj * x[i] for wj, i in zip(w, T)) / sum(w) if T else Q(0)


def tourn_scores(x, g):
    m = len(x)
    S = {}
    for i in range(m):
        for j in range(m):
            if i != j:
                T = [u for u in range(m) if u not in (i, j)]
                S[(i, j)] = g(fit_ordered(x, T), x[i], tuple(sorted((x[i], x[j]))))
    return S


def accepts_pairwise_py(k, S, u, m, delta=Q(0)):
    stranger = sum(1 for j in range(m) if j != u and S[(u, j)] + delta < S[(j, u)])
    return stranger + k < m


def g_bayesish(b, p, pair):
    # conformity: closeness to the fit, with the pair's spread as an add-one-in style weight
    return -abs(p - b) * (1 + Q(abs(pair[0] - pair[1]), 5))


n_tour, worst = 0, Q(1)
for m in (3, 4, 5):
    for z in multisets(range(4), m):
        for k in range(0, m // 2 + 1):
            G = list(itertools.permutations(range(m)))
            acc = Q(sum(1 for p in G if accepts_pairwise_py(k, tourn_scores(tuple(z[q] for q in p), g_bayesish), 0, m)),
                    len(G))
            bound = 1 - Q(max(2 * k - 1, 0), m)
            assert acc >= bound
            worst = min(worst, acc - bound)
            n_tour += 1


# (c) The pair fit must not know which unit is the test: one that marks the test down (here it
# always scores the test unit, unit 0, lower) refuses it on every data set.
m, k = 4, 1
refused = 0
for z in multisets(range(3), m):
    x = tuple(Q(v) for v in z)
    S = {}
    for i in range(m):
        for j in range(m):
            if i != j:
                T = [u for u in range(m) if u not in (i, j)]
                S[(i, j)] = -abs(x[i] - fit_ordered(x, T)) - (10 if i == 0 else 0)
    refused += not accepts_pairwise_py(k, S, 0, m)
assert refused == len(multisets(range(3), m))

# (d) Stability, per data set and with any reference vector R: the 2-eps-inflated rule refuses at
# most (k + a) + 2F/(a + 1) units for every a, F the ordered pairs with |S(u, j) - R(u)| > eps.
n_stab = 0
for trial in range(3000):
    m = rr.randint(2, 6)
    S = {(u, j): Q(rr.randint(-6, 6), 4) for u in range(m) for j in range(m) if u != j}
    R = [Q(rr.randint(-6, 6), 4) for _ in range(m)]
    eps = Q(rr.randint(0, 4), 8)
    F = sum(1 for (u, j), v in S.items() if abs(v - R[u]) > eps)
    for k in range(0, m + 1):
        refused = sum(1 for u in range(m) if not accepts_pairwise_py(k, S, u, m, 2 * eps))
        for a in range(0, m + 1):
            assert refused <= (k + a) + Q(2 * F, a + 1)
    n_stab += 1
# their constant against the optimized one: miss <= k/m + min_a [a/m + 2 (m - 1) nu/(a + 1)]
rows = []
for m, nu in ((20, 0.01), (100, 0.001), (100, 0.01), (1000, 0.0001)):
    best = min(a / m + 2 * (m - 1) * nu / (a + 1) for a in range(0, m + 1))
    rows.append(f"m {m}, nu {nu}: {best:.4f} against 3 sqrt(nu) = {3 * math.sqrt(nu):.4f}")
print(f"Slice 13: the index identity on {n_idx} random (m, alpha, count); the tournament covers >= "
      f"1 - (2k - 1)/m over every ordering on {n_tour} (data, k) cases with an order-reading pair fit "
      f"(least slack {worst}); a pair fit that knows the test refuses it on every data set; the "
      f"inflated rule's count bound on {n_stab} random matrices with arbitrary references; optimized "
      f"stability cost: " + "; ".join(rows))


print()
print("All assertions hold.")
