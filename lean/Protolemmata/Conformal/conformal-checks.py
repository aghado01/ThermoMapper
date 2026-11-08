"""Numeric checks of conformal-lemmas.md (CB-1 to CB-4, M-1 and M-2), labelled by item, before any Lean.

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
print("All assertions hold.")
