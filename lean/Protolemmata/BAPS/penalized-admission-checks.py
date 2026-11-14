"""Numeric checks for penalized-admission-lemmas.md (PA-5 to PA-7), labelled by item.

Checks, not proofs. Exact rational arithmetic (fractions) for the quadratic-form lines and the kernel claims; seeded
random families with rank-deficient S_k, with and without a common null vector; small banded families built from a
cubic B-spline Gram and second-difference penalties. Each line prints PASS or FAIL with its residual or reading; the run
ends with a count. Needs numpy and scipy.

Run: python -I -X utf8 penalized-admission-checks.py
"""
import math
from fractions import Fraction as Fr

import numpy as np
from scipy import optimize
from scipy.interpolate import BSpline

rng = np.random.default_rng(20261114)
RESULTS = []


def report(label, reading, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {label}: {reading}")
    RESULTS.append((label, bool(ok)))


# ---------------------------------------------------------------- exact rational matrices (lists of lists)

def zeros(n, m):
    return [[Fr(0)] * m for _ in range(n)]


def T(a):
    return [list(r) for r in zip(*a)]


def mm(a, b):
    bt = T(b)
    return [[sum((x * y for x, y in zip(r, c)), Fr(0)) for c in bt] for r in a]


def mv(a, x):
    return [sum((y * z for y, z in zip(r, x)), Fr(0)) for r in a]


def qf(a, x):
    return sum((xi * yi for xi, yi in zip(x, mv(a, x))), Fr(0))


def comb(c, S):
    n = len(S[0])
    A = zeros(n, n)
    for ck, Sk in zip(c, S):
        for i in range(n):
            for j in range(n):
                A[i][j] += ck * Sk[i][j]
    return A


def gram(G):
    return mm(T(G), G)


def rank_null(a):
    """Exact rank and a basis of the null space, by reduced row echelon form."""
    m = [list(r) for r in a]
    rows, cols = len(m), len(m[0])
    piv, r = [], 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        pv = m[r][c]
        m[r] = [v / pv for v in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [vi - f * vr for vi, vr in zip(m[i], m[r])]
        piv.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in piv]
    basis = []
    for f in free:
        v = [Fr(0)] * cols
        v[f] = Fr(1)
        for i, pc in enumerate(piv):
            v[pc] = -m[i][f]
        basis.append(v)
    return len(piv), basis


def is_pd(a):
    """Exact test: every pivot of the symmetric elimination without pivoting is positive."""
    m = [list(r) for r in a]
    n = len(m)
    for k in range(n):
        if m[k][k] <= 0:
            return False
        for i in range(k + 1, n):
            f = m[i][k] / m[k][k]
            for j in range(k, n):
                m[i][j] -= f * m[k][j]
    return True


def rand_fr(lo=1, hi=40):
    return Fr(int(rng.integers(lo, hi)), int(rng.integers(lo, hi)))


def rand_vec(n, lo=-6, hi=7):
    return [Fr(int(v)) for v in rng.integers(lo, hi, n)]


def family_exact(n, K, common):
    """K rank-deficient PSD S_k = G_k^T G_k with integer G_k of r_k < n rows. With `common`, every row is made
    orthogonal to one integer vector v, which then lies in every kernel; without, the stack has full column rank."""
    while True:
        v = rand_vec(n, -3, 4)
        if all(x == 0 for x in v):
            continue
        Gs = []
        for _ in range(K):
            r = int(rng.integers(1, n))
            G = [rand_vec(n, -3, 4) for _ in range(r)]
            if common:
                vv = sum(x * x for x in v)
                G = [[vv * g - sum(a * b for a, b in zip(row, v)) * vi for g, vi in zip(row, v)] for row in G]
            Gs.append(G)
        S = [gram(G) for G in Gs]
        stack_rank, _ = rank_null([row for Sk in S for row in Sk])
        if all(rank_null(Sk)[0] < n for Sk in S) and (stack_rank < n) == common:
            return S, (v if common else None)


# ---------------------------------------------------------------- PA-5: identification does not read the coefficients

print("== PA-5: identification does not read the coefficients")
DIAG = [[[Fr(1), Fr(0)], [Fr(0), Fr(0)]], [[Fr(0), Fr(0)], [Fr(0), Fr(1)]]]
fams = [(DIAG, None, "diag(1,0), diag(0,1)")]
for t in range(8):
    n, K = int(rng.integers(3, 6)), int(rng.integers(2, 4))
    S, v = family_exact(n, K, common=(t % 2 == 0))
    fams.append((S, v, f"n={n}, K={K}, {'common null' if v else 'no common null'}"))

bad_split = bad_kernel = bad_zero = bad_pos = 0
pd_rows = []
for S, v, name in fams:
    n = len(S[0])
    stack_rank, common_null = rank_null([row for Sk in S for row in Sk])
    flags = []
    for _ in range(5):
        c = [rand_fr() for _ in S]
        A = comb(c, S)
        # (i) the form splits as a nonnegative sum, and vanishes exactly on the common kernel
        for _ in range(4):
            x = rand_vec(n)
            if qf(A, x) != sum((ck * qf(Sk, x) for ck, Sk in zip(c, S)), Fr(0)):
                bad_split += 1
            in_common = all(all(y == 0 for y in mv(Sk, x)) for Sk in S)
            if (qf(A, x) == 0) != in_common:
                bad_pos += 1
        rkA, nullA = rank_null(A)
        if rkA != stack_rank or any(any(y != 0 for y in mv(Sk, w)) for Sk in S for w in nullA):
            bad_kernel += 1
        for w in common_null:
            if qf(A, w) != 0 or any(y != 0 for y in mv(A, w)):
                bad_zero += 1
        flags.append(is_pd(A))
    pd_rows.append((name, flags, stack_rank == n))
report("PA-5(i) xᵀA(c)x = Σ cₖ xᵀSₖx, exact, 9 families × 5 c × 4 x", f"{bad_split} mismatches", bad_split == 0)
report("PA-5(i) xᵀA(c)x = 0 exactly when Sₖx = 0 for every k, 180 draws", f"{bad_pos} mismatches", bad_pos == 0)
report("PA-5(i) ker A(c) = ⋂ₖ ker Sₖ: equal dimension and Sₖ kills ker A(c), 45 (family, c)", f"{bad_kernel} mismatches",
       bad_kernel == 0)
report("PA-5(i) the common null vectors are null for A(c) at every c", f"{bad_zero} mismatches", bad_zero == 0)
agree = all(len(set(f)) == 1 and f[0] == full for _, f, full in pd_rows)
report("PA-5(ii) A(c) PD at one c ⇔ at all 5 c ⇔ no common null vector",
       "; ".join(f"{nm}: PD {f[0]}" for nm, f, _ in pd_rows), agree)


# PA-5(iii): BAPS's instance with exact hat-function designs and second-difference penalties
def hat_design(points, n):
    return [[max(Fr(0), 1 - abs(x - j)) for j in range(n)] for x in points]


def d2(n):
    return [[Fr(1) if j == i else Fr(-2) if j == i + 1 else Fr(1) if j == i + 2 else Fr(0) for j in range(n)]
            for i in range(n - 2)]


n = 7
D = d2(n)
P = gram(D)
P1 = gram([row for row in D if all(row[j] == 0 for j in range(4, n))])
P2 = gram([row for row in D if all(row[j] == 0 for j in range(0, 3))])
rows_ok = True
cases = []
for pts, nm in [([Fr(5, 2)] * 3, "data at one point"), ([Fr(1, 2), Fr(9, 2)], "two points"),
                ([Fr(1, 3), Fr(5, 4)], "two points in the first block")]:
    Z = hat_design(pts, n)
    ZZ = gram(Z)
    for pen_name, Ps in [("ZᵀZ + λP", [P]), ("ZᵀZ + λ₁P₁ + λ₂P₂", [P1, P2]),
                         ("ZᵀZ + Σᵢ λᵢdᵢdᵢᵀ", [gram([row]) for row in D])]:
        stack = Z + [row for Pj in Ps for row in Pj]
        rk, _ = rank_null(stack)
        common = rk < n
        sing = [not is_pd(comb([Fr(1)] + [rand_fr() for _ in Ps], [ZZ] + Ps)) for _ in range(4)]
        ok = len(set(sing)) == 1 and sing[0] == common
        rows_ok &= ok
        cases.append(f"{nm}, {pen_name}: singular {sing[0]}, common null {common}")
report("PA-5(iii) ZᵀZ + Σⱼ λⱼPⱼ singular ⇔ some x ≠ 0 has Zx = 0 and Pⱼx = 0, hat design, n = 7",
       "; ".join(cases), rows_ok)

# ---------------------------------------------------------------- PA-6: the bound

print("== PA-6: the bound")
bad = 0
draws = 0
for S, v, name in fams[1:]:
    for _ in range(4):
        c = [rand_fr() for _ in S]
        cp = [rand_fr() for _ in S]
        rho = [a / b for a, b in zip(c, cp)]
        m_, M_ = min(rho), max(rho)
        A, Ap = comb(c, S), comb(cp, S)
        for _ in range(5):
            x = rand_vec(len(S[0]))
            draws += 1
            if not (m_ * qf(Ap, x) <= qf(A, x) <= M_ * qf(Ap, x)):
                bad += 1
report(f"PA-6(i) m·xᵀA(c′)x ≤ xᵀA(c)x ≤ M·xᵀA(c′)x, exact, {draws} draws (rank-deficient Sₖ, both kinds)",
       f"{bad} violations", bad == 0)


def lam(A):
    w = np.linalg.eigvalsh(A)
    return w[0], w[-1]


worst = {"min": 0.0, "max": 0.0, "kappa": 0.0}
viol = 0
for t in range(300):
    n = int(rng.integers(4, 16))
    K = int(rng.integers(2, 6))
    Ss = []
    while True:
        Ss = []
        for _ in range(K):
            r = int(rng.integers(1, n))
            G = rng.standard_normal((r, n))
            Ss.append(G.T @ G)
        if np.linalg.matrix_rank(np.vstack(Ss)) == n:
            break
    c = np.exp(rng.uniform(-7, 7, K))
    cp = np.exp(rng.uniform(-7, 7, K))
    rho = c / cp
    m_, M_ = rho.min(), rho.max()
    A = sum(ck * Sk for ck, Sk in zip(c, Ss))
    Ap = sum(ck * Sk for ck, Sk in zip(cp, Ss))
    lmin, lmax = lam(A)
    lpmin, lpmax = lam(Ap)
    if lpmin <= 0:
        continue
    r1 = (m_ * lpmin) / lmin
    r2 = lmax / (M_ * lpmax)
    r3 = (lmax / lmin) / ((M_ / m_) * (lpmax / lpmin))
    worst["min"] = max(worst["min"], r1)
    worst["max"] = max(worst["max"], r2)
    worst["kappa"] = max(worst["kappa"], r3)
    viol += (r1 > 1 + 1e-6) + (r2 > 1 + 1e-9) + (r3 > 1 + 1e-6)
report("PA-6(ii) λ_min(A(c)) ≥ m·λ_min(A(c′)), λ_max ≤ M·λ_max, κ₂ ≤ (M/m)·κ₂′, 300 seeded families",
       f"largest readings (bound = 1): λ_min {worst['min']:.6f}, λ_max {worst['max']:.6f}, κ₂ {worst['kappa']:.6f}",
       viol == 0)

ok = True
reads = []
for c1, c2 in [(Fr(3), Fr(1, 5)), (Fr(1, 7), Fr(11, 2)), (Fr(2), Fr(2))]:
    A = comb([c1, c2], DIAG)
    diagonal = A[0][1] == 0 and A[1][0] == 0  # so A's eigenvalues are its diagonal entries
    kappa = max(A[0][0], A[1][1]) / min(A[0][0], A[1][1])
    rho = [c1 / 1, c2 / 1]
    ok &= diagonal and kappa == max(rho) / min(rho)
    reads.append(f"c = ({c1}, {c2}): κ₂ = {kappa}")
report("PA-6(iii) M/m attained: S₁ = diag(1,0), S₂ = diag(0,1), c′ = (1,1), κ₂(A(c′)) = 1, κ₂(A(c)) = M/m",
       "; ".join(reads), ok)

S, _ = family_exact(4, 3, common=False)
c = [Fr(1), rand_fr(), rand_fr()]
cp = [Fr(1), rand_fr(), rand_fr()]
rho = [a / b for a, b in zip(c, cp)]
report("PA-6(iv) a coefficient held at 1 in both enters with ratio 1, so m ≤ 1 ≤ M",
       f"ρ = ({', '.join(str(r) for r in rho)})", min(rho) <= 1 <= max(rho))

# ---------------------------------------------------------------- PA-7: from the bound to what the gate reads

print("== PA-7: from the bound to what the gate reads")


def unit_scale(A):
    d = 1 / np.sqrt(np.diag(A))
    return A * np.outer(d, d)


def norm1(A):
    return np.abs(A).sum(axis=0).max()


def kappa2(A):
    lmin, lmax = lam(A)
    return lmax / lmin


def kappa1(A):
    return norm1(A) * norm1(np.linalg.inv(A))


def nnz_rows(A):
    return int((A != 0).sum(axis=1).max())


def banded_pd(n, w, delta):
    Tm = np.zeros((n, n))
    for k in range(w + 1):
        vals = rng.standard_normal(n - k)
        Tm += np.diag(vals, k)
        if k:
            Tm += np.diag(vals, -k)
    Tm += (delta - np.linalg.eigvalsh(Tm)[0]) * np.eye(n)
    return unit_scale(Tm)


steps_bad = 0
worst_ratio = 0.0
for t in range(400):
    n = int(rng.integers(3, 14))
    w = int(rng.integers(0, min(4, n)))
    B = banded_pd(n, w, 10.0 ** rng.uniform(-6, 0))
    q = nnz_rows(B)
    f = np.exp(rng.uniform(-7, 7, n))
    FBF = B * np.outer(f, f)
    bmin, bmax = lam(B)
    fmin, fmax = lam(FBF)
    steps_bad += np.abs(B).max() > 1 + 1e-12
    steps_bad += bmax > q * (1 + 1e-12)
    steps_bad += fmax < (f ** 2).max() * (1 - 1e-9)
    steps_bad += bmin < fmin / (f ** 2).max() * (1 - 1e-6)
    ratio = kappa2(B) / (q * kappa2(FBF))
    worst_ratio = max(worst_ratio, ratio)
report("PA-7(i) sketch: |Bᵢⱼ| ≤ 1, λ_max(B) ≤ q, λ_max(FBF) ≥ max fᵢ², λ_min(B) ≥ λ_min(FBF)/max fᵢ², 400 draws",
       f"{steps_bad} violations", steps_bad == 0)
report("PA-7(i) κ₂(B) ≤ q·κ₂(FBF), unit-diagonal banded PD B, F log-uniform over e^±7",
       f"largest κ₂(B)/(q·κ₂(FBF)) = {worst_ratio:.4f}", worst_ratio <= 1 + 1e-6)

# how much of q the unit diagonal leaves unused: κ₂(B) against the best diagonal scaling found by local search
reads = []
for n, w in [(6, 1), (8, 2), (10, 3)]:
    B = banded_pd(n, w, 1e-3)
    q = nnz_rows(B)
    best = optimize.minimize(lambda g: np.log(kappa2(B * np.outer(np.exp(g), np.exp(g)))), np.zeros(n),
                             method="Nelder-Mead", options={"maxiter": 4000, "xatol": 1e-8, "fatol": 1e-10})
    reads.append(f"n={n}, q={q}: κ₂(B)/min_F κ₂(FBF) ≈ {kappa2(B) / math.exp(best.fun):.3f}")
report("PA-7(i) reading: the unit diagonal against a local search over F (≤ q by van der Sluis 4.3)",
       "; ".join(reads), True)

bad = [0, 0, 0]
for t in range(300):
    n = int(rng.integers(2, 20))
    A = rng.standard_normal((n, n))
    As = A + A.T
    bad[0] += np.linalg.norm(As, 2) > norm1(As) * (1 + 1e-12)
    bad[2] += norm1(A) > math.sqrt(n) * np.linalg.norm(A, 2) * (1 + 1e-12)
    q = int(rng.integers(1, n + 1))
    Sp = np.zeros((n, n))
    for j in range(n):
        idx = rng.choice(n, size=q, replace=False)
        Sp[idx, j] = rng.standard_normal(q)
    bad[1] += norm1(Sp) > math.sqrt(q) * np.linalg.norm(Sp, 2) * (1 + 1e-12)
report("PA-7(ii) ‖A‖₂ ≤ ‖A‖₁ symmetric; ‖A‖₁ ≤ √q‖A‖₂ with q nonzeros per column; ‖A‖₁ ≤ √n‖A‖₂, 300 each",
       f"violations {[int(b) for b in bad]}", sum(bad) == 0)
col = np.zeros((9, 9))
col[:, 0] = 1
report("PA-7(ii) √n attained by a column of ones (n = 9)",
       f"‖A‖₁/‖A‖₂ = {norm1(col) / np.linalg.norm(col, 2):.6f}", abs(norm1(col) / np.linalg.norm(col, 2) - 3) < 1e-12)


# composed, on BAPS's banded families
def bspline_design(x, n, k=3):
    t = np.r_[[0.0] * k, np.linspace(0, 1, n - k + 1), [1.0] * k]
    return BSpline.design_matrix(x, t, k).toarray()


def d2f(n):
    D = np.zeros((n - 2, n))
    for i in range(n - 2):
        D[i, i:i + 3] = [1, -2, 1]
    return D


chain_bad = 0
worst = {"posed": 0.0, "sharp": 0.0, "inv": 0.0, "plain": 0.0}
counts = {"ZᵀZ + λP": 0, "ZᵀZ + Σⱼ λⱼPⱼ": 0, "ZᵀZ + Σᵢ λᵢdᵢdᵢᵀ": 0}
for t in range(240):
    n = int(rng.integers(6, 30))
    mpts = int(rng.integers(2, 3 * n))
    x = np.sort(rng.uniform(0, 1, mpts) ** rng.uniform(0.3, 3))
    Z = bspline_design(x, n)
    D = d2f(n)
    kind = list(counts)[t % 3]
    if kind == "ZᵀZ + λP":
        S = [Z.T @ Z, D.T @ D]
    elif kind == "ZᵀZ + Σⱼ λⱼPⱼ":
        J = int(rng.integers(2, 4))
        cuts = np.array_split(np.arange(n - 2), J)
        S = [Z.T @ Z] + [D[cut].T @ D[cut] for cut in cuts]
    else:
        S = [Z.T @ Z] + [np.outer(D[i], D[i]) for i in range(n - 2)]
    K = len(S)
    c = np.r_[1.0, np.exp(rng.uniform(-7, 7, K - 1))]
    cr = np.r_[1.0, np.exp(rng.uniform(-7, 7, K - 1))]
    A = sum(ck * Sk for ck, Sk in zip(c, S))
    Ar = sum(ck * Sk for ck, Sk in zip(cr, S))
    if lam(Ar)[0] <= 1e-12 * lam(Ar)[1] or lam(A)[0] <= 1e-12 * lam(A)[1]:
        continue
    counts[kind] += 1
    pattern = np.zeros((n, n), bool)
    for Sk in S:
        pattern |= Sk != 0
    q = int(pattern.sum(axis=1).max())
    rho = c / cr
    Mm = rho.max() / rho.min()
    E, Er = unit_scale(A), unit_scale(Ar)
    k1, k1r = kappa1(E), kappa1(Er)
    inv_r = norm1(np.linalg.inv(Er))
    posed = q * math.sqrt(q * n) * Mm * k1r
    sharp = q * math.sqrt(n) * Mm * k1r
    via_inv = q * math.sqrt(n) * Mm * inv_r
    chain_bad += (k1 > posed * (1 + 1e-6)) + (k1 > sharp * (1 + 1e-6)) + (k1 > via_inv * (1 + 1e-6))
    # the pieces the sharper chain uses
    chain_bad += norm1(E) > q * (1 + 1e-12)
    chain_bad += norm1(np.linalg.inv(E)) > math.sqrt(n) / lam(E)[0] * (1 + 1e-6)
    chain_bad += lam(E)[0] < lam(Er)[0] / Mm * (1 - 1e-6)
    chain_bad += (np.diag(A) / np.diag(Ar)).max() > rho.max() * (1 + 1e-12)
    worst["posed"] = max(worst["posed"], k1 / posed)
    worst["sharp"] = max(worst["sharp"], k1 / sharp)
    worst["inv"] = max(worst["inv"], k1 / via_inv)
    worst["plain"] = max(worst["plain"], k1 / (Mm * k1r))
report(f"PA-7(iii) κ₁(E(c)) ≤ q·√(qn)·(M/m)·κ₁(E(c_ref)) and the sharper forms, B-spline Gram + D₂ penalties {counts}",
       f"{chain_bad} violations; largest κ₁(E)/bound: posed {worst['posed']:.2e}, q·√n {worst['sharp']:.2e}, "
       f"q·√n·‖E_ref⁻¹‖₁ {worst['inv']:.2e}", chain_bad == 0)
report("PA-7(iii) reading: largest κ₁(E(c)) / ((M/m)·κ₁(E(c_ref))), the share of the proved factor in use",
       f"{worst['plain']:.3f}", True)

# ---------------------------------------------------------------- summary
fails = [lbl for lbl, ok in RESULTS if not ok]
print(f"\n{len(RESULTS) - len(fails)} of {len(RESULTS)} PASS")
for lbl in fails:
    print("FAIL:", lbl)
