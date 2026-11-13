"""Numeric checks for active-subspace-lemmas.md (AS-1 to AS-4), labelled by item.

Checks, not proofs. Exact rational arithmetic (fractions) where the claim is algebra on rationals, small cases
enumerated, and seeded random cases beside chosen ones. Each line prints PASS or FAIL with its residual or reading;
the run ends with a count. Needs numpy and scipy.

Run: python -I -X utf8 active-subspace-checks.py
"""
import itertools
import math
from fractions import Fraction as Fr

import numpy as np
from numpy.polynomial import legendre as leg
from scipy import special, stats

rng = np.random.default_rng(20261113)
RESULTS = []


def report(label, reading, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {label}: {reading}")
    RESULTS.append((label, bool(ok)))


# ---------------------------------------------------------------- exact rational matrices (lists of lists)

def mat(rows):
    return [[Fr(x) for x in r] for r in rows]


def T(a):
    return [list(r) for r in zip(*a)]


def mm(a, b):
    bt = T(b)
    return [[sum((x * y for x, y in zip(r, c)), Fr(0)) for c in bt] for r in a]


def add(a, b, s=Fr(1)):
    return [[x + s * y for x, y in zip(r, q)] for r, q in zip(a, b)]


def scale(a, s):
    return [[s * x for x in r] for r in a]


def eye(n):
    return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]


def zeros(n, m=None):
    return [[Fr(0)] * (n if m is None else m) for _ in range(n)]


def tr(a):
    return sum((a[i][i] for i in range(len(a))), Fr(0))


def inv(a):
    n = len(a)
    m = [list(r) + e for r, e in zip(a, eye(n))]
    for c in range(n):
        p = next((r for r in range(c, n) if m[r][c] != 0), None)
        if p is None:
            return None
        m[c], m[p] = m[p], m[c]
        pv = m[c][c]
        m[c] = [x / pv for x in m[c]]
        for r in range(n):
            if r != c and m[r][c] != 0:
                f = m[r][c]
                m[r] = [x - f * y for x, y in zip(m[r], m[c])]
    return [r[n:] for r in m]


def det(a):
    n = len(a)
    m = [list(r) for r in a]
    d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if m[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            m[c], m[p] = m[p], m[c]
            d = -d
        d *= m[c][c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            m[r] = [x - f * y for x, y in zip(m[r], m[c])]
    return d


def frob2(a):
    return sum((x * x for r in a for x in r), Fr(0))


def projection_exact(basis):
    """P = A (AᵀA)⁻¹ Aᵀ for an integer basis A (p × k), or the zero matrix at k = 0."""
    p = len(basis)
    if not basis or not basis[0]:
        return zeros(p)
    a = mat(basis)
    g = inv(mm(T(a), a))
    return None if g is None else mm(mm(a, g), T(a))


def random_projection_exact(p, k):
    while True:
        a = [[int(rng.integers(-3, 4)) for _ in range(k)] for _ in range(p)]
        P = projection_exact(a)
        if P is not None:
            return P


def random_weights_exact(s):
    w = [Fr(int(rng.integers(1, 9))) for _ in range(s)]
    tot = sum(w)
    return [x / tot for x in w]


def to_np(a):
    return np.array([[float(x) for x in r] for r in a])


def orth_basis(p, k):
    q, _ = np.linalg.qr(rng.standard_normal((p, max(k, 1))))
    return q[:, :k]


def principal_cosines(U, V):
    if U.shape[1] == 0 or V.shape[1] == 0:
        return np.zeros(0)
    return np.clip(np.linalg.svd(U.T @ V, compute_uv=False), 0.0, 1.0)


# ================================================================ AS-1: P̄ over mixed ranks

def as1():
    print("\n== AS-1: P̄ over mixed ranks")
    # (a) exact: projections, P̄'s symmetry and trace, over 200 seeded posteriors of mixed ranks
    worst_sym = worst_idem = 0
    trace_ok = True
    loss_ok = True
    eig_lo, eig_hi = 1.0, 0.0
    for _ in range(200):
        p = int(rng.integers(2, 6))
        S = int(rng.integers(1, 5))
        ks = [int(rng.integers(0, p + 1)) for _ in range(S)]
        Ps = [random_projection_exact(p, k) for k in ks]
        w = random_weights_exact(S)
        for P, k in zip(Ps, ks):
            worst_sym = max(worst_sym, frob2(add(P, T(P), Fr(-1))))
            worst_idem = max(worst_idem, frob2(add(mm(P, P), P, Fr(-1))))
            trace_ok &= tr(P) == k
        Pbar = zeros(p)
        for wi, P in zip(w, Ps):
            Pbar = add(Pbar, P, wi)
        worst_sym = max(worst_sym, frob2(add(Pbar, T(Pbar), Fr(-1))))
        trace_ok &= tr(Pbar) == sum(wi * k for wi, k in zip(w, ks))
        ev = np.linalg.eigvalsh(to_np(Pbar))
        eig_lo, eig_hi = min(eig_lo, ev.min()), max(eig_hi, ev.max())
        # (b) the loss identity, exactly, at a projection of every rank
        for k in range(p + 1):
            P = random_projection_exact(p, k)
            lhs = sum((wi * frob2(add(P, Ps_i, Fr(-1))) for wi, Ps_i in zip(w, Ps)), Fr(0))
            rhs = tr(P) + sum(wi * k_i for wi, k_i in zip(w, ks)) - 2 * tr(mm(P, Pbar))
            loss_ok &= lhs == rhs
    report("AS-1a projections exact: symmetric, idempotent", f"worst ‖P−Pᵀ‖², ‖P²−P‖² = {worst_sym}, {worst_idem}",
           worst_sym == 0 and worst_idem == 0)
    report("AS-1a tr P_s = k_s and tr P̄ = Σ w_s k_s, exactly", "200 posteriors", trace_ok)
    report("AS-1a 0 ≼ P̄ ≼ I", f"eigenvalues in [{eig_lo:.3e}, {eig_hi:.15f}]", eig_lo > -1e-12 and eig_hi < 1 + 1e-12)
    report("AS-1b Σ w‖P−P_s‖² = tr P + Σ w k_s − 2 tr(P P̄), exactly", "every rank of P, 200 posteriors", loss_ok)

    # (c) the minimizer: keep P̄'s eigenvectors with eigenvalue above ½
    worst_margin = math.inf
    for _ in range(150):
        p = int(rng.integers(2, 7))
        S = int(rng.integers(1, 6))
        ks = [int(rng.integers(0, p + 1)) for _ in range(S)]
        bases = [orth_basis(p, k) for k in ks]
        w = rng.dirichlet(np.ones(S))
        Ps = [B @ B.T for B in bases]
        Pbar = sum(wi * P for wi, P in zip(w, Ps))
        mu, V = np.linalg.eigh(Pbar)

        def loss(P):
            return sum(wi * np.sum((P - Q) ** 2) for wi, Q in zip(w, Ps))

        Vs = V[:, mu > 0.5]
        best = loss(Vs @ Vs.T)
        for r in range(p + 1):
            for sub in itertools.combinations(range(p), r):
                Vsub = V[:, list(sub)]
                worst_margin = min(worst_margin, loss(Vsub @ Vsub.T) - best)
        for _ in range(60):
            B = orth_basis(p, int(rng.integers(0, p + 1)))
            worst_margin = min(worst_margin, loss(B @ B.T) - best)
    report("AS-1b the rank summary is least over projections of every rank",
           f"least margin over eigen-subsets and random projections = {worst_margin:.2e}", worst_margin > -1e-12)

    # (c') an eigenvalue of exactly ½ is indifferent: half the mass holds e1, half holds nothing
    e1 = projection_exact([[1], [0], [0]])
    Ps, w = [e1, zeros(3)], [Fr(1, 2), Fr(1, 2)]
    l_in = sum(wi * frob2(add(e1, Q, Fr(-1))) for wi, Q in zip(w, Ps))
    l_out = sum(wi * frob2(Q) for wi, Q in zip(w, Ps))
    report("AS-1b eigenvalue exactly ½: keeping it or not, the same loss", f"{l_in} and {l_out}", l_in == l_out)

    # (c'') coordinate subspaces: Frobenius loss is Hamming loss, P̄'s diagonal the inclusion probabilities
    ok = True
    for _ in range(100):
        p = int(rng.integers(2, 7))
        S = int(rng.integers(1, 6))
        gam = rng.integers(0, 2, size=(S, p))
        w = rng.dirichlet(np.ones(S))
        incl = w @ gam
        Pbar = sum(wi * np.diag(g) for wi, g in zip(w, gam))
        ok &= np.allclose(np.diag(Pbar), incl) and np.allclose(Pbar, np.diag(incl))
        g = rng.integers(0, 2, size=p)
        ok &= abs(sum(wi * np.sum((np.diag(g) - np.diag(q)) ** 2) for wi, q in zip(w, gam))
                  - sum(wi * np.sum(np.abs(g - q)) for wi, q in zip(w, gam))) < 1e-13
    report("AS-1b coordinate subspaces: P̄ = diag(inclusion probabilities), Frobenius loss = Hamming loss", "100 cases", ok)

    # (d) across ranks: ½‖P_s − P‖² = |k_s − k|/2 + Σ sin²θ over the min(k_s, k) angles
    worst = worst_yl = 0.0
    for _ in range(300):
        p = int(rng.integers(2, 9))
        k1, k2 = int(rng.integers(0, p + 1)), int(rng.integers(0, p + 1))
        U, V = orth_basis(p, k1), orth_basis(p, k2)
        c = principal_cosines(U, V)
        half = 0.5 * np.sum((U @ U.T - V @ V.T) ** 2)
        formula = abs(k1 - k2) / 2 + np.sum(1 - c ** 2)
        assert len(c) == min(k1, k2)
        worst = max(worst, abs(half - formula))
        ye_lim = abs(k1 - k2) + np.sum(1 - c ** 2)  # Ye and Lim's chordal metric, squared (their Table 3)
        worst_yl = max(worst_yl, abs((ye_lim - half) - abs(k1 - k2) / 2))
    report("AS-1c ½‖P_s−P‖²_F = |k_s−k|/2 + Σ sin²θ_i over min(k_s,k) angles", f"worst {worst:.2e} over 300 pairs",
           worst < 1e-12)
    report("AS-1c Ye–Lim's chordal metric² exceeds ½‖P_s−P‖²_F by exactly |k_s−k|/2", f"worst {worst_yl:.2e}",
           worst_yl < 1e-12)


# ================================================================ AS-2: the matrix Poincaré check

def sup_ratio(A, B):
    """sup vᵀAv / vᵀBv for symmetric A and positive definite B, by whitening B: a law packed against an end makes
    B badly conditioned, and solve(B, A) then reads spurious ratios above the truth (seen at 1.9 here)."""
    s, U = np.linalg.eigh(B)
    keep = s > s.max() * 1e-13
    Wm = U[:, keep] / np.sqrt(s[keep])
    return float(np.linalg.eigvalsh(Wm.T @ A @ Wm).max())


def law_product(spaces):
    """All points of a finite product law with their weights."""
    pts = list(itertools.product(*[range(len(w)) for _, w in spaces]))
    wts = [math.prod((spaces[i][1][x[i]] for i in range(len(x))), start=Fr(1)) for x in pts]
    return pts, wts


def var_exact(vals, wts):
    m = sum((w * v for v, w in zip(vals, wts)), Fr(0))
    return sum((w * (v - m) ** 2 for v, w in zip(vals, wts)), Fr(0))


def as2():
    print("\n== AS-2: the matrix Poincaré check")
    # (a) the reduction: vᵀ Cov(f) v = Var(vᵀf) and vᵀ E[J D Jᵀ] v = Σ_i D_i E[(∂_i vᵀf)²], exactly
    ok = True
    for _ in range(100):
        n, m, p = int(rng.integers(2, 6)), int(rng.integers(1, 4)), int(rng.integers(1, 4))
        wts = random_weights_exact(n)
        f = [[Fr(int(rng.integers(-5, 6))) for _ in range(m)] for _ in range(n)]
        J = [[[Fr(int(rng.integers(-5, 6))) for _ in range(p)] for _ in range(m)] for _ in range(n)]
        D = [Fr(int(rng.integers(0, 6)), int(rng.integers(1, 4))) for _ in range(p)]
        mean = [sum(w * f[a][j] for a, w in enumerate(wts)) for j in range(m)]
        cov = [[sum(w * (f[a][i] - mean[i]) * (f[a][j] - mean[j]) for a, w in enumerate(wts)) for j in range(m)]
               for i in range(m)]
        G = [[sum(w * sum(J[a][i][q] * D[q] * J[a][j][q] for q in range(p)) for a, w in enumerate(wts))
              for j in range(m)] for i in range(m)]
        v = [Fr(int(rng.integers(-4, 5))) for _ in range(m)]
        qa = sum(v[i] * cov[i][j] * v[j] for i in range(m) for j in range(m))
        qb = sum(v[i] * G[i][j] * v[j] for i in range(m) for j in range(m))
        ok &= qa == var_exact([sum(v[j] * f[a][j] for j in range(m)) for a in range(n)], wts)
        ok &= qb == sum(D[q] * sum(w * sum(v[j] * J[a][j][q] for j in range(m)) ** 2 for a, w in enumerate(wts))
                        for q in range(p))
    report("AS-2a the two quadratic forms, exactly", "100 finite laws", ok)

    # (b) Efron–Stein on finite product spaces, exactly; equality for additive g
    es_ok, add_ok, strict = True, True, 0
    for _ in range(400):
        d = int(rng.integers(1, 5))
        spaces = []
        for _ in range(d):
            sz = int(rng.integers(2, 4))
            spaces.append((None, random_weights_exact(sz)))
        pts, wts = law_product(spaces)
        additive = rng.random() < 0.3
        if additive:
            parts = [[Fr(int(rng.integers(-5, 6))) for _ in range(len(w))] for _, w in spaces]
            g = {x: sum(parts[i][x[i]] for i in range(d)) for x in pts}
        else:
            g = {x: Fr(int(rng.integers(-9, 10))) for x in pts}
        lhs = var_exact([g[x] for x in pts], wts)
        rhs = Fr(0)
        for i in range(d):
            wi = spaces[i][1]
            for x, w in zip(pts, wts):
                vals = [g[x[:i] + (y,) + x[i + 1:]] for y in range(len(wi))]
                rhs += w * var_exact(vals, wi)
        es_ok &= lhs <= rhs
        if additive:
            add_ok &= lhs == rhs
        elif lhs < rhs:
            strict += 1
    report("AS-2b Efron–Stein Var g ≤ Σ_i E[Var_i g], exactly", "400 product laws, 1 to 4 coordinates", es_ok)
    report("AS-2b equality for additive g, exactly", f"strict on {strict} non-additive cases", add_ok)

    # (c) the uniform's constant (b − a)²/π²: Galerkin supremum, equality at cos, random smooth functions
    xg, wg = leg.leggauss(200)
    worst, eq_resid, sup_gap = 0.0, 0.0, 0.0
    for case in range(300):
        a = rng.uniform(-3, 3)
        b = a + rng.uniform(0.1, 5)
        x = a + (b - a) * (xg + 1) / 2
        w = wg / 2  # the uniform law's weights
        c = (b - a) ** 2 / math.pi ** 2
        if case == 0:
            h, dh = np.cos(math.pi * (x - a) / (b - a)), -math.pi / (b - a) * np.sin(math.pi * (x - a) / (b - a))
            eq_resid = abs(np.sum(w * (h - np.sum(w * h)) ** 2) / (c * np.sum(w * dh ** 2)) - 1)
        deg = int(rng.integers(1, 15))
        coef = rng.standard_normal(deg + 1) / (1 + np.arange(deg + 1))
        t = (x - a) / (b - a) * 2 - 1
        fr = rng.standard_normal()
        om = rng.uniform(0, 9)
        h = leg.legval(t, coef) + fr * np.sin(om * x)
        dh = leg.legval(t, leg.legder(coef)) * 2 / (b - a) + fr * om * np.cos(om * x)
        var = np.sum(w * (h - np.sum(w * h)) ** 2)
        en = np.sum(w * dh ** 2)
        if en > 0:
            worst = max(worst, var / (c * en))
    # the Galerkin supremum over polynomials of degree ≤ 30 approaches 1 from below
    t = xg
    V = np.stack([leg.legval(t, np.eye(31)[j]) for j in range(1, 31)], 1)
    dV = np.stack([leg.legval(t, leg.legder(np.eye(31)[j])) for j in range(1, 31)], 1) * 2  # d/dx on [0, 1]
    w = wg / 2
    Vc = V - (w @ V)
    A = (Vc * w[:, None]).T @ Vc
    B = (dV * w[:, None]).T @ dV
    sup = sup_ratio(A, B) * math.pi ** 2
    sup_gap = abs(sup - 1)
    report("AS-2c the uniform: equality at cos(π(x−a)/(b−a))", f"|ratio − 1| = {eq_resid:.2e}", eq_resid < 1e-12)
    report("AS-2c the uniform: Var h ≤ (b−a)²/π² E[h′²] on 300 random smooth h", f"largest ratio {worst:.6f}",
           worst <= 1 + 1e-12)
    report("AS-2c the uniform: Galerkin supremum over degree ≤ 30 is 1", f"|sup·π²/(b−a)² − 1| = {sup_gap:.2e}",
           sup_gap < 1e-10)

    # (d) the matrix check under a product uniform law, by tensor quadrature
    def matrix_check(funcs, grads, boxes, D):
        nodes, weights = leg.leggauss(40)
        axes = [(lo + (hi - lo) * (nodes + 1) / 2, weights / 2) for lo, hi in boxes]
        grid = np.meshgrid(*[ax[0] for ax in axes], indexing="ij")
        W = np.ones_like(grid[0])
        for i, ax in enumerate(axes):
            shape = [1] * len(axes)
            shape[i] = -1
            W = W * ax[1].reshape(shape)
        Wf = W.ravel()
        F = np.stack([np.broadcast_to(f(*grid), W.shape).ravel() for f in funcs])
        Fm = F - (F @ Wf)[:, None]
        cov = (Fm * Wf) @ Fm.T
        Jg = np.stack([np.stack([np.broadcast_to(g(*grid), W.shape).ravel() for g in row]) for row in grads])  # m × p × n
        G = np.einsum("aik,i,bik,k->ab", Jg, np.asarray(D), Jg, Wf)
        return sup_ratio(cov, G)

    pi = math.pi
    one = matrix_check([lambda x, y: np.cos(pi * x), lambda x, y: np.cos(pi * y)],
                       [[lambda x, y: -pi * np.sin(pi * x), lambda x, y: 0 * x],
                        [lambda x, y: 0 * x, lambda x, y: -pi * np.sin(pi * y)]],
                       [(0, 1), (0, 1)], [1 / pi ** 2, 1 / pi ** 2])
    half = matrix_check([lambda x, y: np.cos(pi * x), lambda x, y: np.cos(pi * y)],
                        [[lambda x, y: -pi * np.sin(pi * x), lambda x, y: 0 * x],
                         [lambda x, y: 0 * x, lambda x, y: -pi * np.sin(pi * y)]],
                        [(0, 1), (0, 1)], [0.5 / pi ** 2, 0.5 / pi ** 2])
    report("AS-2d (cos πx₁, cos πx₂) on [0,1]²: largest generalized eigenvalue 1; D halved: 2", f"{one:.15f}, {half:.12f}",
           abs(one - 1) < 1e-12 and abs(half - 2) < 1e-10)

    worst_diag, gain = 0.0, []
    for _ in range(60):
        boxes = []
        for _ in range(2):
            lo = rng.uniform(-2, 2)
            boxes.append((lo, lo + rng.uniform(0.2, 4)))
        D = [(hi - lo) ** 2 / pi ** 2 for lo, hi in boxes]
        fr = rng.uniform(0.3, 3, size=(3, 2))
        ph = rng.uniform(0, 6, size=(3, 2))
        cc = rng.standard_normal((3, 3))
        funcs = [
            (lambda x, y, r=r: cc[r, 0] * np.sin(fr[r, 0] * x + ph[r, 0]) + cc[r, 1] * np.cos(fr[r, 1] * y + ph[r, 1])
             + cc[r, 2] * np.sin(fr[r, 0] * x + ph[r, 0]) * np.cos(fr[r, 1] * y + ph[r, 1])) for r in range(3)]
        grads = [[
            (lambda x, y, r=r: cc[r, 0] * fr[r, 0] * np.cos(fr[r, 0] * x + ph[r, 0])
             + cc[r, 2] * fr[r, 0] * np.cos(fr[r, 0] * x + ph[r, 0]) * np.cos(fr[r, 1] * y + ph[r, 1])),
            (lambda x, y, r=r: -cc[r, 1] * fr[r, 1] * np.sin(fr[r, 1] * y + ph[r, 1])
             - cc[r, 2] * fr[r, 1] * np.sin(fr[r, 0] * x + ph[r, 0]) * np.sin(fr[r, 1] * y + ph[r, 1]))]
            for r in range(3)]
        rd = matrix_check(funcs, grads, boxes, D)
        rs = matrix_check(funcs, grads, boxes, [max(D)] * 2)
        worst_diag = max(worst_diag, rd)
        gain.append(rd / rs)
    report("AS-2d matrix check, 60 random three-output fields on random boxes: largest generalized eigenvalue ≤ 1",
           f"worst {worst_diag:.6f}", worst_diag <= 1 + 1e-10)
    report("AS-2d the diagonal D is never looser than the scalar max(D)",
           f"diagonal reading / scalar reading in [{min(gain):.3f}, {max(gain):.3f}]", min(gain) >= 1 - 1e-12)

    # (e) truncated normals: the Galerkin Poincaré constant against min(σ², (b−a)²/π²); a bimodal law beats both
    def poincare_constant(logdens, a, b, deg=40):
        x01, w01 = leg.leggauss(400)
        x = a + (b - a) * (x01 + 1) / 2
        lw = logdens(x)
        w = w01 * np.exp(lw - lw.max())
        w = w / w.sum()
        t = x01
        V = np.stack([leg.legval(t, np.eye(deg + 1)[j]) for j in range(1, deg + 1)], 1)
        dV = np.stack([leg.legval(t, leg.legder(np.eye(deg + 1)[j])) for j in range(1, deg + 1)], 1) * 2 / (b - a)
        Vc = V - (w @ V)
        A = (Vc * w[:, None]).T @ Vc
        B = (dV * w[:, None]).T @ dV
        return sup_ratio(A, B)

    worst_tn = 0.0
    for _ in range(200):
        a = rng.uniform(-2, 2)
        b = a + rng.uniform(0.2, 4)
        mu = rng.uniform(a - 1, b + 1)
        sd = rng.uniform(0.05, 3)
        cp = poincare_constant(lambda x: -0.5 * ((x - mu) / sd) ** 2, a, b)
        worst_tn = max(worst_tn, cp / min(sd ** 2, (b - a) ** 2 / pi ** 2))
    report("AS-2e truncated normals: Poincaré constant ≤ min(σ², (b−a)²/π²) (Bakry–Émery, Payne–Weinberger)",
           f"largest constant / bound over 200 = {worst_tn:.6f}", worst_tn <= 1 + 1e-9)
    a, b, sd = 0.0, 1.0, 0.06
    cp = poincare_constant(lambda x: np.logaddexp(-0.5 * ((x - 0.2) / sd) ** 2, -0.5 * ((x - 0.8) / sd) ** 2), a, b)
    report("AS-2e two-normal mixture on [0,1], sds 0.06: constant above both σ² and 1/π² (so no defect reading there)",
           f"constant {cp:.4f} against σ² {sd ** 2:.4f}, 1/π² {1 / pi ** 2:.4f}", cp > max(sd ** 2, 1 / pi ** 2))


# ================================================================ AS-3: E36's direction moves

def ps_log_normalizer(a, kappa):
    """log of the power spherical normalizer on S^{a−1} (De Cao and Aziz 2020, eq. for N_X), a ≥ 2."""
    al, be = (a - 1) / 2 + kappa, (a - 1) / 2
    return (al + be) * math.log(2) + be * math.log(math.pi) + special.gammaln(al) - special.gammaln(al + be)


def as3():
    print("\n== AS-3: E36's direction moves")
    # (a) the normalizer integrates the density to one, via the marginal of t = μᵀx
    worst = 0.0
    for a in (2, 3, 4, 7):
        for kappa in (0.0, 0.5, 1.0, 10.0, 1000.0):
            be = (a - 1) / 2
            # ∫_{S^{a−1}} (1 + μᵀx)^κ dσ = |S^{a−2}| ∫ (1+t)^κ (1−t²)^{(a−3)/2} dt = |S^{a−2}| 2^{κ+a−2} B(κ+β, β)
            log_area = math.log(2) + be * math.log(math.pi) - special.gammaln(be)
            log_int = (kappa + a - 2) * math.log(2) + special.betaln(kappa + be, be)
            worst = max(worst, abs(log_area + log_int - ps_log_normalizer(a, kappa)))
    report("AS-3a the power spherical normalizer depends on (κ, a) only and integrates to one", f"worst log gap {worst:.2e}",
           worst < 1e-10)
    # Monte Carlo of the paper's sampler: the cosine's mean against (α − β)/(α + β)
    worst_z = 0.0
    for a in (2, 3, 5):
        for kappa in (1.0, 50.0, 1000.0):
            al, be = (a - 1) / 2 + kappa, (a - 1) / 2
            z = rng.beta(al, be, size=200000)
            t = 2 * z - 1
            se = t.std() / math.sqrt(len(t))
            worst_z = max(worst_z, abs(t.mean() - (al - be) / (al + be)) / se)
    report("AS-3a the sampler's cosine mean against (α−β)/(α+β)", f"worst {worst_z:.2f} standard errors", worst_z < 4.5)
    # symmetry: the density is a function of μᵀx, and the normalizer does not read μ
    worst_sym = 0.0
    for a in (2, 3, 5):
        for _ in range(100):
            u, v = rng.standard_normal(a), rng.standard_normal(a)
            u, v = u / np.linalg.norm(u), v / np.linalg.norm(v)
            kappa = rng.uniform(0, 2000)
            q_uv = kappa * math.log1p(float(u @ v)) - ps_log_normalizer(a, kappa)
            q_vu = kappa * math.log1p(float(v @ u)) - ps_log_normalizer(a, kappa)
            worst_sym = max(worst_sym, abs(q_uv - q_vu))
    report("AS-3a q(θ′|θ) = q(θ|θ′) on one sphere", f"worst log gap {worst_sym:.2e}", worst_sym < 1e-9)
    report("AS-3a at a = 1 the law on S⁰ = {±μ} puts (1 − 1)^κ = 0 on −μ for κ > 0: the change move never flips a sign",
           f"(1 + μᵀ(−μ))^1000 = {(1 + (-1.0)) ** 1000}", (1 + (-1.0)) ** 1000 == 0.0)

    # (b) the change move on a finite toy: equally spaced directions on a circle, a symmetric kernel, t₀ drawn from
    #     its prior given the new direction; the ratio is the likelihood's alone. Detailed balance over all pairs.
    def change_kernel(thetas, kappa, tprior, lik, t_from_new=True):
        n, m = len(thetas), tprior.shape[1]
        k = np.array([[(1 + math.cos(a - b)) ** kappa for b in thetas] for a in thetas])
        k = k / k.sum(1, keepdims=True)
        S = n * m
        P = np.zeros((S, S))
        for i in range(n):
            for s in range(m):
                for j in range(n):
                    for r in range(m):
                        src = j if t_from_new else i
                        q = k[i, j] * tprior[src, r]
                        P[i * m + s, j * m + r] += q * min(1.0, lik[j, r] / lik[i, s])
                P[i * m + s, i * m + s] += 1 - P[i * m + s].sum()
        pi_ = (tprior * lik).ravel()  # uniform prior on θ
        return P, pi_ / pi_.sum(), k

    n, m = 9, 4
    tprior = rng.dirichlet(np.ones(m), size=n)  # π(t₀ | θ), depending on θ as the quantile bounds do
    lik = rng.uniform(0.1, 3, size=(n, m))
    grid = np.arange(n) * 2 * math.pi / n
    P, pi_, k = change_kernel(grid, 6.0, tprior, lik)
    db = np.abs(pi_[:, None] * P - (pi_[:, None] * P).T).max()
    report("AS-3b change move, equally spaced circle: detailed balance with the likelihood ratio alone",
           f"max |π P − (π P)ᵀ| = {db:.2e}; kernel symmetric to {np.abs(k - k.T).max():.1e}", db < 1e-15)
    P, pi_, _ = change_kernel(grid, 6.0, tprior, lik, t_from_new=False)
    db_bad = np.abs(pi_[:, None] * P - (pi_[:, None] * P).T).max()
    report("AS-3b t₀ drawn from its prior given the old direction breaks balance", f"max gap {db_bad:.2e}", db_bad > 1e-4)
    uneven = np.sort(rng.uniform(0, 2 * math.pi, size=n))
    P, pi_, k = change_kernel(uneven, 6.0, tprior, lik)
    db_un = np.abs(pi_[:, None] * P - (pi_[:, None] * P).T).max()
    report("AS-3b uneven circle: the row normalizer reads θ, the kernel is not symmetric, balance fails",
           f"kernel asymmetry {np.abs(k - k.T).max():.2e}, balance gap {db_un:.2e}", db_un > 1e-4)

    # (c) birth and death on a finite toy, exactly: multisets of ridge types up to M_max under a truncated
    #     Poisson count, a ridge prior π_r, a proposal q_r ≠ π_r. A_birth = LR · λ/(M+1) · P_D(M+1)/P_B(M) · π_r/q_r.
    R = 4
    M_max = 3
    lam = Fr(3, 2)
    pr = random_weights_exact(R)
    qr = random_weights_exact(R)
    states = [tuple(sorted(c)) for M in range(M_max + 1) for c in itertools.combinations_with_replacement(range(R), M)]
    L = {s: Fr(int(rng.integers(1, 20)), int(rng.integers(1, 7))) for s in states}

    def PB(M):
        return Fr(1, 2) if M < M_max else Fr(0)

    def PD(M):
        return Fr(1, 2) if M > 0 else Fr(0)

    def target(s):
        M = len(s)
        mult = math.prod(math.factorial(s.count(r)) for r in set(s))
        return L[s] * lam ** M / math.factorial(M) * Fr(math.factorial(M), mult) * math.prod((pr[r] for r in s), start=Fr(1))

    def birth_ratio(s, r, extra=Fr(1)):
        t = tuple(sorted(s + (r,)))
        M = len(s)
        return L[t] / L[s] * lam / (M + 1) * PD(M + 1) / PB(M) * pr[r] / qr[r] * extra

    def kernel(extra=Fr(1)):
        K = {s: {} for s in states}
        for s in states:
            M = len(s)
            if PB(M) > 0:
                for r in range(R):
                    t = tuple(sorted(s + (r,)))
                    K[s][t] = K[s].get(t, Fr(0)) + PB(M) * qr[r] * min(Fr(1), birth_ratio(s, r, extra))
            if PD(M) > 0:
                for r in set(s):
                    t = list(s)
                    t.remove(r)
                    t = tuple(t)
                    K[s][t] = K[s].get(t, Fr(0)) + PD(M) * Fr(s.count(r), M) * min(Fr(1), 1 / birth_ratio(t, r, extra))
            K[s][s] = K[s].get(s, Fr(0)) + 1 - sum(K[s].values())
        return K

    def balance_gap(K):
        return max(abs(target(s) * K[s].get(t, Fr(0)) - target(t) * K[t].get(s, Fr(0))) for s in states for t in states)

    K = kernel()
    gap = balance_gap(K)
    stat = max(abs(sum(target(s) * K[s].get(t, Fr(0)) for s in states) - target(t)) for t in states)
    report("AS-3c birth/death, multisets, exactly: detailed balance and stationarity", f"gaps {gap}, {stat}",
           gap == 0 and stat == 0)
    Kd = {s: {} for s in states}  # the doubled count: an extra 1/(M+1) in the birth ratio
    for s in states:
        M = len(s)
        if PB(M) > 0:
            for r in range(R):
                t = tuple(sorted(s + (r,)))
                Kd[s][t] = Kd[s].get(t, Fr(0)) + PB(M) * qr[r] * min(Fr(1), birth_ratio(s, r, Fr(1, M + 1)))
        if PD(M) > 0:
            for r in set(s):
                t = list(s)
                t.remove(r)
                t = tuple(t)
                Kd[s][t] = Kd[s].get(t, Fr(0)) + PD(M) * Fr(s.count(r), M) * min(Fr(1), 1 / birth_ratio(t, r, Fr(1, M)))
        Kd[s][s] = Kd[s].get(s, Fr(0)) + 1 - sum(Kd[s].values())
    gap_dbl = balance_gap(Kd)
    report("AS-3c an extra 1/(M+1) in the birth ratio (the count read twice) breaks balance", f"gap {float(gap_dbl):.3e}",
           gap_dbl > 0)

    # the same ratio on ordered tuples, birth inserting at a uniform position, death deleting a uniform position
    tuples = [c for M in range(M_max + 1) for c in itertools.product(range(R), repeat=M)]
    Lt = {c: L[tuple(sorted(c))] for c in tuples}

    def target_t(c):
        return Lt[c] * lam ** len(c) / math.factorial(len(c)) * math.prod((pr[r] for r in c), start=Fr(1))

    Kt = {c: {} for c in tuples}
    for c in tuples:
        M = len(c)
        if PB(M) > 0:
            for pos in range(M + 1):
                for r in range(R):
                    t = c[:pos] + (r,) + c[pos:]
                    ratio = Lt[t] / Lt[c] * lam / (M + 1) * PD(M + 1) / PB(M) * pr[r] / qr[r]
                    Kt[c][t] = Kt[c].get(t, Fr(0)) + PB(M) * Fr(1, M + 1) * qr[r] * min(Fr(1), ratio)
        if PD(M) > 0:
            for pos in range(M):
                t = c[:pos] + c[pos + 1:]
                r = c[pos]
                ratio = Lt[c] / Lt[t] * lam / M * PD(M) / PB(M - 1) * pr[r] / qr[r]  # the birth t → c
                Kt[c][t] = Kt[c].get(t, Fr(0)) + PD(M) * Fr(1, M) * min(Fr(1), 1 / ratio)
        Kt[c][c] = Kt[c].get(c, Fr(0)) + 1 - sum(Kt[c].values())
    gap_t = max(abs(target_t(s) * Kt[s].get(t, Fr(0)) - target_t(t) * Kt[t].get(s, Fr(0))) for s in tuples for t in tuples)
    report("AS-3c the same ratio on ordered tuples (insert and delete at uniform positions), exactly", f"gap {gap_t}",
           gap_t == 0)

    # (d) dimension matching on a continuous toy: M ∈ {0, 1}, a ridge (θ on S¹, t₀ on [0, 1]); births draw from the
    #     prior, the Jacobian is 1; a change move with the power spherical on S¹. Against the exact posterior.
    lam_c, kap0, th0 = 1.5, 2.0, 0.7

    def h(th, t):
        return math.exp(kap0 * math.cos(th - th0)) * (1 + t)

    eh = special.i0(kap0) * 1.5
    p1_exact = lam_c * eh / (1 + lam_c * eh)
    cos_exact = math.cos(th0) * special.i1(kap0) / special.i0(kap0)
    M, th, t0 = 0, 0.0, 0.5
    n1 = 0
    cs = 0.0
    N = 400000
    u = rng.random((N, 4))
    zs = rng.beta(20 + 0.5, 0.5, size=N)
    sg = np.where(rng.random(N) < 0.5, -1.0, 1.0)
    for it in range(N):
        mv = u[it, 0]
        if mv < 1 / 3:  # birth (from M = 0) or a refused birth (at M = 1, the cap)
            if M == 0:
                thn, tn = 2 * math.pi * u[it, 1], u[it, 2]
                if u[it, 3] < min(1.0, h(thn, tn) * lam_c * (1 / 3) / (1 / 3)):
                    M, th, t0 = 1, thn, tn
        elif mv < 2 / 3:  # death
            if M == 1 and u[it, 3] < min(1.0, 1 / (h(th, t0) * lam_c)):
                M = 0
        elif M == 1:  # change: θ′ from the power spherical about θ, κ = 20, t₀′ from its prior
            thn = th + sg[it] * math.acos(2 * zs[it] - 1)
            tn = u[it, 2]
            if u[it, 3] < min(1.0, h(thn, tn) / h(th, t0)):
                th, t0 = thn, tn
        if M == 1:
            n1 += 1
            cs += math.cos(th)
    p1 = n1 / N
    report("AS-3d continuous toy: P(M = 1) against the exact posterior", f"{p1:.4f} against {p1_exact:.4f}",
           abs(p1 - p1_exact) < 0.01)
    report("AS-3d continuous toy: E[cos θ | M = 1] against the exact", f"{cs / max(n1, 1):.4f} against {cos_exact:.4f}",
           abs(cs / max(n1, 1) - cos_exact) < 0.02)


# ================================================================ AS-4: rotate-and-refit's non-invariance

def as4():
    print("\n== AS-4: rotate-and-refit's non-invariance")
    # (a) the frame: a hinge along a·u₁ + b·u₂ has full rank on a designed n × n grid, so no sum of M < n
    #     separable terms equals it; a tensor-product hinge sum of M terms has grid rank at most M
    full = True
    for _ in range(200):
        a = Fr(int(rng.integers(1, 9)) * (1 if rng.random() < 0.5 else -1), int(rng.integers(1, 5)))
        b = Fr(int(rng.integers(1, 9)) * (1 if rng.random() < 0.5 else -1), int(rng.integers(1, 5)))
        t = Fr(int(rng.integers(-5, 6)), int(rng.integers(1, 4)))
        n = int(rng.integers(2, 8))
        xs = [Fr(i) / a for i in range(n)]
        ys = [(t - j + Fr(1, 2)) / b for j in range(n)]
        G = [[max(Fr(0), a * x + b * y - t) for y in ys] for x in xs]
        full &= det(G) != 0 and all(G[i][i] == Fr(1, 2) for i in range(n))
    report("AS-4a (a u₁ + b u₂ − t)₊ on the designed grid: lower triangular, diagonal ½, nonzero determinant",
           "200 random a, b ≠ 0, t, n ≤ 7, exactly", full)
    worst_rank = 0
    for _ in range(100):
        M = int(rng.integers(1, 8))
        ug = rng.uniform(0, 1, 40)
        vg = rng.uniform(0, 1, 40)
        Z = np.zeros((40, 40))
        for _ in range(M):
            kind = rng.integers(0, 3)
            fu = (np.maximum(0, rng.choice([-1, 1]) * (ug - rng.uniform())) if kind != 1 else np.ones(40))
            fv = (np.maximum(0, rng.choice([-1, 1]) * (vg - rng.uniform())) if kind != 0 else np.ones(40))
            Z += rng.standard_normal() * np.outer(fu, fv)
        worst_rank = max(worst_rank, np.linalg.matrix_rank(Z, tol=1e-10) - M)
    report("AS-4a an M-term tensor-product hinge sum has grid rank ≤ M", f"largest rank − M = {worst_rank}", worst_rank <= 0)
    ranks = {}
    ug = np.linspace(-1, 1, 30)  # the kink crosses the grid at every angle
    for deg in (0, 10, 45, 90, 135, 180, 270):
        c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
        Z = np.maximum(0, c * ug[:, None] - s * ug[None, :] - 0.1)
        ranks[deg] = int(np.linalg.matrix_rank(Z, tol=1e-9))
    ok = all((ranks[d] <= 1) == (d % 90 == 0) for d in ranks)
    report("AS-4a hinge in u₁ after rotating by φ: rank 1 at multiples of 90°, large otherwise (30 × 30 grid)",
           str(ranks), ok)

    # (b) the scale: C_z = LᵀC_xL; orthogonal L carries eigenpairs; a diagonal L moves the leading direction;
    #     the inactive subspace ker C maps by L⁻¹ under every invertible L, so a rank-k C's range is equivariant
    C = mat([[1, 0], [0, Fr(1, 2)]])
    L = mat([[1, 0], [0, 2]])
    Cz = mm(mm(T(L), C), L)
    report("AS-4b C = diag(1, ½), L = diag(1, 2): LᵀCL = diag(1, 2), the leading direction moves from x₁ to x₂",
           f"LᵀCL = {[[str(x) for x in r] for r in Cz]}", Cz == mat([[1, 0], [0, 2]]))
    worst_orth = worst_ker = worst_rng = 0.0
    moved = 0
    for _ in range(200):
        p = int(rng.integers(2, 7))
        k = int(rng.integers(1, p))
        X = rng.standard_normal((p, p))
        Cf = X @ X.T
        Q, _ = np.linalg.qr(rng.standard_normal((p, p)))
        ev1, V1 = np.linalg.eigh(Cf)
        ev2, V2 = np.linalg.eigh(Q.T @ Cf @ Q)
        worst_orth = max(worst_orth, np.abs(ev1 - ev2).max(),
                         1 - principal_cosines(Q.T @ V1[:, -k:], V2[:, -k:]).min())
        Lg = rng.standard_normal((p, p)) + 3 * np.eye(p)
        B = rng.standard_normal((p, k))
        Ck = B @ B.T  # rank k: an exact k-index function's C
        Czk = Lg.T @ Ck @ Lg
        ns = np.linalg.svd(Ck)[2][k:].T  # ker C
        worst_ker = max(worst_ker, np.abs(Czk @ np.linalg.solve(Lg, ns)).max() / np.abs(Czk).max())
        top = np.linalg.eigh(Czk)[1][:, -k:]
        worst_rng = max(worst_rng, 1 - principal_cosines(top, np.linalg.qr(Lg.T @ B)[0]).min())
        # a non-conformal L moves the leading direction of some full-rank C (built from L's SVD)
        U, sv, Vt = np.linalg.svd(Lg)
        i, j = int(np.argmax(sv)), int(np.argmin(sv))
        c = np.full(p, 0.1)
        c[j] = 1.0
        c[i] = 0.9  # c_j > c_i, but σ_i² c_i > σ_j² c_j whenever σ_i² / σ_j² > 1/0.9
        if sv[i] ** 2 * 0.9 > sv[j] ** 2 * 1.0:
            Cx = U @ np.diag(c) @ U.T
            lead_x = np.linalg.eigh(Cx)[1][:, -1]
            lead_z = np.linalg.eigh(Lg.T @ Cx @ Lg)[1][:, -1]
            moved += int(1 - principal_cosines((Lg.T @ lead_x)[:, None] / np.linalg.norm(Lg.T @ lead_x),
                                               lead_z[:, None]).min() > 1e-3)
    report("AS-4b orthogonal L: eigenvalues kept, leading subspaces carried by Lᵀ", f"worst {worst_orth:.2e}",
           worst_orth < 1e-10)
    report("AS-4b every invertible L: ker(LᵀCL) = L⁻¹ ker C", f"worst {worst_ker:.2e}", worst_ker < 1e-10)
    report("AS-4b rank-k C: the leading k of LᵀCL is Lᵀ range C, whatever L", f"worst 1 − cos {worst_rng:.2e}",
           worst_rng < 1e-9)
    report("AS-4b a non-conformal L moves the leading direction of a C built on its singular vectors",
           f"moved in {moved} of the cases built", moved > 150)

    # (c) the g-form's hat matrix is invariant to B ↦ B T (T invertible), exactly; the raw ridge's marginal is not
    hat_ok = True
    for _ in range(60):
        n, k = int(rng.integers(3, 7)), int(rng.integers(1, 3))
        Bm = mat([[int(rng.integers(-4, 5)) for _ in range(k)] for _ in range(n)])
        g = inv(mm(T(Bm), Bm))
        if g is None:
            continue
        Tm = mat([[int(rng.integers(-3, 4)) for _ in range(k)] for _ in range(k)])
        if det(Tm) == 0:
            continue
        BT = mm(Bm, Tm)
        H1 = mm(mm(Bm, g), T(Bm))
        H2 = mm(mm(BT, inv(mm(T(BT), BT))), T(BT))
        hat_ok &= H1 == H2
    report("AS-4c the g-form's hat B(BᵀB)⁻¹Bᵀ under B ↦ BT, T invertible (so R² and the marginal): unchanged, exactly",
           "60 designs", hat_ok)
    y = rng.standard_normal(30)
    b = rng.uniform(0, 1, 30)

    def raw_ridge_logml(col, tau=1.0, s2=0.3):
        S = s2 * np.eye(30) + tau * np.outer(col, col)
        return -0.5 * (np.linalg.slogdet(S)[1] + y @ np.linalg.solve(S, y))

    gap = abs(raw_ridge_logml(b) - raw_ridge_logml(4 * b))
    report("AS-4c raw ridge: the log marginal moves when a column is scaled by 4", f"moves by {gap:.4f}", gap > 1e-3)

    # (d) the training-range rule: the evaluated hinge does not read the declared box; κ does, by (W′/W)^q
    ok_h, ok_k = True, True
    for _ in range(300):
        q = int(rng.integers(1, 4))
        nf = int(rng.integers(1, 4))
        val_W = val_W2 = Fr(1)
        kap_W = kap_W2 = Fr(1)
        ratio = Fr(1)
        for _ in range(nf):
            lo = Fr(int(rng.integers(-9, 0)), 2)
            hi = lo + Fr(int(rng.integers(2, 12)), 2)
            t = lo + (hi - lo) * Fr(int(rng.integers(1, 7)), 8)
            x = lo + (hi - lo) * Fr(int(rng.integers(-4, 13)), 8)  # inside and beyond the training range
            W = (hi - lo) * Fr(int(rng.integers(1, 9)), 1)  # a declared box's width, ≥ the training width
            W2 = W * Fr(int(rng.integers(3, 13)), 2)  # another declared box, wider by 1.5 to 6
            up = rng.random() < 0.5
            raw = max(Fr(0), (x - t) if up else (t - x))
            peak = (hi - t) if up else (t - lo)
            # in the box's units: hinge (raw/W)^q, peak (peak/W)^q; normalized: their ratio
            val_W *= (raw / W) ** q / (peak / W) ** q
            val_W2 *= (raw / W2) ** q / (peak / W2) ** q
            kap_W *= (peak / W) ** q
            kap_W2 *= (peak / W2) ** q
            ratio *= (W2 / W) ** q
        ok_h &= val_W == val_W2
        ok_k &= kap_W / kap_W2 == ratio and ratio != 1
    report("AS-4d training-range rule: the evaluated hinge product is the same under any declared box, exactly",
           "300 products of 1 to 3 factors, degrees 1 to 3", ok_h)
    report("AS-4d κ moves with the box by the product of (W′/W)^q", "300 products", ok_k)


if __name__ == "__main__":
    as1()
    as2()
    as3()
    as4()
    bad = [lab for lab, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(bad)} of {len(RESULTS)} checks pass")
    for lab in bad:
        print("  FAILED:", lab)
