"""Searches for penalized-admission-lemmas.md, PA-7: how much of the proved q·((1 + √n)/2) the whole chain can use.

A family with m·c′ₖ ≤ cₖ ≤ M·c′ₖ gives m·A(c′) ≼ A(c) ≼ M·A(c′) (PA-6(a) in Loewner form), and every such pair is a
two-term family: S₁ = (A − mA′)/(M − m), S₂ = (MA′ − A)/(M − m) with c′ = (1, 1), c = (M, m). So a search over pairs
searches over every family. With m = 1: A′ = LLᵀ, A = A′ + (M − 1)·A′^{1/2} Q A′^{1/2}, 0 ≼ Q ≼ I.

Readings, not checks: Powell from seeded random starts, so a reading is a lower bound on the supremum.
- T = ‖E(c)⁻¹‖₁ / ((M/m)·‖E(c′)⁻¹‖₁), the inverse part, proved at most (1 + √n)/2;
- R = κ₁(E(c)) / ((M/m)·‖E(c′)⁻¹‖₁), dense (q = n), proved at most q·(1 + √n)/2.

Run: python -I -X utf8 penalized-admission-searches.py [seed]   (the recorded readings are seed 3, the default)
     python -I -X utf8 penalized-admission-searches.py banded   (banded two-term families, seed 5: does T grow with n
     at fixed q? The dense spike family of finding 6 shows it does when q = n.)
"""
import math
import sys
import numpy as np
from scipy import optimize, linalg

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 3)


def unit_scale(A):
    d = 1 / np.sqrt(np.diag(A))
    return A * np.outer(d, d)


def norm1(A):
    return np.abs(A).sum(axis=0).max()


def pair(p, n):
    i = 0
    L = np.tril(p[i:i + n * n].reshape(n, n)); i += n * n
    V = np.linalg.qr(p[i:i + n * n].reshape(n, n))[0]; i += n * n
    s = 1 / (1 + np.exp(-p[i:i + n])); i += n
    M = 1 + math.exp(min(p[i], 30))
    Ap = L @ L.T + 1e-9 * np.eye(n)
    R = linalg.sqrtm(Ap).real
    Q = V @ np.diag(s) @ V.T
    A = Ap + (M - 1) * R @ Q @ R
    return A, Ap, M


def readings(p, n):
    try:
        A, Ap, M = pair(p, n)
        if np.linalg.cond(Ap) > 1e13 or np.linalg.cond(A) > 1e13:
            return 0.0, 0.0
        E, Ep = unit_scale(A), unit_scale(Ap)
        Ei = np.linalg.inv(E)
        den = M * norm1(np.linalg.inv(Ep))
        return norm1(Ei) / den, norm1(E) * norm1(Ei) / den
    except Exception:
        return 0.0, 0.0


if "banded" not in sys.argv:
    for n, trials in [(3, 16), (4, 16), (6, 10), (8, 8), (12, 5), (16, 4)]:
        bt = br = 0
        for which in (0, 1):
            for t in range(trials):
                x0 = np.r_[rng.standard_normal(2 * n * n), rng.uniform(-4, 4, n), rng.uniform(0, 8)]
                res = optimize.minimize(lambda p: -readings(p, n)[which], x0, method="Powell",
                                        options={"maxiter": 40000, "maxfev": 40000, "xtol": 1e-8, "ftol": 1e-11})
                val = -res.fun
                if which == 0:
                    bt = max(bt, val)
                else:
                    br = max(br, val)
        g = (1 + math.sqrt(n)) / 2
        print(f"n={n}: T {bt:.4f} (gamma {g:.4f}), R {br:.4f} (q(1+sqrt n)/2 = {n * g:.4f}, R/q = {br / n:.4f})", flush=True)


# ---------------------------------------------------------------- banded two-term families

def band_L(p, n, w):
    L = np.zeros((n, n))
    i = 0
    for k in range(w + 1):
        L += np.diag(p[i:i + n - k], -k)
        i += n - k
    return L, i


def T_banded(p, n, w, K):
    S, i = [], 0
    for _ in range(K):
        L, used = band_L(p[i:], n, w)
        i += used
        S.append(L @ L.T)
    lc = p[i:i + K]
    c = np.exp(lc)
    cp = np.ones(K)
    A = sum(ck * Sk for ck, Sk in zip(c, S))
    Ap = sum(S)
    try:
        if np.linalg.cond(Ap) > 1e12 or np.linalg.cond(A) > 1e12:
            return 0.0
        E, Ep = unit_scale(A), unit_scale(Ap)
        rho = c / cp
        return norm1(np.linalg.inv(E)) / (rho.max() / rho.min() * norm1(np.linalg.inv(Ep)))
    except Exception:
        return 0.0


if "banded" in sys.argv:
    rng = np.random.default_rng(5)
    for w in (1, 2):
        for n in (8, 16, 32, 64):
            K = 2
            dim = K * sum(n - k for k in range(w + 1)) + K
            best = 0
            for t in range(4):
                x0 = np.r_[rng.standard_normal(dim - K), rng.uniform(-3, 3, K)]
                res = optimize.minimize(lambda p: -T_banded(p, n, w, K), x0, method="Powell",
                                        options={"maxiter": 30000, "maxfev": 30000, "xtol": 1e-8, "ftol": 1e-11})
                best = max(best, -res.fun)
            print(f"half-width {w} (q={2 * w + 1}), n={n}: T {best:.4f}  ((1+sqrt n)/2 = {(1 + math.sqrt(n)) / 2:.3f})", flush=True)
