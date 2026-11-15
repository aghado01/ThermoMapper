"""Rebuild of inverse_part_unbounded_tridiagonal's instance, apart from the chip's checks.

Written from the Lean definitions at 6f6bb27 (hubEdge, pathDeg, pathFamily, hubRho), not from
penalized-admission-checks.py. A path on N + 2 nodes; edge (j, j + 1) has weight 1 for j < N and
N for j = N. S1 = L, S2 = |L|, c' = (1, 1), c = (rho, 1), rho = 2(N + 1)^2 + 1, m = 1, M = rho.

Reads, in 50-digit arithmetic by tridiagonal solves:
  ||E(c)^-1||_1                      against the Lean bound (rho + 1) sqrt(N + 1) / 20,
  T = ||E(c)^-1||_1 / (M ||E(c')^-1||_1)   against sqrt(n),
  kappa_1(E(c))                      against the proved ceiling 3 (1 + sqrt n)/2 (M/m) ||E(c')^-1||_1.
A float64 dense inverse cross-checks the solves, and the hypotheses are read on the matrices.
"""
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 50


def edges(N):
    return [mp.mpf(1) if j < N else mp.mpf(N) for j in range(N + 1)]


def degrees(N):
    e, n = edges(N), N + 2
    return [(e[i - 1] if i >= 1 else 0) + (e[i] if i + 1 < n else 0) for i in range(n)]


def comb(N, c0, c1):
    """Diagonal and off-diagonal of c0 L + c1 |L|."""
    e, d = edges(N), degrees(N)
    return [(c0 + c1) * x for x in d], [(c1 - c0) * x for x in e]


def unit_scale(diag, off):
    s = [mp.sqrt(x) for x in diag]
    return [mp.mpf(1)] * len(diag), [off[j] / (s[j] * s[j + 1]) for j in range(len(off))]


def inverse_norm1(diag, off):
    """Largest absolute column sum of the inverse of a symmetric tridiagonal matrix."""
    n = len(diag)
    piv, mult = [diag[0]], []
    for i in range(1, n):
        mult.append(off[i - 1] / piv[i - 1])
        piv.append(diag[i] - mult[-1] * off[i - 1])
    best, arg = mp.mpf(0), -1
    for j in range(n):
        y = [mp.mpf(0)] * n
        y[j] = mp.mpf(1)
        for i in range(j + 1, n):
            y[i] = -mult[i - 1] * y[i - 1]
        x = [mp.mpf(0)] * n
        x[n - 1] = y[n - 1] / piv[n - 1]
        for i in range(n - 2, -1, -1):
            x[i] = (y[i] - off[i] * x[i + 1]) / piv[i]
        total = mp.fsum(abs(v) for v in x)
        if total > best:
            best, arg = total, j
    return best, arg


def dense(diag, off):
    n = len(diag)
    a = np.zeros((n, n))
    for i in range(n):
        a[i, i] = float(diag[i])
    for j in range(n - 1):
        a[j, j + 1] = a[j + 1, j] = float(off[j])
    return a


def main(sizes):
    print(f"{'N':>5} {'n':>5} {'rho':>10} {'||E(c)^-1||_1':>16} {'Lean bound':>14} {'ratio':>7} "
          f"{'T':>8} {'T/sqrt n':>9} {'kappa_1/ceiling':>16} {'col':>5} {'f64 rel':>9}")
    ok = True
    for N in sizes:
        n = N + 2
        rho = 2 * mp.mpf(N + 1) ** 2 + 1
        # Hypotheses, read on the matrices.
        d1, o1 = comb(N, mp.mpf(1), mp.mpf(0))      # S1 = L
        d2, o2 = comb(N, mp.mpf(0), mp.mpf(1))      # S2 = |L|
        if n <= 258:
            psd = min(np.linalg.eigvalsh(dense(d1, o1)).min(), np.linalg.eigvalsh(dense(d2, o2)).min())
            ok &= psd > -1e-9 * float(N + 1)
        dref, oref = comb(N, mp.mpf(1), mp.mpf(1))  # A(c') = 2 diag(deg)
        ok &= all(x == 0 for x in oref) and all(x > 0 for x in dref)
        ref_inv = mp.mpf(1)                         # E(c') = I
        # The reading.
        dc, oc = comb(N, rho, mp.mpf(1))
        de, oe = unit_scale(dc, oc)
        inv, col = inverse_norm1(de, oe)
        norm = max((abs(oe[j - 1]) if j >= 1 else 0) + 1 + (abs(oe[j]) if j + 1 < n else 0)
                   for j in range(n))
        bound = (rho + 1) * mp.sqrt(N + 1) / 20
        T = inv / (rho * ref_inv)
        ceiling = 3 * (1 + mp.sqrt(n)) / 2 * rho * ref_inv
        rel = mp.mpf(0)
        if n <= 1026:
            f = np.abs(np.linalg.inv(dense(de, oe))).sum(axis=0).max()
            rel = abs(mp.mpf(float(f)) - inv) / inv
        ok &= inv >= bound and norm * inv <= ceiling and rel < mp.mpf(10) ** -6
        print(f"{N:>5} {n:>5} {mp.nstr(rho, 8):>10} {mp.nstr(inv, 10):>16} {mp.nstr(bound, 8):>14} "
              f"{mp.nstr(inv / bound, 4):>7} {mp.nstr(T, 5):>8} {mp.nstr(T / mp.sqrt(n), 4):>9} "
              f"{mp.nstr(norm * inv / ceiling, 4):>16} {col:>5} {mp.nstr(rel, 2):>9}")
    print("all readings consistent" if ok else "A READING FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main([1, 4, 16, 64, 256, 1024]))
