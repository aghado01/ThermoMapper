"""c_M = E[det(Theta' Theta)] under BPPR's sparse direction law, against the closed form p!/((p - M)! p^M)
(coordinating session, slice 4, 3ef7dc4). Written from the law's definition, apart from the engine:
a uniform on 1..A, J uniform among the subsets of size a, theta_J uniform on its subsphere, M columns independent.
By Cauchy-Binet the mean depends on E[theta theta'] = I/p alone, which this law has; the Monte Carlo below reads it.
"""
import math
import sys

import numpy as np


def draw(rng, p, cap, m, count):
    theta = np.zeros((count, p, m))
    for k in range(m):
        a = rng.integers(1, cap + 1, size=count)
        for c in range(count):
            j = rng.choice(p, size=a[c], replace=False)
            w = rng.standard_normal(a[c])
            theta[c, j, k] = w / np.linalg.norm(w)
    return theta


def main():
    rng = np.random.default_rng(20261115)
    count = 60_000
    ok = True
    print(f"{'p':>3} {'A':>2} {'M':>2} {'closed form':>12} {'Monte Carlo':>12} {'z':>6}   E[theta theta'] against I/p, worst entry")
    for p, cap, m in [(6, 3, 2), (6, 3, 3), (6, 3, 6), (10, 3, 5), (4, 2, 4), (6, 1, 3), (8, 3, 9)]:
        theta = draw(rng, p, cap, m, count)
        gram = np.einsum('cik,cil->ckl', theta, theta)
        det = np.linalg.det(gram)
        closed = math.factorial(p) / (math.factorial(p - m) * p ** m) if m <= p else 0.0
        se = det.std(ddof=1) / math.sqrt(count)
        # Past M = p the determinant is exactly 0, and what is read is rounding: its size is the check there, not a z-score.
        z = (det.mean() - closed) / se if m <= p and se > 0 else 0.0
        if m > p:
            ok &= np.abs(det).max() < 1e-9
        second = np.einsum('ci,cj->ij', theta[:, :, 0], theta[:, :, 0]) / count
        worst = np.abs(second - np.eye(p) / p).max()
        ok &= abs(z) < 4.0 and worst < 0.01
        print(f"{p:>3} {cap:>2} {m:>2} {closed:>12.6f} {det.mean():>12.6f} {z:>6.2f}   {worst:.4f}")
    print("the closed form holds for the sparse law" if ok else "A CASE DISAGREES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
