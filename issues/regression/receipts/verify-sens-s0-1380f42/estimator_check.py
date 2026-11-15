"""The sampled-sensitivity chip's estimator comparison, rebuilt apart from the engine (coordinating session, 1380f42).

Sobol's g-function with a = (0, 1, 4.5, 9, 99, 99) under the uniform law on [0, 1]^6, whose main variances are
V_i = 1 / (3 (1 + a_i)^2) and V = prod(1 + V_i) - 1 in closed form. For each input, over R replicates at N rows:
  Saltelli's 2010 main variance, centred: cov(f(B), f(A_B^i) - f(A));        N (p + 2) values
  the two-difference form:               mean((f(B) - f(B_A^i)) (f(A_B^i) - f(A))) / 2;   N (2p + 2) values
  Jansen's total:                        mean((f(A) - f(A_B^i))^2) / 2.
Reported: each estimator's mean against the truth in units of its standard error (bias), its relative spread, and the
two-difference form's gain at equal cost, (spread ratio)^2 * (p + 2) / (2p + 2).
"""
import sys

import numpy as np

A_COEF = np.array([0.0, 1.0, 4.5, 9.0, 99.0, 99.0])


def g(x):
    return np.prod((np.abs(4.0 * x - 2.0) + A_COEF) / (1.0 + A_COEF), axis=-1)


def main(n=4096, reps=300, seed=20261115):
    rng = np.random.default_rng(seed)
    p = len(A_COEF)
    v_main = 1.0 / (3.0 * (1.0 + A_COEF) ** 2)
    v = np.prod(1.0 + v_main) - 1.0
    v_total = v_main * np.prod(1.0 + v_main) / (1.0 + v_main)   # V_i * prod_{j != i} (1 + V_j)
    salt = np.empty((reps, p))
    two = np.empty((reps, p))
    jan = np.empty((reps, p))
    for r in range(reps):
        a = rng.random((n, p))
        b = rng.random((n, p))
        ya, yb = g(a), g(b)
        for i in range(p):
            ab = a.copy()
            ab[:, i] = b[:, i]          # A with input i from B
            ba = b.copy()
            ba[:, i] = a[:, i]          # B with input i from A
            yab, yba = g(ab), g(ba)
            d = yab - ya
            salt[r, i] = np.sum((yb - yb.mean()) * (d - d.mean())) / (n - 1)
            two[r, i] = 0.5 * np.mean((yb - yba) * (yab - ya))
            jan[r, i] = 0.5 * np.mean((ya - yab) ** 2)
    print(f"g-function, N = {n}, {reps} replicates; V = {v:.6f}")
    print(f"{'a_i':>5} {'S_i':>9} {'Saltelli rel':>13} {'two-diff rel':>13} {'gain':>8} "
          f"{'bias z: S':>10} {'T':>6} {'J':>6} {'Jansen rel':>11}")
    ok = True
    for i in range(p):
        s_sd, t_sd, j_sd = salt[:, i].std(ddof=1), two[:, i].std(ddof=1), jan[:, i].std(ddof=1)
        zs = (salt[:, i].mean() - v_main[i]) / (s_sd / np.sqrt(reps))
        zt = (two[:, i].mean() - v_main[i]) / (t_sd / np.sqrt(reps))
        zj = (jan[:, i].mean() - v_total[i]) / (j_sd / np.sqrt(reps))
        gain = (s_sd / t_sd) ** 2 * (p + 2) / (2 * p + 2)
        ok &= max(abs(zs), abs(zt), abs(zj)) < 4.0
        print(f"{A_COEF[i]:>5g} {v_main[i] / v:>9.2e} {s_sd / v_main[i]:>12.1%} {t_sd / v_main[i]:>12.1%} {gain:>8.1f} "
              f"{zs:>10.2f} {zt:>6.2f} {zj:>6.2f} {j_sd / v_total[i]:>10.1%}")
    print("no bias beyond 4 standard errors" if ok else "A BIAS BEYOND 4 STANDARD ERRORS")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
