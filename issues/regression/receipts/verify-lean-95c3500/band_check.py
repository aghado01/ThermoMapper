# Independent rebuild of the tridiagonal family, in 50-digit arithmetic (mpmath), with no dense inverse in doubles.
# Path graph on n nodes, unit edges except the two at node j = n//2, of weight W.
# S1 = L (weighted Laplacian), S2 = 2 diag(L) - L (signless Laplacian). c' = (1, 1): A(c') = 2 diag(L), E(c') = I.
# c = (Mm, 1): A(c) = (Mm - 1) L + 2 diag(L), ratios (Mm, 1), so M/m = Mm.
# T = ||E(c)^-1||_1 / ((M/m) ||E(c')^-1||_1) = ||E(c)^-1||_1 / Mm.
import math
from mpmath import mp, mpf, sqrt as msqrt
mp.dps = 50

def family(n, W):
    j = n // 2
    w = [mpf(1)] * (n - 1)
    w[j - 1] = w[j] = mpf(W)
    deg = [(w[k - 1] if k > 0 else 0) + (w[k] if k < n - 1 else 0) for k in range(n)]
    return j, w, deg

def column_sum(n, W, Mm, col):
    """Sum of |x| for E x = e_col, E the unit-diagonal scaling of A(c); Thomas algorithm, E is SPD tridiagonal."""
    j, w, deg = family(n, W)
    Mm = mpf(Mm)
    off = [-(Mm - 1) * w[k] / ((Mm + 1) * msqrt(deg[k] * deg[k + 1])) for k in range(n - 1)]
    d = [mpf(1)] * n
    rhs = [mpf(0)] * n
    rhs[col] = mpf(1)
    for k in range(1, n):
        f = off[k - 1] / d[k - 1]
        d[k] -= f * off[k - 1]
        rhs[k] -= f * rhs[k - 1]
    x = [mpf(0)] * n
    x[n - 1] = rhs[n - 1] / d[n - 1]
    for k in range(n - 2, -1, -1):
        x[k] = (rhs[k] - off[k] * x[k + 1]) / d[k]
    return sum(abs(v) for v in x)

def lead(n, W, Mm):
    """The rank-one part alone: E's near-null vector is u = sqrt(deg), with eigenvalue delta = 2/(Mm + 1)."""
    j, w, deg = family(n, W)
    u = [msqrt(v) for v in deg]
    return (mpf(Mm) + 1) / (2 * mpf(Mm)) * u[j] * sum(u) / sum(v * v for v in u)

def best_W(n, Mm):
    grid = [math.exp(t * math.log(4 * n) / 39) for t in range(40)]
    return max(grid, key=lambda W: float(column_sum(n, W, Mm, n // 2)))

print("M/m = 1e12:   n      W    T at the hub column   over all columns   rank-one part   T/sqrt(n)   (1+sqrt n)/2")
for n in (16, 64, 256, 1024):
    W = best_W(n, 1e12)
    hub = column_sum(n, W, 1e12, n // 2) / mpf(1e12)
    allc = max(column_sum(n, W, 1e12, c) for c in range(n)) / mpf(1e12) if n <= 256 else None
    print(f"{n:17d} {W:7.1f} {float(hub):12.4f} {'' if allc is None else format(float(allc), '18.4f'):>18} {float(lead(n, W, 1e12)):14.4f} {float(hub)/math.sqrt(n):10.3f} {(1+math.sqrt(n))/2:12.2f}")
print("n = 256, by M/m:")
for Mm in (1e2, 1e4, 1e6):
    W = best_W(256, Mm)
    print(f"   M/m = {Mm:.0e}: W = {W:.1f}, T at the hub column = {float(column_sum(256, W, Mm, 128) / mpf(Mm)):.4f}")
