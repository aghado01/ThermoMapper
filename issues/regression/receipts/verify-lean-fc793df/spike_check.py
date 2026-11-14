# Independent rebuild of the spike family: S1 = I - vv^T, S2 = vv^T, v = (3/5, 4/(5r), ..., 4/(5r)) on r^2+1 coordinates.
# Exact rationals for A and A^{-1}; floats only for the square roots of the diagonal.
from fractions import Fraction as F
import math

def run(r):
    n = r * r + 1
    v = [F(3, 5)] + [F(4, 5 * r)] * (n - 1)
    assert sum(x * x for x in v) == 1
    # A(c) at c = (2, 1): 2(I - vv^T) + vv^T = 2I - vv^T.  A(c') at c' = (1, 1) is I.
    A = [[(2 if i == j else 0) - v[i] * v[j] for j in range(n)] for i in range(n)]
    # Claimed inverse (Sherman-Morrison): (I + vv^T)/2. Check A * Ainv = I exactly on two rows (cost), all rows for small n.
    Ainv = [[((1 if i == j else 0) + v[i] * v[j]) / 2 for j in range(n)] for i in range(n)]
    rows = range(n) if n <= 50 else (0, 1, n - 1)
    for i in rows:
        for j in range(n):
            assert sum(A[i][k] * Ainv[k][j] for k in range(n)) == (1 if i == j else 0)
    d = [math.sqrt(A[i][i]) for i in range(n)]
    # E = D^{-1/2} A D^{-1/2};  E^{-1} = D^{1/2} A^{-1} D^{1/2}
    inv1 = max(sum(abs(float(Ainv[i][j])) * d[i] * d[j] for i in range(n)) for j in (0, 1))
    e1 = max(sum(abs(float(A[i][j])) / (d[i] * d[j]) for i in range(n)) for j in (0, 1))
    return n, inv1, e1

print(" r     n   ||E(c)^-1||_1   6r/25   (M/m)(1+sqrt n)/2   kappa_1(E(c))   kappa_1 / ((M/m)*kappa_1(E(c')))")
for r in (1, 2, 5, 10, 20, 40):
    n, inv1, e1 = run(r)
    print(f"{r:2d} {n:5d} {inv1:14.4f} {6*r/25:8.3f} {2*(1+math.sqrt(n))/2:14.3f} {inv1*e1:16.3f} {inv1*e1/2:18.3f}")
