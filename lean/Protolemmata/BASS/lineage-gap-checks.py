"""Numeric checks of the statements posed in this folder's protolemmata, labelled by item.

Checks, not proofs: each line prints a residual or an estimate to compare by eye. Needs numpy.
"""
import itertools
import math

import numpy as np
from numpy.polynomial import legendre as leg

rng = np.random.default_rng(20261101)


def expect(f, mus):
    out = f
    for ax in reversed(range(len(mus))):
        out = np.tensordot(out, mus[ax], axes=([ax], [0]))
    return float(out)


def coord_mean(f, mus, v):
    m = np.tensordot(f, mus[v], axes=([v], [0]))
    return np.broadcast_to(np.expand_dims(m, v), f.shape).copy()


def hoeffding(f, mus):
    d = len(mus)
    parts = {}
    for k in range(d + 1):
        for u in itertools.combinations(range(d), k):
            g = f.copy()
            for v in range(d):
                if v not in u:
                    g = coord_mean(g, mus, v)
            for k2 in range(k):
                for w in itertools.combinations(u, k2):
                    g = g - parts[w]
            parts[u] = g
    return parts


# A sum of product terms on a finite product space
sizes = (3, 4, 3)
mus = [rng.dirichlet(np.ones(s)) for s in sizes]
terms = []
f = np.zeros(sizes)
for _ in range(5):
    kappa = rng.normal()
    fac = [rng.normal(size=s) for s in sizes]
    terms.append((kappa, fac))
    f += kappa * np.einsum("i,j,k->ijk", *fac)
parts = hoeffding(f, mus)
V = {u: expect(p * p, mus) for u, p in parts.items() if u}
Vtot = sum(V.values())

print("SG-1 components sum to f:", np.abs(sum(parts.values()) - f).max())
print("SG-1 Var f = sum V_u:", abs(Vtot - (expect(f * f, mus) - expect(f, mus) ** 2)))
worst = max(
    np.abs(sum(coord_mean(p, mus, v) - p for v in range(3)) + len(u) * p).max() for u, p in parts.items()
)
print("SG-1 L f_u = -|u| f_u, worst:", worst)
energy = sum(expect(f * (f - coord_mean(f, mus, v)), mus) for v in range(3))
print("SG-1 sum of totals = sum |u| V_u:", abs(energy - sum(len(u) * x for u, x in V.items())))
print("SG-1 Efron-Stein slack (>= 0):", energy - Vtot)

# The order spectrum from one polynomial per term pair, with no subset enumeration
poly = np.zeros(4)
for (ka, fa), (kb, fb) in itertools.product(terms, terms):
    p = np.array([ka * kb])
    for v in range(3):
        mean_product = float(np.sum(fa[v] * mus[v]) * np.sum(fb[v] * mus[v]))
        covariance = float(np.sum(fa[v] * fb[v] * mus[v])) - mean_product
        p = np.convolve(p, np.array([mean_product, covariance]))
    poly[: len(p)] += p
by_order = np.zeros(4)
for u, x in V.items():
    by_order[len(u)] += x
print("SG-1 V_k from pair polynomials vs enumeration:", np.abs(poly[1:] - by_order[1:]).max())

# Under dependence: Gibbs averages stop commuting, and the gap sets how far totals undercount
joint = rng.dirichlet(np.ones(9)).reshape(3, 3) + 2.0 * np.diag(rng.uniform(1, 2, 3))
joint /= joint.sum()
g2 = rng.normal(size=(3, 3))


def gibbs(axis):
    """The conditional-mean operator of a coordinate, as a matrix on vectorized functions."""
    out = np.zeros((9, 9))
    for a, b in itertools.product(range(3), range(3)):
        row = a * 3 + b
        if axis == 0:
            w = joint[:, b] / joint[:, b].sum()
            for c in range(3):
                out[row, c * 3 + b] = w[c]
        else:
            w = joint[a, :] / joint[a, :].sum()
            for c in range(3):
                out[row, a * 3 + c] = w[c]
    return out


E0, E1 = gibbs(0), gibbs(1)
pi = joint.ravel()
D = np.diag(np.sqrt(pi))
Di = np.diag(1 / np.sqrt(pi))
sym = D @ (E0 + E1 - 2 * np.eye(9)) @ Di
gap = np.sort(-np.linalg.eigvalsh((sym + sym.T) / 2))[1]
vec = g2.ravel()
var = pi @ vec**2 - (pi @ vec) ** 2
totals = sum(pi @ (vec * (vec - E @ vec)) for E in (E0, E1))
print("SG-1 dependent: commutator norm", np.abs(E0 @ E1 - E1 @ E0).max(), " gap", gap)
print("SG-1 dependent: Var", var, "<= totals/gap", totals / gap, " (totals", totals, ")")

# SG-2: sqrt(V_u) is 1-Lipschitz in L2(mu)
g = f + rng.normal(size=sizes)
pg = hoeffding(g, mus)
excess = max(
    abs(np.sqrt(expect(parts[u] ** 2, mus)) - np.sqrt(expect(pg[u] ** 2, mus))) - np.sqrt(expect((f - g) ** 2, mus))
    for u in parts
)
print("SG-2 Lipschitz excess (<= 0):", excess)

# SG-2: a triangle wave of amplitude h with n teeth keeps its L2 size and grows its gradient
xs = (np.arange(200000) + 0.5) / 200000
for teeth in (1, 10, 100):
    tri = 0.01 * (2 * np.abs((xs * teeth) % 1 - 0.5))
    slope = np.gradient(tri, xs)
    print(f"SG-2 teeth {teeth:4d}: L2 {np.sqrt(np.mean(tri**2)):.5f}  gradient L2 {np.sqrt(np.mean(slope**2)):.3f}")

# SG-3: a joint mixture of products
K = 4
weights = rng.dirichlet(np.ones(K))
comp = [[rng.dirichlet(np.ones(s)) for s in sizes] for _ in range(K)]
mixture = sum(w * np.einsum("i,j,k->ijk", *c) for w, c in zip(weights, comp))
(ka, fa), (kb, fb) = terms[0], terms[1]
brute = float(np.sum(mixture * ka * kb * np.einsum("i,j,k->ijk", *fa) * np.einsum("i,j,k->ijk", *fb)))
closed = sum(
    w * ka * kb * np.prod([np.sum(fa[v] * fb[v] * c[v]) for v in range(3)]) for w, c in zip(weights, comp)
)
marg = [mixture.sum(axis=tuple(a for a in range(3) if a != v)) for v in range(3)]
print("SG-3 pair integral vs the joint:", abs(brute - closed))
print("SG-3 the joint is dependent, by:", np.abs(mixture - np.einsum("i,j,k->ijk", *marg)).max())

# SG-4: shifted Legendre
nodes, wts = leg.leggauss(64)
x = (nodes + 1) / 2
w = wts / 2


def psi(n, z):
    c = np.zeros(n + 1)
    c[n] = 1
    return np.sqrt(2 * n + 1) * leg.legval(2 * z - 1, c)


def dpsi(n, z):
    c = np.zeros(n + 1)
    c[n] = 1
    return np.sqrt(2 * n + 1) * 2 * leg.legval(2 * z - 1, leg.legder(c))


beta = rng.normal(size=7)
f1 = sum(b * psi(n, x) for n, b in enumerate(beta))
df1 = sum(b * dpsi(n, x) for n, b in enumerate(beta))
weighted = np.sum(w * x * (1 - x) * df1**2)
print("SG-4 E[x(1-x)f'^2] vs sum n(n+1)a^2:", abs(weighted - sum(n * (n + 1) * b**2 for n, b in enumerate(beta))))
print("SG-4 Var f", np.sum(w * f1**2) - np.sum(w * f1) ** 2, "<= half of it", weighted / 2)
for alpha in [(5, 1), (3, 3)]:
    q, d = sum(a > 0 for a in alpha), sum(alpha)
    print("SG-4", alpha, "KHAOS's 1+q(d+q-2) =", 1 + q * (d + q - 2), " sum a(a+1) =", sum(a * (a + 1) for a in alpha))

# SG-4(b): the sum over a tensor basis, under the product uniform measure on [0, 1]^d.
# f = sum_alpha a_alpha prod_v psi_{alpha_v}(x_v). Per variable,
# E[x_v(1-x_v)(d_v f)^2] = sum_alpha a_alpha^2 alpha_v(alpha_v+1); summed over v, sum_alpha a^2 sum_v
# alpha_v(alpha_v+1). total(v) = E[(f - E_v f)^2], with E_v integrating out x_v, is
# sum_{alpha_v >= 1} a_alpha^2 <= half the per-variable sum, equal when every alpha_v <= 1.
def tensor_check(d, support, a, nq=12):
    nodes1, wts1 = leg.leggauss(nq)
    x1, w1 = (nodes1 + 1) / 2, wts1 / 2
    grids = np.meshgrid(*([x1] * d), indexing="ij")
    W = np.ones([nq] * d)
    for v in range(d):
        shape = [1] * d
        shape[v] = nq
        W = W * w1.reshape(shape)
    f = sum(c * np.prod([psi(al[v], grids[v]) for v in range(d)], axis=0) for al, c in zip(support, a))
    rows = []
    for v in range(d):
        df = sum(c * dpsi(al[v], grids[v]) * np.prod([psi(al[u], grids[u]) for u in range(d) if u != v], axis=0)
                 for al, c in zip(support, a))
        weighted = np.sum(W * grids[v] * (1 - grids[v]) * df**2)
        Evf = np.sum(f * w1.reshape([nq if u == v else 1 for u in range(d)]), axis=v, keepdims=True)
        total = np.sum(W * (f - Evf) ** 2)
        rows.append((weighted, sum(c**2 * al[v] * (al[v] + 1) for al, c in zip(support, a)),
                     total, sum(c**2 for al, c in zip(support, a) if al[v] >= 1)))
    return rows


worst = [0.0, 0.0, 0.0]
for trial in range(300):
    d = int(rng.integers(1, 4))
    deg = int(rng.integers(1, 4))
    allidx = list(itertools.product(range(deg + 1), repeat=d))
    k = int(rng.integers(1, len(allidx) + 1))
    support = [allidx[i] for i in rng.choice(len(allidx), size=k, replace=False)]
    a = rng.normal(size=k)
    rows = tensor_check(d, support, a)
    total_sum = sum(r[0] for r in rows)
    worst[0] = max(worst[0], abs(total_sum - sum(c**2 * sum(al_v * (al_v + 1) for al_v in al)
                                                   for al, c in zip(support, a))))
    for weighted, coef, total, tcoef in rows:
        worst[1] = max(worst[1], abs(weighted - coef), abs(total - tcoef))
        worst[2] = max(worst[2], total - weighted / 2)
print(f"SG-4(b) tensor, 300 seeded draws (d <= 3, degree <= 3 per variable, random sparse supports): "
      f"sum over v against sum a^2 sum a_v(a_v+1), worst {worst[0]:.1e}; per variable, both identities, "
      f"worst {worst[1]:.1e}; total(v) - half the weighted form, largest {worst[2]:.1e} (rounding: the bound holds, "
      f"with equality at degree 1)")
assert worst[0] < 1e-9 and worst[1] < 1e-9 and worst[2] < 1e-12
sharp = tensor_check(3, [(1, 0, 1), (0, 1, 1), (1, 1, 0), (0, 0, 0)], rng.normal(size=4))
print("SG-4(b) degree <= 1 in every variable: total(v) against half the weighted form:",
      [f"{t:.6f} = {w / 2:.6f}" for w, _, t, _ in sharp])

# SG-5: E[C]'s top direction against the mean projector's
a_dir = np.array([1.0, 0.0])
b_dir = np.array([np.cos(1.2), np.sin(1.2)])
draws = [(1.0, a_dir), (1.0, a_dir), (100.0, b_dir)]
EC = sum(lam * np.outer(v, v) for lam, v in draws) / 3
EP = sum(np.outer(v, v) for _, v in draws) / 3
top = lambda A: np.linalg.eigh(A)[1][:, -1]
print("SG-5 degrees between the two:", np.degrees(np.arccos(min(1.0, abs(top(EC) @ top(EP))))))

# SC-1: KHAOS's printed birth and death ratios against the Metropolis-Hastings form, on one term
q_max, p_dim, s_q, s_d, M, P_B, P_D = 3, 5, 1.0, 1.0, 4, 0.3, 0.3
chi = np.array([1, 1, 0, 0, 0])
q_order, degree = int(chi.sum()), 4
partitions = 3.0  # binom(degree - 1, q_order - 1)
base = rng.dirichlet(np.ones(p_dim))
eta = {q0: np.minimum(0.95, q0 * base) for q0 in range(1, q_max + 1)}
bern = {q0: np.prod(e**chi * (1 - e) ** (1 - chi)) for q0, e in eta.items()}
c_q = sum(q0**-s_q for q0 in range(1, q_max + 1))
c_d = sum(float(d) ** -s_d for d in range(q_order, p_dim + 1))
marginal = sum(q0**-s_q * bern[q0] for q0 in range(1, q_max + 1)) / c_q
printed_birth = lambda m, q0: P_D * partitions * bern[q0] / (P_B * (m + 1) * q0**s_q * c_q * degree**s_d * c_d)
printed_death = lambda m: P_B * degree**s_d * c_d / (P_D * m * partitions) * marginal
forward = P_B * marginal * degree**-s_d / c_d / partitions
right_birth = lambda m: (P_D / (m + 1)) / forward
right_death = lambda m: forward / (P_D / m)
for q0 in (1, 2):
    print(f"SC-1 printed A_B(M, q0={q0}) * A_D(M+1):", printed_birth(M, q0) * printed_death(M + 1))
print("SC-1 printed product's closed form at q0=1:", bern[1] * marginal / ((M + 1) ** 2 * c_q))
print("SC-1 Metropolis-Hastings A_B(M) * A_D(M+1):", right_birth(M) * right_death(M + 1))

# SC-2: a reversible kernel's gap against any set's flow
n_states = 30
target = rng.dirichlet(np.ones(n_states))
prop = rng.dirichlet(np.ones(n_states), size=n_states)
prop = (prop + prop.T) / 2
prop /= prop.sum(axis=1).max()
kernel = np.zeros((n_states, n_states))
for s, t in itertools.product(range(n_states), range(n_states)):
    if s != t:
        kernel[s, t] = prop[s, t] * min(1.0, target[t] / target[s])
kernel[np.diag_indices(n_states)] = 1 - kernel.sum(axis=1)
root = np.diag(np.sqrt(target))
rooti = np.diag(1 / np.sqrt(target))
sym = root @ kernel @ rooti
true_gap = 1 - np.sort(np.linalg.eigvalsh((sym + sym.T) / 2))[-2]
worst_slack = np.inf
for _ in range(20000):
    A = rng.random(n_states) < rng.uniform(0.05, 0.95)
    if A.all() or not A.any():
        continue
    mass = target[A].sum()
    flow = (target[A, None] * kernel[A][:, ~A]).sum()
    worst_slack = min(worst_slack, flow / (mass * (1 - mass)) - true_gap)
print("SC-2 detailed balance:", np.abs(target[:, None] * kernel - (target[:, None] * kernel).T).max())
print("SC-2 gap", true_gap, " smallest flow/variance - gap over random sets (>= 0):", worst_slack)

# SC-3: KHAOS's posterior precision
Psi = rng.normal(size=(40, 6))
gm = rng.uniform(0.2, 1.0, 6)
g02 = 2.3
Gram = Psi.T @ Psi
G = (g02 * np.outer(gm, gm) + 1) / (g02 * np.outer(gm, gm))
S0 = g02 * np.diag(gm) @ np.linalg.inv(Gram) @ np.diag(gm)
print("SC-3 G o Gram vs Gram + S0^-1:", np.abs(G * Gram - (Gram + np.linalg.inv(S0))).max())

# GE-1, GE-2: the cross-fit difference estimator
mean_fn = lambda z: 3.0 + np.sin(3 * z)


def ridge_fit(xtr, ytr):
    """Well behaved: its mean residual is about the same on and off its training points."""
    B = lambda z: np.vander(z, 9, increasing=True)
    A = B(xtr)
    coef = np.linalg.solve(A.T @ A + np.eye(9), A.T @ ytr)
    return lambda z: B(z) @ coef


def spike_fit(xtr, ytr, h=0.004):
    """Built to expose the bias: it returns to zero between its training points."""
    ker = lambda a: np.exp(-((a[:, None] - xtr[None, :]) ** 2) / (2 * h**2))
    coef = np.linalg.solve(ker(xtr) + 1e-3 * np.eye(len(xtr)), ytr)
    return lambda z: ker(z) @ coef


reps = 4000
for name, fit in (("ridge", ridge_fit), ("spike", spike_fit)):
    for N, n in [(60, 30), (60, 12), (600, 12), (6000, 12)]:
        pool = rng.uniform(size=N)
        goal = mean_fn(pool).mean()
        errs, formula = [], []
        for _ in range(reps):
            idx = rng.permutation(N)[:n]
            y = mean_fn(pool[idx]) + 0.3 * rng.normal(size=n)
            one, two, y1, y2 = idx[: n // 2], idx[n // 2 :], y[: n // 2], y[n // 2 :]
            est, pred = 0.0, 0.0
            for tr, ytr, te, yte in ((one, y1, two, y2), (two, y2, one, y1)):
                fitted = fit(pool[tr], ytr)(pool)
                est += fitted.mean() + (yte - fitted[te]).mean()
                r = mean_fn(pool) - fitted
                out = np.ones(N, bool)
                out[tr] = False
                pred += (len(tr) / N) * (r[out].mean() - r[tr].mean())
            errs.append(est / 2 - goal)
            formula.append(pred / 2)
        e = np.array(errs)
        print(
            f"GE-1 {name} N={N:5d} n={n:3d}: bias {e.mean():+.4f} +- {e.std() / np.sqrt(reps):.4f}"
            f"  formula {np.mean(formula):+.4f}"
        )

gx, gw = leg.leggauss(4000)
gx, gw = (gx + 1) / 2, gw / 2
truth = np.sum(gw * mean_fn(gx))
for name, fit in (("ridge", ridge_fit), ("spike", spike_fit)):
    for n in (12, 30):
        errs = []
        for _ in range(reps):
            zs = rng.uniform(size=n)
            ys = mean_fn(zs) + 0.3 * rng.normal(size=n)
            one, two = np.arange(n // 2), np.arange(n // 2, n)
            est = 0.0
            for tr, te in ((one, two), (two, one)):
                model = fit(zs[tr], ys[tr])
                est += np.sum(gw * model(gx)) + (ys[te] - model(zs[te])).mean()
            errs.append(est / 2 - truth)
        e = np.array(errs)
        print(f"GE-2 {name} iid with the integral, n={n:3d}: bias {e.mean():+.4f} +- {e.std() / np.sqrt(reps):.4f}")

# GE-3: telescoping levels, the tuned weight, and known sampling probabilities
level0 = lambda z: 2.5 + 0.8 * np.sin(3 * z)
level1 = lambda z: mean_fn(z) + 0.2 * np.cos(5 * z)
goal = np.sum(gw * mean_fn(gx))
pool_n, runs, sub = 400, 40, 10
tele = []
for _ in range(20000):
    pool = rng.uniform(size=pool_n)
    zs = rng.uniform(size=runs)
    y1 = level1(zs) + 0.2 * rng.normal(size=runs)
    y2 = mean_fn(zs[:sub]) + 0.3 * rng.normal(size=sub)
    tele.append(level0(pool).mean() + (y1 - level0(zs)).mean() + (y2 - y1[:sub]).mean())
tele = np.array(tele)
print(f"GE-3 telescoping bias {tele.mean() - goal:+.4f} +- {tele.std() / np.sqrt(len(tele)):.4f}")

pool_n, runs = 200, 40
zs = rng.uniform(size=400000)
ys = mean_fn(zs) + 0.3 * rng.normal(size=zs.size)
cov = np.cov(ys, level0(zs))[0, 1]
var_f = level0(zs).var()
best = cov / ((1 + runs / pool_n) * var_f)
variance = lambda lam: lam**2 * var_f / pool_n + (ys - lam * level0(zs)).var() / runs
grid = np.linspace(0, 2, 2001)
print("GE-3 tuned weight", best, " grid minimizer", grid[np.argmin([variance(v) for v in grid])])
print("GE-3 variance at weight 0, 1, tuned:", variance(0.0), variance(1.0), variance(best))

pool = rng.uniform(size=300)
goal = mean_fn(pool).mean()
inclusion = np.clip(0.03 + 0.25 * np.abs(mean_fn(pool) - level0(pool)), 0.02, 1.0)
ipw = []
for _ in range(40000):
    taken = rng.random(300) < inclusion
    y = mean_fn(pool[taken]) + 0.3 * rng.normal(size=taken.sum())
    ipw.append(level0(pool).mean() + np.sum((y - level0(pool[taken])) / inclusion[taken]) / 300)
ipw = np.array(ipw)
print(f"GE-3 known probabilities bias {ipw.mean() - goal:+.4f} +- {ipw.std() / np.sqrt(len(ipw)):.4f}")


# SC-4: the mean-field Potts barrier at alpha_c, per site, along the symmetric path (Woodard et al., eq. 10-11).
def potts_f(x, q, alpha):
    y = (1 - x) / (q - 1)
    return alpha / 2 * (x * x + (q - 1) * y * y) - x * np.log(x) - (q - 1) * y * np.log(y)


for q in (3, 4):
    alpha_c = 2 * (q - 1) * np.log(q - 1) / (q - 2)
    xs = np.linspace(1 / q + 1e-6, (q - 1) / q - 1e-6, 200001)
    c = potts_f(1 / q, q, alpha_c) - potts_f(xs, q, alpha_c).min()
    n = 10.0
    for _ in range(200):
        n = max(10.0, 3 * np.log(n) / c)
    print(f"SC-4 Potts q={q}: barrier per site {c:.5f}, across N=64..2048 e^{c * (2048 - 64):.2f} against"
          f" N^3's e^{3 * np.log(32):.2f}, e^(cN) overtakes N^3 near N={n:,.0f}")

# SC-5 under the proper sigma^2 prior IG(c, c s^2), BASS's default (Sigma2Prior.ScaledToResponse(1e-3)),
# with the g-prior. Given sigma^2 the centred response is N(0, sigma^2 (I + g P_X)), so
# |I + g P| = (1 + g)^M and y'(I + gP)^(-1) y = S = TSS (1 + g(1 - R^2))/(1 + g). Integrating sigma^2
# against IG(c, b) gives log p = -(M/2) log(1 + g) - (nu/2 + c) log(b + S/2) + const(nu, c, b).
def ig_log_marginal_direct(Xc, yc, g, c, b, nu):
    """The log marginal from the design itself: the Gaussian's determinant and quadratic form, and
    the sigma^2 integral by quadrature in log sigma^2."""
    P = Xc @ np.linalg.pinv(Xc)
    A = np.eye(len(yc)) + g * P
    _, logdet = np.linalg.slogdet(A)
    S = yc @ np.linalg.solve(A, yc)
    t = np.linspace(-12, 12, 20001)
    s2 = np.exp(t)
    logint = (-nu / 2 * t - S / (2 * s2) + c * np.log(b) - math.lgamma(c) - (c + 1) * t - b / s2) + t
    m = logint.max()
    return -logdet / 2 + m + np.log(np.sum(np.exp(logint - m)) * (t[1] - t[0])), S


def ig_profile(g, nu, c, kappa, R2):
    return -(nu / 2 + c) * np.log(kappa + (1 + g * (1 - R2)) / (1 + g))


worst_closed, worst_const, least_slack = 0.0, 0.0, np.inf
for trial in range(400):
    N = int(rng.integers(8, 30))
    M = int(rng.integers(1, 5))
    k = int(rng.integers(1, 4))
    X = rng.normal(size=(N, M + k))
    y = X[:, :M] @ rng.normal(size=M) * rng.uniform(0, 2) + rng.normal(size=N) * rng.uniform(0.1, 3)
    yc = y - y.mean()
    Xc = X - X.mean(axis=0)
    tss = yc @ yc
    s2 = tss / N
    c = [1e-3, 0.5, 2.0][trial % 3]
    b = c * s2
    g = float(np.exp(rng.uniform(-2, 7)))
    nu = N - 1
    L1, S1 = ig_log_marginal_direct(Xc[:, :M], yc, g, c, b, nu)
    L2, S2 = ig_log_marginal_direct(Xc, yc, g, c, b, nu)
    R1 = 1 - (Xc[:, :M] @ np.linalg.lstsq(Xc[:, :M], yc, rcond=None)[0] - yc) @ (
        Xc[:, :M] @ np.linalg.lstsq(Xc[:, :M], yc, rcond=None)[0] - yc) / tss
    R2_ = 1 - (Xc @ np.linalg.lstsq(Xc, yc, rcond=None)[0] - yc) @ (Xc @ np.linalg.lstsq(Xc, yc, rcond=None)[0] - yc) / tss
    kappa = 2 * b / tss
    # closed form: -(M/2) log(1+g) + profile + a constant the same for both structures
    const1 = L1 - (-(M / 2) * np.log1p(g) + ig_profile(g, nu, c, kappa, R1))
    const2 = L2 - (-((M + k) / 2) * np.log1p(g) + ig_profile(g, nu, c, kappa, R2_))
    worst_const = max(worst_const, abs(const1 - const2))
    worst_closed = max(worst_closed, abs(S1 - tss * (1 + g * (1 - R1)) / (1 + g)))
    least_slack = min(least_slack, L2 - (L1 - k * np.log1p(g) / 2))
assert worst_const < 1e-6 and worst_closed < 1e-6 * 1e3 and least_slack >= -1e-9
print(f"SC-5 under IG(c, c s^2), 400 seeded nested designs (N 8-29, M 1-4, k 1-3, c in {{1e-3, 1/2, 2}}, "
      f"g from e^-2 to e^7): the log marginal from the design and sigma^2 quadrature is -(M/2) log(1+g) plus "
      f"the profile up to one constant, worst gap {worst_const:.1e}; the superset step's least slack "
      f"{least_slack:.2e} >= 0")
print(f"SC-5 at c = kappa = 0 the IG profile is the Jeffreys profile: "
      f"{abs(ig_profile(3.0, 20, 0, 0, 0.4) - 10 * (np.log(4.0) - np.log(1 + 3.0 * 0.6))):.1e}")
