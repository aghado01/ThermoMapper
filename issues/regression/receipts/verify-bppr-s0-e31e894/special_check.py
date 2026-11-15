"""Checks of BPPR slice 0's primitives from outside the engine (coordinating session, e31e894).

1. The power spherical normalizer N = 2^(a+b) pi^b Gamma(a)/Gamma(a+b), a = (d-1)/2 + kappa, b = (d-1)/2,
   against quadrature of the surface integral |S^(d-2)| * int (1+t)^kappa (1-t^2)^((d-3)/2) dt.
2. LogGammaRatio's algorithm, transcribed into doubles from GammaRatio.cs, against mpmath's loggamma:
   the error by the size of s, and the normalizer it gives at kappa up to 1e5.
3. Cin(2 pi) and Friedman's three expectations, against quadrature.
"""
import math
import random
import sys

import mpmath as mp

mp.mp.dps = 60
STIRLING = [1.0 / 12.0, -1.0 / 360.0, 1.0 / 1260.0, -1.0 / 1680.0, 1.0 / 1188.0, -691.0 / 360360.0]


def log1p_goldberg(x):
    if x == math.inf:
        return x
    u = 1.0 + x
    if u == 1.0:
        return x
    try:
        corrected = math.log(u) * x / (u - 1.0)
    except OverflowError:
        corrected = math.inf
    return math.log(u) if corrected == math.inf else corrected


def stirling_remainder(y):
    inverse = 1.0 / y
    square = inverse * inverse
    total = 0.0
    for c in reversed(STIRLING):
        total = total * square + c
    return total * inverse


def remainder_difference(x, y, s):
    """R(y) - R(x) term by term: q^m - p^m = (q - p) * sum_j q^j p^(m-1-j), with q - p = -s/(x y)."""
    p, q = 1.0 / x, 1.0 / y
    total = 0.0
    for k, c in enumerate(STIRLING):
        m = 2 * k + 1
        total += c * sum(q ** j * p ** (m - 1 - j) for j in range(m))
    return -s / (x * y) * total


def log_gamma_ratio(x, s, differenced=False):
    """The engine's algorithm; returns the value and the largest magnitude it adds, the scale of its rounding."""
    if s == 0.0:
        return 0.0, 0.0
    steps, size = 0.0, 0.0
    while x < 16.0:
        ratio = s / x
        step = (math.log(s) - math.log(x)) if ratio == math.inf else log1p_goldberg(ratio)
        steps -= step
        size = max(size, abs(step), abs(steps))
        x += 1.0
    y = x + s
    first, second = (x - 0.5) * log1p_goldberg(s / x), s * (math.log(y) - 1.0)
    tail = remainder_difference(x, y, s) if differenced else stirling_remainder(y) - stirling_remainder(x)
    return steps + first + second + tail, max(size, abs(first), abs(second))


def truth_ratio(x, s):
    return mp.loggamma(mp.mpf(x) + mp.mpf(s)) - mp.loggamma(mp.mpf(x))


def main():
    ok = True

    print("1. normalizer formula against quadrature of the surface integral")
    worst = mp.mpf(0)
    for d in (2, 3, 4, 7, 16):
        for kappa in (0, 0.5, 3, 30):
            b = mp.mpf(d - 1) / 2
            a = b + kappa
            formula = (a + b) * mp.log(2) + b * mp.log(mp.pi) + mp.loggamma(a) - mp.loggamma(a + b)
            sphere = 2 * mp.pi ** b / mp.gamma(b)   # |S^(d-2)|
            if d == 2:
                integral = mp.quad(lambda th: (1 + mp.cos(th)) ** kappa, [0, mp.pi / 2, mp.pi]) * 2
            else:
                integral = sphere * mp.quad(lambda t: (1 + t) ** kappa * (1 - t * t) ** ((mp.mpf(d) - 3) / 2),
                                            [-1, -0.5, 0, 0.5, 1])
            worst = max(worst, abs(formula - mp.log(integral)))
    print(f"   20 cases, d in 2..16, kappa in 0..30: worst |ln N - ln quadrature| = {mp.nstr(worst, 3)}")
    ok &= worst < mp.mpf(10) ** -20

    print("2a. LogGammaRatio's algorithm in doubles against loggamma, by the size of s (400 x in [1e-6, 1e9] each)")
    print("    'scaled' is the error over the largest term the algorithm adds; 'rel' is over the result, away from its zeros")
    print(f"   {'s':>8} {'worst abs':>10} {'scaled':>9} {'rel, as built':>14} {'rel, differenced':>17}")
    for s in (1e-12, 1e-9, 1e-6, 1e-3, 0.5, 1.0, 7.5, 1e3, 1e6):
        rng = random.Random(20261114)
        wa = ws = wr = wd = 0.0
        for _ in range(400):
            x = 10.0 ** rng.uniform(-6, 9)
            want = truth_ratio(x, s)
            got, size = log_gamma_ratio(x, s)
            fixed, _ = log_gamma_ratio(x, s, differenced=True)
            ea = float(abs(mp.mpf(got) - want))
            wa, ws = max(wa, ea), max(ws, ea / size)
            if abs(want) > 0.05 * size:          # away from a zero of the ratio
                wr = max(wr, float(abs(mp.mpf(got) - want) / abs(want)))
                wd = max(wd, float(abs(mp.mpf(fixed) - want) / abs(want)))
        print(f"   {s:>8.0e} {wa:>10.2e} {ws:>9.2e} {wr:>14.2e} {wd:>17.2e}")
        ok &= wd < 2e-14 and (ws < 2e-14 or s < 1e-3)
    print("    finding: as built, the two remainders are subtracted as rounded values, so the error floors near 2e-18")
    print("    and the relative error grows as s falls below 1e-3; differenced term by term it holds at every s.")
    print("2b. ln N through that algorithm, d = 2..40: absolute error, as built and by two lgamma calls")
    for kappa in (0.0, 0.5, 30.0, 1e3, 1e5):
        built = plain = 0.0
        for d in range(2, 41):
            b = 0.5 * (d - 1)
            a = b + kappa
            want = (mp.mpf(a) + b) * mp.log(2) + b * mp.log(mp.pi) + mp.loggamma(a) - mp.loggamma(mp.mpf(a) + b)
            got = (a + b) * math.log(2.0) + b * math.log(math.pi) - log_gamma_ratio(a, b)[0]
            two = (a + b) * math.log(2.0) + b * math.log(math.pi) + math.lgamma(a) - math.lgamma(a + b)
            built, plain = max(built, float(abs(mp.mpf(got) - want))), max(plain, float(abs(mp.mpf(two) - want)))
        print(f"   kappa {kappa:>8g}: as built {built:.2e}, two lgamma {plain:.2e}")
        ok &= built < 1e-10
    print("3. Cin(2 pi) and Friedman's expectations")
    cin = mp.euler + mp.log(2 * mp.pi) - mp.ci(2 * mp.pi)
    print(f"   Cin(2 pi) = {mp.nstr(cin, 17)} (the report: 2.437653393057224)")
    ok &= abs(cin - mp.mpf("2.437653393057224")) < mp.mpf(10) ** -15
    mp.mp.dps = 25
    c = lambda x1, x2: mp.cos(mp.pi * x1 * x2)
    e1 = mp.quad(lambda x1, x2: x2 * c(x1, x2), [0, 1], [0, 1])
    e2 = mp.quad(lambda x1, x2: x2 ** 2 * c(x1, x2) ** 2, [0, 1], [0, 1])
    e3 = mp.quad(lambda x1, x2: x1 * x2 * c(x1, x2) ** 2, [0, 1], [0, 1])
    d1 = abs(e1 - 2 / mp.pi ** 2)
    d2 = abs(e2 - (mp.mpf(1) / 6 - 1 / (8 * mp.pi ** 2)))
    d3 = abs(e3 - (mp.mpf(1) / 8 - cin / (8 * mp.pi ** 2)))
    print(f"   E[x2 c] - 2/pi^2 = {mp.nstr(d1, 2)}; E[x2^2 c^2] - (1/6 - 1/(8 pi^2)) = {mp.nstr(d2, 2)}; "
          f"E[x1 x2 c^2] - (1/8 - Cin/(8 pi^2)) = {mp.nstr(d3, 2)}")
    ok &= max(d1, d2, d3) < mp.mpf(10) ** -15

    print("all consistent" if ok else "A CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
