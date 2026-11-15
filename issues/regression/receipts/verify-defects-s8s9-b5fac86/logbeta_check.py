"""LogBeta's two forms (slice 9, b8418c4), transcribed into doubles and read against mpmath, apart from the chip's
400-digit decimal reference (coordinating session).

  s < 10:   ln Gamma(s) - LogGammaRatio(l, s)                      (GammaRatio.cs as fixed in 117230e)
  s >= 10:  -ln(l)/2 + ln(2 pi)/2 + mu(s) + mu(l) - mu(s + l) + (s - 1/2) ln(w) + l log1p(-w),  w = s/(s + l)
against the three-term form it replaces, lgamma(a) + lgamma(b) - lgamma(a + b). Python's lgamma stands in for the
engine's LogGamma, so a difference of an ulp or two in the small-shape form is the stand-in's and not the algorithm's.
The error is read over max(1, |ln B|), as the engine's docstring states its bound.
"""
import math
import random
import sys

import mpmath as mp

mp.mp.dps = 420   # ln Gamma(1e300) is near 7e302: a difference of order 1 needs more than 302 digits
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


def remainder_difference(x, y, s):
    p, q, h, q_power, total = 1.0 / x, 1.0 / y, 1.0, 1.0, STIRLING[0]
    for k in range(1, len(STIRLING)):
        q_power *= q
        h = p * h + q_power
        q_power *= q
        h = p * h + q_power
        total += STIRLING[k] * h
    return -(s / y) / x * total


def log_gamma_ratio(x, s):
    if s == 0.0:
        return 0.0
    steps = 0.0
    while x < 16.0:
        ratio = s / x
        steps -= (math.log(s) - math.log(x)) if ratio == math.inf else log1p_goldberg(ratio)
        x += 1.0
    y = x + s
    return steps + (x - 0.5) * log1p_goldberg(s / x) + s * (math.log(y) - 1.0) + remainder_difference(x, y, s)


def stirling_correction(x):
    f = 1.0 / (x * x)
    return (1.0 / 12 - f * (1.0 / 360 - f * (1.0 / 1260 - f * (1.0 / 1680 - f * (1.0 / 1188 - f * (691.0 / 360360 - f * (
        1.0 / 156 - f * 3617.0 / 122400))))))) / x


def log_beta_new(a, b):
    s, l = min(a, b), max(a, b)
    if s < 10.0:
        return math.lgamma(s) - log_gamma_ratio(l, s)
    total = s + l
    w = 1.0 / (1.0 + l / s) if total == math.inf else s / total
    return (-0.5 * math.log(l) + 0.5 * math.log(2.0 * math.pi) + stirling_correction(s) + stirling_correction(l)
            - stirling_correction(total) + (s - 0.5) * math.log(w) + l * log1p_goldberg(-w))


def log_beta_old(a, b):
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def truth(a, b):
    a, b = mp.mpf(a), mp.mpf(b)
    return mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)


def main():
    rng = random.Random(20261115)
    regimes = {
        "both under 10": lambda: (10 ** rng.uniform(-6, 1), 10 ** rng.uniform(-6, 1)),
        "one under 10, one to 1e15": lambda: (10 ** rng.uniform(-6, 1), 10 ** rng.uniform(1, 15)),
        "one under 10, one to 1e300": lambda: (10 ** rng.uniform(-12, 1), 10 ** rng.uniform(15, 300)),
        "both 10 to 1e6": lambda: (10 ** rng.uniform(1, 6), 10 ** rng.uniform(1, 6)),
        "balanced, 1e6 to 1e15": lambda: (lambda x: (x, x * 10 ** rng.uniform(-1, 1)))(10 ** rng.uniform(6, 15)),
        "balanced, 1e15 to 1e290": lambda: (lambda x: (x, x * 10 ** rng.uniform(-1, 1)))(10 ** rng.uniform(15, 290)),
        "both at least 10, far apart": lambda: (10 ** rng.uniform(1, 3), 10 ** rng.uniform(10, 290)),
    }
    print(f"{'regime':<30} {'new, worst':>11} {'old, worst':>11}   error over max(1, |ln B|), 300 points each")
    ok = True
    for name, draw in regimes.items():
        new_worst = old_worst = 0.0
        for _ in range(300):
            a, b = draw()
            want = truth(a, b)
            scale = max(mp.mpf(1), abs(want))
            new_worst = max(new_worst, float(abs(mp.mpf(log_beta_new(a, b)) - want) / scale))
            try:
                old = log_beta_old(a, b)
            except OverflowError:
                old = math.inf
            old_worst = max(old_worst, float(abs(mp.mpf(old) - want) / scale)) if math.isfinite(old) else math.inf
        print(f"{name:<30} {new_worst:>11.2e} {old_worst:>11.2e}")
        ok &= new_worst < 5e-15
    print("ln B(x, 1) = -ln x:")
    for x in (1e5, 1e16, 1e100):
        print(f"   x = {x:g}: new {log_beta_new(x, 1.0)!r}, old {log_beta_old(x, 1.0)!r}, -ln x {-math.log(x)!r}")
    print(f"ln B(1/2, 1/2) = ln pi = {math.log(math.pi)!r}: new {log_beta_new(0.5, 0.5)!r}, old {log_beta_old(0.5, 0.5)!r}")
    print("new form within 5e-15 everywhere" if ok else "THE NEW FORM EXCEEDS 5e-15 SOMEWHERE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
