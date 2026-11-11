Lengths 1/1 of S10's, 3 design(s) per benchmark, Release, 18 logical processors, server GC False; 60.4 minutes in all.

### The verdict

| arm | score | clears |
| --- | --- | --- |
| A1, g, τ held at 1/N, λ (10, 10) | NaN | no |
| A2, ridge, τ sampled, λ (10, 10) | NaN | no |
| A3, ridge, τ held at 1, λ (10, 10) | NaN | no |
| A4, g, Zellner–Siow, λ (1, 1) | NaN | no |
| A5, g, τ held at 1/N, λ (1, 1) | NaN | no |
| A6, ridge, τ sampled, λ (1, 1) | NaN | no |
| A7, ridge, τ held at 1, λ (1, 1) | NaN | no |

No arm clears: no profile is built, every option stays an option, and F-R1 closes as read.

### The gate (noise-free)

Fits whose predictive R̂ exceeds 1.05, each leaving its (benchmark, design) out for every arm:

- Friedman, k = 0, A0: R̂ 2.368 (sigma2)
- Friedman, k = 0, A1: R̂ 1.625 (f3)
- Friedman, k = 0, A2: R̂ 3.069 (sigma2)
- Friedman, k = 0, A3: R̂ 2.144 (sigma2)
- Friedman, k = 0, A4: R̂ 3.376 (sigma2)
- Friedman, k = 0, A5: R̂ 1.357 (f3)
- Friedman, k = 0, A6: R̂ 3.215 (sigma2)
- Friedman, k = 0, A7: R̂ 2.166 (sigma2)
- Friedman, k = 1, A0: R̂ 2.047 (f3)
- Friedman, k = 1, A1: R̂ 1.983 (f3)
- Friedman, k = 1, A2: R̂ 2.823 (sigma2)
- Friedman, k = 1, A3: R̂ 1.844 (sigma2)
- Friedman, k = 1, A4: R̂ 3.338 (sigma2)
- Friedman, k = 1, A5: R̂ 1.688 (f2)
- Friedman, k = 1, A6: R̂ 2.696 (sigma2)
- Friedman, k = 1, A7: R̂ 2.104 (sigma2)
- Friedman, k = 2, A0: R̂ 3.350 (sigma2)
- Friedman, k = 2, A1: R̂ 2.782 (f2)
- Friedman, k = 2, A2: R̂ 3.438 (sigma2)
- Friedman, k = 2, A3: R̂ 2.660 (sigma2)
- Friedman, k = 2, A4: R̂ 2.726 (f1)
- Friedman, k = 2, A5: R̂ 1.724 (f2)
- Friedman, k = 2, A6: R̂ 2.054 (f1)
- Friedman, k = 2, A7: R̂ 1.949 (sigma2)
- Ishigami, k = 0, A0: R̂ 2.161 (sigma2)
- Ishigami, k = 0, A1: R̂ 1.823 (f1)
- Ishigami, k = 0, A2: R̂ 3.156 (sigma2)
- Ishigami, k = 0, A3: R̂ 2.207 (sigma2)
- Ishigami, k = 0, A4: R̂ 2.610 (sigma2)
- Ishigami, k = 0, A5: R̂ 1.755 (f3)
- Ishigami, k = 0, A6: R̂ 2.762 (sigma2)
- Ishigami, k = 0, A7: R̂ 2.124 (sigma2)
- Ishigami, k = 1, A0: R̂ 2.374 (sigma2)
- Ishigami, k = 1, A1: R̂ 1.906 (sigma2)
- Ishigami, k = 1, A2: R̂ 3.195 (sigma2)
- Ishigami, k = 1, A3: R̂ 1.554 (sigma2)
- Ishigami, k = 1, A4: R̂ 2.804 (sigma2)
- Ishigami, k = 1, A5: R̂ 1.823 (sigma2)
- Ishigami, k = 1, A6: R̂ 2.141 (sigma2)
- Ishigami, k = 1, A7: R̂ 1.721 (sigma2)
- Ishigami, k = 2, A0: R̂ 2.017 (f3)
- Ishigami, k = 2, A1: R̂ 2.138 (sigma2)
- Ishigami, k = 2, A2: R̂ 2.455 (sigma2)
- Ishigami, k = 2, A3: R̂ 1.557 (f1)
- Ishigami, k = 2, A4: R̂ 2.156 (f1)
- Ishigami, k = 2, A5: R̂ 1.902 (sigma2)
- Ishigami, k = 2, A6: R̂ 2.656 (sigma2)
- Ishigami, k = 2, A7: R̂ 1.651 (f1)
- borehole, k = 0, A0: R̂ 2.056 (f3)
- borehole, k = 0, A1: R̂ 2.681 (f2)
- borehole, k = 0, A2: R̂ 2.094 (f3)
- borehole, k = 0, A3: R̂ 1.719 (f3)
- borehole, k = 0, A4: R̂ 2.036 (f2)
- borehole, k = 0, A5: R̂ 1.742 (f3)
- borehole, k = 0, A6: R̂ 2.493 (f1)
- borehole, k = 0, A7: R̂ 1.607 (f2)
- borehole, k = 1, A0: R̂ 2.662 (f2)
- borehole, k = 1, A1: R̂ 2.072 (f1)
- borehole, k = 1, A2: R̂ 2.391 (sigma2)
- borehole, k = 1, A3: R̂ 2.208 (f3)
- borehole, k = 1, A4: R̂ 2.534 (f2)
- borehole, k = 1, A5: R̂ 1.493 (f3)
- borehole, k = 1, A6: R̂ 2.899 (sigma2)
- borehole, k = 1, A7: R̂ 1.709 (sigma2)
- borehole, k = 2, A0: R̂ 2.230 (sigma2)
- borehole, k = 2, A1: R̂ 1.341 (f1)
- borehole, k = 2, A2: R̂ 2.043 (f2)
- borehole, k = 2, A3: R̂ 2.007 (sigma2)
- borehole, k = 2, A4: R̂ 2.410 (f1)
- borehole, k = 2, A5: R̂ 1.725 (f2)
- borehole, k = 2, A6: R̂ 2.158 (sigma2)
- borehole, k = 2, A7: R̂ 2.145 (sigma2)
- functional, k = 0, A0: R̂ 1.913 (f2)
- functional, k = 0, A1: R̂ 1.904 (f2)
- functional, k = 0, A2: R̂ 1.878 (sigma2)
- functional, k = 0, A3: R̂ 2.123 (sigma2)
- functional, k = 0, A4: R̂ 1.674 (sigma2)
- functional, k = 0, A5: R̂ 1.863 (f2)
- functional, k = 0, A6: R̂ 1.786 (sigma2)
- functional, k = 0, A7: R̂ 2.068 (f2)
- functional, k = 1, A0: R̂ 2.703 (f2)
- functional, k = 1, A1: R̂ 2.295 (f2)
- functional, k = 1, A2: R̂ 2.505 (sigma2)
- functional, k = 1, A3: R̂ 2.441 (sigma2)
- functional, k = 1, A4: R̂ 3.022 (sigma2)
- functional, k = 1, A5: R̂ 2.431 (f2)
- functional, k = 1, A6: R̂ 2.477 (sigma2)
- functional, k = 1, A7: R̂ 1.815 (sigma2)
- functional, k = 2, A0: R̂ 3.165 (sigma2)
- functional, k = 2, A1: R̂ 1.841 (f1)
- functional, k = 2, A2: R̂ 2.222 (sigma2)
- functional, k = 2, A3: R̂ 1.599 (sigma2)
- functional, k = 2, A4: R̂ 2.728 (sigma2)
- functional, k = 2, A5: R̂ 2.177 (f1)
- functional, k = 2, A6: R̂ 2.963 (sigma2)
- functional, k = 2, A7: R̂ 2.013 (sigma2)

Scored: none. Unscored: Friedman, Ishigami, borehole, functional.

### Context for the verdict

Each arm's per-design ln(NRMSE/NRMSE_A0), and their sd:

| arm |
| --- |
| A1, g, τ held at 1/N, λ (10, 10) |
| A2, ridge, τ sampled, λ (10, 10) |
| A3, ridge, τ held at 1, λ (10, 10) |
| A4, g, Zellner–Siow, λ (1, 1) |
| A5, g, τ held at 1/N, λ (1, 1) |
| A6, ridge, τ sampled, λ (1, 1) |
| A7, ridge, τ held at 1, λ (1, 1) |

Each arm's score with one benchmark left out:

| arm |
| --- |
| A1, g, τ held at 1/N, λ (10, 10) |
| A2, ridge, τ sampled, λ (10, 10) |
| A3, ridge, τ held at 1, λ (10, 10) |
| A4, g, Zellner–Siow, λ (1, 1) |
| A5, g, τ held at 1/N, λ (1, 1) |
| A6, ridge, τ sampled, λ (1, 1) |
| A7, ridge, τ held at 1, λ (1, 1) |

### The noisy controls (read, not ruled)

| arm | score | loses to A0 (score > 1) |
| --- | --- | --- |
| A1, g, τ held at 1/N, λ (10, 10) | NaN | no |
| A2, ridge, τ sampled, λ (10, 10) | NaN | no |
| A3, ridge, τ held at 1, λ (10, 10) | NaN | no |
| A4, g, Zellner–Siow, λ (1, 1) | NaN | no |
| A5, g, τ held at 1/N, λ (1, 1) | NaN | no |
| A6, ridge, τ sampled, λ (1, 1) | NaN | no |
| A7, ridge, τ held at 1, λ (1, 1) | NaN | no |

### The gate (noisy controls)

Fits whose predictive R̂ exceeds 1.05, each leaving its (benchmark, design) out for every arm:

- Friedman, σ = 1, k = 0, A0: R̂ 1.681 (f1)
- Friedman, σ = 1, k = 0, A1: R̂ 1.482 (f1)
- Friedman, σ = 1, k = 0, A2: R̂ 1.921 (f3)
- Friedman, σ = 1, k = 0, A3: R̂ 1.556 (f3)
- Friedman, σ = 1, k = 0, A4: R̂ 1.589 (f1)
- Friedman, σ = 1, k = 0, A5: R̂ 1.572 (f2)
- Friedman, σ = 1, k = 0, A6: R̂ 1.583 (f3)
- Friedman, σ = 1, k = 0, A7: R̂ 1.517 (f1)
- Friedman, σ = 1, k = 1, A0: R̂ 1.798 (f3)
- Friedman, σ = 1, k = 1, A1: R̂ 1.930 (f3)
- Friedman, σ = 1, k = 1, A2: R̂ 1.473 (f3)
- Friedman, σ = 1, k = 1, A3: R̂ 1.660 (f3)
- Friedman, σ = 1, k = 1, A4: R̂ 1.546 (f3)
- Friedman, σ = 1, k = 1, A5: R̂ 1.796 (f3)
- Friedman, σ = 1, k = 1, A6: R̂ 1.458 (f3)
- Friedman, σ = 1, k = 1, A7: R̂ 1.695 (f3)
- Friedman, σ = 1, k = 2, A0: R̂ 1.758 (f2)
- Friedman, σ = 1, k = 2, A1: R̂ 1.992 (f2)
- Friedman, σ = 1, k = 2, A2: R̂ 1.231 (f2)
- Friedman, σ = 1, k = 2, A3: R̂ 2.134 (f2)
- Friedman, σ = 1, k = 2, A4: R̂ 2.308 (f2)
- Friedman, σ = 1, k = 2, A5: R̂ 2.141 (f2)
- Friedman, σ = 1, k = 2, A6: R̂ 1.726 (f2)
- Friedman, σ = 1, k = 2, A7: R̂ 1.437 (f2)
- Ishigami, σ = 0.37, k = 0, A0: R̂ 2.151 (f1)
- Ishigami, σ = 0.37, k = 0, A1: R̂ 1.534 (f1)
- Ishigami, σ = 0.37, k = 0, A2: R̂ 1.668 (f3)
- Ishigami, σ = 0.37, k = 0, A3: R̂ 1.853 (f3)
- Ishigami, σ = 0.37, k = 0, A4: R̂ 2.457 (f3)
- Ishigami, σ = 0.37, k = 0, A5: R̂ 1.768 (sigma2)
- Ishigami, σ = 0.37, k = 0, A6: R̂ 1.495 (f3)
- Ishigami, σ = 0.37, k = 0, A7: R̂ 1.408 (f3)
- Ishigami, σ = 0.37, k = 1, A0: R̂ 2.199 (f1)
- Ishigami, σ = 0.37, k = 1, A1: R̂ 2.562 (f3)
- Ishigami, σ = 0.37, k = 1, A2: R̂ 2.169 (f1)
- Ishigami, σ = 0.37, k = 1, A3: R̂ 1.657 (f3)
- Ishigami, σ = 0.37, k = 1, A4: R̂ 1.864 (f1)
- Ishigami, σ = 0.37, k = 1, A5: R̂ 1.584 (f1)
- Ishigami, σ = 0.37, k = 1, A6: R̂ 1.847 (sigma2)
- Ishigami, σ = 0.37, k = 1, A7: R̂ 1.290 (f3)
- Ishigami, σ = 0.37, k = 2, A0: R̂ 1.786 (f1)
- Ishigami, σ = 0.37, k = 2, A1: R̂ 1.704 (sigma2)
- Ishigami, σ = 0.37, k = 2, A2: R̂ 1.827 (f3)
- Ishigami, σ = 0.37, k = 2, A3: R̂ 1.802 (f3)
- Ishigami, σ = 0.37, k = 2, A4: R̂ 1.802 (f1)
- Ishigami, σ = 0.37, k = 2, A5: R̂ 1.951 (f1)
- Ishigami, σ = 0.37, k = 2, A6: R̂ 2.412 (f1)
- Ishigami, σ = 0.37, k = 2, A7: R̂ 1.693 (f1)

Scored: none. Unscored: Friedman, σ = 1, Ishigami, σ = 0.37.

### Per fit

NRMSE and coverage on the holdout; the predictive R̂ (the quantity setting it) and least ESS; D10's verdict (run; predictive, level, structural); M's median [95%] and largest, against the cap; σ's median; τ/τ₀'s median; the fit's seconds; the functional residual r₁'s mean against its expectation:

| benchmark | k | arm | NRMSE | coverage | pred. R̂ | pred. ESS | D10 | M | max M / cap | σ med. | τ/τ₀ med. | s | r₁ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman | 0 | A0 | 0.009887 | 0.990 | 2.368 (sigma2) | 5 | N (P N, L N, S N) | 24 [19, 31] | 32 / 499 | 0.0908 | 0.00863 | 1.3 |  |
| Friedman | 0 | A1 | 0.01742 | 0.984 | 1.625 (f3) | 7 | N (P N, L N, S N) | 16 [14, 19] | 20 / 499 | 0.257 | 1 | 0.8 |  |
| Friedman | 0 | A2 | 0.009746 | 0.989 | 3.069 (sigma2) | 5 | N (P N, L N, S N) | 28 [21, 31] | 32 / 499 | 0.078 | 0.00379 | 3.8 |  |
| Friedman | 0 | A3 | 0.01203 | 0.997 | 2.144 (sigma2) | 5 | N (P N, L N, S N) | 38 [29, 49] | 51 / 499 | 0.17 | 1 | 3.0 |  |
| Friedman | 0 | A4 | 0.008511 | 0.978 | 3.376 (sigma2) | 4 | N (P N, L N, S N) | 28 [22, 32] | 33 / 499 | 0.0656 | 0.00551 | 1.3 |  |
| Friedman | 0 | A5 | 0.01269 | 0.996 | 1.357 (f3) | 9 | N (P N, L N, S N) | 19 [17, 22] | 24 / 499 | 0.249 | 1 | 1.0 |  |
| Friedman | 0 | A6 | 0.007742 | 0.996 | 3.215 (sigma2) | 4 | N (P N, L N, S N) | 31 [28, 44] | 44 / 499 | 0.0605 | 0.00245 | 2.9 |  |
| Friedman | 0 | A7 | 0.01004 | 0.993 | 2.166 (sigma2) | 5 | N (P N, L N, S N) | 58 [43, 69] | 71 / 499 | 0.138 | 1 | 10.7 |  |
| Ishigami | 0 | A0 | 0.05982 | 0.949 | 2.161 (sigma2) | 5 | N (P N, L N, S N) | 38 [30, 42] | 45 / 499 | 0.216 | 0.114 | 1.2 |  |
| Ishigami | 0 | A1 | 0.06596 | 0.941 | 1.823 (f1) | 6 | N (P N, L N, S N) | 34 [29, 37] | 38 / 499 | 0.275 | 1 | 0.9 |  |
| Ishigami | 0 | A2 | 0.06996 | 0.934 | 3.156 (sigma2) | 4 | N (P N, L N, S N) | 39 [33, 51] | 52 / 499 | 0.21 | 0.00845 | 4.3 |  |
| Ishigami | 0 | A3 | 0.08171 | 0.950 | 2.207 (sigma2) | 5 | N (P N, L N, S N) | 49 [41, 55] | 58 / 499 | 0.37 | 1 | 7.9 |  |
| Ishigami | 0 | A4 | 0.0478 | 0.963 | 2.610 (sigma2) | 5 | N (P N, L N, S N) | 49 [42, 56] | 58 / 499 | 0.156 | 0.0849 | 1.5 |  |
| Ishigami | 0 | A5 | 0.06251 | 0.941 | 1.755 (f3) | 6 | N (P N, L N, S N) | 40 [37, 46] | 47 / 499 | 0.243 | 1 | 1.3 |  |
| Ishigami | 0 | A6 | 0.05628 | 0.932 | 2.762 (sigma2) | 5 | N (P N, L N, S N) | 58 [43, 68] | 73 / 499 | 0.148 | 0.00625 | 14.4 |  |
| Ishigami | 0 | A7 | 0.06773 | 0.967 | 2.124 (sigma2) | 5 | N (P N, L N, S N) | 84 [72, 101] | 108 / 499 | 0.289 | 1 | 68.1 |  |
| borehole | 0 | A0 | 0.008627 | 0.960 | 2.056 (f3) | 5 | N (P N, L N, S N) | 32 [29, 38] | 40 / 499 | 0.512 | 0.00393 | 0.8 |  |
| borehole | 0 | A1 | 0.01698 | 0.978 | 2.681 (f2) | 5 | N (P N, L N, S N) | 16 [16, 19] | 20 / 499 | 2.4 | 1 | 0.5 |  |
| borehole | 0 | A2 | 0.008156 | 0.958 | 2.094 (f3) | 5 | N (P N, L N, S N) | 40 [28, 44] | 46 / 499 | 0.421 | 0.0117 | 2.1 |  |
| borehole | 0 | A3 | 0.008628 | 0.981 | 1.719 (f3) | 6 | N (P N, L N, S N) | 43 [33, 47] | 50 / 499 | 0.808 | 1 | 4.6 |  |
| borehole | 0 | A4 | 0.007648 | 0.978 | 2.036 (f2) | 5 | N (P N, L N, S N) | 36 [32, 42] | 43 / 499 | 0.442 | 0.00335 | 1.0 |  |
| borehole | 0 | A5 | 0.01557 | 0.984 | 1.742 (f3) | 6 | N (P N, L N, S N) | 17 [15, 22] | 23 / 499 | 2.33 | 1 | 0.5 |  |
| borehole | 0 | A6 | 0.006266 | 0.963 | 2.493 (f1) | 5 | N (P N, L N, S N) | 46 [38, 49] | 51 / 499 | 0.343 | 0.00726 | 6.2 |  |
| borehole | 0 | A7 | 0.007811 | 0.987 | 1.607 (f2) | 7 | N (P N, L N, S N) | 70 [59, 76] | 79 / 499 | 0.681 | 1 | 19.6 |  |
| functional | 0 | A0 | 0.008497 | 0.958 | 1.913 (f2) | 6 | N (P N, L N, S N) | 31 [25, 38] | 39 / 1000 | 0.0016 | 0.00011 | 1.0 | 0.225 (-0.050) |
| functional | 0 | A1 | 0.004722 | 0.986 | 1.904 (f2) | 6 | N (P N, L N, S N) | 11 [9, 13] | 13 / 1000 | 0.0258 | 1 | 0.6 | 0.500 (-0.050) |
| functional | 0 | A2 | 0.003963 | 0.986 | 1.878 (sigma2) | 6 | N (P N, L N, S N) | 32 [30, 36] | 38 / 1000 | 0.00131 | 5.96E-05 | 4.3 | 0.120 (-0.050) |
| functional | 0 | A3 | 0.004438 | 0.905 | 2.123 (sigma2) | 5 | N (P N, L N, S N) | 76 [65, 85] | 89 / 1000 | 0.00798 | 1 | 24.3 | 0.217 (-0.050) |
| functional | 0 | A4 | 0.001133 | 0.966 | 1.674 (sigma2) | 6 | N (P N, L N, S N) | 36 [32, 45] | 46 / 1000 | 0.00145 | 0.000107 | 1.0 | 0.337 (-0.050) |
| functional | 0 | A5 | 0.007118 | 0.997 | 1.863 (f2) | 6 | N (P N, L N, S N) | 13 [11, 17] | 19 / 1000 | 0.0257 | 1 | 0.6 | 0.267 (-0.050) |
| functional | 0 | A6 | 0.000412 | 0.998 | 1.786 (sigma2) | 6 | N (P N, L N, S N) | 31 [27, 39] | 40 / 1000 | 0.00128 | 3.84E-05 | 6.2 | 0.141 (-0.050) |
| functional | 0 | A7 | 0.003932 | 0.916 | 2.068 (f2) | 5 | N (P N, L N, S N) | 108 [91, 126] | 128 / 1000 | 0.00715 | 1 | 92.0 | 0.302 (-0.050) |
| Friedman, σ = 1 | 0 | A0 | 0.06756 | 0.921 | 1.681 (f1) | 6 | N (P N, L N, S N) | 10 [8, 11] | 13 / 499 | 1.07 | 0.501 | 0.5 |  |
| Friedman, σ = 1 | 0 | A1 | 0.06647 | 0.921 | 1.482 (f1) | 8 | N (P N, L N, S N) | 9 [7, 13] | 14 / 499 | 1.08 | 1 | 0.5 |  |
| Friedman, σ = 1 | 0 | A2 | 0.06934 | 0.911 | 1.921 (f3) | 6 | N (P N, L N, S N) | 10 [9, 12] | 13 / 499 | 1.06 | 0.23 | 1.0 |  |
| Friedman, σ = 1 | 0 | A3 | 0.06148 | 0.960 | 1.556 (f3) | 7 | N (P N, L N, S N) | 10 [9, 13] | 16 / 499 | 1.09 | 1 | 0.8 |  |
| Friedman, σ = 1 | 0 | A4 | 0.06439 | 0.930 | 1.589 (f1) | 7 | N (P N, L N, S N) | 10 [9, 12] | 13 / 499 | 1.05 | 0.519 | 0.8 |  |
| Friedman, σ = 1 | 0 | A5 | 0.06283 | 0.931 | 1.572 (f2) | 7 | N (P N, L N, S N) | 11 [10, 13] | 15 / 499 | 1.07 | 1 | 0.4 |  |
| Friedman, σ = 1 | 0 | A6 | 0.05867 | 0.955 | 1.583 (f3) | 7 | N (P N, L N, S N) | 11 [10, 14] | 16 / 499 | 1.04 | 0.197 | 0.8 |  |
| Friedman, σ = 1 | 0 | A7 | 0.065 | 0.925 | 1.517 (f1) | 7 | N (P N, L N, S N) | 11 [10, 14] | 16 / 499 | 1.08 | 1 | 0.7 |  |
| Ishigami, σ = 0.37 | 0 | A0 | 0.1012 | 0.921 | 2.151 (f1) | 5 | N (P N, L N, S N) | 27 [23, 30] | 31 / 499 | 0.478 | 0.441 | 0.7 |  |
| Ishigami, σ = 0.37 | 0 | A1 | 0.09656 | 0.889 | 1.534 (f1) | 7 | N (P N, L N, S N) | 27 [25, 30] | 33 / 499 | 0.48 | 1 | 1.1 |  |
| Ishigami, σ = 0.37 | 0 | A2 | 0.104 | 0.880 | 1.668 (f3) | 6 | N (P N, L N, S N) | 24 [22, 29] | 30 / 499 | 0.47 | 0.0199 | 1.5 |  |
| Ishigami, σ = 0.37 | 0 | A3 | 0.1091 | 0.913 | 1.853 (f3) | 6 | N (P N, L N, S N) | 33 [30, 37] | 39 / 499 | 0.564 | 1 | 4.1 |  |
| Ishigami, σ = 0.37 | 0 | A4 | 0.1015 | 0.920 | 2.457 (f3) | 5 | N (P N, L N, S N) | 35 [24, 42] | 44 / 499 | 0.442 | 0.483 | 1.3 |  |
| Ishigami, σ = 0.37 | 0 | A5 | 0.09697 | 0.914 | 1.768 (sigma2) | 6 | N (P N, L N, S N) | 31 [27, 34] | 35 / 499 | 0.463 | 1 | 0.9 |  |
| Ishigami, σ = 0.37 | 0 | A6 | 0.09216 | 0.905 | 1.495 (f3) | 7 | N (P N, L N, S N) | 34 [29, 39] | 43 / 499 | 0.44 | 0.0412 | 7.5 |  |
| Ishigami, σ = 0.37 | 0 | A7 | 0.09496 | 0.933 | 1.408 (f3) | 8 | N (P N, L N, S N) | 49 [40, 68] | 72 / 499 | 0.508 | 1 | 22.4 |  |
| Friedman | 1 | A0 | 0.01235 | 0.975 | 2.047 (f3) | 5 | N (P N, L N, S N) | 24 [22, 30] | 31 / 499 | 0.0991 | 0.00935 | 0.8 |  |
| Friedman | 1 | A1 | 0.01907 | 0.967 | 1.983 (f3) | 5 | N (P N, L N, S N) | 16 [15, 19] | 19 / 499 | 0.266 | 1 | 0.7 |  |
| Friedman | 1 | A2 | 0.007664 | 0.990 | 2.823 (sigma2) | 5 | N (P N, L N, S N) | 29 [23, 35] | 36 / 499 | 0.0646 | 0.00347 | 1.6 |  |
| Friedman | 1 | A3 | 0.01212 | 0.981 | 1.844 (sigma2) | 6 | N (P N, L N, S N) | 39 [35, 48] | 50 / 499 | 0.165 | 1 | 5.2 |  |
| Friedman | 1 | A4 | 0.01103 | 0.985 | 3.338 (sigma2) | 4 | N (P N, L N, S N) | 27 [21, 32] | 33 / 499 | 0.0848 | 0.00767 | 1.2 |  |
| Friedman | 1 | A5 | 0.02004 | 0.946 | 1.688 (f2) | 6 | N (P N, L N, S N) | 17 [15, 20] | 21 / 499 | 0.258 | 1 | 0.8 |  |
| Friedman | 1 | A6 | 0.006869 | 0.987 | 2.696 (sigma2) | 5 | N (P N, L N, S N) | 39 [34, 49] | 51 / 499 | 0.0471 | 0.00243 | 5.4 |  |
| Friedman | 1 | A7 | 0.01082 | 0.995 | 2.104 (sigma2) | 5 | N (P N, L N, S N) | 52 [41, 65] | 68 / 499 | 0.152 | 1 | 15.3 |  |
| Ishigami | 1 | A0 | 0.05832 | 0.946 | 2.374 (sigma2) | 5 | N (P N, L N, S N) | 38 [28, 47] | 49 / 499 | 0.25 | 0.148 | 1.0 |  |
| Ishigami | 1 | A1 | 0.0652 | 0.943 | 1.906 (sigma2) | 6 | N (P N, L N, S N) | 33 [29, 36] | 37 / 499 | 0.314 | 1 | 1.3 |  |
| Ishigami | 1 | A2 | 0.06038 | 0.936 | 3.195 (sigma2) | 4 | N (P N, L N, S N) | 39 [32, 50] | 52 / 499 | 0.243 | 0.0115 | 3.0 |  |
| Ishigami | 1 | A3 | 0.06534 | 0.946 | 1.554 (sigma2) | 7 | N (P N, L N, S N) | 47 [40, 53] | 55 / 499 | 0.384 | 1 | 7.1 |  |
| Ishigami | 1 | A4 | 0.06182 | 0.945 | 2.804 (sigma2) | 5 | N (P N, L N, S N) | 48 [42, 56] | 58 / 499 | 0.203 | 0.151 | 0.8 |  |
| Ishigami | 1 | A5 | 0.06807 | 0.926 | 1.823 (sigma2) | 6 | N (P N, L N, S N) | 39 [33, 45] | 46 / 499 | 0.271 | 1 | 1.2 |  |
| Ishigami | 1 | A6 | 0.04538 | 0.963 | 2.141 (sigma2) | 5 | N (P N, L N, S N) | 54 [46, 65] | 67 / 499 | 0.162 | 0.0057 | 12.7 |  |
| Ishigami | 1 | A7 | 0.06268 | 0.948 | 1.721 (sigma2) | 6 | N (P N, L N, S N) | 74 [66, 83] | 88 / 499 | 0.314 | 1 | 43.9 |  |
| borehole | 1 | A0 | 0.008847 | 0.957 | 2.662 (f2) | 5 | N (P N, L N, S N) | 32 [28, 36] | 36 / 499 | 0.486 | 0.00369 | 0.9 |  |
| borehole | 1 | A1 | 0.019 | 0.985 | 2.072 (f1) | 5 | N (P N, L N, S N) | 15 [13, 18] | 20 / 499 | 2.34 | 1 | 0.6 |  |
| borehole | 1 | A2 | 0.007081 | 0.966 | 2.391 (sigma2) | 5 | N (P N, L N, S N) | 37 [33, 42] | 44 / 499 | 0.399 | 0.00828 | 3.5 |  |
| borehole | 1 | A3 | 0.009051 | 0.972 | 2.208 (f3) | 5 | N (P N, L N, S N) | 38 [34, 44] | 47 / 499 | 0.787 | 1 | 5.3 |  |
| borehole | 1 | A4 | 0.00751 | 0.963 | 2.534 (f2) | 5 | N (P N, L N, S N) | 37 [31, 43] | 44 / 499 | 0.41 | 0.00311 | 1.1 |  |
| borehole | 1 | A5 | 0.01575 | 0.984 | 1.493 (f3) | 7 | N (P N, L N, S N) | 17 [15, 21] | 22 / 499 | 2.29 | 1 | 0.8 |  |
| borehole | 1 | A6 | 0.006475 | 0.967 | 2.899 (sigma2) | 5 | N (P N, L N, S N) | 46 [39, 51] | 53 / 499 | 0.362 | 0.0117 | 6.1 |  |
| borehole | 1 | A7 | 0.006732 | 0.989 | 1.709 (sigma2) | 6 | N (P N, L N, S N) | 59 [51, 71] | 75 / 499 | 0.632 | 1 | 18.6 |  |
| functional | 1 | A0 | 0.00157 | 0.969 | 2.703 (f2) | 5 | N (P N, L N, S N) | 33.5 [26, 43] | 44 / 1000 | 0.00182 | 0.000156 | 1.6 | 0.165 (-0.050) |
| functional | 1 | A1 | 0.005021 | 0.999 | 2.295 (f2) | 5 | N (P N, L N, S N) | 10 [9, 17] | 18 / 1000 | 0.0263 | 1 | 0.5 | 0.376 (-0.050) |
| functional | 1 | A2 | 0.001921 | 0.960 | 2.505 (sigma2) | 5 | N (P N, L N, S N) | 31 [26, 36] | 38 / 1000 | 0.00168 | 0.000115 | 5.2 | 0.124 (-0.050) |
| functional | 1 | A3 | 0.004385 | 0.947 | 2.441 (sigma2) | 5 | N (P N, L N, S N) | 69 [64, 80] | 82 / 1000 | 0.00878 | 1 | 23.8 | 0.270 (-0.050) |
| functional | 1 | A4 | 0.004284 | 0.960 | 3.022 (sigma2) | 5 | N (P N, L N, S N) | 39 [36, 48] | 50 / 1000 | 0.00169 | 0.000158 | 1.4 | 0.144 (-0.050) |
| functional | 1 | A5 | 0.004441 | 0.999 | 2.431 (f2) | 5 | N (P N, L N, S N) | 11 [10, 19.02] | 21 / 1000 | 0.0261 | 1 | 0.9 | 0.474 (-0.050) |
| functional | 1 | A6 | 0.001081 | 0.988 | 2.477 (sigma2) | 5 | N (P N, L N, S N) | 37 [28, 55] | 56 / 1000 | 0.0014 | 0.00011 | 9.9 | 0.216 (-0.050) |
| functional | 1 | A7 | 0.003838 | 0.959 | 1.815 (sigma2) | 6 | N (P N, L N, S N) | 105 [97, 115] | 119 / 1000 | 0.00768 | 1 | 113.6 | 0.291 (-0.050) |
| Friedman, σ = 1 | 1 | A0 | 0.05098 | 0.949 | 1.798 (f3) | 6 | N (P N, L N, S N) | 11 [8, 12] | 12 / 499 | 0.948 | 0.39 | 0.6 |  |
| Friedman, σ = 1 | 1 | A1 | 0.0517 | 0.961 | 1.930 (f3) | 6 | N (P N, L N, S N) | 10 [9, 12] | 14 / 499 | 0.974 | 1 | 0.8 |  |
| Friedman, σ = 1 | 1 | A2 | 0.06131 | 0.915 | 1.473 (f3) | 8 | N (P N, L N, S N) | 10 [9, 11] | 13 / 499 | 0.972 | 0.236 | 1.2 |  |
| Friedman, σ = 1 | 1 | A3 | 0.06155 | 0.929 | 1.660 (f3) | 6 | N (P N, L N, S N) | 11 [10, 13] | 15 / 499 | 1 | 1 | 0.8 |  |
| Friedman, σ = 1 | 1 | A4 | 0.05893 | 0.949 | 1.546 (f3) | 7 | N (P N, L N, S N) | 12 [8, 15] | 17 / 499 | 0.968 | 0.445 | 0.9 |  |
| Friedman, σ = 1 | 1 | A5 | 0.05634 | 0.948 | 1.796 (f3) | 6 | N (P N, L N, S N) | 11 [10, 13] | 14 / 499 | 0.975 | 1 | 1.1 |  |
| Friedman, σ = 1 | 1 | A6 | 0.05372 | 0.937 | 1.458 (f3) | 8 | N (P N, L N, S N) | 12 [9, 15] | 18 / 499 | 0.95 | 0.221 | 1.7 |  |
| Friedman, σ = 1 | 1 | A7 | 0.05054 | 0.969 | 1.695 (f3) | 6 | N (P N, L N, S N) | 13 [10, 16] | 18 / 499 | 0.978 | 1 | 1.3 |  |
| Ishigami, σ = 0.37 | 1 | A0 | 0.08831 | 0.929 | 2.199 (f1) | 5 | N (P N, L N, S N) | 28 [21, 30] | 31 / 499 | 0.492 | 0.473 | 1.6 |  |
| Ishigami, σ = 0.37 | 1 | A1 | 0.08074 | 0.941 | 2.562 (f3) | 5 | N (P N, L N, S N) | 27 [24, 30] | 31 / 499 | 0.467 | 1 | 1.6 |  |
| Ishigami, σ = 0.37 | 1 | A2 | 0.0913 | 0.912 | 2.169 (f1) | 5 | N (P N, L N, S N) | 29 [27, 32] | 35 / 499 | 0.46 | 0.0276 | 3.2 |  |
| Ishigami, σ = 0.37 | 1 | A3 | 0.09481 | 0.897 | 1.657 (f3) | 6 | N (P N, L N, S N) | 32 [29, 38] | 41 / 499 | 0.567 | 1 | 8.2 |  |
| Ishigami, σ = 0.37 | 1 | A4 | 0.0727 | 0.935 | 1.864 (f1) | 6 | N (P N, L N, S N) | 29 [26, 32] | 34 / 499 | 0.428 | 0.397 | 0.9 |  |
| Ishigami, σ = 0.37 | 1 | A5 | 0.07999 | 0.918 | 1.584 (f1) | 7 | N (P N, L N, S N) | 33 [26, 37] | 39 / 499 | 0.444 | 1 | 2.2 |  |
| Ishigami, σ = 0.37 | 1 | A6 | 0.07715 | 0.925 | 1.847 (sigma2) | 6 | N (P N, L N, S N) | 35 [26, 41] | 45 / 499 | 0.429 | 0.0279 | 6.3 |  |
| Ishigami, σ = 0.37 | 1 | A7 | 0.07777 | 0.935 | 1.290 (f3) | 10 | N (P N, L N, S N) | 51 [45, 60] | 64 / 499 | 0.499 | 1 | 27.2 |  |
| Friedman | 2 | A0 | 0.01333 | 0.967 | 3.350 (sigma2) | 4 | N (P N, L N, S N) | 24 [19, 30] | 30 / 499 | 0.0823 | 0.00741 | 1.3 |  |
| Friedman | 2 | A1 | 0.01565 | 0.989 | 2.782 (f2) | 5 | N (P N, L N, S N) | 17 [16, 19] | 20 / 499 | 0.258 | 1 | 1.2 |  |
| Friedman | 2 | A2 | 0.009781 | 0.978 | 3.438 (sigma2) | 4 | N (P N, L N, S N) | 28 [20, 37] | 38 / 499 | 0.0675 | 0.00363 | 4.1 |  |
| Friedman | 2 | A3 | 0.01524 | 0.981 | 2.660 (sigma2) | 5 | N (P N, L N, S N) | 35 [28, 43] | 46 / 499 | 0.2 | 1 | 5.8 |  |
| Friedman | 2 | A4 | 0.007037 | 0.991 | 2.726 (f1) | 5 | N (P N, L N, S N) | 32 [26, 36] | 39 / 499 | 0.06 | 0.00439 | 1.4 |  |
| Friedman | 2 | A5 | 0.01316 | 0.996 | 1.724 (f2) | 6 | N (P N, L N, S N) | 18 [16, 20] | 21 / 499 | 0.252 | 1 | 1.1 |  |
| Friedman | 2 | A6 | 0.007253 | 0.996 | 2.054 (f1) | 5 | N (P N, L N, S N) | 31 [24, 36] | 38 / 499 | 0.0667 | 0.00324 | 5.2 |  |
| Friedman | 2 | A7 | 0.01112 | 0.990 | 1.949 (sigma2) | 5 | N (P N, L N, S N) | 56 [44, 69] | 74 / 499 | 0.148 | 1 | 21.1 |  |
| Ishigami | 2 | A0 | 0.05699 | 0.954 | 2.017 (f3) | 5 | N (P N, L N, S N) | 42 [34, 50] | 52 / 499 | 0.176 | 0.0952 | 2.1 |  |
| Ishigami | 2 | A1 | 0.08852 | 0.960 | 2.138 (sigma2) | 5 | N (P N, L N, S N) | 34 [30, 40] | 41 / 499 | 0.288 | 1 | 1.5 |  |
| Ishigami | 2 | A2 | 0.06484 | 0.947 | 2.455 (sigma2) | 5 | N (P N, L N, S N) | 41 [34, 53] | 55 / 499 | 0.168 | 0.00444 | 7.9 |  |
| Ishigami | 2 | A3 | 0.09014 | 0.952 | 1.557 (f1) | 7 | N (P N, L N, S N) | 45 [40, 53] | 57 / 499 | 0.381 | 1 | 12.2 |  |
| Ishigami | 2 | A4 | 0.06296 | 0.940 | 2.156 (f1) | 5 | N (P N, L N, S N) | 46 [39, 56] | 58 / 499 | 0.173 | 0.103 | 2.7 |  |
| Ishigami | 2 | A5 | 0.07329 | 0.956 | 1.902 (sigma2) | 6 | N (P N, L N, S N) | 42 [35, 47] | 49 / 499 | 0.265 | 1 | 2.1 |  |
| Ishigami | 2 | A6 | 0.06537 | 0.945 | 2.656 (sigma2) | 5 | N (P N, L N, S N) | 58 [47, 70] | 74 / 499 | 0.149 | 0.0089 | 26.9 |  |
| Ishigami | 2 | A7 | 0.0776 | 0.964 | 1.651 (f1) | 6 | N (P N, L N, S N) | 77.5 [69, 90] | 93 / 499 | 0.303 | 1 | 99.9 |  |
| borehole | 2 | A0 | 0.007734 | 0.971 | 2.230 (sigma2) | 5 | N (P N, L N, S N) | 33 [29, 38] | 39 / 499 | 0.458 | 0.00347 | 1.8 |  |
| borehole | 2 | A1 | 0.01558 | 0.976 | 1.341 (f1) | 9 | N (P N, L N, S N) | 16 [15, 17] | 20 / 499 | 2.28 | 1 | 0.9 |  |
| borehole | 2 | A2 | 0.006756 | 0.965 | 2.043 (f2) | 5 | N (P N, L N, S N) | 36 [33, 41] | 42 / 499 | 0.389 | 0.00841 | 3.8 |  |
| borehole | 2 | A3 | 0.009652 | 0.977 | 2.007 (sigma2) | 5 | N (P N, L N, S N) | 40 [34, 45] | 46 / 499 | 0.816 | 1 | 9.0 |  |
| borehole | 2 | A4 | 0.008006 | 0.967 | 2.410 (f1) | 5 | N (P N, L N, S N) | 37 [32, 40] | 41 / 499 | 0.449 | 0.00349 | 1.8 |  |
| borehole | 2 | A5 | 0.01787 | 0.987 | 1.725 (f2) | 6 | N (P N, L N, S N) | 18 [16, 21] | 23 / 499 | 2.32 | 1 | 1.2 |  |
| borehole | 2 | A6 | 0.00663 | 0.969 | 2.158 (sigma2) | 5 | N (P N, L N, S N) | 43 [39, 56] | 59 / 499 | 0.366 | 0.0103 | 11.7 |  |
| borehole | 2 | A7 | 0.00786 | 0.984 | 2.145 (sigma2) | 5 | N (P N, L N, S N) | 63 [55, 74] | 76 / 499 | 0.656 | 1 | 36.9 |  |
| functional | 2 | A0 | 0.001472 | 0.983 | 3.165 (sigma2) | 4 | N (P N, L N, S N) | 29 [26, 32] | 34 / 1000 | 0.00203 | 0.000129 | 1.6 | 0.298 (-0.050) |
| functional | 2 | A1 | 0.006637 | 0.920 | 1.841 (f1) | 6 | N (P N, L N, S N) | 9 [7, 11] | 12 / 1000 | 0.0292 | 1 | 0.9 | 0.523 (-0.050) |
| functional | 2 | A2 | 0.001182 | 0.989 | 2.222 (sigma2) | 5 | N (P N, L N, S N) | 31 [26, 45] | 46 / 1000 | 0.00176 | 0.000112 | 5.2 | 0.069 (-0.050) |
| functional | 2 | A3 | 0.004443 | 0.940 | 1.599 (sigma2) | 7 | N (P N, L N, S N) | 69 [58, 74] | 77 / 1000 | 0.0095 | 1 | 30.3 | 0.091 (-0.050) |
| functional | 2 | A4 | 0.001352 | 0.980 | 2.728 (sigma2) | 5 | N (P N, L N, S N) | 32 [25, 42] | 42 / 1000 | 0.00201 | 0.000134 | 2.5 | 0.180 (-0.050) |
| functional | 2 | A5 | 0.004814 | 0.987 | 2.177 (f1) | 5 | N (P N, L N, S N) | 11 [9, 16] | 17 / 1000 | 0.0291 | 1 | 0.7 | 0.548 (-0.050) |
| functional | 2 | A6 | 0.001022 | 0.991 | 2.963 (sigma2) | 5 | N (P N, L N, S N) | 41.5 [31, 52] | 54 / 1000 | 0.00162 | 0.000124 | 13.5 | 0.127 (-0.050) |
| functional | 2 | A7 | 0.004394 | 0.945 | 2.013 (sigma2) | 5 | N (P N, L N, S N) | 95 [72, 109] | 114 / 1000 | 0.00871 | 1 | 136.5 | 0.079 (-0.050) |
| Friedman, σ = 1 | 2 | A0 | 0.06321 | 0.950 | 1.758 (f2) | 6 | N (P N, L N, S N) | 11 [9, 13] | 15 / 499 | 1.03 | 0.49 | 0.9 |  |
| Friedman, σ = 1 | 2 | A1 | 0.06987 | 0.936 | 1.992 (f2) | 5 | N (P N, L N, S N) | 11 [10, 12] | 13 / 499 | 1.05 | 1 | 0.8 |  |
| Friedman, σ = 1 | 2 | A2 | 0.07087 | 0.891 | 1.231 (f2) | 12 | N (P N, L N, S N) | 11 [9, 13] | 14 / 499 | 1.03 | 0.232 | 1.3 |  |
| Friedman, σ = 1 | 2 | A3 | 0.06418 | 0.956 | 2.134 (f2) | 5 | N (P N, L N, S N) | 11 [10, 13] | 14 / 499 | 1.08 | 1 | 1.1 |  |
| Friedman, σ = 1 | 2 | A4 | 0.07291 | 0.914 | 2.308 (f2) | 5 | N (P N, L N, S N) | 11 [8, 13] | 15 / 499 | 1.03 | 0.469 | 0.8 |  |
| Friedman, σ = 1 | 2 | A5 | 0.06185 | 0.950 | 2.141 (f2) | 5 | N (P N, L N, S N) | 11 [9, 15] | 17 / 499 | 1.05 | 1 | 0.7 |  |
| Friedman, σ = 1 | 2 | A6 | 0.05767 | 0.965 | 1.726 (f2) | 6 | N (P N, L N, S N) | 10 [8, 12] | 13 / 499 | 1.04 | 0.164 | 1.2 |  |
| Friedman, σ = 1 | 2 | A7 | 0.06824 | 0.949 | 1.437 (f2) | 8 | N (P N, L N, S N) | 13 [10, 16] | 19 / 499 | 1.05 | 1 | 1.3 |  |
| Ishigami, σ = 0.37 | 2 | A0 | 0.09225 | 0.923 | 1.786 (f1) | 6 | N (P N, L N, S N) | 26 [24, 29] | 29 / 499 | 0.458 | 0.409 | 1.2 |  |
| Ishigami, σ = 0.37 | 2 | A1 | 0.08969 | 0.921 | 1.704 (sigma2) | 6 | N (P N, L N, S N) | 26 [23, 30] | 30 / 499 | 0.481 | 1 | 1.5 |  |
| Ishigami, σ = 0.37 | 2 | A2 | 0.1182 | 0.847 | 1.827 (f3) | 6 | N (P N, L N, S N) | 25 [22, 32] | 33 / 499 | 0.48 | 0.0228 | 2.8 |  |
| Ishigami, σ = 0.37 | 2 | A3 | 0.1038 | 0.907 | 1.802 (f3) | 6 | N (P N, L N, S N) | 33 [30, 38] | 41 / 499 | 0.564 | 1 | 5.8 |  |
| Ishigami, σ = 0.37 | 2 | A4 | 0.08361 | 0.929 | 1.802 (f1) | 6 | N (P N, L N, S N) | 30 [26, 34] | 35 / 499 | 0.43 | 0.401 | 1.9 |  |
| Ishigami, σ = 0.37 | 2 | A5 | 0.08696 | 0.937 | 1.951 (f1) | 5 | N (P N, L N, S N) | 30 [28, 34] | 36 / 499 | 0.459 | 1 | 1.2 |  |
| Ishigami, σ = 0.37 | 2 | A6 | 0.09771 | 0.918 | 2.412 (f1) | 5 | N (P N, L N, S N) | 32 [29, 39] | 42 / 499 | 0.427 | 0.0219 | 8.7 |  |
| Ishigami, σ = 0.37 | 2 | A7 | 0.09029 | 0.954 | 1.693 (f1) | 6 | N (P N, L N, S N) | 50 [43, 59] | 63 / 499 | 0.499 | 1 | 26.2 |  |

### Times

The whole fit's seconds, the mean over designs, and the reading's (prediction on the holdout) beside it:

| arm | Friedman | Ishigami | borehole | functional | Friedman, σ = 1 | Ishigami, σ = 0.37 |
| --- | --- | --- | --- | --- | --- | --- |
| A0, g, Zellner–Siow, λ (10, 10) | 1.1 (18.5) | 1.4 (20.3) | 1.2 (16.6) | 1.4 (6.1) | 0.7 (15.5) | 1.2 (19.4) |
| A1, g, τ held at 1/N, λ (10, 10) | 0.9 (16.6) | 1.2 (20.2) | 0.7 (13.8) | 0.7 (6.1) | 0.7 (15.4) | 1.4 (22.1) |
| A2, ridge, τ sampled, λ (10, 10) | 3.1 (19.1) | 5.1 (19.0) | 3.1 (17.0) | 4.9 (6.3) | 1.1 (15.6) | 2.5 (16.5) |
| A3, ridge, τ held at 1, λ (10, 10) | 4.6 (22.8) | 9.1 (20.8) | 6.3 (17.6) | 26.1 (7.5) | 0.9 (14.2) | 6.0 (22.0) |
| A4, g, Zellner–Siow, λ (1, 1) | 1.3 (17.5) | 1.7 (19.1) | 1.3 (18.6) | 1.6 (5.9) | 0.9 (14.8) | 1.4 (17.9) |
| A5, g, τ held at 1/N, λ (1, 1) | 1.0 (16.5) | 1.5 (17.0) | 0.9 (15.8) | 0.8 (5.9) | 0.7 (15.9) | 1.5 (22.3) |
| A6, ridge, τ sampled, λ (1, 1) | 4.5 (20.8) | 18.0 (22.1) | 8.0 (18.3) | 9.9 (7.8) | 1.2 (14.4) | 7.5 (19.8) |
| A7, ridge, τ held at 1, λ (1, 1) | 15.7 (19.8) | 70.6 (24.9) | 25.0 (20.6) | 114.0 (9.4) | 1.1 (14.9) | 25.3 (20.5) |

The fits: 20.0 minutes timed; the readings 39.5 minutes; CPU over wall 0.770.

Fits with a draw at the basis cap: none.
