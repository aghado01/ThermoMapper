# Flat-histogram sampler, Phase 1 (N17): hand-over, 2026-10-29

Full report: `issues/regression/bass/flat-histogram-phase1-brief.md`, under Reports ("Phase 1 built").

- Commits on `main`, none pushed: `d7f4ab8` (weights + target wrapper), `3627a13` (SAMC learning + run), `54f417b` (reweighting + ln Z(b)), `2fa3bc0` (WHAM warm start, spec, resume), `f6412f1` (report). 19 new files under `src/maths/samplers/flathistogram/` and `tests/samplers/`; no existing file changed; `Temperature` untouched.
- After the report: `4c34930` (WHAM's Ferrenberg–Swendsen 1/g weights), `ffe6bb5` (report), `27fd047` (proper-prior correction), `9f561d0` (FixedPoint's Newton start: log-linear per cell, refuses unflat solves; the old one failed on ln g spanning >745 nats or steep against the spacing), `001357b` (mean-field Potts toy + reading behind FLAT_HISTOGRAM_MEANFIELD_READING), `3c95d9b` (its report).
- Mean-field reading (q = 4, N to 2048, equal proposals): flat histogram α ≈ 2.2 power law and the only sampler still crossing at N = 2048 (ln peak/dip 12.2); fixed b_t exponential (power law predicted 17.7 and 14.1 round trips where none came); tempering's replicas stop crossing from N = 512 (not distinguished, 3 sizes); q = 3 too weak to separate.
- Suites at `2fa3bc0`: Samplers 251/251 (45 new), Regression 1 303/1 303.
- SAMC update with hats: θ_k += γ_t·(1/C)·Σ_c 1_in·(φ_k − π_k); outside visits inert; limit = flat conditional hat occupancy (convex G). `SamcLearning.FixedPoint` gives that limit exactly for a known g (the tests' reference).
- Toys: double well (quadrature + closed-form ln g), 6-state Potts 3×3 (6^9 enumerated; bimodal at b_c, peak/dip 2.37). Learned ln w to θ* 0.043 (double well, 64k) / 0.040 (Potts); ln Z(b) and laws within the chains' errors (12-seed calibration RMS z 1.1–1.4 vs t₇ 1.18).
- For Phase 2: the ratio estimator needs a proper prior only, its normalizer cancelling, so BASS's unnormalized structure prior needs no correction (corrected 2026-10-30; adding ln C + ln Z_F would be a ~0.76-nat error); the pilot must be on the linear ℓ_b = b·ℓ₁ form (`FromLadders` refuses others); BASS needs its own `IMulticanonicalChain` (blocks must target π_w); `GainStart` default (pilot draws / chains) is regime-dependent; ComparativeSweep fold of 3 claims parked in `tests/samplers/FlatHistogramSweep.cs`.
- Exploration: on the 3×3 Potts, adapted PT round trips (1 640 / 10⁶ proposals) beat the flat histogram's (760) at equal work.
- `artifacts\chip-fh\` deleted.
