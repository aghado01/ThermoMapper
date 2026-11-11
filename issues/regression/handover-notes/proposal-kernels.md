# E4's remaining proposal kernels (KHAOS chip, from 2026-11-08)

- Plan: `issues/regression/bass/proposal-kernels-brief.md`, committed 702d427. Waiting for the coordinating session's verification before slice 1.
- Shape: a birth mixture through F-X4, q = (1 − ω_H − ω_R)·q_base + ω_H·q_H + ω_R·q_R. The base keeps a positive share, so support is kept. Deaths are uniform and changes are the base kernel's. Off by default, so hashes hold.
- Slice 1, the history birth (Nott and Kohn 2005):
  - weights max(ĥ_u, 0.01), where ĥ_u is u's inclusion frequency over the chain's past;
  - they replace z_u (hinge, for every order, K = 1 included) or make_weights' (aⱼ^α + ε) (Legendre);
  - refreshed every 100 sweeps in burn-in and frozen at its end (D10), per chain and per rung, and carried in checkpoints;
  - under SMC, the weights come from each stage's population. Nott and Kohn's best linear predictor is recorded, not built.
- Slice 2, the residual birth:
  - a parent uniform over the eligible terms plus the intercept, then a new factor ∝ exp(½·b·Δℓ);
  - Δℓ = −(n/2)ln(1 − (cᵀr)²/(‖c‖²‖r‖²)), with c = B_p∘φ_f;
  - hinge knots are scored by cumulative sums;
  - the exact q sums over t's sub-terms in S, and a death's reverse needs r(S′) from the candidate (an F-X4 context change).
- Gates: Σq = 1 and χ² on draws; exact detailed balance against the enumeration; the chain's law; Geweke; resume.
- Measurements:
  - ESS per second on Friedman (p 10), Ishigami (both priors), and Friedman at p 50, which the plan adds;
  - for the residual birth, the trap fixture's reach (engine schedule and the 1.25 cap) and MCMC's climb from empty.

## Slice 1 landed (2026-11-08/09): b1df14e (code), 6454c5e (docs, measurement, slow tier)
- Built as planned, except one history per chain under a ladder (the b = 1 rung's, shared by all rungs; survives resizes).
- Gates all pass:
  - exact detailed balance, worst 5.7e-14;
  - kernel law: Σq = 1, ln q vs statement 9e-16, χ² z −1.39..0.29;
  - chain law after a freeze: |z| 2.07 / 2.21;
  - Geweke: 2.32 / 1.45;
  - ladder, SMC (−72.2912 ± 0.0137 vs −72.2910), flat histogram, resume bit for bit;
  - Commit scope 1,874/1,874; whole assembly 2,676/2,676.
- Measurement (pre-registered, 4 seeds × 4 chains, 5k + 5k, ESS per CPU s): no gain across D10's groups.
  - Friedman p10: M 0.73× (z −4.2).
  - Friedman p50: f0 1.50× (z 2.1).
  - Ishigami: within the seeds' spread.
  - Result: opt-in, no default candidate.
- Trap to remember: tiered JIT polluted the first timed fit even after a short warm-up; measurements run with DOTNET_TieredCompilation=0.
- Slow tier rewritten at 16 s (profile ran ~1.6× slow under load): 80 methods.
