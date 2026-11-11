# Stage-4 shared additions: hand-over for S3, S5 and S8 (2026-10-28)

From the chip running `issues/regression/bass/stage4-shared-additions-brief.md`. Full report: that brief, "Reports". Commits: `a2034ea` (A, Wallenius), `e0995be` (B, the flat intercept), `5ce7453` (C, per-replica blocks).

- **S3, Wallenius.** `Maths.Distributions.Wallenius.LogProbability(z, set)`: z the p weights, set the K variable indices in any order, K ≥ 1. `Wallenius.Sample(rng, z, set)` draws the law if the birth wants it.
  - **Deviation from D01 and the brief:** no inclusion–exclusion. The default is the sum over draw orders, done by a recursion over the 2^K subsets. Every term is positive, so it is exact at any weights, and it costs 0.1–0.5 µs at K ≤ 5, against 40–85 µs for the quadrature.
  - The quadrature stays selectable, and Auto takes it for K > 10. D01's bullet ("inclusion–exclusion up to K = 3 and quadrature above") no longer describes the code; the fitting segment is the main session's to amend.
- **S5, the flat intercept.** `CollapsedGPrior.FromCentredLeastSquares(design.Factor, β̂_c, design.ResidualSquares, design.FittedSquares, N, design.ResponseMean, config.Intercept)` takes Mean, Zero and Flat. Its arguments map one to one onto `BassDesign`'s members.
  - The coefficients are (α, β_c), the intercept first, so `Columns` = M + 1.
  - Under Flat the marginal is relative, so it gives no evidence (`block.Intercept == GPriorCentre.Flat`), and b = 0 throws: the linear form there is b·ℓ₁.
  - Pass the statistics of ỹ under every convention: Zero adds Nȳ² itself.
- **S8, per-replica blocks.** τ's slice is an `IBlock<ReversibleJumpChain<BassState>>`. It reads `chain.Current`, `chain.CurrentEvaluation` and `chain.Temperature`, leaves `Temperature.LogDensity` invariant (with log τ's Jacobian if it slices on log τ), and calls `chain.SetState(next, evaluation)`.
  - Pass `blocksForReplica` to the ladder's constructor and to `Restore`.
  - Untempered, the same block runs in `BlockScan<ReversibleJumpChain<BassState>>([new MoveStep<BassState>(), tau], chain, rng, c => c)`, or in the engine's run object, as BARS does.
