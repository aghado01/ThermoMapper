# E34 slice 5: the PCA mode's joint option (chip, 2026-11-10)

Commits: `4135eb9` (code and `MultiOutputBasisTests`), `cb4849b` (report in `multi-output-brief.md`, decision 2 recorded,
notation ledger). The decision 2 log hint is its own commit after them (see the report to the main session).

## For other chips

- `BassComponents.FitJoint`/`AssembleJoint`: a decomposition's r scores as r outputs of one multi-output fit.
  `BassBasisFit.Joint` is that fit; `Components` are views of it (scalar `BassFit` records: column k of the
  coefficients, σ² = Σ_kk, config without `OutputCovariance`, own fingerprint on `ComponentData(k)`, no evidence or
  checkpoint). Every component shares one structure per draw, and `JointDraws.Draw(k, s)` is the same for every k.
- Readers of the components' surfaces need nothing new. Readers of the components' noise must read Σ whole when
  `Joint` is set: `BassBasisFit.NoiseCovariance(jointDraw)` / `NoiseRoot(jointDraw)` (internal). Prediction, the
  pointwise density and the calibration emulator do; a new noise reader on basis fits should too.
- Coefficient-space pointwise units must hold whole cases under the joint option.
- Rank-1 joint fit = separate fit draw for draw under MCMC; not under SMC (several outputs' SMC takes its own path at
  K = 1; follow-up recorded for slice 2's owner).
- Draw-for-draw comparisons across mixed bases break at exact ties (equal-span hinges on gridded inputs): rounding
  decides whether `ReversibleJumpChain` draws a uniform. Gate on continuous inputs or read the target per state.
