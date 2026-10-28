# BASS stage 4: the shared additions (chip brief)

Date: 2026-10-28. A chip of the [fitting segment's stage-4 plan](port-fitting-brief.md#stage-4-plan-decided-2026-10-28), which the user decided on 2026-10-28. The plan names three additions to shared code that BASS's steps need, each generic, tested in its own host, and landed before its first BASS use. This chip builds them in the order the main session needs them: Wallenius for S3, the flat-intercept form for S5, the per-replica block for S8. The main session runs S4 meanwhile.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md) (zero-dependency C#, net10.0, lean docstrings, no compatibility shims, tests assert corrected references), then the fitting segment's D01, D03, D04 and D05.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch.
- **Build and test into your own artifacts root.** The main session builds the default root at the same time, so add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-shared\` to every `dotnet build` and `dotnet test`, and do not run `scripts/fact-harness.ps1`, which builds the default root. Delete `artifacts\chip-shared\` when you finish.
- **Touch only your tasks' files** and the new test files. The main session edits `src/maths/regression/spline/bass/**`, `tests/regression/Bass*.cs`, `tests/oracle/Bass*.cs`, the fitting segment and its checklist. The oracle chip edits `r/oracles/**`, `r/scripts/**` and its own brief.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per task, in the repo's message style, ending with the `Co-Authored-By` line. Do not push.
- **Line endings:** `docs/Maths/*.md` are CRLF; briefs and C# are LF. Check `git diff --stat` after a scripted edit.
- **Report** under "Reports" at the end of this brief, newest first: what landed, the commits, the tests and suite counts, and anything the main session must know (an API choice, a deviation, an open question).

## Task A: Wallenius' variable-set probability (Maths.Distributions)

- **What.** The probability that weighted sampling without replacement draws the set S of K distinct variables, each draw choosing among the remaining variables in proportion to its weight z: Wallenius' multivariate noncentral hypergeometric with one ball per variable, P(S) = ∫₀¹ Πᵢ∈S (1 − t^{αᵢ}) dt with αᵢ = zᵢ / Σ_{k∉S} zₖ (D01's "variable-set probability"). BASS's births need log P(S), for K ≥ 2.
- **How.** Inclusion–exclusion gives P(S) = Σ_{T⊆S} (−1)^{|T|} / (1 + Σᵢ∈T αᵢ), exact in exact arithmetic; R's `dmwnchBass` evaluates it so. The alternating sum loses digits under skewed weights: the reference-fixtures report measured 1·10⁻⁷ at K = 2 and 3·10⁻⁵ at K = 3 with weights spanning five decades, and recommended choosing by an error estimate rather than by K alone. So use inclusion–exclusion where its cancellation error estimate is small, and quadrature of the positive integrand otherwise (the plan's rule, refined by that finding). Handle S holding every variable (P = 1) and K = 1 (zᵢ / Σz; `dmwnchBass` doubles it there, a path R never takes).
- **Where.** `src/maths/distributions/Wallenius.cs`; tests under `tests/distributions/`.
- **Checks.** Against brute force over the K! draw orders (Σ over orders of Π z / remaining weight), at R's default weights (w₂ = 5 plus counts) and at weights spanning five decades, K = 1 … 5; the probabilities over every K-subset sum to one; and F-O5's values (`r/oracles/bass_wallenius_oracle.R`).

## Task B: the flat-intercept form in the shared g-prior block (Maths.Regression)

- **What.** `GPriorCentre.Flat`, Liang et al.'s (2008) convention, exists since F-C2 and BARS refuses it (`BandedGPrior.RefuseFlat`). BASS selects it as an option (D03), so `CollapsedGPrior` needs its closed forms: the intercept takes a flat prior and is integrated out, the other columns are centred and take the g-prior, and k = M. The intercept spends one degree of freedom, so σ²'s shape reads N − 1 where level centring reads N; under the Jeffreys limit the marginal is (1 + g)^{(N−1−M)/2} (1 + g(1 − R²))^{−(N−1)/2}, g = 1/τ (D03's propriety bullet).
- **Tempering.** On the likelihood path at inverse temperature b, the intercept's integral contributes (2πσ²/(bN))^{1/2}, so the shape reads (Nb − 1)/2 and b = 0 is out of reach on the integrated form; the linear form reaches it (D05). The absolute marginal is undefined under the flat prior, so the form reports no evidence; within one fit, differences between structures are well defined.
- **Conditionals.** D04 composes a, σ² and λ at retention, so the flat form also needs σ²'s conditional and a's draw: with centred columns orthogonal to 1, the intercept's conditional is N(ȳ, σ²/(Nb)) given σ², apart from the other coefficients.
- **Shape.** BASS's S5 will build blocks through `CollapsedGPrior.FromLeastSquares` from bordered factors of the centred basis. Choose how the flat form is selected (a flag, a form parameter, a separate factory) with that use in mind, record why, and keep every existing path bit-identical.
- **Where.** `src/maths/regression/gprior/CollapsedGPrior.cs` and its tests under `tests/regression/`; the block section of `docs/Maths/BayesianRegression.md` (CRLF), including its "Left open: a flat intercept" note.
- **Checks.** F-O9's identities (`r/oracles/bass_intercept_oracle.R`): the level-against-Liang identity on its explicit design, Liang's Bayes factors on F-O8's grid, and invariance to shifting y; the marginal against quadrature over g and over σ²; the tempered form at several b against quadrature; BARS still refuses `Flat`; the Regression suite unchanged elsewhere.

## Task C: a per-replica block in parallel tempering (Maths.Samplers)

- **What.** `ParallelTempering` takes a reversible-jump move palette per replica (`movesForReplica`), and nothing else. BASS's chain state is (structure, τ), and D04 updates τ by slice sampling on log τ, a block that is not a reversible-jump move. Under tempering each replica must run that block at its own b (D05: τ's tempered conditional). Add a generic way to run, on each replica, a kernel that leaves that replica's tempered target π_b invariant, between its moves and the swaps.
- **Constraints.**
  - Generic: no BASS, BARS or spline types in Samplers.
  - Its draws come from a further prefix-stable `SeedTree` stream, so a ladder without a block keeps every stream bit for bit, as the exact reference rung did (`92540de`).
  - Checkpoints and snapshots carry what the block needs, so a resumed run is bit-identical to an uninterrupted one.
  - Consider the untempered case too: `BlockScan` composes `IBlock`s for a lone chain, and one block contract serving both would be best. Say what you chose and why.
- **Where.** `src/maths/samplers/rjmcmc/ParallelTempering.cs` and what it needs; tests under `tests/samplers/`.
- **Checks.**
  - On an analytic target with a discrete index moved by reversible-jump moves and a continuous coordinate moved by the block, each rung's law matches its tempered target, by enumeration and quadrature.
  - Swaps keep detailed balance.
  - Resume is bit-identical.
  - Every existing parallel-tempering test passes unchanged, and BARS stays bit-identical: the Samplers and Regression suites, with the seed sweep's per-seed values unchanged (`SEED_SWEEP=50`, compared against `artifacts/test-runs/SeedSweep/20261028_015051357…` and `20261028_015056249…`).

## Reports
