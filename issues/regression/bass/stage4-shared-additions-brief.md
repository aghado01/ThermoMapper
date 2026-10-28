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

### The three additions (2026-10-28)

Commits: `a2034ea` (Task A), `e0995be` (Task B), `5ce7453` (Task C). Suites at `5ce7453`: Distributions 119/119, Samplers 144/144, Regression 499/499. The Regression count includes the main session's uncommitted `BassDesignTests`, which compiled in the shared tree.

**BARS bit-identical.** The seed sweep (`SEED_SWEEP=50`) matches the references byte for byte, table and summary:
- Regression at `e0995be`: `20261028_035946474…` against `015051357…`, 1,370 rows.
- Samplers at `5ce7453`: `20261028_040411118…` against `015056249…`, 154 rows.
- Regression at `5ce7453`: `20261028_041239296…` against `015051357…`, 1,370 rows.

**Task A: Wallenius (`Maths.Distributions.Wallenius`).**
- *API.* `LogProbability(weights, set, method = Auto)`, `Probability(…)` and `Sample(rng, weights, set)`, with `WalleniusMethod` {`Auto`, `DrawOrders`, `Quadrature`}. The weights are all p of them; the set holds K distinct indices, in any order, K ≥ 1.
- *A deviation from the brief and from D01: no inclusion–exclusion.* The default evaluation is the brute-force sum over the K! draw orders, accumulated by a recursion over the 2^K subsets of S: G(U) = Σ_{i∈U} zᵢ·G(U∖i)/(R + Z_U), with R the weight left out of S.
  - Every term and every denominator is a sum of positives, so it is exact to a few ulps at any weights, in O(K·2^K). No error estimate is needed, and the choice the brief asked for goes.
  - Scaled by Π min(zᵢ/R, 1) it stays within [1/(K+1)!, K!], so log P keeps its digits below double range (tested at P ≈ 10⁻³⁵⁰).
  - *Why.* The brief's fallback, quadrature, costs 40–85 µs a call against 0.1–0.5 µs for the recursion at K ≤ 5 (Release). Inclusion–exclusion's own estimate, ε·Σ|terms|/P, already rejects at p = 50 with equal weights at K = 2 (about 10⁻¹²). So births at moderate p would have paid the quadrature on every proposal and on every death's reverse.
- *R's alternating sum is not ported.* It loses up to 7.7·10⁻³ under F-O5's five decades, which the tests assert corrected. Where R is right (R's weights, K = 2 and 3) the recursion agrees with it to its own error.
- *Quadrature.* Adaptive G7–K15 of Chesson's positive integrand on v = ln u, relative 10⁻¹⁴, in log space. It is selectable, and `Auto` takes it above K = 10. The explicit recursion runs to K = 20.
- *Boundaries.* K = 1 is zᵢ/Σz, by −log1p(R/zᵢ); R's `dmwnchBass` doubles it there. The whole set has probability 1.
- *Tests.* `WalleniusTests`, 26 (the commit message says 24):
  - brute force over the draw orders at R's default weights and at five decades, K = 1…5, both evaluations;
  - sums to one over every K-subset;
  - 1/C(p, K) at p = 50 and at K = 11;
  - F-O5's values, pinned from `bass_wallenius_oracle.R` runs;
  - the two evaluations against each other at twelve decades;
  - P near 1 and below underflow;
  - the sampler's law;
  - refusals.
- *For S3.* `Wallenius.LogProbability(z, set)` with z = w₂ + appearances is the variable-set factor of a birth, and of a death's reverse. `Wallenius.Sample` draws that law, if the birth wants it from here. D01's bullet ("inclusion–exclusion up to K = 3 and quadrature above") no longer describes the code; amending it is the fitting segment's business.

**Task B: the flat intercept, in the block's centred form.**
- *Shape: a factory over the centred design, selected by the convention.* `CollapsedGPrior.FromCentredLeastSquares(centredGram, centredCoefficients, residualSumOfSquares, centredFittedSumOfSquares, rows, responseMean, intercept)`.
  - *Why this shape.* S4's `BassDesign` holds the factor of G_c alone, with ȳ apart. Its `Factor`, `Coefficients`, `ResidualSquares`, `FittedSquares`, `Observations` and `ResponseMean` map one to one onto the arguments, and the config carries the convention as a `GPriorCentre`. `FromLeastSquares` would have needed the intercept bordered into the factor.
  - *Beyond the brief.* The factory takes `Mean` and `Zero` as well as `Flat`, so S5 builds every convention's block alike. Through it, `Mean` equals `FromLeastSquares` on the explicit [1, Z_c] with μ₀ = ȳ·e₀ (10⁻¹²). `Zero` adds Nȳ² to Q itself, so pass ỹ's statistics under every convention. `Median` is refused.
  - The coefficients are (α, β_c), the intercept first (`Columns` = M + 1). `Intercept` reports the convention; it is null from the existing factories.
- *The flat form.* k = M, and ν = Nb − 1, which `GPriorMarginal.DegreesOfFreedom` carries into σ²'s shape and the 2π term. The determinant term is (M/2)·ln(τ/(τ+b)) − ½·ln(Nb).
  - b = 0 throws `ArgumentOutOfRangeException`. The linear form there is the consumer's: b times the b = 1 marginal.
  - The improper σ² prior needs ν > 0, and a proper one needs g₁ + ν/2 > 0; both refuse otherwise.
  - The marginal is relative, so read `block.Intercept == GPriorCentre.Flat` as "no evidence".
- *Conditionals.* σ² | y ~ IG(g₁ + ν/2, g₂ + S/2).
  - `Sample` draws σ², then β_c through the factor, then α ~ N(ȳ, σ²/(Nb)), apart from β_c. `InterceptCovarianceScale` is c₀ = 1/b under Flat and 1/(τ+b) inside the g-prior.
  - `LeaveOutLogDensity` and `ConditionalVariance` take the design [1, Z_c], column 0 the intercept's, and carry its variance. The leave-out's Student t has ν = N − 1 under Flat.
- *Existing paths.* `FromLeastSquares` and `FromWeighted` are unchanged bit for bit; the sweep and the existing bit-for-bit test say so.
- *BARS.* BARS still refuses `Flat`. Its reason in `BandedGPrior.cs` is reworded, keeping the phrase `NormalCentreTests` reads, so that test is unchanged.
- *Docs.* `BayesianRegression.md`: the block section and its checks, the "Left open" note and the alternatives bullet.
- *Tests.* `CollapsedGPriorInterceptTests`, 28:
  - dense algebra on [1, Z_c] at four (τ, b), and on the uncentred [1, H] for the two g-prior conventions;
  - level centring against the explicit column;
  - the flat marginal and E[σ²] against quadrature over α and σ², at b = 1, 0.5 and 0.2, under both priors;
  - Liang's and level centring's Bayes factors under Zellner–Siow against F-O9's oracle on F-O8's grid (values pinned from a grid run), to 10⁻⁸;
  - on a hinge design, shift invariance at offsets 10 and −250 (the zero centre moves), with the g-forms;
  - draws, leave-out densities and readout variances against dense counterparts;
  - `MarginalAt` bit for bit, without allocating;
  - refusals.

**Task C: per-replica blocks in `ParallelTempering`.**
- *The contract is `IBlock` itself, over the replica chain: `IBlock<ReversibleJumpChain<TState>>`.*
  - A block reads `chain.Current`, `chain.CurrentEvaluation` and `chain.Temperature` (its rung now, after any adaptation). It updates the state leaving `Temperature.LogDensity` invariant, with log τ's Jacobian if it slices on log τ, and hands the state back with its evaluation through `SetState`, so the swaps read a current evaluation.
  - *Why over the chain, not the state.* A reversible-jump state is immutable and shared, by swaps and by the common start, and the chain caches its evaluation; so the block works through the chain. `IBlock` over a mutable model state (BAPS, HPS) is unchanged.
  - *One contract for both cases.* `MoveStep<T>`, the chain's own step as a block, lets a lone chain run the same blocks in `BlockScan<ReversibleJumpChain<T>>([new MoveStep<T>(), tau], chain, rng, c => c)`. BARS's pattern, the engine's run object holding the chain, serves too.
- *API.*
  - The constructor and `Restore` take `blocksForReplica`, called once per replica, as `movesForReplica` is. `Restore` refuses a mismatch.
  - `HasBlocks` and `BlockStats(replica)` report them.
  - `LadderSnapshot` gains `BlockRngs`, `BlockAttempts` and `BlockAccepts`, all null without blocks.
- *Order and streams.*
  - Each scan runs every replica's moves, then its blocks in order on their cadence (counted by scan from 1), then the swaps. The exact reference rung runs none.
  - The streams are `SeedTree.Derive(masterSeed, 2L + 2)`: the first L + 2 as before, and replica i's block stream at L + 2 + i. A ladder without blocks derives only L + 2.
  - Mutable block state, such as an adaptive width, is the engine's to save, as for `BlockScan`.
- *Docs.* `Samplers.md`: the API rows, a "Per-replica blocks" paragraph, the step order and the source map.
- *Tests.* `ParallelTemperingBlockTests`, 7. The state is (k, θ): an index moved only by reversible-jump moves, and θ only by a slice block.
  - Each rung's P(k) and E[θ | k] against enumeration and quadrature, on both paths.
  - Each pair's swap rejection against its expectation under the product of the rung laws, to within 10⁻³ of it (tolerance 0.01). Stale evaluations or an unbalanced swap would move it.
  - Bit-identical resume, and the mismatch refusals.
  - A block that draws nothing leaves the trajectory bit for bit; it keeps its cadence, and does not run on the reference rung.
  - The lone chain through `BlockScan` samples the posterior.
  - Validation.

**Notes.**
- Each seed sweep writes a new timestamped directory under `artifacts/test-runs/SeedSweep/`. Those runs (`035946474`, `040411118`, `041239296`) were left there as evidence.
- A hand-over note for S3, S5 and S8 sits in the session group's shared notes folder (`stage4-shared-additions-handover.md`).
