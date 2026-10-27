# BASS port: shared numerics and sampler contracts: checklist

Companion to [port-numerics-brief.md](port-numerics-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Decisions

- [ ] **N-D10** Adaptation policy for adaptive kernels (note): drafted 2026-10-25 as a proposal; settles at the user's design pass.

### Build (stages 2–3)

- [x] **N1** Bordered insertion on the strict dense factor, with candidate work private until acceptance: `BorderedCholesky`.
- [x] **N2** Givens deletion, and a changed term as a staged delete and insert: `BorderedCholesky`.
- [x] **N3** Strict rebuild: periodic and on failure, with drift checks against fresh factorizations: `BorderedCholesky.Rebuild`, `FactorDrift`.
- [x] **N4** diag(G⁻¹) maintained through births and deaths, for D06's order-free gate: `BorderedCholesky.LeaveOneOutPivot`.
- [x] **N5** Stable thin SVD with numerical-rank handling in LinAlg: `ThinSvd`.
- [x] **N6** Diagonal-plus-low-rank solves and log-determinants, without m × m inverses: `DiagonalPlusLowRank`.
- [x] **N7** Face-splitting design contraction with the Cartesian (Hadamard) shortcut: `KhatriRao`, both intercept conventions.
- [x] **N8** Generic SMC core in `samplers/smc` (optional; specified in the calibration segment): `SmcRun`, `SmcSampler`.
- [x] **N9** Qualified baseline rerun after the shared changes, before the first BASS integration tests: requalified at `3a04601`.
- [x] **N10** Cross-products summed pairwise or with compensation: `PairwiseSum`. F-O7's near-constant design is the acceptance case: its left-to-right error is 10⁻³ relative in s.
- [ ] **N11** Evidence over a tempering ladder: path sampling and stepping-stone from per-rung traces, with Monte Carlo error (D03, D05). Chip: [shared-infrastructure-brief.md](shared-infrastructure-brief.md).
- [ ] **N12** The g-prior's log marginal at (τ, b) as one scalar, without the mean vector (D04's τ slice, D05's swaps). Same chip.
- [ ] **N13** External reference for the rank-normalized and PSIS diagnostics against R `posterior` and `loo` (OQ9; D10's gate). Same chip; needs the user's approval for the download.

### Room for extensions

- [ ] **N-X1** The face-splitting contraction covers tensor designs (ABMRS).
- [ ] **N-X2** `ISpdFactor` admits a wide-band factor beside `DenseCholesky`.
- [ ] **N-X3** `ConjugateUpdates` expressible as limits of a GIG update.

## Reports

### Baseline requalified (N9) — 2026-10-27

**Qualified revision:** `3a04601`, the step 3.5 baseline (`2efa134`) plus 84 engine commits. The toolchain is step 3.5's: .NET SDK 10.0.201, `net10.0`, Release, on the same machine.

**Suites** (`scripts/fact-harness.ps1 -Project`, each case its own `dotnet test` process):

| Suite | Cases | Passed | Run |
|---|---|---|---|
| `Maths.Regression.Tests` | 406 | 406 | `20261027_020935664` |
| `Maths.Samplers.Tests` | 119 | 119 | `20261027_021721150` |
| `Maths.LinAlg.Tests` | 118 | 118 | `20261027_021924460` |
| `Maths.Distributions.Tests` | 93 | 93 | `20261027_022127437` |

**Downstream consumers:**
- `Maths.Geometry.Tests` passes 51 of 51 and `Synthetic.Tests` 38 of 38, under plain `dotnet test`.
- `Clustering.Graphical.SPC` builds.
- Earlier the same day the Pca migration also ran Mapper 28/28, Oracle 5/5, VizCore 240 (2 skipped by design) and DimReduction 83/83.

**The drift record.** Every one of the 84 engine commits since `2efa134` (in `regression`, `samplers`, `linalg`, `distributions` and `rng`) has a recorded account of what it changes:
- 76 are named in the preparation brief's or this segment's reports, each with its probe or its declared path change;
- 7 are this segment's N-series, recorded by commit range above, which add new code beside the engines;
- `87397fb` changed only a doc comment in `Pointwise.cs`; its report is in the preparation brief under its title.

**Seed sweep** (`SEED_SWEEP=50`):
- **Hold counts.** Regression's `20261027_023326165` matches the last sweep (`20261022_182007796`) claim by claim, except the adaptive-kernel comparison on the data axis, which went from 49 to 50 of 50. Samplers' `20261027_023344333` is identical to `20261021_195902630`.
- **Per-seed values.** To trace those, the sweep was rerun on the source of `6bc5134^` (`bc51934`, extracted with `git archive`, scratch) and split there:
  - *Oct 22 → `bc51934`:* every `HierarchicalPSpline` row changed, and nothing else. That is `d2db820`'s declared path change, the ensemble shell with overdispersed starts, whose agreement in distribution with the single chain is in its report.
  - *`bc51934` → HEAD:* 20–56% of the rows of the five BARS claims changed, and nothing else: adaptive kernel, EBIC, Gibbs loss, Student-t and contaminated-normal weights. That is `6bc5134`, which sheds knots from a start the posterior refuses. The later commits in that span change a doc comment, a PSIS readout and code outside the BARS sampler.
  - *Unchanged throughout:* BAPS, the anisotropic and tensor splines, the penalty and REML claims and the changepoint claim.
- **The 49 → 50.** It is seed 46 (data seed 372476911). On Oct 22 the too-tight fixed kernel reported a peak-location bulk ESS of 4,586, the signature of a chain held at a refused start. It now reports 9.7, against the adaptive kernel's 359.

**Unlike step 3.5,** this rerun did not rebuild a field-by-field output probe across revisions. Its drift evidence is each step's own recorded probe, plus the seed sweep's per-seed table split at `6bc5134^`.

**Gate status:** closed at `3a04601`. The first BASS integration tests measure against this baseline.

### The generic SMC core (N8) — 2026-10-27

`3a04601`, in `src/maths/samplers/smc/`. No existing code changed. The calibration segment's list, item by item:
- **Contracts.**
  - *`SmcStage<TState>`* is a target at a temperature, so a tempered run walks one target from b = 0 to 1 and a cut run walks different conditional targets.
  - *`ISmcParticles<TState>`* draws from the first stage, which must be normalized, and forks a particle for resampling's offspring.
  - *`ISmcKernel<TState>`* is a stage-invariant move. It gets a `Prepare` call with the population before the stage's moves, so it can tune from the past while staying fixed during the stage, as D10's policy asks.
  - Densities come through `ISamplingTarget`: moving to another temperature of the same target queries each particle's evaluation instead of repeating its work.
- **Weights** (`SmcWeights`): log-sum-exp normalization, ESS, the conditional ESS (Zhou, Johansen and Aston 2016), and the log-normalizer increment.
- **Resampling:** multinomial, stratified, systematic and residual. A zero-weight particle is never chosen.
- **Driver** (`SmcRun`, stage by stage): reweight, resample when ESS/N falls below a fraction, then mutate, each particle on its own `SeedTree` stream. `SmcSampler` runs a fixed sequence, or tempers adaptively by choosing each next b through a bisection on the conditional ESS.
- **Diagnostics, kept apart from chain R̂ and autocorrelation ESS:** per stage, the ESS, conditional ESS, whether the population resampled, distinct ancestors, surviving families (distinct initial lineages), acceptance, and the log-normalizer increment.
- **Checkpoints** resume the uninterrupted run bit for bit, because every stream is derived by position.
- **`RandomWalkKernel`:** Gaussian random-walk Metropolis for real vectors, scaled per coordinate from the population's spread.

**Tests** (`tests/samplers/SmcTests.cs`, 10):
- the weights at extreme scales;
- each resampling scheme's offspring averaging N·W_i over 20,000 repetitions, with systematic and residual keeping their deterministic parts;
- a conjugate Gaussian: the posterior moments and the closed-form evidence. Over 20 seeds (scratch) the log-evidence error had mean −0.0003 and sd 0.039;
- a bimodal likelihood, whose modes' masses (0.3 and 0.7) are recovered;
- without resampling or moves, the evidence telescoping to plain importance sampling exactly;
- results independent of parallelism;
- checkpoint resume.

**Not done:** BASS's fitting and cut-calibration adapters (stage 8), the cut SMC's bridges, OQ2's exact reference rung for parallel tempering (and the OQ1 rerun it enables), and SPC population annealing.

### N6, N7, and the SVD callers moved — 2026-10-27

Commits `09572c6`..`28fec55`.

- **`DiagonalPlusLowRank`** (N6): Σ = D + U·diag(s)·Uᵀ.
  - The factor of I + WᵀW (W = D^(−½)·U·diag(√s)) is taken as the thin SVD of W. That makes bᵀΣ⁻¹b and ln|Σ| sums of nonnegative terms, and nothing m × m is formed.
  - Solves and draws from N(0, Σ) and N(0, Σ⁻¹) each cost O(m·r).
  - Along a heavy direction (weight 10¹⁶) it keeps the quadratic form to 10⁻¹⁰ relative. The Woodbury subtraction loses all of it; the test checks both.
- **`KhatriRao`** (N7): the direct functional design X = A ⊙ H, through its factors.
  - It provides the Hadamard Gram, birth borders and response cross-products, for R's convention and centred for the Liang convention (the three-term form), so N7 does not wait on the intercept call.
  - It also accumulates over scattered (case, location) samples through a pairwise cascade.
  - *Terminology:* BASS's design is the column-wise Khatri–Rao product. Its rows are Hadamard products of factor rows, and its transpose is the face-splitting product. ABMRS's rows are face-splitting products themselves. The brief's "Room for extensions" paragraph now says so.
- **`ThinSvd` fixed.** The Grassmann exp map gave it two parallel columns of norms 1.8·10⁻¹⁶¹ and 0.115. ζ² overflowed, the rotation became the identity, and the Jacobi never converged. Three changes fix it:
  - rotation angles are now computed through the ratio of the norms, and norms and cosines are rescaled outside the normal range;
  - columns below ε·‖R‖_F are left alone, and their left vectors are completed;
  - non-finite input is refused.

  Regression tests cover that case and extreme column scales. `PairwiseSum`'s accumulators are now vector lanes, bit-identical to the scalar form.
- **Callers moved onto `ThinSvd`:**
  - *`Pca`* now reads the SVD of the centred data, at 8.0 s against 22.3 s on a random 6,238 × 617 matrix, with eigenvalues agreeing to 2·10⁻¹³. Wide inputs are slower: 0.53 s against 0.32 s at 300 × 2,000. Components are now signed so that their largest entry is positive.
  - *`StiefelManifold` and `GrassmannManifold`* use it in place of their private Gram-route SVDs. They keep their rule of zero left vectors where σ < 10⁻¹²: at the cut locus that direction's sign is arbitrary, and the corrupted-block tests of MoMPCA and SPRED depend on the two choices cancelling.
- **Suites:**
  - LinAlg 118/118, Geometry 49/49, Mapper 28/28, Oracle 5/5, VizCore 240 passed with its 2 by-design skips, DimReduction 83/83.
  - Before the fix, the manifold change broke 8 tests: 2 MoMPCA, 4 DistributedSpred and 2 SpredCylinder.
  - SPRED's ISOLET pilot reruns its warm start to 10⁻¹⁴. Its annealed objectives move at 10⁻⁵ relative, because the annealer's trajectory diverges from bit-level differences. The pilot JSONs that the suite rewrites in `issues/spred/pilot/` were restored to their committed versions.
- **Noticed:** Maths.Oracle failed 2 of 5 once at baseline and then passed on a rerun with no change in between.

### Decision-free stage 2 numerics (N1–N5, N10) — 2026-10-26

Three new LinAlg types, in `src/maths/linalg/`. No existing code changed, and no caller is migrated yet.

- **`PairwiseSum`** (N10): `Of` and `Dot` sum pairwise. The tree runs over blocks of 128 terms, and each block is summed in four interleaved accumulators. It uses plain double arithmetic in a fixed order, so the bits are the same on every machine. `MatrixOps.Dot` and `ColumnGramMatrix` keep their SIMD left-to-right sums for their existing callers.
- **`BorderedCholesky`** (N1–N4), an `ISpdFactor` that keeps A and diag(A⁻¹) beside L.
  - *Births* border the factor with a new last row, in `DenseCholesky`'s arithmetic. A factor grown one row at a time is bit-identical to the strict one.
  - *Deaths* drop a row of L and restore the triangle with Givens rotations.
  - *Changes* are a death and then a birth, staged together on a copy. **The new column goes last**, so the BASS design layer reorders its terms to match.
  - *Proposals* (`ProposeBirth`, `ProposeDeath`, `ProposeChange`) return a `FactorProposal` without touching the accepted factor. It carries positive definiteness, the relative pivot, the log-determinant change and the smallest leave-one-out pivot. `ProposedInverseQuadraticForm` gives bᵀA′⁻¹b for the candidate, and `Accept` commits it. Deaths and changes are staged on a second copy and swapped in, which costs an O(k²) copy per proposal.
  - *diag(A⁻¹)* gains w·wᵀ/d on a birth. On a death it loses v·vᵀ/v_m, and any entry that keeps less than 10⁻³ of itself is recomputed from the new factor, because that subtraction cancels when a near-duplicate dies. `LeaveOneOutPivot(j)` = 1/(A_jj·(A⁻¹)_jj) is D06's s_j.
  - *Rebuilds* (`Rebuild`, and automatically every `rebuildInterval` accepted updates) refactor the kept A strictly and recompute diag(A⁻¹). They report the drift removed as a `FactorDrift`, and keep the incremental factor if the strict one fails. The default interval of 100 is a placeholder, not a measured choice.
- **`ThinSvd`** (N5): Householder QR, then one-sided Jacobi rotations on R, without forming AᵀA. U is completed where σ is exactly zero. The rank tolerance is max(m, n)·ε·σ₁, the same as the O-O4 oracle's.

**Tests** (`tests/linalg/`, 33 new): `Maths.LinAlg.Tests` 98 of 98 under the fact harness (`20261026_222717257`).
- *`PairwiseSumTests`*: errors are measured against exact BigInteger sums and stay within (32 + ⌈log₂(N/128)⌉ + 3)·u·Σ|terms|, including on F-O7's near-constant column.
- *`BorderedCholeskyTests`*:
  - births are bit-identical to `DenseCholesky`, condition estimate included;
  - over 400 random births, deaths and changes with no rebuild, the factor tracks a fresh factorization to 10⁻¹²·√A_ii, diag(A⁻¹) to 10⁻¹¹ relative, and the log-determinant, the proposals' predictions and bᵀA⁻¹b likewise;
  - the final rebuild finds drift ≤ 10⁻¹²;
  - rejected proposals leave every bit of the accepted factor unchanged;
  - leave-one-out pivots match the squared sine against the others, computed through `ThinSvd`, to 10⁻⁸;
  - D06's order counterexample: the newest pivot moves between 10⁻⁶ and 10⁻¹² with insertion order, while the leave-one-out vector does not;
  - after a near-duplicate's death, diag(A⁻¹) holds to 10⁻¹³. With the cancellation guard disabled, this test fails.
- *`ThinSvdTests`*:
  - reconstruction and orthonormality, for tall, wide and square matrices;
  - a constructed spectrum from 1 to 10⁻¹¹, recovered to 10·m·ε absolute;
  - rank 7 from a known product;
  - centred curves lose one rank;
  - zero matrices, and a zero column completed in the middle.

**Not done here:**
- N6 (low-rank solves), N7 (face-splitting contraction, which waits on the intercept call), N8 and N9;
- moving `Pca` and the geometry manifolds' private Gram-route SVDs (`StiefelManifold`, `GrassmannManifold`) onto `ThinSvd`;
- benchmarks.

### Cross-product summation — 2026-10-26

New item N10, from D06's fixture F-O7. On a design with a near-constant column, left-to-right summation of the Gram gave 7·10⁻⁴ to 2·10⁻³ relative error in s ≈ 10⁻⁹. Pairwise summation gave 2·10⁻⁷ to 3·10⁻⁷, which is ε·κ. Nothing in `src/maths/linalg` sums with compensation yet.

### Adaptation policy drafted — 2026-10-25

D10's note written as a proposal (`53a4876`, carried into this segment at the split).

