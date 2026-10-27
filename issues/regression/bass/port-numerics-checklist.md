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
- [ ] **N6** Diagonal-plus-low-rank solves and log-determinants, without m × m inverses.
- [ ] **N7** Face-splitting design contraction with the Cartesian (Hadamard) shortcut.
- [ ] **N8** Generic SMC core in `samplers/smc` (optional; specified in the calibration segment).
- [ ] **N9** Qualified baseline rerun after the shared changes, before the first BASS integration tests.
- [x] **N10** Cross-products summed pairwise or with compensation: `PairwiseSum`. F-O7's near-constant design is the acceptance case: its left-to-right error is 10⁻³ relative in s.

### Room for extensions

- [ ] **N-X1** The face-splitting contraction covers tensor designs (ABMRS).
- [ ] **N-X2** `ISpdFactor` admits a wide-band factor beside `DenseCholesky`.
- [ ] **N-X3** `ConjugateUpdates` expressible as limits of a GIG update.

## Reports

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

