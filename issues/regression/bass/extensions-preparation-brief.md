# BASS's extension seams (chip brief)

Date: 2026-10-31. The port's segments each listed "Room for extensions": seams to put in place so that the [extensions](extensions-brief.md) plug in without reworking the port. This chip builds the ones still open, as contracts with today's members only.
- *Its rule:* every seam is a refactor that changes no draw. BASS, BARS and every engine sharing the touched code stay bit-identical.
- *What it builds* (the checklists' ids):
  - F-X3, term types behind one design-kernel contract, and O-X4, term derivatives available to sensitivity through it;
  - F-X1, the observation-model seam with an offset;
  - F-X2, the sweep's latent slot;
  - F-X6, the declared output-side inner product;
  - N-X1, the face-splitting contraction open to tensor designs;
  - N-X2, `ISpdFactor` admitting a wide-band factor;
  - N-X3, `ConjugateUpdates` expressible as limits of a GIG update.
- *Not in scope:*
  - F-X4, the proposal-kernel contract, which was decided on 2026-10-29 to be built with KHAOS's Bernoulli-inclusion proposal, not before it;
  - any extension itself (no generalized-hyperbolic noise, no tensor B-splines, no GIG sampler, no correlated residual);
  - any default change;
  - `docs/Maths` and the checklists, which the main session updates after verifying.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md), whose design principles govern every contract here. Then:
- the "Room for extensions" sections of [port-fitting-brief.md](port-fitting-brief.md), [port-numerics-brief.md](port-numerics-brief.md) and [port-outputs-brief.md](port-outputs-brief.md), which state each seam;
- the extensions each seam serves, in [extensions-brief.md](extensions-brief.md):
  - generalized-hyperbolic noise (GBMARS) for F-X1, F-X2 and N-X3;
  - free-knot tensor-product surfaces (ABMRS) for F-X3, N-X1 and N-X2;
  - active-subspace analysis for O-X4;
  - residuals correlated along the output index for F-X6;
- D04 in the fitting brief (the sweep's order, the robust order's latent slots);
- the derivations review §7, item 6, which the fitting brief cites for the term-type contract.

## The seams, as this chip builds them

Each seam gets the smallest contract that its extension needs and today's code fills. Where a shape is a real choice, make it, name it in the report with its alternative, and keep it easy to revise; the user reviews the shapes.

1. **F-X3, the design-kernel contract.**
   - The term types (hinge products, categorical subsets, the functional variable) sit behind one contract. Each supplies its column evaluation, its support count, its moments and its derivatives.
   - Tensor-product B-splines and directional hinges will join later as members. Write the contract so they can without changing it, and say how each would fill it.
   - *O-X4:* term derivatives reach sensitivity through it. Check each term type's derivative against central differences.
2. **F-X1, the observation-model seam.** BASS's observation model hands the Gaussian core a working response and weights, (y − offset, W), not weights alone. Gaussian BASS passes a zero offset and unit weights, so its marginal and draws are unchanged.
3. **F-X2, the sweep's latent slot.** D04's robust order (structure, τ, a, latent scales, then σ² and the skew) becomes the sweep's general form, with the latent slots empty for Gaussian BASS. An empty slot draws nothing and consumes no random number.
4. **F-X6, the output-side inner product.** The direct contractions are written G = (AᵀA)∘(HᵀQH) and h_j = A_jᵀYQH_j, with Q = I declared. A later Q = R⁻¹ then changes only the output-side factors, plus Σᵢ yᵢᵀQyᵢ and a log|R| term per case in the collapsed marginal. Q = I takes today's arithmetic path exactly, so that a declared identity costs nothing.
5. **N-X1, the face-splitting contraction.** `KhatriRao` covers tensor designs: a row that is itself a face-splitting (row-wise Kronecker) product, as ABMRS's rows are. Check it against the explicit product on small tensor designs.
6. **N-X2, a wide-band factor.** `ISpdFactor` admits a wide-band Cholesky beside `DenseCholesky`, chosen by size.
   - Build the factor and its selection rule, and check it against the dense factor on block-banded Grams: its solves, its log-determinant, its draws, and its status at near-singularity.
   - BASS keeps the dense path. Nothing selects the wide-band factor until a caller declares a banded Gram.
7. **N-X3, conjugate updates as GIG limits.** `ConjugateUpdates`' Gamma and inverse-gamma updates are parameterized as the limits of a GIG update, so GBMARS's latent-scale block generalizes them. No GIG sampler is built.
   - The work is the parameterization and its documentation in the code.
   - Add a test that the GIG density's limits match each update's density.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Other chips work beside you.**
  - The flat-histogram chip owns `src/maths/samplers/flathistogram/`, `BassFlatHistogramRun.cs`, `BassMulticanonical.cs` and the flat-histogram tests.
  - The F-R1 chip owns `tests/regression/BassDeterministicProfileTests.cs`.
  - Change no signature their files call. If a seam needs one changed, finish that seam after their reports, and name it in yours. If a build fails mid-change in a file you do not own, wait and retry; never edit it.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-ext\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-ext\` when you finish.
- **Refactor by writing the new beside the old,** moving callers over, then deleting the superseded code in a final cleanup step: no shims, no `[Obsolete]`, no compatibility layers (the repo is pre-release).
- **Bit-identity is the gate for every step:**
  - `BassBitIdentityTests`' hashes hold;
  - BARS's hash holds;
  - the Regression, Samplers, LinAlg and BASS oracle suites pass.
  - Where a seam touches shared code (`KhatriRao`, `ISpdFactor`, `ConjugateUpdates`), every engine using it is bit-identical. Find them by search, and say which you checked.
- **Commit your own files only, by explicit path,** one commit per seam, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's untracked notes under `issues/regression/` and other sessions' edits stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.

## The steps

In this order: F-X3 with O-X4, then F-X1, F-X2, F-X6, N-X1, N-X2 and N-X3, then the cleanup. Each lands with its checks and ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next. Name the passes and their fixes in the report.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - no draw changed anywhere, with the hashes named;
  - each contract carries what its extension needs, read against the extensions brief;
  - nothing built ahead of its extension beyond the contract and today's members.
- *Coverage:*
  - each term type through the contract, including a functional fit;
  - the derivatives at a hinge's knot, where the derivative jumps;
  - Q = I on a direct fit;
  - the wide-band factor at the band's edge and near singularity.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- each seam's contract, its shape calls with their alternatives, and how its extension will fill it;
- the checks with their tolerances, and the bit-identity evidence;
- the suites with their counts;
- the review passes and their fixes;
- anything the main session must know, any seam deferred for another chip's files included.

## Reports

Newest first.

### The extension seams, built — 2026-10-31

All seven seams landed as contracts with today's members only, one commit each, and no draw changed. F-X4 and the extensions themselves were not touched. Nothing is deferred for another chip's files: no signature the flat-histogram or F-R1 chips call was changed (`BassModel.Create(config, data[, form])`, `BassChain`'s constructor and `Restore`, `BassStructureMoves` are as they were; the new `Create` overload and the optional parameters are additions).

**Commits** (on `main`, in order):

| Seam | Commit |
|---|---|
| F-X3 with O-X4 | `9b0c189` |
| F-X1 | `ebe2308` |
| F-X2 | `27b86b9` |
| F-X6 | `0131c78` |
| N-X1 | `df00275` |
| N-X2 | `aa90524` |
| N-X3 | `83d3d09` |
| Review fix (the ladder test) | `cbacbdb` |
| Cleanup (descriptions) | `d62f596` |

#### The contracts, their shape calls, and how each extension fills them

**F-X3, the design-kernel contract** ([DesignKernel.cs](../../../src/maths/regression/spline/bass/DesignKernel.cs), [DesignKernels.cs](../../../src/maths/regression/spline/bass/DesignKernels.cs)).
- *Shape.* `IDesignKernel<TPoints>` is one term type's part of a term, per side of the Khatri–Rao column (`CaseInputs` for the case side, `OutputInputs` for the output side). It supplies `Reads`, `Validate`, `Multiply` (the column evaluation, applied in order to a column of ones), `MultiplyDerivative` (one variable, for a `Differentiable` part), `Scale`, and `Factors`, the separable `DesignFactor`s for the moment service. Today's members: the continuous hinge product, the categorical subsets, the functional hinge product, held by `DesignKernels`. `BassBasis` composes columns from them and `BassMoments` composes factors from them; every caller of `BassBasis` is unchanged.
- *Calls and alternatives.*
  1. *The contract is per part, generic in the points type; `BassTerm` stays closed.* Alternatives: one non-generic interface with a side enum and both inputs; or `BassTerm` holding a list of polymorphic parts. The second touches every consumer of `Continuous`, `Categorical` and `Functional`, so a new part kind adds a typed list to `BassTerm` (the extension's change), not a member to the contract.
  2. *Support is not a kernel member.* The admission gate reads the support of the product (hinge × indicator), which only the composed column has, so `BassBasis.Support` stays a function of that column. Alternative: a per-kernel count with intersection semantics, which needs the joint column anyway.
  3. *A derivative at a kink is the mean of the one-sided limits*: s/2 for a degree-one hinge at its knot, 0 from degree two up. That is what a central difference returns, and the knots are the training inputs, so every case sits on one. Alternatives: the right derivative, or the Clarke interval.
  4. *A variable belongs to one part*, so a derivative is the differentiable part's derivative times the other parts' values, with no general product rule. A categorical part is not differentiable and refuses (`NotSupportedException`); the alternative, zero, would claim a derivative in a level index.
  5. *Moments are the separable factors* (`HingeDesignFactor`, `SubsetDesignFactor`, value records, so `BassMoments` caches a factor's mean and a pair's moments by the factor itself). Alternative: kernel-level `Mean(term)` and `Cross(a, b)`, which would move the ANOVA bookkeeping into every kernel.
  6. *`DesignKernels` is a fixed set,* with no registration hook; a hook is the extension's.
  7. *O-X4 reaches sensitivity as* `BassMoments.FactorDerivative` (beside `FactorValue`) and `BassBasis.CaseDerivative` / `OutputDerivative`. Alternative: a service producing ∇f per draw, which the active-subspace work writes from these.
- *How the extensions fill it.* A tensor-product B-spline part (ABMRS) is a case-side kernel whose `Multiply` writes the product of one-variable B-splines, whose derivative differentiates one factor, and whose factors are one-variable piecewise polynomials that `PiecewisePolynomial` integrates as it does a hinge; `BassTerm` gains its list. A directional hinge is a case-side kernel for (wᵀx − t)₊, differentiable in every variable of w. It is not a product of one-variable factors, so it declares none and supplies its moments through a further interface that the moment service tests for. A test implements `IDesignKernel<CaseInputs>` for it, derivative against central differences, to show that the signatures carry it.

**F-X1, the observation-model seam** ([BassObservationModel.cs](../../../src/maths/regression/spline/bass/BassObservationModel.cs)).
- *Shape.* `IBassObservationModel.Working(response)` returns a `WorkingResponse`: `Values` = y − offset (the response itself, uncopied, when there is no offset), `Weights`, and `Output` (F-X6). `GaussianObservation` is the one member: zero offset, unit weights. `BassModel` carries the model and its working response; every design fits it, and `BassResolved.Resolve` reads its level and spread, so an offset is the shifted data in the resolved record, the target and the design, bit for bit.
- *Calls and alternatives.*
  1. *The core refuses weights* (`NotSupportedException`, naming F-X1) rather than forming a weighted Gram, which is GBMARS's. Alternative: build the weighted Gram, centring and marginal now, which is the extension.
  2. *The observation model is a runtime argument:* a new `BassModel.Create(config, data, observation, form)` overload. Alternative: a `BassConfig` field, which would change the serialized config and the run fingerprint.
  3. *`BassDesign` takes the `WorkingResponse` as an optional parameter,* so its callers are unchanged.
- *How GBMARS fills it.* Its observation model returns y − β√w·v and weights 1/(w·v) at the current latent state, and the core adds the weighted Gram. *For the main session:* the readers of the raw response, the pointwise densities, the training fit, the residual-dependence check and prediction, still read y, not y − offset. Putting the offset back is the extension's.

**F-X2, the sweep's latent slots** (same file; [BassSweep.cs](../../../src/maths/regression/spline/bass/BassSweep.cs), [BassLadder.cs](../../../src/maths/regression/spline/bass/BassLadder.cs), [BassSmc.cs](../../../src/maths/regression/spline/bass/BassSmc.cs)).
- *Shape.* `BassSweepSlot` is D04's order (`Structure`, `Tau`, `Coefficients`, `Latent`, `Scale`), written once in `BassSweepSlots`. The chain, the tempering replicas and the SMC kernel all assemble their blocks from it: BASS's move step and τ block, then the observation model's `Blocks(slot)` for the later three. Gaussian BASS returns none, so its scan has the two blocks it had, and an empty slot draws nothing and consumes no random number.
- *Calls and alternatives.* Blocks are the existing `IBlock<ReversibleJumpChain<BassState>>`, so the latent state lives with the observation model, not in `BassState`; alternative: `BassState` gains a latent field, which changes every snapshot and so waits for an extension that needs it. The observation model fills only the three later slots.
- *How GBMARS fills it.* A block for a, one for the latent scales, one for σ² and the skew, in `Blocks(slot)`.

**F-X6, the output-side inner product** ([InnerProduct.cs](../../../src/maths/linalg/InnerProduct.cs), [KhatriRao.cs](../../../src/maths/linalg/KhatriRao.cs)).
- *Shape.* `IInnerProduct` is Q on the output index: `Dot` over stacked vectors (k cases of length m end to end), so one member gives the Gram entries hᵀQh and the energy Σᵢ yᵢᵀQyᵢ; `LogDeterminant` for the collapsed marginal's log|R| term (n·ln|R| = −n·ln|Q|); `IsIdentity`; `Dimension`. `IdentityInnerProduct.Dot` is `PairwiseSum.Dot`. `KhatriRao.Gram`, `Border` and `CrossProducts` take an optional product (null is Q = I), written G = (AᵀA)∘(HᵀQH) and h_j = a_jᵀ·Y·Q·h_j. `WorkingResponse.Output` carries it, and `BassDesign` routes ‖ỹ‖², the residual's sum of squares, the output-side Gram entries of a border and the cross product through it.
- *Calls and alternatives.*
  1. *The centred forms stay Q = I.* Under Q the intercept's orthogonal complement is the Q-centred one, and the three-term form carries over with m replaced by 1ᵀQ1, the means by Q-weighted means and the inner products by hᵀQh (documented on `KhatriRao`); that centring, the level and the σ² prior's scale are the extension's, so `BassDesign` refuses a Q other than the identity, with the reason. Alternative: write the Q-centred forms now, untested by any user.
  2. *Q rides the `WorkingResponse`,* with the offset and weights, as the noise's part of what the core fits. Alternatives: a `BassConfig` field (serialized form), or a `BassData` member.
  3. *`Dot` alone,* with no `Apply` (Qh): a tridiagonal Q keeps each product linear in m, and `Apply` would only save recomputation.
- *How correlated residuals fill it.* A tridiagonal `IInnerProduct` for the Markov residual, plus the Q-centred forms.

**N-X1, the face-splitting contraction** ([KhatriRao.cs](../../../src/maths/linalg/KhatriRao.cs)).
- *Shape.* `FaceSplittingGram(axes)` and `FaceSplittingCrossProducts(axes, y)` for X = B₁ ∗ … ∗ B_d (row i = the row-wise Kronecker product of the axes' rows, columns lexicographic, axis 0 slowest), accumulated over each row's nonzeros in the block cascade the scattered Khatri–Rao Gram already used; `FaceSplittingHalfBandwidth(axes)` for the factor's selection; `KroneckerGram(axes)` for a full grid. The scattered Gram's loop moved into the shared `AccumulateGram` and is bit-identical (a pinned hash).
- *Calls and alternatives.* Dense output, bounded by 46,340 columns like `DenseCholesky`; alternative, LAPACK band output feeding a `DecomposeBanded`, which ABMRS wants at real sizes and which nothing needs yet. Gridded cross products are not built.

**N-X2, the wide-band factor** ([WideBandCholesky.cs](../../../src/maths/linalg/WideBandCholesky.cs)).
- *Shape.* `WideBandCholesky : ISpdFactor`: strict Cholesky of a Gram zero beyond a declared half-bandwidth, row-packed so that every inner product reads two contiguous runs; O(n·p²) to factor, O(n·p) per solve. `SpdFactors.UsesWideBand(order, halfBandwidth)` is the rule (order ≥ 64 and 2(p + 1) ≤ n) and `SpdFactors.ForBand(gram, halfBandwidth)` applies it. BASS keeps the dense path; nothing selects the band.
- *Calls and alternatives.*
  1. *A new class beside `BandCholesky`,* which already exists and handles a band, because ABMRS's Gram needs the dense factor's conventions: strict, no floor, and the draw L⁻ᵀz, so one z gives one draw under either factor. `BandCholesky` draws with LDLᵀ's D^(−½) and floors pivots. On a Gram that is zero beyond the band the two factors do the same operations on the same nonzero terms in the same order, so everything is bit-identical to the dense factor's. Alternative: `BandCholesky(strict, Ldlt)`.
  2. *The rule is the operation count's,* not a measured crossover, and it errs toward the dense factor: timed, the band is already 1.5× faster at 2(p + 1) = n and 16× faster at p + 1 = n/8.
  3. *The band is the caller's declaration,* a precondition, as `BandCholesky`'s is.

**N-X3, conjugate updates as GIG limits** ([GigParameters.cs](../../../src/maths/samplers/blocks/GigParameters.cs), [ConjugateUpdates.cs](../../../src/maths/samplers/blocks/ConjugateUpdates.cs)).
- *Shape.* `GigParameters(P, A, B)` is GIG(p, a, b) ∝ x^(p−1)·exp(−(a·x + b/x)/2), a parameterization with `OfGamma`, `OfInverseGamma`, `ToGamma`, `ToInverseGamma`, `IsGamma`, `IsInverseGamma` and `LogKernel`. `ConjugateUpdates.NormalVariance` and `NormalPrecision` gain GIG forms, GIG(p − n/2, a, b + S) and GIG(p + n/2, a + S, b), of which the inverse-gamma and gamma updates are the limits a = 0 and b = 0. No GIG sampler, distribution or Bessel K.
- *Calls and alternatives.* The two existing updates keep their formulas, tested equal to the GIG forms' limits bit for bit; alternative, routing them through the GIG form, which gives the same bits but puts a conversion in HPS's and BAPS's hot paths.

#### Checks

Every tolerance below is the test's own; "exact" is `Assert.Equal` on doubles.
- *F-X3, O-X4.* The factors with the kernels' scale reproduce the column to 1e-13·(1 + |x|) on the scalar fixture's seven terms (hinge, categorical, mixed, a hinge facing its end, where c = 0 and R leaves it at 1) at degrees 1–3, to the same on the functional fixture's terms on both sides, and to 1e-12 on every term of a direct functional fit. Derivatives against central differences: interior points, h = 1e-6, 1e-6·(1 + |d|), degrees 1–3, scalar terms; the training inputs, which are the knots, h = 1e-7, 5e-5·(1 + |d|), where a degree-one derivative jumps; the output side the same on the functional terms and on every term with an output part in the functional fit; `FactorDerivative`, h = 1e-7, 1e-5·(1 + |d|), at the knot and 0.03 either side. A categorical part refuses a derivative; a term without the variable is constant in it. A directional hinge, a kernel reading two variables, implements the contract and matches central differences to six digits.
- *F-X1.* A design on y − offset equals a design on the shifted data, exact: the births' proposals, ‖ỹ‖², the residual and fitted sums of squares, the coefficients and the residual, scalar and direct functional. The resolved record's level, mean square and σ² prior are exact. Two 150-step chains, one on the shifted data and one on the offset, agree state by state (M, ln τ, both sums of squares), exact. Weights are refused; a working response of another response is refused.
- *F-X2.* The order is D04's; Gaussian BASS's sweep is [move step, τ block], one block with τ held; a model that fills the later slots runs them after τ in slot order, on the scan's stream (the first sweep agrees with Gaussian's, later ones diverge once a block draws), and a ladder of three rungs runs them at every rung. An observation model that fills every slot with nothing walks Gaussian's chain, scan stream, move stream and block counters word for word over 100 steps.
- *F-X6.* Q = I: `Gram`, `Border` and `CrossProducts` with the identity equal those without, exact; the identity's `Dot` is `PairwiseSum.Dot`; a direct design with a declared identity equals the default one through four births and a rebuild, exact. A tridiagonal Markov Q: G to 1e-12, the border to 1e-12, h_j to 1e-11 (relative) against the materialized design and I ⊗ Q; and −½Σᵢ yᵢᵀQyᵢ + (n/2)·ln|Q| equals the Gaussian log density with the dense (nm)² covariance I ⊗ R, to 1e-10. A dimension mismatch is refused; the core refuses Q ≠ I with the reason, a scalar response refuses a Q, and the working response refuses a mismatched one.
- *N-X1.* `FaceSplittingGram` against XᵀX of the materialized row-wise Kronecker product, 1e-13 relative, two axes (7 × 6 hats, 700 scattered samples) and three (4 × 3 × 5, 300); Xᵀy to 1e-12; the half-bandwidth of an 8 × 7 hat design is p·ν₂ + p = 8 at p = 1, no Gram entry lies beyond it and some lie on it; `KroneckerGram` on a grid against the materialized product, 1e-13, and equal to the face-splitting Gram over the same grid. The scattered Khatri–Rao Gram is unchanged bit for bit: a hash of its entries (2,123 samples, zeros in the factors) was taken with the old loop, `6e4c55f9a3d70b21`, and pinned.
- *N-X2.* On the 72 × 72 block-banded Gram of a 9 × 8 tensor design (half-bandwidth 9), everything equals the dense factor's, exact: the factor's entries, relative pivots, `LogDet`, `Status`, `Solve`, `InverseQuadraticForm`, `SampleInnovation` from one z, `EquilibratedCondition`. L·Lᵀ reconstructs the Gram to 1e-12, and A·x = b to 1e-10. At the band's edge: the declared bandwidth is the Gram's (an entry on it is nonzero, none beyond), a declaration above it by 5, or by 1000 (clamped to the full band), changes nothing, a diagonal is a band of 0, and an entry placed beyond the band is no part of what the factor sees. Near singularity, a pentadiagonal order-80 matrix with a pair of rows 1 − 10⁻⁸ alike reports the dense relative pivot (2·10⁻⁸) and a pivot ratio below 10⁻⁷, all else exact; an exactly dependent pair and an indefinite one stop at row 41 with the dense factor's `FailedAt` and `Status`, `LogDet` NaN, and every solve, draw, condition estimate and `Factor` refused. The rule's boundaries are tested and what `ForBand` returns. Timed at n = 1600, p = 80, 6.6 s dense against 0.13 s band, under load.
- *N-X3.* The GIG update is Bayes' rule for a prior with a, b > 0, to 1e-12, for the variance and the precision forms. Over 200 random updates the limits give the gamma and inverse-gamma updates' shapes, rates and scales, exact. The kernel at b = 0, or a = 0, differs from the Gamma, or inverse-gamma, log density by a constant, to 1e-12. The normalized GIG density, its constant by the integral form of Bessel K (checked against K_{1/2} and K_{3/2} in closed form, 9 digits), converges to each update's density with error ≤ 4·b (or 4·a·3) at b, a = 10⁻², 10⁻⁴, 10⁻⁶, 10⁻⁸, monotonically, below 10⁻⁷ at the last.

#### The suites and bit-identity

*Bit-identity.* `BassBitIdentityTests` (eight: six BASS fits, `ee7d60a2738b56d8`, `414e30ea867f9058`, `137dad3cbeb33acb`, `72fda364ee08f1f6`, `3e5f43ed894f60bc`, `e57e8cbe0c39d339`; Friedman tempered, `0d067c924599ed93`; the tempered BARS fit, `7edb59818f141da9`) passed at the baseline, in the Regression run of each seam's snapshot, in the cumulative run, and on the final committed tree. Engines sharing the touched code, found by search:
- `KhatriRao`: only `BassDesign` (direct fits, the three functional hashes) and its tests; its scattered Gram is pinned.
- `ISpdFactor`: unchanged; `WideBandCholesky` is new and nothing selects it.
- `ConjugateUpdates`: BARS (`StudentTWeights`), BAPS (`Baps`), the hierarchical P-spline, BGP's sampler and the collapsed g-prior's σ² posterior; the existing four-argument formulas are unchanged and only overloads were added, so the BARS hash and the BAPS, HPS and BGP suites stand. One cref in `HierarchicalPSpline.cs` became ambiguous with the overload and was made explicit.
- `BassBasis` and `BassMoments`: the posterior surface, pointwise densities, prediction, calibration emulators and sensitivity, all through the same `CaseColumn` and `OutputColumn`; the sensitivity suites (moments against tensor quadrature, Sobol, effects, basis sensitivity) pass.

*Suites,* from the baseline at `f14e8ea` to the final tree. Each seam's snapshot ran Regression and Oracle (N-X1 and the cumulative run also LinAlg and Samplers), from a copy of its build:
| Suite | Baseline | After each seam | Final (cumulative) |
|---|---|---|---|
| LinAlg | 120 pass | 123 (F-X6), 129 (N-X1) | 139 pass (+19: inner product 3, face-splitting 6, wide-band 10) |
| Samplers | 265 pass, 2 fail | 271 (N-X1) | 277 pass (+6 GIG limits) |
| Regression | 1,342 pass, 3 fail | 1,366 of 1,367 (F-X3), 1,375 (F-X1), 1,380 (F-X2), 1,381 (F-X6), 1,383 (N-X1) | 1,383 pass |
| Oracle | 101 pass, 7 fail | 106 of 108 (F-X3, F-X6, N-X1); 102 (F-X1) and 103 (F-X2) | 106 pass, 2 fail |

Every failure is not mine. The baseline's Samplers failures were the flat-histogram chip's two `WhamTests` warm-start tests and its Regression failures three `BassFlatHistogramTests` (one, a different seed, was still failing in the F-X3 snapshot); both cleared as that chip's work landed. The Oracle suite's failures are all `TimeoutException`s of the R oracles at 120 s (7 of 7 at the baseline, 6 of 6 at F-X1 and 5 of 5 at F-X2, when five snapshot runs shared the machine); the two that remain in the cumulative run are `DiagramMetricsTransportParityTests`, TDA's transport oracle, which fails the same way at the baseline. The BASS oracle tests (the chain, tempering, SMC and fit parity) pass in the F-X3, F-X6, N-X1 and cumulative runs.

The cumulative run is of the tree before the ladder test (`cbacbdb`) and the doc-only edits that followed; the final committed tree was rebuilt and ran the eight hashes and all the seams' tests.

#### The review passes

Each seam ended with a read of its own code and tests for correctness and coverage; the cumulative run closed them. What the passes found:
- *F-X3.* The first run of the knot test failed on its own premise: it required a degree-one knot under a case at every degree, and was made to require a knot under a case. No code change.
- *F-X1.* The first cut read the offset in the design but left `BassResolved` on the raw response, so the zero centre's level and the scaled σ² prior would have disagreed with the design's under an offset. `Resolve` now takes the working response, and the test compares the resolved record, the design and a chain.
- *F-X2.* A stub that implemented only `Working` stopped compiling when `Blocks` joined the interface (caught by the build). The ladder's and the SMC kernel's reading of the slots was covered only through the shared assembly and the Gaussian hashes; a test now runs a filled model's slots at every rung of a ladder (`cbacbdb`). The SMC kernel's reading stays covered by its Gaussian parity tests and the shared assembly.
- *F-X6.* Giving the centred forms a Q parameter was considered and dropped, since only the identity could use it. They stay at Q = I, the design refuses Q ≠ I with the reason, and the Q-centred form is recorded on `KhatriRao`.
- *N-X1.* Moving the scattered Gram's loop into the shared accumulation was gated by a hash taken with the old loop, before the change.
- *N-X2.* The first near-singularity test asserted a condition estimate above 10⁷ for a 2·10⁻⁸ pivot; Hager's estimate found a smaller one, which `AssertSameFactor` shows is the dense factor's too. The assertion became the pivot ratio.
- *N-X3.* The new overloads made one cref in `HierarchicalPSpline.cs` ambiguous (a compiler warning); it names the four-argument form now.
- *Cleanup (`d62f596`).* The superseded code went with each seam, because what each seam replaced was private to its class or file: the basis's hinge loop and delegate, the moment service's hinge caches and `TermFactor`, the chain's own block assembly, the ladder's and SMC kernel's `BassTauBlock` construction, the scattered Gram's inline accumulation. So the final sweep found nothing standing beside a replacement; it corrected the two descriptions that still said "moves, then τ". This departs from "write beside, delete in a final sweep" in letter: no caller of a replaced member exists outside the file that replaced it, and the hashes were the A/B comparison.

#### For the main session

- *Nothing is deferred.* No other chip's file or signature was touched. The working tree still holds the user's notes (`opus-bayesian-evals.md`, `perplexity-bart-etc.md` and two untracked `perplexity-*` notes) and the other chips' edits, unstaged.
- *For the checklists and `docs/Maths`:* F-X1, F-X2, F-X3, F-X6, O-X4, N-X1, N-X2 and N-X3 can be ticked. F-X4 stays open.
- *What each extension still has to write,* beyond filling the contract:
  - GBMARS: the weighted Gram, centring and marginal (the core refuses weights), the latent state and its blocks, putting the offset back in prediction and in the pointwise densities, the GIG sampler and Bessel K, the decision on `BassState`.
  - ABMRS: a `BassTerm` slot or its own term type for the tensor part, the knot moves, band output from the face-splitting accumulation with a `DecomposeBanded` on the wide-band factor, and the factor's use by its target.
  - Correlated residuals: the Q-centred forms, the Q-aware level and σ² prior scale, the tridiagonal `IInnerProduct` and its θ updates.
  - Active subspaces: the gradient integrals from `FactorDerivative`, and the moments of a non-separable kernel.
- *The shape calls most worth the user's review:* F-X3's per-part generic contract with `BassTerm` closed; the kink convention (the mean of the one-sided limits); F-X1's refusal of weights and the observation model as a runtime argument; F-X6's Q riding the working response, with the centred forms left at Q = I; N-X2's new class beside `BandCholesky`.
- *The wide-band factor is faster than the rule needs:* 1.5× at 2(p + 1) = n. The rule could move toward the band; it was left on the operation count.
- *Test artefacts:* `artifacts\chip-ext\` is deleted.
