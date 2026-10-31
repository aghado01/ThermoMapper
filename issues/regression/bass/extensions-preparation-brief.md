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
