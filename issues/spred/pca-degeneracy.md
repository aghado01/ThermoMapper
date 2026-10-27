I'll build N6 from the thin SVD of D^(−½)·U·S^(½) so its quadratic forms and log-determinant stay as sums of nonnegative terms, avoiding the cancellation issue with Woodbury. For N7 I'll build both intercept conventions since neither depends on your choice, and before touching `Pca` (which has nine downstream consumers), I'll optimize the SVD's inner loops and benchmark the results.

Vector lanes give the same bits as the scalar accumulators, and all 36 tests pass. Now N6.

Moving to N7: I've identified the design as the column-wise Khatri–Rao product A ⊙ H (column j = a_j ⊗ h_j), whose transpose is the face-splitting product referenced in the brief's title.

Both build. Now their tests: N6 against the dense covariance, N7 against the design materialized explicitly.

All 14 tests pass, confirming the Woodbury subtraction loses every digit along the heavy direction. Next I'll update the callers, starting with timing `Pca` against a `ThinSvd` route at ISOLET size, since VizCore and SPRED run it on 6,238 × 617 data.

ThinSvd is much faster than Pca at ISOLET scale (8.0s vs 22.3s) with eigenvalues matching to 2·10⁻¹³, and it's only slower for wide inputs (300×2,000). Since this is a net win, I'll baseline the affected suites and migrate Pca to ThinSvd.

While the baseline runs, I'll leave the sources alone so its build sees the unmigrated code.

The baseline run of the four affected suites is still going in the background. The migration script is ready, and I'll apply it once the baseline reports.

Baseline: Geometry 49/49, DimReduction 83/83 (18 minutes), Mapper 28/28, and Maths.Oracle already fails 2 of 5 before any change. Before migrating I need to know which two fail.

This time all five oracle tests pass, so the two failures at baseline were intermittent. The relevant two, PCA parity and the Grassmann median, pass either way. Applying the migration.

Build succeeds, with the only warning being a pre-existing one in `Spherical.cs`. I'll now run the full test suite in the background—LinAlg, Geometry, Mapper, Oracle, VizCore, and DimReduction—which should take about 30 minutes.

While that runs, I'll update the numerics docs for N6 and N7. The migration's report waits for its results.

The baseline run modified four tracked files under `issues/spred/pilot/`. I'll back up these baseline versions to compare against the migration's changes later, then restore them afterward.

The baseline run reproduced the committed objective values exactly; only the date and timing fields changed. After the migrated run I'll compare objectives against these copies. Waiting for the after-run.

The migration broke eight tests—two `MoMPCATests` in Geometry and six `DistributedSpredTests`/`SpredCylinderTests` in DimReduction, all related to corrupted blocks or the Grassmann median—that passed at baseline. I'll read the failure messages to determine whether the migration or the tests themselves are at fault.

Found two issues: one test's Jacobi solver fails to converge within the Grassmann exponential map, and another now recovers a subspace 0.37 radians off (versus the previous 0.1 tolerance). I'll start by identifying the input causing the Jacobi solver to stall.

I traced the Jacobi stall to subnormal column norms causing ζ² to overflow and the rotation to degenerate—I'll fix it with overflow-free scaling and an ε·‖R‖_F negligibility cutoff. The second failure is semantic, not numerical: the old manifold helper zeroes U columns when σ < 10⁻¹² at the cut locus, which the corrupted-block tests depend on, so I'll preserve that rule in the adapter.

LinAlg (115/115) and Geometry (49/49) suites pass. Now I'll add regression tests for the stalled shape, an extreme-scale case, and the non-finite guard, then rerun the remaining suites.

The previously failing SPRED classes now pass, 34 of 34. Next I'll delete the helpers the migration left dead, then do the final run of all six suites on that state in the background.

They have no callers left. I'll delete them in a separate commit after the migration commit; removing dead code doesn't change behavior, so the full run now covers both states.

About 30 minutes, mostly DimReduction. I'll write the report once the results are in.

All six test suites pass in the final state, and I've restored the pilot JSONs the suite rewrote; two unrelated files in `issues/regression/` are left untouched. Now committing the changes in pieces.

Now the cleanup commit: the dead helpers in `Pca` and `GrassmannManifold`.

Last step: the numerics checklist and the `Pca` status row.

N6, N7 and the move onto the new SVD are all done and committed (`09572c6` to `e2b1474`). All six affected test suites pass. Of the numerics items, only N8 (the optional SMC core) and N9 (the baseline rerun) remain.

- **N6, [DiagonalPlusLowRank](src/maths/linalg/DiagonalPlusLowRank.cs):** the calibration covariance Σ = D + U·diag(s)·Uᵀ. It is built on the thin SVD instead of the usual Woodbury identity. That way the quadratic form and log-determinant come out as sums of nonnegative terms, and no m × m matrix is ever formed. The test shows the difference: along a heavily weighted direction it keeps the quadratic form to 10⁻¹⁰ relative, and Woodbury's subtraction loses every digit.

- **N7, [KhatriRao](src/maths/linalg/KhatriRao.cs):** the direct functional design's Gram matrix, the border a birth needs, and the cross-products with the response, all computed through A and H without building the full design.
  - It has both R's form and the centred Liang form, so it doesn't wait on your intercept decision.
  - It also handles observations at scattered (case, location) pairs.
  - The brief's name, "face-splitting", isn't quite right for BASS's design. BASS's design is the column-wise Khatri–Rao product (its transpose is the face-splitting product); ABMRS's rows really are face-splitting products. I corrected the brief's wording.

- **Moving the callers turned up two problems:**
  - **A real SVD bug.** The Grassmann manifold handed `ThinSvd` two parallel columns with norms 1.8·10⁻¹⁶¹ and 0.115. An intermediate value overflowed and the algorithm looped until it gave up. It now computes the rotation in a way that can't overflow, ignores columns too small to matter, and rejects NaN or infinite input. Regression tests cover these cases.
  - **A behavior the manifold tests relied on.** The old manifold SVD set some vectors to zero when the subspaces sit at right angles. At that point the log map's direction is ambiguous, and zeroing makes the two choices cancel. Two MoMPCA tests and six SPRED tests depend on that, so the manifolds keep the rule.

- **`Pca` is now faster.** On a 6,238 × 617 matrix (ISOLET size) it takes 8.0 s instead of 22.3 s, with the same eigenvalues to 2·10⁻¹³. Only wide inputs are slower: 0.53 s against 0.32 s at 300 × 2,000.

  Component signs now follow a fixed rule (largest entry positive). Anything that depended on the old arbitrary signs could flip, though no current test did.

- **Final suite results:**

  | Suite | Result |
  |---|---|
  | LinAlg | 118 of 118 |
  | Geometry | 49 of 49 |
  | Mapper | 28 of 28 |
  | Oracle | 5 of 5 |
  | VizCore | 240 passed, 2 skipped by design |
  | DimReduction | 83 of 83 |

  The SPRED ISOLET pilot reproduces its starting objective to 10⁻¹⁴. Its final objectives shift at the 10⁻⁵ level, because its random search follows a different path once the bits differ.

**Two side effects you should know about:
- **The DimReduction suite rewrites committed files.** It takes 18 minutes and writes its pilot results into the tracked files in `issues/spred/pilot/`. I restored them after each run.
- **The Oracle suite was flaky once.** Before any of these changes it failed 2 of 5 once, then passed on a rerun.**

I left `issues/regression/clumped-problems.md`, `error-capture-proposal.md` and `issues/spc/opus-gvdv-ideas.md` uncommitted, since they aren't mine.

# Build the extension: rotate MoMPCA's basis to the principal axes of the median of the blocks' covariances projected into the robust subspace, ordered and signed like Pca, with those robust variances in the result.

Builds. Now the tests: a fixture where two of five blocks inflate the variance along y. All-data PCA would then make y the top axis, while the robust axes keep x on top with the clean variances.

The test expectation was too tight—the geometric median's balance point actually lands around 0.56-0.58 for y, matching hand calculations. I'll widen the bounds accordingly and add a direct comparison against the all-data PCA value (~76).

All Geometry (51/51) and LinAlg (118/118) tests pass. I moved Pca's sign rule into a shared public Pca.Orient so other consumers get identical behavior; now updating the dev-sequence note and committing.

The extension is built and committed (`a997a10`). The Geometry suite passes 51 of 51 and LinAlg 118 of 118. `MoMPcaResult` is now `(Mean, Basis, Variances, Alpha)`, and its `Basis` holds robust principal axes:

1. Each block's covariance is taken about its own mean, with the same n denominator as `Pca`, and projected into the robust subspace.
2. Those projected covariances are combined by a geometric median under the Frobenius norm, the same median routine `MoMPCA` already uses.
3. The median's eigenvectors become the axes, in descending order of eigenvalue and signed like `Pca`'s. The eigenvalues are the `Variances`.

The span is unchanged, so the paper's estimate is untouched, and all the subspace tests still pass. With a single block the result is exactly the full-data PCA, eigenvalues included. `Pca`'s sign rule is now shared as a public `Pca.Orient`, with the same logic, so `Pca`'s results don't change.

**The new test.** Five blocks are ellipses in the xy plane. Three clean ones have variances near (2, 0.5); in the two corrupted ones the y-variance is near 200. Over all the data that makes y the top axis, and the test confirms plain `Pca` puts y first. `MoMPCA` keeps x first, with variances of about 2.0 and 0.56.

**One thing to know about the robust variances.** I first expected about 0.5 for y, and the test failed at 0.561. That turned out to be correct behaviour, so I widened the tolerance rather than changing the code. A geometric median is pulled toward outliers by roughly the clean blocks' own spread, not by how extreme the outliers are. Here the clean x-variances range from 1.9 to 2.1, which puts y at about 0.56; solving the median's balance condition by hand gives about the same. All the data together give 76 for the same quantity.

When the median's eigenvalues tie, the axes between them are arbitrary. That's the same instability the paper gives as its reason for not estimating a basis. Both points, and the fact that this goes beyond the paper, are written up in [dev-sequence.md](issues/spred/dev-sequence.md).

# is there a way to tie break?

Not in a way that comes from the data. Under an exact tie, the median covariance is the same in every direction within the tied block, so every orthonormal basis of that block fits equally well. Any rule that picks one either breaks rotation equivariance or brings in information from outside the variances. That's a property of the question rather than the method. It applies equally to plain `Pca`, and it's why the paper won't estimate a basis.

What you can do:

1. **Report the tied group as a subspace, not as axes.** Eigenvectors inside a tied group are unstable. The projector onto the whole group is stable as long as the group is well separated from the other eigenvalues (Davis–Kahan). The open question is when two eigenvalues count as tied, and `MoMPCA` already holds the data to answer it. Each block's projected covariance gives that block's variance along each robust axis. If the gap between two robust variances is smaller than their spread across blocks (a MAD, say), treat them as one group. The paper's node bootstrap (§6.3) is the heavier version of the same test.

2. **Fix a basis within a group by convention.** If a consumer really needs axes:
   - **Coordinate-aligned:** the first axis is the coordinate direction with the largest projection into the group, then the next, orthonormalized. It's deterministic but tied to the coordinate system.
   - **Procrustes to a reference:** align to a reference basis, such as the full-data PCA or a previous run's axes. This is what the paper does for its embedding figures. It changes smoothly with the data as long as the reference stays fixed.

3. **Break the tie with real information, from a secondary criterion inside the tied group:**
   - **Non-robust full-data variance:** this lets the outliers back in, but only among directions the robust estimate considers equal.
   - **Higher moments, as ICA does:** these resolve genuinely non-Gaussian structure. Data that really are rotationally symmetric, like the circle blocks in the tests, stay tied, because there is nothing to find.

A convention fixes the exact tie but not the real problem, which is near-ties. When two variances are close, the axes swing a lot under small changes to the data, and no convention stops that. Any rule for choosing a frame also has to jump somewhere, because no continuous choice exists over all subspaces.

I'd do 1 and 2 together:
- group eigenvalues whose gap is inside the block spread;
- return the variances with group labels, so consumers can tell which axes are trustworthy individually;
- inside each group, fix the basis by Procrustes to the full-data PCA, falling back to coordinate alignment.

The tie test only works in `MoMPCA`, which has blocks. For `Pca` the honest options are to report the gaps or to add a bootstrap.

# this reminds of something from one of the voroninski papers where there's some kind of algorithmic jump that happens cyclically that avoids a degenerate condition

You're probably thinking of the negation check. It's in Heckel, Huang, Hand & Voroninski, *Deep Denoising* (arXiv 1805.08855, Algorithm 1), and grows out of the loss-landscape result in Hand & Voroninski 2017 (arXiv 1705.07576). Both are in your gauntlet library.

- **The degenerate condition.** The loss has its global minimum at x₀, plus a critical region near a negative multiple −ρx₀. There the curvature is positive along ±x₀ and zero in every orthogonal direction, and the loss is higher than near x₀. Plain gradient descent can stall in that flat region.
- **The jump.** Before every iteration, the algorithm checks whether f(−xᵢ) < f(xᵢ), and if so negates the iterate. The degeneracy sits at a known, discrete, symmetry-related spot, and the objective can tell the two candidates apart. So the algorithm compares them and jumps every step, instead of trying to descend out of the flat region.

**It doesn't apply to eigenvalue ties.** The trick needs two conditions: the alternatives are a discrete set, and the objective differs between them. An exact eigenvalue tie fails both. Every rotation within the tied block is an alternative, and the variance is identical under all of them, so comparing objective values can't pick one. `Pca.Orient` already covers the one discrete case, the sign of an axis, and there too the objective can't tell the sides apart, so a convention is all that's possible.

**It does match a different degeneracy we just discussed: the Grassmann median at the cut locus.** When a block's subspace sits at a principal angle of exactly π/2 from the current iterate, the log map has two directions, ±.

- **What we do now.** The manifold adapter leaves that U column at zero, so the two directions cancel. That picks the zero element of the subgradient set. It's a valid subgradient, but it gives up on descent: at the cut locus the distance to that block is at a local maximum, and moving either way reduces it. The zero direction is exactly where Weiszfeld can stall.
- **What the negation check would do instead.** At such a direction, take a trial step both ways and keep whichever lowers Σₖ d(U, Uₖ). The other blocks break the symmetry, so unlike an eigenvalue tie, there is a real answer.
- **What we already have in that spirit.** `MoMPCA`'s medoid warm start is a one-shot version of the same idea: it avoids starting on a degenerate point. The per-iteration check would also cover iterates that drift onto a cut locus later.

The open question is whether that change improves robustness. With the cancelling convention, a corrupted block at the cut locus exerts no pull in that direction, which is convenient. With the objective-based choice it pulls one way, though only with Weiszfeld weight about 1/(π/2), against much larger weights from clean blocks close to the iterate. I'd expect the median to stay with the clean majority either way. That's a hypothesis, though. It would need testing on the cut-locus fixtures, plus a case where cancelling actually stalls.