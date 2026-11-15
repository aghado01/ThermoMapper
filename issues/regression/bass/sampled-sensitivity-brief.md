# Sobol indices and effects for terms without factors: E25's first build

Date: 2026-11-15. Status: **approved by the user on 2026-11-15** at its recommendations ("Yes, start it", to starting a chip on this brief at its recommendations). Drafted by the second coordinating session at the user's word that E36's dependencies come first, from [E25's entry](extensions-brief.md#sensitivity-under-dependent-inputs), the [lineage review](lineage-review.md#separability-is-what-the-closed-forms-rest-on)'s note on separability, [E36's brief](bppr-brief.md) (decision 10) and a reading of the sensitivity readouts' public surface. The estimators are cited from memory until slice 0 reads them. Reports go under "Reports" below, newest first.

## Why

- **A BPPR fit has no Sobol indices and no effects.** E36's directional hinge (`9e249dc`) reads several inputs through one projection, so it has no one-variable factors and the moment service refuses it. `Sensitivity`, `Effects` and the basis fits' versions refuse any fit holding one. Prediction, the gradient and E8's C already read such a fit.
- **The direction is settled.** The lineage review: a kernel without factors is sampled under any measure, with exact values per draw, so only the integral over x is approximate; Sobol, C and effects sit above both routes and know neither. E8 built that shape for C, two routes chosen per query. This brief gives Sobol and the effects the same shape.
- **Scope: the kernel side.** A product law, and terms without factors. E25's measure side (dependent inputs: accumulated local effects, Shapley effects, G4's mixture of products) stays in its entry for a later brief, and nothing here closes it off.

## Sources

- **In the repository:** E25's entry and [G4](lineage-gaps-register.md#g4-one-moment-service-across-dependent-inputs); E8's [build brief](active-subspace-brief.md), for the two routes, the shared points and the error rows; E36's brief and its slice 0 and slice 1 reports.
- **The estimators, cited from memory and not known to be on disk:** Sobol' 2001 and Saltelli et al. 2010 for the pick-freeze forms; Jansen 1999 for the total index; Janon et al. 2014 for their asymptotic variances; Owen 2013 for small indices; Sobol' and Kucherenko 2009 for the derivative bound on a total index. Slice 0 reads those the user supplies and marks the rest.
- **Oracles.** R's `BASS` gives Sobol indices for separable fits only, which the exact route already meets. Whether `BayesPPR` 0.2.0 has a sensitivity readout is for slice 0 to read from the pinned package.

## What is already built

- **The readouts' shapes:** `SobolAnalysis` over per-draw `SobolDecomposition`s (variance, main variances by subset, total variances by input, the rest, and a status), `EffectAnalysis` over grids, the basis fits' `BasisSobolAnalysis` and `BasisEffectAnalysis`, and `SensitivitySummary`.
- **The exact route:** `BassMoments` under a `ProductMeasure`, and `DesignKernels.Separable`, which says when it applies.
- **A sampled route to copy:** E8's. N points drawn once from the product law and shared by every draw, so the draws' spread is the posterior's; each entry carries its Monte Carlo standard error; `Auto` takes the exact route when every term has factors.
- **Truth from outside the model:** E36's references (`ProjectionPursuitReferences`), E8's `ReferenceFunction` family and Friedman's function.

## Decisions

Each was open with a recommendation, and the user approved the brief at its recommendations, so the recommended form is the one to build. The alternatives stay recorded for the page each decision lands on, with its trade-offs, when built.

1. **The route.**
   - *Recommended:* `SensitivitySpec` gains a route as E8's spec has: `Auto`, `Exact`, `Sampled`. `Auto` is exact when every term has factors and sampled otherwise. `Exact` refuses a fit without them. `Sampled` runs on any fit, so a separable fit checks the estimator against the closed forms.
2. **What is sampled.**
   - *Recommended, staged:* first the whole surface. Later, and only if mixed bases are built (E28), the split f = f_sep + f_ns: the separable part's pieces exact from the moment service, and only what involves f_ns sampled. That is a control variate, the closed forms supplying efficiency.
   - *Alternative:* quadrature per ridge. A direction has at most three active coordinates, so a ridge's conditional means are integrals in at most two dimensions and a pair's moment in at most six. It has no Monte Carlo error, but it holds for sparse directions and a product law only, and the hinge's kink sits across the grid.
3. **The estimator.**
   - *Recommended:* pick-freeze on two base point sets A and B shared by every draw. A main variance by Saltelli's 2010 form, a total by Jansen's, a requested subset by one more block. First-order and total indices cost N(p + 2) surface values per draw; columns are computed once per distinct term, as the term table shares them across draws.
   - *Reported with it:* each variance piece's standard error, and each index's by the delta method, as E8 reports C's.
   - *Alternatives:* a given-data estimator from one sample (cheaper, biased, with a bandwidth to choose); the derivative bound on a total index from E8's C alone (free where C exists, a bound and not an index).
   - *Points:* the engine's generator from the spec's seed, as E8's route. A low-discrepancy set is an option to take up only if the error rows call for it.
4. **A piece that comes out below zero.**
   - *Recommended:* report the estimate and its error as they are, with a status flag where the estimate is below zero by more than its error allows. Nothing is clamped: a clamp hides how large the error is.
5. **The effects.**
   - *Recommended:* E[f | x_S] on the grid by averaging over the shared base points with x_S pinned, each grid value with its standard error, in the `EffectAnalysis` shape. The error is null on the exact route.
6. **The measure.**
   - *Recommended:* a product law only. A joint sample is refused with its reason: under dependent inputs Sobol's decomposition is not defined (G4: Hoeffding is lost), and E25's other readouts are the route.
7. **The Poincaré check.**
   - *Recommended:* restore it for a fit holding a directional hinge. E36's slice 1 left it out because its totals came from the moment route; the sampled totals supply them, with the bound read within their error.
8. **The basis fits, functional forms and several outputs.**
   - *Recommended:* after the scalar fit, in a slice of their own, each carrying the route and the error rows through.

## The slices

One commit per concern; tests with the code they test; each slice ends with its gates and a fresh-eyes review before the next starts.

| Slice | Delivers | Its own gates |
| --- | --- | --- |
| 0. Sources, references and the estimators | The estimators read; `BayesPPR` read for a sensitivity readout; reference functions with closed-form Sobol indices and effects from outside the model (E36's references, Friedman's function, Ishigami's function, Sobol's g-function); the pick-freeze primitives over any function | Each estimator against the closed forms on chosen cases and a seeded sweep: no bias beyond its stated error, the error falling as N^(−1/2), small indices and a zero-variance function included |
| 1. Sobol on the sampled route, a scalar fit | Decisions 1, 3, 4, 6 and 7; the decomposition with its error rows; `SensitivitySummary` | A separable fit's sampled indices against its exact ones within the stated error, on every existing reference; a hand-built directional fit against its closed form; the exact route's results bit for bit; the Poincaré check on a directional fit |
| 2. Effects on the sampled route | Decision 5 | Sampled effects against exact ones on separable fits; a directional fit's against the closed form; the grid's mean against the surface's mean |
| 3. Basis fits, functional forms and several outputs | Decision 8 | Each readout's sampled route against its exact one on separable fits; E34's and D07's references |
| 4. The split, if mixed bases are built | Decision 2's second stage | The split against the whole-surface estimate within both errors; the error's reduction stated; the cost row decides its default |

## Gates common to every slice

- **The exact route untouched.** With the route at `Auto`, every fit built before E36 reads bit for bit as before: the pinned tier at each commit, the whole assembly at the slice's end.
- **Truth from outside the model.** Reference indices and effects come from closed forms or from a method apart from the estimator under test.
- **Primitives.** Chosen cases and a seeded random sweep, far ranges included.
- **The cost row.** Per draw against the exact route on a separable fit, and against prediction: surface values, wall clock, memory. Indicative, with the load beside it.
- **Docs.** Chips leave `docs/` alone; each report lists the doc changes owed.

## Not in this brief

- **Dependent inputs:** accumulated local effects, Shapley effects and G4's mixture of products, which are E25's measure side.
- **E24's guarded estimators.** Their guarantees cover linear quantities; Sobol shares depend on the emulator (E24's own limits).
- **Mixed bases themselves** (E28), and a ridge times a categorical factor (E36's decision 7), which are fitting.
- **Shapley effects under a product law.**

## Working rules

As [E36's brief](bppr-brief.md#working-rules) has them, with `artifacts\chip-sensitivity\` for builds. "The main session" is the second coordinating session, which verifies each slice from a clean archive. Files: `src/maths/regression/spline/bass/sensitivity/**` and their tests. The BPPR chip writes in BASS's structure chain and has edited three of the active-subspace files; E8's chip owns them. A change inside either's files is agreed by message first.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; decisions taken beyond this brief; what was left unbuilt and why; follow-ups with the decision each would change; the doc changes owed.

## Reports

### Slice 2 verified from a clean archive (coordinating session, 2026-11-15)

- **Commit:** `9d96c34`, on a snapshot of `84fa2e6`, through `scripts/verify-csharp.ps1` over its five changed files, in the two projects that select.
- **Result:** Regression 2,760 of 2,760 and Oracle 57 of 57 at commit scope, `SampledEffectsTests` 11 and `SampledSensitivityTests` 18 among them, and E8's `SampledPoincareCheckTests`, which calls this route.
- **The gates' cuts, checked and let stand; the user is told they loosen.**
  - *The bias cut by Bonferroni at a family-wise 10⁻³.* At 5,910 cells a flat 4.5 fails by chance about 3.5% of the time per fit, 5,910 × 2Φ(−4.53), so the 4.53 that tripped it is the size that cut should have expected. The cell reads −0.53 over 200 fresh seeds, with spread 1.011: no bias.
  - *Its cost:* at 24 seeds the smallest bias a fit-level cell can show moves from about 0.92 of a standard error to 1.07. The primitives' gates at 300 replicates read 0.23.
  - *The pooled spread's band, [0.85, 1.15],* is pre-emptive: none of the eight tables was outside [0.9, 1.1]. It rests on those tables' own range, 0.909 to 1.082, which makes ±0.1 about a two-sigma band.
- **Its truths are outside the estimator:** the exact route on separable fits; tensor quadrature on the hand-built directional fit; and quadrature of a chain-built K = 4 ridge fit's own predicted surface, read at this session's request.
- **Receipts:** `artifacts/verify-sens-s2-2b951c1/`, with the changed list, the counters and the snapshot tracked in [receipts](../receipts/verify-sens-s2-2b951c1/).

### Slice 2: effects on the sampled route (2026-11-15)

**Landed** in `9d96c34`, which carries its tests. It does not depend on decision 3's form: an effect has no main-variance form.
- *Decision 5, effects by pinning.* `BassSensitivity.Effects` resolves the route as `Of` does (decision 1), and the sampled route reads a scalar fit.
  - *The cell:* each cell is the effect's ANOVA component, the mean over the 2N shared base rows of the per-row alternating sum Σ_{v⊆u} (−1)^{|u|−|v|} f(g_v, x_{−v}). The base rows are A and B, the points the Sobol route draws from the same seed.
  - *The estimator:* unbiased under the product law, its error the sum's sample deviation over √2N, `PinnedEffects`' estimator from slice 0. The signs add to zero, so the intercept drops from the sums, and each draw's mean takes it back.
  - *The order of evaluation:* the base block and each input's blocks at its grid values are kept per draw. A pair's cells are evaluated in chunks and reduced against them at once.
  - *The memory:*
    - draws go in batches whose kept blocks stay under 8M values;
    - chunks go under the column budget slice 1 set;
    - each chunk's term columns are read through prediction's `BassSurface`.
- *What changed in the records:*
  - `Effect` gains `Errors` (draws × cells); `EffectFunctions` gains `MeanError`; `EffectAnalysis` gains `Route` and `Points`. Each is null or 0 on the exact route.
  - `Points` is N, the rows in each base set, as on `SobolAnalysis`; the effects average over both sets.
- *`EffectSpec.EvaluationBudget`,* 2·10⁹ term values by default, refuses up front a request whose pinned blocks would evaluate more, naming the knobs.
- *The request, shared by both routes:* `EffectRequest` (the effects asked for, or the default list of mains and the pairs some term holds; the grids; the cell budget) moved out of `Effects` and is called by both. It finds a term's variables through slice 1's `TermVariables`, so the default pairs come from a directional term's or a ridge column's active coordinates too. The exact route's code is otherwise as it was.
- *The gates' cuts, in both slices' fit-level tallies:*
  - The bias cut is Bonferroni's at a family-wise level of 10⁻³ (`BiasCut`), so it grows with the cells read: 4.2 at 42 cells, 5.2 at 6,000.
    - *Why:* at 5,910 cells a flat 4.5 fails by chance about 3.5% of the time per fit (5,910 × 2Φ(−4.53), the coordinating session's figure).
    - *What it costs:* at 24 seeds the smallest bias a fit-level cell can show moves from about 0.92 of a standard error (4.5/√24) to 1.07 (5.23/√24). The primitives' gates, at 300 replicates, still read 0.23.
  - The pooled spread of z is held within [0.85, 1.15]. The cells of one seed share its points, so that spread has about the seeds' few degrees of freedom. Slices 1 and 2's eight separable tables read 0.909 to 1.082, which makes ±0.1 about a two-sigma band. None of them was outside [0.9, 1.1], so the wider band is pre-emptive. The primitives' gates (slice 0, 300 replicates) read calibration sharply.
  - The trigger: under the old flat cut of 4.5, one cell of the 5,910 in "hinges, the boxes" read 4.53. Over 200 fresh seeds the same cell reads √R·mean z −0.53 and sd of z 1.011 (probe `artifacts/chip-sensitivity/probes/EffectCellProbe.cs`): chance, the extreme of many correlated cells.
  - Slice 1's thresholds move with them, and slice 1's numbers stand under both.

**Gates.** From a snapshot of `HEAD` (`f62b4ab`) with the five files, unless a row says otherwise.
- *Separable fits, sampled against exact.* Every draw's mean and every cell of every default effect (mains, and the pairs the terms hold, categorical levels among them), on 7-point grids, over 24 seeds at N = 1024 (2,048 base rows). The fits are slice 1's four.

  | Fit | Cells | Exact zeros | Worst √R·\|mean z\| (cut) | Pooled sd of z | Per-cell sd | Beyond 3 | Worst z |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | hinges, mixed laws | 5,910 | 630 | 1.57 (5.23) | 0.909 | 0.61 to 1.31 | 0.02% | −3.35 |
  | hinges, the boxes | 5,910 | 630 | 4.53 (5.23) | 0.963 | 0.67 to 1.27 | 0.27% | 4.34 |
  | degree 2, empirical and Student | 6,540 | 0 | 2.89 (5.25) | 0.952 | 0.65 to 1.28 | 0.18% | 3.31 |
  | Legendre | 10,776 | 2,744 | 2.58 (5.34) | 1.002 | 0.51 to 1.37 | 0.39% | 4.03 |

  The exact zeros are inputs a draw's terms do not read, zero on both routes. The two routes list the same effects on the same grids.
- *Two fits without factors, against quadrature,* a method apart from the estimator.
  - *Slice 1's hand-built directional fit:* tensor quadrature of its surface, 48 panels. Its six default effects are three mains on 6-point grids and the three pairs its hinges hold. Over 40 seeds at N = 2048:

    | Degree | Cells | Worst √R·\|mean z\| | Pooled sd of z | Per-cell sd | Beyond 3 | Worst z |
    | --- | --- | --- | --- | --- | --- | --- |
    | 1 | 254 | 2.58 | 1.015 | 0.74 to 1.24 | 0.30% | 3.70 |
    | 2 | 254 | 1.79 | 1.032 | 0.73 to 1.22 | 0.37% | 4.10 |

  - *A ridge-block fit from the chain,* E36's slice 3, at the coordinating session's request.
    - `Structure.Directional = new DirectionalBasis { Columns = 4 }` on 120 cases of exp(0.6x₀ + 0.8x₁) + sin 3x₀ with noise sd 0.05, on [0, 1]².
    - Six draws, every term a ridge column, up to 8 columns (two blocks).
    - The truth is quadrature of the fit's own predicted surface, 3-node Gauss–Legendre on 64 panels per input; halving the panels moves it by 1.1·10⁻⁵ against effects up to 1.28.
    - Over 40 seeds at N = 2048: 384 cells (mains on 7 points and the pair), worst √R·|mean z| 2.26, pooled sd 0.920, per-cell sd 0.70 to 1.15, none beyond 3, worst z 2.82.
- *The grid's mean against the surface's mean.* On 25-point grids, under the uniform law on the boxes, a main effect's mean over its grid is 0 up to the grid's trapezoid error; under the level probabilities it is 0 exactly.
  - The sampled grid mean equals the truth's within its weighted error, a bound on its own error: at worst 1.60 of it on the separable fit, 0.63 on the directional fit (degree 2, quadrature truth).
  - The categorical input's exact grid mean is 0 to 10⁻¹².
  - Each draw's sampled mean sits within four of its errors of the exact mean.
- *The route:* `Auto` on a separable fit is the exact route bit for bit (values and mean), with null errors, `Route = Exact` and `Points = 0`.
- *The shared points:* the seed reproduces the analysis bit for bit, and another seed moves it. A draw read alone is its row, cells, errors and mean, bit for bit. The effects are named in that test, since the default pairs follow the selected draws' terms.
- *The refusals:* `Exact` names the term without factors; `Points` below 2; the evaluation budget, with its count; a functional fit, naming slice 3.
- *The exact route's own tests:* `BassSensitivityEffectTests` (12 of 12) after the request moved into `EffectRequest`.
- *The cost row* (wall clock and allocation per draw, the machine at 1 to 28% load, so indicative; probe `artifacts/chip-sensitivity/probes/EffectsCostProbe.cs`, not committed):

  | Per draw | Exact | Sampled, N = 1024 | Sampled, N = 4096 |
  | --- | --- | --- | --- |
  | Separable fit (60 draws, 16 distinct terms), mains on 100 points | 0.020 ms | 15 ms, 6.3 MB | 70 ms, 54 MB |
  | The same, mains and pairs on 25 points | 0.029 ms | 46 ms, 6.3 MB | 191 ms, 64 MB |
  | Directional fit (2 draws), mains on 100 points | — | | 128 ms, 188 MB |
  | The same, mains and pairs on 25 points | — | | 808 ms, 1.08 GB |

  - *What drives it:* pinning costs cells × rows × terms, against the exact route's cells × terms, so effects run 3 to 4 orders of magnitude above the exact route. The evaluation budget bounds a request.
  - *The allocation column* is churn, not peak: each chunk rebuilds its stacked cases. Peak memory stays within the store and column budgets.
  - *A posterior with many distinct terms* (a chain over continuous directions) shares few columns across a batch. Its cost per draw approaches that of evaluating every pinned block for every draw alone.
- *Fitting untouched:*
  - *Selected at commit scope:* test-select over the five files, per test project from disk, selects Regression 252 of 347 classes and Oracle 13 of 40. No other project holds the files in its closure.
  - *Regression at commit scope (with the pinned tier):* 2,756 of 2,756.
  - *Oracle at commit scope:* 57 of 57.
  - *The whole assembly:* Regression 3,672 of 3,672, in 19 minutes at below-normal priority.
  - *After `HEAD` moved:* E8's `7cca408` (decision 7, calling this route) landed during the runs, so the commit went onto `8236e0f`. On a snapshot of it, the classes this slice touched passed 77 of 77, with E8's `SampledPoincareCheckTests`, `DirectionalHingeTests`, `RidgeColumnTests`, `BassSensitivityEffectTests` and `PickFreezeTests` among them.

**Decisions taken beyond the brief.**
- *The base rows are A and B,* the Sobol route's 2N points from the same seed, so a fit's indices and effects at one seed read the same points.
- *The evaluation budget.* The exact route's cell budget bounds the values held; the sampled route's cost is set by the rows and terms behind each cell, so a second budget bounds it.
- *The request shared by both routes,* rather than copied.
- *The gates' cuts* by Bonferroni and by the pooled spread's degrees of freedom, as above.

**Unbuilt.**
- *Basis fits, functional forms and several outputs:* slice 3.
- *Reusing a separable posterior's term columns across all its draws* within the store budget. Today a batch of draws shares them, and the batch is sized by the stored blocks. Worth building only if separable fits are read on the sampled route for more than a check: `Auto` reads them exactly.

**Follow-ups.**
- *A smaller default grid on the sampled route.* The default (100 points, all pairs the terms hold) evaluates about 10⁹ term values per batch on a three-input fit, a minute or more. Decision: whether `EffectSpec` takes a route-dependent default (say 25 points), or the budget's refusal suffices. Either keeps the exact route's default.
- *Reusing the stacked cases across batches* would remove most of the allocation churn without moving a bit. Decide against the cost row if effects of long posteriors are read often.

**Doc changes owed.**
- `docs/Maths/Regression/Sensitivity.md`:
  - decision 5's estimator and its derivation;
  - the base rows shared with the Sobol route;
  - the evaluation budget, and the order of evaluation (kept blocks, chunked pairs);
  - the error rows on effects;
  - the cost against the exact route.
- The sensitivity page's testing notes: the gates' Bonferroni cut and the pooled spread's degrees of freedom, as a convention for stochastic gates over many correlated cells.

### Slice 1 verified from a clean archive (coordinating session, 2026-11-15)

- **Commits:** `5f9cdef` and `c5cea88`, on a snapshot of `07f3c2a`, through `scripts/verify-csharp.ps1` over their ten changed files, in the two projects that select.
- **Result:** Regression 2,771 of 2,771 and Oracle 57 of 57 at commit scope. The slice's classes are among them: `SampledSensitivityTests` 18, `PickFreezeTests` 23, and the BPPR chip's `DirectionalHingeTests` 12 and `RidgeColumnTests` 6 with their changed assertions.
- **Read.**
  - *The dispatch:* `OfCore` resolves the route and hands a sampled one on; `FromCoefficients` is as the BPPR chip left it.
  - *`Auto`'s rule* asks the kernels (`DesignKernels.Separable`) over the terms the selected draws hold, and `Exact`'s refusal names the term it meets and no kind of term.
  - *The two-difference form for a subset:* E[D^B_u D^A_u] = 2^|u|·(−1)^|u|·V_u, each d ⊆ u arising as v Δ w in 2^|u| ways.
  - *The stated errors' calibration:* the z-scores' pooled spread (0.95 to 1.08) and their share beyond 3 (at most 0.57%) are what about 1,200 cells over 24 seeds give under calibrated errors; the per-cell spreads' range is the noise of a deviation read from 24 replicates.
- **Its truths are outside the estimator:** the exact route on separable fits, and tensor quadrature on the directional fit, with the quadrature's own error read beside it.
- **Open with the user:** decision 3's default form. Saltelli's stands until they rule, and both forms are gated.
- **Passed on:** decision 7's change is in E8's readout. The first coordinating session has recorded it in E8's brief; E8 takes it after its slice 8a, starting from this note.
- **Receipts:** `artifacts/verify-sens-s1-07f3c2a/`, with the changed list, the counters and the snapshot tracked in [receipts](../receipts/verify-sens-s1-07f3c2a/).

### Slice 1: Sobol on the sampled route, a scalar fit (2026-11-15)

The coordinating session's note on starting slice 1: build what decision 3's form does not touch. The form of a main variance is with the user, so it stays a choice made in one place, `PickFreezeDesign`, and slice 1's gates read both forms. Until the user rules, the default is Saltelli's centred form, the brief's approved decision 3.

**Landed** (each commit carries its own tests):
- `5f9cdef`: the form as a choice (`MainVarianceForm`). It also registers slice 0's centring comparison in `ComparativeSweep`.
  - `PickFreezeDesign` takes the form, Saltelli's by default.
  - The two-difference form reads each block's reverse hybrid B_A^u (`Reverse`), N(2p + 2) values for first-order and total indices.
  - `PickFreeze.Estimate` takes the reverse hybrids, required by that form and refused by Saltelli's.
  - *The form:* (−1)^{|u|}·mean(D^B_u D^A_u)/2^{|u|}. f(B_A^v) and f(A_B^w) share x on v Δ w, so E[D^B_u D^A_u] = Σ_{v,w⊆u} s(v)s(w)τ_{vΔw} = 2^{|u|}(−1)^{|u|}V_u, since each d ⊆ u arises as v Δ w in 2^{|u|} ways. It is unbiased and shift-invariant.
  - *Its error:* sd/√N of its per-row term, which is a product of two differences, each small when the piece is.
  - *The registration:* the comparison moves into `PickFreezeTests.CentringSpread(shift, seed)`, as `CONTRIBUTING.md`'s rule on comparative assertions asks. Slice 0's test asserted it inline; my review found that.
- `c5cea88`: the sampled route for a scalar fit.
  - *Decision 1, the route:* `SensitivitySpec` gains `Route` (`Auto`, `Exact`, `Sampled`), `Points` (N = 4096) and `Seed` (1).
    - `Auto` is exact when every term the selected draws hold has factors (`DesignKernels.Separable`), sampled otherwise.
    - `Exact` refuses a term without factors up front, naming the term it meets: "is not a product of one-variable factors". It names no kind of term, so the refusal holds for directional hinges and ridge columns alike.
    - `Sampled` runs on any scalar fit.
  - *The points:* 2N rows from E8's sampler (`BassActiveSubspaces.Points`) from the seed. A is the first N rows and B the next N, so the rows are independent and every draw reads the same points.
  - *The evaluation:* A, B and each block's hybrid (and reverse) are stacked as one set of cases.
    - Each distinct term's column is read once per batch of draws through prediction's `BassSurface`, the batch sized so the columns stay under 8M values.
    - Each draw's surface is evaluated with the intercept set to 0: every piece is shift-invariant, and the zero-variance rule then reads the same magnitude the exact route reads.
    - A singleton block on an input no term reads takes A's values (B's for its reverse), bit for bit, and is not evaluated.
  - *The subsets:* as the exact route reports them, those inside some term up to the order. A term without factors reads, as its variables, the continuous ones a case kernel says it reads (`IDesignKernel.ReadsVariable`), so ridge columns will be read the same way once they land.
  - *Decision 4, the error rows:* `SobolDecomposition` gains `VarianceError`, `MainVarianceError`, `TotalVarianceError`, `MainIndexError` and `TotalIndexError`, all null on the exact route.
    - `DecompositionStatus` gains `BelowZeroBeyondError`.
    - Nothing is clamped on the sampled route.
    - `SobolAnalysis` gains `Route` and `Points`.
  - *`SensitivitySummary`:* each piece's mean standard error over its draws under the weights, with the flagged draws' count, and `Route` and `Points` on the analysis summary. The draws share their points, so their errors correlate. The mean of the errors bounds the error of the posterior mean from above (the triangle inequality on standard deviations). It is not the posterior mean's own error.
  - *Decision 6, the measure:* `SensitivitySpec` holds a product law only. It has no field for a joint sample, so there is nothing to refuse.
  - *Out of this slice:* a functional fit and several outputs are refused on the sampled route, naming slice 3.
  - *Tests I don't own:* the BPPR chip's `DirectionalHingeTests` and `RidgeColumnTests` asserted that `fit.Sensitivity` refuses a fit holding a directional hinge or a ridge column. Each now asserts both halves of decision 1, at the coordinating session's direction: `Exact` refuses up front naming the term, and `Auto` returns an analysis on the sampled route. The other refusals (`Effects`, `FromCoefficients`) keep their assertions.
  - *Shared files:* the BPPR chip landed its hunks in `BassSensitivity.cs` and both tests (`332692b`) before this commit, so the commit holds only mine: `partial` and the route's dispatch line in `OfCore`, which leave `FromCoefficients` alone, and the two assertions.

**Gates.** All from a snapshot of `HEAD` (`72d3941`) with the slice's files and only my hunks in the shared files, unless a row says otherwise.
- *Separable fits, sampled against exact, both forms.* Every draw's every piece (variance, mains and main indices for every reported subset, totals and total indices), over 24 seeds at N = 2048. z = (sampled − exact)/stated error, per cell (draw and piece).
  - *The fits:* `SensitivityFixture`'s scalar hinge fit at degree 1 under the mixed laws (a two-normal mixture, a sub-range uniform, unequal level probabilities) and under the boxes; its degree-2 fit under an empirical marginal and a Student-t law; and a KHAOS Legendre fit.

  | Fit | Form | Cells | Exact zeros | Worst √R·\|mean z\| | Pooled sd of z | Per-cell sd | Beyond 3 | Worst z |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | hinges, mixed laws | Saltelli | 1,140 | 120 | 2.44 | 0.947 | 0.67 to 1.20 | 0.04% | 3.28 |
  | hinges, mixed laws | two-difference | 1,140 | 120 | 3.06 | 0.967 | 0.67 to 1.44 | 0.42% | −3.92 |
  | hinges, the boxes | Saltelli | 1,140 | 120 | 3.13 | 1.082 | 0.76 to 1.42 | 0.57% | 4.12 |
  | hinges, the boxes | two-difference | 1,140 | 120 | 2.93 | 0.973 | 0.75 to 1.31 | 0.23% | −3.76 |
  | degree 2, empirical and Student | Saltelli | 1,200 | 60 | 2.29 | 1.016 | 0.69 to 1.35 | 0.25% | −3.30 |
  | degree 2, empirical and Student | two-difference | 1,200 | 60 | 1.76 | 1.006 | 0.61 to 1.37 | 0.24% | −3.43 |
  | Legendre | Saltelli | 1,408 | 112 | 2.22 | 0.960 | 0.61 to 1.27 | 0.17% | 3.60 |
  | Legendre | two-difference | 1,408 | 112 | 3.49 | 1.027 | 0.73 to 1.35 | 0.29% | −4.32 |

  The cuts are: √R·|mean z| ≤ 4.5, the largest of about 1,200 cells near normal; pooled sd in [0.9, 1.1]; at most 1% beyond 3. Zero-variance draws agree on both routes, and the reported subsets are the same list.
- *A hand-built directional fit against tensor quadrature,* a method apart from the estimator. The fit has two draws: (Dense, −0.5) and (Pair, 0.3), and (Pair, 0.3), (Single, 0) and (Dense, 0.7), on DirectionalHingeTests' three inputs with unequal boxes, x₂'s declared. Quadrature uses 3-node Gauss–Legendre on 48 panels per input. Halving the panels moves the main variances by 9.9·10⁻⁶ of V at degree 1 and 1.3·10⁻⁸ at degree 2. Run over 60 seeds at N = 4096, 42 cells:

  | Degree | Form | Worst √R·\|mean z\| | Pooled sd of z | Per-cell sd | Beyond 3 | Worst z |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1 | Saltelli | 1.44 | 1.076 | 0.89 to 1.24 | 0.36% | −3.83 |
  | 1 | two-difference | 0.81 | 1.059 | 0.89 to 1.26 | 0.40% | −3.77 |
  | 2 | Saltelli | 1.38 | 1.044 | 0.86 to 1.24 | 0.60% | −3.90 |
  | 2 | two-difference | 1.17 | 1.045 | 0.91 to 1.18 | 0.40% | 3.79 |

  Brief's wording: it says the directional fit is checked against its closed form. A directional hinge under a product uniform law has one only through the projection's piecewise-polynomial law, so quadrature stands in, with its own error read beside it.
- *The exact route bit for bit:*
  - `Auto` on a separable fit is `Exact`'s arrays bit for bit, with `Route = Exact`, `Points = 0` and null error rows.
  - The fitting code is untouched; the commit scope's pinned tier and the whole assembly are below.
- *The points from the seed, shared by every draw:*
  - The same seed reproduces the analysis bit for bit, and another seed moves it.
  - A draw read alone (`Draws = [d]`) is its row in the full analysis bit for bit: its variance, totals and their errors.
- *The summary's rows:* each mean error is the weighted mean of the decomposition's own; the flagged count matches; on the exact route they are null.
- *The primitives under the two-difference form* (slice 0's gates, extended):
  - every reference, 300 replicates at N = 2048: worst bias 2.59 of its error, stated error / spread from 0.912 to 1.133;
  - the 300-case sweep: 5,198 z-scores, mean −0.082, sd 0.985, 0.31% beyond 3, worst −4.24;
  - the g-function's S = 7.2·10⁻⁵ input read to a relative spread of 3.1% (Saltelli's: 324%), and its pair piece 1.1·10⁻⁹ to 7.3% (Saltelli's spread is 1.3·10⁻⁶, 1,200 times the piece).
- *Decision 7, the Poincaré check, read here before E8's readout takes it.* Each draw's sampled total(v), at N = 8192, against C_P(v)·C_vv, with C from E8's sampled route at its own seed and C_P from `InputPoincareConstant.Of`. At degrees 1 and 2: no violation within three errors of the two independent estimates. The largest total/bound is 0.802 at degree 1 and 0.768 at degree 2; the largest (total − bound)/error is −11.2 and −7.5.
- *The cost row* (wall clock, `GC.GetTotalAllocatedBytes`; the machine at 9 to 31% load, so indicative; probe `artifacts/chip-sensitivity/probes/SampledCostProbe.cs`, not committed). The separable fit has 60 draws, 16 distinct terms and 3 inputs, one of them categorical, with subsets to order 3 (7 blocks).

  | Per draw | Values | ms | MB |
  | --- | --- | --- | --- |
  | Exact route | — | 0.059 | 0.004 |
  | Sampled, N = 1024, Saltelli | 9,216 | 0.91 | 0.074 |
  | Prediction at 9,216 cases | | 0.14 | 0.094 |
  | Sampled, N = 1024, two-difference | 16,384 | 1.00 | 0.099 |
  | Prediction at 16,384 cases | | 0.22 | 0.166 |
  | Sampled, N = 4096, Saltelli | 36,864 | 2.09 | 0.29 |
  | Prediction at 36,864 cases | | 0.51 | 0.37 |
  | Sampled, N = 4096, two-difference | 65,536 | 3.88 | 0.39 |
  | Prediction at 65,536 cases | | 1.24 | 0.66 |
  | Directional fit, N = 4096, Saltelli (degree 1; degree 2) | | 5.1; 4.3 | 3.1 |
  | Directional fit, N = 4096, two-difference (degree 1; degree 2) | | 6.6; 7.1 | 5.2 |

  The sampled route costs 35 to 65 times the exact route at N = 4096, and 3 to 4 times prediction at as many cases. The difference is the case preparation, the stacking and the reduction. The directional fit's two draws carry the fixed costs of drawing and preparing the points.
- *Fitting untouched:*
  - *Selected at commit scope:* test-select over the slice's ten files, per test project from disk, selects Regression 252 of 345 classes and Oracle 13 of 40. No other project holds the files in its closure.
  - *Regression at commit scope (with the pinned tier):* 2,771 of 2,771.
  - *Oracle at commit scope:* 57 of 57.
  - *The whole assembly:* Regression 3,651 of 3,651, in 17 minutes at below-normal priority. The commits went onto `d4b7b8e`, two documentation commits past `72d3941`, which touch no source or test file.
  - *An earlier full run on `2a6bb6b`:* Regression 2,757 of 2,758 at commit scope, Oracle 57 of 57, and the whole assembly 3,637 of 3,638. The one failure was `RidgeColumnTests`' old refusal, which `332692b` had just added; this slice's hunk replaces it.

**Decisions taken beyond the brief.**
- *The two-difference form built now,* as a mode of the design, so the user's choice of default changes one place. The coordinating session asked for the choice to stay in one place, and slice 1's gates read both forms. The public spec does not expose the form yet; the gates reach it through an internal `BassSensitivity.OfInForm`.
- *The intercept is left out of the sampled values,* as the exact route's magnitudes leave it out. Every piece is shift-invariant, and the zero-variance rule then reads the same scale on both routes.
- *One stacked set of cases per analysis,* each distinct term's column read once per batch of draws, the batch held under 8M column values (64 MB).
- *The draws are not checked on this route.* Prediction's `BassDrawFit.Pool` checks that a draw's term ids ascend; the sampled route reads draws as the exact route does, through `Setup`, which does not. A hand-built draw with unordered ids reads correctly here and is refused by E8's gradient readout; my first directional fixture had one.

**Unbuilt.**
- *Decision 7 inside E8's readout* (`BassActiveSubspace.cs` is E8's chip's): the check is read in a test above, and the change E8 needs is below.
- *Effects (slice 2); basis fits, functional forms and several outputs (slice 3).*

**Follow-ups.**
- *Decision 7, for E8's readout* (it changes E8's `Poincare` for fits holding a term without factors, scalar ones now):
  - Where `Poincare` is null for a fit that is not separable under a product law, read the totals from `BassSensitivity.Of(fit, data, new SensitivitySpec { Measures = spec.Measures, Route = SensitivityRoute.Sampled, MaxOrder = 1, Points = spec.Points, Seed = SeedTree.Derive(spec.Seed, 1)[0], Draws = spec.Draws })`.
  - `TotalVariance[d, v]` and `TotalVarianceError[d, v]`, v the readout's measure index.
  - `MaxOrder = 1` keeps the cost at N(p + 2), and a seed apart from C's keeps the two estimates independent.
  - Count a violation where total − C_P·C_vv > 3·√(se_T² + (C_P·se_C)²), beside the finite-element constant's relative-error allowance.
  - A functional, basis or several-output fit's totals come with slice 3.
- *The form's default* (decision 3), with the user. If they choose the two-difference form, it is the design constructor's default and nothing else. If they choose both as named modes, `SensitivitySpec` gains the field that `OfInForm` stands in for.
- *The posterior mean's own Monte Carlo error:* the summary gives a bound. The exact error needs the per-row terms averaged across draws before the deviation is taken, one more pass. Worth it only if the bound proves loose for a decision.

**Doc changes owed.**
- `docs/Maths/Regression/Sensitivity.md`:
  - the route (decision 1) with its alternatives;
  - the points from the seed, shared by the draws;
  - the stacked evaluation and its column budget;
  - the error rows and the flag (decision 4);
  - the summary's bound on the posterior mean's error;
  - the two-difference form, its derivation and its evidence beside Saltelli's;
  - the intercept left out;
  - decision 6's product law.
- `docs/Maths/Regression/ActiveSubspaces.md`: decision 7's check on a directional fit, once E8 takes it.
- `docs/Maths/Regression/BASS.md`: `Auto`'s rule, now shared by Sobol and E8.

### Slice 0 verified from a clean archive (coordinating session, 2026-11-15)

- **Commits:** `e013006`, `9578185` and `9940ab9`, on a snapshot of `1380f42`, through `scripts/verify-csharp.ps1` over the six files.
- **Result:** Regression 93 of 93 at commit scope, the slice's three classes and the pinned tier. Oracle holds two of the files in its closure and selects no class; no other project selects.
- **Read.** The main variance as Cov(f(B), D_u), by Möbius inversion of the closed variances; Jansen's total; and the two-difference form, whose expectation is 2τ_u on expanding its four products.
- **Checked from outside the engine** ([script and output](../receipts/verify-sens-s0-1380f42/)): the g-function at N = 4096 over 300 replicates, in numpy from the closed forms.
  - No estimator is biased: the largest |z| of a mean against its truth is 2.1.
  - *Relative errors* as S falls from 0.72 to 7.2·10⁻⁵: Saltelli's centred form 2.8%, 6.2%, 18%, 31% and 320% to 344%; the two-difference form 2.1% to 3.3% throughout.
  - *Gains at equal cost:* 1.1, 2.6, 21, 50 and 6,400 to 8,000. They agree with the report's table within the replicates' noise.
  - *Jansen's total* reads the two small inputs to 3.3% at the base cost. So under Saltelli's form an inactive input is still told from its total. What the two-difference form adds is the main variances themselves and the interaction pieces built from them.
- **For the user:** decision 3's form for a main variance, the report's first follow-up.
- **Receipts:** `artifacts/verify-sens-s0-1380f42/`, with the changed list, the counters, the snapshot and the outside check tracked in [receipts](../receipts/verify-sens-s0-1380f42/).

### Slice 0: sources, references and the estimators (2026-11-15)

**Landed** (each commit carries its own tests):
- `e013006` adds the references (`tests/regression/SobolReferences.cs`). Each one is a function from outside BASS's model with its ANOVA decomposition in closed form under the product uniform law on its box: the mean, every closed variance τ_u = Var(E[f | x_u]), and every conditional mean E[f | x_u]. The main variances, totals and effects follow by Möbius inversion, unless the family states them directly to keep their digits.
  - *Ishigami's function* (a = 7, b = 0.1). x₃ has no main effect but a total one, through V₁₃ = 8b²π⁸/225.
  - *Sobol's g-function* (a = 0, 1, 4.5, 9, 99, 99). Its two small indices are S = 7.2·10⁻⁵.
  - *Friedman's function* on [0, 1]⁶, through E[g | x₁] = 10(1 − cos πx₁)/(πx₁). With g = 10 sin(πx₁x₂), its closed forms are:
    - E[g] = 10·Cin(π)/π;
    - E[g²] = 50(1 − Si(2π)/(2π));
    - E[E[g | x₁]²] = (100/π²)(2π Si(π) − π Si(2π) − 4), by parts.

    Si and Cin are summed in decimal.
  - *E36's four ridge references* (`RidgeSobolReference`). The moment generating function about the box's centre is used, so every exponent stays within the projection's own range. τ_u is centred, as Π S(a + b) − Π S(a)S(b) telescoped through sinhc − 1, and each total is in its own closed form, so that small pieces keep their digits.
  - *A constant.*
  - *`SobolQuadrature`*, the check: tensor Gauss–Legendre quadrature of E[f | x_u] and of the centred τ_u. It is a method apart from both the closed forms and the estimators.
- `9578185` adds `PickFreeze` and `PickFreezeDesign` (`src/maths/regression/sensitivity/PickFreeze.cs`, `Maths.Regression.Sensitivity`). They work on any function, from values or from a function. The blocks are every input alone plus the downward closure of the requested subsets: first-order and total indices cost N(p + 2) values, and each requested pair one block more.
- `9940ab9` adds `PinnedEffects` and `PinnedEffectDesign` (`PinnedEffects.cs`): decision 5's estimator over any function. A cell's ANOVA component is the per-row alternating sum over v ⊆ u, with the cells running last input fastest, as `Effect` has them.

**The estimators and their errors.** The sources are not on disk, so the forms below are derived here; the citations are from memory.
- *Setup.* Rows j = 1…N are independent: a_j and b_j are drawn independently from the product law, and c^u_j = (b_{j,u}, a_{j,−u}) is the hybrid A_B^u. Write y^A, y^B and y^u for f's values on A, B and the hybrid. Rows of y^B and y^v share x_v; rows of y^A and y^B share nothing. So:
  - Cov(y^B, y^v) = τ_v and Cov(y^B, y^A) = 0;
  - Cov(y^A, y^v) = τ_{−v} = V − T_v.
- *A main variance* is V_u = Σ_{v⊆u} (−1)^{|u|−|v|} τ_v, the variance of f's ANOVA component on u.
  - Let D_u = Σ_{v⊆u} (−1)^{|u|−|v|} y^v, with y^∅ = y^A. Then E[D_u] = 0, since the signs sum to zero, and Cov(y^B, D_u) = V_u.
  - The estimate is their sample covariance over N − 1: exactly unbiased, and shift-invariant. At u = {i} it is Saltelli's 2010 form, mean(f(B)(f(A_B^i) − f(A))), with f(B) centred on its own mean. For a pair, D is the second difference y^{ij} − y^i − y^j + y^A.
  - Its variance is Var((y^B − μ)D_u)/N + O(N⁻²). The stated error is sd_j((y^B_j − ȳ^B)(D_j − D̄))/√N.
- *A total* is Jansen's 1999 form, T̂_i = Σ (y^A − y^i)²/(2N). It is unbiased, since E(y^A − y^i)² = 2V − 2τ_{−i} = 2T_i. Its stated error is sd(½Δ²)/√N.
- *The variance* is the sample variance over the 2N values of A and B. The per-row term is ψ_j = ½((y^A_j − m)² + (y^B_j − m)²), and the stated error is sd(ψ)/√N.
- *An index* is a piece over the variance. By the delta method, Var(Ŝ) ≈ Var(t_j − Sψ_j)/(N V²), t_j being the piece's per-row term; the stated error is sd(t_j − Ŝψ_j)/(√N V̂). Its bias is O(1/N), below the error's O(N^(−1/2)).
- *An effect cell* is e_j(g) = Σ_{v⊆u} (−1)^{|u|−|v|} f(g_v, x_{j,−v}). Under a product law, E[e_j] = Σ_v (−1)^{|u|−|v|} E[f | x_v = g_v] = f_u(g). The cell is their mean, unbiased, with stated error sd/√N. The sum is shift-invariant, and the base rows' variation outside u cancels wherever f is additive across u.
- *The cases where a quantity is undefined or flagged:*
  - The variance is zero when V̂ ≤ tolerance × the rounding scale, the mean of |y|·|y − ȳ|. The indices are then NaN, never 0 or 1. A variation below about the tolerance (10⁻¹⁰) times the mean's size is zero, and a smaller one above it is read. This keeps the exact route's margin over the rounding floor, about 10⁶ε.
  - A piece below zero by more than 3 standard errors plus the tolerance times that scale sets `BelowZeroBeyondError`. Nothing is clamped (decision 4).
  - The per-row deviations are scaled by their largest term. Without that, terms already in squared units would square again and leave the range beyond about 10^±77.

**Gates.** Unless a row says otherwise, every number comes from a run over the commits' files on a snapshot of `HEAD` (`5923ebc`), in `artifacts/chip-sensitivity/snap`.
- *The references against quadrature.* Every closed variance τ_u of every subset of the varying inputs agrees to at most 1.9·10⁻¹⁴ of V (Ishigami's worst; the ridges reach 2.5·10⁻¹⁵ and Friedman 1.3·10⁻¹⁴). Each finer check holds too:
  - conditional means at seeded points, to 4.2·10⁻¹⁵ of the scale;
  - the stated mains and totals against their Möbius sums, to 3.0·10⁻¹⁶ of V;
  - the main variances summing to V exactly (0, or 2.2·10⁻¹⁶);
  - f₀ + Σ_u f_u rebuilding f, with each main effect averaging to 0, to 1.3·10⁻¹⁵;
  - Si(π), Si(2π), Cin(π) and Cin(2π) against quadrature, within 6.7·10⁻¹⁶;
  - Ishigami's V against a²/8 + bπ⁴/5 + b²π⁸/18 + ½ (cited from memory), to 10⁻¹⁴.
- *No bias beyond the stated error, on chosen cases.* 300 replicates at N = 2048, with every input alone and every pair. Pieces whose truth and estimates are zero within rounding (inert inputs, additive pairs) are counted as exact zeros.

  | Reference | Pieces tested | Exact zeros | Worst bias / its error | Stated error / spread |
  | --- | --- | --- | --- | --- |
  | Ishigami | 16 | 4 | 1.44 | 0.933 to 1.035 |
  | g-function | 56 | 0 | 2.42 | 0.907 to 1.068 |
  | Friedman | 24 | 32 | 2.24 | 0.904 to 1.133 |
  | single index | 12 | 18 | 2.05 | 0.934 to 1.027 |
  | dense single index | 30 | 12 | 2.59 | 0.914 to 1.055 |
  | single index on a skewed box | 12 | 8 | 2.26 | 0.945 to 1.010 |
  | two ridges | 18 | 24 | 1.81 | 0.936 to 1.073 |

  The cuts are 4 and [0.85, 1.15]. The pieces tested include the g-function's small indices (S = 7.2·10⁻⁵) and Ishigami's S₃ = 0.
- *Effects, likewise* (200 replicates at N = 1024; every main effect and two pairs; 5-point grids). Worst bias 2.34 of its error; stated error / spread from 0.900 to 1.127 over 350 cells, with 75 exact zeros.
- *The error falls as N^(−1/2).* N runs from 256 to 16,384, with 100 replicates each.
  - *The spread:* slopes from −0.439 to −0.536 over Ishigami's and the g-function's 40 pieces, and from −0.468 to −0.530 over nine effect cells (cut ±0.08).
  - *The stated error:* slopes from −0.493 to −0.505 for the pieces and −0.498 to −0.502 for the cells (cut ±0.05).
- *Zero variance:*
  - Constants 0, 1, −3.7·10⁵ and 10²⁰⁰ read `ZeroVariance`, with NaN indices and exact zero pieces, errors and effects.
  - 10⁶(1 + 10⁻¹²x) reads zero, and 10⁶(1 + 10⁻⁸x) is read, its total index within its error of 1.
- *Shift invariance and the correction's evidence:*
  - f + c reads as f to 2.8·10⁻⁸ of the error at c = 10⁸ sd (3.0·10⁻¹⁴ at 10 sd), and the effects to 2.5·10⁻⁷.
  - The uncentred form's spread of V₁ on Ishigami is 1.44 times the centred one's at mean 0, 10.5 times at 10 sd and 94.9 times at 100 sd, as 2μ²T_i/N predicts.
- *Seeded sweep over far ranges* (300 cases at N = 4096):
  - *The cases:* g-functions with up to 8 inputs and a_i from 10⁻³ to 10⁴, so indices down to about 10⁻⁸ of V; Ishigami with a ∈ [0, 20] and b ∈ [0, ½], b = 0 a fifth of the time; and exponential ridge sums on boxes 10⁻³ to 10³ wide, centred out to ±10³.
  - *The far ranges:* each case is scaled by 10^±100 and shifted by 10⁻² to 10⁹ standard deviations, with up to three pairs.
  - *Estimator:* 5,274 z-scores and 610 exact zeros; mean −0.035, sd 1.003, 0.28% beyond 3, worst −4.07.
  - *Effects:* 200 cases with three random points per input, 861 z-scores; mean 0.065, sd 1.037, 0.12% beyond 3, worst −3.03.
  - *What the sweep found, before any of it ran clean:*
    - the deviations' overflow, fixed above;
    - a mean-square zero rule that would have refused a well-resolved function with a large mean, replaced by the rounding scale;
    - the ridge reference overflowing far from the origin, and losing small totals to cancellation, both fixed in the reference.
- *The values route against the function route:* bit for bit, for both estimators. The design's blocks, hybrids, cells and refusals are each tested.
- *The cost row* (the reduction alone, from values: wall clock per estimate, the machine at 96 to 99% load, so indicative). No fit or prediction is wired yet; the comparison against the exact route and prediction comes with slice 1. The probe is `artifacts/chip-sensitivity/probes/SmallIndexProbe.cs`, run from the snapshot and not committed.
  - p = 3, N = 4096 (20,480 values): 1.3 ms.
  - p = 6, N = 4096 (32,768 values): 6.5 ms; with all 15 pairs (94,208 values), 9.5 ms. Each allocates 161 KB.
  - p = 6, N = 65,536: 32 ms and 2.5 MB.
  - An effect of 100 cells: 15 ms; a 5 × 5 pair: 7.4 ms.
- *Fitting untouched.* The code is new, and nothing outside the slice calls it.
  - *Selected at commit scope:* test-select, run per test project from disk with the six files, selects Regression 3 of 339 classes (`PickFreezeTests`, `PinnedEffectsTests`, `SobolReferenceTests`). Oracle holds the two engine files in its closure and selects 0 of 39. No other project holds them.
  - *Regression at commit scope (with the pinned tier):* 93 of 93.
  - *The whole assembly:* Regression 3,585 of 3,585, in 22 minutes, at below-normal priority.

**`BayesPPR` 0.2.0** has no sensitivity readout. Read from the deparsed namespace (`artifacts/chip-bppr/oracle-source/`), its exports are `bppr`, `bppr_pca`, `bppr_resume`, `predict`, `plot`, `print` and `summary`, and none reads a variance decomposition or an effect. A BPPR fit's indices therefore have no oracle, and slice 1's directional checks rest on the closed forms above.

**Sources.**
- None of the six estimator papers is on disk; the coordinating session has the list. The forms above are derived here, and their citations stay from memory.
- On disk and read:
  - Roustant et al., 1612.03689. It states the derivative bound S_i^T ≤ C(μ_i)ν_i/D (Lamboni et al. 2013; Sobol' and Kucherenko 2009), which is decision 7's check, and gives the Poincaré constants on intervals.
  - Francom et al., 2506.11471, a review with no estimator formulas.
- Ishigami's and the g-function's closed forms are cited from memory and checked against quadrature above.

**Decisions taken beyond the brief.**
- *The home.* The estimators read any function, so they live in `src/maths/regression/sensitivity/` (`Maths.Regression.Sensitivity`), beside E8's generic `Maths.Regression.Subspaces`, rather than inside the BASS folder. Slice 1 calls them from there.
- *Saltelli's form, centred.* The uncentred form is not shift-invariant: its variance carries 2μ²T_u/N, 95 times the spread at a mean of 100 sd. Centring on f(B)'s own mean, with the covariance over N − 1, keeps it exactly unbiased at no cost. This is a correction, not an option.
- *The zero-variance rule.* It reads the values' rounding scale, the mean of |y|·|y − ȳ|, not their mean square, so a large mean alone does not make a function zero. Slice 1 passes values without the intercept, as the exact route's magnitudes leave it out, so the two routes judge zero alike.
- *The effect estimator* is built here, beside the indices, as one of "the estimators". Slice 2 wires it to fits.
- *Exact zeros in the gates.* A piece whose truth and every estimate are zero within rounding is counted apart rather than z-scored.

**Unbuilt.**
- *The route in fits:* slices 1 and 2.
- *Decision 3's points:* the engine's generator from the spec's seed, which is slice 1's, through E8's `Points`.

**Follow-ups.**
- **The main variance's form, for the user before slice 1 (it changes decision 3).** The two-difference form, ½·mean((f(B) − f(B_A^u))(f(A_B^u) − f(A))) after Owen 2013 (cited from memory), costs N(2p + 2) values against N(p + 2), and its error scales with the piece itself. From the probe (300 replicates at N = 4096):

  | Piece | S | Saltelli's relative error | Two-difference relative error | Gain at equal cost |
  | --- | --- | --- | --- | --- |
  | g-function, largest | 0.72 | 2.7% | 2.1% | 0.90 |
  | g-function | 0.18 | 6.9% | 2.7% | 3.8 |
  | g-function | 0.024 | 18% | 3.3% | 16 |
  | g-function | 0.0072 | 30% | 3.3% | 48 |
  | g-function, small | 7.2·10⁻⁵ | 320% | 3.2% | ≈5,000 |
  | Ishigami, S₁ | 0.31 | 5.4% | 4.0% | 1.2 |
  | Ishigami, S₂ | 0.44 | 3.1% | 1.75% | 2.0 |
  | Ishigami, S₃ | 0 | spread 0.20 | spread 0.099 | 2.6 |

  - At N = 4096, Saltelli's form cannot tell a small index from zero; the two-difference form reads it to 3%.
  - It loses at equal cost only on the dominant piece, and there by 10%.
  - *Recommended:* the two-difference form for main variances, with Jansen's totals, which share its blocks, unchanged.
  - *If taken,* it is one more block per input and a change inside `PickFreeze.Estimate`, before slice 1 binds. *If declined,* Saltelli's centred form stays, and small indices carry the errors above.
- *Low-discrepancy points* (decision 3's option). Nothing in slice 0 calls for them: the stated errors are calibrated to within 15%.

**Doc changes owed.**
- `docs/Maths/Regression/Sensitivity.md` (or the sensitivity page):
  - the sampled route's estimators and their derivations, as above: the main variance as a covariance with D_u, Jansen's total, the pooled variance, and the delta-method errors;
  - the correction to Saltelli's form, with the evidence;
  - the effect estimator;
  - the zero-variance rule's rounding scale;
  - the follow-up's numbers, once decided.
- `docs/Maths/Regression/Sources&Changes.md`: the uncentred Saltelli form's shift dependence, a source corrected; the six estimator papers, cited from memory.
- `docs/Maths/Regression/ActiveSubspaces.md`, or wherever E36's references are documented: the Sobol references beside them.
