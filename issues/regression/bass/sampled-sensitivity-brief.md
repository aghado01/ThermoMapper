# Sobol indices and effects for terms without factors: E25's first build

Date: 2026-11-15. Status: **approved by the user on 2026-11-15** at its recommendations ("Yes, start it", to starting a chip on this brief at its recommendations). Drafted by the second coordinating session at the user's word that E36's dependencies come first, from [E25's entry](extensions-brief.md#sensitivity-under-dependent-inputs), the [lineage review](lineage-review.md#separability-is-what-the-closed-forms-rest-on)'s note on separability, [E36's brief](bppr-brief.md) (decision 10) and a reading of the sensitivity readouts' public surface. The estimators are cited from memory until slice 0 reads them. Reports go under "Reports" below, newest first.

## Why

- **A BPPR fit has no Sobol indices and no effects.** E36's directional hinge (`9e249dc`) reads several inputs through one projection, so it has no one-variable factors and the moment service refuses it. `Sensitivity`, `Effects` and the basis fits' versions refuse any fit holding one. Prediction, the gradient and E8's C already read such a fit.
- **The direction is settled.** The lineage review: a kernel without factors is sampled under any measure, with exact values per draw, so only the integral over x is approximate; Sobol, C and effects sit above both routes and know neither. E8 built that shape for C, two routes chosen per query. This brief gives Sobol and the effects the same shape.
- **Scope: the kernel side.** A product law, and terms without factors. E25's measure side (dependent inputs: accumulated local effects, Shapley effects, G4's mixture of products) stays in its entry for a later brief, and nothing here closes it off.

## Sources

- **In the repository:** E25's entry and [G4](lineage-gaps-register.md#g4-one-moment-service-across-dependent-inputs); E8's [build brief](active-subspace-brief.md), for the two routes, the shared points and the error rows; E36's brief and its slice 0 and slice 1 reports.
- **The estimators, read on 2026-11-15** (the report "The estimator papers, read against slices 0 to 3" below): Sobol' 2001 and Saltelli et al. 2010 for the pick-freeze forms; Jansen 1999 for the total index; Janon et al. 2014 for their asymptotic variances; Owen 2013 for small indices; Sobol' and Kucherenko 2009 for the derivative bound on a total index.
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
   - *Default decided 2026-11-15, pending the public readout change.* Use `MainVarianceForm.TwoDifference` for sampled main variances and expose the existing centred `MainVarianceForm.Saltelli` as an explicit choice on `SensitivitySpec`. Both use pick-freeze on two independent base sets A and B shared by every draw, and Jansen's estimator for total variances. First-order and total indices cost N(2p + 2) surface values per draw with the reverse hybrids and N(p + 2) with Saltelli. This choice applies wherever the Sobol readout takes the sampled route, including direct, named-output and basis fits. The exact route uses its analytic moments.
   - *Acceptance.* With the same seed and points, an omitted form must reproduce explicit `TwoDifference`, and explicit `Saltelli` must reproduce the previous sampled results. Run the existing closed-form and Monte Carlo error gates for both modes through the public readouts. The weak-index precision and equal-cost trade-off are measured in [slice 0's report](#slice-0-sources-references-and-the-estimators-2026-11-15).
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

### Decision 3: public sampled estimator choice — in progress (2026-11-15)

- Owner: the public sampled-estimator chat. Base snapshot `a83db0c`; other chats' structure-chain changes and document moves are excluded.
- `SensitivitySpec.MainVarianceForm` defaults to `TwoDifference`, carried through scalar/named-output, direct-functional and basis sampled designs. The primitive's default stays unchanged. The existing two-form public-readout gates are migrated off the temporary internal overrides; those overrides and their nullable plumbing will be removed in a separate cleanup commit.
- Baseline receipt: `artifacts/chip-sensitivity-public-form/baseline-results/baseline.trx`, 7/7. Public Saltelli tables at the base commit are fingerprinted for scalar, directional, named-output, integrated, conditioned, partially conditioned and basis views at matched seed/points/draws. The new pinned tests assert those fingerprints with explicit Saltelli.
- Scoped gates: **78/78** in `scoped-results/scoped.trx` (2 min 25 s), including both public forms against exact moments/quadrature and both primitive forms' closed references, far-range sweeps and small-index checks. The first capture run precedes the final pin assertions and expanded Auto checks; those are in the commit run. Calibration tables: `scoped-calibration.txt`.
- Running: the commit selection (264/359 Regression classes, plus pinned), under `artifacts/chip-sensitivity-public-form/`. Two processors, two test threads, single-node builds, below-normal MSBuild priority. The first parallel-node restore stalled and was cancelled; single-node restore/build proceeds, with NuGet's existing access-denied vulnerability-cache warning. Oracle selects 13/40 classes and will run separately.
- Next: final pins and route checks, commit selection and whole-assembly gates, cleanup of the superseded overrides, then the completed handoff here. Mixed-basis control variates remain deferred.

### The term-variables rule and its test verified from a clean archive (coordinating session, 2026-11-15)

- **Commit:** `5e99c83`, on a snapshot of `0802553`, through `scripts/verify-csharp.ps1` over its two changed files, in the two projects that select.
- **Result:** Regression 2,847 of 2,847 and Oracle 57 of 57 at commit scope, `SampledTermVariablesTests`' three cases among them.
- **Read.** The source change is the comment alone, which now states the rule as a completeness requirement. The test asserts, for each kind of term, that at least one input is left out, so no kind is tested vacuously.
- **Why it was asked for.** The BPPR chip's slice 5 found by one fixture that a category gating a ridge was never varied; the report shows this test failing on exactly the two gated kinds with that fix removed. The next term type whose reads are understated fails here.
- **Receipts:** `artifacts/verify-sens-vars-0802553/`, with the changed list, the counters and the snapshot tracked in [receipts](../receipts/verify-sens-vars-0802553/).

### A term's variables: the rule stated and tested (2026-11-15)

**Landed** in `5e99c83`, at the coordinating session's request after E36's slice 5 (`f11c8cc`) added a non-separable term's categorical factors to `TermVariablesOf`.
- *The comment* states the rule as a completeness requirement: every variable a term's value depends on. The route skips the block of an input no held term reads, giving it A's values bit for bit, so a variable left out reads zero. A non-separable term's categorical factors come from the term itself, since `ReadsVariable` answers for continuous variables only.
- *The test,* `SampledTermVariablesTests`. For each term and each input the rule leaves out, the term's column on A's rows equals its column on the route's own hybrid moving that input alone (`Stacked`), bit for bit. An output variable is moved alone at its locations.
  - *Hand terms* on the scalar fixture, at degrees 1 and 2, eleven kinds: a hinge, a hinge pair, a hinge with a categorical factor, a categorical term, Legendre alone and paired, a directional hinge on two inputs and on one, a gated directional hinge, a ridge column and a gated ridge column. Each leaves an input out (16 pairs per degree), and each column varies over A's rows, so it reads something.
  - *The fixtures' fitted terms:* the scalar fit (16 held terms, 21 pairs), the functional fit with t integrated (15, 20) and the two-output fit with both integrated (25, 46). The commit message's "61 terms" should read 56.
  - *Its power:* with slice 5's line removed in a snapshot, built into a fresh folder, the gated directional hinge and the gated ridge column fail at both degrees, and nothing else does.
  - *The basis route* (`BassBasisSampledSensitivity`, `BassBasisEffects`) calls the same function, so it needs no test of its own.
- *Gates,* on a snapshot of `d93c812` with the two files, every file touched, built into a fresh folder (`artifacts/chip-sensitivity/b-d93c812/`). The commit landed on `e96222c`, which changes no code since `d93c812`.
  - *Selected at commit scope,* per test project from disk: Regression 262 of 357 classes and Oracle 13 of 40. No other project holds the files.
  - *Regression at commit scope* (with the pinned tier): 2,847 of 2,847. *Oracle:* 57 of 57.
  - *This chip's seven classes in full:* 105 of 105.

### The estimator papers, read against slices 0 to 3 (2026-11-15)

The user supplied the six papers the brief cited from memory: Jansen 1999, Sobol' 2001, Sobol' and Kucherenko 2009 and Saltelli et al. 2010 in `supellex-private/PDFs/sensitivity/` (by DOI), and Owen 2013 (`1204.4763v1`) and Janon et al. (`1303.6451v1`, the preprint of the 2014 paper) in `supellex-private/gauntlet/stats/`. Each is read for its forms, their credits and their variances. Every form the engine computes is as the reports state it, and no code changes. Two credits are corrected below, and the remaining citations are now read.

**What the papers state, against the reports.**
- *Saltelli et al. 2010, Table 2.*
  - Formula (b), the paper's own, is mean(f(B)(f(A_B^(i)) − f(A))) on the triplet A, B, A_B^(i), where A_B^(i) is A with column i from B, as the engine's hybrid. This is the reports' "Saltelli's 2010 form".
  - Formula (f), credited to Jansen 1999, is ½·mean((f(A) − f(A_B^(i)))²), its best practice for a total. Together they cost N(k + 2), as decision 3 has them.
- *Jansen 1999, the alternating design.* The all-effects variance (the total) is estimated without bias by the mean of ½(Y_{i,1} − Y_{i,2})², the two values differing only in the group's draws. The engine's total is this form. Jansen's consecutive terms are correlated, so he estimates their mean's variance as a moving average's; the engine's rows are independent, so sd/√N stands.
- *Sobol' 2001.*
  - The total as half the mean squared difference (its (8)). Its Theorem 4 orders the variances: for a total the squared-difference form's is the smaller, for a main variance the product form's.
  - Its remark 3 (p. 277): the variances lose accuracy when f₀ is large, so subtract a crude c₀ ≈ f₀, after its reference [12].
- *Owen 2013.*
  - "Correlation 1", f(x)(f(x_u:y_−u) − f(y)), is credited to Saltelli 2002 and Mauntz 2002 and called correlated sampling by Sobol' and Myshetskaya 2007. It is Saltelli 2010's (b) with x = B and y = A.
  - Sobol' and Myshetskaya 2007 (read through Owen) found its variance grows with μ, and replaced f(x) by f(x) − c with c near μ: Owen's "Oracle 1" at c = μ.
  - "Correlation 2", (f(x) − f(z_u:x_−u))(f(x_u:y_−u) − f(y)), takes a third independent vector z. Its expectation is the closed variance τ_u (Owen's τ²_u), at four values per sample. On the g-function its efficiency over Correlation 1 is largest on the smallest index (4,256 in his Table 2).
- *Janon et al. 2014.*
  - S_N = Cov(Y, Y^X)/Var(Y) on one pick-freeze pair (after Homma and Saltelli 1996); T_N pools the mean and variance over Y and Y^X (after Monod et al. 2006).
  - Both are asymptotically normal when E Y⁴ < ∞ (their Proposition 2.2), with σ²_T = Var((Y − μ)(Y^X − μ) − ½S((Y − μ)² + (Y^X − μ)²))/V².
  - σ_T ≤ σ_S, with equality only at S = 0 or 1, and T_N is efficient among regular estimators that are functions of the pairs.
- *Sobol' and Kucherenko 2009.*
  - Theorem 2: S_i^T ≤ ν_i/(π²D), with ν_i the mean of (∂f/∂x_i)², for inputs uniform on the unit cube.
  - Theorem 3: σ_i²c²/D ≤ S_i^T ≤ σ_i²C²/D when c ≤ |∂f/∂x_i| ≤ C, under any law of variance σ_i².
  - For a normal input, S_i^T ≤ σ_i²ν_i/D.

**Corrections and precisions.**
- *The centring is credited, not new.* Slice 0 called the centred form "a correction, not an option".
  - The decision stands, but the idea is Sobol' 2001's c₀ and Sobol' and Myshetskaya 2007's c near μ. Saltelli 2010's (b) as printed is uncentred.
  - What the engine adds is f(B)'s sample mean in c's place, with the covariance over N − 1, which keeps the estimate exactly unbiased with no oracle.
  - The growth's closed form, 2μ²T_u/N, and its evidence (95 times the spread at a mean of 100 sd) stay as reported.
- *The two-difference form, precisely.* The reports' "after Owen 2013" reads: Owen's Correlation 2 with z's u-columns taken from A, so it needs no third base set.
  - In ½·mean((f(B) − f(B_A^u))(f(A_B^u) − f(A))), the cross term f(B_A^u)·f(A) shares a_u, so it adds τ_u where Correlation 2's f(z_u:x_−u)·f(y) adds nothing. The expectation is 2τ_u, hence the ½.
  - The cost is the same: N(2p + 2) for every input alone.
  - On an additive piece the per-row term is ½(g_i(b_i) − g_i(a_i))². For a normal piece of variance σ², its variance is 2σ⁴, against Correlation 2's 5σ⁴, so at equal cost the two-difference form is 2.5 times as efficient there. This is derived here and checked by Monte Carlo over 10⁶ rows (2.005 and 5.016). The two are not compared on interaction pieces.
  - Owen's family (his Theorem 1) needs z independent of y, so the form lies outside it. The form for |u| > 1, the product of the two Möbius differences over 2^|u|, estimates the main variance V_u and is the engine's own; Owen's forms estimate the closed τ_u.
  - Owen's Table 2 agrees in direction with slice 0's probe, the gain largest on the smallest index. The numbers do not compare: the g-functions and the baselines differ.
- *The index errors* are Janon et al.'s delta method in σ²_T's form: the piece's per-row term less S times the pooled per-row term ψ.
  - Two differences: ψ pools over A and B, which are independent, not over one pick-freeze pair; and the numerator is Cov(y^B, D_u), with f(A) in it. So Janon's efficiency result does not cover the engine's estimator, and the reports claim none.
  - The condition E f⁴ < ∞ holds for every law the engine takes: uniform laws, truncated mixtures and empirical laws on a bounded box, where a fit's surface is bounded.
- *Decision 7's bound.* Roustant et al.'s constant C(μ_i), 1/π² on [0, 1] and σ² for a normal input, is Sobol' and Kucherenko's Theorem 2 and their normal case; on [a, b], (b − a)²/π² by scaling.
- *The g-function's closed forms* are Saltelli 2010's Appendix A: V_i = (1/3)/(1 + a_i)², a closed variance Π_{j∈u}(1 + V_j) − 1, a total V_i·Π_{j≠i}(1 + V_j). The reference's values, checked against quadrature in slice 0, are these forms.
- *Ishigami's function.* Janon et al. state it with a = 7 and b = 0.1 on [−π, π]³, the reference's parameters, and its indices S₁ = 0.3139, S₂ = 0.4424 and S₃ = 0. None of the six states its variance. The form from memory, a²/8 + bπ⁴/5 + b²π⁸/18 + ½ = 13.845, gives their indices to four places (V₁ = 4.346, V₂ = 6.125), beside slice 0's quadrature check to 10⁻¹⁴.

**Candidate readings, recorded, not queued.**
- *Correlation 2 on interaction pieces.* It would enter decision 3's choice only by beating the two-difference form there by more than its 2.5 times loss on additive pieces. Not proposed.
- *The variance pooled over every block,* as T_N pools over its pair: each hybrid's values are draws of f under the law. The decision is the variance's estimator: adopt it if the indices' errors fall at the same N on the references; otherwise nothing changes.

**Doc changes owed** (refining slice 0's):
- `docs/Maths/Regression/Sensitivity.md`:
  - the forms' credits as read: Saltelli 2010's (b), and Jansen 1999's total, which Saltelli 2010 lists as (f);
  - the centring credited to Sobol' 2001 and to Sobol' and Myshetskaya 2007, with the engine's sample-mean form;
  - the two-difference form as Correlation 2 with z's u-columns from A, its 2τ_u, and 2σ⁴ against 5σ⁴ on an additive piece;
  - the index errors as Janon et al.'s σ²_T form pooled over A and B, and E f⁴ < ∞ on a bounded box.
- `docs/Maths/Regression/Sources&Changes.md`: the six papers read, with their paths, in place of slice 0's "cited from memory"; Saltelli 2010's (b) as printed, centred after Sobol' 2001 and Sobol' and Myshetskaya 2007, in place of "a source corrected".
- `docs/Maths/Regression/ActiveSubspaces.md`: decision 7's bound credited to Sobol' and Kucherenko 2009 (Theorem 2 and the normal case) beside Roustant et al.

### Slice 3 verified from a clean archive; slices 0 to 3 are all landed (coordinating session, 2026-11-15)

- **Commit:** `986814d`, on a snapshot of `e8446eb`, through `scripts/verify-csharp.ps1` over its nine changed files, in the two projects that select.
- **Result:** Regression 2,779 of 2,779 and Oracle 57 of 57 at commit scope, `SampledOutputsTests` 13, `SampledEffectsTests` 11 and `SampledSensitivityTests` 18 among them.
- **The two points the chip raised, checked.**
  - *The basis Sobol gate at 80 seeds.* At 24 seeds the two-difference form read 1.09% beyond 3 against the 1% cut; over 120 fresh seeds the same cells read 0%, with spread 1.026. More seeds under the same cut is the conservative answer, and it stands.
  - *The component accessors refuse on the sampled route,* with their reason: they are the exact route's quantities.
- **One reading asked for here and added** (`57d3cfa`): the basis effects' pooled spread of 0.882 passed only under slice 2's widened band. Over 120 fresh seeds every piece kind returns to about 1, spread over stated error 0.970 to 1.008 by median and the z-scores' spread 0.958 to 0.994. The stated errors are calibrated there, and the 0.882 was a shared fluctuation of correlated cells.
- **Where this brief stands.** Slices 0 to 3 are landed and verified, which is its "done". Slice 4, the split for mixed bases, waits on mixed bases being built.
  - With the user: decision 3's default form for a main variance; the gates' Bonferroni cut and wider spread band, let stand here.
- **Receipts:** `artifacts/verify-sens-s3-e8446eb/`, with the changed list, the counters and the snapshot tracked in [receipts](../receipts/verify-sens-s3-e8446eb/).

### Slice 3: basis fits, functional forms and several outputs (2026-11-15)

**Landed** in `986814d`, which carries its tests. Decision 8 carries the route and the error rows through each readout. The main-variance form stays the one choice in `PickFreezeDesign`, Saltelli's by default, and the gates read both.
- *Several outputs (E34):*
  - `BassMultiOutput.Sensitivity` takes the route as `BassSensitivity.Of` does. The declared output functional rides `Setup`'s projected coefficients, so the sampled route reads the projected surface as it reads a scalar fit's.
  - E34 has no effects readout, so there is nothing more to carry.
  - The chain fits a directional basis on several outputs (37 directional term reads in the gate's fit), and its Sobol route resolves sampled.
- *Basis fits (D07): Sobol.* `BassBasisSensitivity.Of` resolves the route over the components' held terms; `Exact` refuses the term it meets.
  - *A location's view:* the components' surfaces are read at the shared points, and each location's reconstructed surface Σₖ g_ℓk·wₖ is reduced by `PickFreeze`. The pieces are bilinear in the components, so this is the exact route's view read at the points, with exact errors.
  - *The integrated form:* the location is one more input, drawn uniformly over the training locations beside each row from a stream apart from the points' (`SeedTree.Derive(seed, 1)`), and y_m(ℓ) + Σₖ g_ℓk(wₖ + aₖ₀) is read at each row's own location. It is unbiased with exact errors, at some variance over integrating ℓ exactly.
  - *What changed in the record:* `BasisSobolAnalysis` gains a sampled constructor, `Route` and `Points`, and `BasisSobolSummary` gains `Route` and `Points`.
  - *What the sampled route refuses:* the ℓ-free component covariances (`ComponentCovariance`, `ComponentMainCovariance`, `ComponentTotalCovariance`, `ComponentMeans`) are the exact route's and refuse with their reason there.
- *Basis fits: effects.* The components are evaluated together at the pinned points, so each view Σₖ gₖ·eₖ carries its exact error √(gᵀΣg). Σ is the components' cell estimates' covariance, their per-row alternating sums' sample covariance over 2N.
  - *The core:* `EffectCore` serves r components; the direct route is its r = 1 case, and slice 2's tables read the same.
  - *The views:* `ComponentEffects`, `AtLocation` and `Integrated` (with the location's main effect and each input's two-way effect with the location) all carry errors. `BasisEffectAnalysis` gains `Route` and `Points`.
- *A direct functional fit's forms,* Sobol and effects. The integrated form draws the output variables with the points, from a stream apart (`SeedTree.Derive(seed, 2)`), and the blocks swap and pin them as they swap and pin inputs. The conditioned form holds them, or those `ConditionOn` names, at each coordinate, a decomposition or a set of effects per coordinate.
  - *A row's value:* Σⱼ aⱼ·Bⱼ(x_r)·hⱼ(t_r) at the row's own output location, each term's case and output factors read once per surface.
  - *Its use:* a direct functional response refuses a directional basis (E36), so on these forms the sampled route serves as a check on the exact route.
- *Shared across the readouts:*
  - a term's variables through the kernel contract (`TermVariablesOf`), so the default subsets and pairs come from a directional term's or a ridge column's active coordinates;
  - the effect request (`EffectRequest`) taking the measure and transform, not a direct fit's setup;
  - the decomposition tables (`SampledTables`).
- *Two of my own tests* asserted the slice-1 refusal of a functional fit on the sampled route. They now read its sampled analysis and form.

**Gates.** From a snapshot of `HEAD` (`418638a`) with the nine files. The bias cut is `BiasCut` (Bonferroni, 10⁻³); the pooled spread's band is [0.85, 1.15]; at most 1% may lie beyond 3.

| Readout | Fit, truth | Seeds, N | Cells (exact zeros) | Worst √R·\|mean z\| | Pooled sd | Per-cell sd | Beyond 3 | Worst z |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Several outputs, Saltelli | E34's `Fixed()`, two functionals, exact | 24, 2048 | 880 (0) | 2.85 | 1.009 | 0.64 to 1.34 | 0.30% | −3.48 |
| Several outputs, two-difference | the same | 24, 2048 | 880 (0) | 1.94 | 1.043 | 0.63 to 1.44 | 0.35% | −3.81 |
| Basis Sobol, Saltelli | D07's `FitOf()` (2 components), 5 locations and integrated, exact | 80, 1024 | 1,520 (0) | 3.46 | 0.999 | 0.86 to 1.17 | 0.32% | −3.66 |
| Basis Sobol, two-difference | the same | 80, 1024 | 1,520 (0) | 3.46 | 1.001 | 0.78 to 1.17 | 0.48% | −4.16 |
| Basis Sobol | directional components from the chain, 4 locations and integrated, quadrature of the predicted output | 30, 2048 | 195 (0) | 1.89 | 0.918 | 0.68 to 1.14 | 0.05% | −3.52 |
| Basis effects | `FitOf()`, components, 3 locations and integrated, exact | 24, 1024 | 9,060 (720) | 1.84 | 0.882 | 0.71 to 1.02 | 0.00% | 2.60 |
| Basis effects | directional components, 4 locations, quadrature | 30, 2048 | 432 (0) | 3.08 | 0.985 | 0.83 to 1.08 | 0.00% | −2.77 |
| Functional Sobol, Saltelli | integrated, conditioned on t (7 coordinates), and two outputs conditioned on t₀ with t₁ integrated, exact | 24, 1024 | 3,200 (40) | 3.28 | 1.030 | 0.60 to 1.32 | 0.77% | −3.79 |
| Functional Sobol, two-difference | the same | 24, 1024 | 3,200 (40) | 3.33 | 1.044 | 0.74 to 1.46 | 0.24% | 3.44 |
| Functional effects | the same three forms, 5-point grids, exact | 24, 1024 | 9,640 (0) | 2.96 | 0.943 | 0.58 to 1.39 | 0.07% | 3.74 |

- *The several-output identity, bit for bit.* On the directional fit from the chain and on `Fixed()`, the sampled analysis on a declared functional equals the scalar fit's whose draws carry the projected coefficients, at the same seed: main and total variances and their errors.
- *The basis Sobol gate's seed count.*
  - At 24 seeds and N = 2048, the two-difference form read a pooled sd of 1.134 with 1.09% beyond 3. A probe over the same views per piece kind found the location pair pieces at spread / RMS stated error 1.197.
  - Over 120 fresh seeds the same cells read 1.024, sd of z 1.026 and 0% beyond 3 (probe `artifacts/chip-sensitivity/probes/BasisCalibrationProbe.cs`).
  - The views are bilinear in two components, so their cells move together and the pooled statistics have about the seeds' few degrees of freedom. The gate reads 80 seeds on the same master.
  - *The basis effects' pooled sd of 0.882,* on the band's low side, read over 120 fresh seeds at the gate's views, grid and N (probe `artifacts/chip-sensitivity/probes/BasisEffectsCalibrationProbe.cs`, at the coordinating session's request). Each piece kind returns to about 1, so the stated error is calibrated, not conservative; the 0.882 was the same shared fluctuation as the Sobol cells'. Spread over RMS stated error, median (10% to 90%), with z's sd:
    - per component: mains 0.995 (0.872 to 1.021), z sd 0.973; pairs 0.970 (0.926 to 1.018), 0.974; means 1.008, 0.973;
    - per location: mains 0.986 (0.872 to 1.007), 0.958; pairs 0.970, 0.976; means 0.994, 0.979;
    - integrated: input mains 0.988, 0.964; pairs 0.970, 0.974; the location's main effect 1.001, 0.991; each input with the location 0.995 (0.882 to 1.008), 0.974; means 0.997, 0.994.
- *The quadrature's own error:* halving the panels moves the directional basis fit's integrated variance by 2.0·10⁻⁵ of itself.
- *The routes:* `Auto` on a separable basis fit is the exact route, bit for bit, for Sobol and effects, with null errors; a directional basis fit resolves sampled; `Exact` refuses its term, naming it.
- *The cost row* (the machine at 5 to 18% load, so indicative; probe `artifacts/chip-sensitivity/probes/OutputsCostProbe.cs`, not committed):

  | Per draw | Exact | Sampled, N = 4096 |
  | --- | --- | --- |
  | Basis Sobol, 2 components, every one of 15 locations' views | 0.116 ms | 13.3 ms, 2.6 MB |
  | Basis effects on 25 points, every location | 0.242 ms | 126 ms, 65 MB |
  | Functional Sobol, integrated | 0.055 ms | 1.42 ms, 0.35 MB |
  | Functional Sobol, conditioned on 7 coordinates | 0.016 ms | 5.32 ms, 1.69 MB |

- *Fitting untouched:*
  - *Selected at commit scope:* test-select over the nine files, per test project from disk, selects Regression 255 of 350 classes and Oracle 13 of 40. No other project holds them.
  - *Regression at commit scope (with the pinned tier):* 2,779 of 2,779.
  - *Oracle at commit scope:* 57 of 57.
  - *The whole assembly:* Regression 3,695 of 3,695, in 19 minutes at below-normal priority. `HEAD` did not move between the snapshot and the commit, and `output/BassComponents.cs`, which E8's slice 8c will edit, is unchanged.

**Decisions taken beyond the brief.**
- *The basis integrated form samples the location* beside each row, rather than integrating it exactly over an estimated component state. It is simpler and unbiased with exact errors; the alternative is the follow-up below.
- *The sampled basis analysis keeps no component covariances.* The views read the reconstructed surfaces directly, and the component accessors refuse there with their reason.
- *Basis effects evaluate the components together,* so their views' errors are exact rather than bounded. `EffectCore` serves r components, and the direct route is r = 1.
- *Separate streams for the integrated output variables and the basis locations,* derived from the spec's seed by `SeedTree`.
- *The basis Sobol gate at 80 seeds,* for the reason above.

**Unbuilt.**
- *The sampled route's component covariances for a basis fit.* They are estimable as symmetrized cross pick-freeze pieces, at r(r + 1)/2 pieces per block.
- *An effects readout for several outputs (E34),* which the engine does not have on either route.
- *Slice 4 (the split f_sep + f_ns),* which waits on mixed bases (E28) and is not this chip's to start.

**Follow-ups.**
- *The basis integrated form, exact in ℓ.* A family estimator could weigh every location at every row, taking the per-row terms of each location's pieces, at about m times the reduction. Decision: whether the variance that sampling ℓ adds matters. The gate's cells show no bias; a variance comparison at equal cost would decide.
- *Basis effects' cost* (126 ms per draw at 25 points, N = 4096). Reusing the stacked cases across batches and chunks would cut its allocation churn. Decide against the cost row if long basis posteriors' effects are read often.

**Doc changes owed.**
- `docs/Maths/Regression/Sensitivity.md`:
  - decision 8 per readout: the projected functional, the basis views and their bilinearity, the integrated form's sampled location, the basis effects' √(gᵀΣg), and the functional forms' row evaluation and output-variable draws;
  - the streams derived from the seed.
- `docs/Maths/Regression/BASS.md` (basis fits): the component accessors exact-route only.

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
