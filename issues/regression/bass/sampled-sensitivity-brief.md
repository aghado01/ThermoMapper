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
