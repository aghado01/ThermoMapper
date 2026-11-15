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
