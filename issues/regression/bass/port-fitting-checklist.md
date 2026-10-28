# BASS port: fitting: checklist

Companion to [port-fitting-brief.md](port-fitting-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Decisions (stage 1)

- [x] **F-D01** Term representation and prior measure: decided 2026-10-27 as drafted; R's categorical subset measure by default, the paper's selectable.
- [x] **F-D02** Count prior at the capacity limit: decided 2026-10-27 as drafted (global truncation, λ integrated out, h₁ = h₂ = 10).
- [x] **F-D03** Coefficient-prior modes and the τ hyperprior: decided 2026-10-27 with proper defaults, σ² scaled to the response and the intercept level-centred (F-O9); Jeffreys σ², Liang and the zero centre selectable.
- [x] **F-D04** Collapsed and conditional update order: decided 2026-10-27 as drafted, with k = M + 1 under level centring.
- [x] **F-D05** Tempering and SMC target paths: decided 2026-10-27; the integrated form by default, the linear form selectable.
- [x] **F-D06** Numerical model admission: decided 2026-10-27; the count as the support default, the centred angles under every convention.
- [x] **F-C1** Contract shapes: data, config and result, decided 2026-10-28 as recommended.
- [x] **F-C2** The shared-code changes the contracts need, before stage 4's first fit: `Sigma2Prior` as a serializable record, shared with BARS's config; `GPriorCentre.Flat` for the Liang convention, refused by BARS until the shared block supports it.
- [ ] **F-R1** Deterministic-simulator profile: drafted 2026-10-26 as a proposal; its coefficient prior is set by a pre-registered comparison.
- [ ] **F-M1** Move schedule: equal over available moves (D01) against the prior-cancelling option, in a pre-registered comparison.

### Reference fixtures (stage 1)

- [x] **F-O1** Propriety of the collapsed marginal over τ: `r/oracles/bass_marginal_oracle.R`.
- [x] **F-O2** Order-free admission and R's `getQf` rule: `r/oracles/bass_admission_oracle.R`.
- [x] **F-O3** Count and order laws under a constant likelihood, by enumeration (D01, D02): `bass_enumeration_oracle.R`'s `laws` mode. Checks: D01's kernel and the prior-cancelling schedule against the target's state, count and order laws; the order tilt at M = 1 against the admitted prior mass; λ at retention against the prior tilted by F.
- [x] **F-O4** Boundary ratios π(1)/π(0) and π(M_max)/π(M_max − 1), and R's fixed thirds as a defect fixture: the same mode. Checks: both schedules exact; R's thirds 3× at the bottom, 2/3 at the top, right in the middle.
- [x] **F-O5** Wallenius probabilities against brute-force sums over draw orders: `r/oracles/bass_wallenius_oracle.R`.
- [x] **F-O6** Defect fixtures: R's tempering and swap density, R's Jeffreys mode, KHAOS's double count factor: `r/oracles/bass_enumeration_oracle.R` (modes `swap` and `count`) and `r/oracles/bass_jeffreys_oracle.R`.
- [x] **F-O9** Intercept conventions: Liang's flat intercept against level centring, with R's zero centre beside them, `r/oracles/bass_intercept_oracle.R`. Checks: on F-O8's grid, level centring's extra credit against Liang (positive, 0.4–0.5 nats for one term, 0.2–0.3 between structures one term apart, flat in N); on an explicit 40-point design, the level marginal against the g-form identity, and the three conventions' Bayes factors at offsets 0, 10 and −250.
- [x] **F-O8** τ's prior: the unit-information Gamma family, Zellner–Siow against Francom et al. 2018's shape 1/N, `r/oracles/bass_tau_prior_oracle.R`. Checks: the log-moments by quadrature against ψ and ψ′; Zellner–Siow's offset from unit information is the same at every N; the null-exit penalty and its growth; the near-cancellation between structures that both carry τ.
- [x] **F-O7** Centred admission angles for the Liang default, as the `centred` and `hadamard` modes of `r/oracles/bass_admission_oracle.R`. Checks: uncentred s_j = centred s_j · ρ_j; the downdate's accuracy on a near-constant column; and the centred Hadamard form against explicit centring.

### Build (stage 4)

Sequenced in the fitting segment's [stage-4 plan](port-fitting-brief.md#stage-4-plan-decided-2026-10-28), decided 2026-10-28: steps S1–S10, with three shared additions (Wallenius in Distributions, a per-replica block in `ParallelTempering`, the flat-intercept form in `CollapsedGPrior`).

- [x] **F1** Term representation: continuous, categorical and functional parts, with interned terms (S2).
- [x] **F2** Design layer and the direct functional contractions: the terms' factors (S2) and `BassDesign` (S4).
- [ ] **F3** Collapsed g-prior target under both intercept conventions, with the fixed-g mode.
- [ ] **F4** Birth, death and change moves with exact boundary probabilities.
- [ ] **F5** D04's sweep, with composition draws at retention.
- [ ] **F6** Parallel tempering on the likelihood path, with the exact reference rung (D03's proper defaults supply it).
- [x] **F7** D06's two gates and refusal reasons (S4).
- [ ] **F8** Retention, checkpoints and exact restart.
- [ ] **F9** Pointwise density producer on the OQ17 seam, with its chains for PSIS's relative efficiency.
- [ ] **F10** Chains start on the support.
- [ ] **F11** Direct fits report the residual-dependence check along the output coordinates (proposal).

### Room for extensions

- [ ] **F-X1** Observation-model seam carries an offset with the weights.
- [ ] **F-X2** The sweep has a latent slot.
- [ ] **F-X3** Term types behind one design-kernel contract, with derivatives.
- [ ] **F-X4** Proposal kernels behind one contract with exact forward and reverse probabilities.
- [ ] **F-X5** Coefficient priors as a declared set.
- [ ] **F-X6** The output-side inner product is declared (Q = I now).

## Reports

### S4: the design layer and admission — 2026-10-28

`BassDesign` (`src/maths/regression/spline/bass/BassDesign.cs`), one per chain, with `BassAdmission` (the rules, resolved from the config), `BassRefusal` (cap, support, numerical, variance inflation) and `DesignProposal` (a candidate's account).
- **The centred basis.** The design holds Z_c and ỹ = y − ȳ under every intercept convention (D03's computation): scalar fits keep explicitly centred columns, and direct fits keep the case-side and output-side factors with Gram entries from the centred Hadamard form. Every entry and inner product is formed without cancellation and summed pairwise (N10). Its statistics are ‖r‖² and Q_c = b_cᵀG_c⁻¹b_c; each convention adds its own pieces in S5.
- **Moves.** Births, deaths and changes go through `BorderedCholesky`'s private proposals; a change's new column goes last. `Load` builds a structure wholesale and strictly, for starts and restores (S7).
- **Admission (D06),** in order: the cap; a repeated term (a duplicated column, s = 0); the support threshold per part, by count or, as the option, by mass; the order-free gate, every leave-one-out angle at least s_min; the declared variance inflation cap. A death needs no gate. When the bordered factor fails the gate or comes within a factor of two of s_min, a strict factorization of the candidate's Gram decides, and an admitted candidate is accepted from it.
- **Statistics.** A candidate's ‖r′‖² updates the current residual's: a birth removes (zᵀr)²/δ, a death adds β̂_m²/(G⁻¹)_mm, and a change does both through the reduced model, with zᵀr_reduced = zᵀr + β̂_m·uᵀG⁻¹e_m/(G⁻¹)_mm. These lose digits only when one move explains nearly all the residual (below 10⁻⁸ of it), where the candidate's residual is formed explicitly. Taking ‖ỹ‖² − Q′ instead would lose them on every good fit. On each acceptance β̂ is solved and the residual recomputed explicitly, r = ỹ − Z_cβ̂, so it cannot drift. That replaces D06's update-then-cadence wording with something at least as accurate and no costlier, since acceptances are rare next to proposals. The factor still rebuilds strictly on its cadence.
- **Tests.**
  - `BassDesignTests` (10) check against a reference that shares nothing with the bordered factor: the explicitly centred design, its least squares and its leave-one-out angles from a thin SVD.
    - Random sequences of 400 proposals on scalar and direct data, at s_min = 10⁻¹⁰ and at a binding 3·10⁻³: about 300 proposals compared per run, every decision equal (except within 5% of the floor), ‖r′‖² and Q′ within 10⁻⁸, the angles within 10⁻⁶, and each accepted state's residual, ‖r‖² and Q.
    - The binding floor refused 15 and 6 proposals on the angle gate and admitted 24 within a factor of two of it, where the strict factorization decided.
    - F-O2's order-freedom with three hinges at neighbouring knots: the smallest angle is 5.5·10⁻⁵; h₂ after h₁ shows 1.6·10⁻³, and h₃ added last 2.2·10⁻⁴. At a floor of 1.1·10⁻⁴ every one of the six orders refuses the third birth, where a newest-column rule would admit the set in four of them.
    - One case per refusal reason, including a categorical subset and its complement, which span one space with the intercept (D01's complement pairs).
    - Column mass as the option; `Load` against the incremental path; a near-exact fit through the explicit path (‖r′‖² ≤ 10⁻²⁴‖ỹ‖²).
    - 1,500 accepted moves on direct data: the residual within 10⁻⁹, then a rebuild finding a factor drift of 3·10⁻¹⁵ at M = 11.
  - `BassAdmissionParityTests` (3, `Maths.Oracle.Tests`): the design's angles against `bass_admission_oracle.R`'s QR of the explicitly centred design, within 10⁻⁶ relative (F-O7): random scalar terms to 5·10⁻¹²; direct terms by the three-term form against the materialized design to 3·10⁻¹³; and F-O7's hard case to 3·10⁻⁹ at an angle of 10⁻⁶. That case is an indicator true for all but one of 10⁵ cases, nearly collinear with a hinge once centred.
  - Suites: Regression 471/471; Oracle 22/22, two TDA transport tests passing on a rerun after timing out under the chips' load (their R oracle runs near the bridge's 120-second limit alone).
- **For S10's cost baseline.** `KhatriRao.CentredBorder` re-centres every current factor on each call, so a direct proposal allocates O(M(n + m)). That is correct, but a candidate for keeping the centred factors in the design.

### S2: terms and the basis — 2026-10-28

In `src/maths/regression/spline/bass/`:
- **Terms** (`BassTerm`). A continuous, a categorical and a functional part, each of distinct variables, with at least one factor. Terms are canonical, each part sorted by variable, and compare by value, so a model's terms form a set (D01's M! count). `HingeFactor(Variable, Sign, Knot)`; `CategoricalFactor` holds its subset as sorted level indices, so a factor of any number of levels is one value (the hub's edge case of more than 64 levels). `BassTermTable` interns terms in the order first seen, one writer at a time; a run keeps one per chain and merges them in chain order (S7).
- **Inputs** (`BassInputs`). The scaled columns through the transform; each variable's knots, its distinct scaled values (R's `unique` on the scaled column); the categorical training-level indices, −1 outside them; each predictor's place in the data. `CaseInputs` and `OutputInputs` also carry prediction's new inputs. The S1 review's note to take knots from the raw values is withdrawn: a hinge reads only the scaled knot, so raw values that scale to one double are one knot location, and counting them twice would bend D01's knot measure.
- **The basis** (`BassBasis`), R's construction. The case side is Π hinges^q / c times the indicators, with c the hinge product's peak on [0, 1], and c = 1 where that is 0 (R's guard for a hinge facing its end of the range); the output side is the functional hinges, normalized the same way. Powers follow R's `^`: x·x at 2, the library power above. `Validate` refuses an index outside the inputs, a subset holding every level, and a functional factor on a scalar response. `Support` counts each part's nonzero values and mass.
- **Tests.**
  - `BassBasisTests` (10): canonical terms; the table; knots and level indices; hinges by hand at degrees 1 and 2; the zero-constant guard; interactions peaking at one; categorical indicators with an unseen level; support per part; and at degrees 1–3, a direct functional design's separable form against an explicit nm × k design from an independent evaluator, with the Khatri–Rao Gram against the explicit one.
  - `BassBasisParityTests` (3, `Maths.Oracle.Tests`): 60 random terms mixing every part, at degrees 1, 2 and 3, on tied predictors, against R BASS's own `makeBasis`, `makeBasisCat` and `scale_range` at the pinned commit (`r/oracles/bass_basis_oracle.R`, new): within 10⁻¹⁴, and the scaled columns within 10⁻¹⁵. A hand-checked two-term fixture confirms the oracle runs and gives the expected values.
  - Suites: Regression 461/461, Oracle 19/19.

### S1: the contracts — 2026-10-28

In `src/maths/regression/spline/bass/` (namespace `Maths.Regression.Spline.Bass`):
- **`BassData`** with `ContinuousPredictor`, `CategoricalPredictor`, `ScalarResponse` and `FunctionalResponse`. Each copies its values and refuses its defects with a reason: no predictor, a name used twice (output variables included), mismatched lengths, a constant continuous column (no knot), a non-finite value, a code outside the levels, a repeated level, a factor taking fewer than two levels.
- **`BassInputTransform`**: R's range scaling to [0, 1] with exact ends and linear beyond, `Contains` for prediction's flag, and each factor's training levels (the declared levels some case takes, by name).
- **`BassConfig`** and its records, all values, with the decided defaults. `ToJson` and `FromJson` use canonical options (enums by name), and `Validate` names the field in every refusal; a manifest is checked on load. `SampledG`/`FixedG` and `EqualOverAvailable`/`PriorCancelling` are polymorphic records with a `Kind` discriminator.
- **`Bass.Fingerprint`**: `RunFingerprint` over the canonical JSON and the data's exact bits.
- **`BassResolved`**: N = n_x·n_r, the basis cap min(1000, N − 1), the interaction caps by the predictors of each kind, the support thresholds min(20, 0.1·n) per part, the centre's level, s², σ²'s prior as resolved, τ's prior Gamma(a, rate a·N) or the fixed τ (1/N by default), and the automation log of every auto-pick with its reason. It refuses a zero response under the zero centre, a constant one under level centring or the flat intercept, and a support threshold no term could meet.

The refinements to the decided sketch are recorded under the fitting segment's contract choices ("Refined in S1"). Tests: `BassContractTests` (12), each refusal, the copies, the transform's edges, JSON round trips by value, the fingerprint's sensitivity to a setting and to one ulp of data, and the resolved values against direct computation. Regression suite 450/450 (438 before, 12 new).

### Stage 4 planned (a proposal) — 2026-10-28

After the baseline was requalified at `cf37584`, the fitting segment gained a stage-4 plan: F1–F11 in ten steps, with what each reuses, the three shared additions it needs and the fixtures that check each. Four choices are put to the user: the ridge priors deferred, τ under tempering through a per-replica block, starts from exact prior draws, and where the work runs.

### F-C2: the shared-code changes — 2026-10-28

- **`Sigma2Prior` is a value.** A sealed record of `Sigma2PriorKind` (`ScaledToResponse`, `Jeffreys`, `InverseGamma`) and that kind's parameters, the kind-plus-parameters shape CONTRIBUTING gives sub-configs. A parameter the kind does not take is null, not NaN, so the JSON conventions that refuse named floats read it. The constructor checks the parameters against the kind, so a prior read from JSON is checked as a built one is. `IsProper` and `ReadsResponse` are read from the kind and not serialized. The factories and `ToString` are unchanged, so BARS's call sites and its checkpoint fingerprint are too. Its first `System.Text.Json` attribute (`JsonIgnore`) enters the maths layer; BASS's config unions will need the same namespace.
- **`GPriorCentre.Flat`,** appended after `Zero` so the existing values keep their numbers. BARS refuses it at `WeightedNormalModel`'s constructor, both `BandedGPrior.Fit` overloads and `BandedGPrior.Level`, before any work, with the reason: B-splines do not split the constant out of the basis, and the block has no flat-intercept form yet.
- **Tests:** `Sigma2PriorTests` (JSON round trip under plain and snake-case, string-enum options; value equality; refusals built and read back; the fingerprint's spelling; resolution) and `NormalCentreTests.FlatCentre_IsRefused_WithItsReason`. Regression suite 438/438 (432 before, 6 new). Samplers does not reference Regression and was not rerun.

### Contract shapes decided (F-C1); stage 1 complete — 2026-10-28

The user took the five recommendations: values throughout; `GPriorCentre` gains `Flat`; `Sigma2Prior` becomes a shared serializable record; R's 9,000 and 1,000 as provisional sampling defaults; unseen levels refused by default. The two shared-code changes are F-C2, before stage 4's first fit. With D01–D06 decided and every stage-1 fixture built, stage 1 is complete; the baseline rerun opens stage 4.

### Contract shapes drafted (F-C1) — 2026-10-28

At the user's request, stage 1's last item: the data, config and result contracts for fitting, in the segment's "Contract shapes" section, consuming D01–D06 as decided.
- *Data:* named predictor columns, continuous and categorical; scalar and functional responses, with output variables for the coordinates and the covariates known at them; refusals at entry; one input transform fitted from the training data.
- *Config:* values only, so it serializes and fingerprints whole, unlike `BarsConfig`'s instances (OQ10); the decisions' closed sets as records; R's defaults where the decisions leave a number open, the sampling lengths provisional.
- *Result:* every draw's generative state in per-chain draw stores over an interned term table, SMC weights, diagnostics including D01's and D02's tilts and D06's refusals by reason, the evidence with its priors, the checkpoint beside it; fitted values, the pointwise density, prediction and sensitivity as views.
- *Open for the user:* values throughout; extending `GPriorCentre` with `Flat`; `Sigma2Prior` as a shared serializable record; the sampling defaults; the unseen-level policy.

### Count, order and λ laws and the boundary ratios (F-O3, F-O4) — 2026-10-27

`bass_enumeration_oracle.R` gains a `laws` mode for D01 and D02 as decided: a constant likelihood; F = the support rule and D06's centred angle gate; the count prior with λ integrated out (h₁ = h₂ = 10), so a birth from M carries (M + h₁)/(h₂ + 1). Three kernels, each from its own transition matrix: D01's equal-over-available move probabilities; the prior-cancelling schedule (c = 0.4), the option kept for F-M1; and D01's selection with R's fixed thirds in the ratios. Run by `bass_stage1_checks.R`; every check passes, F-O6's included, since the refactor keeps the earlier modes' defaults.
- *Fixtures.* p = 3 on {0, 1}³, M_max = 3, npart = 2: 60 terms, 18 supported, 913 states in F, 75 refused by the gate (the Wallenius proposal active, and a middle state). p = 2 on a four-level grid, M_max = 2, npart = 3: 80 terms, 36 supported, 665 states, 2 refused.
- *F-O3.* D01's kernel and the cancelling schedule match the target's state, count and order laws to 10⁻¹² (TV). The tilt is visible: at M = 1 the order law is 2/3 and 1/3 against the unconstrained ½ and ½, equal to the admitted prior mass per order (0.50 and 0.25 on p = 3; 0.75 and 0.375 on p = 2). λ at retention, Gamma(h₁ + M, h₂ + 1) mixed over the count law, equals the prior tilted by F to 2·10⁻¹⁴: D02's "exact under global truncation".
- *F-O4.* Both schedules give the boundary ratios exactly. R's thirds give π(1)/π(0) three times the target's and π(M_max)/π(M_max − 1) two-thirds of it, and are right in the middle (π(2)/π(1) at M_max = 3); their state law is 0.25–0.26 off in total variation.
- *The cancelling schedule's move probabilities* are proper at every M.
- **Stage 1 now has the contract shapes left.**

### The design pass: D01–D06 decided — 2026-10-27

The user's design pass took every recommendation, and adopted D10's adaptation note.
- **D03, proper defaults.** σ² ~ IG(10⁻³, 10⁻³·s²) scaled to the response, and the intercept level-centred on the mean, so every default prior is proper: evidence by default, exact fits proper, the reference rung and SMC from the prior on either form, OQ19 settled for BASS. Jeffreys σ², the Liang convention and R's zero centre stay selectable. The cost is about 0.23 nats of extra credit per added term against Liang (F-O9).
- **A correction made in the pass.** The marginal depends only on the column space of [1, Z] under every intercept convention, so the centred basis [1, Z_c] computes it exactly. The alternative's knock-on for D06 (a stricter gate on R's uncentred Gram) and its claim to spare the centred downdate were wrong: D06 reads the centred angles under every convention and no longer depends on the intercept choice.
- **D01:** as drafted; R's categorical subset measure by default, the paper's uniform selectable, recorded in the register as a divergence between paper and code. **D02, D04:** as drafted (k = M + 1 under level centring). **D05:** the integrated form by default, the linear form selectable. **D06:** the count as the support default, column mass an option.
- **Applied:** D01–D06 and the regime note in this segment, the calibration brief's SMC paragraph, D10's note in the numerics segment, the hub's status, decision table and register.
- **What remains of stage 1:** F-O3 and F-O4 (now unblocked), and the contract shapes. A baseline rerun then opens stage 4.

### Intercept conventions fixture (F-O9) — 2026-10-27

At the user's request, so the D03 alternative's numbers reproduce from the repository: `r/oracles/bass_intercept_oracle.R`, run by `bass_stage1_checks.R` (every check passes).
- *Grid mode.* Liang's and level centring's Bayes factors against the intercept-only structure on F-O8's grid, converged over twice the range to 10⁻⁸. Level centring's extra credit is 0.43–0.47 nats for one term and 0.68–0.69 for two, and between structures one term apart 0.230, 0.237 and 0.241 at N = 40, 100 and 1000.
- *Design mode,* on a 40-point design with one and two hinges:
  - the explicit level marginal (R's τ form on y − ȳ over the uncentred Gram) equals the g-form identity, N in place of N − 1 in one exponent, to 10⁻⁸;
  - level centring (9.023 and 8.065 nats) and Liang (8.559 and 7.369) do not move at offsets 0, 10 and −250;
  - R's zero centre does: the two-hinge structure's Bayes factor goes from 5.56 to 2.16 and −3.87 nats. This also serves the shift fixture the current D03 draft lists, read on Bayes factors rather than on τ's posterior.
- *A numerical note.* The first run failed in `integrate()` on the zero centre at offset −250, where y'y and Q agree to six digits; the oracle now takes the residual directly, as D03 warns.

### D03's proper-defaults alternative — 2026-10-27

At the user's request after the readiness review, drafted beside D03 for the design pass, not decided: `Sigma2Prior.ScaledToResponse()` for σ² and level centring (`GPriorCentre.Mean`) for the intercept, so that every default prior is proper. The Liang convention, the zero centre and Jeffreys σ² stay selectable.
- *The cost, on F-O8's grid (fixture F-O9 since).* Level centring's marginal is Liang's with N in place of N − 1 in one exponent. Against Liang it favours a structure over the intercept-only one by 0.43–0.47 nats at M = 1 and 1.06–1.08 at M = 5, and two structures one term apart by 0.23–0.24, at N = 40, 100 and 1000.
- *What it buys:* evidence by default, proper exact fits, the reference rung and SMC from the prior without a separate mode, OQ19 settled for BASS, and R's uncentred numerics.
- *Knock-on edits* to D03–D06, the calibration brief and the regime note are listed in the alternative and not applied.
- *Found by the review, independent of the choice.* The calibration brief says the improper default cannot start a prior-to-posterior SMC run, while D05's linear form starts one from the prior over (structure, τ), which is proper under the improper σ² and intercept priors. The alternative removes the paragraph; under the current draft it needs the linear-form exception.

### D05's evidence inputs from the shared infrastructure — 2026-10-27

At the user's request, the shared-infrastructure session's four points are folded into the D05 draft's evidence bullet as inputs for the design pass, not decisions: the estimator follows the tempered form; the linear form makes path sampling free; BASS's priors record adds τ's prior, the term measure, M's count prior and the interaction cap; and SMC's error under-covers at one sweep per stage. Their evidence is in the numerics checklist's reports ("The linear form's path sampling, and the evidence's priors"; "Follow-ups"). D03's evidence bullet now points at D05's rather than naming path sampling.

### τ's prior corrected (F-O8) — 2026-10-27

The user's decision: every source is checked, papers included, and a design shown flawed is corrected (the hub's validation section; the register gains a **Corrected** outcome). Francom et al. 2018's Gamma(1/N, 1) is the first entry. D03 now states the τ prior as the unit-information Gamma family Gamma(a, rate a·n), default a = ½ (Zellner–Siow), with the shape a declared constant; Francom's prior is the member a = 1/N.

`bass_tau_prior_oracle.R`, run by `bass_stage1_checks.R` (all checks pass):
- *Log-moments.* Quadrature matches ψ and ψ′ to 10⁻⁸. E[ln τ] − ln(1/N) is −1.27 at every N under Zellner–Siow, and −36.9, −96.0 and −993.7 at N = 40, 100 and 1000 under Francom's shape; sd(ln τ) is 2.22 against 40, 100 and 1000. Francom's median is e⁻²⁸, e⁻⁷⁰ and e⁻⁶⁹⁴.
- *Leaving the intercept-only structure* (Liang convention, one term): Francom's prior costs 1.37–1.96 nats more at N = 40, 2.05–2.82 at N = 100 and 4.24–5.08 at N = 1000, over R² from 0.05 to 0.5.
- *Between two structures that both carry τ* (two terms at R² = 0.25 against one at 0.20): the priors differ by 0.17–0.22 nats.
- *Limits.* The Bayes factors are quadrature only, checked for convergence over twice the range; the cases are a grid, not a study. R's convention, where every structure carries τ, is covered only through the near-cancellation case.

### The parity rule — 2026-10-27

The user's decision: the C# preserves no R defect to align with a fixture (the hub's validation section). The `BassLiteral` options are gone from D03 (R's τ rate 2/n) and D06 (R's `getQf` admission), and from the declared set of coefficient priors. R's zero-centred intercept stays selectable as DiMatteo, Genovese and Kass's prior, `GPriorCentre.Zero`, not as a parity option. Earlier reports that mention `BassLiteral` record the drafts as they stood.

### F-O6 at three predictors, and reversibility — 2026-10-27

From the reference-fixtures chip's follow-up, `f7acf2a`. The enumeration oracle now reports each kernel's detailed balance with its own stationary law. The checks add p = 3 to the swap and count fixtures, so the Wallenius factor of R's proposal varies with the state. All checks pass.
- **Swap at p = 3.** N = 27 on a 3 × 3 × 3 grid, a support rule of 4 nonzero values: 126 terms, 48 supported, 1174 states, 14,160 moving transitions.
  - *Validity.* R's `getQf` agrees with the strict rank on every state, and R's `alpha` is b times its b = 1 value. At b = 1 R's kernel, with the state's move probabilities, is the posterior to 7·10⁻¹⁵. With its fixed thirds it is 0.100 off in total variation (F-O4's defect, larger here than at p = 2).
  - *Targets.* At b = 0.5 and 0.2, R's kernel is 0.077 and 0.143 in total variation from D05's likelihood path, 0.055 and 0.081 from the prior times p_1^b, and 0.41 to 0.65 from the two powered paths.
  - *Swaps.* R's swap misses balance against its own rung laws by spreads of 1.88 and 1.14, and against the candidate families by 0.88 to 6.18. D05's kernels reach their target to 7·10⁻¹⁶, and its swap balances on them to 1·10⁻¹⁴.
- **R's tempered kernel is not reversible at p = 3.** Its gap from detailed balance with its own stationary law is 5·10⁻⁷ at b = 0.5 and 1·10⁻³ at b = 0.2, against 2·10⁻¹⁶ at b = 1 and at most 7·10⁻¹⁶ for D05's kernels. At p = 2, where every interaction draws the one pair, R's kernel stays reversible at every rung (at most 5·10⁻¹⁶).
  - *Why.* R's acceptance is the untempered ratio to the power b, so the proposal ratio enters as (q_rev/q_fwd)^b. The kernel then has a detailed-balance target only if that factor's (1 − b) power is a ratio of state functions. At p = 3 the Wallenius factor's is not.
  - *So below b = 1* R's within-rung moves balance against no law at all, before any swap. That sharpens D05's "no single target".
- **Count law at p = 3.** A constant likelihood, two knots per variable, M_max = 2: 60 terms and 1831 states. D01's ratio gives the stated law (0.202, 0.404, 0.394) to 6·10⁻¹⁴. KHAOS's gives the stated law times 1/M! (0.252, 0.503, 0.245) to 4·10⁻¹⁴, 0.149 from the stated law.
- **Not covered:** p ≥ 4, and K_max = 3, where the Wallenius factor takes three-variable sets. Both grow the state space past dense GTH.

### Wallenius and defect fixtures (F-O5, F-O6) — 2026-10-27

From the reference-fixtures chip ([brief](reference-fixtures-brief.md)), `a6960f0`. Three oracles, run by `r/scripts/bass_stage1_checks.R`, all checks passing. They source the pinned R BASS checkout (`2671014`) through `_common.R`'s `bass_source`, and evaluate R's expressions as written, extracted from the function bodies, rather than transcriptions of them.
- **F-O5,** `r/oracles/bass_wallenius_oracle.R`. The reference is brute force over the K! draw orders. It sums to 1 over every enumerated order and equals 1/C(p, K) under equal weights, to 3·10⁻¹⁶. Its denominators sum the weights left in the urn: subtracting the drawn ones from the total, as a first version did, cost 9·10⁻¹³ with weights spanning five decades. R's own `sample(prob = z, replace = FALSE)` draws that law (largest |z| 2.8 over 20,000 draws per order).
  - *R's range.* `dmwnchBass` is exact at K = 2, 3: 5·10⁻¹⁵ with equal weights (p = 6), 1.4·10⁻¹³ with R-like counts (p = 8), 1.6·10⁻¹² with p = 50.
  - *Beyond it* the alternating sum loses digits: 5·10⁻¹³ at K = 4, 5 (p = 8) and 4·10⁻⁹ at K = 5 (p = 50). With weights spanning five decades it loses them from K = 2: 1·10⁻⁷ at K = 2, 3·10⁻⁵ at K = 3, 3·10⁻⁴ at K = 4, 8·10⁻³ at K = 5. Quadrature of the positive integrand (on t = e^{−u}) stays within 2·10⁻¹⁵ of brute force everywhere.
  - *At K = 1* `dmwnchBass` returns twice the probability, because its loop runs over `1:0`. R never calls it there (`logProbChangeMod` scores the uniform 1/p it samples), so this is a trap for a literal port, not an R defect.
- **F-O6 (a),** `r/oracles/bass_enumeration_oracle.R`, mode `swap`. The structure space is p = 2 with N = 24 on a 4 × 3 grid, K_max = 2 and M_max = 2: 62 terms, 30 with at least 3 nonzero values, 464 states and 6930 moving transitions. The continuous parameters are fixed at σ² = 0.1, τ = 1/N, λ = 2, the slice on which R's structure moves act. Stationary laws come from each kernel's transition matrix by GTH elimination, which keeps every component's relative accuracy (a dense solve left the corrected swap's spread at 10⁻⁶).
  - *Validity.* R's `getQf` agrees with the strict rank in both column orders on every state. R's `alpha` is exactly b times its b = 1 value. At b = 1, R's kernel with the state's move probabilities in place of its fixed thirds is the posterior to 2·10⁻¹⁴. With the thirds it is 0.020 off in total variation: F-O4's defect, measured in passing.
  - *R's tempered kernel targets none of the candidate laws.* At b = 0.5 and 0.2 its total variation is 0.086 and 0.145 from D05's likelihood path (prior times p_b), 0.036 and 0.078 from the prior times p_1^b, 0.26 and 0.51 from the powered collapsed path, and 0.32 and 0.58 from the joint over (structure, a) powered and then integrated.
  - *Swaps.* R's swap fails balance against its own rung laws: the spread of log s₁ − log s₂ − (b₁ − b₂)·lp is 1.78 (1 ↔ 0.5) and 1.07 (0.5 ↔ 0.2), where balance needs 0. Against the four candidate families it is 0.82 to 5.26. D05's corrected kernels reach their target to 4·10⁻¹⁶, and D05's swap balances on them to 5·10⁻¹⁵.
- **F-O6 (b),** `r/oracles/bass_jeffreys_oracle.R`. N = 40, an intercept and two hinges, σ² = 0.04. A tensor Gauss–Legendre rule for the integral over a matches the Gaussian closed form to 2·10⁻¹⁴. That confirms D03's coherent factors per birth, c(2πσ²)^{1/2}(|G′|/|G|)^{−1/2} (flat) and (2π)^{1/2} (Jeffreys), and its σ² shapes, g₁ + (N − k)/2 and g₁ + N/2.
  - *R's factors,* evaluated as written: the births give 2.667 and 2.272 on the log scale, against −0.967 and −0.179 (flat) and 0.919 (Jeffreys). The change move gives −0.047 against 0.093 and 0. The σ² shape is 19.5 against 19 and 20.
  - *The implied priors are confirmed by quadrature.* σ^{−2k}|G|^{3/4} reproduces R's birth, death and change factors, and σ^{−(k−1)} reproduces its σ² shape, each to 1·10⁻¹⁴. Two different priors, so no single target, as D03 says.
- **F-O6 (c),** `bass_enumeration_oracle.R`, mode `count`. A constant likelihood, p = 2 with two knots per variable, K_max = 2, M_max = 3, λ = 2: 24 terms and 2325 states. D01's birth ratio gives the stated count law (0.165, 0.330, 0.315, 0.190 for M = 0–3) to 6·10⁻¹⁴ in total variation. KHAOS's double count gives the stated law times 1/M! (0.241, 0.482, 0.230, 0.046) to 6·10⁻¹⁴, 0.229 from the stated one.
- **Findings for the drafts** (the briefs are not edited):
  - *D01, the inclusion–exclusion cut.* The cut at K ≤ 3 holds at R's default w₂ = 5, where the weights stay within a small factor of each other. With a small w₂ or skewed counts, the sum loses digits from K = 2 (above). Choosing by an error estimate, about ε·2^K/P(S), rather than by K alone would cover both. A port that routes K = 1 through the Wallenius formula must not copy R's loop.
  - *D05, R's lp.* `lp` departs from the collapsed log posterior by more than the term priors. It is the log joint density at a's conditional mode, so it also carries ½ log|G| − (k/2) log(2πσ²/(1 + τ)), which depends on the structure and on R's column scaling (`const`). The spreads above include both.
- **Not covered:**
  - The enumeration has p = 2, so every interaction draws the one pair and the Wallenius factor is 1.
  - R's tempered draws of λ, σ² and τ are held fixed on the slice.
  - Categorical and functional terms, and the robust modes.
  - The Jeffreys fixture is one design at one σ².

### Centred admission fixture (F-O7) — 2026-10-26

`r/oracles/bass_admission_oracle.R` gains two modes:
- `centred` reads the same design file as the default mode, with the intercept first, and drops the intercept.
- `hadamard` reads A and H for direct fits and builds the materialized design itself, for fixtures only.

The default mode's values are unchanged, apart from a new `mode` field. The QR of the explicitly centred design is the reference; it uses tol = 0, because R's LINPACK default drops nearly dependent columns. Cross-products are summed in plain double, pairwise or left to right, so that the numbers do not depend on the platform's long double. All checks are run by `r/scripts/bass_stage1_checks.R` and pass:
- *The identity.* On a hinge design, uncentred s = centred s · ρ to 3·10⁻¹⁵. The default mode agrees with the QR angles to 1·10⁻¹².
- *Near-constant column* (N = 10⁵, an indicator true for all but one case, ρ = 10⁻⁵).
  - Centred s is 1.0·10⁻⁹, which the Liang gate admits; uncentred s is 1.0·10⁻¹⁴, which R's convention refuses.
  - Relative errors in s: 3.4·10⁻³ by the downdate and 3.2·10⁻⁷ explicitly centred, both summed pairwise (the check asks for at least 100× apart, inside the band). Explicitly centred but summed left to right, 6.9·10⁻⁴.
  - With a thousand zeros (ρ = 10⁻²) the downdate and explicit forms agree, at 4·10⁻⁸ and 2·10⁻⁷. Left to right gives 1.8·10⁻³.
- *Centred Hadamard form* (n = 60, m = 80, means large against spread, with one design-only and one functional-only term). It matches explicit centring to 2·10⁻¹⁶ in the Gram and 3·10⁻¹⁶ in s; the downdate gives 4·10⁻¹³ and 1·10⁻¹³.

This corrects the scratch reading in the report below. Most of that error came from the summation order on near-constant columns, not from accumulating over the rows in general. D06's numerics bullet is revised to match, and the numerics segment gains N10 (pairwise or compensated cross-products). No C# side yet.

### D06's Gram under the Liang convention — 2026-10-26

From the user's review. D06 now says the gate reads the Gram that D03's intercept convention maintains. Under the Liang default that is the centred Gram, with s_j the reciprocal of the usual VIF. Under R's convention it is the uncentred Gram, with the intercept as a column. The two are related by uncentred s_j = centred s_j · ρ_j. `BassLiteral` reads `getQf`'s pivots rather than s_j.

A scratch check (not a fixture) confirmed the identity to 3·10⁻¹⁵ and found two numerical points, now in D06's numerics bullet:
- Accumulating the Gram over N rows, not ε alone, sets the error of s near the gate: up to 2·10⁻² at N = 10⁵.
- D03's rank-one downdate loses about ε/ρ_j on near-constant columns: 4·10⁻² against 8·10⁻⁴ at ρ = 10⁻⁵.

Both stay inside the factor-of-two confirmation band. The strict rebuild under the Liang convention uses explicitly centred columns, or the three-term centred Hadamard form in direct fits. New fixture item: F-O7.

### Andros et al. read — 2026-10-26

From arXiv 2406.18751 and 2602.09351. The direct functional algebra records that covariates known at the output locations enter as further functional variables (R's `xx.func` already allows it), and D01's knot wording now covers several functional variables. A residual-dependence check for direct fits is proposed (F11), and the output-side inner product is declared as prep for correlated residuals (F-X6). The deterministic-simulator profile reads that check.

### Lineage insights folded in — 2026-10-26

D01 gains the prior-cancelling schedule as an option and two "considered, not adopted" notes (KHAOS's delayed rejection, a uniform prior over basis functions). The deterministic-simulator profile is drafted. The hub's register ("What the lineage taught") records every insight and where it landed.

### First reference fixtures — 2026-10-25

`33b99bc`: `bass_marginal_oracle.R` and `bass_admission_oracle.R`, with `r/scripts/bass_stage1_checks.R` reproducing the D03 propriety and D06 admission numbers. No C# side yet; parity tests come with F3 and F7.

### D01–D06 drafted — 2026-10-25

Proposals, not decisions: `097e25b` (D01–D03), `53a4876` (D04, D06), `0844bdc` (D05), `a32c701` (categorical terms), `77274fb` (Jeffreys mode, propriety, order-free gate, s_min).

