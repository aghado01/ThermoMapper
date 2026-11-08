# The lineage gaps' Lean leads — chip brief

Date: 2026-11-01. Status: done, report below. Written at the user's request from the Francom-lineage review.

## Context

The [lineage-gaps register](../../../issues/regression/bass/lineage-gaps-register.md) records where the Francom lineage and the literature around it stop short. This folder's protolemmata pose the statements behind it. They have passed numeric checks ([`lineage-gap-checks.py`](lineage-gap-checks.py)) and nothing more: no statement has been compiled, and every paper-dependent claim comes from an abstract or a search snippet.

**This chip's job:** turn the self-contained statements into compiling, proved enthymemata, then run the proof-walks now that the papers are on disk, and correct the record with what it finds. A statement that turns out false, or true only with a further hypothesis, is a finding worth more than a proof. The harness exists for those.

**Read first:**
- [lean/README.md](../../README.md): the taxonomy, meta-CI's rules, import discipline.
- This folder's four protolemmata, in full: [sensitivity-generator-lemmas.md](sensitivity-generator-lemmas.md) (SG), [structure-chain-lemmas.md](structure-chain-lemmas.md) (SC), [guarded-estimator-lemmas.md](guarded-estimator-lemmas.md) (GE), [contraction-walk.md](contraction-walk.md).
- [docs/Maths/PosteriorContraction.md](../../../docs/Maths/PosteriorContraction.md), and `lean/Enthymemata/BgpCovering.lean` for the namespace, the import style and how cited theorems enter as labelled axioms.

## Environment

- **Toolchain** Lean 4.32.2 and Mathlib v4.32.2, as pinned. Don't bump them.
- **Mathlib's build output is absent** from `lean/.lake` (`.lake/packages/mathlib/.lake/build` is empty; the compressed cache in `~/.cache/mathlib` is on C:), and **D: is 98% full, with about 6.7 GB free.**
  - Fetch by scoped module only: `lake exe cache get Mathlib.RingTheory.Polynomial.ShiftedLegendre`. A bare `lake exe cache get` restores everything; leave that out.
  - Check `df -h /d` before each fetch and after. **If a fetch would leave under 3 GB free, stop and report.** Keep to what the modules need, and delete nothing to make room: disk decisions are the user's.
  - Whether the rest of `lean/` still has its build output is unknown. Build only this chip's modules (`lake build Enthymemata.BassStructureChain`), and run `scripts/meta-ci.ps1 -NoBuild` for the taxonomy and ledger gates. Run the full `lake build` only if the disk allows, and say in the report which of these ran.
- **Iterate fast:** call `lean` directly with `LEAN_PATH` from `lake env printenv LEAN_PATH`. Linters run only under `lake build <Module>`, so finish with that.
- **Quirks:** `Π` is a reserved token; identifiers with a combining macron don't parse; `omit … in` and `open … in` go before the docstring; `lean/` is LF, so write files from Python with `newline='\n'`; meta-CI's ledger counts the word for an apology in comments too, so keep it out of docstrings. `π` as a variable name is fine unless `Real` is open.
- **Papers:** arXiv sources in `D:/aghado01/codex-scientiae/supellex-private/gauntlet/bayesian/<id>/<id>-tex/`; PDFs in `…/PDFs/Bayesian/`; Ghosal & van der Vaart 2017 by chapter in `…/PDFs/GhosalVanDerVaart2017/`. `pdftotext` works, and `pypdf` keeps math glyphs.
- **Python:** `python` is a pyenv `.bat` shim that breaks on `|` in arguments, so write scripts to files. It has no numpy; `C:/Users/azrie/PDenv/PyenvPython/pyenv-win/versions/3.12.6/python.exe` does.
- **The checkout is shared with other sessions.** Commit by pathspec (`git add -- paths`, then `git commit -- paths`), leave nothing staged, work on main, push nothing. End each commit message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Never touch `issues/regression/opus-bayesian-evals.md`, `perplexity-*.md` or anything outside the files named here.

## Scope

New modules, flat in `lean/Enthymemata/`, namespace `BASS`, each added to `Enthymemata.lean`: `BassStructureChain`, `GuardedEstimator`, `BassSensitivity`, `BassContractionCover`. Take the statements in the protolemmata as a first draft. Their Lean has never been compiled, so expect to adjust instances and hypotheses.

### Part 1: self-contained statements, in this order

Easiest and most load-bearing first, so an early stop still leaves value.

1. **SC-1:** `mhRatio_mul_swap`, `mh_detailed_balance`.
2. **SC-2:** `gap_le_cut`, through the two computations for the indicator of A.
3. **SC-3:** `khaos_precision_entry`, then the matrix form with `Matrix.hadamard` and positive definiteness of the sum. `Matrix.PosSemidef.hadamard` exists, though the sum form should not need it.
4. **GE-1:** `crossfit_gap`, then the conditional expectation it comes from, over simple random samples of the pool outside A, on a finite pool.
5. **GE-3:** telescoping unbiasedness, λ* as the minimizer of the stated quadratic, and Horvitz–Thompson unbiasedness under independent inclusion, all as finite statements. GE-2 follows the same pattern with an abstract expectation; take it if it comes cheaply.
6. **SG-1:**
   - `coordMean` idempotent, self-adjoint and commuting;
   - `generator_eigen`;
   - the components P_u as a complete orthogonal family;
   - Σ_v total(v) = Σ_u |u|·V_u;
   - `efron_stein`;
   - (e), the factorization of U_ρ on product-form functions.

   Prefer proving (b)–(d) over an abstract commuting family of self-adjoint idempotents, then instantiating, so the continuous case can follow.
7. **SG-2(a):** the 1-Lipschitz bound on √V_u.
8. **SG-3:** `mixture_product_integral`; then the continuous form from `MeasureTheory.integral_fintype_prod_eq_prod`, if cheap.
9. **SG-4:** `shiftedLegendre_jacobi`, from `coeff_shiftedLegendre`. Then orthogonality by integration by parts and the weighted identity, as far as they go without Legendre normalization.
10. **SG-5:** the 2 × 2 counterexample, with r = 1 by Rayleigh quotient.
11. **Contraction cover:** a union of finitely many subspaces' balls, in an M-dimensional Euclidean norm, covered at radius ξr/4 by a count free of r and of any norm bound. This is the walk's central claim about local entropy, and it is checkable with no paper. Build on `BgpCovering`'s lemmas where they fit.

### Part 2: proof-walks, with the papers

Method: step through the proof against BASS's setting. The first step that doesn't transport is the lemma to pose. Ask whether the proof goes through, not only whether the statement's hypotheses hold, and add new items to the protolemma with the walk's table.

1. **Ronen et al. (`2210.09352v1`) onto BASS's structure chain,** at fixed τ (SC-0, SC-4). What is the bottleneck set, which step uses trees, and what would the barrier's growth in N have to be? Note what Kim & Roćková (`2306.00126v1`) add on the upper-bound side.
2. **Woodard, Schmidler & Huber (`Woodard2009.pdf`).** State persistence exactly, and what measuring it on BASS would take from a flat-histogram run's ln g. Resolve the question SC-4 raises: the Potts toy fits tempering's replicas at N³ for q = 3 while the theorem, as reported in abstracts, says torpid for q ≥ 3. Read the toy's tests and the docs to tell whether the sizes are too short or the settings differ. Report; change no toy.
3. **The contraction walk against Ghosal & van der Vaart 2017, Chapter 8.** Confirm or correct:
   - that Theorem 8.19(b) as printed asks only for local entropy;
   - that the Bgp fixed-constant axiom's global cover is a strengthening, so a local restatement would stay a consequence of the printed theorem (propose it; don't edit `BgpUnknownVarianceTheorem`);
   - each row of the walk's table.

   *MARS via LASSO* is not on disk yet; if it arrives during the chip, read it for an M-term approximation rate. Otherwise leave prior mass conditional, as written.
4. **Prior-art checks for posed statements.** Read Roustant et al. (`1906.09883v1`, `1612.03689v1`) for whether they state SG-1's frame; Yoo & Ghosal (`PDFs/Bayesian/Ghosal2015-1411.6716v3.pdf`) and Odin et al. (`2403.03540v2`) for SG-2's limits; Bürkner et al. (`2204.06043v2`) for SG-4. Correct each "not found" and each abstract-level summary that reading changes, in the protolemma and the register, with what the paper actually says.

### Out of scope

- Any change to `src/` or `tests/`, including the exact-gap measurement in the register's G1.
- G6 (orthogonal discrepancy), which is a derivation for later, and the KHAOS package's source.
- Promotion to `Lemmas/`: proof-closed enthymemata wait for semantic review.
- New design. Where a walk suggests one, record it as a finding.

## Self-review

After each module, and after each walk, take a fresh-eyes pass before moving on, for correctness and for coverage gaps, and fix what it finds. For a module: does the statement say what the protolemma claims, with no hypothesis doing the work? Run `#print axioms` on the main results. For a walk: does each verdict rest on a line you read? Name the passes and their fixes in the report.

## Reporting

Append a report below, newest first, under "Reports". Update each protolemma's status line as its items land. It should say:
- what is proved, by declaration name, and what remains apologizing;
- every statement that changed, why, and any that proved false;
- each walk's table, and the lemma each walk posed;
- which build gates ran, and the disk's free space at the end;
- corrections made to the register, and any that change what a gap claims;
- what stopped you, if anything did.

## The 2026-11-08 queue: SG-4(b), SC-5's inverse-gamma form, CW-1

Three open items this chip left. The user asked for them on 2026-11-08, relayed by the BASS extensions session, to be taken in this order with a report after each. The method is this chip's: numeric checks first in `lineage-gap-checks.py`, standard axioms only, the module and barrel builds green, meta-CI passing, and D:'s free space watched.
1. **SG-4(b), the sum over a tensor basis.**
   - The statement: Σ_v E[x_v(1−x_v)(∂_v f)²] = Σ_α a_α²·Σ_v α_v(α_v+1) under the product uniform law, with total(v) ≤ ½·E[x_v(1−x_v)(∂_v f)²], sharp at degree 1.
   - It backs KHAOS slice 4's readouts from coefficients and the generator-eigenvalue prior.
   - Record what Roustant et al. already state.
2. **SC-5's proper IG(c, c·s²) form** of `collapsedLogMarginal_superset_ge`.
   - Its counterpart under D03's Jeffreys form is `jeffreys_superset_ge`.
   - Check which σ² prior is BASS's default, and whether the path profile and G1 readings need a caveat.
3. **CW-1, the small lemma** (`contraction-walk.md`).
   - The statement: with ‖f − f*‖ₙ ≤ σ₀ε/6 and |σ − σ₀| ≤ σ₀ε/6, (f, σ) lies in the KL ball of radius ε, through UV-1's per-observation bounds.
   - Only CW-1. Not the walk's other contraction items.

## Reports

### 2026-11-08: CW-1 proved, with ε ≤ 3 added

**Proved,** in `lean/Enthymemata/BgpGaussianKL.lean` (now 27 declarations, PROOF-CLOSED; `#print axioms` shows the standard three only).
- `klDiv_le_and_evariance_le_linear`: UV-1's per-observation bounds K ≤ 10D + 2Δ and V ≤ 18D + 16Δ, under |σ − σ₀| ≤ σ₀/2.
- `klBall_of_designNorm`: CW-1. With Σ_i (m_i − m₀_i)² ≤ n(σ₀ε/6)², that is ‖f − f*‖ₙ ≤ σ₀ε/6, and with |σ − σ₀| ≤ σ₀ε/6 and ε ≤ 3, the summed KL divergences and log-ratio variances are at most nε².
- `klBall_needs_eps_le`: without a bound on ε the statement is false.

**What changed in the statement: ε ≤ 3 is added.** The per-observation bounds need r ∈ [½, 2], which the slab gives only for ε ≤ 3, as in `klDiv_le_and_evariance_le_of_slab`.
- *Without it,* ε = 5.9 lets σ = σ₀/60 and the KL divergence is about 1795, against ε² = 34.8.
- *The route is unaffected:* the contraction route takes ε → 0, so the hypothesis costs it nothing.
- *The averages are tighter than stated:* at most ε²/3 and 17ε²/18.
- *The hypothesis 0 ≤ ε is not needed,* since a negative ε empties the conditions.

**Checks first** (`lineage-gap-checks.py`, CW-1): 5,000 seeded draws, with n up to 39 and ε ∈ [0, 3].
- The per-observation bounds hold in every draw.
- The averages reach at most 0.12ε² and 0.69ε².
- The ε = 5.9 failure is reproduced.

**What it decides.** The prior-mass step can be stated in the empirical norm BASS's fits use, with no sup-norm slab: CW-1 sits between `BgpGaussianKL` and the contraction route.
- The walk's other items, which were not in scope, stay open: the local `fixed_constant_contraction`, σ²'s upper tail by one test, and the fast rate at K ≥ 2.

**Gates.**
- `lake build Enthymemata.BgpGaussianKL`: clean, no warnings.
- `lake build Enthymemata`: green.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- D: had 40 GB free.

### 2026-11-08: SC-5's proper inverse-gamma form proved

**Proved,** in `lean/Enthymemata/BassContractionCover.lean` (now 31 declarations, PROOF-CLOSED; `#print axioms` shows the standard three only).
- `invGammaProfile`, the collapsed marginal's profile under IG(c, c·s²): −(ν/2 + c)·log(κ + (1 + g(1 − R²))/(1 + g)), with κ = 2c·s²/TSS.
- `invGammaProfile_zero`: at c = κ = 0 it is the Jeffreys profile.
- `invGammaProfile_mono`: it increases in R² for g > 0, κ ≥ 0 and ν/2 + c ≥ 0.
- `invGamma_superset_ge`: the superset step falls by at most (k/2) log(1 + g), as under the Jeffreys form.

**What changed in the statement.** Nothing. "c added to the exponent and c·s² to the sum" is exactly the profile above.
- *Derivation.* Given σ², the centred response is N(0, σ²(I + gP)), with |I + gP| = (1 + g)^M and quadratic form S = TSS·(1 + g(1 − R²))/(1 + g). Integrating σ² against IG(c, b) gives −(M/2) log(1 + g) − (ν/2 + c) log(b + S/2) plus a constant common to every structure.
- *Scope.* s² reads the response alone, so κ is the same for every structure. The bound needs the g-form, whose marginal is a function of R².

**Checks first** (`lineage-gap-checks.py`, SC-5 under IG).
- *The setup:* 400 seeded nested designs, with N from 8 to 29, M from 1 to 4, k from 1 to 3, c in {10⁻³, ½, 2}, and g from e⁻² to e⁷.
- *The closed form:* the log marginal, from the design's determinant and quadratic form and σ² by quadrature, equals −(M/2) log(1 + g) plus the profile, up to one constant within 9e-13.
- *The superset step* never fails.

**The decision: SC-5's bound holds under the σ² prior BASS fits actually run.** The default (`docs/Maths/BayesianRegression.md`, D03) is the g-form with level centring (ν = N) and `Sigma2Prior.ScaledToResponse(10⁻³)`, which is IG(10⁻³, 10⁻³·s²). The Jeffreys option was already covered.
- *The readings need no caveat on the g-form.* The path profile's superset reading (32,577,747 edges, no violation) and the G1 readings that cite SC-5 are now backed by a proof under either σ² prior.
- *The path profile's report does not name its σ² prior.* Either one is covered, so this changes nothing, but the brief could say which it ran.
- *A caveat applies only to the other coefficient forms* (`RidgeForm`, `RawRidge`, `ModifiedG`). Their marginals are not functions of R² alone, so SC-5 makes no claim there. No reading is known to have used them.

**Gates.**
- `lake build Enthymemata.BassContractionCover`: clean, no warnings.
- `lake build Enthymemata`: green.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- D: had 40 GB free.

### 2026-11-08: SG-4(b)'s tensor form proved

**Proved,** in `lean/Enthymemata/BassIntegrals.lean` (now 44 declarations, PROOF-CLOSED; `#print axioms` on the eight new results shows `propext`, `Classical.choice` and `Quot.sound` only).
- `integral_jacobi_tensorPartial`: per variable.
- `sum_integral_jacobi_tensorPartial`: summed over variables.
- `tensorTotal_eq`: the total numerator from the coefficients.
- `tensorTotal_le` and `tensorTotal_eq_half`: the bound and its sharpness.
- Beside them: `hasDerivAt_tensorExpansion` (`tensorPartial` is ∂_v f), `integral_tensor_sq` (Parseval on the tensor basis), `tensorCondMean_eq` (E_v f keeps the terms constant in x_v), and `integral_cube_prod` (Fubini on products).

**What changed in the statement.**
- Nothing.
- The coefficients are on the unnormalized P_n, as in the one-variable lemmas: a_α² N_α, with N_α = Π_v ⟨P_{α_v}, P_{α_v}⟩, is the squared orthonormal coefficient. So the statements read as the protolemma's.
- total(v) is defined as E[(f − E_v f)²], with E_v integrating out x_v, and proved equal to Σ_{α_v≥1} a_α² N_α. So the bound is about the Sobol total numerator itself, not a coefficient formula standing in for it.

**Checks first** (`lineage-gap-checks.py`, SG-4(b) tensor):
- 300 seeded sparse expansions, with d ≤ 3 and degree ≤ 3 per variable;
- both identities agree to 6e-13;
- the bound holds in every draw, equal at degree 1.

**Theirs and ours** (read on 2026-11-08):
- *Theirs: Roustant et al. 1906.09883v1.* They state the total index as a sum of squared coefficients with ℓ_v ≥ 1, and the per-variable derivative rewriting. The weighted form is given only as "adapted accordingly".
- *Theirs: Lamboni et al. 2013,* as cited in 1612.03689v1. The unweighted DGSM bound S^T ≤ C·ν/D.
- *Ours:*
  - the summed identity as the tensor generator's eigenvalue Σ_v α_v(α_v+1);
  - the weighted constant ½ on the Jacobi basis, with sharpness;
  - the formal proofs.

**What it fixes for KHAOS slice 4.** Under the product uniform law on [0, 1]^d and KHAOS's shifted Legendre basis, with b_α the orthonormal coefficients:
- each total numerator is Σ_{α_v≥1} b_α²;
- each weighted derivative measure is Σ_α b_α² α_v(α_v+1);
- their sum over v is Σ_α b_α² Σ_v α_v(α_v+1).

All of them are exact from the coefficients, with no quadrature, and total(v) ≤ ½ of the weighted measure. The readouts hold for that law and basis only. An input law other than the product uniform needs its own basis.

**Gates.**
- `lake build Enthymemata.BassIntegrals`: clean, no warnings.
- `lake build Enthymemata`: green.
- `scripts/meta-ci.ps1 -NoBuild`: PASS.
- D: had 40 GB free, and nothing was fetched.

### 2026-11-01: Parts 1 and 2 done, nothing false, three first-draft claims corrected

**Proved, with nothing apologizing.** All eleven of Part 1's items, and SG-3's measure form and SG-4(b)'s identities besides. `#print axioms` on 43 main results reads `propext, Classical.choice, Quot.sound` on every one: no cited axiom, no `sorryAx`.

| Module | Main declarations |
|---|---|
| `BassStructureChain` | SC-1: `mhRatio_mul_swap`, `mh_detailed_balance`, `not_isMhPair_of_mul_ne_one`. SC-2: `varianceOf_indicator`, `dirichlet_indicator`, `gap_le_cut`, `relaxation_ge_of_cut`, `gap_le_cut_of_stationary`. SC-3: `khaos_precision_entry`, `khaos_precision_eq`, `khaos_precision_eq_add_inv`, `khaos_precision_posDef`, `khaos_precision_of_const`, `khaos_precision_of_diagonal`. SC-5: `cutFlow_le_separator`, `gap_le_separator`, `collapsedLogMarginal_superset_ge` |
| `GuardedEstimator` | GE-1: `crossfit_gap`, `srs_mean`, `crossfit_bias`, `srs_difference_unbiased`. GE-2: `guarded_fold_unbiased`. GE-3: `telescope_srs`, `telescope_iid`, `tunedVariance_sub_tuned`, `tunedVariance_tuned_le`, `tunedVariance_one_gt_zero_iff`, `horvitz_thompson_unbiased` |
| `BassSensitivity` | SG-1 over `CommutingMeans`: `sum_component`, `component_component_of_ne`, `component_idem`, `generator_component`, `parseval`, `total_eq_sum`, `sum_total`, `dirichlet_eq_sum_total`, `efron_stein_form`, `efron_stein_form_eq_iff`, `noise_apply`, `form_noise`. On the finite product: `coordMean_idem`, `coordMean_comm`, `coordMean_selfAdjoint`, `efron_stein`, `efron_stein_eq_iff_additive`, `noise_prod_form`. `generator_eigen`. SG-2(a): `sqrt_compVar_sub_le`, `sqrt_compVar_sub_le_form`. SG-3: `mixture_product_integral`, `mixture_pair_integral`. SG-4(a): `shiftedLegendre_jacobi`. SG-5: `rayleigh_max_is_eigen`, `no_common_eigenvector`, `top_eigenspace_ne_subspace_mean` |
| `BassIntegrals` | SG-3: `mixture_product_integral_pi`. SG-4(b): `l2Form_derivative_jacobi`, `jacobiForm_shiftedLegendre`, `shiftedLegendre_orthogonal`, `variance_le_half_jacobiForm`, `variance_eq_half_jacobiForm_of_degree_one` |
| `BassContractionCover` | `card_le_of_pairwise_dist`, `exists_cover_subspace_ball`, `exists_local_cover_iUnion`, `exists_local_cover_with_scale`, `card_finsets_le`, `exists_dictionary_local_cover`, `log_dictionary_count_le`; `exists_mTerm_approx`; `sum_powersetCard_prod_le`; `jeffreysProfile_mono`, `rSquaredOf_mono`, `jeffreys_superset_ge` |

`BassIntegrals` is a fifth module beyond the brief's four, so `BassSensitivity`'s imports stay free of measure theory. Lean items beyond the brief's list came from the walks: SC-5, the M-term bound and the normalizer.

**Not done.**
- SG-4(b)'s sum over a tensor basis, Σ_v E[x_v(1−x_v)(∂_v f)²] = Σ_α a_α²Σ_v α_v(α_v+1). The one-variable identity is proved, and the sum is that identity per variable by Fubini.
- CW-1, the ‖·‖ₙ form of the KL ball, is posed, since it sits against `BgpGaussianKL`'s internals.
- The local form of `fixed_constant_contraction` is proposed in the walk. `BgpUnknownVarianceTheorem` is unchanged.
- GE-3's quadratic is taken as stated, per the brief, and not derived from independence.

**Statements that changed.** None proved false.
- `gap_le_cut` drops `∑ π = 1`, which the inequality never uses.
- `noise_prod_form` takes μ nonnegative, since it goes through the form in `productMeans`.
- `efron_stein`'s equality clause needed full support to read as "additive", which `efron_stein_eq_iff_additive` adds. Over the abstract family it reads V_u = 0 for |u| ≥ 2.
- `mixture_product_integral` and `generator_eigen` lose `DecidableEq` binders nothing reads.
- SG-5's example has rational entries (a, b at 53°) and proves more: the two matrices share no eigenvector.

**What the first draft had wrong**, as claims rather than as Lean statements.
1. *SC-4's linear barrier.* The expected linear growth in N, and so "the same exponential form", came from looking only at deaths. See the first walk.
2. *SC-3's "a third form unless g ∝ 1"* needs a Gram with off-diagonal entries: on a diagonal Gram it is diagonal for any g.
3. *Bhatnagar & Randall's q ≥ 3.* They proved q = 3 at α_c; Woodard et al. state q ≥ 3, α ≥ α_c, with some proofs in a thesis.

**Walk 1: Ronen et al. onto the structure chain.** Full table in [SC-4](structure-chain-lemmas.md).

| Step | Verdict |
|---|---|
| Finite reversible chain | Transports. F grows with N where theirs does not |
| Mixing from conductance | Transports |
| **The bottleneck: one edge through the trivial tree** | **Fails, and this is what uses trees.** BASS has no root, so a term can be removed from any state |
| Flow through one state | Posed as a separator: `gap_le_separator` |
| Likelihood ratio for one fixed tree | Transports only if a separator exists, and then uniformly over it |
| Two trees exact on f₀ | Changes: both sides need large explained variance |

*Lemma posed: SC-5.* A state can hold a term and its substitute, so S → S ∪ T → T loses at most k·½ ln(1 + g), not N times a variance (`jeffreys_superset_ge`, no hypothesis left on the profile). So the barrier's expected growth is polynomial in N, of degree growing with the k terms exchanged at once, and superpolynomial only if k or M grows with N. F-M2's births (−½ ln N at three N) against its deaths (median −234 nats) fit this. Kim & Roćková's lower bound has the same shape: a per-node Occam factor to the power of the nodes grown through noise. As printed it is n^{Θ(log n)} at L ~ log n, which the paper calls exponential. Their twiggy move suggests a several-term exchange, recorded as a finding. *Not verified:* this is a reasoned expectation about mixing, not a mixing bound. The [spectral-gap brief](../../../issues/regression/bass/spectral-gap-brief.md) pre-registered the same rival and can decide it.

**Walk 2: Woodard, Schmidler & Huber.** Persistence is γ(A, φ) = min{1, φ[A]/π[A]}. The bound and what measuring it takes are in SC-4: a basin label and a per-level tally, since ln g alone does not hold a set's share. *Potts toy: resolved, the sizes are too short.* The settings match the theorem's. The barrier per site is 0.00113 nats at q = 3 (computed, and matching the test's own header), so e^{cN} moves by about e^{2.2} across the sweep where N³ moves by e^{10.4}, and it overtakes N³ only near N ≈ 27,000. No toy changed. *To your attention:* `docs/Maths/Samplers.md` says the reading "agrees with" the theorem, which holds at q = 4 only. That file is outside this chip.

**Walk 3: the contraction chain against the book's Chapter 8.**

| Check | Result |
|---|---|
| Theorem 8.19(b) asks only local entropy | **Confirmed:** (8.21) |
| The fixed-constant axiom's global cover is a strengthening | **Confirmed:** the book names it, (8.22)–(8.23). A local restatement stays a consequence; proposed in the walk, not edited |
| UV-1 | Transports with one small lemma, CW-1: `slab_subset_klBall` is sup-norm, and the prior mass works in ‖·‖ₙ |
| UV-5 | Transports. Lemma 8.27(ii) needs less: σ bounded above only |
| UV-2 | New. **"Fixed g" cannot be a constant** (CW-2): the coefficient mass costs N‖f_M‖ₙ²/(2gσ²), so g has to be of order N |
| UV-3 | New and simple, as expected. The default cap min(1000, N − 1) has to grow for the asymptotics |
| UV-4 | New, and **proved:** `exists_dictionary_local_cover` |
| The normalizer | **No obstacle:** W_F lies in [1, 2.59] at the defaults, free of N and the design |
| σ²'s upper tail | Looks removable by one test rather than a sieve. An expectation, not checked in Lean |

*MARS via LASSO* arrived during the chip. **It states no M-term rate,** and "logarithmic in dimension" means n's exponent is free of d, with constants that depend on d, so it is a fixed-d result. *Lemma posed, CW-3:* `exists_mTerm_approx` (Maurey) gives ε(M)² = V²R²/M in ‖·‖ₙ, reading neither the dictionary nor the design. So prior mass is no longer wholly conditional: ε_n² ≈ V·(log(pN)/N)^{1/2} for every interaction order. That is slower than least squares' n^{−4/5}. The fast rate holds by hand at K = 1 and is open at K ≥ 2.

**Walk 4: prior art.** Corrections are in SG's prior art and the register.
- Roustant et al. (`1906.09883v1`) state SG-1(b) and SG-4(b)'s identity outright. What they lack is a resampling generator over subsets and a prior from the eigenvalues.
- Σ_v total(v) = Σ_u |u|·V_u is, from memory and unread, mean dimension and total influence.
- Yoo & Ghosal's quasi-uniform knots are exactly SG-2's missing bound.
- Odin et al. hold the subspace as a parameter, so that does not transport, but it shows a route from L² through detectability.
- Bürkner et al.'s R2D2 is exchangeable across polynomials.

**Register corrections that change what a gap claims.**
- *G1:* expected barrier growth is polynomial, not linear, and the Potts question is closed.
- *G2:* several expectations became proved or corrected, listed there; the prior-mass step has a slow rate; "logarithmic in dimension" is corrected.
- *G5:* what is likely new narrows to the variance-by-order readout and the prior.

The reading list now says what was read.

**Self-review passes and their fixes.**
- *`BassStructureChain`:* dropped the unused `hπ`; removed `hP`, unused in the separator bound.
- *SC-5:* the pass found `h R2 ≤ h R2'` doing the work, so `rSquaredOf_mono` and `jeffreys_superset_ge` now discharge it for the Jeffreys form. The proper IG form is argued, not written out.
- *`BassSensitivity`:* `form_algebra_adj`'s first proof skipped the fact that a product of self-adjoints needs commuting. `noise_prod_form`'s first form ran over a list fold rather than the operator, and now reads `F.noise`. "Additive" was stated as V_u = 0, and is now proved literally.
- *All modules:* `DecidableEq` binders nothing reads were removed. `set_option linter.unusedFintypeInType false` stays on four lemmas that need finiteness in the proof only: `eq_of_update_eq`, `eq_of_update_eq_of_ne`, `nonempty_of_prob`, `khaos_precision_posDef`.
- *Walk 1:* the two-trees row was too loose, since Lemma 3 needs both sides.
- *Walk 2:* a `pypdf` re-read corrected Corollary 3.1, where plain-text extraction had garbled γ's argument.
- *Walk 3:* the g-form bound was re-derived in general g, which produced CW-2.

**Build gates.**
- *On `lean/.lake` (D:):* `lake build` of this chip's five modules, 2,704 jobs, 0 warnings. `scripts/meta-ci.ps1 -NoBuild` passes, and its ledger reads all five PROOF-CLOSED.
- *Full `lake build`:* in a scratch workspace on C:, with the same sources, toolchain and Mathlib: 3,150 jobs, 0 errors. Its 8 warnings are the apologies already in `BARS` (6), `Stability` and `CumulativeField/Local`.
- *Not run:* the full build on D:. It would only repeat the C: result, so `lean/.lake` holds output for this chip's closure alone. The first full build there will fetch about 0.45 GB more.

**Disk.**
- D: started at 6.6 GB free. This chip's scoped fetches took about 2.6 GB, including the first closure and the cache tool, which left a low of 3,387 MB, above the 3 GB floor. It reads **14,269 MB free** at the end, because something outside this chip released about 11 GB. Nothing was deleted to make room.
- Iteration ran in `C:\Users\azrie\AppData\Local\Temp\bw`, shared-object clones of `lean/.lake/packages` with Mathlib unpacked. It is safe to delete.

**What stopped me.** Nothing. The user's reboot paused the chip once, and the compile restarted from disk.
