# BGP-4 open pieces (UV-1c, UV-4a, UV-6) — chip brief

Date: 2026-10-24. Status: open, for a separate thread. Written at the user's request from the main regression thread.

## Context

BGP-4 extends TWCD2025's GP contraction theorem to an unknown noise variance, the case the implementation already runs
(`GpNoiseModel.UnknownScale`). The protolemma [bgp-unknown-variance-lemmas.md](bgp-unknown-variance-lemmas.md) splits
the extension into cited pieces, pieces that need a proof, and pieces worth formalizing.

**Already done**, in `lean/Enthymemata/BgpUnknownVariance.lean`, which compiles with no `sorry` (PROOF-CLOSED in
meta-CI's ledger): UV-1a, UV-1b, UV-1d, UV-2a, UV-2b and UV-2b′. Their names are listed in the protolemma's
"Lean-worthy pieces" table.

**This chip owns the three open rows of that table:**

| Piece | Statement | Protolemma section |
|---|---|---|
| UV-1c | The Gaussian KL and log-ratio variance closed forms | UV-1 |
| UV-4a | Covering a scaled family and an interval | UV-4 |
| UV-6 | The composed theorem, as an enthymema | UV-6 (with UV-5 as its forward boundary) |

Read first:
- [lean/README.md](../../README.md): the taxonomy (protolemmata → enthymemata → lemmas → theorems), meta-CI's rules and
  import discipline.
- The protolemma, in full. Its citation boundary (TWCD2025 Theorems 3.1–3.3 and (59)–(61); Ghosal & van der Vaart
  2007, Theorem 4; Hellinger tests) is settled: cite those, don't re-prove them.
- `lean/Enthymemata/BgpUnknownVariance.lean`: the proved pieces, the namespace (`BGP`) and the import style to match.
- `lean/Enthymemata/Stability.lean`: the house pattern for a cited theorem entering Lean as a labelled `axiom` with its
  reference, which UV-6 will need for GvdV07 Theorem 4.

## Environment

- **Toolchain** `leanprover/lean4:v4.32.2`, with Mathlib pinned in `lean/lake-manifest.json` (rev `905b958…`). Work at
  that pin; don't bump it.
- **Build and check:** `pwsh lean/scripts/meta-ci.ps1` (build plus taxonomy checks; `-Validate` for strict ledger mode).
  Compile a single module directly while iterating.
- **Cache:** fetch only the modules a new import needs, for example `lake exe cache get Mathlib.Topology.MetricSpace.CoveringNumbers`.
  A bare `lake exe cache get` restores all of Mathlib, so avoid it.
- **Imports:** scoped modules only; the `import Mathlib` umbrella fails meta-CI.

## Scope

### UV-1c — the Gaussian KL and log-ratio variance closed forms

For one observation, with δ = μ − μ₀ and r = σ₀/σ:

```text
KL(N(μ₀,σ₀²) ‖ N(μ,σ²)) = −log r + (r² − 1)/2 + δ² r²/(2σ₀²)
Var_{N(μ₀,σ₀²)}[log p₀/p] = (r² − 1)²/2 + δ² r⁴/σ₀²
```

- **What the pinned Mathlib has** (checked 2026-10-24):
  - `gaussianReal μ v` with `integral_id_gaussianReal`, `variance_id_gaussianReal`, `memLp_id_gaussianReal` (all moments
    finite), and the affine map lemmas (`gaussianReal_map_const_mul`, `…_add_const`).
  - `InformationTheory.klDiv` (`Mathlib/InformationTheory/KullbackLeibler/Basic.lean`), with `klDiv_of_ac_of_integrable`
    and the log-likelihood ratio `llr`.
  - `ProbabilityTheory.centralMoment`.
  - **Not found:** a closed-form Gaussian KL, or the Gaussian fourth moment.
- **The route.** Write the log ratio as a quadratic in x. Its mean needs E z and E z². Its variance also needs E z³ = 0 and
  E z⁴ = 3 (Var z² = 2). Derive those from the moment API, or from the map lemmas and symmetry. If the fourth moment is
  out of reach, state the variance identity with it as an explicit hypothesis and record the gap.
- **Done when:** both identities are stated against Mathlib's `klDiv` and variance, proved or with only a recorded,
  named gap. The protolemma states them in the `(μ₀, σ₀, μ, σ)` parameterization; match it, so UV-1a/b/d compose with no
  algebra at the seam.

### UV-4a — covering a scaled family and an interval

- **The statement.** A cover of B at radius ε/(2b) and a cover of [a, b] at radius ε/(2R), with R bounding B, cover
  {s • g | s ∈ [a, b], g ∈ B} at radius ε. So the covering numbers multiply. Also state the interval's own count:
  [a, b] at radius ε needs at most ⌈(b − a)/(2ε)⌉ + 1.
- **The proof idea:** s•g − s′•g′ = s•(g − g′) + (s − s′)•g′.
- **What the pinned Mathlib has:** `Mathlib/Topology/MetricSpace/CoveringNumbers.lean`, with `Metric.coveringNumber`,
  `externalCoveringNumber`, `IsCover`, `IsCover.coveringNumber_le_encard`, monotonicity (`coveringNumber_anti`,
  `…_mono_set`). Covering numbers are `ℕ∞` and radii `ℝ≥0`.
- **Choose deliberately:** internal or external covering numbers. The product step is easiest with external covers built
  from the two given covers. If the protolemma's (61′) needs the internal one, the factor-2 conversion is standard.
- **Done when:** the product bound and the interval bound compile without `sorry`, in a normed space over ℝ general
  enough to take the GP sieve's sup-norm space.

### UV-6 — the composed theorem, as an enthymema

- **The statement.** Under Assumptions 1–2 and the proper prior [P1] with σ₀ ∈ (σ_min, σ_max), the posterior puts
  vanishing mass outside {‖f − f*‖ₙ + |σ − σ₀| ≤ C ε̄ₙ}. See the protolemma's UV-6 for the rate and the dependence of C.
- **The shape.** An implication assembled from:
  - the prior-mass condition (UV-1, UV-2);
  - the sieve and remaining mass (UV-3);
  - entropy (UV-4, with UV-4a);
  - the tests (UV-5).

  The general theorem (GvdV07 Theorem 4) enters as a cited `axiom`, following `Stability.lean`'s pattern, with its
  reference in the docstring.
- **Where the `sorry`s sit.** By design, exactly at:
  - UV-5, the tests for ‖f − f′‖ₙ + |σ − σ′|;
  - UV-4's "to confirm against the paper" check (Rₙ at most polynomial in n);
  - UV-3's sieve, if not discharged along the way.

  Each `sorry` gets a comment naming its protolemma section.
- **The goal is the statement, not the proof.** The composed statement must compile, so the dependency structure is
  checked, with the forward boundaries visible in meta-CI's ledger as `apologizing(n)`.
- **Where it lives:** a new enthymema module, for example `Enthymemata/BgpUnknownVarianceTheorem.lean`, registered in
  `Enthymemata.lean`. That way the PROOF-CLOSED bookkeeping module stays closed and promotable on its own.

## Working rules

- **The math is settled.** Keep the protolemma's statements and citation boundary. If a statement turns out wrong or
  mis-scoped while formalizing, stop and record it in the protolemma, rather than adjusting it silently to fit what proves.
- **Don't promote.** Promotion to `Lemmas/` follows a semantic review, which is the user's call.
- **Commits:** targeted commits on `main`, one concern each, ending with the project's co-author trailer.
- **Keep the protolemma current.** Update its "Lean-worthy pieces" table and the per-piece "Tier" notes as pieces land.

## Done criteria

1. UV-1c and UV-4a compile without `sorry`, or each has its gap recorded as a named hypothesis and the reason noted in
   the protolemma.
2. UV-6's statement compiles as an enthymema, with `sorry`s only at the forward boundaries listed above.
3. `pwsh lean/scripts/meta-ci.ps1` is green, and the ledger shows the new module's apologies, with nothing new in the
   verified stages.
4. The protolemma's table reflects the new status.
5. A report is appended below.

## Report

*(To be appended by the thread that picks this up: what landed with commit hashes, the ledger state, the gaps recorded
and why, and anything in the protolemma that needed correcting.)*

### 2026-10-24 — chip thread (UV-1c, UV-4a, UV-6)

**What landed.** These are targeted commits on `main`, in this order:
- `f08fabe` UV-1c in `Enthymemata/BgpGaussianKL.lean`, PROOF-CLOSED.
- `7230805` UV-4a in `Enthymemata/BgpCovering.lean`, PROOF-CLOSED.
- `dcbe6a3` UV-6 in `Enthymemata/BgpUnknownVarianceTheorem.lean`, `apologizing(1)`.
- The protolemma update and this report, in the commit that adds this text.

Each code commit registers its module in `Enthymemata.lean`, and each builds on its own.

**Ledger.** `pwsh lean/scripts/meta-ci.ps1`: PASS.

| Module | Decls | Sorries | State |
|---|---|---|---|
| BgpCovering | 7 | 0 | PROOF-CLOSED |
| BgpGaussianKL | 24 | 0 | PROOF-CLOSED |
| BgpUnknownVariance | 9 | 0 | PROOF-CLOSED (unchanged) |
| BgpUnknownVarianceTheorem | 39 | 1 | apologizing(1) |

- Nothing new is in `Lemmas/` or `Theorems/`, and nothing was promoted.
- The new modules build with no linter warnings beyond the one `sorry`.
- `#print axioms BGP.posterior_contraction_unknown_variance` gives `sorryAx` (only through `exists_tests`) and
  `BGP.gvdv07_contraction`. Every supporting lemma depends only on `propext`, `Classical.choice` and `Quot.sound`.

**UV-1c: no gap.**
- `klDiv_gaussianReal` states the KL against Mathlib's `klDiv`. `variance_llr_gaussianReal` and
  `evariance_llr_gaussianReal` state the variance against `variance` and `evariance` of `llr`.
- Both use the protolemma's `(μ₀, σ₀, μ, σ)`, with the variance written as `(σ ^ 2).toNNReal`.
- The pinned Mathlib has no Gaussian fourth moment. `moments_gaussianReal_zero` reads it off the MGF, through
  `iteratedDeriv_mgf_zero` and four explicit derivatives of exp(v t²/2).
- `klDiv_le_and_evariance_le_of_slab` composes UV-1a, UV-1b, UV-1c and UV-1d into the protolemma's slab inclusion,
  with c = σ₀/6.

**UV-4a: no gap.**
- The cover-level lemma is `isCover_image2_smul`. It works with internal covering numbers, since products of
  internal covers stay inside the scaled family.
- Also landed: the brief's radii (`…_le_half`), the external form through Mathlib's factor-2 comparison, and the
  interval bound ⌈(b − a)/(2ε)⌉ + 1.
- All of it holds over any seminormed space on ℝ.

**UV-6: the statement compiles, and the proof is the composition.** `posterior_contraction_unknown_variance` feeds
the cited axiom `gvdv07_contraction` from:
- UV-1 and UV-2: `eventually_prior_mass`, using `exists_amplitude_window`, `noise_slab_mass`, `prior_slab_mass`,
  UV-2a and UV-2b′;
- UV-3: `prior_compl_sieve_le`, proved;
- UV-4 with UV-4a: `eventually_sieve_cover` and `closure_sieve_subset_cover`;
- UV-5: `exists_tests`, the one `sorry`, commented with its section.

Where this departs from the brief's list of `sorry`s, and why:
- **UV-4's Rₙ check is the named hypothesis `hRadius`, not a `sorry`.** The same goes for UV-2's "to confirm" check,
  which is `h59`. Both are properties of the paper's prior and sieve, which Lean does not construct and the theorem
  takes as hypotheses. A `sorry` would have asserted them for an arbitrary sieve, where they are false. So it would
  not be a deferred true claim.
- **UV-3's sieve is discharged**, so it carries no `sorry`.
- **[P1]'s priors on s and η enter by the properties the proof uses.** These are: support in the box with positive
  mass on subintervals, and at least the log-uniform mass on subintervals, which the log-uniform law meets with
  equality. The theorem therefore covers a class of priors that contains [P1].
- **The axiom's form.** It has two rates, a fixed M, a testing hypothesis for a general semimetric, and hypotheses
  stated on the strong side (measurable sieves and neighbourhood parts, tests at every radius, nεₙ² → ∞). That keeps
  the axiom a consequence of the cited theorem.
  - **To confirm at review:** whether the printed Theorem 4 states the fixed-M form or only "every Mₙ → ∞".
  - The KL neighbourhood is taken in `ℝ≥0∞`, so a parameter with an infinite divergence is not counted as 0.

**Corrections to the protolemma** (recorded there, not applied silently):
- **The exponent 4 in (60′) is not enough by itself** (under UV-3). Composed with the general theorem, the
  remaining-mass exponent must beat (59′)'s prior-mass exponent by a σ₀-dependent margin. In units of the slab
  radius that is C₂ > C₁ + 144/σ₀² (`hC₂`), and (59′)'s exponent is not 1.
  - The sieve can be enlarged to meet the margin.
  - **To confirm against the paper:** that TWCD2025's sieve gives (60) at any exponent with (61) keeping its order.
  - This is the one place the formal statement departs from the protolemma's text, and it does so through an
    explicit hypothesis.
- **A refinement of UV-2 step 4.** The slab factor costs an arbitrarily small fraction of nεₙ², not "one more
  nεₙ²".
- **A stale note.** `Pow.Asymptotics` is now in the local cache.

**Environment notes.**
- Caches were fetched by module: `lake exe cache get Mathlib.Probability.Distributions.Gaussian.Real
  Mathlib.InformationTheory.KullbackLeibler.Basic Mathlib.Topology.MetricSpace.CoveringNumbers`. The import closure
  is 3103 files, and `Pow.Asymptotics` and `PiL2` came with it.
- `Π` is a reserved token, and identifiers with a combining macron (`C̄`) do not parse. The Lean names are `μW`, `Cb`
  and `kb`.

**Open after this chip:**
- the UV-5 tests;
- the three "to confirm" checks against the papers: `h59`, `hRadius`, and (60) at any exponent;
- the GvdV07 fixed-M citation detail;
- semantic review of BgpCovering, BgpGaussianKL and BgpUnknownVariance for promotion. BgpGaussianKL imports
  BgpUnknownVariance.

### 2026-10-25 — follow-up (UV-5 tests, the paper checks, the GvdV07 check)

**What landed.** These are commits on `main`:
- `d971edb` UV-5's Gaussian tails in `Enthymemata/BgpGaussianTails.lean`, PROOF-CLOSED and axiom-free.
- `17aae2c` in `Enthymemata/BgpUnknownVarianceTheorem.lean`:
  - UV-5 proved (`exists_tests`);
  - the general theorem re-cited, in two routes;
  - the composed theorems.
- This commit: the protolemma update and this report.

**Ledger.** `pwsh lean/scripts/meta-ci.ps1`: PASS. The new BgpGaussianTails is PROOF-CLOSED, and
BgpUnknownVarianceTheorem goes from `apologizing(1)` to PROOF-CLOSED.
- No `sorry` is left in the BGP-4 modules, and nothing was promoted.
- `#print axioms` shows each composed theorem resting only on its own cited axioms plus Lean's standard ones.
- BgpUnknownVarianceTheorem is PROOF-CLOSED in the ledger, but its theorems are conditional on cited axioms. Whether
  that module can be promoted is a review question.

**UV-5, the tests.**
- `exists_tests` holds with K = σ_min²/(1024σ_max⁴) and ξ = 1/8, split into three case lemmas:
  - the likelihood-ratio mean test;
  - the scale-up test centred at f*;
  - the scale-down test centred at f₁.
- The tails come from completing the square under `gaussianReal`, independence under `Measure.pi` (`iIndepFun_pi`)
  and Chernoff, all as Lebesgue integrals, so with no integrability side conditions.
- **A correction to the protolemma's sketch.** Both scale tests cannot be centred at f₁. Under the truth that leaves
  a mean shift of order r, which can cancel the scale gap when r ~ σ₀. Recorded under UV-5.

**The checks against TWCD2025** (local copy, `graveyard/bibliotecha/.../TWCD2025.md`, Appendix A.1):
- **`h59`, (59) uniform over the amplitude box: confirmed.** Part I uses f* only through (A2)'s constants, and (A2)
  is monotone. So f*/s satisfies (A2) with (ν₁/s_min, ν₂/s_min²) for every s in the box, and Part I run once with
  those constants gives the bound uniformly. That is stronger than the continuity the protolemma asked for.
- **(60) at any exponent: confirmed.** The exponent 4 comes only from the constants in eq. (80). Scaling rₙ^ϱ and Nₙ²
  gives any C₂, and Part III's entropy then changes only through C̄₂. That is the tunable form, and
  `posterior_contraction_unknown_variance_of_tunable_sieve` states it, removing the margin hypothesis.
- **`hRadius`, Rₙ polynomial: confirmed.** Rₙ ≤ Nₙ(rₙ/δₙ)^{D/2} + εₙ, since the kernel has h(0) = 1, and this is
  polynomial in n by eq. (80) and δₙ.
- `h59` and `hRadius` stay hypotheses, because they are facts about objects Lean does not construct. Their docstrings
  now say "confirmed".

**The GvdV07 check** (arXiv 0708.0491, read from the PDF): **the printed Theorem 4 is not the fixed-M form.**
- It has one rate and concludes for every Mₙ → ∞. A remark after it allows any semimetric with Lemma 2's test
  property.
- The fixed-M two-rate form that the protolemma, TWCD2025 Theorem 3.1 and vdVvZ08b Theorem 3.3 state is Ghosal &
  van der Vaart (2001, Theorem 2.1). vdVvZ09 §3 says so, calling it "a slight generalization of the main result in
  [GGvdV00], formulated as Theorem 2.1 in [GvdV01]".
- For independent, non-identically distributed observations that form is not printed as a theorem.

What changed in Lean as a result:
- **The old axiom is renamed and re-cited.** `gvdv07_contraction` is now `fixed_constant_contraction`, with the
  citation above. It still gives `posterior_contraction_unknown_variance`, the protolemma's fixed-C statement.
- **The printed route is added.** `gvdv07_theorem1` and `gvdv07_lemma1` are the printed Theorem 1 and Lemma 1.
  `posterior_contraction_unknown_variance_gvdv07` composes them in Lean: Theorem 1 on the sieve at Lε̄ₙ, Lemma 1 off
  it at 6εₙ/σ₀, glued by `posterior_le_add`. The result is the printed conclusion, every Mₙ → ∞, with the smaller
  margin C₂ > C₁ + 72/σ₀².
- **The citation boundary is unchanged, and the finding is recorded.** A dated "Citation check" note sits under it in
  the protolemma.

**Prior art found on the way.** van der Vaart & van Zanten (2008, Theorem 3.3) is the unknown-σ, fixed-design result
for a single GP prior with σ on [a, b]. Its proof is a two-line sketch that leaves the unknown-scale tests to the
reader. BGP-4 is its extension to the hierarchical bandwidth prior and the amplitude/noise-ratio parameterization,
with those tests now written out.

**Environment note.** This session ran under heavy memory pressure: 63 GB committed on 31.5 GB of RAM. A Lean check
that loads Mathlib took 8–10 minutes instead of 1–2, so iteration was slow.

**Open after this:**
- semantic review and promotion: the four PROOF-CLOSED bookkeeping modules, and the axiom-conditional theorem module;
- whether to keep both routes of the general theorem, or only the printed one;
- the non-Lean items of UV-6: random design, the posterior mean, the σ Bernstein–von Mises.
