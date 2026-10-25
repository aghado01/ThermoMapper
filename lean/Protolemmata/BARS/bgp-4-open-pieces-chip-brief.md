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
