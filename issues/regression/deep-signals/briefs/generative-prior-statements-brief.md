# Generative-prior statements: the Lean chip's fourth brief

Date: 2026-11-16. Status: **posed at the user's request on 2026-11-16** ("Yes, pose both", to Lean and numerical briefs on the compressed-sensing results). Drafted by the coordinating session. Its sibling is the numerical [generative-prior probes brief](generative-prior-probes-brief.md). The direction it serves is in the [approach brief](approach-brief.md#the-compressed-sensing-guarantees-used-by-analogy-not-yet-tested).

## Why

The compressed-sensing theorems for generative priors are the certificate the approach leans on: a code of dimension k is recoverable from a number of measurements linear in k, by descent, because the objective's landscape is benign. The work has used them by analogy only. Their proofs are long concentration arguments, but their deterministic skeleton is short and exact, and it is what the numerical probes will read against: the angle map a ReLU layer induces, the matrix the weight condition names, the exact zero set of the deterministic comparison field h, and what the isometry condition buys.

## Sources

LaTeX under `D:\aghado01\codex-scientiae\supellex-private\gauntlet\voroninski`.

- **`1705.07576v3/1705.07576v3-tex/deep_cs_arxiv.tex`** (Hand and Voroninski, *Global Guarantees for Enforcing Deep Generative Priors by Empirical Risk*):
  - "Main Results": both theorems, the two conditions, and the footnote that gives Q as an integral factoring in polar coordinates;
  - `\label{defn-g}` in "Proofs": the angle map;
  - "Approximate angle contraction property", for the properties of g the paper uses;
  - "Control of the zeros of h", Lemma `lemma:Seps`: the deterministic comparison field h and the constant ρ_d.
- **`2011.00288v2`** (`dense_pr.tex`), "Proof Technique": the same expectation for the phaseless case, Φ_{x,y}.
- **`1805.08855v2`**, the autoencoder proposition: the k/n noise reduction.

## The statements

Numbering continues. Each is the coordinating session's posing from the sections above; take the paper's exact forms where they differ and say what changed.

**S17. The angle map.** g(θ) = arccos(((π − θ)·cos θ + sin θ)/π) on [0, π].
- Its argument lies in [0, 1], so g maps [0, π] into [0, π/2]; g(0) = 0 and g(π) = π/2.
- g(θ) < θ on (0, π]. The heart is sin θ > θ·cos θ there.
- g is increasing, so the iterates θᵢ = g(θᵢ₋₁) decrease to 0 from any start.
- Add whichever further properties of g the paper's angle-contraction section states and uses.

**S18. The weight condition's matrix as an expectation.**
- *The angular identity.* With u(φ) = (cos φ, sin φ) and M(θ) = [[cos θ, sin θ], [sin θ, −cos θ]],
  (1/π)·∫_{θ−π/2}^{π/2} u(φ)·u(φ)ᵀ dφ = ((π − θ)/2π)·I + (sin θ/2π)·M(θ).
  The interval is the set of directions with positive inner product with both e₁ and (cos θ, sin θ).
- *The Gaussian form,* if Mathlib's Gaussian measure makes the polar factorization workable: for a standard Gaussian w in the plane, E[1{w·x > 0}·1{w·y > 0}·wwᵀ] is that matrix for x = e₁ and y at angle θ. Otherwise carry the factorization as a named premise.
- *The link to S17.* From Q's form, x̂ᵀQ_{x,y}ŷ = ((π − θ)·cos θ + sin θ)/2π and x̂ᵀQ_{x,x}x̂ = ½. So the cosine of the angle between the two outputs of a layer that meets the condition exactly is cos g(θ).

**S19. What the isometry condition buys.** Suppose |⟨A(G x₁ − G x₂), A(G x₃ − G x₄)⟩ − ⟨G x₁ − G x₂, G x₃ − G x₄⟩| ≤ ε·‖G x₁ − G x₂‖·‖G x₃ − G x₄‖ for all codes.
- Then (1 − ε)·‖G x₁ − G x₂‖² ≤ ‖A(G x₁ − G x₂)‖² ≤ (1 + ε)·‖G x₁ − G x₂‖².
- For ε < 1, A∘G is injective wherever G is.
- The compressed objective is sandwiched by the uncompressed one: (1 − ε)·f_I ≤ f_A ≤ (1 + ε)·f_I, with f_A(x) = ½‖A(G x − G x₀)‖².

**S20. The zeros of the explicit deterministic comparison field.** For nonzero x and x₀ with angle θ̄₀, put θ̄ᵢ = g(θ̄ᵢ₋₁) and

h/‖x₀‖ = −2^{−d}·(Π_{i<d} (π − θ̄ᵢ)/π)·x̂₀ + 2^{−d}·[‖x‖/‖x₀‖ − Σ_{i<d} (sin θ̄ᵢ/π)·Π_{i<j<d} (π − θ̄ⱼ)/π]·x̂.

- For d ≥ 2, h = 0 exactly when x = x₀ or x = −ρ_d·x₀, where ρ_d is the bracket's sum started from θ̌₀ = π.
- For 0 < θ̄₀ < π the two directions are independent and the product in the leading term is positive, since every later angle is at most π/2 (S17). At θ̄₀ = 0 every angle is 0, so the bracket is ‖x‖/‖x₀‖. The leading term supplies the −1 when the two coefficients of x̂ = x̂₀ are combined. At θ̄₀ = π the first coefficient vanishes and the bracket is ‖x‖/‖x₀‖ − ρ_d.
- ρ₂ = 1/π, and ρ_d > 0 for d ≥ 2. That ρ_d → 1 is the paper's claim; prove it if it comes cheaply.

**S21. A k-dimensional prior keeps k/n of white noise.** For an orthogonal projection P of rank k on ℝⁿ and a vector η with second moments E[ηηᵀ] = σ²·I, E‖Pη‖² = σ²·k, against E‖η‖² = σ²·n. State it on second moments.

**S22, optional and last. Sign patterns.** n hyperplanes through the origin of ℝᵏ cut it into at most 2·Σ_{j<k} C(n − 1, j) regions. This is the count behind the log(n₁⋯n_d) in the sample complexity: a bias-free ReLU layer is linear on each region.

## Scope and order

- **Order:** S19, S21, S17, S20 first; then S18; S22 last.
- **A witness for every statement** and `#print axioms` for every closed result. A named premise gets an instance that meets it wherever one is cheap.
- **Placement** in the harness's stages as before; name modules by what they state.
- **Left to the numerical chip:** every concentration claim, and whether anything of ours meets either condition.
- **Left alone:** `src/`, `tests/`, `docs/`. The Q extension's build still waits on the user's approval and is not part of this brief.

## Done

- S17, S19, S20 and S21 proved without `sorry`, each with its witness. S18's angular identity proved, its Gaussian form proved or compiling with its premise named.
- A report appended to this brief: per statement, the module, status, axioms and each change from the posed form. It ends in a **Methods** section that documents the logic of the proofs concisely, with exposition and the key mathematical points, like a dense methods paper. It carries no pointers into `docs/`.

## Working rules

- Commits land on `main` by pathspec, each one informative about its work, with no co-author trailer; nothing stays staged.
- `lake build` runs at below-normal priority with capped threads.
- Each step ends with a fresh-eyes review of its own output, fixed before it is reported.
- Messages are terse events; results travel in the report.

## Round 1 report — 2026-11-16

S19 and S21 are proof-closed in `Enthymemata` for semantic review. S17, S20, S18 and optional S22 remain queued in that order. The S20 posing above now identifies the deterministic comparison field and locates the −1 in the leading term at angle zero.

### S19 — generator secants

Module [IsometryCondition.lean](../../../../lean/Enthymemata/GenerativePrior/IsometryCondition.lean), commit `b7706b6`. All names below are in `GenerativePrior`.

- `RangeRestrictedIsometry` quantifies over two arbitrary generator secants. `rric_sandwich` bounds their diagonal squared norms. `rric_injective` proves separation of any two distinct generated outputs for ε < 1, and `rric_injOn` transfers injectivity from G to A ∘ G on any specified set of codes.
- `rric_objective_sandwich` bounds `compressedObjective` by 1 ± ε times `uncompressedObjective`. Both use the same reference code and the factor ½. The statement allows an arbitrary code type, real inner-product spaces and a linear measurement map.
- The checkout already contained an untracked draft of this module. This round retained its RRIC, diagonal and objective proofs, pairwise separation proof, and `rric_scaled_id` witness. It added `rric_injOn`, fixed the ambiguous linear-map type in the draft's final example, named it `rric_scaled_id_separates`, registered the module in the harness, and added and checked all six axiom prints. With G the identity on ℝ, A = √(3/2) I meets RRIC at ε = ½ and separates codes 1 and 2.

### S21 — projected second moments

Module [IsotropicProjection.lean](../../../../lean/Enthymemata/GenerativePrior/IsotropicProjection.lean), commit `8590038`. The projection P is fixed before sampling. `IsotropicSecondMoments` requires a probability measure μ, integrability of every product ηᵢηⱼ, and the actual integral ∫ ηᵢηⱼ dμ = σ² δᵢⱼ. Its scalar parameter `v` denotes the common raw diagonal moment σ². This premise permits dependent coordinates and requires no Gaussian assumption or separate centering hypothesis.

- `IsotropicSecondMoments.integral_quadratic` proves E[ηᵀMη] = σ² tr(M) for every fixed real matrix M. `ambient_energy` and `projected_energy` in that namespace give E‖η‖² = σ²n and E‖Pη‖² = σ²k for k = rank(P). `projected_fraction` gives E‖Pη‖² = (k/n) E‖η‖² when n > 0, including the zero-energy case.
- `SignedAxes.law`, `noise` and `proj` specify four equally likely outcomes ±2e₁, ±2e₂ and projection onto e₁. `SignedAxes.law_isProbabilityMeasure`, `isotropic`, `orthogonal`, `rank_one` and `energies` verify the premises and the expected energies. Each draw has energy 4 and its projection has energy 4 or 0. Their means are 4 and 2, respectively.
- This proves the fixed-projection expectation requested in S21. Proposition `prop:autodec` in the cited source concerns a Gaussian probability bound across piecewise-linear autoencoder regions. That concentration result has different hypotheses and remains outside this round.

### Replay

Both statement modules, their import module and the `Enthymemata` aggregate compiled in the shared checkout. The taxonomy check passed. The two proof files carry all 15 `#print axioms` commands. Each lists exactly `propext`, `Classical.choice` and `Quot.sound`. The aggregate replayed existing `sorry` warnings in `BARS`, `CumulativeField.Local` and `Stability`.

From the repository root, use PowerShell with the portable Elan launcher on PATH.

```powershell
(Get-Process -Id $PID).PriorityClass = 'BelowNormal'
$env:LEAN_NUM_THREADS = '2'
$env:ELAN_HOME = Split-Path -Parent (Split-Path -Parent (Get-Command lake.exe).Source)
Set-Location -LiteralPath lean
lake build +Enthymemata.GenerativePrior.IsometryCondition +Enthymemata.GenerativePrior.IsotropicProjection +Enthymemata.GenerativePrior +Enthymemata
lake env lean Enthymemata/GenerativePrior/IsometryCondition.lean
lake env lean Enthymemata/GenerativePrior/IsotropicProjection.lean
pwsh -NoProfile -File scripts/meta-ci.ps1 -NoBuild
```

### Doc changes owed

At the next documentation pass, state the secant hypothesis alongside the objective bound and identify k/n as a fraction of expected energy under a fixed projection.

### Methods

For S19, use the same secant u in both slots of RRIC and replace each self-inner-product by a squared norm. Splitting the absolute value yields both bounds. If Au vanishes and ε < 1, the positive factor 1 − ε forces u = 0. Linearity identifies Au with the difference of the measured outputs. Choosing the reference code as the second endpoint and dividing both inequalities by two gives the objective bounds.

For S21, expand ηᵀMη as a finite sum of Mᵢⱼηᵢηⱼ. Product integrability justifies moving the integral through both sums. Substituting the given moments removes the off-diagonal terms and leaves σ² tr(M). For M = I this trace is n. Symmetry and idempotence give PᵀP = P, so the projected energy is the quadratic form for M = P. The harness's existing `ActiveSubspace.IsOrthProj.rank_eq_trace` identifies its trace with k.

## Round 2 report — S17

S17 is proved in [AngleMap.lean](../../../../lean/Enthymemata/GenerativePrior/AngleMap.lean), commit `0b3ced3` on main. In namespace `GenerativePrior`, `reluAngleCos` is the argument of arccos and `reluAngle` is the angle map g. The depth-n comparison angle for S20 is `reluAngle^[n] θ`, with depth zero equal to the starting angle.

| Statement | Declarations |
| --- | --- |
| The argument decreases strictly from 1 to 0 on [0, π], and cos(g(θ)) equals that argument. | `hasDerivAt_reluAngleCos`, `strictAntiOn_reluAngleCos`, `reluAngleCos_zero`, `reluAngleCos_pi`, `reluAngleCos_mem`, `cos_reluAngle` |
| g is continuous and strictly increasing, with g(0) = 0, g(π) = π/2 and exact image [0, π/2]. | `continuous_reluAngle`, `strictMonoOn_reluAngle`, `reluAngle_zero`, `reluAngle_pi`, `reluAngle_mem`, `reluAngle_mapsTo`, `reluAngle_image` |
| θ cos θ < sin θ and g(θ) < θ throughout (0, π]. Every allowed orbit decreases weakly to 0, strictly for a positive start. Each positive-depth angle is in [0, π/2]. | `mul_cos_lt_sin`, `reluAngle_lt`, `reluAngle_iterate_antitone`, `reluAngle_iterate_strictAnti`, `reluAngle_iterate_tendsto_zero`, `reluAngle_iterate_succ_mem` |
| On (0, π), g has derivative in [0, 1]. Its Lipschitz constant on the closed interval is at most 1. Step errors bounded by δ accumulate to at most nδ. | `hasDerivAt_reluAngle`, `deriv_reluAngle_mem`, `lipschitzOnWith_reluAngle`, `reluAngle_iterate_error` |
| The antipodal witness follows π ↦ π/2 ↦ arccos(1/π). It attains the common upper bound for all orbits of depth at least two. | `reluAngle_iterate_two_pi`, `reluAngle_iterate_le_arccos_inv_pi` |

**Source coverage.** In `deep_cs_arxiv.tex`, the angle-contraction section (lines 546–624) uses monotonicity and comparison with the identity. Its layer lemma, `lemma:control-thetatwo`, bounds the one-step angular error by 4√ε under the weight distribution condition with 0 < ε < 0.1. The next section (lines 707–720) uses nonexpansiveness to propagate this error and uses the two-step angular bound. The Lean error theorem accepts measured angles in [0, π] and the one-step error as hypotheses. Substitution of δ = 4√ε yields the paper's d-layer estimate. The derivative theorem is stated on the open interval, while the finite-difference theorem includes both endpoints. The WDC-to-layer estimate and the separate reciprocal rate bounds in the later zero-set analysis (line 1136) are outside this proof.

### Replay

The module, its import module and the `Enthymemata` aggregate built successfully. The aggregate replayed the existing warnings listed in round 1. Direct elaboration printed exactly `propext`, `Classical.choice` and `Quot.sound` for all 25 theorems, including the antipodal witness. The taxonomy check passed and marked the module proof-closed. Its declaration counter reported 23 of the 27 declarations (25 theorems and two definitions). The matcher at `lean/scripts/meta-ci.ps1` line 168 skips the four theorem lines beginning with `@[simp]`.

With the portable Elan launcher on PATH, run from the repository root in PowerShell.

```powershell
(Get-Process -Id $PID).PriorityClass = 'BelowNormal'
$env:LEAN_NUM_THREADS = '2'
$env:ELAN_HOME = Split-Path -Parent (Split-Path -Parent (Get-Command lake.exe).Source)
Set-Location -LiteralPath lean
lake build +Enthymemata.GenerativePrior.AngleMap +Enthymemata.GenerativePrior +Enthymemata
lake env lean Enthymemata/GenerativePrior/AngleMap.lean
pwsh -NoProfile -File scripts/meta-ci.ps1 -NoBuild
```

### Doc changes owed

Record the exact angular range, the antipodal two-step bound, and the hypothesis on per-layer errors alongside their use in the comparison field. S20, S18 and optional S22 remain queued.

### Methods

Write q for the arccos argument. Its derivative −(π − θ) sin θ / π is negative on (0, π). The endpoint values bound q, and the two order reversals in arccos ∘ q make g strictly increasing. The derivative of sin θ − θ cos θ is θ sin θ, positive in the interior. Comparison with its value at zero gives q(θ) > cos θ for 0 < θ ≤ π, hence g(θ) < θ. Each orbit is bounded below by zero. Continuity makes its limit a fixed point, and strict decrease excludes every positive fixed point.

For the Lipschitz bound, put a = (π − θ)/π and b = sin θ/π. Both are nonnegative and a + b ≤ 1 because sin θ ≤ θ. Since q = a cos θ + b, expansion gives q² + (a sin θ)² ≤ (a + b)² ≤ 1. On the interior this bounds g′ = a sin θ / √(1 − q²) by 1. The mean value theorem and continuity at the endpoints give the closed-interval estimate. The triangle inequality adds at most δ per layer. Monotonicity makes the antipodal start maximal at depth two, and each subsequent step can only decrease the angle.
