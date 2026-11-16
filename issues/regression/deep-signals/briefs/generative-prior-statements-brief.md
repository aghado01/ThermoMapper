# Generative-prior statements: the Lean chip's fourth brief

Date: 2026-11-16. Status: **posed at the user's request on 2026-11-16** ("Yes, pose both", to Lean and numerical briefs on the compressed-sensing results). Drafted by the coordinating session. Its sibling is the numerical [generative-prior probes brief](generative-prior-probes-brief.md). The direction it serves is in the [approach brief](approach-brief.md#the-compressed-sensing-guarantees-used-by-analogy-not-yet-tested).

## Why

The compressed-sensing theorems for generative priors are the certificate the approach leans on: a code of dimension k is recoverable from a number of measurements linear in k, by descent, because the objective's landscape is benign. The work has used them by analogy only. Their proofs are long concentration arguments, but their deterministic skeleton is short and exact, and it is what the numerical probes will read against: the angle map a ReLU layer induces, the matrix the weight condition names, the exact zero set of the expected gradient, and what the isometry condition buys.

## Sources

LaTeX under `D:\aghado01\codex-scientiae\supellex-private\gauntlet\voroninski`.

- **`1705.07576v3/1705.07576v3-tex/deep_cs_arxiv.tex`** (Hand and Voroninski, *Global Guarantees for Enforcing Deep Generative Priors by Empirical Risk*):
  - "Main Results": both theorems, the two conditions, and the footnote that gives Q as an integral factoring in polar coordinates;
  - `\label{defn-g}` in "Proofs": the angle map;
  - "Approximate angle contraction property", for the properties of g the paper uses;
  - "Control of the zeros of h", Lemma `lemma:Seps`: the expected gradient h and the constant ρ_d.
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

**S20. The zeros of the expected gradient.** For nonzero x and x₀ with angle θ̄₀, put θ̄ᵢ = g(θ̄ᵢ₋₁) and

h/‖x₀‖ = −2^{−d}·(Π_{i<d} (π − θ̄ᵢ)/π)·x̂₀ + 2^{−d}·[‖x‖/‖x₀‖ − Σ_{i<d} (sin θ̄ᵢ/π)·Π_{i<j<d} (π − θ̄ⱼ)/π]·x̂.

- For d ≥ 2, h = 0 exactly when x = x₀ or x = −ρ_d·x₀, where ρ_d is the bracket's sum started from θ̌₀ = π.
- The argument the coordinating session expects: for 0 < θ̄₀ < π the two directions are independent and the first coefficient is positive, since every later angle is at most π/2 (S17). At θ̄₀ = 0 every angle is 0 and the bracket is ‖x‖/‖x₀‖ − 1. At θ̄₀ = π the first coefficient vanishes and the bracket is ‖x‖/‖x₀‖ − ρ_d.
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
