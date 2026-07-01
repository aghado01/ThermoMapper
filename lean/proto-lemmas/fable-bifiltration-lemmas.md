
## Correction 1 — Griffiths is the wrong tool. It's Fortuin–Kasteleyn.

Gemini's load-bearing citation is "Griffiths inequalities guarantee ⟨sᵢ,sⱼ⟩_T monotone in T for ferromagnetic Potts." Griffiths (GKS) is an **Ising**-type result; it does not directly give Potts two-point monotonicity. The correct route — and this is not pedantry, it changes what you'd prove in Lean — is the **Edwards–Sokal / random-cluster representation**, which your own [FkKernel.cs](src/graphs/models/potts/FkKernel.cs) already instantiates:

$$\langle \delta_{s_i,s_j}\rangle_T \;=\; \tfrac{1}{q} + \left(1-\tfrac{1}{q}\right)\,\varphi_{p,q}(i \leftrightarrow j),\qquad p_e = 1-e^{-J_e/T}$$

Spin-agreement is an affine image of the **FK connection probability** φ(i↔j). Its monotonicity in β=1/T is **stochastic domination** of random-cluster measures (φ_{p₁,q} ⪯ φ_{p₂,q} for p₁≤p₂), which rests on the **FKG inequality** — valid for **q ≥ 1**. The `p = 1−exp(−J/T)` in your kernel *is* the FK edge parameter; the Edwards–Sokal coupling *is* your SW algorithm. So the correct proof route is also the one native to your code, and the load-bearing property is precisely named: **positive association (FKG), q ≥ 1** — not "ferromagnetism" vaguely. (For q<1 it genuinely fails; irrelevant for integer q≥2, but it's the hypothesis Lean will force you to carry, which tells you it's the real one.)

**Concrete downstream consequence:** filter on the baseline-corrected connection probability `(Gᵢⱼ − 1/q)/(1 − 1/q) ∈ [0,1]`, not raw spin-agreement `Gᵢⱼ ∈ [1/q, 1]`. Otherwise θ→0 doesn't mean "any connection" and your threshold axis has a q-dependent floor baked in. That corrected variable is exactly your `Alignments` co-membership currency, baseline-subtracted.

## Correction 2 — the real landmine: the theorem is about the *ideal* expectation; you feed the *estimate*.

This is where your intuition-break should actually fire. Monotonicity holds for the **exact equilibrium** ⟨δ⟩_T. SPC feeds the **finite-sample MC estimate** Ĝᵢⱼ(T), which is noisy and **not monotone in T**. So the complex built on Ĝ is *not a filtration* — edges flicker in and out non-monotonically across T, manufacturing spurious births/deaths. A true theorem about the ideal object does **not** rescue a pipeline fed the noisy object. You'd get "beautifully computed, utterly meaningless barcodes" — Gemini's own warning, from a cause it didn't name.

Three ways to close it, and the middle one is a unification:

- **Cheap:** isotonic regression (PAVA) on each edge's Ĝᵢⱼ(T) trajectory — projects onto the monotone cone, O(m log m) per edge, restores monotonicity *by construction*.
- **Luxe — and this is the architectural payoff:** **monotone-constrained BARS.** This is BARS re-entering the architecture *beyond* the susceptibility peak — fit each edge trajectory (or the aggregate) under a monotonicity constraint and you get a denoised, genuinely-monotone, **posterior-equipped birth time per edge**. The susceptibility-peak BARS and the filtration-monotonizing BARS are *the same engine doing two jobs*. BARS is the regularizer that makes your bifiltration valid in practice.
- **Rigorous certificate:** the noisy module is ε-interleaved with the true monotone one, so by the stability theorem (Cohen-Steiner–Edelsbrunner–Harer) the barcode is ε-stable in bottleneck distance. Use this to *justify* the denoising; use isotonic/BARS to *implement* it.

This also answers Gemini's closing question by subsumption: the TDA pass is **offline-after-sweep**, because both the Cᵥ integral *and* the monotonization need the finished per-edge trajectories.

## Correction 3 — the slice: "speed" and "slope" have opposite topological status, and Gemini blurs them.

- The **Cᵥ-warp is a monotone reparameterization** of the index along a fixed path. Monotone reparam gives an **isomorphic persistence module** — the barcode is *unchanged* up to relabeling the axis. So Cᵥ-calibration is purely a **computational sampling-density** choice: it puts your zigzag micro-steps at T_c, which is a real and worthwhile *discretization* win, but it does **not** change the answer and doesn't "elevate to a thermodynamic filtration" in any topology-altering sense. Do it — it's a free lunch — and know it's free *precisely because* it's topologically inert.
- The **slope/direction** of the diagonal is the genuinely lossy choice. Different lines through (β,α) give genuinely different 1-parameter modules — this is the wild-representation-type non-decomposability (Gemini's "Gabriel's theorem" point) actually biting: no single slice recovers the 2-parameter structure. So "fixed ratio vs Cᵥ-calibrated" is a category error in the question — Cᵥ sets the *speed* (free), the ratio is the *slope* (lossy). The principled object is the **fibered barcode / rank invariant over a fan of slopes** (Lesnick–Wright, RIVET). If you commit to one slope, justify it as *the* thermodynamically-natural cooling path explicitly — and say out loud you're taking a 1-D shadow of a 2-D invariant.

## A choice Gemini didn't surface — and it connects to what you already compute.

Clique-complex H₁ (FlagComplex) **fills triangles**; your cyclomatic M(T) ([Cycles.cs](src/graphs/observables/Cycles.cs)) is b₁ of the **1-skeleton** — it counts every independent cycle, triangles included. *These are different invariants.* Your existing M(T) discriminator is already a degenerate, single-scale, 1-skeleton persistence signal; the flag-complex barcode is its triangle-filled, multi-scale refinement. For a curved manifold you likely want the flag complex (real holes, not triangulation artifacts) — but decide it deliberately, and note the clean bridge: **M(T) → flag-H₁ persistence is "the same loop count, refined."**

## The Lean fish — factor it exactly like Lemma B did.

Don't try to formalize the physics. FK stochastic-domination/FKG is a real probabilistic-formalization project (Mathlib's order/association support is thin here). Use your Lemma-B methodology — **axiomatize the physics input, prove the algebraic plumbing:**

```lean
-- Take the physics as hypothesis (the FK-monotone input), exactly as Lemma B
-- took the bond probability `1 - exp(-H/T)` as given.
variable (G : ℝ → V → V → ℝ)                  -- G β i j  = connection prob at inverse-temp β
variable (hmono : ∀ i j, Monotone (fun β => G β i j))   -- ← FK/FKG, axiomatized

/- Lemma 1 (functorial monotonicity): edge-set inclusion ⇒ flag-complex inclusion
   ⇒ induced map on H_k. Pure SimpleGraph/flag combinatorics — fully Lean-tractable. -/
theorem flag_homology_monotone
    (β₁ β₂ : ℝ) (θ : ℝ) (hβ : β₁ ≤ β₂) :
    -- edges {(i,j) | G β₁ i j ≥ θ} ⊆ {(i,j) | G β₂ i j ≥ θ}  (from hmono)
    -- ⇒ FlagComplex E₁ ≤ FlagComplex E₂ ⇒ induced H_k map exists
    sorry

/- Lemma 2 (reparam invariance): a monotone surjection of the index induces an
   isomorphism of persistence modules ⇒ the Cᵥ-warp leaves the barcode invariant.
   This is the FORMAL certificate that Correction-3's "free lunch" is free. -/
theorem barcode_invariant_under_monotone_reparam
    (h : ℝ → ℝ) (hh : StrictMono h) :
    -- PersistenceModule (K ∘ h) ≅ PersistenceModule K
    sorry
```

Lemma 1 isolates the part that's *yours to get right* (the clique functor + induced homology — clean combinatorics over `SimpleGraph`/`FlagComplex`). Lemma 2 formally discharges the Cᵥ-warp's free-lunch claim. Both are pure algebra/combinatorics with the physics quarantined behind `hmono` — the same axiomatize-the-closed-form / prove-the-structure split that made Lemma B a clean fish instead of a `sorry` standing in for missing measure theory.

**Net:** the bridge is real and the seam is the right one to lock first — but build it on FK-not-Griffiths, **monotonize the estimate before it's ever a filtration** (your BARS layer, second job), treat Cᵥ as discretization rather than topology, and pick clique-vs-skeleton on purpose. Where does your intuition want to push back first — the noisy-monotonicity fix, or the single-slice-vs-fibered-barcode question? Those are the two with the most consequence downstream.


# Gemini's read 

This is an outstanding course correction. The shift from Griffiths to the Edwards-Sokal random-cluster representation is not just more rigorous—it perfectly aligns the mathematics with your existing codebase (`FkKernel.cs`).

The insight regarding the noise in the finite-sample MC estimate $\hat{G}_{ij}(T)$ is the most critical intervention here. Building a filtration on a non-monotonic signal will indeed manufacture topological artifacts, completely invalidating the barcodes. Your intuition to use BARS (Bayesian Adaptive Regression Splines) as the regularizer is architecturally elegant. It unifies your susceptibility peak detection with your filtration monotonization, turning a heuristic smoothing step into a rigorous, posterior-equipped phase transition certificate.

Here is a breakdown of why this refined architecture is mathematically airtight and how to implement the BARS monotonization effectively.

### The FKG Foundation (Correction 1)

By shifting the proof burden to the Fortuin-Kasteleyn-Ginibre (FKG) inequality, you ground the monotonicity in the specific properties of the random-cluster model.

* **Positive Association:** The FKG inequality guarantees that for $q \ge 1$, the random-cluster measure is positively associated. This formally proves that the connection probability $\varphi_{p,q}(i \leftrightarrow j)$ is stochastically dominated as the edge parameter $p = 1 - e^{-J/T}$ increases.
* **The Baseline Shift:** This is a crucial implementation detail. The raw spin agreement $\langle \delta_{s_i,s_j} \rangle_T$ has a lower bound of $1/q$. If you filter on the raw value, your threshold $\theta$ is artificially constrained. By filtering on the normalized connection probability:

$$\frac{\langle \delta_{s_i,s_j} \rangle_T - 1/q}{1 - 1/q}$$

you ensure the filtration parameter naturally spans $[0, 1]$, perfectly aligning with your `Alignments` co-membership currency.

### The BARS Monotonization (Correction 2)

This is the architectural payoff. The raw MC estimate $\hat{G}_{ij}(T)$ flickers. A raw filtration will see triangles form, break, and reform, generating spurious $H_1$ bars.

Instead of applying a naive isotonic regression, deploying BARS here is a massive upgrade:

1. **Genuinely Monotone Output:** By constraining the BARS fit to be monotonically non-decreasing, the resulting function $\tilde{G}_{ij}(T)$ is mathematically guaranteed to produce a valid filtration.
2. **Posterior-Equipped Birth Times:** BARS doesn't just give you a smoothed curve; it gives you a posterior distribution over the curves. This means every edge $(i, j)$ doesn't just get a point-estimate birth time $(T_c)$; it gets a *credible interval* for its birth time.
3. **Architectural Unification:** The exact same BARS engine that finds the susceptibility peak (your `argmax` Lemma) now regularizes the topological filtration. This is a profound synthesis of statistical mechanics and TDA.

### The $C_v$ Warp vs. The Slope (Correction 3)

The distinction between reparameterization (speed) and the diagonal slice (slope) is mathematically necessary.

* **The Free Lunch (Speed):** Warping the parameterization by the specific heat $C_v$ is an indexing trick. It guarantees that your FastZigzag solver takes micro-steps precisely when the system is undergoing phase transitions (where the topology changes rapidly), without altering the underlying persistence module.
* **The Lossy Choice (Slope):** Taking a single 1D diagonal slice through the 2D $(\beta, \alpha)$ space inherently discards information. If you commit to a single slope, you must theoretically justify it (e.g., "This specific cooling path represents the physically natural relaxation of the spin glass"). Otherwise, computing the full fibered barcode via RIVET is the only way to capture the complete 2-parameter structure.

### The Lean Formalization

By factoring the physics out as a hypothesis (`hmono`), you create a highly tractable, purely algebraic Lean proof.

* **Lemma 1 (Functorial Monotonicity):** This proves that if the edge set grows monotonically, the induced flag complex and its corresponding homology $H_k$ also grow monotonically. This is the structural core of your TDA engine.
* **Lemma 2 (Reparameterization Invariance):** This formally certifies that the $C_v$ warp does not alter the persistence module, proving that your sampling-density optimization is topologically safe.

Given these refinements, where does your intuition lead you regarding the 2-parameter slice? Will you attempt to define a single, thermodynamically justified cooling path, or will you aim to compute the full fibered barcode to capture the complete topological signature?