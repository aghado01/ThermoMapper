Yes — Grimmett 2006 is *the* canonical source for exactly the FKG/stochastic-domination machinery I leaned on, and consulting it directly is how we de-risk the axiomatized step rather than hand-waving it. Let me navigate the contents.

The Grimmett machinery is exactly what grounds our axiom. Let me pin the precise theorem chain — Chapter 3 has the random-cluster comparison inequalities and Chapter 1 the Potts coupling.

Grimmett doesn't just *support* the axiom — it lets us **shrink** it. That's the real payoff of consulting it: the precise theorem chain shows most of what I lumped into `hmono` is actually Lean-provable plumbing, and only a small probabilistic core needs citing. Here's the exact chain, with the numbers pinned.

## The theorem chain (each step now has a citation or a proof)

The estimand we need monotone is spin-agreement ⟨δ_{sᵢ,sⱼ}⟩_T. Decompose via Edwards–Sokal and walk it:

1. **`p_e(β) = 1 − e^{−βJ_e}` is non-decreasing in β** for J_e ≥ 0. `dp_e/dT = −e^{−J_e/T}·J_e/T² < 0`, so every edge parameter rises as you cool. *Lean-trivial* (`Real.exp` monotonicity).
2. **{i ↔ j} (open-subgraph connectivity) is an increasing event** — adding edges never breaks a connection. *Lean-tractable* (`SimpleGraph.Reachable.mono` over the edge set).
3. **Comparison inequality (Grimmett Thm 3.21):** `φ_{p₁,q} ≤_st φ_{p₂,q}` when `p₁ ≤ p₂`, `q ≥ 1`. This is the genuinely hard probabilistic core — it rests on **positive association (Thm 3.8)** → **FKG lattice condition** → **Holley (Thm 2.1)**. *This is the one to axiomatize.*
4. **Apply ≤_st to the increasing event:** `φ_{p₁,q}(i↔j) ≤ φ_{p₂,q}(i↔j)`. *Lean-trivial* — it's literally the definition of `≤_st` (Grimmett §2.1: µ₁ ≤_st µ₂ ⟺ µ₁(X) ≤ µ₂(X) for all increasing X).
5. **Edwards–Sokal two-point identity:** `⟨δ_{sᵢ,sⱼ}⟩_T = 1/q + (1−1/q)·φ_{p(β),q}(i↔j)` — from the cluster-spin coupling (Grimmett §11.2, the "sample ω from φ, paint each open cluster a uniform colour" recipe). *Cite as the second input.*
6. **Compose:** affine with positive slope (`1−1/q > 0` for q ≥ 2) ∘ monotone ⟹ ⟨δ⟩_T non-decreasing in β = `hmono`. *Lean: composition.*

So the refined factoring **axiomatizes exactly two Grimmett results** (Thm 3.21 + the §11.2 coupling) and **proves the other four steps**. That's a much tighter quarantine than "take all of `hmono` on faith" — it isolates the irreducibly-deep theorem (FKG/Holley on the random-cluster measure) from the calculus and combinatorics that are yours to discharge.

## Two things Grimmett makes precise that you'd want to know

**The `q ≥ 1` hypothesis is load-bearing, and Grimmett flags exactly where it fails.** §3.3 notes outright "the failure of stochastic monotonicity when q < 1." Your Potts has integer q ≥ 2, so you're safely interior — but this tells you *which* property carries the bridge: not "ferromagnetism" loosely, but **positive association of the random-cluster measure, which holds iff q ≥ 1**. That's the hypothesis your Lean axiom should name, and it's the honest minimal one.

**You're in the cleanest possible regime.** Thm 3.8 carries the explicit caveat (line 130) that it "applies only to finite graphs G" — and your kNN coupling graph *is* finite, free boundary. So you need **none** of the Chapter 4 infinite-volume machinery, none of the phase-coexistence/non-uniqueness subtleties. The comparison inequality applies directly. The finite graph is the easy case, not a compromise.

## The epistemic upgrade — and it retroactively justifies last turn's monotonization

This is the part worth sitting with. Grimmett **proves the estimand ⟨δ⟩_T is genuinely monotone** (3.8 + 3.21). So when you isotonic-regress or monotone-BARS the noisy SW estimate Ĝ(T), you are **not** imposing a convenient fiction to make PH behave — you're projecting the estimator back onto a property the random-cluster theory *guarantees the true target has*, which finite sampling corrupted. That flips the monotonization from "a hack" to "principled regularization toward a theorem-certified structure." And it's clean under your own validation-independence statute: the monotone constraint originates **outside** the estimator (Grimmett's theorem, not the data it's fitting), so enforcing it isn't circular. The monotonizer recovers signal the theory promises is there.

## The one adversarial check — and it's specific

Theorem 3.21 is stated for a **scalar** p (homogeneous edges). Your couplings are **heterogeneous** — `J_e = kernel(d_g(i,j))` varies per edge, the whole point of curved-metric coupling. Cooling moves the *vector* **p**(T) = (p_e(T))_e, each coordinate by a different amount. So the result you actually need is the **inhomogeneous comparison inequality**: `φ_{**p₁**,q} ≤_st φ_{**p₂**,q}` whenever `**p₁** ≤ **p₂**` coordinatewise. This almost certainly holds — the random-cluster measure with per-edge p_e is standard (Grimmett's Thm 3.1 conditional probabilities are already stated edge-wise), and the Holley argument is coordinate-local — but it's the precise statement to confirm in his framework rather than assume, because it's the version your code instantiates. Step 1 already gives you what you need to invoke it: all p_e move the *same direction* as T drops (J_e ≥ 0), so **p**(T) is coordinatewise-monotone and the inhomogeneous ≤_st (once confirmed) closes the chain.

Net: Grimmett converts the bridge from "cite a vague Potts-Griffiths folklore" to a four-step proof resting on two named, finite-graph, q≥1 theorems — and tells you the exact hypothesis (positive association) and the exact extension to verify (inhomogeneous comparison). Want me to draft the refined Lean skeleton with the axiom narrowed to just `comparison_inequality` + `edwards_sokal_twopoint`, and the four plumbing lemmas as real (provable) obligations? That's the version where the `sorry` count actually reflects the math you'd cite versus the math you'd prove.