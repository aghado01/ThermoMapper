Drafting it into the lean dir alongside your existing lemma files, matching the convention (`.md` with embedded Lean + the axiom/obligation accounting).

Drafted to [bifiltration-bridge-lemmas.md](.discussion/issues/lean/bifiltration-bridge-lemmas.md). Here's the tour.

**The factoring delivered.** The file quarantines exactly the two Grimmett results as `axiom`s and discharges everything else:

- **A1 `comparison_inequality`** (Thm 3.21) — stated for a *vector* `p : V → V → ℝ`, deliberately, because that's the inhomogeneous version SPC's heterogeneous couplings need. Carries `1 ≤ q`.
- **A2 `edwards_sokal_twopoint`** (§11.2) — the single algebraic identity `⟨δ⟩ = 1/q + (1−1/q)·φ(i↔j)`.
- **P1–P5 + T1 are fully written proofs, not `sorry`** — the FK parameter's β-monotonicity (`Real.exp_le_exp`), the increasing-event combinatorics (`Reachable.mono`), the comparison fired on that event, the affine composition, and the cut-graph edge-set growth. The physics→edge-set spine closes.
- **T2/T3 are the only other `sorry`s, and they're construction obligations** in your TDA layer (flag-complex homology functor; reparam invariance), not Mathlib gaps — they wait on `FlagComplex`/`PersistentHomology` getting Lean definitions, then they're real proofs.

So the honest count is **4 `sorry`s, two kinds**: two you *cite* (A1, A2), two you *build* (T2, T3). Everything between physics and topology is green.

**Three design choices worth your eye:**

1. **Config = `SimpleGraph V` ordered by its own lattice `≤`.** This is what makes "increasing event" and "adding edges preserves connectivity" land on Mathlib's existing API (`Reachable.mono`) instead of a bespoke Bool-vector order. The connection event `Conn i j` and `Increasing` predicate fall out cleanly.

2. **`q > 1` vs `q ≥ 1` is tracked separately and on purpose** — the comparison inequality needs `q ≥ 1`, the affine slope `1−1/q > 0` needs `q > 1`. The proof of `spin_agreement_mono` actually uses both, and I split the hypotheses so the dependency is visible (and `q = 1` percolation is documented as the degenerate-but-valid case).

3. **The `IsRandomCluster` / `IsPottsTwoPoint` opaque specs** are the stand-ins for the unformalized measure-theoretic definitions — the axioms are stated *under* them, exactly the way Lemma B was stated under the closed-form bond probability. This is what lets you cite Grimmett without formalizing the random-cluster σ-algebra.

**Caveats for the proof pass** (where I'm guessing Mathlib's moving API): `Real.exp_le_exp` and `SimpleGraph.Reachable.mono` are the two names I'd verify first — the shapes are right but the exact identifiers drift across Mathlib versions. And `opaque ... : Prop` may want to be `axiom IsRandomCluster : (...) → Prop` instead if Lean complains about inhabitation — same effect, declares the predicate without a body.

**One thing I deliberately did *not* discharge:** the inhomogeneous-comparison check is *inside* A1's statement (the vector `p`), flagged in the notes — Grimmett writes 3.21 scalar, and confirming the per-edge version against his §3.1 conditional-probability machinery is the one citation-step that isn't quite turnkey. I'd treat A1 as "cited modulo that check" rather than fully closed.

Want me to tighten anything — narrow A1 to a cleaner per-edge statement, or stub the `PersistenceModule` interface so T2/T3 have a real codomain to prove against instead of the `True` placeholder?