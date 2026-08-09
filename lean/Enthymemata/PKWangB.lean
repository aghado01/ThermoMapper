import Mathlib.Data.Fintype.Basic
import Mathlib.Data.Real.Basic

/-!
# PKWang — the open obligations

The apologizing remainder of the PKWang formalization; the proved part lives
in `Lemmas.PKWang` (threshold reduction + deterministic cut graph).

* **Lemma A** — the MeanField global ranking is the doubly-pooled limit of
  the density-adaptive LocalField per-site ranking.
* **Lemma B, part 2** — single-linkage equivalence of the deterministic cut.

Protolemma and provenance: `lean/Protolemmata/spc-pkwang-lemmas.md`.
-/

namespace PKWang

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- Local field: cumulative energy at edge `(u, v)`, ranked ONLY against the
site-local neighborhood of `u` — the sum of `J u w` over all `w` with
`J u w ≤ J u v`. Obligation: the `Finset.sort`-based cumulative sum. -/
noncomputable def localHCum (J : V → V → ℝ) (u v : V) : ℝ := sorry

/-- Mean field: cumulative energy at edge `(u, v)`, ranked against the global
edge pool — the sum of `J x y` over all pairs `(x, y)` with `J x y ≤ J u v`.
Obligation: same machinery as `localHCum`, global index set. -/
noncomputable def globalHCum (J : V → V → ℝ) (u v : V) : ℝ := sorry

/-!
### TODO (Lemma A statement)

The pooling relationship `localHCum → globalHCum`: state the algebraic limit
that maps the LocalField arrays onto the MeanField array as local constraints
relax to the global pool. Statement deferred until the two defs above are
real — a sorried statement over sorried defs apologizes twice and says
nothing.

### TODO (Lemma B, part 2): single-linkage equivalence

The connected components of `Lemmas.PKWang.deterministicCutGraph` are exactly
the single-linkage clusters at energy threshold `T * log 2`. Needs a
formalized `singleLinkageClustering` to state against — likely via
`SimpleGraph.ConnectedComponent` on the threshold subgraph of the complete
weighted graph.
-/

end PKWang
