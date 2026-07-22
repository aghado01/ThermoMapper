---
name: feedback_docstrings_lean
description: "docstrings = lean what-it-does + credit + key fact (critiques/derivations/prior-art → design notes); ALSO weave the function-algebra vocabulary as *locating narrative* + cref discipline so the parallel structure is legible"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 164b2936-306f-4ef1-bfdf-52132c232a52
---

Docstrings must be lean operational references — what the code does, source attribution, and the one or two behavioral facts that matter. Do NOT put academic critique, full derivations, lemma proofs, "this paper was wrong" commentary, or prior-art quibbling in a docstring.

**Why:** a docstring serves the next reader of the code, not a reviewer of the literature. Mixing in critique/derivation bloats it and misplaces the content — "doc strings should be doc strings, not academic reviews."

**How to apply:** operational summary + attribution (e.g. "mean-field formulation due to Wang 2020") + the single functional fact that matters (e.g. "computes the closed form instead of Monte Carlo"). Park reductions, lemmas, reduction-to-single-linkage critiques, and prior-art in `ThermoMapper/issues/` design notes. Concrete instance: [[project_wang2020_spc]] PKWang docstring stays lean; the single-linkage reduction and MC critique live in the spc-samplers plan appendix.

**Surface the parallel structure (the enrichment).** Beyond bare what-it-does, weave the project's common
function-algebra vocabulary ([[project_unifying_vocabulary]]: *field / accumulate / reduce / observable /
currency / sampler / solver / affinity / form-degree*) as **locating narrative** — one clause showing where
the thing sits in the recurring shape, e.g. *"accumulates the per-draw bond events and reduces them to the
`Affinities` currency."* The vocabulary is load-bearing prose (it *locates*), not ornament: it's what lets a
reader scanning docstrings **see** the parallel structure Azriel takes care to retain — the same verbs/nouns
recurring across tiers. Keep it to a clause, not a paragraph (still lean).

**cref discipline.** `<see cref="..."/>` the vocabulary-nouns that are real types (`Affinities`, `Accumulator`,
`IReduction`, …) so the parallel structure is *navigable*, not just narrated. The two rules compose: lean
exposition + a locating clause in the shared vocabulary + crefs on the real nouns.
