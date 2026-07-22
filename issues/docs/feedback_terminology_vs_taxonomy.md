---
name: feedback-terminology-vs-taxonomy
description: "draws vocab promiscuously across fields by fitness (no anchor; cybernetician); math liked for generality + as reasoning scaffolding but NOT reified (functional → IReduction, no IFunctional); field terms first-class & often exact (observables, Betti)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 8b422e2b-4a68-4698-962d-25a53e4b8384
---

The user draws vocabulary and taxonomy **promiscuously across fields** — math, physics, topology, CS,
cybernetics — at discretion, picking the fittest / most-precise term per spot. He's a **cybernetician**,
fluid across boundaries by design; **no single field is the anchor.** This is a data-science project, not
a pure-math reframing — do **not** subordinate field terms to math.

What he likes about **math vocabulary** (functionals, maps, fields-over-indexes, scalar vs vector fields)
is its **generality** — posed correctly it's "not easily wrong" and survives domain crossings. But it's
valued as **reasoning scaffolding**, not a privileged anchor, and **not reified**: *"I like talking in
terms of functionals on the way to writing `IReduction` — but I don't want an `IFunctional`."* So: reason
in the general term, **build the concrete one** (named from whatever field fits best).

**Field terms are first-class and often exact — not lossy overlays.** *"observables"* (physics) was the
**first move** of this whole reorg (`graphs/observables`, `graphs/models/*/observables`) — a cornerstone,
not a metaphor. *"Betti"* (topology, `b₁ = E−V+C`) is a precise invariant that "shan't be ignored."
Likewise *currency*, *kernel*.

**Fitness test for a borrowed term** (his criterion, articulated via *observables*): **intuitive** +
**native to a serious community that studies the object** (physics for graphs — data-adjacent, and matching
*this* project's thermodynamic register: χ/Cv/M are literally physics observables) + **carries its weight
without loss of generality**. *Observable* is the ideal case — the physics term and the algebra's Observable
*kind* are the **same object** (`field → value`), zero gap, so the borrow loses nothing; it even spans the
deterministic (cyclomatic = reduce a partition) and expectation (χ = reduce an MC ensemble) cases unchanged.
The test cuts both ways: intuitive-but-narrowing (drags baggage) **or** general-but-inert both fail.

**Meta — I've over-corrected here twice; don't.** Don't over-systematize his terminology into a fixed
hierarchy (neither "naming is soft, defer it" nor "math is the anchor, field terms are shell"). He's
eclectic *by intent*. The one durable rule is **reason-in-general-terms vs reify-concrete-ones**
(functional → `IReduction`; no `IFunctional`). When naming defers (e.g. the genus), it's because it isn't
yet *forced* — not because naming is categorically soft. Relates to [[project_unifying_vocabulary]],
[[user_reads_for_architecture]], [[feedback_adversarial_partner]].
