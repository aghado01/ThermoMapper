---
name: Don't oversplit in design discussions
description: When user proposes a naming/structure pattern, match it as-is without introducing asymmetries or extra shapes they didn't ask for.
type: feedback
originSessionId: 20f9fe2d-c568-48c1-ba67-ecd6949b9fc8
---
When the user sketches a structural pattern (namespaces, file layouts, API shapes) and offers two examples to illustrate it, do NOT invent asymmetric treatments for the two examples just because their current implementations differ.

**Concrete instance (2026-04-18, ps.core.pwshspc):** User asked for `Spc.Synthetic.Generate.TwoMoons` and `Spc.Hashish.Compute.SimHash`. I responded with two *different* shapes — "verb-as-class, noun-as-method" for synthetic vs. "verb-as-namespace, noun-as-class" for hashish — justifying the split by pointing out that synthetic generators are pure functions while hash algorithms carry state. User called this "oversplitting" and stepped away to think.

**Why:** Design conversations benefit from a single unified pattern the user can evaluate. Introducing asymmetry multiplies the decisions they have to make and obscures the pattern they were sketching. Even if the implementation details suggest two shapes, propose the unified shape first and let the user surface the asymmetry if they want it.

**How to apply:** When the user gives examples of a structural pattern, treat the examples as hints about the *shape they want*, not as constraints that must be derived from existing implementation state. Propose the simplest unified rendering of that shape. If an implementation detail genuinely doesn't fit (e.g. stateful class vs static function), surface it as a single question ("Hashish classes carry state — OK with that under the same pattern?"), not as a fork in the proposal.
