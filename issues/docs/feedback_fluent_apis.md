---
name: feedback-fluent-apis
description: "When a static-method API has grown a long parameter list that branches across discrete pipeline stages, refactor to a fluent builder with one method per stage rather than adding more parameters"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ef8e192e-720a-4222-addb-e184edf1d3ba
---

When a static method's parameter list grows to encode discrete pipeline
stages (e.g. topology → filter → repair → refine → scale, each picked
by an enum), prefer refactoring to a fluent builder with one method per
stage over adding more parameters or splitting into more sibling static
methods.

**Why:** The user chose the fluent-compiler surface for the graph
construction refactor (5-stage pipeline) over keeping the static
`GraphBuilder.BuildResult` facade. Reason: the fluent shape "encourages
stage-level thinking" and matches the per-stage decision manifest the
auto-mode compiler emits. A 9-parameter static method hides the staged
nature; a `.WithTopology(...).WithFilter(...).WithRepair(...)` chain
makes the stages first-class.

**How to apply:**
- When a static API has ≥ 6 parameters where several are enum/strategy
  selectors for discrete pipeline phases, that's the signal to migrate
- Keep the old static surface as an `[Obsolete]` shim that internally
  constructs the equivalent fluent chain — gives callers a migration
  window, doesn't force a big-bang rename
- The fluent chain should make per-stage overrides explicit and let
  auto/diagnostic modes drive defaults when a stage isn't overridden
- Persist per-stage decisions in a manifest so auto/explicit choices
  remain auditable end-to-end (the SPC/Mapper output ships with the
  graph-construction manifest as a sibling artifact)

This pairs with [[feedback_design_before_build]] for project layout
expectations.
