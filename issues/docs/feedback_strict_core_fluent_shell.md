---
name: feedback-strict-core-fluent-shell
description: "For multi-stage pipelines with order-dependent interrupts and overrides, prefer a pure declarative DTO consumed by the backend; put any fluent / chained API at the CLI / REPL boundary where it translates user intent into the DTO"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ef8e192e-720a-4222-addb-e184edf1d3ba
---

When building a multi-stage pipeline backend in C# (graphs construction,
clustering, mapper, etc.) that has interdependent logic — interrupts,
overrides, auto-picks, conditional skips — keep the **backend purely
declarative**: a single deeply-immutable `XxxConfig` record passed to a
`Build(config)` entry point. The engine looks at the whole config
holistically before executing.

**Never put a fluent builder in the backend.** The fluent API lives at
the CLI / REPL boundary, where it translates messy human input into the
clean config DTO.

**Why** (from the user, condensed):

- A fluent chain implies order: `.WithGlobalBandwidth(0.5).WithLmpRescaling()`
  reads like both happen in sequence. But LMP fundamentally renders
  global bandwidth dead — the fluent API misrepresents the real semantics.
- The user's Matlab-style instinct: "I want to define the state of the
  universe and hit run." Declarative config matches that mental model.
- 1:1 reproducibility is trivial when the manifest IS the config record.
  Serialize the DTO, deserialize later, get the same manifold.
- The fluent CLI layer can validate / route / warn / throw on
  catastrophic input before the engine ever sees it. Engine stays a
  deterministic mathematical kernel.
- Python `notebook/mvp/bench.py` is already structured this way —
  config-driven, CLI translates to args, engine runs.

**How to apply:**

- Top-level config: one `XxxConfig` record with `init`-only properties,
  `required` for the few non-defaultable fields, all else nullable or
  with sensible defaults
- Per-stage sub-configs: nested records carrying strategy enum (`Auto`,
  `Explicit1`, `Explicit2`, ...) + that strategy's parameters
- Backend entry point: `public static XxxResult Build(XxxConfig config,
  ...other immutable inputs like distance functions)` — no instance
  state, no method chaining
- Engine evaluates the whole config first, resolves interrupts /
  auto-picks / fallbacks internally, then executes the linear pipeline
- CLI fluent layer is its own class living under the REPL — e.g.
  `GraphCommandBuilder` in `UserRepl.Commands` — that has `.UseKnn()`
  etc. and emits the strict config via `.BuildConfiguration()`
- Manifest emission is straightforward: snapshot the config + diagnostic
  reports + automation log (which auto-picks fired and why)

**What this rules out:**

- Builder patterns that mutate internal state across method calls in
  the backend
- Sequence-dependent semantics in the backend API
- Splitting "fluent" and "explicit" entry points across the backend —
  there's only one, the declarative `Build(config)`

This pairs with [[feedback_no_compat_bandaids]] (clean breaks at the
end-state) and [[feedback_defer_deletes_during_refactor]] (keep old
shapes alive during the rebuild, sweep at the end).
