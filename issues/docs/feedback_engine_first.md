---
name: feedback_engine_first
description: "deliberate strategy — build engine/primitives from first principles to soundness FIRST, defer applications/consumers, circle back on corrections apps expose; do NOT push \"build a consumer to validate\" mid-engine"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: f4c32996-1e21-444a-b832-be1c336459d1
---

Azriel's standing strategy: **build the engine from first principles, make it sound, build applications later, circle back on any corrections the applications expose.** This is a STRATEGY, not procrastination — the project's goals cut across many disciplines and applications emerge as new capabilities are built, so locking primitives to a premature consumer is a known failure mode.

**Why:** he has repeatedly built into applications too quickly, then incurred painful rebuilding/reorganizing once "another consumer" turned up needing primitives exposed differently / more completely. Engine-first avoids that churn. He is also actively developing the underlying THEORY in parallel — sound bedrock makes it easier to propose meaningful applications without speculation.

**How to apply:** when he's mid-engine, do NOT recommend "wire up a consumer to validate it" or otherwise push application-ward — validate with oracles / brute-force / property checks instead (the engine's own correctness spec, not a downstream app). Sequence further engine/primitive depth (representatives, paper-faithful extensions, exact bounds) ahead of consumers. Surface integration goals when *he* raises them ([[user_looks_ahead]]), then defer. Consumers/applications come after the bedrock is sound, by his explicit call. Aligns with [[project_faithfulness]] and [[feedback_execute_settled_decisions]]; this is the *engine ≫ application* ordering made explicit. Incident: after Z5c I recommended wiring persistent-Mapper-over-T first; he corrected to "keep building engine (representatives, then §5 codim-one), discuss consumers later."
