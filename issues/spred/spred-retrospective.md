Convention matches (d×k column-major both sides). Now the tests — the distributed suite is the heart of the recent work.

Grassmann exp/log are proper Edelman forms with an orthonormalize guard — sound. Now the remaining test files and the annealer test additions.

Two more checks: static mutable state in the graph/PH pipeline (the shell-lift thread-safety concern from the `Parallel.For` wrap), then the ISOLET brief.

Review complete — I read all the SPRED source ([Spred.cs](src/tda/dim-reduction/Spred.cs), [DistributedSpred.cs](src/tda/dim-reduction/DistributedSpred.cs), [PersistenceObjective.cs](src/tda/dim-reduction/PersistenceObjective.cs), [PersistenceObjectiveConfig.cs](src/tda/dim-reduction/PersistenceObjectiveConfig.cs), [SubspaceAnnealer.cs](src/maths/geometry/dim-reduction/SubspaceAnnealer.cs)), all five test files, the aggregation primitives they lean on (Grassmann manifold, IRLS/GeometricMedian), and the three issues docs.

**Overall verdict: this is good work.** The architecture matches the design docs faithfully, the determinism and cancellation plumbing is done properly, and I specifically verified several load-bearing details that could have been silent landmines and were not. No correctness-critical bug found. The findings below are: two things I verified safe (worth knowing), four code issues (real but non-catastrophic), five test blind spots, and two doc-level gaps — roughly in priority order within each group.

## Verified sound (the scary ones)

1. **The Grassmann median genuinely iterates.** [DistributedSpred.cs:319](src/tda/dim-reduction/DistributedSpred.cs:319) calls `GeometricMedian.Compute` with `opts = default`, and `default(IrlsOptions)` is all-zeros — which would have meant _zero iterations_, making every robustness test vacuously return its warm start (which is always the clean block in your fixtures). But [irls.cs:45](src/maths/geometry/solver/irls.cs:45) has a sentinel: `MaxIterations == 0` → `IrlsOptions.Default` (200 iterations, OptimalityCheck singularity policy). So the aggregation tests are not vacuous on that axis.
2. **Conventions and invariants line up end-to-end.** The k×d row-major projection ↔ d×k column-major Grassmann frame packing is consistent in `ProjectionToFrame`/`FrameToProjection` and the tests' `PackFrame`; the Grassmann exp/log/distance are proper Edelman forms with an orthonormalize guard, so the aggregate stays orthonormal; `Prune` keeps essential bars because `Bar.Persistence` is `+∞` for them; and the `maxDegreeOfParallelism = 1` default is genuinely load-bearing — the graph pipeline has its own unbounded inner `Parallel.For`s (DirectedKnn, EpsilonBall, PathNeighborRefiner), exactly as the dev-sequence claims. The `AggregateProjections_DuplicateSubspaceWinsMedian` test is a nice touch — it validates the Grassmann quotient invariance (rotated frame = same point), which is the whole design premise.

## Code findings

1. **Validation admits blocks the pipeline can't survive, and the failure surface differs serial vs parallel.** `ValidateInputs` allows `blockCount` up to `data.Length` ([DistributedSpred.cs:231](src/tda/dim-reduction/DistributedSpred.cs:231)), but a block needs enough rows for the kNN recipe (K+1) and a meaningful PCA warm start. A too-small block throws from the `PersistenceObjective` _constructor_ — `GraphPathologyException` is caught only inside `Evaluate` ([PersistenceObjective.cs:58](src/tda/dim-reduction/PersistenceObjective.cs:58)), not around reference-barcode construction ([PersistenceObjective.cs:42](src/tda/dim-reduction/PersistenceObjective.cs:42)). On the serial path that escapes raw; on the parallel path `Parallel.For` wraps it in `AggregateException`. So "serial == parallel" holds for results but not failure modes. Either tighten validation (or document the block-size floor) and normalize the exception surface in `RunBlocks`.
2. **Redundant full objective evaluations in the diagnostics path.** The single-block path evaluates the _identical_ projection three times — `LocalObjective` in `RunBlock`, then `AggregateObjective` and the full-data objective at [DistributedSpred.cs:126-132](src/tda/dim-reduction/DistributedSpred.cs:126) — each a complete graph+PH+Wasserstein pass (seconds at scale). One value serves all three. Relatedly, `RunBlock` re-evaluates the annealer's result ([DistributedSpred.cs:297](src/tda/dim-reduction/DistributedSpred.cs:297)) because `SubspaceAnnealer.Compute` doesn't return `bestValue` it already holds — having the engine return (projection, value) would save one full eval per block everywhere.
3. **No config validation on `PersistenceObjectiveConfig`.** A `Dimensions` entry with `Dim >= MaxDimension` silently compares empty diagrams — the objective can become identically 0 and the SA optimizes nothing, with no error. Negative `MinPersistence` and non-positive `WassersteinOrder` are also unchecked. A cheap `Validate()` at objective construction would close this; it matters more now that the ISOLET brief has people hand-writing configs.
4. **Minor: block-seed aliasing.** `seed + 1009 * block` means (seed 1, block 1) and (seed 1010, block 0) share a stream. Deterministic and provenance-recorded — consistent with your RNG-hygiene stance — but a splitmix-style derivation would remove aliasing across runs. Also `Spred.Compute` accepts ragged rows (fails deep in `Project` with `IndexOutOfRangeException`) while `DistributedSpred` validates them — inconsistent surfaces.

## Test blind spots

1. **Most distributed tests never anneal.** Everything except `ParallelMatchesSeededSerialRun` (3 iters) and `AdversarialLowIteration` (6 iters) runs `maxIters: 0`, so the corrupted-block, tradeoff, and matches-global tests are really facts about _PCA-per-block + Grassmann median_. That's a legitimate fast-contract-test choice, but the names and the dev-sequence phrasing ("clean-majority/corrupted-block robustness fixture") imply the SPRED estimator was exercised. Worth a one-line comment in the fixtures stating what's actually under test.
2. **The median warm start is never stressed.** Aggregation initializes at `frames[0]`, and every fixture puts clean blocks first — so the median always _starts at the right answer_. One permutation with the corrupted block at index 0 would actually test Weiszfeld convergence away from the initialization. Cheap and high-value.
3. No test for the small-block degeneracy or the `AggregateException` asymmetry (finding 1).
4. Mid-run cancellation is tested at the engine level, but the parallel `Parallel.For` cancellation path only sees a pre-canceled token.
5. **Thread-safety rests on one green run.** Serial==parallel equality on 3 blocks is real evidence, but given the `GpRegression.Fit` shared-scratch incident in this repo's history, the honest closure item for P4 is a static-mutable-scratch audit of the graphs/TDA.Ph stack (I checked for `[ThreadStatic]` — none — but didn't do a full static-field sweep).

## Docs

1. **The ISOLET brief's Phase 2 is compute-infeasible as specified.** The H0-only config prescribes `MinPersistence = 0`, keeping all ~n finite H0 bars, and the exact Hungarian Wasserstein is the O(n³) cost your own P0 profile measured (165 ms at n=200). Extrapolating to 8 blocks of ~975 rows: roughly 19 s _per evaluation per block_, times maxIters, times seeds, times three target dimensions. Phase 2 silently depends on the un-landed P1 approximation lever (Sinkhorn/sliced OT) or P2 subsampling — the brief should state that gate explicitly, since it otherwise reads as runnable today. (The brief's contiguous-block/speaker-order hazard, by contrast, is correctly identified, and the shuffled-partition deliverable is real work — the code has no shuffle option yet.)
2. **Memory hygiene.** The brief formally supersedes the "ISOLET PCA wall" interpretation, so I patched my `project_isolet_pca_wall` memory note with the correction (published target is 93% purity / 65% coverage + hierarchy, wall demoted to historical control).
3. Trivially: the cylinder test comment cites the paper's "(λ, 1−λ)" next to weights (0.5, 1.0) that don't sum to 1 — deliberate generalization, slightly misleading comment.

Nothing here blocks the P4 closure — the two things I'd act on first are the small-block/exception-surface issue (only sharp edge a caller can actually cut themselves on) and the corrupted-block-first median test (cheapest way to convert the robustness story from "plausible" to "demonstrated"). Happy to chip any of these off into a follow-up task if you want to delegate them.

All chips are queued and the in-session items are done. Here's the full disposition ledger:

## Memory (addressed first, per your correction)

- New memory [feedback_markbrain_flow_direction.md](C:\Users\azrie.claude\projects\D--aghado01-ThermoMapper\memory\feedback_markbrain_flow_direction.md): MarkBrain is repository/archive; project memories are the proximal, current record — edit directly, they migrate to the vault eventually.
- Corrected the MEMORY.md header itself (it said "vault canonical, do not hand-edit" — that inverted framing is what produced my bogus "will be overwritten on sync" caveat) and added the index line.
- Consequence: the `project_isolet_pca_wall` supersession edit from the review stands as-is; nothing to propagate.

## Edited here (doc-level observations)

- **ISOLET brief** — added an "H0 matching-cost gate" subsection after the H0-first config in [isolet-benchmark-brief.md](issues/spred/isolet-benchmark-brief.md): sizes the exact-Hungarian cost at block scale (~15–20 s/eval at 975 rows, extrapolated from your own P0 profile), states full-data exact H0 is infeasible, and gates Phase 2 on Sinkhorn/sliced OT, declared subsampling, or a pilot-derived wall-clock budget. Also added a matching row to the risks table.
- **dev-sequence** — the P1 "Approximate the Wasserstein" bullet now cross-references that gate in [dev-sequence.md](issues/spred/dev-sequence.md), so the two docs agree on what blocks the ISOLET screen.

## Chipped off (five chips pending)

1. **Harden DistributedSpred block-size and error surface** — the small-block/unguarded-constructor issue plus normalizing the serial-vs-parallel exception surface; includes the deterministic mid-run parallel cancellation test (via an internal `RunBlocks` seam rather than timing races).
2. **Stress Grassmann median with corrupted-first blocks** — the warm-start blind spot, plus the clarifying comments on which tests exercise annealing vs. PCA+median aggregation.
3. **Validate PersistenceObjective config and inputs** — construction-time `Validate()` (the silent `Dim >= MaxDimension` constant-objective trap being the headline), ragged-row consistency, and the cylinder `(λ, 1−λ)` comment nit as a ride-along.
4. **Deduplicate DistributedSpred diagnostic evaluations** — the single-block triple-eval fix, with the annealer-returns-best-value API change as a clearly-scoped optional stretch.
5. **Audit graph/PH pipeline static state for parallelism** — the shell-lift-style static-scratch sweep; if clean, it adds a cheap repeated serial-vs-parallel stress fact.

Two practical notes: chips 1, 2, and 5 all touch `DistributedSpredTests.cs` and chips 1 and 4 both touch `DistributedSpred.cs`, so run them one at a time (or merge deliberately) rather than concurrently. And the chip launcher defaults to a fresh worktree — given the long-path worktree failures recorded for this monorepo, start them in the main tree if the worktree option misbehaves.

## Remaining open item (your call, no chip)

**Block-seed derivation** — `seed + 1009 * block` aliases across runs (seed 1/block 1 ≡ seed 1010/block 0). It's deterministic and the derived seed is recorded in diagnostics, so provenance holds; a splitmix-style derivation would remove the aliasing but churns every seeded distributed test expectation. My recommendation: leave it and fold it into the already-deferred repo-wide RNG-hygiene audit, where it belongs with its siblings, rather than spending a churn cycle on it now.

That accounts for every finding from the review: two verified-safe items needed no action, four code findings and five test blind spots are covered by chips 1–5, both doc gaps are edited in place, and the seed question is the one deliberate deferral.

Let me look at the chip round-up, the current working-tree state, and your updated guidance.

The pilot-test diff came back empty — inconsistent with the earlier status. Checking whether it got staged or reverted in the meantime, and reading the two pool sites:

Tied off. The ledger is appended to [chips-open-ends.md](issues/spred/chips-open-ends.md) and committed as `be2b09f` (that file only, on top of the concurrent lane's `2c84cac`).

**The headline: nothing from the chip flurry was actually lost.** All five work products were re-landed by the reconciliation commits — `10d7939` carries chips 1–4 wholesale (error-surface hardening, corrupted-first tests with the maxIters:0 clarifying comments, the validation battery, `SubspaceAnnealerResult` dedup) and the real prize: the bug chip 2 exposed got a proper root-cause fix. The corrupted leading block sits at principal angle π/2 from the clean majority — the Grassmann **cut locus** — where the log map degenerates and Weiszfeld froze at the `frames[0]` initialization. The fix is `GeometricMedian.MedoidIndex` as the shared warm-start selector, with MoMPCA re-seated on it too (where the same stall was silently poisoning scale calibration). Suites went 47/49 → 57/57.

**What I closed this turn:**

- Chip 5's one unaudited residual — the `Spectral.cs`/`GraphLaplacian.cs` ArrayPool sites, reachable under `Parallel.For` via the `EffectiveResistance` filtration knob — verified by direct read: rent/return in `try/finally`, call-local buffers, no static state. The clean bill now covers that path.
- The round-up file is now the closed coordination ledger, per your updated brief-centric dispatch conventions.

**Genuinely open, and deliberately not touched here:**

1. **S0 annealer mobility** — `2c84cac` landed mid-turn from a concurrent session: the pilot is blocked on the SA's fixed schedule (`0.99^iter` cooling, `0.1·temp` step), not compute. That's dev-sequence P3 (adaptive cooling, larger steps, restart ensembles), and it's being actively probed in that lane, so I stayed off it.
2. **Block-seed derivation aliasing** — still parked for the repo-wide RNG-hygiene audit, unchanged.

One connective observation for whoever picks up P3: the step budget problem is geometric. On Gr(30, 617) the initial 0.1-radian geodesic step is tiny relative to typical inter-subspace distances (which scale like √k·(π/2) ≈ 8.6 radians there), so the anneal can only ever explore a small neighborhood of the PCA warm start — consistent with "blocked on mobility." A step schedule scaled to the manifold's diameter (or acceptance-rate-adaptive) is probably the first lever, ahead of restarts.

Now the commit (scoped to the three spred docs) and the single bundled chip:

Triage done, committed as `161117b`, one bundled chip issued. The disposition:

**Execute (chipped as one unit):** candidates 1+2 — two-plane Givens proposals with the acceptance-adaptive step as their controller. I pressure-tested the chip's recommendation before adopting it and it holds, with one sharpening now recorded in the brief: adaptive step _alone_ would degenerate (acceptance-targeting buys its 25% by shrinking steps toward a diffusion crawl), so the two-plane move carries the fix — it restores O(1)-dimensional moves whose improving fraction doesn't collapse at intrinsic dimension 17,610. One design point worth noting: the Givens move is itself a Grassmann geodesic (rank-1 horizontal tangent, single σ = θ), so the design.md "geodesic SA on Grassmann" deviation story survives intact — this is a proposal-distribution change, not a manifold-story change.

**Track:** gradient-informed Riemannian search (revisit only if 1+2 stall on S0), and restart-ensembles/plateau-early-stop as cheap later add-ons once the options surface exists.

What landed where:

- [annealer-mobility-brief.md](issues/spred/annealer-mobility-brief.md) — the settled spec: proposal design, `SubspaceAnnealerOptions` DTO, determinism contract, and the validation plan. The key test requirement is a synthetic high-codimension mobility fact (quadratic-on-Grassmann, old proposals stall / new ones descend) so mobility is provable in seconds without touching ISOLET.
- [dev-sequence.md](issues/spred/dev-sequence.md) — P3 re-pointed: "proposal-scale-bound, not budget-bound," with the probe numbers and the brief link.
- [chip-grassman-median.md](issues/spred/chip-grassman-median.md) — committed with a disposition section closing the loop, per the report-back convention.
- Chip **"Rebuild SubspaceAnnealer proposals per mobility brief"** — scoped to the annealer + geometry tests, explicitly fenced off `SpredIsoletPilotTests`/`pilot/` so it can't collide with the concurrent pilot lane; its report goes into the brief's Report-back section.

One alternative you have: the probe session that wrote chip-grassman-median.md offered to build this on your word and already holds the S0 context — if that session is still alive, relaying the brief to it instead of starting the chip is equally valid; the brief was written to be executable by either. The S0 acceptance re-probe stays with the pilot lane in either case.

Read it. Honest verdict: yes, there's one genuinely transferable idea in there for the annealing walk, but it's not where the doc points it — and it's worth being clear about both halves.

**What doesn't transfer.** The doc's native shape is _candidate-set pruning_: many discrete basins/clusters exist simultaneously, you score them cheaply, keep K, refine survivors. The annealer isn't a candidate-set problem — it's a sequential walk, and its S0 failure was a proposal-geometry defect (vanishing improving fraction of isotropic moves), which the two-plane redesign already attacks directly.

Top-K doesn't compete with that fix. Also, side note on provenance: that's a Perplexity riff whose citations are half keyword collisions — crochet cardigan patterns and "Cascade Series SPC" _vinyl flooring_ — so I'd treat it as an ideation prompt, not literature. The ideas that survive have proper names elsewhere.

**What does transfer: the cascade economics.** The core heuristic underneath cascade top-K is "never spend the expensive oracle on an unscreened candidate." Your annealer is now the textbook case for it: one PH evaluation costs ~0.092 s at S0 scale while generating a proposal costs microseconds — five to six orders of magnitude of asymmetry. The current design (and the redesign in the brief) spends the full oracle on _every_ proposal, one per iteration. Two concrete instantiations:

1. **Screened multi-proposal steps — best-of-M per iteration.** Generate M candidate two-plane moves, score them with a _cheap_ stage, and spend the exact evaluation only on the surrogate-best. The two-plane structure makes this unusually cheap: a Givens candidate (i, v, θ) changes only coordinate i of every projected point, so with the current projection cached, screening a candidate is O(n) after one O(nd) dot-product pass per fresh v. And the cheap stage doesn't need new theory — the fidelity ladder already exists and is oracle-validated: sliced/Sinkhorn `DistanceKind` on a landmark subsample.

The principled version of this is **multiple-try Metropolis** (Liu–Liang–Wong), which corrects the acceptance ratio if reversibility ever matters; for optimization-mode SA with best-so-far tracking, the greedy version is fine. One honest gate before committing: measure surrogate fidelity first (rank correlation between the cheap score and true Δobjective on the cylinder plus one S0 block) — if the cheap stage misranks, the cascade just adds cost.

2. **Fidelity-cascade restarts.** The doc's stage-0-coarse → top-K → refine pattern maps exactly onto the parked "restart ensembles" remainder of P3, upgraded to **successive halving / Hyperband shape**: many seeds × short budget × cheap fidelity (sliced + subsampled + coarse `MinPersistence`), keep top-K by objective, graduate survivors to exact Hungarian on full data. Zero new math — pure orchestration over knobs that already exist in the config surface.

**Where the doc's native idea actually belongs.** Basin-level top-K with persistence guardrails is really an SPC-side concept, and there it's more than an analogy — your ISOLET resolver work is already implicitly top-K-shaped (selective clustering _is_ a purity-coverage gate; lineage persistence _is_ a survival criterion across the temperature filtration). If cascade pruning earns a home in this codebase, I'd bet it's in the Track B resolver/hierarchy machinery before the annealer. The "thermal filtration" framing in the doc is a nice articulation of what your lineage resolver already does.

So: not a rival to the two-plane redesign, but a legitimate phase-2 layer on top of it (item 1), plus a concrete upgrade to the tracked restart-ensembles item (item 2). If you want, I'll fold both into the mobility brief's "Tracked, not executed" section with the surrogate-fidelity measurement as the explicit gate — say the word and it's a two-minute edit.
