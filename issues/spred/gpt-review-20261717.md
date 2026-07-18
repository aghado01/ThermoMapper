Overall, the recent work is strong: SPRED is now a coherent research implementation with sound layering, useful diagnostics, scalable execution, and credible component-level validation. It is not yet benchmark-mature—the post-redesign ISOLET gate remains open, and I found two material semantic/integration gaps to fix first.
Review findings
[P1] Equal-count essential bars ignore their birth locations.
[DiagramMetrics.cs (line 65)](D:/aghado01/ThermoMapper/src/tda/ph/DiagramMetrics.cs:65) counts essential bars, removes them from matching, and only charges for count surplus. Consequently, diagrams containing essential bars (1, ∞) and (100, ∞) have distance zero if their counts agree. Sliced and Sinkhorn inherit the same policy.
This particularly matters after MinPersistence pruning: the cylinder report says pruning leaves exactly one essential H1 bar, so the H1 objective becomes a Betti-count gate rather than persistence-scale fidelity. That materially qualifies the “full persistence” and “near-exact” language in [dev-sequence.md (line 149)](D:/aghado01/ThermoMapper/issues/spred/dev-sequence.md:149).
I would either match equal-count essential bars by birth coordinate and penalize only surplus bars, or explicitly name/document this as count-only semantics. Add a same-count/different-birth test across all three backends.

[P1] The new annealer controls stop at the engine boundary.
SubspaceAnnealer.Compute accepts SubspaceAnnealerOptions at [SubspaceAnnealer.cs (line 51)](D:/aghado01/ThermoMapper/src/maths/geometry/dim-reduction/SubspaceAnnealer.cs:51), but neither [Spred.Compute (line 28)](D:/aghado01/ThermoMapper/src/tda/dim-reduction/Spred.cs:28) nor either public [DistributedSpred (line 72)](D:/aghado01/ThermoMapper/src/tda/dim-reduction/DistributedSpred.cs:72) entry point exposes them.
That blocks the prescribed S0 re-probe: the mobility report explicitly says to calibrate InitialTemperature, while the ISOLET pilot calls DistributedSpred. Today that requires bypassing the public driver. Threading the pure options DTO through these APIs preserves the layering perfectly.

[P2] Sinkhorn’s forbidden assignment cells are not actually forbidden.
BuildCost represents impossible diagonal matches with a finite big sentinel at [DiagramMetrics.cs (line 271)](D:/aghado01/ThermoMapper/src/tda/ph/DiagramMetrics.cs:271). SinkhornAssignment then includes every cell in its log-sum-exp calculations, relying on those weights to underflow at [DiagramMetrics.cs (line 292)](D:/aghado01/ThermoMapper/src/tda/ph/DiagramMetrics.cs:292).
This holds in the tested small-ε regime, but the public API accepts any finite positive ε. With a sufficiently large ε, forbidden cells receive substantial mass and the method no longer solves transport over the declared diagonal-augmented geometry. Mask sentinel cells as log-kernel −∞ and omit them from the final transport sum.

[P2] Numerical validation remains incomplete around the newly hardened surfaces.
Examples:
SubspaceAnnealerOptions relational checks do not reject NaN or positive infinity at [SubspaceAnnealer.cs (line 149)](D:/aghado01/ThermoMapper/src/maths/geometry/dim-reduction/SubspaceAnnealer.cs:149).
VarianceRegularizer is deliberately signed, but should still be finite; [PersistenceObjective.ValidateConfig (line 56)](D:/aghado01/ThermoMapper/src/tda/dim-reduction/PersistenceObjective.cs:56) does not check it.
EssentialPolicy.FinitePenalty accepts negative or non-finite values at [DiagramMetrics.cs (line 46)](D:/aghado01/ThermoMapper/src/tda/ph/DiagramMetrics.cs:46).
ValidateP allows NaN at [DiagramMetrics.cs (line 226)](D:/aghado01/ThermoMapper/src/tda/ph/DiagramMetrics.cs:226).
An invalid DiagramDistanceKind silently falls back to exact Wasserstein instead of failing.
These are inexpensive hardening fixes, but important because non-finite objectives can silently poison the annealer.

[P3] The documentation has fallen behind the implementation again.
[design.md (line 46)](D:/aghado01/ThermoMapper/issues/spred/design.md:46) still specifies the retired isotropic, cooling-tied proposal.
[annealer-mobility-brief.md (line 3)](D:/aghado01/ThermoMapper/issues/spred/annealer-mobility-brief.md:3) says “ready to execute” despite its implementation report.
[dev-sequence.md (line 179)](D:/aghado01/ThermoMapper/issues/spred/dev-sequence.md:179) still calls P3 active/executing rather than “implementation landed; S0 validation pending.”
The historical probe in [isolet-benchmark-brief.md (line 397)](D:/aghado01/ThermoMapper/issues/spred/isolet-benchmark-brief.md:397) is dated July 17 despite the July 16 commit and describes the old annealer in present tense.
The PH sequence still says P1b’s commit is pending at [dev-sequencing.md (line 45)](D:/aghado01/ThermoMapper/issues/ph/dev-sequencing.md:45), although it landed.
issues/spred/spred-retrospective.md remains untracked.

SPRED arc status
Area Current state
Architecture Sound and appropriately layered
Generic search engine Two-plane Grassmann moves, adaptive per-column scales, deterministic RNG, cancellation
PH objective Weighted H0/H1, graph-restricted Rips, exact/sliced/Sinkhorn distances, pruning, variance term
P0 profiling Closed; diagram matching—not graph construction or PH reduction—was correctly identified as the wall
P1 evaluation cost Largely closed through pruning and sliced screening; semantic issues above remain
P2 problem-size reduction Deferred: landmarks/subsampling/witness constructions
P3 mobility Redesign implemented and synthetic high-codimension fact passes; real S0 re-probe not run
P4 distributed SPRED Substantially complete: bounded concurrency, deterministic block seeds, Grassmann median, cancellation and diagnostics
Scientific validation Cylinder succeeds qualitatively; ISOLET feasibility established, effectiveness not yet demonstrated

The most important distinction is that the mobility redesign has passed an engine-level synthetic test—two-plane proposals descend where isotropic proposals remain frozen—but has not yet passed the actual ISOLET S0 gate. The remaining concerns in the mobility report are real: a flat PCA eigentail may require paired rank-2 moves, and the default temperature may melt a good warm start.
There is therefore no evidence yet that SPRED beats PCA or random projections on ISOLET, or that its projection improves downstream SPC. That is the boundary between “working research implementation” and “demonstrated method.”
Place in the broader PH stack
The placement is exactly right:
Maths.Geometry/SubspaceAnnealer
↑
TDA.DimReduction/PersistenceObjective + Spred
↙ ↘
Graphs.Proximity TDA.Ph
SPRED is not another PH reducer. It is the first substantial consumer that turns PH output into an optimization signal. That is an important maturation step for the repository: topology is becoming something the system can optimize against, not merely compute afterward.
The underlying PH platform is already broad and comparatively mature. The complete local PH suite passed 328/328 and includes ordinary persistence/cohomology reducers, graph and full-Rips H0/H1, cycle representatives, vineyards, multiple zigzag engines, collapse machinery, and nerves.
The active conditioned-persistence arc is currently:
P0 ConditionedFiltration: landed.
P1a residual prior: landed.
P1b reach-axis slices: landed.
P1c conditioned non-monotone sweep into graph zigzag: scoped, not implemented.
Filled-H1 zigzag and directed/flux work: later gates.
Full Rips is likewise at a clear intermediate point: threshold-bounded complete Euclidean H0/H1 is available, while alternate metrics, enclosing-radius defaults, higher-simplex enumeration, and dense performance routing remain open.
SPRED currently consumes the stable ordinary-PH path. Future conditioned, zigzag, or directed-topology objectives should become additional objective providers in TDA.DimReduction; none of that belongs in the generic annealer.
Recommended next sequence
Fix essential-bar semantics and expose SubspaceAnnealerOptions through both SPRED drivers.
Mask Sinkhorn’s forbidden support and close the finite-value validation holes.
Reconcile the SPRED and PH status documents.
Run the S0 eigengap pre-check, calibrate temperature from measured proposal deltas, then rerun the warm-start versus redesigned annealer comparison.
If S0 still crawls, implement the paired rank-2 proposal before spending on the full ISOLET screen.
Only then run the PCA/random/SPRED and downstream-SPC benchmark matrix.
For the visualization work, SPRED is already a promising first-class client: projected/reference barcodes, per-dimension objective terms, block subspaces and aggregation diagnostics are all useful visual currencies. The missing piece is a neutral per-iteration trace/observer contract—owned by the engine or consumer layer and merely consumed by viz.
Validation performed here: 12 focused geometry tests, 66 non-benchmark SPRED tests, and 328 PH tests all passed. I did not rerun the opt-in ISOLET pilots or the expensive external R transport oracle. No source files were changed.

---

## Parent-thread disposition (2026-07-16)

Every finding independently verified against source before acting — **all six confirmed, zero false
positives.** F1/F3/F4 are latent or boundary issues, not live bugs at current defaults, but all real.

| # | Sev | Verified | Disposition |
|---|---|---|---|
| Essential bars matched by count, births ignored | P1 | `DiagramMetrics` charges only `surplus·perBar^p`; equal-count essentials never compare births — a preserved H1 loop scores 0 at any scale. Material: makes the pruned H1 objective a Betti gate. | **Chip A** — match equal-count essentials by birth coordinate (all 3 backends), penalize only surplus; same-count/different-birth test ×3. |
| Annealer options stop at engine boundary | P1 | `Spred.Compute`/`DistributedSpred` expose no `SubspaceAnnealerOptions`; S0 re-probe can't set `InitialTemperature` via the public driver. | **Chip B** — thread the pure options DTO through both drivers. **Prerequisite for the pilot-lane S0 re-probe; collides with the pilot lane on `DistributedSpred.cs` — sequence deliberately.** |
| Sinkhorn forbidden cells not masked | P2 | LSE loops include forbidden cells; only the final sum is guarded at −700. Safe at ε=0.01, breaks at large ε + small diagrams. | **Chip A** — mask sentinel cells as log-kernel −∞ in the iterations, omit from the transport sum. |
| Finite-value validation holes | P2 | `ValidateOptions` (NaN/+∞ slip past relational checks), `FinitePenalty`, `ValidateP` (NaN), `VarianceRegularizer` (finite-but-signed), `_ => Wasserstein` enum fallback — all unguarded. | **Chip A** (`ValidateP`, `FinitePenalty`) + **Chip B** (`ValidateOptions`, `VarianceRegularizer`, reject invalid `DiagramDistanceKind`). |
| Docs behind implementation | P3 | SPRED docs confirmed stale. | **Done in-thread:** design.md proposal section, mobility-brief status, dev-sequence P3 reconciled to two-plane/landed. **Flagged, not touched:** isolet-brief §probe present-tense/date (pilot lane); ph/dev-sequencing.md P1b "*Commit pending*" is stale (`66d6f24` landed) — PH track, flag to its owner; `spred-retrospective.md` untracked (user's call: commit or gitignore). |

**Chip split rationale (complementary real estate, per dispatch discipline):**
- **Chip A — PH diagram-distance semantics + hardening** (F1, F3, F4-`ValidateP`/`FinitePenalty`):
  `DiagramMetrics.cs` + `TDA.Ph.Tests` only. No SPRED-driver overlap, no pilot-lane collision — safe
  to fire anytime. Carries the one genuine correctness fix (F1).
- **Chip B — SPRED consumer surface: options threading + validation** (F2, F4-options/regularizer/kind):
  `SubspaceAnnealer.cs`, `Spred.cs`, `DistributedSpred.cs`, `PersistenceObjective.cs` + geometry/
  TDA.DimReduction tests. Touches `DistributedSpred.cs`, which the active pilot lane depends on, and
  is the S0 re-probe's prerequisite → **the pilot lane should own or gate it**, not a parallel chip.

Accepted forward-looking notes (no action now): SPRED as the first PH-consumer that turns PH into an
optimization signal (placement affirmed); a neutral per-iteration trace/observer contract for viz is
a real seed → parked to the viz track, not the annealer.

### Chip A result (2026-07-17)

F1/F3/F4 landed in `DiagramMetrics.cs` + `TDA.Ph.Tests` — 343/343 green, 15 new facts. Essential
bars now match by birth in all three backends (their only finite coordinate; death = ∞ on both
sides collapses the L∞ ground metric to |Δbirth|), solved with the existing KMN `MinAssignment` on
a birth-only matrix whose surplus rows/columns stay at zero — sorted-order pairing is exact for
equal counts but cannot choose which surplus bars go unmatched under `FinitePenalty`, so one
assignment path serves both, and essential counts are Betti-scale so O(k³) is free. Sliced adds the
term once, slice-independent (no finite death to project; the transport is already 1-D); Sinkhorn
adds it exactly, never smoothed — cross-backend equality is pinned by test. Surplus stays charged
at perBar^p and `InfiniteOnMismatch` still returns +∞ on count mismatch. Sinkhorn's sentinel cells
are now masked as log-kernel −∞ inside both LSE sweeps and skipped in the transport sum; a
forced-support fixture (two bars vs empty ⇒ the constrained plan is the identity at any ε) verifies
exact agreement at ε = 1 and ε = 10, where the unmasked kernel drifted to ≈ 138 / ≈ 323 against the
exact 101. `ValidateP` rejects NaN; `FinitePenalty` requires finite perBar ≥ 0 (zero deliberately
disables the surplus charge, negative would reward essential-count mismatch). One caveat vs the
chip brief: "close to exact W at large ε" is a theorem only on forced-support fixtures — on general
diagrams large ε smears mass across *admissible* cells (ordinary entropic bias), so the general
large-ε test asserts LP feasibility (finite, ≥ exact optimum) instead.

### Chip B result (2026-07-17, parent thread)

F2 + F4-remainder landed in `5f9fc2e` — executed in the parent thread per the "pilot lane should
own or gate it" call (the interrupted session had picked B up; this session finished it).
`Spred.Compute` and both `DistributedSpred` entry points take an optional `SubspaceAnnealerOptions`
forwarded to every engine call; the parameter sits **after** `maxDegreeOfParallelism` so the pilot
suite's positional call sites compile unchanged (no collision — `SpredIsoletPilotTests` untouched).
Validation moved onto `SubspaceAnnealerOptions.Validate()` (public, NaN-proof pattern checks) so
drivers reject a bad record *before* the reference-barcode build, with the offending property as
`ParamName`, never wrapped in block context. `ValidateConfig` rejects non-finite
`VarianceRegularizer` (negative stays legal — PCA-spirit) and undeclared `DiagramDistanceKind`;
the `_ => Wasserstein` silent fallback is now an explicit case plus a throwing arm. Forwarding
facts live on the cylinder fixture deliberately — only where the anneal moves is a dropped
parameter detectable (best-so-far tracking makes optimal-warm-start fixtures option-invariant).
Maths.Geometry.Tests 41/41; TDA.DimReduction.Tests 76/76. **The S0 re-probe is unblocked** —
`InitialTemperature` (mobility finding 2) is now reachable from the public drivers; that re-probe
stays with the pilot lane.
