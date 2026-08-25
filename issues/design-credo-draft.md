# Design Credo — candidate articulations for ratification

> Draft 2026-07-23, from the primer-governance session. Format follows the MSCE
> decision-guidance shape: each principle = **name** + credo line + *lean*
> (arbitration rule for when it bites) + anti-pattern it exists to prevent.
> Status: candidate pool. Strike, rewrite, or promote entries through use —
> this document is L2, not scripture.

## I. Epistemic posture — how claims earn trust

### 1. Evidence over authority (the replication gate)
A method from a paper is not integrated until its reported results reproduce in
my harness. Citation is a lead, not a license.
*Lean:* when replication is expensive, replicate the smallest result that would
falsify the method if it failed.
**Anti-pattern:** integration by citation — trusting the abstract, inheriting the bug.

### 2. External calibration (differential oracles)
Correctness is established by side-by-side comparison against independent
implementations and sources. Agreement builds trust; divergence localizes error.
*Lean:* prefer an imperfect external oracle over a perfect self-referential test —
testing my code against my own understanding of it proves only consistency.
**Anti-pattern:** closed-loop validation.

### 3. Triangulation
Consult multiple independent perspectives (models, sources, implementations);
I am the synthesis layer, not any single oracle.
*Lean:* diversity of angle beats redundancy of the same angle.
**Anti-pattern:** polling one source three times and calling it consensus.

### 4. Stratified verification
Each architectural layer is tested at its own level of abstraction: property and
unit tests at primitives, differential and replication tests at engines,
scenario tests at applications. Test layers grow with the work, and patterns
that emerge in testing get named and kept.
*Lean:* a bug caught at the wrong layer is a missing test at the right one.
**Anti-pattern:** testing everything through the top layer.

## II. Architectural values — where things live

### 5. Mechanism before policy (engines first)
Build the general mechanism from the problem's invariants before the specific
application. Engines compute; applications decide. Reusability is the default
posture, not an afterthought.
*Lean:* if a decision could differ between two plausible callers, it is policy —
lift it out of the engine.
**Anti-pattern:** use-case logic baked into computation.

### 6. Anticipatory generality (resolves engine-first vs. rule-of-three)
Generalize on **invariants** — structure the problem provably has — and trust
deep invariants to compound: they carry large adjacent possibles, and the
engine will be exapted for purposes invisible at design time (evidence:
ThermoMapper, SPC clustering → persistent homology, mapper, graph theory).
Engines are convex bets — bounded downside (build cost), unbounded upside
(option value on unforeseeable futures).
*Lean:* restraint governs features, not foundations. Build nothing for a
specific imagined caller; foreclose nothing the invariants already permit.
Surface area (features, options, configuration) still waits for recurrence.
**Anti-pattern:** speculative generality (features for imagined callers) — and
its dual, foreclosure (corner-cutting that seals off the space the mathematics
held open).

### 7. The sinking rule
Functionality lives at the lowest stratum where it remains fully general;
specificity rises, generality sinks. Dependencies point downward only.
*Lean:* when unsure, keep it high — code sinks easily once proven general, but
a primitive that turns out to be policy contaminates everything above it.
**Anti-pattern:** premature sinking (false generality in the foundation).

### 8. Architectural annealing (the tetris feeling)
A new use case that doesn't fit the current structure is a perturbation, not an
exception. Two lawful responses: a new engine crystallizes behind it, or an
existing engine re-anneals — restructures into a generalization where old and
new functionality both fit at lower total complexity.
*Lean:* re-anneal when the new case shares invariants with an existing engine;
crystallize a new engine when the invariants differ. Never force the piece by
deforming the application layer around it.
**Anti-pattern:** bolt-on accretion — each addition locally cheap, the structure
globally ratcheting toward illegibility.

## III. The loop — how work proceeds

### 9. Human-gated crystallization
Lessons, conventions, and this document itself follow the pipeline:
capture → candidate pool → ratification → primer. Automation handles the floor
(capture, recurrence-counting); judgment holds the ceiling (promotion).
*Lean:* nothing becomes a standing rule on one occurrence.
**Anti-pattern:** skill-ifying a one-off; equally, re-deriving a lesson that
already recurred three times.

### 10. The deployment gate (holds are lawful)
Error costs are asymmetric across artifact classes. Object-level artifacts
(code, tests) fail locally and cheaply — draft them loose and sharpen through
use. Governing artifacts (primers, credos, constitutions) are upstream
amplifiers: an embedded flaw propagates into all downstream work and warps the
very judgment used to evaluate it. So the pipeline splits: draft freely into
the **candidate pool** (inert — safe to be wrong in; this is where the
blank-page problem gets solved), but **deploy** into active priming only at
satisfaction. "Hold" is a lawful verdict at the gate, distinct from paralysis:
an unwritten rule costs linear re-derivation; a wrong deployed rule costs
compounding drift.
*Lean:* when unsure whether an articulation is ready, leave it in the pool and
let another episode of evidence accumulate.
**Anti-pattern:** deploying guidance to feel finished; equally, mistaking a
deliberate hold for procrastination and forcing a draft live.
*(Boundary evidence: the ~2024 attempt where prematurely deployed guidance
took projects off the rails.)*

### 11. Full relaxation (no quenching)
A perturbation's reconfiguration ripples across dimensions that cannot all be
seen at once. The work is to follow the cascade — adjust, discover the next
strain, adjust again — recursively until the ripples quiet. In annealing
terms: cool slowly. Stopping mid-relaxation is quenching: it freezes defects
into the lattice as hidden strain (half-migrated seams, quiet inconsistencies)
that later perturbations amplify. The god-file is this principle's inverse —
accretion without relaxation.
*Lean:* the stubbornness is load-bearing. Accept a perturbation only with the
full cascade budgeted, or defer it entirely (a lawful hold under §10); never
accept it and stop halfway.
**Anti-pattern:** the half-annealed refactor — worse than none, because the
strain is now invisible.

---

*Standard vocabulary pointers, for searching the literature:* differential /
oracle-based testing, N-version programming, mechanism–policy separation
(classic systems), stratified design (Abelson & Sussman), stable-dependencies
principle, rule of three / refactoring to generality (Fowler), YAGNI vs.
first-principles design, adjacent possible (Kauffman), exaptation
(Gould & Vrba), real options / convex payoffs, annealing vs. quenching,
relaxation cascades / self-organized criticality (Bak).
