# Brief — P1a: residual-prior conditioning (the prior generalizes, monotone δ)

**Status:** build ticket — **local `issues/ph` (transient; archive to MarkBrain when P1 closes).** First
increment of the synthesis §6 **P1**. **Scope in one line:** lift P0's *fixed-weight* content edges to an
**edge prior + residual admission** `r_ij = |t_j − (t_i + τ_ij)| < δ` — undirected, monotone, read through P0's
existing barcode path — proving Li's residual band subsumes P0's fixed-similarity the way P0 proved it subsumes
SIFTS. **No Δ, no directedness, no zigzag, no gauge.**

**Background (MarkBrain — canonical design records):**
- synthesis §2 (Li construction, neurons stripped) + §6 P1 — [opus-brief-conditioned-persistence-synthesis.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/opus-brief-conditioned-persistence-synthesis.md)
- the phase this builds on — [p0-conditioned-filtration-brief.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/p0-conditioned-filtration-brief.md)
- the Φ-half seed — [backbone-conditioned-persistence.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/persistent-homology/backbone-conditioned-persistence.md)
- Li 2025 — arXiv 2508.11646 (`codex-scientiae: compendia/ph/2508.11646v1.md`), read as PH groundwork (Defs 1–2).

**Sequencing:** [dev-sequencing.md](dev-sequencing.md) — P1a is the main-line head after P0 (`9e3aa71`).

---

## 1. The one new idea — residual weight against a prior

P0 admitted a content (similarity) edge at its raw distance `d`. P1 admits it at its **residual against a
prior**:

- an observation field `{t_i}` (scalar to start; vector + metric later, like P0's metric swap),
- an edge prior `τ_ij` — the prediction that `t_j ≈ t_i + τ_ij`,
- the content edge `(i,j)` carries filtration value `r_ij = |t_j − (t_i + τ_ij)|` (undirected: symmetrized, §3),
- sweeping the slack `δ` admits edges in **residual order** — that sweep *is* the δ-filtration.

A **return relative to the prior** is a chord that closes a loop; its **birth is its residual**. The prior sets
how *surprising* a return is: a predicted return has residual → 0 (born early / cheap), an unpredicted one has
high residual (born late — a persistent surprise). SIFTS / Li / thermal are three choices of `(t, τ)`, one
mechanism.

## 2. The load-bearing distinction (caught in scoping — don't skip)

**The fixed backbone stays a separate anchor; the *residual band generalizes the content edges only.***
- P0's **backbone** (reading-order path at ε₀) → still passed as the fixed prior `K₀` at ε₀, unchanged.
- **`τ≡0`** does **not** reproduce P0's backbone — it reproduces P0's **similarity**: `r_ij = |t_j − t_i|`, the
  raw distance in `t`-space. (This is exactly the synthesis's "SIFTS = reading-order backbone **+** `τ≡0`
  content" — the backbone is not the τ≡0 part.)
- So P1a keeps `ConditionedFiltration`'s two-argument shape: `backbone` (fixed K₀) + `content` (now
  residual-weighted). A naive "τ≡0 recovers the whole P0 fixture" is wrong — τ≡0 recovers the *chord*, the
  backbone is still handed in fixed.

## 3. Why it's reuse-maximal — a producer feeding P0 unchanged

The residual is just a **weight**. P1a is a thin **upstream producer** —
`(observations, prior) → residual-weighted content edges (i,j,r)` — feeding the **existing**
`ConditionedFiltration.BuildGraph` / `ComputeBarcode`. Union convention, Rips, reducer, barcode: all P0. **The
only new code is the residual computation**, mirroring how P0's only new code was the union.

- **Placement:** `src/tda/ph/ResidualPrior.cs`, namespace `TDA.Ph` (companion to `ConditionedFiltration`).
- **Vocabulary:** functional identifiers (`residual`, `prior`, `observations`, `predicted`); the `t / τ / δ`
  symbols stay in prose and docstrings, not the type system (project terminology discipline).

## 4. The one design knob — undirected symmetrization

Li's residual is **directed** (`i→j`: predict `j` from `i`). P1a is undirected (directedness is P2), so it must
fold `r_ij` and `r_ji`:
- **`min`** — admit as soon as *either* direction's prediction is satisfied. **Recommended:** matches "the loop
  closes when the earliest supporting prediction enters," keeps the filtration monotone and clean.
- `max` (both directions must agree) / `mean` (average deviation) — record as alternatives; defer unless a test
  motivates one.

`τ≡0` makes `r_ij = |t_j − t_i|` already symmetric, so the knob is **inert on the subsumption test** — it cannot
disturb the P0-recovery claim.

## 5. API (shape, not gospel)

```csharp
public enum ResidualSymmetry { Min, Max, Mean }

public static class ResidualPrior
{
    // Residual-weighted content edges from an observation field + edge prior.
    // Each prior entry (i, j, tau) is a candidate edge; r_ij = |t_j - (t_i + tau)|, symmetrized.
    public static IReadOnlyList<(int i, int j, double r)> ResidualEdges(
        double[] observations,                              // {t_i}
        IReadOnlyList<(int i, int j, double tau)> prior,    // predicted t_j = t_i + tau
        ResidualSymmetry symmetry = ResidualSymmetry.Min);
}
```

Use: `ConditionedFiltration.ComputeBarcode(n, backbone, ResidualPrior.ResidualEdges(t, prior))`.

## 6. Tests / exit (`TDA.Ph.Tests`)

1. **`τ≡0` subsumes P0's content (the load-bearing claim).** Backbone = P0's path at ε₀ (unchanged); content =
   `ResidualEdges` with `τ≡0` and observations placing the chord `(0,5)` at residual `d`. Barcode = **P0's exact
   result** — one essential H₁ born at `d`. τ≡0 residual = raw distance = P0's similarity; backbone untouched.
2. **The prior shifts the return's birth (the P1 payoff P0 can't express).** Same backbone + one content chord.
   Prior A (`τ≡0`, unpredicted) → chord residual = its raw distance → return born **late** (persistent surprise).
   Prior B (predicts the chord, `τ` = the observed gap) → residual → 0 → the **same** return born at **0** (the
   prior expected it). Same loop (`β₁ = 1` both), **birth/persistence set by the prior**.
3. **Residual ordering = the δ-filtration.** Two content chords at different residuals are born in residual order;
   each bar's birth equals its residual.
4. **Green** in `TDA.Ph.Tests`; barcode routes through the existing reducer unchanged — no new reduction code.

## 7. Out of scope (P1b and beyond — do not creep)

- **`Δ` reach axis** (the prior-magnitude bound `τ_ij ≤ Δ`, the second filtration parameter) → **P1b**.
- **Non-monotone / zigzag slices** (trajectories, sliding windows, `T`-sweeps) reusing the built zigzag stack →
  **P1b** ([zigzag-frontier.md](../../../MarkBrain/ThermoMapper/issues/tda-purification/zigzag-engine/zigzag-frontier.md)).
- **Full multiparameter module** (`K_{δ,Δ}` rank invariant) → **Z6**.
- **Directedness / magnetic flux / directed flag** → **P2**. Gauge `q`, sheaf, `λ_q` → **P3–P4**.

**P1a = symmetric residual band, monotone δ, read as ordinary PH.**

> If P1a needs more than `ResidualEdges` + the P0 call + a test file, something is being rebuilt — the residual
> is a weight, and the engine already eats weights.
