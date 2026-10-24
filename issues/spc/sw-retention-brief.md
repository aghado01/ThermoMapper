# SW retention brief

Status: side brief, 2026-10-24. Parked for later and not scheduled. Written at the user's request, so that the Swendsen–Wang items of the BASS preparation brief's retention inventory (`issues/regression/bass-preparation-brief.md`, "Retention inventory", SPC samplers) have a home outside the regression track.

The principle, from the user: keep each draw's generative state, derive views from it, and where alternatives keep different information, offer both as options.

## Where things stand

- **The cycle** (`src/clustering/graphical/spc/runtime/core/sampler/SwendsenWang.cs`). Every sweep computes the spins, the FK bonds and their clusters (union-find), the energy, the magnetization from the largest colour, Σ|c|², and the largest cluster's size.
- **The accumulator** (`Accumulator.cs`). It keeps:
  - running sums of E, E², m, m², Σ|c|² and Σ|c|² without the largest cluster;
  - the log-binned cluster-size histogram;
  - the final spins and RNG words, for resume;
  - the optional per-edge counts (`AccumulationSpec`: affinities, alignments, co-membership) and per-node sums (the cluster-size and giant-participation landscapes).
- **The runner** (`SwRunner.cs`). It discards `RunBudget.BurnIn` cycles, then accumulates `Cycles`. A sweep keeps every run's accumulator, per temperature and replica (`SweepResult.SweepRuns`).
- **The serialization.** The accumulator is the resume state and is serialized (`SwRunSerialization.cs`, `AccumulatorSerializer`), so any retained series is part of the binary format and has to append across a resume.

## Items

### S1. Per-sweep scalar series

- **What:** one record per measured sweep holding E, m, Σ|c|², Σ|c|² without the largest cluster, and the largest cluster's size.
- **What it unlocks:**
  - autocorrelation times, so χ and C_v from a single run get a Monte Carlo error bar (today only the spread across replicates gives one);
  - Ferrenberg–Swendsen reweighting, single- and multi-histogram, which needs per-sweep energies paired with the observables;
  - equilibration checks on the measured stretch.
- **Cost:** five doubles per sweep, about 0.4 MB per run at 10⁴ sweeps. Over a sweep that is the temperature count × the replica count × that.
- **Reuse:** `Maths.Samplers` already has rank-normalized ESS from Geyer's autocorrelation sum (`src/maths/samplers/ensemble/ConvergenceDiagnostics.cs`, over a `ChainTrace`). SW is not on the ensemble shell, so the series would feed those estimators directly.

### S2. Higher moments

- m⁴, for the Binder cumulant U = 1 − ⟨m⁴⟩/(3⟨m²⟩²), and E³ and E⁴, for the error on C_v and the energy distribution's shape.
- With S1 these are views over the series. If S1 ends up optional, running sums of them are worth keeping regardless: they cost three doubles per run.

### S3. The largest cluster per sweep

- The cycle finds max |c| every sweep, and today it survives only folded into Σ|c|² without the largest cluster (seam 7 of the retention inventory).
- Its mean over N is the percolation order parameter. It belongs in S1's record.

### S4. The burn-in trace

- The same per-sweep record, kept for the burn-in cycles, so the burn-in length can be checked, not assumed.
- An option, since the trace's purpose is diagnostic.

### S5. Per-sweep partitions, an option with thinning

- **What:** the cluster labels of each kept sweep, N integers per sweep.
- **What it unlocks:** partition posteriors, draw-level lineage across T, soft memberships (the soft-clustering unification), and a process trace of the clusters (the process-visualization ambition). Today only co-membership counts survive.
- **Cost:** about 40 MB per run at N = 1,000 and 10⁴ sweeps, so gigabytes over a sweep's temperatures and replicas. It is an option with thinning, never a default.
- **Persistence:** at that size it belongs behind Archivory's chunk engine, which is not built. An in-memory option with heavy thinning could come first.

### S6. Replica spread in `SweepCurves`

- `SweepCurves.ByTemperature` (`src/clustering/graphical/spc/profiling/SweepCurves.cs`) averages each frame's reduction over the replicas at each T and drops their spread.
- The runs still carry it, so this is a view change: return the spread (or a standard error) beside each mean.

## Sequencing

1. **S1 with S3:** one per-sweep record, with the serialization change and resume appending.
2. **S2:** views over S1, or running sums if S1 is optional.
3. **S6:** independent and small; it can go any time.
4. **S4:** follows S1's record shape.
5. **S5:** waits on the persistence decision; in-memory thinning first if it is wanted sooner.

**Working rules** (from the regression track):
- Retention must not change the RNG stream. Check each change against the old code with a probe that diffs the accumulator's outputs bit for bit.
- One concern per commit.

## Checks worth having

Validation ground truth from outside the estimator's own model:
- **Reweighting:** Ferrenberg–Swendsen estimates at a neighbouring temperature against a direct simulation there.
- **Binder crossing:** on the 2D q-state Potts lattice, the crossing against the exact T_c = 1/ln(1 + √q) (J = 1, k_B = 1).
- **Autocorrelation:** τ_int from the series against batch means on the same run.

## Open questions

| ID | Question | Notes |
|---|---|---|
| SQ1 | Is S1 on by default? | Its size argues yes (0.4 MB per run); the binary format change argues for a version bump either way |
| SQ2 | Where do the series live in the result? | On `Accumulator` (travels with resume), or beside it in `SwRunResult` (keeps the resume state lean) |
| SQ3 | Storage format for persisted series and partitions | The interop tier's membrane doctrine points at a standard columnar spec at the crossing |
| SQ4 | Does SW join a shared diagnostics layer? | The MCMC engine family sketch has reversible engines sharing mixing diagnostics; S1 is the precondition |

## Report

Nothing yet.
