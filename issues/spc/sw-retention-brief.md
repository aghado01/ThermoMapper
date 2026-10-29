# SW retention brief

Status: side brief, 2026-10-24. Parked for later and not scheduled. Written at the user's request, so that the Swendsen–Wang items of the BASS preparation brief's retention inventory (`issues/regression/bass-preparation-brief.md`, "Retention inventory", SPC samplers) have a home outside the regression track.

Revised 2026-10-28 after the spin-glass threads discussion (`issues/spc/opus-spin-glass-threads.md`): the `Maths.Samplers` containers S1 and S5 can reuse, R-hat across replicas, temperature integration as an S1 consumer, the state S5 keeps (SQ5), and the checks that follow. Nothing in S1–S7 had landed by then.

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

- **What:** one record per measured sweep holding E, m, Σ|c|² and the largest cluster's size c_max. Σ|c|² without the largest cluster is Σ|c|² − c_max², a view, so the record carries four values.
- **What it unlocks:**
  - autocorrelation times, so χ and C_v from a single run get a Monte Carlo error bar (today only the spread across replicates gives one);
  - R-hat across the replicas at each temperature. Each replica starts from its own random spins (seeded per replica through `SpcSeedHelper.Derive`), so the replicas are overdispersed chains. Replicas that disagree mark the transitions; SPC has no convergence readout today;
  - Ferrenberg–Swendsen reweighting, single- and multi-histogram, which needs per-sweep energies paired with the observables. It reweights the scalar curves (χ, C_v), for example to place the χ peak between grid points. Reweighting the per-edge co-membership readout needs per-sweep partitions (S5);
  - the free energy by temperature integration (Farhi, Hed, Bon, Caticha, Mak and Domany 2013): ln Z(β) − N ln q = −∫₀^β ⟨H⟩ dβ′, with the grid's own temperatures as the intermediates.
    - For Potts, Z(0) = qᴺ for every coupling graph on the same nodes, and H is bounded, so the integral gives the absolute free energy and needs no cutoff. The cycle records H = −Σ J δ once per undirected edge, the Hamiltonian it samples.
    - A trapezoid over the running ⟨E⟩ sums works today. Stepping-stone, and path sampling with a Monte Carlo error, need the per-sweep energies.
    - The integral starts at β = 0, which the grid does not reach. The closure is a reference rung, as BASS has: iid uniform spins are exact β = 0 draws, and stepping-stone from them up through a few hot rungs closes the integral (threads note, Thread 5);
  - equilibration checks on the measured stretch.
- **Cost:** four doubles per sweep, about 0.3 MB per run at 10⁴ sweeps. Over a sweep that is the temperature count × the replica count × that, about 60 MB at 48 temperatures and 4 replicas.
- **Reuse:**
  - `ChainTrace` (`src/maths/samplers/ensemble/ChainTrace.cs`) holds a set of functionals per draw, which is S1's record. `ChainDiagnostics.Diagnose` (`ConvergenceDiagnostics.cs`) gives rank-normalized split R-hat (bulk and tail) and ESS from Geyer's autocorrelation sum over a set of traces, with the replicas at one temperature as the chains.
  - The energy column feeds `RungDraws.FromChains` and `LadderEvidence.PathSampling` / `SteppingStone` (`src/maths/samplers/evidence/LadderEvidence.cs`).
  - SW is not on the ensemble shell and need not be: it reaches the diagnostics through these data shapes (SQ4).

### S2. Higher moments

- m⁴, for the Binder cumulant U = 1 − ⟨m⁴⟩/(3⟨m²⟩²), and E³ and E⁴, for the error on C_v and the energy distribution's shape.
- With S1 these are views over the series. If S1 ends up optional, running sums of them are worth keeping regardless: they cost three doubles per run.

### S3. The largest cluster per sweep

- The cycle finds max |c| every sweep, and today it survives only folded into Σ|c|² without the largest cluster (seam 7 of the retention inventory).
- Its mean over N is the percolation order parameter. It is S1's c_max.

### S4. The burn-in trace

- The same per-sweep record, kept for the burn-in cycles, so the burn-in length can be checked, not assumed.
- An option, since the trace's purpose is diagnostic.

### S5. Per-sweep partitions, an option with thinning

- **What:** the cluster labels of each kept sweep, N integers per sweep.
- **What it unlocks:** partition posteriors, draw-level lineage across T, soft memberships (the soft-clustering unification), and a process trace of the clusters (the process-visualization ambition). Today only co-membership counts survive.
  - Replica link overlaps: `Graphs.Observables.LinkOverlap` between two replicas at one temperature, over the kept sweeps, gives P(Q) per temperature. This works for ferromagnetic SPC, with no frustration involved: a bimodal P(Q) near a transition says the replicas settle into different partitions, which the averaged co-membership hides.
  - Clustering the sampled states by a partition distance (threads note, Threads 3 and 4).
- **Which state (SQ5):** SW's Markov state is the spins; the FK clusters are auxiliary. The spins do not determine the clusters (the bonds are random), while the clusters with each cluster's new colour determine the next spins (the recolouring pass of `SwendsenWang.Draw`). They keep different information, so the options are the spins, the FK labels, or both.
  - Spins take one byte per node while q ≤ 255, a quarter of int labels.
  - The link overlap reads either. On spins it is the Potts link overlap; on FK labels it compares co-memberships, which is what SPC reads.
- **Cost:** about 40 MB per run at N = 1,000 and 10⁴ sweeps as int labels (10 MB as byte spins), so gigabytes over a sweep's temperatures and replicas. It is an option with thinning, never a default.
- **Container:** `DrawStore<TDraw>` (`src/maths/samplers/ensemble/DrawStore.cs`) keeps every `every`-th draw, and `AppendFrom` covers resume.
- **Persistence:** at that size it belongs behind Archivory's chunk engine, which is not built. An in-memory option with heavy thinning could come first.

### S6. Replica spread in `SweepCurves`

- `SweepCurves.ByTemperature` (`src/clustering/graphical/spc/profiling/SweepCurves.cs`) averages each frame's reduction over the replicas at each T and drops their spread.
- The runs still carry it, so this is a view change: return the spread (or a standard error) beside each mean.

### S7. An exact cluster-size histogram

- **What:** counts of clusters at each exact size 1..N, summed over the measured sweeps, beside (or in place of) the
  log-binned histogram, whose bins [2^k, 2^(k+1)) hold counts only.
- **What it unlocks:** the exact mean of every cluster functional Σ_k g(c_k), since each is linear in those counts.
  - S₂ = ⟨Σ_k c_k²⟩/N², which is already recorded through Σ|c|².
  - The size-floored count F_n = ⟨Σ_k (c_k − n)₊⟩/n, at any floor n, chosen after the run.
  - The expected partition entropy, through ⟨Σ_k c_k ln c_k⟩.
  - For any superadditive g these means are monotone in β, which makes them endpoint-only certificates for the regime-A
    stopping rule (`issues/spc-BARS/stopping-theory-primer.md`, §12.4, §12.7 and §12.13). The log-binned histogram only
    bounds them.
- **Cost:** N 64-bit counts per run, about 26 KB at N = 3,200, far less than S1. At high T nearly every node is a singleton each sweep, so a 32-bit count of size 1 overflows after about 2³¹/N sweeps (670,000 at N = 3,200); the log-binned histogram's int bins share that limit.
- **Relation to S5:** per-sweep partitions give the same means and much more, at gigabytes. S7 is the cheap route to the
  means alone.

## Sequencing

1. **S1 with S3:** one per-sweep record, with the serialization change and resume appending.
2. **S2:** views over S1, or running sums if S1 is optional.
3. **S6 and S7:** independent and small; either can go any time.
4. **S4:** follows S1's record shape.
5. **S5:** waits on the persistence decision; in-memory thinning first if it is wanted sooner.

**Working rules** (from the regression track):
- Retention must not change the RNG stream. Check each change against the old code with a probe that diffs the accumulator's outputs bit for bit.
- One concern per commit.

## Checks worth having

Validation ground truth from outside the estimator's own model:
- **Reweighting:** Ferrenberg–Swendsen estimates at a neighbouring temperature against a direct simulation there.
- **Binder crossing:** on the 2D q-state Potts lattice with q ≤ 4 (q = 2 or 3), the crossing against the exact T_c = 1/ln(1 + √q) (J = 1, k_B = 1). For q > 4 the transition is first-order and the cumulant dips instead of crossing; SPC's default q = 20 is in that regime.
- **Autocorrelation:** τ_int from the series against batch means on the same run.
- **Temperature integration:**
  - On any tree graph, ln Z is exact for every q: Z = qᴺ Πₑ (1 + (e^{βJₑ} − 1)/q). It checks the reference-rung closure.
  - For q = 2 on an L × L torus, Kaufman's exact finite-lattice Ising partition function (Potts J corresponds to Ising J/2, plus a constant −J/2 per edge in H, so ln Z shifts by βJ/2 per edge).

## Open questions

| ID | Question | Notes |
|---|---|---|
| SQ1 | Is S1 on by default? | Its size argues yes (0.3 MB per run); the binary format change argues for a version bump either way. Leaning yes (2026-10-28): keep each draw's state, and about 60 MB per sweep in memory |
| SQ2 | Where do the series live in the result? | On `Accumulator` (travels with resume), or beside it in `SwRunResult` (keeps the resume state lean). Leaning `Accumulator` (2026-10-28): resume has to append the series, and the running sums are its reductions, so the series is the resume state and the sums are views; a side file would break resume |
| SQ3 | Storage format for persisted series and partitions | The interop tier's membrane doctrine points at a standard columnar spec at the crossing. The SPCX checkpoint stays internal; S5's persistence waits on Archivory |
| SQ4 | Does SW join a shared diagnostics layer? | The MCMC engine family sketch has reversible engines sharing mixing diagnostics; S1 is the precondition. Proposed (2026-10-28): through the data shapes (`ChainTrace`, `ChainDiagnostics`, `RungDraws`), not the ensemble shell |
| SQ5 | Which state does S5 keep: spins, FK labels, or both? | Spins are the Markov state and cost a byte per node; FK labels are what SPC reads; both is the full generative state of a sweep. They keep different information, so all three are candidate options |

## Report

Nothing yet.
