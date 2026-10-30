# SW retention brief

Status: side brief, 2026-10-24. Parked for later and not scheduled. Written at the user's request, so that the Swendsen–Wang items of the BASS preparation brief's retention inventory (`issues/regression/bass-preparation-brief.md`, "Retention inventory", SPC samplers) have a home outside the regression track.

Revised 2026-10-28 after the spin-glass threads discussion (`issues/spc/opus-spin-glass-threads.md`): the `Maths.Samplers` containers S1 and S5 can reuse, R-hat across replicas, temperature integration as an S1 consumer, the state S5 keeps (SQ5), and the checks that follow. Nothing in S1–S7 had landed by then.

Corrected 2026-10-29: SPC has no mid-run resume (the SPCX file caches finished tasks), so the series has nothing to append across; SQ2's leaning keeps its conclusion on the cache's grounds. The log-binned histogram counts in 64 bits since `1c42741` (SPCX v8).

Chipped 2026-10-29: the items that stay inside SPC (S1–S4, S6, S7, and the Binder check) go to [the chip brief](sw-retention-chip-brief.md), under the user's frame that accumulation is the run's state and the sums its views. S5 and the consumers in `Maths.Samplers` (R-hat, ESS, the verdict, temperature integration, reweighting) stay here.

Built 2026-10-29: S1–S4, S6 and S7 landed through the chip (`6de243a`–`d7269ea`, review `8c4c1ae`; [its report](sw-retention-chip-brief.md#reports)). What remains here is S5 and the samplers-side consumers. S2's energy moments, the energy cumulant and C_v's Monte Carlo error followed in `c50b273`.

Built 2026-10-30: S5, as kept partitions and spins, every k-th draw (`a777534`, SPCX v12). Thinned to about one draw per autocorrelation time, it needs no Archivory. What remains here is the samplers-side consumers, and the partitions' own consumers (VI-loss cuts, credible balls, the robust window, replica P(Q)), which sit with the Bayesian clustering session.

The principle, from the user: keep each draw's generative state, derive views from it, and where alternatives keep different information, offer both as options.

## Where things stand

- **The cycle** (`src/clustering/graphical/spc/runtime/core/sampler/SwendsenWang.cs`). Every sweep computes the spins, the FK bonds and their clusters (union-find), the energy, the magnetization from the largest colour, Σ|c|², and the largest cluster's size.
- **The accumulator** (`Accumulator.cs`). It keeps:
  - running sums of E, E², m, m², Σ|c|² and Σ|c|² without the largest cluster;
  - the log-binned cluster-size histogram;
  - the final spins and RNG words, which the engine's `Restore` would resume from (nothing calls it);
  - the optional per-edge counts (`AccumulationSpec`: affinities, alignments, co-membership) and per-node sums (the cluster-size and giant-participation landscapes);
  - since the chip: the per-draw record, always kept (`EnergyPerDraw`, `MagPerDraw`, `SumSqClusterSizesPerDraw`, `MaxClusterSizePerDraw`), and two options, the burn-in trace (`BurnInTrace`) and the exact cluster-size histogram (`ExactClusterSizeHistogram`).
- **The runner** (`SwRunner.cs`). It discards `RunBudget.BurnIn` cycles, then accumulates `Cycles`. A sweep keeps every run's accumulator, per temperature and replica (`SweepResult.SweepRuns`).
- **The serialization.** The accumulator is serialized as SPCX v11 (`SwRunSerialization.cs`, `AccumulatorSerializer`; the record and both options since v9–v11) as a cache of finished tasks: `SpcExecutor` loads a completed task's file (`SpcxDiskFrameSink.TryLoad`) instead of rerunning it, and the loaded result carries the accumulator alone. There is no mid-run resume. A retained series kept on the accumulator becomes part of the binary format and returns with a cached task; one kept anywhere else is lost on that path.

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
- **Cost:** as built, 28 bytes per sweep (E and m as doubles, Σ|c|² as a long, c_max as an int), about 0.28 MB per run at 10⁴ sweeps. Over a sweep that is the temperature count × the replica count × that, about 54 MB at 48 temperatures and 4 replicas.
- **Reuse:**
  - `ChainTrace` (`src/maths/samplers/ensemble/ChainTrace.cs`) holds a set of functionals per draw, which is S1's record. `ChainDiagnostics.Diagnose` (`ConvergenceDiagnostics.cs`) gives rank-normalized split R-hat (bulk and tail) and ESS from Geyer's autocorrelation sum over a set of traces, with the replicas at one temperature as the chains.
  - The energy column feeds `RungDraws.FromChains` and `LadderEvidence.PathSampling` / `SteppingStone` (`src/maths/samplers/evidence/LadderEvidence.cs`).
  - SW is not on the ensemble shell and need not be: it reaches the diagnostics through these data shapes (SQ4).

### S2. Higher moments

- m⁴, for the Binder cumulant U = 1 − ⟨m⁴⟩/(3⟨m²⟩²), and E³ and E⁴, for the error on C_v and the energy distribution's shape.
- With S1 these are views over the series, and S1 is always kept. The Binder cumulant is built (`BinderCumulant`, `f542951`). So are ⟨E⟩ to ⟨E⁴⟩ and the energy cumulant V = 1 − ⟨E⁴⟩/(3⟨E²⟩²) (`EnergyCumulant`, `c50b273`), the first-order diagnostic for SPC's default q = 20.
- C_v's Monte Carlo error is built without E³ and E⁴: C_v is the mean of y = (E − Ē)²/T², so its error is sd(y)/√ESS(y), with the ESS from `ChainDiagnostics.MeanEss` over the replicas (`c50b273`). SPC's project already referenced `Maths.Samplers`; the call is on plain arrays, the user's choice over a port (2026-10-30). It matches Gaussian theory (iid 0.98, AR(1) 1.03) and the spread of C_v across SW replicas (1.04–1.05 at 128 replicas).
- Generalized the same day to `MonteCarloError` (`OfMean`, `OfVariance`), for the two curves the stopping rule reads (`issues/spc-BARS/stopping-theory-primer.md`, §12.13 and §12.15): the detection curve Var(m), and χ_FK, whose error over N is S₂'s within-run SE. Against the spread across 32 SW replicas: Var(m) 1.15, χ_FK 1.05–1.07; against theory, the mean form 1.00 (iid) and 1.03 (AR(1)).

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
- **As built (2026-10-30, `a777534`).** The gigabyte figures above assume every sweep kept. The consumers need about as many partitions as there are effectively independent draws, which the per-draw record's ESS now measures; at a few hundred per run (0.8 MB at N = 1,000 and 200 draws) they fit in memory and the SPCX cache, as S7 did.
  - `AccumulationSpec.PartitionEvery` keeps the FK partition of every k-th draw (draw t when t mod k = 0, as `DrawStore` thins), canonically labelled by each cluster's smallest node index. `SpinsEvery` keeps the leaving spins as bytes (q ≤ 256). Both are `Accumulator` arrays, flattened draw-major (`KeptPartitions`, `KeptSpins`), not `DrawStore`, so SPC holds no `Maths.Samplers` type, as with S1.
  - Checked against what the run already accumulates: from k = 1 partitions, Σ|c|², c_max, the exact histogram and the co-membership counts rebuild exactly; kept spins colour each cluster uniformly and give the record's m bit for bit; thinning keeps exactly those draws of the same stream; the accumulator probe stays bit-identical.
  - Not built: thinning by the run's own measured τ (keep all, thin at the end), compression of the canonical labels (`System.IO.Compression`, base library), and Archivory for keeping every draw at scale.

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
- **Cost:** N 64-bit counts per run, about 26 KB at N = 3,200, far less than S1. At high T nearly every node is a singleton each sweep, so the count of size 1 gains about N a sweep, and a 32-bit count would wrap after about 2³¹/N sweeps (670,000 at N = 3,200). The counts are 64-bit, as the log-binned histogram's have been since `1c42741`.
- **Relation to S5:** per-sweep partitions give the same means and much more, at gigabytes. S7 is the cheap route to the
  means alone.

## Sequencing

1. **S1 with S3:** one per-sweep record, with the serialization change. Appending across a resume comes only if mid-run resume is ever wired. *Done (`bc70bf2`).*
2. **S2:** views over S1. *Done: Binder (`f542951`); energy moments, energy cumulant and C_v's error (`c50b273`).*
3. **S6 and S7:** independent and small. *Done (`96de69d`, `d7269ea`).*
4. **S4:** follows S1's record shape. *Done (`1322e47`).*
5. **S5:** waits on the persistence decision; in-memory thinning first if it is wanted sooner. *Done, thinned, in the SPCX cache (`a777534`).*

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
| SQ1 | Is S1 on by default? | Its size argues yes (0.3 MB per run); the binary format change argues for a version bump either way. Leaning yes (2026-10-28): keep each draw's state, and about 60 MB per sweep in memory. **Decided 2026-10-29 (user):** always collected, as the sums and the histogram are, with no switch and no thinning |
| SQ2 | Where do the series live in the result? | On `Accumulator` (travels with resume), or beside it in `SwRunResult` (keeps the resume state lean). Leaning `Accumulator` (2026-10-28): the running sums are the series' reductions, so the series is the state and the sums are views. The resume argument given then does not hold (corrected 2026-10-29: no mid-run resume exists). The cache does: a task loaded from its SPCX file carries the accumulator alone, so a series kept anywhere else is lost for cached tasks. **Decided 2026-10-29 (user):** on `Accumulator`; accumulation is the run's state, and the sums are its views |
| SQ3 | Storage format for persisted series and partitions | The interop tier's membrane doctrine points at a standard columnar spec at the crossing. The SPCX checkpoint stays internal; S5's persistence waits on Archivory |
| SQ4 | Does SW join a shared diagnostics layer? | The MCMC engine family sketch has reversible engines sharing mixing diagnostics; S1 is the precondition. Proposed (2026-10-28): through the data shapes (`ChainTrace`, `ChainDiagnostics`, `RungDraws`), not the ensemble shell |
| SQ5 | Which state does S5 keep: spins, FK labels, or both? | Spins are the Markov state and cost a byte per node; FK labels are what SPC reads; both is the full generative state of a sweep. They keep different information, so all three are candidate options. **Decided 2026-10-30 (user):** both, as independent options (`PartitionEvery`, `SpinsEvery`); FK labels are what the Bayesian consumers read |

## Report

The SPC-side items' report, with the tests, the Binder crossing and the reversible defaults left for the user, is in [the chip brief](sw-retention-chip-brief.md#reports).
