# SMC's adaptive mutation count: report

**Status (2026-11-08).** The option and the readout are built and gated, the measurement has run, and the rule is recorded
in `Samplers.md`, and BASS's fit takes the option as `SmcSpec.AdaptiveMutation`. The adaptive count **stays an option**;
the fixed count stays the default (reading below), confirmed by the user on 2026-11-08 ("Keep fixed") against the
rule's letter. `BayesianRegression.md` names the option (`03b6d13`).

**To rerun the measurement:** build `Maths.Regression.Tests` in Release, then run `dotnet test` filtered to
`BassSmcMutationQualificationTests`, at below-normal priority. Set `BASS_QUALIFICATION=1` and `BASS_COST_CLOCK=cycles`; the
`SMC_MUTATION_*` settings in the class's summary size the run. GBASS took about 23 minutes and KHAOS about 9, at roughly
1.7 cores under load.

The chip's assignment (2026-11-08, from the GBASS/KHAOS coordinating session): let the engine's SMC choose how many
mutation sweeps each stage takes, and measure whether that beats the fixed count at equal cost. Found in the Kimi review
read against the code (`40113e4`): `SmcOptions.MutationSteps` was fixed for a whole run, and where one reversible-jump
sweep moves little a stage leaves its particles near where resampling put them. The known cases are KHAOS's Legendre fit on
Ishigami (about 2% of births and almost no swaps accepted) and GBASS's skewed heavy-tailed noise (a ridge between the
skew and the latents).

## The measurement's decision, fixed before it ran

**Arms.** The fixed count at 2, 5 and 10 sweeps per stage, the cost frontier. The adaptive count with a budget of 5 sweeps
per stage, at c = 0.1 and c = 0.01.

**Fixtures.**
- GBASS slice 5's evidence fixtures: t(5) and the skewed generalized hyperbolic, n 12. Each has a plain Monte Carlo ln p(y)
  over 2,000,000 prior draws.
- KHAOS's Ishigami fit under the raw ridge: degree 8, order 3, n 120, σ² ~ IG(3, 0.5) so the normalizer is an evidence.
  It has no exact ln p(y), so its replicates' spread stands for the error.
- Posterior means of M, σ² and the skew, each against eight long MCMC chains.

**Cost.** Seconds on the thread's own cycle counter (indicative on the shared machine), with sweeps per particle beside them.
Equal cost is read through the work-normalized error, MSE × seconds, which a change of N leaves alone.

**The rule.** Adaptive becomes `SmcSpec`'s default, at the better c, when both of these hold:
- on at least one fixture, its ln Ẑ work ratio to the fixed 5 has a 95% interval wholly below 1;
- on no fixture is any of its work ratios (ln Ẑ, M, σ², skew) wholly above 1 on the same interval.

Otherwise it stays an option, with this reading recorded.

**Amended after the pilot, before the measurement** (N 100, a few replicates, references cut short):
- *What the pilot showed.*
  - On both GBASS fixtures every stage wanted at least 5 sweeps, so the adaptive runs were the fixed runs bit for bit.
  - On KHAOS most replicates of every arm collapsed to near-empty structures (M ≈ 0.2, σ² ≈ 12.5, ln Ẑ ≈ −349.7). One adaptive
    replicate grew the posterior's structure (M 9.2, σ² 0.34) and read ln Ẑ 43 nats higher.
- *Why the rule changed.* A replicates' SD mixes trapped and found runs, so it is no error measure there.
- *The amendments.*
  - KHAOS's ln p(y) reference is the pooled ln of the mean Ẑ over every replicate of every arm. Ẑ is unbiased in Z, so that
    pool leans on the replicates that found the posterior; it is a floor if none did.
  - Each arm reports how many replicates reached the posterior: weighted σ² within twice the chains'.
  - The interval on a work ratio comes from a bootstrap over replicates, paired, because every arm's replicate r starts from
    one seed.
  - A budget of 20 joins the arms: fixed 20, and adaptive 20 at c = 0.01.
  - Each adaptive arm is read against the fixed arm at its own budget.
- *The decision rule itself* is unchanged.

## What landed

| Commit | What |
|---|---|
| `8ec1bc8` | `SmcOptions.AdaptiveMutation` (`src/maths/samplers/smc/SmcMutation.cs`), `SmcStageReport.Steps` and `MovedFraction`, `SmcStageReadout`; core tests |
| `9261c9f` | `BassSmcReplicate.StageReadout`, the per-replicate readout |
| `ebe5757` | The oracle law test under the adaptive count |
| `7d0869a` | The measurement harness (`BassSmcMutationQualificationTests`), its tables (`smc-adaptive-mutation-results.csv`, `-replicates.csv`, `-stages.csv`) and this report |
| `11436c4` | `Samplers.md`: the rule, its alternatives and trade-offs, the library's holdings, the measured reading |
| `7dee840` | `SmcSpec.AdaptiveMutation` (validation, serialized only when set), `BassSmcRun`'s pass-through and progress, `BassSmcTests` |

**The rule.** Stage t reads stage t − 1's report: k sweeps left a fraction 1 − m of the live particles where they were (a
particle moved when any of its moves was accepted). So stage t wants ⌈k·ln c / ln(1 − m)⌉, at least 1. It takes that
within its credit, k̄·t − Σ_{s<t} k_s with k̄ = `MutationSteps`. Easy stages pass unspent sweeps to hard ones, and a run
never sweeps more than the fixed count over the same stages. Why the moved fraction, why a credit, the invariance
argument and the alternatives not taken are in `Samplers.md`, "SMC's sweeps per stage".

**The readout.** `BassSmcReplicate.StageReadout` lays out, stage by stage: b, sweeps, the structure move's acceptance,
the fraction of particles whose structure changed, ESS/N, resampling and families. `LeastMoved` names the worst stage.
On KHAOS it reads 0.02–0.06 at the last stages, which is the trap below.

**Gates.**
- *The fixed count's draws are unchanged.* The regression assembly's commit scope passed, pinned tier included: 1,755 of
  1,755 with the core, and 1,541 of 1,541 with the readout.
- *The laws under the adaptive count.*
  - Core (`SmcTests`): held at b = 0 the population keeps the prior; tempered to b = 1 it meets the conjugate posterior and
    evidence, and stays there.
  - BASS (`BassSmcParityTests`, against the R enumeration oracle): held at b = 0 it keeps the prior's count law, E[τ] and
    inclusions. At b = 1 it meets the exact posterior, and ln Ẑ pools to −28.600 ± 0.029 against −28.589. In both, the
    count moved off the budget.
- *Resume bit for bit.*
  - Core: resume, and independence from the parallelism.
  - BASS: an adaptive fit whose count moves resumes bit for bit in 3-stage pieces. The wiring's commit scope passed
    1,597 of 1,597, pinned tier included.

## The measurement

Sweeps per particle and seconds are means per replicate. Seconds are on the thread's cycle counter, indicative only, with
the machine at about 65% load from other sessions. The work ratio is MSE × s against the fixed arm at the same budget,
with a paired bootstrap's 95% interval. Full columns are in `smc-adaptive-mutation-results.csv`.

| Fixture | Arm | Sweeps | s | Reached | ln Ẑ mean ± sd | RMSE ln Ẑ | RMSE M | RMSE σ² | ln Ẑ work ratio |
|---|---|---|---|---|---|---|---|---|---|
| t(5), ref −17.257 ± 0.010 | fixed 2 | 32 | 1.09 | 32/32 | −17.239 ± 0.030 | 0.034 | 0.051 | 0.011 | |
| | fixed 5 | 80 | 2.14 | 32/32 | −17.242 ± 0.024 | 0.028 | 0.046 | 0.012 | |
| | fixed 10 | 160 | 3.94 | 32/32 | −17.238 ± 0.021 | 0.028 | 0.037 | 0.010 | |
| | fixed 20 | 320 | 7.42 | 32/32 | −17.244 ± 0.020 | 0.023 | 0.036 | 0.009 | |
| | adaptive 5, c 0.1 and 0.01 | 80 | 2.08 | 32/32 | identical to fixed 5 | | | | clock noise only |
| | adaptive 20, c 0.01 | 307 | 7.16 | 32/32 | −17.242 ± 0.016 | 0.022 | 0.037 | 0.009 | 0.86 [0.44, 1.54] |
| skewed GH, ref −17.980 ± 0.030 | fixed 2 | 45 | 1.70 | 32/32 | −18.004 ± 0.027 | 0.036 | 0.046 | 0.013 | |
| | fixed 5 | 114 | 3.67 | 32/32 | −18.013 ± 0.030 | 0.044 | 0.052 | 0.011 | |
| | fixed 10 | 227 | 6.68 | 32/32 | −18.011 ± 0.024 | 0.039 | 0.044 | 0.012 | |
| | fixed 20 | 456 | 12.98 | 32/32 | −18.013 ± 0.025 | 0.042 | 0.045 | 0.011 | |
| | adaptive 5, c 0.1 and 0.01 | 114 | 3.70 | 32/32 | identical to fixed 5 | | | | clock noise only |
| | adaptive 20, c 0.01 | 448 | 12.70 | 32/32 | −18.006 ± 0.025 | 0.036 | 0.039 | 0.007 | 0.74 [0.44, 1.25] |
| KHAOS, floor −287.78 | fixed 2 | 122 | 4.11 | 0/16 | −347.3 ± 5.5 | 59.8 | 9.1 | 10.8 | |
| | fixed 5 | 264 | 7.93 | 1/16 | −347.9 ± 7.6 | 60.6 | 9.6 | 11.6 | |
| | fixed 10 | 516 | 15.27 | 1/16 | −347.9 ± 7.4 | 60.6 | 9.5 | 11.3 | |
| | fixed 20 | 1,129 | 28.94 | 2/16 | −342.0 ± 21.5 | 58.1 | 9.0 | 10.6 | |
| | adaptive 5, c 0.1 | 239 | 5.46 | 0/16 | −349.83 ± 0.09 | 62.1 | 9.9 | 12.1 | 0.72 [0.69, 0.80] |
| | adaptive 5, c 0.01 | 243 | 6.59 | 0/16 | −349.80 ± 0.20 | 62.0 | 9.9 | 11.9 | 0.87 [0.83, 0.97] |
| | adaptive 20, c 0.01 | 964 | 12.80 | 0/16 | −349.81 ± 0.12 | 62.0 | 9.6 | 11.3 | 0.50 [0.44, 0.64] |

N 500 throughout; 32 replicates on GBASS, 16 on KHAOS. The KHAOS chains read M 10.146 ± 0.007 and σ² 0.2245 ± 0.0002.

**The reading.**
- *GBASS: the count never moves at a budget of 5.* The structure move accepts 30–55% at every b, and some particles sit
  stuck at the structure's cap. So the moved fraction at 5 sweeps never reaches the 0.944 at which a stage would want
  fewer, and the adaptive run is the fixed run bit for bit.
  - At a budget of 20, the count saves 2–4% of the sweeps. Its ln Ẑ work ratios, 0.86 and 0.74, have intervals well
    across 1.
  - The fixed frontier favours few sweeps here: by MSE × s, fixed 2 is the cheapest per unit error on ln Ẑ.
  - The skewed GH's references (±0.030) limit its RMSE column. Read its SDs there.
  - The structure's moved fraction cannot see the skew–latent ridge. τ's slice step and the latents always count as
    moved, and the ridge lies in them.
- *KHAOS: every arm is trapped.*
  - The schedule spends about 40 stages at b < 0.001, then the moved fraction collapses over the last eight, from 0.9 to
    0.02.
  - Most replicates end at near-empty structures: M ≈ 0.1–0.8, σ² ≈ 12, ln Ẑ ≈ −349.8.
  - The fixed arms grew the structure in 6 of 64 replicates: 4 reached M ≈ 10, and fixed 2's two reached M ≈ 7. Those took
    108–124 stages and read ln Ẑ from −334 to −283. Even the found runs disagree by 37 nats.
  - The adaptive arms grew it in none of 48.
  - The allocation shows why the count fails here. It saves early: 3 sweeps while 94% of particles move, at b ≈ 0. It
    then spends its credit on the moderately hard stages just before the collapse (6 → 20 sweeps at a budget of 5; up to
    143 at a budget of 20). Reading only the previous stage, it reaches the last, hardest stages with nothing left above
    the budget.
  - Fewer early sweeps also leave less diversity for selection to find.

**The decision: adaptive stays an option; the fixed count stays `SmcSpec`'s default.**
- On GBASS the rule's first condition fails: at a budget of 5 the arms are identical, and at 20 the intervals cross 1.
- On KHAOS the letter of the rule holds: the ln Ẑ work ratios are 0.72 [0.69, 0.80] and 0.50 [0.44, 0.64]. Its premise
  does not. MSE × s compares at equal cost only where the error is variance, which falls as N grows. On KHAOS the error is
  a trap's bias, about 60 nats in every arm.
  - The adaptive arms read lower only because they cost less: a trapped run ends in about 48 stages, an escaping one in
    about 115.
  - They are no less wrong: 0 of 48 against 6 of 64.
- I did not apply the rule's letter, and say so here: this departs from the recorded rule. Should the user read it the
  other way, the change is one default in `SmcSpec`, and with it the re-baselining of the pinned SMC hash.
- **Confirmed by the user, 2026-11-08: "Keep fixed."** The fixed count stays `SmcSpec`'s default. A later rule of this
  kind compares work ratios only among runs that reach the posterior, or counts the runs that do beside them.

## Open, each with the decision it would change

- **KHAOS's SMC is not to be trusted at these budgets.** Up to 20 sweeps per stage at N 500, at most 2 of 16 replicates
  reach the posterior. The rest report evidence about 60 nats low, with genealogy errors under 1.
  - *Decision:* whether KHAOS's SMC fit and its evidence are refused, flagged (say on `LeastMoved` near 0), or given a
    start or kernel that can find the structure: a reference nearer the posterior, or a birth that proposes from the
    data's residual. KHAOS slice 2 should know before its SMC path is relied on.
  - **Decided 2026-11-08: flagged.** The coordinating session's call, the user unsure; reversible. KHAOS slice 5 measured
    the mechanism: a near-empty population reweights alike, so the schedule jumps to b = 1 and ln Ẑ is read there.
    `LeastMoved` did not separate it: it read under 0.1 on reached replicates too. BASS's SMC fit now flags a jump over a
    collapsed population and replicates whose ln Ẑ disagree beyond their errors. The remedies are recorded, not built
    ([the KHAOS brief](khaos-brief.md#reports)).
- **A count that sees continuous coordinates.** The moved fraction cannot see GBASS's skew–latent ridge. The
  decorrelation of ℓ across a stage's sweeps can (the `Samplers.md` section's alternative).
  - *Decision:* whether an adaptive count can help GBASS's skewed horseshoe, where the ridge fails Geweke at b = ½.
  - Measure it there; this measurement's skewed GH does not show the ridge.
- **Front-loaded allocation.** On search-limited targets the early, cheap sweeps are where structure is found, so a
  credit that saves them works against the run.
  - *Decision:* whether KHAOS's SMC can be made to work by allocation alone.
  - Worth measuring only once the trap above has a remedy to compare against.
- **`MutationSteps` 5 against 2.** On the small GBASS fixtures, 2 sweeps per stage is cheaper per unit ln Ẑ error than 5.
  - *Decision:* `SmcSpec.MutationSteps`'s default.
  - That needs the stage-8 Friedman fixtures, where R3 of the Friedman band brief already compares 20 against 5. This
    reading alone does not move it.
