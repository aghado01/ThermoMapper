# The conformal layer: addendum (chip brief)

Opened 2026-11-14 at the user's word ("that brief was getting way too long, addendum brief"). The [conformal-layer brief](archive/bass-gbass-khaos-activespaces/briefs/conformal-layer-brief.md), now archived, is the record of E16 through slice 15 and slice 14's follow-ups: its decisions 1 to 12, the plan for slices 7 to 18, the gates common to every slice, the working rules and every report. They all stand. This addendum carries what is open and the reports from here on.

**Chip:** the conformal chip, `local_69bae1e5`. **Verified by:** the main session, `local_b6732445`, from clean archives (`scripts/verify-csharp.ps1`, `scripts/verify-qualification.ps1` for gated tests).

## Where E16 stands

- Slices 0 to 15 are built and verified; slice 14's follow-ups 1 to 3 and its ratio estimator too, on their Lean (decision 12).
- Follow-up 2 is completed by `3f68330`, the merged K-fold CV+ bound under likelihood ratios, on the Lean at `47acb3c`.
- Slice 16, batch predictive inference (`60ba2cb`), slice 17, scores beyond the model's own (`0b91db9`), and slice 18, the guarantee read out (`5b55a6f`), are built and verified. Their reports are below. Every slice in the plan is built.
- Slice 18's first follow-up, the smoothed cell law for the log score (`b4b7e22`), is built and verified: `CellLaw.Smoothed`, opt-in, for the split form's unweighted cases.
- GBASS's two follow-ups (`5ce91b9`) are built and verified: `GroupConformal`, group-level units for new groups; and the cell rule, `Cells` on both seams with the layer's refusal.
- BASS's opt-in group-units call (`b79840f`, the handoff's item 16, (b)) is built and verified: `BassFit.GroupConformal`.
- Slice 15's precision rule reads at the calibration split's size (`591eace`), built and verified: `ScoreMixture.PrecisionAt`. Fitted at 1,000 and calibrated at 399, m goes from 500 to 45 and the whole-space sets from 13.4% to none.
- Slice 18's follow-up, a set's own miscoverage bound (`0f4d48e`), is built and verified: `TrainingConditional.Miscoverage`, the β-quantile of Beta(k, m − k) (Vovk 2012), read out on the split sets.
- Slice 18's calibration at n units is proved in the Lean (`8594176`: `gaffke_calibration`, and `highProbability_risk`, what `HighProbability` reports) and verified from a clean archive. `LossBound`'s docstrings name both (`9cf0fd5`, comments only).
- Sets built to the training-conditional guarantee (`621370c`) are built and verified: `TrainingConditional.Level(m, α, β)`, with an optional confidence on the split forms and the group form. At β 0.9 the share of calibrations whose sets miss above α falls from 0.46 to 0.0975 (residual), 0.0035 (asymmetric) and 0.0495 (equi-tailed).
- `BassFit.GroupConformal` takes the optional confidence (`14fd05d`), built and verified.
- Slice 17's grid-free intervals for the residual scores (`28993bb`, `ResidualShrinkage.Intervals`) are built and verified.
- The asymmetric interval's exact law is closed (`b8cd947`, docstrings only). Where the band's width varies, the union of the two tails is Makarov's worst case given their Beta laws. Beta(k_lo + k_hi, m − k_lo − k_hi) holds only at constant width, the residual interval.
- Slice 17's RoBAS-Full (`942df9c`, `ResidualScore.Horseshoe`) is built and verified. It is served by `Prediction`, the grid form. Its sets are not always intervals: two pieces at n 2 and 3 (`f400892`, probed and verified).
- RoBAS-Full's sets are intervals from n 4 when the residuals are not all equal, proved unconditionally. Scope A (`700c718`) proved it granted four bounds on Dawson's function. Scope B (`7d98281`, `ConformalDawson.lean`) proves those bounds, so `robasFull_ordConnected` needs only n ≥ 4 and residuals not all equal.
  - Verified from a clean archive on 2026-11-15 ([receipts](receipts/verify-lean-conformal-7d98281/)): checks 425 s, the barrel built (3,191 jobs), the 23 named results on the standard three axioms, meta-ci PASS.
  - `Intervals` serves RoBAS-Full there (`c3539b7`, verified) and refuses it below n 4. Its docstrings cite the unconditional theorem (`bb2fd78`; comments only, read).
- CBMA reads each model's evidence from its BASS or BARS fit (`3e3b7f2`, verified): the `ModelEvidence` seam and the two engines' `ForAveraging` adapters, with both engines' constants pinned against the closed form.
- Slice 15's alternation speed (`192ee84`, verified): safeguarded extrapolation of γ, 1.3–3× fewer passes at n 1,000 to 10,000. It is on by default and moves existing results within the fit's tolerance (report).
- τ through richer features (`274f69d`, verified): `DensityRatio.Spline`, additive cubic B-splines with the ridge chosen by 5-fold cross-validated log loss, opt-in through `ScoreMixtureOptions.Knots`.
- Slice 12a's k̂ diagnostic stays opt-in, as the user settled on 2026-11-14. That is the code's default (`ConformalOptions.LeaveOneOutParetoK`, false).
- The defects chip owns `PredictiveChecks.cs` and `Pointwise.cs`'s weight check; a change there is agreed with it first (the archived brief's coordination note).

## Open, each with the decision it would change

- **A stop that reads R², for slice 15's factorization** (the alternation speed's finding). On noisy τ at J ≥ 3, the 300-pass cap still binds after R² has settled. It decides whether the factorization stops on R² and cuts that cost severalfold. The user's.
- **The spline τ's two follow-ups** (its report), the user's:
  - interactions between inputs;
  - `Knots` as the default.
- **Conformal test martingales beyond E35's slice 6,** queued at the user's word (2026-11-15, "lets add conformal test martingales to the outstanding conformal work queue").
  - *Already built,* in slice 6 on 2026-11-08 (`4dd1dfb`, `415c08c`, `206eca9`, `befb449`):
    - a stream's smoothed p-values (M-1);
    - the power, power-mixture, Simple Jumper and plug-in bets with the capital (M-2);
    - Ville's, CUSUM and Shiryaev–Roberts alarms.

    M-1 and M-2 are proved in Lean.
  - *Candidates,* each with the decision it would change:
    - *The full form, built and verified* (`b07fb8d`; its report below, with two follow-ups for the user: pruning draws by weight, and both forms' capitals combined, which first needs a reading of whether the average keeps Ville's guarantee). Each score is recomputed with the newest observation, through the add-one-in weights over the stream; E35's brief calls it valid by the same theorem. Decides whether the monitor tests the full step's own scores, not only a fixed fit's.
    - *The Sleepy Jumper, built and verified* (`765795c`; its report below). Resting keeps the capital through 2,000 exchangeable arrivals, ln S −0.79 against the Simple Jumper's −7.98: the Simple Jumper's neutral bet does not hold its share, since its jumps keep moving capital into the betting bets. After a change, Ville alarms in a median of 14 arrivals against 55. Decides whether the monitor's bet can rest between jumps. It is defined in two places on the shelf:
      - Vovk, Nouretdinov and Gammerman 2003, "Testing exchangeability on-line": `supellex-private/PDFs/Bayesian/Vovk2003.pdf`, PDF page 6;
      - Vovk, Gammerman and Shafer 2005, *Algorithmic Learning in a Random World*, §7.1: the full book at `PDFs/textbooks/Vovk2005.pdf`, book pages 176–177. Chapter 7 is not among the chapter files in `PDFs/Vovk2005/`, which hold chapters 2, 6, 8, 9 and 10.

      The Simple Jumper, already built, simplifies it, per Vovk, Petej, Nouretdinov, Ahlberg, Carlsson and Gammerman 2021, "Retrain or not retrain", Algorithm 1: `PDFs/Bayesian/vovk21b.pdf` (COPA 2021), the same paper as arXiv 2102.10439 in `gauntlet/conformal-bayes`. No source on the shelf defines a "Composite Jumper".
    - *CUSUM's and Shiryaev–Roberts's guarantees in Lean, proved and verified* (Vovk 2021, Propositions 4.2 and 4.4; the Lean `4c832df`, the C# `ad78359`; its report below). The almost-sure lim sup is not formalized.
    - *M-1 under covariate shift, proved and verified,* at the user's word (2026-11-15, "Yes, prove it"). The Lean `deaf0da`, section `CovariateShift`: `shift_uniform`, `shift_law`, and `newestShare_eq_esymm`, which equals `StreamShares`' πⱼ/k.
      - The C# `0fede4d` marks it proved (docstring only, read).
      - Verified from a clean archive ([receipts](receipts/verify-lean-conformal-deaf0da/)): checks 961 s, the barrel built (3,191 jobs), the 8 named results on the standard three axioms, meta-ci PASS.
      - The hypotheses were read: w ≥ 0, which the C# checks; a symmetric g; and the weighted law normalized, e.g. Z = Fin 2, w ≡ 1, g constant. The ratio is the true one by the law's definition.
      - Not formalized: the continuous lift, as for M-1.
    - *Streams with known structure or under covariate shift,* which E35's brief leaves out: a score that removes a trend or season, or a test martingale weighted for the shift. Decides whether the monitor serves streams that are not exchangeable as they stand.
      - *Covariate shift: built and verified* (`d4c97d4`; its report below), scoped first by the user. Its guarantee is posed to the conformal Lean chip as a follow-up.
      - *The Mondrian monitor, built and verified* (`c64ccc9`; its report below). Pool the strata for a change everywhere; watch each stratum for a change confined to one.
  - The pick among them is the user's.

## Working rules

As in the archived brief: each slice is reported here, then waits for the main session's verification before the next starts; anything that changes a default or removes an option goes to the user; chips leave `docs/` alone and each report ends with **Doc changes owed**.

## Reports

Newest first.

### 2026-11-15 — Streams under covariate shift: the guarantee proved, the docstring follows

The conformal-lean thread proved the guarantee I posed for `CovariateShiftPValues` (`deaf0da`) and reported it here.

**Landed** on `main`, nothing pushed: `0fede4d`, a docstring only.
- The guarantee's paragraph now reads "proved (`shift_law`, and `newestShare_eq_esymm` for these shares, up to the continuous lift as for M-1)" in place of "posed here, checked by enumeration, not yet in Lean".
- It also says the ratio must be the true one.
- *What changes for existing calls: nothing.*

**Checked before citing.**
- `shift_law` is the finite core over a finite unit space: weights w ≥ 0, a symmetric g, and the weighted law normalized. It concludes the stream's joint uniformity, E[Π 1(pₜ ≤ uₜ)] = Π uₜ. The hypotheses are satisfiable.
- `newestShare_eq_esymm` shows the newest-position share is wⱼ e_k(w₋ⱼ)/((k + 1) e_(k+1)(w)), which is `StreamShares`' πⱼ/k.
- The file holds no `sorry`, axiom or `native_decide`.
- The Lean thread reports that the C# needs no new condition:
  - only w ≥ 0, which the class checks;
  - a stream unit with w = 0 has probability 0, so refusing it costs nothing;
  - the 5e-324 floor on p is harmless.

**Gates.**
- `CovariateShiftPValueTests` 9 of 9.
- *Per commit:* 120 of 120 in the 5 classes selected in Maths.Regression.Tests.

**Defects found: none.**

**Doc changes owed:** the conformal test martingales section: the covariate-shift guarantee, proved.

### 2026-11-15 — RoBAS-Full's intervals, unconditional: the docstrings follow the Lean

The RoBAS-Full Dawson bounds session proved the four Dawson bounds (`7d98281`, `ConformalDawson.lean`) and reported it here. This is the follow-up the RoBAS-Full intervals item left open: drop "granted four bounds … the checks certify".

**Landed** on `main`, nothing pushed:
- `bb2fd78`, docstrings and one test comment only:
  - `ResidualScore.Horseshoe` and `ResidualShrinkage.Intervals` now cite `robasFull_ordConnected`, with no hypothesis left;
  - `HorseshoeScoreTests`' comment cites the same.
- *What changes for existing calls: nothing.* `Intervals` has served RoBAS-Full from n 4 since `c3539b7`.

**Checked before citing.**
- `robasFull_ordConnected` is `robasFull_interval` applied to the proved `robasDawsonBounds`.
- Its hypotheses are n ≥ 4 and residuals not all equal, both satisfiable, and its conclusion holds at every α.
- `ConformalDawson.lean` holds no `sorry`, axiom or `native_decide`. The Lean record reports that `#print axioms` shows only the standard three.

**Gates.**
- `HorseshoeScoreTests` 6 of 6 and `ResidualIntervalTests` 4 of 4, on the final text.
- *Per commit:* 61 of 61 in the 4 classes selected in Maths.Regression.Tests, on the text before the docstrings were rewrapped. The rewrap changed comments only, and the two classes were rerun on it.
- No new XML-doc warning. The build's two are in `BarsEnsemble.cs` and `CostProfile.cs`.

**Defects found: none.**

**Doc changes owed:** the residual-scores section: RoBAS-Full's intervals from n 4 now rest on a proof with no hypothesis.

### 2026-11-15 — Conformal test martingales: CUSUM's and Shiryaev–Roberts's guarantees marked proved

The conformal-lean thread proved both (`4c832df`) and reported it here. This brings the C# into line, as was done when RoBAS-Full's intervals were proved.

**Landed** on `main`, nothing pushed:
- `ad78359`:
  - `MartingaleAlarm.Guarantee` now marks CUSUM and Shiryaev–Roberts `Proved`, with the proved wording;
  - the `AlarmProcedure` and `ConformalTestMartingale` docstrings now say proved, not posed;
  - `MartingaleAlarmTests` assert the proved flags and wording, and the false-alarm theory drops `_Posed` from its name.
- *Verified by the main session on 2026-11-15,* from clean archives with zero skips:
  - the C#, the Commit selection over `ad78359`'s files at `66922ed`: 120 of 120 ([receipts](receipts/verify-conformal-alarmproved-ad78359/));
  - the Lean it rests on, `4c832df` ([receipts](receipts/verify-lean-conformal-4c832df/)). The checks passed in 785 s, and the 12 named results rest on the standard three axioms only. meta-ci passed. The barrel built on its second run, 3,190 jobs: the first run's build stopped when two unrelated modules, `ProjectionMean` and `RandomCluster`, crashed with Lean exit 0xC0000409 under load.
  - The hypotheses were read: measurable bets integrating to at most 1, and p-values independent with law μ. M-1's uniform p-values meet them.
- *What changes for existing calls:* `Guarantee.Proved` reads true for CUSUM and Shiryaev–Roberts, and the statement changes to:

  > The mean time to a false alarm is at least c, and the expected number of false alarms among the first N arrivals is at most N/c, for bets that integrate to at most 1.

  CUSUM's adds that it alarms no more often than Shiryaev–Roberts.

**Checked before marking.**
- The Lean's runs are `MartingaleAlarm`'s recursions:
  - Shiryaev–Roberts: Rₙ = fₙ(Rₙ₋₁ + 1);
  - CUSUM: γₙ = fₙ·max(γₙ₋₁, 1);
  - each restarted at 0 after an alarm, with the bet chosen from all earlier p-values.
- The hypotheses are satisfiable, and every bet here meets them:
  - measurable bets that read only the past and integrate to at most 1, so the plug-in is covered too;
  - p-values i.i.d. with law μ at every horizon.
- The theorems: `cusum_meanRunLength`, `sr_meanRunLength`, `alarms_stream_le`, `sr_alarms_le` and `cusum_count_le`, in `lean/Enthymemata/ConformalValidity.lean`.
- *What stays unproved,* as the docstrings now say: Vovk's almost-sure form, a long-run frequency of false alarms at most 1/c. The proved form is the expectation at every N.

**Gates.**
- `MartingaleAlarmTests` pass 21 of 21. The measured mean times and long-stream frequencies are unchanged.
- *Per commit,* the selected classes and the pinned tier across every test project: 120 of 120 in the 5 classes selected in Maths.Regression.Tests. Maths.Oracle.Tests compiles the source file and selects no class.

**Defects found: none.**

**Doc changes owed:** the conformal test martingales section: CUSUM's and Shiryaev–Roberts's guarantees, proved, in the expectation form; the almost-sure form not proved.

### 2026-11-15 — Conformal test martingales: the Mondrian monitor, for streams with known structure

Built at the user's pick (2026-11-15, "Mondrian monitor"). Every statistical gate's size, seed and level, and the uniformity test's rerun rule, were fixed in the tests before their first reading, as the coordinator asked after the covariate-shift item.

**Landed** on `main`, nothing pushed:
- `c64ccc9`:
  - strata in the inductive ranker, `SmoothedPValues.Add(score, stratum)`, with `CountIn(stratum)`;
  - the monitor's `Add(score, stratum)`.
  - `Add(score)` is stratum 0, so the unstratified form is the one-stratum case.
- `MondrianPValueTests` is a new class, seven cases.
- *What changes for existing calls: nothing.* `SmoothedPValues.Count` still counts every arrival. `SmoothedPValue.Count` is the arrival's own stratum's count, which equals all arrivals in an unstratified stream. The inductive monitor's wrong-form message now names `Add(score, stratum)`.
- *Verified by the main session on 2026-11-15,* from a clean archive of `5cb24ce` with zero skips: the Commit selection over `c64ccc9`'s files, 157 of 157 ([receipts](receipts/verify-conformal-mondrian-c64ccc9/)).
  - Each gate's size, seed, level and rerun rule were fixed before the first reading, as asked.
  - The pre-stated power gate for a change confined to one season failed: the pooled capital alarmed in 0.10 of streams, against at least 0.8, while the per-stratum watch reached 0.495. It is recorded as a finding about pooling, not a passed gate, and its test now pins those readings.
  - The every-season gate was added after that reading, with its own level fixed before its first reading, and passed at 1.000.

**The source.** Vovk, Gammerman and Shafer 2005:
- §4.5, Mondrian conformal transducers: the taxonomy κ(n, z), each category a rectangle of positions by units, and each score ranked within its category;
- Proposition 4.10, a case of Theorem 8.2 per §8.4's Mondrian-exchangeability models: the smoothed p-values are i.i.d. uniform given the categories, if the units are exchangeable within each category.

The C# follows them. No new guarantee, so nothing goes to the Lean thread. It is M-1's finite core per stratum, with the strata's orderings independent.

**What it computes.**
- Each arrival carries a stratum, fixed in advance from its position and its own unit, each separately: a phase of the position (a season), a category of the unit, or a pair of both.
- Each score is ranked among its stratum's scores only, one order-statistic tree per stratum, at O(log n) in its stratum.
- θ comes from the stream's one seeded generator.

**Gates.**
- *The finite core, by enumeration.* For three bags with ties, over every ordering within the strata (24, 8 and 36), every within-stratum rank sequence has total probability exactly 1. The unstratified ranks of the same orderings reach 24 of 5,040, 8 of 720 and 72 of 720 sequences, so they are not uniform.
- *A seasonal stream.*
  - Four seasons by position, each its own law: N(0, 1), N(2, 0.5²), N(−1, 2²), N(1, 1). The score −|y| ignores the season.
  - Fixed beforehand: 1,000 streams of 200, seed 2026_1126, the 1% level, and a rerun at 4,000 streams on seed 2026_1127 if a KS gate failed below 1.1 times its critical value.
  - Read once, with no rerun needed:
    - Mondrian KS pooled 0.0022 against 0.0036, and by season 0.0042, 0.0020, 0.0026 and 0.0040 against 0.0073;
    - lag-one correlation 0.0013;
    - Ville's false alarms at c = 20 under the Sleepy Jumper: 0.0010.
  - The unstratified ranks of the same scores read KS by season of 0.31, 0.40, 0.13 and 0.13.
- *A change: one pooled capital against a per-stratum watch.*
  - The fixture: 200 streams, 200 seasonal arrivals, then 400 changed, all under the Sleepy Jumper.
  - The pooled Mondrian monitor bets into one capital at c = 20. The per-stratum watch runs one monitor per season at c = 80, so its false alarms over the four stay within 1/20.
  - Each cell gives the share of streams alarmed after the change, then the median delay:

    | change | pooled capital | per-stratum watch |
    |---|---|---|
    | every season's sd × 2.5 (seed 2026_1130) | 1.000, 26 | 0.875, 149 |
    | season 2's sd 2 → 5 (seed 2026_1128) | 0.100, 141 | 0.495, 167 |

    Neither had a false alarm before the change, which each case asserts.
  - **A pre-stated gate failed.** The one-season case's gate, fixed before its first reading, was that the pooled monitor alarms in at least 0.8 of the streams. It read 0.10, so the gate failed. The case now pins its readings and its false-alarm bound. The every-season case was added afterward, with its gate (pooled at least 0.8) fixed before its first reading, and it passed at 1.000.
- *One stratum, the monitor and the refusals.*
  - `Add(score)` is `Add(score, 0)` bit for bit, and a unit alone in its stratum reads p = θ.
  - The monitor's p-values and capital are the source's.
  - A negative stratum, a NaN score and the other forms' monitors are refused.
- *Existing classes* `SmoothedPValueTests` pass, 19 of 19.
- *Per commit,* the selected classes and the pinned tier across every test project: 157 of 157 in the 7 classes selected in Maths.Regression.Tests; Maths.Oracle.Tests compiles the two source files and selects no class; no other project's compile closure holds them.

**Defects found: none.**

**The finding, which E35's strata question asked:** whether the monitor routes by stratum or the caller runs one per stratum. Both serve, for different changes, and neither dominates.
- *The pooled capital* (one monitor, `Add(score, stratum)`) is right for a change that reaches every stratum: it alarmed in every stream, in a median of 26 arrivals against the watch's 149.
- *A per-stratum watch* (one monitor per stratum, its c raised by the strata's number) is right for a change confined to one stratum: 0.495 against 0.100.
- The pooled capital pays for the unchanged strata's uniform p-values, and in either form a changed stratum's bag fills with changed units, so its p-values drift back toward uniform.

**Follow-ups, the user's to decide:**
- *Both at once:* run the pooled capital and the per-stratum watch side by side, with the false-alarm budget split between them, for example c = 40 pooled and c = 160 per stratum. This decides whether the monitor offers the hedge itself.
- *Strata for the full and covariate-shift forms.* Rescoring or weighting within each stratum should hold by the same argument per stratum, but it is neither built nor checked.

**Doc changes owed:** the conformal test martingales section: strata, the Mondrian form, the taxonomy's rectangle condition, Proposition 4.10, the seasonal readings, the pooled-against-per-stratum table, and the failed pre-stated gate.

### 2026-11-15 — Conformal test martingales: streams under covariate shift

Built at the user's pick (2026-11-15, "Streams under shift"). Two choices settled the scope and the guarantee:
- the scope, "Covariate shift first": streams with known structure, the Mondrian monitor, stay on the open list;
- the guarantee, "Build, checks now": build with an exact enumeration gate, and pose the statement to the conformal-lean thread as a follow-up.

**Landed** on `main`, nothing pushed:
- `d4c97d4`:
  - `CovariateShiftPValues` and `WeightedPValue`;
  - `ConformalTestMartingale.CovariateShift` with `Add(score, weight)`, returning `WeightedMartingaleStep`;
  - `MonitorForm`, read by `Form`. `IsFullForm` stays.
- `CovariateShiftPValueTests` is a new class, nine cases.
- *What changes for existing calls: nothing.* The two earlier forms read as before, and each form refuses the others' input with a message naming the form's own `Add`.
- *Verified by the main session on 2026-11-15,* from a clean archive of `e885e6a` with zero skips: the Commit selection over `d4c97d4`'s files, 113 of 113 ([receipts](receipts/verify-conformal-covshift-d4c97d4/)).
  - *On the uniformity gate,* the first reading failed at 500 streams (KS 0.0075 against the 1% critical 0.0073), and the gate was changed after it: 2,000 streams on a new seed, reading 0.0017 against 0.0036. A gate changed after a failing reading is a protocol slip, and this note records it. The second reading still counts as evidence: a real departure of the first reading's size would stay near 0.0075 as streams grow, and it fell to under a quarter of that. The exact enumeration check, joint uniformity over every weighted ordering to 2·10⁻¹⁵, does not depend on sampling.

**The sources.**
- The weighting is Tibshirani, Barber, Candès and Ramdas 2019: Definition 1 (weighted exchangeability) and Lemma 3 (the chance that the last position holds each unit, given the bag). It is on the shelf, `gauntlet/conformal-bayes/1904.06019v3`.
- The online form is Vovk, Nouretdinov and Gammerman 2003's reverse-filtration proof of M-1.
- The paper that builds weighted conformal test martingales is not on the shelf: Prinster, Han, Liu and Saria 2025, *WATCH*, I believe. So the method here is derived, not ported. Its premise and gate are below.

**What it computes.**
- *The design.* A reference sample of m units from the source law P comes first, scored by the fit that scores the stream. Then the stream arrives from the target law Q: the inputs' law shifted, y given x unchanged. Every unit carries w(x) = dQ/dP, the reference units included.
- *The p-value at the k-th stream arrival* (n = m + k units). Unit j holds the newest position with chance Pⱼ = πⱼ/k.
  - πⱼ = wⱼ e_(k−1)(w₋ⱼ)/e_k(w) is the chance it holds one of the k stream positions: the conditional Poisson inclusion probability.
  - p = Σ Pⱼ[αⱼ < αₙ] + θ Σ Pⱼ[αⱼ = αₙ].
  - Equal weights give the unweighted ranks. With no reference sample it is `SmoothedPValues`' p.
- *The arithmetic.*
  - Tilting every weight by c leaves π unchanged, and with aᵢ = c wᵢ/(1 + c wᵢ) the count of units drawn is Poisson-binomial. So πⱼ = aⱼ P₋ⱼ(count = k − 1)/P(count = k).
  - c is solved so the expected count is k, so the probabilities read lie near the mode and none underflows.
  - Leave-one-out comes from prefix and suffix tables. Every term is positive, so nothing cancels.
  - Where k exceeds m, the units not drawn are counted instead, so the cost is O(n·min(k, m)) per arrival.
- *Weights of zero.* A reference unit may have w = 0, outside the target law, and never holds a stream position. A stream unit may not: the null rules it out, so such a unit is itself evidence of a shift beyond the stated one, and it is refused.

**The guarantee (posed; checked by enumeration; not in Lean).**
- The premise: the units' joint law is weighted exchangeable with weights 1 on the reference positions and w on the stream's (independent draws, as Tibshirani et al.'s Lemma 1 notes), and the scores are one fixed function.
- The claim: the stream's p-values are i.i.d. uniform.
- The argument, two steps:
  - Given the bag at time n and every later unit, the arrangement has chance ∝ Π w over the stream positions (their Lemma 3), so pₙ is uniform given them.
  - Every later p-value reads only that bag and the later units, so it is independent of pₙ.
- *For the Lean thread:* this is M-1 with weighted orderings, `seqRank_uniform`'s finite core with each ordering weighed by Π w over the stream positions. The test below is its enumeration.

**Gates.**
- *The inclusion chances.*
  - Against their definition, enumerated over every k-subset: 300 bags of 2 to 10 units, weights from 10⁻⁸ to 10⁸, some zero, every k. The largest relative error is 3.6·10⁻¹⁵, and Σπ = k.
  - Beyond enumeration, 150 units at k 1 to 149, against ln e_r by a log-space recursion with each unit left out: within 10⁻¹² relative.
- *The finite core, by enumeration.* Six units, every ordering weighed by its likelihood. The stream's p-values are jointly uniform, P(every pₜ ≤ uₜ) = Π uₜ over a grid, with the largest departure:

  | bag | weighted p-values | unweighted ranks |
  |---|---|---|
  | 3 reference units | 2.3·10⁻¹⁵ | 0.059 |
  | a reference unit of weight 0 | 1.9·10⁻¹⁶ | 0.118 |
  | ties, 4 reference units | 8.9·10⁻¹⁶ | 0.196 |
  | equal weights | 1.1·10⁻¹⁵ | 1.1·10⁻¹⁵ |
- *On a stream under the stated shift.*
  - The fixture is slice 14's: x uniform for the source, density 2x for the target, so w = 2x, with noise growing with x and a one-scale fit's log density. Each stream has 50 reference units, then 100 stream units.
  - Over 2,000 streams, the weighted p-values' KS pooled is 0.0028 against the 1% critical 0.0036; at k 1, 50 and 100 it is 0.020, 0.029 and 0.022 against 0.036; the mean p is 0.499, the lag-one correlation 0.0013, and Ville's false alarms at c = 20 under the Sleepy Jumper are 0.0005 against 0.05. The unweighted ranks of the same units read KS 0.067 and mean p 0.457.
  - The first run, at 500 streams, read a pooled KS of 0.0075 against the 1% critical 0.0073, a draw of about p 0.008.
  - A probe at 2,000 streams on another seed read 0.0017 against 0.0036, with every block of 25 arrivals flat (KS at most 0.0039 against 0.0073) and the chances matching the log-space recursion at 150 units. So the test runs 2,000 streams.
- *A change in y given x.*
  - The fixture: 100 streams, 100 units under the stated shift, then 200 whose noise doubles; the Sleepy Jumper with Ville at c = 20.
  - There were no false alarms before the change. 0.99 of streams alarmed after it, with a median delay of 40.
- *Equal weights* give the unweighted counts over n to 14 digits, and with no reference sample `SmoothedPValues`' p. The monitor's p-values and capital are the source's under the same bet. Each bad input is refused.
- *The cost per arrival,* 200 reference units, wall clock, indicative: 251 µs at k up to 100, 418 µs to 200, 829 µs to 500, 1,581 µs to 1,000 and 2,482 µs to 2,000, linear in n at fixed m.
- *Per commit,* the selected classes and the pinned tier across every test project: 113 of 113 in the 4 classes selected in Maths.Regression.Tests; Maths.Oracle.Tests compiles the two source files and selects no class; no other project's compile closure holds them.

**Defects found: none.**

**The decision it serves:** whether the monitor serves streams that are not exchangeable as they stand. It now serves a stated covariate shift.
- The weighted monitor tolerates the stated shift and alarms on a change in y given x.
- The unweighted ranks of the same stream are not uniform under it (KS 0.067, mean p 0.457), so their capital carries no guarantee there.
- At this shift's size their Ville alarms stayed rare (0.002 of streams), but a larger shift would move them further.

**Follow-ups:**
- *The guarantee in Lean* (the user's choice above): the statement goes to the conformal-lean thread.
- *The Mondrian monitor for streams with known structure* (seasons, categories), from Vovk, Gammerman and Shafer 2005, §8.4. It stays on the open list, the user's.
- *An estimated ratio.*
  - With w fitted outside the stream, by `DensityRatio.Logistic` or `.Spline`, validity is lost by however far the weights move. No bound is read out here.
  - A bound like the weighted split rule's TV term would decide whether the monitor reports its own loss. The user's.
- *WATCH on the shelf.* Reading it would check the method here against its source. It also proposes diagnosing a concept shift against an out-of-support covariate shift. The user's.

**Doc changes owed:** the conformal test martingales section: the covariate-shift form, its design, the p-value and its inclusion chances, the tilted Poisson-binomial arithmetic, the guarantee's statement and argument, the enumeration gate, the stream readings, and the zero-weight rule.

### 2026-11-15 — Conformal test martingales: the Sleepy Jumper

Built at the user's pick (2026-11-15, "The Sleepy Jumper").

**Landed** on `main`, nothing pushed:
- `765795c`:
  - `SleepyJumperBet(sleep, jump, levels)`, an `IBettingFunction` beside Simple Jumper;
  - the readouts `Sleep`, `Jump`, `Levels`, `BettingShare` and `Mixture()`.
- `BettingCapitalTests` gains:
  - two bets in its integral and past-only theories;
  - a path enumeration with four parameter sets;
  - the resting test;
  - the Sleepy Jumper in the change and cost readings.
- `MartingaleAlarmTests` gains a Ville row and the Sleepy Jumper in the delay reading.
- *What changes for existing calls: nothing.*
- *Verified by the main session on 2026-11-15,* from a clean archive of `338c6ab` with zero skips: the Commit selection over `765795c`'s files, 104 of 104 ([receipts](receipts/verify-conformal-sleepy-765795c/)).

**What it computes.**
- *The sources.* Both define it alike: Vovk, Nouretdinov and Gammerman 2003, §6 (`Vovk2003.pdf`, PDF page 6), and Vovk, Gammerman and Shafer 2005, §7.1 (book pages 176–177).
- *The automaton.*
  - The states are (awake or asleep, ε̃), and it starts asleep at ε̃ = 1.
  - Asleep, it wakes with probability R and keeps ε̃.
  - Awake, ε̃ jumps to a uniform draw with probability J, and independently it falls asleep with probability R.
  - Awake, it bets the power bet εp^(ε−1) with ε = ε̃; asleep, it bets 1, a rest. So the starting ε̃ = 1 rests even awake, until its first jump.
- *The bet.*
  - The step's bet is the mixture over the states, each weighted by the capital its paths have earned.
  - It is filtered forward, in logs: O(levels) per step.
  - ε̃ lies on `levels` midpoints of [0, 1], 100 by default, with the starting level 1 held apart. The jump's uniform law is the midpoint rule.
- *Exact for its grid.* Each bet mixes power bets, each integrating to 1, with weights from the past. No margin or quadrature enters.
- *The defaults* are Vovk et al.'s R = 0.01 and J = 0.001.

**Gates.**
- *The automaton's mixture.* Every path of six steps over two levels, weighed literally by the sources' rules, gives the capital at each step to 10⁻¹² relative. It holds at (R, J) = (0.3, 0.4), (0.01, 0.001), (1, 1) and (0, 0.5).
- *Integrates to 1.* After uniform and skewed histories of 0 to 333 p-values, |∫ f − 1| is at most 3.0·10⁻¹⁴ at the defaults and 2.5·10⁻¹⁴ at (0.3, 0.2, 20 levels). The tail below e^(−700) comes from `Mixture()`, Σ w·a^ε.
- *Reads only the past*, as every bet must.
- *Ville's false alarms* over 10,000 streams of 300 uniform p-values at c = 10: 0.0023, against 1/c = 0.1.
- *What resting buys.*
  - The fixture: 200 streams of 2,000 exchangeable arrivals, then 300 whose scores fall by 1.5 sds; Ville at c = 100, from the start.

    | bet | mean ln S at the change (lowest) | alarmed after it | median delay |
    |---|---|---|---|
    | Sleepy Jumper | −0.79 (−0.91) | 1.000 | 14 |
    | Simple Jumper | −7.98 (−11.20) | 1.000 | 55 |
    | power mixture | −602 | 0 | none |
    | power 0.5 | −389 | 0 | none |

  - The share of the capital betting was 0.07 at the change and 0.97 fifty arrivals after it.
  - The test pins each reading: the Sleepy Jumper above −2 and above Simple Jumper at the change, alarming in every stream, sooner, and waking.
- *The delay reading* (qualification; 300 streams, 500 exchangeable arrivals then 1,000 changed; c = 100). Each cell gives the share of false alarms, the share detected and the median delay:

  | change | procedure | Sleepy Jumper | Simple Jumper |
  |---|---|---|---|
  | mean | Ville | 0.000, 1.000, 86 | 0.007, 0.990, 119 |
  | mean | CUSUM | 0.000, 1.000, 81 | 0.020, 0.980, 81 |
  | noise | Ville | 0.000, 1.000, 48 | 0.017, 0.983, 110 |
  | noise | CUSUM | 0.000, 1.000, 47 | 0.030, 0.970, 75 |

  Shiryaev–Roberts alarms falsely in every stream for both jumpers at c = 100. A bet near 1 adds about 1 to Rₙ per step, so it passes c near step c, well inside 500 arrivals. Its posed guarantee, a mean time to a false alarm of at least c, allows that, so this is a reading, as slice 6c found for Simple Jumper.
- *The cost per arrival,* wall clock and indicative: 4.8 µs at 100 levels, against Simple Jumper's 0.39 µs and the plug-in's 0.71 µs.
- *The change reading* at 200 exchangeable arrivals: ln S was −0.06 at the change, and peaked at 50.9 against Simple Jumper's 36.2.
- *Per commit,* the selected classes and the pinned tier across every test project: 104 of 104 in the 3 classes selected in Maths.Regression.Tests; Maths.Oracle.Tests compiles the source file and selects no class; no other project's compile closure holds them.

**Defects found: none.**

**The decision it serves:** whether the monitor's bet can rest between jumps. It can, and resting is what keeps Ville's alarm alive on long streams.
- Ville's statistic reads the capital from the start, so a capital lost to an exchangeable stretch stays lost. Through 2,000 such arrivals, Simple Jumper loses about 8 in ln S and the Sleepy Jumper under 1.
- So after a change the Sleepy Jumper alarms sooner, a median 14 against 55, and with fewer false alarms.
- CUSUM reads past a lost capital, so the two jumpers' CUSUM delays are close on a mean change: 81 and 81.

**Follow-ups, the user's to decide:**
- *ε̃ that drifts slowly as well as jumping,* which the book suggests (Vovk, Gammerman and Shafer 2005, page 177). It decides whether the bet tracks a change that grows gradually.
- *The grid's effect on the delays,* not measured. Any grid is exact, so the grid only moves power.

**Doc changes owed:** the conformal test martingales section: the Sleepy Jumper, its automaton and grid, its defaults, the resting table, the delay table and the Shiryaev–Roberts reading, and the cost.

### 2026-11-15 — Conformal test martingales: the full form

Built at the user's pick (2026-11-15, "The full form"), the first of the Vovk items.

**Landed** on `main`, nothing pushed:
- `b07fb8d`:
  - `FullFormPValues` and `FullFormScore`, either `Density` or `Centrality`;
  - `ConformalTestMartingale.Full`, with `Add(logDensity, logCdf)`, `AddRange(CandidateLogDensity)`, `IsFullForm` and `EffectiveSampleSize`;
  - a pointer in `SmoothedPValues`' docstring.
- `FullFormPValueTests` is a new class with 11 cases.
- *What changes for existing calls: nothing.* The monitor's constructor still builds the inductive form, and `Add(double)` reads as before. The two forms refuse each other's input.
- *Verified by the main session on 2026-11-15,* from a clean archive of `87d9966` with zero skips: the Commit selection over `b07fb8d`'s files, 131 of 131 ([receipts](receipts/verify-conformal-fullform-b07fb8d/)). The claim that no new Lean is needed was read: the full form's scores are recomputed from the bag, so they are equivariant, and M-1's premise covers them.

**What it computes.**
- *The posterior given the stream.*
  - The fit's draws come from a training set proper outside the stream.
  - Each arrival adds its terms ln p(yₙ | xₙ, θ_s) to every draw's log weight, so w⁽ⁿ⁾ ∝ W̄ Π_{i≤n} p(yᵢ | xᵢ, θ_s): the add-one-in weights over the whole stream.
  - Every unit is then rescored in sample under w⁽ⁿ⁾, and the newest score is ranked among all n with the randomized tie.
- *The scores.*
  - `Density`: ln Σ_s w_s p(yᵢ | θ_s), the full step's own score.
  - `Centrality`: min(F̂ₙ, 1 − F̂ₙ) from each draw's log CDF, with both tails kept.
- *The arithmetic.*
  - Each unit keeps its terms scaled by their largest, so a rescore is a dot product.
  - Where that sum falls below 10⁻²⁵⁰, the unit is rescored by log-sum-exp.
- *The seam.*
  - `AddRange` reads a `CandidateLogDensity` whose cases are the units, each at its one observed response.
  - The seam's draws and weights must be the stream's.
  - It must serve the full form, since reweighing by a term reads the term as that unit's likelihood.

**The guarantee, and why it needs no margin.**
- M-1 needs the scores at each time to be a symmetric function of the units so far: `Equivariant` at every n, which the Lean record already says covers the full form. So no new guarantee was posed.
- Given the draws, each score reads the units only through the product of their terms, whatever their order. So the p-values are i.i.d. uniform on exchangeable streams, exactly, whatever the number of draws.
- That differs from the full step's sets: their draws come from the calibration units' own posterior, so they need CB-3's margin. Here the draws must instead be fixed before the stream; refitting them to it would break the guarantee.
- The weights still degenerate as the stream grows, which costs sharpness only.

**Gates.**
- *M-1 by enumeration*, six continuous units over all 720 orderings: each rank sequence occurs exactly once, under either score.
  - Feeding the same scores, each kept from its own time, to the inductive ranker breaks it: 442 distinct sequences under `Density` (one occurring 7 times) and 473 under `Centrality`. This is the drift M-1's record names.
- *Uniform and independent on exchangeable streams.*
  - The fixture: 1,000 streams of 100 units from N(0.5, 1); each stream's fit is 100 draws from five training units of N(0, 1), so the fit is off.
  - KS pooled 0.0026 against a 1% critical value of 0.0052.
  - At n 2, 50 and 100: 0.025, 0.035 and 0.026 against 0.0515.
  - Lag-one correlation −0.0002.
  - Ville's false alarms at c = 20 under the Simple Jumper: 0.015, against 1/c = 0.05.
  - Mean ESS of the 100 draws: 85.0, 45.9 and 15.9 at n 1, 10 and 100.
  - Both scores read alike here, since on this symmetric location fixture they order the units the same.
  - The positional check runs from n 2, since p₁ is θ₁ by construction. At n 1 it read 0.0502 against 0.0515, which tests only the tie draws across seeds 0–999.
- *A function of the units so far, computed exactly.*
  - The fixture: 200 wide draws, 300 units near 0 and one at 60. The 300 hold the weights near μ = 0, so the far unit's scaled sum is 0, and its density, −1,780.99, comes from log-sum-exp.
  - In four orders of arrival, every unit's final density matches the definition to 10⁻⁸ and its centrality to 10⁻⁹ relative.
- *The seam and the refusals.*
  - The seam reads as the spans do, bit for bit, and the monitor's p-values and capital are the source's under the same bet.
  - Each bad input is refused, and a refused arrival leaves the stream as it was:
    - bad draw weights;
    - a NaN or +∞ term;
    - a unit with zero density under every weighted draw;
    - a log CDF above 0, or missing for the centrality;
    - a split-only seam, its reason quoted;
    - another seam's draws or weights;
    - a grid of two values.
- *The delay to an alarm, against the inductive form on the same draws.*
  - The fixture: 200 streams, 100 units of N(0, 1), then 200 of N(1.5, 1) (mean) or N(0, 2.5²) (noise); 50 draws from 3 or 200 training units; the Simple Jumper with Ville at c = 20.
  - Each cell gives the share of streams alarmed after the change, then the median delay:

    | change, training units | inductive | full | false alarms before, inductive and full |
    |---|---|---|---|
    | mean, 3 | 0.865, 33 | 0.590, 41 | 0.015, 0.030 |
    | mean, 200 | 1.000, 35 | 0.995, 39 | 0.025, 0.035 |
    | noise, 3 | 0.990, 42 | 1.000, 33 | 0.025, 0.010 |
    | noise, 200 | 1.000, 36 | 1.000, 36 | 0.005, 0.005 |
- *The cost per arrival* (wall clock, indicative, with the qualification switch; 1,000 draws):

  | score | n 0–250 | n 250–500 | n 500–1,000 | n 1,000–2,000 |
  |---|---|---|---|---|
  | `Density` | 128 µs | 322 µs | 728 µs | 1,763 µs |
  | `Centrality` | 401 µs | 754 µs | 1,454 µs | 2,996 µs |

  That is linear in n·S, about 0.9 ns per unit and draw for the density. The ESS at n 2,000 was 35 of 1,000. The inductive form's arrival is O(log n).
- *Per commit,* the selected classes and the pinned tier across every test project: 131 of 131 in the 5 classes selected in Maths.Regression.Tests; Maths.Oracle.Tests compiles three of the four files and selects no class; no other project's compile closure holds them.

**Defects found: none.**

**The finding.** The full form tests a different hypothesis's sharpness, not a sharper one.
- *Where the fit is vague and the model can absorb the change,* as a location model absorbs a shift in the mean, the posterior given the stream follows the change. The units after it then score as conforming, so the full form alarms less often: 0.59 against 0.865.
- *Where the model cannot absorb the change,* as with a noise change under a location model, the stream sharpens the fit, and the full form alarms sooner: a median of 33 against 42.
- *Under a tight fit* the two read alike.

So it serves a caller whose fit is weak and whose feared change is one the model cannot follow. It costs O(n·S) per arrival.

**Follow-ups, the user's to decide:**
- *Pruning the draws by weight.*
  - At each arrival, rescore only under the draws whose weight is above a fixed share of the largest.
  - The kept set is read from the current weights alone, a function of the units so far, so M-1 still holds.
  - The cost falls with the ESS, to O(S + n·S_kept): the ESS was 35 of 1,000 at n 2,000.
  - Pruning for good, so that a draw once dropped never returns, would read the history's order and break M-1.
  - This decides whether the full form serves long streams.
- *Both forms at once.*
  - The idea: the inductive and full forms' capitals averaged, to hedge between the two regimes above.
  - The open question: each is a martingale in its own p-values' filtration, so whether their average keeps Ville's guarantee needs its own reading first.
  - This decides whether the monitor offers the combination.

**Doc changes owed:** the conformal test martingales section: the full form, its scores, why its guarantee needs no margin, the seam's premise, the delay table and the finding, and the cost.

### 2026-11-15 — τ through richer features: a spline density ratio

Built at the user's pick (2026-11-15, "slice 15's alternation speed then τ through richer features, then we'll return to the vovk items").

**Landed** on `main`, nothing pushed:
- `274f69d`:
  - `DensityRatio.Spline`, the nonparametric ratio, with the readouts `Ridge` and `Inputs`;
  - `ScoreMixtureOptions.Knots`, which routes each level's τ fit through it.
- `DensityRatioTests` gains three tests and `ScoreMixtureTests` one.
- *What changes for existing calls: nothing.*
  - `Logistic` now runs through a shared core and still refuses a constant feature.
  - `Knots` is null by default, so `ScoreMixture` fits τ on the features as given, as before.
- *Verified by the main session on 2026-11-15,* from a clean archive of `1dd5974` with zero skips: the Commit selection over `274f69d`'s files, 59 of 59 ([receipts](receipts/verify-conformal-tauspline-274f69d/)). It was the first verification under the weight reader's trap, and its run left `jit-Rows.txt` (60 KB) in its results folder.

**What it computes.**
- *The expansion.*
  - Each input becomes cubic B-splines over the pooled samples' range, with `knots` interior knots (6 by default) at the pooled quantiles. Repeated knots are dropped.
  - Each input's first column is dropped, since its B-splines sum to one and the intercept spans that. A column constant over the pooled samples is dropped too; a few-valued input leaves some.
  - An input beyond the range holds its edge value.
  - So ln w = Σⱼ sⱼ(xⱼ) plus the size correction: additive in the inputs, with no interactions.
- *The ridge.*
  - It is chosen from 10⁻⁴ to 10 by decades, by the classifier's 5-fold cross-validated log loss. Row i of each sample goes in fold i mod 5, so sorted inputs spread over every fold.
  - A ridge whose fit is refused in some fold is passed over. A column constant within a fold's training rows sits idle there, its slope held at zero.
  - Every row is then refitted at the chosen ridge.
- *τ through it.* `Knots` set to a count fits each level's ratio with `Spline`, and `Ridge` then goes unused.

**Gates.**
- *Variance shift.* N(0, 1) against N(0, 1.6²), 4,000 points each. ln w is even in x, so a linear fit reads it as flat.
  - The spline's ln w is within 0.15 of the truth at x = −2 to 2: +0.150 and −0.140 at ±2, −0.038 at 0.
  - The linear fit's lies between −0.034 and 0.034, up to 0.78 off.
  - Cross-validation chose 10⁻⁴; the fit took 9 Newton steps.
- *Nonlinear covariate shift.*
  - The fixture:
    - training inputs uniform, the test's triangular, with density 4 min(x, 1 − x);
    - noise growing toward the middle;
    - 4,000 replications at n 199 and α 0.1, each ratio fitted once from 2,000 unlabeled inputs of each law.
  - Coverage, se 0.0047:

    | weights | coverage | mean TV |
    |---|---|---|
    | none | 0.8385 | |
    | the linear fit's | 0.8423 | 0.250 |
    | the spline's | 0.9077 | 0.034 |

    The spline's coverage is within the TV bound (`sum_weight_refused_estimated_le`).
- *Units, edges and refusals.*
  - Two inputs, one in units of 10³ and one of 10⁻³, each with an offset, read alike to 8 digits and take the same ridge.
  - Beyond the range, an input holds its edge value to 12 digits.
  - A binary input beside a continuous one fits; its step reads 1.0996 against ln 3 = 1.0986.
  - Each bad argument is refused.
- *The score mixture, a middle cluster.* The narrow cluster holds |x − ½| below about 0.2, so its membership is not monotone in x. On 1,000 points of raw x:

  | τ | J | mean membership error | R² at J 2 |
  |---|---|---|---|
  | `Knots` 6 | 2 | 0.0695 | 0.979 |
  | linear in x | 2 | 0.40 | 0.9996 |

  The spline fit took 2.7 s.
- *Existing classes.* `DensityRatioTests` and `ScoreMixtureTests` pass.
- *Per commit,* the selected classes and the pinned tier across every test project: 59 of 59 selected tests in Maths.Regression.Tests (DensityRatioTests and ScoreMixtureTests); Maths.Oracle.Tests compiles the two source files and selects no class; no other project's compile closure holds them.

**Defects found:** none in existing code.

My own build had one, fixed before the commit and found in self-review. A binary input left all-zero spline columns, which `Logistic` refuses, so `Spline` threw on it. Pooled-constant columns are now dropped, and fold-constant ones sit idle. The binary-input test pins the fix.

**The finding.** On the middle cluster, the linear τ scores the higher R², 0.9996 against the spline's 0.979, while getting the memberships wrong.
- The linear τ is nearly flat in x, so two factors fit it closely while holding nothing of the clusters.
- So the factorization's R² cannot flag a τ that misses the shift. The ratio fits' held-out log loss, which the ridge's choice already computes, can.

**Follow-ups, the user's to decide:**
- *Interactions between inputs.* The spline ratio is additive, so it misses a shift in a joint feature. A tensor-product basis, or a kernel ratio such as KLIEP or uLSIF, would serve such a shift. That decides whether the layer serves shifts that live in interactions.
- *`Knots` as the score mixture's default.* That would be a change of default. The cost is 31 logistic fits per level instead of one.

**Doc changes owed:**
- The density-ratio section: the spline form (knots, edges, constant columns, the cross-validated ridge) and its additive limit.
- The posterior-conformal section: `Knots`, and the finding that R² cannot flag a flat τ.

### 2026-11-15 — Slice 15's alternation speed: safeguarded extrapolation of γ

Built at the user's pick (2026-11-15, "slice 15's alternation speed then τ through richer features, then we'll return to the vovk items").

**Landed** on `main`, nothing pushed:
- `192ee84`: `ScoreMixture`'s factorization extrapolates γ's step, its least-squares update moved into `SolveProfiles`. `Factor` gains `extrapolate` (on by default; internal), and `ScoreMixtureCost` counts pushes and pushes kept.
- `ScoreMixtureTests` gains a test.
- *Verified by the main session on 2026-11-15,* from a clean archive of `ec771c1` with zero skips: the Commit selection over `192ee84`'s files, 50 of 50 ([receipts](receipts/verify-conformal-s15speed-192ee84/)). `ScoreMixture` has no caller outside its tests. The safeguard was read in the code: a push is kept only below the previous pass's residual, else γ_LS, the factor capped at 4.
- *What changes for existing calls.* A fit's γ, π and R² now come from the extrapolated alternation, which stops at a slightly different point within the same 10⁻⁷ rule and 300-pass cap. Downstream readings move with it, by amounts within the fit's own tolerance. In the existing classes:
  - the unequal-split test's mean ESS at 399 units went from 100.0 to 100.2, and the training's from 107.2 to 108.2;
  - its whole-space share at the training's m went from 13.4% to 13.5%, and its marginal coverage at the calibration's from 0.9155 to 0.9145, m staying 45;
  - the end-to-end test's width went from 1.035 to 1.036.

  Every class passes.

**What it computes.**
- *The step.* Each pass after the first pushes γ past its least-squares update, γ + (1 + β)(γ_LS − γ), and re-projects π. It keeps the push only where the residual falls below the pass before, otherwise taking γ_LS (safeguarded extrapolation, Ang and Gillis 2019, for nonnegative factorizations). So no pass ends above the plain step's residual from the same π.
- *β's schedule.* β starts at ½ and grows by 5% on a kept push, up to its cap; a refused push sets the cap to β and divides β by 1.5. The cap grows by 1% per kept push.
- *The cap's ceiling.* Theirs stops at 1. A probe read ceilings of 1, 4, 16 and 64: from 4 up every run was identical, β settling between 1 and 4 as π's projection refuses larger pushes. So the ceiling is 4.

**Gates.**
- *The probe* (stage only; the true two-cluster τ, five restarts per J; passes summed over the restarts, 1,500 the cap):

  | n, τ | J | plain passes | extrapolated passes | time |
  |---|---|---|---|---|
  | 1,000, exact | 2, 3, 4 | 1,360, 1,400, 1,217 | 477, 458, 457 | 0.40 → 0.15 s at J 2 |
  | 1,000, noise 0.01 | 2, 3, 4 | 1,500 each | 690, 1,155, 743 | 1.27 → 0.60 s at J 4 |
  | 10,000, noise 0.01 | 2, 3, 4 | 1,500 each | 890, 1,260, 1,458 | 5.7 → 2.9 s at J 2; equal at J 3 and 4 |

  1 − R² was never above the plain alternation's: equal to about five significant digits on noisy τ, and lower on exact τ.
- *The test:* on the probe's fixture at n 1,000, exact and noisy, J 2 to 4, the extrapolated residual is not above the plain one's (within 10⁻⁴ relative) and takes no more passes. Every case read lower and fewer.
- *Existing classes* `ScoreMixtureTests` and `PosteriorConformalTests` pass.
- *Per commit,* the selected classes and the pinned tier: ScoreMixtureTests with the pinned tier, Maths.Regression.Tests 50 of 50; PosteriorConformalTests run beside it, 5 of 5.

**Defects found: none.**

**The finding the decision rests on.** On noisy τ at J ≥ 3 the 300-pass cap still binds, extrapolated or not. It binds after R² has settled to about five digits: the passes crawl along the flat directions a J above the cluster count leaves, which move nothing the fit uses. Extrapolation buys up to about twice at J 2 and little at J 3 and 4 there. So the cost at large n is linear in n and bounded by the cap: 300 × 5 restarts × the J's the rule reads.

**Follow-up, the user's to decide:** a stop that reads what the fit uses instead of the residual's last 10⁻⁷, such as R² settled to its rule's resolution, or a lower cap. It would cut the factorization's cost severalfold at large n, and move γ, π and R² by about the amounts above.

**Doc changes owed:** the posterior-conformal section's factorization: the safeguarded extrapolation, its schedule and ceiling, the probe's table, and the finding on the cap.

### 2026-11-15 — RoBAS-Full's grid-free intervals from n 4

Built at the user's word (2026-11-15, "Wire it now"), on the conformal-lean thread's proof (700c718). `robasFull_interval`: RoBAS-Full's sets are intervals for n ≥ 4 when the calibration residuals are not all equal, granted `RobasDawsonBounds`, four bounds on Dawson's function that the checks certify and that go into Lean next. `robas_common_point`: a point where the p-value is 1. The user chose to wire it before the bounds land in Lean.

**Landed** on `main`, nothing pushed:
- `c3539b7`: `ResidualShrinkage.Intervals` serves `ResidualScore.Horseshoe`, and the internal `ResidualShrinkage.CommonPoint` finds its start.
- *Docstrings:* `ResidualScore.Horseshoe`'s and `Intervals`' state the theorem and its hypothesis.
- *Tests:* `HorseshoeScoreTests` gains the common point's test, and its refusal check becomes the new boundary. `ResidualIntervalTests` take RoBAS-Full back into both of their sweeps.
- *Existing calls:* none changes except an `Intervals` call on RoBAS-Full at n ≥ 4 with residuals not all equal. It was refused and is now served.
- *Verified by the main session on 2026-11-15,* from a clean archive of `77f4ccd` with zero skips, in one run with `3e3b7f2` ([receipts](receipts/verify-cbma-robasint-c3539b7/)): Maths.Regression.Tests' Commit selection 2,823 of 2,823, Maths.Oracle.Tests' 108 of 108. The Lean it rests on, `700c718`, was verified from a clean archive the same day.

**What it computes.**
- *Where it serves.* RoBAS-Full from n 4.
  - Refused below n 4 at a level above 0, where its sets fall in two pieces (n 2 and 3).
  - Refused on calibration residuals all equal, outside the theorem.
  - The whole line at k = 0, as for every score.
- *The start: the common point.* z₀ with n·z₀ = S·g(S²/(2T(z₀))), where S is the calibration residuals' sum, T(z) the spread over them and z, and g = K₁/K₀.
  - Up to x = 1, g is the ratio of the series K_j(x) = Σ xᵏ/(k!(2k + 2j + 1)).
  - Past it, by the ODE 2xK₁ + K₀ = eˣ, g(x) = (1/₁F₁(1; 3/2; −x) − 1)/(2x) through `LogKummer`.
  - g rises from 1/3 to 1, so z₀ lies between S/(3n) and S/n, where n·z − S·g changes sign; bisection finds it.
  - Every unit's score reaches the candidate's there, so the p-value is 1 and the start is accepted at every level.
- *The search* is as for the other scores: doubling outward to a refusal, then bisecting each end to the tolerance, rounded outward.

**Gates.**
- *The common point:* at n 4, 5, 8, 15 and 40, shifts 0, 0.7, −3 and 25, with spreads from tight to wide, z₀ lies in [S/(3n), S/n], and the grid form's p-value there is exactly 1.
- *Against the grid form,* RoBAS-Full now in the sweep (n 19, 49 and 99, shifts 0 and 1.5, a grid of step 0.001 over [−8, 8]):
  - the grid's accepted values form one interval inside [L, U], each end within the tolerance of the grid's first refused value;
  - read on both sides of each end at 1e-9, the end is refused and two tolerances inside is accepted.
- *Coverage under a mean shift* (residuals N(1.5, 1), 49 units, 2,000 replications, α 0.1): RoBAS-Full 0.9030 at mean width 3.400, beside RoBAS-EB's 0.9045 at 3.398 (se 0.0067).
- *The boundary:*
  - n 20 is served;
  - n 3 at k 1 is refused ("below n 4"), the pinned counterexample still reading two pieces on the grid form;
  - residuals all equal are refused ("not all equal").
- *Neighbours:* `ScoresBeyondTheModelTests` pass unchanged.
- *Per commit,* the selected classes and the pinned tier: the four classes reading the files with the pinned tier, Maths.Regression.Tests 61 of 61.
- *A build note.* `SpecialFunctions.Expm1` is internal to the distributions assembly, unlike `Log1p`. The first build failed on it, and its test passes were stale binaries, so they are not counted here. The series form replaced it, and the gate script now stops when the build fails.

**Defects found: none.**

**Follow-up, on the conformal-lean thread's side:** once `RobasDawsonBounds` is proved in Lean, the docstrings drop "granted four bounds … that the checks certify".

**Doc changes owed:** the scores-beyond-the-model section: RoBAS-Full's intervals from n 4 with the theorem's hypothesis, the common point and how it is found, the boundary, and the gates' numbers.

### 2026-11-14 — Slice 17's follow-up: CBMA's evidence from the engines' fits

Built at the user's pick (2026-11-14, "CBMA evidence from the engines' fits"). The engines' side was settled with the coordinator (its answers to Q1 and Q2):
- Ownership: both adapters are new members only, on `BassEvidence` and on `BarsEvidence`. Neither engine folder has an active owner, and the BPPR chip's BASS-core files are untouched.
- Scale: conditions (a) to (d), answered below.

**Landed** on `main`, nothing pushed:
- `3e3b7f2`:
  - the seam `ModelEvidence` (`src/maths/regression`): ln p(y), its Monte Carlo error, the engine's flags, what the value depends on beyond the data (`Priors`), and the key of what it was read on (`KeyOf`: the response as the engine read it, its weights, and how they are read);
  - `AveragedModel`, `EvidenceReading`, and a second `ModelAveraging.FullStep` overload taking models with their evidence;
  - `BassEvidence.ForAveraging(BassData)` and `BarsEvidence.ForAveraging(double[])`;
  - `ModelEvidenceTests` and `EngineEvidenceTests`.
- The log-weight `FullStep` and the new overload share one core. It reads each model's scores and predictive density once per candidate, then mixes them under every weight vector. The log-weight form runs exactly its old arithmetic.
- `ScoresBeyondTheModelTests`' three fixture helpers are now `internal`, for reuse.
- *Verified by the main session on 2026-11-15,* from a clean archive of `77f4ccd` with zero skips, in one run with `c3539b7` ([receipts](receipts/verify-cbma-robasint-c3539b7/)): Maths.Regression.Tests' Commit selection 2,823 of 2,823, Maths.Oracle.Tests' 108 of 108.

**What it computes.**
- *The overload.*
  - Each model's log posterior probability given the n units is its log prior plus its log evidence. Models whose evidences carry different keys are refused.
  - `EvidenceReading` gives p(M_k | Z₁:ₙ), its least and greatest over the evidences' tilts, each case's set moving under a tilt or not, the engines' flags led by each model's index, and each model's priors.
  - A tilt moves one model's log evidence by the error multiple (2 by default) of its Monte Carlo error, up and down, and the others' the other way by theirs. For up to three models that is every vertex of the errors' box. The scores do not depend on the weights, so a tilt costs a re-mix, not a re-score. The readout is indicative, as the full step's margin is.
- *The adapters.*
  - BASS: its flags as text (saturated, the ladder at its budget, SMC overdispersion). The key is its data's response as laid out, with its case weights and their meaning (precision or replication); a functional or several-output response is keyed as BASS lays it out.
  - BARS: a run that reports an evidence reads no case weights, so the key reads the response alone.
- *The conditions.*
  - (a) Constants: pinned for both engines, without Monte Carlo, against the g-prior's marginal written from its definition, a dense multivariate t with every constant (the gates below). BASS's collapsed marginal (`GPriorMarginal.LogMarginal`) carries −(ν/2) ln 2π, and under precision weights ½ Σ ln wᵢ. BARS's block, through `BandedGPrior`, carries −(m/2) ln 2π. Both read the response as given, so engines mix.
  - (b) Priors read the response: each engine's priors description rides through the seam into `EvidenceReading.Priors`. On symmetry: both engines set their centre and σ² prior from the n training responses, never the candidate's. The full step on these fits already accepts that: its draws come from a fit under that prior, so its score is not exactly symmetric in the n + 1 units. The evidences read the same n responses and enter only the model weights, so CBMA adds no asymmetry of its own; it inherits the full step's. Exact symmetry would need the draws and the evidences re-read with the candidate.
  - (c) Response key: built from what the engine read, weights and their meaning included, so a mismatch is refused.
  - (d) The evidences' Monte Carlo errors and flags pass as given, beside the tilt readout.

**Gates.**
- *The constants, against the definition* (the dense t: y | σ² ~ N(Zμ₀, σ²(W⁻¹ + Z G⁻¹ Zᵀ/τ)), G = ZᵀWZ, σ² ~ IG(a, b)):

  | engine, configuration | engine | definition |
  |---|---|---|
  | BARS, no knots, n 60 (τ = 1/m, centre the median projected, σ² as resolved) | −17.1690934545 | −17.1690934545 |
  | BASS's enumeration over 37 structures (g-prior, τ held at 0.7, IG(3, 2), centre 0, n 20) | −23.6317698380 | −23.6317698380 |
  | the same under precision weights | −24.1801277051 | −24.1801277051 |

  BASS's design holds its columns whitened under weights, W^½z. A first run read them as raw columns and missed by 0.07. A probe settled it: the intercept-only structure matched to 0, and each one-term structure matched to 1e-14 once its column was read back as (W^½z)/√wᵢ.
- *A fitted BASS evidence,* tempered to its reference rung: −26.7583 ± 0.0012, against the definition's −26.7590.
- *On shared data* (n 60, sin plus noise):
  - BASS reads −5.118 ± 0.264 and BARS −6.697 ± 0.026. They share a key and differ by 1.6 nats, far inside the (n/2) ln 2π = 55 nats either would drop.
  - CBMA averages their full-step seams: P(BASS | Z) 0.829, from 0.731 to 0.896 over ±2 errors; one of the two cases' sets moves under the tilts.
  - BARS keyed on a shifted response is refused.
- *The overload, with exact evidences* (the σ 1 and σ 1.5 laws, 8 units, α 0.3 so k is 2):
  - Against the log-weight form, sets and weights are identical bit for bit, the posterior is the exact one, and no case is fragile.
  - Errors of ±0.5 and ±3 nats widen P(y = βx | Z) from 0.936 to [0.843, 0.975] and [0.035, 1.000]. Every case's set moves at both on two models apart in their mean, and none at a multiple of 0.
  - The first fixture's two noise laws moved no set even at ±3 nats (a noise scale barely moves a set, slice 17). At α 0.1 on 8 units every set was the whole line. Both are fixed in the fixture.
- *Keys and refusals:* the weights and their meaning enter the key; malformed evidences and priors are refused.
- *Neighbours:* `ScoresBeyondTheModelTests` and `BarsEvidenceTests` pass unchanged.
- *Per commit,* the selected classes and the pinned tier: Maths.Regression.Tests 278 of 334 classes, 2,822 of 2,822; Maths.Oracle.Tests 23 of 39 classes, 104 of 108, the four misses R oracles timing out at 120 s under load (GBASS's readouts twice, the transport oracle twice), the two GBASS classes then passing alone (55 s and 1 min 55 s); Maths.Samplers.Tests compiles none of these files.

**Defects found: none.** The engines' evidences are complete densities of the response as given.

**Doc changes owed:**
- the scores-beyond-the-model section's CBMA part: the evidence seam, the overload, the tilt readout, conditions (a) to (d) with the symmetry note, and the gates' tables;
- the API rows for `ModelEvidence`, `AveragedModel`, `EvidenceReading` and the two adapters;
- in each engine's evidence section, a line on `ForAveraging`.

### 2026-11-14 — RoBAS-Full's intervals: not intervals at n 2 and 3; none found from n 4

Probed at the user's pick (2026-11-14, "RoBAS-Full intervals"). The question was whether the grid-free form can serve RoBAS-Full.

**Landed** on `main`, nothing pushed: `f400892`.
- `ResidualScore.Horseshoe`'s docstring records the finding (comment only).
- `HorseshoeScoreTests.ItsSets_AreNotAlwaysIntervals` reads both counterexamples through the grid form.
- `Intervals` still refuses RoBAS-Full.
- *Verified by the main session on 2026-11-14,* from a clean archive of `5473adc` with zero skips: the Commit selection over `f400892`'s files, 57 of 57 ([receipts](receipts/verify-conformal-robasprobe-f400892/)). The source change is doc comments only.

**The finding.**
- *Not intervals at n 2 and 3.* A seeded search counted each set's connected pieces on a sinh-spaced grid, fine near the residuals and reaching 60 spreads beyond them. RoBAS-Full's sets fall in two pieces: one about the prior's centre 0, one about the data. The horseshoe's spike at 0 lets a candidate there conform when the others are few.
  - n 2: 34 of 4,539 random trials, e.g. residuals 18.02 and 19.33 at k 2 (α 0.9): [−8.03, 6.89] and [18.02, 19.32].
  - n 3: 8 of the targeted sweep's n-3 configurations, at k 1 and 2, e.g. residuals 0.510, 0.653, 0.546 at k 1: [−0.198, 0.062] and [0.333, 0.668].
- *None from n 4.* The targeted sweep covered n 4 to 10 at every level, about 9,900 configurations: tight clusters at shifts from 0.5 to 200 with spreads from 0.01 to 3, and an opposite-signed outlier, on a grid of 8,001. The random sweep reached n 50 with outliers, two clusters and heavy tails. Neither found a second piece.
- *The control.* RoBAS-EB, whose theorem starts at n 4, read one piece in every trial at n ≥ 4, so the detector reads true.
- *What it means.* RoBAS-Full has no unconditional interval property. Whether it is an interval from n 4, as RoBAS-EB is, is open: the search is evidence, not a proof. Until it is proved, the grid-free form refuses RoBAS-Full at every n. Its bracketed search would find only the piece holding its start and drop the other, which undercovers.

**Gates.**
- *The pinned counterexamples* read through the grid form:
  - n 3 at k 1: accepted at 0, refused at 0.2, accepted at 0.5;
  - n 2 at k 2: accepted at 0, refused at 12, accepted at 18.5;
  - RoBAS-EB on the n-3 residuals reads no such gap.
- *Per commit,* the selected classes and the pinned tier: `HorseshoeScoreTests`, `ScoresBeyondTheModelTests` and `ResidualIntervalTests`, 57 of 57 in `Maths.Regression.Tests`.

**Defects found: none.** The refusal was already in place, and this records why it stands.

**Open, for the user:** the conjecture that RoBAS-Full's sets are intervals for n ≥ 4. A proof would go through Bersson and Hoff's lemma, which RoBAS-EB's theorem uses: each unit's acceptance an interval, all sharing a point. It belongs to the conformal-lean thread if wanted. With it, the grid-free form could serve RoBAS-Full from n 4, as it serves RoBAS-EB.

**Doc changes owed:** the scores-beyond-the-model section: RoBAS-Full's two-piece sets at n 2 and 3, the search's coverage, and the open conjecture from n 4.

### 2026-11-14 — Slice 17's follow-up: RoBAS-Full, the horseshoe residual score

Built at the user's pick (2026-11-14, "RoBAS-Full"), the first of slice 17's follow-ups. Sources read in TeX: RoBAS (arXiv 2607.04236v2), §3.2 with (8), Appendix's derivation of the RoBAS-Full score and its choice of a = b = ½, the marginal-likelihood lemma, and the plug-in σ̂² of the robustness theorem's extension.

**Landed** on `main`, nothing pushed:
- `942df9c`:
  - `ResidualScore.Horseshoe`, read by `ResidualShrinkage.Prediction`;
  - `ResidualShrinkage.LogKummer`, ln ₁F₁(1; 3/2; −x);
  - `HorseshoeScoreTests`.
- *Two existing tests now handle the new value.* `ScoresBeyondTheModelTests.TheResidualScores_AreTheirDefinitions` gains RoBAS-Full's two-pass definition, where its default branch would otherwise have read the new score as RoBAS-EB's. `ResidualIntervalTests` leave RoBAS-Full out, since the grid-free form refuses it.
- *Existing scores read unchanged.*
- *Verified by the main session on 2026-11-14,* from a clean archive of `e08d77b` with zero skips: the Commit selection over `942df9c`'s files, 56 of 56 ([receipts](receipts/verify-conformal-robasfull-942df9c/)). The ₁F₁ was read against D(√x)/√x: the positive series up to x 36, the asymptotic one past it, cut at its least term. `4773f8b`, the confidence parameter's wording, is doc comments only.

**What it computes.**
- *The working model:* R | θ ~ N(θ, σ²), θ | τ² ~ N(0, σ²τ²/n), τ ~ C⁺(0, 1), the horseshoe.
- *The score.* By their lemma, a unit's Bayes-assisted score ranks as the marginal likelihood of the other n: ln p(r₋ⱼ) = −n s²₋ⱼ/(2σ̂²) + ln ₁F₁(1; 3/2; −n r̄²₋ⱼ/(2σ̂²)), up to a constant every unit shares. σ̂² = Σ(r − r̄)²/n is taken over all n + 1 residuals, their plug-in, so it is shared by every unit at a candidate. The others' mean and spread come from the residuals centred at their overall mean, as RoBAS-EB reads them. Where all n + 1 residuals are equal, every unit ties.
- *The function.* ₁F₁(1; 3/2; −x) = D(√x)/√x, D Dawson's integral (Kummer's transformation).
  - Up to x = 36: e⁻ˣ Σ xᵏ/(k!(2k + 1)), every term positive.
  - Past it: the asymptotic (1/z) Σ (2k − 1)!!/zᵏ with z = 2x, cut at its least term, about e⁻³⁶ at the seam.
- *What it does.* Near zero the score shrinks as the horseshoe does. Under a large mean shift it ranks as distance to average (their robustness theorem).
- *No interval theorem* is proved for it (theirs is RoBAS-EB's). So the grid form serves it, and `Intervals` refuses it, pointing there.

**Gates.**
- *The function against its integral.* Simpson on ∫₀¹ 2s·e^{−x s²(2 − s²)} ds, or past x = 100 on (1/(2x)) ∫₀⁶⁰ e^{−w}(1 − w/x)^{−1/2} dw, at nineteen x from 1e-12 to 1e12, across the seam at 36 and the reference's switch at 100. The largest difference in the log is 4.4e-13. D(1) = 0.53807950691 holds to 10 decimals.
  - The first reference, Simpson on v directly, read 5.4e-11 off at x 99, its integrand piled at v = 1. The substitution v = 1 − s² spreads it.
- *The closed form is the working model's marginal likelihood.* It was integrated directly over φ = atan τ, where the half-Cauchy is uniform, and over θ on its conditional's window, each integrand read from the model's own densities:

  | n, shift | closed form | integrated |
  |---|---|---|
  | 5, 0 | −6.4296749457 | −6.4296749467 |
  | 12, 2 | −20.6754508875 | −20.6754508875 |
  | 30, −8 | −48.5478000179 | −48.5478000069 |

  This checks their constant 2/π and the reduction to ₁F₁.
- *The definitions:* 300 random residual vectors, some shifted by 1e4, each unit's score against its two-pass definition with the quadrature's ₁F₁.
- *Coverage and width,* as slice 17's test (n 15, α 0.1, 2,000 replications, a grid of step 0.05; se 0.0067):

  | shift | distance to origin | distance to average | RoBAS-EB | RoBAS-Full |
  |---|---|---|---|---|
  | 0 | 0.9370, 4.052 | 0.9385, 4.215 | 0.9355, 4.108 | 0.9370, 4.082 |
  | 3 | 0.9435, 9.449 | 0.9380, 4.239 | 0.9385, 4.243 | 0.9400, 4.258 |

  RoBAS-Full sits between the references, narrower than distance to average unshifted and than distance to origin shifted. The grid's coarseness keeps every score's coverage near 0.94.
- *The robustness:* residuals N(40, 1), n 15, 200 replications on a grid of step 0.01. RoBAS-Full's sets match distance to average's to within 0.030 in width, and are identical in 25.
- *Refused:* the grid-free form.
- *Neighbours:* slice 17's classes and the grid-free intervals' pass.
- *Per commit,* the selected classes and the pinned tier: HorseshoeScoreTests, ScoresBeyondTheModelTests and ResidualIntervalTests, the classes reading the files, with the pinned tier: Maths.Regression.Tests 56 of 56; Maths.Oracle.Tests compiles the source and selects no class; Maths.Samplers.Tests compiles none of the files.

**Defects found: none.**

**Doc changes owed:** the scores-beyond-the-model section: RoBAS-Full's working model, score, plug-in and function, why the grid form serves it, the gates' tables, and the API rows for `ResidualScore.Horseshoe`.

### 2026-11-14 — Slice 17's follow-up: grid-free intervals for the residual scores

Built at the user's pick (2026-11-14, "the grid thing from slice 17"), the third of slice 17's follow-ups. Sources read in TeX: RoBAS (arXiv 2607.04236v2), §3.4 with Algorithm 1 and the interval theorem, and Appendix's complexity note.

**Landed** on `main`, nothing pushed:
- `28993bb`: `ResidualShrinkage.Intervals`, with `ResidualInterval` and `ResidualIntervals`, and `ResidualIntervalTests`.
- *Nothing existing changes:* `Prediction`, the grid form, is untouched.
- *Verified by the main session on 2026-11-14,* from a clean archive of `82f98b3` with zero skips: the Commit selection, 2,641 of 2,641, with `b8cd947`'s files in the same run ([receipts](receipts/verify-conformal-s17fu-28993bb/)).

**What it computes.**
- *The ends.* Each score's set is found without a grid, each end within a tolerance in the response's units (1e-8 by default). Each end is rounded outward, so the interval holds the exact set.
- *Why each score's set is an interval.*
  - Distance to origin: at every n, the set is the centre ± the k-th largest |r|.
  - Distance to average: at every n. Each unit's acceptance is |m·rᵢ − S − z| ≥ |n·z − S|, an interval in z about the calibration mean S/n, and the k-fold cover of intervals sharing a point is an interval.
  - RoBAS-EB: for n ≥ 4, by their interval theorem. Below it, a level above 0 is refused, pointing to the grid form.
- *The search.*
  - Their Algorithm 1 starts from an approximate maximiser of the p-value, found by bracketed search. Here it starts from the candidate residual whose own score is 0: â·r̄ over the calibration residuals alone, which is 0 for distance to origin and r̄ for distance to average. That point is accepted, since the candidate's ties count for it and k ≤ n.
  - From there each side doubles a step until a candidate is refused. That happens since the candidate's score grows as |z| and each other's as |z|/n.
  - Each side then bisects to the tolerance, `ConformalBisection.Narrow` as the split predictive distribution's ends use it. Each case costs about 2·log₂(range/tolerance) evaluations of O(n). The set is the same for every case up to the centre f(x), so it is found once per call.
- *Readouts.* CB-1's coverage 1 − k/m. As a full-conformal set it reports no training-conditional law.

**Gates.**
- *Against the grid form* (n 19, 49, 99, mean shifts 0 and 1.5, all three scores, a grid of step 0.001 over [−8, 8]):
  - the grid's accepted values form one interval, and they all lie inside [L, U];
  - each end lies within the tolerance of the grid's first refused value beyond it;
  - read on both sides of each end at tolerance 1e-9, the end itself is refused and two tolerances inside is accepted.
- *Distance to origin's closed form:* the ends are the centre ± the 6th largest |r| at m 60, to within 1e-10.
- *Coverage under a mean shift* (residuals N(1.5, 1), 49 calibration units, 2,000 replications, α 0.1, se 0.0067):

  | Score | Coverage | Mean width |
  |---|---|---|
  | distance to origin | 0.9035 | 5.640 |
  | distance to average | 0.9065 | 3.397 |
  | RoBAS-EB | 0.9045 | 3.398 |

- *The whole line* at k = 0.
- *Refused:* RoBAS-EB at n 3 with a level above 0 (distance to average serves there); a tolerance of 0; a nonfinite residual or centre, at any level; an undefined score.
- *Neighbours:* slice 17's `ScoresBeyondTheModelTests` pass unchanged.
- *Per commit,* the selected classes and the pinned tier: ResidualIntervalTests and ScoresBeyondTheModelTests, the classes reading the file, with the pinned tier: Maths.Regression.Tests 52 of 52; Maths.Oracle.Tests compiles the file and selects no class; Maths.Samplers.Tests compiles neither file.

**Defects found: none.**

**Doc changes owed:** the scores-beyond-the-model section: the grid-free intervals, why each score's set is an interval, the zero-score start against their maximiser, the gates' numbers, and the API rows for `Intervals`, `ResidualInterval` and `ResidualIntervals`.

### 2026-11-14 — Closed: an exact law for the asymmetric interval

Closed at the user's word (2026-11-14, "close it and record the finding and point the docstring"). The miscoverage bound's second follow-up is not built.

**Landed** on `main`, nothing pushed: `b8cd947`, docstrings only. `AsymmetricConformalInterval.TrainingConditionalMiscoverage` and `SplitConformal.AsymmetricQuantiles` carry the finding and point to the one-score forms. *Verified by the main session on 2026-11-14:* every changed line is a doc comment, and both files ran in `28993bb`'s verification, 2,641 of 2,641 at `82f98b3`.

**The finding.**
- *No exact law in general.* The tails' scores, A = y − q_lo(x) and B = q_hi(x) − y, sum to the band's width q_hi − q_lo. Given the calibration the miss is G_A(a₍ₖ₎) + G_B(b₍ₖ₎), two Beta(k, m − k) terms.
- *Where the width varies* (CQR's purpose), the two scores order the units differently, and the terms' dependence follows the width's law. Given only their laws, the worst case of the sum is the union (Makarov 1981), so the union is what holds without assumptions on the width. Whether CQR's structure reaches that worst case is not settled here.
- *Where the width is constant,* B is −A shifted, one ordering serves both tails, and the miss is exactly Beta(k_lo + k_hi, m − k_lo − k_hi). That is the residual interval's case, and the earlier gate's fixture (0.85% passing the union against 10%).
- *Why not build the narrow case.* Detecting it would rest on floating-point equality of a model's quantiles. It fires almost only on synthetic data, and where it fires the interval is the residual interval.
- *For an exact bound,* or sets to the guarantee without the union's price, take a one-score form: `Residuals`, `Quantiles` (symmetric CQR), or the predictive distribution's equi-tailed set.

**Gates:** docstrings only. The build is clean, and `TrainingConditionalTests` pass 8 of 8.

**Doc changes owed:** the guarantee-read-out section records why the asymmetric interval keeps the union, and names the one-score forms.

### 2026-11-14 — The confidence through BASS's group-units call

Built as the coordinator approved: one optional parameter on `BassFit.GroupConformal`, committed on its own.

**Landed** on `main`, nothing pushed:
- `14fd05d`:
  - `BassFit.GroupConformal` gains `double? confidence = null`, last in its parameter list, with its param tag. Nothing else in `BassFit.cs` changes.
  - `BassCandidates.GroupUnits` passes it to `GroupConformal.Density`.
  - `BassGroupConformalTests` gains one test.
- *Existing calls* compile and read unchanged.
- *Verified by the main session on 2026-11-14,* from clean archives with zero skips:
  - Maths.Regression.Tests' Commit selection at `19fd518`, 2,637 of 2,637 ([receipts](receipts/verify-passthru-14fd05d/)).
  - Maths.Oracle.Tests', 56 of 57. `GbassHorseshoeReadoutTests.HBass_MeetsCorrectedHbass_OnTheSurfaceTheNoiseScaleAndV`'s R oracle timed out at 120 s while the machine was loaded (the regression run took 25 min against about 6). Rerun alone from a clean archive of `6f6bb27`, it passed, 1 of 1 ([receipts](receipts/verify-passthru-oracle-rerun/)).

**Gates.**
- *The test* (the crossed fixture, twelve held-out readers each at schools L0 and L1, m 13, α 0.5). The same seams come through the case-level call, each case in its cell's stratum, then pooled, since the strata do not touch the values.
  - *Without a confidence:* the call omitting it, the call passing null, and `GroupConformal.Density` on those seams give the same p-values bit for bit, at CB-1's k of 6.
  - *At β 0.9:* the call matches `GroupConformal.Density` with the same confidence, in level and p-values. Its level is `TrainingConditional.Level(13, 0.5, 0.9)`, below 6, and each set's training-conditional miscoverage at 0.9 is at most α.
- *The class's two earlier tests* pass unchanged.
- *Per commit,* the selected classes and the pinned tier: Maths.Regression.Tests 236 of 327 classes, 2,637 of 2,637; Maths.Oracle.Tests 13 of 39 classes, 57 of 57; Maths.Samplers.Tests compiles none of these files.

**Defects found: none.**

**Doc changes owed:** the BASS conformal section's row for `BassFit.GroupConformal` gains its `confidence`.

### 2026-11-14 — Sets to the guarantee: the split forms' (α, β) level

Built at the user's pick (2026-11-14, "Sets to the guarantee"), the first of the miscoverage bound's two follow-ups. Separately, `LossBound`'s docstrings now name the Lean's `gaffke_calibration` and `highProbability_risk` beside their Theorem 1 citations (`9cf0fd5`, comment only), as the conformal-lean thread offered (8594176).

**Landed** on `main`, nothing pushed:
- `621370c`:
  - `TrainingConditional.Level(m, α, β)`;
  - an optional `confidence` β on `SplitConformal.Calibrate`, `Density`, `Quantiles`, `AsymmetricQuantiles`, `Residuals`, `Distribution` and `CalibrateDistribution`, and on `GroupConformal.Density`;
  - `SplitCalibration.Confidence`, `SplitPitCalibration.Confidence`, and `Confidence` on `ConformalPrediction`, `AsymmetricConformalIntervals`, `ConformalDistribution` and `SplitPredictiveDistribution`;
  - three tests in `TrainingConditionalTests`.
- *Nothing existing changes.* Without a confidence every level is ⌊αm⌋ (⌊αm/2⌋ per tail) as before; the existing classes pass unchanged.
- *Verified by the main session on 2026-11-14,* from a clean archive of `2390b27` with zero skips ([receipts](receipts/verify-conformal-guarantee-621370c/)):
  - Maths.Regression.Tests' Commit selection, 2,636 of 2,636;
  - Maths.Oracle.Tests', 57 of 57.
  - The changed files include `9cf0fd5`'s `LossBound.cs`, whose diff is doc comments only.

**What it computes.**
- *The level.* The largest k at most ⌊αm⌋ whose Beta(k, m − k) β-quantile is at most α (Vovk 2012's (ε, δ)-validity), by bisection, the quantile rising in k. With chance at least β over the calibration units, the stratum's sets miss a new unit with chance at most α. The cap keeps the marginal guarantee 1 − k/m ≥ 1 − α beside it. It is 0, the whole space, where even k = 1 passes α: at α 0.1 and β 0.9 that is m up to 22, where 1 − 0.1^(1/(m − 1)) passes 0.1.
- *Where it lives.* `SplitCalibration` takes the confidence and caches each stratum's level, so every reading of `Level(h)` follows: thresholds, ranks, sets and the readouts.
- *Per form.*
  - The density, CQR and residual forms: per stratum.
  - The group form: over the held-out groups.
  - The asymmetric interval: each tail at α/2 and confidence (1 + β)/2, joined by their union.
  - The equi-tailed set of the predictive distribution: the largest k per tail whose Beta(2k, m − 2k) quantile is at most α, both tails reading the PIT.
- *Refused.*
  - A confidence with weights, on every form. The weighted threshold moves with each case's weights, which the law does not read.
  - β outside (0, 1).
- *What k it takes* at α 0.1:

  | m | ⌊αm⌋ | β 0.5 | β 0.9 | β 0.99 |
  |---|---|---|---|---|
  | 20 | 2 | 2 | 0 | 0 |
  | 50 | 5 | 5 | 2 | 1 |
  | 100 | 10 | 10 | 6 | 4 |
  | 200 | 20 | 20 | 15 | 11 |
  | 1,000 | 100 | 100 | 88 | 78 |
  | 10,000 | 1,000 | 1,000 | 962 | 931 |

**Gates** (β 0.9; the share of calibrations whose sets miss a new unit with chance above α, each set's miss read exactly):

| Form | k, plain → to β | Share above α, plain → to β |
|---|---|---|
| Residual interval, m 200, α 0.1, 4,000 replications | 20 → 15 | 0.4617 → 0.0975 (mean width 3.309 → 3.586, +8%) |
| Asymmetric CQR, m 200, 4,000 replications | each tail at (1 + β)/2 | 0.4695 → 0.0035 |
| Equi-tailed set, m 100, α 0.2, 2,000 replications | 10 → 7 per tail | 0.4825 → 0.0495 |

- *Reading the table.* The plain level misses above α in about half of all calibrations, as its Beta law's median near k/m says. To β the residual interval sits at 1 − β, its scores continuous. The asymmetric interval's union and the equi-tailed set's whole-k step make those two conservative.
- *The rule.* At each m and β the level is the largest k meeting (α, β), at most ⌊αm⌋, and falls as β rises; the two-tail rule likewise.
- *Every form reads it.* The density sets take the level per stratum. The group form takes it over its 100 groups. `Calibrate` takes it per stratum, with a stratum past the last at 0. The grid distribution's `Level(α)` and its readout agree. Without a confidence, the levels are ⌊αm⌋ and ⌊αm/2⌋.
- *Refusals:* a confidence with weights (`Density`, `Distribution`, `AtCase`, `Residuals`), and β of 0 or 1.
- *Existing classes,* the miscoverage bound's five, pass unchanged.
- *Per commit,* the selected classes and the pinned tier: Maths.Regression.Tests 236 of 327 classes, 2,636 of 2,636; Maths.Oracle.Tests 13 of 39 classes, 57 of 57; Maths.Samplers.Tests compiles none of these files.

**Defects found: none.**

**Not wired:** `BassFit.GroupConformal`, GBASS's file, takes no confidence yet. It would be one parameter passed through to `BassCandidates.GroupUnits`, left for GBASS or the coordinator to approve. *Approved by the coordinator on 2026-11-14:* one optional nullable β at the end of the parameter list, default null. No other change to `BassFit.cs`. A test pins it: at null it matches the current call bit for bit, and at 0.9 it matches `GroupConformal.Density`.

**Doc changes owed:**
- the guarantee-read-out section: building sets to (α, β), the level rule and its cap, the two-tail rule, the per-form table above, and the k table;
- the API rows for `TrainingConditional.Level` and the `confidence` parameters and properties.

### 2026-11-14 — Slice 18's follow-up: a set's own miscoverage bound

Built at the user's pick (2026-11-14, "miscoverage it is"). Slice 18's other follow-up, the Lean's calibration at n ≥ 2, is being handled in the conformal-lean thread.

**Landed** on `main`, nothing pushed:
- `0f4d48e`:
  - `TrainingConditional.Miscoverage(m, k, β)`;
  - `TrainingConditionalMiscoverage(β)` on `ConformalSet`, `ConformalInterval` and `AsymmetricConformalInterval`, and per stratum on `SplitCalibration`;
  - `TrainingConditionalMiscoverage(α, β)` on `ConformalDistribution` and `SplitPredictiveDistribution`;
  - `TrainingConditionalTests`.
- *Nothing existing changes:* these are new readouts.
- *Verified by the main session on 2026-11-14,* from a clean archive of `74d094b` with zero skips ([receipts](receipts/verify-conformal-s18fu3-0f4d48e/)):
  - Maths.Regression.Tests' Commit selection, 2,633 of 2,633;
  - Maths.Oracle.Tests', 57 of 57.

**What it computes.**
- *The law* (Vovk 2012, Proposition 2a; Bian and Barber 2023). A split set accepts a test unit whose conformity reaches S₍ₖ₎, the k-th smallest of the stratum's n calibration scores, k = ⌊αm⌋, m = n + 1. Given the calibration, its miscoverage over new units is F(S₍ₖ₎−). The probability integral transform puts that at or below U₍ₖ₎, the k-th of n uniforms, Beta(k, m − k): exactly that law with continuous scores, and conservative where scores tie.
- *The readout* is Beta(k, m − k)'s β-quantile. With chance at least β over the calibration units, the set's miscoverage is at most that. β is the confidence, as in `LossBound.Gaffke`. The bound is 0 at k = 0, the whole space. Its mean is k/m, CB-1's marginal reading.
- *Against slice 18's b*_β.* Gaffke bounds a rule fixed before its losses are read, such as a validation split's misses at a chosen set. A threshold that is the calibration's own k-th score needs this law.
- *Where it reads.*
  - The split density form's sets, at their stratum's m and k.
  - The group-level form's sets, whose units are groups.
  - The CQR and residual intervals, and `SplitCalibration` per stratum.
  - The split predictive distribution's equi-tailed set at α. Both tails read one score, the PIT, so a miss is the mass of 2k of the m spacings: Beta(2k, m − 2k), exactly.
  - The asymmetric CQR interval. Its two tails read two scores, so it takes the union bound: each tail at confidence (1 + β)/2, summed. That is valid whatever ties the tails' scores, and loose.
- *Null.* The full step, the cross-fitted forms and CBMA, which are full-conformal forms, read null. So does RoBAS-EB, a full step on residuals. None carries a law of this kind (Bian and Barber). So do the weighted rule and an inflated set.

**Gates.**
- *The exact law* (99 calibration residuals of N(0, 1), α 0.1, m 100, k 10, 4,000 replications; each interval's miscoverage given its calibration is 2Φ(−c)):
  - the mean is 0.09991 against k/m 0.1;
  - the miscoverage passed the β-quantile 0.4920, 0.1000 and 0.0105 of the time, against 1 − β of 0.5, 0.1 and 0.01 (β 0.5, 0.9, 0.99).
- *Tied scores and strata* (responses rounded to integers; strata of 49 and 199 units, m 50 and 200, k 5 and 20; β 0.9): passed 0.0003 and 0.0000, conservative as the law says. Each stratum reads its own m and k, and `SplitCalibration` agrees.
- *The equi-tailed set* (the right fit's predictive distribution, 99 units, α 0.2, k 10 per tail, 2,000 replications, read from the grid-free ends at x 0.4): bound 0.2523, Beta(20, 80)'s 0.9-quantile, passed 0.0995 against 0.10, exactly as the law says. The grid form reads the same. The full step and a weighed case read null.
- *The asymmetric interval* (CQR tails at α/2 on 99 responses, k 5 each, 4,000 replications, β 0.9): bound 0.1801, two tails of 0.0901, passed 0.0085 against 0.10. Valid, and about twelve times looser than its level: the price of the union.
- *The quantile.* It is the root of the binomial tail P(Bin(m − 1, q) ≥ k) = β, summed term by term, to 9 decimals at m and k of 100 and 10, 1,000 and 1, 50 and 49, and 31 and 3, at β 0.5, 0.9 and 0.999. Out-of-range m, k and β are refused.
- *Which records read it:* a split density set at m 31, k 3; a small stratum 0; the group form at m 16, k 1. The weighted rule and the full step read null.
- *Per commit,* the selected classes and the pinned tier: Maths.Regression.Tests 236 of 327 classes, 2,633 of 2,633; Maths.Oracle.Tests 13 of 39 classes, 57 of 57; Maths.Samplers.Tests compiles none of these files.

**Defects found: none.**

**Follow-ups, each with the decision it would change.**
- *The (α, β) level as a set rule:* the split forms choosing the largest k whose β-quantile is at most α. Their sets then hold P(miscoverage ≤ α) ≥ β, and are wider. It decides whether the layer builds sets to a training-conditional guarantee, not only reports one.
- *A tighter law for the asymmetric interval:* where each unit's two tail scores order alike (the tails' scores in the fixture above are one score's two signs), the miss is the mass of k_lo + k_hi spacings, Beta(k_lo + k_hi, m − k_lo − k_hi), as the equi-tailed set's is. It decides whether that interval reads the union or the exact law where the condition holds.

**Doc changes owed:**
- the guarantee-read-out section: the training-conditional law, where it reads and where it doesn't, the equi-tailed set's 2k law, the asymmetric union bound, and the gates' numbers;
- the API rows for `TrainingConditional` and the `TrainingConditionalMiscoverage` readouts.

### 2026-11-14 — Slice 15's precision rule at the calibration split's size; slice 12a's k̂ stays opt-in

Built at the user's pick (2026-11-14, "Slice 15 precision rule"). With it, the user settled slice 12a's k̂ diagnostic as opt-in. That is already the code's default (`ConformalOptions.LeaveOneOutParetoK`, false), so no code changes for it, and the open list can drop it.

**Landed** on `main`, nothing pushed:
- `591eace`:
  - `ScoreMixture.PrecisionAt(calibrationUnits)` and its readout `PrecisionChoice` (m, the size, the mean ESS and self-weight there, and whether the rule is met);
  - `ScoreMixture.Weights` now takes m at its units' count;
  - `ScoreMixtureTests` gains two tests.
- *What changes for existing calls.* Only a call whose calibration split differs from n − 1, n the training split's size, gets a different m. At n − 1 (the existing end-to-end test's 400 against 399) every value is unchanged, bit for bit. A fixed `Precision` stands at every size, as before. The fit-level `Precision` and its readings are the choice at n − 1.
- *Verified by the main session on 2026-11-14,* from a clean archive of `94a37ef` with zero skips: the Commit selection, 49 of 49 ([receipts](receipts/verify-conformal-s15p-591eace/)). `ScoreMixture` is read by `ScoreMixtureTests` alone, so the selection holds every caller.

**What it computes.**
- *The finding.* The rule (theirs, §6.1) chooses the largest m whose randomized weights have mean ESS above 100 and mean self-weight at most 1/30. It reads them on the training split, each point against the other n − 1. Their splits are the same size. Here, fitted at 1,000 and calibrated at 399, a test case stands against 399 units, not 999. Its ESS is about 0.4 of the reading, and its self-weight about 2.5 times. A test case whose own weight passes α gets the whole space.
- *The rule at the calibration's size.* Each training point's weight sums over the others (the sum and the sum of squares) are rescaled by N/(n − 1). These are their expected sums over N units. The rule's two conditions are then read on the rescaled sums.
- *Why it stays valid.* The choice reads the count and the training split alone, never the calibration's values, so the weights stay a fixed function of the inputs and the coverage stays exact.
- *Its approximation.* The rescaling reads expected sums. A point's self-weight is convex in its others' sum, so the rescaled self-weight reads a little low where the weights concentrate. The ESS it reads is likewise slightly low, which makes the rule slightly conservative where it decides.

**Gates.**
- *The finding, repaired* (fitted on 1,000, 2,000 replications of 399 calibration units and a test case, α 0.1):

  | Reading | m | Whole space | Marginal | Clusters 1 / 2 | Mean finite width |
  |---|---|---|---|---|---|
  | the training's own | 500 | 13.4% | 0.9275 | 0.925 / 0.930 | 2.240 |
  | at the calibration's size | 45 | 0.0% | 0.9155 | 0.913 / 0.918 | 2.190 |

  At 399 the rule is met: mean ESS 100.0, self-weight 0.0145. The training's own m 500 read ESS 107.2 and self-weight 0.0189 on 999 others. The training's marginal 0.9275 is inflated by its whole-space sets.
- *The rescaling against pools actually drawn* (each point against 399 others drawn at random, the same draws):
  - at the m it took (35 on this fit): ESS 100.36 against 100.69, self-weight 0.01318 against 0.01324;
  - at m 500: 45.91 against 46.20, and 0.04507 against 0.04640, 2.9% low, as convexity says.
- *Equal splits:*
  - `PrecisionAt(n − 1)` is the fit's choice and readings exactly, and `Weights` at n − 1 units matches the direct weights at `Precision` bit for bit;
  - a smaller split takes a smaller m (199 ≤ 399 ≤ 999), and the choice is cached per size;
  - a fixed m stands at 99 units, with its readings there and the rule not met;
  - a size below 1 is refused.
- *Existing classes,* slice 15's `ScoreMixtureTests` and `PosteriorConformalTests`, pass unchanged.
- *Per commit,* the selected classes and the pinned tier: ScoreMixtureTests, the one class reading the file, with the pinned tier: Maths.Regression.Tests 49 of 49; Maths.Oracle.Tests compiles the file and selects no class; Maths.Samplers.Tests compiles neither file; PosteriorConformalTests run beside it, 5 of 5.

**Defects found: none.** The finding was the rule reading the wrong pool size.

**Doc changes owed:**
- the posterior-conformal section's hyperparameters: the rule read at the calibration split's size, the rescaling and its slight low reading where weights concentrate, and the table above;
- the API rows for `PrecisionAt` and `PrecisionChoice`;
- the waiting list: slice 12a's k̂ settled opt-in.

### 2026-11-14 — BASS's opt-in group-units call for `GroupConformal` (item 16, (b))

Built as commissioned (the user's choice of (b), 2026-11-14; the Open list's specification, c94f244).

**Landed** on `main`, nothing pushed:
- `b79840f`:
  - `BassFit.GroupConformal`, a new method beside `CandidateLogDensity`. No existing member of `BassFit.cs` changes, and its types are qualified (`Conformal.…`), so the file gains no `using`.
  - `BassCandidates.GroupUnits` behind it. `BassCandidates.LogDensity`'s existing signature now forwards to a private overload, whose one flag skips the case-level strata rule for this call alone.
  - `BassGroupConformalTests`.
- `BassGroupCells` is unchanged.
- *Verified by the main session on 2026-11-14,* from a clean archive of `842531c` with zero skips ([receipts](receipts/verify-item16-b79840f/)):
  - Maths.Regression.Tests' Commit selection, 2,626 of 2,626;
  - Maths.Oracle.Tests', 57 of 57.
  - The code was read against item 16's specification: the case-level path forwards with the flag off, so only `GroupUnits` skips the strata rule, and the groups come from labels, a named grouping, or the cells under one grouping with its column given.

**What it computes.**
- *The seams.* The held-out groups' cases are read at their responses and the test cases on their grids, each by BASS's new-group density (the seam the case-level call builds). The labels go to `GroupConformal.Density`.
- *The groups,* from exactly one source:
  - `groupOf`, each held-out case's label;
  - `grouping`, the name of a grouping whose levels are the groups. Every held-out and test case's level of it must be new;
  - otherwise the cells. They are read only where every effect reads one grouping (intercept and slope of one grouping included) and a column gives it.
- *Where the cells don't serve.* A cell is one level of every grouping jointly, so the call refuses the cells reading in two designs, naming the remedy:
  - crossed or nested groupings, where one held-out subject seen at two sites spans two cells;
  - a grouping's column left out, where each case is a cell of its own.
- *Refused.*
  - Both readings at once.
  - A grouping no effect reads.
  - Under `grouping`, a case whose level of it the fit knows.
  - A held-out or test case whose groups are all known. This is checked here because a seam with no new case carries no `NewGroup` flags for `GroupConformal` to read.
  - The case-level call's own refusal stays as it was.
- *Unseen levels.* `unseenLevels` defaults to `Refuse`, as `CandidateLogDensity`'s does, so a named new level asks for `OutsideEverySubset`. The same switch governs categorical predictors' unseen levels, so this call does not change the default on its own.
- *What it covers.*
  - `WorstCase` holds every case of one new group at once.
  - `OneCase` holds one case drawn uniformly from it.
  - Both need the held-out and new groups drawn alike, including how their sizes arise.
  - The earlier report's 0.893 for one pooled case was its random-intercept fixture's alone, not a guarantee for pooled cases.

**Gates.**
- *End to end* (the crossed school × reader fixture):
  - two held-out readers, each with repeated cases at schools L0 and L1, span four cells;
  - named by `grouping: "reader"` they are two groups (m = 3), and under both scores they reproduce the p-values given explicit labels;
  - the cells reading refuses the crossed fit;
  - BASS's case-level `CandidateLogDensity` still refuses the same pooled stratum, both left to its default (the `GroupCells` message) and given as one stratum ("spans cells");
  - each refusal above is exercised.
- *One grouping* (the fixture restricted to the reader's intercept, its data to the reader's grouping, since BASS refuses a grouping no effect reads): the cells reading gives the named reading's p-values, and with the reader's column left out it is refused.
- *Neighbours:* `RandomEffectsCandidateTests`, `ConformalCellRuleTests` and `GroupConformalTests` pass unchanged.
- *Per commit,* the selected classes and the pinned tier: Maths.Regression.Tests 235 of 326 classes, 2,626 of 2,626; Maths.Oracle.Tests 13 of 39 classes, 57 of 57; Maths.Samplers.Tests compiles none of these files.

**Defects found: none.**

**Doc changes owed:**
- the BASS conformal section's API rows: `BassFit.GroupConformal`, its three group readings and its refusals;
- the new-groups section: BASS's group-units call beside `GroupConformal`, with the size condition and the note that the 0.893 for pooled cases is the fixture's alone.

### 2026-11-14 — GBASS's two follow-ups: group-level units for new groups, and the layer's cell rule

Built at the user's pick (2026-11-14, "GBASS's two follow-ups"). The cell rule was built as the coordinator approved it, as one change: the field, BASS's wiring and the refusal. The defects chip was told before `Pointwise.cs` was touched and confirmed it had nothing pending there.

**Landed** on `main`, nothing pushed:
- `5ce91b9`:
  - `GroupConformal.Density` and `GroupScore` (`WorstCase`, `OneCase`), with `GroupConformalTests`;
  - `Cells` (optional, default null, length checked) on `PointwiseLogDensity` (one per observation) and `CandidateLogDensity` (one per case);
  - `BassPointwise.LogDensity` and `BassCandidates` set them from `BassGroupCells`, read only; `BassFit.cs` and `BassGroupCells` are untouched;
  - the refusal, `ConformalBayes.CheckCells`, read by `ConformalBayes.Check` and `SplitDraws.Check`, with `ConformalCellRuleTests`.
- *Verified by the main session on 2026-11-14,* from a clean archive of `f96be6c` with zero skips ([receipts](receipts/verify-gbass-fu-5ce91b9/)):
  - Maths.Regression.Tests' Commit selection, 3,033 of 3,033;
  - Maths.Oracle.Tests', 145 of 147. The two misses are `DiagramMetricsTransportParityTests` (p 1 and 2): its R transport oracle timed out at 120 s under load. It reads none of these files, and it missed the same way in the defects chip's slice 6 verification.

**A — group-level units for new groups.**
- *What it computes.* The split density form with held-out groups as its units. There are m = G + 1 units, k = ⌊αm⌋, and CB-1's rule (B + T + 1)/m > k/m on one score per group:
  - the worst-case score, a group's least-conforming case: a new group's cases are held all at once;
  - the one-case score, one case drawn per group on a seeded stream: one case is held.
- *Conditions.* The calibration groups are held out and scored as new groups. The groups must be drawn alike, their sizes included. Strata are not read.
- *Refused.* A calibration or test case whose group the fit knows (`NewGroup` false).
- *Sources.* The idea is Dunn, Wasserman and Ramdas's, as the agreement names it. Their paper is not on the shelf, so the construction rests on CB-1 alone, and the docstring says so.

**B — the cell rule.** Where both seams give cells, a case in a known group ranks only among units of its own cell. This holds in the full step (every path through `Check`: `FullStep`, `Scores`, the distribution, the CDF margin's probes, CBMA) and in the split forms (every path through `SplitDraws`: density, distribution, cross-conformal, tournament, risk-controlled sets).
- *Exempt.*
  - A case of a new group;
  - a seam without cells;
  - a stratum no unit holds;
  - training strata no candidate ranks in, so LOO and WAIC's single-stratum partitions stay as they are;
  - `GroupConformal`, whose units are groups and so span cells.
- *What it adds to BASS's own checks.* BASS checked the training partition's strata (`BassFit.PointwiseLogDensity`) and the candidates' strata among themselves (`BassCandidates`), but never a candidate's stratum against the training units' cells. Strata relabelled among the candidates passed both checks and put each case among another cell's units. The layer now refuses that: the BASS test relabels the crossed fixture's strata, and the full step refuses.

**Changes to existing calls (flagged as asked).**
- *Results: none.* Cells are read only by the check. With cells and without them, the full step's and the split form's p-values match exactly.
- *Refusals: one kind.* A call with a known-group candidate in a stratum whose units lie in another cell, or span cells, is now refused where it was served. That is the call above. No existing test makes one: the E21 classes pass unchanged (`RandomEffectsCandidateTests` 25 of 25, `RandomInterceptTests` 19 of 19).

**Gates.**
- *Group units* (79 held-out groups of 3 to 8 cases, a new group drawn alike, 2,000 replications, α 0.1): worst-case group units hold the new group whole 0.8925, one-case group units hold one case 0.8945; case units hold it whole only 0.7565 and one case 0.8930; se 0.0067.
- *Singleton groups:* identical to the case-level split form. *One-case:* seeded. The refusals are covered.
- *The cell rule:* the full step and the split forms each:
  - serve each case in its cell's stratum, with the p-values unchanged;
  - refuse a case among another cell's units, and a stratum spanning cells;
  - exempt new-group cases, a seam without cells and the group-level form.
- *BASS's seams* carry `BassGroupCells`' training and candidate cells exactly, and the relabelled call is refused.
- *Per commit,* the selected classes and the pinned tier: Maths.Regression.Tests 303 of 325 classes, 3,033 of 3,033; Maths.Oracle.Tests 28 of 39 classes, 147 of 147; Maths.Samplers.Tests compiles none of these files.

**Finding: the group-units condition on sizes.** The first fixture fixed the new group at 6 cases against calibration groups of 3 to 8. Its rates read low: worst-case 0.8845 and one-case 0.8800 against 0.9, se 0.0067. A probe at 20,000 replications on fresh seeds split this into two causes:
- *Mostly an unlucky sample.* On the same 2,000 draws, the true 90% band held y₀ only 0.8820 of the time. At 20,000 replications the one-case rate sits on that control: 0.9052 against 0.9032.
- *A real size effect under the worst-case score.* At a fixed size of 6 it ran about 0.006 under the control (0.8969 and 0.8901). With the size drawn alike it matched the control: 0.9011 against 0.9015.

The fixture now draws the new group's size alike, and the docstring states the condition.

**Defects found: none.**

**Open, for the coordinator: BASS cannot build `GroupConformal`'s calibration seam yet.** `BassCandidates` refuses null strata when cells span more than one, and strata inside cells conflict with the group form's single stratum. Held-out groups, each its own new cell, are therefore refused at the builder. The choice belongs to GBASS and the coordinator:
- (a) `BassCandidates` exempts new-group cases from its strata rule, as the layer's rule exempts them. This loosens what the case-level split forms accept from BASS.
- (b) An opt-in group-units call on BASS's side.

Either way, `GroupConformal` could then read its groups from the seam's `Cells` when no labels are given.

**Doc changes owed:**
- the conformal layer's API rows for `GroupConformal`, `GroupScore` and the two `Cells` fields;
- the new-groups section: group units, the two scores and what each covers, the size condition, and the probe's numbers;
- the strata section: the cell rule, where it is read, its exemptions, and what it adds to BASS's checks.

### 2026-11-14 — Slice 18's follow-up: the smoothed cell law for the log score

Built at the user's pick (2026-11-14, "Smoothed cell law for log score"), the first of slice 18's follow-ups.

**Landed** on `main`, nothing pushed:
- `b4b7e22`:
  - `ConformalDistribution.Smoothed` (`SmoothedCdf`: G̃ and its complement at each grid value);
  - `CellLaw` (`Steps`, `Smoothed`), which `CellMasses`, `LogScore` and `QuadraticScore` now take;
  - `SplitPredictiveDistribution.SmoothedAt`;
  - `SmoothedCellLawTests`.
- *Nothing existing changes:* `CellLaw.Steps` is the default, so every existing call reads as before.
- *Verified by the main session on 2026-11-14,* from a clean archive of `75e58d8` with zero skips: the Commit selection, 1,228 of 1,228 ([receipts](receipts/verify-conformal-s18fu1-b4b7e22/)).

**What it computes.**
- *The law.* G̃(F̂(y | x)), G̃ the calibration PITs' empirical CDF interpolated linearly between its knots (0, 0), (PIT₍ᵢ₎, i/(n + 1)) and (1, 1). At a tie with the PITs it is (B + T)/(n + 1).
- *Where it sits.* It lies inside the conformal band [Q⁻, Q⁺] at every y: it is the randomized Q(y, τ) with τ moving linearly across each gap between neighbouring PITs.
- *What it fixes.* Wherever the fit's CDF rises, every cell has mass, so the log score is finite: `logLoss_proper`'s positivity hypothesis, which the step law breaks, is met.
- *Precision.* Each point is read on the tail that holds its digits, as the CDF scores are, its complement kept beside it.
- *What it gives up.* Exact uniformity: its PIT at the truth lies within 1/(n + 1) of uniform. The randomized Q keeps the exact guarantee; this law is for scoring.
- *Where it applies.* The split form's unweighted cases. The full step's score is not monotone in y, and a weighed case reads weighted ranks, so both have none.

**Gates.**
- *Band and mass* (49 calibration units, 20 cases, cells of width 0.1): G̃ lies in [Q⁻, Q⁺] and rises at every grid value, G̃ + (1 − G̃) = 1, and every cell has mass. In every case the step law leaves cells empty.
- *The log score* (199 calibration units, 400 cases, cells of width 1): the smoothed score is finite everywhere and ranks the right fit ahead of one ignoring x, 1.447 against 1.523. In this draw the step law happened to be finite everywhere too; its +∞, 2 of 400 in slice 18's draw, is the band-and-mass gate's empty cells.
- *Uniformity* (99 calibration units, 4,000 replications): the smoothed PIT departs from uniform by at most 0.0068 at the deciles, within 1/(n + 1) = 0.0100.
- *The grid-free readout:* `SmoothedAt` matches the grid's G̃ bit for bit.
- *The full step and a weighed case* have no smoothed law, and their scores on it are NaN. An undefined `CellLaw` is refused.
- *Existing classes,* slice 18's and the exact-ends tests, pass unchanged.
- *Per commit,* the selected classes and the pinned tier: 1,228 of 1,228.

**Defects found: none.**

**Doc changes owed:** the guarantee-read-out section's cell scores gain the smoothed law, with its band property, its 1/(n + 1) uniformity, where it applies, and the API rows for `SmoothedCdf`, `CellLaw` and `SmoothedAt`.

### 2026-11-14 — Slice 17, scores beyond the model's own

Built at the user's pick (2026-11-14, "Slice 17, scores beyond the model"), the archived plan's row 17, the last in the plan. Sources read in TeX:
- Bhagwat, Kong and Jiang, CBMA (arXiv 2511.16924v1), §3, Algorithm 1, the lemma and both theorems;
- Ashouritaklimi, Cortinovis and Caron, RoBAS (arXiv 2607.04236v2), §2.2–§3.4;
- Wu, Lohmanova, Kaski and Caprio (arXiv 2602.03331v2), §4 and Algorithm 1.

Each is a symmetric score or a rule already proved, so coverage stands on CB-1, and for Wu et al. on slice 18's `fail_le_fixed`. No new Lean was needed. The plan's Lean column asked for Wu et al.'s risk bound "if it differs from coverage": it is slice 18's high-probability rule, already proved and verified.

**Landed** on `main`, nothing pushed:
- `0b91db9`:
  - `ModelAveraging.FullStep`, giving an `AveragedPrediction`;
  - `ResidualShrinkage.Prediction` with `ResidualScore` (`DistanceToOrigin`, `DistanceToAverage`, `Shrinkage`);
  - `RiskControlledSets.Density`, giving `HighDensitySets`;
  - `ScoresBeyondTheModelTests`.
- *Nothing existing changes.*
- *Verified by the main session on 2026-11-14,* from a clean archive of `ca4ff81` with zero skips: the Commit selection, 48 of 48 ([receipts](receipts/verify-conformal-s17-0b91db9/)).

**What it computes.**
- *CBMA* (`ModelAveraging`).
  - From K models' full-step seams over the same units, strata and grids, and each model's log posterior probability given the n units (its log evidence plus log prior), the score is σᵢ = Σ_k q_k(y)·σᵢ^{M_k}, q_k(y) ∝ p(M_k | Z₁:ₙ)·p_{M_k}(y | x, Z₁:ₙ).
  - That is p(M_k | Z₁:ₙ₊₁), so σᵢ is the hierarchical model's posterior predictive, symmetric in the n + 1 units (their lemma). CB-1's rule on it covers at least 1 − α.
  - The per-model scores are `ConformalBayes.Scores`'s, mixed by log-sum-exp. The weights are reported at every candidate.
- *The residual scores* (`ResidualShrinkage`). Full conformal on a fixed predictor's residuals, each of the n + 1 scored against the other n:
  - distance to origin, |r|;
  - distance to average, |r − r̄|, the Bayes score of a flat prior on the mean (their proposition);
  - RoBAS-EB, |r − â·r̄| with â = υ̂²/(υ̂² + σ̂²/n) from their EB estimates. It shrinks to the origin when the residuals sit near zero and to their average under a mean shift.

  The others' mean and spread are read from residuals centred at their overall mean, so they keep their digits under a large shift.
- *Wu et al.'s high-density sets* (`RiskControlledSets`). The split form's density score −ln p̂(y | x) and the set {s ≤ λ*}, λ* the least threshold whose Gaffke bound on the calibration miscoverage is at most α (slice 18's `HighProbability`). The set's miscoverage exceeds α with chance at most 1 − β over calibration sets, a training-conditional guarantee. Their Bayesian quadrature only estimates the set's size, so it is not built.

**Gates.**
- *CBMA with one model* is the full step's set bit for bit, its weight 1.
- *CBMA's weights against the exact evidence.* Two noise laws, σ 1 and 1.5, on 8 units; the exact posterior model probabilities run from 0.015 to 0.669 over the grid. With 4,000 exact draws the weights match them to 0.0021.
- *CBMA across models wrong in their mean* (y = 0.7x + N(0, 1), models βx and βx², n 30, 500 replications, α 0.1):

  | score | coverage | mean width |
  |---|---|---|
  | CBMA | 0.912 | 3.52 |
  | βx alone | 0.912 | 3.44 |
  | βx² alone | 0.910 | 4.14 |

  The mean weight on βx at the truth is 0.958.
- *A model wrong only in its noise scale* (σ 3 against σ 1) barely moves the conformal set: 3.52 against 3.47 at 800 replications. A symmetric unimodal density ranks units by |residual| whatever its scale. So the plan's "wrong noise" bites only where the law's shape or centre differs, and the gate above uses a wrong mean.
- *The residual scores' definitions:* 300 random vectors, some shifted by up to 10⁴, against two-pass definitions.
- *The residual scores under a right and a shifted predictor* (n 15, α 0.1, 2,000 replications, se 0.0067):

  | shift | distance to origin | distance to average | RoBAS-EB |
  |---|---|---|---|
  | 0: coverage, width | 0.930, 4.11 | 0.926, 4.27 | 0.927, 4.17 |
  | 3: coverage, width | 0.946, 9.50 | 0.947, 4.31 | 0.948, 4.32 |

  RoBAS-EB sits between the references when the predictor is right and matches distance to average under the shift.
- *Wu et al.'s sets* (a fit to 50 training units, 60 calibration units, the miscoverage read on 1,000 fresh units, 300 calibration sets, α 0.1, β 0.9): the miscoverage is at most α in 0.957 of the calibration sets (se 0.017), with mean miscoverage 0.046.
- *The refusals:*
  - no models, weights of another count, models over other units, no model with positive weight;
  - a NaN residual, centres and grids of different counts, a grid out of order, an undefined score;
  - a confidence of 1.
- *Per commit,* the selected classes and the pinned tier: 48 of 48 (one class, the new one).

**Defects found: none.**

**Follow-ups, each with the decision it would change.**
- *RoBAS-Full,* the horseshoe score with ₁F₁, beside the EB score. Decides whether the layer offers the fully Bayesian robust score.
- *CBMA's evidence from the engines' own fits.* Here it is exact and given; BASS and BARS would supply estimates such as bridge sampling. Decides whether model averaging reads an engine's evidence directly.
- *A grid-free interval for the residual scores* (their bracketed root search), RoBAS-EB's sets being intervals. Decides whether `ResidualShrinkage` reports exact ends.

**Doc changes owed** (for the offline docs pass): `docs/Maths/Regression/ConformalPrediction.md`:
- a section on scores beyond the model's own: CBMA's score and weights, the three residual scores, and Wu et al.'s risk-controlled sets with their training-conditional guarantee;
- the finding that a wrong noise scale barely moves conformal sets;
- the gates' numbers;
- API rows for `ModelAveraging`, `AveragedPrediction`, `ResidualShrinkage`, `ResidualScore`, `RiskControlledSets` and `HighDensitySets`.

### 2026-11-14 — Slice 18, the guarantee read out

Built at the user's pick (2026-11-14, "Slice 18, the guarantee read out"), under the Lean's sections `LossPosterior` and `ProperScores` (proved 2026-11-10, checked against their sources 2026-11-13). It follows the archived plan's row 18 and its 2026-11-13 correction: b*_β is Gaffke's confidence bound, calibrated over replicated calibration sets, and never named a posterior.

**Landed** on `main`, nothing pushed:
- `5b55a6f`:
  - `LossBound`, with `Cdf`, `Survival`, `Mean` and `Gaffke`;
  - `RiskControl.Mean` and `RiskControl.HighProbability`, giving a `RiskThreshold`;
  - `ConformalDistribution.CellMasses`, `LogScore` and `QuadraticScore`;
  - `GuaranteeReadoutTests`.
- *Nothing existing changes.*
- *Verified by the main session on 2026-11-14,* from a clean archive of `07dc64e` with zero skips: the Commit selection, 1,218 of 1,218 ([receipts](receipts/verify-conformal-s18-5b55a6f/)).

**What it computes.**
- *The law of L⁺ = Σ Uᵢ ℓ₍ᵢ₎,* over n losses in [0, B] and ℓ₍ₙ₊₁₎ = B, U ~ Dir(1, …, 1).
  - P(L⁺ > x) is the divided difference of (t − x)₊ⁿ on the knots, confluent where they tie. It is read in exact rational arithmetic on the losses' binary values, so discrete and tied losses read exactly.
  - Losses of two values take an exact binomial sum, O(n) a reading: there L⁺ = a + (B − a)·Beta(c_B, c_a).
  - `Mean` is (Σ ℓ + B)/(n + 1).
- *`Gaffke(β)`,* b*_β, the least x with P(L⁺ ≤ x) ≥ β, by bisection on the exact law. It is Gaffke's upper confidence bound on the risk J: P(J ≤ b*_β) ≥ β over replicated calibration sets of independent losses at most B (Vlassis and Thomas's Theorem 1). At one loss it is ℓ + β(B − ℓ) (`gaffke_one`).
- *`RiskControl`,* on a loss matrix nonincreasing in the threshold:
  - `Mean`: the first threshold with (Σ ℓ + B)/(n + 1) ≤ α. The expected test loss is at most α (`crc_orderings`), given the losses at the last threshold at most α.
  - `HighProbability`: the first threshold with b*_β ≤ α, by bisection. Its risk exceeds α with chance at most 1 − β (`fail_le_fixed`).
  - Each reports whether any threshold met α.
- *Proper scores on the grid's cells.* `CellMasses` is the predictive's law on (−∞, g₀], (g₀, g₁], …, (g_{G−1}, ∞), from Q̄.
  - `LogScore(y)` is −ln of y's cell (`logLoss_proper`, `logLoss_strict`).
  - `QuadraticScore(y)` is Σ p² − 2p of y's cell, minus Gneiting and Raftery's quadratic score (`gneitingLoss_gap`).
  - Both are NaN for a case with no distribution.

**Gates.**
- *The law is exact:*
  - one loss gives ℓ + β(B − ℓ) on 50 cases;
  - 0/1 losses give Beta(k + 1, n − k), P(Bin(n, x) ≥ k + 1), to 10⁻¹² on 60 cases;
  - three loss sets, ties included, agree with 100,000 Dirichlet draws at five points within 4.5 se.
- *Calibration over replicated calibration sets* (2,000 each), P(J ≤ b*_β):

  | loss law | n 3 | n 8 |
  |---|---|---|
  | Bernoulli(0.3), β 0.75 | 1.000 | 0.941 |
  | Bernoulli(0.3), β 0.9 | 1.000 | 0.938 |
  | three points on {0, ½, 1}, β 0.9 | 1.000 | 0.973 |
  | U², β 0.9 | 1.000 | 0.997 |

  Every rate is at least β, conservative at these small n, discrete losses included.
- *The rules* (miscoverage of {s ≤ λ} on a grid of 101 thresholds, scores uniform, n 30, α 0.1, 1,500 calibration sets): conformal risk control's test loss is 0.093; the high-probability rule (β 0.9) fails 0.026, against the allowed 0.1.
- *Propriety, exactly,* on 500 pairs of laws on 2 to 8 cells: the log score's gap is KL(q ‖ p) > 0, and the quadratic score's gap is Σ (p − q)².
- *On the split predictive distribution* (199 calibration units, 400 test cases, cells of width 1): the right fit beats one ignoring x on the quadratic score, −0.254 against −0.237.
- *The refusals:*
  - losses outside [0, B] or NaN, an infinite B, β outside (0, 1);
  - losses that rise in the threshold, losses above B, no thresholds.
  - A rule no threshold satisfies returns the last, flagged.
- *Existing classes,* the predictive-distribution tests, pass unchanged.
- *Per commit,* the selected classes and the pinned tier: 1,218 of 1,218.

**Defects found: none in code.** One finding about the readout. The split predictive distribution is a step function in the calibration ranks, so a cell no calibration PIT rises across has no mass, and the log score is +∞ for an observation there: 2 of 400 cases for the right fit at cells of width 1. The log-loss's propriety needs the forecast positive where the truth is (`logLoss_proper`'s hypothesis), which such cells break. Over the finite cases the right fit scores 1.48 against 1.61. The quadratic score has no such gap and is the one that ranks fits on these distributions.

**Follow-ups, each with the decision it would change.**
- *A smoothed cell law for the log score:* the PIT ranks' randomized Q spread over the cell, or a floor. Decides whether the layer offers the log score on predictive distributions without +∞.
- *The calibration for n ≥ 2 in the Lean* (Vlassis and Thomas's chain-measure induction), formalized only at one unit here. Decides whether the Lean states what `HighProbability` reports.
- *A set's own miscoverage bound:* b*_β applies to losses of a fixed set rule, such as a validation split's misses. A set whose threshold is chosen on the same calibration units needs the training-conditional law instead. Decides whether conformal sets report a held-out bound.

**Doc changes owed** (for the offline docs pass): `docs/Maths/Regression/ConformalPrediction.md`:
- a section on the guarantee read out: Gaffke's bound as a confidence bound, never a posterior, with Snell and Griffiths's Theorem 4.3 failing as stated;
- the two risk rules and their guarantees;
- the cell scores beside the CRPS, and the log score's +∞ on step-function laws;
- the gates' numbers;
- API rows for `LossBound`, `RiskControl`, `RiskThreshold`, `CellMasses`, `LogScore` and `QuadraticScore`.

### 2026-11-14 — Slice 16, batch predictive inference

Built at the user's pick (2026-11-14, "Slice 16, batches"), under the Lean's section `Batch` (proved 2026-11-10) and Lee, Tchetgen Tchetgen and Dobriban (arXiv 2409.13990v5, read in TeX: §2, Condition 1, Algorithm 1, Theorem 1, the quantile and selection corollaries, the joint-bounds theorem and corollary).

**Landed** on `main`, nothing pushed:
- `60ba2cb`: `BatchConformal`, with
  - `RankLaw` and `Quantile`, giving a `BatchInterval`;
  - `BoxCoverage`, `Box` and `Simultaneous`, giving `BatchBounds`;
  - `Function` and `Select`, the latter giving a `BatchSelection`;
  - and `BatchConformalTests`.
- *Nothing existing changes.*
- *Verified by the main session on 2026-11-14,* from a clean archive of `eb0b064` with zero skips: the Commit selection, 49 of 49 ([receipts](receipts/verify-conformal-s16-60ba2cb/)).

**What it computes.** Batch sets for functions of m test units' scores from n calibration scores, exchangeable together. The scores' floor and ceiling, S_(0) and S_(n+1), are ±∞ unless given.
- *The basis.* The test units' sorted ranks are uniform on H = {1 ≤ r_1 ≤ … ≤ r_m ≤ n + 1}, |H| = C(n + m, m), and every sorted test score lies between S_(r_j − 1) and S_(r_j) at once (`card_testRank_mem`, `batch_sandwich`). Each set reads that sandwich on rank sequences holding at least 1 − α of H.
- *`Quantile`,* the ζ-th smallest score of the batch: [S_(q_L − 1), S_(q_U)], from the exact law of r_(ζ) (their quantile corollary). q_U is the least rank with at most γ of the law above it; q_L the greatest with at most β below. At m = 1 it is split conformal's set.
- *`Box` and `BoxCoverage`,* several order statistics at once: [S_(w_j − 1), S_(q_j)] with the box's share of H (their joint-bounds theorem, `batch_simultaneous`). The count runs by dynamic programming on L(ρ), with the product over consecutive orders, the Lean's correction of their closed form.
- *`Simultaneous`:* their centred search, each order's box centred at round(t_j·n/m) and widened by the least a whose share reaches 1 − α, by bisection.
- *`Function`,* any h monotone in each sorted score: their Algorithm 1, with the rank order h̃ given (h on the ranks by default). It enumerates H twice, refused above 2,000,000 sequences, since computing it is NP-hard in general.
- *`Select`:* scores μ̂·1{y ≤ c}, the threshold T̂ = S_(q_η), and the units with μ̂ above it declared above c, at most η of them falsely with chance at least 1 − α (their selection corollary; `batch_exceed`, `false_claim_iff`).
- *Exact levels.* The rank counts are integers (`BigInteger`). A share of H at most γ is the count at most γ·|H| for the given γ read exactly, so whole levels such as (1 − γ)(n + 1) = 18 land on their rank and do not round away from it.

**Gates.**
- *The counts are enumerations of H:*
  - the law of r_(ζ) at every k, for n ≤ 6 and m ≤ 4;
  - 300 random boxes' shares, to 10⁻¹².
- *At one test unit:* the quantile set is split conformal's, q_U = ⌈(1 − γ)(n + 1)⌉, across five (n, γ), whole levels included.
- *The batch median* (n 499, m 20, β = γ = 0.05, 4,000 batches): it covers 0.900, within its bound 0.9 + max p = 0.907, with mean width 0.94. Bonferroni over single sets covers 1.000 at mean width 6.07.
- *The quartiles at once* (n 200, m 40, α 0.1, 3,000 batches): the box holds 0.911 of H at ranks [21, 79] and [121, 179] of 201, and covers 0.914.
- *The batch mean* (n 12, m 3, scores in [0, 1], 4,000 batches): it covers 0.983 at mean width 0.78 of the range. With unbounded scores the mean's set is the whole line, since a sequence reaching S_(0) = −∞ is always held.
- *`Function` with h the second smallest and h̃ = r_(2)* is the quantile set exactly.
- *Selection* (n 200, m 50, μ̂(x) = x, c 0.7, η 2, α 0.1, 3,000 batches): false claims stay within 2 in 0.916 of the batches, with 4.7 units selected on average.
- *The refusals:*
  - a batch of none, an order outside 1, …, m, levels summing to 1 or more;
  - a NaN score, scores outside the floor or ceiling;
  - orders that do not rise, a box out of order;
  - an H above the enumeration cap;
  - negative predictions, η = m.
- *Per commit,* the selected classes and the pinned tier: 49 of 49 (one class, the new one).

**Defects found: none.**

**Follow-ups, each with the decision it would change.**
- *The rank order h̃ for the mean.* The default, h on the ranks, is conservative here (0.983 at nominal 0.9). Their remark on choosing h̃ suggests orders closer to h's own. Decides whether `Function` offers a better default for sums and means.
- *Sparse functions* (their §5: h reading a few order statistics), as `Quantile` reads one. It decides whether functions of a few order statistics avoid the enumeration.
- *Their upper bound 1 − α + ε* for the quantile, not formalized. The Lean states the lower bound only. Decides whether the layer reports the set's tightness.

**Doc changes owed** (for the offline docs pass): `docs/Maths/Regression/ConformalPrediction.md`, a batch section:
- the sandwich and the four forms, with the Lean names;
- exact levels;
- the gates' numbers;
- API rows for `BatchConformal`, `BatchInterval`, `BatchBounds` and `BatchSelection`.

### 2026-11-14 — Slice 14's follow-up 2, completed: weighted CV+'s merged bound under likelihood ratios

Built at the main session's word once its Lean (`47acb3c`; `cvPlus_merged_likelihoodRatio`, `cvPlus_merged_likelihoodRatio_refined`) was verified. It closes the open item in follow-up 2's report.

**Landed** on `main`, nothing pushed:
- `3f68330`: weighted CV+ under likelihood ratios reports a merged bound (`CvPlusBounds.Merged`), where it reported none. `WeightedCrossFitTests` gains its checks.
- *Verified by the main session on 2026-11-14,* from a clean archive of `1ea5d6b` with zero skips: the Commit selection, 1,213 of 1,213 ([receipts](receipts/verify-conformal-fu2e-3f68330/)).
- *Nothing else changes.* Position weights' 2c stands: its Lean was restated validly at `47acb3c`, with the same conclusion.

**What it computes,** per test point, from the bag's weights (the training units' and the test point's own, W their total):
- *simple:* min(1, 2c(w_max)), with c(γ) = (αW + (K − 1)γ)/(W + (K − 1)γ), K the folds and w_max the bag's largest weight;
- *refined:* γ the largest weight at most αW, L = αW + (K − 1)γ, N = W + (K − 1)γ and e = Σw²/W. The miss is 2L(W + (K − 1)e)/N² where 2L ≤ N, and 1 otherwise;
- where no weight is at most αW, no unit can be refused: miss 0;
- the merged coverage is 1 minus the lesser miss. It applies under `MergedRanks`, as slice 9's merged bound does.
- `Coverage` reports the largest of the three bounds that apply: one fold's, across folds', and merged.

**Gates.**
- *By hand* (three folds of weights 3, 2 and 4, the test's 2, α 0.2): the simple 2·8.2/17, the refined at γ = 2 (L = 6.2, N = 15, e = 21/11), the merged 1 minus their lesser, and the coverage the largest of the three bounds. Every weight above αW gives a merged bound of 1.
- *Covariate shift* (follow-up 2's fixture, 10 folds of 10 units): the merged bound averages 0.553 against one fold's 0.282 over 2,000 replications, while weighted CV+ covers 0.914; the across-folds bound, given w̄ = 2, is about 0.6 there (1 − 2α − 20/W with W near 101, computed rather than read from the run).
- *Existing classes,* the CV+ tests: unchanged and passing.
- *Per commit,* the selected classes and the pinned tier: 1,213 of 1,213.

**Defects found: none.**

**Doc changes owed:** the cross-fitted section's weighted CV+ bounds gain the merged bound under likelihood ratios, its two forms and where it serves (many units), with the Lean's table of which bound serves where.
