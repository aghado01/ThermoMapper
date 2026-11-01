# The flat-histogram sampler's loose ends (chip brief)

Date: 2026-10-31. Three items Phase 2 left and one the warm start's report did, for the chip that made the warm start robust ([its brief](flat-histogram-warmstart-brief.md)), to run after that brief is reported. Each lands with its checks and ends with a fresh-eyes review, fixed before the next. Its report goes under "Reports" below, newest first. `docs/Maths` and the checklists are the main session's.

Read first: the design's "Phase 1, as built" and "Phase 2, as built" ([flat-histogram-design.md](flat-histogram-design.md)), Phase 2's "For the main session" ([its report](flat-histogram-phase2-brief.md#reports)), and the warm start's ([its report](flat-histogram-warmstart-brief.md#reports)).

## The items

1. **The window's top.** A pilot that never found Friedman's posterior basin reads the window's top below the posterior's ℓ₁ (−778.3 on seed 8017).
   - *What it costs today:* production put 1–13% of each chain's draws above the window, where the slope-1 extension holds the posterior's own law and the b = 1 weights read them exactly. But the learning spent 8–29% of its visits there, which SAMC's update leaves inert.
   - *What the replications add* (2026-10-31): on seeds 8217 and 8417 the pilots' tops sat higher (−773.6 and −771.7), learning spent 1.9–3.5% of its visits above the window, and the b = 1 Kish ESS fell to 22–32% of seed 8017's. So the top trades learning's inert visits against production's b = 1 draws; read both.
   - *Build* an option, not a default: the window's top read again from the learning's own visits, extended at the same node spacing while the gain is still large, within burn-in (D10). A wider margin above the pilot's b = 1 quantile is the simpler alternative; build it too if cheap.
   - *Check:* on the toys, an option left off is bit for bit today's run; on Friedman seed 8017, read each option against today's window: the learning's visits above the window, production's round trips, the b = 1 Kish ESS, ln Ẑ(1), a1's share and the time. One run each, a reading, recommending no default.
2. **The seed-sweep claims.** Fold `tests/samplers/FlatHistogramSweep.cs` and `tests/regression/BassFlatHistogramSweep.cs` into `ComparativeSweep`, as Phase 1's report proposed.
   - *Check:* the registered claims run under `SEED_SWEEP` as before, and nothing else changes.
3. **SMC's lost mass below the lower split** (a reading). On the linear form, fitting SMC from the prior read ln Z(0.10) about 5.5 nats low against the flat histogram's curve, losing mass between b = 0.05 and 0.10 as its lineages fell to 172–220 of 4,000.
   - D05's SMC path is the integrated form, so read the same stages there: ln Z at b = 0.03, 0.05, 0.07 and 0.10 on both forms (4 replicates, 4,000 particles), against the flat histogram's curve on seed 8017 and exact importance sampling from the prior.
   - Report the families at each b, and whether more mutation steps per stage (5 and 20) close the gap.
   - Recommend no default.

4. **The warm start's readout and adapter** (the warm start's "For the main session").
   - *A `TiedNodes` field* on `BassFlatHistogramPilotReadout`: the nodes the warm start tied, empty where the fixed point accepted ĝ, so a fit's readout says whether its start was tied. Correct `BassDiagnostics`' doc of `NodeDraws`, which says only that an empty node is interpolated.
   - *The adapter's `BassFlatHistogramRun.NodeDraws`* duplicates `Wham.NodeDraws`: call the engine's.
   - *The refusal's message names the failing segment:* move the hull reading (`FlatHistogramLimit.Infeasibility`, test code) into the engine as a diagnostic, and have a refusal say at which segment flat occupancy leaves the hull, and by how much, rather than only how far Newton strayed.
   - *Check:* `BassBitIdentityTests`' hashes and every pinned start hold; the field on a tied start (seed 8217's pilot, or a synthetic edge pilot) and on an accepted one; the message on the gap pilots and on a pilot too thin to tie.

## Working rules

As the warm-start brief's, with `artifacts\chip-fh4\` as the artifacts root. Your files are those of the warm-start brief, plus `tests/samplers/ComparativeSweep.cs` and the two sweep files above, and for item 4 `src/maths/regression/spline/bass/BassDiagnostics.cs` (the readout's field and doc) and `BassFlatHistogramRun.cs` (the count and the field's fill).

## Reports

Newest first.

### The flat-histogram sampler's loose ends (chip, 2026-10-31)

**What landed**, on `main`, nothing pushed:
- `9d15376`: today's toy runs (the double well from flat and from a pilot, Potts from WHAM) and the gap pilots' fits (7803 accepted, 7802 tied), hashed bit for bit on the unchanged engine.
- `72788a2`, item 1: the window top's two options, their checks, and seed 8017's reading, written before its runs.
- `0454bf8`, item 1's review fix: `ExtendTopUntil`'s doc says the raised top can climb into a sparse tail.
- `25e5124`, item 2: the sweeps folded into `ComparativeSweep`.
- `f093987`, item 3: its reading's test, written before its run, committed against the tree before item 4.
- `0ec16d3`, item 4: the readout's tied nodes, the adapter's count, and the refusal's segment. Item 3's run then ran on the committed tree.

**1. The window's top.**
- *Built*, two options of `FlatHistogramSpec`, both off by default:
  - `ExtendTopUntil`: in the learning's first iterations, a visit at least a node spacing above the window raises its top to the highest node at or below the visit. The new nodes sit at the same spacing, ln w on the slope-1 extension, so w(ℓ) is unchanged and only SAMC's reach grows. `MulticanonicalWeights.ExtendedAbove`, `SamcLearning.ExtendAbove` and the learning histogram's `Extend` carry it. Production's window is frozen with its weights (D10), and a checkpoint resumes the raised window.
  - `WindowTopMargin`: a pilot window's top at the b = 1 rung's (1 − q)-quantile plus that many of its sd, not held inside the pilot's highest draw (`Wham.Window`'s new argument; the adapter's start shape follows).
- *Checks, the options off:* today's toy runs and gap fits hold their hashes, as do every earlier start pin, `BassBitIdentityTests` and seed 8017's run below (ln Ẑ(1) −841.36 and a1 65.9%, as qualified).
- *Checks, the options on:*
  - The double well over [−30, −4], raised over the first 4,000 of 32,000 iterations, ends at [−30, 0], four nodes higher; its chains start at x = 1, where ℓ₁ = 0. Against the fixed window, the learning's visits above the window fall from 86.8% to 0.0%. Its learned ln w sits 0.040 from θ* over the wider window, against 0.062 for a run declared over it, and ln Z(0.5) and ln Z(1) reweigh within 0.5 errors.
  - The same run in five pieces, cut within the extension, at its end and in production, ends bit for bit where the whole run does.
  - Fixture A, τ held, from a top three spacings low, rises three nodes to F's highest ℓ₁, −15.80. Its learned ln w sits 0.102 from θ*, within 0.15, and ln Z(b) matches the enumeration within its errors.
  - The margin's window arithmetic holds on the double well. On fixture A the adapter's start is the engine's margin window (from a flat start, below).
  - *Read in exploration:* on fixture A with τ held, the gap pilot's highest draw is F's highest ℓ₁, so any margin adds nodes no state reaches. There SAMC's limit does not exist, and the WHAM start is refused even with ties. The option's doc says so.

*Seed 8017, arm (i), one run each, written before the runs* (the qualification's lengths, Thin 200, 4 chains):

| | today's window | `ExtendTopUntil` 50,000 | `WindowTopMargin` 3 |
| --- | --- | --- | --- |
| window, start → end; K | [−1551.0, −778.29], fixed; 204 | → [−1551.0, −751.65] by iteration 4,000; 204 → 211 | [−1551.0, −767.73], fixed; 206 |
| the start | accepted as it stood | the same start | tied (level below node 3, slope 1 above node 202: the top nodes held no pilot draw) |
| learning: visits above the window | 29.2% | 0.1% | 1.1% |
| learning: round trips per chain; flatness | 5, 9, 5, 4; 2.80 | 1, 1, 3, 3; 0.90 | 3, 1, 3, 7; 1.77 |
| production: round trips per chain (all) | 17, 9, 12, 9 (47) | 0, 1, 0, 2 (3) | 4, 5, 6, 9 (24) |
| production: flatness; above the window | 0.39; 6.4% | 0.57; 0.0% | 1.06; 0.7% |
| b = 1 Kish ESS of 500,000 (per chain) | 34,256 (9,186; 5,242; 9,487; 10,341) | 7,462 (1,816; 2,618; 1,099; 1,930) | 14,543 (2,590; 3,271; 2,305; 6,381) |
| ln Ẑ(1) ± error (chains' SD) | −841.36 ± 0.17 (0.34) | −841.45 ± 0.19 (0.39) | −841.11 ± 0.27 (0.54) |
| a1, pooled (χ² p) | 65.9% (0.917) | 69.6% (0.984) | 65.5% (0.999) |
| a1 + x6·x8, pooled (χ² p) | 10.6% (0.197) | 8.5% (0.518) | 10.9% (< 0.001) |
| D10's verdict on π_w | Converged | Converged | Converged |
| time: pilot, learning, production, all | 75, 377, 440, 892 s | 83, 405, 457, 946 s | 76, 375, 455, 906 s |

- The margin's pilot, read before its run: b = 1 rung's 99% quantile −779.19, sd 3.82, highest draw −778.29, so today's top is the cap; with 3 sd the top is −767.73.
- *What it shows* (kept apart from the table):
  - *Raising the top empties the learning's inert visits,* from 29% to 0.1% (extended) and 1.1% (margin), and costs production's b = 1 draws, to 22% and 42% of today's ESS. That is the trade the brief foresaw from the replications.
  - *The extension as built is greedy.* Once the raised nodes flatten, the chains reach higher ℓ₁, which raises the top again. Within 4,000 iterations it climbed 26.6 nats, to −751.65, into a tail of rare states. Production then rarely reached the top bins (two chains' top four nodes at 0.0–0.3 of flat), so round trips fell to 3, two chains with none, failing the qualification's first criterion. The posterior's readings still agree with today's: ln Ẑ(1) within 0.1 nat, a1 within 4 points.
  - *The margin tunnels* (24 round trips, fewest 4) at 42% of today's ESS, its ln Ẑ(1) within 0.25 nats. But a1 + x6·x8 parts across its chains (6.6–14.8%, p < 0.001), as on seed 8417. Its start was tied: its top nodes held no pilot draw and the fixed point refused ĝ.
  - Neither option buys a reading today's window lacks on this seed. The times sit within 6% of one another on a shared machine (builds overlapped the runs).
  - A gentler extension is the open variant, not built: a quantile of the learning's visits rather than the highest, a cap on the rise, or a stop once the visits above fall below a share.
- No default is recommended.

**2. The seed-sweep claims.**
- `FlatHistogramSweep`'s three claims (SAMC's double well; the warm start on the double well and on Potts) are now in `tests/samplers/ComparativeSweep.cs`, and BASS's warm-against-flat claim in `tests/regression/ComparativeSweep.cs`. Both claims are unchanged, and the two sweep files are deleted. The samplers host does not reference the regression assembly, so BASS's claim had to go in the regression registry, a file the brief did not list.
- *Check:* under `SEED_SWEEP` = 1, each claim's rows (own seeds and one chain-axis seed: 8 in the samplers, 2 in the regression) read bit for bit as the separate sweeps' did, at full precision, data and chain seeds included; the registries' other claims run beside them. Their tables now land under each host's suite name.
- *A reading the check surfaced, unchanged by the fold:* two of the folded claims fail at their first chain-axis seed, the same before and after.
  - The warm start's 16k against the flat start's 32k on the double well: ratio 1.254.
  - BASS's warm start's 4k against the flat start's 8k: ratio 1.188.
  - The Potts and SAMC claims hold there (ratios 0.30–0.72).

**3. SMC below the lower split, on both forms** (a reading; `f093987`, run at `0ec16d3` in 20 min 32 s).
- *What it reads:* ln Z(b) at b = 0.03, 0.05, 0.07 and 0.10 on each form, by four methods:
  - SMC from the prior with a stage forced at each b (4 replicates × 4,000 particles, C4's defaults), at 5 sweeps a stage (C4's default, Phase 2's reading) and at 20, with each replicate's families at each b;
  - importance sampling from 200,000 exact prior draws (Phase 2's streams, reproducing its linear-form column to the digit);
  - seed 8017's arm (i), rerun bit for bit (ln Ẑ(1) −841.36).
- *The integrated form has no flat-histogram curve of its own:* the two forms' Z(b) differ between b = 0 and 1, the integrated form's the lower by Hölder. There the flat histogram's 500,000 production draws are reweighted by that form's tempered marginal, Ẑ(b) = Σ m_b(S, τ)/w(ℓ₁) over Σ 1/w(ℓ₁), each draw's m_b from its structure and τ (its rebuilt ℓ₁ checked against its level). On the linear form the same code reproduces the curve to 6 digits.

*The linear form:*

| b | importance sampling | flat histogram (Kish ESS) | SMC, 5 sweeps: pooled; replicates | families | SMC, 20 sweeps: pooled; replicates | families |
| --- | --- | --- | --- | --- | --- | --- |
| 0.03 | −45.00 ± 0.01 | −45.01 ± 0.00 (26,593) | −45.02 ± 0.01; −45.00 to −45.04 | 2,242–4,000 | −45.02 ± 0.00; all −45.02 | 4,000 |
| 0.05 | −73.02 ± 0.21 | −73.02 ± 0.01 (165,678) | −73.08 ± 0.05; −72.99 to −73.21 | 1,393–1,478 | −73.04 ± 0.01; −73.02 to −73.05 | 2,169–2,228 |
| 0.07 | −98.84 ± 0.50 | −97.61 ± 0.01 (136,703) | −97.89 ± 0.12; −97.59 to −98.19 | 478–524 | −97.73 ± 0.03; −97.63 to −97.77 | 881–1,316 |
| 0.10 | −136.52 ± 0.85 | −125.93 ± 0.06 (49,007) | −131.51 ± 0.22; −130.89 to −131.91 | 188–217 | −129.83 ± 0.99; −126.96, −131.03, −130.11, −131.22 | 69, 519, 496, 639 |

*The integrated form (D05's path):*

| b | importance sampling | flat histogram, reweighted (Kish ESS) | SMC, 5 sweeps: pooled; replicates | families | SMC, 20 sweeps: pooled; replicates | families |
| --- | --- | --- | --- | --- | --- | --- |
| 0.03 | −54.29 ± 0.00 | −54.29 ± 0.01 (2,991) | −54.29 ± 0.00; −54.29 to −54.31 | 4,000 | −54.29 ± 0.00; −54.28 to −54.29 | 4,000 |
| 0.05 | −84.48 ± 0.02 | −84.47 ± 0.01 (7,128) | −84.47 ± 0.01; −84.45 to −84.49 | 2,154–4,000 | −84.47 ± 0.00; −84.47 to −84.48 | 4,000 |
| 0.07 | −113.19 ± 0.10 | −112.95 ± 0.01 (32,146) | −112.95 ± 0.04; −112.87 to −113.05 | 1,359–1,428 | −112.94 ± 0.01; −112.92 to −112.97 | 2,081–2,149 |
| 0.10 | −153.66 ± 0.33 | −150.17 ± 0.02 (33,452) | −150.59 ± 0.14; −150.26 to −150.89 | 442–576 | −150.58 ± 0.03; −150.51 to −150.63 | 823–855 |

- SMC's stages: 63–67 (5 sweeps) and 67–75 (20) on the linear form, 64–67 on the integrated; 58 s, 116 s, 50 s and 92 s.
- *What they show* (kept apart from the tables):
  - *On the linear form the loss is SMC's and grows past b = 0.05:* 0.28 nats at 0.07 and 5.6 at 0.10, with 188–217 families left at 5 sweeps. 20 sweeps narrow it to 3.9 nats pooled but do not close it. One replicate came 1.0 nat short of the flat histogram while its families fell to 69; the other three stayed 4.2–5.3 low. That pattern is rare high-likelihood mass that a population reaches by chance, not mixing at a fixed b that more sweeps repair.
  - *On the integrated form, D05's path, SMC loses about 0.4 nats by b = 0.10* (−150.59 ± 0.14 and −150.58 ± 0.03 against −150.17 ± 0.02), keeping 442–855 families. At 0.07 and below it agrees with the flat histogram within 0.01, and 20 sweeps tighten its replicates without moving the 0.4 nats.
  - The three estimators fall in the order Phase 2 read: importance sampling reads low first, SMC next. The integrated form's flat-histogram reading rests on one seed's draws, reweighted.
  - Importance sampling from the prior holds to b = 0.05 on both forms and reads low beyond, as before.
  - So the 5.5 nats Phase 2 read belong to the linear form; on the path D05's SMC takes, the same stages lose about 0.4 nats.
- No default is recommended.

**4. The warm start's readout and adapter.**
- *`BassFlatHistogramPilotReadout.TiedNodes`:* the nodes the warm start tied, empty where SAMC's fixed point accepted ĝ. The engine reports them: `MulticanonicalResult.TiedNodes` and `MulticanonicalCheckpoint.TiedNodes`, from `Wham.WarmStart`'s out list, so the adapter fills the field without a second WHAM solve, and a resumed fit reports the checkpoint's.
- *`BassFlatHistogramRun.NodeDraws`* is gone; the adapter calls `Wham.NodeDraws`, the same sum, asserted equal to the readout's.
- *The refusal names the failing segment:*
  - `SamcLearning.Infeasibility` is the hull reading, moved from test code into the engine as a public diagnostic.
  - `FixedPoint`'s refusal and the warm start's (both kinds) now read "SAMC's limit does not exist: flat occupancy leaves the hull of the measure's hats at segment m (ℓ … to …), by g of its mass, and the solve stopped with a node's occupancy d from flat". Inside the hull they say so and give d.
  - The adapter's wrapper (window, K, b_lo, empty nodes) is unchanged.
- *Docs:* `BassDiagnostics`' readout says the window's top may have been raised and that a WHAM start may be tied; its `NodeDraws` no longer says only that an empty node is interpolated; `TiedNodes` is documented.
- *Checks:*
  - `BassBitIdentityTests`, every pinned start (seed 8017's included, `BASS_QUALIFICATION`) and today's runs and fits hold.
  - The field: gap pilot 7802 tied [0, 7, 10, 11, 13, 14, 15, 16], equal to the warm start's own list and kept through a cut and resumed fit; 7803 tied [], accepted.
  - The message:
    - 7801 at segment 13 (ℓ −22.504 to −21.3867), by 0.0736;
    - 7802 at segment 16 (ℓ −19.1521 to −18.0348), by 0.127;
    - the pilot too thin to tie at segment 3 (ℓ −32.5261 to −31.5102), by 0.126, "and fewer than two nodes hold 10 of the pilot's draws to tie the rest to";
    - a three-node toy at segment 0, by 1.5.

**Suites** (Release, at `0ec16d3`): Samplers 282 passed (42 s); Regression 1,409 passed (5 min 35 s); the BASS oracle classes (`Maths.Oracle.Tests.Bass*`) 90 passed (5 min 42 s); seed 8017's pinned start (`BASS_QUALIFICATION`) passed (1 min 16 s). Item 2 folded one sweep test into each host's registry, item 3 added one test to Regression and item 4 one to Samplers, so the counts read as after item 1 (Samplers 282, Regression 1,409 at `72788a2`).

**Review passes and their fixes:**
- *Item 1, before its runs:* off is bit for bit (pins); a resume across the extension; non-finite levels never raise the top.
  - The double well's top landing at exactly 0 was traced to its chains' start at ℓ₁ = 0, not a fault.
  - The margin's first adapter test (WHAM, fixture A, τ held) was refused for a real reason: nodes past F's highest ℓ₁. The option's doc says so, and the test reads the adapter from a flat start.
- *Item 1, after its runs:* the extension's climb is documented on the option (`0454bf8`).
- *Item 2:* the claims compared row for row; no reference to the deleted files remains.
- *Item 3:*
  - Its checks held: the flat histogram's rerun reproduced −841.36; every one of the 500,000 draws' rebuilt ℓ₁ equalled its level; the linear-form reweighting equalled the curve; and the prior draws reproduced Phase 2's importance-sampling column.
  - What it cannot settle: whether the integrated form's 0.4 nats at b = 0.10 is SMC's loss or a seed's, its reference being one flat-histogram run reweighted. Its Kish ESS at b = 0.03 is 2,991, the lowest of the readings, where it agrees with importance sampling.
  - The forced marks are fewer than Phase 2's, so its 5-sweep linear values differ from that run's by up to 0.1.
- *Item 4:*
  - The qualification sweep's printout parsed the refusal's last word as the departure; it now reads the hull reading, asserting the message names the same segment.
  - `ReadWarmStart` prints the message whole.
  - A gap read as 1.5000000000000004 is compared to 12 places.
- *Process:* one compile, after an edit script failed, ran with an empty artifacts root and wrote `bin/` and `obj/` under four `projects/` folders and `tests/Maths.Samplers.Tests` (21:07). They were removed; nothing else was touched.

**For the main session.**
- *Two registered warm-start claims fail at their first swept chain seed* (item 2's reading), the double well's (1.254) and BASS's (1.188). They predate this chip; more seeds under `SEED_SWEEP` would say whether "half the iterations from the warm start beat the flat start's full length" holds for the sampler or for its seeds.
- *The extension rule as built climbs into the sparse tail* on Friedman; a quantile or capped rule is the open variant (item 1).
- *SMC's loss below the split is mostly the linear form's* (item 3): D05's integrated path loses about 0.4 nats at b = 0.10 on this reading, against the linear form's 5.6.
- *Engine API:* `MulticanonicalResult` and `MulticanonicalCheckpoint` gained `TiedNodes`; `FlatHistogramSpec` gained `WindowTopMargin` and `ExtendTopUntil`, which enter every run's fingerprint and a flat-histogram config's JSON, so checkpoints taken before `72788a2` do not resume.
- `docs/Maths` (Samplers.md, BayesianRegression.md) says nothing of the two options, `TiedNodes` or `Infeasibility`, and the checklists are untouched.
- `artifacts\chip-fh4\` is deleted.

Recommend no default.
