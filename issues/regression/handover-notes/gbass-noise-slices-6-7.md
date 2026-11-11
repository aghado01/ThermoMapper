# GBASS noise, slices 6 and 7: BARS under the noise, and interweaving (2026-11-08)

Landed: slice 6 in `957629b` (code) and `2b2f56c` (docs, report); slice 7 in `211d768` (code) and `c52c170` (docs, ledger, checklist, report). Stopped after slice 7, as asked. E1's brief is built; the checklist marks it.

- **Slice 6, BARS (OQ20):**
  - `StudentTWeights` now draws through `StudentTMixing`, number for number, and takes the law's marginal.
  - `BarsNoise` covers quantile (AL, skew held), skewed (β drawn) and pure scale mixtures on any law, on `BarsConfig.Noise`.
  - Under a skew σ² is held on the chain, the knot step reads the marginal given σ², and the priors read the data once (`WeightedNormalModel.ReadingFrom`).
  - BASS and BARS now share the σ² and β conditionals (`VarianceMeanConditionals`); BASS's arithmetic is unchanged.
  - Gates: Geweke on four noises; quantile recovery (10.5% / 50.2% / 89.5% of the data below the q = 0.1 / 0.5 / 0.9 curves); skew recovery (three of four data sets covered); bit-for-bit resume. The whole assembly ran once at the slice's end: 2,189 of 2,189.
- **Slice 7, interweaving:**
  - `BassNoise.Interweave` (opt-in) adds `BassInterweavingBlock`: z = σ·v held, σ by one slice step on ln σ, then v = z/σ. Each law gains `IMixingLaw.LogPrior`.
  - The step is valid: alone (σ moved only by it) Geweke reads 1.75 (skewed horseshoe) and 2.31 (GH), and the interwoven sweep passes for t, GH and the quantile noise.
  - **The main session's gate:** the skewed horseshoe at b = ½ still fails on seed 20261107 (4.96, and 10.12 at 4×) and passes on 20261108 (3.07). It also fails on that seed without the step (4.27), so it's the horseshoe's v mixing slowly. The β·v joint move stays a follow-up, with these readings as its measurement.
  - **What it buys:** σ²'s ESS per second doesn't improve (0.04–1.15×; about 20× worse under the quantile noise). mean ln v improves 3–6× under the skewed GH and horseshoe. A fit costs 1.3–2.5×. It stays opt-in. The follow-up is z′ = σ²·v under pure scale mixtures, where the likelihood is free of σ.
  - Slice-end run: the whole assembly, 2,241 of 2,241.
- **Tests:** from slice 7 I followed the new rule in `CONTRIBUTING.md` (`scripts/test-select.ps1`), run with `-SourceRoot` on a snapshot of HEAD plus my files, because other sessions' work in progress sat in the tree.
- **Harness lessons:**
  - The first BARS cost row was wrong twice: the test thread's clock missed the pool's chains, and the process's CPU time counted the pool's spinning. The row now steps one chain on the test thread (`BarsEnsemble.ChainOn`).
  - `git archive` stamps files with the commit's time, so a snapshot over an old build output silently kept a stale DLL. A fresh snapshot now drops its build output.
- **Open for the user:** the flat histogram under a noise (slice 5's question), still refused.

- **Addendum (2026-11-08): the mixture's tail table past its ends**, `a670046`, its own commit after slice 7.
  - **Defect:** `TailTable.Read` read the side near 1 by its own integral past the 10⁻¹² ends. That lost 1 − F's digits and read up to +3.0·10⁻¹⁵ (generalized hyperbolic) or +6.3·10⁻¹⁵ (Normal-Wald).
  - **Fix:** that side is now ln(1 − e^(far tail)), the far tail integrated on its own side.
  - **Checked** on four families: no log tail above 0. Distributions' selected class 31 of 31; regression selection with the pinned tier 1,724 of 1,724.
- **Addendum (2026-11-08): the move along β·v**, `3be066e` (code) and `49b0cd3` (docs, report).
  - `BassNoise.SkewRescale` (opt-in, drawn skew only) adds `BassSkewRescaleBlock`: (β, v) → (β/κ, κv), the Jacobian κ^(N−1) in A, ln κ ~ N(0, h(x)²) with h from the state's curvature (Hastings ratio in A).
  - **Valid:** its orbit's exact law holds on a frozen state at b = 1 and ½ (|z| ≤ 0.87, 44% acceptance), and Geweke passes under GH, Normal-Wald and skewed t, b = 1 and ½.
  - **Measurement:** it doesn't pass where the scan failed. The skewed horseshoe on seed 20261107 reads 11.13 at b = ½, and 4.04–4.69 at 4×; the scan alone also fails at b = 1 at 4× (4.13). It doesn't raise σ²'s or β's ESS per second, and the ladder's cold pair swaps 0 either way. **So it stays opt-in, and the readings stand as the skewed horseshoe's mixing limit.**
  - **Tests:** per-commit selection 1,558 of 1,558; slice-end whole assembly 2,317 of 2,317.

- **Addendum (2026-11-08): the move along β·v, corrected**, `1a8459c` (code, tests, gates) and `717de68` (page, report, ledger anchors). This supersedes the previous addendum's readings.
  - **Defect (the main session's review):** after an accepted move, σ²'s, the shape's and interweaving's blocks read coefficients drawn before it. **A second one beside it:** κ's scale read the coefficients through R, which is not exact on the a-marginal.
  - **Fix:** the coefficients are drawn afresh after an accepted move, and the scale's likelihood term is b·N/2, from x alone. The other exact form, given a at O(N), was built, measured and removed: it mixed v's mean 1.7× slower per second under the skewed horseshoe and failed its Geweke at b = 1 (4.01, 4.21).
  - **Skewed horseshoe Geweke, kept form:** passes at b = 1 (1.93, 1.43 at 4×, 2.41) and ¾ (2.15, 2.69, 2.59), and these are now gates. It fails at b = ½ (9.71, 4.29, 5.18). The scan alone failed b = 1 at 4× (4.13).
  - **Stays opt-in:** per second it buys v's mean, 3.5×, and not ln σ², β or M, at 1.2–1.6× a run's time.
  - **Interweaving** reads a only conditionally (confirmed).
  - **Tests:** per-commit 1,618 of 1,618; slice-end 2,320 of 2,320.

- **Addendum (2026-11-08): the move along β·v is the skewed horseshoe's default**, `bf21578` (code, tests) and `de4d965` (page, report, GBMARS row, ledger anchors). This follows the user's "Default there".
  - `BassNoise.SkewRescale` is `bool?`. Null takes `SkewRescaleByDefault`, which is GBP(1, ½, ½) with β drawn; true or false declares it; true under a held skew is refused. `RescalesSkew` is what the sweep runs.
  - Null isn't serialized, so earlier configs keep their JSON and fingerprint. An undeclared skewed-horseshoe config now draws differently under the same fingerprint.
  - `BassResolved.AutoPicks` logs `Noise.SkewRescale` on or off with its reason wherever a drawn skew is left to the default.
  - **Resume:** neither a refusal nor a restore was built. BASS checkpoints are in-memory live state that is never serialized, so none from before the change can reach the new build. Within a build the setting resolves from the fingerprinted config.
  - No pinned hash moved. The reach Geweke's plain forms and interweaving's measurement declare the move off, so their readings stand. The GBP skewed-horseshoe Geweke under the default reads 1.99.
  - **Tests:** per-commit 1,627 of 1,627; slice-end 2,329 of 2,329.

- **Addendum (2026-11-08): the mixture's tail table across a cusp at 0**, `019738d` (fix, sweep, conformal test unskipped) and `d0728d6` (docs, report). The conformal chip found it in its slice 5.
  - **Cause:** the table's grid near 0. A mixing density ~ v^(α−1) at 0 gives the noise a |x|^(2α−1) singularity at 0. Only poles (α ≤ ½) stepped into 0. A cusp or kink (α > ½, e.g. GBP(1, 1, 1.5)) took ordinary steps, so a drifted median left one cubic cell spanning 0.
  - **Fix:** every finite left exponent takes the pole's walk, with 0 as a node. Only a pole's cells are still read by the integral.
  - **What moves:** GBP tables at pa > ½ and variance-gamma tables at p > ½. Unchanged bit for bit: the horseshoe, the Strawderman–Berger and GIG mixtures with b > 0.
  - **Sweep (40 laws):** worst 5.8·10⁻² before, 7.9·10⁻⁶ after. Commit selections: distributions 32/32, regression 1,952/1,952.

- **Addendum (2026-11-08): E21's ridge blocks, plan and slice 0.** Brief `2d788b2` (`issues/regression/bass/ridge-blocks-brief.md`); slice 0 `7cf1761` (code), `98c39b5` (docs: new `docs/Maths/RandomEffects.md`, report).
  - `WorkingResponse.Whitening` carries the block; BassDesign whitens by the reading it's on; ProposeFollow at another scale refits as for moved weights.
  - Fixed: hinge prefix tables were read under a block (wrong border), at n ≥ 500 cases.
  - Held-s fits bit for bit (three hash pins). Per-commit 1,845/1,845; slice-end 2,640/2,640.
  - **Slice 1 must settle:** level centring makes priors read the working response, so ProposeFollow refuses; a drawn s needs priors frozen on the unwhitened data.
  - **Left as found:** `BassSweep.cs` holds the KHAOS chip's uncommitted edit.

- **Addendum (2026-11-09): E21 slice 1 (drawn s), IN PROGRESS, uncommitted** (usage limit reached mid-slice).
  - **Source written, compiles:**
    - `BassGroupScalePrior` (half-Cauchy on s) and `BassRandomIntercept.ScalePrior`; `RandomInterceptBlock.At(s)` and `LogScaleInformation`.
    - `BassGroupScaleBlock`: MH on ln s via ProposeFollow, in τ's slot.
    - `BassModel`: `ScaleDrawn`, `ReadingMoves`, `BaseReading`, `BlockAt(draw)`; under drawn s the reading carries the priors' frozen reference (`BassPriorReference.Moments`, or `About(0)` under the zero centre).
    - `BassDraw.GroupScale`, with `WorkingOn` rebuilding the reading at the draw's s; the effects drawn from the reading's block.
    - `BassPriorSampler.DrawPrior` draws s; SMC and the ladder adopt the sampler's design under `ReadingMoves`.
    - Config refusals (A <= 0; flat histogram with drawn s); the resolution log; prediction at the draw's s.
  - **Files (all uncommitted):** BassRandomIntercepts, BassModel, BassObservationModel, BassSweep, BassState, BassPriorSampler, BassSmc, BassLadder, BassConfig, BassResolved, posterior/BassPrediction, tests/regression/RandomEffectsTests.cs.
  - **Next:**
    - Build, then run RandomEffectsTests, which has the slice 1 gates written but not yet run: switch, orbit law, draws, resume, ladder law with ln s, SMC evidence.
    - Then the held-s pins, KHAOS's RandomInterceptTests, the Commit scope and a commit.
    - Then Geweke, recovery and the cost row (qualification), the docs (RandomEffects.md section) and the report.
