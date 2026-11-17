# Generative-prior probes: the numerical chip's third brief

Date: 2026-11-16. Status: **G1–G5 complete at the registered finite scope**; k = 16 calibration was conditionally skipped by its original gate. Originally posed at the user's request on 2026-11-16 ("Yes, pose both", to Lean and numerical briefs on the compressed-sensing results). Drafted by the coordinating session. Its sibling is the Lean chip's [generative-prior statements brief](generative-prior-statements-brief.md). The direction it serves is in the [approach brief](approach-brief.md#the-compressed-sensing-guarantees-used-by-analogy-not-yet-tested).

## Why

The approach leans on theorems that say a code of dimension k is recoverable from a number of measurements linear in k, by descent, when two conditions hold: each generator layer acts on pairs of inputs through their angle (the weight distribution condition), and the measurement map is a near-isometry on the generator's range (the range restricted isometry condition). The pilot tested a Bayesian analogue on a parametric synthesizer. Nothing has yet run the papers' own object, measured either condition for anything of ours, or put the engine's inversion beside the papers' algorithm. These probes do that, on the engine as built.

## Sources

- **The papers,** LaTeX under `D:\aghado01\codex-scientiae\supellex-private\gauntlet\voroninski`:
  - `1705.07576v3` (`deep_cs_arxiv.tex`), "Main Results": the landscape theorem, both conditions, the sizes under which they hold;
  - `1812.04176v1`: Algorithm 1 (subgradient descent with a sign check), its step size, and its experiments (two layers of 250 and 600, m = 150, k from 10 to 100, success at relative error below 10⁻³);
  - `2008.10579v1`: the same from magnitudes;
  - `1805.08855v2`: the k/n noise reduction.
- **The Lean chip's statements** as they land: the angle map g, the matrix Q_{x,y}, the exact zero set of the expected gradient, and ρ₂ = 1/π for the impostor's place at depth 2.
- **The engine:** `src/maths/regression/spline/bass/calibration/` (the emulator contract, the nested sampler, the tempering ladder, cut SMC), `DirectionalHinge.cs` and `DirectionalTermSpace.cs` for the direction prior.
- **Your own pilot:** its `SynthesizerEmulator`, its resonance synthesizer, and its receipts.

## The probes

**G1. The papers' regime.**
- A random expansive generator without bias and a Gaussian measurement matrix at the paper's sizes; Algorithm 1 as printed.
- Success against k at fixed m, and against m/k; depth 2 and 3; the magnitude objective beside the linear one.
- Without the sign check, where descent from random starts ends: near x₀, or near −ρ_d·x₀. At depth two compare the linear objective with ρ₂ = 1/π and the magnitude objective with its own ρ₂ ≈ 0.66915 (`dpr.tex`, definition of ρ_d).
- *Decisions:* the baseline curve every later probe is read against; whether the impostor sits where the theory puts it.

**G2. The same inverse problem through the calibration stack.**
- The generator behind the emulator contract, as the pilot's reference inversion put the synthesizer there; the code on a declared box; k from 2 upward as far as the samplers go; m at several multiples of k.
- The nested sampler, the tempered ladder and cut SMC. Read recovery error, 90% coverage, convergence, mass near −ρ_d·x₀, and wall time against k.
- At retained states, the cosine between −∇f and the direction to x₀.
- *Decisions:* how the engine's inversion scales in k against descent; whether the impostor traps chains and tempering clears it; whether a gradient-informed move on the code is warranted.

**G3. The weight condition's residual.**
- For a layer of n rows on ℝᵏ, the largest spectral gap between Σᵢ 1{wᵢ·x > 0}·1{wᵢ·y > 0}·wᵢwᵢᵀ and Q_{x,y} over many sampled pairs across angles, against n/k.
- Layers: Gaussian rows; rows uniform on a sphere; rows drawn from the engine's direction prior (sparse unit directions, at most three active coordinates), without and with thresholds. Scale each so that Σᵢ wᵢwᵢᵀ is near the identity.
- Include x = y = e₁ and read the second diagonal alongside the sampled maximum. Under the threshold-free sparse prior the population gap is at least 3/7 at k = 8. For thresholds, declare the input radii and training design used to draw knots. [Derivation and assumptions](approach-brief.md#analytical-controls-for-the-generative-prior-probes)
- *Decisions:* how wide a layer must be for a given ε; whether a layer of directional hinges under its prior is near enough to angle-preserving for an inversion through it to inherit the benign landscape, and where its residual floors.

**G4. The isometry constant for structured views.**
- On the pilot's resonance manifold and on G1's generator range, sample secant pairs and estimate the smallest ε for which each measurement map satisfies the isometry condition, after its best global scale, against m.
- Maps: Gaussian rows; rows of a linear short-time transform at one window length, as real and imaginary parts; the same across three window lengths; each of those whitened by its atoms' Gram; the pilot's kept-bin subsets.
- Count real rows and rank. For each pair of secants normalize the inner-product distortion by the two original norms, and record how the global scale was chosen. Check linear Gram whitening against the row-space projector. For kept log powers, take the difference of the two transformed signals and report the nonlinear distortion separately. A finite sample maximum bounds the uniform worst case from below.
- *Decisions:* whether structured views are near-isometric on these manifolds at the m where recovery worked; whether whitening lowers ε; how far they sit from Gaussian rows.

**G5. The noise reduction, from receipts already held.** In the pilot's reference inversions, examine squared reconstruction error against the clean held-back log powers by k/m and waveform noise level. S21 gives k/n for a fixed projector under isotropic noise in the same coordinates. The pilot adds noise to the waveform and scores log powers. No new runs. *Decision:* whether the receipts support a k/m trend in this prediction task.

## Scope

- **Probes live under `artifacts/generative-prior/`,** with small settings, probe sources and result summaries tracked at the matching path under `issues/regression/receipts/`.
- **The engine is read and called as built.** `src/`, `tests/` and `docs/` stay as they are. A gap comes back as a primitive proposed in the report, stated for any application.
- **Ground truth comes from outside the estimator:** the code that made each signal, a rank or norm computed independently.
- **Gates come first.** Before each probe's first reading, write its gate into the report: sizes, seeds, tolerances and the rerun rule for a marginal failure. A gate changed after a failure is recorded as a change.
- **Order:** G1, G3, G4, then G2, then G5.

## Done

- G1 to G5 each run against its gate, with its decisions answered or the reason one could not be.
- A report appended to this brief: the gates as pre-stated, the readings, what each decision came to, the primitives proposed, and what was left unread or unrun. It ends in an **Exposition** section of self-contained segments the user can splice into docs later, technical and concise. It carries no pointers into `docs/`, and `docs/` is not opened.

## Working rules

- Commits land on `main` by pathspec, each one informative about its work, with no co-author trailer; nothing stays staged.
- Builds set `ArtifactsRoot` inline on each command. Heavy runs go at below-normal priority with capped threads.
- Every run is seeded and reproducible from the tracked settings on a pinned commit.
- Each step ends with a fresh-eyes review of its own output, fixed before it is reported.
- Messages are terse events; results travel in the report.

## G1 report — first checkpoint, 2026-11-16

### Protocol registered before the first numerical reading

**Published comparison.** The noiseless success plot has a value of 1 at k = 10 in 30 runs. The commented table in the same source confirms the count. Huang, Hand, Heckel, and Voroninski's experiment uses **Algorithm 2** (practical gradient descent). Algorithm 1 admits a Clarke subgradient at a kink. Algorithm 2 always takes the masked derivative with a strict-positive ReLU gate. Away from kinks the two steps coincide. The authors report no exactly zero preactivations in their experiments, and qualify the convergence argument at kinks. See `1812.04176v1-tex/CSDL_PH10_TechRep.tex`, lines 488–508, 665–693, 790–805 and 815–824.

**One cell.** Use independent entries, real binary64 arithmetic, k = 10, n₁ = 250, n₂ = 600, m = 150. Both layers are bias-free ReLUs. The entries of W₁, W₂ and A are centered normals with variances 1/250, 1/600 and 1/150, respectively. Draw the true code x★ from N(0,I₁₀) and set y = A ReLU(W₂ ReLU(W₁x★)) with zero noise. Algorithm 2 specifies an arbitrary nonzero start. Here x₀ is a fresh, independent N(0,I₁₀) draw for each problem.

**Seeds.** Use `SeedTree.Derive(181204176, 30)` in repetition order 1 through 30. For each repetition seed, derive five child seeds in order W₁, W₂, A, x★, x₀. Each initializes its own `Xoshiro256PlusPlus` stream sampled by `Normal.SampleStandard`. Fill matrices in row-major order. Redraw all five objects every repetition. The receipt will include the resolved seeds.

**Iteration and scoring.** For f(x) = ½‖AG(x) − y‖² compare f(−x) with f(x) *every iteration*. Negate only on strict improvement. At the chosen sign z, compute the gradient by backpropagation using each preactivation > 0. Stop and return z when ‖∇f(z)‖₂ < 2⁻⁵² (the binary64 spacing at one). Otherwise take the full step z − ν∇f(z), ν = 2ᵈ/d² = 1. At 50,000 updates return the last updated code. Test relative Euclidean code error against the generating x★ only after termination. Success requires error strictly below 10⁻³.

**Gates.** Reproduce the observed 30/30 count. Score every repetition including iteration-limit stops. Independently check both sign objectives and the gradient at the selected initial sign using scalar forward sums and the fixed-mask directional derivative along each coordinate. Require absolute discrepancies ≤ 10⁻¹⁰ × max(1, |reference|) per scalar. Repeat the objective and gradient checks at each returned code. Require finite readings and zero exactly zero layer preactivations on the selected gradients for the claim about equivalence on this run. Count sign flips, the two stop reasons and exact zeros separately from recovery.

**Rerun rule.** A 29/30 outcome is marginal. Replay the same 30 seeded problems once and inspect the failed endpoint, stopping gradient, sign choices and objective checks. For a larger discrepancy or any numerical check failure, diagnose implementation, sources and failed rows first. Retain the original count and gate. A replay tests reproducibility and adds no independent trials. No new seeds or parameter tuning in this checkpoint.

### Readings and decision

**Passed the preregistered 30/30 recovery gate.** All 30 final codes were within 9.24 × 10⁻¹⁶ relative Euclidean error of their independently drawn generating code. The paper's k = 10 point is reproduced with fresh matrices, codes and initializations. [Primary summary](../../../../artifacts/generative-prior/g1-k10-primary/results/summary.json) · [Primary per-repetition rows](../../../../artifacts/generative-prior/g1-k10-primary/results/trials.jsonl)

| Reading | Observed |
| --- | ---: |
| Successful final codes at error < 10⁻³ | 30/30 |
| Relative code error, min / median / max | 2.5074 × 10⁻¹⁶ / 5.7276 × 10⁻¹⁶ / 9.2344 × 10⁻¹⁶ |
| Gradient stops | 29, after 189–417 updates |
| Update-limit stops | 1, at 50,000 updates |
| Median updates, including the capped run | 291 |
| Strict sign flips / comparisons | 17 / 58,659 |
| Exact zero preactivations at selected signs | 0 |
| Largest scaled discrepancy in the objective and gradient checks | 7.7716 × 10⁻¹⁶ (gate 10⁻¹⁰) |

All 17 sign changes occurred at the first comparison, in distinct repetitions. Every selected sign was inside a differentiable region of the generator, so Algorithm 2 took the direction specified by Algorithm 1 throughout this cohort. The effect of the sign check on recovery remains an open experimental question.

**Repetition 13 reached the update limit.** Its seed was −1738237470, its relative code error was 6.8430 × 10⁻¹⁶ and its final objective was 5.7565 × 10⁻³¹. The final gradient norm was 2.617852353400183 × 10⁻¹⁶, or 1.179 times the stopping threshold. Regenerating its network and observation from the recorded streams and applying one signed update to the returned code reproduced all ten coordinate bit patterns. The update had rounded to a fixed point. This repetition counts as a recovery under the original gate.

**Checks.** A separate calculation from the recorded endpoint vectors confirmed all 30 relative errors and success flags. The stored first-step sign decisions agreed with the two signed objectives, each gradient stop satisfied the strict norm bound, and all 150 component seeds were distinct. A fresh build and replay from the pin produced byte-identical primary and replay summaries and trial rows on this machine. [Replay summary](../../../../artifacts/generative-prior/g1-k10-replay/results/summary.json) · [Replay trial rows](../../../../artifacts/generative-prior/g1-k10-replay/results/trials.jsonl)

### Replay

The protocol landed in `dfb8db4`. The code, settings and engine used for both runs are pinned at `3842563b259d32ca9a997af6890151d0aac15526`. Run the following from the repository root in PowerShell with the .NET 10 SDK. Use a fresh run label on each invocation.

```powershell
./issues/regression/receipts/generative-prior/replay.ps1 -Revision 3842563b259d32ca9a997af6890151d0aac15526 -RunName g1-k10-review
```

The command archives the pinned sources under `artifacts/generative-prior/g1-k10-review/source/`, builds in Release and writes `results/summary.json` and `results/trials.jsonl` under that run. The numerical loop is serial and runs below normal priority with one logical processor exposed to .NET. The observed runtime was .NET 10.0.5, x64. Each trial row contains its seeds in registered order, truth, start, estimate, sign-flip iterations (one-based), stopping reason and check readings. The tracked copies use repository line endings. [Settings](../../receipts/generative-prior/settings/g1-k10.json) · [Tracked summary](../../receipts/generative-prior/results/g1-k10/summary.json) · [Tracked trial rows](../../receipts/generative-prior/results/g1-k10/trials.jsonl)

### Review boundary

This supplies the first linear baseline at m/k = 15. The remaining k and m/k values, depth three, sign ablation, magnitude objective and G2–G5 are queued. For the future negative-ray comparison at depth two, the linear prediction is ρ₂ = 1/π ≈ 0.31831. In the magnitude paper's `dpr.tex`, lines 623–625, the angles starting at π give θ₁ = π/2 and θ₂ = arccos(1/π). Substitution gives ρ₂ = 2 sin(θ₂)/π + (π − 2θ₂)/π² ≈ 0.6691522942. This signed linear run did not measure either negative-ray location. No engine primitive was needed for this checkpoint.

### Doc changes owed

The numerical baseline, stopping-rule observation and distinction between the linear and magnitude radii are available for the next documentation batch.

## G1 report — linear curve and paired sign comparison

### Protocol registered before the first new numerical reading

**Finite-realization radial control.** Before running descent, evaluate y = AG(x★) and v = AG(−x★). Record vᵀy, vᵀv and a★ = max(0, vᵀy/vᵀv) when vᵀv > 0, or mark the ray flat when v = 0. Positive homogeneity makes a★ the exact minimizer of f(−a x★) over a ≥ 0 for that realization. Keep both the population-radius deviation from 1/π and deviation from a★ alongside each endpoint's orientation.

**Domain and reference.** Keep the first checkpoint and its gate as recorded above. Add independent signed cohorts at k = 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, with 30 problems per k. Reuse the first signed k = 10 cohort. Keep widths 250/600, m = 150, zero noise, the published variances, independent Gaussian truth and start, unit step, strict objective comparison, relative-error success threshold 10⁻³, gradient threshold 2⁻⁵² and 50,000-update limit. The paper's noiseless reference counts, in increasing k including 10, are 30, 30, 30, 30, 30, 30, 30, 30, 30, 28, 28, 25, 17, 7, 4 out of 30 (`CSDL_PH10_TechRep.tex`, commented table at lines 815–824, corroborated by `PSR_versus_k.eps`).

**Seeds and pairing.** The original k = 10 seeds remain `SeedTree.Derive(181204176,30)`. Derive 14 cell masters by `SeedTree.Derive(181204177,14)`, in the added-k order above, then 30 repetition seeds from each cell master. Derive W₁, W₂, A, truth and start streams from each repetition seed exactly as in checkpoint one. Each cell redraws every object. Pair sign-on and sign-off at **k = 10, 80, 100 only**, using identical matrices, truth and start in each pair. For k = 10 run only the new sign-off arm and reuse the archived signed endpoints. In sign-off descent always take the gradient at the current code and apply the same stopping rule. Initial opposite-sign and endpoint objectives remain diagnostics. This paired domain is fixed before reading outcomes.

**Curve comparison gate.** Report m/k, both cohort counts, rates and 95% Wilson intervals. Compare each signed count with the paper's independent 30-run cohort using the two-sided Fisher exact test (sum of hypergeometric probabilities no larger than the observed table's probability). A cell is discrepant at p < 0.05/15, controlling the family-wise rate across the published curve. Curve agreement requires no discrepant cells. No monotonicity or equality of empirical counts is required. The original k = 10 gate retains its own meaning. Check all new rows for finite readings, independent scalar objective/gradient discrepancies ≤ 10⁻¹⁰ and exact-zero selected preactivations. Record zero counts separately and restrict any equivalence claim accordingly.

**Paired comparison and endpoint classes.** Report both-success, on-only, off-only and neither counts at each paired k, the difference in success fractions, and the two-sided exact McNemar p = min(1, 2 P[Binomial(b+c,½) ≤ min(b,c)]). Use p < 0.05/3 for any claim of a paired recovery difference. No effect size or sign advantage is required to complete this comparison. Classify each returned endpoint in order as recovered (relative code error < 10⁻³), negative neighbourhood (angle to −x★ ≤ 20 degrees AND |‖x‖/‖x★‖ − 1/π| ≤ (1/π)/4), or elsewhere. Keep raw cosine to x★, angle to −x★, radius ratio, signed axial projection and radius deviation, plus separate orientation/radius pass flags. A zero endpoint receives no negative orientation. These neighbourhoods are finite-width diagnostics. The magnitude radius ≈ 0.66915 stays queued.

**Repeated-state acceleration and replay.** Detect only bit-identical code states within a bounded recent history. Given fixed matrices/data, a fixed deterministic mode and a pure gradient evaluation, repeating the full code repeats all later transitions and stop decisions. If a cycle contains no gradient stop, skip whole copies of that observed cycle toward update 50,000. Add skipped comparisons, zero counts and sign-flip iteration indices exactly and leave the unskipped remainder to the original loop. Record evaluated updates, skipped updates and cycle length. Before the curve run, compare accelerated signed k = 10 against all 30 original trial fields and endpoints, and replay a full 50,000-update capped path with acceleration disabled. Validate any newly observed cycle length on at least one full unaccelerated path. Run serially at below-normal priority with one .NET processor exposed.

**Rerun and diagnosis.** Keep every original row. For a curve discrepancy, numerical check failure, a paired count anomaly or failed shortcut equivalence, replay the affected cell with identical seeds, inspect endpoint geometry and stopping/accounting, and diagnose against these original gates. Add no fresh seeds and tune no step, stopping threshold or neighbourhood after observing the results. A fresh pinned replay must exercise both sign modes and a skipped cycle. Preserve endpoints and seeds for independent checking. Stop after this linear comparison.

### Linear curve readings

**The curve passed the registered comparison gate.** All 15 two-sided Fisher p values exceed 0.05/15. There were 420 new independent problems, with the original 30 signed k = 10 problems reused. The signed k = 10 diagnostic replay matched every original trial field and endpoint bit and supplied the additional geometry for that cohort. The intervals below describe the new cohort's recovery probability. Both cohorts' Wilson intervals and the unrounded p values are retained in the analysis receipt. [Actual analysis](../../../../artifacts/generative-prior/g1-linear-primary/analysis/analysis.json) · [Tracked analysis](../../receipts/generative-prior/results/g1-linear-curve/analysis.json) · [Curve CSV](../../receipts/generative-prior/results/g1-linear-curve/curve.csv)

| k | m/k | Signed successes / 30 | Paper / 30 | Wilson 95% interval | Fisher p |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 15 | 30 | 30 | 0.88649–1 | 1 |
| 15 | 10 | 30 | 30 | 0.88649–1 | 1 |
| 20 | 7.5 | 30 | 30 | 0.88649–1 | 1 |
| 25 | 6 | 30 | 30 | 0.88649–1 | 1 |
| 30 | 5 | 30 | 30 | 0.88649–1 | 1 |
| 35 | 4.28571 | 30 | 30 | 0.88649–1 | 1 |
| 40 | 3.75 | 30 | 30 | 0.88649–1 | 1 |
| 45 | 3.33333 | 30 | 30 | 0.88649–1 | 1 |
| 50 | 3 | 30 | 30 | 0.88649–1 | 1 |
| 55 | 2.72727 | 29 | 28 | 0.83330–0.99409 | 1 |
| 60 | 2.5 | 30 | 28 | 0.88649–1 | 0.49153 |
| 70 | 2.14286 | 29 | 25 | 0.83330–0.99409 | 0.19451 |
| 80 | 1.875 | 20 | 17 | 0.48780–0.80770 | 0.59591 |
| 90 | 1.66667 | 9 | 7 | 0.16665–0.47876 | 0.77102 |
| 100 | 1.5 | 1 | 4 | 0.00591–0.16670 | 0.35326 |

### Paired sign comparison

The three paired cells contain 90 matched problems. Each pair shares all five RNG streams, matrices, observations, truth and start. Independent checks confirmed the pairing. The success advantage at k = 80 was 4/30 = 13.33 percentage points, with two pairs recovering only when the sign check was disabled. No paired test crossed the registered 0.05/3 threshold.

| k | Sign-on / 30 | Sign-off / 30 | Both | On only | Off only | Neither | Exact McNemar p |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 30 | 30 | 30 | 0 | 0 | 0 | 1 |
| 80 | 20 | 16 | 14 | 6 | 2 | 8 | 0.2890625 |
| 100 | 1 | 1 | 1 | 0 | 0 | 29 | 1 |

### Negative-ray readings and stopping

All 105 failed endpoints across the 540 recorded arm endpoints were classified **elsewhere**. None passed the 20-degree negative-orientation condition or the population-radius band on its own. The failed endpoints' angles to −x★ ranged from 105.254° to 170.758°, and their relative radii ranged from 0.962260 to 1.759443. Their finite-realization negative-ray optima a★ ranged from 0.125326 to 0.608903. Every ray was nonflat. The raw control scalars, endpoint vectors and both radius deviations remain in each trial row.

| k | Mode | Recovered | Negative neighbourhood | Elsewhere |
| ---: | --- | ---: | ---: | ---: |
| 10 | on | 30 | 0 | 0 |
| 10 | off | 30 | 0 | 0 |
| 80 | on | 20 | 0 | 10 |
| 80 | off | 16 | 0 | 14 |
| 100 | on | 1 | 0 | 29 |
| 100 | off | 1 | 0 | 29 |

The cohort therefore supplies no observed endpoint in the predicted negative neighbourhood at this budget. All failed runs reached 50,000 updates, with final gradient norms between 0.0015763 and 0.184191. These readings describe budget-limited endpoints. Across all arms, 343 runs stopped on the gradient norm and 197 reached the update cap, including 92 successful recoveries. The loop accounted for 10,765,396 updates, evaluating 10,567,773 and skipping 197,623 proven repeated updates in four period-one cycles. The maximal scalar-reference discrepancy was 1.272 × 10⁻¹⁵ against the registered 10⁻¹⁰ bound, all readings were finite, and the selected-preactivation zero count was zero. All 2,250 component seeds across the 450 independent problems were distinct.

### Pins, receipts and replay

The protocol is pinned at `4116728`, the numerical implementation/settings at `8a77d0006ac8e50ca3539f3e3fd26a0097d57bb2`, and the independent audit at `968fed6`. Each primary cell has `trials.jsonl` and `summary.json` under `artifacts/generative-prior/g1-linear-primary/results/k<K>-on/`, with matching `-off/` directories at the three paired k values. The k = 10 signed diagnostic rows originated in `g1-linear-validation/results/k10-on/` and are copied into the primary directory for the combined audit. Its original first-checkpoint rows remain the baseline. Small tracked copies of all cell outputs sit under `issues/regression/receipts/generative-prior/results/g1-linear-curve/`. [Settings](../../receipts/generative-prior/settings/g1-linear-curve.json) · [Example actual paired endpoints](../../../../artifacts/generative-prior/g1-linear-primary/results/k80-off/trials.jsonl)

The signed k = 10 repetition 13 with acceleration disabled executed all 50,000 updates and agreed with the accelerated result in every numerical field and endpoint bit. A fresh build replayed repetition 13 in both modes and reproduced the primary accelerated rows, including skipped-update accounting. A full unaccelerated sign-off run also matched the accelerated numerical fields and endpoint bits. Both controls matched comparison counts, zero counts and flip indices. All observed skipped cycles had period one. The verification manifest records the checks, source pins and SHA-256 hashes of the 42 retained output files. These receipts use LF line endings and occupy about 2.2 MB, plus the manifest. [Verification manifest](../../receipts/generative-prior/results/g1-linear-curve/verification.json) · [Actual signed full control](../../../../artifacts/generative-prior/g1-linear-full-control/results/k10-on/trials.jsonl) · [Actual sign-off full control](../../../../artifacts/generative-prior/g1-linear-replay/results-full/k10-off/trials.jsonl)

The shortcut's equivalence follows from the state it caches. No RNG is consumed inside descent. Fixed matrices, data, sign mode and code bits determine the next selected mask, gradient, sign choice, stop decision and code bits. On observing a repeated code, every transition in that cycle has already run without a gradient stop. Repeating a whole cycle returns the same code. Its zero counts and one-based flip indices are translated and accumulated before the remaining original iterations run. The original stopping rule and returned endpoint are preserved.

To replay the whole linear round, run the following from the repository root with the .NET 10 SDK. It also regenerates the signed k = 10 diagnostic rows on the same fixed problems.

```powershell
./issues/regression/receipts/generative-prior/replay.ps1 -Revision 8a77d00 -RunName g1-linear-review -Experiment linear
./artifacts/evidence-contrasts/venv/Scripts/python.exe ./issues/regression/receipts/generative-prior/audit-linear.py --results artifacts/generative-prior/g1-linear-review/results --baseline issues/regression/receipts/generative-prior/results/g1-k10/trials.jsonl --output artifacts/generative-prior/g1-linear-review/analysis
```

The audit uses only Python's standard library, so another Python 3 interpreter can run the second command. A short replay exercising both sign modes and period-one acceleration is available with `-Cells 'k10-on,k10-off' -TrialIndices 13`. Add `-FullIterations` with a fresh run name to execute all 50,000 updates in each selected arm. All heavy execution ran serially at below-normal priority with `DOTNET_PROCESSOR_COUNT=1`, on .NET 10.0.5 x64.

### Round decision and doc changes owed

The linear curve agrees with the published experiment under the predeclared independent-cohort comparison. The paired data show a measured advantage of 4/30 at k = 80 and equal recovery counts at k = 10 and 100. Their exact paired tests supply insufficient evidence for a recovery advantage across these cells. The predicted negative-ray neighbourhood was unoccupied by returned endpoints in either mode. The widths, finite-realization controls and stopping budget accompany that observation. No engine primitive was needed.

The complete curve, paired outcomes, endpoint geometry, exact ray control and stopping observations are available for the next documentation batch. Depth three, magnitudes, the remaining G1 extensions and G2–G5 remain queued.

## G1 report — fixed-k measurement count

### Protocol registered before the first measurement-count readings

**Domain and pairing.** Fix k = 80, widths 250/600, zero noise and signed descent. Compare m = 80, 120, 150, 240, giving m/k = 1, 1.5, 1.875, 3. Use the same 30 problems as the preceding signed k = 80 cohort, whose cell master is 1600914390 (child 11 of `SeedTree.Derive(181204177,14)`). Derive repetition and five component seeds as before. W₁, W₂, truth and Gaussian start remain paired across every m. Draw a single row-major standard Gaussian array U with 240 rows and 600 columns per repetition from its existing measurement stream. A_m is its first m rows multiplied entrywise by 1/√m, in the existing multiplication order. The first 150 rows reproduce the preceding measurement realization. Reuse the verified m = 150 endpoints and stop/event counts.

**Binary64 and data checks.** Keep the existing C# reference, pairwise sums, masked gradient, strict sign comparison, unit step, success threshold < 10⁻³, gradient stop < 2⁻⁵², 50,000-update cap and exact repeated-state shortcut. Before new descent, verify every scaled measurement entry bit for bit against the corresponding U prefix, compare generator hashes across m, and regenerate the archived m = 150 truth/start, initial signed objectives, endpoint objective, final gradient norm and negative-ray control. Require bit-identical scalar readings and vectors for this preservation check. Record SHA-256 hashes of W₁, W₂, each standard prefix and each scaled A_m, plus mismatch counts and resolved seeds. Independent scalar objective and gradient checks retain the 10⁻¹⁰ × max(1, |reference|) tolerance. Require finite readings and record selected-preactivation zeros. A zero count qualifies the equivalence to the paper's subgradient algorithm at the selected iterates.

**Named primary contrast.** The primary question is whether this fixed finite solver recovers more of the paired problems at m = 240 than at m = 80. Use the two-sided exact McNemar test on discordant pairs at α = 0.05. Evidence for increased recovery requires a positive paired success difference and p < 0.05. This is the sole confirmatory contrast. Report all six pairwise four-cell tables, discordant repetition indices, success differences and exact p values as descriptive readings, alongside Wilson 95% marginal intervals. A nonsignificant primary contrast remains a completed experiment with that reading.

**Endpoint and stopping readings.** Preserve every returned code, seed, truth, start, relative error, sign event and stop reason. Reuse the recovered / negative neighbourhood / elsewhere definitions from the linear round, including the 20-degree orientation bound and ±25% population-radius band around 1/π. Keep the exact finite-realization negative-ray optimum and raw geometry. Report the four-outcome success pattern for each repetition and every reversal in which a lower m succeeds while a higher m fails, with both endpoint classes, errors, final gradients and stop reasons. In exact arithmetic the nested noiseless row constraints have nested feasible sets. Recovery by the specified finite solver is a separate observation, and a recovered code by itself does not establish global identifiability.

**Rerun and execution rules.** Freeze the original results, seeds, thresholds and neighbourhoods. If a data-preservation or numerical gate fails, diagnose and replay the affected repetition with the same settings before interpreting recovery. Investigate individual reversals against the stored endpoints and stop reasons, without changing the solver or extending its cap. Replay a changed-m path on a fresh pinned build and compare its numerical fields and endpoint bits. A new skipped-cycle length receives a full unaccelerated control as in the preceding round. All heavy execution remains serial, below normal priority, with one .NET processor exposed. This round ends at the fixed-k comparison.

### Fixed-k readings and primary decision

**The named primary contrast supports increased recovery by the fixed solver at m = 240 versus m = 80.** Twenty-nine paired problems changed from failure to recovery, none changed from recovery to failure, and one failed in both cells. The estimated increase is 29/30 = 96.67 percentage points. The two-sided exact McNemar p value is 3.725290298461914 × 10⁻⁹ against the preregistered α = 0.05. The experiment contains 30 underlying paired problems, with 90 new arm runs and 30 archived m = 150 endpoints. [Actual analysis](../../../../artifacts/generative-prior/g1-measurements-primary/analysis/analysis.json) · [Tracked analysis](../../receipts/generative-prior/results/g1-measurements/analysis.json)

| m | m/k | Recovered / 30 | Wilson 95% interval | Negative neighbourhood | Elsewhere | Gradient stops | Update-limit stops |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 80 | 1 | 0 | 0–0.11351 | 0 | 30 | 0 | 30 |
| 120 | 1.5 | 4 | 0.05310–0.29681 | 0 | 26 | 0 | 30 |
| 150 | 1.875 | 20 | 0.48780–0.80770 | 0 | 10 | 0 | 30 |
| 240 | 3 | 29 | 0.83330–0.99409 | 0 | 1 | 1 | 29 |

The reused m = 150 trial file agrees with its preceding receipt after line-ending normalization, and every original field and endpoint bit is preserved. All 30 regenerated m = 150 initial signed objectives, endpoint objectives, final gradient norms, ray controls, truths and starts agreed bit for bit. The four generators' hashes agreed within every problem, and every scaled Gaussian prefix entry passed its bit check. The data-check rows contain the five component seeds and the hashes of the generator matrices, truth/start vectors, standard Gaussian prefixes and scaled measurement matrices. [Data checks](../../receipts/generative-prior/results/g1-measurements/data-checks.jsonl)

### Raw paired outcomes and individual readings

All six comparisons have zero success reversals. The p values outside the named 80-versus-240 contrast are descriptive. The full lists of gained and lost repetition indices and all 30 four-bit success patterns appear in the analysis receipt.

| Lower m | Higher m | Both recover | Lower only | Higher only | Neither | Exact McNemar p |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 80 | 120 | 0 | 0 | 4 | 26 | 0.125 |
| 80 | 150 | 0 | 0 | 20 | 10 | 1.90735 × 10⁻⁶ |
| 80 | 240 | 0 | 0 | 29 | 1 | 3.72529 × 10⁻⁹ |
| 120 | 150 | 4 | 0 | 16 | 10 | 3.05176 × 10⁻⁵ |
| 120 | 240 | 4 | 0 | 25 | 1 | 5.96046 × 10⁻⁸ |
| 150 | 240 | 20 | 0 | 9 | 1 | 0.00390625 |

The m = 120 recovered repetition indices were 11, 12, 19 and 21. The nine gains from 150 to 240 were 1, 4, 6, 7, 14, 17, 20, 24 and 28. Repetition 27 failed at every m. These partitions also account for the smaller-to-larger comparisons.

Relative error can increase while a pair remains in the same success class. As a descriptive check, the audit retains every pairwise relative-error increase, including variations between already recovered codes. Increases involving a failed endpoint occurred for repetitions 1 and 24 from 80 to 120, repetition 24 from 80 to 150 and from 120 to 150, and repetition 27 from 150 to 240. The corresponding raw errors are below.

| Repetition | Error at m = 80 | Error at m = 120 | Error at m = 150 | Error at m = 240 |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 1.45199011 | 1.47888388 | 0.80638013 | 5.2948 × 10⁻¹⁶ |
| 24 | 1.04620676 | 1.06026218 | 1.07439788 | 4.3045 × 10⁻¹⁶ |
| 27 | 1.13111584 | 0.85671663 | 0.77663723 | 0.83223969 |

Repetition 27's seed was 74883976. At m = 240 its final objective was 0.35075344 and gradient norm 0.08706143 after 50,000 updates. Its angle to −x★ was 131.772° and radius ratio 1.035498, placing it in the elsewhere class. All 67 failed arm endpoints reached the update limit with nonzero gradients, ranging from 0.00380094 to 0.15750543. Of the 53 successful endpoints, 52 also reached the limit. The sole gradient stop was repetition 1 at m = 240, after 13,046 updates.

### Finite solver and identifiability

For an exact real-valued standard Gaussian row array U, the noiseless code-feasibility sets are S_m = {x | U[1:m]G(x) = U[1:m]G(x★)}. Their nesting gives S₂₄₀ ⊆ S₁₅₀ ⊆ S₁₂₀ ⊆ S₈₀, since the positive 1/√m scale cancels from each exact equation. Global identifiability asks whether S_m contains only x★.

This round measures a prescribed start and at most 50,000 unit-step signed updates on each finite binary64 objective. Increasing m changes that objective and its gradient field. The observed 0 → 4 → 20 → 29 recovery counts support improved performance of this solver in the paired cohort. A cap failure leaves global identifiability unresolved, and one recovered code supplies a recovery witness without enumerating all feasible codes. Repetition 27's nonzero residual objective and gradient at the cap document its unfinished trajectory.

### Numerical checks, pins and replay

All data-preservation and numerical gates passed. The largest independent scalar-reference discrepancy across the 120 endpoints was 1.2872 × 10⁻¹⁵ against the 10⁻¹⁰ bound, and the selected-preactivation zero count was zero. The new arms evaluated 4,463,046 updates, with 1,500,000 archived updates reused for m = 150. No exact cycle was skipped in these cells.

The protocol is pinned at `f61f907`, and the numerical experiment/settings at `c8912cc9009cbf638836a12428aed623b8d9d855`. The audit began at `3968003`; its final receipt commit adds the descriptive error-increase lists. The verification manifest pins the exact audit script by SHA-256. The original `Network.cs`, `Reference.cs` and `Recover` method body agree with the preceding pinned experiment, with their hashes recorded in that manifest. All new runs were serial, below normal priority, with `DOTNET_PROCESSOR_COUNT=1`, on .NET 10.0.5 x64.

Primary outputs are `artifacts/generative-prior/g1-measurements-primary/results/m<M>/trials.jsonl` and `summary.json`, with `data-checks.jsonl` beside the cell directories. A fresh pinned build replayed repetition 1 at both m = 80 and m = 240, reproducing every trial field, endpoint bit and data-check field. The compact tracked receipts retain all 120 endpoints and occupy about 759 kB plus their verification manifest. [Verification manifest](../../receipts/generative-prior/results/g1-measurements/verification.json) · [Actual m = 240 endpoints](../../../../artifacts/generative-prior/g1-measurements-primary/results/m240/trials.jsonl) · [Actual changed-m replay](../../../../artifacts/generative-prior/g1-measurements-replay/results/m240/trials.jsonl)

From the repository root, replay the whole fixed-k round with a fresh run name, then run the standard-library Python audit.

```powershell
./issues/regression/receipts/generative-prior/replay.ps1 -Revision c8912cc -RunName g1-measurements-review -Experiment measurements
./artifacts/evidence-contrasts/venv/Scripts/python.exe ./issues/regression/receipts/generative-prior/audit-measurements.py --results artifacts/generative-prior/g1-measurements-review/results --output artifacts/generative-prior/g1-measurements-review/analysis
```

Another Python 3 interpreter can run the audit. For the short changed-m replay, add `-Cells 'm80,m240' -TrialIndices 1` to the first command and use a fresh run name. The full-cohort audit expects all four measurement cells and 30 repetitions each.

### Round decision and doc changes owed

The paired fixed-k experiment supports the named increase in recovery from m = 80 to 240 for the frozen solver, with all original numerical gates met. There were no success reversals, and the error-increase diagnostics and unfinished repetition 27 are retained for review. The nested-row construction, paired counts, finite-solver interpretation and replay are available for the next documentation batch. This closes the fixed-k comparison.

## Remaining numerical work — active sequence

The completed execution sequence was G1, then G3, G4, G2 and G5. Each block received committed settings and criteria before readings. Source inspection, implementation controls, finite-sample findings and engine gaps are recorded in this brief. The six Lean witnesses are mathematical controls under their stated hypotheses.

## G1 report — depth and magnitude extensions

### Status

G1 is complete at the declared finite experimental scope. Its bridge, all extension cells and all full-cycle controls passed. The remaining probes proceed in order below.

### Protocol registered before extension readings

**Cells and pairing.** Use 30 repetitions per cell. Add magnitude-only observations at depth two, widths 250/600, k = 10, 80, 100, m = 150. Add depth-three widths 250/400/600 at k = 10 and 80, m = 150, for both linear and magnitude observations. Run every new cell with and without the per-iteration sign check. Reuse the original cell masters and repetition streams at the same k, so truth, start, W₁ and A remain paired across depths and objectives. Within each depth/objective all weights and observations are shared by the two sign arms. W₂ changes with the declared intermediate width at depth three. The third layer uses child 5 of `SeedTree.Derive(repetitionSeed,6)`; children 0–4 retain their original roles. Gaussian row variances remain 1/nᵢ and 1/m, and there are no biases.

**Algorithm and source distinction.** Apply the published signed gradient recurrence to f(x) = ½‖AG(x) − y‖² or f(x) = ½‖|AG(x)| − b‖², with y = AG(x★) and b = |AG(x★)|. Use strict objective improvement for negation, strict-positive ReLU masks, sgn(0) = 0 for the magnitude derivative, step 2ᵈ/d² (1 at depth two and 8/9 at depth three), norm stop < 2⁻⁵² and 50,000 updates. This tests the stated DPR algorithm on random Gaussian generators. The magnitude paper's MNIST experiments instead used Adam with opposite starts (`dpr.tex`, lines 567–569). Score ‖x − x★‖/‖x★‖ < 10⁻³ against the generating latent code. A ReLU generator has no latent odd symmetry that would justify replacing this score with min over ±x★.

**Controls and interpretation gates.** Before new cells, bridge the extended code on all 30 original signed k = 10 problems and require every original trial field and endpoint bit to agree. At each new initial selected sign and endpoint, compare scalar forward losses and fixed-mask coordinate derivatives with the existing independent reference at tolerance 10⁻¹⁰ × max(1, |reference|). Require finite readings. Count exact ReLU and magnitude preactivation zeros separately. Report recovery counts, Wilson 95% intervals, stop reasons, sign events, endpoint classes and raw returned vectors. Compare matched sign arms with exact two-sided McNemar tests. Treat seven sign comparisons as one family at α = 0.05/7. Recovery counts remain descriptive whenever that threshold is unmet. No unpublished target success rate is imposed.

**Negative-ray controls.** Starting at θ₀ = π, iterate g(θ) = arccos(((π − θ)cosθ + sinθ)/π). The linear radius is cosθ_d. The magnitude radius is 2 sinθ_d/π + ((π − 2θ_d)/π)cosθ_d, giving the depth-two checks 1/π and 0.6691522942. Use the same 20-degree angular bound and ±25% band around the objective/depth-specific radius, with recovered endpoints classified first and an elsewhere class retained. Record the exact finite-ray optimum a★ = max(0, ⟨v,y⟩/‖v‖²) for linear data and max(0, ⟨|v|,b⟩/‖v‖²) for magnitudes, where v = AG(−x★), marking a flat ray if v = 0. Population-radius and finite-ray deviations are separate readings.

**Reruns.** Preserve original rows and gates. A bridge, reference or finiteness failure triggers same-seed diagnosis and replay before interpretation. A sign-comparison significance finding is reported at the registered family threshold, with no tuning or added seeds. Validate any newly skipped deterministic cycle length on a full unaccelerated path. Replay at least one new magnitude and one new depth-three path on a fresh pinned build. Keep heavy execution serial, capped and below normal priority.

### Extension readings and decisions

All 420 extension endpoints passed the scalar-reference, finiteness and kink-count gates. There were 153 recoveries and 267 elsewhere endpoints, with no endpoint in its objective/depth-specific negative neighbourhood. No tested sign comparison crossed the registered 0.05/7 family threshold. The magnitude k = 10 point improved by 7/30 with the sign check at each depth, with eight on-only recoveries and one off-only recovery. The raw p = 0.0390625 in each case remains above the family threshold. [Actual analysis](../../../../artifacts/generative-prior/g1-extensions-primary/analysis/analysis.json) · [Tracked analysis](../../receipts/generative-prior/results/g1-extensions/analysis.json)

| Depth | Widths | Objective | k | Sign-on / 30 | Sign-off / 30 | Both | On only | Off only | Neither | Exact paired p |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 250/600 | magnitude | 10 | 26 | 19 | 18 | 8 | 1 | 3 | 0.0390625 |
| 2 | 250/600 | magnitude | 80 | 0 | 0 | 0 | 0 | 0 | 30 | 1 |
| 2 | 250/600 | magnitude | 100 | 0 | 0 | 0 | 0 | 0 | 30 | 1 |
| 3 | 250/400/600 | linear | 10 | 28 | 30 | 28 | 0 | 2 | 0 | 0.5 |
| 3 | 250/400/600 | linear | 80 | 3 | 2 | 2 | 1 | 0 | 27 | 1 |
| 3 | 250/400/600 | magnitude | 10 | 26 | 19 | 18 | 8 | 1 | 3 | 0.0390625 |
| 3 | 250/400/600 | magnitude | 80 | 0 | 0 | 0 | 0 | 0 | 30 | 1 |

The matched depth-two signed linear counts were 30, 20 and 1 at k = 10,80,100 in the preceding receipts. Their magnitude counts here were 26,0,0. These are separate forward observation models on the paired generators. The Gaussian DPR experiment does not reproduce the paper's MNIST/Adam experiment. The stronger depth-three losses also coincide with a changed intermediate architecture and the declared 8/9 step. These finite readings establish no universal depth or objective effect.

The radii used in the classifier were 0.3183098862 for depth-two linear, 0.6691522942 for depth-two magnitude, 0.4937310902 for depth-three linear and 0.7159203500 for depth-three magnitude. No endpoint passed the 20-degree negative-orientation gate. Twenty-three passed a population-radius band alone. Thus a radial match would have misclassified directions if used by itself. The full negative stationary-point location remains unobserved in this cohort.

The exact finite-ray controls gave the following descriptive k = 10 readings, shared by sign arms. They minimize the loss on the fixed negative ray, with full-space stationarity left to the actual gradient.

| Depth/objective | Population radius | Finite-ray minimum | Finite-ray mean | Finite-ray maximum |
| --- | ---: | ---: | ---: | ---: |
| 2 / magnitude | 0.66915229 | 0.48347940 | 0.68162601 | 0.84042577 |
| 3 / linear | 0.49373109 | 0.29656075 | 0.51133998 | 0.68627838 |
| 3 / magnitude | 0.71592035 | 0.44494995 | 0.72807827 | 0.96363856 |

### Numerical controls and interpretation for later probes

There were 119 gradient stops and 301 update-limit stops, including 34 successful recoveries. Thirty-five trajectories reached exactly repeated states. Their periods were 1,2,9,38,48,51,96 and 106. Every period received a full unaccelerated 50,000-update control on a fresh pinned build, agreeing in numerical fields, stop/event counts and endpoint bits. The shortcut evaluated 13,406,886 updates and accounted for 1,717,139 skipped updates. The largest scalar-reference discrepancy was 1.4433 × 10⁻¹⁵ against the 10⁻¹⁰ bound. Selected ReLU and magnitude preactivation zero counts were both zero.

The failures therefore include finite-step cycles with nonzero gradients. Their endpoints are finite-optimizer readings. The population negative-ray predictions and the papers' convergence claims require their stated distribution/isometry and step assumptions. The present outcomes supply no capture at those negative locations. G2 keeps its declared linear noisy likelihood and reads posterior mass and convergence separately from these optimization trajectories. G3 and G4 measure their own sampled condition residuals, with no causal attribution of G1 failures to an unmeasured condition.

### Extension pins, receipts and replay

The protocol is `d579c14`, numerical source/settings `0fe9045`, independent audit `b51874a`, and full-cycle driver `73c7b82`. Primary cells reside at `artifacts/generative-prior/g1-extensions-primary/results/<cell>-<on|off>/`. All 420 endpoints, summaries, the bridge rows, generated full-control settings and control outputs are tracked under `receipts/generative-prior/results/g1-extensions/`. The manifest links each copy to its actual output and SHA-256. [Verification manifest](../../receipts/generative-prior/results/g1-extensions/verification.json) · [Actual full controls](../../../../artifacts/generative-prior/g1-extension-cycle-controls/control-plan.json)

```powershell
./issues/regression/receipts/generative-prior/replay.ps1 -Revision 0fe9045 -RunName g1-extensions-review -Experiment extensions
./artifacts/evidence-contrasts/venv/Scripts/python.exe ./issues/regression/receipts/generative-prior/audit-extensions.py --results artifacts/generative-prior/g1-extensions-review/results --output artifacts/generative-prior/g1-extensions-review/analysis
./issues/regression/receipts/generative-prior/verify-extension-cycles.ps1 -RunName g1-extension-controls-review -Analysis artifacts/generative-prior/g1-extensions-review/analysis/analysis.json
```

Use fresh run names. The original signed k = 10 bridge can be selected with `-Cells bridge-k10`. All numerical runs remained capped and below normal priority. Doc changes owed include the depth/magnitude table, the four distinct radii, the absence of negative-direction capture and the distinction between finite cycles and posterior trapping.

## G3 report — layer weight conditions

### Protocol frozen before G3 execution

Execute after the G1 extension report. Use k = 2, 8, 16 and n/k = 8, 32, 128, 512, with 10 independent repetitions. Derive k masters by `SeedTree.Derive(181204180,3)`, then ten repetition seeds each. Each repetition derives separate streams for Gaussian rows, spherical rows, sparse directions, knots, the training design and pair orientations. Widths use nested raw-row prefixes rescaled by n, so repeated widths within a repetition are paired.

Compare Gaussian rows N(0,I/n), uniform-sphere rows scaled by √(k/n), and the engine's independent sparse directions with active count uniform on 1..min(3,k), uniform coordinate subsets and `DirectionalTermSpace.DrawDirection`. Use the same sparse directions without thresholds and with engine-drawn knots. Train the knot law on 512 independent standard Gaussian cases, transformed by the engine's fitted standardization, with q = 1 − 20/512 and local probability 2/3. Check the engine knot bounds against sorted scalar projections and a type-7 quantile at tolerance 10⁻¹⁰. Scale knots by the same √(k/n) factor as sparse rows.

For each matrix, measure ‖WᵀW − I‖₂ and the maximum sampled ‖Σ 1{w·x > t}1{w·y > t}wwᵀ − Q(x,y)‖₂. Use angles 0, π/12, π/4, π/2, 3π/4 and π, with four independent random orientations each. Include x = y = e₁, x = e₁ and y = e₂, and opposite e₁ explicitly. Threshold-free inputs have unit radius. Thresholded inputs use equal radii 0.5, 1 and 2 in standardized coordinates. Q uses the Gaussian angle formula, and parallel/opposite controls are I/2 and zero.

Check dense symmetric eigenvalues against the scalar Jacobi oracle on the explicit controls for every cell, to 10⁻⁹ × max(1, operator norm), with eigenpair residuals at the same tolerance. Report ten-repetition medians and ranges. The empirical width decision is the first tested n/k whose median sampled maximum is ≤ 0.1. A sampled maximum is a lower bound on the uniform supremum, so meeting that reading gives sampled evidence only.

At the widest threshold-free layer, compare the pooled e₁/e₁ second-diagonal mean with 1/2 for Gaussian/spherical rows and (A − 1)/(4(k − 1)) for sparse rows, A = min(3,k). Require agreement within four repetition-based standard errors plus 10⁻¹⁰. At k = 8 the sparse prediction is 1/14, giving the 3/7 population gap from I/2. Record this witness beside every sampled maximum. No finite-width sample is required to meet an unmeasured uniform WDC.

On a numerical-control failure replay that same seeded case and diagnose before interpreting the layer law. A lone witness deviation between four and five standard errors gets one replay on the same seeds to check reproducibility and is reported against its original gate. No new widths, seeds or tolerances are substituted. Preserve pair matrices for the widest k = 8 controls, scalar witness rows for every cell, and compact summaries. Execution remains serial, capped and below normal priority.

### G3 readings and decision

G3 is complete. All 720 scalar rows, 180 explicit k = 8 matrix controls, eigenvalue/eigenpair checks and knot checks passed. All nine population-witness comparisons passed their four-standard-error gate. A fresh build reproduced all three numerical output files byte for byte. [Tracked summary](../../receipts/generative-prior/results/g3-weights/summary.json) · [Independent audit](../../receipts/generative-prior/results/g3-weights/audit.json) · [Output/replay hashes](../../receipts/generative-prior/results/g3-weights/verification.json)

| k | Gaussian: first sampled n/k | Sphere: first sampled n/k | Sparse: widest median gap | Thresholded sparse: widest median gaps, radii 0.5 / 1 / 2 |
| ---: | ---: | ---: | ---: | --- |
| 2 | 512 | 128 | 0.23936 | 0.76097 / 0.74532 / 0.71366 |
| 8 | 512 | 512 | 0.44174 | 0.81536 / 0.81193 / 0.80173 |
| 16 | 512 | 512 | 0.47688 | 0.82946 / 0.82817 / 0.82520 |

Neither sparse law reached the sampled 0.1 target at any tested width. The widest Gaussian medians were 0.06741, 0.07567 and 0.07083; spherical medians were 0.05080, 0.06593 and 0.06510. These first tested ratios depend on the finite angle/orientation cohort and median rule. They provide no uniform WDC certificate or claim of a minimal sufficient width.

At k = 8, the threshold-free sparse e₁/e₁ second-diagonal mean was 0.0690655 (SE 0.00223364), agreeing with 1/14. Its implied mean diagonal gap from 1/2 was 0.4309345, near the analytic 3/7 population obstruction. Gaussian and spherical means were 0.5020101 (SE 0.0045040) and 0.5018654 (SE 0.0050285). At k = 2 and 16 the sparse witness means were 0.2573584 and 0.0351154, agreeing with 1/4 and 1/30. Meanwhile the widest full-Gram median gaps for sparse directions were only 0.03525, 0.07740 and 0.08456. Thus approximate second-moment isotropy coexists with a substantial activation-conditioned defect. Increasing width concentrates this declared sparse law around a different kernel. Its Gaussian WDC premise cannot be justified by width alone. Spherical rows remain a plausible dense alternative on this sampled evidence.

Thresholds further change the activation geometry and make it depend on radius. The table records the tested shifts without attributing each residual to a single mechanism. Negative knots can jointly activate opposite inputs; the threshold-free opposite-input control remains exactly zero. These are engine knot-law probes in standardized coordinates, with the training design and fitted scaling recorded in the source. They provide no claim that BASS's own inference requires Gaussian WDC.

### G3 pins, replay and doc changes owed

Protocol `0094720`; numerical source/settings `5de814a95db95da35fbceeb8ef18367a48b3a7b6`; settings SHA-256 `312054542a54dfc2fa9c1b76170cb7502088306f2fb7d8c480d81be490abbecd`. Actual primary and fresh replay directories are `artifacts/generative-prior/g3-weights-primary/` and `g3-weights-replay/`. The tracked receipts retain every scalar row and widest k = 8 matrix. Both builds completed with zero errors and 45 existing engine warnings. Execution was serial, capped and below normal priority.

```powershell
./issues/regression/receipts/generative-prior/lab-replay.ps1 -Revision 5de814a -Experiment g3 -RunName g3-review
./artifacts/evidence-contrasts/venv/Scripts/python.exe ./issues/regression/receipts/generative-prior/audit-weights.py --results artifacts/generative-prior/g3-review/results --output artifacts/generative-prior/g3-review/analysis
```

Doc changes owed: record the sparse population floor, distinguish full-Gram isotropy from activation-conditioned isotropy, and qualify the Gaussian/spherical width table as finite sampled evidence. No production change is commissioned by G3.

## G4 report — structured measurement geometry

### Protocol frozen before G4 execution

Execute after G3. Compare the resonance pilot's k = 3 waveform manifold (512 samples, 8 kHz, its declared box) with the depth-two k = 80 Gaussian generator from the first preceding k = 80 repetition. Treat its ordered 600 output coordinates as a 600-sample vector for the measurement test. This latter ordering is an experimental coordinate choice. Use five independently seeded measurement/pair repetitions from `SeedTree.Derive(181204181,5)`, each split into input, Gaussian-row and kept-bin-order streams. Draw 64 secants per manifold/repetition from independent pairs of latent codes, using the pilot's uniform box and standard Gaussian codes for G1. Test all 64 same-secant pairs plus 32 disjoint cross-secant pairs.

Use real row counts 6, 12, 24, 48, 96, 150 and 240. Compare Gaussian rows N(0,I/m), a centered length-256 periodic unit-energy Hann transform, and three centered transforms of lengths 128/256/512. Frequencies span 150–3600 Hz on a logarithmic grid at each count, equally split over the three windows. One frequency uses the geometric midpoint. Real and imaginary rows count separately. All maps within a repetition share the secants. At the pilot's kept-bin counts 3, 6, 12 and 24, use a fixed nested random subset of its 64 analysis frequencies and the full-frame unit-energy Hann window. Read both the complex linear rows (6,12,24,48 real rows) and the original nonlinear log-power differences with floor 10⁻⁴ (3,6,12,24 real observables).

Whiten the one-window, three-window and kept-complex maps through the positive eigenvalues of their row Gram matrices. Define numerical rank by λ > 10⁻⁸ λ_max and report raw rows, rank and the retained condition number. For nonlinear log power, report its Jacobian rank at the first secant's first signal and identify that local reference, alongside the underlying complex-atom rank. Independently orthonormalize the retained unwhitened row combinations by two-pass Gram–Schmidt. Require whitened row orthonormality and the normalized-secant inner products against that row-space projector to agree within 10⁻⁵. Record eigenpair residuals and the cutoff. The projector is the retained numerical row space.

For each pair of secants normalize distortion by their two original signal norms. For nonlinear log power, transform the two endpoints before differencing, then use those same original norms. Choose one nonnegative global squared scale γ per finite test set by minimizing max_i |γ⟨B s_i,B t_i⟩/(‖s_i‖‖t_i‖) − ⟨s_i,t_i⟩/(‖s_i‖‖t_i‖)|. Solve the one-dimensional feasibility intervals by bisection and check representative fits against all nonnegative intersections of the corresponding signed affine lines. Require scale-objective agreement within 10⁻⁸.

Report ε, fitted γ and its five-repetition range at each real row count, plus the matching Gaussian reading. A median ε ≤ 0.2 supplies sampled near-isometry evidence at that count. A larger sampled optimum supplies a lower bound ruling out a uniform ε ≤ 0.2 for that realization under any global scale. A passing finite search supplies no uniform certificate. Read the nonlinear log-power family separately, with its observable count and nonlinearity explicit.

All source generators, settings, secant endpoints, normalized inner-product tables and rank/whitening controls are retained or reproducible from recorded seeds. On a numerical control failure, replay the affected map and diagnose before its geometric interpretation. Keep original scales, samples, cutoff and gates, adding no frequencies or centering changes after readings. Execute serially, capped and below normal priority.

### G4 readings and decision

G4 is complete at the registered maps and samples. All 470 normalized-inner-product tables and 640 generating secants passed their controls. Maximum whitening/projector discrepancy was 1.27480 × 10⁻⁸ against 10⁻⁵, maximum scale-reference discrepancy 2.22045 × 10⁻¹⁶ against 10⁻⁸, and maximum stored eigenpair residual 3.35835 × 10⁻⁹. An independent audit checked the saved minimax objectives, feasible scale intervals and infeasibility immediately below each optimum. [Summary for every count and scale range](../../receipts/generative-prior/results/g4-audio/summary.json) · [Independent audit](../../receipts/generative-prior/results/g4-audio/audit.json)

| Map at 240 real rows | Retained row rank | Resonance median ε (range) | G1 median ε (range) |
| --- | ---: | --- | --- |
| Gaussian | 240 | 0.19214 (0.14536–0.22248) | 0.23373 (0.18439–0.27045) |
| Single window, raw | 158 | 0.99063 (0.98731–0.99636) | 0.42379 (0.39127–0.50412) |
| Single window, whitened | 158 | 0.17232 (0.06857–0.23357) | 0.29028 (0.23087–0.33946) |
| Three windows, raw | 202 | 0.74291 (0.69585–0.79581) | 0.44899 (0.36976–0.50867) |
| Three windows, whitened | 202 | 0.22146 (0.15290–0.26330) | 0.20454 (0.16819–0.23124) |

The only median decisions reaching ε ≤ 0.2 were resonance/Gaussian and resonance/single-white, both first at 240 rows. Every other registered map/manifold family failed the median target at every tested count. The near value 0.20454 remains above the original gate. Individual passing realizations and failing realizations coexist in several widest cohorts; the full ranges are retained.

| Kept-bin map at 24 bins | Real observables | Rank | Resonance median ε (range) | G1 median ε (range) |
| --- | ---: | ---: | --- | --- |
| Complex, raw | 48 | 48 | 0.83114 (0.75424–0.98484) | 0.47359 (0.41820–0.52471) |
| Complex, whitened | 48 | 48 | 0.83442 (0.76206–0.98487) | 0.41040 (0.40532–0.46611) |
| Nonlinear log power | 24 | 24, local Jacobian | 0.96339 (0.94952–0.96931) | 0.66577 (0.61443–0.73106) |

Whitening removes correlations within the retained numerical row space. Its agreement with the independent orthonormal projector control confirms that interpretation at the declared cutoff. It leaves the retained row space fixed: at the widest count the single-window and three-window ranks were 158 and 202, and their largest retained Gram condition numbers across this study were 8.55144 × 10⁷ and 4.34908 × 10⁷. The observed differences between manifolds and maps therefore call for manifold-specific geometric checks. Whitened kept-complex rows still had large residuals, and the nonlinear log-power map had a separate, substantially larger finite distortion. The local Jacobian rank describes that map at the declared reference signal; it supplies no global injectivity statement.

Each ε above is minimized over one nonnegative global squared scale for the entire finite constraint set. Whenever it exceeds 0.2, these sampled secants already rule out uniform ε ≤ 0.2 under every such scale for that realization. Passing a finite set supplies sampled evidence only. These measurements concern the declared secant products; they do not establish the complete range-restricted isometry premise on all generator differences. Treating the G1 output coordinates as ordered samples is also part of this experiment's scope.

### G4 pins, replay and doc changes owed

Protocol `8e527e3`; numerical source/settings `488f75c15ce573f69ec61663eb8fa27d903115d8`; settings SHA-256 `9e1f52cf8a56d1f5735dba7003b43354857eaebf503669da29392fae68db4bcc`. Actual primary output: `artifacts/generative-prior/g4-audio-primary/results/`. A fresh process using that isolated build reproduced all numerical output bytes in `g4-audio-replay/results/`. The tracked receipt includes all secant codes, norms and inner-product tables, with separate actual-file and LF-normalized copy hashes. [Verification](../../receipts/generative-prior/results/g4-audio/verification.json)

```powershell
./issues/regression/receipts/generative-prior/lab-replay.ps1 -Revision 488f75c -Experiment g4 -RunName g4-review
./artifacts/evidence-contrasts/venv/Scripts/python.exe ./issues/regression/receipts/generative-prior/audit-audio.py --results artifacts/generative-prior/g4-review/results --output artifacts/generative-prior/g4-review/analysis
```

The isolated build had zero errors and 45 existing engine warnings. Execution remained serial, capped and below normal priority. Doc changes owed: keep real row counts, numerical rank, global scale and nonlinear observable counts explicit; explain whitening as the retained projector; record the finite counterexamples and the two median passing cases. No transform or engine change is commissioned by G4.

## G2 report — calibration through the exact Gaussian generator

### Protocol frozen before G2 execution

Execute after G4. Use bias-free depth-two widths 250/600 and linear Gaussian measurements, with k = 2, 4, 8, m/k = 2, 4, 8 and six paired problems per cell. Derive four k masters from `SeedTree.Derive(181204182,4)`, with the fourth reserved for the conditional k = 16 expansion, six repetition seeds per k and ten branches per repetition. Branches 0–4 retain G1's weight, measurement, truth and start roles. Draw truth and start from standard normals conditioned on the declared box [−3,3]^k. Branch 5 supplies a shared standard Gaussian measurement-noise prefix. Use A_m = U[1:m]/√m and noise 0.1·z[1:m]/√m, keeping expected total noise energy 0.01 while m changes. Branches 6–8 supply nested, tempered and SMC seeds, and branch 9 is reserved for replay. Pair networks, truth, noise prefixes and initial descent codes across m and samplers.

Wrap the exact generator in the existing emulator contract, with one draw, zero emulator variance and no discrepancy. Use the engine's Gaussian variance group with an inverse-gamma prior of six degrees of freedom and deviation scale 0.1/√m. The actual noise deviation is fixed at that scale. The calibration prior is uniform on the declared box. The measured coverage is conditional on this truth-generating distribution and model specification.

Run the nested sampler with three chains, burn-in 500, decorrelation every 100 scans, batches of 250 up to 1000 retained samples per chain, R̂ target 1.05 and minimum bulk/tail ESS per chain 50. Run the same sampler on the ladder [1,0.6,0.35,0.2,0.1,0.05]. Run cut SMC with 256 particles, one outer draw, target conditional ESS 0.95, five mutation scans per stage, resampling below 0.5 and at most 1000 stages per leg. Use three independently seeded SMC populations per problem, derived from branch 8. Report R̂/ESS convergence for chain methods and ESS, lineage and between-population agreement for SMC. A population's final normalized weights alone supply no independence claim.

The inversion decision reads median relative code error ≤ 0.1 at ratios 4 and 8, plus at least five of six converged chain runs per cell. For SMC, require final ESS ≥ 0.25N in each population and coordinate-mean ranges across the three populations ≤ 0.5 pooled posterior standard deviations, in at least five of six problems. Report surviving initial lineages as a separate reading. Report 90% coordinate intervals and coverage, raw coordinate errors, noise estimates, negative-neighbourhood mass, wall time and all convergence failures. Pool unconverged chains only into labelled diagnostic readings. Report each method's frontier at these fixed budgets. Attempt k = 16 with the same six-repetition design only if every method meets its error and diagnostic gates at k = 8 and its median calibration time there is ≤ 120 seconds. A stopped expansion records this budget-defined frontier.

Compare with the G1 signed descent recurrence on the same noisy forward problem, using the same start, unit step and 50,000-update cap. Retain the unconstrained descent's returned code and recovery score, with the calibration prior box explicit. For at most 200 deterministically spaced posterior states per run, compute the cosine between the direction to truth and the negative gradient of the noisy and clean quadratic losses. Separate states above and below relative error 0.1. Record the objective-specific negative-neighbourhood mass using the original 20-degree/±25% linear-radius definition.

At k = 2, repetition 1, ratio 4, independently integrate out the variance analytically and evaluate the posterior on 201² and 401² midpoint grids over the box. Require grid posterior means to agree within 0.005 physical units and the 5%/95% quantiles within 0.02. Compare each sampler's posterior means with the 401² grid to 0.04 physical units, keeping the original statistical reading if a sampler fails. Verify emulator predictions and gradients with the scalar reference to 10⁻¹⁰ × max(1,|reference|), and retain raw posterior states locally with small tracked summaries.

On a numerical-control failure replay the affected seeded problem before interpretation. A convergence or recovery failure is a result at the registered budget and is retained. No new seeds or tuning replace a failing cohort. A fresh replay checks one k = 2 calibration per method, with timing excluded from numerical equality. Heavy execution remains capped, serial and below normal priority.

Before G2 readings, the replay selector is fixed as `-Selection k2-r1-ratio4`. It filters the declared full experiment after seed derivation, retaining all three methods and the grid control for that original problem. It changes no settings or streams. Timing fields are excluded from replay equality.

### G2 readings and budget frontier

The primary experiment completed all 54 paired problems and 162 calibrations. Independent reconstruction from all 342,222 retained states passed the mean, weighted-quantile, coverage, variance, negative-mass and SMC ESS/population-spread checks. All scalar generator controls passed. The grid means differed by at most 2.44908 × 10⁻⁶ and its 5%/95% quantiles by 0.0138460, within the original 0.005/0.02 gates. Sampler mean discrepancies from the fine grid were 0.00612350 (nested), 0.00280571 (tempered) and 0.00104226 (SMC), all below 0.04. The nested run's grid agreement coexists with a failed convergence verdict. [Full cell summary](../../receipts/generative-prior/results/g2-calibration/summary.json) · [State audit and each failed diagnostic](../../receipts/generative-prior/results/g2-calibration/audit.json)

| k | m/k | Signed descent median error | Nested median / diagnostic passes | Tempered median / diagnostic passes | SMC median / diagnostic passes | SMC 90% coordinate coverage |
| ---: | ---: | ---: | --- | --- | --- | ---: |
| 2 | 4 | 0.08840 | 0.03717 / 0 of 6 | 0.03065 / 6 of 6 | 0.03404 / 6 of 6 | 12/12 |
| 2 | 8 | 0.02632 | 0.02670 / 0 of 6 | 0.02618 / 6 of 6 | 0.02654 / 6 of 6 | 11/12 |
| 4 | 4 | 0.04901 | 0.05703 / 0 of 6 | 0.04296 / 0 of 6 | 0.04756 / 6 of 6 | 20/24 |
| 4 | 8 | 0.03755 | 0.04438 / 0 of 6 | 0.03833 / 0 of 6 | 0.03927 / 6 of 6 | 21/24 |
| 8 | 4 | 0.04971 | 0.09872 / 0 of 6 | 0.06792 / 0 of 6 | 0.05265 / 6 of 6 | 43/48 |
| 8 | 8 | 0.02742 | 0.04482 / 0 of 6 | 0.03478 / 0 of 6 | 0.02549 / 6 of 6 | 43/48 |

At both decision ratios, SMC met the median-error and five-of-six diagnostic gates through k = 8. Tempering met them at k = 2 and failed the diagnostic requirement at k = 4 and 8. Ordinary nested sampling failed the diagnostic requirement already at k = 2. These are frontiers at the declared budgets; they establish no impossibility at larger budgets or an algorithm's universal dimension limit. All-method qualification at k = 8 failed, so the registered conditional k = 16 expansion was skipped. It was the only unrun G2 dimension.

The ratio-2 cells are retained too. SMC's median errors were 0.10690, 0.08824 and 0.11469 at k = 2,4,8; its diagnostic passes were 6,6,4 of 6. The two k = 8 failures had between-population spreads 0.51621 and 0.53770, above 0.5. Every problem's minimum final particle ESS exceeded 128, comfortably above the 64-particle gate. At k = 8 the minimum surviving initial-family count ranged from 2 to 11 out of 256; at lower dimensions it ranged from 8 to 89. Adequate final weights and mean agreement do not erase genealogy loss.

Nested and tempered columns include explicitly labelled diagnostic pools whenever the original engine verdict failed. Their raw interval coverage is retained beside diagnostic-qualified coverage in the summary; in cells with no qualified runs that latter reading is null. The verdict gates coordinates, log variance and log likelihood. The saved `MaximumRHat`/`MinimumEssPerChain` extrema cover coordinates and log variance; they are narrower summaries than the verdict. Sampling used thinning 1, with the separate decorrelation move every 100 scans. Particle ESS and chain ESS have separate fields and units.

The coverage fractions are small-cohort, coordinate-wise readings conditional on the Gaussian truth law restricted to the box, this exact generator, the fixed noise deviations and the declared inference prior. The inference prior on the code is uniform on the box and the variance has an inverse-gamma prior; the truth/noise cohort is not drawn from that full prior predictive distribution. No general 90% coverage guarantee follows. The audit retains posterior noise estimates against the actual 0.1/√m deviation for every cell. At decision ratios 4 and 8, SMC's median posterior RMS deviations ranged from 0.9943 to 1.0708 times their actual deviations.

| k | Descent median seconds, ratios 4 / 8 | Nested | Tempered | SMC, three populations |
| ---: | --- | --- | --- | --- |
| 2 | 0.115 / 0.105 | 0.253 / 0.291 | 0.863 / 0.826 | 9.916 / 11.947 |
| 4 | 0.079 / 0.055 | 0.190 / 0.197 | 0.941 / 0.937 | 10.061 / 11.180 |
| 8 | 0.195 / 0.093 | 0.227 / 0.285 | 1.197 / 1.326 | 13.632 / 17.497 |

Times describe these serial, one-processor, below-normal-priority runs, including their finite stopping rules and runtime warm-up. They exclude grid construction and receipt analysis. The fast chain timings often end at the retained-sample budget with failed diagnostics. Descent used the same noisy data and start, with an unconstrained optimization path; all returned codes nevertheless stayed within the calibration prior box. At ratio 2 its median errors were 0.18392, 0.34571 and 0.31285, and the k = 8 median time was 9.825 seconds. The posterior mean and the finite loss minimizer answer different inferential questions.

### G2 negative neighbourhood, alignment and proposed primitive

Only the k = 2, ratio-2, repetition-3 problem had retained mass in the declared negative neighbourhood: 0.0040000 for tempered sampling and 0.00141031 for SMC. Both passed their diagnostics. Every ordinary nested reading and every other problem/method had zero recorded mass there. This cohort supplies no evidence that ordinary chains were trapped at that neighbourhood or that tempering cleared such a trap. Its convergence failures require their own diagnosis.

The 200-state deterministic diagnostic subsets per calibration retained 32,400 states in total. Their unweighted cosine summaries are separated at relative code error 0.1; they are state diagnostics, including unequally weighted SMC particles and unconverged chain states. Above that threshold the median clean cosine across the k/method groups ranged from 0.67771 to 0.85186 and the median noisy cosine from 0.55826 to 0.78676. For k = 8 SMC, the above-threshold subset had 1,326 states with medians 0.67771 (clean) and 0.56163 (noisy); below the threshold its 2,274 states had medians 0.89545 and 0.66092. The saved rows carry individual weights and both cosines. These local readings indicate useful first-order information on the sampled states; they prove no global inward-gradient or posterior-convergence claim.

The existing `ICalibrationEmulator` accepted the generator without an engine change. It exposes predictions and noise, with no derivative contract. The calibration code uses adaptive random walks, decorrelation and the configured tempering/SMC machinery; the general sampler library already has `NewtonMetropolis`, which requests both a gradient and a Hessian. A proposed general primitive is an optional prediction-adjoint/first-derivative contract, propagated through measurement likelihood and parameter coordinates into a log-target score, with an explicitly Metropolis-corrected first-order move and defined bounded-domain behaviour. Qualification should compare equal computational budgets against the current kernels on independent targets and preserve the existing diagnostic gates. The alignment and chain failures warrant that investigation; this experiment has not tested its benefit. No engine implementation is authorized here.

### G2 pins, receipts and replay

G2 is complete. Protocol `3b3889f`, replay selection `ae10b7e`; numerical source/settings `3296cb48ac848ab14529cd1dd22ddc0f8f604966`, settings SHA-256 `4914fd80ebb43c8717aae3ba66c6d2d43e7f84630b213df5199eef942950b7b7`. Actual primary output is `artifacts/generative-prior/g2-calibration-primary/results/`. Full posterior and alignment states remain there; their hashes, problem seeds, raw descent endpoints, compact calibrations and grid/summary/audit receipts are tracked. The fixed fresh-build replay at `artifacts/generative-prior/g2-calibration-replay/` agreed exactly after excluding timing: all three calibration rows, the problem row, all 6,018 posterior states, all 600 alignment states and the grid control. Replay group decisions from a single problem have no cohort interpretation. Both isolated builds had zero errors and 45 existing engine warnings. [Primary verification](../../receipts/generative-prior/results/g2-calibration/verification.json) · [Replay receipts](../../receipts/generative-prior/results/g2-calibration/replay/verification.json)

```powershell
./issues/regression/receipts/generative-prior/lab-replay.ps1 -Revision 3296cb4 -Experiment g2 -RunName g2-review
./issues/regression/receipts/generative-prior/lab-replay.ps1 -Revision 3296cb4 -Experiment g2 -RunName g2-review-replay -Selection k2-r1-ratio4
./artifacts/evidence-contrasts/venv/Scripts/python.exe ./issues/regression/receipts/generative-prior/audit-calibration.py --results artifacts/generative-prior/g2-review/results --replay artifacts/generative-prior/g2-review-replay/results --output artifacts/generative-prior/g2-review/analysis
```

Doc changes owed: record the budget-defined sampler frontiers, diagnostic qualification of posterior summaries, genealogy loss, the measured gradient alignments and the absence of observed impostor trapping. The derivative seam remains a proposed primitive. The noise/prior and coordinate assumptions accompany every coverage statement.

## G5 report — noise trend in existing reference receipts

### Protocol frozen before G5 readings

Execute after G2, using only the existing `artifacts/latent-recovery/results/stage1/s1-k<K>-replicates.csv` for k = 1..5 and `stage2/s2-k<K>-replicates.csv` for k = 2,3. Select emulator = oracle, method = nested and discrepancy = 0. Hash every original file and save the selected original fields and source row indices as a portable input receipt. No new inversions, synthesizer draws or noise realizations are added.

Use per-repetition squared reconstruction RMSE at clean held-back log-power bins as the response. Report its median and quartiles against k/m and waveform noise variance in every cell, retaining unconverged original runs and reporting their fraction. Stage 1 has three waveform noise levels and four measurement ratios. Stage 2 has one waveform noise level, so it supplies no noise-level slope.

Within each stage fit log MSE with a separate intercept per k, a log(k/m) slope and, where identifiable, a log waveform-noise-variance slope. Estimate 95% percentile intervals from 2000 cluster bootstrap resamples, drawing whole repetition clusters within k so paired ratios and noise levels remain together. Stage seeds are children of `SeedTree.Derive(181204183,2)`; bootstrap seeds derive structurally from each stage seed, using the pinned SplitMix64/xoshiro256++ streams. Verify those integer streams against the C# implementation.

A positive 95% lower bound for the log(k/m) slope supports a decreasing-error trend as m grows in this saved cohort. An interval containing one is compatible with a proportional k/m shape within the fitted model. These are descriptive, specification-dependent readings. Also fit the converged-only rows as a sensitivity check, with their counts explicit.

The stored RMSEs are rounded to four decimal places. Retain lower and upper squared-error bounds from max(0,RMSE − 0.00005) and RMSE + 0.00005, and repeat the point slope fit at those bounds. Zero rounded RMSE is left-censored and is omitted from the primary log fit, with its count recorded; the upper-bound sensitivity retains it. Do not substitute an arbitrary logarithm floor.

S21 gives a fixed rank-k projector's expected k/n energy fraction under isotropic noise in the same coordinates. Here noise was added to the waveform and reconstruction is scored in nonlinear log-power coordinates, with a fitted variance and a finite posterior approximation. Interpret a measured trend in those coordinates, retaining that distinction from the projection identity. A data-shape, stream or arithmetic-control failure triggers same-input diagnosis and replay. Preserve all gates and original rows, and report rounding sensitivity or missing receipts rather than generating replacement data.

### G5 readings and decision

G5 is complete using the original seven input files only. Data-shape gates passed for 1,800 stage-1 and 160 stage-2 rows, with no duplicate paired keys and no zero rounded RMSEs. Integer controls matched ten seed children, ten 64-bit outputs and ten bounded integers for each of three positive/negative seed cases against the isolated C# assembly. An analytic response with known slopes 1.25 and 0.75 reproduced both within 5 × 10⁻¹⁶, and the single-noise control correctly left its noise slope unidentified. A pinned portable-row replay reproduced every stage result and every portable CSV byte. [Full cells and source-file hashes](../../receipts/generative-prior/results/g5-noise/summary.json) · [Arithmetic/stream verification](../../receipts/generative-prior/results/g5-noise/stream-verification.json) · [Independent cell and replay audit](../../receipts/generative-prior/results/g5-noise/audit.json)

| Saved cohort | Rows / wholly converged rows | Waveform noise variances | m/k values | log(k/m) slope, bootstrap 95% interval | log noise-variance slope, bootstrap 95% interval |
| --- | --- | --- | --- | --- | --- |
| Stage 1, k = 1..5 | 1,800 / 936 (52%) | 10⁻⁶, 10⁻⁵, 10⁻⁴ | 1,2,4,8 | 2.54197 [2.39159, 2.68489] | 0.45027 [0.41519, 0.48345] |
| Stage 2, k = 2,3 | 160 / 12 (7.5%) | 10⁻³ | 2,4,8,16 | 1.65977 [1.25489, 2.05488] | unidentified |

Both log(k/m) intervals have positive lower bounds, meeting the registered descriptive trend criterion. Neither interval includes one, so a proportional k/m shape is outside these fitted intervals. Stage 1 also has a positive fitted noise-variance slope. The model contains one intercept per k and common slopes, with no interaction terms. These are readings of finite stored reconstructions under that specification and bootstrap design, conditional on its repetition clusters.

The stage-1 converged-only point slopes were 2.60991 for log(k/m) and 0.43415 for noise variance, retaining the same trend direction. The 12 converged stage-2 rows all have m/k = 2, so its converged-only ratio slope is unidentified after the k intercepts. Stage 2's primary trend therefore remains a diagnostic reading dominated by unconverged saved inversions. The analysis generated no replacements for them and supplies no additional qualification of those posteriors.

| Example cell: k = 3 | Waveform variance | m/k | Median squared held-back log-power RMSE |
| --- | ---: | ---: | ---: |
| Stage 1 | 10⁻⁵ | 1 | 0.17413153 |
| Stage 1 | 10⁻⁵ | 2 | 0.03542068 |
| Stage 1 | 10⁻⁵ | 4 | 0.00561065 |
| Stage 1 | 10⁻⁵ | 8 | 0.00085328 |
| Stage 2 | 10⁻³ | 2 | 3.49145992 |
| Stage 2 | 10⁻³ | 4 | 2.64490346 |
| Stage 2 | 10⁻³ | 8 | 0.48613730 |
| Stage 2 | 10⁻³ | 16 | 0.36088969 |

All 68 cells retain quartiles, convergence counts and lower/upper rounded-error median bounds. Point fits at the two rounding bounds moved the stage-1 ratio slope over 2.53877–2.54521 and the noise slope over 0.44945–0.45111; the stage-2 ratio slope moved over 1.65946–1.66007. These two simultaneous bound fits are sensitivity checks, with no claim that they bound every possible coefficient over independently rounded rows. No censoring rule was activated because no selected RMSE rounded to zero.

S21 concerns a fixed orthogonal projector and isotropic noise in those same coordinates. The saved pilot adds waveform noise, observes a nonlinear log power with a floor, fits a noise variance and approximates a posterior, then scores clean held-back log powers. Its fitted exponents provide evidence of decreasing prediction error with more measurements in these cohorts. They do not verify S21's k/n identity, a benign generator landscape or an asymptotic noise-reduction rate.

### G5 pins, replay and doc changes owed

Protocol/settings `9512c90`, input-shape controls `711d34b`, source-hash field `1f5fc7d`, stream controls `ef1b489` and analytic slope control `db5e882`. Numerical source/settings pin: `1667f0b`. Primary settings SHA-256 `6ff30860d58acf3f11ccd85e45ad574dada3b0768e989c6553ded39d356288e6`; analysis-source SHA-256 `9b9c7d78fc364ca29b253b908715e977a6a97974590afa9f1cbd749560d06502`; portable CSV SHA-256 `21f0dfc25eafe72cf151ad4128afe95e17ffcde9ca316c6a7bee05035d58d05f`. All three hashes also matched the replay. Actual outputs are `artifacts/generative-prior/g5-noise-primary/results/` and `g5-noise-replay/results/`; stream controls are in `g5-noise-controls/`. [Primary verification](../../receipts/generative-prior/results/g5-noise/verification.json) · [Replay verification](../../receipts/generative-prior/results/g5-noise/replay/verification.json)

The initial bootstrap launch did not explicitly verify its process priority. The pinned replay used one worker and verified inherited BelowNormal priority inside Python before analysis; its complete numerical agreement is recorded. This was an execution-control correction, with no data, gate or statistical setting changed.

```powershell
./issues/regression/receipts/generative-prior/verify-noise-streams.ps1 -RngAssembly artifacts/generative-prior/g2-calibration-primary/build/bin/GenerativePriorLab/Release/net10.0/Maths.Rng.dll -Python ./artifacts/evidence-contrasts/venv/Scripts/python.exe -Output artifacts/generative-prior/g5-stream-review
./issues/regression/receipts/generative-prior/noise-replay.ps1 -Revision 1667f0b -RunName g5-review -Python ./artifacts/evidence-contrasts/venv/Scripts/python.exe -Rows issues/regression/receipts/generative-prior/results/g5-noise/oracle-rows.csv
```

Doc changes owed: describe the coordinate-specific trend and its non-unit fitted slopes, rounding sensitivity, the stage-1 converged-only point check and the unidentified stage-2 sensitivity. Keep the distinction from fixed-projector noise energy explicit. No new inversion or primitive implementation was needed for G5.

## Completion and remaining limits

All five probes have their declared decisions, raw/compact receipts, source pins and controls above. No required numerical block remains open. The conditional k = 16 expansion was skipped as specified. The negative stationary neighbourhood was not observed as an optimizer destination, and no ordinary-chain impostor trap was established. Finite pair searches leave uniform WDC/RRIC qualification open; sampler frontiers and coverage remain budget/model-specific; G5's stage-2 convergence limitation remains unresolved by its saved data. The proposed derivative seam is future work requiring its own qualification. Production `src/`, `tests/` and `docs/` were not edited. The Exposition below supplies self-contained segments for later document work.

## Exposition

In noiseless recovery through a bias-free Gaussian ReLU generator of widths 10 → 250 → 600, 150 Gaussian measurements sufficed to recover the generating code in all 30 independent repetitions of the signed descent experiment. At each iteration the method evaluated the quadratic measurement loss at both signs and took a unit gradient step from the lower-loss sign. Seventeen starts changed sign, all on the first comparison. No selected iterate had an exactly zero layer preactivation, so the practical gradient agreed with the paper's subgradient algorithm throughout these runs.

Recovery and numerical termination are separate measurements. All relative code errors were below 9.24 × 10⁻¹⁶. In one run the signed update rounded back to the same binary64 vector while the gradient norm remained above the spacing at one. That run satisfied the experimental recovery criterion at the 50,000-update limit.

Across the published depth-two linear-recovery grid, with widths 250/600 and m = 150, signed descent recovered 30/30 codes through k = 50, 29/30 at k = 55, 30/30 at k = 60, 29/30 at k = 70, 20/30 at k = 80, 9/30 at k = 90 and 1/30 at k = 100. All 15 comparisons with the paper's independent 30-run cohorts passed the preregistered family-wise Fisher-test gate.

With networks, truths and starts paired, disabling the sign check changed recovery from 20/30 to 16/30 at k = 80. Six pairs recovered only with the check and two only without it. Both modes recovered 30/30 at k = 10 and the same single code at k = 100. None of the three exact paired tests reached the registered significance threshold.

For each finite generator and measurement matrix, positive homogeneity gives the exact negative-ray optimum a★ = max(0, ⟨AG(−x★), AG(x★)⟩/‖AG(−x★)‖²) when the denominator is positive. This control distinguishes a realization's radial shift from a returned code pointing in another direction. Every failed endpoint in these runs reached the update limit and lay outside the declared negative neighbourhood around −x★/π.

At fixed k = 80 and widths 250/600, nested Gaussian measurement rows gave 0, 4, 20 and 29 successful code recoveries among the same 30 paired problems at m = 80, 120, 150 and 240. The prescribed signed unit-step solver retained its 50,000-update cap and 10⁻³ relative-code-error threshold. The named 240-versus-80 contrast had 29 gains, zero losses and one joint failure, with two-sided exact McNemar p = 3.72529 × 10⁻⁹. Every individual success pattern was monotone across the four m values.

Exact noiseless feasible sets shrink under nested measurement rows. The finite optimization trajectory also depends on how those rows reshape the objective. The fixed-budget recovery counts measure that solver's performance, while global identifiability concerns uniqueness over the entire feasible set. Relative code errors increased on some individual failed trajectories despite the monotone binary success patterns.

The depth and magnitude extensions used the published sign-comparison recurrence on paired Gaussian generators. At depth two, magnitude recovery counted 26/30 at k = 10 and 0/30 at k = 80 and 100. At depth three, with widths 250/400/600, linear recovery counted 28/30 at k = 10 and 3/30 at k = 80; magnitude recovery counted 26/30 and 0/30. None of seven paired sign comparisons reached the registered family threshold. Every endpoint lay outside its objective-specific negative neighbourhood. Exact repeated-state controls established that some failed trajectories were finite-step cycles with nonzero gradients. Their locations do not identify posterior modes or validate the population negative-ray stationary-point predictions.

Width affects sampling fluctuations around a layer law, while the law also fixes the population activation kernel. Under the engine's threshold-free sparse directions, with active count uniform from one to three at k = 8, the e₁/e₁ second-diagonal expectation is 1/14. The empirical mean at n/k = 512 was 0.06907 with repetition-based standard error 0.00223. Its separation from the Gaussian value 1/2 agrees with the analytic 3/7 obstruction, even though the full-Gram median defect was only 0.07740. Gaussian and spherical rows met the sampled 0.1 activation-kernel target at sufficiently wide tested layers; neither sparse law did. These finite width decisions apply to the sampled angles and orientations, with the analytic sparse witness providing the population distinction.

Structured audio measurements showed a separate row-space effect. At 240 real rows on the resonance manifold, the single-window median finite distortion fell from 0.99063 to 0.17232 after whitening, with retained numerical rank 158. The whitened products agreed with an independent orthonormal basis for that retained row space. At 24 kept frequencies, whitened complex measurements still had median distortion 0.83442, and nonlinear log power had median distortion 0.96339 across 24 real observables. These optima allow one global scale for each entire sampled constraint set. A finite optimum above a proposed uniform tolerance supplies a counterexample within that set; a finite optimum below it supplies sampled evidence whose extension to all generator secants remains unproved.

The exact Gaussian generator also served the existing Bayesian calibration contract. For k = 2,4,8, six problems at each m/k = 2,4,8 shared Gaussian measurement/noise prefixes and used a uniform code prior on [−3,3]ᵏ, an inferred scalar noise variance and no emulator uncertainty or discrepancy. At ratios 4 and 8, three independently seeded 256-particle SMC populations passed the registered mean-error and between-population diagnostic gates through k = 8. Tempered chains passed at k = 2; ordinary nested chains failed their convergence gates already there. All methods' medians could be close to truth while their diagnostic verdicts differed. These outcomes describe fixed computational budgets and this model, with no universal sampler dimension limit or coverage guarantee.

The low-dimensional independent grid agreed with every sampler's mean within 0.00613 physical units. Nevertheless its ordinary nested run failed convergence, illustrating why agreement of one posterior functional does not certify the entire retained sample. The largest recorded mass near the declared negative neighbourhood was 0.004 in one tempered k = 2, ratio-2 run, with 0.00141 for SMC on the same problem. No ordinary nested run recorded mass there. These data support no attribution of chain failure to trapping at the predicted impostor. At retained diagnostic states, clean and noisy loss gradients usually pointed toward truth; the measured alignment motivates qualification of a general derivative-aware proposal seam, whose actual sampling benefit remains untested.

Saved resonance-pilot reconstructions showed decreasing clean held-back log-power error as measurement count grew. A log-squared-error model with one intercept per code dimension gave log(k/m) slopes 2.54197 [2.39159, 2.68489] in the first cohort and 1.65977 [1.25489, 2.05488] in the second, using 2,000 whole-repetition cluster bootstrap resamples. Both intervals excluded zero and one. The first cohort's converged-only point slope was 2.60991; the second had only 12 wholly converged rows, all at one ratio, leaving that sensitivity unidentified. The trend concerns nonlinear log-power prediction under the saved finite inversions. It supplies no verification of the k/n noise-energy identity for a fixed projector in isotropic-noise coordinates.
