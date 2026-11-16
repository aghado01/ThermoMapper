# Generative-prior probes: the numerical chip's third brief

Date: 2026-11-16. Status: **posed at the user's request on 2026-11-16** ("Yes, pose both", to Lean and numerical briefs on the compressed-sensing results). Drafted by the coordinating session. Its sibling is the Lean chip's [generative-prior statements brief](generative-prior-statements-brief.md). The direction it serves is in the [approach brief](approach-brief.md#the-compressed-sensing-guarantees-used-by-analogy-not-yet-tested).

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

## Exposition

In noiseless recovery through a bias-free Gaussian ReLU generator of widths 10 → 250 → 600, 150 Gaussian measurements sufficed to recover the generating code in all 30 independent repetitions of the signed descent experiment. At each iteration the method evaluated the quadratic measurement loss at both signs and took a unit gradient step from the lower-loss sign. Seventeen starts changed sign, all on the first comparison. No selected iterate had an exactly zero layer preactivation, so the practical gradient agreed with the paper's subgradient algorithm throughout these runs.

Recovery and numerical termination are separate measurements. All relative code errors were below 9.24 × 10⁻¹⁶. In one run the signed update rounded back to the same binary64 vector while the gradient norm remained above the spacing at one. That run satisfied the experimental recovery criterion at the 50,000-update limit.

Across the published depth-two linear-recovery grid, with widths 250/600 and m = 150, signed descent recovered 30/30 codes through k = 50, 29/30 at k = 55, 30/30 at k = 60, 29/30 at k = 70, 20/30 at k = 80, 9/30 at k = 90 and 1/30 at k = 100. All 15 comparisons with the paper's independent 30-run cohorts passed the preregistered family-wise Fisher-test gate.

With networks, truths and starts paired, disabling the sign check changed recovery from 20/30 to 16/30 at k = 80. Six pairs recovered only with the check and two only without it. Both modes recovered 30/30 at k = 10 and the same single code at k = 100. None of the three exact paired tests reached the registered significance threshold.

For each finite generator and measurement matrix, positive homogeneity gives the exact negative-ray optimum a★ = max(0, ⟨AG(−x★), AG(x★)⟩/‖AG(−x★)‖²) when the denominator is positive. This control distinguishes a realization's radial shift from a returned code pointing in another direction. Every failed endpoint in these runs reached the update limit and lay outside the declared negative neighbourhood around −x★/π.
