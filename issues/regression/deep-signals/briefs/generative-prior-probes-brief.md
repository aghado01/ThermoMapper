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
- *Decisions:* how wide a layer must be for a given ε; whether a layer of directional hinges under its prior is near enough to angle-preserving for an inversion through it to inherit the benign landscape, and where its residual floors.

**G4. The isometry constant for structured views.**
- On the pilot's resonance manifold and on G1's generator range, sample secant pairs and estimate the smallest ε for which each measurement map satisfies the isometry condition, after its best global scale, against m.
- Maps: Gaussian rows; rows of a linear short-time transform at one window length, as real and imaginary parts; the same across three window lengths; each of those whitened by its atoms' Gram; the pilot's kept-bin subsets.
- *Decisions:* whether structured views are near-isometric on these manifolds at the m where recovery worked; whether whitening lowers ε; how far they sit from Gaussian rows.

**G5. The noise reduction, from receipts already held.** In the pilot's reference inversions, reconstruction error at the held-back bins against the noise level and k/m. No new runs. *Decision:* whether the k/m rate shows at these sizes.

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
