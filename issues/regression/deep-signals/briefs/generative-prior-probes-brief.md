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
- Without the sign check, where descent from random starts ends: near x₀, or near −ρ_d·x₀. Read ρ₂ against 1/π.
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
