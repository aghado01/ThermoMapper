# Lineage gaps: register

Where the Francom lineage and the literature around it stop short, and what closing each gap would take with this engine. The engine is built for diverse applications, so gaps are stated by problem class. Nothing here is scheduled: an entry starts on the user's go-ahead. Lineage findings that land keep their record in the port hub's [register](../bass-full-port-brief.md#what-the-lineage-taught); extension work keeps its entry in the [extensions brief](extensions-brief.md).

**How to read it.**
- *Literature status* comes from web searches on 2026-11-01. Unless a paper is marked as read, it is known from its abstract or a search snippet. "Not found" means those searches did not turn it up, not that it does not exist. Most references were procured the same day and await reading (see the [reading list](#reading-list)).
- *Statements* are posed in [`lean/Protolemmata/BASS/`](../../../lean/Protolemmata/BASS/), and [`lineage-gap-checks.py`](../../../lean/Protolemmata/BASS/lineage-gap-checks.py) checks them numerically. Nothing is proved yet.
- *Measurements* cited below are this engine's own, with their reports linked.

| Gap | What the literature leaves | What this engine holds | Kind of work | Statements |
|---|---|---|---|---|
| [G1](#g1-does-the-sampler-reach-the-posterior) | Whether Bayesian MARS's sampler reaches its posterior | F-R1's failures, the flat histogram's a1, enumerated kernels | Measurement, then a proof-walk | [SC-0, SC-2, SC-4](../../../lean/Protolemmata/BASS/structure-chain-lemmas.md) |
| [G2](#g2-posterior-contraction-for-bass) | Any contraction result for Bayesian MARS | The BGP proof chain, in Lean | Proof-walk, then new lemmas | [The walk](../../../lean/Protolemmata/BASS/contraction-walk.md) |
| [G3](#g3-sobol-indices-and-active-subspaces-have-different-standing) | What a posterior's active subspace is worth | The moment service, term derivatives | A lemma, a recovery check, a caveat on E8 | [SG-2, SG-5](../../../lean/Protolemmata/BASS/sensitivity-generator-lemmas.md) |
| [G4](#g4-one-moment-service-across-dependent-inputs) | Closed forms under dependent inputs | `BassMoments`, `IIntervalMoments` | A new measure | [SG-3](../../../lean/Protolemmata/BASS/sensitivity-generator-lemmas.md) |
| [G5](#g5-sensitivity-as-a-markov-generator) | One frame for Sobol, DGSM and active subspaces; a smoothness-matched PCE prior | The mixing layer, graph Laplacians, F-X3's per-part contract | Formalization, a readout, a prior option | [SG-1, SG-4](../../../lean/Protolemmata/BASS/sensitivity-generator-lemmas.md) |
| [G6](#g6-orthogonal-discrepancy-for-cut-calibration) | Cut calibration with a discrepancy | D09's cut posterior, stage 7's D·v discrepancy | A derivation, then a calibration option | Not posed yet |

Two other things sit beside these: a [hypothesis](#h1-the-steep-band-as-a-transition-in-sample-size) to falsify, and [what was checked and is not open](#checked-and-not-open).

---

## G1: Does the sampler reach the posterior?

**What the literature leaves.** The lineage judges fits by prediction and coverage. BART has mixing theory: an exponential lower bound for a simplified model (Ronen et al., arXiv 2210.09352), upper bounds for Bayesian CART (Kim & Roćková, arXiv 2306.00126). Nothing was found for Bayesian MARS.

**What this engine holds.**
- All 144 of F-R1's untempered fits fail the predictive R̂ ≤ 1.05 gate, at medians 2.14 and 1.74 ([report](deterministic-profile-brief.md)).
- On S10's Friedman fixture (N = 500) the flat histogram finds a1, 4 × x1·x2 + 2 × x3 + x4 + x5, at 66–68%, with ln Ẑ(1) of −841.36 ± 0.17, −841.55 ± 0.07 and −841.14 ± 0.32 over three seeds. That is about 45 nats above the reference arms, whose basins hold at most 1.3% ([docs](../../../docs/Maths/BayesianRegression.md); [report](flat-histogram-warmstart-brief.md)). These are exploration runs, one per seed, on one fixture.
- F-M2 found local modes tens of nats apart, every term load-bearing, and births lost to an Occam cost near ½ ln N ([measurements](local-move-measurements-brief.md)).

**How to close it.**
1. *Exact gaps on enumerated kernels.* At fixed τ the structure chain is a finite reversible chain (SC-0), and `BassKernelEnumeration` already builds its kernels. Compute 1 − λ₂ per move mix and SC-2's lower bound over the basins, then ask whether R̂ and the flat histogram's claims agree with them. Nothing in samplers computes a gap today. An N sweep reads how the gap scales. The fixtures are small, so this shows scaling rather than Friedman's barrier.
2. *A proof-walk of Ronen et al. onto BASS* (SC-4). The barrier argument suggests linear growth in N, which is what an exponential bound needs.
3. *Persistence* (Woodard, Schmidler & Huber 2009) as the reading of the steep ladder band, measured from a flat-histogram run's ln g. Their result and Bhatnagar & Randall's mean-field Potts q ≥ 3 sit near the [Potts toy](../../../docs/Maths/Samplers.md)'s q = 4 failure, but the toy's fitted N³ at q = 3 looks polynomial where the theorem says torpid. Resolve that first.

**Would land in:** the mixing layer shared across engines; F-M2's record. *Unlikely to help:* locally balanced proposals, an idea that went in here and came out. F-M2 measured every term load-bearing, and its use-weighted births gained only ×1.07, so there is little sign of an attractive neighbour to scan for. That is a reading of those measurements, not a test of these proposals.

---

## G2: Posterior contraction for BASS

**What the literature leaves.** Nothing found for Bayesian MARS. *MARS via LASSO* (Ki, Fang & Guntuboyina, *AoS* 2024) gives a rate, logarithmic in dimension, for least squares under a variation constraint.

**What this engine holds.** [PosteriorContraction](../../../docs/Maths/PosteriorContraction.md): Gaussian KL, tests and covers in Lean, with the general theorem cited as axioms.

**How to close it.** [The walk](../../../lean/Protolemmata/BASS/contraction-walk.md), with its reading of each step. The main expectations:
- The Gaussian pieces transport as they stand.
- Prior mass needs an M-term approximation rate, and the g-form makes the coefficient mass easy because it is Gaussian in the same norm as the metric.
- Heavy-tailed τ may cost nothing in entropy, unlike the GP case. Each structure is a finite-dimensional subspace, so a local cover needs no norm bound.
- The work sits in the data-dependent priors (σ²'s scale, the level-centred intercept) and the structure prior's normalization.

Start at fixed g, compact σ and a normalized structure prior.

**Would land in:** docs/Maths beside the BGP result; `lean/`. **Pays off in:** G3's Sobol reading, which inherits any L²(μ) rate.

---

## G3: Sobol indices and active subspaces have different standing

**What the literature leaves.** Yoo & Ghosal (arXiv 1411.6716) get derivative rates for tensor-product B-splines, and Odin, Bachoc & Lagnoux (arXiv 2403.03540) subspace recovery for GP priors. Nothing was found for gradient outer products under an adaptive hinge posterior.

**The finding (SG-2, SG-5).**
- √V_u is 1-Lipschitz in L²(μ), so a rate for f passes straight to Sobol. C = E[∇f∇fᵀ] has no such bound: a small triangle wave with many teeth is small in L² and large in gradient.
- E[C]'s top eigenspace weights draws by their eigenvalues, so it is not the mean of the draws' subspaces. They sit 54° apart in a three-draw example.
- Under `EmpiricalMeasure` the atoms sit on knots, so C reads F-X3's kink convention.

**How to close it.** Prove the Lipschitz bound. Before E8's directions reduce anything, add a recovery check for C on functions where it is known, at growing N. Offer E[C] and the subspace mean as separate readouts, and say which measure a C was computed under.

**Would land in:** [E8](extensions-brief.md#active-subspace-analysis).

---

## G4: One moment service across dependent inputs

**What the literature leaves.** Shapley effects are the usual reading under dependence (Iooss & Prieur, arXiv 1707.01334), with analytic forms in Gaussian linear cases and otherwise a double Monte Carlo loop. Nothing was found on mixtures of products here.

**The finding (SG-3).** A joint mixture of products keeps every pair integral closed. One component is `ProductMeasure`; N components with delta marginals are a joint empirical sample. So this replaces a split between a closed-form tier and a Monte Carlo tier with one formula, and latent-class inputs, product-kernel smoothing and possibly Bernstein copulas sit between the two.

**Limits.** Hoeffding does not apply under dependence, so what is exact is the moments and C. Shapley drops from a double loop to a single one, not to a closed form. Separated components weaken the Poincaré constant behind active-subspace bounds. Cost is linear in the components, and the delta case keeps its cheaper direct path.

**Would land in:** `InputMeasures`, `BassMoments`, `SensitivitySpec`; a new extension entry on the go-ahead. It also widens the [guarded estimators](../../../lean/Protolemmata/BASS/guarded-estimator-lemmas.md)' exact case.

---

## G5: Sensitivity as a Markov generator

**What the literature leaves.** The pieces exist separately (see [SG's prior art](../../../lean/Protolemmata/BASS/sensitivity-generator-lemmas.md)). The searches found no sensitivity paper stating the frame, though it may be folklore. R2D2 on polynomial chaos exists (Bürkner et al., arXiv 2204.06043), but nothing connecting coefficient priors to the Jacobi generator's spectrum turned up.

**The finding (SG-1, SG-4).**
- Hoeffding components are eigenfunctions of resampling dynamics, with eigenvalue −|u|. So Σ_v total(v) = Σ_u |u|·V_u is a Dirichlet energy and Efron–Stein a spectral gap of 1.
- Variance by interaction order comes from one polynomial per term pair, from factors `BassMoments` already caches.
- DGSM and active subspaces are the same construction under the Langevin generator, and a generator per variable covers categorical and ordinal inputs.
- Under dependence the spectrum is lost, but the Dirichlet form remains, with Gibbs dynamics' gap measuring how far totals undercount.
- Shifted Legendre is the Jacobi generator's eigenbasis, so a weighted DGSM comes straight from KHAOS's coefficients, and a prior can use the eigenvalues: (5,1) and (3,3) both score 13 under KHAOS's penalty, and 32 and 24 here.

**Connections.** This is the same language as the MCMC family's spectral gaps and the graph Laplacians, which is why SC-2 and SG-1 share a Dirichlet form.

**How to close it.** Formalize SG-1 on finite products, then as an abstract statement. Build the variance-by-order readout, and test the Efron–Stein identity against `SobolDecomposition` as a property check. The prior is an option beside KHAOS's, to earn its place predictively.

**Would land in:** the sensitivity stage; [E5](extensions-brief.md#khaoss-modified-g-prior).

---

## G6: Orthogonal discrepancy for cut calibration

**What the literature leaves.** Gopalan & Francom's note says the cut approach still needs hybridizing with discrepancy modelling. Plumlee (2017) keeps θ identified by making the discrepancy orthogonal to the emulator's gradient, through GP kernels; Xie & Xu's projected calibration avoids the GP, and arXiv 2605.30492 treats constraints too. Nothing was found on a basis discrepancy with a linear constraint.

**The candidate.** Stage 7's discrepancy is δ = D·v with v Gaussian. Orthogonality to ∂_θ f at θ is a linear constraint on v, C(θ)·v = 0, whose entries come from the moment service with `FactorDerivative` when D's columns are piecewise polynomial, or a finite sum on a grid. Conditioning a Gaussian on a linear constraint stays conjugate, so the collapsed likelihood becomes a low-rank update at each θ.

**Open, and not small.** The constraint moves with θ, and a hinge emulator's ∂_θ f is piecewise constant, so the constrained prior jumps as θ crosses a knot. What that does to θ's posterior, and to the cut's modular reading, has not been derived. So this one is a derivation to attempt before anything is posed for Lean.

**Would land in:** D09's options; E14's neighbourhood.

---

## H1: The steep band as a transition in sample size

**Hypothesis, untested.** Under the integrated form, a power posterior at b has the likelihood scale and the Occam term of a full posterior on bN cases. So Friedman's steep band at b ≈ 0.16–0.32 ([readings](friedman-band-brief.md)) predicts that fits on nested subsamples of about 80–160 of the 500 cases switch basin. If they do, the band is a Bayesian phase transition in sample size, in singular learning theory's sense, and the ladder reads where more data changes the model's structure.

**Falsifying it.** Flat-histogram fits on nested subsamples, read by basin. The linear form's Occam term scales differently, so the prediction is for the integrated form only.

---

## Checked and not open

- **Multi-fidelity BASS exists in the lineage's orbit.** Tucker, Collins, Huerta & Brown (arXiv 2609.27174, read to the end of its introduction) use BASS emulators for multi-fidelity mapping and fusion, with elastic alignment. An AR1 model through F-X1's offset with a stage-8 cut may still differ; compare after reading it.
- **Closed-form Sobol and Shapley for BART** are in Statistica Sinica SS-2024-0318, so that is not a gap on the tree side.
- **A variance-share prior on polynomial chaos** is R2D2 (Bürkner et al.). The open piece in G5 is the smoothness-matched one.

## Manuscript readings to confirm

KHAOS's printed birth, death and Mutate1 ratios invert some proposal factors, beyond the double count the hub records ([SC-1](../../../lean/Protolemmata/BASS/structure-chain-lemmas.md)). That is a reading of the manuscript only; the `khaos` package's source and the supplement have not been read. A KHAOS port derives its ratios and tests their reciprocity either way, under F-X4's contract.

## Reading list

**On disk, 2026-11-01, unread unless noted.**
- In `gauntlet/bayesian`: Ronen et al. (2210.09352v1); Kim & Roćková (2306.00126v1); Odin, Bachoc & Lagnoux (2403.03540v2); Tucker et al. (2609.27174v1, read to the end of its introduction); Roustant et al. (1612.03689v1, 1906.09883v1); Bürkner et al. (2204.06043v2).
- In `PDFs/Bayesian`: Yoo & Ghosal (`Ghosal2015-1411.6716v3.pdf`); Woodard, Schmidler & Huber (`Woodard2009.pdf`); Plumlee (`Plumlee2016.pdf`).

**Being tracked down.** Ki, Fang & Guntuboyina, *MARS via LASSO*: arXiv 2111.11694, DOI 10.1214/24-AOS2384.

**Not found on disk by name or arXiv id.** Bhatnagar & Randall (2004); Woodard et al.'s AAP companion, if `Woodard2009.pdf` is the EJP paper; Shen & Ghosal (1403.0625); Xie & Xu (1803.01231); arXiv 2605.30492; Lamboni et al. (1202.0943); arXiv 2412.04918 and 2605.30979; Iooss & Prieur (1707.01334); O'Donnell's *Analysis of Boolean Functions*; Sancetta & Satchell (2004) on Bernstein copulas; Zrnic & Candès (2024); Angelopoulos, Duchi & Zrnic (2023); Giles (2008); Peherstorfer, Willcox & Gunzburger (2016, 2018).
