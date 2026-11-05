# The Francom lineage, read against the engine

A synthesis of the Francom group's papers and what sits around them, read against the extension seams. The engine is built for diverse applications, so capabilities are stated by the class of problem they serve. Nothing here is scheduled. Two companions carry the detail: the [extensions brief](extensions-brief.md) holds each piece of work as an entry (E1–E30), and the [gaps register](lineage-gaps-register.md) holds where the literature stops short (G1–G6, H1).

**What rests on what.** Read here: GBMARS, KHAOS, the co-active subspaces paper, elastic calibration and Francom and Nachtsheim's sensitivity review. Everything else is marked: *unread* means known from an abstract or a search result, and *from memory* means not checked against a source. [Reading leads](#reading-leads) says which papers are on disk.

## One engine with five axes

Every paper in the lineage reuses one core: a reversible-jump sampler over sums of products of one-variable factors, with a conjugate Gaussian core collapsed by a g-prior. Each paper changes one axis.

| Axis | Lineage papers | Seam |
|---|---|---|
| Basis: the term types a model builds from | Hinges (BASS 2018, 2020), Legendre polynomials (KHAOS 2025), ridge functions (BPPR 2024) | F-X3, the design-kernel contract |
| Noise | Gaussian, then generalized hyperbolic (GBMARS 2024) | F-X1's offset, F-X2's latent slot, N-X3 |
| Output | PCA and direct functional (2018), landmark warping (2022), elastic warping (2023), correlated residuals (Andros 2024, 2026) | The outputs segment, F-X6 |
| Functionals: quantities computed from the fitted surface | Sobol (2018), active subspaces (2024), co-active subspaces (2025), Sobol by quantile (GBMARS) | The moment service, O-X4 |
| Use | Modular and cut calibration (2019, Mathews 2024, Gopalan and Francom) | D09, stage 8 |

So the seams cut along the same lines as the authors' own research, which is decent evidence that they are the right seams. The sampler is not an axis in these papers. It is where this engine has gone past them (the ladders, the flat histogram) and where G1 sits.

Reading the basis axis as a choice of inductive bias: hinges for kinks and thresholds, polynomials for smooth surfaces, ridge functions for low-dimensional structure among many inputs. Holding all three behind one contract is the engine's value. [E19](extensions-brief.md#polynomial-chaos-terms) and [E28](extensions-brief.md#bases-that-mix-term-types) follow from that reading.

## Separability is what the closed forms rest on

Every closed form here comes from two facts together: a term is κ·Πᵥ φᵥ(xᵥ) ([`DesignKernel.cs`](../../../src/maths/regression/spline/bass/DesignKernel.cs)), and the inputs follow a product measure. The mean, variance and Sobol indices, main effects, the active-subspace matrix and its diagonal, co-active subspaces and the functional active subspace then all reduce to one-variable integrals.

Two things break it, from different places:
- **A non-separable kernel,** a modelling choice: a directional hinge, a BPPR ridge function, a GP term.
- **A non-product measure,** an analysis choice: dependent inputs, observational covariates, a calibrated posterior such as D09's θ.

Both fall back to the same thing, Monte Carlo over the measure with exact values and derivatives per draw, so only the integral over x is approximate.

**Open, to settle before E8, which is the first consumer that cannot assume a product measure.** A fit made only of separable kernels still loses its closed forms under a dependent measure, so the choice between exact and sampled depends on the kernel and the measure together, per query. `Factors` takes a `ProductMeasure` by type, and that is worth keeping. G4's joint mixture of products then settles the measure side with one formula. The kernel side remains: a kernel with no factors has to be sampled under any measure. Sobol, C and effects would sit above both, knowing neither.

**The same test reaches trees (a reading here, untested).** A hard-split leaf is a product of one-variable indicators, so it is separable and piecewise polynomial of degree zero, and the moment service could read a forest. That fits the closed-form Sobol and Shapley results for BART the register records. Soft routing stays separable, but its integrals are not polynomial. Oblique splits are not separable. This concerns moments only: the [note on trees](../perplexity-bart-etc.md) argues against forcing trees through the spline fitting core, since leaf sufficient statistics do the work there.

## F-X1 is a data-augmentation seam

F-X1 was built for noise, but it is more general. Many observation models turn Gaussian once an extra latent variable is conditioned on, and each then hands the core a working response and weights.

| Observation model | Latent | Working response, weights | Beyond the seam |
|---|---|---|---|
| Generalized hyperbolic (GBMARS) | v | y − β√w·v, 1/(w·v) | Weighted Gram, GIG sampler, Bessel K |
| Probit (Albert and Chib) | z, truncated | z, 1 | σ² held at 1, a truncated-normal draw |
| Ordinal (Lawrence et al. 2008) | z, cutpoints | z, 1 | As probit, plus cutpoints |
| Logit (Pólya–Gamma; Polson, Scott and Windle 2013) | ω | (y − ½)/ω, ω | Weighted Gram, σ² held at 1, a Pólya–Gamma sampler |
| Logistic noise, or logit with a latent z (Stefanski 1991; Holmes and Held 2006, both unread) | v, and z for binary data | y or z, 1/v | Weighted Gram, a draw of v given the residual, which is neither GIG nor generalized beta prime |
| Scale shared across a declared group | One scale per group | y, the group's weight | Weighted Gram, indexing by group |
| Dirichlet-process mixture of normals (DPMBART) | μᵢ, σᵢ | y − μᵢ, 1/σᵢ² | Weighted Gram, the DP's state and moves |

The weighted Gram unlocks most of these, so it is the most useful single piece here. It lets GBMARS split in two: first the weighted Gram, an engine change to centring and the marginal, then the GIG and GH pieces, which are model-specific. [E20](extensions-brief.md#data-augmented-observation-models) has the rest, including what probit needs that was first called free. The weighted Gram and σ² held are built since 2026-11-02 ([the upgrades](extensions-brief.md#upgrades-to-existing-code)).

A scale that varies with x, as in heteroscedastic GPs or NFT BART, is outside this seam. It needs a second modelled surface, not a latent per observation.

## Two seams that turn out to be shared

- **Coefficient priors.** KHAOS's modified g-prior needs the collapsed core to accept a prior precision other than ΨᵀΨ/g. SC-3 shows its G∘ΨᵀΨ is ΨᵀΨ + S₀⁻¹, the ordinary conjugate form. If that becomes a declared builder, ridge blocks are a second kind of prior it can build, and random-effect intercepts on any grouping variable are ridge blocks ([E21](extensions-brief.md#a-prior-precision-builder-with-ridge-blocks)).
- **Output-side inner products.** F-X6 has the identity, and E12 plans the tridiagonal Markov residual. A low-rank-plus-diagonal Q is a third member, for correlated responses ([E22](extensions-brief.md#correlated-outputs-and-coupling-within-a-run)).

## Stochastic simulators: a latent through the noise or through the basis

A simulator whose output varies at a fixed input can have a latent variable in either of two places, and the latent slot can hold both, though they are different models.
- **Through the noise.** GBMARS's latent scales give heavy tails and skew, and a unimodal output law at each x.
- **Through the basis.** Zhu and Sudret's stochastic polynomial chaos gives the emulator an extra latent input, which captures multimodal output laws without replicate runs ([E26](extensions-brief.md#stochastic-emulation-through-latent-inputs)).

Item response theory is the parametric case of the second: an emulator with latent inputs, ability and difficulty.

## What "simulator" means in this literature

A simulator is a program encoding a mechanistic theory. One run can take hours, so a study gets perhaps 100p–1000p runs for p inputs. The input space is designed, not observed. Rerunning a deterministic code returns the same answer, so its error is model error; agent-based and Monte Carlo codes give a distribution at each input. Inputs come in two kinds, conditions you control and parameters θ you want to infer.

The workflow (Sacks et al. 1989; Kennedy and O'Hagan 2001; Higdon et al. 2008, from memory), with where each step lands here:
1. Design runs that fill the input space. Outside the engine today.
2. Emulate with a cheap stand-in that carries its own uncertainty. BASS, BGP and the entries on the basis axis.
3. Analyse sensitivity on the emulator. The moment service, E8, E25.
4. Calibrate against observations, y = f(x, θ) + δ(x) + ε. D09 and stage 7.
5. Choose the next runs. No seam yet ([E30](extensions-brief.md#choosing-the-next-runs)).

That keeps four kinds of uncertainty apart: which parameters (θ), interpolation (the emulator saw only some runs), stochasticity (ε, or the simulator's own randomness) and discrepancy (δ).

**Discrepancy and covariate shift are different failures.** Reweighting to another input measure handles covariate shift: the same function under a different input distribution. E25's measures and the moment service do that. Discrepancy is the function itself differing at the same x, which no reweighting reaches, and it is learned from a small anchored set where both are observed. Without it, Brynjarsdóttir and O'Hagan (2014, from memory) show that more data gives tighter intervals around the wrong quantity. So a sensitivity or calibration reading says which of the two it has handled. G6 is the open piece on the cut side.

**Where the framing stops.** Sobol indices mean only what the declared measure means, so a domain with no parameter box has to supply one or live with dependent inputs. Active subspaces need continuous, differentiable inputs. Emulators pay only when runs are expensive.

## Build order

**GBMARS first, then KHAOS: the user's call, 2026-11-04,** with CRAN's GBASS pinned as the oracle before any of it ([the build brief](gbass-noise-brief.md)). This reverses the order first proposed here, KHAOS, E8, GBMARS, which called GBMARS the heaviest of the three by counting the weighted Gram. That was built on 2026-11-02, so what remains is a GIG sampler, Bessel K, the latent state and the noise law's readers. E8 keeps its place after KHAOS.
1. **GBMARS (E1, with E20).** The noise axis: a working response at a latent state, on seams built for it. The structure chain stays as it is.
2. **KHAOS (E19, with E4, E5, E21).** The same lineage on other axes: the basis, the proposals and the coefficient prior. It tests F-X3 with a member that is not a hinge and builds F-X4 as already decided. It is deterministic and conjugate, and its only new distribution is the Poisson-binomial. It pays for grouped-data problems through E21. GBMARS builds none of this, so it is a diff of its own, with CRAN's `khaos` as its oracle. Its method carries over, and so does the product: a polynomial-chaos basis under GBMARS's noise needs nothing more.
3. **E8.** It only consumes what is built, and it has strong internal checks. The consumer shape above gets settled here.

E24 is small, sits beside all of them and can land at any point.

## Reading leads

| Lead | Bears on | Status |
|---|---|---|
| Rumsey, Gibson, Francom and Morris, *All Emulators are Wrong…*, arXiv 2512.09060, with `duqling` | E27 | Unread; not on disk |
| Zhu and Sudret, stochastic polynomial chaos, arXiv 2202.03344 | E26; the strongest adjacent lead | Unread; not on disk |
| Hutchings, Bingham, Lawrence (and Rumsey in v2), arXiv 2405.16298: fast emulation, modular calibration and active learning for functional output | Output and use axes, E30 | Unread; not on disk |
| Francom, Sansó and Kupresanin, landmark-warped emulators, *SIAM/ASA JUQ* 10 (2022) | E11 | Unread; not on disk |
| Zahm et al., arXiv 1801.07922, vector-valued dimension reduction | E8 | Unread; not on disk |
| Stefanski, a normal scale mixture representation of the logistic, *Statistics & Probability Letters* 11 (1991); Holmes and Held, *Bayesian Analysis* 1 (2006), from memory | E20; the mixing law is in GBASS | Unread; not on disk |
| Sun, Kong and Pal, the modified half-normal, *Communications in Statistics* (2021) | E1's global scale; `rMHN` is in GBASS | Unread; not on disk |
| Wycoff, Binois and Wild, arXiv 1907.11572, closed-form C for GPs | E8, for BGP | Unread; not on disk |
| Collins, Francom and Rumsey, BPPR, arXiv 2210.09181 | E8's end point | Unread; not on disk |
| Rumsey et al., *Discovering Active Subspaces*, arXiv 2307.11241 | E8 | Unread; not on disk |
| Zhang et al. 2003, *Identification of Linear Directions in Multivariate Adaptive Spline Models* | E8's frequentist precursor | On disk (`PDFs/Bayesian/Zhang2003_StatAssc_98-462.pdf`), unread |
| Tucker, Collins, Huerta and Brown, arXiv 2609.27174 | E23 | On disk, read to the end of its introduction |
| Francom and Nachtsheim, arXiv 2506.11471 | E25: it sorts methods by "enough runs for an emulator?" and "independent inputs?", and the engine covers one quadrant | Read |
| Leclercq 2018, BOLFI (`PDFs/Bayesian`) | Calibration with no tractable likelihood | On disk, unread |

## Against the notes in `issues/regression/`

The [Perplexity notes](../perplexity-bart-etc.md) are about trees and an integrated framework, and were read against this review on 2026-11-02. They agree with it more than they conflict, and two of their points check out against the code.
- **An engine boundary they name, which holds for BASS too.** BARS's `IObservationModel` integrates out coefficients and scale, while BART's tree update integrates one tree's leaves given the others and a shared variance. Gaussian BASS is on BARS's side: D04 integrates a and σ² out of the structure step. So an additive model with BASS as one component needs the structure step's marginal given σ², which D04's generalized-hyperbolic mode reads when σ² is not integrated and E20's probit needs too. F-X1's offset is already the partial residual (a reading here, untested). Their warning stands: such a sum is a new model, with centring to declare and component recovery to test.
- **The component multiplier.** Their a_i·g(x) for BCF's treatment forest is a weighted Gram in other clothes: weights a_i², working response r_i/a_i, and a_i = 0 carries no information. So a varying-coefficient z·τ(x) is one more consumer of E20's core (a reading here, untested).
- **Checked against the code.** `Truncated<T>` cannot draw, which E20's probit needs. And [ProbabilityAndInformation.md:147](../../../docs/Maths/ProbabilityAndInformation.md) says no distribution has a quantile, while `Normal.Quantile` exists and line 51 lists it. That line was corrected on 2026-11-02 (`95c7a7c`), and the draws and the other quantiles were built the same day ([the report](extensions-checklist.md#reports)).
- **Already absorbed.** Their tempering distinction, ∫p(y | s, η)^b against (∫p(y | s, η))^b, is the integrated and linear forms. Their naming point, that GBMARS's "generalized" is not generalized Bayes, is in E13.
- **Not taken up here, and open.** A capability inventory for the regression families, narrow capability descriptors that refuse an unsupported plan before a run, a composition layer, and four separate verdicts: target correctness, computational adequacy, predictive adequacy and evidence qualification. Their order puts this after the spline families and trees at the end, as a test that the engine is basis-agnostic. That does not conflict with the order above, and this review's axes would seed the inventory. These are now E33, with the tree family's engine points in E31 and the GP's in E32 (2026-11-02). The tree family's reading list (SBART, DART, GP-BART, BCF, survival and NFT BART, XBART) stays in the notes, along with the question the note asks at BCF, of whether its decomposition is a direction for studying bias in machine learning.

The [evals note](../opus-bayesian-evals.md) applies these ideas to AI evaluation as a stochastic computer experiment. It is one consumer. The engine pieces it asks for are entries here, stated by problem class. Where the two differ, this review and the register are the later reading: multi-fidelity BASS is in the lineage's orbit after all (Tucker et al.), and probit needs more than the note's "fits now".
