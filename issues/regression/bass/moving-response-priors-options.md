# Priors that stop reading a moving response: the options

The third row of [Upgrades to existing code](extensions-brief.md#upgrades-to-existing-code), written out for the user's call. **The user ruled on 2026-11-02: build steps 1 and 2 of [the order below](#what-to-build-in-order).** Those are built ([the forms](../../../docs/Maths/BayesianRegression.md#bass-priors-that-stop-reading-a-moving-response), [the report](extensions-checklist.md#reports)); steps 3 and 4 are open. Defaults stay the user's, per engine. The row beside it, [a design whose working response moves](extensions-brief.md#upgrades-to-existing-code), is separate work: this note is only about what the priors read.

## What reads the response today

| Read | Where | What for | Kind |
|---|---|---|---|
| The level c | `BassResolved.Resolve`; `BassDesign.PriorLevel` | α's prior centre under level centring | A prior |
| s² | `BassResolved.Resolve` | σ²'s IG(c, c·s²) | A prior |
| ȳ_w | `Resolve`, the target; the design | Least-squares statistics, the zero centre's energy | A statistic: it should move with the design |
| Constant or zero response | `Resolve` | Refusals (D03) | Gates on the priors above |

`Flat` with `Jeffreys` or `Held` reads neither prior.

**What the weighted core has already settled.** Under `WeightedPriorGram.Unweighted`, D03's robust row, the level, s² and the coefficient prior read no weights. So when only the weights move and y stays fixed (a scale mixture at β = 0, a scale per group), this row is already answered. Sabanés Bové and Held's GLM g-prior makes the same choice: its Gram reads the model's known weights, never the moving working weights. What remains is a moving *value*: y − β√w·v (GBMARS), the latent z (probit), (y − ½)/ω (Pólya–Gamma), y − μᵢ (a mixture), y − g(x) (E32), y − ρf_lo (E23).

## The problem

**A moving working response is not data.** y − β√w·v is a function of β, w and v, so a prior that reads it is a prior conditioned on parameters. Two things are coherent: read something that does not move, or carry the conditional into every block. An uncorrected re-read does neither, since its conditionals share no joint, and CONTRIBUTING's rule is that such a design is corrected rather than kept. BARS's `WeightedNormalModel` re-reads at every evaluation, but its y is fixed and the read ignores weights, so it gets the same value each time.

**The constraint.** D03 reads y so that the prior is proper and the fit is equivariant. Noise with a location-scale group, y → ky + m with α → kα + m, a → ka, √w → k√w and v, β unchanged, covers Gaussian, t, Laplace, asymmetric Laplace and Normal-Wald. Under that group a prior can have at most two of: proper, reading no data, equivariant, because a proper measure cannot be invariant under a non-compact group. Each route below gives up one. Binary responses have no such group, which is why a declared prior costs them nothing.

## What GBMARS does

Read from `knrumsey/GBASS` and the paper's §3, 2026-11-02.

- **No prior reads the working response.** z and V⁻¹ are rebuilt each sweep once w, β and v have moved, and only the structure step and the coefficients' draw read them. The priors are a | w, τ ~ N(0, wτI), τ ~ IG(½, N/2), β ~ N(m, s²), and w's below.
- **The structure step is the marginal given w,** the form row 2 built. Coefficients are integrated out; w is a Gibbs block. Under skew w sits in z, so it could not be collapsed in any case.
- **The ridge is a stated choice.** The paper considers the g-prior on BᵀV⁻¹B (citing Alhamzawi and Yu) and declines it because it complicates the local scales' conditional.
- **The intercept is in the ridge,** uncentred and shrunk to 0 with the shared τ.
- **Every read of y is of the raw y, once.** y is not standardized.

| Entry | w's prior | Scale-free | Reads y |
|---|---|---|---|
| `gbass`, `tbass` | GIG(0, 0, 0), so 1/w | Yes: a's prior scales with w | A floor, w ≥ var(y)/N |
| `qbass`, `nwbass` | GIG(−0.1, 0, 0.1), so IG(0.1, 0.05) | No: raw units | The same floor |
| `hbass` | GBP(1, ½, ½), so half-Cauchy on √w | Intended: `scale = var(y)` | That scale, and the same floor |

So GBMARS is the invariant, declared and frozen routes below, by entry point, and never has this problem. Starts and proposal widths also read var(y), which does not touch the target.

**What it leaves open** (read from the code, not run; fixtures if GBMARS parity is ever claimed):
- *Shift.* A shift k in y lands in a₀, so in `pen`, so in τ, and each birth carries −½ ln τ. On τ's conditional mean at N = 1000, M = 10, that is about 1.2 nats per term at k = 100√w and 3.5 at k = 1000√w.
- *Small units under `qbass` and `nwbass`.* b = 0.1 outweighs the residuals once N·w is well under 0.1, and w then sits near 0.1/N.
- *The floor.* Only the β = 0 path applies it; the modified half-normal draw skips it, so `qbass` has it at q = ½ alone. Under `hbass` w is in units of var(y) but the floor stays var(y)/N, so it passes w's natural scale once var(y) exceeds about N.

## The routes

| Route | Proper | Reads | Equivariant | Here |
|---|---|---|---|---|
| Invariant prior | No | Nothing | Yes | Built: `Flat`, `Jeffreys` |
| Frozen read | Yes | Observed data, once | Yes | Built: `BassPriorReference` |
| Declared | Yes | Nothing | No | Built: `InverseGamma`, `Zero`, `Held`, a reference's stated level and κ |
| Hierarchical | Per hyperprior | Per hyperprior | Per hyperprior | Built: `Sigma2Prior.Hierarchical`, on the state and the `Scale` slot |
| Conditional on the latents | Yes | Parameters | Yes | Gram side: `Weighted`; its latent correction is new |
| Pilot read | Yes | Posterior expectation of the moving values, then frozen | Yes | New |

**Invariant prior.** Flat on α, 1/σ² on σ², and the g-prior scaled by σ² on the rest. It reads nothing, so nothing can move under it.
- *Why it is sound:* improper constants on common parameters cancel across models sharing the group (Berger, Pericchi and Varshavsky 1998; Bayarri et al. 2012). That covers structures and also noise laws, so Gaussian against t against Normal-Wald can be compared this way.
- *Propriety:* Fernández and Steel (2000) find the posterior under 1/σ exists for scale mixtures whatever the design or mixing, while moments can fail. Skew is outside what was read. F's exact-fit refusal still has to reach a moving response.
- *Cost:* no b = 0 rung, no SMC from the prior, no absolute evidence; the fractional rung below restores an evidence.

**Frozen read.** A fixed function of the observed y gives a fixed target whatever moves later. This is D03's reading. Richardson and Green (1997) do this with latent allocations in an RJMCMC, and BART and rstanarm do it for scale.
- *Cost:* the double use of the data, and an evidence that stays empirical-Bayes.
- *The read must mean something per model,* so the model names it:
  - Moments of y, as today.
  - Median and MAD under heavy tails.
  - The sample q-quantile under the asymmetric Laplace, where f is the q-quantile and ȳ is not its level.
  - Under skew, E[y] = α + √w·β·μ_v, so ȳ centres that sum, not α. Centring α on ȳ − √w·β·μ_v reads parameters rather than latents; it adds 1/w and 1/√w terms to w's conditional, which stay in the modified half-normal family, and a quadratic to β's, while v's conditional is untouched (derived here, not run).
- *On the Gram:* a frozen W₀, the information at a fixed point (Li and Clyde 2018). `PriorPencil` already takes any prior Gram.
- *Standardizing y at entry* gives the same numbers with every output to map back, so it is not a separate route.

**Declared.** Right where the group is absent. For binary responses, Sabanés Bové and Held's prior has a flat intercept and covariance g·c·(XᵀWX)⁻¹ with c = v(h(0))/h′(0)²: 4 for logit, π/2 for probit, 1 for Poisson log. It is anchored at an imaginary sample, so it reads no data. D03's unweighted Gram is this prior's Gram at unit weights; c is the piece a reference has to carry for a σ² held at 1. Where the group is present, it suits a response whose units are known.

**Hierarchical.** The reference becomes a parameter, and this route connects the others. IG(a, b) mixed over 1/b gives exactly 1/σ² for any a (checked numerically), and a flat centre gives the flat intercept, so the invariant prior is this route's scale-free end and the frozen read is its plug-in. A proper hyperprior needs a scale, read once from y: Richardson and Green's rate on the variance prior, or `hbass`'s half-Cauchy on var(y).
- *For:* proper and equivariant, with heavy tails, so a poor read hurts less than under a plug-in. Half-Cauchy has a conjugate two-block form (Wand et al. 2011).
- *Cost:* a block per sweep, and the marginal is then given the hyperparameter, as with τ.

**Conditional on the latents.** p(φ | z) kept as a conditional prior, with its factor carried into the latent block, whose closed form becomes a Metropolis proposal. It is coherent and it is a different model: integrating v out no longer gives the stated noise law.
- *On the Gram* it has a precedent (Alhamzawi and Yu), and `WeightedPriorGram.Weighted` under a moving W is this. v's conditional gains τ(Bᵢa)²/σ² in its GIG's b, which is conjugate, and ½ln|BᵀV⁻¹B|, which needs the correction: a rank-one determinant update per observation.
- *On the level and s²* I know of no precedent, and it adds nothing the frozen read lacks, since equivariance is already met.

**Pilot read.** The one sound way to read the moving values: as an estimate of a constant. Iterate φ̂ = E[statistic(z) | y, φ̂] over pilot chains and freeze it (Casella 2001). It takes the model's own notion of level and scale, so skew and groups need no formula. Cost: a pilot stage.

**A fractional rung** is not a prior, but it is how the invariant prior gets an evidence. O'Hagan's (1995) fraction f of the likelihood makes the base proper, so ln Z(1) − ln Z(f) is defined, on a ladder that stops at f rather than 0. It reads data only through the likelihood, so latents raise no reading question, though the integrated and linear forms differ once they exist. Not yet checked: whether the ladder can report this.

## By consumer

| Consumer | Group on y | What moves | Routes that fit |
|---|---|---|---|
| Probit, logit (E20) | None | z; ω | Declared with c; σ² held at 1 (built) |
| t, Laplace, variance-gamma, a scale per group | Location-scale | W only | Invariant, frozen or hierarchical; already free of W |
| Skewed generalized hyperbolic (E1) | Location-scale | y − β√w·v, W | Invariant, or frozen with the level on the mean; σ² given (built) |
| Asymmetric Laplace (E1) | Location-scale | y − β√w·v, W | Invariant, or frozen with the sample quantile; σ² given |
| A mixture (E20) | Per component | y − μᵢ, W | `Zero`, since μᵢ carries the level; held |
| E32 | From outside | y − g(x) | `Zero` or declared; given from outside (built) |
| E23 with ρ drawn | Location-scale | y − ρf_lo | Frozen from y, or invariant |

## What to build, in order

D03 stays Gaussian BASS's default. Defaults elsewhere are the user's, per engine, with the rest selectable.

1. **Built: this row, on the seam, needing no latent blocks.**
   - The level's energy Σwᵢ·(ȳ_w − c)² in `CentredFit`; with no reference, no arithmetic is added.
   - `BassPriorReference` on `WorkingResponse`: a frozen level, mean square and coefficient-prior constant (κ, Sabanés Bové and Held's c), with the rule named by the observation model. This changed the seam's shape, as approved.
   - A frozen W₀ on `WeightedPriorGram`.
   - Tests of the routes under a moving response: chains on y and on ky + m visit the same structures. A design cannot follow a response yet, so the response moves by a model per snapshot, not through `ReplaceTarget` as first written here, which swaps only the target.
2. **Built: the `Scale` slot.** The hierarchical forms: a hyperprior on the IG's rate, and half-Cauchy, with their scale from the reference. The rate went on the state, not the target, and a slice step replaced the two-block form as the one built ([why](../../../docs/Maths/BayesianRegression.md#bass-priors-that-stop-reading-a-moving-response)).
3. **With the latent blocks.** The corrected conditional prior on the Gram; the level on the mean under skew; the pilot read.
4. **At the run.** The fractional rung.

## Decided with it, or left open

- **F under a moving W.** The gate reads the weighted Gram, as this upgrade asks. So a latent move can put the current structure below s_min, outside F, where its posterior is 0. A latent block must then treat that draw as off the support. The alternative is to gate on the prior's columns alone, which keeps F fixed but admits a numerically singular weighted Gram. Left to the row beside this one, flagged here since the same decision reaches it.
- **Refusals.** Under a frozen read the constant and zero refusals read the reference, not the moving values, as built. A held σ² lifts them.
- **The evidence.** Under a frozen read it stays empirical-Bayes in the level and scale, as `BassEvidencePriors` already records for D03's default.

## How far each source was read

- **Read:** GBASS `R/gbass.R` and `R/wrappers.R`; `R/nwbass.R` in part, with a search for its floor and w's draw; GBMARS §§2–3 as text, formulas stripped (its TeX was checked earlier for the [brief](extensions-brief.md#generalized-hyperbolic-noise-gbmars)); Sabanés Bové and Held 2011 §2.1 with Table 1, whose values the formula for c reproduces.
- **Abstract or search result only:** Berger, Pericchi and Varshavsky 1998; Bayarri et al. 2012; Fernández and Steel 2000.
- **From memory:** Richardson and Green 1997; BART and rstanarm's defaults; Casella 2001; O'Hagan 1995; Wand et al. 2011; Li and Clyde 2018. Alhamzawi and Yu is known only through GBMARS's citation.
