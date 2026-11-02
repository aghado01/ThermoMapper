# Priors that stop reading a moving response: the options

The third row of [Upgrades to existing code](extensions-brief.md#upgrades-to-existing-code), written out for the user's call before it is built. Nothing here is decided. The row beside it, [a design whose working response moves](extensions-brief.md#upgrades-to-existing-code), is separate work: this note is only about what the priors read.

## What reads the response today

| Read | Where | What for | Kind |
|---|---|---|---|
| The level c | `BassResolved.Resolve`; `BassDesign.PriorLevel` | α's prior centre under level centring | A prior |
| s² | `BassResolved.Resolve` | σ²'s IG(c, c·s²) | A prior |
| ȳ_w | `Resolve`, the target; the design | Least-squares statistics, the zero centre's energy | A statistic: it should move with the design |
| Constant or zero response | `Resolve` | Refusals (D03) | Gates on the priors above |

**What the weighted core has already settled.** Under `WeightedPriorGram.Unweighted`, D03's robust row, the level, s² and the coefficient prior read no weights. So when only the weights move and y stays fixed (a scale mixture at β = 0, a scale per group), this row is already answered. What remains is a moving *value*: y − β√w·v (GBMARS), the latent z (probit), (y − ½)/ω (Pólya–Gamma), y − μᵢ (a mixture), y − g(x) (E32), y − ρf_lo (E23).

**Why it cannot simply re-read.** A prior read from z is p(α | z) and p(σ² | z). The joint then has those factors, and z's full conditional gains them. The latent block's closed form (a truncated normal, a GIG, a Pólya–Gamma) is then no longer the conditional, so an uncorrected re-read samples something other than the stated model. BARS's `WeightedNormalModel` re-reads at every evaluation, but its y is fixed and the read ignores weights, so it gets the same value each time. It would meet this problem only with a response per sweep.

**What every option but re-reading needs.** With c frozen and ȳ_w moving, the marginal must carry the level's energy, Σwᵢ·(ȳ_w − c)². Today's `CentredFit` takes c = ȳ under `Mean`, so its energy is zero. The general form adds exactly 0 for Gaussian BASS, so the hashes should hold. `PriorPencil` already carries it.

## The options

**A. Declared constants only.** σ² takes `InverseGamma(a, b)`, `Jeffreys` or `Held`, all of which exist. The level takes `Zero` or `Flat`, which exist, or a new declared level. A model whose values move says so, and `Resolve` refuses `ScaledToResponse` and `Mean` with the reason.
- *For:* the least code. It uses no data, so an evidence under it is a plain evidence.
- *Against:* it gives up D03's scale-free default exactly where the raw units are large or unknown (OQ19). The zero centre under-fits a response far from 0, and the flat intercept loses b = 0 and the evidence.
- *Fits:* probit and logit, where σ² = 1 fixes the scale and N(0, 1/(τN)) on α is natural; a mixture; E32. *Fits badly:* GBMARS on raw units.

**B. A frozen reference.** The observation model names what the priors read once, at `Resolve`, and they never read it again: the response y itself, or declared numbers (a level and a mean square). Gaussian BASS and any fixed offset keep today's reading, the working response.
- *For:* it keeps D03's scale-free, shift-equivariant default for GBMARS, group scales and E23. It is D03's own empirical-Bayes reading, stated once, and A is its special case with declared numbers. The chain's target is a fixed model.
- *Against:* the reference must mean something per model. y's level is not the level of y − β√w·v under skew, and y ∈ {0, 1} says nothing about (y − ½)/ω, so those declare numbers. It keeps the double use of the data that D03 accepted. For E32, y includes the GP's share, so s² overstates the BASS part's noise, though only weakly at c = 10⁻³.
- *Seam:* `WorkingResponse` gains an optional reference, such as `BassPriorReference(Level, MeanSquare)`. This changes the seam's shape, which has no recorded ruling.

**C. Sample them.** A flat hyperprior on α's centre is the flat intercept, which exists. σ² ~ IG(a, b) with b ~ Gamma on a declared scale is conjugate and would go in the `Scale` slot.
- *For:* coherent with no data read, and the scale adapts.
- *Against:* it adds a block and still states a scale in the response's units, unless the hyperprior is improper, which loses b = 0 and the evidence. No consumer asks for it yet.

**D. Standardize once, at entry.** y → (y − m)/s with m and s fixed, and unit-scale constants on the result. OQ19 settled the scale-sensitive τ this way.
- *For:* the priors are plain constants, and the numbers are B's.
- *Against:* every output maps back: σ², the coefficients, predictions, the densities' Jacobian and the evidence. That is a much larger change than B for the same frozen numbers.

**Ruled out: re-reading each sweep.** Uncorrected, it fails its stated intent, as above, and CONTRIBUTING's rule is that such a design is corrected rather than kept as an option. Corrected, by a Metropolis step on the latents with p(α | z)p(σ² | z), it gives up the closed-form latent blocks that make these models cheap.

## Recommendation

**B, with declared numbers as its explicit form.** It keeps the default the user chose in D03 wherever a reference means something, and A's constants remain available where it doesn't, without a second mechanism. C can join later as a prior kind, since it does not conflict. D is B's numbers at a higher price.

| Consumer | What moves | Level under B | σ² |
|---|---|---|---|
| Probit, logit (E20) | z; ω | Declared, 0 or Φ⁻¹(ȳ) | Held at 1 (built) |
| GBMARS (E1) | y − β√w·v, W | Frozen from y, or declared | Per E1's u-prior; given (built) when drawn apart |
| A scale per group | W only | Already free of W | Already free of W |
| A mixture (E20) | y − μᵢ, W | Declared, or `Zero`, since μᵢ carries the level | Held at 1 |
| E32 | y − g(x) | `Zero` or declared, as one component carries the level | Given from outside (built) |
| E23 with ρ drawn | y − ρf_lo | Frozen from y, or declared | Frozen from y, or declared |

## Decided with it, or left open

- **F under a moving W.** The gate reads the weighted Gram, as this upgrade asks. So a latent move can put the current structure below s_min, outside F, where its posterior is 0. A latent block must then treat that draw as off the support. The alternative is to gate on the prior's columns alone, which keeps F fixed but admits a numerically singular weighted Gram. Left to the row beside this one, flagged here since the same decision reaches it.
- **Refusals.** Under B the constant and zero refusals would read the reference, not the moving values. A held σ² already lifts them.
- **The evidence.** Under B it stays empirical-Bayes in the level and scale, as `BassEvidencePriors` already records for D03's default.
