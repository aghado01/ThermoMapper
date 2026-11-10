# Independent random slopes — slice 3 qualification

E21, 2026-11-09. Commission: [ridge-blocks-brief.md](ridge-blocks-brief.md), decision 6 and slice 3. The user's three-tier rescaling contract and all 94 combined focused checks pass. Commit selection passed 2,124/2,124 with no skips. The full Slice run will follow on this immutable implementation commit. Independent sign-off and every later slice remain held.

## Implementation

- `BassRandomSlope` names a continuous input and explicitly selects `Unit` or `Standardized` scaling. JSON must contain the scaling field. The config owns an ordered value copy, with equality and hashing, as it does for intercepts.
- The training contribution is `z_i b_g`. Level-space constraints have rows `a_g = Σ_{i in g} z_i` and `c_g = Σ_{i in g} z_i²`. A joint thin SVD reveals their rank; an orthonormal complement defines the isotropic constrained prior. This avoids independently magnifying a nearly zero constraint row. No default or correlated intercept/slope covariance is introduced.
- Intercepts and slopes share `IBassRandomEffect`'s grouping/scale/prior contract and the existing joint row whitening. The one-intercept closed form remains in place. Retained records resolve by `(Grouping, Input)`, with null input for intercepts; null input is omitted from JSON. Each slope has its own held or drawn scale.
- Training pointwise evaluation, known-level prediction, population spread and conditional mean surfaces multiply effects by the fitted input scale. Existing nested-pair validation remains in the common group resolution. Other noise families, multi-output readout/sensitivity and conformal paths are outside this slice.

`Standardized` uses the training root mean square deviation (divisor n), after the existing unit transform. Positive changes of physical units are the rescaling cases below. Input reversal is a different model under `Unit` scaling because it changes the origin of an independent random slope.

## Evidence so far

Builds use isolated archives plus only this slice's files. [The receipts](../../../artifacts/e21-ridge-s3-qualification/) ([tracked inputs](../receipts/e21-ridge-s3-qualification/)) record `baseline.txt` for the historical absolute-tolerance trials and `accepted-contract-baseline.txt` for the final contract, including the newer readout/conformal commits. The archives exclude R. Runs use Release, below-normal priority, `DOTNET_PROCESSOR_COUNT=2`, one xUnit thread and disabled build servers. The final private build root is `artifacts/chip-ridge-blocks/slice3/accepted`; a pre-existing neighboring snapshot was left untouched.

The historical blocker handoff, `handoff-focused.trx`, had 86 passing cases and two absolute-tolerance failures, with no skips. The revised prediction/rescaling batch, `accepted-rescaling.trx`, passes all 24 cases under the user-approved contract:

| Check | Result |
| --- | --- |
| Full dense case covariance, whitening, log determinant and Gaussian marginal, both input scales, with/without intercept | Passed |
| Dense conditional effect means and covariance, including posterior intercept/slope cross covariance | Passed |
| Both centring constraints in conditional and retained draws, to 10⁻¹² | Passed |
| Dependent constraint rows, zero-dimensional slopes, several named slope inputs and mixed predictor positions | Passed |
| Value equality/hash, caller mutation isolation, required JSON scaling and invalid declarations | Passed |
| Plain and ladder two-chain resume | Passed bit for bit |
| Known/population predictions, fitted training scaling, pointwise and dense conditional mean surface | Passed |
| Existing valid/invalid nesting cases, within-parent geometry, unequal-scale whitening and named-group resume | Passed |
| Dense enumeration of the structures and both scales, 64 versus 96 probability-coordinate quadrature nodes | Passed for both input scales; Δln Z ≤ 1.5×10⁻¹⁴, ΔE[ln s] ≤ 3.5×10⁻⁴ |
| Geweke, independent intercept/slope, both effect scales and τ drawn | Four cases of 10,000 pairs passed; worst standardized paired difference 2.491 |
| Positive power-of-two input/domain rescaling (`c=4` and `c=8`), held and drawn scales | Full fits identical bit for bit |
| Nonexact input/domain rescaling (`c=3.7`), held scales | Relative error ≤ 8.73×10⁻¹⁵, below 10⁻¹² |
| Nonexact rescaling, primitive whitened columns/response, marginals and 1,000 conditional draws at a fixed state | Relative error ≤ 1.63×10⁻¹⁴, below 10⁻¹² |
| Nonexact rescaling, both scales drawn | Every decision identical; full continuous-state/draw relative error ≤ 3.70×10⁻¹¹, below 10⁻⁹ |

The separate `legacy-pins-and-resume.trx` run passed all five cases: the three old one-group draw hashes, both ungrouped fingerprints (one test), and one-group resume. Earlier failed compilation receipts record two test-call migrations, not engine failures.

### Recovery and the scope of its comparison

The controlled recovery benchmark has a linear population surface, held coefficient shrinkage and an unknown slope scale. Eight unbalanced-group datasets use two-point or shifted-exponential effects projected into the convention, rather than effects sampled from its centred normal. Rao–Blackwellized mean surfaces measure held-out and population RMSE. This isolates the specified level/slope convention; Geweke separately exercises drawn coefficient shrinkage.

| Reading | Unit | Standardized |
| --- | --- | --- |
| Pooled 90% effect-interval coverage | 29/32 | 31/32 |
| Mean standardized effect error squared | 0.857 | 0.719 |
| Realized contrast RMS/σ in the 90% scale interval | 4/4 | 4/4 |
| σ in its 90% interval | 4/4 | 4/4 |
| Worst level/slope standardized error | 1.822 | 1.831 |
| Largest held-out known-group RMSE | 0.0842 | 0.0679 |

Population RMSE matches the no-block control to the reported precision. The eight comparison claims are registered in `ComparativeSweep`; `SEED_SWEEP=3` passed all own-seed, chain-axis and data-axis comparisons (`recovery-seed-sweep.trx` and its linked tables). This is a controlled convention check, not a claim of general population-risk dominance.

Exploratory trials with a weak quadratic population term and drawn coefficient shrinkage failed the per-dataset population no-harm check, despite recovering slope effects and the level/slope in the cases reached (`recovery-initial.trx`). A linear trial with shrinkage drawn also failed that comparison (`linear-recovery-and-rescaling.trx`, `linear-rb-recovery.trx`). These receipts predate the joint constraint-rank correction and are preserved as bounded preliminary observations. The controlled benchmark does not supersede or qualify those regimes.

### Adopted rescaling contract and complete decision path

The user set the three tiers in decision 6 after reviewing the original failures. Both declared input scales are tested at ×4 and ×8 with held and drawn scales. These exact-arithmetic controls compare every sampled parameter's double bits, including signed zeros. At ×3.7 the same fit parameters have retained absolute differences up to 1.4448×10⁻¹¹ under `Unit` and 4.4067×10⁻¹² under `Standardized`; the new contract bounds their relative differences rather than treating accumulated adaptive-trajectory rounding as an absolute 10⁻¹² failure. No input rounding was introduced.

For paired vectors a,b of length m, define R² = (Σaᵢ² + Σbᵢ²)/(2m). The symmetric relative error is maxᵢ |aᵢ−bᵢ| / max(|aᵢ|, |bᵢ|, R). The compared family has common units: a whitened column, a draw's coefficient vector, one block's effects, or the three enumerated structure marginals. R supplies the family's measured amplitude for exact or near-zero coordinates; multiplying the units multiplies both the numerator and R. Both vectors identically zero have zero error. Positive scalar parameters (σ², τ, λ and each effect scale) instead use |a−b|/max(a,b). There is no absolute floor of 1. RMS is evaluated after scaling by the largest magnitude to avoid square underflow. A test changes tiny vectors' units by 10⁻¹⁵⁰ and verifies the same relative error, alongside zero and deliberately out-of-tolerance controls.

An internal, non-fingerprinted observation hook reads the actual untempered `Bass.Fit` chain after every sweep. The rescaling fixture uses scalar Gaussian noise, raw ridge with **τ held at 0.7**, and one intercept plus one independent slope on grouping `g`. The sampled-scale cases draw **both group scales**, starting at 0.9 and 0.6; the held-scale controls keep those values. There is one chain, 200 burn-in sweeps and 100 retained sweeps, with thinning 1. These are distinct from the Geweke checks above, which draw τ as well as both group scales.

The hook records all 300 sweeps, not just retained structures: each move's and block's cumulative attempt/accept counts, null/refusal tallies, empty births, structure, both RNG states, scales and τ. The rescaling fixture has exactly one structure update plus one update per drawn scale per sweep; explicit attempt-count checks establish that each 0/1 accept-counter increment is one individual decision. Thus rejected scale proposals cannot be hidden by a later accepted move. Both scales must have rejected proposals in burn-in and sampling. Every discrete reading and RNG state agrees, and a separate observed/unobserved fit comparison confirms the hook consumes no randomness and changes no draws.

Historical absolute failures remain in `handoff-focused.trx` and earlier receipts. `sources-before-contract.tar`, `source-manifest-before-contract.json` and the original baseline preserve the prior handoff source. The final receipts also hold code/test hashes checked against the frozen build source, the slice's source bundle, explicit paths and live working-tree state. The private build snapshot remains until the immutable-commit Slice run finishes.
