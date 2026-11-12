
## Status and caveats

- **Against the completion contract** (2026-11-01; [the completion report](../../../issues/regression/bass-completion-report.md), an audit from R's source up, checked by the main session). Every capability P01–P15 exists and its central behaviour is verified, and the audit found no wrong number. Every item of R's surface is now built or left out by a recorded decision (2026-11-02). Two benchmark measures stand between it and the contract: statistical error per unit time and the sensitivity moment work alone, both in the statistical half, not run.
- **This section describes BASS through stage 4** (S1–S10, 2026-10-28), **with stages 5 and 6** (2026-10-29). Calibration (stage 7) is built to D09, decided on 2026-10-29, and SMC (stage 8) is built for fitting (C4) and cut calibration (C5). The flat-histogram sampler is built for fitting (N17 Phase 2, 2026-10-30). Stage 9 builds the presentation data; rendering and text formatting wait for a consumer above maths.
  - *Stage 5 and 6's choices left as reversible defaults* are listed in their chips' reports (`issues/regression/bass/stage5-outputs-brief.md`, `stage6-sensitivity-brief.md`). Among them: basis sensitivity's integrated form includes the mean curve y_m in the location's effect (confirmed by the user on 2026-10-29; the centred anomaly y − y_m would be a later option); `Predict`, `Sensitivity` and `Effects` take the training data, since the fit keeps neither its knots nor its inputs; a basis fit's intervals are weighted quantiles of its draws.
  - *Open from the qualification:* the defaults mix slowly between structures. R̂ of M read 1.5–2.5 on the recovery surfaces at 5,000 + 2,000 sweeps, as R BASS's chains do, with births and deaths accepting about 1.5%; tempering lowers it, but the local moves bound it. Friedman's σ reads high untempered for that reason and for the fit's bias.
  - *Built:* the contracts and their refusals, the terms and the basis, the structure prior and the proposals, the design and admission, the collapsed target, the sweep, the run with its diagnostics and checkpoints, the self-sizing tempered ladder and the evidence, the pointwise density and the residual check, the qualification runs, and F-M2's local-move options (Phase 1).
  - *Built since:* stages 5–8 as above, and the coefficient prior's two axes, with the normalized ridge and τ's shared family (2026-10-29). The extension seams are built as contracts that no extension fills yet (2026-10-31; [above](#bass-the-extension-seams)).
  - *Not yet built:* the robust coefficient prior, KHAOS's modified g-prior, and BARS on τ's shared family.
  - *Meanwhile,* the decisions and their evidence are in the [fitting segment](../../../issues/regression/bass/port-fitting-brief.md) and its [checklist's reports](../../../issues/regression/bass/port-fitting-checklist.md#reports), and the scope and the lineage register in the [port hub](../../../issues/regression/bass-full-port-brief.md).
- **BASS's efficiency defaults wait on measurement.** The sampling lengths (R's 9,000 burn-in sweeps and 1,000 kept draws) are provisional (the contract choices), as BARS's are (OQ14). The equal move schedule and the untempered default stand until BASS's own pre-registered comparisons say otherwise (F-M1; D05); S10's cost baseline and recovery runs (the qualification section) are where those comparisons start. They found the local moves bound the defaults' mixing (births and deaths accept about 1.5%), so that is the next design item, F-M2, decided 2026-10-28; the run lengths and the ladder's sizing, its censored steep band included, follow it. The deterministic-simulator profile (a wider λ prior, the ridge as the alternative coefficient prior, the residual-dependence check read before intervals) is a proposal, F-R1; the ridge it would compare against the g-form is built (2026-10-29).
- **`EbicPrior` reproduces He, Yang and Kang's count law, not their location measure.** Their knots are subsets of a discrete candidate grid, while BARS's are continuous, so a fit with the same γ and C is the same prior on k, not the same model. For γ > 1 the law rises again past C/2: keep C ≥ 2·k_max.
- **The BAPS `Collapsed` scheme's posterior of λ is heavy-tailed.** Under the uniform-edf prior the density of ln λ decays like $e^{-u}$, so $E[\lambda]$ is infinite; `BapsResult` therefore reports the median. $\ell_p$ is computed as $y^\top y-\hat\beta^\top Z^\top y$ and loses precision above about λ = 10¹⁰ (NaN near 10¹⁴), but the prior keeps the chain far below that range. The λ modes target different posteriors.
- **The efficiency defaults are provisional** (preparation brief, OQ14). The move schedule (constant 0.4/0.4/0.2), the relocation kernel (fixed local Beta, τ = 50), the hot end (0.05) and the run length (2,000 draws per chain) are as they were. The measurements below found regime-dependent winners, so the defaults wait for one designed pass that covers both directions of prior–likelihood conflict, with a criterion fixed beforehand. The adaptive relocation kernel (`adaptKernelTau`) passed the confirmation and is available as a setting.
- **Knot count and peak location mix slowly** (OQ8). Measured at the defaults, with the centred Normal model (`BirthBenchmark`, 10 seeds × 4 chains × 10,000 draws, 2026-10-22):

  | Fit | Knot-count τ (median, range) | Peak-location τ | Draws per chain for 100 effective knot counts |
  |---|---|---|---|
  | Data gap, n = 120 | 8 (7–9) | 1 | 770 |
  | Sine, n = 100 | 14 (12–16) | 19 | 1,440 |
  | Jump, n = 150 | 35 (33–43) | 55 | 3,530 |
  | Sine with outliers, n = 120 | 46 (39–52) | 39 | 4,600 |
  | Skewed design, n = 150 | 55 (42–79) | 248 | 5,520 |
  | Narrow peak, n = 150 | 61 (39–160) | 333 | 6,150 |
  | Oscillation, about 13 knots, n = 200 | 298 (174–776) | 32 | 29,850 |

  - **At the defaults** (2,000 draws per chain, no early stop) a run gets 7–260 effective knot counts per chain, and on peaked fits 6–8 effective peak locations per chain.
  - **Judge a run by ESS,** not by draw count: `BarsResult.Knots` and `Peak.LocationConvergence` report it. An R̂ stop's ESS gate waits on the knot count only (OQ7), which on peaked fits is several times faster than the peak location.
  - **Before step 3.9,** under local births, the sine fit's knot-count τ was about 450.
- **Local births strand isolated knots** (OQ11). This applies when ε < 1 and births are local.
  - **The mechanism.** A death's acceptance carries the reverse birth's density at the removed knot. Under local births, a knot far from every other is almost never removed: from knots at 0.39 and 0.96 that density is $e^{-37}$ or $e^{-61}$, against a posterior gain of $e^{3.8}$ to $e^{5.3}$.
  - **Seen in practice.** In `EbicTests`' runs, chains started that way held their knots for about 1,500 steps.
  - **The default avoids it.** Since step 3.9, births come from the uniform base (ε = 1, above). The knot count's autocorrelation time on the adaptive-stop fit fell from about 390 to 14, so an R̂ stop's ESS gate needs about 1,400 draws per chain there, not 45,000.
- **Units** (preparation brief, OQ19). The variance priors behave differently by engine.
  - **Raw-unit:** BAPS's isotropic Gibbs schemes and `HierarchicalPSpline` apply IG(10⁻³, 10⁻³) to their variances in the response's own units, so a fit can depend on them. Bach & Klein's reference code standardizes y first.
  - **Scale-free:** BAPS's anisotropic scheme (above) and the collapsed scheme.
  - **Until OQ19 is settled,** fit those engines on a response of order 1, or pass `noisePrior: Jeffreys` and scaled smoothing priors.
- **Reporting limits.**
  - Acceptance rates count the sampling phase only, since step 3.4; BGP's included burn-in before.
  - A single chain gets the monitored functionals' rank-normalized R̂ and ESS from its two halves; the classic R̂ of the other functionals, such as the fit map, needs two chains and is unavailable.
