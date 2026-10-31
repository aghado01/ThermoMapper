
Emulation is Francom's application, but it isn't what your synthesis is for. What you've built is inference on *features* of an unknown function: where it peaks, where it crosses a threshold, where it changes regime, which inputs interact. Each comes with a posterior, a check that the answer survives a change of smoothness prior (BARS vs BAPS vs BGP), and a calibration layer on top. Emulation is the special case where evaluations are expensive and deterministic. A lot of empirical science asks feature questions about curves, and usually answers them with point estimates.

**Mechanistic interpretability is emulation in disguise.** A forward pass at temperature zero is a deterministic simulator of the interventions you control, so Francom's framing applies as is:
- Dose-response curves over ablation strength get BARS's threshold and peak readouts.
- Layer profiles get the exact changepoint posterior.
- A factorial partial-lesion design gets BASS's Sobol decomposition, where the gap between total and first-order indices measures redundancy.
- Your original ThermoMapper idea carries over too. Fit sparse, then spend new evaluations near the located feature. On rented GPUs, that's sampling dose levels near a threshold instead of on a uniform grid.

**Neurophysiology is BARS's original home.** DiMatteo, Genovese and Kass built it for firing-rate curves. Peak latency, peak rate and comparisons between conditions were the motivating questions. Your old experimental designs map onto the rest:
- 2D receptive or tuning fields are ABMRS surfaces.
- Response onset is a changepoint.
- A muscimol inactivation design is a factorial lesion, so its analysis is the same Sobol question as the ablation one.

**Payment integrity.** This is where the calibration layer matters most:
- **Flagging with a controlled false-positive rate.** Conformal p-values give finite-sample guarantees. Bates, Candès, Lei, Romano and Sesia (2023, *Annals of Statistics*) showed that Benjamini–Hochberg on them controls the false discovery rate. False flags cost provider goodwill, so "at most 5% of flagged claims are false alarms" is a contractual-grade statement.
- **Expected-payment bands.** GBMARS's asymmetric-Laplace route gives quantile surfaces, with the calibration caveat from earlier.
- **Billing shifts.** A changepoint posterior says *when* a provider's pattern shifted and how sure that is, not only that it did.
- **Scale.** Claims volume is large and MCMC doesn't scale to it, so this works on per-provider or aggregated curves, not on raw claim lines.

**Studying ML training runs.**
- Capability against scale or compute gets a posterior on where an ability turns on.
- Hyperparameter importance becomes Sobol indices with uncertainty, where the usual approach (functional ANOVA on random forests) gives point estimates.
- Training seeds are exactly the exchangeable whole-curve units that E17 and E18 need. Conformal bands for a new seed's learning curve, and posteriors on a curve's turning points robust to the noise model, both work there without extra assumptions. SPC sweeps need you to budget for replicates; seeds come free.

**Where it won't compete.** For raw prediction at large n, or for black-box optimization, GPs, gradient boosting and Bayesian optimization stay the default. BASS's hinges also follow the input axes, so structure along rotated directions needs the projection work from E8. Its distinctive value is the feature readouts with calibrated uncertainty, plus the check across priors.

Mech interp and neurophysiology share the factorial-lesion structure, so one worked example could show the method in both. The neurophysiology half would also connect your published work to the new tooling.

Building them together pays off because they share a small set of primitives. The unlocks are those primitives, more than any one engine, and most of the cross-pollination is the same primitive running under different priors.

## The primitives

1. **Exact readouts are roots and integrals of piecewise polynomials.** BARS's peak and span readouts are one case of something more general. Any feature defined by a polynomial equation in x, f, f′ and f″ has an exact posterior per draw. So do features defined between two fits, since f − g on merged knots is still piecewise polynomial. That covers:
   - peaks, onsets and inflection points (f″ = 0);
   - time above a threshold and area under the curve;
   - the point where two curves cross;
   - break-even points;
   - any condition that clears denominators into a polynomial. For example, revenue p·D(p) peaks where elasticity equals −1, which is the polynomial condition D + pD′ = 0.
2. **Variance decomposition with uncertainty.** Sobol indices and active subspaces work on any tensor-product representation, including a surrogate fitted to a black box.
3. **A noise family.** GBMARS gives heavy tails, skew and quantiles.
4. **Calibration and comparison.** Conformal sets over retained draws, LOO, evidence and proper scores sit on top of every engine. That makes "does the conclusion survive a change of prior" a routine check.
5. **Curves as data.** Functional outputs, elastic alignment, whole-curve bands, and replicates as exchangeable units.
6. **Sequential design aimed at features.** You sample where the feature is uncertain, rather than where the function is. Bayesian optimization is one special case.

## Domains

- **Pharmacokinetics.** Noncompartmental analysis reads Cmax, Tmax, AUC and time above the minimum inhibitory concentration from a trapezoid rule on sparse samples. Each of those is an exact readout here, with a posterior. Continuous glucose monitoring's time-in-range (70–180 mg/dL) is the same span readout. Sparse sampling per subject means it wants a hierarchical layer.
- **Pricing and marketing.**
  - Revenue- or profit-maximizing price is an exact peak readout of (p − c)·D(p).
  - Marginal return on ad spend is a derivative readout.
  - The saturation point of a response curve is an inflection or threshold.
  - Bayesian marketing-mix modeling is already standard practice, but it mostly uses fixed parametric saturation curves.
- **Building energy.** ASHRAE's change-point models for energy baselines are literally hinges: heating and cooling energy against temperature, with a balance point at the knot. BASS's knot posterior is a posterior on the balance point. Extra inputs such as humidity and occupancy get Sobol indices. Measurement-and-verification needs uncertainty on the savings counterfactual, which is the calibration layer's job, with days or weeks as the units because hourly residuals are autocorrelated (E12).
- **Agriculture and climate.** Schlenker and Roberts (2009) found U.S. crop yields fall steeply above roughly 29–32°C depending on the crop. That threshold is a span readout, and its posterior is what an adaptation study actually needs.
- **Analytical chemistry.** Spectra and chromatograms are peak-finding problems where the number of peaks is itself uncertain, which is what your multi-peak lemmas cover. Retention-time drift between runs is a warping problem, so elastic alignment plus peak readouts is a chromatography pipeline.
- **Genomics.** Copy-number segmentation is changepoint detection; circular binary segmentation is the classic method. The exact changepoint engine gives posteriors on breakpoint location and count, and the Cauchy segment model handles outlier probes.
- **Operations and site reliability.**
  - Latency targets are upper quantiles, so p99 latency against load is an asymmetric-Laplace GBMARS fit.
  - The knee, where queueing blows up, is a threshold readout.
  - Deploy regressions are changepoints.
  - Alerting on conformal p-values controls false discoveries.
- **Causal inference.**
  - Heterogeneous treatment effects as a surface over covariates get Sobol indices that say which covariates moderate the effect. Bayesian causal forests (BART) are the comparison point.
  - The point where one treatment overtakes another is a crossing readout.
  - Regression discontinuity and kink designs are knots at a known or unknown cutoff.
- **Interpretable ML.** Any trained model is a deterministic simulator, so emulating it gives analytic global importance with interactions, plus active subspaces. Functional ANOVA on random forests (Hutter, Hoos and Leyton-Brown 2014) is the incumbent. The caveat: Sobol under a product input measure is an interventional question. On correlated tabular features it extrapolates off the data, which is the same problem interventional SHAP has.
- **Agent-based and stochastic simulators.** These are used for epidemics, traffic and economies. Calibrating them is Francom's problem with heteroskedastic noise, which GBMARS, the cut posterior and semi-modular inference cover.
- **Finance.** Yield curves are a classic free-knot spline problem. Implied-volatility surfaces over strike and maturity are tensor surfaces with convexity constraints. GBMARS's normal-inverse-Gaussian and variance-gamma noise families originated as models of asset returns.
- **Sports analytics.** Aging curves are a peak readout (peak age, with a posterior) on misaligned careers, which is an alignment problem.
- **Social science methodology.** Multiverse and specification-curve analyses ask whether a finding survives analytic choices. Your cross-engine check is a Bayesian version of that for smoothness priors.

## Cross-pollination worth exploiting

- **Step or ramp?** Fit the changepoint engine and BARS to the same data and compare them by LOO. That answers whether an apparent discontinuity is a jump or a steep smooth rise. The question comes up in regression discontinuity, regime shifts and dose thresholds, and few toolkits can pose it.
- **One readout layer, several engines, deliberate disagreement.**
  - BASS's degree-1 maxima sit at knots.
  - ABMRS gives smooth extrema.
  - BGP needs a grid.
  - When peak posteriors disagree across engines, the disagreement shows which part of the answer comes from the prior.
- **Three independent Sobol computations.** KHAOS's sums of squared coefficients, BASS's hinge integrals and pick-freeze Monte Carlo compute the same quantity three ways. Agreement among them is a correctness test for each.
- **Conformal as engine selection.** All engines meet the same coverage guarantee, so interval width at fixed coverage is a fair efficiency comparison.
- **The thermodynamic layer applies to everything.** The tempering ladder, free-energy curve and real log canonical threshold came from SPC but apply to every engine. The threshold measures complexity in a way that compares across very different priors.