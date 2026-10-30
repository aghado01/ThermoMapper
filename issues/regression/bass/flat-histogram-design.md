# A flat-histogram sampler, shared (design draft)

Date: 2026-10-29. Status: a draft for the user's design pass. The user decided to build it whatever the reference-path test reads ([that test's brief](friedman-reference-brief.md)); what is left to decide is how, which "The decisions" below lays out with recommendations.

## Why

- **BASS: Friedman.** The likelihood-tempering path crosses a first-order transition at b ≈ 0.23–0.28. There ℓ₁ is bimodal at fixed b within each chain, its phases 120–190 nats apart at a fixed M ([the measurement](friedman-band-brief.md#reports)).
  - Parallel tempering completes no round trip there.
  - SMC's lineages collapse through it at every particle count read, and its replicates' ln Ẑ disagrees by 6.5–8.4 nats up to 16,000 particles.
  - Tempering along b is the path that has the transition.
- **SPC: Potts.** Potts models with q > 4 have first-order transitions, where Swendsen–Wang slows sharply and fixed-temperature sweeps straddle phases. Frustrated SPC (`issues/spc/opus-spin-glass-threads.md`, Thread 4) needs a sampler that crosses such barriers.
- **One readout for both.** The density of states gives the free-energy curve ln Z(b) at every b with no path discretization. That is E13's curve for BASS, and the spin-glass threads' Thread 5 plank 3 for SPC.
- **The standard remedy** for first-order transitions is a flat-histogram method: sample the energy evenly, so a chain passes between phases through the disordered region instead of tunnelling at a fixed temperature (Berg and Neuhaus 1992; Wang and Landau 2001).

## The idea for a Bayesian target

- *The energy and the base measure.* Take the untempered log-likelihood ℓ₁(x) as the energy, and the prior p(x) as the base measure.
- *The density of states* g(ℓ) is then the prior's mass at likelihood level ℓ, nested sampling's object (Skilling 2006).
- *The multicanonical target* is π_w(x) ∝ p(x)·w(ℓ₁(x)), with w ≈ 1/g over a window [ℓ_lo, ℓ_hi], so that ℓ₁'s histogram under π_w is flat there.
  - Outside the window, w continues as e^{b_lo·ℓ} below and e^{ℓ} above, so the chain neither wanders into the prior's far tail nor past the posterior.
- *From one run on the frozen weights:*
  - the posterior, by importance weights e^{ℓ₁}/w(ℓ₁);
  - every tempered posterior on the linear path for b in [b_lo, 1], by e^{b·ℓ₁}/w(ℓ₁);
  - the evidence curve Z(b) = E_prior[e^{b·ℓ₁}] = E_w[e^{b·ℓ₁}/w]/E_w[1/w], a ratio estimator with no rungs, which needs the prior normalized (D03's defaults are).
- *The energy is generic.* Every target exposes `LogPrior` and `LogLikelihood(1)` (`ISamplingTarget`), so π_w wraps any engine's target. BASS's ℓ₁ is its collapsed log marginal at b = 1. D05's integrated form plays no part: b = 1 is the same posterior on either form, and the reweighted intermediate b are the linear form's.

## The decisions (for the user's pass, with recommendations)

1. **The algorithm.** Learn the weights by stochastic approximation during burn-in only, then freeze them. The production run is plain MCMC on π_w, under D10's adaptation policy.
   - *The learning:* SAMC (Liang, Liu and Carroll 2007, *JASA* 102, 305–320), with a gain decaying as 1/t. It has a convergence theorem, and its estimates keep improving. Belardinelli and Pereyra's 1/t Wang–Landau (*Phys. Rev. E* 75, 046701, 2007) is the physics form of the same fix.
   - *Not adopted:* the original Wang–Landau, which halves its modification factor on a flatness check. Its error saturates, and it has no convergence theorem. Its flatness check stays as a diagnostic.
2. **Chains.** During learning, the chains share one weight function, as PAWL's interacting walkers do (Bornn, Jacob, Del Moral and Doucet 2013, *JCGS* 22(3), with a Bayesian variable-selection example). Production then runs independent chains on the frozen weights, so they share one target and D10's verdict applies to them. Belardinelli and Pereyra (2016) find that walkers beyond a critical count no longer reduce the learning error, so the count stays modest (the ensemble's 4).
3. **The warm start.** Initialize ln g from a short self-sized tempering pilot. Its rungs' ℓ₁ histograms combine by multiple-histogram reweighting (Ferrenberg and Swendsen 1989; WHAM) into an estimate of ln g over the pilot's range. The machinery exists (`TraceRungs`, the self-sized ladder), and learning from that start corrects the pilot instead of discovering g from flat.
4. **The window and the weights.**
   - The window runs from ℓ₁'s low quantile at the pilot's hot end b_lo to its high quantile at b = 1, with margins.
   - ln w is piecewise linear in ℓ₁ over K bins, continuous so bin edges add no artefact. K defaults from the window's width over the posterior's ℓ₁ sd, and is selectable.
   - Adaptive binning, PAWL's, is a later option.
5. **The readouts.**
   - The posterior's weighted draws, through `DrawWeights`, which every view already reads (stage 8).
   - The ln Z(b) curve, with its Monte Carlo error from the chains' spread.
   - The mixing measure: tunnelling, the round trips between the window's ends.
   - The b = 1 weights' Kish ESS, and the histogram's flatness.
   - D10's verdict on the production chains, for π_w's own functionals. The posterior's errors come from the chains' spread of the reweighted estimates, as SMC's replicates do.
6. **Where it lives.**
   - A generic engine in `Maths.Samplers.FlatHistogram`: the weight function, SAMC's learning, the multicanonical stage density, the run, the reweighting and Z(b).
   - *One seam, shared:* the stage density. Parallel tempering, SMC and this engine all evaluate "log prior + a function of ℓ" through one interface (today's `Temperature`, generalized). BASS's structure moves, already run at any SMC stage or ladder rung, then run on π_w unchanged.
   - *Adapters:*
     - BASS first (`BassSampling.FlatHistogram`, a spec, null by default);
     - SPC later, in its own track: the Potts energy H with single-site moves, and the multibondic cluster algorithm (Janke and Kappler, *Phys. Rev. Lett.* 74, 212, 1995) as the cluster variant. It is cautioned against at large q.
7. **Validation.**
   - *In the Samplers host* (no Regression dependency):
     - a double well with a known g;
     - a small q-state Potts lattice whose density of states enumerates exactly (q^N up to about 10⁷), at a q and size where the energy's distribution is bimodal near the transition: the learned ln g against the exact one, the reweighted laws, and Z(b) at every b. The harness picks the size and reports the bimodality;
     - the seam's change leaving parallel tempering and SMC bit-identical.
   - *BASS:*
     - fixture A's enumerable F, with τ held, whose exact g comes from each structure's ℓ₁ and prior mass. Check the learned weights, the reweighted posterior and Z(b) against it, then sampled τ by quadrature;
     - then the Friedman qualification: whether the production chains' reweighted posteriors and ln Z(1) agree to about a nat, the basins' masses, tunnelling counts, and the cost against SMC.
8. **The phases.**
   - Phase 1: the engine, the seam and the Samplers validation.
   - Phase 2: the BASS adapter, its enumeration checks and the Friedman qualification. That qualification reads; it sets no default.
   - Phase 3: SPC's adapter, in SPC's track.
9. **The name.** "Flat-histogram" is the umbrella across fields, "multicanonical" the target (Berg), SAMC and Wang–Landau the learning algorithms, and "generalized ensemble" the physics umbrella. The recommendation is `FlatHistogram` for the namespace, with `MulticanonicalWeights`, `SamcLearning` and `MulticanonicalRun` inside, so each name stays true of its method.

## What to expect, honestly

- *The overhead is the window's width.* Production spends equal time per ℓ bin, so the share near b = 1 is about the posterior's ℓ₁ spread over the window's width. On Friedman that is roughly 10 nats in about 400, so the b = 1 ESS per sweep falls by a factor of about 40 against an ideal sampler. The bet is that crossing between basins costs far more than that today.
- *Learning can be slow.* The window is wide, and BASS's moves can jump hundreds of nats in ℓ₁, a load-bearing death among them. The warm start is the mitigation; the Friedman qualification reads whether it is enough.
- *First-order barriers fall from exponential to polynomial in the Potts case,* the classic result. For BASS's structural phases nobody has measured it; Phase 2 does.
- *Not adopted here:* nested sampling. It integrates the same g, but needs draws from the prior above a likelihood threshold, which BASS's reversible-jump moves give only by a constrained chain. The flat-histogram family uses the kernels this repo has.

## Open questions for the user

The nine decisions above, chiefly:
- 1 (SAMC with a 1/t gain, learning in burn-in only);
- 3 (the WHAM warm start from a tempering pilot);
- 6 (the shared stage-density seam);
- 8 (the phases).
