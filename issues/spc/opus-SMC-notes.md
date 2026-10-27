Of the four, only OQ2 fits now, and only partly.

**Not now: the BASS fitting adapter, the cut-calibration adapter, and the cut-SMC bridges.** They need a BASS model to adapt, which is stage 4, plus a settled tempering path (D05) and calibration target (D09). Once those exist the adapters are thin, because N8 already holds the hard parts.

**OQ2, the exact reference rung for parallel tempering: yes, for BARS.** The idea is that the coldest-from-the-data rung, b = 0, draws a fresh independent sample from the reference each time instead of taking an MCMC step. Every trip down the ladder then restarts from an exact draw, which is the setting of Syed et al.'s round-trip analysis. The brief assigned it to the SMC stage because it needs what N8 just built: exact draws from a normalized reference. For BARS:
- **Powered path:** buildable now. The reference is k ~ Poisson(1) with uniform knot locations. It's defined even with Normal BARS's improper σ² prior, because powering never re-integrates σ².
- **Likelihood path:** buildable now for the linear-form models such as Poisson, where the reference is the knot prior restricted to admitted configurations, drawn by rejection.
- **Normal model on the likelihood path:** needs a proper σ² prior, which is OQ19 and part of your design pass.
- **The OQ1 rerun:** the brief pairs this with rerunning OQ1, the comparison that set BARS's default tempering path. That rerun could change the default, so it would come back to you as a decision.

For BASS the rung matters less. D05's draft gives BASS's powered path no normalizable reference, so only its likelihood path in the proper-prior mode gets the rung.

**SPC population annealing: what it was about.** It came up in the SMC architecture discussion ([astra-smc-in-thermomapper.md](issues/regression/astra-smc-in-thermomapper.md)).
- **Today:** SPC's Swendsen–Wang sampler runs a fixed grid of temperatures as independent jobs, and each one equilibrates from scratch.
- **The proposal:** run SW through the N8 core instead, with each particle a whole Potts spin configuration. Stepping from β to β′ reweights each particle by exp(−(β′ − β)·H(s)), then resamples, runs SW sweeps at β′ and measures.
- **What it buys:**
  - each temperature starts from the previous temperature's equilibrium instead of from scratch;
  - adaptive spacing crowds temperatures where the distribution changes fast, which is exactly SPC's transitions: susceptibility peaks and cluster splits;
  - free-energy and partition-function estimates come with an exact reference, Z(0) = qⁿ at β = 0.
- **What it costs:** resampled particles are correlated copies, so you need the family and ancestry diagnostics N8 now reports. It also needs two changes on the SPC side:
  - a way to move a spin state to a new temperature with fresh RNG and accumulators, since `SwendsenWang` refuses to restore a checkpoint at a different temperature;
  - weights that use the energy of the configuration actually transferred, since today the energy is measured before recolouring.
- **First validation:** tiny Potts graphs whose equilibrium can be enumerated exactly, then a comparison with independent SW runs at equal cost.

It's SPC-track work, outside the BASS contract, so it doesn't fit your "don't digress" for now. It is the most natural first real consumer of the SMC core, though, if you want one later.