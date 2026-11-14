# GBASS noise, slice 5: the engine's reach (2026-11-07)

Landed: `37137f2` (code and tests) and `8e42f54` (docs, ledger, the brief's report). Stopped after slice 5, as asked. Slices 6 and 7 wait for the user, beside KHAOS's build.

- **What runs now under a noise:**
  - Tempered ladders and SMC, on the likelihood path in the integrated form, both the configs' defaults. Each rung tempers the likelihood given the latents, and every noise block reads b (`BassCoefficientsBlock.Power`).
  - An exact reference rung, and SMC's start. These draw the latents from their prior (`BassPriorSampler.DrawPrior`, `VarianceMeanMixture.PriorLatents`).
  - Evidence wherever the values are complete.
  - A direct functional response, with a latent per observation.
  - Categorical predictors.
- **Refused under a noise, each with its reason:**
  - the linear form and the powered path;
  - the flat histogram;
  - a Jeffreys σ² under a skew, wherever b = 0 draws it.
- **What MCMC fits keep:** a chain's start keeps v = 1, so their draws are bit for bit what they were.
- **Ladder fix found by the resume test:** under a noise, the replica that takes the reference rung's draw now adopts the sampler's design (`BassLadder.HandOver`). Before, a resumed ladder threw at its first swap from the reference, because the replica read the response as the old model had. Gaussian ladders are unchanged.
- **Open reading:** the skewed horseshoe fails Geweke at b = ½, with |z| 5.3 to 9.2. The evidence points at the scan crawling along a β·v ridge:
  - the step stays on its conditional under a drift at b = ½;
  - the symmetric horseshoe at ½ passes;
  - the skewed horseshoe at ¾ passes.
  
  A joint move along that ridge is a follow-up.
- **Decision for the user:** whether the flat histogram takes a noise. The brief's report gives two routes, (a) and (b), each a sampler of its own, and (c) keep it refused.
- **Cost:** a rung at b = 1 costs what an MCMC chain costs, and slice 1's t sweep did not move (962 µs against 961). The reach cost row now runs an untimed warm-up pass, because a cold JIT had read 1.5× slow.
- **Boundaries:**
  - I touched nothing of the conformal layer's.
  - Its candidate density refuses every direct functional fit, Gaussian included, because whole curves are E17's. So no follow-up for it here.
  - `r/` is untouched.
- **Gates:**
  - The whole regression assembly passes 2,121 of 2,121 on `e43f209`, and 2,160 of 2,160 on `b408bf4`, each with my diff. The pinned Gaussian hashes are unmoved.
  - Distributions 347, samplers 285 and the GBASS oracles 27, all passing.
  - Geweke: seven rungs, worst |z| 2.62.
  - Evidence against brute force: z at most 2.90.
- **Stale row handed off:** the GBMARS row in BayesianRegression.md's "What each extension still writes" is stale for slices 4 and 5. I asked khaos-slices to carry the update, because that file holds its uncommitted edits. If it declines, the text is in my message to it.
- **Lesson:** gate.ps1's `-Qualification` set BASS_QUALIFICATION in the PowerShell session, and a later "full suite" in the same session ran every qualification gate. That took 5+ hours instead of the usual 8–18 minutes. The gate now clears the variable when the switch is off.
