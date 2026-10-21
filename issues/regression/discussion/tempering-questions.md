First, the scope. Q5 only applies when `temperLevels` is set, and it never changes the posterior, because the cold rung is the same on both paths. What it settles is what the hot rungs do and how well the ladder mixes. It is also the precedent BASS, SMC and evidence estimation will inherit. So it matters little for BARS today and a lot as architecture.

## The case for the likelihood path

**1. The powered path's hot rungs aren't a property of the model.** Step 3.2 showed powered hot rungs drifting to Poisson(λ^β) knots, but that number depends on how the knots are stored, not just on the model. I checked the same BARS posterior under two equally valid representations:

| Rung | Prior mean | Powered, ordered knots (ours) | Powered, labeled i.i.d. knots |
|---|---|---|---|
| λ = 4, β = 0.25 | 4 | 1.41 knots | 5.67 knots |
| λ = 5, β = 0.05 | 5 | 1.08 knots | 14.17 knots |

The same model and the same β give hot rungs with about 1 knot in one representation and about 14 in the other. The likelihood path's hot rungs approach p(k) in every representation. That is the real force of the reference-measure point: under the powered path, the scaffolding's behaviour is an accident of bookkeeping.

**2. Coherence downstream.**
- **BASS:** on the likelihood path, the λ and τ Gibbs conditionals are the same at every rung. On the powered path they would have to be tempered too, which is what R's `rgammaTemper` does, and R's mix has no single target.
- **SMC:** it walks this same path from prior to posterior and reuses the ℓ_b machinery.
- **Evidence:** thermodynamic integration over this ladder gives log p(y). On varying dimension the powered path has no usable β = 0 endpoint.

**3. Swaps are cheap.** Each one is O(1) from cached block statistics.

## The case against, or at least for caution

**1. No efficiency advantage measured.** Between-seed spread of mean knots was about the same for both ladders (SE 0.017 against 0.018). That isn't an efficiency comparison, but it also gives no empirical win. The case is structural.

**2. The hot rungs sit where the prior is.** With a diffuse count prior, such as `UniformPrior(50)`, likelihood hot rungs wander to around 25 knots. That makes proposals expensive, and the prior is far from the posterior, so the ladder may need more rungs to bridge them. Powered hot rungs (in our ordered bookkeeping) sit at few knots: cheap, and arguably closer to a small-k posterior. Point 1 above says that closeness is accidental, but accidents can still be convenient. Step 3.3's round-trip measurements are where this gets settled.

**3. The evidence benefit isn't available to BARS as specified.** The improper 1/σ² prior refuses b = 0, and evidence with an improper prior on σ² isn't defined anyway. Getting it would need a proper IG(g₁, g₂) prior. That moves BARS slightly off the DiMatteo–Genovese–Kass marginal and belongs with BASS's prior decisions.

**4. Only the Normal model at unit weights can use it.** Poisson, Whittle and Gibbs-loss would need tempered Laplace marginals (the IRLS objective with the log-likelihood scaled by b). That's feasible but not built.

**5. `temperBetaMin` = 0.05 was never chosen for this path.** At b = 0.05 with m = 120, the hot rung still carries about 6 observations' worth of data. The ladder shape is step 3.3's job.

## My recommendation

Decide Q5 on principle now: likelihood path as BARS's default, mainly because of point 1 and BASS coherence. Record a revisit trigger with it: if 3.3's round-trip measurements show the powered ladder winning clearly per unit of work on the validation fits, reopen the default.

Three things I'd like your take on:

1. **Default or required?** The sampler has no default path, but BARS does. Your strict-core instinct might argue for making `temperPath` required whenever `temperLevels` is set, which puts the declaration where it's used.
2. **Proper σ² prior.** Is it worth flagging now as a BASS-stage question, since it decides whether tempered runs can ever produce evidence? Or leave it parked with the coefficient-prior decisions?
3. **Tempered marginals for Poisson and Whittle.** Worth putting on the list, or is the powered path enough for those models?

---
