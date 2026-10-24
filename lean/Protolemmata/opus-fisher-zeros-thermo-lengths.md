That framing is right: the question is how concentrated the posterior is given everything gathered, and volley-to-volley change is a learning curve to watch, not the criterion. In posterior terms, the stopping condition is concentration of the features that define the regions of interest, computed per draw as the readout already does:
- the posterior probability of the significant-peak count, P(K = k̂) ≥ 1 − α;
- each peak location's credible width below a target, set naturally by the window the dense phase will sample;
- the same for the span endpoints.

That uses all the data and discards nothing.

**On "conformal":** I think you mean conformal prediction. It earns a place, but as a check on the posterior rather than as the stopping rule itself.
- **What it gives:** distribution-free coverage for intervals around new observations. It's valid even when the model is wrong.
- **Why that matters here:** BARS's credible bands are honest only if the model is right, meaning the Normal noise and the spline prior. Agreement between conformal intervals and the credible bands says the posterior's confidence is calibrated, and disagreement says the confidence can't be trusted however concentrated it looks.
- **Conformal Bayes** (Fong and Holmes, 2021) fits your constraint exactly. It computes full conformal intervals from the existing posterior draws, by importance-weighting them with each point added back in: no data held out, no refits.
- **Your replicates are the right calibration source.** They measure the noise at each temperature independently of any curve model, and they handle the heteroscedasticity near the transition that plain conformal assumes away.
- **The limit:** conformal covers predictions of χ at a temperature. It doesn't directly cover derived features like the peak location, so the feature criterion stays Bayesian and conformal checks the layer underneath.

**Pushing on "accumulating evidence can only improve the estimate":** it holds only under conditions, and they shape the rule.
1. **Resolution.** The posterior can only be confident about what the grid can resolve. A narrow peak hiding between grid points is invisible, and the posterior will be confidently wrong about the peak count until a volley lands near it. So the rule needs a resolution clause: the grid spacing h_v small against the narrowest feature you care about, or each span containing several grid points.
2. **Concentration isn't monotone in practice.** A later volley can reveal a second peak and widen the intervals. That's the method being honest, so the rule shouldn't treat widening as a failure.
3. **Misspecification.** Under a wrong model, more data concentrates the posterior confidently on the wrong answer. That's exactly what the conformal check guards against.
4. **Monte Carlo error doesn't shrink with data.** Each fit needs enough feature ESS, which is cheap here: just run long.

So the pieces compose:
- **concentration** on the features as the criterion;
- **resolution and calibration** as the conditions that make the concentration trustworthy;
- **the learning curve** across volleys as a report.

I've added this framing to the showcase memory as discussion, not a decision.

# Yea conformal bayes sounds like it would be useful. and stopping condition need not be about the features. I've wondered if i might be able to derive something like a minimal sample density for a uniform grid on a susceptibility curve from the potts model. kind of like a nyquist frequency on temperature that would yield some kind of bound on how many grid points should be suficient to recover the true curve

This is a missing piece of your "approximation economy" pillar. `fable-BARS-lemma.md` parks the Jackson bound because nothing supplies the curve's smoothness, and the Potts model can.

**Why not Nyquist literally.** Nyquist needs a band-limited signal. A finite system's χ(T) isn't band-limited, but it is analytic, and for analytic functions the right analogue is the width of the strip of analyticity around the real axis. Two physical quantities set the scale. They're related, since both blow up near the transition.

**1. Fisher zeros set how smooth the curve is.**
- For N spins, Z(β) is a finite sum of exponentials, so its zeros in the complex β-plane (Fisher zeros, Fisher 1965) lie off the real axis. Thermodynamic curves are analytic within a strip whose half-width is set by the nearest zero, Im β₁.
- **The bound.** Cauchy estimates bound every derivative of χ inside that strip, roughly |χ⁽ⁿ⁾| ≲ n!·M/(Im β₁)ⁿ. Feed that into a cubic spline's error bound (about h⁴·max|f⁗|), and the spacing needed for accuracy ε is h ∝ Im β₁ · ε^{1/4}.
- **So the minimal density scales as 1/Im β₁.** Finite-size scaling says that distance shrinks with system size, Im β₁ ~ L^{−1/ν} on lattices (Itzykson, Pearson and Zuber 1983). Bigger systems have sharper peaks and need denser grids.

**2. Thermodynamic length says how many temperatures are distinguishable.**
- The Fisher information of the Gibbs family in β is Var(E), so ds = σ_E dβ, and the thermodynamic length L = ∫ σ_E dβ counts the statistically distinguishable temperatures along the path (Salamon and Berry 1983; Crooks 2007).
- **Its practical form is the histogram-overlap rule** from Ferrenberg–Swendsen reweighting and replica exchange: neighbouring temperatures need overlapping energy distributions, Δβ ≲ c/σ_E. This one is a direct link to your information-geometry interests.

**Why this fits your constraint.** Both quantities can be estimated from data you already collect. Every SW run at a grid point gives energy samples, so Var(E) and its higher cumulants come for free. Those cumulants are also exactly what the β-derivatives of thermodynamic averages are, so they estimate the smoothness too. That gives a stopping rule without features: stop densifying when every grid interval satisfies Δβ ≤ c/σ̂_E, or when h is at most a fraction of the estimated strip width. It's a sufficiency certificate for the design, checked against data, and it composes with the rest:
- the density rule covers approximation, the bias side;
- the replicates and BARS cover noise;
- conformal Bayes checks calibration.

**Caveats to pressure-test:**
- **Uniform T against β.** dT = −T² dβ, so strip widths in T vary along the axis. A uniform T-grid has to satisfy the worst point, near the transition.
- **Hierarchical SPC transitions** give several near-real zeros, and the narrowest one sets the density.
- **Constants matter.** These are scaling laws with constants, not sharp "n points suffice" theorems, unless the constants are tracked.
- **Data graphs aren't lattices.** The finite-size exponents are unknown there, which is why estimating from data beats prescribing from theory.

There's a Lean angle too. "β-derivatives of averages are cumulants" and bounded-variable cumulant bounds are tractable, and Mathlib has Cauchy estimates. So the boring, load-bearing half of this bound looks provable in your harness. I've recorded it all in the showcase memory as discussion.


# is it possible to derive some kind of initial estimate for the number of points in my starting grid before i start bisecting? kind of like domany calculating the nominal guage to constrain the effective window that gets normalized to 0,1

Checking the existing bracket write-up so the answer builds on what's already settled, then I'll answer.

Yes, and most of it can be had before any SW run, from the same kind of graph-only inputs the gauge uses. But I'd change one thing in the framing first: because each volley's grid contains the last, what's worth estimating up front is the **finest spacing you need before you trust a stop**, not the starting count.

## Why a spacing floor rather than n₀

With midpoint volleys, n_v − 1 = 2^v (n₀ − 1). The fits are independent and nothing is thrown away. So starting at 10 points and refusing to stop before volley v* ends with exactly the same data as starting at n_{v*}. Starting coarse costs nothing, and the early fits give you the learning curve for free.

n₀ still decides two things:
- **Which lattice you're on.** Starting from n₀ − 1 = 9 or 8 gives lattices that never share interior points.
- **What the first fit can express.** Cubic BARS with k interior knots needs at least k + 4 points to pass the admission gate.

The real danger of a coarse grid is the stopping rule firing on an aliased curve. A bump that falls entirely between grid points is invisible to posterior confidence. It's invisible to conformal calibration too, since those points sit on the same grid. So the up-front estimate belongs in the stopping rule, as its resolution clause: a spacing h*, or equivalently a minimum volley v* = ⌈log₂((n* − 1)/(n₀ − 1))⌉.

## Where h* comes from: one small calculation

This is my own derivation, so check it. Take a finite system switching between two states, with an energy gap ΔE and an entropy gap ΔS = ln g (g = the ratio of state counts). Then:
- the switch sits at T* = ΔE/ΔS;
- its 25–75% width is 2 ln 3/ΔE in β (β = 1/T). That is **about 2.2/ΔS in ln T**, whatever ΔE is;
- the nearest Fisher zero is at Im β₁ = π/ΔE;
- its thermodynamic length ∫σ_E dβ is exactly π, for any ΔE and ΔS. The integral reduces to B(½,½), which is the Fisher–Rao distance between two non-overlapping distributions.

So the two Nyquist-style ideas from last time give the same answer for a switch. Every switch costs the same statistical distance, and its width depends on the entropy jump alone.

Now apply it to SPC's two kinds of feature, treating each cluster as one rigid giant spin:

- **A split** (A∪B breaks into A and B; J_AB = the total coupling between them). Aligned, there are q states at energy 0. Misaligned, there are q(q−1) states at energy J_AB. So ΔS = ln(q−1), T* = J_AB/ln(q−1) (the mean-field cut temperature), and the width is **about 2.2/ln(q−1), roughly 0.75 in ln T at q = 20, whatever the cluster sizes.**
  - The Potts fluctuations themselves smear χ on that scale. A cluster that lives less than about 0.75 in ln T never makes its own bump; it shows up as a shoulder.
  - That is the band limit you were after, and it depends only on q. Resolution improves only logarithmically in q.
- **A melt** (a cluster of n_c spins going free). Here ΔS ≈ (n_c − 1) ln q, so it is very sharp. But it is a *step*: the drop in χ at the hot edge of the superparamagnetic plateau. A coarse grid brackets a step to one interval and can't lose it. Pinning it down exactly is the dense phase's job.

So the aliasing risk is short-lived bumps, and their lifetimes can be read off the graph for free:
- **Read lifetimes from PKWang's closed-form ladder** (cut temperature per edge T_e = Hcum_e/ln 2). `PKWang.EstimateBracket` already uses its ends to set the window; the interior holds the lifetimes.
- **The lifetime of a cluster** with at least n_min points (the same size floor the readout's prominence gate implies) is ℓ_c = ln(T_breakup/T_separation).
- **Uncertain constants drop out in ln T.** ln 2 versus ln(q−1), the 1/K̂ gauge, and any multiplicative mean-field shift all cancel, because a lifetime is a ratio of temperatures.

That gives:

  h* = min over clusters of max(ℓ_c, 2.2/ln(q−1)) / c_BARS

Here c_BARS is how many points per feature width BARS needs to register a bump. **I don't know that number; it has to be measured.**

**Worked example.** `AutoGridFixedSweep`'s default bracket spans 100× in T, which is 4.6 e-folds. The smear floor alone gives h* ≤ 0.75/c_BARS. With c_BARS ≈ 3 that's about 0.25, so **n* ≈ 20 is the worst case for any graph at q = 20.** That's your 10 points plus one volley.

## What this means for the [0, 1] normalization

Feature widths are constant in ln T, not in T.
- **Linear mapping starves the cold end.** If the bracket is mapped linearly onto [0, 1], a uniform grid over a 100× bracket puts only about 6% of its points in the bottom two e-folds (0.05–0.37 T_c). Most points land in the flat hot tail.
- **Normalizing ln T fixes this.** Each volley then takes geometric midpoints, and equal spacing means equal feature width. Every split has the same width there, so ln T already is the thermodynamic-length coordinate everywhere except at the melts. It also gives BARS a roughly constant feature scale, which suits a single knot prior.
- **The exact equal-length coordinate needs the switch locations,** which is what discovery finds. Its after-the-fact twin: the SW energy statistics from discovery give σ_E(T), so the dense phase can place points at equal thermodynamic length within each region of interest (the Ferrenberg–Swendsen overlap condition).

## Where this could be wrong

- **Clusters aren't rigid.** Mixed clusters peel apart gradually. That only widens features, which makes the cap safer.
- **Simultaneous k-way splits** have ΔS = (k−1) ln(q−1), so they are k−1 times narrower.
- **Shifts aren't purely multiplicative.** Real fluctuations move transitions relative to the mean-field ladder in ways that don't fully cancel.
- **The 2.2 is a 25–75% width.** A bump is detectable over a wider footprint, since its tails decay exponentially and replicates lower the noise. So the estimate errs conservative.

Two checks would settle this without new engine work:
- **Measure c_BARS** from synthetic bumps of known width, at known spacing and replicate counts, BirthBenchmark-style.
- **Test the 0.75 floor** against the existing dense χ(T) sweeps (the Domany parity battery and the BWD toy) by measuring the narrowest distinct bump in ln T.

# yea but i mean just estimating based on the potts model and the data something like a minimum number of points to sample, a priori, before any BARS or SW sweeps

Yes. Using only the graph and q, no BARS and no SW, the answer is about **14 points** at q = 20 over a 100× temperature window, for any graph, provided the grid is evenly spaced in ln T. This is my derivation from idealized physics and nothing has been measured, so check it.

There are three different minimum counts, depending on what the points have to guarantee. Only the first is the discovery count.

## 1. Seeing every feature: the discovery count

Everything here is a calculation on the graph. Treat each cluster as one rigid giant spin:

- **A split** of A∪B into A and B happens at about T ≈ J_AB/ln(q−1), where J_AB is the total coupling between the two halves.
  - The split's edge is about **w_q = 2.2/ln(q−1) wide in ln T**, which is about 0.75 at q = 20.
  - That width depends only on q, so no bump in χ can be narrower than it.
- **A melt** of cluster C happens at about T ≈ (1−1/q)·d̄_C/(2 ln q), where d̄_C is C's mean internal weighted degree. That's the same form as the gauge-style T_c. It is sharp, but it's a step, and a coarse grid can't miss a step: it just narrows it down to one interval.

From these:

- **The window.** It runs from the first split (the weakest cut between two significant clusters) to the last melt (the most tightly coupled cluster). Its size in e-folds is E = ln(T_last melt / T_first split), plus some padding.
- **The narrowest feature.** This is the shortest-lived significant cluster, measured in ln T from when it separates to when it breaks up. The split edge sets a floor on any bump's width, so w_min = min over clusters of max(lifetime, w_q).
- **The sampling constant is plain geometry.** A grid with spacing h ≤ w/2 guarantees at least two samples above half-maximum on any bump of half-max width w. That's enough to tell a bump from a one-point blip.

On a grid evenly spaced in ln T, that gives:

  **n_min = 1 + ⌈2E / w_min⌉ ≤ 1 + ⌈2E·ln(q−1)/2.2⌉**

The right-hand bound holds for any graph, since w_min ≥ w_q. With E = 4.6 (a 100× window) that's the 14 points above. If the cluster tree has only long-lived clusters, you need fewer.

This replaces my earlier figure of about 20, which divided by a BARS constant; the constant here is geometric. Measurement noise is a separate matter for the k replicates and doesn't change the count.

**The normalization matters a lot.** Evenly spaced in T instead of ln T, the cold end forces h ≤ T_lo·w_min/2, which means about **260 points** for the same window. So the evenly spaced grid on [0,1] should be evenly spaced in ln T, with each volley taking geometric midpoints.

## 2. Recovering the whole curve to accuracy ε: the dense count

This is the Fisher-zero version of the Nyquist idea, and it gives a very different number.

- χ is smooth (analytic) out to the nearest Fisher zero in the complex plane. For a melt that distance is about π/((n_c−1) ln q) in ln T.
- Combining that with the cubic-spline error bound gives a spacing of h ≈ a·(16ε/5)^{1/4}, where a is that distance.
- For a 200-point cluster at ε = 1%, that's h ≈ 0.002, or roughly 2,000 points over the window (order of magnitude).

That's a dense sweep. But the requirement is only global if you sample evenly everywhere; it really applies only near the melt edges. The ratio of count 2 to count 1 is roughly the size of the largest cluster. That ratio is exactly why your approach pays: discovery pays count 1, which doesn't depend on cluster size, and the dense phase pays count 2 only inside the regions of interest.

## 3. Reweighting to any temperature: the replica-exchange count

This is the thermodynamic-length count: neighbouring temperatures' energy histograms have to overlap, so n ≈ ∫σ_E dβ. It grows like √N, because it promises the whole ensemble at every temperature, which is far more than the χ curve needs. It looks like the natural "Potts Nyquist", but it's the wrong count for discovery. It's the right tool for placing the dense phase's points so they support histogram reweighting.

## What's idealized

- **Rigid giant spins.** Real clusters peel apart gradually, which only widens features. So count 1 stays an upper bound.
- **Simultaneous k-way splits** have a narrower edge, by a factor of k−1.
- **Mean-field temperatures shift under fluctuations.** Lifetimes measured in ln T largely cancel the unknown constants, but the window E carries that uncertainty, so pad it.
- **You need to set the minimum size of a "significant" cluster.** Without it, the weakest cut just isolates one loosely attached point and E blows up.

The cheapest check needs no new engine work: measure the narrowest distinct bump in ln T in the existing dense χ(T) sweeps (the Domany parity battery and the BWD toy) and compare it with the 0.75 floor.