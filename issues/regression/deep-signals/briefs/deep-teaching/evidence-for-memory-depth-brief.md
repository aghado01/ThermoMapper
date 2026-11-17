# How much evidence supports memory at depth d?

Research brief, 2026-11-17. **Commissioned Lean investigation; section 1 in progress.** For a given observation design, how many independent records, spanning how many events, suffice to detect a dependence at lag \(d\), recover the minimal context tree, or predict the next event to a specified tolerance? The design decision for audio is whether to spend a recording budget on longer spans, different windows, or retaining the pairing between measurements.

The depth-four work supplies one case where a measurement design can be studied exactly. This follow-up starts from the already proved collision at **every** \(d\ge4\). Its first two questions can proceed independently of the depth-four completion proof. [Population cutoff and witness](../context-tree-memory-identifiability-brief.md#verified-result--2026-11-16) · [Depth-four measurement design](01-depth-four-completion-brief.md) · [Predictive ambiguity](02-predictive-ambiguity-brief.md)

## Experiment and recovery targets

Use stationary binary sources invariant under global sign reversal, with positive word masses and memory at most a known \(d\). The *true depth* is the maximum leaf depth of the minimal suffix tree. Start with independent length-\((d+1)\) blocks drawn from the stationary law. Charge for the number of blocks, the span in ordered events, and the measurements retained per block. In the separate-view experiment, each observation comes from a fresh block, with its view chosen in advance. In the paired experiment the learner also knows which measurements came from the same block.

| Decision | Quantity to establish |
| --- | --- |
| Detect a specified deep alternative | Smallest sample size giving both false-positive and false-negative probabilities at most \(\delta\) against a named shallower source. |
| Recover the tree across a class | A uniform bound on the probability of returning the wrong minimal tree, with the frequency of relevant histories and the smallest transition difference specified. Exact depth alone asks only for the maximum leaf depth. |
| Predict from a signed deployment history \(H\) | A bound on \(\mathbb E_H[(\widehat q(H)-q(H))^2]\), where \(q(h)=P(X_{\rm next}=+1\mid H=h)\), with confidence \(1-\delta\) over training. This permits small errors on rare branches. |

## 1. Finite evidence for the existing all-depth pair

Fix \(d\ge4\), a known strength \(0<a<1/2\), and \(0<\delta<1/2\). Compare the fair IID source \(P_0\) with the exact-depth-\(d\) stationary source \(P_a\) in the preceding proof. Every separate rectangular power law through width \(d+1\) agrees. For the taper on sites \(\{0,1,2,d\}\), retain the bit that its DC power equals 16. Its success probability is \(p_0=1/8\) under \(P_0\) and \(p_a=(1+a)/8\) under \(P_a\). These population facts and the connection to the actual Fourier observations have Lean declarations. [Fourier collision](../../../../../lean/Enthymemata/TimeFrequency/ContextTreeFourier.lean) · [Asymmetric event masses](../../../../../lean/Enthymemata/TimeFrequency/ContextTreeCollision.lean)

**First proof target.** For a fixed allocation with \(n\) taper bits and any number of independent rectangular observations, prove the following finite-experiment conclusions together.

1. The likelihood ratio for \(P_a\) against \(P_0\) depends only on the number of successful taper bits. With zero taper bits it equals one on the common support. Under \(P_a\), expected log likelihood ratio in bits equals \(n\,\mathrm{KL}_2(\mathrm{Bern}(p_a)\Vert\mathrm{Bern}(p_0))\). Prior odds between these two specified sources multiply by the realized likelihood ratio.
2. Classify by comparing the observed success frequency with the midpoint of \(p_0,p_a\). Their gap is \(a/8\). One-sided Hoeffding bounds suggest that \(n\ge128\,a^{-2}\ln(1/\delta)\) suffices for each error probability to be at most \(\delta\).
3. Obtain a necessary budget for **any** test of the same data. Data processing to the binary decision gives \(n\,\mathrm{kl}(p_a,p_0)\ge\mathrm{kl}(1-\delta,\delta)\), where \(\mathrm{kl}\) is Bernoulli divergence in natural logarithms. Using \(\mathrm{kl}(p_a,p_0)\le a^2/7\) gives \(n\ge7\,\mathrm{kl}(1-\delta,\delta)/a^2\). Check constants, the direction of the likelihood ratio and integer rounding in the formal statement.

This would give matching order \(a^{-2}\log(1/\delta)\) for small \(\delta\) on the specified pair. At fixed \(a\), the bit budget is independent of \(d\). The selected sites span \(d+1\) events. With \(a\) tending to zero, the required sample size diverges even at fixed depth. A simple-vs-simple likelihood ratio measures evidence for these two sources. Bayesian evidence for competing *depth classes* requires priors on their trees and transition parameters and marginal likelihoods of the **observed powers**.

## 2. Give a constructive recovery bound for arbitrary depth

**Second proof target.** Allow paired adjacent-cell powers on the same length-\((d+1)\) block. Retain the \(d\) unnormalized powers \(Y_i=(X_i+X_{i+1})^2\). Since \(X_i\in\{-1,+1\}\), each \(Y_i=2+2X_iX_{i+1}\) specifies one relative sign. Propagating along the path recovers the ordered block up to its global sign. Under source sign symmetry, choosing either lift with an independent fair coin reproduces the signed block law exactly. Prove this equivalence of experiments and use it to construct an estimator. The cost is \(d\) paired two-site measurements per independent block spanning \(d+1\) events. [Two-sample Fourier recovery](../phaseless-audio/01-two-sample-recovery-brief.md)

Let every length-\(d\) history have probability at least \(\pi_{\min}>0\). For every suffix whose full-history extensions have differing next-positive probabilities, assume the range of those probabilities is at least \(\gamma>0\). Estimate the transitions by counts of histories and extensions in the lifted blocks. A simultaneous error strictly below \(\gamma/4\) lets the estimator distinguish zero range from range at least \(\gamma\), by splitting at empirical range greater than \(\gamma/2\). The shortest suffix with constant transitions determines each leaf.

Derive an explicit sufficient block count. A candidate order is

\[
  \frac{\log(2^d/\delta)}{\pi_{\min}\gamma^2}.
\]

Bound the probability of too few visits to any history, then the transition errors given its visit count. State the constants and the estimator's behavior for unseen histories. Retain \(\pi_{\min}\) in the answer. Even uniform histories have mass \(2^{-d}\). A complementary lower-bound question is whether a rare, predictive split can be hidden at this cost under the same stationarity and symmetry assumptions. The pair in question 1 establishes a lower bound for its discrimination task. Tightness of this whole-tree upper bound requires its own construction.

## 3. Decide which depths a finite sample supports

Fix a candidate maximum depth \(D\) **before** examining the records. Construct a simultaneous confidence set of source trees from uncertainty in the observable frequencies. Project the set to the possible true depths. Its minimum would give a lower confidence bound on memory: reporting “at least \(k\)” requires excluding every shallower explanation. Specify the outcome for an empty confidence set. For a binary split, report both the probability of falsely retaining it and the conditions under which a true split will be retained. These address different errors.

Use the depth-four measurements as the first unpaired case for this rule. Propagate simultaneous histogram errors through the actual reconstruction, including division by rare history masses. Compare its evidence requirements with paired adjacent powers at the same block budget and state the measurement costs of both designs. For a design admitting the exact collision, show why the corresponding depths must remain indistinguishable. For prediction, ask when a shallow suffix already attains a given excess-loss tolerance despite ambiguity about deeper leaves. [Deep Teaching research questions](../deep-teaching-audio-research-brief.md#proposed-lean-questions-in-dependency-order)

## Extension to recordings and proof acceptance

Moving to one ordered recording changes the statistical experiment. Specify its length, window overlap, retained temporal order and source dependence. Prove concentration for that sampling scheme or give a bound in independent recordings with an explicitly defined statistic per recording. Temporal pairing can change the observable law as well as the variance. An elapsed-time statement for speech also needs the mapping from event sites to time.

Kontoyiannis gives bounds on predictive likelihood and posterior consistency for observed discrete time series. Galves and Leonardi give nonuniform exponential bounds for empirical context trees. Use these as comparisons when stating what assumptions the power observations introduce. [Bayesian context trees](https://arxiv.org/abs/2211.02676) · [Empirical context-tree bounds](https://arxiv.org/abs/0710.5900)

For each completed question, record the experiment, source class, estimator or test, error event, and a sufficient or necessary observation budget with its units. Link declarations and report module builds and axiom checks under the [Lean harness conventions](../../../../../lean/README.md). The formulas above are proposed proof targets. The first deliverable is the finite-sample upper and lower bound instantiated on the physical all-depth collision and its distinguishing measurement.

## Commissioning and in-progress record — 2026-11-17

The user, through the research coordinator, commissioned this as the next Lean assignment after depth-four completion and predictive ambiguity. Work begins with section 1, then the paired-observation and supported-depth questions. No depth-evidence theorem is yet reported as verified.

The first implementation will separate the finite observed likelihood experiment, the Bernoulli testing budgets, and the link to the existing physical collision. Rectangular observations are independent ancillary observations with the same law under both hypotheses. A fixed count of independently sampled taper bits has probability (1+a)/8 under the alternative and 1/8 under the fair source. Randomized decisions must be included in the lower bound, or any restriction to deterministic tests must be explicit.

Mathlib supplies bounded-variable Hoeffding through `ProbabilityTheory.hasSubgaussianMGF_of_mem_Icc` and `ProbabilityTheory.HasSubgaussianMGF.measure_sum_ge_le_of_iIndepFun`. The lower bound will require proved finite data processing and the binary decision divergence bound, not an assumed information inequality. The logarithm inequality log x ≤ x−1 gives the proposed per-bit upper bound a²/7. A candidate integer sufficiency statement uses ceil(128 a⁻² log(1/δ)); necessity uses ceil(7 kl(1−δ,δ)/a²).

Proofs, axiom checks, corrections, and sample units will be recorded below as each slice closes. Later sections remain unverified.
