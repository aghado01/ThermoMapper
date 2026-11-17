# How much evidence supports memory at depth d?

Research brief, 2026-11-17. **Commissioned Lean investigation; sections 1–3 proved.** For a given observation design, how many independent records, spanning how many events, suffice to detect a dependence at lag \(d\), recover the minimal context tree, or predict the next event to a specified tolerance? The design decision for audio is whether to spend a recording budget on longer spans, different windows, or retaining the pairing between measurements.

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

## Commissioning and sequence — 2026-11-17

The user, through the research coordinator, commissioned section 1, then the paired-observation and supported-depth questions. After this brief's sequence, the user has queued [predictive information in audio](../predictive-information-in-audio-brief.md), starting with its one-bit obstruction and three-target joint-activity retention results.

## Section 1 result — 2026-11-17

**Proved finite-experiment bounds.** Fix depth d=m+4, known strength 0<a<1/2, and error target 0<δ<1/2. Draw each recorded view from a fresh stationary length-(d+1) block. Choose all rectangular widths and offsets in advance; retain the complete Fourier power vector within each rectangular view. In addition, record n independent bits indicating DC power 16 on taper sites {0,1,2,d}. No temporal overlap or data-adaptive view allocation is covered.

For k successful taper bits, the likelihood ratio of the specified deep source against the fair source is

\[
 (1+a)^k(1-a/7)^{n-k}.
\]

The independently recorded rectangular vectors have a common joint law under both hypotheses and cancel from the likelihood ratio. With n=0, that ratio is one on the common support; no test has both error probabilities below 1/2. Expected log evidence under the alternative is n kl((1+a)/8,1/8) nats, or this quantity divided by ln 2 in bits. For positive prior weights u,v on these two sources, posterior odds are (u/v) times the realized ratio.

| Guarantee | Number n of retained taper bits |
| --- | --- |
| Necessary for any randomized test with both errors ≤δ | n ≥ ceil(7 kl(1−δ,δ)/a²) |
| Sufficient for the midpoint count test | n ≥ ceil(128 ln(1/δ)/a²) |

The midpoint test accepts the deep source when k≥n(2+a)/16. Each error is at most exp(−na²/128). Necessity uses proved finite data processing, including randomized decisions t(data)∈[0,1], and kl((1+a)/8,1/8)≤a²/7. Also kl(1−δ,δ)=(1−2δ)ln((1−δ)/δ)>0. Thus the bounds have the same small-δ order a⁻²log(1/δ); they are not an exact optimal constant.

**Units and interpretation.** Each taper bit costs one independent block, four selected sites, and a span of d+1 ordered events. At fixed a, the bit budgets do not depend on d. Additional rectangular records cost blocks and measurements but carry no discrimination information for this pair. At fixed δ, the necessary bound diverges as a tends to zero, so exact-depth detection cannot have a uniform finite budget over arbitrarily weak members of this family. This is evidence between two known sources; it does not supply Bayesian marginal likelihoods for entire depth classes. The existing stationary finite-state construction supplies the block experiment; no infinite path-space measure or single-stream concentration claim is added.

**Declarations.** All names are under `TimeFrequency.ContextTree.DepthEvidence`.

- [MemoryDepthLikelihood.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthLikelihood.lean): `taper_likelihood`, `zero_tapers`, `expected_log_evidence`, `expected_bits`, `posterior_odds`.
- [MemoryDepthTesting.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthTesting.lean): `binary_data_processing`, `testing_divergence`, `taper_divergence_upper`, `necessary_budget_rounded`, `decision_divergence_pos`.
- [MemoryDepthUpper.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthUpper.lean): finite-law/measure interface, `bit_upper_tail`, `bit_lower_tail`, `midpoint_errors`, `sufficient_budget_rounded`.
- [MemoryDepthObservations.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthObservations.lean): actual raw-record pushforwards `raw_taper_records`, `raw_rectangular_records`; retained-law equality `rectangular_experiment_eq`; `taper_sites_and_span`, `every_depth`.
- [MemoryDepthEvidence.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthEvidence.lean): `common_only_obstruction`, `physical_budgets`, instantiating both rounded budgets with the actual independent rectangular observations.

**Verification.** All five modules compile in the aggregate. The 49 public `#print axioms` outputs contain only subsets of `[propext, Classical.choice, Quot.sound]`; the new modules have no proof holes or postulated axioms. They remain in Enthymemata pending semantic review and promotion. Full `lake build` passed (3804 jobs), using BelowNormal priority and `LEAN_NUM_THREADS=2`. The new modules emitted no warnings; existing unrelated notices remain. The taxonomy gate also passed when section 2 was registered.

## Section 2 result — 2026-11-17

**Paired experiment.** For any positive sign-symmetric law on ordered length-(d+1) binary blocks, retain d adjacent powers (X_i+X_(i+1))² from each block. Their values 0 and 4 specify the relative signs. Propagation along the path gives exactly two signed lifts. Choosing their initial sign with an independent fair coin reproduces the original signed block law. This statement needs sign symmetry, but not stationarity. Stationarity and the memory bound enter when interpreting the recovered block kernel as the source's transition law and minimal tree.

The proofs connect actual raw paired-power records, their relative-sign law, randomized lifted datasets, and history-response records. The observation cost is **n independent blocks, each spanning d+1 events and retaining d paired two-site powers**. Fair coins are internal randomization, not additional recordings.

**Estimator and assumptions.** Let π_h be the probability of length-d history h and q_h its next-positive probability. Assume π_h≥π_min>0. For every suffix whose extensions have nonconstant q_h, assume some two extensions differ by at least γ, with 0<γ≤1. For n reconstructed blocks, let N_h count visits and S_h count positive responses. Set qhat_h=S_h/N_h when N_h>0 and qhat_h=1/2 otherwise. Treat a suffix as nonsplitting exactly when every pair of its extension estimates differs by at most γ/2; return the shortest such suffixes.

For 0<t≤1/2 the proved bounds are

\[
\begin{aligned}
 P(\exists h:\ |\widehat q_h-q_h|\ge t)
 &\le 2^{d+1}\exp(-n\pi_{\min}t^2),\\
 P(\exists h:\ N_h\le n\pi_h/2)
 &\le 2^d\exp(-n\pi_{\min}/8).
\end{aligned}
\]

The first bound includes unseen histories. It follows from a direct exponential-moment bound on the masked count-ratio numerator, rather than the brief's proposed intermediate conditioning on the visit count. The separate low-visit bound is also proved. No conditional-on-count theorem is claimed. On the simultaneous event |qhat_h−q_h|<γ/4, the empirical nonsplitting predicate equals true transition constancy for every suffix, hence the returned tree equals the minimal tree.

Consequently a sufficient number of blocks for tree error at most δ is

\[
 n\ge\left\lceil
 \frac{16}{\pi_{\min}\gamma^2}
 \ln\frac{2^{d+1}}{\delta}
 \right\rceil .
\]

This is an explicit sufficient bound, not an optimal sample complexity. Uniform histories already have π_min=2⁻ᵈ. Neither the π_min dependence nor the separation assumption can be omitted from this guarantee. The whole-tree lower-bound construction remains open, distinct from section 1's proved simple-pair testing necessity. Overlapping windows from one recording remain outside the experiment.

**Declarations.** Namespace `TimeFrequency.ContextTree.DepthEvidence`.

- [MemoryDepthPaired.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthPaired.lean): `paired_reconstruction`, `paired_record_law`, `randomized_lift_law`, `randomized_lift_records`, `history_response_records`.
- [MemoryDepthConcentration.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthConcentration.lean): `finite_chernoff`, `bounded_mgf`, finite union bounds.
- [MemoryDepthCounts.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthCounts.lean): `estimate`, `masked_sum`, `ratio_mgf`, positive conditional-law parameters.
- [MemoryDepthBounds.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthBounds.lean): `transition_error`, `simultaneous_transition_error`, `transition_budget`, `too_few_visits`, `any_too_few_visits`.
- [MemoryDepthTree.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthTree.lean): `conditional_kernel`, `threshold_contexts`, `threshold_tree`, `paired_tree_budget`.

**Verification.** Full `lake build` passed (3809 jobs), BelowNormal priority with `LEAN_NUM_THREADS=2`. The 42 public theorem axiom listings in section 2, and 91 across sections 1–2, contain only subsets of `[propext, Classical.choice, Quot.sound]`. The new modules emit no warnings. The taxonomy gate `lean/scripts/meta-ci.ps1 -NoBuild` passes for sections 1–2, with all ten new modules PROOF-CLOSED. They remain Enthymemata pending semantic review and promotion.

## Section 3 result — 2026-11-17

**Supported depths.** Fix the maximum history length D and candidate class before observing data. Form a simultaneous confidence set of valid source laws, then project it to the maximum leaf length of their minimal suffix trees. The minimum supported depth is a lower confidence bound; an empty set returns `None`, not a depth claim. Whenever the set contains the truth, this minimum cannot exceed the true block-kernel depth. Overstatement therefore has probability at most the noncoverage probability.

Here “depth” is formally the depth of the D-history conditional kernel. The current `Source` structure expresses stationary sign-symmetric cylinder laws, not a Markov memory bound. Interpreting this finite-kernel depth as the entire process's memory requires the declared restriction to sources of memory at most D. The confidence-set construction does not infer that restriction from the data.

A suffix is retained as a split only when the confidence set is nonempty and every candidate has nonconstant transitions there. For a truly constant suffix, false retention has probability at most noncoverage. In a transition band of radius t, a true difference greater than 4t forces retention on coverage; missed retention then has probability at most noncoverage. These are distinct error statements. Section 2 supplies the paired experiment's simultaneous coverage budget.

**Unpaired depth-four experiment.** Record n independent length-five blocks for each of six preselected views: rectangular DC powers of widths 2, 3, 4 and 5, and the two four-site tapers omitting sites 3 and 2. Each block supplies one scalar power. Normalize the four rectangular powers by 4, 9, 16 and 25. The width-four and width-five groups also supply their top-bin indicators; the two taper groups supply their top-bin indicators. These eight statistics lie in [0,1]; the two statistics sharing each rectangular group need not be independent.

Write m for their eight empirical means. Their expectations are the actual eight source summaries used by the explicit depth-four inversion. Simultaneously,

\[
 P\{\max_j|m_j-\mathbb E m_j|\ge\varepsilon\}
 \le16e^{-2n\varepsilon^2}.
\]

Thus n ≥ ceil(ln(16/δ)/(2ε²)) suffices for noncoverage at most δ, for ε>0 and 0<δ<1. The charge is **6n independent blocks, span five events per block, 6n retained scalar powers**. Prefix views use only their specified sites. No overlapping-recording or adaptive-allocation bound is asserted.

**Stability through the actual inversion.** The eight summaries recover all 32 five-symbol atom masses by an explicit affine map. A kernel-checked integer coefficient certificate bounds every row's absolute sum by C=385/16. Hence summary error at most ε gives atom error at most Cε. For two valid sources with every four-symbol history mass at least η>0, atom error at most E gives transition error at most E/η. This ratio bound uses the validity of the source masses and the history floor.

Restrict the candidate class to that common history floor and retain sources whose summaries lie within ε of m. On coverage, every retained candidate Q satisfies

\[
 \max_h |q_Q(h)-q_P(h)|\le 2C\varepsilon/\eta.
\]

This uses valid candidate laws; it does not treat an unconstrained affine estimate with potentially negative atoms as a probability law. The confidence set is a mathematical feasibility set, not an implemented numerical solver. A true split with gap greater than 4Cε/η is retained on coverage. If every nonconstant suffix has a gap at least γ and 2Cε/η<γ/4, thresholding any retained candidate's transitions at range γ/2 recovers the true minimal tree. Its excess Brier prediction loss is at most (2Cε/η)².

**Budget comparison.** For a target transition radius u, substituting ε=ηu/(2C) gives the following conservative sufficient budgets. The paired formula is for 0<u≤1/2; both use 0<δ<1.

| Depth-four design | Sufficient independent blocks | Retained measurements |
| --- | --- | --- |
| Paired adjacent powers | N ≥ ceil(ln(32/δ)/(ηu²)) | 4N scalar powers, grouped by block |
| Six unpaired views | n ≥ ceil(2C² ln(16/δ)/(η²u²)) per view, N=6n | N scalar powers, one per block |

At a common total block budget N, the unpaired design can allocate floor(N/6) blocks per view; the paired design uses all N blocks. Both require a maximum span of five events. The unpaired bound loses additional factors through the conservative inversion estimate and division by rare histories. These are guarantees for the specified estimators, not matching minimax bounds or a proof of optimal measurement allocation.

**Collision and prediction tolerance.** Under a common observed law, if a supported-depth set covers each of depths 0 and d with error at most δ, it contains both with probability at least 1−2δ. Applied to the existing all-depth rectangular collision, additional independent records of those same views cannot resolve the depth ambiguity uniformly.

Depth ambiguity need not prevent accurate prediction. For every positive tolerance τ, the existing depth-four family contains a positive exact-depth-four alternative with the same rectangular Fourier laws as the fair source, whose three-symbol suffix predictor is 1/2 and whose excess Brier loss is at most τ. The exact loss along the family parameter (a,0) is a²/2. This is an existence result for arbitrarily weak deep dependence, not a uniform shallow-predictor guarantee over all deep sources. Brier loss is in squared-probability units, not bits.

**Relation to observed-sequence results.** Kontoyiannis studies prior predictive likelihood and posterior consistency for observed discrete time series; Galves and Leonardi give nonuniform exponential bounds and consistency for an empirical Context estimator, allowing unbounded trees. Their abstract-level scope differs from the independent, transformed block experiment proved here. Here the power observation must first identify the relevant law, and reconstruction incurs its own stability and history-mass factors. The present results do not transfer either paper's single-series guarantees to overlapping power windows. [Kontoyiannis](https://arxiv.org/abs/2211.02676v4) · [Galves and Leonardi](https://arxiv.org/abs/0710.5900v3)

**Declarations.** Namespace `TimeFrequency.ContextTree.DepthEvidence`.

- [MemoryDepthConfidence.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthConfidence.lean): supported depths, empty outcome, coverage, false/missed split events, common-law ambiguity.
- [MemoryDepthReconstruction.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthReconstruction.lean): the eight-summary inversion and exact source reconstruction.
- [MemoryDepthStability.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthStability.lean): integer row certificate, affine stability, ratio and source-kernel bounds.
- [MemoryDepthEmpirical.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthEmpirical.lean): bounded-statistic concentration and separate-view rounded budget.
- [MemoryDepthUnpaired.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthUnpaired.lean): physical statistics, expectations, six-view allocation and simultaneous error.
- [MemoryDepthPrediction.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthPrediction.lean): transition-error prediction bound and arbitrarily weak deep example.
- [MemoryDepthPowerConfidence.lean](../../../../../lean/Enthymemata/TimeFrequency/MemoryDepthPowerConfidence.lean): valid-source confidence sets, split errors, candidate tree and prediction guarantees.

**Verification.** Full `lake build` passed (3816 jobs), BelowNormal priority with `LEAN_NUM_THREADS=2`. Section 3 has 52 public axiom listings; all 143 across sections 1–3 use only subsets of `[propext, Classical.choice, Quot.sound]`. The new modules emit no warnings and contain no proof holes or postulated axioms. The repository-wide taxonomy check `lean/scripts/meta-ci.ps1 -NoBuild` passes after registration of the queued predictive-information modules. All 17 memory-depth modules are PROOF-CLOSED. That structural check does not qualify the queued proofs; this build receipt covers the memory-depth sequence. All modules remain Enthymemata pending semantic review and promotion.

## Doc changes owed

The research overview now records observation budgets, sampling assumptions and units, and the distinction between pair discrimination, tree recovery and prediction tolerance. Its exact-operator-norm question remains explicitly open; the proved reconstruction constant is conservative. No engine or documentation-site changes are included in this slice.
