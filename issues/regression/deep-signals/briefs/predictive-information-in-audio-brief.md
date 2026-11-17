# What predictive information does a learned audio representation retain?

Research proposal and commissioned first-slice proof report · 2026-11-17

**For a specified future observable in conversational speech, how much does earlier audio improve prediction once the recent audio is known? Which self-supervised objectives make a learned representation preserve that improvement, and how many independent conversations does it take to learn and demonstrate it?**

This question guides a choice of context, measurements, and training targets for a native audio model. Fix the prediction task before choosing an encoder. An instructive first target is the *joint* future activity of two speakers. A second target could be the quantized energy in a held-back acoustic band. Both have outcomes one can read from a recording under a declared measurement rule. A foundation model would need a family of such targets and horizons. A guarantee applies to the family on which it is established.

## An experiment with a clock and an observable

For the joint-activity example, take separate speaker channels and choose a time \(t\). Let \(R\) be the last 200 ms of both channels, \(O\) the preceding two seconds, and \(Y=(A,B)\) two bits indicating whether each channel meets a fixed activity criterion in the 200 ms interval starting 500 ms after \(t\). These durations illustrate a question one could ask. Changing the prediction horizon, history span, channel access or activity rule defines a different experiment. If an activity detector supplies the bits, its output defines the observable. On mixed audio, attributing activity to individual speakers requires an additional observation model.

Write \(X=(R,O)\). Train an encoder on other recordings, freeze it, and form \(Z=f(X)\) on a new conversation. The deployment predictor in this experiment receives \((R,Z)\). This makes the contribution of the older segment measurable while giving every predictor the same recent audio. For an encoder that must replace the *whole* history, evaluate a separate predictor receiving \(Z\) alone. Use causal feature extraction, including window alignment, so every encoder input is available by \(t\).

## Three quantities in bits per prediction

For a finite target, use base-two log loss under the deployment distribution. Conditional entropy is the loss of the best possible probabilistic predictor with the specified inputs. With finite expectations and a fixed encoder,

| Quantity | Population expression | Interpretation |
| --- | --- | --- |
| Gain available in older audio | \(G_X=H(Y\mid R)-H(Y\mid X)=I(Y;O\mid R)\) | How much observing the older segment can improve prediction. |
| Gain preserved by the encoder | \(G_Z=H(Y\mid R)-H(Y\mid R,Z)=I(Y;Z\mid R)\) | How much of that improvement a predictor can obtain from the representation and recent audio. |
| Information discarded | \(G_X-G_Z=I(Y;O\mid R,Z)\) | The expected divergence in bits between the future laws given \(X\) and given \((R,Z)\). |

A representation is predictively sufficient for this target exactly when \(P(Y\mid X)=P(Y\mid R,Z)\) almost surely. An approximate theorem should bound the discarded bits by a stated tolerance. That theorem concerns the conditional law, so it permits multiple latent parameterizations and works for a recurrent state, a codec, or a generator-based inference procedure. State a resource constraint, such as a finite code alphabet or transmission budget, whenever claiming that compression preserves the gain. Copying the entire history gives an exact unconstrained baseline.

A fitted decoder \(q\) incurs a further cost. Its excess log loss over the full-history optimum decomposes as

\[
L(q)-H(Y\mid X)
= (G_X-G_Z)
+\mathbb E\,\mathrm{KL}_2\!\left(P(Y\mid R,Z)\,\Vert\,q(\cdot\mid R,Z)\right).
\]

The first term measures the encoder. The second measures how well the decoder reads it. This distinction matters when a finite model performs worse with longer context even though additional history can only improve the population optimum. Bialek, Nemenman and Tishby study information shared by past and future. Conditioning on \(R\) and fixing \(Y\) gives the task-specific increment here. [Predictability, complexity and learning](https://arxiv.org/abs/physics/0007070)

## A context tree gives one exact case

Take \(Y\) to be the next sign, \(R\) the most recent \(d-1\) signs, and \(O\) the preceding sign. For every \(d\ge4\), the project has constructed a fair IID source and a positive stationary sign-symmetric binary source of exact memory \(d\) whose separate rectangular-power laws agree through width \(d+1\). The latter has two positive-probability histories with the same \(R\) and different next-sign probabilities. The direct information-theoretic consequence is \(G_X=0\) for the fair source and \(G_X>0\) for the other. Thus even exact knowledge of these training summaries leaves the *presence* of incremental predictive information undecidable across the pair. The existing Lean statements establish the collision and transition difference. The first-slice result below formalizes this zero/positive distinction in bits. [Stationary sources and Fourier observations](context-tree-memory-identifiability-brief.md#verified-result--2026-11-16)

At depth four, the completed proofs go further. Two specified four-site taper laws recover the source throughout the stated class. Every single four-site indicator taper admits a collision. With signed deployment histories, the sharp supremal minimax excess Brier loss on the constructed family matching the fair observations is \(1/8\), including after one taper. Brier loss has squared-probability units. The corresponding value in bits requires its own calculation. The training summaries in these proofs describe populations, whereas \(Z\) describes one deployed history. They enter different parts of a learning theorem. [Completion report](deep-teaching/01-depth-four-completion-brief.md#result-and-proof-report--2026-11-17) · [Prediction report](deep-teaching/02-predictive-ambiguity-brief.md#result-and-proof-report--2026-11-17) · [Finite-sample evidence follow-up](deep-teaching/evidence-for-memory-depth-brief.md)

## First representation theorem

**Can low self-supervised training loss certify that *any encoder achieving it* preserves the joint future activity law?** Use \(Y=(A,B)\) from the audio experiment. Compare heads predicting two observable bits, \(A\) and \(B\), with three heads predicting \(A,B,AB\). The product is observed on the same future interval and records simultaneous activity. All heads receive \((R,Z)\). This comparison changes one training target while holding the desired future event and deployment inputs fixed.

**Obstruction.** Construct a finite law with constant \(R\) and a fair older bit \(O\). Conditional on \(O=0\), the future is equally likely to be \(00\) or \(11\). Conditional on \(O=1\), it is equally likely to be \(01\) or \(10\). Both individual activity probabilities are always one half. A constant encoder therefore attains full-history Bayes risk on both marginal training heads. Yet the joint future determines \(O\), so this encoder discards exactly one bit per prediction. This is the lower-bound witness for the objective formalized in the first-slice report. It specifies a possible dependency between observable events, without assuming a context tree.

The witness admits a waveform realization. Encode \(O\) by silence or a fixed pulse in the older segment, give every example the same recent segment, and render each future bit as silence or a pulse on its channel. Choose the energy threshold strictly between their energies. A two-value code storing \(O\) attains full-history prediction of the joint law, so a one-bit representation budget suffices. The marginal objective also admits the constant optimum under that same budget. This isolates the choice of training target from the encoder's capacity. All three labels can be computed from the same recorded interval.

**Positive statement.** Write \(u=P(A=1\mid X)\), \(v=P(B=1\mid X)\) and \(w=P(AB=1\mid X)\). The four probabilities of \(00,10,01,11\) are \((1-u-v+w,u-w,v-w,w)\). The same inversion holds given \((R,Z)\). Prove that at full-history Bayes risk on the three heads, the learned representation retains the entire conditional law of \(Y\).

For an approximate version, let \(\eta\) be the **sum of the three population squared-loss excesses above their respective full-history Bayes risks**. Conditional-expectation projection bounds the sum of mean squared errors between the three full and encoded conditional means by \(\eta\), even for imperfect heads. The linear part of the displayed inversion has squared Frobenius norm eight. The first-slice proof establishes the bounds

\[
\mathbb E\|P(Y\mid X)-P(Y\mid R,Z)\|_2^2\le8\eta,\qquad
G_X-G_Z\le\frac{8\eta}{\beta\ln 2}
\]

for the second inequality when all four encoded conditional probabilities are at least \(\beta>0\) almost surely. The log bound uses the encoded conditional law in the denominator. The squared bound also applies on the simplex boundary. These constants are sufficient targets. Report any improvement or counterexample with its hypotheses. The result states which extra self-supervision forces a predictive distinction to survive, with a quantitative penalty for imperfect learning.

## From an objective to a learned audio representation

1. **Which training views determine the target?** For a finite future alphabet, each observable head measures a linear functional of its conditional law. Characterize the differences between admissible future laws that all heads erase. Give a pair of histories with the same recent context realizing one erased difference, or prove a lower bound on how much the head means change for each admissible difference. That lower bound controls the analogue of the factor eight. Apply the criterion to an acoustic target, such as joint occupancy of two future frequency bands. This asks which *events* the measurement objective should retain.
2. **Can the chosen encoder attain the required excess under a budget?** Put the encoder and heads in a specified hypothesis class. Bound its best three-head excess at the chosen code size and add training and estimation errors. For a generator-based teacher, specify the predictions available from causal audio and bound their conditional bias or mismatch against the observable head targets. With known recording transformations, an invariant code is justified for a target only if collapsing the corresponding histories preserves its conditional law given \(R\). A collision with different future laws supplies a lower bound on the cost of that invariance.
3. **How much unlabeled data suffices?** Start with independent conversations, a declared sampling rule for \(t\), bounded heads, and empirical risk averaged per conversation. If the chosen class contains an encoder and heads within \(\eta_{\rm app}\) of full-history Bayes risk, a uniform empirical-to-population error \(u_n\) and empirical optimization error \(\eta_{\rm opt}\) give a candidate population excess \(\eta_{\rm app}+2u_n+\eta_{\rm opt}\). Derive \(u_n\) for an explicit class and substitute into the retention bound. Count data used to learn the prior or teacher. Compare against direct prediction of the same joint observable with the same recordings, causal inputs and representation budget. A claimed saving needs a smaller sufficient budget or a proved risk advantage for a named comparison.
4. **What can be demonstrated on new conversations?** Freeze the encoder and decoders, measure the paired difference in log scores on held-out conversations, and state the confidence level. Specify smoothing or integrability for the log scores and the split unit when speakers recur. This estimates the gain of the fitted predictors. Inferring the population information gain also requires control of the decoders' excess losses. Repeat with fixed targets at different horizons to learn which history and timing resolution each target uses.

## Relation to speech research and Lean scope

Laskowski and Shriberg found that older activity could be represented at coarser timing resolution, with activity seven seconds back still useful for predicting speech onset on their two conversational corpora. Robertson and Dunbar found that phone discriminability peaked around 40 ms of preceding context in their contrastive pretraining setup. Ekstedt and Skantze study the dependencies among future voice activity events with a self-supervised target. These studies motivate specifying the observable, horizon and representation each time a context length is reported. [Conversational history compression](https://www.cs.cmu.edu/~kornel/pubs/0004937.pdf) · [Acoustic context in pretraining](https://arxiv.org/abs/2312.01515) · [Voice Activity Projection](https://arxiv.org/abs/2205.09812)

For a first Lean slice, take finite histories and targets, prove the one-bit obstruction and the three-head retention bound, then instantiate the criterion on a declared measurement of a future audio block. Quantify over eligible encoders and heads. State which conclusions hold for every encoder, which depend on an optimization guarantee, and which still require evidence from recordings. The entropy identities supply the bookkeeping. The deliverable that guides model design is the change of observable targets that forces preservation, together with a bound at the loss the learner can attain. [Deep Teaching research brief](deep-teaching-audio-research-brief.md) · [Lean harness](../../../../lean/README.md)

## First-slice result and proof report — 2026-11-17

**Model and quantifiers.** Use a finite retained history support X with a strictly positive history law, and an arbitrary nonnegative conditional distribution p_x on the four future outcomes 00, 10, 01, 11. Zero future probabilities are allowed. A frozen deployment feature f maps X into a finite used alphabet and is surjective onto that alphabet. For the proposed audio experiment, instantiate f(x)=(R(x),Z(x)); the theorem does not require that a representation replace the recent audio. Surjectivity only removes unused code values. The encoded conditional law is defined by the actual weighted fibers of f, rather than supplied as a premise.

The results quantify over every such frozen encoder and every three real-valued decoder heads. There is no premise that the heads are Bayes-optimal, no learning algorithm and no optimization guarantee. The code-size comparison below is specific to the one-bit witness; the general theorem is conditional on the excess the chosen encoder and heads attain.

**A marginal optimum can discard one bit.** Let the older bit O be fair and recent context constant. Given O=0, the future is 00 or 11 with probability 1/2 each; given O=1, it is 10 or 01 with probability 1/2 each. Both marginal activity means equal 1/2 at every history. A constant encoder with constant heads 1/2 therefore attains the full-history Bayes risk on both marginal tasks, pointwise in O. Its encoded joint law is uniform on the four outcomes.

Its expected divergence from the full-history joint law is exactly **one bit per prediction**. The conditional-entropy difference is also proved to be one. A two-value identity code storing O reproduces the full conditional joint law. Both this code and the constant code meet an alphabet-size-at-most-two constraint. Thus the marginal objective admits a maximally lossy solution for the available one-bit joint distinction even though that distinction fits the stated budget.

**Observable targets and waveform realization.** The three heads predict A, B and A AND B, all on the same future outcome; their means are u=p10+p11, v=p01+p11 and w=p11. Their equality to the expectations of the three explicit finite observables is proved. The inversion is

\[
 (p_{00},p_{10},p_{01},p_{11})=(1-u-v+w,u-w,v-w,w).
\]

The waveform witness uses finite sampled blocks: silence or a fixed pulse encodes the older bit and each future channel's activity; the recent block is common silence. Squared sample energy and a threshold strictly between zero and pulse energy recover each bit exactly. A one-sample unit pulse with threshold 1/2 is an explicit realization. Past, recent and future blocks have ordered slot timestamps −2, −1 and +1 relative to cutoff zero. These are dimensionless event slots, not a calibrated speech-time or detector model. The two future channels share the future slot, so their conjunction needs no additional recording or future interval. The construction exhibits a possible audio dependency; it is not a claim about conversational speech frequencies.

**Three heads force joint retention.** Define η as the sum of the three population squared-loss excesses above their respective full-history Bayes risks. The encoded conditional means are weighted projections onto the fibers of f. The proved projection identity separates representation error from the additional error of the supplied heads. It implies

\[
 \mathbb E\|p_X-q_{f(X)}\|_2^2\le8\eta .
\]

The constant eight is a sufficient bound, not an optimal operator norm. At η=0, the entire conditional joint law agrees at every retained history. The squared bound remains valid on the probability-simplex boundary and accommodates imperfect heads.

If every encoded probability q_z(y) is at least β>0, then

\[
 \mathbb E\,\mathrm{KL}_2(p_X\Vert q_{f(X)})
 = H(Y\mid f(X))-H(Y\mid X)
 \le \frac{8\eta}{\beta\ln2}.
\]

The floor is on **q**, the encoded law in the divergence denominator, not on p. Full-history probabilities may vanish. The formal definitions use natural logarithms internally and divide by ln 2. The entropy identity ties the divergence bound to discarded bits. For f=(R,Z), this is precisely the difference between the older-audio gain available from X and the gain retained by (R,Z). No floor-free bound in bits is claimed.

**All-depth context-tree case.** With uniform histories from the existing stationary construction at every depth d=m+4, take recent context to be the tail of the history and the target to be the next sign. Its transition probabilities lie strictly between zero and one. The existing exact-memory witness provides two histories with the same recent context and different next-sign probabilities. A generic finite projection theorem proves that such a distinction yields strictly positive expected conditional divergence in bits. The fair source has zero gain.

The combined declaration retains the actual equality of all separate rectangular Fourier-power laws through width d+1 alongside the zero/positive information conclusions. It therefore upgrades the earlier transition-based obstruction to an information statement. It does not calculate a universal numerical bit gap, and the earlier sharp Brier value 1/8 is not reused as a value in bits.

**Declarations.** Namespace `TimeFrequency.PredictiveActivity`.

- [PredictiveActivityCore.lean](../../../../lean/Enthymemata/TimeFrequency/PredictiveActivityCore.lean): actual observable means, exact inversion, injectivity and the factor-eight squared bound.
- [PredictiveActivityProjection.lean](../../../../lean/Enthymemata/TimeFrequency/PredictiveActivityProjection.lean): weighted code fibers, conditional means, orthogonality and the projection/error identities.
- [PredictiveActivityRetention.lean](../../../../lean/Enthymemata/TimeFrequency/PredictiveActivityRetention.lean): actual encoded law, head excess, approximate joint retention and exact retention at zero excess.
- [PredictiveActivityBits.lean](../../../../lean/Enthymemata/TimeFrequency/PredictiveActivityBits.lean): scalar divergence bound, encoded probability floor and discarded-bit bound.
- [PredictiveActivityEntropy.lean](../../../../lean/Enthymemata/TimeFrequency/PredictiveActivityEntropy.lean): weighted log identity and conditional entropy interpretation.
- [PredictiveActivityObstruction.lean](../../../../lean/Enthymemata/TimeFrequency/PredictiveActivityObstruction.lean): constant marginal optimum, exact one-bit loss, identity-code retention, code sizes and pulse realization.
- [PredictiveInformationPositive.lean](../../../../lean/Enthymemata/TimeFrequency/PredictiveInformationPositive.lean): divergence positivity and information lost by merging differing conditional probabilities.
- [PredictiveInformationContext.lean](../../../../lean/Enthymemata/TimeFrequency/PredictiveInformationContext.lean): actual all-depth transition law, fair/deep information and the rectangular-observation obstruction.

**Verification.** Full `lake build` passed (3824 jobs), BelowNormal priority with `LEAN_NUM_THREADS=2`. All 52 public axiom listings in these eight modules, and all 195 across this slice and the memory-depth evidence sequence, use only subsets of `[propext, Classical.choice, Quot.sound]` (including one empty axiom set). The new modules emit no warnings and contain no proof holes or postulated axioms. The taxonomy gate `lean/scripts/meta-ci.ps1 -NoBuild` passes, with all eight new modules PROOF-CLOSED and reached by the aggregate import closure. Modules remain Enthymemata pending semantic review and promotion.

**Scope remaining.** A finite-sample theorem for learning a constrained encoder, real-recording evidence, longer acoustic-band targets, continuous-time waveform analysis and optimal constants remain unproved. This first slice establishes what the specified population objective forces every eligible encoder to retain, conditional on its attained population loss.

## Doc changes owed

Update the approach's target-selection and teacher claims in the next documentation batch: individual future marginals can erase predictive dependence; adding their same-interval conjunction identifies a binary joint law; small three-head population excess gives the stated squared and floored bit guarantees. Preserve the distinction between a theorem about attained population loss and a guarantee that training attains it.
