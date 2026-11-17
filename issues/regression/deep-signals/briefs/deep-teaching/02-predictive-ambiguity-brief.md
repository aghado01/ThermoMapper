# Predictive consequences of the depth-four measurement completion

Commissioned 2026-11-17 as question 2 of the [Deep Teaching research brief](../deep-teaching-audio-research-brief.md), following [measurement completion](01-depth-four-completion-brief.md).

## Information and loss

Training supplies exact, separately labeled population laws of the stated windows. Deployment supplies the signed length-four history h. The source has the uniform history law in the collision family. The target is the conditional probability of the next sign being positive. Evaluate the history-averaged excess squared probability loss, E_i(t)=Σ_h 2⁻⁴(t(h)−q_i(h))². For a Bernoulli target this is Brier excess when the score is (t−Y)².

## Proof obligations

1. Derive the conditional probability and exact squared prediction distance for two parameters in the stationary family proved in slice 1. Quantify the residual directions under each single-taper augmentation.
2. Prove the sharp minimax excess for two specified observationally identical sources, with the midpoint predictor attaining it. Instantiate the constants on the physical Fourier collision.
3. Characterize the family’s linear predictive functionals invariant under the omitted directions. Include average next-positive probability and prediction given shorter suffixes. Specify deployment history and training information separately.
4. Prove that both added tapers remove the population ambiguity throughout the positive depth-four class. Explain the resulting zero ambiguity term under exact summaries. Finite-data estimation and stability require further results.
5. Report declarations, builds, axioms and Doc changes owed. Formulate concrete subsequent stability, teacher-transfer and sampling questions from the established constants.

## Acceptance

The main result quantifies loss due to the actual stationary Fourier collision and what the completed design repairs. The generic two-point identity supports that result. Retaining ordered training recordings or pairing observations changes the information contract.
