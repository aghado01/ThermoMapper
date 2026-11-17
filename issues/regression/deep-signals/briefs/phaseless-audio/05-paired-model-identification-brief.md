# Slice 5 — identification from paired views of the same waveform

**Question.** If an utterance contributes a labeled tuple of window views, what can its joint population law determine about the waveform law? The [approach](../approach-brief.md) proposes same-signal views. [Tachella et al.](README.md) identify model sets under conditions on separate linear images. Here the datum is a joint law of fixed nonlinear measurements.

1. **State the population experiment.** Let C be a compact real waveform class closed under sign, K = C/{±1}, and M the joint continuous map induced on K by the labeled powers. The sign quotient can equivalently be represented by the compact set of outer products xxᵀ. Show that if M is injective on K, then it has a continuous inverse on its compact image in the Euclidean measurement space. Deduce that M preserves distinction between compact signal sets and that equality of pushed-forward Borel probability laws implies equality of the laws on K. Make the support or class assumption common to both candidate laws.
2. **Converse and ambiguity.** From a single collision at two distinct points of K, construct unequal point masses with equal joint observation law. On the signed waveform space, exhibit the ambiguity between point masses at x and −x whenever x ≠ 0. The target here is the quotient waveform law. Identifying a generator parameterization would require additional assumptions.
3. **Pairing witness.** Use x = (a,b), the labeled views M₁ = (|a|,|b|) and M₂ = (|a+b|,|a−b|), and two equiprobable models
   - X supported on {(1,2),(2,−1)}
   - Y supported on {(1,−2),(2,1)}.

   Compute all four labeled tuples. Each individual view has the same distribution under X and Y. Their joint laws differ, as do their quotient supports. Give an event whose joint probability is 1/2 under one and 0 under the other. Squaring the coordinates gives an equivalent observable witness.

**Done when** the paired-law implication, collision converse and finite witness are verified, with the observation regime stated. Record whether a particular pilot or proposed multiwindow family has met the injectivity premise using [slice 3](03-finite-window-identifiability-brief.md). The law theorem itself can close even when that application is open. Finite training data and a particular optimizer are separate questions.
