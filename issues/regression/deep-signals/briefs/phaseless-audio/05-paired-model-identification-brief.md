# Slice 5 — identification from paired views of the same waveform

**Question.** If an utterance contributes a labeled tuple of window views, what can its joint population law determine about the waveform law? The [approach](../approach-brief.md) proposes same-signal views. [Tachella et al.](README.md) identify model sets under conditions on separate linear images. Here the datum is a joint law of fixed nonlinear measurements.

1. **State the population experiment.** Let C be a compact real waveform class closed under sign, K = C/{±1}, and M the joint continuous map induced on K by the labeled powers. The sign quotient can equivalently be represented by the compact set of outer products xxᵀ. Show that if M is injective on K, then it has a continuous inverse on its compact image in the Euclidean measurement space. Deduce that M preserves distinction between compact signal sets and that equality of pushed-forward Borel probability laws implies equality of the laws on K. Make the support or class assumption common to both candidate laws.
2. **Converse and ambiguity.** From a single collision at two distinct points of K, construct unequal point masses with equal joint observation law. On the signed waveform space, exhibit the ambiguity between point masses at x and −x whenever x ≠ 0. The target here is the quotient waveform law. Identifying a generator parameterization would require additional assumptions.
3. **Pairing witness.** Use x = (a,b), the labeled views M₁ = (|a|,|b|) and M₂ = (|a+b|,|a−b|), and two equiprobable models
   - X supported on {(1,2),(2,−1)}
   - Y supported on {(1,−2),(2,1)}.

   Compute all four labeled tuples. Each individual view has the same distribution under X and Y. Their joint laws differ, as do their quotient supports. Give an event whose joint probability is 1/2 under one and 0 under the other. Squaring the coordinates gives an equivalent observable witness.

**Done when** the paired-law implication, collision converse and finite witness are verified, with the observation regime stated. Record whether a particular pilot or proposed multiwindow family has met the injectivity premise using [slice 3](03-finite-window-identifiability-brief.md). The law theorem itself can close even when that application is open. Finite training data and a particular optimizer are separate questions.

## Verified result — 2026-11-16

Let both candidate laws live on the same compact class of real waveforms modulo global sign. Represent this class by its outer products and observe the labeled tuple of powers from the same waveform. `outerClass_compact` proves compactness. `lifted_injective_iff` expresses injectivity of the induced continuous map as pair separation on the waveform class. Under that premise, `imageHomeomorph` and `inverse_recovers` construct a continuous inverse on the measured image. `image_sets_eq_iff` distinguishes all subsets, including compact subsets. `paired_power_laws` proves equality of the two pushed-forward Borel measures exactly when the two laws on the outer-product class agree. This includes probability measures. [PairedViewLaws.lean](../../../../../lean/Enthymemata/TimeFrequency/PairedViewLaws.lean)

A collision at distinct quotient points gives unequal Dirac laws with identical paired measurements (`collision_diracs`). On signed waveforms, `signed_diracs` gives the two laws concentrated at x and −x for any nonzero x. The full-grid translated Hann meets the separation premise on any common compact class whose nonzero members have connected support in its lag graph. Each single-frame pilot fails separation on unrestricted waveforms, and separation on the two restricted synthesis families is open. Generator parameterization and finite-sample estimation require further assumptions. [Slice 3](03-finite-window-identifiability-brief.md)

| Model | Waveform | Singleton magnitudes | DC and Nyquist magnitudes |
| --- | --- | --- | --- |
| X | (1,2) | (1,2) | (3,1) |
| X | (2,−1) | (2,1) | (1,3) |
| Y | (1,−2) | (1,2) | (1,3) |
| Y | (2,1) | (2,1) | (3,1) |

Each listed atom has mass 1/2 within its model. The event that the paired tuple equals ((1,2),(3,1)) has probabilities 1/2 and 0. Each labeled view separately has the same probability law under the two models. The quotient supports differ because neither sign of either representative in Y equals (1,2). `views_squared` connects these magnitudes to the actual singleton and length-two Fourier powers. `four_tuples`, `equal_single_views`, `different_joint_laws` and `different_quotient_supports` check the witness.

### Verification

`lake build +Enthymemata.TimeFrequency.PairedViewLaws` passed. All 15 public theorem axiom listings used only standard Lean axioms.

### Doc changes owed

State the common signal class and the pairing regime alongside a population-law claim. Include the four tuples as an example where losing the pairing changes what can be learned.
