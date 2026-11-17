# Conformal Bayes for conversation audio — proof sequence

Posed 2026-11-16. The [audio approach](../approach-brief.md) proposes using one view of a recording to predict another. Start with a coverage statement that names the observable and the sampling unit. Use it to study how ambiguity limits the precision of a prediction. The last three slices address transport, decisions and monitoring, each with a different probability statement.

| Order | Assignment | Depends on |
| --- | --- | --- |
| 1 | [A held-back spectral band on a new conversation](01-conversation-band-brief.md). Prove simultaneous bin coverage and distinguish a randomly sampled turn from every turn. | Existing split ranks and group-score contracts |
| 2 | [Set size under a known phase ambiguity](02-ambiguity-and-set-size-brief.md). Compute a sharp lower bound and exhibit the gain from one informative view. | The observable coverage target in 1 and the verified [gap witness](../phaseless-audio/01-two-sample-recovery-brief.md) |
| 3 | [Transporting synthetic latent coverage](03-latent-transport-brief.md). Derive the weighted guarantee and break it by changing the latent law at the same observed input. | Sampling unit from 1, witness from 2, existing weighted ranks |
| 4 | [Early starts on recorded turns](04-recorded-turn-risk-brief.md). Bound the expected fraction of eligible turns interrupted in a new conversation. | Conversation sampling from 1 and the existing conformal risk rule |
| 5 | [An anytime alarm on a conversation stream](05-conversation-alarm-brief.md). Connect randomized sequential ranks to the exact product-law premise of Ville's inequality. | Permutation averaging from 1 and the existing sequential rank and capital results |

The order after slice 2 sets reading priority. Slices 4 and 5 can use the conversation result directly. All five briefs pose statements for proof. Keep empirical claims about real conversations, simulators and shifts identifiable as hypotheses in the reports.

**Proof and reporting.** Reuse the declarations in [ConformalValidity.lean](../../../../../lean/Enthymemata/ConformalValidity.lean) and follow the [Lean harness conventions](../../../../../lean/README.md). Append to each brief the exact conclusion and hypotheses, its witness or counterexample, the declarations and commit, and the verification commands and results. If a proposed implication fails, record the counterexample and the corrected claim. List any documentation work under **Doc changes owed**.

**Local reading.** The versioned papers below have TeX sources under `codex-scientiae/supellex-private/gauntlet/conformal-bayes`. Each brief points to the paper that owns its statistical argument. The project-specific witnesses use the phaseless audio proofs.
