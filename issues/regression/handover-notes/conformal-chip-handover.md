# Conformal chip (CB-1 to CB-4): hand-over

Done in commit 99b464f on main. The record is the report under "Reports" in `lean/Protolemmata/Conformal/conformal-chip-brief.md`. The statements and decisions are in `conformal-lemmas.md` beside it, and the proofs in `lean/Enthymemata/ConformalValidity.lean` (PROOF-CLOSED).

For the GBASS sessions, about the pointwise contract (CB-4):
- The conformal layer reads every unit, the candidate included, by one density.
- The default is the density with every observation-level latent integrated. That is OQ17's scope rule, which GBASS's pointwise content already follows.
- Data units at their retained latents with the candidate integrated undercovers. It never covers more than the symmetric latent score, and it reaches coverage 0 in a two-unit counterexample. So the layer should refuse a latent-conditioned matrix for a latent-noise model.

E16's Hoff bullet was corrected on 2026-11-06 at the user's call. It now separates Hoff's leave-one-out score from the in-sample add-one-in step.
