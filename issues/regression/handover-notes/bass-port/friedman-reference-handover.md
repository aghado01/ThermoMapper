# Friedman: path from a pilot reference — hand-over, 2026-10-29

Full report: `issues/regression/bass/friedman-reference-brief.md`, under Reports ("The reference path measured").

- Commits on `main`, none pushed: `abdcd4d` (harness `tests/regression/BassFriedmanReferenceTests.cs`, behind BASS_QUALIFICATION), `99b9cb1` (report). No production change; no seam needed (`ReversibleJumpChain` and `BassTauBlock` take any target, so the kernel is C4's sweep with γ_b in the target's place).
- Checks pass: kernel vs γ_b by enumeration on fixture A (b 0.3, 0.7); q's draws vs density; exact ln Z on B and C by this harness's enumeration (−28.5886, −30.9026, matching the oracle), SMC from q within 4 errors; ε = 1 vs C4's from-prior SMC on B (gaps 0.001–0.006).
- `BassTarget`'s structure prior is D02's count weight w, unnormalized both ways: from a non-prior start, ln p(y) = normalizer − ln C − ln Z_F (ln C 0.9531, ln Z_F −0.194 on Friedman).
- Friedman: at 4,000 particles every arm (all pilots / no SMC pilot × ε 0.05 / 0.5) agrees within a nat (sd 0.26–0.58), on −886.4 to −887.5. That is the evidence of e's own structures (−887.15 with no move), and 22–23 nats below the 99.9% lower bound (−864.6) that R2's from-prior replicates put on ln p(y) via Markov on unbiased Ẑ. Population holds 100% of its weight on e's structures through b = 0.75 (e's height over ε·p₀ ~120 nats); leaves only in the last stages; no families collapse (20–174 at b = 1); A at 65–98%; no a1/a2.
- From-prior pilots reproduced R2 bit-for-bit on the current code (ridge chip's commits included).
- `artifacts\chip-reference\` deleted.
