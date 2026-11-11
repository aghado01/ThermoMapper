# BGP-4 Lean chip (UV-1c, UV-4a, UV-6, then UV-5 and the citation checks): status

**2026-10-24.** Commits on main: f08fabe (UV-1c), 7230805 (UV-4a), dcbe6a3 (UV-6), 9f9df85 (protolemma + report).

**2026-10-25 follow-up.** The commits are in the chip brief's report; full detail is there, in
lean/Protolemmata/BARS/bgp-4-open-pieces-chip-brief.md under "Report".
- **UV-5 is proved.** The new BgpGaussianTails.lean holds the Gaussian tails, and exists_tests covers three cases.
  No `sorry` is left in the BGP-4 modules.
- **The TWCD2025 checks are confirmed against Appendix A.1:** (59) uniform over the amplitude box (Part I is monotone
  in (A2)'s constants), (60) at any exponent (the constants in eq. 80), and Rₙ polynomial (eq. 75 with h(0) = 1).
- **The GvdV07 check:** the printed Theorem 4 has one rate and concludes for every Mₙ → ∞, not for a fixed M.
  - The fixed-M two-rate form is Ghosal & van der Vaart (2001, Thm 2.1), as vdVvZ09 says.
  - In Lean, the axiom gvdv07_contraction is renamed fixed_constant_contraction and re-cited.
  - The printed Theorem 1 and Lemma 1 are new axioms, giving a second composed theorem for every Mₙ → ∞.
- **Prior art:** vdVvZ08b Theorem 3.3 is the unknown-σ fixed-design GP result. Its proof is a two-line sketch.
- **Environment:** the machine was paging (63 GB committed on 31.5 GB of RAM), so Lean loads took 8–10 minutes.
- **Touched only lean/.** issues/gemini-pkwang-cleanup.md (deleted) and issues/spc/ are someone else's.

Committed 2026-10-25: d971edb (tails), 17aae2c (theorem module), 6ead543 (protolemma + brief report). Full meta-ci PASS; all five Bgp modules PROOF-CLOSED.
