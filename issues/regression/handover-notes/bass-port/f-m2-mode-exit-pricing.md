# F-M2: pricing the mode exits, P1–P4 (chip report pointer)

Date: 2026-10-29. The full report is in `issues/regression/bass/mode-exit-pricing-brief.md`, section "Reports" ("Mode exits priced: P1–P4").

Commits on `main`:
- `17f8b5f`: the helper additions in `BassMixing.cs`, and P1.
- `108cefe`: P2. `b735551`: P3. `bc1e6da`: P4. `811360e`: P4's header fix.
- `f90c3cf`: the report.

There is no production change and no internal accessor, and `BassSweep.cs` was not touched. Regression is 838/838 with the switch off. M2–M4 reread unchanged.

The readings, in short (the detail and each reading's strength are in the report):
- *Fixed states.* Friedman chain 3's burned state lies below its whole sampling phase, so it was read at its most frequent structure. Phase 0's M2 readings for that chain came from that edge state.
- *P1:* no pair death in any of the 12 states costs under 57 nats, and every pair costs more than its cheaper single. Co-adapted pairs are substitutes, not jointly redundant.
- *P2:*
  - The prior and proposal ratios cancel to 1e-16.
  - No cheap fixed-M exit leaves Friedman's worse modes 0 and 3.
  - The mixed worse chain 2 and the functional worse chain 3 each have a near-free two-step path: a prune or variable change, then a death, ending at the better modes' level.
- *P3:* against the full state's own knot optimum, which is the needed control, Friedman's x1·x2 block re-represents in every chain, and so do the functional fit's worse chains. The endpoints are 137–1,237 knot steps away on Friedman.
- *P4:* the band's rungs are frozen for deaths. The median term's death α is at most 1e-3 on 30–31 of 31 rungs per chain, and Δℓ_b tracks b·Δℓ₁.
- *Design points left open for the user:* the part and move selection under several parts, the "part emptied" prune, and the variable change with no unused variable.
