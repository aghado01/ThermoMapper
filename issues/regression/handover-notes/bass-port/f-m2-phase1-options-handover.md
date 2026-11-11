# F-M2 Phase 1: the options chip's hand-over (2026-10-29)

The full report is prepended under "Reports" in `issues/regression/bass/local-move-options-brief.md` (repo `D:\aghado01\ThermoMapper`).

- Commits on `main`, none pushed: `cf431c0` (Step 0), `7a08c37` (Step 1), `afafe03` (Step 2), `368d7c3` (Step 3), `a99c02f` (Step 4), `c167af3` (Step 5), `aaa8034` (review fixes), `856ecfd` (the smoke harness), `8217eba` (a pinned JSON shape), `9a0d277` (the report). Suites at `8217eba`: Regression 1,030, Samplers 194, BASS oracle 54, all passing. `174b8db`, interleaved, is the user's own docs commit.
- The options are each null by default:
  - `BassSampling.LocalSteps(0.5, 0.02)`;
  - `BassSampling.FactorMoves(1, 1, 1)`, thirds of the change's slot, with move keys `prune`, `grow` and `variable`;
  - `TemperingSpec.Split(0.9)`, self-sizing ladders only;
  - each echoed by `BassResolved` (`StepWindows`, `FactorShares`, `SplitThreshold`).
- Exactness: detailed balance by enumeration on |F| = 2,263, 2,918 (categorical, both measures) and 593 (functional). The worst log gap was 5.7e-14 against a tolerance of 1e-9. Negative controls fail as they should. Visit frequencies and the constant-likelihood laws pass with every option on.
- A part's K = 1 → 0 cancels exactly, because each part's order is uniform on {0..cap}. The cancellation is asserted over every term of the multi-part fixtures.
- Bit-identity holds with the options off: `BassBitIdentityTests` (8 hashes recorded before any production edit), plus the older S8 and resize pins.
- Choices for the main session to confirm:
  - A null move is a ratio of −∞: attempted and rejected, drawing the chain's uniform, with the options on only.
  - A prune of a term's only factor is refused with no `BassRefusal` reason.
  - 4b splits in place of sizing while a pair is steep, and the last round still equalizes.
  - `ResetTallies()` keeps its signature.
  - With the family on, `Support` counts the family's candidate terms.
- Smoke reading (one seed, read only): the local steps accept 0.83 of admitted proposals on Friedman and mixed, and 0.40 on functional. The family accepts at most 0.03. R̂ and the chains' bands move both ways by option. With 4b the ideal round-trip rate is 0.0091 against 0.0014; there are no observed round trips in either run, and steep pairs total 5 in both.
