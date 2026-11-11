# BASS completion audit: hand-over (2026-11-01)

Report: `issues/regression/bass-completion-report.md`. Matrix: `issues/regression/bass/completion-audit-matrix.md`. Probes: `tests/regression/BassCompletionAuditProbes.cs`. The brief's Reports section has the note for the main session.

- **Verdict:** not yet fulfilled in full. 215 rows: 183 Fulfilled, 12 limited, 17 left out by decision, 3 Gap. No wrong number found.
- **Contract gaps, ranked:** G-1 `other` (P05 names it), G-2 `func.var` on several output variables (skipped probe), G-3 `mcmc.use` on Sobol, G-4 `temp.ladder` as a list, G-5 `start.temper`, G-6 prediction's two switches, G-7 and G-8 the two open benchmark measures.
- **Open for the user:** three rulings in the report's "What the brief does not rule on"; they decide whether G-3 to G-6 stay gaps.
- **Suites:** Samplers 284/284, Regression 1,507 (1 skipped, G-2), BASS oracle 90/90, in the main checkout with `cb0400b` (another session's flat-histogram refusal) in the tree.
- **Re-trace:** 18 rows, one failed (`npart.func`), widened to seven moved rows, now covered by a probe.
- **Not run here:** LinAlg and Distributions suites; anything behind `BASS_QUALIFICATION` or `SEED_SWEEP`.

## Gaps closed (2026-11-02)

G-1 `0ab1e4d`, G-2 `4c34948`, G-3 `1364cbb` + `b844f53`, G-4 `3a83d0a`, G-6 `26befff`; G-5 left out as ruled. Closing docs `971c84b`. Verdicts F 192, FL 5, D 18, G 0. For what to check, read the brief's Reports section (`issues/regression/bass/completion-audit-brief.md`); G-7 and G-8 stay open.
