We're close to the end of the BASS port itself. Everything R's BASS can do is now built at the maths layer, and the last mathematical stage, SMC, is running. What remains is checking the work against the port's contract and measuring its performance, rather than building new capabilities.

**R's public surface (the completion contract, P01–P15)**

| R function | Status |
|---|---|
| `bass` (scalar and direct functional fitting, tempering, restart) | Built: stage 4, plus the F-M2 options and N14's verdict |
| `bassPCA`, `bassBasis` | Built: stage 5 |
| `sobol`, `sobolBasis` (with effect functions) | Built: stage 6, plus effect functions for basis fits, which R lacks |
| `calibrate.bassBasis` | Built: stage 7, with R's defects corrected |
| `predict.bass`, `predict.bassBasis` | Built: stage 5 |
| The three plot methods, two print and two summary methods | The data is built (stage 9). Rendering and text formatting sit above maths and wait for a consumer and viz v2's contracts |

**Beyond R:**
- tempering with a self-sizing ladder, and the evidence;
- the convergence verdict;
- the local-move options;
- SMC for fitting (C4) and for cut calibration (C5), both running now. R has no SMC. The hub tracks the SMC backends apart from capability completion.

**What closes the port after stage 8:**
1. **A completion audit.** Check what's built against the hub's own contract: P01–P15, its list of edge cases, and its benchmark list. Some edge cases I can't vouch for from memory. One is the categorical representation above 64 levels, which the hub's modernization table calls for and which I haven't checked.
2. **Performance qualification.** The hub asks for benchmarks by case count, predictors, active terms, interaction order, response rank and category count, plus memory and error per unit time. That work isn't done. It folds in:
   - the parked efficiency pass;
   - N15, the `DiagonalPlusLowRank` cost;
   - the run-length defaults, which D10 made a budget but which still have no converged reading to size them from.
3. **Rendering and formatting above maths.** These are the top half of the plot and print methods. They're outside maths and wait for a consumer, so I'd treat them as part of viz rather than as a port blocker.

**Open, but not blocking completion:**
- Friedman's worse modes. C4's qualification may answer this: if SMC populations hold the better representation where MCMC chains don't, a later pre-registered pass would weigh SMC as the default for structured fits.
- F-R1, the deterministic-simulator profile, which needs the ridge prior built first.
- OQ2, exact reference rungs on long ladders.

**After BASS:**
- the notation unification pass, deferred until after BASS and its extensions;
- the extensions brief, E1–E15: GBMARS noise, KHAOS, ABMRS, active subspaces, cross-engine evidence, semi-modular calibration, calibrated Gibbs intervals;
- the BARS revisits you mentioned (E15, and E1's BARS side);
- on the SPC side, Thread 7, and the Domany-battery comparison of temperature rules.

My recommendation for after stage 8 is the completion audit. It's read-only, mostly mechanical, and turns "close to done" into a concrete list.