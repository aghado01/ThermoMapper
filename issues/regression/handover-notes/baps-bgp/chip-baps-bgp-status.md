# BAPS/BGP chip (B04, B05): working notes

Started 2026-10-19 at f5afb7a.

- While this chip runs, there may be uncommitted edits under src/maths/regression/spline/baps/**,
  src/maths/regression/bgp/**, and new test files tests/regression/{Baps,Bgp,PenalizedSpline,HierarchicalPSpline}*.cs.
  Each check lands in the same commit as its fix. Please don't stage these paths.
- Scratch probes run from a console project in the chip's own scratchpad, not in the shared test project.

## Status 2026-10-19 (chip)

Landed on main: 0ea0329 (R03 fix), 4745e00 (R10 support rejection), c2f0a1d (R06 protection test),
9f9a0b3 (R10 wording), 3d57bab (item 6: BAPS R-hat moments form; the sums-based ChainDiagnostics.RHat
overload now has no callers), plus the decision memo in the chip brief.
Waiting on the user: R02 target (memo recommends a uniform prior on edf) and R10 transform
(memo recommends instrument-owned similarity + response centering). The BGP transform would change
the EmpiricalBayesBandwidthPrior constructor (T_n read from the GP's own distances).
Interim suite (fact-harness, Maths.Regression.Tests, Release) at code revision 3d57bab: 105 passed, 0 failed
(artifacts/test-runs/TestHarness.Runner/20261019_183307417).

## Final 2026-10-19 (chip complete)

User decisions: R02 = A1 (uniform edf prior) + LambdaMedian for all modes; R10 = T2 + M1 (+None) + response centering.
Further commits: 1d1cada (LambdaMedian), 5fe826e (R02 edf prior), 1f10ba0 (BGP fitted transform; prior ctor is now (gp, ...)),
ed60885 (enum theory data), d0464ac (report appended to the chip brief).
Suite at ed60885: 117/117 (artifacts/test-runs/TestHarness.Runner/20261019_190012452). BGP_BENCH benchmark passes.
Requests for the main track are in the brief's report: delete the sums-based ChainDiagnostics.RHat (no callers left, and its doc
line is stale); a new BandCholesky call site in PenalizedSpline.PenaltySpectrum; BAPS/BGP baselines change by design;
fact-harness fails string-argument theory cases (MSB4177).
No uncommitted chip files remain.

## R10 follow-up 2026-10-19
f24d9d7: GP constant mean integrated out under a flat prior (GLS mean in GpFit.Mean; GpTransform maps inputs only).
fbe6d0d: BgpResult.Bandwidth / BandwidthOriginal (BandwidthSummary), alongside InputScale.
Suite at fbe6d0d: 121/121 (fact harness and plain dotnet test); BGP_BENCH passes. Report updated (Item 2a).
