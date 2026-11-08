# Targeted test runs — report (2026-11-07)

Routine runs of `tests/Maths.Regression.Tests` were the whole assembly at every commit: 2,160 tests, 10 m 32 s on the loaded machine at the last full run, plus the build. The cause was the common gate in the active briefs ("the eight pinned hashes and every existing fingerprint hold"), with the fingerprint tests untagged, so only the whole assembly could check it. This work makes the routine run targeted: a profile, a pinned tier, a slow tier, selection by change, and the rule in `CONTRIBUTING.md`.

## What landed

| Commit | What |
|---|---|
| `39c1008` | The pinned tier: `[Trait("Tier", "Pinned")]` on 9 methods or classes, 31 cases |
| `42de39d` | `scripts/test-select.ps1` and its name graph `scripts/test-select.cs` |
| `82e91eb` | The slow tier, `tests/Maths.Regression.Tests/slow-tests.txt`, from the profile |
| `1d8ffa1` | The rule in `CONTRIBUTING.md`'s testing standards, `AGENTS.md`'s script list, and the common-gate and working-rule lines of the GBASS, KHAOS and conformal briefs |
| The commit carrying this report | `CONTRIBUTING.md`'s section wrapped to the file's width, with `-SourceRoot` for snapshot builds |

Nothing pushed; each commit by pathspec. The conformal brief carried the conformal session's uncommitted line 134, so it went in as `HEAD` plus my two lines, and their line stays in the working tree.

## The profile

The whole assembly at `d8b53ea`, from a `git archive` snapshot built into this work's own artifacts root, run under `dotnet test` with a trx logger at below-normal priority, 6 xUnit threads, build servers off, with two other chips' test hosts running.

- 2,172 of 2,172 passed in 7 m 17 s of wall time (build 27 s). 2,572 s summed over the cases.
- The time is concentrated:

  | Cases at or above | Cases | Share of the time |
  |---|---|---|
  | 10 s | 70 (3%) | 69% |
  | 5 s | 101 | 78% |
  | 2 s | 185 | 88% |
  | 1 s | 295 | 94% |

  The top case is 98 s (`BassCalibrationCutSmcTests.AgreesWithTheNestedReference_OnBassFits(basis: True)`). By class, `BassFlatHistogramTests` holds 435 s, `BassPriorReferenceTests` 330 s, `BassCalibrationCutSmcTests` 260 s, `BassMovingResponseTests` 212 s and `BassTemperingTests` 208 s: five classes, 56%.
- xUnit runs a class's tests one after another, so the assembly's wall time cannot fall below its longest class: `BassFlatHistogramTests`' 435 s (7 m 15 s) sets the 7 m 17 s. More threads would not shorten it.
- **Decision taken: a slow tier.** The methods whose cases sum to 10 s or more, 46 of 1,174 methods, hold 2,068 s (80%). The rest is 504 s, about 85 s of thread time per worker at 6 threads, and its longest class is 32 s. A 5 s cut would add 17 methods for 5 points more, so 10 s stands.

## The pinned tier

Each candidate was read, starting from `rg` over `fingerprint|pinned|bit for bit|hash` and the hex literals; membership is an assertion against an expectation stored from an earlier commit.

| Test | What it pins | Cases | Profile |
|---|---|---|---|
| `BassBitIdentityTests` (class) | The eight pinned hashes | 8 | 3.0 s |
| `BassTemperingSizingTests.FixedLevels_RunAsS8Did` | S8's seven fixed-ladder hashes | 7 | 0.3 s |
| `BassFlatHistogramTests.Pilot_PoolsItsLaddersRungs_…` | The pilot's accepted starts | 1 | 15.7 s |
| `BassFlatHistogramTests.WarmStart_OnFortyEightFixtureAPilots_…` | The 48 pilots' accepted starts | 1 | 45.9 s |
| `BassFlatHistogramTests.Fit_FromAPilotWithAGap_…` | A gap pilot's accepted start | 2 | 2.2 s |
| `BassFlatHistogramTests.TodaysFits_FromTheGapPilots_ArePinned` | Two fits from the gap pilots | 1 | 0.9 s |
| `BassFlatHistogramQualificationTests.Friedman_ArmOne_WarmStart_IsPinned` | Seed 8017's start (behind `BASS_QUALIFICATION`) | 1 | 0 s |
| `BassMovePassTests.Friedman_AtN500_IsS10s` | S10's data, with the default config's fingerprint | 1 | 0 s |
| `BassColumnBitsTests` (class) | The column kernels against verbatim copies of their old forms | 9 | 0.7 s |

31 cases, 69 s summed; 65 s of it runs serially in `BassFlatHistogramTests`, so the tier alone takes about a minute of wall time. Checked in a run: a `Tier=Pinned` filter took exactly these 31 cases.

Read and left out:
- the bit-for-bit comparisons between two runs of today's code: resume against the uninterrupted run, an offset against shifted data, the memo against the direct form (`BassMarginalTermsTests`), one draw by two paths. They are properties the selection covers;
- `BassPrimitiveCostTests`' old-against-new bits, a qualification cost harness whose assertions `BassColumnBitsTests` already holds;
- `EbicTests` and `SpanCrossingTests`, "pinned" against closed forms within a tolerance.

**Finding: one fingerprint pin.** No test pins a config's fingerprint except through `Friedman_AtN500_IsS10s`, which hashes the default config's fingerprint on S10's data. The gate's "every existing fingerprint" for a config with an option set (tempering, SMC, the flat histogram, ridge, a noise, Legendre terms) rests on each option's own test that its absent key stays out of the JSON (for example `KhaosChainTests.Config_SerializesTheBasisAndTheRawRidge_AndAHingeConfigAsItWas`). See open items.

## Selection by change

`src/repo-audit` gave no type-reference graph to reuse: its semantic compile yields diagnostics and assembly-level used references, and its syntax walker's per-file names sit in internal records. The graph is its own, in `scripts/test-select.cs`:

- It covers the test project's compile closure, read from the csproj files, which keeps it inside the assembly's sources.
- A changed file reaches each file that names one of its top-level types, delegates or extension methods, to the fixed point. A test class is selected when its file is reached.
- Comments and strings are blanked; interpolation holes are kept.
- A name read as a member or a call is not a type's, unless its dotted chain starts at one of the closure's namespaces (`Maths.Distributions.Normal.SampleStandard`).

The first version read every identifier and selected every class for any change: the private test-nested records `Fit` and `Run` matched the methods `Bass.Fit(…)` and `.Run(…)`.

The script finds the same 204 test classes the profile ran.

**The check on recent commits.** Each commit's selection was computed from its `src/` changes alone, in a snapshot of its own tree; the test classes the commit touched must be among those selected.

| Commit | What | Classes selected | Touched classes selected | Thread-seconds (profile), without the slow tier |
|---|---|---|---|---|
| `4e62cf6` | Efficiency lever 1, which moved two pinned hashes | 82 of 177 | 3 of 3, `BassBitIdentityTests` among them | 236 |
| `37137f2` | GBASS slice 5 | 104 of 202 | 3 of 4 | 277 |
| `baf69f9` | KHAOS slice 1, raw ridge | 106 of 204 | 2 of 2 | 278 |
| `b408bf4` | KHAOS slice 1, Legendre part | 104 of 200 | 2 of 2 | 270 |
| `6dfa192` | Conformal split form | 109 of 198 | 3 of 3 | 281 |
| `fa05048` | Two distributions primitives, not yet used by the assembly | 0 of 198 | none touched | 0 |

The one class not selected, `GbpMixingTests` in `37137f2`, gained two cases that exercise an unchanged primitive at a new point: `GbpMixing` lives in `MixingLaw.cs`, last changed in `4662639`. No `src/` change in that commit could move it. A chip editing that test file selects it by the file's own change.

Breadth by where the change is (HEAD, 211 classes): a conformal file 10, `Cholesky` 15, a BASS core file or a KHAOS primitive about 112, `Normal` 208, the RNG 210. The graph follows references, not execution, so a hub (`BassModel`, `BassFit`) passes a leaf's change to everything that names the hub.

## Timings, before and after

Measured on a snapshot of `957629b` with the pinned traits, built once (35 s), each run `--no-build`, on the shared machine at below-normal priority, 6 threads.

| Run | Tests | Test time | End to end |
|---|---|---|---|
| Before: the whole assembly (the profile) | 2,172 | 7 m 17 s | 7 m 53 s with the build |
| Before: the last reported full run, default threads, loaded | 2,160 | 10 m 32 s | plus the build |
| Edit, one BASS core file (107 classes, no slow tier) | 1,476 | 1 m 47 s | 114 s |
| Commit, one conformal file (8 classes and the pinned tier) | 88 | 1 m 6 s | 86 s |
| Commit, one BASS core file (107 classes and the pinned tier) | 1,487 | 4 m 11 s | 258 s, loaded (see below) |
| Edit, every class without the slow tier (a build file changed) | 2,042 | 2 m 57 s | 185 s |

The selection itself takes about 7 s (the graph's compile and a read of 475 files).

These are single readings on a machine whose load moved under them: the Commit run on the BASS core file ran while other chips' test hosts held the cores, and its common tests summed 1,156 s against 632 s in the Edit run two minutes before (`BassFactorMoveTests` 79 s against 6.5 s in the profile). The Edit run's own sum is about twice what the profile predicts for its classes (about 280 s). The profile's per-class shares are the steadier guide: a selection without the slow tier is 9 to 11% of the assembly's thread time for a BASS core change, and under 1% for a leaf (15 s, 0.6%, for a conformal file).

## Open items, each with the decision it would change

1. **Pin the fingerprints of option-on configs.** A theory over the bit-identity fixtures' configs and one config per option (tempering, SMC, flat histogram, ridge, raw ridge with Legendre terms, a noise, local steps, factor moves), each pinning `Bass.Fingerprint`, would run in milliseconds and join the pinned tier. *Decision:* whether "every existing fingerprint holds" is checked by a test at each commit. If it is added, the gate is mechanical. If not, that half of the gate rests on each option's JSON test, which catches an option's own key appearing but not a changed default elsewhere in the config.
2. **`BassFlatHistogramTests`' serial 435 s sets the slice-end wall time.** *Decision:* whether to split its enumeration and end-to-end law checks into a class of their own, so xUnit runs them beside the rest. If split, the whole-assembly wall time falls toward the next class (`BassPriorReferenceTests`, 330 s), about 1.5 min less. If left, a slice's end costs about 7 min, as now. The flat-histogram chip owns that file.
3. **The pinned tier's minute is mostly one test.** `WarmStart_OnFortyEightFixtureAPilots_…` (46 s) and `Pilot_Pools…` (16 s) share the 48 pilots' fixture. *Decision:* whether a commit pays about 65 s for them. Keep them per commit while flat-histogram code moves; if commits rarely touch that code, moving the two to the slice end would leave a pinned tier of about 5 s.
4. **Narrower selection for hub changes.** A per-class coverage profile (which `src/` files each class runs) would select by execution rather than reference. *Decision:* worth building only if per-edit runs on core files (about 2 min now) remain the bottleneck. Its cost is a coverage run per refresh, and a stale profile misses new paths.
5. **The slow list drifts.** It was measured under load at `d8b53ea`. A new profile (`-WriteSlowList <trx>`) rewrites it; nothing waits on it until a slow method grows or a fast one becomes slow.
6. **Other assemblies.** `Maths.Distributions.Tests` (347) and `Maths.Samplers.Tests` (285) have no slow list. The script takes them with `-Project`, and their selection works as is. Nothing waits on profiling them while they run in a few minutes.
