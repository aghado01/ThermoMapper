# fact-harness: theory cases (MSB4177) — fixed, committed 1f96cb4

2026-10-19. Files: src/test-harness/{FactRunner,Program,PathUtil}.cs, scripts/fact-harness.ps1. Committed on main as 1f96cb4 (user approved).

## Root cause
The harness took each `--list-tests` line (a display name) and ran `--filter FullyQualifiedName=<display name>`.
- A theory display name carries its arguments, but xUnit's FullyQualifiedName does not, so a theory case could never match.
- Parentheses are grouping syntax in the VSTest filter. Enum and int cases hit "Incorrect format", matched 0 tests, and exited 0, so they **passed without running**. That was 27 of the 117 at ed60885 (every theory case).
- String cases: the `"` toggles quoting in the MSBuild `-property:VSTestTestCaseFilter=...` that `dotnet test` generates, which gives MSB4177 (or MSB1006 for `;`).
- Two-layer escaping (VSTest grammar `\` plus MSBuild `%XX`) fixes everything except a string argument that itself contains `"`: VSTest's MSBuild task re-quotes backslash-quote pairs wrongly before calling vstest.console.

## Fix
- Each case gets `NNN.runsettings` in the run dir, containing `RunConfiguration/TestCaseFilter = DisplayName=<VSTest-escaped name>` and `TreatNoTestsAsError=true`, passed with `--settings`. The XML route avoids both command-line layers; the ordinal-only filename keeps `%`/`;` out of the MSBuild property.
- A filter that matches nothing now FAILS (exit 1).
- `SuiteNameFromFact` drops the argument list before splitting on '.'. Before this, `(x: 0.5)` and `(bad: 9.99E+199)` produced fake suites. Example: in the ed60885 run, `NonFiniteSupport_IsRejected(bad: 9` and `KernelAffinity_ScalesAsIntrinsicDimension(d: 1, expectedSlope: 0` were filed as their own suites.
- summary.json and suite-manifest format are unchanged. The only new files are the per-case `NNN.runsettings` in the run dir.

## Verification
- Throwaway project in this session's scratchpad (not in the repo), 19 cases: plain strings, `; , & | = ! ~ \ % $() @()`, an embedded `"`, parentheses, dots, an empty string, a mixed double/string/enum case, a fact, a custom-DisplayName fact, and one deliberately failing case.
  - Before the fix: 14 FAIL (MSB4177/MSB1006) and 3 PASS that ran nothing.
  - After: 18 PASS and 1 FAIL (the deliberate one). Every log shows `Total: 1`.
- Maths.Regression.Tests at fbe6d0d (main moved during this work): 121/121 passed, all 29 theory cases included, every log `Total: 1`, no "No test matches" (artifacts/test-runs/TestHarness.Runner/20261019_202408907).

## Follow-ups (not done)
- Any earlier harness "pass" for a theory case is void. The enum workaround in ed60885 can return to string InlineData if wanted (main-track call; those test files are not in this change).
- `--common-arg --settings/--filter` would override the per-case selection; a guard is possible.
- VSTest `=` is case-insensitive, so two cases that differ only in letter case would each run both.
- A project with its own RunSettingsFilePath would have it replaced by the per-case file. None exist in the repo today.
