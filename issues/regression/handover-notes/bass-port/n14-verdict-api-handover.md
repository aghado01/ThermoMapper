# N14: the shared verdict landed (e8547d3, 2026-10-29)

Update: N14 complete. BASS `e88b33e` (log posterior, groups predictive/level/structural, `BassSampling.StoppingScope`, `BassDiagnostics.Verdict` + `ModeSeparation`), clients `1c26dcd` (BARS groups predictive/structural; BGP/BAPS/HPS results carry `Verdict`), report `6fe653a` in the brief. Every BASS config fingerprint moved (the JSON now holds `StoppingScope`); pre-N14 checkpoints will not resume.

For stage 7's C3. All in `Maths.Samplers.Ensemble` (`src/maths/samplers/ensemble/ConvergenceVerdict.cs`, `ChainEnsemble.cs`).

- `ChainEnsemble.ReferenceRHatTarget = 1.01`, beside `MinEssPerChain`.
- `EnsembleRun.Verdict` : `ConvergenceVerdict(ConvergenceOutcome Outcome, double Target, double MinEssPerChain, IReadOnlyList<GroupVerdict> Groups)`, with `Group(string name)` lookup (throws `KeyNotFoundException`).
- `GroupVerdict(string Name, ConvergenceOutcome Outcome, bool GatesStop, IReadOnlyList<FunctionalReading> Failing)`.
- `FunctionalReading(int Index, string Name, RHatDiagnostic RHat, EssDiagnostic? Ess)`: `Ess` is null for an unmonitored functional (classic R-hat only).
- `ConvergenceOutcome { Converged, NotConverged, Undetermined }`.
- `IEnsembleModel<TDraw>.FunctionalNames` (default "f0", "f1", ...) and `.Groups` (`IReadOnlyList<FunctionalGroup>`, default one group "all" = `FunctionalGroup.AllName`). `FunctionalGroup(string Name, IReadOnlyList<int> Functionals)`; each functional in exactly one group, checked at `Run` entry.
- `ChainEnsemble.Run(..., ParallelOptions? parallel = null, IReadOnlyCollection<string>? stoppingScope = null)`: the scope is the LAST parameter (so no positional caller breaks); pass it by name. Null = every group; empty or unknown names are refused (`ArgumentException`, param "stoppingScope").
- Rules: frozen and unavailable readings FAIL (NotConverged), as B01 and D10 say; only one chain, or every gating functional constant, is Undetermined. (The brief's API bullet said "constant or frozen" for Undetermined; the rules section and D10 say frozen fails, which is what is built.)
- A run with `rHatTarget = 0` still gets a verdict, read against 1.01.
