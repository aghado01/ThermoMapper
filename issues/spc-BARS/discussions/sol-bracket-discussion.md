
# an open question is how to pose adjustments to the thermal curve's window when it seems like the bracketing might be too narrow. this question is intertwined with the upcoming BARS and graph compiler work because in some documentation somewhere i have isolated the mechanism for how this gets computed and connected it to the lineage work from Domany and other and their approach to defining the bracket, and my goal is to make the bracket a configurable aspect of doing a run for research purposes

<details><summary>6 previous messages</summary>

> I’ll locate that bracket discussion and reconstruct the existing mechanism before proposing anything. I’ll trace the terms around thermal windows, clipping, Domany temperature schedules, normalized \([0,1]\) coordinates, lineage persistence, and graph-derived temperature scaling, then separate what belongs in graph compilation from what belongs in the run configuration.
>
> <details><summary>Ran 6 commands</summary>
>
> - Ran `rg -n -i -g '*.md' -g '*.cs' -g '*.json' "thermal window|temperature window|window boundary|bracket|bracketing|clip(ped|ping)?|normalized thermal|\[0,\s*1\]|Domany|Blatt|lineage|Tmin|Tmax|temperature range|AutoGrid|thermal scale" issues src tests presets lean README.md CONTRIBUTING.md TODO.md`
> - Ran `rg -n -i -g '*.md' -g '*.cs' "EstimateBracket|T_ps|T_fs|pseudo.?transition|first plateau|hot end|cold end|coldOvershoot|hotOvershoot|bracket adequacy|bracket-adequacy|schedule.*bracket|bracket.*schedule" issues src tests lean`
> - Ran `rg -l -i -g '*.md' "bracket|T_fs|T_ps|coldOvershoot|hotOvershoot|AutoGridFixedSweep|clip-semantics|75240f4c" issues | Sort-Object`
> - Ran `Get-ChildItem -LiteralPath '.' -Force | Select-Object Name,Mode; rg --hidden -l -i -g '*.md' -g '!artifacts/**' -g '!.git/**' "bracket|T_fs|T_ps|coldOvershoot|hotOvershoot|AutoGridFixedSweep|clip-semantics|75240f4c" . | Sort-Object`
> - Ran ``$files = @('src\clustering\graphical\spc\runtime\scheduling\SpcScheduleHelpers.cs','src\clustering\graphical\spc\runtime\scheduling\AutoGridFixedSweep.cs','src\graphs\models\potts\BwdPottsCriticalEstimate.cs','src\clustering\graphical\spc\runtime\core\solver\PKWang.cs','src\clustering\graphical\spc\profiling\SpcProfileAnalysis.cs'); foreach ($f in $files) { "`n===== $f =====`n"; Get-Content -LiteralPath $f -Raw }``
> - Ran `git log --all --oneline -- '.discussion/issues/spc-parity/domany-parity-plan.md' '.discussion/issues/**' | Select-Object -First 40; git log --all --name-only --pretty=format: | rg -i "domany-parity-plan|spc-samplers|bracket|thermal.*window|bars" | Sort-Object -Unique | Select-Object -First 200`
>
> </details>
>
> The current implementation exposes three different bracket semantics already: a graph-statistic heuristic, the BWD \(T\_{ps}(q)\) anchor valid only under the normalized coupling, and a solver-native PKWang bracket from cumulative-energy cut quantiles. I’m looking for the design note that reconciles those with the lineage rule, because they should not be collapsed into one generic “auto” mode.
>
> <details><summary>Ran 10 commands</summary>
>
> - Ran `rg -l -i -g '*.md' -g '*.txt' -g '*.jsonl' -g '!**/node_modules/**' -g '!**/.git/**' -g '!**/artifacts/**' "bracket-adequacy|T_ps\(q\)|domany-parity-plan|Scheduling scope|coldOvershoot|hotOvershoot|SP-plateau lower bound" 'D:\aghado01'`
> - Ran ``$files = @(
'D:\aghado01\MarkBrain\ThermoMapper\issues\spc-e2e-maturity\spc-parity\temperatures-auto-brief.md',
'D:\aghado01\MarkBrain\ThermoMapper\issues\spc-e2e-maturity\spc-parity\domany-parity-plan.md',
'D:\aghado01\MarkBrain\ThermoMapper\issues\architecture-overhaul\dendrogram-integration-map.md',
'D:\aghado01\MarkBrain\ThermoMapper\issues\graph-engine-expansion\coupling-normalization.md',
'D:\aghado01\MarkBrain\ThermoMapper\issues\bars\opus-bars-asides.md',
'D:\aghado01\MarkBrain\ThermoMapper\claude\project_thermomapper_architecture.md'
); foreach ($f in $files) { "`n===== $f =====`n"; rg -n -i -C 8 "bracket|window|T_ps|T_fs|schedule|clip|lineage|BARS|temperature" $f }``
> - Ran ``Get-Content -LiteralPath 'D:\aghado01\MarkBrain\ThermoMapper\issues\graph-engine-expansion\coupling-normalization.md' -Raw; "`n===== dendrogram bracket/lineage excerpts =====`n"; rg -n -i -C 12 "bracket|T_fs|T_ps|regime border|lineage|temperature range|grid-spacing|schedule" 'D:\aghado01\MarkBrain\ThermoMapper\issues\architecture-overhaul\dendrogram-integration-map.md'``
> - Ran `rg -n "^#|^##|^###" 'D:\aghado01\MarkBrain\ThermoMapper\issues\architecture-overhaul\dendrogram-integration-map.md'; Get-Content -LiteralPath 'D:\aghado01\MarkBrain\ThermoMapper\issues\architecture-overhaul\dendrogram-integration-map.md' -TotalCount 260`
> - Ran `rg -n -i "temperature.*(estimate|scale|calibr|bracket)|TpsUpperBound|CouplingGauge|MeanDegree|TemperaturesResolved|TemperatureGridSpec|requested.*resolved" src/graphs src/user-repl src/clustering tests | Select-Object -First 300`
> - Ran ``$p='src\user-repl\commands\SpcCommand.cs'; $lines=Get-Content -LiteralPath $p; for($i=440;$i -le 510;$i++){ '{0,4}: {1}' -f ($i+1),$lines[$i] }; "`n===== Graph config/calibration clues =====`n"; rg -n -C 5 "MeanEdgeDistance|Coordination|1/K|mean-degree|Gauge|Temperature" src/graphs/pipeline/scalers/GlobalBandwidthScaler.cs src/graphs/GraphCompilerConfig.cs src/graphs/GraphBuildResult.cs src/graphs/*.cs``
> - Ran ``$files=@('src\graphs\GraphBuildResult.cs','src\graphs\GraphCompilerConfig.cs','src\user-repl\commands\RunManifest.cs','src\user-repl\commands\TemperatureGridSpec.cs'); foreach($f in $files){ if(Test-Path -LiteralPath $f){"`n===== $f =====`n";Get-Content -LiteralPath $f -Raw}}``
> - Ran `rg -n -C 12 "SignificantPeakSpans|Clipped|clip|SpanCoverage|PeakSpan" src/maths/regression/spline/bars/SplineExtrema.cs src/maths/regression/spline/bars/BarsEnsemble.cs tests | Select-Object -First 400`
> - Ran `rg -n -i -C 10 "bracket|T_fs|T_ps|lineage|window|schedule|overshoot" 'D:\aghado01\MarkBrain\ThermoMapper\claude\project_domany1999_parity.md' 'D:\aghado01\MarkBrain\ThermoMapper\claude\project_arch_overhaul_progress.md' 'D:\aghado01\MarkBrain\ThermoMapper\claude\MEMORY.md' 'D:\aghado01\MarkBrain\ThermoMapper\issues\thread-049cafab-32c9-4038-afd7-2f38d69a6f41.md'`
> - Ran `$p='src\clustering\graphical\spc\partitions\hierarchical\LineagePersistence.cs'; $lines=Get-Content -LiteralPath $p; for($i=0;$i -le 140;$i++){ '{0,4}: {1}' -f ($i+1),$lines[$i] }; for($i=250;$i -le 310;$i++){ '{0,4}: {1}' -f ($i+1),$lines[$i] }`
>
> </details>
>
> I found the design thread. The key distinction it exposes is that three things are currently called a “window”: the outer acquisition bracket, the BWD \(T*{fs}\!-\!T*{ps}\) plateau, and the resolver’s lineage window. They must become separate contracts. In particular, the parity test already records that clipping the hot side at the first \(T\_{ps}\) hides later lineage splits.
>
> <details><summary>Plan</summary>
>
> Located the original Domany parity, auto-temperature, coupling-normalization, BARS span, and lineage notes; reconstructed the existing implementations and the missing posterior boundary evidence.
>
> - [x] Locate thermal-window and bracket documentation across issues, source, and snapshots
> - [x] Reconstruct how graph statistics, Domany-style schedules, and lineage behavior currently define the bracket
> - [x] Analyze adaptive expansion criteria under the joint BARS curve posterior
> - [x] Propose a configurable research-run bracket contract and provenance model
>
> </details>

</details>

I found the relevant design thread. The existing documentation already contains most of the conceptual pieces, but they are distributed across the Domany parity plan, auto-temperature brief, coupling-normalization note, BARS span discussion, and lineage implementation.

The clean resolution is to separate three objects currently called a “window”:

1. **Acquisition bracket** \([T_{\min},T_{\max}]\): the physical domain normalized to \([0,1]\) and fitted by BARS.
2. **Thermodynamic feature regions**: \(T*{fs}\!-\!T*{ps}\), posterior peak spans, PH components, plateaus, and other analytically extracted regions.
3. **Resolver window**: the subset of temperatures a lineage or hierarchy resolver is allowed to consume.

Conflating them is the main source of ambiguity.

## What the earlier design established

The [Domany parity plan](D:/aghado01/MarkBrain/ThermoMapper/issues/spc-e2e-maturity/spc-parity/domany-parity-plan.md) records that the classical BWD/Domany method used a generously bracketed fixed sweep. The bracket was estimated before sampling; \(T*{fs}\) and \(T*{ps}\) were subsequently extracted as boundaries of the first superparamagnetic plateau.

The [auto-temperature brief](D:/aghado01/MarkBrain/ThermoMapper/issues/spc-e2e-maturity/spc-parity/temperatures-auto-brief.md) then made bracket selection explicit and provenance-bearing:

- manual temperatures remain possible;
- `auto` is requested explicitly;
- both requested and resolved schedules are recorded;
- the bracket depends on the coupling construction.

The [coupling-normalization note](D:/aghado01/MarkBrain/ThermoMapper/issues/graph-engine-expansion/coupling-normalization.md) contains the important graph-engine connection:

> The coupling construction declares its temperature calibration.

For the BWD \(1/\widehat K\)-normalized coupling, the \(q\)-only \(T\_{ps}(q)\) estimate is valid. For an unnormalized coupling, the scale remains graph-dependent. For PKWang, temperature is a single-linkage cut height and the bracket comes from the cumulative-energy ladder instead.

Those are different calibration families, not alternate implementations of one formula.

## Bracket adjustment should be an outer inference loop

The BARS reconnaissance loop you described sits inside a larger bracket-validation loop:

\[
\text{resolve initial bracket}
\rightarrow
\text{fit joint curve posterior}
\rightarrow
\text{evaluate boundary adequacy}
\rightarrow
\begin{cases}
\text{expand and refit},\\
\text{accept and extract features}.
\end{cases}
\]

Within a bracket epoch \(B*r=[T*{\min}^{(r)},T\_{\max}^{(r)}]\), you perform the endpoint-inclusive complementary schedules and recursively fit the union of measurements.

If the bracket is rejected, expand directionally—for example on a log-temperature axis,

\[
T*{\min}^{(r+1)}=T*{\min}^{(r)}/\rho*L,\qquad
T*{\max}^{(r+1)}=\rho*R T*{\max}^{(r)}.
\]

All old observations remain valid. They are remapped into the new normalized coordinate and the joint curve is refitted using the enlarged union. Nothing should be extrapolated from the old spline into the newly exposed physical interval.

This gives two nested stopping conditions:

- the curve is sufficiently resolved within the current bracket;
- the bracket itself has sufficient posterior evidence of containing the relevant thermal behavior.

Bracket adequacy is therefore part of “confidence reached,” not an unrelated precondition.

## Boundary adequacy should be posterior evidence

The current BARS engine already computes directional clipping for every draw in `PeakSpan.LeftClipped` and `RightClipped` in [SplineExtrema.cs](./src/maths/regression/spline/bars/SplineExtrema.cs:12). But those flags are discarded when `BarsEnsemble` pools the draws. `BarsResult` retains span coverage but not the posterior probability of directional clipping.

That is the nearest concrete engine gap.

A useful posterior boundary report would include:

- probability of any significant left-clipped span;
- probability of any significant right-clipped span;
- expected number of clipped spans on each side;
- probability of a significant boundary maximum;
- probability that the endpoint derivative points outward:
  - \(f'(0)<0\) on the cold side;
  - \(f'(1)>0\) on the hot side;
- posterior endpoint height relative to the inferred stable baseline;
- probability that a persistent superlevel component reaches either boundary.

These should be probabilities or expected counts—not a single clip boolean.

`SpanCoverage[0]` and `SpanCoverage[^1]` are insufficient because they do not distinguish an ordinary wide span from a genuinely truncated one.

There is also an identifiability limit worth stating plainly:

> An in-window curve cannot conclusively distinguish a genuine boundary transition from a transition truncated by a narrow physical bracket.

The only decisive experiment is to sample beyond the boundary. Consequently, posterior boundary evidence should trigger a controlled exploratory expansion, not a claim that the bracket is definitely wrong.

## Domany landmarks should not become outer boundaries

\(T*{fs}\) and \(T*{ps}\) are curve features, not acquisition-bracket definitions.

The lineage parity test already encodes this lesson in [BwdParityTests.cs](./tests/VizCore.Tests/BwdParityTests.cs:381): it excludes the cold ferromagnetic regime below \(T*{fs}\), but deliberately does not clip the hot side at the first \(T*{ps}\), because doing so hides the hotter versicolor/virginica split.

Therefore, useful research resolver modes might include:

- full accepted acquisition bracket;
- first Domany plateau \([T_{fs},T_{ps}]\);
- post-ferromagnetic discovery range \([T_{fs},T_{\max}]\);
- union of posterior significant spans;
- explicit physical interval.

These belong to the resolver configuration, separately from bracket acquisition.

## Lineages can provide a second adequacy signal

Once targeted clustering has produced a temperature stack, lineages themselves can reveal boundary censoring:

- a lineage present at \(T\_{\min}\) has a left-censored birth;
- a lineage surviving to \(T\_{\max}\) has a right-censored death;
- several distinct clusters already present at the cold boundary may indicate that the cold coherent regime was missed;
- substantial coherent lineages surviving at the hot boundary may indicate that fragmentation has not been observed.

These are structural follow-up signals. They can recommend another bracket epoch or a subsequent run, but they should not replace the analytic curve’s initial boundary diagnostics.

There is a related issue in the current lineage score. [LineagePersistence.cs](./src/clustering/graphical/spc/partitions/hierarchical/LineagePersistence.cs:24) uses

\[
\text{LevelCount}\times\text{Decidedness}.
\]

The comment calls step count grid-agnostic, but it is actually sampling-density-dependent. Once BARS creates complementary refinements and targeted dense regions, a lineage can score higher merely because more temperatures were sampled in its neighborhood.

A BARS-compatible lineage measure should integrate survival against a declared thermal measure:

\[
S(L)=
\overline{\delta}_L
\int_{I_L} w(x)\,dx,
\]

where \(x\in[0,1]\) is the bracket’s declared normalized coordinate. Numerically, this can use interval widths or quadrature weights. It then remains invariant when the same curve is sampled more densely, while still allowing research comparisons between linear-T, log-T, or thermodynamically warped measures.

## Proposed configurable contract

The bracket configuration should factor three independent choices:

```text
ThermalWindow
├── Initial source
│   ├── Manual
│   ├── Construction calibration
│   └── Solver-native calibration
├── Coordinate
│   ├── Linear temperature
│   └── Log temperature
├── Boundary policy
│   ├── Fixed
│   ├── Diagnose only
│   └── Posterior expansion
└── Expansion limits
    ├── Directional factor
    ├── Maximum epochs
    ├── Absolute temperature bounds
    └── Evidence threshold
```

Named research presets could then express:

- **Domany fixed**: construction-calibrated, generous fixed bracket;
- **manual fixed**: exact investigator-specified range;
- **BARS diagnostic**: report bracket inadequacy without mutating the run;
- **BARS expandable**: directional posterior-triggered expansion.

This preserves the classical method as an experimental baseline rather than silently replacing it.

## Graph-compiler boundary

The graph compiler should not own adaptive BARS policy. It should declare the coupling’s scale semantics and calibration evidence; SPC should decide how that calibration becomes a run bracket.

The current CLI branches directly on `BandwidthStrategy` in [SpcCommand.cs](./src/user-repl/commands/SpcCommand.cs:477). That is brittle for precompiled graphs and future projections because it infers calibration from the requested option rather than reading it from the emitted graph artifact.

The graph-engine renovation could instead emit a serializable calibration descriptor containing:

- coupling gauge and normalization;
- calibration family/key;
- resolved scale statistics;
- whether the BWD \(q\)-only anchor is admissible;
- relevant assumptions and provenance.

The SPC bracket resolver then combines that descriptor with \(q\), the requested bracket policy, and any manual padding. Thus:

> The projection declares temperature-scale semantics; the run declares bracket policy.

Finally, provenance should record more than the final range:

```text
requested policy
initial calibration and resolved bracket
coordinate transform
each expansion epoch
boundary evidence causing expansion
new temperatures measured
final accepted bracket
termination reason
resolver window selected afterward
```

That makes bracket adjustment a legitimate experimental axis rather than hidden scheduler behavior.

The most immediate enabling addition is therefore not the expansion scheduler itself. It is retaining and pooling the directional boundary evidence that `PeakSpan` already computes. Once that posterior artifact exists, fixed, diagnostic, and adaptive bracket policies can all consume the same evidence without contaminating the general BARS engine with SPC-specific decisions.
