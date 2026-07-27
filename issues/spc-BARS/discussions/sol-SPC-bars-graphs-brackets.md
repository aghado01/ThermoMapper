# Sheaf literature is already in queue for synthesis review, the magnetic laplacian was the first step in that direction, and i have more pressing priorities before i develop on that such as finishing my zigzag work, SPRED, the graph engine renovation and the next-gen SPC with BARS thermal curve estimation

That sequencing makes sense. I’ll treat the sheaf/groupoid material as a queued synthesis horizon, with the magnetic Laplacian understood as an intentional first foothold—not as a recommendation for near-term development.

The practical priority stack I infer is:

1. Finish the zigzag foundations and establish trustworthy nonmonotone persistence semantics.
2. Advance SPRED on its own architectural track.
3. Complete the graph-engine renovation and preserve its neutral artifact boundary.
4. Build next-generation SPC around posterior thermal-curve estimation with BARS.
5. Only afterward revisit higher-order structural equivalence, sheaves, and HOPE-inspired consolidation.

Within that sequence, two findings from my review are immediately relevant:

- `NerveDiff.MatchNodes` appears to assume a point belongs to only one Mapper node. Overlapping Mapper covers violate that assumption, so this may matter directly to the zigzag work.
- BARS already exposes posterior peak intensity, span coverage, peak modes, and convergence diagnostics. Next-gen SPC therefore needs an integration and statistical-design layer—not a new peak-intensity primitive.

For SPC+BARS, the important near-term design questions seem to be:

- which `SweepProfile` channels receive posterior fits;
- how Monte Carlo uncertainty and autocorrelation become observation weights;
- whether curves such as co-membership require monotonic constraints;
- whether node/edge curves are modeled individually, through shared knots, or through latent response modes;
- what typed posterior artifact crosses from the BARS engine into SPC hierarchy and partition strategies.

I’ll keep the broader HOPE synthesis parked as a design lens: useful for deciding what response artifacts should be preserved now, but not a justification for pulling future sheaf or meta-clustering work ahead of the engine roadmap.

# review the `lean` directory protolemmas, enthymemes and lemmas for context on the BARS x SPC theoretical work

<details><summary>3 previous messages</summary>

> I’ll read the Lean project as mathematical context, not merely scan theorem names. I’ll start from the project structure and supporting notes, then trace the protolemmas, enthymemes, and proved lemmas that constrain the intended BARS×SPC theory.
>
> The formal queue reveals a deliberate chain rather than a single “BARS finds the SPC peak” claim: exact draw-wise spline reductions, thermodynamic monotonicity, a \((\beta,\theta)\) bifiltration, and finally a confidence pushforward from curve uncertainty to barcode uncertainty. I’m checking each link for whether it is proved, merely stated, or still only described informally.
>
> The tier labels need a careful correction: several markdown briefs describe P1–P5/T1 as “proved,” but that code has not been graduated into the compiled Lean files. I’m running the project’s own meta-CI now so the final review distinguishes mathematical proof sketches from declarations that Lean has actually accepted.

</details>

The validation is still installing the pinned Lean 4.30 toolchain in this shell context; meanwhile, the source audit already establishes that the no-apologies `Lemmas` library is empty. The only completed proof body is the PKWang threshold reduction, currently retained in the enthymeme tier; the BARS declarations all still contain `sorry`.

> no need to litigate the Lean machinery. i was pointing you to read the concepts related to BARS and SPC to understand what my original work is aiming at in the upcoming SPC updates

Understood—I over-indexed on the formalization status. The conceptual material makes the upcoming SPC direction much clearer.

The intended BARS×SPC work is not simply “replace the susceptibility peak detector with a smoother.” It is a redesign of the thermal sweep as a posterior inference process.

BARS has four distinct roles:

1. **Continuous thermal reconstruction.** Sparse Swendsen–Wang observations become posterior draws of continuous response curves rather than disconnected values at a temperature grid.

2. **Draw-wise landmark extraction.** Peaks, peak counts, prominence spans, and widths are reduced within each posterior draw before aggregation. This preserves uncertainty that is destroyed by detecting features only on the posterior mean.

3. **Adaptive experimental design.** The initial sweep discovers candidate transition regions; posterior peak intensity, span coverage, and uncertainty then determine where additional SW sampling is worthwhile.

4. **Structural regularization.** For quantities whose ideal estimands are theorem-backed as monotone—particularly FK connection or normalized spin-agreement trajectories—constrained curve estimation repairs finite-Monte-Carlo flicker and yields posterior edge-birth temperatures.

That last role connects BARS to the \((\beta,\theta)\) construction:

\[
K*{\beta,\theta}
=
\operatorname{Flag}
\left\{(i,j):G*{ij}(\beta)\geq\theta\right\}.
\]

Cooling along \(\beta\) strengthens ideal FK connectivity, while lowering \(\theta\) admits more edges. The result is inherently two-parameter. A path through this surface produces a conventional persistence module; varying the path reveals which conclusions depend on the chosen thermodynamic slice.

The documents also make an important distinction between three kinds of uncertainty that next-generation SPC should preserve:

- **Curve multimodality:** one thermal response contains several transitions.
- **Posterior multimodality:** draws disagree about where a transition occurs.
- **Structural nonmonotonicity:** cluster systems merge and split across temperature.

Those require different instruments. Peak-count posteriors address the first, tempering and location posteriors address the second, and zigzag persistence addresses the third. They should not be collapsed into one generic notion of “multiple peaks.”

A natural next-generation SPC flow therefore looks like:

\[
\text{sparse SW sweep}
\rightarrow
\text{posterior thermal fields}
\rightarrow
\text{adaptive resampling}
\rightarrow
\text{posterior landmarks and edge births}
\rightarrow
\begin{cases}
\text{SPC hierarchy/lineages}\\
(\beta,\theta)\text{ topology}\\
\text{structural zigzag}
\end{cases}
\]

The existing `BarsResult.PeakIntensity` and `SpanCoverage` fit this architecture directly. Their significance is not merely visualization:

- peak intensity is the matching-free first moment of the random set of transitions;
- span coverage describes posterior support for structurally broad transition regions;
- the narrowest credible structural span can serve as a post-hoc audit of whether the first temperature grid had enough resolution;
- phase-two sampling can target a combination of intensity, uncertainty, and unresolved span width.

This explains why BARS and zigzag are complementary. BARS should monotonize only quantities whose ideal target is actually monotone. It should not erase genuine merge/split behavior in partitions or force susceptibility itself to be monotone. The denoised monotone edge field supports the bifiltration; the nonmonotone evolution of clusters remains material for zigzag analysis.

A few conceptual guards in the notes are worth retaining:

- Structural peak width and epistemic uncertainty in peak location are different objects.
- Exact spline extrema remove numerical search error, but do not establish that the fitted posterior is calibrated.
- The existence of examples where a mean curve loses peaks proves that mean-first reduction is unsafe; it does not prove averaging always undercounts.
- Phase-one coverage is logically prior to exact feature extraction. Exact roots of an under-resolved fitted curve are exact answers about the wrong curve.
- Sup-norm-to-bottleneck confidence propagation requires the relevant filtration stability map. Exact peak and span computation helps landmark uncertainty, but does not by itself certify an entire persistence diagram.
- Edgewise independent BARS fits may be neither computationally practical nor statistically optimal. Shared-knot, hierarchical, or latent-field estimation is likely the eventual form if the full edge × temperature surface is modeled.

So I now understand the original contribution as a **Bayesian thermodynamic resolution engine**: SPC generates stochastic multiscale evidence, BARS turns it into posterior thermal geometry and allocates further computation, and PH/zigzag turn the resulting relational evolution into structural invariants. The susceptibility peak becomes only one observable in that system—not the defining output.

# the workflow also inverts the feature extraction and clustering. I iteratively estimate the joint posterior whole curve via uniform temperature schedules starting with the initial pass that includes measurement at the normalized thermal window boundary points [0,1], do successive complementary schedule sweeps in between previous points, and recursively fit the union of all points measured up to the current iteration of this stage. then when confidence is reached, i analytically extract features using things i've already built in the bars engine or things that i can develop further there or in SPCs feature extractor utilities. and its after this that targeted expensive SW sweeps can be made in the identified regions of interest of the given thermal curve and perform clustering

Yes—that is a more consequential inversion than I had captured.

Conventional SPC effectively does:

\[
\text{dense thermal sweep}
\rightarrow
\text{cluster at every temperature}
\rightarrow
\text{inspect aggregate curves}
\rightarrow
\text{select a partition}.
\]

Your next-generation workflow instead does:

\[
\text{global curve inference}
\rightarrow
\text{analytic feature extraction}
\rightarrow
\text{targeted high-cost SW}
\rightarrow
\text{clustering}.
\]

More precisely:

### 1. Coverage-preserving thermal reconnaissance

Begin with the normalized thermal window boundaries:

\[
S_0=\{0,1\}.
\]

Each subsequent schedule complements the existing design by sampling between previously measured temperatures—for example, midpoint refinement—while retaining all observations:

\[
U*r=\bigcup*{j=0}^{r}S_j.
\]

At iteration \(r\), BARS is recursively refit to the entire accumulated dataset on \(U_r\), not merely to the newest sweep. This produces an increasingly resolved **joint posterior over the complete thermal curve**. The nested design preserves global coverage, prevents early peak guesses from creating blind regions, and makes the refinement history reproducible.

The key point is that this stage is not yet trying to cluster at candidate temperatures. It is learning the global thermodynamic landscape sufficiently well to determine where clustering effort has scientific value.

### 2. Stop on functional confidence

The stopping decision should concern the posterior curve as a whole, not merely stabilization of one estimated \(T_c\). Possible evidence includes:

- convergence of the BARS chains and functionals;
- stabilization of posterior peak count and peak intensity;
- stabilization of significant spans;
- sufficiently narrow uncertainty over the function or relevant derivatives;
- absence of unresolved spans near the current schedule spacing;
- adequate behavior at both normalized-window boundaries;
- agreement of features across successive complementary refinements.

This is where the span machinery can act as a self-audit: a posterior-supported feature whose structural width is near or below the current sampling gap signals that the reconnaissance schedule is still inadequate.

### 3. Analytic feature extraction

Only once the whole-curve posterior is trusted do you analytically extract features from posterior draws:

- global and local extrema;
- peak counts and locations;
- prominence;
- significant-peak spans;
- crossings and widths;
- plateaus, cliffs, shoulders, or derivative landmarks;
- posterior peak intensity and span coverage;
- potentially relationships among several thermal channels.

This preserves the crucial reduction order:

\[
\text{feature posterior}
=
\text{aggregate}\bigl(\text{analytic feature per curve draw}\bigr),
\]

rather than extracting features from the posterior mean curve.

The feature extractor therefore converts a curve posterior into a posterior collection of thermal regions of interest. Those regions may overlap, vary in number across draws, and carry different kinds of evidence.

### 4. Targeted production sweeps and clustering

The expensive SW work then occurs only in the identified regions:

- increase equilibration, draws, or replicas;
- accumulate richer edge and node currencies;
- resolve competing transitions separately;
- construct partitions, dense temperature stacks, lineages, and dendrograms;
- compare clustering behavior across posterior-supported thermal regions.

Thus feature extraction is not a readout from already-computed clusters. It becomes the mechanism that decides **where clustering should be performed and at what resolution**.

This makes BARS an upstream experimental-design engine, not a downstream smoothing accessory.

The separation also suggests two distinct SW configurations:

- a reconnaissance configuration optimized for economical, globally representative thermal measurements;
- a production configuration optimized for reliable affinities, alignments, co-membership, landscapes, and final clustering within selected regions.

The resulting architecture is closer to sequential Bayesian experimental design:

\[
\boxed{
\text{uniform complementary exploration}
\rightarrow
\text{global posterior confidence}
\rightarrow
\text{analytic ROI discovery}
\rightarrow
\text{conditional computational exploitation}
}
\]

That is the real inversion: ThermoMapper first learns the thermodynamic coordinate system in which clustering is meaningful, then spends the clustering budget within that learned geometry.

# right but the analytic curve itself also provides a novel substrate and connections to PH and opens new ways to use SPC to examine thermal landscapes, precisely because of the theoretical guarantees afforded by estimating the joint posterior of the analytic curve rather than piece wise disjoint sampling heuristics of the old paradigm

Exactly. I was still treating the analytic curve primarily as an upstream scheduler. More accurately, it is a new mathematical data object at the center of SPC.

The old paradigm produces a collection of separate estimates:

\[
\{(T*j,\widehat y_j,\operatorname{SE}\_j)\}*{j=1}^m.
\]

Even if those points are interpolated afterward, that does not automatically define a coherent probability law over functions. Cross-temperature dependence, feature ordering, derivative uncertainty, and the probability of entire curve shapes remain unspecified.

Your workflow instead estimates

\[
\Pi(df\mid\mathcal D),
\]

a joint posterior on the complete thermal response \(f:[0,1]\to\mathbb R\). Each posterior draw is a coherent, globally defined thermal landscape. That makes pathwise questions well-posed:

- How many transitions does this draw contain?
- Which peaks coexist?
- In what order do they occur?
- Which shoulders resolve into separate modes?
- How wide is each transition region?
- Does a level-set component reach the thermal-window boundary?
- What topology does the complete thermal landscape induce?

These cannot be recovered correctly from independent pointwise intervals or from a pooled mean curve.

## The curve as a PH substrate

For a scalar thermal curve, the natural first construction is superlevel-set persistence. As the level \(\ell\) descends,

\[
X\_\ell(f)=\{T\in[0,1]:f(T)\geq\ell\},
\]

components are born at local maxima and merge at intervening minima. Consequently:

- peaks become births;
- valleys determine mergers;
- bar persistence becomes a global prominence measure;
- level-set components supply thermal support intervals;
- the merge tree gives a hierarchy of transitions.

This is a principled generalization of local peak prominence. A small shoulder that never forms a persistent component is treated differently from a secondary transition separated by a deep valley.

Because a BARS draw is piecewise polynomial, much of this can be constructed analytically:

1. enumerate exact derivative roots;
2. classify extrema;
3. order critical values;
4. solve exact level crossings;
5. build the one-dimensional merge tree or \(H_0\) persistence diagram without imposing a dense evaluation grid.

The existing extrema and span machinery therefore contains much of the computational foundation for an analytic curve-PH engine.

There is also an especially clean correspondence:

\[
\text{peak set}
\longrightarrow
\text{marked peak process}
\longrightarrow
\text{persistence diagram}.
\]

A peak is marked only by location, height, or local prominence in the simpler BARS readout. PH upgrades it to a globally paired object with birth, death, persistence, and potentially a support interval. The posterior peak-intensity field could therefore be complemented by a **posterior bar-intensity measure** over thermal location × persistence.

## PH must also be draw-wise

The same noncommutation principle applies:

\[
\operatorname{Dgm}(\mathbb E[f])
\neq
\mathbb E[\operatorname{Dgm}(f)].
\]

The right-hand side is schematic because diagrams are not vectors, but the point is that PH must be computed on every posterior curve draw. Those results can then be summarized through:

- posterior probability that a bar exists;
- birth/death/persistence distributions;
- persistence landscapes;
- diagram-space densities;
- bottleneck credible regions;
- intensity measures for random collections of bars;
- posterior merge-tree summaries.

This gives an uncertainty-aware topology of the thermal landscape, rather than PH performed on a single smoothed estimate.

The theoretical advantage of the joint curve posterior is that these are genuine pushforwards of one probability law. The dependence among different temperatures is already encoded in each draw. A disjoint sweep cannot supply this without introducing another model afterward.

Under the relevant tameness and stability assumptions, a functional credible band also has a direct topological interpretation:

\[
\lVert f-\widehat f\rVert\_\infty\leq\varepsilon
\quad\Longrightarrow\quad
d_B\!\left(\operatorname{Dgm}(f),
\operatorname{Dgm}(\widehat f)\right)\leq\varepsilon.
\]

Thus curve uncertainty becomes a bound on uncertainty in the topology of the inferred thermal landscape.

## Window boundaries become topological information

Including \(T=0\) and \(T=1\) from the first schedule is important beyond uniform coverage. It makes boundary behavior visible throughout inference.

A feature whose superlevel component remains attached to a boundary may represent:

- a transition outside the chosen thermal window;
- a truncated peak;
- a persistent ordered or disordered regime;
- an inadequately bracketed feature.

That can be represented explicitly using boundary-censoring semantics, relative persistence, or extended persistence, rather than being handled only as a peak-detector “clip” flag. Subsequent complementary schedules refine the interior without ever losing knowledge of the enclosing domain.

## The broader SPC thermal landscape

The scalar analytic curve is the first substrate, not the last. ThermoMapper already produces several potential functions and fields:

\[
\chi(T),\quad C*v(T),\quad H*{\mathrm{label}}(T),
\]

along with

\[
A_e(T),\quad G_e(T),\quad C_e(T),\quad L_v(T)
\]

for edge and node quantities.

This creates several progressively richer topological objects:

- **Scalar curve topology:** persistent extrema and transition hierarchy.
- **Multi-channel thermal trajectory:** relationships among susceptibility, energy, entropy, and structural signals.
- **Edge/node field topology:** functions on graph × temperature rather than on temperature alone.
- **Thermodynamic bifiltration:** temperature and edge threshold as separate ordered axes.
- **Structural zigzag:** actual clusters and Mapper nerves merging and splitting across temperature.

The analytic scalar curve can organize these larger objects without reducing to a scheduler. Its persistent components define thermodynamically meaningful regions in which the graph-valued fields can be examined more deeply. Targeted SW clustering then becomes a conditional structural interrogation of particular posterior-supported thermal regimes.

So the inversion is stronger than “feature extraction before clustering”:

> First infer a probability distribution over complete thermal landscapes; then use analytic geometry and topology to determine what structures, regimes, and transitions exist; finally deploy SPC clustering as a targeted instrument for examining those inferred structures.

That is a genuine departure from the old grid-sweep heuristic. The curve is not merely a better interpolation of the sweep—it is the posterior geometric substrate on which the next generation of SPC analysis is built.

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
