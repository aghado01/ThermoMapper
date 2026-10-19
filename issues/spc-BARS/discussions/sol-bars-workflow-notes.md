
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
