# Ive found myself wondering recently about what I could actually do with ThermoMapper other than an open source project about scientific computing. Its a bit out of place in terms of its language C# and the community that does things like this. I wonder about products, like some kind of website that might draw viewers or users since c# and perhaps c#/f# later . I've got all this work on TDA and clustering and a vis engine that's already in web developer languages. The main differentiator is the rare capabilities from different domains integrated into a unified framework, even if its incomplete.

Thought for 38s

The viewer is already a website in miniature. What is missing is whether the docs ever named a user for it. I am checking SPRED, persistence, and the MarkBrain notes before treating any product shape as new.

Add to chat

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

54 sources

Interrupted

# [interruption] please focus on reading the materials for what they are which is an unfinished project with an uncertain future beyond completing and validating implementation details. I didnt say it was a product and I am trying to think about what to do with here. Don't be pedantic

The browser half is already a website in pieces. VizCoreSmoke writes self-contained HTML: EyeTorus, the crescent-ellipsoid, the Möbius-ellipsoid, data embedded, Three.js from a CDN, orbit and zoom, a photo mode. VizApi already hosts a live editor at /eye, where the geometry is regenerated in C# and the browser only draws, lights, and drags. Hosting those files, plus the editor, is a site. Nobody installs .NET to look at it. The gallery today shows generator labels, not a clustering run, so it is a geometry cabinet rather than a demonstration of the stack. That is still the shortest thing that can draw a viewer, because the objects are unusual and the classical reference is real: EyeTorus is a toroidal reading of the Blatt, Wiseman and Domany 1996 toy.

The notes in issues/viz already describe the next surface, and it is the one that uses the integration. A film runs a real producer, keeps provenance, and can be opened as a demonstration, a diagnosis, or a fork. The order is fixed and matches what is actually in the tree: the EyeTorus graph laboratory, then a pinned SPC sweep on the same fixture, then a barcode linked to its cycles, then Mapper as a membership relation rather than a second point cloud, then the statistical landscapes. src/viz contracts exist and deliberately do not yet call the science modules. viz-core is the reference until those slices land. That is the unfinished part. It is also the part worth putting on a site, one film at a time, only after the producer behind it has been checked.

That sequence is where the rare combination shows up, and it is not a generic TDA page. The concepts doc treats an edge as a temperature-dependent currency: a coupling, then an affinity or alignment, then a cut, a landscape, or a filtration value. Almost nobody else has an interactive sweep of that. A public "heating the eye" page, temperature tied to assignments, susceptibility, and the hierarchy, with the Domany and Iris cases as calibration rather than decoration, is a thing people in that line would open. There is essentially no modern interactive demo of superparamagnetic clustering. The barcode film is the other obvious draw: selecting a bar highlights the cycle in the point cloud, and selecting the cycle highlights the bar. Persistence diagrams alone are common. Persistence with a body is what the viz program is for, and it is what most library HTML exports do not do. Mapper as a relation is the same idea: cover element, local cluster, members, nerve, kept as linked artifacts instead of flattened into anonymous nodes.

SPRED is a different kind of page, and it should stay later. It is implemented as Grassmann annealing of a projection whose barcode stays close to the data's own, with a block variant aggregated by a Grassmann median. A small validated fixture, original cloud beside the projection, both barcodes, the subspace moving, is a viewer magnet for a narrow audience and a poor front door. The same is true of the conditioned-filtration path. The named scientific use is a task-event backbone on a neural trial cloud, not a plot of four location means. That wants data and a finished producer. It is a later film, not the first public page.

C# only matters if the offer is "install this." The academic TDA tools people already use are Python libraries with HTML on the side. They will not adopt a .NET package, and a later F# migration does not change that. It also does not change the site. The language stays behind the host. A NuGet package, or an F# one, is for the person who wants the engine after they have seen a film. The site is for everyone else.

So the concrete things, in the order the unfinished work actually supports:

- Publish the gallery and the eye editor. No new science. It draws people who like strange geometry and gives the project a URL that is not a README.
- When the graph laboratory is validated, make that the first interactive page: the same EyeTorus under Euclidean, mutual-kNN, geodesic, and the curved realizations, with the false bridge visible. That is the question the film already asks.
- When SPC on that fixture is validated, add the temperature sweep as the second page. That is the differentiator a Python Mapper demo does not have.
- Add the barcode-with-body and the Mapper relation only after those producers are checked. Each is a page, not a platform.
- Leave upload-anything, SPRED, and the neural-manifold backbone until the fixture films exist and the corresponding code is no longer partial.

A site of that shape does not pretend the framework is finished. It is the workbench the viz notes already describe, published one validated slice at a time, with the science remaining in C# or later F#.

# please do have a closer look at the TDA capabilities and things under maths. SPC was always a centerpiece and the project was originally created as an SPC implementation idea because its basically non-existent in modern libraries or research despite what i believe remains a uniquely powerful tool. with PH in particular, I have alot more than just persistence diagrams. I used quite a bit of context engineering with research paper lineages and other persistence libraries to build an engine with zigzag persistence as well as starting to build some original idea of conditioned persistence that build on existing literature. There are some bugs here and there for sure but those can be fixed. i also think my bayesian regression work is useful and powerful as well, with a useful connection to SPC that i have in mind. but realworld applications of all these more exotic things

Thought for 1m 8s

SPC is the centerpiece in the code, not only in the origin story, and the exotic pieces around it are real objects with a narrower set of uses than “TDA in general.” The bugs in the status sections do not change that map.

What is built is a Potts model on a coupling graph, after Blatt, Wiseman and Domany. Swendsen–Wang samples Fortuin–Kasteleyn bonds, pe=1−e−Je/Tp_e = 1 - e^{-J_e/T}pe​=1−e−Je​/T, or PKWang evaluates a closed-form cumulative-energy ladder. A sweep produces three edge currencies rather than a label vector: affinity as bond survival, alignment as spin agreement, co-membership when it is asked for. The observables are the ones the phase is defined by: FK susceptibility χ\\chiχ, energy, specific heat, magnetization. The hierarchies are also built: Blatt phases between magnetization-variance peaks, a dense temperature stack, an excess-of-mass walk on a thermal landscape, lineage tracking, and a tree from the mean-cluster-size curve δˉ(T)\\bar\\delta(T)δˉ(T). The chosen temperature today is the heuristic in SpcAnalysis: midpoint of the χ\\chiχ plateau after a 3-point smooth and a 0.8 drop. That heuristic is the seam, not the method.

The reason this stayed out of modern libraries is cost and fashion, not a refutation. Domany’s own line used it on yeast cell-cycle expression, where a susceptibility peak marked a split and the number of clusters fell out of the phase rather than being supplied. The same group put it on coexpression networks, promoter k-mer tables, satellite images, speech, textures, and named collaborations on fMRI and protein structures. HDBSCAN and Leiden took the niche because they are cheaper and ship in Python. They do not give the object SPC gives: one sweep, a stable band in which dense regions are internally aligned and mutually uncorrelated, three currencies, and a hierarchy whose cost axis is temperature. No kkk. On a million-cell graph the sampler still loses to Leiden. On a heterogeneous table where the question is the stable resolution, it has not been replaced.

The Bayesian regression connection is designed and not wired into SpcAnalysis. The regression module is four conjugate families: free-knot BARS by reversible-jump MCMC, with coefficients integrated out; penalized BAPS, including hierarchical and bivariate splines; a Bayesian GP bandwidth; an exact changepoint by dynamic programming. The visible domain use is a curve on a grid: where it peaks, how many peaks, where it changes level, with a posterior on location and span rather than a point estimate. That is the DiMatteo–Genovese–Kass use, and it is why Kass used BARS on neural rate curves, to put an interval on the time of peak firing.

The SPC inversion written down in the design notes is sharper than “smooth χ(T)\\chi(T)χ(T).” Do not run Swendsen–Wang at every temperature. Fit the joint posterior of the curve, stop when the functional confidence is enough, read extrema and spans off the analytic curve, and spend the sampler only inside those windows.

A second construction in the same notes is the flag bifiltration Kβ,θ=Flag{Gij(β)≥θ}K\_{\\beta,\\theta} = \\mathrm{Flag}\\{G\_{ij}(\\beta) \\ge \\theta\\}Kβ,θ​=Flag{Gij​(β)≥θ}, plus superlevel persistence of the curve computed draw by draw, so the diagram of the mean curve is not treated as the mean diagram. BARS should monotonize only quantities that are ideally monotone, such as the FK connection. Susceptibility and the partitions are not monotone, and they belong to zigzag. Spectra and dose-response can use BARS with no Potts model at all.

Persistent homology here is three engines, not a diagram plot. Standard reduction, cohomology with clearing, and involuted homology in the Čufar–Virk style, which returns cycle representatives. Flag and graph Rips, effective-resistance weights, exact, sliced, and entropic Wasserstein. H1CycleEdges turns those representatives into a protect set for graph construction, so a load-bearing loop edge is not rescaled away. The documented bugs are local: the union-find path inverts the elder rule when vertex births differ, LazyRips is partial, bottleneck distance is a stub, and the cocycle payload is a reduced coboundary rather than a true cocycle.

Conditioned persistence is the original construction, and it is in the tree. ConditionedFiltration puts a backbone K0K_0K0​ in from ε0\\varepsilon_0ε0​ and lets similarity edges enter at their distances. An H1H_1H1​ bar is then a chord that closes a loop along that backbone, not a loop of the point cloud in the abstract. ResidualPrior builds the content edges as prediction residuals, r=φ(∣tj−(ti+τ)∣,∣ti−(tj+τ)∣)r = \\varphi(|t_j - (t_i+\\tau)|, |t_i - (t_j+\\tau)|)r=φ(∣tj​−(ti​+τ)∣,∣ti​−(tj​+τ)∣), with a reach bound ∣τ∣≤Δ|\\tau| \\le \\Delta∣τ∣≤Δ. A return the prior expected is born early. A surprising return is born late. The design notes treat a fixed backbone with τ≡0\\tau \\equiv 0τ≡0 as subsuming the SIFTS reading, and the residual band as the prior that shifts birth. Directed, gauge, and sheaf phases are not implemented. This is a different question from ordinary Rips. Ordinary Rips asks what the cloud is. Conditioned persistence asks what remains once a known axis is already there.

Zigzag is built for the case that axis is not a filtration. Fast zigzag cones the sequence and answers it with one standard persistence run, following Dey and Hou. The graph engines track H0H_0H0​ and H1H_1H1​ of a graph that gains and loses edges, with a dense oracle beside them. Codimension-one complexes reduce through Alexander duality to a dual-graph H0H_0H0​. Strong collapse is a reference. Reflection only deletes the newest cell, and the vineyard port is unwired. The missing driver is the one that would read a non-monotone conditioned sweep as components over temperature, merges and splits both. That is the persistent-Mapper-over-TTT end of the SPC story, and it is not built. Graph zigzag H1H_1H1​ is unfilled cycle space.

The applications that match these objects, rather than the TDA brand, are these.

Heterogeneous tables with no known kkk. Expression matrices, promoter or sequence-feature tables, chemical or spectral libraries, the image and texture sets Domany already ran. The deliverable is the temperature hierarchy and the three currencies, not a single cut. BARS sits in front of the sweep: posterior spans on χ(T)\\chi(T)χ(T) or specific heat, then Swendsen–Wang only inside them. That is the practical form of the connection, and it is also how you report the transition instead of the 0.8 cliff.

A known axis plus residual similarity. Reading order through a document or a sequence alignment: the backbone is the order, the chords are embedding similarities, an H1H_1H1​ bar is a return off the line. Trial structure in a neural population cloud: the backbone is the task-event sequence, and an H1H_1H1​ bar is recurrence that the event order does not explain. That is the flagship already named in the MarkBrain notes, on the Bichot-style PFC recordings. Absolute homology of the cloud and homology relative to the events are different questions. The relative time-warp half of that is not built. Temperature itself is the same kind of axis: components and nerves over TTT, which is SPC read topologically rather than by a cut.

Things that gain and lose edges. Sliding windows, contact graphs, windowed functional connectivity, and SPC partitions that split as well as merge. An ordinary dendrogram cannot see a split. Zigzag can. Susceptibility along the sweep belongs here, not in a monotone filtration.

Curves with an unknown number of peaks, with no clustering at all. Spectra, dose-response, a peristimulus rate, a susceptibility trace from someone else’s sampler. BARS, BAPS, or the changepoint dynamic program are the tool, and the posterior on peak location is the product. Neuroscience already used BARS for exactly that.

SPRED is not an application of this stack. It consumes barcode distances and anneals a projection on the Grassmannian. It does not belong in the SPC or backbone story, and it should not be run on a handful of location means.

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

31 sources

so long as the observables are continuous, one can transform the curve to get monotonous i think. the idea is somewhere under lean protolemata or otherwise

Thought for 44s

The earlier split was wrong. The notes already record the transform, and susceptibility does not have to be handed to zigzag just because the raw Monte Carlo curve wiggles.

It is in lean/Protolemmata/fable-bifiltration-lemmas.md, Corrections 2 and 3, and sharpened in bifiltration-bridge-lemmas.md. Not a compiled lemma. Two different transforms, and they do different jobs.

The continuous nonnegative observable is the free one. Specific heat is clipped nonnegative in the engine. Its cumulative is a strictly monotone reparameterization of the index along a fixed path. A strictly monotone reindex induces an isomorphic persistence module, so the barcode is unchanged up to axis labels. The notes call that the CvC_vCv​-warp: it spends sample resolution at TcT_cTc​, and it is topologically inert. Speed is free. The slope through (β,θ)(\\beta, \\theta)(β,θ) is the lossy choice, because one line is a one-dimensional shadow of the bifiltration. That is the continuous-to-monotone step, and it does not require the observable itself to have been monotone.

The other transform is for an estimand the Potts theory says is monotone and a sample that is not. Ferromagnetic spin agreement, J≥0J \\ge 0J≥0, q>1q > 1q>1, is non-decreasing as the system cools. That is spin_agreement_mono. The finite-sample G^ij(T)\\hat G\_{ij}(T)G^ij​(T) is not, so a complex built on the raw estimate is not a filtration: edges flicker and manufacture births. The close is to project each edge trajectory back onto the monotone cone. Isotonic regression is the cheap projection. Monotone-constrained BARS is the one the note calls the architectural payoff, because it is the same engine as the susceptibility-peak readout, and it returns a posterior birth time per edge. The stability certificate is that the noisy module is ε\\varepsilonε-interleaved with the true monotone one, so the barcode moves by at most ε\\varepsilonε in bottleneck distance. Once spin_agreement_mono is the estimand, that projection recovers a property the theory already has, rather than inventing one. Both the cumulative and the per-edge monotonization need the finished trajectories, so the persistence pass sits after the sweep.

What still is not a filtration is a partition that splits as well as merges. The curve does not have to be abandoned because the samples are non-monotone. The reindex and the projection are what make it one.

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

23 sources

# right, so there are some potentially interesting results and fertile ground. Many things ive hinted at in the lean work came out of implementing things i didn't realize would connect, especially the bifiltration and persistent homology ideas with bars, and i've only barely dug in that work because ive been more concerned with validation as well as some important architectural overhaul to be sure things can be trusted. SPC legacy approach has been validated and replicated many of the key findings of the original lineage and i'm very eager to see what kind of speed up and performance gains or lack thereof I get with the most straight forward application of BARS and intelligent sampling. I've wondered about applications of SPC to things that haven't been explored, like perhaps mechanistic interpretability

Thought for 50s

The order you kept is the one the tree shows. The bifiltration and the BARS connection are protolemmata because they fell out of building the engines, and the engines had to be trusted first. That is not a delay on the interesting part. It is why the interesting part is usable.

The connection itself is already written, in lean/Protolemmata/fable-bifiltration-lemmas.md and the SPC–BARS discussion, and it is the kind that only shows up once both sides exist. A filtration on the raw Monte Carlo edge field manufactures births, because G^ij(T)\\hat G\_{ij}(T)G^ij​(T) flickers. The same BARS engine that reads the susceptibility peak, under a monotone constraint, is what puts that field back on the cone the Potts theory already guarantees for the estimand, and it returns a posterior on edge birth rather than a point TcT_cTc​. The CvC_vCv​ cumulative is the other half: a continuous nonnegative observable reindexes the axis for free, and the barcode does not move. Superlevel persistence of a piecewise-polynomial draw is then closed form, derivative roots and level crossings, taken draw by draw so the diagram of the mean curve is not treated as the mean diagram. None of that is dug in. SpcAnalysis still picks TclusT\_{\\text{clus}}Tclus​ with the smoothed χ\\chiχ peak and the 0.8 drop. The Lean ledger still has the obligations. That matches what you said.

The legacy path is the validated object. The fixtures are the lineage ones: the 1995 toy, the 1996 iris, the 1999 iris, with the cuts the papers actually used, spin agreement, co-membership, peripheral capture, magnetization-variance peaks, the superparamagnetic plateau. The notes treat that replication as ground truth for the thermal curve, not only for the partition. That is the right baseline for the experiment you want, because the question is whether a different acquisition of the same curve spends less.

That experiment is also already written, and it is the straightforward one, not the bifiltration. Coverage-preserving reconnaissance on the normalized window: measure the endpoints, then the complementary midpoints, refit BARS on the union each round, stop when the joint posterior of the whole curve is confident, with the claim kept honest as “no feature narrower than this width” rather than “nothing unresolved.” Extrema and spans come off the draws analytically. Only then does the expensive Swendsen–Wang run, inside those windows, with the rich currencies. Two sampler configs, cheap reconnaissance and production. The first comparison is the sparse posterior against the dense classical curve you already trust: coverage, miss rate, draws. Whether that is a speedup is exactly the measurement, and the notes do not predict one. PKWang is a different speed path. It is a cumulative-energy ladder, not the same temperature axis, so it does not answer this question. PeakIntensity and SpanCoverage already exist on the BARS result. The missing piece is the integration layer that feeds the sweep and stops it.

Mechanistic interpretability is not an experiment in the notes. It is a reasonable later cloud, and the reason it fits is the object SPC returns, not the word clustering. A residual stream or a set of contexts is a point cloud with no known kkk, heterogeneous density, and a resolution question: at what coupling scale do features cohere, and where does a susceptibility peak mark a split. That is the gene-expression question Domany actually used SPC for, moved onto activations. The three currencies matter more than the labels. Affinity and alignment on a token or feature graph are co-activation structure. A hierarchy whose cost is temperature is a splitting story, which is closer to polysemanticity and feature splitting than a flat clustering is. An SAE wants monosemantic directions. SPC does not. It wants the stable band. Those can sit on the same cloud without being the same method.

The other half is already named, and it is the conditioned-filtration primitive rather than a new clustering idea. Token position, or layer order, is a backbone. Residual similarity is the chord. A return the prior expected is born early. A surprising return is born late. That is a probe of trajectory structure in a contextualized embedding, not a feature dictionary, and it is the same shape as the task-event backbone on a neural trial cloud. The Grassmann comparison you already built is the third piece, and it is the one that matches coding-subspace questions. A Pearson of a few mean vectors is the wrong geometry for that, which is the point of the Li et al. reading.

None of that has been run, and it should not be the next run. The next run is the BARS acquisition against the replicated fixtures. If the sparse posterior covers the dense curve, the same acquisition is what you would point at an activation graph. If it misses narrow peaks, you learn that before spending the interpretability cloud.

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

![](https://www.google.com/s2/favicons?domain=github.com&sz=256)

22 sources

There are two formulas, and they are not the same estimate.

Order the directed scores on an edge so \(a\le b\). The tuning parameter is the weight on the larger one:

\[
H(\alpha)=(1-\alpha)\,a+\alpha\,b=a+\alpha(b-a),\qquad \alpha\in[0,1].
\]

Mutual is \(\alpha=0\), balanced is \(\alpha=1/2\), inclusive is \(\alpha=1\).

At a fixed score cut \(\theta\), that edge is born at

\[
\alpha_e=
\begin{cases}
0 & \text{if }a\ge\theta,\\[4pt]
\dfrac{\theta-a}{b-a} & \text{if }a<\theta\le b,\\[10pt]
\text{no birth} & \text{if }b<\theta.
\end{cases}
\]

It is present for every \(\alpha\ge\alpha*e\). The distinct \(\alpha_e\) are the only places the partition can change. Sort them and walk the components. Let \([\alpha*{\mathrm{lo}},\alpha\_{\mathrm{hi}}]\) be the stable interval of that merge tree, the excess-of-mass interval if you walk a height, otherwise the longest run of \(\alpha\) on which the selected partition does not change. The analytical estimate is the midpoint

\[
\hat\alpha=\frac{\alpha*{\mathrm{lo}}+\alpha*{\mathrm{hi}}}{2}.
\]

The partition is constant on the whole interval, so the midpoint is a label for it, not a sharper cut. If the interval is a point, nothing between mutual and inclusive is identified. Report \(\alpha=0\) and keep \(b-a\) as the one-sided field.

The other formula does not estimate an \(\alpha\) for the graph. If the two directed scores are noisy readings of one strength, pool them with weights \(w_u,w_v\), incident degree or local mass:

\[
H=\frac{w*u H*{u\to v}+w*v H*{v\to u}}{w_u+w_v}.
\]

In the ordered basis that is an edge-wise weight \(\alpha*{uv}=w*{\max}/(w_u+w_v)\), different on every edge. It is closed form and it is not \(\hat\alpha\).
