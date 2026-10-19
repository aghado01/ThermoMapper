
# If I'm right about SPC and BARS playing well together then i should see a meaningful change in the economics of SPC applications. everything ive seen about SPC since domany's work have bemoaned the speed and doing my sparse BARS fitting would open up new use cases for an underutilized and powerful method, which is really what i want to be able to do with this separation

That cut is the right SPC package, and it is the experiment the separation is for. It is not yet a measured change in the economics.

The engine already takes a coupling graph, not a geometry. Swendsen–Wang and PKWang are alternative producers of the same affinity currency. PKWang is a deterministic cumulative-energy ladder: affinities only, temperature in $H^{\mathrm{cum}}$ units, not on the Swendsen–Wang scale. Leaving it in the research repo is correct for this package.
[User note: PKWang path is going to be pulled out of SPC's clustering module and become its own module with some additional development. Depending on how it works out in testing, it maybe included as a feature of SPC python/numba implementation or perhaps its own small package separately]

A BARS-on-$\chi$ demo that also offers a closed-form surrogate confounds the thing you want to measure, which is the cost of the stochastic sweep. Riemannian geometry is likewise not an SPC dependency. It lives in the graph compiler’s non-Euclidean bandwidth calibration.
[user note: Not quite true, an original motivation for ThermoMapper project has always been to extend SPC to work in non-euclidean curved spaces and manifold-valued data. But this is the research and may not be in the V1 python port of SPC.]

What the simplified managed package actually ships, from docs/Graphs/GraphConstruction.md and docs/Clustering/SpcEngine.md:

Euclidean graph compiler. GraphCompiler.Build is a distance-over-indices pipeline: topology, filter, repair, refinement, projection, out to a symmetric weighted CsrGraph. The Euclidean path is k-nearest neighbors, MST-min repair (Borůvka bridges only), a Gaussian or Laplacian kernel, bandwidth from the per-node k-th neighbor, and the Blatt–Wiseman–Domany gauge $\kappa = 1/\hat K$ so $J = \hat K^{-1}\exp(-d^2/2a^2)$. That gauge is the original method, and it should be the default rather than a hidden option.

Swendsen–Wang on that graph. Fortuin–Kasteleyn bond probability $p_e = 1-e^{-J_e/T}$, q-state Potts (their default is the original $q=20$), affinities plus alignments and co-membership, a fixed or auto temperature grid, and the susceptibility profile that locates the superparamagnetic band.

The small graph substrate that compiler needs: CSR, the Euclidean kernels, the observables. Not Maths.Geometry, not heat-kernel bandwidth on hyperbolic or spherical manifolds, not the H1 protect-set, not Mapper.

“A better Euclidean compiler” is not specified as a defect in those docs. The Euclidean path already exists. Unless you have a named failure — hubness, the bandwidth estimator, the MST repair adding bridges you do not want — the package should ship that path, not a second compiler invented for the split.

The speed complaint since Domany is real, and it is about Monte Carlo times the temperature grid, not about a missing C kernel. Blatt, Wiseman, and Domany (PRL 1996; Neural Computation 1997) claimed the opposite of a slow method: Swendsen–Wang equilibration scales with $N$ and is independent of embedding dimension, because the graph is sparse. Later work treated the sweep as the bottleneck anyway. O

tt and Stoop (Phys. Rev. E, 2006) proposed belief propagation specifically to accelerate the averages against Monte Carlo, and then had to deal with critical slowing down. Spike-sorting practice often swapped in Wolff for simplicity. Yelibi and Gebbie’s fast SPC (arXiv 1810.02529, Physica A 2020) exists as a ground-state likelihood surrogate because finite-temperature annealing is expensive. The PyPI package spclustering is a fixed-bracket wrapper of Domany’s C. None of those is sparse thermal experimental design. That is the gap, and it is unoccupied.

Your own cost formula is what BARS can change. Stage 1 is $O(|\mathcal T|\cdot R\cdot(B+C)\cdot(N+\mathrm{nnz}))$, spread over workers. Stage 2 is one equilibrium pass at the chosen temperature, and at the documented defaults — on the order of 6,000 equilibrium draws against about 1,200 per sweep task — that pass can dominate wall time. Graph construction is paid once per dataset and is not what that literature complains about. A BARS fit on a one-dimensional thermal curve does not make a draw cheaper. It can only cut $|\mathcal T|$ and decide where the equilibrium budget is spent: a sparse reconnaissance, a posterior curve, per-draw peaks and spans, then production draws inside the posterior bracket, expanded only where clipping or endpoint mass says the window is wrong. issues/spc-BARS already states that inversion, and it already records the unproven parts. Expansion still costs extra Swendsen–Wang evaluations. Edgewise independent fits may be neither practical nor optimal. The lean approximation-economy pillar — sparse grid plus spline versus dense grid — is explicitly not formalized.

So the two packages are the instrument, not the result. Same Euclidean coupling graph, two schedules: a dense fixed bracket, and a sparse grid plus BARS spans plus refined draws. Report temperatures evaluated, equilibrium draws spent, wall time, and whether $T_{fs}$, $T_{ps}$, and the resulting partition moved. If the chosen temperature is stable and the sweep count drops, that is the economic change, and it is the use case the Domany line never shipped. If the posterior bracket keeps expanding back to the dense grid, the claim fails on that data and you will know it from the demo rather than from the package description. Policy stays in that demo. The regression package remains the full module, with its own compiled backend, and does not get named for this application.59 sources

# yea but im not satisfied with the implementation of the current compiler. its unwieldy and messy. It needs some rewrite and for this port, some streamlining for it to be build more to purpose rather than trying to be the nexus of the much more complex thermomapper project and its ambitions. My lift for SPC would be more focused with software engineering and more mature design

The current compiler should not be ported. It is the nexus, and the mess is structural, not a matter of trimming a few options at the call site.

GraphCompiler.Build is a five-stage pipeline with a strategy interface per stage, driven by one immutable config that also has to serve Mapper, SPRED, HDBSCAN’s planned distance graph, and the persistent-homology protect set. The Euclidean Potts path is one projection at the end of that. Along the way the config carries manifold inference (Poincaré, spherical, Fisher–Rao), a fidelity enum whose Intrinsic branch is the heat kernel plus Van Vleck, spherical intrinsic modes, geometry-dependent bandwidth consistency factors, local mutual proximity, and protected edges that LMP must not reweight. AffinityProjection throws. Distance versus coupling is a fork in the same build. Pathology interrupts are control flow, not a report. That is why it feels unwieldy: every stage is written so a later ThermoMapper ambition can plug in, and the SPC caller inherits the plug points.

It also has defects a purpose-built compiler should not reproduce. The OR filter can leave rows directed, and the CSR builders then drop an edge nominated only by its higher-index endpoint. MST-min re-symmetrizes and sets EnsureConnectedApplied even when it added no bridge. Zero bandwidth is not guarded, so MAD of zero on quantized data yields an empty graph. Bandwidth labels in the health report are inverted. MaxNearZeroEdgeRatio fails on LMP builds. Directed kNN and the path-neighbor refiner are both Θ(n²) distance work, which is the wrong default cost for a package whose selling point is a cheaper sweep.

A compiler built for this package is a function, not a nexus. Points in Euclidean space, or an explicit Euclidean distance, go in. A symmetric CSR of non-negative couplings comes out. The pipeline is fixed and the choices are arguments:

k-nearest neighbors, with mutual or union chosen by the caller, not by a hubness gate that swaps the filter underneath them.
Explicit symmetrization before anything else is built, so the directed-row bug cannot recur.

Optional MST-min repair, reported as the bridges actually added, not as a flag that means “the stage ran.”
One kernel, Gaussian, one bandwidth rule, and the Blatt–Wiseman–Domany gauge $\kappa = 1/\hat K$ so $J = \hat K^{-1}\exp(-d^2/2a^2)$. Reject zero bandwidth instead of emitting an empty graph.

No metric inference, no fidelity enum, no LMP, no protect set, no distance-projection fork, no pathology exception as the way the build is steered. A health report can exist beside the build. It should not be woven through it.

Do not keep a strategy interface per stage when v1 has one implementation of each. That pattern is what let the research compiler become the place every other module negotiates. The contract with Swendsen–Wang is the CSR and the gauge. Tests for the package compiler are symmetry, non-negative $J$, connectedness when repair is on, the $\kappa$ normalization, and a hard failure on degenerate bandwidth.

ThermoMapper keeps the five-stage compiler. The package gets a smaller one written against the SPC use, and the research compiler is a source of lessons and of the bugs not to copy, not a base class to wrap. That is the software-engineering lift, and it is the right place to spend it: the sweep engine is already graph-in, and the compiler is what currently makes “run SPC on Euclidean data” a tour of the rest of the project.18 sources
