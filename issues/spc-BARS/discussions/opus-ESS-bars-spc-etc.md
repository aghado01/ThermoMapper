
# i wonder how the MCMC data augmentation might bear on ThermoMapper's SPC sampler conceptually


SPC is Superparamagnetic Clustering — Blatt–Wiseman–Domany, with Swendsen–Wang cluster updates on a Potts model. That makes the DTA connection much more concrete than I expected. Let me check how the couplings are built.
Confirmed: `J = exp(-d²/2δ²)` with a **single global bandwidth**. Let me check whether the scaler already supports per-point bandwidths.
Good question — and the answer is sharper than I expected, because ThermoMapper already has the *obvious* version and DTA's actual move is the non-obvious one.

**What's already there.** `LocalMutualProximityScaler` re-weights couplings from local neighborhood counts instead of a global MAD bandwidth, precisely so long MST-bridge edges don't collapse under a global Gaussian. That's static homogenization: a function of observed geometry, applied once as a pipeline stage before any sampling. It's the Zelnik-Manor/mutual-reachability family, and it's the right fix for *density* heteroscedasticity.

**What DTA adds beyond that** is the part in its abstract that's easy to skim past: the augmented data may be a deterministic function of *latent* data and *unknown parameters*, not just observed data. LMP can't be that — it's computed before the chain runs. So the DTA transposition isn't "scale locally," it's "let the homogenizing transformation depend on the current state of the sampler."

**And for SPC there's a specific heteroscedasticity that LMP structurally cannot reach.** Geometric local scaling equalizes coupling magnitudes, but the quantity that actually varies across a real dataset is each region's *local critical temperature*. A single global `T` puts some regions below their `T_c` (already frozen) while others sit above theirs (still disordered) — so the superparamagnetic window isn't one window, it's a different window per region. That's a property of the model at temperature, i.e. of the latent configuration, which is exactly DTA's domain and exactly not LMP's.

**The mechanism-level payoff would be mixing, not just cluster quality.** Swendsen–Wang's efficiency lives in the bond representation, `p_ij = 1 − exp(−J_ij/T)` between aligned spins. Heteroscedastic `J` gives pathological bond probabilities — clusters that either span the graph or never nucleate — which is the classic SW slowdown. Homogenizing `J` in the sampler's own coordinates is directly DTA's advertised benefit (faster convergence of Gibbs/EM), transposed to cluster MC. And ThermoMapper already measures the needed signal: `SusceptibilityKind`, `SweepCurves`, `LabelEntropyCurve`.

**Here's the trap, and it's one we already have a theorem for.** An adaptive scheme that lets local scale move jointly with the configuration is precisely what You's scale-selection paper (`2605.08001`) proves degenerate — the scale runs to the boundary and the objective collapses to a marginal. In SPC terms that's the failure everyone recognizes: adaptive temperature melts everything into one cluster or shatters it into singletons. Same theorem, different clothes.

**Which yields the genuinely interesting conclusion: SPC's existing design is already the safe alternative.** A temperature *sweep* that characterizes behavior across scales, rather than optimizing for a best scale, is exactly alternative #1 in that paper — the sensitivity path, with its uniform consistency and functional CLT. So the paper retroactively justifies why SPC sweeps instead of fits, and the DTA-flavored adaptive version is the tempting move that needs one of the paper's other two constructions (marginal radial calibration, or the bounded balance equation) to be safe.

There's also a tidy structural symmetry worth naming: DTA transforms the *data* to homogenize variance; `ScaledManifold`/`H_α` transforms the *metric* to the same end. Those are one operation from opposite sides, and ThermoMapper has both in the tree already.

**The zero-stretch application, if you want something immediately usable:** `2605.03266`'s intrinsic ESS applies to the SPC sampler with no assumption-bending at all, because SPC *is* MCMC SPC's output summary is a spin–spin correlation matrix, so "coordinate-wise ESS is not canonical under reparameterization" is the exact pathology the paper names. 

The caveat about geodesic Gaussian kernels failing positive-definiteness on curved spaces would bite if the kernel were built on manifold-valued state.

---

# yea lets capture the thermomapper-direction asides for later, this point about ESS, the discussion about scatter/covariance/weiszfeld etc, and others as distinct topic briefs

Now cross-linking the two tracks from the roadmap.
Four briefs written to `D:\aipithicus\aipithicus-issues\thermomapper\briefs\`, matching the house format from `quantized-abliteration-surgery.md` (status / date / evidence checkout / primary question / executive finding), with Obsidian wikilinks between them and a cross-reference table added to the roadmap.

| Brief | What it holds |
|---|---|
| `intrinsic-ess-for-spc` | The ESS point — the one result that applies with zero assumption-stretching, since SPC *is* MCMC. Includes why it's not a generic diagnostic (spin–spin correlation is a matrix, so no canonical scalar to take an ESS of) and the geodesic-Gaussian positive-definiteness trap. |
| `weiszfeld-scatter-vs-sandwich-covariance` | The scatter/covariance distinction with the radial-moment table, the `WeiszfeldScatter.cs:5` breakdown-point correction, naming discipline, and the three-rung metric hierarchy. |
| `scale-selection-degeneracy` | The theorem that relative scale isn't identifiable from a median loss, its three sanctioned alternatives, and the two places it binds — SPC's adaptive temperature and ABLITERALUS R7. |
| `dta-adaptive-homogenization-spc` | The DTA transposition: what LMP already does vs. what DTA adds, why the real heteroscedasticity is per-region `T_c` rather than density, and that the payoff is mixing speed. |

**Two placement decisions worth flagging.** I created a new top-level `thermomapper/` area rather than filing these under `abliteralus/`, since you framed them as a separate track — easy to move if you'd rather they sat with the abliteralus briefs. And I put the Weiszfeld and scale-selection briefs there too even though both span the projects, because both originated in the ThermoMapper review; the roadmap's R2 and R7 rows carry the actionable versions and now link across.

**One substantive thing surfaced while writing them up.** The scale-selection and median-of-means papers together imply there is no single correct `α`: the MoM paper offers *two* calibration rules, a robust block-scale one and an inference-optimal one. So "which scale is right" depends on whether you want robustness or inference, and that choice is a prerequisite to picking among the three alternatives recorded in the scale-selection brief and it's a decision R7 will have to make explicitly rather than inherit.

The vault isn't a git repo (it's an Obsidian vault), so nothing to commit there; the roadmap cross-link is committed in abliteralus as `b6f309b`.

---

# I'm curious to cross-examine our cross-examination of kisungyou's work and my SPC implementation given my planned SPC x BARS work, which resembles adaptive, but is decidedly not, and instead my application of BARS to SPC thermal curves is about iterative RJMCMC fitting on increasingly dense evenly-spaced SW simulations in the normalized bracket. see thermomapper/lean protolemmas, enthymemes and lemmas on the subject

Important correction to register before I look: what I flagged as the degeneracy trap was scale moving *jointly with the latent configuration*. A fixed evenly-spaced temperature design with RJMCMC curve-fitting on top is a different object entirely. Let me read your Lean material.
`BARS.lean` says to work from the protolemmata, not memory. Let me find them.
Now the bracket and densification discussion — that's the part your framing turns on.
The bracket loop is expansion-only with an explicit identifiability limit — that structurally forecloses the degeneracy I raised. One last check on my strongest remaining claim.
This substantially revises two of the briefs I wrote last turn. Let me capture the cross-examination and correct them.
Read `BARS.lean`, the readout obligations, the brackets discussion, and the `bars/` implementation. You're right that it isn't adaptive, and the design is better defended than "not adaptive" — I'd have kept worrying about the outer bracket loop, and it turns out that's closed too.

**Where my caution actually fails.** Temperature is a design variable on a deterministic nested schedule (`S₀ = {0,1}`, complementary midpoint refinement, refit to the accumulated union `U_r`), and the RJMCMC is over spline knots — a *representation* of the curve, not a scale competing with the estimate. Nothing is jointly minimized over location and scale, so the degeneracy mechanism has no purchase. 

The one place a scale genuinely does get adapted on evidence from its own fit is the bracket loop, and the degenerate direction there would be *narrowing until everything looks resolved*. **Expansion-only forecloses it** — monotone expansion on a log axis, all prior observations retained and remapped, nothing extrapolated into newly exposed interval, plus the identifiability limit stated outright. That removes the degenerate direction which is stronger than any of You's three alternatives requires. I've corrected `scale-selection-degeneracy.md`.

**My DTA proposal narrows to almost nothing — for a good reason.** 

- The design already homogenizes twice over: the normalized bracket *is* a deterministic reparameterization of the temperature axis.
- free-knot adaptivity absorbs residual structure. Making the representation adaptive is cleaner than transforming couplings and needs no latent-state feedback. What survives is physical and narrow: **critical slowing down**. 
- SW autocorrelation time rises near the transition, so observation noise on the thermal curve is heteroscedastic in `T` and worst exactly at the feature being localized. That heteroscedasticity lives in the observation model, not the coupling field. Narrowed in the DTA brief.

**The ESS claim gets substantially stronger, and lands on a seam you already built.** It's not a diagnostic — it's the **weight**. `IObservationModel.LogMarginalLikelihood(design, y, weights)` takes per-point weights, and `WeightedNormalModel` documents `W ≡ I` as the *non-robust* case, with the marginal likelihood explicitly conditional on `W`. Feed nominal precision from raw sweep count and you understate uncertainty precisely where autocorrelation is highest — so the posterior is most overconfident about peak location exactly where the science depends on it, and it fails invisibly, because a confident-looking peak is what you expect to see. The coordinate-free property is load-bearing here rather than decorative: susceptibility comes from a spin–spin correlation matrix, so there's no canonical scalar whose classical ESS you'd take.

**Two things I'd raise that aren't in the notes.**

The stopping rule references "absence of unresolved spans near the current schedule spacing." After a bracket expansion the accumulated union is deliberately non-uniform in normalized coordinates — prior points compress toward the interior, new sampling goes to the newly exposed edges. That's correct behaviour, but it makes "the current schedule spacing" ambiguous, and it sits directly on a stopping criterion. The self-audit should use the **local** gap at the feature's location: otherwise a sparse edge inflates the global gap and a resolved interior feature reads as unresolved, or a dense interior deflates it and an edge feature reads as resolved when it isn't.

And MP-2 (`count_of_mean_undercounts`) and SP-3 (`fwhm_of_mean_ne_mean_of_fwhm`) are the same theorem class as the scatter-versus-sandwich distinction — a nonlinear functional not commuting with expectation. 

You've formalized the general principle and enforced it at the readout. 

The ABLITERALUS roadmap sits on the same principle at R4 and R8 without saying so, and its per-layer treatment is vindicated by it rather than just convenient. Want me to write that discipline into the roadmap explicitly? It's a small edit and it'd give R4/R8 a stated reason rather than an implicit one.

New brief at `thermomapper/briefs/spc-bars-cross-examination.md`, with the two corrections applied in place.

---

# its important to discern the uniform iterative sampling, refitting on the net set of sampled points in the curve fitting process from the dense simulations that run on analytically identified regions or features of interest in order to do the subsequent clustering. because SPC x BARS subverts the classical "sweep" approach by using its fixed-grid strategy to identify the analytic areas of interst, and then densely sample those

That distinction is load-bearing and I flattened it. Let me get it into the brief properly, because it changes where two of my findings actually land.
Registered, and it relocates two of my findings.

**The two regimes are doing unrelated jobs.** Regime A is uniform nested grid → refit BARS on the accumulated union `U_r` → joint posterior over the *whole* curve, on an economical SW config. Regime B is dense, expensive SW at analytically identified ROIs, and it's the only place clustering happens. "Densify" means grid refinement in A and replica/equilibration depth in B — unrelated operations, and an argument that slides between them is wrong wherever it lands. Mine was sliding.

**What I hadn't connected: the closed-form readout is what licenses the sparse grid.** Because peaks and spans come analytically from posterior draws — candidates exactly `{endpoints} ∪ {knots} ∪ ⋃ᵢ {interior derivative roots}`, argmax in a finite closed-form set — feature location is limited by the *posterior*, not by grid spacing. So O1 and O2 aren't tidiness results; they're the premise the entire cost inversion rests on. Remove the zero-scan-slop certificate and the sparse grid stops being sufficient.

**My ESS framing was subtly wrong and is now stated correctly.** I'd written that noise is worst "exactly where the peak is being localized," which smuggles in targeting. The regime-A grid is uniform and doesn't know where `T_c` is. The right statement: whichever grid points happen to land near the transition are the *least reliable measurements on the curve*, and they're the ones determining where the peak gets inferred. Uniform design, non-uniform observation quality, degrading toward the feature the fit exists to localize.

**And there's a second argument for per-point weights that doesn't need critical slowing down at all.** The design already plans two SW configurations. The moment observations of different provenance are pooled — most obviously if regime-B measurements at resolved ROIs are ever folded back to sharpen the curve — precision differs by orders of magnitude, and a few excellent points get drowned by many economical ones under `W ≡ I`. Weights are what make heterogeneous-provenance observations poolable at all. Whether B should feed back into the curve is a real open design question; the weighting machinery is a precondition either way.

**The bigger version of my §4 gap, now that the regimes are separate.** Regime A has *two* knobs — grid density and per-point run length — and one stopping rule spanning both. Two different inadequacies present identically as a wide posterior: the grid is too coarse to resolve a narrow feature, or the runs at existing points are too short to be worth fitting. The remedies are opposite, and choosing wrong burns the reconnaissance budget in the way that most resembles progress. 

ESS separates them cleanly: low ESS at existing points means the measurements are inadequate; adequate ESS with a wide posterior *between* points means the grid is. So the self-audit wants a local gap **and** a local ESS — and that's the same quantity from §3 doing double duty, which is the tidiest argument in the pass for computing it at all.

Brief updated in place at `thermomapper/briefs/spc-bars-cross-examination.md`, with each finding now tagged to its regime.

---

# you're still missing the point. the iterative grid samples evenly every iteration, and each fit cycle fits on the superset of all previous points rather than composing previous fits, e.,g. SW [0 1 2 3 4 5] -> fit (F(0,1,2,3,4,5)) -> SW at [0.5, 1.5, 2.5, 3.5, 4.5] -> fit(F(1, 1.5, 2, 2.5...)) -> SW [0.25, 0.75.....]. The number of points on the curve sampled can change while spacing remains even for any given iteration and the superset becomes increasingly dense and confidence on the joint posterior gets better with each iteration. you're also not accounting for replicates where i can run the same simulation N times per point on the curve to improve confidence. The stopping condition is based on joint posterior confidence on the curve fit, and after that threshold is reached, we have analytic features or spans of interest extracted from the final joint posterior estimate/fit which is what informs the second regime. but i could also be misunderstanding, my knowledge of the ESS paper is very limited

You're right on both, and the second one substantially deflates my ESS claim. Let me correct the brief rather than defend it.
Both corrections are mine, not misunderstandings on your part.

**The spacing finding is withdrawn entirely.** I asserted the accumulated union was non-uniform. It isn't — each iteration adds an evenly spaced complementary set, so the union is a dyadic refinement that stays uniform at every iteration with spacing halving each round: `{0,1,2,3,4,5}` → `+{0.5,1.5,…}` → `+{0.25,0.75,…}`. "Current schedule spacing" is exactly well-defined and the span self-audit against it is sound as written. I invented the premise.

**Replicates deflate the ESS claim, and this is the part worth your attention given you said your knowledge there is thin.** `N` independent runs per grid point give you the observation variance *empirically and model-free* — the sample variance across replicates **is** the precision estimate, with no autocorrelation correction needed between them, because they're independent by construction. That's a better input to `W` than anything derived from an autocorrelation model, and it needs no kernel, no lag window, no mixing assumption. So: per-point `W` is still the right vehicle; replicate variance, not ESS, is the right thing to put in it. My "ESS is the correct weight" was overreach.

**What genuinely survives is narrow but real, and it's the one thing replicates can't do.** Replicate spread measures *precision, not accuracy*. If every replicate at a point is under-equilibrated, they're all biased the same direction — toward the initial condition — and their spread is **small**. The point then presents as tightly measured, gets weighted heavily, and is systematically wrong. No statistic computed across replicates can see this, because it's a bias, not a variance. Two things make it worth guarding rather than hypothetical: equilibration time is longest near the transition, so the bias is worst exactly at the feature the fit exists to localize; and it fails in the direction that looks like success. Per-replicate autocorrelation is what catches it, because it asks whether each chain ran long enough for its measurement to be a draw from the target rather than from its start. That's a burn-in adequacy check on individual runs — much smaller than what I claimed, and the only piece replicates don't already answer better.

**One thing in the design's favour I should name**, since refitting on the superset rather than composing fits removes a failure I'd otherwise have gone looking for: no compounding of approximation error across iterations, early knot configurations don't constrain later ones, and successive fits are directly comparable. That last part is what makes "agreement of features across successive complementary refinements" a real convergence signal rather than a self-fulfilling one — composed fits would have inherited their predecessors' structure and agreed with themselves.

The two-knobs point survives in weaker form and no longer needs ESS: large replicate spread at measured points means the measurements are inadequate; small spread with a wide posterior *between* points means the grid is. Both readable from the fit plus replicate variances.

Brief updated — §3 rewritten, §4 marked retracted with the reasoning kept so the withdrawal is legible rather than silently deleted.

---

# yea but i dont think there's anything we can do better for accuracy given that these are SW MCMC simulations being run at each point. in SPC we have a latent thermodynamical system determined by a potts model. i think the only things that approach improving accuracy are things like how the normalized bracket is determined in the first place. accuracy is also something we strive for in regime B after coarsely identified regions of interest *induced by curve fitting and analysis of the fitted curve* guides the final schedule of deep simulations with more expensive observables. i dont think we can get much more accurate than that,

You're right, and it exposes something I've been carrying wrongly for several turns: **SW is itself the standard fix for critical slowing down.** Cluster updates drop the dynamic exponent by roughly an order of magnitude versus single-spin-flip — that's the reason to use it. I've been invoking critical slowing down as a live concern about an algorithm chosen specifically to suppress it. Let me correct that and the briefs that lean on it.
Conceded — and the DTA thread is now fully collapsed rather than narrowed, since the residual I retreated to doesn't survive either.

**The correction that matters:** I kept invoking critical slowing down as the physical justification for heteroscedastic observation noise, across the DTA brief, the ESS brief, and two turns of discussion. Swendsen–Wang exists precisely to suppress that effect — cluster updates cut the dynamic exponent by about an order of magnitude against single-spin-flip dynamics. The residual isn't literally zero, SW's exponent is small but positive, but with independent replicates that residual is *measured* rather than modelled, so it needs no special handling. My equilibration-bias point sat on top of the same weak foundation and goes with it.

**And your framing of where accuracy lives is right.** Regime A is deliberately cheap because its job is *location*, not measurement; the ROIs it produces are regions with posterior width, not point estimates, so a small displacement is absorbed by the region rather than propagated. Regime B is where accuracy is bought, with longer runs and more expensive observables, spent only where the curve says it's worth spending. There's no better measurement to be had at a regime-A point than running the appropriate sampler on the latent Potts system.

**One lever I'd add to your list, upstream of the bracket and larger than it.** The bracket selects a window on a fixed system. The **graph** *defines* the system — kernel choice, bandwidth strategy, neighbourhood structure, global-versus-LMP scaling. So the bracket governs accuracy relative to the model, while the graph governs whether the model is about the data's cluster structure at all. You can estimate a thermal curve to arbitrary posterior precision on a graph encoding the wrong affinities, and every diagnostic in both regimes will look healthy: the curve is well-resolved, the spans are prominent, the clustering is stable. It's the one place where a wrong choice yields a confident, internally consistent answer to the wrong question, and neither regime can detect it, because both take the graph as given.

That's presumably part of what the graph compiler renovation is for. Observable choice — susceptibility versus label entropy versus edge currencies — is a smaller third lever, since they differ in signal-to-noise for locating a transition.

Briefs corrected: §2 now records the collapse and why, with an accuracy-levers section, and the DTA brief is marked as retained only so the idea isn't re-derived.

---

# right, and this is why i have a graph compiler engine and work to do on that in addition to my integration with BARS. but first things first, i have already replicated most of the work of BWD on the classical sweeping approach used in SPC, and that has fixed graph settings and known results on known datasets, so i can test my method exhaustively under those constraints

That's the right sequencing — fix the graph at known-good settings and vary one thing at a time. BWD replication gives you something rarer than a baseline, though: ground truth for the *curve itself*, not just the final partition. That enables a validation you can't otherwise get, and it bears on one specific weakness.

Right, and that ordering also makes the BARS validation much stronger than a partition-agreement check, because the dense sweep gives you ground truth for the **curve**. So you can test coverage directly: does the sparse-grid BARS posterior contain the densely measured classical curve at its nominal rate? That's a real calibration check against real ground truth — the thing that's normally unavailable, and the one that certifies the inferential apparatus rather than just one of its outputs. Iris is the good case for feature agreement, since the parity test already documents the hotter versicolor/virginica split that a first-`T_ps` clip hides — MP-2 in the field rather than on paper.

**The thing I'd aim it at specifically.** The stopping rule's self-audit catches features that are *partially* resolved — a posterior-supported span narrower than the current spacing. It cannot catch one that's **entirely unsampled**. A transition of width `w` falling wholly between grid points at spacing `h > w` puts no grid point on it, the neighbours read as baseline, and BARS fits a smooth curve through them. There's no narrow span to flag because there's no span at all: the audit stays silent and the posterior is confidently featureless precisely where a feature exists. 

**Dyadic refinement resolves it eventually since `h` halves each round, but the stopping rule may fire first, because everything visible looks fine.**

No internal criterion can close that — a sparse design can't know what it didn't sample. Which is exactly why the BWD replication is the right place to attack it: with known feature inventories you can measure the miss rate against spacing and convert an unbounded risk into a characterized operating curve. And it suggests the stopping rule should certify **"no features narrower than `w_min`"** rather than "no unresolved features" — honest, falsifiable, directly calibrated by that first measurement, and it turns the audit from a self-report into a guarantee with a stated failure boundary.

Captured in the cross-examination brief. That also gives the graph-compiler work a cleaner starting point later: once the method is characterized at fixed graph settings, anything that shifts when you start varying kernels and scalers is attributable to the graph rather than confounded with the estimator.


# User notes 

1. is there a way to estimate w_min apriori based on the Potts model or with graph diagnostics? maybe something else? Since the susceptibility and other curves are continuous and modeling a physical system when assumptions are met, this seems like it should be possible. This would be akin to estimating something like a Nyquist frequency for latent curves. I'm thinking that perhaps something in the spectral density of the graph or in the shape of the susceptibility curve itself could give us a clue as to the typical lengthscale of features. 

2. Need to separate issues raised about the current design versus planned changes not seen. The current scheduler will probably be updated. The graph compiler will be updated to make the SPC run path more clear and not confounded by ill-concieved features of the graph compiler design. 