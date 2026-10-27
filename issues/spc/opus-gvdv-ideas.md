This is a discussion question. First I'm reading what the project already records about SPC's resolution spine and its entropy readout, and checking whether Chapter 14 covers geometry-dependent partitions at all.

Mostly it would show that SPC already has a partition prior. At each temperature, Swendsen–Wang draws partitions from a well-defined law, and that law is a product partition model. The useful questions are what a Bayesian reading lets you do with SPC's draws, and whether adding a likelihood is worth it.

## SPC's own partition law

Swendsen–Wang samples the Fortuin–Kasteleyn (FK) random-cluster measure. Its open-bond pattern ω on the graph has weight Π_{open} p_e Π_{closed}(1−p_e) q^{k(ω)}, where k(ω) counts components, singletons included. On a bond, p_e = 1 − e^{−J_e/T}.

Summing over the bond patterns whose components equal a given partition π, then absorbing the factor for the cut edges, gives:

**P_T(π) ∝ Π_{blocks B} c(B),  c(B) = q · Σ_{connected spanning subgraphs H of G[B]} Π_{e∈H} (e^{J_e/T} − 1)**

That is the book's Example 14.2, a product partition model (Hartigan's), with one difference:
- **In the book**, the cohesion depends only on block size, which makes the prior exchangeable.
- **In SPC**, the cohesion is the block's own connectivity, weighted by its couplings. So the prior is geometric and not exchangeable.

In this picture:
- **q** is the reward per cluster, the role the DP or Ewens concentration α plays.
- **T** sets how cheap it is to cut an edge.

Two facts tie this to code you already have:
- **The baseline-corrected Alignments field is this law's co-clustering matrix.** Given the FK clusters, spins are coloured uniformly at random, one colour per cluster. So ⟨δ(sᵢ,sⱼ)⟩ = 1/q + (1−1/q)·P(i↔j), and your `(G−1/q)/(1−1/q)` is exactly P(i↔j). In Bayesian-clustering terms that is the posterior similarity matrix.
- **On a path graph, the cohesion reduces to a changepoint prior.** Connected blocks are contiguous segments, the cohesion is q·Π w_e, and the law becomes independent Bernoulli changepoints with geometry-dependent rates. That is the same product-partition family Barry & Hartigan used for changepoints. So SPC and your changepoint engine share prior machinery at the base.

What SPC lacks is a likelihood. The data enter only through J. So SPC's "posterior" is a covariate-driven prior with nothing scoring the data given the partition.

## What the framing buys without changing the sampler

Decision theory over the draws you already make.
- **A better point estimate.** Blatt's crisp step thresholds co-membership at ½ and closes by connected components. Your soft-clustering note already names that closure as single linkage, the source of chaining. Bayesian clustering instead picks the partition that minimizes expected loss:
  - **Binder's loss** with equal costs says "link i and j when P(i↔j) > ½", like Blatt's rule, but solves it under transitivity instead of by closure. Same draws, a resolver that doesn't chain.
  - **Variation of information**, the VI loss (Wade & Ghahramani; the SALSO search of Dahl, Johnson & Müller), is the usual alternative.
- **What each loss needs.** Binder's expected loss depends only on the pairwise P(i↔j), so you can compute it from the accumulated currency today. The VI loss needs the partitions themselves, one per draw. That is your retain-don't-discard rule with a concrete payoff.
- **Uncertainty per temperature.** A credible ball of partitions, in VI distance around the point estimate, gives a posterior radius at each T. It is a candidate stability axis for the T-stack resolver, next to the fuzzy-entropy plateau you already track.

The limit to be honest about: each T gives a different prior. Without a likelihood, the Bayesian view can summarize each slice of the T-stack, but it can't choose between them. Choosing T stays a stability criterion, or needs the next step.

## Adding a likelihood

The posterior would be P(π | x) ∝ P_T(π) · Π_B m(x_B), where m(x_B) is a block marginal likelihood.
- **What it buys:** T, q and the number of clusters get learned through hyperpriors and the marginal likelihood, and the SW→f-SPC handoff (log L_c = log L − EN) becomes one coherent posterior instead of two stages.
- **Neighbouring models:** PPMx (Müller, Quintana & Rosner), the distance-dependent CRP (Blei & Frazier), whose partitions are components of random links and are structurally closest to SW bonds, and Dahl's EPA distribution. None is in this book. Chapter 14 covers the exchangeable family plus covariate-dependent species-sampling processes.

Pushing back on doing this:
- **Data counted twice.** In unsupervised SPC the geometry that builds J is the same data the likelihood scores, so the posterior would be overconfident. PPMx and the ddCRP avoid this because their covariates and responses are different variables. A likelihood earns its place only with a second view of the data: the graph built on one set of features, the likelihood on another.
- **Shape assumptions return.** A Gaussian block likelihood brings back the convex-cluster bias whose absence is why you use SPC. A weak or nonparametric block likelihood avoids that but is harder to make conjugate.
- **The sampler breaks.** SW is fast because q^K comes free from colouring components with integer q. A likelihood term breaks the cluster flip, so you'd need Metropolis-within-SW or split-merge moves (Jain & Neal), and low-T mixing is where that hurts. Replacing q^K with Pitman–Yor-style weights for power-law cluster sizes loses SW the same way.

## Two other Chapter 14 pieces

- **The Indian buffet process** at the end of the chapter is the Bayesian model for overlapping membership. It handles feature allocations, not partitions, and it is where to look if overlap becomes first-class beyond soft weights.
- **Gibbs-type and Pitman–Yor priors** are the catalogue for the cluster-count term, if q^K turns out to be the wrong reward for your data. Spike sorting, with a few large units and many small ones, is the case where it could matter.

## My call

The first step is the decision-theory resolver over existing draws. It is cheap, needs no sampler change (Binder runs on the accumulated Alignments), and is falsifiable. Run Blatt's rule and the Binder-optimal and VI-optimal partitions on the same SW draws, on the BWD toy and on ISOLET. If chaining is what costs SPC on ISOLET, this shows it. The likelihood model is worth it only once there's a second view of the data; without one it double-counts.

The FK-to-product-partition derivation and the Alignments identity are standard and I checked them by hand here. The papers cited are from memory, not checked against the texts.

Before answering, I'm checking what the arch-overhaul corpus already settled about resolving the T-stack and the bifiltration, so this builds on it rather than re-deriving.

The FK view from last turn gives the framing: the geometry is the prior. At each T, SPC's partition law is a product partition model whose cohesion is the graph's own connectivity. So Bayesian ideas break the geometry only when they bring in a likelihood with shape assumptions (Gaussian blocks), or a prior that ignores the graph (CRP and Pitman–Yor). Without a likelihood there are three places to put Bayesian machinery, and none of them touches J.

## 1. The T axis: trajectory engines

This part is already in your corpus. `resolution-and-bars.md` §3 has BARS on χ(T), giving a posterior over T* that is the exact pushforward of the curve posterior. The bifiltration note has BARS fitting each edge's Ĝᵢⱼ(T) as a monotone curve, which gives each edge a birth temperature with a posterior attached. One natural extension is per-cluster:
- Run the changepoint engine on a candidate cluster's trajectory across T: its size, or its mean internal connection probability.
- That gives a posterior over the cluster's split temperatures.
- So a split comes with a credible interval instead of an eyeballed jump.

## 2. The decision layer: losses instead of thresholds

This is the new piece, and I think it is the main one for automation.

Binder's loss charges a for a false merge and b for a false split. The pairwise Bayes rule links i and j when P(i↔j) > a/(a+b). Two consequences:
- **Your bifiltration already has this axis.** The θ axis of the (β, θ) bifiltration, read off the baseline-corrected Alignments, is the cost ratio of false merges to false splits. So the bifiltration is (prior resolution) × (loss asymmetry), which gives θ a meaning beyond "a threshold".
- **Blatt's crisp step is a heuristic, not the optimum.** It thresholds at θ and closes by connected components; your soft-clustering note already calls that closure single linkage. The actual Binder-optimal partition respects transitivity without that closure, so it doesn't chain. That makes it a new member of the flat resolver family, next to argmax, threshold-to-abstain and closure.
- **Abstention gets a cost.** Add an abstain cost c and a point goes to `Unassigned` when every assignment's expected loss exceeds c. That is Chow's reject rule carried over to clustering.

That leaves three stated costs (a, b, c) in place of tuned thresholds. This fits your "calibration is not API" rule: costs carry across datasets in a way per-dataset thresholds don't.

Two engineering facts:
- **Binder needs only pairwise P(i↔j)**, which the accumulated Alignments already provide. You could build it today.
- **The VI loss, credible balls and "does block B exist at T" need the partitions themselves**, one per draw. SW already computes component labels every sweep in union-find, so keeping them costs storage, not compute. Thin them, though: N × draws × T-grid adds up fast. For example, 10⁴ points × 1000 draws × 50 temperatures is about 2 GB, versus about 200 MB at 100 draws.

## 3. The prior class: robust Bayes over the window

The T-scan is a family of geometric priors. Robust Bayesian analysis studies quantities that stay stable across a class of priors, and that is a formal version of the superparamagnetic window:
- A feature F (a block, a pair linked, a split) is robust on a window W if P_T(F) ≥ 1−α for every T in W, with Monte Carlo error bands.
- The window itself becomes the largest T-interval on which the loss-optimal partition stays inside its credible ball.

## Pushing back

- **Robustness has four axes:** the choice of T, sampler noise, graph construction (k, local scale, the GraphCompiler) and data perturbation. These tools cover the first two well.
- **The graph is upstream of all of it.** If the graph bridges two clusters, the FK prior reports the bridge faithfully, and every layer above it becomes calibrated about the wrong object. Your corpus already says resolution burden tracks graph quality, and this doesn't change that.
- **Without a likelihood, Bayesian language can't make the choice of T true.** The robust window is a stability criterion, with the same epistemic standing as SPC's own window. What it adds is explicitness, uncertainty and stated costs; it doesn't add ground truth.
- **The likelihood side already has a home:** f-SPC, the inverse estimator. The layers above are all forward-side and stay geometric.

## What I'd build, in order

Each step is useful on its own:
1. **A Binder resolver** on the baseline-corrected Alignments, one partition per T, scored against Blatt's closure on the Domany battery and ISOLET. This tests the chaining hypothesis directly.
2. **Retain thinned per-draw partitions.** That enables VI-optimal estimates, credible balls, and probabilities that a cluster exists at a given T, which together define the robust window.
3. **Cluster split temperatures** from the changepoint engine, feeding the `PartitionHierarchy`/wave_clus resolver that is still unbuilt.

One thing to keep apart: your multiparameter note records the warped path γ through (β, θ) as a working frame, not a settled one. Reading θ as a cost ratio gives one option for it: a fixed-θ slice uses the same loss at every resolution. Letting θ vary along γ would mean the merge/split costs change with scale, and whether that is ever desirable is open.

If you want, I can write this up as a note in `issues/spc-e2e-maturity/` next to `resolution-and-bars.md`.