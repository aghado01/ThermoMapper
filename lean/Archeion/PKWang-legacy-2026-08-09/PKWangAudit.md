# PKWang audit — correction of the 2020 parallel-SPC construction

Status: active protolemma specification. Paper-specific provenance belongs
here and in `PKWangAudit.lean`; the recovered construction is specified under
`Protolemmata/CumulativeField/`.

Source: [PKW2020.md](D:/aghado01/graveyard/codex-scientiae/bibliotecha/compendia/clustering/PKW2020.md),
*Parallel architecture to accelerate superparamagnetic clustering algorithm*.
The superseded first Lean register is preserved under
`Archeion/PKWang-legacy-2026-08-09/`.

## Scope

This audit separates two questions that the paper conflates:

1. Does its parallel construction sample, or otherwise reproduce observables
   of, the inhomogeneous ferromagnetic Potts model used by SPC?
2. Given the surrogate construction exactly as written, how much of its
   randomized mask-generation pipeline survives analytic reduction?

The first answer is no: the surrogate does not sample the claimed model. The
second is a collapse: at finite `M`, all masks and their averaging reduce to a
single empirical quantile; at population level, even that randomness becomes
an analytic cut. The legitimate Swendsen--Wang chain is not being called
spurious. The paper's replacement consists of iid scalar-energy draws and does
not target that chain.

## Published construction, stated faithfully

For a finite collection of directed nearest-neighbour couplings, sort

\[
J_1 \le \cdots \le J_L,
\qquad
c_n=\sum_{m\le n}J_m.
\]

At temperature `T`, draw iid scalar budgets

\[
E_k=-T\log(1-r_k), \qquad r_k\sim\operatorname{Uniform}(0,1),
\]

match each `E_k` to the *closest* cumulative total `c_n`, mark the prefix
edges unequal-spin and the remaining edges equal-spin, average those masks,
and retain edges whose average exceeds `0.5`.

The literal closest-match language matters. If

\[
b_1=0,
\qquad
b_r=\frac{c_{r-1}+c_r}{2}\quad(r\ge2),
\]

then the rank-`r` mask is active precisely when `E_k < b_r` (up to null
midpoint ties). A floor-prefix reinterpretation instead uses `c_r`; that is a
cleaner surrogate, but it is not the published matching rule.

## Role of the reported simulations

The `M` iid scalar draws are used in the clustering path, not merely for a
numerical verification: each draw selects a lookup mask, the masks are averaged
as the claimed spin--spin correlation, and the `0.5` cut determines clusters.
Separately, Figure 1a compares their histogram with the asserted energy law.

These replacement draws are inverse-transform Monte Carlo, not MCMC
iterations. Thus the precise redundancy claim is: at finite `M`, the mask
generation and averaging collapse to one empirical order statistic; at
population level, the whole randomized surrogate collapses to a closed-form
cut. This says nothing against the genuine sequential Swendsen--Wang chain.

## Audit theorem slate

### PKA-1 — the sampled energy law is not the Potts energy law

For a Potts configuration `S`,

\[
P_T(S)=Z^{-1}e^{-H(S)/T}.
\]

Its induced energy law contains the density of states

\[
P_T(H=h)=\Omega_q(h)e^{-h/T}/Z.
\]

It is discrete, graph- and `q`-dependent, and supported below a finite
Hamiltonian bound. The paper's exponential surrogate is continuous,
`q`-free, and has positive mass above every finite bound:

\[
P(E>H_{\max})=e^{-H_{\max}/T}>0.
\]

The bounded-support contradiction is the preferred first formal theorem; it
avoids unnecessary density-of-states infrastructure while decisively showing
that equation (6) cannot be the pushforward of the finite Potts measure.

### PKA-2 — prefix masks need not be spin-realizable

Same-spin indicators must be transitive. A three-vertex prefix can declare
`01` and `12` equal-spin while declaring `02` unequal-spin. No labeling can
realize that mask. Consequently the construction does not generate the spin
states whose correlations equation (5) denotes.

### PKA-3 — finite simulation collapses to an empirical cut

For samples `E_1,...,E_M`, the reported affinity at boundary `b` is only the
empirical CDF

\[
\widehat G_M(b)=M^{-1}\sum_k \mathbf 1[E_k<b].
\]

For `0<theta<1`,

\[
\widehat G_M(b)>\theta
\iff
b>E_{(\lfloor M\theta\rfloor+1)}
\]

outside sample ties. At the paper's `M=300`, `theta=0.5`, the complete final
partition depends on the 151st order statistic. The `M` generated pairwise
masks and their averaging are unnecessary even at finite `M`.

### PKA-4 — the population cut

As `M` grows, the empirical quantile converges to

\[
T\kappa(\theta),
\qquad
\kappa(\theta)=-\log(1-\theta).
\]

The paper hard-codes the half threshold, not `log 2`; `T log 2` is the analytic
consequence the paper does not state. Under the literal construction it is
compared with `b_r`; under the corrected floor surrogate it is compared with
`c_r`.

### PKA-5 — global pooling is a linkage reparameterization

For nonnegative couplings, after resolving equal values as whole tie classes,
the global cumulative score is monotone in the original coupling. Its threshold
graphs are therefore the same graph family as similarity-oriented,
coupling-threshold single linkage, with only the cut coordinate changed. The
Potts spins, `q`, Monte Carlo, and proposed parallel mask hardware do not
contribute to that hierarchy.

### PKA-6 — directed multiplicity and ties are underspecified

The paper sorts `N*K` directed neighbour entries. Reciprocal neighbours are
duplicated, unilateral neighbours are not, and equal reciprocal couplings can
receive different sequential cumulative totals. A faithful implementation
must state its directed topology, matching tie convention, and zero-prefix
convention. A repaired undirected control must instead count each edge once
and filter complete tie classes.

## Verified landing and remaining work

`Lemmas/PKWangAudit.lean` currently contains:

- the positive exponential-tail fact;
- the non-transitive-mask obstruction;
- the empirical-CDF count reduction;
- the population half-cut corollary.

Still to formalize:

- a finite Potts Hamiltonian and its support bound;
- the probability-law non-equivalence corollary;
- closest-prefix Voronoi cells and the midpoint formula;
- the empirical-count/order-statistic equivalence;
- the tie-corrected global-field/single-linkage equivalence.

## Interpretation boundary

The audit does not claim that exponential cuts or cumulative edge fields are
useless. It establishes that the paper's stated Potts justification fails and
that its spin/mask simulation is unnecessary. Finite-`M` randomness is not
declared nonexistent: it is exactly a scalar empirical-quantile perturbation,
not an MCMC trajectory or a sample of Potts spin states. The independently
specified construction recovered from the cumulative-field idea lives in
`CumulativeField` and owes its own calibration, filtration, stability, and
empirical validation story.
