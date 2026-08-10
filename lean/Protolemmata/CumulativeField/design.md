# CumulativeField — recovered local cumulative-linkage construction

Status: design-stage protolemmata. This is a ThermoMapper-native construction,
not an SPC sampler and not a corrected claim of Potts correspondence.

## Mathematical object

Let `G=(V,E)` be a finite base graph carrying nonnegative symmetric couplings
`J`. For a vertex `u`, define its directed local cumulative curve

\[
L_u(t)=\sum_{v:\,uv\in E,\;J_{uv}\le t}J_{uv}.
\]

The directed score of edge `uv` from endpoint `u` is

\[
H_{u\to v}=L_u(J_{uv}).
\]

Filtering by `<=` rather than sequentially sorting is essential: every edge in
one coupling tie class receives the same score, so the construction is
permutation invariant.

The live method has three explicit stages:

\[
\text{local calibration}
\longrightarrow
\text{endpoint symmetrization}
\longrightarrow
\text{graph cut and }\pi_0.
\]

The exponential affinity link is an optional coordinate map on the final
score; it does not supply a Potts interpretation.

## Calibration remains an open design choice

Raw cumulative mass is the first concrete candidate:

\[
H^{\mathrm{raw}}_{u\to v}=L_u(J_{uv}).
\]

It preserves coupling magnitude but also carries local total-mass and degree
effects. Locality alone therefore does not establish density invariance.
Alternatives to compare before fixing the live semantics include:

- cumulative mass fraction `L_u(J_uv) / L_u(+infinity)`;
- local rank or empirical-quantile fraction;
- a locally scaled coupling followed by cumulative mass;
- a typed calibration supplied by the graph artifact.

The graph compiler declares coupling units, gauge, and calibration evidence.
The clustering consumer chooses the cumulative calibration; it must not infer
semantics from a graph-construction option name.

## Endpoint symmetrization

Once the two directed scores are commensurate, form one undirected score before
applying a nonlinear affinity link:

\[
H_{\mathrm{mut}}=\min(H_{u\to v},H_{v\to u}),
\]

\[
H_{\mathrm{bal}}=\frac{H_{u\to v}+H_{v\to u}}2,
\]

\[
H_{\mathrm{inc}}=\max(H_{u\to v},H_{v\to u}).
\]

The balanced rule is an arithmetic mean in cumulative-hazard space. Under
`p_T(H)=1-exp(-H/T)` it becomes

\[
p_{\mathrm{bal}}
=1-\sqrt{(1-p_{u\to v})(1-p_{v\to u})}.
\]

All three are parameter-free and satisfy

\[
H_{\mathrm{mut}}\le H_{\mathrm{bal}}\le H_{\mathrm{inc}}.
\]

Consequently their threshold graphs form a diagnostic sandwich. Mutual edges
are bilateral core evidence; Inclusive-only edges expose one-sided local
evidence. A continuously tuned interpolation is deferred until an estimand
requires the extra degree of freedom.

Arithmetic averaging *after* the exponential link is a distinct research
mode. Its effective score is the temperature-dependent soft minimum

\[
-T\log\frac{e^{-a/T}+e^{-b/T}}2,
\]

so edge rankings may change with temperature. It must never hide behind the
unqualified name `Mean`.

## Cut coordinate and filtration

For `0<theta<1`, define

\[
\kappa(\theta)=-\log(1-\theta),
\qquad
\tau=T\kappa(\theta).
\]

Then

\[
1-e^{-H/T}>\theta \iff H>\tau.
\]

For every temperature-independent score-space symmetrizer, `T` and `theta`
therefore parameterize one cut axis rather than a substantive bifiltration.
Raising `tau` removes edges, so connected components yield the associated
merge hierarchy directly.

For the deterministic solver, the exact critical ladder can be computed once;
an acquisition bracket is then an observation/readout window rather than part
of the underlying hierarchy. Full stochastic SPC retains its separate physical
temperature bracket and BARS boundary-adequacy loop.

## Global control and pooling

At a common level `t`, summing local curves gives the directed global pool
exactly:

\[
\sum_u L_u(t)=G_{\mathrm{directed}}(t).
\]

On a symmetric graph with symmetric couplings, the directed pool is twice the
unique-edge pool. Under homogeneous local cumulative curves, global and local
scores differ by a fixed scale and induce the same hierarchy after cut-axis
reparameterization.

This is a comparison/control theorem, not a limit assertion and not a
vindication of the 2020 paper. Heterogeneous local curves are precisely the
regime in which the recovered method can differ from global linkage.

## Lean theorem slate

### Verified core

- `affinity_gt_iff`: arbitrary-threshold exponential cut reduction;
- `cutScale_half`: the conventional `log 2` slice;
- `cutGraph_antitone`: higher cuts remove edges;
- `cutGraph_mono_score`: larger fields add edges;
- Mutual/HazardMean/Inclusive symmetry and score sandwich.

### Active obligations

- monotonicity of each local cumulative curve under nonnegative couplings;
- exact directed-pooling identity and the factor-two undirected corollary;
- homogeneous-field hierarchy equivalence;
- graph sandwich and induced component-partition refinement;
- equivalence between the cut-graph component hierarchy and weighted
  single-linkage;
- calibration covariance under positive coupling rescaling;
- bracket covariance in the native `tau` coordinate;
- counterexample showing post-link arithmetic mean can reverse edge order as
  temperature changes.

## Naming and provenance

The mathematical namespace is `CumulativeField`. A finalized clustering
implementation may be named `LocalCumulativeLinkage` once the calibration
contract is settled. The paper name is confined to `PKWangAudit`; constructive
modules may cite the audit as provenance in design documentation but do not
import it or inherit its claimed SPC semantics.
