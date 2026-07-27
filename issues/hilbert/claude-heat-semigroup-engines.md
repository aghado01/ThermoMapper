# Heat-Semigroup Engines — affordability, the subgraph reframing, correspondence

**Status:** engines-facing analysis and forward sketch; discussion-grade, nothing
scheduled. Third doc of the arc — companions:
[claude-HOPE-paper-analysis.md](claude-HOPE-paper-analysis.md) (first read),
[claude-hilbert-synthesis.md](claude-hilbert-synthesis.md) (HOPE × NLD × stack).

**Seeds (Azriel):** (1) the NLD matrix exponential is too expensive for large dense
point-cloud graphs and the affordability question was unresolved; (2) the real
interest is not network-type data but graphs on point clouds; (3) the conceptual
extension — systematically examining *subgraphs within one large graph* rather than
side-by-side whole-network comparison; (4) if the data is self-similar
("fractillitude"), correspondence-*measure* approaches may open inside the subgraph
formulation. These four seeds structure the three movements below.

---

## I. Affordability — the cost anatomy of NLD-type computation

**The exponential is a solved problem at the action level; the probe structure is the
real cost.** NLD as published eigendecomposes each Laplacian (`O(N³)`) and probes with
all `N` basis vectors. Neither is necessary.

**Actions, not matrices.** Everything needed is `exp(−tL)·v` for specific probes:

- **Chebyshev expansion** of the heat kernel in `L`: the basis vectors `T_j(L̃)v` are
  `t`-independent — only the coefficients (modified Bessel functions of `t·λ_max`)
  carry `t`. One basis per probe evaluates the **entire time grid** as reweighted
  sums: NLD's `∫dt` — its whole multi-scale content — costs barely more than one
  flow. Degree needed grows like `√(t·λ_max)`.
- **Krylov/Lanczos** actions: `O(m·nnz)` per probe, near-best polynomial accuracy.
- **Rational/contour methods** (Hale–Higham–Trefethen) + fast Laplacian solvers:
  ~8–16 shifted solves `(L+σI)⁻¹v`, uniform accuracy across `t` — the high-accuracy
  route when polynomials strain.

**The t-axis splits onto owned machinery.** Small-to-moderate `t`: polynomial methods
— and for *localized* probes, push-style heat-kernel algorithms (Kloster–Gleich
family) cost proportional to the **flow's support, independent of N**. Large `t`:
only the bottom spectrum survives `e^{−tλ}` — a LOBPCG job (`Spectral.cs`, already
owned). Crossover governed by the spectral gap.

**Probe reduction — the L1→L2 relaxation.** NLD computes the entrywise-L1 norm of
heat-kernel differences over the time grid; entrywise L1 of an implicit matrix has no
cheap unbiased estimator. The L2 sibling is trace-shaped:

```
‖e^{−tL₁} − e^{−tL₂}‖_F² = tr(e^{−2tL₁}) + tr(e^{−2tL₂}) − 2·tr(e^{−tL₁}e^{−tL₂})
```

and traces of matrix functions are Hutchinson/Hutch++ territory: ~10–100 stochastic
probes instead of `N`. **Declare the rung**: L1/TV is where NLD's bridge-sensitivity
claims live; L2/trace is where affordability lives. A fidelity choice, not a silent
swap. (`tr(e^{−tL})` is the heat trace — the spectral invariant that returns in
Movement III as the dimension screen.)

| quantity | estimator | cost driver |
|---|---|---|
| one flow, full t-grid | Chebyshev shared basis | `O(m·nnz)` |
| localized flow, small t | push (support-adaptive) | independent of N |
| Frobenius-relaxed NLD | Hutchinson + actions | ~10–100 probes |
| large-t anything | bottom-r LOBPCG | already owned |
| exact entrywise-L1 NLD | no cheap estimator | the honest price of TV |

## II. The subgraph reframing — from global metric to local field theory

Wild point-cloud data hands you **one big graph** (or a construction-filtration family
of them), not a population of networks. The extension dissolves NLD's frame: the heat
semigroup becomes a *localized family of diffusion observables on one graph*; global
pairwise distance becomes the degenerate special case. Locality then stops being an
economy and becomes semantics: a probe seeded in `S` at small `t` *is* a local
computation.

**The boundary-condition menu is the mathematical heart.** A subgraph has no unique
Laplacian; the choice is a declared coupling of subsystem to environment — the
(gauge, measure, metric) audit extended by a boundary axis:

- **Island** (induced, degrees recomputed): `S` as a closed system.
- **Dirichlet** (delete complement rows/cols): boundary grounded — absorbing bath.
  Heat *content* `Q_S(t)` under this choice measures how well `S` traps diffusion;
  small-t asymptotics carry boundary-size (perimeter) information; conductance is its
  `t→0` shadow.
- **Schur complement / Kron reduction** (eliminate the complement): the environment
  folded into effective edges — integrating out degrees of freedom, a decimation/RG
  move; exactly preserves effective resistances among retained nodes. Densifies, but
  dense-on-`|S|` is fine when `S ≪ G`.

Closed system / absorbing bath / integrated-out environment: thermodynamic boundary
coupling in graph terms.

**Observable families this unlocks** (all compositions over Movement-I primitives):

- **HKS field** `[e^{−tL}]_ii` over (node × t) — multi-scale behavioral node
  descriptor; clustering nodes by HKS profiles = behavioral identity again.
- **Heat content** `Q_S(t)` — community quality *across scales*, not at one
  resolution.
- **Subgraph autonomy** — flow in `S`-with-environment vs `S`-as-island: how much a
  region's diffusion identity depends on its surroundings. Modules = dynamically
  closed subsystems; a principled community criterion, not a combinatorial one.

**Placement against the literature and the stack.** Markov stability
(Delvenne–Barahona line) already treats communities as diffusion-trapping subgraphs
with `t` as resolution and plateau-detection as model selection — the diffusion twin
of the SPC temperature sweep (honest asymmetry: `t` generates a true semigroup, `T`
does not; the plateau logic transfers, the algebra doesn't fully). Underneath sits
the cleaner bridge: **`e^{−tL}` is the Gibbs propagator of the free (harmonic) field
on the graph; Potts correlations are the interacting cousin.** Diffusion machinery is
the Gaussian sector of the theory the stack already runs — not a bolt-on. (The
precise formal reach of this free↔interacting bridge is an open question worth its
own note someday.)

## III. Correspondence gated by self-similarity

Cross-subgraph comparison (`S₁` vs `S₂`, different node sets) needs correspondence.
The seed: if the data is self-similar, correspondence-*measure* approaches open. The
established name for the object is a **coupling** (Gromov–Wasserstein): optimize
`π ∈ Π(μ₁, μ₂)` minimizing pairwise-distance distortion; the optimal `π` *is* the
soft correspondence, the residual is the distance. Ingredients per subgraph: a metric
(diffusion distance at scale `t` — Movement I) and a declared node measure — the
metric–measure thesis was already GW-shaped; mm-spaces are Gromov's native objects,
and Mémoli's spectral GW compares them through heat kernels specifically.

**The correspondence-fidelity ladder:**

- **Rung 0 — invariant descriptors** (heat traces, spectral measures, HKS
  distributions): correspondence-free, cheap, always available.
- **Rung 1 — functional maps**: correspondences as small `r×r` matrices between
  truncated Laplacian eigenbases, constrained by HKS-type descriptors. Built entirely
  from LOBPCG + Chebyshev primitives. Correspondence-cheap.
- **Rung 2 — entropic GW couplings** (Sinkhorn-solved): full correspondence-measure,
  committal, tractable at subgraph sizes.

**The gate is load-bearing.** GW's failure mode: it *always returns a coupling* —
without an external check, motif discovery hallucinates. Two disciplines:

1. **Screen before coupling** (cheap, trace-based necessary conditions):
   - spectral dimension from the heat-trace slope `tr(e^{−tL}) ~ t^{−d_s/2}`;
     *stability across a t-window* = scale-invariance;
   - entropy susceptibility `−dS/d log t` of `ρ_t = e^{−tL}/Z` (the Laplacian-RG
     literature's diagnostic: flat = scale-invariant window, peaks = characteristic
     scales);
   - **local dimension field**: `[e^{−tL}]_ii ~ t^{−d_s(i)/2}` — the HKS log-slope
     gives pointwise spectral dimension, a field over (node × scale) that maps
     *where* self-similarity holds before any Sinkhorn iteration is spent.
2. **Null ensembles**: GW distortion means "similar" only against a reference
   (rewired / matched random-geometric subgraphs) — validation ground truth from
   outside the estimator, per the independence principle.

**False-fractal hazard (flagged hard):** on kNN graphs over nonuniform samples,
density artifacts manufacture power-law-looking windows — the same trap that produced
collapsed "scale-free network" claims. Power-law statistical hygiene applies, and the
DTM/α measure-correction work sits *upstream* of any fractal claim. Load-bearing, not
polish.

**Three question types, one gate structure (screen → couple):**

- **Local** — nested balls `B(x,r)` vs `B(x,cr)`, same center across scale: local
  self-similarity; the nesting is a filtration, so this sits beside persistence.
- **Regional** — different centers, matched scale: motif discovery / recurring local
  geometry. Mapper cover cells are a natural subgraph family; GW-clustering cover
  cells would enrich Mapper nodes with intrinsic-shape labels.
- **Global** — the coarsening ladder vs itself: the RG fixed-point question.

**Self-similarity as accelerant, not just object.** If coarse structure matches,
match coarse first and refine within matched blocks: hierarchical/multiscale GW,
using the same Schur/coarsening ladder as Movement II — the fractal structure pays
for its own measurement. Ancestry worth crediting: fractal *image compression*
encodes an image by range–domain block correspondences under contraction; this is
its graph-Laplacian descendant. And the HOPE echo: redundancy across neurons (HOPE)
vs redundancy across scale (here) — both measure structure as "how much of the
object is explained by the rest of itself"; compression as structure detection.

**Carrier note:** a coupling is mathematically the stack's soft-assignment currency
(doubly-stochastic-flavored membership matrix) — the carrier type exists. New someday:
one solver family (Sinkhorn/GW = Inference-fit) and the screens as Observables.

## IV. Engine factorization — "Hilbert first-class" in concrete terms

Four primitives; everything above is composition:

1. **Matrix-function actions** `f(L)·v`, `f ∈ {exp, resolvent, polynomial filter,
   spectral projector}`, backend by regime (Chebyshev / Krylov / push / LOBPCG /
   contour).
2. **Declared measure** on graph signals — weighted `L²(V, μ)`,
   `μ ∈ {uniform, degree, density}` — making the sym/rw/1-K̂ menu code-real as
   inner-product choices.
3. **Sub-operator constructor** with the boundary menu (island / Dirichlet / Schur).
4. **Stochastic estimators** — Hutchinson/Hutch++ traces, KPM spectral densities.

Compositions that fall out at roughly a page each: NLD and its L2 sibling, HKS
fields, heat content, Markov-stability-type sweeps, subgraph autonomy, heat-trace
intrinsic dimension (Weyl slope — feeds the DR track), local dimension fields, von
Neumann graph entropy and its susceptibility, spectral graph wavelets, heat-kernel
local clustering, and the metric+measure ingredients for functional maps / GW.

Placement by the location rule throughout; the primitives are concrete types where
they land, and the compositions are consumers. Code-graph analysis remains a
calibration toy; the target is wild point-cloud data — the reframing in Movement II
is what aligns the machinery with that target.

## Open edges

- Quantify what the L1→L2 relaxation loses on the cases that motivated NLD
  (bridge-sensitivity) — is TV's advantage material at point-cloud scale?
- Which boundary-menu entry serves which observable — autonomy wants
  with-environment vs island; what does motif discovery want?
- Screen statistics: what constitutes a "stable" dimension window (power-law hygiene,
  window-fitting discipline).
- GW practical envelope: subgraph sizes, entropic-regularization strength,
  hierarchical refinement depth.
- The free↔interacting bridge (`e^{−tL}` ↔ Potts/FK correlations): how far does the
  formal statement reach beyond the structural analogy.
- Directed / asymmetric variants (non-self-adjoint generators) are out of scope until
  a use case forces them.
