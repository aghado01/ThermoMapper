---
provider: codex
format: Structural
exported_at: 2026-08-10T00:29:46.0908055Z
thread_id: 019fe387-d599-7403-bee7-440e8d2a6416
exchanges: 27
user_label: Aipithicus
models:
  - gpt-5.6-sol
exclude:
  - thinking
  - commentary
  - tool-calls
  - tool-results
  - subagents
  - synthetic
  - timestamps
  - session-markers
  - exchange-markers
---

# so i was just thinking through some of my lean enthymemata. i haven't discharged any of them because i suspected they would change. i want to review what i have against this discussion and examine further D:\aghado01\ThermoMapper\lean\proto-lemmas\Fable5-BARS-revisions.md

You were right not to discharge these yet. My recommendation is to preserve `PKWangA`, retain a few generic analytic lemmas, and revise the remaining interfaces before investing in proofs.

The discussion’s durable insight is that many constructions reuse connected components across parameter changes. Its strongest claim—“one `π₀` functor unifies merge trees, Mapper, zigzags, and bifiltrations”—compresses several genuinely different objects too far.

## Main architectural correction

I would formalize four separate layers:

1. `GraphFiltration I V`: an ordered family of graphs with inclusion maps.
2. `ComponentDiagram`: connected components plus their induced maps.
3. `MergeForest`: the display/event structure derived from a constructible one-parameter component diagram.
4. `Mapper`: pull back a cover, take connected components, then form a nerve.

Mapper can reuse the component machinery, but it is not itself merely `π₀`.

In particular, the proposed “Mapper on a nested sublevel cover degenerates to the merge tree” statement is false for ordinary Mapper. If three connected cover elements satisfy \(U_1\subset U_2\subset U_3\), every finite intersection is nonempty, so their nerve is a filled triangle. The merge tree is a chain. Recovering the merge tree requires the component display poset, an adjacent interval cover/multinerve construction, or an explicit comparison theorem—not the ordinary nerve of an arbitrary nested cover. This agrees with the standard definition of Mapper as the nerve of connected pullback components and with the more careful cover hypotheses used in convergence results: [Carrière–Michel–Oudot](https://jmlr.org/papers/volume19/17-291/17-291.pdf) and the [categorical treatment by Munch–Wang](https://arxiv.org/abs/1512.04108).

The existing C# implementation also supports this separation: [PersistentMapper.cs](/D:/aghado01/ThermoMapper/src/tda/mapper/PersistentMapper.cs) produces framewise nerves, while [NerveDiff.cs](/D:/aghado01/ThermoMapper/src/tda/ph/nerves/NerveDiff.cs:123) uses overlap-based matching and permits both merges and splits. Those matches are correspondences, not automatically functorial inclusion maps.

## Assessment by enthymeme

| Area | Recommendation | Principal reason |
|---|---|---|
| `PKWangA` | Keep | Already proved and structurally sound |
| `PKWangB` | Rewrite definitions first | Current definitions omit graph adjacency and do not match tie/symmetrization behavior |
| BARS v2 | Keep generic calculus lemmas; redesign observables | Improved logically, but still does not model the engine’s prominence-relative span |
| Potts/FK | Redefine probability model | Parameters belong to finite edges with values in \([0,1]\), not arbitrary `V → V → ℝ` |
| Bifiltration | Fix indexing and functor interfaces | Threshold direction is dual; reparameterization needs an order equivalence |
| Stability | Split into several adapters | Current file conflates scalar curves, graph filtrations, and Mapper persistence—and does not compile |
| Ascent | Weaken/restate billing theorem | Current proposed crossing result is mostly a generic partition fact, not a valley theorem |
| Mapper | Define independently from merge forests | Components are shared infrastructure; the nerve construction is additional data |

### PKWang

The useful generalization of the existing threshold lemma is

\[
1-e^{-H/T}>\theta
\quad\Longleftrightarrow\quad
H>T\log\!\frac1{1-\theta},
\]

under \(T>0\) and \(0\le\theta<1\).

The vague “local field tends to global field under pooling” statement can be replaced by an exact finite identity. Define cumulative profiles over actual graph edges:

\[
H_u(c)=\sum_{\{u,v\}\in E,\;J_{uv}\le c}J_{uv},
\qquad
H_G(c)=\sum_{\{u,v\}\in E,\;J_{uv}\le c}J_{uv}.
\]

Then

\[
\sum_u H_u(c)=2H_G(c).
\]

That is a precise double-counting theorem, with any homogeneity approximation stated separately.

Several implementation decisions must be settled first:

- [PKWangB.lean](/D:/aghado01/ThermoMapper/lean/enthymemes/PKWangB.lean:20) currently sums over all vertex pairs rather than edges.
- The implementation uses sequential, rank-inclusive sums after sorting; tied couplings can therefore receive different cumulative values ([MeanField.cs](/D:/aghado01/ThermoMapper/src/clustering/graphical/spc/runtime/core/solver/MeanField.cs:57), [LocalField.cs](/D:/aghado01/ThermoMapper/src/clustering/graphical/spc/runtime/core/solver/LocalField.cs:54)). A filtered-sum definition that includes an entire tie class would be a specification change.
- Mutual and Inclusive symmetrization are `min` and `max`, but Mean averages the probabilities ([EdgeFieldSymmetrization.cs](/D:/aghado01/ThermoMapper/src/graphs/primitives/EdgeFieldSymmetrization.cs:17)). Only `min` and `max` commute cleanly with the increasing exponential kernel; Mean does not admit a single fixed averaged `Hcum`.
- Clustering uses `g > θ`, while the thermal dendrogram searches with `g ≥ θ` ([AffinityThreshold.cs](/D:/aghado01/ThermoMapper/src/clustering/graphical/spc/partitions/strategies/AffinityThreshold.cs:38), [ThermalDendrogram.cs](/D:/aghado01/ThermoMapper/src/clustering/dendrogram/ThermalDendrogram.cs:76)).

The suggested relationship to robust single linkage/HDBSCAN should remain labelled an analogy until an order- or distortion-equivalence theorem exists. The classical robust single-linkage result assumes sampled metric-density geometry, not arbitrary cumulative edge fields: [Chaudhuri–Dasgupta](https://papers.nips.cc/paper_files/paper/2010/file/b534ba68236ba543ae44b22bd110a1d6-Paper.pdf).

### BARS

[BARS.v2.lean](/D:/aghado01/ThermoMapper/lean/enthymemes/BARS.v2.lean:54) makes three worthwhile corrections:

- interval-free Fermat;
- the antitone shoulder theorem;
- separation of global superlevel occupancy from per-peak width.

But it still needs revision before proof work:

- Define both full-topology `interiorPeak` and interval-relative `windowPeak`; the engine counts boundary candidates.
- Globalize spline candidates as derivative roots plus knots and endpoints. A knot maximum need not be a local maximum of either polynomial extension outside its span.
- Replace global half-maximum FWHM with a `PeakBracket` using the operational prominence/drop level:
  \[
  \ell=\text{mergeLevel}+\rho(\text{peakHeight}-\text{mergeLevel}).
  \]
- State count results for prominence-gated peaks, not only weak `IsLocalMax`; plateaus make the current `ncard` contract junk to zero.
- Monotonicity alone does not provide a unique crossing when plateaus exist. Use `StrictMono`, or define canonical left/right generalized inverses.
- Describe closed-form roots as removing scan/iteration error, not all numerical error.

There is also an ordinary draft issue: the Mathlib 4.30 root-finiteness theorem is `Polynomial.finite_setOf_isRoot`, not the name currently used at [line 71](/D:/aghado01/ThermoMapper/lean/enthymemes/BARS.v2.lean:71).

### Potts, bifiltration, and stability

The random-cluster layer should be defined over the finite edge set with \(p_e\in[0,1]\), together with the boundary condition and probability measure. That is the native model in Grimmett’s [Random-Cluster Model](https://www.statslab.cam.ac.uk/~grg/books/rcm1-1.pdf). Raw spin agreement should also be normalized if it is to act like a connection probability with threshold in \([0,1]\).

For the product filtration, make the order explicit as something like \((\beta,\theta^{op})\). A `StrictMono` change of coordinates is not enough for a filtration isomorphism; use an order isomorphism or restrict the result to the image.

Stability should then be split into:

- scalar-function/sublevel persistence;
- fixed-complex edge-birth stability;
- Mapper or zigzag stability.

The current [Stability.lean](/D:/aghado01/ThermoMapper/lean/enthymemes/Stability.lean:23) does not compile directly: `opaque Dgm` and `opaque landscape` require inhabitants for their opaque result types. These should temporarily be `axiom`/`constant` declarations or real structures. Persistence-landscape 1-Lipschitz stability is a reasonable cited assumption: [Bubenik, *Statistical Topological Data Analysis using Persistence Landscapes*](https://www.jmlr.org/papers/v16/bubenik15a.html).

## Build coverage is presently misleading

The aggregate Lean build succeeds, but [Enthymemes.lean](/D:/aghado01/ThermoMapper/lean/Enthymemes.lean:1):

- imports the old `BARS.lean`, not `BARS.v2.lean`;
- omits `Stability.lean`;
- therefore does not expose the failures above.

There are also stale `lean/prelemmas/...` references even though the directory is now `lean/proto-lemmas`, and the `lean/enthymemes` versus `Enthymemes.*` casing is likely to fail on Linux. Status automation should distinguish `PROVED`, trusted/cited axioms, obligations, and unimplemented targets rather than counting only `sorry`.

## Recommended sequence

1. Repair aggregate imports, paths, casing, and status checks.
2. Define weighted graph fields, edge tie policy, graph filtrations, component diagrams, and merge forests.
3. Finish PKWang at arbitrary \(\theta\), including symmetrization-specific theorems and the exact pooling identity.
4. Redefine BARS peaks, brackets, prominence levels, and spline candidate sets; then discharge the generic calculus lemmas.
5. Build the bounded-edge Potts/FK model and monotonicity.
6. Define the ordered bifiltration and initially take only \(H_0\).
7. Split and repair stability.
8. Define Mapper as pullback components plus nerve, then prove comparison results under explicit cover conditions.
9. Return to Ascent once “basin,” “valley,” and the comparison target are formal objects.

So the immediate move should be a definition-and-contract pass, not a proof-discharge pass. No files were changed during this review.

---

# ive actually been meaning to update the organization and taxonomy a bit. i now have 'protolemmata', 'enthymemata', 'lemmas' and 'archeion' , maintaining the Greek register and the contemporary landing zone of 'lemmas' for verified deliverables, and archeion for the greek register archive where stale protolemmata and other cleanup items can go

That taxonomy is much cleaner, and it sharpens the lifecycle:

```text
protolemmata → enthymemata → lemmas
      ↘             ↘
             archeion
```

`archeion` should be a side exit, not a fourth maturity stage: material goes there because it is superseded or no longer active, not necessarily because it was wrong.

I’d give the directories these contracts:

- `Protolemmata/`: conjectures, design arguments, candidate statements, counterexamples, and source notes. Expected to change freely.
- `Enthymemata/`: compiling Lean statements whose definitions and hypotheses are accepted, but whose proofs or explicitly tracked assumptions remain incomplete.
- `Lemmas/`: verified deliverables—no `sorry`, and any trusted external theorem exposed clearly rather than hidden behind an opaque declaration.
- `Archeion/`: superseded drafts, abandoned proof shapes, transcripts, and cleanup artifacts. Excluded from default builds, but indexed with an archive reason and successor link.

Applied to the material we just reviewed:

- `Fable5-BARS-revisions.md` remains a protolemma while its architecture is being corrected.
- `BARS.v2.lean` becomes the active `Enthymemata/BARS.lean` only after its observable definitions are settled.
- The old `BARS.lean` goes to `Archeion/` once superseded.
- `PKWangA.lean`, being proved, belongs in `Lemmas/`.
- `PKWangB`, `Stability`, and genuinely stated-but-undischarged results belong in `Enthymemata/`.
- Header-only Potts, bifiltration, and Ascent shells are still protolemmata in substance; they should not acquire enthymema status until they contain accepted Lean statements.
- The prospective component-diagram/merge-forest architecture should begin in `Protolemmata/`, rather than being prematurely installed as an enthymema.

For filesystem and module names, I recommend consistent Lean-style capitalization everywhere:

```text
lean/
  Protolemmata/
  Enthymemata/
  Lemmas/
  Archeion/
```

with imports such as `Enthymemata.BARS`. That also eliminates the current Windows-tolerated, Linux-sensitive `enthymemes` versus `Enthymemes` mismatch.

Each archived item only needs a small provenance header:

```yaml
status: archived
archived-on: 2026-08-08
reason: superseded
successor: ../Protolemmata/...
```

This taxonomy also suggests a better build policy: build `Lemmas` and `Enthymemata`, scan the latter for declared obligations, and exclude `Protolemmata` and `Archeion` entirely. It gives your Greek register real semantic force rather than making it merely decorative.

---

# ok we can rename with capitalization. I would also add a Theorems folder but i dont have any of those at the moment, at least nothing of consequence meriting that label

---

# why not fix the stale tool path issue and other staleness problems, because they aren't going to un-stale themselves

---

# which portable powershell profile are you talking about, the pwsh profile.ps1?

---

# all of the portable environment bindings are in the user scope, not machine

---

# fyi, all i do is set the persistent values in user scope nowadays, i don't use app launchers. its just user scope + shell instance specific profile.ps1 for console history basically

---

# oh i guess i haven't set anything up for lean yet . i have a placeholder `D:\aghado01\command-center\portable-env\env-user\set-lean-user.ps1` where i was planning to set my ambient lean environment so that the lean harness can be launched and integrated into project similarly to how R is handled

Done. Lean now follows your actual environment model: persistent user scope only; no profile, launcher, or machine-scope changes.

- [set-lean-user.ps1](/D:/aghado01/command-center/portable-env/env-user/set-lean-user.ps1) now prepends `PDenv\elan\bin` to the user PATH and persists `ELAN_HOME`.
- The real user registry was verified. A PATH backup was written to `C:\Users\azrie\AppData\Local\portable-env\path-backups\PATH-User-20260808-165158.txt`.
- [meta-ci.ps1](/D:/aghado01/ThermoMapper/lean/scripts/meta-ci.ps1) now tolerates stale inherited PATHs by resolving Lake through `PORTABLE_ROOT`; it passed without manual PATH injection.
- The taxonomy is now `Protolemmata → Enthymemata → Lemmas`, with `Theorems` and the `Archeion` side exit, all capitalized and documented in [lean/README.md](/D:/aghado01/ThermoMapper/lean/README.md).
- The proved PKWang core was promoted to [Lemmas/PKWang.lean](/D:/aghado01/ThermoMapper/lean/Lemmas/PKWang.lean).
- Stability’s abstract constructors were repaired and the module can no longer evade the aggregate build.
- The noncompiling BARS revision was preserved honestly as [Protolemmata/BARS.v2.lean.md](/D:/aghado01/ThermoMapper/lean/Protolemmata/BARS.v2.lean.md), with its current failures recorded.
- Full validation passed: 8,487 jobs under Lean 4.30.0. Only the pre-existing deliberate obligations warn.
- The accidental Lean 4.32.2 installation was removed; only 4.30.0 and 4.31.0-rc1 remain.
- No stale active taxonomy paths remain. No files were staged or committed, and unrelated user changes were preserved.

The substantive Fable 5 revisions—tameness/windowed stability, generalized PKWang thresholding, and the BARS v2 statement changes—remain deliberately protolemmatic; I did not disguise semantic revision as cleanup.

---

# i dont mind upgrade to 4.32.2 and disposing 4.31.0

---

# i had to reboot application; i've set the ELAN env variable user scope and added to path. project-scoped remains project-scoped

---

# sandbox unlocked

Upgrade complete.

- Project pinned to Lean/Mathlib 4.32.2 in [lean-toolchain](/D:/aghado01/ThermoMapper/lean/lean-toolchain) and [lakefile.toml](/D:/aghado01/ThermoMapper/lean/lakefile.toml).
- Mathlib manifest refreshed to commit `905b95818eb32af7874a58b427f50c1711a5e96c`.
- Full meta-CI passed: 8,667 jobs.
- Updated PKWang for Lean 4.32’s `Std.Symm` structure.
- Removed Lean 4.31.0-rc1. Installed toolchains are now 4.30.0 and 4.32.2.
- Confirmed only `ELAN_HOME` is ambient; all `LEAN_*` compiler bindings remain unset/project-scoped.
- Retested [BARS.v2.lean.md](/D:/aghado01/ThermoMapper/lean/Protolemmata/BARS.v2.lean.md) under 4.32.2; its same four known failures remain documented.
- Updated [README.md](/D:/aghado01/ThermoMapper/lean/README.md) with the new versions and cold-build costs.

No changes were committed or staged. A future worthwhile optimization is replacing the broad `import Mathlib` directives—the cold elaborations took 5–12 minutes.

---

# ok lets update the broad imports and document scoping imports as a best practice, updating them as needed when new dependencies arise

Done. Every `import Mathlib` outside `.lake` has been replaced with scoped imports.

- Active modules compile successfully with explicit dependencies.
- Protolemma examples were updated too.
- [README.md](D:/aghado01/ThermoMapper/lean/README.md) and [AGENTS.md](D:/aghado01/ThermoMapper/AGENTS.md) now document the maintenance practice.
- [meta-ci.ps1](D:/aghado01/ThermoMapper/lean/scripts/meta-ci.ps1) rejects future umbrella imports in active source.
- Full meta-CI passes: 2,512 jobs, down from 8,667 with the umbrella imports.
- `BARS.v2` retains exactly its four documented proof/name gaps; no dependency errors remain.
- `git diff --check` passes. No changes were staged or committed.

---

# ok can we clear dead weight that's not imported or a dependency without breaking things. i'd rather re-cache selectively later with the latest than just have the bulk sitting dormant

Cleanup complete without breaking the Lean harness.

- Footprint reduced from **7.43 GiB to 2.00 GiB**, reclaiming **5.43 GiB**.
- Mathlib artifacts reduced from 8,333 modules to 2,227; auxiliary dependencies retain only the active 2,473-module closure.
- Removed the redundant 822 MiB compressed user cache at `C:\Users\azrie\.cache\mathlib`.
- Package sources, manifest, and toolchain remain intact.
- All 2,512 build jobs are up to date and meta-CI passes.
- [README.md](D:/aghado01/ThermoMapper/lean/README.md) now documents module-specific recaching and warns that bare `lake exe cache get` restores all of Mathlib.

Everything removed is regenerable through scoped `lake exe cache get Mathlib.…` calls after future dependency or version updates. Nothing was staged or committed.

---

# so circling back on the actual lean work we were trying to do before this infra digression, what's the plan?

The short version: repair BARS semantically, extract its proof-closed core, then build the shared `π₀` foundation before expanding into Mapper. The current [Enthymemata/BARS.lean](D:/aghado01/ThermoMapper/lean/Enthymemata/BARS.lean) is stale and should not be proved as written.

1. Reconcile the BARS statement surface

Use [BARS.v2.lean.md](D:/aghado01/ThermoMapper/lean/Protolemmata/BARS.v2.lean.md) as revision source and incorporate the settled corrections:

- Make `local_max_is_critical` interval-free.
- Add the antitone crossing lemma.
- Treat `levelWidth` honestly as global occupancy.
- Add the bracketed per-peak width result.
- Split elementary non-affinity from genuine draw-relative FWHM noncommutation.
- Add the overcount witness alongside the undercount witness.
- Preserve the documented plateau, boundary, and `ncard` contracts.

This produces a compiling enthymema surface before we invest in proofs.

2. Discharge the BARS core

The immediate proof tranche is:

- `local_max_is_critical`
- `finite_local_max`
- `level_set_finite`
- both unique-crossing lemmas
- `superlevel_of_peak_bracket`
- `peak_width_eq_bracket_roots`
- `levelWidth_not_affine`

The obsolete root theorem is now known exactly: `Polynomial.finite_setOf_isRoot`. The two measure failures need explicit set rewriting rather than more imports.

Expected result: only three genuinely substantial obligations remain:

- `count_of_mean_undercounts`
- `count_of_mean_overcounts`
- draw-relative `fwhm_of_mean_ne_mean_of_fwhm`

3. Promote the stable core selectively

Move the reviewed, proof-closed definitions and results into `Lemmas`, while leaving the three witness constructions in `Enthymemata`. We should not wait for the entire BARS file to close before promoting reusable polynomial and bracket machinery.

4. Design the `π₀` foundation

Then create a protolemma/design pass for the connected-components construction over graph filtrations. This becomes the shared formal vocabulary for:

- PKWang threshold chains and dendrograms
- merge trees
- the H₀ bifiltration
- Mapper pullback components
- eventually Ascent basin comparisons

We should first decide whether the Lean surface needs a literal category-theoretic functor or a lighter order-indexed family of component maps.

5. Resume the cross-estate build order

After `π₀`:

- Replace PKWangB’s sorted cumulative sums with deterministic filtered sums, matching the implementation’s diagonal/tie convention.
- Generalize the proved PKWang threshold reduction from `½` to arbitrary `θ ∈ (0,1)`.
- Recast `PottsGriffiths` as the FK-based monotonicity chain.
- Retrofit Stability with tameness and a windowed sup norm; do not pretend its opaque landscape obligation is currently provable.
- Build Bifiltration from FK monotonicity plus `π₀`.
- Define Mapper from cover pullbacks, induced graphs, components, and their nerve.
- Leave Ascent until that shared component vocabulary exists.

The large extensions from the discussion—total `(m,n,k)` peak-count independence, raw-count instability, and the carrier-generic merge tree—remain `Protolemmata` until `π₀` and the gluing semantics are settled. No paper-specific concepts enter the ThermoMapper Lean modules, and nothing currently merits `Theorems`.

My recommended next action is phases 1–2: synchronize BARS and close its reusable core before opening the `π₀` design.

---

# well i think bars_v2 is the new one, bars v1 only needs to be archived i think

Exactly. I mischaracterized it as a reconciliation; this is a replacement.

Both [Enthymemata/BARS.lean](D:/aghado01/ThermoMapper/lean/Enthymemata/BARS.lean) and [Protolemmata/BARS.lean.md](D:/aghado01/ThermoMapper/lean/Protolemmata/BARS.lean.md) are semantically v1—the differences are formatting only. The correct sequence is:

1. Fix the four elaboration gaps in [BARS.v2.lean.md](D:/aghado01/ThermoMapper/lean/Protolemmata/BARS.v2.lean.md) while leaving its genuine obligations as `sorry`.
2. Confirm v2 compiles as an enthymema.
3. Atomically archive both v1 representations under `Archeion/`.
4. Move v2 to `Enthymemata/BARS.lean`; no merging or backporting from v1.
5. Discharge and selectively promote v2’s reusable core.
6. Then proceed to the `π₀` foundation and the broader estate sequence.

That avoids ever leaving the aggregate broken and makes v2 the sole active BARS authority.

---

# i think it also needs a slightly more descriptive name. its not just "BARS at-large" , its about an application of BARS in SPC. but not sure about what to call it

I’d name the new module **`Spc.BarsReadout`**:

- File: `lean/Enthymemata/Spc/BarsReadout.lean`
- Namespace: `Spc.BarsReadout`
- Document title: **“SPC × BARS: posterior thermal-feature readout”**

“Readout” fits because these results do not establish BARS itself; they establish what SPC may validly derive from BARS posterior curves—peaks, counts, brackets, widths, occupancy, and draw-wise summaries.

I’d avoid:

- `Spc.Bars`: still too broad.
- `Spc.BarsPeaks`: too narrow for level sets, widths, and future merge-tree results.
- `Spc.BarsFeatureExtraction`: sounds like implementation machinery.
- `Spc.BarsPosteriorReadout`: precise, but unnecessarily restrictive at the module level; the title and module documentation can carry “posterior.”

The generic polynomial facts can initially live inside this module because that is their present proof role. If they acquire independent consumers later, we can extract them into a neutral analysis lemma module without changing the SPC-facing theory.

The old v1 artifacts can retain their historical `BARS` name when moved into `Archeion`; only v2 becomes `Spc.BarsReadout`.

---

# is there anything that might have been left behind in the v1 BARS.lean of importance?

Literal answer: no important Lean declaration or proof is stranded in v1. All six v1 theorems are still `sorry`, and v2 preserves or strengthens every statement.

- `local_max_is_critical` is strictly stronger in v2 because the unnecessary interval hypothesis is removed.
- Five declarations are unchanged.
- The old fixed-level “FWHM” statement survives exactly as `levelWidth_not_affine`.
- The name `fwhm_of_mean_ne_mean_of_fwhm` is reassigned to the genuine draw-relative result. That rename is the only migration hazard.
- No other Lean module uses these declarations; [Enthymemata.lean](/D:/aghado01/ThermoMapper/lean/Enthymemata.lean:2) only imports the module.

The audit did uncover three important obligations that were already absent from v1 rather than lost by v2:

- Recover `argmax_in_closed_form_set` from the referenced [single-peak protolemma](/D:/aghado01/ThermoMapper/lean/Protolemmata/fable-BARS-lemma.md:33) as a useful boundary-aware corollary.
- Eventually formalize spline globalization and the critical-neighbour test that connects the polynomial facts to `SplineExtrema`.
- Model the actual SPC readout: prominence-gated peaks, configurable boundary treatment, clipping, and prominence-relative per-branch spans. Current v2’s raw `peakCount` and global `halfMax` are scaffolding, not yet an exact model of [SplineExtrema.cs](/D:/aghado01/ThermoMapper/src/maths/regression/spline/bars/SplineExtrema.cs:37).

So I would archive v1 without merging it into v2, while adding those three items to the new `Spc.BarsReadout` ledger. The old cubic-only warning should not be carried forward—the implementation now handles higher-degree spans.

---

# ok lets capture these new obligations either in lean comments, or new protolemmata(s)

Captured in both places:

- Added [spc-bars-readout-obligations.md](/D:/aghado01/ThermoMapper/lean/Protolemmata/spc-bars-readout-obligations.md), covering:
  - boundary-aware closed-form argmax;
  - spline globalization and neighbour-test faithfulness;
  - the actual SPC prominence, boundary, clipping, and branch-relative span contract;
  - intended `Pi0` ownership and discharge order.
- Updated [BARS.v2.lean.md](/D:/aghado01/ThermoMapper/lean/Protolemmata/BARS.v2.lean.md:20) so its comments no longer conflate global half-maximum occupancy with the implemented prominence-relative span.

I also recorded that the older `argmax_expectation_noncommute` sketch is defective as written and that the cubic-only guard is now stale.

All relative links resolve, code fences balance, and the diff passes whitespace checks. No active `.lean` source changed, so a Lean rebuild was unnecessary.

---

# ok lets review the pkwang lemmas. i know we had discussed earlier some things, i believe one item was expanding theta coverage beyond 1/2 and implications of that

Yes—the earlier θ-generalization is valid, and it reveals more structure than merely supporting another parameter value.

For

\[
\kappa(\theta)=-\log(1-\theta),
\]

the general reduction is

\[
1-\exp(-H/T)>\theta
\quad\Longleftrightarrow\quad
H>T\,\kappa(\theta).
\]

I verified this statement against Lean 4.32.2. Algebraically it needs only `T > 0` and `θ < 1`; application semantics add `0 ≤ θ`. At `θ = 1/2`, `κ(θ)=log 2`, recovering the existing theorem in [PKWang.lean](/D:/aghado01/ThermoMapper/lean/Lemmas/PKWang.lean:25).

### The main implication

For MeanField, LocalField-Mutual, and LocalField-Inclusive, the graph depends only on the scalar cut coordinate

\[
\tau=T\,\kappa(\theta).
\]

Thus:

- Higher `T` or higher `θ` removes edges.
- Any two `(T,θ)` pairs with equal `Tκ(θ)` produce exactly the same graph and connected components.
- `θ=1/2` is simply the exponential median convention—not a special physical constant.
- The apparent solver bifiltration `(T,θ)` factors through a one-dimensional chain. It is a degenerate bifiltration or gauge-redundant presentation, not two independent structural axes.

The earlier discussion was therefore directionally right but overstated the independence of the axes. In ordinary order the graph is antitone in both `T` and `θ`; a Lean filtration should use the corresponding order duals.

There is one important exception: LocalField with arithmetic-mean symmetrization averages probabilities after the nonlinear exponential transform. It cannot generally be reduced to a temperature-independent effective `Hcum`. Edge orderings can change with `θ`, making that mode genuinely two-parameter. Mutual and Inclusive remain clean because `min` and `max` commute with the monotone exponential kernel.

### Current lemma audit

| Item | Verdict | Revision |
| --- | --- | --- |
| `pk_wang_closed_form_reduction` | Correct and proved | Generalize to `θ`; retain the half-threshold result only as the median corollary |
| `deterministicCutGraph` | Correct for symmetric effective energies | Factor into a generic scalar cut graph plus `Tκ(θ)` specialization |
| `localHCum` / `globalHCum` | Honest obligations, but currently underspecified | Define filtered sums over an explicit base graph, not every vertex pair |
| Lemma A pooling | “Limit” is too vague | Replace with an exact double-counting identity and a homogeneity corollary |
| Single-linkage result | Right destination, weak current formulation | State the threshold chain through `Pi0`; prove preservation of the coupling order |

Lemma A has a particularly clean exact formulation. If

\[
L_u(t)=\sum_{\substack{e\ni u\\J_e\le t}}J_e,
\qquad
G(t)=\sum_{\substack{e\in E\\J_e\le t}}J_e,
\]

then on a finite undirected graph,

\[
\sum_u L_u(t)=2G(t).
\]

Under homogeneous local cumulative fields, `G` is therefore a fixed positive multiple of each `Lᵤ`. Local and global hierarchies agree after a temperature rescaling—not necessarily at the same numerical temperature. That is a better Mapper-backend control theorem than an unspecified asymptotic limit.

### Live implementation consequences

The Lean review exposed three current code mismatches:

- [PKWang.EstimateBracket](/D:/aghado01/ThermoMapper/src/clustering/graphical/spc/runtime/core/solver/PKWang.cs:113) always divides by `log 2`, even though clustering already accepts arbitrary `theta`. For general θ it must use `κ(θ)`; the endpoints `θ=0` and `θ=1` have degenerate sweeps and need separate handling.
- [MeanField.cs](/D:/aghado01/ThermoMapper/src/clustering/graphical/spc/runtime/core/solver/MeanField.cs:57) and LocalField assign partial cumulative sums within tied coupling classes. The filtered-sum specification assigns the entire tie class one value. Current tests deliberately use distinct couplings, so this discrepancy is untested.
- Partitioning uses strict `G > θ`, while [ThermalDendrogram.cs](/D:/aghado01/ThermoMapper/src/clustering/dendrogram/ThermalDendrogram.cs:76) uses `G ≥ θ`. Their components can disagree exactly at a crossing.

I would revise PKWang in this order:

1. General θ reduction, `κ`, cut-coordinate factorization, and monotonicity—all proof-closed.
2. Correct graph/tie semantics and θ-aware bracketing.
3. Filtered local/global cumulative fields and the exact pooling identity.
4. Mutual/Inclusive commutation theorems plus an explicit Mean counterexample.
5. `Pi0` threshold-chain/single-linkage equivalence.

The genuinely consequential deliverable is not the elementary logarithm identity by itself; it is the theorem that the default PKWang hierarchy is θ-gauge-invariant up to temperature rescaling.

---

# well, the live implementations are certainly subject to change, with the original global hcum formulation likely being relegated to a didactic case study on the error of the pkwang paper, since local field case in some form is more density robust, but the question is how interesting is the choice of symmetrization operation (mean vs anything else), if so then that means there are two parameters because of your point about arithmetic mean and otherwise is there a satisfying alternative op that doesn't induce the extra free parameter, is there any practical interest in parametrizing inclusive vs mutual; bracketing is something i want to expose as a hyperparameter of SPC since outside of the nominal gaussian coupling and for more general graph compiler semantics and configurations, this case this becomes important

Yes—the symmetrization choice is scientifically meaningful, but ordinary arithmetic averaging of the affinities should not be the unnamed default.

Let \(p_T(h)=1-e^{-h/T}\), with directed local values \(a=H_{u\to v}\) and \(b=H_{v\to u}\). The clean choices are:

| Policy | Symmetric field \(H_{\mathrm{eff}}\) | Equivalent affinity | Meaning |
|---|---:|---:|---|
| Mutual | \(\min(a,b)\) | \(\min(p_T(a),p_T(b))\) | both endpoints must endorse |
| Balanced | \((a+b)/2\) | \(1-\sqrt{(1-p_T(a))(1-p_T(b))}\) | average cumulative hazard |
| Inclusive | \(\max(a,b)\) | \(\max(p_T(a),p_T(b))\) | either endpoint may endorse |

The balanced expression is, I think, the satisfying parameter-free alternative you were looking for. It averages in `Hcum`/cumulative-hazard space and therefore geometrically averages the two survival probabilities. It is symmetric, monotone, idempotent, lies between Mutual and Inclusive, and does not invent another structural axis.

For all three,

\[
p_T(H_{\mathrm{eff}})>\theta
\iff
H_{\mathrm{eff}}>T\,\kappa(\theta),
\qquad
\kappa(\theta)=-\log(1-\theta).
\]

So \(T\) and \(\theta\) reduce to the single effective cut

\[
\tau=T\kappa(\theta).
\]

That revises my earlier bifiltration claim: this is formally a two-indexed monotone family, but it is a degenerate one that factors through \(\tau\), not a genuinely two-dimensional solver family.

The current post-link arithmetic mean is different:

\[
\frac{p_T(a)+p_T(b)}2
=
1-\frac{e^{-a/T}+e^{-b/T}}2.
\]

Converted back into an energy, it becomes

\[
H_{\mathrm{soft}}(T)
=
-T\log\!\left(\frac{e^{-a/T}+e^{-b/T}}2\right).
\]

That is a legitimate temperature-dependent soft minimum: it approaches \(\min(a,b)\) in the cold limit and \((a+b)/2\) in the hot limit. Consequently, edge orderings can change with temperature and \(T,\theta\) become genuinely independent parameters.

That is potentially interesting—it amounts to saying reciprocity should be strict when cold and soften as temperature rises—but it is a specific model, not a neutral “mean.” I would retain it only under an explicit name such as `AffinityMean` or `FreeEnergyMean`; the unqualified balanced rule should be `HazardMean` or `FieldMean`.

On Mutual versus Inclusive: yes, there is practical interest, but I would expose it categorically rather than add a continuous mixing parameter. Mutual and union/symmetric neighborhood graphs are established distinct constructions with different cluster-isolation behavior, not cosmetic alternatives; the mutual construction is especially protective against another region linking inward unilaterally. [Maier, Hein, and von Luxburg](https://arxiv.org/abs/0912.3408) analyze precisely this distinction.

For Mapper I would use:

- `Mutual` as the conservative default, especially since false bridges can merge pullback components.
- `HazardMean` as the balanced research alternative.
- `Inclusive` as an explicit recall/periphery or sensitivity mode.

They also produce a useful exact sandwich:

\[
E_{\mathrm{Mutual}}(\tau)
\subseteq
E_{\mathrm{HazardMean}}(\tau)
\subseteq
E_{\mathrm{Inclusive}}(\tau).
\]

Rather than tuning a continuous “mutualness” coefficient, we can run the two envelopes. Edges surviving Mutual are bilateral/core evidence; edges appearing only under Inclusive are one-sided/asymmetry evidence. That difference is itself a density-boundary diagnostic.

There is one prior question to settle, however: the current local fields are raw cumulative sums. Their magnitudes depend on local degree and total coupling mass, so averaging them assumes the endpoint scales are commensurate. “Local field in some form” is exactly the right qualification. I would separate:

1. local calibration—raw cumulative energy, normalized cumulative mass, or rank/quantile;
2. endpoint symmetrization;
3. exponential thermal link.

Local scaling is well motivated for multi-scale and cluttered data, but that literature does not automatically validate this particular raw cumulative-sum calibration. [Zelnik-Manor and Perona](https://papers.nips.cc/paper_files/paper/2004/hash/40173ea48d9567f1f393b20c855bb40b-Abstract.html) supports the local-scale principle, not this exact estimator.

On bracketing: agreed completely. It belongs to the SPC run configuration, with the graph artifact declaring its coupling gauge and calibration evidence. That matches the existing separation in the [bracketing discussion](D:/aghado01/ThermoMapper/issues/spc-BARS/discussions/sol-SPC-bars-graphs-brackets.md:443).

For PKWang specifically, I would make the native bracket live in \(\tau\)-space:

- derive \(H_{\mathrm{eff}}\) after local calibration and symmetrization;
- bracket its critical values or robust quantiles;
- map to \(T\) using \(T=\tau/\kappa(\theta)\);
- record the requested policy, resolved bracket, coordinate, and calibration provenance.

The implementation currently does the opposite: it symmetrizes after applying the exponential in [PKWang.cs](D:/aghado01/ThermoMapper/src/clustering/graphical/spc/runtime/core/solver/PKWang.cs:127), while [EstimateBracket](D:/aghado01/ThermoMapper/src/clustering/graphical/spc/runtime/core/solver/PKWang.cs:82) accepts a symmetrization rule but brackets the unsymmetrized directed values and hardcodes \(\log 2\). Thus the rule currently has no effect on bracket estimation.

For deterministic PKWang, we can go further: construct the complete critical edge ladder once. Then the bracket is only the acquisition/readout window, not something that changes the underlying hierarchy. For stochastic SPC, it remains a genuine experimental hyperparameter controlling where simulation and BARS acquisition occur.

The Lean revision I would make from this is:

- generalize the threshold theorem to arbitrary \(\theta\);
- prove the \((T,\theta)\mapsto\tau\) factorization;
- prove link-commutation for Mutual, HazardMean, and Inclusive;
- prove the graph/partition sandwich;
- record post-link `AffinityMean` separately as a genuinely temperature-dependent soft-min construction.

That gives the global Wang formulation a clean didactic role, while the live local-field theory makes the calibration and symmetrization choices explicit rather than accidentally burying them in `Mean`.

---

# lets review the paper D:\aghado01\graveyard\codex-scientiae\bibliotecha\compendia\clustering\PKW2020.md in contetx here, use mdnav maybe

---

# note: the ln2 is likely hard-coded in the paper and in reading it perhaps you can verify with a fresh look that their MCMC iterations are in fact spurious

Yes. A fresh reading strengthens the conclusion considerably:

> The paper’s replacement \(M\)-draw procedure is computationally spurious even on its own terms—and, independently, it does not sample the Potts model it claims to approximate.

I read the full methodology and results through `mdnav`—98.5% of the document, with only the title and contents unread.

## The fundamental derivation is invalid

The paper begins correctly with a Boltzmann probability over spin configurations:

\[
P_T(S)=Z^{-1}e^{-H(S)/T}.
\]

It then treats this as an exponential probability density over scalar energies and samples

\[
E=-T\log(1-r).
\]

That transition occurs between equations 3 and 6 in [PKW2020.md](/D:/aghado01/graveyard/codex-scientiae/bibliotecha/compendia/clustering/PKW2020.md:40).

But the energy distribution induced by the Potts measure is instead

\[
P_T(H=h)=\frac{\Omega_q(h)e^{-h/T}}{Z},
\]

where \(\Omega_q(h)\) counts the spin configurations having energy \(h\). The omitted density-of-states factor carries the graph structure and \(q\)-dependence.

There are several immediate contradictions:

- A finite Potts Hamiltonian has bounded, discrete support; their exponential draw is continuous and unbounded.
- The Potts energy law depends on \(q\); equation 6 does not.
- They report \(q=20\) in the experiment, but \(q\) has no effect anywhere after the replacement begins.
- Figure 1a merely confirms that inverse-transform samples from an exponential resemble that same exponential. It does not compare them with energies sampled from the Potts system.

The easiest formal refutation is therefore not density-of-states machinery at all:

\[
P_{\mathrm{Potts}}(H>H_{\max})=0,
\qquad
P_{\mathrm{Exp}(T)}(E>H_{\max})=e^{-H_{\max}/T}>0.
\]

So the two laws cannot be equal for any finite nontrivial graph.

## The constructed objects are not spin states

After drawing a scalar \(E\), the paper globally sorts the \(N K\) couplings,

\[
J_1\le\cdots\le J_L,\qquad
c_n=\sum_{m\le n}J_m,
\]

matches \(E\) to the closest \(c_n\), and declares the included edges to have \(\delta_{ij}=0\) and the remaining edges \(\delta_{ij}=1\). This is stated explicitly in [the global-prefix passage](/D:/aghado01/graveyard/codex-scientiae/bibliotecha/compendia/clustering/PKW2020.md:64).

But \(\delta_{ij}=1\) means “\(s_i=s_j\),” which must be an equivalence relation. An arbitrary prefix mask need not be transitive. For example, it can declare

\[
\delta_{12}=1,\qquad
\delta_{23}=1,\qquad
\delta_{13}=0,
\]

which cannot arise from any assignment of spins because \(s_1=s_2=s_3\) would force \(\delta_{13}=1\).

Thus the paper never reconstructs \(S_i\), despite repeatedly calling these masks “spin states.” It constructs thresholded edge masks and later takes their graph closure. That is a graph-clustering procedure, not Potts sampling.

## The \(M\) draws are even more spurious than we had stated

The replacement draws are not MCMC—the paper itself says the new method is non-Markov. They are iid exponential energy budgets. Original SW iterations are meaningful for sampling the actual Potts target; the replacement’s iid draws are unnecessary for its surrogate.

Under the paper’s literal “closest cumulative energy” rule, define

\[
b_1=0,\qquad
b_r=\frac{c_{r-1}+c_r}{2}\quad(r\ge2).
\]

Ignoring probability-zero matching ties, an edge of rank \(r\) is marked same-spin precisely when \(E<b_r\). Therefore its finite-\(M\) reported correlation is merely the empirical CDF:

\[
\widehat G_r
=
\frac1M\sum_{m=1}^M\mathbf 1[E_m<b_r]
=
\widehat F_M(b_r).
\]

Thresholding at arbitrary \(\theta\) gives

\[
\widehat G_r>\theta
\iff
b_r>E_{(\lfloor M\theta\rfloor+1)},
\]

where \(E_{(k)}\) is one order statistic of the samples.

For the paper’s \(M=300\) and \(\theta=0.5\), the entire final clustering depends only on the 151st sorted energy draw. Their proposed \(M\)-way hardware computes \(M\) complete \(\delta\)-matrices, averages them, and thresholds the average—all to recover one random scalar cut.

In the population limit,

\[
G_r=1-e^{-b_r/T},
\]

and hence

\[
G_r>\theta
\iff
b_r>T[-\log(1-\theta)].
\]

At the paper’s hard-coded \(0.5\),

\[
b_r>T\log 2.
\]

So your intuition is right with a refinement: the paper hard-codes \(0.5\), while \(\log 2\) is its latent population-limit cut. At finite \(M=300\), the literal cutoff is the random 151st order statistic, not exactly \(T\log2\).

## Our current closed form is cleaner—but not literal paper parity

ThermoMapper currently uses

\[
E<H_{\mathrm{cum},r}
\]

directly in [WangMonteCarloDemo.cs](/D:/aghado01/ThermoMapper/src/clustering/graphical/spc/demos/WangMonteCarloDemo.cs:39), producing

\[
G_r=1-e^{-c_r/T}.
\]

That corresponds to a floor/prefix-threshold reinterpretation, not the paper’s twice-stated “closest” rule. Literal closest matching produces midpoint \(b_r\), not cumulative endpoint \(c_r\).

Consequently:

- The proved theorem in [PKWang.lean](/D:/aghado01/ThermoMapper/lean/Lemmas/PKWang.lean:25) is mathematically correct.
- It is correct for ThermoMapper’s cleaned exponential-CDF surrogate.
- It is not yet an exact theorem about the published algorithm.
- `ClosedForm_MatchesMonteCarlo` validates our chosen event against its CDF, but does not establish fidelity to the paper’s closest-prefix matching.
- The phrase “Wang 2020 as published” in the demo should be corrected.
- Its reference to nonexistent `issues/lean/{spc-lemmas, wang_2020_scandal}` is also stale.

## What the global formulation actually is

After correcting ties so equal couplings receive the same cumulative value, global \(H_{\mathrm{cum}}\) is a monotone transform of \(J\). Therefore every threshold of \(H_{\mathrm{cum}}\) is simply a threshold of the original couplings.

Its entire hierarchy is ordinary single linkage on the weighted graph, with temperature merely reparameterizing the cut height. The Monte Carlo, \(q\), spins, Potts Hamiltonian, and proposed hardware are all vestigial.

The paper’s empirical results can still show that a thresholded Gaussian-weighted \(k\)-NN graph works on concentric rings and Iris. They do not show equivalence to SPC. The speedup compares a much simpler graph-filtration algorithm with a sampler targeting a different probability law.

## Consequences for the live local-field work

The local field is absent from the paper. It is our exaptation of the only potentially useful residue: cumulative coupling evidence.

That means it should be presented honestly as a new deterministic method:

\[
\text{local calibration}
\rightarrow
\text{endpoint symmetrization}
\rightarrow
\text{thermal monotone link}
\rightarrow
\pi_0.
\]

It is not a repaired Potts sampler or a parallel SPC approximation. Its scientific merits must come from density adaptation, graph-filtration behavior, Mapper suitability, and empirical comparison with robust linkage methods.

This strengthens the earlier symmetrization conclusion:

- Paper-global has no principled symmetrization theory.
- Local fields introduce two endpoint-relative values, so Mutual/HazardMean/Inclusive are genuinely our semantics.
- Symmetrizing calibrated `Hcum` before the exponential preserves the one-dimensional cut \(\tau=T[-\log(1-\theta)]\).
- Raw local cumulative sums still need a calibration decision; locality alone does not guarantee density invariance.

The paper also supplies no principled bracket. It fixes the Gaussian scale \(a\), sweeps temperature, and retrospectively labels Iris \(T\in[1.2,1.4]\) as the SPC phase. Thus our configurable, provenance-bearing bracket design is not deviating from a paper prescription; no general prescription exists there.

## Recommended Lean disposition

I would revise the PKWang line into five precise results:

1. **Generic exponential-CDF cut**  
   The existing theorem, generalized to arbitrary \(\theta\). This is valid mathematics, independent of the paper.

2. **Finite-\(M\) order-statistic collapse**  
   The thresholded empirical CDF depends on one order statistic; the \(M\) masks are unnecessary even before taking \(M\to\infty\).

3. **Literal closest-prefix reduction**  
   Paper fidelity uses midpoint boundaries \(b_r\), distinguishing it from ThermoMapper’s cleaned \(c_r\) formulation.

4. **Prefix mask is not necessarily spin-realizable**  
   A three-vertex counterexample disproves the claim that the constructed \(\delta\)-matrix is a spin state.

5. **Exponential energy law is not the Potts energy law**  
   The bounded-support counterexample gives a particularly clean formal refutation.

Then the tie-corrected global-single-linkage theorem remains as the didactic conclusion, while the local-field definitions and symmetrization lemmas move under their own ThermoMapper method identity.

So I would indeed relegate the global formulation to a case study—but preserve the negative theorems in `Lemmas`, not bury them in `Archeion`. What becomes archival is the superseded proto-exposition and paper-faithful implementation; the verified diagnosis is consequential deliverable-grade mathematics.

---

# were there simulations being used for clustering or was it just a numerical verification

They were used operationally for clustering, not merely as a numerical verification—but they were not genuine spin simulations.

For each temperature, the paper:

1. Draws \(M\) iid exponential energies; the experiment uses \(M=300\).
2. Matches every draw to the closest cumulative-coupling prefix.
3. Produces \(M\) pairwise \(\delta_{ij}\) masks.
4. Averages them into \(\widehat G_{ij}\).
5. Connects pairs with \(\widehat G_{ij}>0.5\).

That procedure is described in [PKW2020.md](/D:/aghado01/graveyard/codex-scientiae/bibliotecha/compendia/clustering/PKW2020.md:64), and the results say the generated samples were used for the reported classifications at [line 70](/D:/aghado01/graveyard/codex-scientiae/bibliotecha/compendia/clustering/PKW2020.md:70).

Figure 1a is only a numerical verification that inverse-transform exponential draws reproduce the exponential distribution they were drawn from. That check is tautological, but the draws themselves subsequently feed the clustering.

The sharper characterization is:

- Conventional SPC: genuine SW Markov-chain spin simulations.
- Paper’s replacement: iid scalar-energy Monte Carlo used to construct clustering masks.
- ThermoMapper closed form: analytic evaluation of what those iid draws estimate.

Moreover, at \(M=300\) and threshold \(0.5\), the final partition depends only on the 151st ordered energy draw. So the simulations were genuinely used by their implementation, but mathematically all 300 masks and their averaging were redundant. In the population limit, that random cut becomes \(T\ln 2\)—or the corresponding midpoint-ladder comparison under the paper’s literal “closest” rule.

---

# so i would have a two lean narratives on this. one is the critical examination of pkwang (lack of correspondence with spc, spurious simulations; effectively a correction against the paper's errors) + the recovered insight and approach that comes from it, tbd

Exactly. They should be two distinct Lean narratives sharing a small neutral mathematical kernel—not “Lemma A and Lemma B” of one purported PKWang method.

1. **Critical examination / correction of PKWang**

   Formalize the published algorithm faithfully, including its closest-prefix rule, then prove:

   - the exponential energy law cannot equal the finite Potts energy law;
   - prefix masks need not correspond to any spin configuration;
   - the replacement uses iid Monte Carlo, not MCMC;
   - finite-\(M\) thresholding collapses to one sample order statistic;
   - the population limit collapses to a deterministic \(\theta\)-quantile cut;
   - the global cumulative ladder is ultimately a reparameterized single-linkage hierarchy.

   This narrative establishes two independent failures: no correspondence with SPC, and unnecessary simulation even for the paper’s surrogate.

2. **Recovered construction / exaptation**

   Start without any claim of Potts equivalence:

   \[
   \text{directed local evidence}
   \rightarrow
   \text{calibration}
   \rightarrow
   \text{symmetrization}
   \rightarrow
   \text{cut filtration}
   \rightarrow
   \pi_0.
   \]

   Its open design questions are genuinely ours:

   - raw cumulative mass versus locally normalized mass/rank;
   - Mutual, HazardMean, or Inclusive symmetrization;
   - density-robustness conditions;
   - the effective cut coordinate \(\tau=T[-\log(1-\theta)]\);
   - bracket calibration and provenance;
   - Mapper pullback-component interpretation.

The shared kernel should be neutral—something like `ThermalCut` or `ExponentialCut`—containing the arbitrary-\(\theta\) CDF equivalence. The audit uses it to expose the collapse; the recovered method uses it constructively. That prevents the latter from inheriting the paper’s invalid statistical-mechanical justification.

A clean eventual organization would be:

- `Lemmas/PKWangAudit.lean` — verified negative/correction results.
- `Archeion/PKWangPublished...` — literal or superseded reproduction material.
- `Lemmas/ExponentialCut.lean` — neutral reusable mathematics.
- `Protolemmata/LocalCumulativeField.md` — recovered method while its semantics remain unsettled.
- Later `Enthymemata/LocalCumulativeField.lean` and ultimately `Lemmas/...`.

One wording discipline matters: call the paper’s replacement draws “spurious iid Monte Carlo,” not “spurious MCMC.” SW’s MCMC is legitimate for actual SPC; PKWang’s error is that its replacement does not sample that target, after which its own Monte Carlo is unnecessary.

The exact local/global pooling identity can become the bridge between the narratives, but not a vindication of the paper: it explains how the failed global pooling relates algebraically to the recovered local construction.

---

# i kind of want to break lemmas up into submodules kind of like my src, with a module for the PKWang-related ones, but renamed to something more fitting, and relegating references to pkwang in the PKWangAudit.lean

Yes—the clean dependency direction is:

\[
\texttt{PKWangAudit} \longrightarrow \texttt{CumulativeField},
\]

with no dependency or naming flowing back from the paper into the recovered mathematics.

`CumulativeField` is my preferred umbrella name. It names the actual mathematical object without prematurely claiming density robustness, SPC equivalence, or physical temperature semantics. The eventual clustering method could receive the more specific name `LocalCumulativeLinkage`.

A sensible structure would be:

```text
lean/
├── Lemmas/
│   ├── CumulativeField.lean
│   └── CumulativeField/
│       ├── Cut.lean
│       ├── CutGraph.lean
│       ├── Symmetrization.lean
│       ├── Pooling.lean
│       └── Linkage.lean
├── Enthymemata/
│   ├── CumulativeField.lean
│   ├── CumulativeField/
│   │   └── ...unfinished portions...
│   └── PKWangAudit.lean
└── Protolemmata/
    ├── CumulativeField/
    │   ├── design.md
    │   ├── calibration.md
    │   └── symmetrization.md
    └── PKWangAudit.md
```

The declarations would live under `namespace CumulativeField`, with neutral names such as:

- `cutScale`
- `exponentialAffinity_gt_iff`
- `cutGraph`
- `localCumulative`
- `globalCumulative`
- `poolingIdentity`
- `mutualField`
- `hazardMeanField`
- `inclusiveField`
- `cutGraph_mono`
- `componentSandwich`

`PKWangAudit.lean` would be the only active Lean source containing names such as `PKWang`, `publishedAlgorithm`, or `Wang2020`. It would import the generic module and prove paper-specific statements:

- published exponential samples are not the Potts energy law;
- published prefix masks need not be spin-realizable;
- finite-\(M\) thresholding equals an order-statistic cut;
- literal closest-prefix matching uses midpoint boundaries;
- the corrected population formulation reduces to `CumulativeField.cutGraph`;
- the global construction is a single-linkage reparameterization.

Because those audit proofs are currently obligations, the file should begin in `Enthymemata/PKWangAudit.lean` and move unchanged to `Lemmas/PKWangAudit.lean` once proof-closed.

The existing files then split cleanly:

- [PKWang.lean](/D:/aghado01/ThermoMapper/lean/Lemmas/PKWang.lean) becomes generic `CumulativeField/Cut.lean` and `CutGraph.lean`.
- [PKWangB.lean](/D:/aghado01/ThermoMapper/lean/Enthymemata/PKWangB.lean) becomes the unfinished `CumulativeField` definitions, pooling, and linkage modules.
- [spc-pkwang-lemmas.md](/D:/aghado01/ThermoMapper/lean/Protolemmata/spc-pkwang-lemmas.md) splits into the audit brief and the constructive `CumulativeField` design material.

One infrastructure detail: meta-CI already discovers nested files recursively, but it currently requires the root stage aggregate to import every leaf directly. For genuine hierarchical module barrels, it should instead verify recursive import closure. Then:

```lean
-- Lemmas.lean
import Lemmas.CumulativeField
import Lemmas.PKWangAudit
```

can remain concise while `Lemmas.CumulativeField` owns its internal imports.

I would also avoid `ThermalLinkage` as the umbrella: once \(T\) and \(\theta\) collapse to the cut coordinate \(\tau\), “thermal” is an optional parameterization rather than the invariant mathematics. `CumulativeField` leaves room for arbitrary graph-compiler score semantics and for the recovered approach to find its final identity.

---

# i'll probably prefer to simply write new documents and lean files rather than editing since we are rewriting by and large, and then archive them. lets do it

Done—the active formal narrative is now cleanly split.

- [CumulativeField](D:/aghado01/ThermoMapper/lean/Lemmas/CumulativeField.lean) contains the paper-independent verified core: arbitrary-θ cut reduction, graph-cut monotonicity, and parameter-free Mutual/HazardMean/Inclusive symmetrizations.
- [Local.lean](D:/aghado01/ThermoMapper/lean/Enthymemata/CumulativeField/Local.lean) holds the provisional local cumulative field. Its sole new `sorry` is the deliberately undis­charged monotonicity obligation.
- [CumulativeField/design.md](D:/aghado01/ThermoMapper/lean/Protolemmata/CumulativeField/design.md) records calibration alternatives, symmetrization semantics, the single effective cut coordinate `τ`, bracketing, pooling, and single-linkage obligations.
- [PKWangAudit.lean](D:/aghado01/ThermoMapper/lean/Lemmas/PKWangAudit.lean) is now the only active formal module carrying the paper name.
- [PKWangAudit.md](D:/aghado01/ThermoMapper/lean/Protolemmata/PKWangAudit.md) contains the critical narrative. It distinguishes genuine SW MCMC from the paper’s iid replacement, notes that the draws actually fed clustering, proves why masks are unnecessary, and clarifies that the paper hard-codes `0.5`—`log 2` is the unstated analytic consequence.
- The old A/B register was moved bit-for-bit into [Archeion/PKWang-legacy-2026-08-09](D:/aghado01/ThermoMapper/lean/Archeion/PKWang-legacy-2026-08-09/README.md).

Hierarchical barrels and recursive import-closure checking are now supported by [meta-ci.ps1](D:/aghado01/ThermoMapper/lean/scripts/meta-ci.ps1). Full verification passed: 2,517 Lean jobs, taxonomy gate green, and `git diff --check` clean.

Your pre-existing BARS changes were left untouched. Nothing was staged or committed.

---

# can we export this chat to `D:\aghado01\ThermoMapper\issues\lean\discussions` ; you can use Thermomapper/.codex/chat-export (leaf dir doesn't exist yet) for the intermediate working dir
