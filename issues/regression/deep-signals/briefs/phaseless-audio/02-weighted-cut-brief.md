# Slice 2 — sharp conditioning for relative signs

**Question.** How much can a sign flip change a real waveform while barely changing its measured edge products? Grohs and Rathmair motivate a weighted cut. Here the weights and the distance follow from [slice 1](01-two-sample-recovery-brief.md). See the [sequence](README.md) for the distinction from the continuous theorem.

Take a finite simple undirected graph, count each edge once, and give vertex v real amplitude xᵥ. The observation on uv is qᵤᵥ(x) = xᵤxᵥ. Work on the nonzero support, with vertex mass μᵥ = xᵥ² and edge weight wᵤᵥ = (xᵤxᵥ)². Write y for x with signs reversed on S and δS for the edges crossing S. Use the Euclidean data norm and waveform distance d±(x,y) = min(‖x−y‖₂, ‖x+y‖₂).

1. Prove ‖q(x)−q(y)‖₂² = 4 Σ_{uv∈δS} wᵤᵥ and d±(x,y)² = 4 min(μ(S), μ(Sᶜ)). Show that every y with |y| = |x| has this form. Complementary cuts give the same quotient.
2. For support with at least two vertices, take h(x) as the minimum crossing weight divided by the lighter side's mass over cuts with both masses positive. Prove the best constant in ‖q(x)−q(y)‖₂² ≥ h(x) d±(x,y)² for fixed amplitudes, including attainment. Characterize h = 0 by disconnection. For empty or singleton support every admissible competitor already agrees up to global sign. State this directly so the empty minimum has no spurious numerical value.
3. Evaluate a quiet-bridge family. On the three-vertex path use x = (1,δ,1), 0 < δ ≤ 1, and flip the last sample. The squared waveform distance is 4 and the squared edge-data change is 4δ². At δ = 0 connect this example to slice 1. For its DC/Nyquist raw power pair, derive the exact conversion factor between the sign-flip power discrepancy and the edge-product discrepancy under slice 1's normalization.

**Done when** the finite identities, the sharp variational statement and the quiet-bridge witness are verified. Report which data norm the constant belongs to. Amplitude variation and general Fourier cells enter in [slice 4](04-measured-stability-brief.md).

## Verified result — 2026-11-16

**Proof revision `1313d60`.** [WeightedSignCut.lean](../../../../../lean/Enthymemata/TimeFrequency/WeightedSignCut.lean) defines the sign-quotient distance using Euclidean waveform norms. The data vector holds one real sample product per undirected edge. A vertex ordering selects one orientation and unused coordinates are zero.

- `exists_flip`, `distance_flip_sq` and `products_flip_sq` show that every fixed-amplitude competitor is a sign flip, with squared discrepancies four times the lighter side’s mass and four times the crossing weight respectively. `complementary_cut` proves invariance of their ratio.
- `sharp_constant_iff` characterizes every valid lower constant for squared distances, and `extremizing_flip` gives an admissible competitor with positive distance that attains the minimum. `properCuts_nonempty` constructs a positive-mass cut from any two nonzero vertices. On the nonzero support, `constant_zero_iff_disconnected` characterizes zero by failure of path connectivity. When at most one vertex is nonzero, `subsingleton_support_sign` proves global-sign agreement directly.
- `quiet_bridge` evaluates the path (1,δ,1) with its last sign reversed. The squared waveform distance is 4 and the squared edge-product difference is 4δ² for every real δ. `raw_pair_factor` gives a factor of eight from squared edge-product discrepancy to the sum of the two squared Fourier power discrepancies at fixed amplitudes. `quiet_bridge_raw` gives 32δ² for the changing frame. At δ=0, `quiet_bridge_zero` reuses slice 1 to exhibit two inequivalent signals with identical one- and two-sample views.

### Verification

`lake build +Enthymemata.TimeFrequency` passed. The module printed 18 axiom checks, each containing only `propext`, `Classical.choice` and `Quot.sound`. `lean/scripts/meta-ci.ps1 -NoBuild` passed on the Lean tree archived from `1313d60`. The concurrent noise-law promotion had included the barrel import in `9cf612b`; `1313d60` contains both the import and its source.

### Doc changes owed

Record the sharp fixed-amplitude squared-distance constant, its support convention, and the factor of eight for a pair of unnormalized length-two Fourier powers. A stability claim for another sampled operator needs its own comparison with these edge products.
