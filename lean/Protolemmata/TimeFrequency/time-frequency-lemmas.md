# Time–frequency and phaseless-recovery lemmas

Candidate statements for [the Lean chip's brief](../../../issues/lean/time-frequency-statements-brief.md), as sharpened for Lean. The coordinating session posed S1–S6 on 2026-11-16 and marked them unchecked. Each entry gives the posing used, the proof route, the module and a status line. The brief's report carries the changes from the posed forms.

## S1. Phase propagation on a connected graph

**Status:** proved, `Enthymemata/TimeFrequency/PhasePropagation.lean`.

- *Posing.* `x, y : V → ℂ`, `x` nowhere zero, `G` a simple graph on `V`. On every edge `y u · conj (y v) = x u · conj (x v)`.
  - From one vertex: if `G` is preconnected and `‖y v₀‖ = ‖x v₀‖` at one vertex, then `y = c · x` with `‖c‖ = 1` (`exists_phase_of_vertex`). As posed, with every magnitude and `G` connected (`exists_phase`).
  - Sharpness: with every magnitude given, the edge products determine every `y` up to a global phase exactly when `G` is preconnected (`exists_phase_iff_preconnected`).
- *Route.* Put `c = y v₀ / x v₀`. One edge carries `y u = c·x u` to `y w = c·x w` when `‖c‖ = 1`: cancel `x u`, conjugate, multiply by `c`. Induct along a walk. For sharpness, take `y = x` on the component of `u` and `y = −x` off it.
- *Polarization.* `a·conj b = ½[(‖a + b‖² − ‖a‖² − ‖b‖²) + i(‖a + i·b‖² − ‖a‖² − ‖b‖²)]`, as posed. With the paper's `a − i·b` the imaginary term changes sign. `exists_phase_of_norms` states S1 on magnitudes alone.
- *Corollary.* On `ZMod n`, `n ≥ 1`, edges `k ~ k + s`: connected exactly when `s` is coprime to `n` (`shiftGraph_connected_iff`). Reachable sets are integer multiples of `s`, and `1` is one exactly when `s` is a unit.

## S2. Window scaling

**Status:** proved, `Enthymemata/TimeFrequency/WindowScaling.lean`.

- *Posing.* Sample `j` at time `tⱼ`, `X(ω) = Σⱼ gⱼ·xⱼ·e^{−iωtⱼ}`. Times are arbitrary, so a hop offset changes nothing.
  - Tone `xⱼ = A·e^{iωtⱼ}` at its own frequency: `|X|² = |A|²·|Σⱼ gⱼ|²`.
  - Impulse `a` at `j₀`: `|X|² = |a|²·|g_{j₀}|²` at every frequency.
  - Noise: on the second-moment matrix `C`, the expected power is `Σⱼₖ uⱼ·C j k·conj uₖ` with `uⱼ = gⱼ·e^{−iωtⱼ}` (`noisePower`). It is the integral of `|X|²` under any measure whose moments `∫ εⱼ·conj εₖ` are `C j k` (`integral_norm_windowed_sq`). White `C` gives `σ²·Σⱼ|gⱼ|²`.
- *Normalization.* `Σⱼ|gⱼ|² = 1`, as the brief fixes it. The rectangular window `1/√w` gives `|A|²·w`, `|a|²/w` and `σ²`.
- *Addition.* Cauchy–Schwarz: no unit-energy window of `w` samples reads more than `|A|²·w` from a tone (`norm_windowed_tone_sq_le`). The rectangular window attains it.

## S3. The flatness bound

**Status:** proved, `Enthymemata/TimeFrequency/FlatnessBound.lean`.

- *Posing.* `‖B‖²_F = Σᵢⱼ Bᵢⱼ²` (`frobSq`).
  - Real symmetric 2×2 `B`: `‖B‖²_F ≥ 2·|det B|`. Equality holds exactly when `B` is a multiple of `I` or has trace zero.
  - Real 2×2 `S`, symmetric `H`: `‖SᵀHS‖²_F ≥ 2·(det S)²·|det H|`.
  - For symmetric `H` and any factor `Σ = S·Sᵀ`, in any dimension: `tr((HΣ)²) = ‖SᵀHS‖²_F`. Hence `tr((HΣ)²) ≥ 2·det Σ·|det H|` in two dimensions.
  - Attainment: for `det H ≠ 0` and `d > 0`, a positive definite `S` with `det S = d` reaches equality.
- *Route.* With `B = [[a, b], [b, c]]`, `‖B‖²_F − 2 det B = (a − c)² + 4b²` and `‖B‖²_F + 2 det B = (a + c)²`.
  - For attainment with `H` definite, take `S ∝ adj H·sgn + √det H·I`, a square root of a multiple of `H⁻¹`. Then `SHS` is scalar.
  - For `H` indefinite, take `S ∝ M + √det M·I` with `M = ‖H‖²_F·I − tr H·H`, which is positive definite. Then `tr(SHS) = 0`.
- *Reading.* `Var(xᵀHx) = 2 tr((HΣ)²)` for `x ~ N(0, Σ)` is the Gaussian fourth-moment identity, not formalized here. A cell centred off zero adds `4μᵀHΣHμ ≥ 0`, so the bound still holds. Given that identity, the standard deviation is at least `2·√det Σ·√|det H|`, and every cell area reaches it. For `det H = 0` the infimum is `0` and no invertible `S` attains it.

## S6. Contraction from a regularity inequality

**Status:** proved, `Enthymemata/TimeFrequency/Contraction.lean`.

- *Posing.* Real normed space, `0 ≤ α ≤ 1`, `4η ≤ r`. On the ball `‖x − x⋆‖ ≤ r`, `‖v x − (x − x⋆)‖ ≤ ½‖x − x⋆‖ + 2η`.
  - Sharpened: from `‖x₀ − x⋆‖ ≤ r`, every iterate stays in the ball, and `‖x_t − x⋆‖ ≤ ρᵗ‖x₀ − x⋆‖ + 4η(1 − ρᵗ)` with `ρ = 1 − α/2`.
  - As posed, with `η ≥ 0`: `‖x_t − x⋆‖ ≤ ρᵗ‖x₀ − x⋆‖ + 4η`.
- *Route.* `x' − x⋆ = (1 − α)(x − x⋆) − α(v x − (x − x⋆))`, so one step gives `ρ‖e‖ + 2αη`. The invariant `ρᵗe₀ + 4η(1 − ρᵗ) ≤ r` carries the ball.

## S5. The Kronecker Gram

**Status:** posed; see the brief's report.

## S4. The whitened marginal

**Status:** posed; see the brief's report.
