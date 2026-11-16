# Time–frequency and phaseless-recovery statements: a Lean chip

Date: 2026-11-16. Status: **commissioned by the user on 2026-11-16** ("one opus chip to start working the lean"). Drafted by the coordinating session. Its sibling is the numerical chip's [latent-recovery pilot](../regression/latent-recovery-pilot-brief.md).

## Why

A discussion with the user on 2026-11-16 explored the regression engines as teacher and measuring stick for unsupervised training of an audio model: several magnitude spectrograms as lossy views of one signal, a low-dimensional generator as the prior, and recovery guarantees in the style of Voroninski's papers. The coordinating session derived six statements along the way and marked every one unchecked. They carry the numerical pilot's conventions and two possible engine extensions, so they are posed here to be proved, corrected or refuted. A statement found false as posed is a result: record the counterexample and the corrected form.

## Sources

- **The harness's rules:** [lean/README.md](../../lean/README.md). Candidate statements start as a note under `lean/Protolemmata/`, compile under `lean/Enthymemata/`, and reach `lean/Lemmas/` when they stop apologizing. Active modules take scoped Mathlib imports. `scripts/verify-lean.ps1` runs the checks from a clean archive.
- **The papers,** LaTeX under `D:\aghado01\codex-scientiae\supellex-private\gauntlet\voroninski` (each folder is a versioned arXiv id; `inventory.jsonl` names the entry `.tex`):
  - `1109.0573v2/PhaseRetrieval.tex`, Theorem `teo:1D`: three structured patterns determine a signal up to global phase when its transform does not vanish and the shift is coprime to n. S1 is its mechanism in general form.
  - `2011.00288v2/dense_pr.tex`, Lemma 1 and Theorem 2: a regularity inequality on subgradients gives linear convergence to the noise level. S6 is its deterministic half.
- **The engine's seams:** the output inner product Q in [BassObservationModel.cs:177](../../src/maths/regression/spline/bass/BassObservationModel.cs) (declared, refused by today's core) for S4; the marginal's form in `src/maths/regression/gprior/CollapsedGPrior.cs` for S4; the product design in the [BAMRS brief](../regression/abmrs-bamrs/briefs-slices/bamrs-brief.md) for S5.

## The statements

Each is the coordinating session's posing. Sharpen hypotheses and conventions as Lean requires, and say what changed.

**S1. Phase propagation on a connected graph.** V finite, G a connected simple graph on V, x, y : V → ℂ with x nowhere zero. If ‖x v‖ = ‖y v‖ for every v and x u · conj(x v) = y u · conj(y v) on every edge, then y = c · x for one c with ‖c‖ = 1.
- *Sharpness.* When V splits into two nonempty parts with no edge between them, some y meets both hypotheses and is no global-phase multiple of x.
- *Polarization,* which supplies the edge products from magnitudes alone: a · conj(b) = ½[(‖a + b‖² − ‖a‖² − ‖b‖²) + i(‖a + i·b‖² − ‖a‖² − ‖b‖²)]. Check the sign convention; the paper uses x − i·Dˢx.
- *Optional corollary.* On ZMod n with edges k ~ k + s, the graph is connected exactly when s is coprime to n, the paper's condition.
- *Why it matters.* It is the reason overlapping or multi-window magnitudes pin a signal, and the reason silence splits a recording into separately recoverable pieces.

**S2. Window scaling.** For a window g on finitely many samples with Σ|g j|² = 1:
- a tone x j = A·e^{iωj} gives |Σ g j · x j · e^{−iωj}|² = |A|²·|Σ g j|²;
- an impulse of amplitude a at j₀ gives |a|²·|g j₀|² at every frequency;
- noise with second moments E[ε i · conj(ε j)] = σ²·δ_ij gives expected power σ² at every frequency. A statement on the second-moment matrix alone avoids measure theory; choose the form that composes best.
- *Corollary for the rectangular window of length w,* g j = 1/√w: the three powers are |A|²·w, |a|²/w and σ².
- *Why it matters.* The pilot classifies content by the slope of log power against log window length (+1, −1, 0). This fixes the normalization those readings assume. Tell the numerical chip if the normalization must change.

**S3. The flatness bound.** For a real symmetric 2×2 matrix B, ‖B‖²_F ≥ 2·|det B|, with its equality case. Hence for a real 2×2 S and symmetric H, ‖Sᵀ H S‖²_F ≥ 2·(det S)²·|det H|.
- *Optional.* For det H ≠ 0 and any d > 0, some symmetric positive definite S with det S = d attains equality.
- *Why it matters.* For x ~ N(0, Σ), Var(xᵀHx) = 2·tr((HΣ)²) = 2·‖Σ^½ H Σ^½‖²_F. So the bound says how flat a quadratic can be made over a cell of fixed area, whatever the cell's shape: the coordinating session's claim that √|det H| times the area is the invariant obstruction. This derivation was flagged as most in need of a second look.

**S4. The whitened marginal.** For y ∈ ℝᴺ, Z ∈ ℝ^{N×p} of full column rank and Ω symmetric positive definite with Ω = LLᵀ, put ỹ = L⁻¹y and Z̃ = L⁻¹Z. Then:
- the residual forms agree: yᵀΩ⁻¹y − yᵀΩ⁻¹Z(ZᵀΩ⁻¹Z)⁻¹ZᵀΩ⁻¹y = ỹᵀỹ − ỹᵀZ̃(Z̃ᵀZ̃)⁻¹Z̃ᵀỹ;
- the Gaussian log density of y under mean Zβ and covariance σ²Ω is that of ỹ under mean Z̃β and covariance σ²I, less ½·log det Ω;
- for Ω = I_n ⊗ R, log det Ω = n·log det R and Ω⁻¹ = I_n ⊗ R⁻¹.
- *Target.* State the consequence in the engine's own terms: the collapsed g-prior marginal under noise covariance σ²(I ⊗ R), with the coefficient prior on the whitened Gram, equals the unit-noise marginal on (ỹ, Z̃) plus the determinant term. Read `CollapsedGPrior`'s marginal and the contract's comment (it names the centring under Q, Σᵢ yᵢᵀQyᵢ and log|R| as the missing terms) and say how the level-centred intercept reads under Q.
- *Why it matters.* It is the specification for the extension that lets the core serve Q ≠ I.

**S5. The Kronecker Gram.** (A ⊗ B)ᵀ(A ⊗ B) = (AᵀA) ⊗ (BᵀB). With both Grams invertible, the projection onto the columns of A ⊗ B is P_A ⊗ P_B, and for y the vectorization of a matrix Y the fitted energy yᵀ(P_A ⊗ P_B)y is ‖P_A Y P_B‖²_F under the matching vec convention. State the convention.
- *Why it matters.* On gridded data a product design's Gram factors by axis. The BAMRS brief assumes scattered rows, where only the row-wise product holds.

**S6. Contraction from a regularity inequality.** In a real normed space, fix x⋆, η ≥ 0, r > 0 with 4η ≤ r, and 0 < α ≤ 1. Suppose ‖v x − (x − x⋆)‖ ≤ ½·‖x − x⋆‖ + 2η whenever ‖x − x⋆‖ ≤ r. Then the iteration x_{t+1} = x_t − α·v(x_t) from ‖x₀ − x⋆‖ ≤ r stays in the ball and satisfies ‖x_t − x⋆‖ ≤ (1 − α/2)ᵗ·‖x₀ − x⋆‖ + 4η.
- The paper's version carries the global sign ambiguity; this one takes a single basin.
- *Why it matters.* It is the reusable shape of "a deterministic condition implies convergence", the corpus's proof style.

## Scope and order

- **Order:** S3, S2, S1 first. They are cheap and they audit the coordinating session. Then S6, S5, S4.
- **A witness for every statement.** Beside each result, an `example` instantiates its hypotheses on a concrete object, so no proof closes on an unsatisfiable premise. Record `#print axioms` for each closed result.
- **Placement** follows the library's pattern: engine-named modules, a subfolder per topic as `ActiveSubspace/` has. S4 and S5 belong beside the regression-family modules. Name modules by what they state. The report lists the placement.
- **Left to the numerical chip:** the dilogarithm covariance of log powers, and every probabilistic claim (concentration of random matrices, the papers' sample-complexity halves).
- **Left alone:** `src/`, `tests/`, `docs/`.

## Done

- S1, S2 and S3 proved without `sorry`, each with its witness, eligible for promotion under the README's rules.
- S4, S5 and S6 at least as compiling statements, with every proof that closes.
- A report appended to this brief: per statement, the module, its status, its axioms, and each change from the posed form with the reason. A **Doc changes owed** heading. An **In progress** entry kept current while work is open.

## Working rules

- Commits land on `main` by pathspec, each one informative about its work, with no co-author trailer. The index is shared with other sessions, so nothing stays staged.
- `lake build` runs at below-normal priority with capped threads; the machine is shared.
- Each step ends with a fresh-eyes review of its own output, fixed before it is reported.
- Messages to the coordinating session or the numerical chip are terse events; results travel in the report.

## Report

### 2026-11-16 — S1 to S6: all six proved, none false as posed

**Landed:** `475c447` (S1, S2, S3, S6), `6324306` (S5, S4 with one step open), and the commit that adds this report (S4 closed, the Khatri–Rao Gram, witnesses, lint).

**Per statement.** Every declaration named below prints only `propext`, `Classical.choice` and `Quot.sound`.

| | Module | Status | Changes from the posed form, and why |
|---|---|---|---|
| S1 | `lean/Enthymemata/TimeFrequency/PhasePropagation.lean` | proved | **V finite** dropped; connectivity alone is used. **One vertex's magnitude** suffices on a preconnected graph (`exists_phase_of_vertex`); the posed form is `exists_phase`. **Sharpness** is an iff: with every magnitude given, the edge products determine every y up to phase exactly when the graph is preconnected (`exists_phase_iff_preconnected`). **Polarization** is correct as posed; the paper's `a − i·b` flips the imaginary term's sign (`mul_conj_eq_polarization_sub`). `exists_phase_of_norms` states S1 on magnitudes alone. **Corollary** `shiftGraph_connected_iff`, for n ≥ 1. |
| S2 | `lean/Enthymemata/TimeFrequency/WindowScaling.lean` | proved | **Sample times** tⱼ are arbitrary, so hop offsets are covered. **The tone's power** is at its own frequency. **Noise** is stated on the second-moment matrix (`noisePower_white`) and as an integral under any measure with those moments (`integral_norm_windowed_sq`). **Added:** no unit-energy window of w samples reads more than \|A\|²·w from a tone (`norm_windowed_tone_sq_le`). **The normalization** Σⱼ\|gⱼ\|² = 1 is unchanged, so the numerical chip was not messaged. |
| S3 | `lean/Enthymemata/TimeFrequency/FlatnessBound.lean` | proved | **Equality** holds exactly when B is scalar or traceless (`frobSq_eq_two_mul_abs_det_iff`). **The variance link** tr((HΣ)²) = ‖SᵀHS‖²_F holds in any dimension and for any factor Σ = SSᵀ (`trace_sq_eq_frobSq`). **The optional attainment** is proved (`exists_posDef_frobSq_eq`). **Not formalized:** the Gaussian identity Var(xᵀHx) = 2 tr((HΣ)²). |
| S4 | `lean/Enthymemata/BassWhitenedMarginal.lean` | proved | **Hypotheses.** Full column rank and positive definiteness are not needed for the residual identity, which holds for any Ω = LLᵀ under Mathlib's inverse (`residual_whiten`). The density identity needs only det L ≠ 0 and σ² ≠ 0 (`gaussLogDensity_whiten`). Full column rank enters only in the closed form. **The marginal covariance** σ²(Ω + τ⁻¹Z(ZᵀΩ⁻¹Z)⁻¹Zᵀ) is taken as the model's marginal; the Gaussian convolution behind it is not derived. **The middle term** is written (ZᵀΩ⁻¹y)ᵀ(ZᵀΩ⁻¹Z)⁻¹(ZᵀΩ⁻¹y), the same scalar. **Added at the user's request:** the Khatri–Rao Gram, `khatriRao_gram` and `khatriRao_gram_one`. |
| S5 | `lean/Enthymemata/BamrsKroneckerGram.lean` | proved | **Invertibility.** The projection identity needs none (`proj_kronecker`). **The convention** is y(i, j) = Y i j with A's index outer, which is Mathlib's `vec Yᵀ` (`rowMajor_energy`). Mathlib's column-major `vec Y` pairs with P_B ⊗ P_A instead (`colMajor_energy`). |
| S6 | `lean/Enthymemata/TimeFrequency/Contraction.lean` | proved | **Hypotheses.** 0 ≤ α suffices (posed 0 < α), r > 0 is not needed, and η ≥ 0 enters only the posed form (`norm_descent_le_posed`). **Sharper bound** (`norm_descent_le`): ‖x_t − x⋆‖ ≤ ρᵗ‖x₀ − x⋆‖ + 4η(1 − ρᵗ). |

**Declarations to print axioms for:**
- *S1:* `exists_phase_of_vertex`, `exists_phase`, `exists_phase_iff_preconnected`, `mul_conj_eq_polarization`, `mul_conj_eq_polarization_sub`, `exists_phase_of_norms`, `shiftGraph_connected_iff`.
- *S2:* `norm_windowed_tone_sq`, `norm_windowed_impulse_sq`, `noisePower_white`, `integral_norm_windowed_sq`, `integral_norm_windowed_sq_white`, `norm_windowed_tone_sq_le`, `rect_tone`, `rect_impulse`, `rect_noise`.
- *S3:* `two_mul_abs_det_le_frobSq`, `frobSq_eq_two_mul_abs_det_iff`, `frobSq_congr_ge`, `trace_sq_eq_frobSq`, `trace_sq_ge`, `exists_posDef_frobSq_eq`.
- *S6:* `norm_descent_step_le`, `norm_descent_le`, `norm_descent_le_posed`.
- *S5:* `BAMRS.kronecker_gram`, `BAMRS.proj_kronecker`, `BAMRS.colMajor_energy`, `BAMRS.rowMajor_energy`.
- *S4:* `BASS.residual_whiten`, `BASS.gaussLogDensity_congr`, `BASS.gaussLogDensity_whiten`, `BASS.gaussLogDensity_gPrior`, `BASS.gPrior_logMarginal_whiten`, `BASS.log_det_one_kronecker`, `BASS.inv_one_kronecker`, `BASS.quadForm_one_kronecker`, `BASS.khatriRao_gram`, `BASS.khatriRao_gram_one`, `BASS.level_whiten`, `BASS.level_one_kronecker`, `BASS.level_eq_mean`.
- S1, S2, S3 and S6 are in namespace `TimeFrequency`.

**Witnesses.** Each module ends with `example`s that instantiate the main results' hypotheses on concrete objects:
- *S1:* the path on ZMod 3 with y = i·x. Its connectivity comes from `shiftGraph_connected_iff` itself. Beside it, the edgeless two-vertex graph as the sharpness counterexample.
- *S2:* rectangular windows of 4 and 2 samples. White noise is realized on two outcomes under the counting measure.
- *S3:* diag(1, 2) for the strict bound, diag(1, −1) for equality and attainment at d = 2.
- *S4:* the correlated pair L = [[1, 0], [1, 1]] with Z = (1, 1)ᵀ, σ² = τ = 1. Q = I for the level.
- *S5:* the one-column design (1, 1)ᵀ on both axes.
- *S6:* v(x) = 3x/2 + 1/10 on ℝ, with η = 1/20, r = 1, α = 1/2.

**Gates:**
- `lean` on each file during drafting.
- `lake build` on `Enthymemata.TimeFrequency`, `Enthymemata.BassWhitenedMarginal` and `Enthymemata.BamrsKroneckerGram` with the linter set: no warnings after fixes.
- `scripts/meta-ci.ps1`: PASS with its build gate, all six modules `PROOF-CLOSED`.
- No cache fetch. D: had 22 GB free at the end.

**Follow-ups:**
- *Promotion review* of the six modules, all `PROOF-CLOSED`.
- *Two identical projection definitions,* `BAMRS.proj` and `BASS.proj`. They merge if one module comes to import the other; nothing waits on it.
- *Left unformalized by scope:* the Gaussian fourth-moment identity behind S3's reading, and the convolution that makes S4's marginal covariance. Both are probabilistic.
- *Notation:* the engine's Q means both the fitted energy and the output inner product. Logged in `issues/notation-ledger.md`.

#### Methods

**S1, phase propagation.**
- *Setup.* Put c = y(v₀)/x(v₀); the magnitude at v₀ gives \|c\| = 1.
- *The edge step* needs nothing more. Suppose y(u) = c·x(u) and y(u)·conj y(w) = x(u)·conj x(w). Cancelling x(u) ≠ 0 gives c·conj y(w) = conj x(w). Conjugating gives conj c·y(w) = x(w). Multiplying by c, with c·conj c = \|c\|² = 1, gives y(w) = c·x(w).
- *Propagation.* Induction on a walk from v₀ carries this to every reachable vertex.
- *Magnitudes propagate too.* The edge products alone force \|c_u\|·\|c_w\| = 1 across every edge. So one vertex's magnitude fixes all of them, and on a bipartite graph one is necessary: y = r·x on one side and x/r on the other passes every edge for r > 0. That remark is not formalized.
- *Sharpness.* If v is unreachable from u, take y = x on u's component R and −x off it. No edge leaves R, so every edge product keeps its sign. A global phase would need c = 1 at u and c = −1 at v.
- *Polarization.* ‖a + b‖² = ‖a‖² + ‖b‖² + 2 Re(a·conj b), and ‖a + ib‖² = ‖a‖² + ‖b‖² + 2 Im(a·conj b), since Re(a·conj(ib)) = Re(−i·a·conj b) = Im(a·conj b). With a − ib the imaginary term changes sign.
- *The shift graph.* Each edge moves by ±s, so every vertex reachable from 0 lies in ℤs. The walk 0, s, 2s, … reaches all of ℕs. So the graph is connected iff 1 ∈ ℤs, iff s is a unit of ZMod n, iff gcd(s, n) = 1.

**S2, window scaling.**
- *Setup.* X(ω) = Σⱼ uⱼxⱼ with uⱼ = gⱼ·e^{−iωtⱼ}, and \|uⱼ\| = \|gⱼ\|.
- *Tone and impulse.* A tone at its own frequency cancels the phases, so X = A·Σgⱼ. An impulse leaves one term.
- *Noise.* \|X\|² = X·conj X = Σⱼₖ uⱼ·conj uₖ·εⱼ·conj εₖ, a finite sum of integrable products. Integrating term by term gives Σⱼₖ uⱼ·C(j, k)·conj uₖ. White C leaves σ²·Σ\|uⱼ\|² = σ²·Σ\|gⱼ\|².
- *The rectangular window* gⱼ = 1/√w has Σgⱼ = √w and \|gⱼ\|² = 1/w. That gives the slopes +1, −1 and 0 exactly.
- *The tone bound* is Cauchy–Schwarz: \|Σgⱼ\|² ≤ (Σ\|gⱼ\|)² ≤ w·Σ\|gⱼ\|².
- *Other window families.* For a family gⱼ ∝ h(j/w), the slopes hold asymptotically, with constants (∫h)²/∫h² for the tone and h(θ)²/∫h² for an impulse at relative position θ. That remark is not formalized.

**S3, the flatness bound.**
- *The bound.* For B = [[a, b], [b, c]], ‖B‖²_F − 2 det B = (a − c)² + 4b² and ‖B‖²_F + 2 det B = (a + c)². So ‖B‖²_F ≥ 2\|det B\|, with equality exactly for a = c, b = 0 (scalar) or for a + c = 0 (traceless).
- *Congruence.* SᵀHS is symmetric with determinant (det S)²·det H.
- *The variance link* is cyclicity: tr((HSSᵀ)²) = tr(SᵀHS·SᵀHS) = Σᵢⱼ Bᵢⱼ² for symmetric B.
- *Attainment* needs SᵀHS scalar (H definite) or traceless (H indefinite). Both survive S ↦ tS, which scales det S by t², so one shape per H suffices.
  - *Definite H, a > 0:* take S₀ = adj H + √(det H)·I, a multiple of the square root of H⁻¹. Cayley–Hamilton, with adj H·H = det H·I and (adj H)² = tr H·adj H − det H·I, gives S₀HS₀ = det H·(tr H + 2√det H)·I. For a < 0 the sign of adj H flips.
  - *Indefinite H:* M = ‖H‖²_F·I − tr H·H has eigenvalues λ₂(λ₂ − λ₁) and λ₁(λ₁ − λ₂), both positive. Its determinant is −det H·((a − c)² + 4b²), and tr(HM) = 0 because tr(H²) = ‖H‖²_F. With S₀ = M + √(det M)·I, S₀² = (tr M + 2√det M)·M, so tr(S₀HS₀) = tr(HS₀²) ∝ tr(HM) = 0.
  - *In Lean,* each entry identity is a polynomial identity modulo δ² = det, closed by `linear_combination`. The coefficient is Hᵢⱼ in the definite cases and tr H in the indefinite one.
  - *Positive definiteness* is Sylvester's criterion in two dimensions, through S₀₀·xᵀSx = (S₀₀x₀ + S₀₁x₁)² + det S·x₁².
- *Reading.* Given the Gaussian identity, SD(xᵀHx) ≥ 2√(det Σ)·√\|det H\| = (2/π)·(area of the 1σ ellipse)·√\|det H\|, sharp at every area.

**S6, contraction.**
- *One step.* With e = x − x⋆, the next error is (1 − α)e − α(v(x) − e). So it is at most (1 − α)‖e‖ + α(‖e‖/2 + 2η) = ρ‖e‖ + 2αη, with ρ = 1 − α/2. α ≤ 1 makes \|1 − α\| = 1 − α.
- *The invariant* ‖e_t‖ ≤ ρᵗe₀ + 4η(1 − ρᵗ) closes because 4ηρ + 2αη = 4η.
- *The ball.* With e₀ ≤ r and 4η ≤ r, the invariant is at most ρᵗr + r(1 − ρᵗ) = r. That keeps every iterate where the hypothesis applies.

**S5, the Kronecker Gram.**
- *The Gram.* (A ⊗ B)ᵀ = Aᵀ ⊗ Bᵀ, and the mixed product (A ⊗ B)(C ⊗ D) = AC ⊗ BD, give the Gram.
- *The projection.* Mathlib's inverse satisfies (X ⊗ Y)⁻¹ = X⁻¹ ⊗ Y⁻¹ unconditionally: both sides vanish when a factor is singular. Hence P_{A⊗B} = P_A ⊗ P_B for any A and B.
- *The energy.* (P_B ⊗ P_A)·vec Y = vec(P_A·Y·P_Bᵀ), and vec Y·vec Z = tr(YᵀZ). Since P is symmetric, and idempotent when its Gram is invertible, tr(Yᵀ·P_A·Y·P_B) = ‖P_A·Y·P_B‖²_F. Transposing Y swaps the factors.

**S4, the whitened marginal.**
- *Whitening.* Ω = LLᵀ gives Ω⁻¹ = (L⁻¹)ᵀL⁻¹; Mathlib's inverse makes both (LLᵀ)⁻¹ = (Lᵀ)⁻¹L⁻¹ and (L⁻¹)ᵀ = (Lᵀ)⁻¹ unconditional. So every Ω⁻¹-form is a Euclidean form of whitened vectors, vᵀΩ⁻¹w = (L⁻¹v)ᵀ(L⁻¹w). The residual, Gram and cross-product identities are three instances.
- *The density.* A covariance congruence C = LMLᵀ has det C = det M·det(LLᵀ) and C⁻¹ = (L⁻¹)ᵀM⁻¹L⁻¹. So log N(y; μ, LMLᵀ) = log N(L⁻¹y; L⁻¹μ, M) − ½ log det(LLᵀ).
- *The marginal.* With Z = LZ̃, Ω + τ⁻¹Z(ZᵀΩ⁻¹Z)⁻¹Zᵀ = L(I + τ⁻¹P̃)Lᵀ. The prior on the whitened Gram is exactly what makes the marginal covariance a congruence of the unit-noise one; a prior on ZᵀZ would not.
- *The closed form.*
  - Weinstein–Aronszajn turns det(I_N + a·Z̃(Z̃ᵀZ̃)⁻¹Z̃ᵀ) into det((1 + a)I_k) = (1 + a)ᵏ.
  - Idempotence gives (I + aP)⁻¹ = I − (a/(1 + a))P, since a − d − ad = 0 at d = a/(1 + a).
  - At a = 1/τ, the log determinant is N log σ² − k log(τ/(τ + 1)), and the quadratic is (ỹᵀỹ − ỹᵀP̃ỹ/(τ + 1))/σ². That is `CollapsedGPrior`'s scale S = ‖r‖² + τ/(τ + 1)·Q at b = 1, over σ², with its determinant term (k/2)·log(τ/(τ + 1)).
- *I ⊗ R.*
  - `det_kronecker` gives det Rⁿ and `inv_kronecker` gives I ⊗ R⁻¹.
  - (I ⊗ L_R)(I ⊗ L_R)ᵀ = I ⊗ L_R·L_Rᵀ, so a factor of R whitens the whole.
  - vᵀ(I ⊗ Q)w = Σᵢ vᵢᵀQwᵢ, from the block-diagonal action ((I ⊗ Q)w)(i, ·) = Q·w(i, ·).
- *The Khatri–Rao Gram.* Column j of A ⊙ H is aⱼ ⊗ hⱼ. Entry (j, j′) of Xᵀ(I ⊗ Q)X is Σᵢ (Aᵢⱼhⱼ)ᵀQ(Aᵢⱼ′hⱼ′) = (Σᵢ AᵢⱼAᵢⱼ′)·hⱼᵀQhⱼ′, which is (AᵀA ∘ HᵀQH)ⱼⱼ′. Hence the whitened Gram of BASS's functional design under Ω = I ⊗ R is (AᵀA)∘(HᵀR⁻¹H), as `KhatriRao.Gram` computes it.
- *The intercept under Q.*
  - The constant column whitens like any other. So the level-centred intercept's level in whitened coordinates is the generalized least-squares level 1ᵀΩ⁻¹y / 1ᵀΩ⁻¹1.
  - Under I_n ⊗ R that is Σᵢ1ᵀQyᵢ / (n·1ᵀQ1) = 1ᵀQȳ / 1ᵀQ1, the Q-weighted level of the mean curve ȳ. Its Gram entry, the mass, is n·1ᵀQ1 in place of n·m.
  - It is the plain grand mean exactly when the constant is a left eigenvector, 1ᵀQ = λ1ᵀ with λ ≠ 0: Q = I, or any Q with equal column sums.
  - Q-centring a column projects out the constant in the W-inner product, x − 1·(1ᵀWx)/(1ᵀW1). That form is not stated in Lean.
