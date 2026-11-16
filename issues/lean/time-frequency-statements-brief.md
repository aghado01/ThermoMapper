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

### In progress

*Updated 2026-11-16.*
- **Landed:** S1, S2, S3 and S6 are proved in `lean/Enthymemata/TimeFrequency/`. Each has a witness, and every closed result prints only `propext`, `Classical.choice` and `Quot.sound`. The candidate note is [time-frequency-lemmas.md](../../lean/Protolemmata/TimeFrequency/time-frequency-lemmas.md).
- **Next:** S5, then S4, beside the regression-family modules.
- **Uncommitted:** nothing.
- **Numerical chip:** not messaged. S2 keeps the energy normalization Σⱼ|gⱼ|² = 1, and nothing posed proved false.
