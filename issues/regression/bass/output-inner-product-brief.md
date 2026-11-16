# The output inner product Q in BASS's core: preparation brief

Date: 2026-11-16. Status: **prepared, not built.** The Lean chip wrote it for Part C of the [noise-law brief](../../lean/noise-law-statements-brief.md). Building starts on the user's go-ahead, which the coordinating session relays. The user said on 2026-11-16 that additive engine work is acceptable so long as it is coherent.

## Why

A direct functional fit treats its noise as independent across output locations (D07). When the misfit or the noise is smooth along the index, a dense grid counts n·m independent observations. Intervals then come out too narrow, and births absorb smooth residual. [extensions-brief.md § Residuals correlated along the output index](extensions-brief.md#residuals-correlated-along-the-output-index) states the model: residual covariance σ²(I_n ⊗ R), Q = R⁻¹.

The working-response contract already declares Q, and today's core refuses it. Two callers ask for it now, both from the numerical chip's [latent-recovery pilot](../latent-recovery-pilot-brief.md#primitives-proposed):
- **Primitive 2,** Q in the fitting core. Its B1 reading supplies Q for spectral noise from the windows alone. The log-power correlation of two noise-dominated cells is Li₂(\|G\|²)/Li₂(1), with the pseudo-Gram entry near DC and Nyquist.
- **Primitive 1,** an emulator's error covariance across outputs. Its Stage 1 fits read residual lag-1 r of 0.70 to 0.83 along ln f.

## The specification: what is proved

Each item below is proved in Lean, with no `sorry` and standard axioms only, and promoted to `lean/Lemmas/`.

- **`BASS.gPrior_logMarginal_whiten`** (`Lemmas/BassWhitenedMarginal.lean`).
  - *Setting:* noise σ²Ω with Ω = LLᵀ, and the coefficient prior β ~ N(0, σ²(τ·ZᵀΩ⁻¹Z)⁻¹) on the *whitened* Gram.
  - *Result:* the collapsed marginal's log density is the unit-noise closed form on (ỹ, Z̃) = (L⁻¹y, L⁻¹Z), less ½·log det Ω:
    `−(N/2)·log(2πσ²) + (k/2)·log(τ/(τ+1)) − S̃/(2σ²) − ½·log det Ω`, with `S̃ = ỹᵀỹ − ỹᵀP̃ỹ/(τ+1)`.
  - That is `CollapsedGPrior`'s determinant term and scale at b = 1 on whitened data.
- **`BASS.marginal_covariance`** (S14). The response's second moments under that prior and noise are σ²(Ω + τ⁻¹Z(ZᵀΩ⁻¹Z)⁻¹Zᵀ), the covariance the marginal above reads.
- **`BASS.log_det_one_kronecker`, `inv_one_kronecker`, `quadForm_one_kronecker`.** For Ω = I_n ⊗ R:
  - log det Ω = n·log det R;
  - Ω⁻¹ = I_n ⊗ R⁻¹;
  - whitening by I ⊗ L_R gives ỹᵀỹ = Σᵢ yᵢᵀQyᵢ.
- **`BASS.khatriRao_gram`.** The functional design X = A ⊙ H, column aⱼ ⊗ hⱼ, case outer, has Xᵀ(I ⊗ Q)X = (AᵀA)∘(HᵀQH). That is `KhatriRao.Gram`'s form with an output product.
- **`BASS.level_whiten`, `level_one_kronecker`, `level_eq_mean`.** The level-centred intercept's level is the generalized least-squares level:
  - 1ᵀΩ⁻¹y / 1ᵀΩ⁻¹1 = Σᵢ1ᵀQyᵢ / (n·1ᵀQ1) = 1ᵀQȳ / 1ᵀQ1, with mass n·1ᵀQ1.
  - It is the plain grand mean exactly when 1ᵀQ = λ1ᵀ.
- **`BASS.residual_whiten`.** The residual forms agree for any L.
- **Tempering,** derived here and not stated in Lean. The congruence gives p(y | β) = \|L\|⁻¹·p_unit(ỹ | β). The power posterior's likelihood raised to b therefore carries b·(−½·log det Ω): the determinant term scales with b, as `GPriorMarginal`'s Σ ln wᵢ term already does.

## What is already built

- **The contract.** `WorkingResponse` takes `IInnerProduct? output` and validates it: a scalar response is refused, and so is a dimension other than the locations (`BassObservationModel.cs:230-236`). Its comment names what is missing: the centring under Q, Σᵢ yᵢᵀQyᵢ and log\|R\| (`BassObservationModel.cs:187`).
- **`IInnerProduct`** (`src/maths/linalg/InnerProduct.cs`). It has `IsIdentity`, `Dimension`, `LogDeterminant` (ln\|Q\|) and `Dot` (Σᵢ uᵢᵀQvᵢ over stacked vectors). The identity is the only implementation.
- **`KhatriRao`'s direct contractions** take `output`: `Gram`, `Border` and `CrossProducts` give (AᵀA)∘(HᵀQH) and aⱼᵀ·Y·Q·hⱼ. Its centred forms serve Q = I only. Its doc (`KhatriRao.cs:29-35`) states the Q-centred generalization: m becomes 1ᵀQ1, means become Q-weighted means, and the factors' products become hᵀQh.
- **`BassDesign` already reads Q in places.** It holds `_output` (`:175`). The response energy (`:384`, `:400`), the residual energy (`:771`, `:1583`, `:1619`), the cross products (`:1524`) and the centred self-energy's output factor (`:1514`) go through `_output.Dot`.
- **Calibration** already has a fixed output correlation in its multivariate formulation (`CalibrationCovariance.cs:61-64`), and the `TruncationTerm.ResidualSet` path that primitive 1 would feed.

## Where the core reads the independent form today

Located by search on 2026-11-16. The building chip confirms each, and adds any the search missed.

| Place | Today | Under Q |
|---|---|---|
| `BassDesign.cs:303` | Refuses `!_output.IsIdentity` | Accepts, within the first slice's scope |
| `BassDesign.CentredSelf` (`:1507-1516`) | Plain means h̄ = Σh/m, the m·h̄² term | Q-weighted mean 1ᵀQh/1ᵀQ1, the (1ᵀQ1)·h̄_Q² term, the `KhatriRao` doc's three-term form |
| `KhatriRao.CentredGram` and the centred border | Q = I only | Q-centred three-term form |
| `BassDesign` response centring (`:375-400`) | y − ȳ with the plain mean; `level` = Σy/N | Centred against the GLS level 1ᵀQȳ/1ᵀQ1 (`level_one_kronecker`) |
| `WorkingResponse.ResponseMean`, `ResponseSquares` (`BassObservationModel.cs:428-449`) | Plain mean and Σ(y − ȳ)² | GLS level ℓ and Σᵢ(yᵢ − ℓ·1)ᵀQ(yᵢ − ℓ·1) |
| `GPriorWeights` / `CentredFit` (`CollapsedGPrior.cs`) | Mass N; the zero centre's energy N·ȳ² | Mass n·1ᵀQ1; energy n·(1ᵀQȳ)²/(1ᵀQ1) |
| The likelihood constant | Σ ln wᵢ and a log-Jacobian (`BassObservationModel.cs:257`) | Adds n·ln\|Q\| through the same `SumLogWeights` route, so it enters at b/2 |
| `BassResolved` mean squares (`:229-265`), `BassTarget` exact-fit floor (`:113`, `:307-309`) | Plain mean squares | Q-weighted mean square, (Σᵢ eᵢᵀQeᵢ)/N |
| `BassPointwise` (`:177`, `:218`, `:232`) | Units are case blocks (`:41-47`), the density a product of per-location normals | The block's density −½·log det(2πσ²R) − ½·eᵀQe/σ² |
| `posterior/BassPrediction` noise | Independent σ·z per location | σ·L_R·z per case |
| `BassResidualDependence` (F11) | Raw residuals against the independent expectation −1/m | Whitened residuals L_R⁻¹eᵢ, with the raw reading beside them |
| `calibration/BassEmulators` (`:10`) | σ²_e·I (D09's amendment) | σ²_e·R for a fit under Q (primitive 2), or the residual set (primitive 1) |

## Decisions

Each carries a recommendation and the slice where it first binds.

1. **Where the coefficient prior's Gram reads under Q.**
   - *Recommended:* on the whitened Gram XᵀΩ⁻¹X = (AᵀA)∘(HᵀQH).
   - *Why:* it makes the marginal covariance a congruence of the unit-noise one, so S4's closed form holds unchanged and the oracle is exact.
   - *The alternative:* the unweighted Gram (AᵀA)∘(HᵀH), as `CollapsedPencilPrior` does for observation weights. It breaks the congruence, needs a pencil, and gives no closed-form oracle.
   - *Binds at:* slice 1.
2. **How Q combines with observation weights.**
   - *Recommended:* refuse weights with Q ≠ I in slice 1, then admit case-level weights in slice 4. Weights of the form wᵢ·1_m keep the Kronecker form, with L = diag(w)^{−½} ⊗ L_R.
   - *Why refuse the rest:* per-location weights give σ²·D^{−½}(I ⊗ R)D^{−½}, which is not a Kronecker product, so the Hadamard shortcut is lost. Today's design already leaves the Khatri–Rao path under weights (`BassDesign.cs:304`).
   - *Binds at:* slices 1 and 4.
3. **Several outputs (E34).**
   - *Recommended:* refuse Q ≠ I with several outputs until E34's Σ and Q are stated together. A per-output Q and a shared Q are different models, and the choice is the user's.
   - *Binds at:* slice 1, as a refusal.
4. **What the level-centred intercept reads.**
   - *Recommended:* the GLS level 1ᵀQȳ/1ᵀQ1, with mass n·1ᵀQ1. Centred columns are W-orthogonal to the constant, x − 1·(1ᵀWx)/(1ᵀW1) with W = I ⊗ Q.
   - *The other conventions:* the zero centre's energy is n·(1ᵀQȳ)²/(1ᵀQ1), and a frozen level's gap is measured in the same mass. The flat intercept integrates α against the whitened constant column L⁻¹1, spending one observation as today.
   - *Note:* the level equals the plain grand mean when Q's column sums are equal.
   - *Binds at:* slice 2.
5. **Tempering and the evidence.**
   - *Recommended:* carry n·ln\|Q\| in `SumLogWeights`. It then enters every tempered marginal at b/2, and every thermodynamic-integration slope at ½.
   - *Why:* the evidence estimators read `GPriorMarginal` and need no change of their own. Evidence comparisons across Q are then honest, since the log\|R\| term is in the likelihood.
   - *Binds at:* slice 3.
6. **What the residual-dependence check reports once Q is declared.**
   - *Recommended:* both readings, with the declared Q named in the report. On whitened residuals L_R⁻¹eᵢ the check tests the declared Q: about −1/m if R is right. On raw residuals it is today's F11 reading.
   - *Binds at:* slice 5.
7. **Where Q comes from.**
   - *Recommended:* the declared `IInnerProduct` only, fixed for the fit. Sampling R's parameters θ is the extension entry's later step: it recomputes the output-side Gram, O(M²m) under a banded Q.
   - *Implementations to add:*
     - dense Q from a given R, by Cholesky;
     - E12's tridiagonal Markov Q, linear in m;
     - B1's table, R_kl = Li₂(\|G_kl\|²)/Li₂(1), dense or banded by window overlap, with the pseudo-Gram entry within two bins of DC and Nyquist.
   - *Binds at:* slice 1 (dense) and slice 5 (the table).

## The slices

Each ends at a gate. The [common gates](#gates-common-to-every-slice) apply to all.

1. **Q in the uncentred core.**
   - *Build:*
     - lift `BassDesign`'s refusal for a direct functional response with no weights and one output;
     - add `DenseInnerProduct`, R given, Cholesky factor, ln\|Q\| from the factor;
     - give `KhatriRao.Gram`, `Border` and `CrossProducts` the declared Q, with the prior on the whitened Gram (decision 1);
     - route n·ln\|Q\| into `SumLogWeights` (decision 5).
   - *Scope:* the zero centre only (no intercept apart).
   - *Gate, the whitening oracle:* for fixed terms, `BassDesign`'s log marginal under Q equals `CollapsedGPrior` built on explicitly whitened (Z̃, ỹ) plus (n/2)·ln\|Q\|, to 1e-12 relative, over a seeded sweep of designs and dense R. This is S4's `gPrior_logMarginal_whiten`.
2. **Centring under Q.**
   - *Build:* the Q-centred three-term Gram and border in `KhatriRao`, the GLS level and mass (decision 4), and `WorkingResponse.ResponseMean` and `ResponseSquares` under Q.
   - *Gate:* the centred block equals an explicit projection off the whitened constant column, to rounding, for each intercept convention. The level matches `level_one_kronecker` on an R without equal column sums, and equals the grand mean on one with them.
3. **Tempering and evidence.**
   - *Build:* nothing new if slice 1 routed the constant. Check the ladder and the evidence estimators read it.
   - *Gate:* at b ∈ {0, 0.3, 1}, the tempered marginal under Q minus the whitened one is (b/2)·n·ln\|Q\|. The evidence estimators agree with each other under Q as they do under Q = I (`evidence-estimator-disagreement-brief.md`'s tolerances).
4. **Case-level weights.**
   - *Build:* weights constant within each case keep the Kronecker form; per-location weights stay refused.
   - *Gate:* the whitening oracle with L = diag(w)^{−½} ⊗ L_R.
5. **The readers downstream.**
   - *Build:* `BassPointwise`'s block density under σ²R; predictive noise σ·L_R·z; F11 on whitened and raw residuals (decision 6); the emulator's σ²_e·R for calibration (primitive 2); B1's table as an `IInnerProduct` (decision 7).
   - *Gate:*
     - a fit to data simulated with residual covariance σ²(I ⊗ R) covers at 90% where the Q = I fit under-covers;
     - the whitened F11 reads about −1/m;
     - the pilot's Stage 1 rerun under B1's Q, read against its R1 and R2 with the gates pre-stated before the run.

## Gates common to every slice

- **The invariant.** With Q = I (`IdentityInnerProduct`, `IsIdentity`), every existing result stays bit for bit. The identity takes today's arithmetic path in every place the table lists. The pinned tier and the slice's `scripts/test-select.ps1 -Scope Commit` run gate it, and the slice's end runs `-Scope Slice`.
- **The oracle** is fixed before its first reading, with its tolerance, its seeds and its rule for a marginal failure.
- **Coverage claims** register in their host's `ComparativeSweep` and read across seeds.

## Not in this brief

- Sampling R's parameters.
- Q with per-location weights.
- Q with several outputs (E34).
- The PCA mode, whose truncation residual already keeps the within-output dependence (D07).
- Primitive 1's residual-set emulator, which needs only `TruncationResiduals()` for direct fits. It is a separate, smaller change.

## Working rules

Those of the BASS slice briefs:
- commits on `main` by pathspec, nothing staged;
- builds below normal priority with capped threads;
- `scripts/test-select.ps1` with an `-ArtifactsRoot` set inline;
- a fresh-eyes review before each report;
- each report ends in an Exposition section the user can splice into docs later; it carries no pointers into `docs/`, and `docs/` is not opened.

## Report

*Not started.*
