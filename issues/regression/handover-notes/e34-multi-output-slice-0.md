# E34 slice 0: the matrix-t references (chip, 2026-11-08)

Brief: `issues/regression/bass/multi-output-brief.md` (report under "Reports"). Landed on `main`, nothing pushed: `d5b5dc6` (references + gates), `b6d5413` (report). Whole regression assembly 2,573/2,573 on 15eb151 + slice. `artifacts/chip-e34/` deleted. Slice 1 waits for the main session's verification note.

For whoever builds slice 1:
- **References live test-side:** `tests/regression/MatrixTReferences.cs`. `MultiOutputCase(Y, X, PriorPrecision Λ₀, PriorMean A₀?, b, FlatIntercept)`; `MatrixTReference.Collapse` → `MatrixTCollapse(R, LogDeterminantTerm, Count = nb − f, …)`; `LogLikelihoodGiven(collapse, Σ)`, `LogMarginal(collapse, ReferenceSigmaPrior)`, `SigmaPosterior`, `MatrixTLogDensity`; `DenseVecReference.LogLikelihoodGiven(case, Σ, KroneckerOrder)`; `ReferenceInverseWishart` (LogDensity, Sample, SamplingRoot/SampleFromRoot); `SigmaQuadrature.Adaptive` (K = 2); `ReferenceSigmaPrior.ScaledToResponse(Y, levels, c)` = decision 2.
- **Test helpers** in `MatrixTReferenceTests`: `Case(y, xc, centre, form, pencilPrior, τ, b)` builds the multi-output case for scalar BASS's three prior forms (g-form Λ₀ = τ·blockdiag(n, G_c); normalized ridge τ·blockdiag(1, diag G_c/n); pencil τD) under Mean/Zero/Flat; `ScalarLogMarginal` reads the engine's scalar block from the centred LS fit.
- **Vec order:** vec stacks A's columns (outputs): Cov(vec A) = Σ ⊗ Λ₀⁻¹.
- **Count under level centring is nb** (intercept inside the prior, as scalar BASS), not b(n − 1); flat is nb − 1 and fails decision 3's condition at nb ≤ K − ν₀ even under a proper prior. S_Y for decision 2: divisor n, about the centre levels (0 under the zero centre). Sent to the main session as items 1–3 of the report; not yet decided.
- **Dense vec route's limit:** its Cholesky works at κ̂ = cond(Σ)·(1 + b·tr XΛ₀⁻¹Xᵀ); past ~10¹³ it reads nothing. The collapse (p × p and K × K) stays exact there.
- **Quadrature lesson:** a box over Σ's Cholesky factor misses the tails (the off-diagonal spreads with L₁₁); Bartlett coordinates of Σ⁻¹ on a moment-matched reference make them independent χ and normal.
- Sources: Brown–Vannucci–Fearn 1998 and Huang–Wand 2013 not on disk; asked the main session. Matrix-t cited from memory (Dickey 1967; Gupta and Nagar 2000).

## Addendum (2026-11-08): the sources read — `b55a827` (gates), `a274708` (docs)
- Papers: `supellex-private\PDFs\Bayesian\jrssb98.pdf`, `HuangWand2013.pdf`. BVF's PDF drops Greek in both pdftotext and pypdf (Type1 math fonts, no ToUnicode). Decoder: dump raw codes per font (`rawtext.py`), map F10/F12 = TeX OML (14 δ, 13 γ, 11 α…), F11 = OT1 capitals (6 Σ, 255 Γ), F7 = OMS (255 −, 24 ∼). pypdf lives in pyenv 3.12.6 only (`C:\Users\azrie\PDenv\PyenvPython\pyenv-win\versions\3.12.6\python.exe -I script.py`); scripts in this session's scratchpad.
- BVF: Σ ∼ IW(δ; Q), Dawid shape δ = ν − q + 1; α ∼ N(α₀, hΣ), h → ∞; (20) exponent (n + δ + q − 1)/2 → count n. Gate: ln p − ln g(γ) constant to 2.8e-14 over 2 × 32 structures; Q_γ = Ψ₀ + R.
- HW: Σ | a ∼ IW(ν + K − 1, 2ν diag(1/a)), a_k ∼ IG(½, 1/A_k²); ρ ∝ (1 − ρ²)^(ν/2 − 1). At K = 1 = `Sigma2Prior.Hierarchical(ν/2, ½)`, A = s. Decision 2's default = HW at ν = 2c → correlations ∝ (1 − ρ²)^(c − 1) (piled at ±1); raised to the main session.
- Reference IW sampler now inverts the Bartlett factor (W was singular in doubles at ν = 0.5).
