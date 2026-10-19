
### 3. Extrema Extraction & Downstream Readout Issues (`SplineExtrema`)

* **Discarded Boundary Clipping Flags**
  * **Source**: [`issues/spc-BARS/discussions/sol-SPC-bars-graphs-brackets.md:509`](file:///d:/aghado01/ThermoMapper/issues/spc-BARS/discussions/sol-SPC-bars-graphs-brackets.md#L509-L514), [`lean/Protolemmata/BARS/spc-bars-readout-obligations.md:124`](file:///d:/aghado01/ThermoMapper/lean/Protolemmata/BARS/spc-bars-readout-obligations.md#L124-L127).
  * **Gap**: In [SplineExtrema.cs](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/SplineExtrema.cs), `PeakSpan.LeftClipped` and `RightClipped` identify when a span was truncated by the $[0,1]$ acquisition boundary. However, [BarsEnsemble](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/BarsEnsemble.cs) discards these flags when pooling draws, preventing the outer bracket expansion loop from detecting boundary clipping.

* **Span Coverage Probability Violations ($\pi(T) > 1$)**
  * **Source**: [`docs/Maths/BayesianRegression.md:48`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L48-L49), [`lean/Protolemmata/BARS/bars-span-lemmas.md:8`](file:///d:/aghado01/ThermoMapper/lean/Protolemmata/BARS/bars-span-lemmas.md#L8-L14).
  * **Observation**: `SpanCoverage` $\pi(T)$ is defined as the posterior expectation of significant spans covering $T$. When a single draw has nested or multiple overlapping peak spans, $\pi(T)$ exceeds 1.0, breaking its interpretation as a valid probability.

* **Non-Commutation of Nonlinear Functionals on the Posterior Mean**
  * **Source**: [`lean/Protolemmata/BARS/bars-multipeak-lemmas.md:108`](file:///d:/aghado01/ThermoMapper/lean/Protolemmata/BARS/bars-multipeak-lemmas.md#L108-L135), [`lean/Protolemmata/BARS/bars-span-lemmas.md:76`](file:///d:/aghado01/ThermoMapper/lean/Protolemmata/BARS/bars-span-lemmas.md#L76-L96).
  * **Observations**:
    1. **Argmax & Expectation**: $\operatorname{argmax} \mathbb{E}[f] \ne \mathbb{E}[\operatorname{argmax} f]$ (two-tent counterexample). Peak locations extracted from the posterior mean curve `r.Fit` are systematically biased.
    2. **Peak Count Undercounting**: `count_of_mean_undercounts` proves that averaging curves fills in intermediate valleys. Extracting peaks from `r.Fit` systematically undercounts transitions. Features must be extracted draw-by-draw and pooled via point-process intensities $\lambda(T)$ or [PeakPosterior](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/SplineExtrema.cs).
    3. **FWHM of Mean $\ne$ Mean of FWHMs**: Span widths must be calculated per draw, never from `r.Fit`.

* **Boundary Maximum Policy & Degree Assumptions**
  * **Source**: [`lean/Protolemmata/BARS/bars-multipeak-lemmas.md:88`](file:///d:/aghado01/ThermoMapper/lean/Protolemmata/BARS/bars-multipeak-lemmas.md#L88-L100), [`lean/Protolemmata/BARS/spc-bars-readout-obligations.md:108`](file:///d:/aghado01/ThermoMapper/lean/Protolemmata/BARS/spc-bars-readout-obligations.md#L108-L113).
  * **Observation**: [SplineExtrema](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/SplineExtrema.cs) unconditionally counts boundary maxima ($i=0, n-1$) as peaks. In thermal curve estimation, boundary peaks are often truncation artifacts that require bracket expansion rather than physical transitions. Needs a configurable policy (`interiorOnly` vs `includeBoundary`).
  * **Degree Assumption**: Span polynomial derivative recovery was originally exact only for cubic splines ($p=3$); higher degrees require explicit Vandermonde interpolation and Durand–Kerner root finding.
