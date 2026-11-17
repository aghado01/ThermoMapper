# 4. What speaker geometry must a predictive encoder preserve?

Proposed follow-up to [profile identification](03-linked-context-profiles-brief.md) and [predictive sufficiency in audio](../predictive-information-in-audio-brief.md). A returning speaker provides a test of representation quality when an earlier snippet helps predict an acoustic measurement on a different turn.

**Experiment.** In the balanced two-source model, let \(X\) be the available reference snippet and \(Z=f(X)\) a frozen code. A fresh same-source segment has measured response \(T\) and held-out context \(V\). Draw \(V\) independently of \(S,X\) with weights \(w_v\). Assume \(\mathbb E[T\mid S,X,V=v]=\mu_v+S a_v\). This explicitly excludes persistent residual channel or content effects after conditioning. The predictor sees \(Z,V\). Define \(m_X=\mathbb E[S\mid X]\), \(m_Z=\mathbb E[S\mid Z]\), and \(A=\sum_v w_v a_v^2\). Begin on finite spaces and restrict statements about conditional laws to positive-probability inputs.

**The objective.** Prove that the best full-snippet predictor is \(\mu_v+a_v m_X\), and the best code-based predictor is \(\mu_v+a_v m_Z\). For any heads \(g_v\), the excess squared loss above the full-snippet Bayes risk should be exactly

\[
A\,\mathbb E[(m_X-m_Z)^2]
+\mathbb E[(g_V(Z)-\mu_V-a_Vm_Z)^2].
\]

For \(A>0\), excess \(\eta\) therefore limits the loss of source information in posterior-mean squared error to \(\eta/A\). If all sampled contexts have zero contrast, construct an uninformative optimal encoder despite an informative reference snippet. Instantiate both cases with actual harmonic blocks. For an ambiguous reference, explicitly allow \(m_X\) strictly between \(-1\) and \(1\).

**Geometry and discrimination.** Give codes the pseudometric obtained from the weighted RMS distance between their *optimal decoded profiles*. Prove it equals \(\sqrt A\,|m_Z(z)-m_Z(z')|\). In this finite binary model the predictive profiles lie on a line segment between the two source profiles. Definite evidence for opposite sources reaches opposite ends when the reference identifies the source. Uncertain snippets give interior points. Use held-out responses to bound the error of a distance computed with fitted heads in mean square over two draws, accounting for both endpoints.

For two independently drawn snippets and independent balanced source draws, the posterior probability of a source match is \((1+m_Xm_{X'})/2\). Derive the corresponding code-based probability. With \(K=\mathbb E[m_X^2]\) and \(L=\mathbb E[m_Z^2]\), the target excess Brier risk for a match is \((K^2-L^2)/4\le \eta/(2A)\) when \(A>0\). This is a verification trial with its own prior on a match. The same-source training pairs obey a different sampling law. Apply the bound to a small reference distribution in which more than one pitch or context occurs per speaker.

**How much listening supplies evidence?** As a follow-up corollary, let the reference contain \(n\) conditionally independent draws from the same speaker, with finite one-draw observation laws \(P_+,P_-\). Include any observed context in the draw. Define their overlap \(\rho=\sum_x\sqrt{P_+(x)P_-(x)}\). Compare the optimal full-reference predictor with an oracle told the speaker. Target a squared-loss gap of at most \(A\rho^n\). Seek a sharp waveform example in which each draw is either shared by the two voices or reveals its source. This counts independently generated segments. Translating it to milliseconds of continuous speech would require a temporal dependence model.

**Which raw shapes can training determine?** Give two codes with identical held-out risk and decoded distances but arbitrarily different Euclidean distances, by jointly changing encoder scale and decoder slope. Also exhibit a nuisance coordinate ignored by the decoder that separates same-source snippets arbitrarily far. If seeking a Euclidean margin, prove a conditional statement using a declared Lipschitz bound on the decoder and uniform prediction error on the tested support. Report when low *mean* risk can still leave an individual pair ambiguous.

**Learning consequence.** The first theorem applies to every encoder attaining a specified population excess. To claim a finite recording budget, fix a constrained encoder/head class and attainable approximation error, bound estimation and optimization errors, and count the recordings used for teacher or pair construction. Name the fixed observation and code budget when comparing a wavelet frontend to a learned frontend. [Related evidence contracts](../deep-teaching/evidence-for-memory-depth-brief.md)

**Closure.** Supply the decomposition, a physical example with informative contexts, the induced distance and a same-speaker prediction bound. Record the label-exchange symmetry and the counterexample to Euclidean shape identification. The model predicts an acoustic response conditional on source and context. Extending the identifiable profile to more speakers, continuous articulation or overlapping sources requires its own observation model.
