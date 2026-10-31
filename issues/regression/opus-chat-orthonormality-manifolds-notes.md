The question that sorts most of this is whether whatever uses the frame is invariant to rotations within its span. If it is, only the subspace matters and you're on the Grassmannian. If it isn't, the frame itself is the object and you're on the Stiefel manifold, up to signs. BASS has both cases, and not always where you'd first expect.

**BASS-PCA lives on the Stiefel manifold.** The reconstruction ΦΦᵀ depends only on the span. But `bassPCA` fits a separate BASS emulator per component, and axis-aligned hinge fits don't transform with a rotation Φ → ΦQ, w → Qᵀw. So the emulators depend on the frame. Two consequences:
- *Eigengaps decide what's reportable.* Davis–Kahan bounds a subspace's rotation by the perturbation over the eigengap. The gap after component r governs the span; the gaps inside the retained block govern the frame. Where retained eigenvalues are close, components rotate into each other between refits. Per-component readouts, such as one component's Sobol indices, mean something only across a clear gap.
- *Functional Sobol indices depend only on the span.* With orthonormal Φ, total and partial output variances are the traces of Cov(w) and Cov(E[w | x_u]), and a rotation leaves both unchanged. That needs Φ orthonormal in the L²(t) inner product. On a non-uniform grid, OQ16's √Δt weighting therefore applies to the SVD itself, not only to the robust centre. With R's unweighted SVD the cross terms don't vanish.

**E6 might not need a new model.** MoMPCA's block PCAs are already subsample replicates of the subspace. Their spread around the median on Gr(r, m), measured in principal angles, is a resampling estimate of basis uncertainty: E18's resampling idea applied to the basis, alongside Sartini et al.'s Bayesian FPCA. Two cautions:
- Each block sees n/K curves, so the angles need shrinking by about √K to represent the full sample, as in subsampling.
- To carry basis uncertainty into the emulators, align each block's frame to a reference by orthogonal Procrustes first. Otherwise component identities get scrambled across draws.

**Active subspaces (E8) are genuinely Grassmannian.** Each posterior draw of f gives C = E[∇f∇fᵀ] and its top-k eigenspace, so the posterior is a distribution on Gr(k, d). A cheap summary is to average the projection matrices WWᵀ over draws:
- The top-k eigenvectors of the average are the extrinsic mean.
- Its eigenvalues lie in [0, 1] and sum to k. All k near 1 means the draws agree on the subspace.

This stops holding one step later. The brief's next step, learning projection directions during the fit, feeds Wᵀx into axis-aligned hinges. W then becomes a Stiefel parameter, and the sampler needs moves on that manifold, such as Cayley or geodesic updates.

**KHAOS's orthonormality is a different kind.** Its polynomials are orthonormal in L²(μ) over the inputs. The design Gram isn't the identity, and there's no frame in Rⁿ to put geometry on. KHAOS meets the PCA side only through the Sobol readout: per-component PCE coefficients give per-component variance shares, and the trace identity above sums them into functional indices.
