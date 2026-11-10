# **Lean 4 Formalization Audit Report**

[https\://gemini.google.com/app/98777d176100e6ff](https://gemini.google.com/app/98777d176100e6ff)

*User prompt: need to audit these latest 'memes*

Response: Across all seven files, the formalizations are proof-closed, non-vacuous, and free of the standard Lean 4 traps (truncated ℕ subtraction, Real.log on nonpositive reals, ENNReal.ofReal clamping negatives to 0, or division-by-zero vacuity). Below is the file-by-file audit covering axiom soundness, edge-case guards, and three specific tightening opportunities.

### **1\. Axiom Audit (\`\`)**

The three axioms (fixed\_constant\_contraction, gvdv07\_theorem1, gvdv07\_lemma1) are sound and accurately reflect Ghosal & van der Vaart (2007, Theorem 1 & Lemma 1; 2017, Theorems 8.19(b) & 8.20):

> * **External vs. internal covers in hentropy:** Both fixed\_constant\_contraction and gvdv07\_theorem1 allow the covering finset N : Finset (E × ℝ) to be external (not restricted to S or Θ n), whereas htests only provides tests around centers θ₁ ∈ S. This is sound because fixed\_constant\_contraction covers at radius ξ \* εbar n / 2 (rather than ξ \* εbar n) and gvdv07\_theorem1 covers at radius ξ \* r / 4 (rather than ξ \* r / 2); by hd\_comm and hd\_triangle, any ball meeting Θ n ⊆ S is contained in a double-radius ball centered at a point in Θ n ⊆ S.  
> * **Omission of d n θ θ \= 0 (hd\_self) in fixed\_constant\_contraction and gvdv07\_theorem1:** The axiom signatures require hd\_comm, hd\_triangle, and hd\_meas, which forces 0 ≤ d n θ θ' (since 0 ≤ d n θ θ ≤ 2 \* d n θ θ'), but they do not explicitly include hd\_self : ∀ n θ, d n θ θ \= 0\. A potential pathological instantiation—setting *dn*​(*θ*,*θ*′)=1 whenever *θ*\=*θ*0​ or *θ*′\=*θ*0​ and 0 otherwise (which satisfies hd\_comm and hd\_triangle)—is blocked only indirectly because hmass and hrem force Θ n ∩ B to be nonempty inside klBall, and htests then demands exponentially powerful tests separating θ₁ ∈ Θ n ∩ B from θ₀, which is impossible on klBall. Adding (hd\_self : ∀ n θ, d n θ θ \= 0\) to both axiom signatures is a zero-cost hygiene improvement (paramDist satisfies it by simp \[paramDist, empiricalNorm\]).  
> * **klBall in ℝ≥0∞:** Defining klBall with klDiv and evariance (llr ...) in ℝ≥0∞ bounded by ENNReal.ofReal (n \* ε ^ 2\) prevents infinite KL divergences or log-likelihood-ratio variances from wrapping to 0 under toReal. Moreover, hS : prior Sᶜ \= 0 and hS\_pos : ∀ θ ∈ S, 0 \< θ.2 ensure gaussianReal has strictly positive variance prior-a.e.

### **2\. Module-by-Module Audit**

> * **BGP-4 Bookkeeping (\`\`):**  
  * **Status:** Clean.  
  * **Verification:** gaussian\_kl\_scale\_part\_le (UV-1a), gaussian\_v\_scale\_part\_le (UV-1b), and scale\_ratio\_bridge (UV-1d) have exact constants on *r*\=*σ*0​/*σ*∈\[1/2,2\]. The Real.log x \= 0 trap for *x*≤0 in log\_inv\_radius\_le\_n\_radius\_sq is explicitly guarded by pairing with eventually\_priorMassRadius\_pos in exp\_neg\_n\_radius\_sq\_le\_radius and log\_inv\_mul\_radius\_le.  
> * **BGP-4 Contraction with Unknown Variance (\`\`):**  
  * **Status:** Clean.  
  * **Verification:** The sieve-complement margin constants are exact: converting from the slab radius *ρn*​\=*εn*​ to the KL radius *εn*′​\=*ρn*​/(*σ*0​/6) multiplies *n*(*εn*′​)2 by 36/*σ*02​, giving 4×36/*σ*02​\=144/*σ*02​ in posterior\_contraction\_unknown\_variance (from (C \+ 4\) n ε\_n^2) and 2×36/*σ*02​\=72/*σ*02​ in posterior\_contraction\_unknown\_variance\_gvdv07 (from gvdv07\_lemma1's *e*2*nεn*2​ factor). In exists\_tests (UV-5), the *n*\=0 edge case is handled cleanly (*φ*≡0, *e*0\=1), and the three test regimes (tests\_mean\_case, tests\_scale\_up\_case, tests\_scale\_down\_case) match the compact noise range \[*σ*min​,*σ*max​\] with *K*\=*σ*min2​/(1024*σ*max4​) and *ξ*\=1/8.  
> * **BASS Guarded Estimators (\`\`):**  
  * **Status:** Clean.  
  * **Verification:** All simple-random-sample and draw lemmas (srs\_mean, crossfit\_gap, crossfit\_bias, drawsExpect\_mean, telescope\_iid, telescope\_srs) carry explicit non-emptiness or cardinality bounds (k \+ 1 ≤ U.card, A.Nonempty, k \+ 1 ≤ Aᶜ.card, 0 \< n, J.Nonempty), guaranteeing strictly positive denominators and binomial coefficients.  
> * **BASS Local Cover, Maurey, and Jeffreys Profile (\`\`):**  
  * **Status:** Clean, with one scaling convention to watch and one optional strengthening.  
  * **Wrinkle 1 (rssOf / rSquaredOf scaling):** rssOf V y is defined as sInf ((fun v \=\> designNormSq (y \- v)) '' V), where designNormSq divides by N (making rssOf the mean squared residual *N*1​RSS, not the unnormalized RSS). Because tss in rSquaredOf V y tss \= 1 \- rssOf V y / tss is an abstract parameter with 0 \< tss, rSquaredOf\_le\_one, rSquaredOf\_mono, and jeffreys\_superset\_ge are completely valid as stated; however, any downstream file that instantiates tss on data y must pass the mean total sum of squares designNormSq (y \- ȳ) rather than ∑*k*​(*yk*​−*y*ˉ​)2.  
  * **Wrinkle 2 (Discarded internal-center property in exists\_local\_cover\_with\_scale):** exists\_cover\_iUnion\_subspace\_balls and BGP.exists\_isCover\_Icc both guarantee that the cover centers lie inside the covered set, but exists\_local\_cover\_with\_scale drops those conjuncts via obtain ⟨Nf, \-, hNfcov, hNfcard⟩ and obtain ⟨Nσ, \-, hNσcov, hNσcard⟩. Keeping (N : Set (E × ℝ)) ⊆ (⋃ k, (V k : Set E) ∩ closedBall f₀ r) ×ˢ Icc (σ₀ \- r) (σ₀ \+ r) is free if an internal cover is needed downstream.  
> * **BASS Sensitivity Integrals & Markov Generator (, ):**  
  * **Status:** Clean.  
  * **Verification:** Lifting hoeffdingOp from a commutative ring into Algebra.adjoin ℝ (Set.range F.mean) (CommutingMeans.algebraOf) and proving self-adjointness by Algebra.adjoin\_induction (form\_algebra\_adj) is a clean architecture. In , \`variance\_le\_half\_jacobiForm\` and \`integral\_expansion\` match \$\\text{Var}(f) \= \\\|f\\\|\_2^2 \- a\_0^2 \\le \\frac{1}{2} J(f, f)\$ on $[0,1]$ without needing the $1/(2n+1)$ normalization. In , the explicit 2×2 counterexample in no\_common\_eigenvector (3\\mathbb{E}\[C\] \= \\begin{psmallmatrix} 38 & 48 \\\\ 48 & 64 \\end{psmallmatrix} and 3\\mathbb{E}\[P\_j\] \= \\frac{1}{25}\\begin{psmallmatrix} 59 & 12 \\\\ 12 & 16 \\end{psmallmatrix} for *a*\=(1,0), *b*\=(3/5,4/5)) and rayleigh\_max\_is\_eigen are exact.  
> * **BASS Structure Chain (\`\`):**  
  * **Status:** Clean, with a natural unification between SC-2 and SC-5.  
  * **Stationary vs. reversible separator bound (cutFlow\_le\_separator / gap\_le\_separator):** SC-2 is proved both for reversible kernels (gap\_le\_cut) and for non-reversible stationary kernels (cutFlow\_compl\_of\_stationary, gap\_le\_cut\_of\_stationary) because a block-scan sweep is stationary without being reversible. However, SC-5 (cutFlow\_le\_separator and gap\_le\_separator) reverts to assuming reversibility (hrev : ∀ s t, π s \* P s t \= π t \* P t s) in order to rewrite ∑*s*∈*A*​∑*t*∈*W*​*π*(*s*)*P*(*s*,*t*)=∑*t*∈*W*​*π*(*t*)∑*s*∈*A*​*P*(*t*,*s*). Because cutFlow\_compl\_of\_stationary already establishes cutFlow(*π*,*P*,*A*)=cutFlow(*π*,*P*,*Ac*)=∑*t*∈*Ac*​∑*s*∈*A*​*π*(*t*)*P*(*t*,*s*) for any stationary kernel, replacing the exit-separator condition hA with the entry-separator condition ∀ t ∉ A, t ∉ W → ∀ s ∈ A, P t s \= 0 (which holds whenever *W* separates *A* and *Ac*∖*W*) yields cutFlow π P A ≤ e \* ∑ t ∈ W, π t and gap\_le\_separator for non-reversible stationary sweeps directly from cutFlow\_compl\_of\_stationary, without hrev.

### **3\. Pointed Action Checklist**

> 1. **\`\` (fixed\_constant\_contraction, \`gvdv07\_theorem1\`):** Optionally add \`(hd\_self : ∀ n θ, d n θ θ \= 0)\` alongside \`hd\_comm\` and \`hd\_triangle\` so \`d n\` is formally a semimetric on its own.  
> 2. **\`\` (rssOf / \`rSquaredOf\`):** Add a docstring note on \`rSquaredOf\` recording that \`tss\` is scaled by 1/*N* (matching \`designNormSq\` in \`rssOf\`).  
> 3. **\`\` (cutFlow\_le\_separator, \`gap\_le\_separator\`):** If SC-5 is intended to apply to the non-reversible block-scan sweeps motivated in SC-2, add \`cutFlow\_le\_separator\_of\_stationary\` and \`gap\_le\_separator\_of\_stationary\` using \`cutFlow\_compl\_of\_stationary\` in place of \`hrev\`.

---

*User prompt: Some more*

Response: All five files are axiom-free, proof-closed, and well-guarded against Lean 4's usual formalization traps. Both sub-stacks—the relative persistent homology trilogy (RelativePersistence → RelativeHomology → IntervalDecomposition) and the Fortuin–Kasteleyn / Potts coupling pair (RandomCluster → EdwardsSokal)—include concrete sanity checks (Fin 3 vertex-edge filtration and Fin 2 single-edge partition/agreement formulas) that rule out vacuity or off-by-one/sign errors in the definitions.

### **1\. Relative Persistent Homology Trilogy (\[cite: 8, 10, 11\])**

> * **Matrix Reduction & Subcomplex Restriction (RelativePersistence.lean):**  
  * **Status:** Clean.  
  * **Verification:** isLow\_iff\_of\_reductions proves via well-founded induction on column index *j* and combination\_top that IsLow R i j is an invariant of the boundary matrix *D*, independent of the upper-triangular *U* in *D*\=*RU*. relative\_isReduction and relative\_barcode then establish that deleting the rows and columns of an initial subcomplex *K* (rel K D) requires no second reduction: finite bars outside *K* are preserved (relative\_isLow\_iff), and relative essential classes outside *K* are exactly absolute essential classes outside *K* plus cells outside *K* that kill a class born in *K* (relative\_essential\_iff).  
> * **Homological Bridge & Simplicial Instantiation (RelativeHomology.lean):**  
  * **Status:** Clean.  
  * **Verification:**  
    * **No ℕ-truncation at *p*\=0:** relCycles K S p is defined as chains F (C.cells S p) ⊓ (chains F K).comap C.boundary.mulVecLin rather than intersecting with C.cells K (p \- 1). Because C.boundary\_deg enforces deg*j*\=deg*i*\+1 whenever ∂*ij*​\=0, the boundary of a *p*\-chain is automatically supported in degree *p*−1 (and vanishes at *p*\=0), avoiding 0 \- 1 \= 0 bugs on ℕ.  
    * **Guarded degree subtraction in SimplicialFiltration.toComplex:** Setting deg := fun i \=\> (X.cell i).card \- 1 is guarded by X.nonempty : ∀ i, (cell i).Nonempty, which ensures (X.cell i).card ≥ 1 before subtracting 1\.  
    * **Quotient honesty & algorithm independence:** relBoundaries\_le\_relCycles proves *Bp*​(*S*,*K*)≤*Zp*​(*S*,*K*) whenever *K*⊆*S* are initial segments, ensuring the .comap (C.relCycles K S p).subtype in RelHomology does not silently clip boundaries outside *Zp*​(*S*,*K*). Although the internal dimension count (finrank\_persistence\_of\_graded) uses a homogeneous reduction *DV*\=*R* (GradedReduction), exists\_gradedReduction constructs one for every FilteredComplex and finrank\_persistence transfers the result to *any* IsReduction C.boundary R U via isLow\_iff\_of\_reductions.  
> * **Interval Module Decomposition & Uniqueness (IntervalDecomposition.lean):**  
  * **Status:** Clean, with one structural subtlety in barcode\_unique worth noting for downstream use.  
  * **Verification:** rep R V i chooses the killing boundary column *R*⋅,*j*​ when *i* is killed by *j* (IsLow R i j), and the cycle column *V*⋅,*i*​ when *i* is never killed—correctly respecting the elder rule so a dying bar's basis element maps to 0 in *Hp*​(*T*,*K*) as soon as *j*∈*T* (persistence\_barBasis).  
  * **Subtlety (barcode\_unique / interval\_count\_eq and invisible intervals):** barcode\_unique proves that any candidate interval family (b : B → ι, d : B → Option ι) matching the persistence ranks has the exact table multiplicities for every forward finite bar (b x \= i, d x \= some j) with ¬*K*(*i*) and *i*\<*j*, and every essential bar (b x \= i, d x \= none) with ¬*K*(*i*). Remarkably, the proof does not even need a well-formedness hypothesis on (b, d). However, because intervalRank b d S T is only tested on initial segments *K*⊆*S*⊆*T*, any degenerate element *x*∈*B* with d x \= some j and *j*≤*b*(*x*) (or with *K*(*b*(*x*)) and *d*(*x*)=some *j*∈*K*) is identically invisible to intervalRank b d S T (since *S*(*b*(*x*))⟹*S*(*j*)⟹*T*(*j*) whenever *j*≤*b*(*x*) and *S* is initial). If a downstream theorem needs a full bijection B ≃ {i // IsBar C R K p i}, pair barcode\_unique with the well-formedness hypothesis ∀ x, ¬ K (b x) ∧ ∀ j, d x \= some j → b x \< j.

### **2\. Random-Cluster & Edwards–Sokal Coupling (\[cite: 9, 12\])**

> * **Random-Cluster Measure & Holley Comparison (RandomCluster.lean):**  
  * **Status:** Clean.  
  * **Verification:**  
    * **Subspace-dimension proof of Grimmett (3.11):** Proving supermodularity of the cluster count *k*(*a*)+*k*(*b*)≤*k*(*a*⊓*b*)+*k*(*a*⊔*b*) (componentCount\_supermodular) by identifying *k*(*G*)=finrankR​(edgeConstant *G*) and applying Submodule.finrank\_sup\_add\_finrank\_inf\_eq with edgeConstant *a*⊔edgeConstant *b*≤edgeConstant(*a*⊓*b*) and edgeConstant *a*⊓edgeConstant *b*≤edgeConstant(*a*⊔*b*) is exceptionally clean—it completely bypasses edge-induction and component-merging casework.  
    * **Positivity & boundary cases:** Because IsEdgeParameter enforces 0≤*pij*​\<1, rcWeight\_bot\_pos guarantees *Z*\>0 (rcPartition\_pos) for all *q*\>0 even when *pe*​\=0 on some or all edges, and Mathlib's holley applies directly without requiring strictly positive weights on all configurations.  
    * **Minor interface note on IsEdgeParameter:** IsEdgeParameter p quantifies nonneg and lt\_one over *all* pairs i j : V, including the diagonal i \= i (even though rcWeight filters out e.IsDiag). For pe J β, pe\_isEdgeParameter satisfies this automatically whenever *Jii*​≥0 and *β*≥0 (in particular for *Jii*​\=0, where *pii*​\=0\<1), so this never obstructs the FK parameterization.  
> * **Potts Model & Two-Point Identity (EdwardsSokal.lean):**  
  * **Status:** Clean.  
  * **Verification:**  
    * **Bijection between subsets of edgePairs V and SimpleGraph V:** prod\_add\_eq\_sum\_config uses SimpleGraph.fromEdgeSet to turn the product over unordered distinct pairs into a sum over Config V, properly discharging looplessness via not\_isDiag\_of\_mem\_edgeSet.  
    * **Cluster-spin counting:** spinsEquiv establishes (ω.ConnectedComponent → Fin n) ≃ {σ : V → Fin n // ω ≤ agreeGraph σ}, giving *nk*(*ω*) compatible spin assignments (card\_compat) and *n*1​*nk*(*ω*) agreeing assignments when ¬*ω*.Reachable *ij* (card\_compat\_agree\_mul), from which pottsPartition\_eq and pottsAgree\_eq follow cleanly.  
    * **External cross-check:** pottsAgree\_one\_edge combines pottsAgree\_eq with rcProb\_one\_edge on Fin 2 to recover P(*σ*0​\=*σ*1​)=*p*01​\+*n*(1−*p*01​)1​, confirming that the edge factors, cluster weights, and normalization constants match across both files.

---

*User prompt: one more fresh one*

Response: This file (Conformal) is axiom-free, proof-closed, and free of formalization traps.

> * **Division-by-zero & ℕ-subtraction guards:** Every coverage theorem takes ht : t ∈ C (or hC : ∀ j, j ∈ C ↔ str j \= str t in coverage\_stratum, which discharges t ∈ C via (hC t).2 rfl), guaranteeing 0 \< C.card, while card\_permsOn\_pos guarantees 0 \< (permsOn C).card via 1 ∈ permsOn C. All ℕ subtractions (C.card \- k in le\_card\_rankLt\_lt, r \- 1 in card\_rank\_eq) are guarded by hk : k ≤ C.card and hr : 1 ≤ r.  
> * **Order-duality & tie-breaking (CB-1):** card\_rankLe\_le proves (C.filter fun u \=\> rankLe C s u ≤ k).card ≤ k cleanly via exists\_max\_image, and le\_card\_rankLt\_lt dualizes it on αᵒᵈ at C.card \- k without redoing the combinatorics. coverage\_tieBreak uses Prod.Lex (α ×ₗ L) to lift any equivariant score with injective labels into coverage\_orderings\_eq.  
> * **Split conformal & CQR interval property (CB-2):** cqrScore negates the Romano–Patterson–Candès nonconformity score into min(*y*−*q*lo​(*x*),*q*hi​(*x*)−*y*), and cqr\_ordConnected establishes that the accepted candidate set is order-connected (Set.OrdConnected, i.e., an interval) without even needing *q*lo​(*x*)≤*q*hi​(*x*).  
> * **Latent vs. integrated vs. mixed scores (CB-4):** integratedScore\_le\_latentScore proves that the latent-posterior score dominates the integrated score via the nonnegativity of the conditional variance ∑*v*​*qθ*​(*v*)(*fθ*​(*v*,*zi*​)−*gθ*​(*zi*​))2≥0, and mixed\_counterexample provides a concrete 2-unit Bernoulli witness where the asymmetric mixed score has 0% coverage (¬ Accepts 1\) on every dataset *z* (since the data unit scores 5/8 and the candidate scores 1/2).  
> * **Error margins & fragile bands (CB-3):** Even though Fragile, AcceptsInflated, and FragileEstimated do not carry a separate 0 ≤ ε j hypothesis, hclose : ∀ j ∈ C, \|ŝ j \- s j\| ≤ ε j forces 0 ≤ ε j on C automatically via abs\_nonneg.

**Two Minor Interface Notes in Margin (CB-3)**

> * **Unit-attached ε in coverage\_close vs. slot-attached ε in coverage\_inflated:** In coverage\_inflated, the error bound ε : Equiv.Perm ι → ι → ℝ can depend on both the permutation σ and the slot j (for example, when the test slot t has a larger approximation error than calibration slots j ≠ t). In coverage\_close, by contrast, ε : ι → ℝ is attached to the underlying data unit (\|ŝ σ j \- S (z ∘ σ) j\| ≤ ε (σ j)), so that Fragile k C (S z) ε (σ t) is a fixed property of the unit σ t and averages out to \|M\| / \|C\| via ratio\_permsOn. When applying coverage\_close to an asymmetric algorithm where the error depends on whether a unit lands in the test slot t or a calibration slot, ε u must be chosen as the envelope across both roles (or a constant bound fun \_ \=\> ε\_max).  
> * **Domain of hε in fragileEstimated\_mono:** fragileEstimated\_mono currently asks for (hε : ∀ j, ε j ≤ ε' j) on all of ι, whereas hclose in fragileEstimated\_of\_fragile and card\_fragile\_le\_estimated is restricted to j ∈ C; weakening hε to (hu : u ∈ C) (hε : ∀ j ∈ C, ε j ≤ ε' j) matches the rest of the section when C is a proper stratum.

---

