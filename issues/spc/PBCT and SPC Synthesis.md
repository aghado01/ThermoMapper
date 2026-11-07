# **PBCT and SPC Synthesis**

[https\://gemini.google.com/app/769b5f8f5313ecb5](https://gemini.google.com/app/769b5f8f5313ecb5)

*User prompt: i wonder if there's som ekind of synthesis between this Approximate learning of parsimonious Bayesian context trees \| Statistics and Computing \| Springer Nature Link and superparamagnetic clustering.*

Response: A synthesis between **Parsimonious Bayesian Context Trees (PBCT)** (Ghani, Heard, & Sanna Passino, *Statistics and Computing*, 2026\) and **Superparamagnetic Clustering (SPC)** exists and is mathematically deep, provided you avoid a naive pairwise mapping that breaks sufficient-statistic additivity.

At its core, PBCT solves a discrete variable-order Markov partitioning problem on the probability simplex, while SPC solves a multi-scale graph partitioning problem via the disordered Potts model. Combining them bridges **Bayesian non-parametric sequence modeling** (Dirichlet-Multinomial Polya urns \+ Chinese Restaurant Process priors), **information geometry** (Fisher-Rao / Jensen-Shannon metrics on predictive distributions), and **computational mechanics** (Crutchfield–Shalizi causal states / *ϵ*\-machines).

### **1\. The Mathematical Isomorphism: PBCT as a *T*→0 Greedy Quench**

In a standard depth-*D* Bayesian Context Tree over a vocabulary *V* (∣*V*∣=*M*), every internal node branches into *M* children, scaling as *MD*. PBCT compresses this by requiring the children of any context node *e* to form a **partition** Π*e*​\={*Ce*,1​,…,*Ce*,*K*​} of *V*. Symbols *u*,*u*′∈*Ce*,*k*​ in the same cluster share a single next-token predictive distribution *ϕe*,*k*​∼Dirichlet(*η*) on the simplex Δ*M*−1.

At node *e*, each candidate preceding symbol *u*∈*V* has an empirical transition count vector x*e*⋅*u*​\=(*Xe*⋅*u*,*v*​)*v*∈*V*​ with total count *Ne*⋅*u*​\=∑*v*∈*V*​*Xe*⋅*u*,*v*​ and empirical distribution *p*^​*e*⋅*u*​\=x*e*⋅*u*​/*Ne*⋅*u*​. With a CRP(*α*) prior on Π*e*​ and conjugate Dirichlet(*η*) prior on *ϕ*, Ghani et al. score the pairwise merge of clusters *Ci*​ and *Cj*​ using the exact posterior Bayes factor:

*sij*​\=*p*(Π∣*α*)*p*(Π(*i*∪*j*)∣*α*)​⋅*B*(x*i*​\+*η*)*B*(x*j*​\+*η*)*B*(x*i*​\+x*j*​\+*η*)*B*(*η*)​

where *B*(⋅) is the multivariate Beta function. Expanding log*sij*​ via Stirling's approximation for moderate-to-large counts *Ni*​,*Nj*​ reveals the hidden statistical mechanics:

log*sij*​≈−Internal Energy Δ*Eij*​ (Information Loss)(*Ni*​\+*Nj*​)*D*JS(*πi*​,*πj*​)​(*p*^​*i*​∥*p*^​*j*​)​​\+Occam / Complexity Entropy Gain Δ*Sij*​2*M*−1​log(*Ni*​\+*Nj*​*Ni*​*Nj*​​)−log*α*\+logΓ(∣*Ci*​∣)Γ(∣*Cj*​∣)Γ(∣*Ci*​∣+∣*Cj*​∣)​​​

where *D*JS(*πi*​,*πj*​)​ is the sample-weighted Jensen-Shannon divergence (whose local quadratic form is the **Fisher-Rao metric** on Δ*M*−1).

**The physical translation:**

> * **PBCT's inference algorithm** (Heller–Ghahramani / Heard model-based agglomerative clustering) starts at *K*\=*M* singletons and greedily merges the pair with maximum log*sij*​ down to *K*\=1. In Potts model terms, this is a **deterministic, zero-temperature (*T*→0) greedy quench** along a single agglomerative trajectory.  
> * **Tree pruning as the ferromagnetic phase:** In PBCT, if the optimal partition at node *e* collapses to the single cluster Π*e*∗​\={*V*}, the context carries no predictive information and recursion terminates (*e* becomes a leaf). In SPC, Π*e*​\={*V*} is the **low-temperature ferromagnetic phase** (*T*\<*Tc*1​), where all spins lock together.  
> * **Full branching as the paramagnetic phase:** The singleton partition Π*e*​\={{*v*1​},…,{*vM*​}} corresponds to the **high-temperature paramagnetic phase** (*T*\>*Tc*2​), where all Fortuin-Kasteleyn (FK) bonds break.  
> * **Parsimonious context clusters as the superparamagnetic regime:** Non-trivial partitions Π*e*​ (1\<*K*\<*M*) correspond directly to the **superparamagnetic grains** (*Tc*1​\<*T*\<*Tc*2​).

### **2\. Calling BS: Two Traps Where a Naive Synthesis Fails**

Before wiring an SPC engine to a PBCT likelihood, two structural incompatibilities must be resolved:

**Trap 1: Pairwise Potts Hamiltonians vs. Non-Additive Coalition Sufficiency** The Potts Hamiltonian *H*(s)=−∑⟨*u*,*v*⟩​*Juv*​*δsu*​,*sv*​​ is strictly a **2-body interaction**. By contrast, the PBCT log marginal likelihood of a cluster *Ck*​ depends on the *pooled* count vector X*Ck*​​\=∑*u*∈*Ck*​​x*e*⋅*u*​ inside log*B*(X*Ck*​​\+*η*). You cannot decompose log*B*(∑*u*∈*Ck*​​x*u*​\+*η*) into a sum of static pairwise couplings ∑*u*,*v*∈*Ck*​​*Juv*​. Worse, if a symbol *u*∈*V* rarely occurs in context *e* (*Ne*⋅*u*​∈{0,1}), its Dirichlet posterior is nearly flat, yielding *suv*​\>1 against *multiple incompatible high-count clusters*. If you naively set *Juv*​∝max(0,log*suv*​), rare symbols act as **percolation bridges** that short-circuit the FK bond graph and artificially merge distinct clusters.

**Trap 2: Top-Down Greedy Myopia on Higher-Order Parity Processes** Ghani et al. infer the tree **top-down** from the root, terminating whenever a node's optimal child partition is {*V*}. If a data stream is generated by a higher-order parity or XOR Markov process (e.g., *Xt*​\=*Xt*−1​⊕*Xt*−2​(mod2) plus noise), the 1st-order transition marginals at depth 1 are completely uniform (*p*^​*u*​\=*p*^​*v*​ for all *u*,*v*∈*V*). PBCT's top-down greedy clustering merges all symbols into {*V*} at depth 1 and **terminates immediately at the root**, completely missing the depth-2 structure.

### **3\. Three Rigorous Architectures for Synthesis**

#### **A. SPC as the Global Proposal Path \+ PBCT Marginal Likelihood as the Plateau Selector**

Both algorithms suffer from complementary weaknesses:

> * PBCT has an exact, parameter-free Bayesian objective LPBCT​(Π*e*​)=log*p*(Π*e*​∣*α*)+∑*k*​log*p*(X*Ck*​​∣*η*), but searches the *BM*​ Bell-number partition space via a myopic *O*(*M*3) greedy merge that cannot backtrack.  
> * SPC explores collective multi-way spin coalitions via Swendsen-Wang / Wolff cluster flips across temperature *T*, reducing an exponential partition space to a 1D thermodynamic path {Π(*T*)}*T*∈\[*T*min​,*T*max​\]​, but relies on heuristic susceptibility thresholds (*χ*(*T*) peaks or correlation cutoffs *G*ˉ*uv*​(*T*)\>0.5) to pick the final partition.

**The synthesis:**

> 1. Map each symbol *u*∈*V* at node *e* to its posterior predictive mean *p*\~​*e*⋅*u*​\=*Ne*⋅*u*​\+*η*0​x*e*⋅*u*​\+*η*​∈Δ*M*−1.  
> 2. Define pairwise couplings using the **sample-weighted Fisher-Rao distance** (or posterior Bhattacharyya affinity) gated by effective sample size *N*eff​(*u*,*v*)=*Ne*⋅*u*​\+*Ne*⋅*v*​*Ne*⋅*u*​*Ne*⋅*v*​​ to eliminate rare-symbol percolation bridges:  
>    *Juv*(*e*)​\=*K*^1​exp(−2*a*2*N*eff​(*u*,*v*)*d*FR2​(*p*\~​*e*⋅*u*​,*p*\~​*e*⋅*v*​)​)  
>    (Isolate symbols with *Ne*⋅*u*​\<*N*min​ during the Potts sweep and assign them post-hoc via maximum posterior predictive responsibility.)  
> 3. Run the SPC temperature sweep to extract the candidate partitions Π(*T*) along the susceptibility curve *χ*(*T*), and score each candidate partition using the **exact multi-way PBCT log marginal posterior** LPBCT​(Π(*T*)).  
> 4. **Result:** If max*T*​LPBCT​(Π(*T*))=LPBCT​({*V*}) (the ferromagnetic *T*→0 state wins), prune the branch into a leaf. Otherwise, cut the tree at the superparamagnetic plateau *T*∗ that maximizes LPBCT​(Π(*T*)).

#### **B. Bottom-Up "Causal State" Percolation (Fixing Top-Down Myopia)**

Instead of running SPC inside each node top-down, run SPC **once globally** over all observed active contexts up to maximum depth *D*, E*D*​\={*e*∈⋃*d*\=1*D*​*Vd*:*Ne*​≥*N*min​}:

> * Embed every observed suffix *e*∈E*D*​ on the predictive simplex Δ*M*−1 via *p*\~​*e*​.  
> * Construct a two-layer coupling graph:  
  1. **Vertical (Tree Hierarchy) Couplings:** Connect each child context *u*⋅*e* to its parent suffix *e* with coupling *Ju*⋅*e*,*e*vert​ proportional to their Dirichlet-Multinomial Bayes factor.  
  2. **Horizontal (Sibling Equivalence) Couplings:** Connect sibling contexts *u*⋅*e* and *v*⋅*e* (or general suffixes *e*,*e*′) with Fisher-Rao coupling *Ju*⋅*e*,*v*⋅*e*horiz​.  
> * Sweeping temperature *T* from high *T* (every depth-*D* suffix distinct) to low *T* (all suffixes collapsed to the root {*V*}) induces a **thermodynamic renormalization flow** over the suffix tree:  
  * If only sibling contexts *u*⋅*e* and *v*⋅*e* lock into the same FK cluster, you recover a **PBCT parsimonious node**.  
  * If child *u*⋅*e* locks onto its parent *e*, the variable-order branch **prunes dynamically from the bottom up**, solving the XOR/parity trap because depth-2 contexts are evaluated before depth-1 collapse can censor them.  
  * If non-sibling contexts *e*,*e*′ are allowed to couple horizontally, the FK clusters generalize PBCT into an **approximate Bayesian *ϵ*\-machine** (partitioning histories into minimal sufficient causal states).

#### **C. Swendsen-Wang Cuts (Barbu–Zhu) for Full Posterior Tree Sampling**

If you want posterior uncertainty over context partitions rather than a single MAP tree, standard Gibbs sampling on CRP partitions mixes poorly because moving one symbol *u*∈*V* at a time across high-count clusters crosses deep likelihood valleys. You can turn SPC into an **exact MCMC proposal engine** via Swendsen-Wang Cuts (Barbu & Zhu, 2005):

> 1. Given current partition Π*e*(*t*)​ at node *e*, activate FK bonds between symbols *u*,*v* in the same cluster with probability *puv*​\=1−exp(−*Juv*(*e*)​/*T*).  
> 2. Pick a connected FK component *Y*⊆*Ck*​ (a collective sub-coalition of symbols) and propose splitting *Y* into a new cluster or merging *Y* into another cluster *Ck*′​.  
> 3. Accept or reject the collective move using the exact ratio of the non-additive PBCT marginal posteriors *p*(Π*e*(*t*)​∣*x*)*p*(Π*e*prop​∣*x*)​ corrected by the boundary bond cut probabilities.

### **Interactive Exploration: Thermodynamic Sweep vs. PBCT Marginal Likelihood**

The visualization below simulates a context node *e* with vocabulary size *M*\=12 partitioned into three true latent equivalence classes on the simplex. Adjust the sequence sample size *N*, the CRP concentration *α*, and the temperature *T* to compare how the **SPC susceptibility *χ*(*T*)** and the **PBCT Log Marginal Posterior LPBCT​(Π(*T*))** co-evolve across the paramagnetic, superparamagnetic, and ferromagnetic regimes.

**Key architectural takeaway:** When sample size *Ne*​ is small (e.g., deep in the tree), the SPC susceptibility *χ*(*T*) will still exhibit geometric phase transitions on the empirical simplex, but the PBCT log Bayes factor ΔLPBCT​(*T*) drops below zero across the entire temperature axis—giving you an exact, hyperparameter-free **thermodynamic stopping rule** to prune the context tree into a leaf.

**Sources**

> 1. [(PDF) Approximate learning of parsimonious Bayesian context trees](https://www.researchgate.net/publication/382655019_Approximate_learning_of_parsimonious_Bayesian_context_trees)

---

*User prompt: This could be interesting. Lets formulate the fisher-rao \+dirichlet multinomial coupling matrix, then design the bottom-up causal state, and outline the swendsen wang cut acceptance ratio for PBCT. Also, something ive been wondering about is like in SPC you run simulations at different temperatures, there is so much information spread across simulations run at different temperatures. like im trying to figure out how to learn across temperatures and use inference to find the best clusters rather than detecting a susceptibility peak doing the cuts at that temperature. like there's gotta be come kind of way to learn from the rich information spread over the temperatures, but i just dont know how and i thought maybe this*

Response: Your intuition about the temperature sweep hits on the single biggest structural flaw in classical Domany–Blatt Superparamagnetic Clustering (SPC): **cutting at a single susceptibility peak *T*∗ discards almost all of the thermodynamic information you just paid to simulate, and fails whenever a dataset has multi-scale density regimes.**

Below is the complete mathematical formulation of all four pieces: (1) the Fisher-Rao \+ Dirichlet-Multinomial coupling matrix with sample-size gating, (2) the bottom-up Causal State suffix algorithm, (3) the exact Barbu–Zhu Swendsen-Wang Cut (SWC) acceptance ratio for PBCT, and (4) the **Multi-Temperature Anti-Chain & Persistence Integration architecture** for learning across the entire temperature spectrum simultaneously.

## **1\. Exact Fisher-Rao \+ Dirichlet-Multinomial Coupling Matrix**

Let *e*∈*V*∗ be a context node and *u*∈*V* (∣*V*∣=*M*) be a candidate preceding symbol with observed next-token transition counts x*eu*​\=(*xeu*,1​,…,*xeu*,*M*​), total count *Neu*​\=∑*m*\=1*M*​*xeu*,*m*​, and conjugate prior *ϕeu*​∼Dirichlet(*η*) with concentration *η*0​\=∑*m*\=1*M*​*ηm*​.

The posterior over the simplex Δ*M*−1 is Dirichlet(*αeu*​), where *αeu*​\=x*eu*​\+*η*, *αeu*,0​\=*Neu*​\+*η*0​, and the posterior predictive mean is *p*\~​*eu*​\=*αeu*​/*αeu*,0​.

### **A. Connecting Fisher-Rao Geodesics to the Dirichlet-Multinomial Bayes Factor**

On the probability simplex Δ*M*−1, the Riemannian metric induced by the Fisher information matrix *gab*​(*p*)=*δab*​/*pa*​ yields the **Fisher-Rao geodesic distance** (via the positive octant sphere embedding *ψm*​(*p*)=*pm*​​):

*d*FR​(*p*,*q*)=2arccos(*m*\=1∑*M*​*pm*​*qm*​​)∈\[0,*π*\]

Now examine the exact PBCT log Bayes factor Λ*uv*(*e*)​ for merging symbols *u* and *v* at context *e*:

Λ*uv*(*e*)​\=log*B*(x*eu*​\+*η*)*B*(x*ev*​\+*η*)*B*(x*eu*​\+x*ev*​\+*η*)*B*(*η*)​

Applying a second-order Laplace/Stirling expansion around the posterior means *p*\~​*eu*​,*p*\~​*ev*​ connects the Dirichlet-Multinomial marginal likelihood directly to the Fisher-Rao geodesic:

Λ*uv*(*e*)​\=−81​*N*eff(*e*)​(*u*,*v*)*d*FR2​(*p*\~​*eu*​,*p*\~​*ev*​)+2*M*−1​log((*αeu*,0​\+*αev*,0​)*η*0​*αeu*,0​*αev*,0​​)\+O(*d*FR4​,*N*eff−1​)

where the **effective pairwise sample size (harmonic reduced mass)** is:

*N*eff(*e*)​(*u*,*v*)=*Neu*​\+*Nev*​*Neu*​*Nev*​​

*(Note: If you want an exact closed-form non-asymptotic affinity on the simplex without Stirling's approximation, the exact Bhattacharyya coefficient between the two Dirichlet posteriors is BDir​(u,v)=∫ΔM−1​Dir(ϕ;αeu​)Dir(ϕ;αev​)​dϕ\=B(αeu​)B(αev​)​B((αeu​\+αev​)/2)​.)*

### **B. Solving the Zero-Count Percolation Trap via Debiased Precision Gating**

In standard geometric SPC, *d*(*u*,*v*)=0⟹*Juv*​\=*J*max​, because zero distance means identical coordinates. In statistical sequence modeling, if two symbols *u*,*v* never appear in context *e* (*Neu*​\=*Nev*​\=0), their posterior means both equal the prior mean (*p*\~​*eu*​\=*p*\~​*ev*​\=*η*/*η*0​), so *d*FR​(*p*\~​*eu*​,*p*\~​*ev*​)=0. If you assign them *J*max​, uninformative symbols form a ferromagnetic clique and short-circuit the Fortuin-Kasteleyn (FK) graph.

Conversely, if two symbols belong to the **same** true causal state (*ϕeu*​\=*ϕev*​) with finite counts, sampling noise inflates their expected squared geodesic distance to:

E\[*d*FR2​(*p*\~​*eu*​,*p*\~​*ev*​)​*ϕeu*​\=*ϕev*​\]≈*N*eff(*e*)​(*u*,*v*)+*η*0​*M*−1​

To separate **manifold geometry** from **sample uncertainty**, factorize the Potts coupling *Juv*(*e*)​ into two distinct terms:

> 1. **Noise-Debiased Squared Geodesic Distance:**  
>    *d*\~FR2​(*u*,*v*)=max(0,*d*FR2​(*p*\~​*eu*​,*p*\~​*ev*​)−*N*eff(*e*)​(*u*,*v*)+*η*0​*M*−1​)  
> 2. **Evidential Precision Gate:**  
>    *wuv*(*e*)​\=*N*eff(*e*)​(*u*,*v*)+*η*0​*N*eff(*e*)​(*u*,*v*)​∈\[0,1)

This yields the **Precision-Gated Fisher-Rao Potts Coupling**:

*Juv*(*e*)​\=*K*^*wuv*(*e*)​​exp(−2*ae*2​*d*\~FR2​(*u*,*v*)​)I(⟨*u*,*v*⟩∈N*k*(*e*)​)

where *ae*​ is the local robust scale (e.g., median mutual *k*\-NN geodesic distance on Δ*M*−1) and *K*^ is the mean weighted degree.

> * **Behavior as *Neu*​→0:** *N*eff(*e*)​(*u*,*v*)→0⟹*wuv*(*e*)​→0⟹*Juv*(*e*)​→0. Data-starved symbols automatically decouple into **paramagnetic spectators** during the Swendsen-Wang sweep and are assigned post-hoc to whichever macro-cluster *Ck*​ maximizes the exact Dirichlet-Multinomial predictive responsibility.  
> * **Behavior as *Neu*​,*Nev*​→∞:** *wuv*(*e*)​→1 and the variance bias vanishes, recovering the exact Riemannian Fisher-Rao kernel on Δ*M*−1.

## **2\. Bottom-Up Causal State SPC Suffix Algorithm**

Top-down PBCT fails on XOR, parity, and delayed-interaction Markov processes because depth-1 transition marginals look identical (*p*\~​*u*​≈*p*\~​*v*​), causing the root to collapse into a single leaf {*V*} before depth-2 contexts are ever inspected.

To fix this and unify PBCT with Crutchfield–Shalizi *ϵ*\-machines, we formulate a **bottom-up multi-layer renormalization over the active suffix lexicon**.

### **Step 1: Active Suffix Extraction & Leave-One-Out Vertical Conditioning**

Extract all contexts up to depth *D* with count at least *N*min​ (Constructed in O(*N*seq​) via a suffix tree):

E≤*D*​\={*e*∈*d*\=0⋃*D*​*Vd*:*Ne*​≥*N*min​}

**Crucial Statistical Catch on Vertical Edges:** In a suffix tree, a parent context's transition count vector is the sum of its children's count vectors: x*e*​\=∑*u*∈*V*​x*ue*​. If you compute a vertical coupling between child x*ue*​ and parent x*e*​, you are correlating x*ue*​ with itself (autogenous bias of order *Nue*​/*Ne*​). The statistically exact test of whether child *ue* splits from parent *e* must compare x*ue*​ against the **complement sibling pool**:

x*e*∖*ue*​\=x*e*​−x*ue*​,*Ne*∖*ue*​\=*Ne*​−*Nue*​

Because x*ue*​ and x*e*∖*ue*​ are disjoint multinomial draws, they are conditionally independent given *ϕe*​.

### **Step 2: Three-Layer Anisotropic Coupling Graph G\=(E≤*D*​,A)**

We construct three distinct edge sets over E≤*D*​:

> 1. **Horizontal Sibling Couplings** ⟨*ue*,*ve*⟩∈Asib​ (same parent suffix *e*): Coupled via *Jue*,*ve*sib​ using the Precision-Gated Fisher-Rao metric between x*ue*​ and x*ve*​. Locking these edges forms **PBCT parsimonious partitions** inside node *e*.  
> 2. **Vertical Parent-Child Couplings** ⟨*ue*,*e*⟩∈Avert​: Coupled via *Jue*,*e*vert​ computed between x*ue*​ and the sibling complement x*e*∖*ue*​. Locking all active children *ue* onto *e* **prunes** node *e* into a leaf.  
> 3. **Cross-Branch Causal Couplings (Optional *ϵ*\-Machine Mode)** ⟨*e*,*e*′⟩∈Acausal​ (*e*,*e*′ do not share a parent): Mutual *k*\-NN Fisher-Rao edges across different branches of the tree. For example, in an XOR process (*Xt*​\=*Xt*−1​⊕*Xt*−2​), histories *e*\=01 and *e*′\=10 live in separate branches of the context tree but have identical predictive distributions and merge into a single **causal state**.

### **Step 3: Bottom-Up Renormalization Pass (*d*\=*D* down to 1)**

Instead of terminating top-down at depth 1:

> 1. Run the multi-temperature SPC partition inference (detailed in Section 4\) at every parent node *e* from the deepest internal level *d*\=*D*−1 up to the root *d*\=0.  
> 2. When evaluating parent node *e* at depth *d*, its children *ue* (at depth *d*\+1) have **already resolved their subtree partitions** Π*ue*∗​.  
> 3. The effective log marginal likelihood of child *ue* is therefore the **subtree-propagated evidence**:  
>    Lsub∗​(*ue*)=max(log*p*DM​(x*ue*​∣*η*),LPBCT​(Π*ue*∗​∣{x*wue*​}*w*∈*V*​))  
> 4. Parent *e* only collapses children *ue* and *ve* together if their pooled leaf likelihood exceeds the sum of their **subtree-propagated** marginal likelihoods Lsub∗​(*ue*)+Lsub∗​(*ve*). In an XOR process, depth-2 children have high subtree marginal likelihoods Lsub∗​(*ue*)≫log*p*DM​(x*ue*​), which blocks the depth-1 root from prematurely collapsing\!

## **3\. Exact Swendsen-Wang Cut (Barbu–Zhu) Acceptance Ratio for PBCT**

If you want to sample the exact posterior over PBCT partitions *P*(Π*e*​∣X*e*​) rather than taking a single point estimate, standard Swendsen-Wang fails because the PBCT log-posterior is **not** a sum of pairwise Potts interactions—it is a non-additive product of multivariate Beta functions over arbitrary-sized coalitions *Ck*​∈Π*e*​:

*P*(Π*e*​∣X*e*​)∝*α*∣Π*e*​∣*C*∈Π*e*​∏​Γ(∣*C*∣)*B*(*η*)*B*(X*C*​\+*η*)​,X*C*​\=*u*∈*C*∑​x*eu*​

**Swendsen-Wang Cuts (SWC)** (Barbu & Zhu, 2005\) turns the pairwise Fisher-Rao FK bond percolation into a **proposal distribution** for collective multi-symbol jumps, corrected by a Metropolis-Hastings acceptance step.

### **A. The SWC Proposal Mechanism at Temperature *T***

Given current partition Π={*C*1​,…,*CK*​} of *V*:

> 1. **Within-Cluster FK Percolation:** For every edge ⟨*u*,*v*⟩ where *u*,*v* belong to the *same* cluster *Ck*​, turn the bond ON independently with probability *quv*​\=1−exp(−*Juv*(*e*)​/*T*). All inter-cluster edges are kept OFF.  
> 2. **Select a Sub-Coalition:** The active bonds partition each *Ck*​ into connected FK components VFK​(Π). Pick one connected component *V*0​⊆*Ck*​ according to a selection distribution *γ*(*V*0​∣VFK​) (e.g., uniform 1/∣VFK​∣).  
> 3. **Propose Target Cluster:** Propose reassigning the entire coalition *V*0​ to a target cluster *Cl*​ (*l*∈{1,…,*K*\+1}, where *CK*\+1​\=∅ represents splitting *V*0​ into a brand-new cluster) with proposal probability *Q*(*l*∣*V*0​,Π), forming candidate partition Π′.

### **B. Exact Cancellation of Internal Bonds & The Acceptance Ratio**

To move from Π→Π′, every edge in the **cut-set** C(*V*0​,*Ck*​∖*V*0​)={⟨*u*,*v*⟩:*u*∈*V*0​,*v*∈*Ck*​∖*V*0​} had to fail to turn ON (probability ∏C(*V*0​,*Ck*​∖*V*0​)​(1−*quv*​)). In the reverse move Π′→Π, *V*0​ is embedded inside *Cl*​∪*V*0​, so every edge in the reverse cut-set C(*V*0​,*Cl*​) must fail to turn ON (probability ∏C(*V*0​,*Cl*​)​(1−*quv*​)), while **every bond internal to *V*0​ has the exact same probability in both directions and cancels out completely\!**

Since 1−*quv*​\=exp(−*Juv*(*e*)​/*T*), the Metropolis-Hastings acceptance probability *α*(Π→Π′) is:

*α*(Π→Π′)=min(1,exp(*TJ*(*V*0​,*Ck*​∖*V*0​)−*J*(*V*0​,*Cl*​)​)⋅*γ*(*V*0​∣Π)*Q*(*l*∣*V*0​,Π)*γ*(*V*0​∣Π′)*Q*(*k*∣*V*0​,Π′)​⋅*P*(Π∣X*e*​)*P*(Π′∣X*e*​)​)

where *J*(*A*,*B*)=∑*u*∈*A*,*v*∈*B*​*Juv*(*e*)​ is the total cut coupling between sets *A* and *B*, and the exact PBCT posterior ratio for transferring *V*0​⊊*Ck*​ into *Cl*​\=∅ is:

*P*(Π∣X*e*​)*P*(Π′∣X*e*​)​\=Γ(∣*Ck*​∣)Γ(∣*Cl*​∣)Γ(∣*Ck*​∣−∣*V*0​∣)Γ(∣*Cl*​∣+∣*V*0​∣)​⋅*B*(X*Ck*​​\+*η*)*B*(X*Cl*​​\+*η*)*B*(X*Ck*​​−X*V*0​​\+*η*)*B*(X*Cl*​​\+X*V*0​​\+*η*)​

*(For a pure split V0​⊊Ck​→CK\+1​\=∅, replace Γ(∣Cl​∣+∣V0​∣)/Γ(∣Cl​∣) with αΓ(∣V0​∣) and B(XCl​​\+XV0​​\+η)/B(XCl​​\+η) with B(XV0​​\+η)/B(η); for a pure merge V0​\=Ck​→Cl​, invert the split ratio.)*

**Why this mixes so fast:** Look at the product of the cut factor and the Bayes factor in log space:

ΔlogAccept≈Exact Multi-Body Energy −Δ*E*PBCT​log*P*(Π∣X*e*​)*P*(Π′∣X*e*​)​​​−Pairwise Surrogate Energy −Δ*E*Potts​*TJ*(*V*0​,*Cl*​)−*J*(*V*0​,*Ck*​∖*V*0​)​​​

Because the pairwise Fisher-Rao couplings approximate the 2-body Taylor expansion of −log*P*(Π∣X*e*​), the **2-body energy barriers cancel out in the exponent**, leaving only 3-body and higher-order Dirichlet-Multinomial coalition synergies to govern rejection.

## **4\. Learning Across Temperatures: Multi-Scale Anti-Chain & Persistence Inference**

Now let us address your core question: **How do we synthesize the rich information spread across simulations at all temperatures {*T*1​,…,*TR*​} instead of slicing at a single susceptibility peak *T*∗?**

### **Why Single-*T* Susceptibility Cuts Are Fundamentally Broken**

Suppose a node (or manifold dataset) contains two macro-branches:

> * **Branch A** has high sample density (or tight geometric curvature) and splits cleanly into clusters *A*1​,*A*2​ on the temperature plateau *T*∈\[0.08,0.16\], then shatters into paramagnetic singletons at *T*\>0.18.  
> * **Branch B** has lower sample density (or broader variance) and remains fused until *T*\=0.22, resolving into *B*1​,*B*2​ on the plateau *T*∈\[0.22,0.34\].

Here, **no single global temperature *T*∗ exists** where {*A*1​,*A*2​,*B*1​,*B*2​} coexist\! Slicing at the first susceptibility peak (*T*\=0.12) under-splits *B*; slicing at the second peak (*T*\=0.28) shatters *A* into dust. Furthermore, slicing at one *T*∗ throws away *R*−1 Markov chain simulations.

There are **three rigorous mathematical frameworks** to integrate across the entire temperature spectrum simultaneously:

### **Framework 1: The Multi-Temperature Grain Hierarchy & Optimal Anti-Chain DP (Best for PBCT \+ SPC)**

Across an evenly spaced temperature grid 0=*T*0​\<*T*1​\<⋯\<*TR*​\=*T*max​, the FK simulations generate an evolving ensemble of candidate clusters (superparamagnetic grains).

> 1. **Extract the Candidate Grain Pool H:** Let H\={*C*⊆*V*} be the set of all macro-clusters that appear across the temperature sweep (either as connected components of the thresholded thermal co-association graph *Guv*​(*Tr*​)≥*τ* across *r*\=1,…,*R*, or directly as high-frequency FK components). By construction, as *T* increases from 0 to *T*max​, grains refine from the root *V* down to singletons {*u*}, forming a **Grain Hierarchy (Component Tree)** H with less than or equal to 2*M*−1 nodes\!  
> 2. **Compute the Thermal Persistence Integral *τ*pers​(*C*) for Every Grain:** For each candidate grain *C*∈H, let *p*intact​(*C*∣*Tr*​) be the empirical probability at temperature *Tr*​ that *C* forms an isolated connected FK component (all internal bonds percolate and all boundary bonds to *V*∖*C* break). Integrating over the temperature axis (which you can evaluate cleanly via your even-grid BARS/spline functional estimator) gives the **Thermodynamic Lifetime** of grain *C*:  
>    *τ*pers​(*C*)=∫*T*min​*T*max​​*p*intact​(*C*∣*T*)*dT*≈*r*\=1∑*R*​*p*intact​(*C*∣*Tr*​)Δ*Tr*​  
> 3. **Assign Each Grain Its Self-Contained Bayesian \+ Thermodynamic Score Ψ(*C*):** Notice a crucial mathematical property of the PBCT marginal posterior:  
>    LPBCT​(Π)=logΓ(*α*\+*M*)Γ(*α*)​\+*C*∈Π∑​*ψ*PBCT​(*C*)  
>    where every candidate cluster *C*∈H has a **strictly local, additive score**:  
>    *ψ*PBCT​(*C*)=log*α*\+logΓ(∣*C*∣)+log*B*(X*C*​\+*η*)−log*B*(*η*)  
>    We can weight each grain by its pure Bayesian marginal log-likelihood *ψ*PBCT​(*C*) (or augment it with a thermodynamic persistence prior *β*pers​log*τ*pers​(*C*) for pure unsupervised manifold SPC where no Dirichlet-Multinomial model exists).  
> 4. **Solve for the Optimal Multi-Temperature Cut via Bottom-Up Dynamic Programming:** A valid partition Π of *V* is simply a **maximal anti-chain (tree cut)** of the multi-temperature grain hierarchy H—every symbol *u*∈*V* is covered by exactly one selected ancestor grain *C*∈H. Because ∑*C*∈Π​Ψ(*C*) is additive over disjoint grains, we find the exact globally optimal **cross-temperature partition Π∗** in O(∣H∣) time via a single bottom-up pass from the singleton leaves to the root *V*:  
>    *V*∗(*C*)=max​Ψ(*C*),*C*′∈children(*C*)∑​*V*∗(*C*′)​  
>    If Ψ(*C*)≥∑*C*′∈children(*C*)​*V*∗(*C*′), we **lock grain *C* at its native birth temperature *T*(*C*)** and prune its higher-temperature fragments; otherwise, we split *C* into its optimal children\! **Result:** Branch A automatically selects its clusters from *T*\=0.12 while Branch B selects its clusters from *T*\=0.28, completely eliminating single-*T* susceptibility cuts.

### **Framework 2: Thermodynamic Integration of the Thermal Co-Association Curve *Guv*​(*T*)**

Even before forming discrete grains, look at the pairwise spin-spin correlation (FK co-occurrence probability) *Guv*​(*T*)=P*T*​(fk(*u*)=fk(*v*)) across *T*∈\[0,*T*max​\].

> * For every pair (*u*,*v*), *Guv*​(0)=1 and lim*T*→∞​*Guv*​(*T*)=0.  
> * If we define the random variable *T*melt​(*u*,*v*) as the thermal scale at which the bond between *u* and *v* breaks, then −*dTdGuv*​(*T*)​ is the probability density of *T*melt​(*u*,*v*).  
> * By integration by parts, the **area under the correlation curve across all simulated temperatures** is the **Expected Critical Melting Temperature** between *u* and *v*:  
>   E\[*T*melt​(*u*,*v*)\]=∫0*T*max​​*Guv*​(*T*)*dT*  
> * Consequently, the **Integrated Thermodynamic Dissimilarity**  
>   *D*thermo​(*u*,*v*)=1−*T*max​1​∫0*T*max​​*Guv*​(*T*)*dT*  
>   combines every single simulation across all temperatures {*T*1​,…,*TR*​} into a noise-suppressed **multi-scale persistence metric**. Because *D*thermo​(*u*,*v*) integrates over all *T*, two points that share a wide superparamagnetic plateau have a huge gap Δ*D*thermo​ relative to points that split early, turning the multi-temperature grain tree H in Framework 1 into a canonical persistence dendrogram (0-dimensional persistent homology over the thermal filtration\!).

### **Framework 3: Bayesian Model Averaging Over Temperature *T***

If you treat temperature *T* as a continuous hyperparameter governing the FK proposal distribution *P*FK​(Π∣*T*), the marginal posterior probability of any cluster assignment *C* integrating out *T* is:

*P*(*C*∈Π∣X)=∫*T*min​*T*max​​*P*FK​(*C*∈Π∣*T*)*P*(*T*∣X)*dT*

where the posterior weight of temperature *Tr*​ is given directly by the PBCT marginal evidence of the partitions sampled at *Tr*​:

*P*(*Tr*​∣X)∝*p*(*Tr*​)*S*1​*s*\=1∑*S*​exp(LPBCT​(Π(*s*)(*Tr*​)))

Rather than hunting for a noisy peak in susceptibility *χ*(*T*), the Bayesian evidence exp(LPBCT​) automatically concentrates the posterior measure *P*(*T*∣X) onto the true superparamagnetic plateau.

## **Interactive Exploration: Multi-Temperature Anti-Chain DP vs. Single-*T* Cuts**

The simulation below demonstrates **Framework 1** on a heterogeneous 12-symbol context node containing two branches with **incompatible critical temperatures**:

> * **Dense Branch A** ({*s*1​,…,*s*6​}) resolves into *A*1​,*A*2​ at **low temperature** (*T*∈\[0.08,0.18\]) and shatters into singletons above *T*\=0.20.  
> * **Sparse Branch B** ({*s*7​,…,*s*12​}) stays fused until *T*\=0.22, resolving into *B*1​,*B*2​ at **high temperature** (*T*∈\[0.24,0.38\]).

Slide the **Single-Temperature Cut *T*** to see why every single-temperature slice fails, then toggle to **Multi-Temperature Anti-Chain DP** to see how integrating across all temperatures locks *A*1​,*A*2​ at *TA*∗​\=0.13 and *B*1​,*B*2​ at *TB*∗​\=0.30 simultaneously.

### **How This Connects to Your Even-Grid Spline / Functional Architecture**

Notice why **Framework 1 \+ Framework 2** pairs naturally with running independent simulations on an **evenly spaced temperature grid** *T*1​,…,*TR*​:

> 1. **Embarrassingly Parallel Simulation across *T*:** Every temperature *Tr*​ runs independently on the precision-gated Fisher-Rao graph.  
> 2. **Smooth Functional Estimation Over *T*:** For each pair ⟨*u*,*v*⟩ (or candidate grain *C*∈H), the empirical survival observations {*G*^*uv*​(*Tr*​)}*r*\=1*R*​ across the even temperature grid can be smoothed or integrated via quadrature / Bayesian spline regression without adaptive temperature-stepping dependencies.  
> 3. **Zero Threshold Tuning:** Instead of hunting for derivatives of global susceptibility *χ*(*T*), the integrated thermal co-association *D*thermo​(*u*,*v*)=1−*T*max​1​∫*Guv*​(*T*)*dT* builds the candidate grain tree H, and the bottom-up recursion *V*∗(*C*)=max(Ψ(*C*),∑*C*′∈ch(*C*)​*V*∗(*C*′)) extracts the multi-scale partition across temperatures in one pass.

---

*User prompt: ok lets write the exact O(M log M) algorithm to build grain hierarchy from G\_uv(T) , then formulate the unsupervised manifold score when no multinomial count exists*

Response: To build the multi-temperature grain hierarchy in O(*M*log*M*) time without building an O(*M*2) dense correlation matrix, we exploit a topological equivalence between **single-linkage filtrations** and **Maximum Spanning Trees (MST)** over the sparse Potts coupling graph E (∣E∣=O(*kM*)). We then formulate three self-contained unsupervised manifold scores Ψ(*C*)—from a zero-parameter thermodynamic persistence functional to a rigorous Riemannian Laplace-Beltrami marginal likelihood—that replace the Dirichlet-Multinomial likelihood when operating on general continuous manifolds.

## **Part 1: Exact O(∣E∣log∣E∣)=O(*M*log*M*) Grain Hierarchy Construction**

### **1\. Why Naive Multi-Temperature Clustering is O(*RM*2)—and How to Collapse It**

In classical SPC, users compute the dense *M*×*M* spin-spin correlation matrix *Guv*​(*Tr*​)=P*Tr*​​(fk(*u*)=fk(*v*)) at each temperature *Tr*​ (*r*\=1,…,*R*) and threshold it at *Guv*​(*Tr*​)≥*θ*. Doing this naively has two fatal flaws:

> 1. **Computational bottleneck:** Tracking all pairs costs O(*RM*2) memory and time.  
> 2. **Non-nested (braided) clusters across *T*:** Because finite-sample MCMC estimates *G*^*uv*​(*Tr*​) fluctuate slightly across independent temperature runs, thresholding at *Tr*​ and *Tr*\+1​ can produce clusters that partially overlap rather than cleanly nesting, breaking the tree structure required for bottom-up dynamic programming.

**The Mathematical Reduction:** Two vertices *u*,*v*∈*V* belong to the same Fortuin-Kasteleyn (FK) connected component in a single sweep if and only if they are connected by a path of active edges *e*∈E. Therefore, the entire multi-scale topology is already completely encoded on the **sparse coupling edge set** E (∣E∣≤*kM*):

> * For each edge *e*\=⟨*u*,*v*⟩∈E at temperature *Tr*​, record the empirical co-cluster probability across *S* Swendsen-Wang sweeps:  
>   *Ge*​(*Tr*​)=*S*1​*s*\=1∑*S*​I(fk(*s*)(*u*,*Tr*​)=fk(*s*)(*v*,*Tr*​))  
> * Across your even temperature grid 0=*T*0​\<*T*1​\<⋯\<*TR*​\=*T*max​ (with spacing Δ*T*\=*T*max​/*R* and boundary conditions *Ge*​(0)=1, *Ge*​(*T*max​)=0), integrate the edge survival curve via composite trapezoidal or Simpson quadrature (or your even-grid spline integral) to get the **Expected Critical Melting Temperature** *T*ˉ*e*​ of edge *e*:  
>   *T*ˉ*e*​\=∫0*T*max​​*Ge*​(*T*)*dT*≈Δ*T*(21​*Ge*​(*T*0​)+*r*\=1∑*R*−1​*Ge*​(*Tr*​)+21​*Ge*​(*TR*​))  
> * Define the **Thermodynamic Ultrametric Dissimilarity** *De*​\=*T*max​−*T*ˉ*e*​∈\[0,*T*max​\].

Because *T*ˉ*e*​ is monotonically related to the thermal scale at which edge *e* breaks, the connected components formed by inserting edges in decreasing order of *T*ˉ*e*​ (from *T*max​ down to 0\) form a strictly nested **Binary Merge Tree (Dendrogram)** T with exactly 2*M*−1 nodes. Furthermore, any cycle in E is redundant for connectivity, so T is determined entirely by the **Maximum Spanning Tree** of (E,*T*ˉ*e*​) via Kruskal's algorithm in O(∣E∣log∣E∣)=O(*M*log*M*) time\!

### **2\. Exact Definition of Thermal Persistence *τ*pers​(*C*) on the Tree**

When edge *ek*​\=⟨*uk*​,*vk*​⟩ merges left child cluster *CL*​ (formed at melting temperature *T*ˉ(*CL*​)) and right child cluster *CR*​ (formed at *T*ˉ(*CR*​)) into parent cluster *CP*​\=*CL*​∪*CR*​, the parent cluster *CP*​ becomes fully connected at temperature:

*T*death​(*CP*​)=*T*ˉ*ek*​​≤min(*T*death​(*CL*​),*T*death​(*CR*​))

For each child *C*∈{*CL*​,*CR*​}, *C* exists as a **coherent, isolated grain** over the temperature window between the temperature where *C* splits from its sibling (*T*birth​(*C*)=*T*death​(*CP*​)=*T*ˉ*ek*​​) and the higher temperature where *C* itself shatters into its own children (*T*death​(*C*)):

> * **Asymptotic / Ultrametric Persistence:**  
>   Δ*T*(*C*)=*T*death​(*C*)−*T*birth​(*C*)≥0  
> * **Exact Finite-*T* Grid Integral (No Ultrametric Approximation):** Want the exact probability *p*intact​(*C*∣*Tr*​) at each grid temperature *Tr*​ without the MST bottleneck approximation? For any node *C*∈T, let *e*in​(*C*) be the weakest internal MST edge that holds *C* together (*e*in​(*C*)=argmin*e*∈MST(*C*)​*T*ˉ*e*​, which is the merge edge of *C*), and let ∂*C*\={*e*∈E:∣*e*∩*C*∣=1} be the cut-set of edges leaving *C*, with strongest boundary edge *e*out​(*C*)=argmax*e*∈∂*C*​*T*ˉ*e*​ (the edge that merges *C* into its parent *CP*​). At temperature *Tr*​, *C* is intact if *e*in​(*C*) is connected and *e*out​(*C*) is broken, giving the exact curve difference:  
>   *τ*pers​(*C*)=∫0*T*max​​max(0,*Ge*in​(*C*)​(*T*)−*Ge*out​(*C*)​(*T*))*dT*\=*T*ˉ*e*in​(*C*)​−*T*ˉ*e*out​(*C*)​\=*T*death​(*C*)−*T*birth​(*C*)  
>   Notice that by linearity of integration and the monotonicity *Ge*in​(*C*)​(*T*)≥*Ge*out​(*C*)​(*T*), the exact integral of the curve difference ∫(*Ge*in​​(*T*)−*Ge*out​​(*T*))*dT* **collapses algebraically to *T*ˉ*e*in​(*C*)​−*T*ˉ*e*out​(*C*)​\!** You don't even need to store the *R* curves during tree construction—integrating *T*ˉ*e*​ per edge first and taking differences on the Kruskal tree gives the exact functional integral\!

*(Degeneracy Compression / ϵ\-Collapse: If multiple internal edges inside a macro-cluster break at essentially the same temperature, ΔT(C)\<ϵT​, we contract child C into parent CP​ so H becomes a k\-ary hierarchy of strictly persistent macro-grains.)*

### **3\. Zero-Allocation O(*M*log*M*) Algorithm (Span\<T\> / Array-Backed)**

Here is the complete algorithm using flat arrays of length 2*M*−1 (leaves 0…*M*−1, internal nodes *M*…2*M*−2):

namespace Spcx.Thermodynamics;

using System;  
using System.Numerics;  
using System.Runtime.CompilerServices;

public readonly record struct CouplingEdge(int U, int V, double TBar, double JCoupling);

public struct GrainNode  
{  
    public int LeftChild;        // \-1 for singleton leaves (0..M-1)  
    public int RightChild;       // \-1 for singleton leaves  
    public int Parent;           // \-1 for root (2M-2)  
    public int Size;             // \|C\|  
    public double TDeath;        // Temperature where C shatters internally (T\_max for leaves)  
    public double TBirth;        // Temperature where C merges into Parent (0 for root)  
    public double Persistence;   // tau\_pers(C) \= TDeath \- TBirth  
    public double InternalCutJ;  // Sum of J\_uv completely inside C  
    public double BoundaryCutJ;  // Sum of J\_uv in cut-set dC \= {(u,v) : u in C, v not in C}  
}

public static class GrainHierarchyBuilder  
{  
    /// \<summary\>  
    /// Builds the exact (2M \- 1)-node Multi-Temperature Grain Hierarchy in O(\|E\| log \|E\|) time.  
    /// Zero heap allocations in the tree construction loop if buffers are pre-allocated.  
    /// \</summary\>  
    public static int BuildHierarchy(  
        int M,  
        double tMax,  
        Span\<CouplingEdge\> edges,         // Length \|E\|, sorted descending by TBar in-place  
        ReadOnlySpan\<double\> vertexDegreeJ, // Length M: sum\_{v} J\_{uv} for each u  
        Span\<GrainNode\> treeNodes,        // Output buffer of length 2M \- 1  
        Span\<int\> dsuParent,              // Scratch buffer of length 2M \- 1  
        Span\<int\> dsuRootNode)            // Scratch buffer of length M  
    {  
        // 1\. Sort sparse edges in descending order of Expected Melting Temperature TBar: O(\|E\| log \|E\|)  
        edges.Sort(static (a, b) \=\> b.TBar.CompareTo(a.TBar));

        // 2\. Initialize M singleton leaves  
        for (int i \= 0; i \< M; i++)  
        {  
            treeNodes\[i\] \= new GrainNode  
            {  
                LeftChild \= \-1,  
                RightChild \= \-1,  
                Parent \= \-1,  
                Size \= 1,  
                TDeath \= tMax,  
                TBirth \= 0.0,  
                Persistence \= tMax,  
                InternalCutJ \= 0.0,  
                BoundaryCutJ \= vertexDegreeJ\[i\]  
            };  
            dsuParent\[i\] \= i;  
            dsuRootNode\[i\] \= i;  
        }

        for (int i \= M; i \< 2 \* M \- 1; i++)  
            dsuParent\[i\] \= i;

        // 3\. Kruskal Canonical Dendrogram Pass: O(\|E\| alpha(M))  
        int nextInternalId \= M;  
        for (int eIdx \= 0; eIdx \< edges.Length; eIdx++)  
        {  
            ref readonly CouplingEdge edge \= ref edges\[eIdx\];  
            int rootU \= Find(edge.U, dsuParent);  
            int rootV \= Find(edge.V, dsuParent);

            if (rootU \== rootV)  
            {  
                // Edge is internal to an already-formed cluster; accumulate internal J if desired  
                continue;  
            }

            // Create new parent grain node C\_P \= C\_L U C\_R at melting temperature edge.TBar  
            int parentId \= nextInternalId++;  
            ref GrainNode left \= ref treeNodes\[rootU\];  
            ref GrainNode right \= ref treeNodes\[rootV\];

            left.Parent \= parentId;  
            left.TBirth \= edge.TBar;  
            left.Persistence \= Math.Max(0.0, left.TDeath \- left.TBirth);

            right.Parent \= parentId;  
            right.TBirth \= edge.TBar;  
            right.Persistence \= Math.Max(0.0, right.TDeath \- right.TBirth);

            double internalJ \= left.InternalCutJ \+ right.InternalCutJ \+ edge.JCoupling;  
            double totalDegreeJ \= (left.BoundaryCutJ \+ 2.0 \* left.InternalCutJ)  
                                \+ (right.BoundaryCutJ \+ 2.0 \* right.InternalCutJ);

            treeNodes\[parentId\] \= new GrainNode  
            {  
                LeftChild \= rootU,  
                RightChild \= rootV,  
                Parent \= \-1,  
                Size \= left.Size \+ right.Size,  
                TDeath \= edge.TBar,  
                TBirth \= 0.0,  
                Persistence \= edge.TBar, // Updated when parentId merges higher up  
                InternalCutJ \= internalJ,  
                BoundaryCutJ \= Math.Max(0.0, totalDegreeJ \- 2.0 \* internalJ)  
            };

            // Union both subtrees into parentId  
            dsuParent\[rootU\] \= parentId;  
            dsuParent\[rootV\] \= parentId;

            if (nextInternalId \== 2 \* M \- 1\)  
                break;  
        }

        return nextInternalId \- 1; // Root index (2M \- 2 if graph is connected)  
    }

    \[MethodImpl(MethodImplOptions.AggressiveInlining)\]  
    private static int Find(int i, Span\<int\> parent)  
    {  
        while (parent\[i\] \!= i)  
        {  
            parent\[i\] \= parent\[parent\[i\]\]; // Path halving  
            i \= parent\[i\];  
        }  
        return i;  
    }  
}

## **Part 2: Formulating the Unsupervised Manifold Score Ψ(*C*) Without Multinomial Counts**

When clustering points {*x*1​,…,*xM*​} on a general d-dimensional Riemannian manifold M—such as subspaces on the Grassmannian Gr(*k*,*n*), orthonormal frames on the Stiefel manifold St(*k*,*n*), hyperbolic embeddings on the Poincaré ball B*d*, or probability vectors where raw counts *Nu*​ don't exist—we no longer have the Dirichlet-Multinomial Beta function log*B*(X*C*​\+*η*).

Instead, we need an **additive cluster quality functional** L(Π)=∑*C*∈Π​Ψ(*C*) so the bottom-up dynamic program

*V*∗(*C*)=max(Ψ(*C*),*V*∗(*CL*​)+*V*∗(*CR*​))

can select the optimal multi-temperature anti-chain across T in O(*M*) time.

Here are **three rigorous formulations** of Ψ(*C*), ordered from pure topology/thermodynamics to full Riemannian Bayesian generative modeling. You can use them individually or combine them.

### **Formulation 1: Pure Thermodynamic Excess Persistence (Zero-Parameter / EOM Analogue)**

If you want a **strictly non-parametric score** that relies *only* on the Potts temperature sweep (no manifold center computation required\!), look at the **Mass-Weighted Thermodynamic Persistence Volume** of grain *C*:

Ψpers​(*C*)=∣*C*∣*γτ*pers​(*C*)=∣*C*∣*γ*∫*T*birth​(*C*)*T*death​(*C*)​*p*intact​(*C*∣*T*)*dT*\=∣*C*∣*γ*(*T*death​(*C*)−*T*birth​(*C*))

where *γ*∈\[0,1\] is a mass-scaling exponent (typically *γ*\=1, analogous to Excess of Mass in HDBSCAN, or *γ*\=1/2 to prevent singletons from dominating).

**Why *γ*\=1 needs a Singleton Cutoff (Calling BS on raw *γ*\=1):** Notice a crucial algebraic identity when *γ*\=1: because ∣*CP*​∣=∣*CL*​∣+∣*CR*​∣ and *T*death​(*CP*​)=*T*birth​(*CL*​)=*T*birth​(*CR*​), the sum of a parent and its children telescope\! Specifically:

∣*CL*​∣*τ*pers​(*CL*​)+∣*CR*​∣*τ*pers​(*CR*​)=∣*CL*​∣*T*death​(*CL*​)+∣*CR*​∣*T*death​(*CR*​)−∣*CP*​∣*T*death​(*CP*​)

If you allow singletons (∣*C*∣=1) or micro-clusters (∣*C*∣\<*m*min​) to compete in the DP with *T*death​({*u*})=*T*max​, the leaves will trivially win every time because *T*max​≥*T*death​(*C*). **The Exact Fix:** Set Ψpers​(*C*)=−∞ for micro-clusters ∣*C*∣\<*m*min​ (e.g., *m*min​\=5), and for valid macro-clusters ∣*C*∣≥*m*min​ define the **Thermodynamic Work Integral**:

Ψthermo​(*C*)=Integrated Cohesive Free Energy*J*in​(*C*)*τ*pers​(*C*)​​−Boundary Isolation Penalty*J*cut​(∂*C*)*T*birth​(*C*)​​

where *J*in​(*C*)=∑⟨*u*,*v*⟩⊂*C*​*Juv*​ is the internal Potts coupling and *J*cut​(∂*C*)=∑*u*∈*C*,*v*∈/*C*​*Juv*​ is the cut-set coupling (both computed in O(1) inside GrainNode during the Kruskal pass above\!).

> * A grain *C* maximizes Ψthermo​(*C*) when it has high internal bond energy *J*in​(*C*), survives over a wide temperature plateau *τ*pers​(*C*), and has minimal boundary leakage *J*cut​(∂*C*).

### **Formulation 2: Riemannian Pseudo-Gaussian / Wrapped Normal Marginal Likelihood (Continuous Conjugate Analogue of PBCT)**

What is the exact continuous Riemannian counterpart of the Dirichlet-Multinomial marginal likelihood log*B*(X*C*​\+*η*)/*B*(*η*)? On a *d*M​\-dimensional Riemannian manifold M with geodesic distance *d*M​(⋅,⋅), suppose points *xu*​∈*C* are drawn from a **Riemannian Normal distribution** centered at the cluster's **Fréchet mean** *μC*​∈M:

*μC*​\=arg*μ*∈Mmin​*u*∈*C*∑​*d*M2​(*xu*​,*μ*)

Let *WC*​\=∑*u*∈*C*​*d*M2​(*xu*​,*μC*​) be the **within-cluster Riemannian sum of squared geodesic dispersion**.

**How to compute *WC*​ in O(1) per tree node without running a Fréchet / Weiszfeld optimizer on every node:** Running an iterative Weiszfeld or Riemannian gradient descent solver on all *M*−1 internal nodes of T would cost O(*M*2). Instead, use one of two exact shortcuts:

> 1. **Riemannian Huygens–Steiner / Pairwise Dispersion Identity:** Around any local patch of M with bounded sectional curvature *κ*, the sum of squared geodesics to the Fréchet mean equals half the average pairwise squared geodesic distance inside *C* (up to O(*κWC*2​) curvature correction):  
>    *WC*​\=*u*∈*C*∑​*d*M2​(*xu*​,*μC*​)≈2∣*C*∣1​*u*,*v*∈*C*∑​*d*M2​(*xu*​,*xv*​)  
> 2. **Chordal Isometric Embedding (Exact O(1) Tree Merge):** For almost every manifold in your stack, the squared geodesic distance is locally equivalent to a **chordal Euclidean inner product** in an ambient embedding Φ:M→R*D*:  
   * **Simplex Δ*M*−1 (Fisher-Rao):** Hellinger sphere map Φ(*p*)=*p*​∈*S*\+*M*−1​, where *d*FR​(*p*,*q*)=2arccos⟨Φ(*p*),Φ(*q*)⟩≈2∥Φ(*p*)−Φ(*q*)∥2​.  
   * **Grassmannian Gr(*k*,*n*):** Projection matrix embedding Φ(*U*)=*UU*⊤∈Sym(*n*), where 21​∥Φ(*U*)−Φ(*V*)∥*F*2​\=∑*i*\=1*k*​sin2*θi*​(*U*,*V*)≈∑*i*\=1*k*​*θi*2​\=*d*Gr2​(*U*,*V*).  
   * **Stiefel St(*k*,*n*):** Ambient *n*×*k* matrix embedding Φ(*U*)=*U*, with chordal metric ∥*U*−*V*∥*F*2​\=2*k*−2Tr(*U*⊤*V*). Because ∥Φ(*xu*​)−Φˉ*C*​∥22​ is quadratic, every node *C*∈T only has to store the **additive first moment vector** S1​(*C*)=∑*u*∈*C*​Φ(*xu*​) and **scalar second moment** *S*2​(*C*)=∑*u*∈*C*​∥Φ(*xu*​)∥2. When merging *CP*​\=*CL*​∪*CR*​:  
     S1​(*CP*​)=S1​(*CL*​)+S1​(*CR*​),*S*2​(*CP*​)=*S*2​(*CL*​)+*S*2​(*CR*​),*WCP*​​\=*S*2​(*CP*​)−∣*CP*​∣∥S1​(*CP*​)∥2​  
     And once the DP selects the final winning clusters *C*∈Π∗ (typically *K*≪*M* clusters), you can polish those *K* cluster centers with 5 iterations of your exact Riemannian Weiszfeld solver\!

**Integrating Out (*μC*​,*σC*2​) to Get the Exact Riemannian Marginal Score ΨRiem​(*C*):** Placing a conjugate Normal-Inverse-Gamma prior on the tangent-space location *μC*​ (with prior dispersion scale *κ*0​\>0 around the global manifold centroid *μ*0​) and radial variance *σC*2​∼Inv-Gamma(*ν*0​/2,*ν*0​*s*02​/2) (where *s*02​\=*a*2 is the local squared Potts scale) and integrating out both (*μC*​,*σC*2​) combined with the CRP prior *α*Γ(∣*C*∣) yields the **Riemannian Student-*t* Log Marginal Likelihood**:

ΨRiem​(*C*)=CRP Partition Priorlog*α*\+logΓ(∣*C*∣)​​−2*d*M​​log(1+*κ*0​∣*C*∣​)\+logΓ(2*ν*0​​)Γ(2*νn*​(*C*)​)​\+2*ν*0​​log(*ν*0​*s*02​)−2*νn*​(*C*)​log(*ν*0​*s*02​\+*WC*​\+*κ*0​\+∣*C*∣*κ*0​∣*C*∣​*d*M2​(*μC*​,*μ*0​))

where *νn*​(*C*)=*ν*0​\+(∣*C*∣−1)*d*M​ are the posterior degrees of freedom.

### **Formulation 3: Unified Thermo-Riemannian Score (Combining Manifold Dispersion \+ Thermal Persistence)**

Why choose between Formulation 1 (which measures **non-linear topological persistence** along curved manifold filaments) and Formulation 2 (which measures **geodesic concentration** around a Fréchet center)? In statistical mechanics, the log prior weight of a macro-state *C* discovered by a thermal filtration is proportional to its **thermodynamic lifetime** *τ*pers​(*C*)/*T*max​. Replacing the flat CRP parameter *α* with the **Thermodynamic Persistence Prior** *α*thermo​(*C*)=*α*0​(*T*max​*τ*pers​(*C*)​)*β* yields the **Unified Thermo-Riemannian Score**:

Ψunified​(*C*)=ΨRiem​(*C*)+*β*log(*T*max​*τ*pers​(*C*)​)−*γ*cut​*J*cut​(∂*C*)

Look at how the three terms cooperate during the bottom-up recursion *V*∗(*CP*​)=max(Ψunified​(*CP*​),*V*∗(*CL*​)+*V*∗(*CR*​)):

> 1. **If two children *CL*​,*CR*​ have distinct Fréchet means** (*d*M​(*μCL*​​,*μCR*​​)≫*s*0​), merging them inflates *WCP*​​ by ∣*CL*​∣+∣*CR*​∣∣*CL*​∣∣*CR*​∣​*d*M2​(*μCL*​​,*μCR*​​), causing ΨRiem​(*CP*​) to drop sharply and **forcing the DP to keep *CL*​ and *CR*​ split**.  
> 2. **If a child *CL*​ is a transient thermal fluctuation** at a single temperature step (*τ*pers​(*CL*​)≈0), its persistence prior *β*log(*τ*pers​(*CL*​)/*T*max​)→−∞ **vetoes the spurious split** and locks *CP*​ intact.  
> 3. **If the entire dataset is homogeneous noise**, the root *V* (*K*\=1) wins over every child split, giving you an exact unsupervised test for whether any cluster structure exists at all.

## **Complete Bottom-Up Anti-Chain DP Implementation in C\#**

Once GrainHierarchyBuilder.BuildHierarchy constructs treeNodes\[0 .. 2M-2\], the bottom-up DP resolves the globally optimal cross-temperature partition Π∗ in a **single linear pass** from index 0 to 2*M*−2 (since Kruskal creates internal nodes *M*…2*M*−2 in strictly topological bottom-up order\!):

public static class MultiTemperatureAntiChainDp  
{  
    /// \<summary\>  
    /// Solves V\*(C) \= max(Psi(C), V\*(C\_L) \+ V\*(C\_R)) in O(M) time without recursion.  
    /// Relies on the invariant that Kruskal assigns parent IDs strictly greater than child IDs.  
    /// \</summary\>  
    public static double SolveOptimalAntiChain(  
        int M,  
        ReadOnlySpan\<GrainNode\> treeNodes,     // Length 2M \- 1  
        ReadOnlySpan\<double\> psiScore,         // Length 2M \- 1: Psi(C) for each node  
        Span\<double\> optimalValue,             // Length 2M \- 1: V\*(C)  
        Span\<bool\> lockGrainHere,              // Length 2M \- 1: true if C is selected in Pi^\*  
        List\<int\> selectedGrainIds)            // Output: node IDs forming the optimal partition Pi^\*  
    {  
        int totalNodes \= 2 \* M \- 1;

        // 1\. Bottom-Up Pass: Leaves (0..M-1) then Internal Nodes (M..2M-2) in increasing index order  
        for (int i \= 0; i \< M; i++)  
        {  
            optimalValue\[i\] \= psiScore\[i\];  
            lockGrainHere\[i\] \= true;  
        }

        for (int i \= M; i \< totalNodes; i++)  
        {  
            ref readonly GrainNode node \= ref treeNodes\[i\];  
            double splitScore \= optimalValue\[node.LeftChild\] \+ optimalValue\[node.RightChild\];  
            double mergeScore \= psiScore\[i\];

            if (mergeScore \>= splitScore)  
            {  
                optimalValue\[i\] \= mergeScore;  
                lockGrainHere\[i\] \= true;  // Lock grain C at its native birth temperature\!  
            }  
            else  
            {  
                optimalValue\[i\] \= splitScore;  
                lockGrainHere\[i\] \= false; // Defer to finer higher-temperature children  
            }  
        }

        // 2\. Top-Down Frontier Extraction from Root (2M \- 2\)  
        selectedGrainIds.Clear();  
        ExtractSelectedGrains(totalNodes \- 1, treeNodes, lockGrainHere, selectedGrainIds);

        return optimalValue\[totalNodes \- 1\];  
    }

    private static void ExtractSelectedGrains(  
        int nodeId,  
        ReadOnlySpan\<GrainNode\> treeNodes,  
        ReadOnlySpan\<bool\> lockGrainHere,  
        List\<int\> selectedGrainIds)  
    {  
        if (lockGrainHere\[nodeId\] \|\| treeNodes\[nodeId\].LeftChild \== \-1)  
        {  
            selectedGrainIds.Add(nodeId);  
            return;  
        }

        ExtractSelectedGrains(treeNodes\[nodeId\].LeftChild, treeNodes, lockGrainHere, selectedGrainIds);  
        ExtractSelectedGrains(treeNodes\[nodeId\].RightChild, treeNodes, lockGrainHere, selectedGrainIds);  
    }  
}

### **Summary of Computational Complexity**

| Pipeline Stage | Naive Multi-*T* SPC | Kruskal \+ Anti-Chain DP Architecture |
| :---- | :---- | :---- |
| **Storage across *R* temperatures** | O(*RM*2) dense matrices | O(∣E∣)=O(*kM*) sparse edge accumulators |
| **Thermal Persistence Integration** | O(*RM*2) | O(*RkM*) SIMD vector sum over edges |
| **Grain Hierarchy Construction** | O(*RM*2) thresholded BFS | O(*kM*log*M*) Kruskal MST pass |
| **Sufficient Statistics & Ψ(*C*)** | Re-computed from scratch | O(*M*⋅*D*) bottom-up moment addition |
| **Cross-Temperature Partition Π∗** | Impossible (single-*T* slice only) | O(*M*) exact bottom-up anti-chain DP |

---

