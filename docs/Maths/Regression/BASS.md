# BASS — Bayesian Adaptive Spline Surfaces

> Adaptive spline surfaces over continuous and categorical predictors, with scalar or direct functional responses, fitted by reversible jump MCMC over products of hinges. Part of `Maths.Regression`.

### BASS: terms and the prior

BASS, Bayesian adaptive spline surfaces (Francom et al. 2018; Francom and Sansó 2020), fits a surface over any number of predictors as a sum of products of hinges, by reversible jump over the collection of terms. It is being ported from R BASS 1.3.1 under the [port hub](../../../issues/regression/bass-full-port-brief.md). Its fitting decisions, D01–D06, were taken in the design pass of 2026-10-27, and the contract shapes and the [stage-4 plan](../../../issues/regression/bass/port-fitting-brief.md#stage-4-plan-decided-2026-10-28) on 2026-10-28; the [fitting segment](../../../issues/regression/bass/port-fitting-brief.md) records each with its derivation. This section and the seven after it describe stage 4 as built: S1–S6 here and in the next three, then the run (S7), tempering and the evidence (S8), the pointwise density and the residual check (S9), and the qualification (S10). Their tests and reviews are in the [checklist's reports](../../../issues/regression/bass/port-fitting-checklist.md#reports). BASS shares BARS's collapsed g-prior block and its reversible-jump engine; what it adds is the term space over many predictors, its prior and proposals, and a design that admits terms in any order.

**The references.** The checks below lean on R oracles under `r/oracles`, run by `r/scripts/bass_stage1_checks.R` and `bass_stage4_checks.R`: stage 1's reference fixtures F-O1–F-O9 (named in the fitting checklist), and the exact posterior over an enumerable structure space, τ integrated by quadrature (`bass_enumeration_oracle.R`'s `posterior` mode). Its fixtures are built by recipe, and the C# tests write their own data:
- *A:* two continuous predictors on {0, ½, 1}, N = 27;
- *B:* two on a four-point grid, N = 32;
- *C:* A's grid crossed with a three-level factor;
- *E:* A's inputs with the exact response 1 + 2(x₁ − ½)₊.

**Terms** (`BassTerm`, D01). A term has up to three parts, each over distinct variables, and at least one factor:
- *continuous:* hinges $(s(x-t))_+^q$ on continuous predictors, with a sign s = ±1, a knot t and the fit's degree q (1 by default, any positive integer);
- *categorical:* indicators 1(z ∈ D) of a proper, nonempty subset D of a factor's training levels;
- *functional:* hinges on the output variables of a direct functional response.

An empty part is the constant 1, and the intercept is not a term. Terms are held canonically, each part sorted by variable, and compare by value, so a model's terms form a set. A `HingeFactor(Variable, Sign, Knot)` indexes its variable's knots. A `CategoricalFactor` holds its subset as sorted level indices, so a factor of any number of levels is one value, with no ceiling at 64. `BassTermTable` interns terms, one id per distinct term in the order first seen, and retained draws hold ids: a storage choice, which leaves the measure unchanged.

**The basis** (`BassBasis`, R's `makeBasis` and `makeBasisCat`). A term evaluates as two factors:
- the case side $a(x,z)=\prod_j(s_j(x_j-t_j))_+^q/c\cdot\prod_l 1(z_l\in D_l)$, with c the hinge product's largest value on [0, 1], or 1 where that is 0 (a knot at the end its hinge faces);
- the output side h(t), the functional hinges normalized the same way.

A scalar response's column is the case side, over its n_x cases. A direct functional response observes each case at the same n_r output locations, N = n_x·n_r observations in all, and its column is $a(x_i)\,h(t_r)$ at case i and location r: a Khatri–Rao product that the design contracts without forming the $N\times M$ design.
- *Inputs.* `BassInputTransform` is fitted once on the training data: R's range scaling to [0, 1] for continuous predictors and output variables, over each column's declared domain when it declares one (`InputDomain`, [Polynomial chaos](../PolynomialChaos.md)) and its training range otherwise, keeping both ends and the training range, and each factor's training levels, the declared levels some case takes (the others are left out and logged). Each variable's knots are its distinct scaled values, so two raw values that scale to one double are one knot.
- *Refused at entry* (`BassData`): a constant continuous column, which has no knot; a factor taking fewer than two levels; a non-finite value; mismatched lengths; a name used twice, output variables included. The data's arrays are copied, so a fit never sees a value change.
- *Checked* (`BassBasisParityTests`): 60 random terms mixing every part, at degrees 1–3 on tied predictors, against R BASS's own functions at the pinned commit (`bass_basis_oracle.R`), within 10⁻¹⁴, and the scaled columns within 10⁻¹⁵. `BassBasisTests` checks a direct design's separable form against an explicit $N\times k$ design.

**The term measure** (`BassStructurePrior`, D01). A model is an unordered collection of M distinct terms. Its measure is the prior on labelled terms times the M! orderings, as in Nott, Kuk and Duc (2005), so a birth from M carries p(M+1)/p(M) once and no further factor of M + 1. A term's prior, for each part the fit has:
- the order K uniform: on {1, …, K_max}, with K_max the part's interaction cap, when the fit has one part; on {0, …, cap} for each part independently when it has several (decided 2026-10-28), the term whose parts are all empty being the intercept, outside the support;
- given K ≥ 1, the variable set uniform over the C(p, K) subsets of the part's p variables;
- each hinge's sign uniform and its knot uniform over its variable's distinct values, n_v of them for variable v;
- each categorical subset from the declared measure.

With one part and no support restriction this is Francom et al.'s constant $c_K=2^{-K}\binom pK^{-1}\prod_k n_{v_k}^{-1}$, times $1/K_{\max}$ for the order. With several parts, the all-empty combination's share of the measure is cut by the global truncation (below), like any term outside F, with no renormalization.
- *Knots on distinct values,* following Francom et al. 2018 ("the number of unique marginal values") and R; a functional variable's knots are the distinct values in its column of the output coordinates. *Not adopted:* Nott, Kuk and Duc's uniform over the n observations, which weights tied values by their multiplicity, and the continuous Unif(0, 1) of the JSS and GBMARS texts, which describes the range rather than the measure the code samples.
- *Categorical subsets.* R's measure is the default (`CategoricalSubsetMeasure.SizeUniform`): the size uniform on {1, …, n_v − 1}, then the subset uniform among those of that size, which mirrors the continuous rule (the order, then the set). Francom et al. 2019 state the uniform over the $2^{n_v}-2$ proper subsets, selectable as `Uniform`. The two agree up to three levels, which covers every factor in the paper's application; from four levels on, the paper's puts more mass on mid-sized subsets. Both are coherent, so this is a divergence between paper and code, not a correction.
- *Repeated terms lie outside the support.* An identical term duplicates a column, which the design refuses (below). Every admitted collection then has distinct members, which keeps the M! count exact.
- *Complement pairs.* With the intercept in every model, a categorical main effect on D and one on its complement span the same space, so the two models are one, and the measure counts it twice. It is a tilt to report, not a defect. A model holding both terms is refused, since their columns and the intercept are collinear.
- *R's scoring of the orders is a defect, not an option.* R samples each part's order on {0, …, cap}, but `logProbChangeMod` scores order 0 by its proposal weight, which moves with the counts, and each order J ≥ 1 at 1/cap, so the orders' masses sum past one ([Corrections](#corrections-to-the-sources)).

**The support F and the count prior** (D01, D02). The prior is restricted globally to the feasible set F: the basis cap M_max on the number of terms, each term's support (below), numerical admission (D06) and, under the Jeffreys σ² option, no exact fit (D03). There are no per-order or per-λ normalizers.
- *λ integrated out* (KHAOS). With λ ~ Gamma(h₁, rate h₂) and M | λ Poisson, one global normalizer keeps λ's integral closed-form. A collection of M terms with its orderings weighs $\Gamma(h_1+M)/(\Gamma(h_1)(h_2+1)^M)$ (`LogCollectionWeight`), and a birth's count factor is $(M+h_1)/((M+1)(h_2+1))$ (`LogCountRatio`). λ leaves the chain, and each retained draw gets λ | M ~ Gamma(h₁ + M, rate h₂ + 1), exact under global truncation: F-O3 matches its mixture over the count law to the prior tilted by F, to 2·10⁻¹⁴.
- *The defaults,* h₁ = h₂ = 10, are R's and Nott, Kuk and Duc's Gamma(10, 10). Francom et al. 2018 name λ's prior as the main control on overfitting.
- *Not adopted: the papers' truncation per λ.* JSS and Francom et al. 2018 write the Poisson truncated and renormalized for each λ. λ's update then needs the cap's normalizer, a regularized incomplete gamma function, which is tractable for the cap alone and intractable once support and admission join F. The global reading is what R and GBMARS sample (GBMARS calls its update approximate because its stated prior is the per-λ one), and one mechanism then covers the cap, support and admission. It changes the model the papers state.
- *Not adopted as the default: per-order normalization,* in which c_K is one over the number of admissible order-K terms, and for which Francom et al. recommend the unconstrained constant "as a conservative proxy". No implementation samples it, and its constant needs a data-dependent count. It stays computable on enumerable fixtures, to measure the tilt.
- *The tilts this leaves.* The order's marginal becomes π(K) ∝ P(admissible | K)/K_max. Products of hinges have more zeros and fail the support count more often, so high orders get less mass than under per-order normalization: on F-O3's fixtures the order law at M = 1 is 2/3 and 1/3 against the unconstrained ½ and ½, in proportion to the admitted prior mass per order (0.50 and 0.25 at p = 3). λ's marginal is its Gamma prior tilted by P_λ(F); at the default cap, min(1000, N − 1), the cap rarely binds. The run is to report the support-rejection rate per order and the posterior mass at the cap (S7).
- *KHAOS's manuscript* keeps both the death's 1/(M + 1) and the labelled prior's, which puts an extra 1/M! on the count law. It is a defect fixture (F-O6): 0.229 in total variation from the stated law on a two-predictor enumeration.

**Proposals** (`BassTermProposal`, D01; Nott, Kuk and Duc, extended by BASS to higher interaction orders and to several parts).
- *Birth.* For each part, its order from $I(K)\propto w_1+(\text{terms of order }K)$; for K = 1 its variable uniformly over the part's p variables; for K ≥ 2 its variables by weighted sampling without replacement, with $z_u\propto w_2+(\text{terms using }u)$; signs, knots and subsets from the prior. A draw that leaves every part empty is the intercept, a failed birth, as in R. The moves propose nothing and count it (`EmptyBirths`), where R counts a rejected birth, so the acceptance tallies differ and the law does not.
- *Death.* A term uniformly. Its reverse is the removed term's birth under the reduced structure's counts.
- *Change.* A term uniformly, then one of its factors uniformly, whose sign and knot, or subset, is redrawn from the prior. For a hinge the proposal is symmetric; for a subset its ratio is the redrawn subset's prior ratio, which cancels the target's. Either way the acceptance is the likelihood ratio alone. A change to the same term is the same state.
- *Within-move ratios,* which the chain adds to its target and selection ratios: for a term t and the current structure S, a birth's $\ln\tfrac1{M+1}-\ln q_b(t\mid S)$, and a death's $\ln q_b(t\mid S\setminus t)+\ln M$, with $q_b$ the birth's probability under a structure's counts. With the prior's ratio they make D02's count factor.
- *The pseudo-counts w₁ and w₂* (`OrderWeight` and `VariableWeight`, 5 each by default) are validated positive, never 0. A death's acceptance carries its reverse birth's probability: with w₁ = 0, an order used only by the removed term would have proposal probability 0 there, and with w₂ = 0 so would a variable used only by a removed interaction, so such deaths could never be accepted. Larger weights pull the births towards uniform orders and variables, smaller ones towards the current counts.
- *The variable set's probability* is Wallenius' noncentral hypergeometric with one ball per variable: for a set V, $P(V)=\int_0^1\prod_{i\in V}(1-t^{\alpha_i})\,dt$ with $\alpha_i=z_i/\sum_{k\notin V}z_k$ (`Maths.Distributions.Wallenius`).
  - *R's evaluation,* `dmwnchBass`, is inclusion–exclusion over the 2^K subsets, an alternating sum that loses digits under skewed weights.
  - *The port's* sums over draw orders instead, by the recursion $G(U)=\sum_{i\in U}z_iG(U\setminus i)/(\bar Z+Z_U)$ from G(∅) = 1, so that P(V) = G(V). Here $\bar Z$ is the weight outside V, Z_U the weight in U, and G(U) the probability that weighted draws from U and that outside weight take U's members first. Every term and denominator is a sum of positives, so it is exact to a few ulps at any weights, in O(K·2^K). Above K = 10 it takes adaptive quadrature of Chesson's positive integrand. D01 was amended to this on 2026-10-28 (`a2034ea`). The recursion costs 0.1–0.5 µs at K ≤ 5, against 40–85 µs for the quadrature.
  - *Checked* (`WalleniusParityTests`): every subset of orders 1–6 against F-O5's brute force over draw orders, within 1.1·10⁻¹⁴ relative, at BASS's default w₂ = 5 with skewed counts and at w₂ = 0.01 with counts in the hundreds. On the same subsets `dmwnchBass` errs by up to 2.2·10⁻⁹ and 2.8·10⁴ relative, at K = 5.
- *Exact discrete draws.* Knots, single variables, subset sizes and shuffles use Lemire's multiply-shift with rejection, which is exactly uniform. The generator's `NextInt` carries a bias of order n/2³² on a draw over n values, which knot draws over thousands of distinct values would feel in the law the acceptance ratio assumes.
- *Considered, not adopted* (the hub's [register](../../../issues/regression/bass-full-port-brief.md#what-the-lineage-taught)):
  - delayed rejection of structurally invalid births (KHAOS), which can condition on its rejected set because that set's probability is computable; BASS's births fail on data-dependent support, whose probability is not, so a failed birth counts as a rejection;
  - a uniform prior over basis functions (KHAOS); BASS's prior, uniform over the order and then over that order's terms, already corrects for how many terms each order has.

  KHAOS's Bernoulli-inclusion variable proposal and its variable swap are extensions.
- *Checked* (`BassProposalTests`, on six term spaces of up to 8,294 terms, which the test enumerates from the definitions): the prior sums to one with the empty combination's mass (2·10⁻¹³); the birth proposal sums to one at count states up to M = 7 (1.9·10⁻¹³); 400,000 births per space match their stated probabilities by Pearson's χ² on 490 to 8,295 cells, every z below 1.6; every death's ratio reverses its birth's within 10⁻¹³; and the change draws and ratios match the exact law.

**Move schedules** (D01). The `ReversibleJumpChain` owns the state-dependent selection ratio, so a schedule is a set of move weights at M.
- *Default: equal over the moves available* (Nott, Kuk and Duc; `EqualOverAvailable`): birth only at M = 0, death and change ½ each at M_max, ⅓ each between.
- *Not ported: R's fixed thirds.* R selects its moves by this rule but puts ⅓ in every ratio. Under a constant likelihood that makes π(1)/π(0) three times the target's and π(M_max)/π(M_max − 1) two-thirds of it (F-O4), and it puts R's own chain 0.020 off its posterior in total variation on F-O6's enumerable space (0.100 at three predictors). Here the thirds have nowhere to enter.
- *Option: prior-cancelling* (`PriorCancelling(C = 0.4)`, C ∈ (0, ½]; Green 1995, as used by BARS's DMGK schedule and by He, Yang and Kang). It sets $b_M=C\min(1,p(M+1)/p(M))$ and $d_M=C\min(1,p(M-1)/p(M))$, with $p(M+1)/p(M)=(M+h_1)/((M+1)(h_2+1))$, and gives the change $1-b_M-d_M$ even at M = 0, where it proposes nothing, so the weights sum to one and the count prior cancels at the boundary too. Its trade-off is BARS's ([DMGK schedule](#bars-target-and-moves)): there the winner depended on the direction of prior–likelihood conflict, so it waited for a pre-registered comparison. F-M2's pass ran it on 2026-10-29, with F-M1 folded in: it did not clear (score 0.864, the mixed fit at 0.793), so equal over available stays the default.
- *Checked* (`BassSweepTests`). Under a constant likelihood, on both schedules, four chains of 60,000 sweeps against the prior on F, which the test enumerates through the design (2,263 states; count law 0.6037, 0.3049, 0.0781, 0.0132): the count law within five batch-means errors, π(1)/π(0) within 5% where R's thirds would triple it, the order law, τ's prior and λ's retained mean. With data, the prior-cancelling schedule matches fixture A's exact posterior at M_max = 2 (`BassChainParityTests`).

### BASS: the design and admission

`BassDesign` (one per chain; D06) holds a structure's columns, the factor of their Gram and the least-squares fit, through births, deaths and changes ([`BorderedCholesky`](../LinearAlgebra.md#bordered-cholesky-borderedcholesky)'s proposals; a change's new column goes last).

**The centred basis.** Under every intercept convention the collapsed marginal depends only on the column space of [1, Z], so the design holds Z_c, the columns centred over all N rows, and ỹ = y − ȳ. The intercept is orthogonal to Z_c, and its Gram entry N stands apart: the Gram is $N\oplus G_c$, with $G_c=Z_c^\top Z_c$. The design's statistics are the least-squares fit's residual ‖r‖² and its fitted sum $Q_c=\tilde y^\top Z_cG_c^{-1}Z_c^\top\tilde y$; each convention adds its own pieces in the target (the zero centre's Nȳ², [below](#bass-the-collapsed-target)).
- *Scalar fits* hold explicitly centred columns. *Direct fits* hold the case-side and output-side factors and form Gram entries by the centred Hadamard form $G_c=(A_c^\top A_c)\circ(H_c^\top H_c)+n_x(\bar a\bar a^\top)\circ(H_c^\top H_c)+n_r(\bar h\bar h^\top)\circ(A_c^\top A_c)$ (`KhatriRao.CentredGram`, `CentredBorder`). It is a sum of three positive semidefinite terms, so the diagonal forms without cancellation, and it keeps the separable structure. Every inner product is a pairwise sum ([`PairwiseSum`](../LinearAlgebra.md#pairwise-summation-pairwisesum), N10).
- *Not used: the rank-one downdate* $B^\top B-N\bar b\bar b^\top$ of the uncentred Gram. For a column whose mean is large against its spread it cancels: the diagonal entry loses about ε/ρ_j relative, with $\rho_j=\lVert b_j-\bar b_j\mathbf 1\rVert^2/\lVert b_j\rVert^2$. On F-O7's indicator true for all but one case in 10⁵ (ρ = 10⁻⁵), s came out with a relative error of 3·10⁻³ by the downdate against 3·10⁻⁷ explicitly centred. On a direct design, the Hadamard form matched explicit centring to 2·10⁻¹⁶ in the Gram, where the downdate gave 4·10⁻¹³.
- *Pairwise sums are part of the decision.* Summed left to right, the Gram of F-O7's designs with a near-constant column gave s ≈ 10⁻⁹ with relative errors of 7·10⁻⁴ to 2·10⁻³, against about 2·10⁻⁷ from ε times the equilibrated condition number: the column's many equal centred values, added after one large one, round the same way every time. Summed pairwise, the errors were 2·10⁻⁷ to 3·10⁻⁷.
- *The design pass's correction* (2026-10-27). The proper-defaults draft claimed that level centring would let the gate read R's uncentred Gram and spare the downdate. The column-space argument makes the centred basis exact under every convention, so admission does not depend on the intercept choice.
- *Cost.* `KhatriRao.CentredBorder` re-centres every current factor on each call, so a direct proposal allocates O(M(n_x + n_r)): correct, and a candidate for S10's cost baseline.

**The numerical gate** (D06, in Q6's two-gate form, as for BARS). A structure is admitted only when every column's angle to the span of all the others and the intercept satisfies $s_j=1/(G_{jj}(G^{-1})_{jj})\ge s_{\min}$, read from the centred Gram. So s_j is the reciprocal of column j's centred variance inflation factor, and the intercept's own s is 1. The floor is $s_{\min}=10^{-10}$ (`AngleFloor`), BARS's `MinRelativePivot` threshold.
- *Why order-free.* A newest-column angle depends on insertion order. With b₁ = e₁, b₂ = e₁ + 10⁻³e₂ and b₃ = e₂ + 10⁻³e₃, the leave-one-out values are (10⁻¹², 10⁻¹², 10⁻⁶): added last, b₃ shows 10⁻⁶ and b₁ shows 10⁻¹². Under a newest-column rule at 10⁻¹⁰, the same set would be admitted or refused by its history, so F would not be a set and the target restricted to it would be undefined. BARS's columns come in a canonical order (sorted knots), so its `MinRelativePivot` is already a property of the set; BASS's terms are unordered.
  - On three hinges at neighbouring knots (`BassDesignTests`) the set's smallest leave-one-out angle is 5.5·10⁻⁵, while a newest-column rule sees 1.6·10⁻³ when h₂ follows h₁ and 2.2·10⁻⁴ when h₃ comes last. At a floor of 1.1·10⁻⁴ all six orders refuse the third birth, where a newest-column rule would admit the set in four of them.
- *Trade-off.* The gate reads every column's angle, not the newest one's. The bordered inverse keeps diag(G⁻¹) in O(M²) per birth or death. The gate is downward closed, since removing a column only widens the others' angles, so a death never fails it. It is at least as strict as a pivot gate in any order.
- *Centred against uncentred.* Both angles start from the same residual, b_j projected off the other columns and 1, and divide it by ‖b_j − b̄_j1‖² or by ‖b_j‖², so uncentred s_j = centred s_j·ρ_j (F-O7, to 3·10⁻¹⁵). The uncentred gate is stricter for columns whose mean is large against their spread: a full hinge on uniform data has ρ = 1/4, and an indicator true for all but a few cases has ρ near zero. F-O7's indicator has a centred s of 10⁻⁹, which the gate admits, and an uncentred s of 10⁻¹⁴, which an uncentred gate would refuse. The marginal reads only the column space, so the centred angle is the one that measures the model.
- *Strict confirmation.* When the bordered factor fails the gate or comes within a factor of two of s_min, a strict factorization of the candidate's Gram decides, and an admitted candidate is accepted from it. A computed s_j carries a relative error of about u_G·κ(E), with u_G the relative error made in forming the Gram and κ(E) the equilibrated condition number, and the band covers it for designs that pass. At a binding floor of 3·10⁻³, `BassDesignTests` saw 24 admitted candidates within the band, decided by the strict factorization.
- *Statistics.* 1/s_min = 10¹⁰ bounds every coefficient's variance inflation. It is a numerical guard, not a modelling restriction: the marginal's Occam factor already discourages redundant terms. *Option:* `MaxVarianceInflation`, a declared cap, a statistical support restriction for stable extrapolation, off by default.
- *Not ported: R's `getQf` rule,* a pivot-ratio gate on the uncentred factor, refusing when max_{j≥2} r_jj / min_j r_jj exceeds 10³ (commented "a hack, otherwise we get huge variance inflation in beta"). It depends on column norms, so it refuses well-conditioned local terms whose columns are small, and more often as n grows. At n = 40,000, a 20-point local term with values 0.01 gives R's ratio 1291, refused, while the design's smallest leave-one-out angle is 0.25. It is a defect fixture (`bass_admission_oracle.R`).

**Support per part** (D01, D06). A term needs a threshold of nonzero values in its case side and, in a direct fit, in its output side, counted apart (R's `npart` and `npart.func`). By default the threshold is min(20, 0.1·n), with n = n_x for the case side and n_r for the output side. It is a prior support restriction inside F, checked on the uncentred factors, not a numerical test.
- *The count, the default:* every BASS source uses it, and Francom et al. 2018 justify it by the edge instability of hinge terms with few nonzero values.
- *Option: column mass* (`SupportKind.ColumnMass`), Σᵢ B(xᵢ) with hinges normalized to a maximum of one, BARS's Q6 rule. It weighs each nonzero value by its size, so a hinge whose nonzero values lie near its knot carries less mass than its count.

**Refusals.** A birth or a change is checked for the cap, a repeated term (a duplicated column, s = 0), the support in each part, the angle gate and the declared variance-inflation cap, in that order; a death needs no gate. A refused candidate carries its reason (`BassRefusal`: `Cap`, `Support`, `Numerical`, `VarianceInflation`, and the target's `ExactFit`), the target scores it off the support, and the moves tally it (`Tally<BassRefusal>`).
- *A ridge is a prior, never a numerical repair.* Under the g-prior a rank-deficient design lies outside F, and a ridge enters only as a declared coefficient prior (D03). ABMRS's code instead drops thin columns silently and adds a 10⁻⁸ ridge; the register corrects that by refusal, as the He, Yang and Kang row of [Corrections](#corrections-to-the-sources) records for BARS.

**The residual.** A candidate's ‖r′‖² updates the current residual's: a birth removes (zᵀr)²/δ, with z the new centred column and δ its pivot; a death adds $\hat\beta_m^2/(G^{-1})_{mm}$; a change does both, through the reduced model. These lose digits only when one move explains nearly all the residual (below 10⁻⁸ of it), and there the candidate's residual is formed explicitly. Taking ‖ỹ‖² − Q′ instead would lose them on every good fit.
- *On every acceptance* β̂ is solved and the residual recomputed, r = ỹ − Z_cβ̂, so it cannot drift; the factor still rebuilds strictly every 100 accepted updates. D06 asked for the residual updated on births and recomputed after deaths and on a cadence. Recomputing on every acceptance is at least as accurate and no costlier, since acceptances are rare next to proposals.
- *Checked* (`BassDesignTests`, against a reference that shares nothing with the bordered factor: the explicitly centred design, its least squares and its angles from a thin SVD):
  - random sequences of 400 proposals on scalar and direct data, at s_min = 10⁻¹⁰ and at a binding 3·10⁻³: every decision equal except within 5% of the floor, ‖r′‖² and Q′ within 10⁻⁸, the angles within 10⁻⁶;
  - 1,500 accepted moves on direct data: the residual within 10⁻⁹, and a rebuild finding a factor drift of 3·10⁻¹⁵ at M = 11;
  - one constructed case per refusal reason.

  `BassAdmissionParityTests` checks the angles against the oracle's QR of the explicitly centred design: to 5·10⁻¹² on scalar terms, to 3·10⁻¹³ on direct terms by the three-term form, and to 3·10⁻⁹ on F-O7's hard case, at an angle of 10⁻⁶.

### BASS: the collapsed target

`BassTarget` (D03, D05) is the posterior over (structure, τ), with the coefficients and σ² integrated out by the shared block's centred form ([above](#the-collapsed-g-prior-block); `CollapsedGPrior.CentredMarginalAt`), an O(1) query on the design's ‖r‖², Q_c, M, N and ȳ; below, Q is Q_c with the convention's additions. A structure's account, `BassStructureFit`, carries those statistics and its terms' ln π, or the design's refusal. Under the g-form, the default ([the ridge](#bass-the-coefficient-priors-two-axes-and-the-ridge) is the other form),

$$y\mid a,\sigma^2\sim N(Ba,\sigma^2I),\qquad a\mid\sigma^2,\tau\sim N(\mu_0,\sigma^2(B^\top B)^{-1}/\tau),$$

with B the intercept and the terms' columns, τ the precision multiplier as in R, JSS and Francom et al., and μ₀ the intercept's centre. That is the form under level centring and the zero centre; the Liang convention takes the intercept out of the g-prior (below).

**Proper defaults** (D03, decided 2026-10-27). σ² is scaled to the response and the intercept level-centred, so every default prior is proper: σ², the intercept, τ under Zellner–Siow, and D02's count prior. For the intercept, shift invariance, a proper prior and no double use of the data cannot all hold, since a proper prior on a location is shift-invariant only if it reads the data:

| Convention | Shift-invariant | Proper | No double use of the data |
|---|---|---|---|
| Liang et al. 2008, a flat intercept (`GPriorCentre.Flat`) | yes | no | yes |
| Level centring (`GPriorCentre.Mean`), the default | yes | yes | no: about 0.23 nats of extra credit per added term |
| R's zero centre (`GPriorCentre.Zero`) | no: up to 9 nats across offsets on F-O9's design | yes | yes |

- *Why level centring.* It is proper; it gives BARS and BASS one convention through the shared block; and it is R's convention applied to y − ȳ, which R's own code advises (a commented-out stop in its σ² update). The intercept stays a column inside the g-prior, so k = M + 1, centred on the level, μ₀ = ȳ·e₀. With the mean, the level's own energy is exactly zero, so τ's conditional reads no level at all. Liang with the scaled σ² prior was the coherent runner-up.
- *What it costs against Liang.* It reads the level from the data at no charge, so σ² gets N degrees of freedom, not N − 1. Under the Jeffreys limit its marginal is Liang's with −N/2 in place of −(N − 1)/2 on (1 + g(1 − R²)), g = 1/τ (F-O9, to 10⁻⁸). Against the intercept-only structure it favours a structure by 0.43–0.47 nats more than Liang at M = 1 and by 1.06–1.08 at M = 5, and between two structures one term apart by 0.23–0.24, at N = 40, 100 and 1000. A term's Occam factor is about ½·ln(1 + g), 1.8–3.5 nats at these N, so the extra credit is a small fraction of it, but systematic and not shrinking with N.
- *What it buys.* An evidence wherever a run can give one (a ladder ending at the reference rung, S8, or the SMC normalizer, stage 8), empirical-Bayes in the level and the scale and so comparable only on the same data; proper exact fits; the exact reference rung and SMC from the prior with no separate mode; and no raw-unit prior among the defaults (OQ19 settled for BASS).
- *The Liang convention* integrates a flat intercept out, with k = M. It gives up the absolute marginal: its values hold up to the flat prior's constant, which cancels from every within-rung ratio and every swap, so it reports no evidence, and on the integrated tempered form (below) it loses b = 0.
- *The zero centre* is DiMatteo, Genovese and Kass's prior, selectable on its own terms, for a response whose level is known to be zero; k = M + 1, and the intercept's fitted energy Nȳ² joins Q. Under it, a y far from zero against its noise pushes τ small, and a small τ makes every term's Occam factor, about (τ/(1+τ))^{1/2}, harsher, so fits come out sparser than the data warrant. On F-O9's 40-point design the two-hinge structure's Bayes factor goes from 5.56 to 2.16 and −3.87 nats at offsets 0, 10 and −250, where level centring and Liang do not move.
- *The median centre* is refused: it is for the noise extensions' weights. For supplied-basis and PCA fits the coefficient responses already have zero mean across cases, so the conventions nearly coincide there; the choice bites for scalar and direct fits.
- *The σ² prior,* `Sigma2Prior.ScaledToResponse(10⁻³)`: IG(c, c·s²), with c = 10⁻³ and s² the response's mean square about the centre's level, as Normal BARS's default. It is scale-free: multiplying y by a constant moves every structure's marginal by exactly −N times its log (checked at 10⁻⁶, 10³ and 10⁶). At c = 10⁻³ the differences between structures sit within 0.05 nats of the Jeffreys prior's on BARS's sine fit. The Jeffreys prior 1/σ², R's default, stays selectable as the reference prior that reads nothing from y, and so does a declared inverse gamma, a prior in y's units, under which rescaling y moves the structures unequally.

**τ's prior under the g-form** (D03): the unit-information Gamma, τ ~ Gamma(a, rate a·N), whose mean is 1/N for every shape a. The default a = ½ is Zellner–Siow, equivalently g = 1/τ ~ IG(½, N/2) (Liang et al. 2008). The shape is a declared constant (`TauPrior.UnitInformationGamma(a)`), never a function of N: E[ln τ] − ln(1/N) = ψ(a) − ln a at every N, −1.27 at a = ½. A larger a tightens τ towards unit information, with the held τ₀ as the limit a → ∞. The family's other members, hyper-g and hyper-g/n, are [below](#bass-the-coefficient-priors-two-axes-and-the-ridge). N counts observations: N = n_x·n_r for a direct functional response. The sources:

| Source | Prior as written | Stated intent | E[τ] as written |
|---|---|---|---|
| Francom et al. 2018 | Gamma(1/N, 1) | "centered over the unit information prior" | 1/N |
| R documentation (`b.tau`) | rate "one over the number of observations", shape 1/2 | "centers … on the unit information prior" | n/2 |
| JSS 2020 | Gamma(1/2, 2/n), "(shape and rate)" | "the Zellner–Siow Cauchy prior" | n/4 |
| GBMARS | g ~ IG(1/2, 2/n) on the variance multiplier | "leads to the Zellner–Siow Cauchy prior" | n/4 |
| R code (the default) | Gamma(1/2, rate 2/n) | none | n/4 |

Every stated intent centres τ on unit information, and two sources name Zellner–Siow, which the written 2/n gives only when read as a scale.
- *Francom et al. 2018's Gamma(1/N, 1) is the member a = 1/N, corrected.* A shape that shrinks with N keeps the mean at 1/N but spreads ln τ over about N nats: E[ln τ] = ψ(1/N) ≈ −N and sd(ln τ) ≈ N, so the median is about 2⁻ᴺ (e⁻²⁸ at N = 40). Under the Liang convention, leaving the intercept-only structure costs 1.4–5.1 nats more than under Zellner–Siow, growing with N, while between two structures that both carry τ the two priors differ by under 0.3 nats (F-O8).
- *R's code rate 2/n is not ported.* It contradicts every stated intent, R's documentation included, and it would buy no parity: τ's prior does not enter the deterministic statistics, which are compared at a fixed τ.
- *On u = ln τ.* The target holds τ's prior as a density on u, the coordinate τ's slice step moves (D04), so the powered path powers τ in log τ, D05's recorded choice. Where eᵘ leaves the doubles, τ is off the support, a prior mass below e⁻³⁷⁰.

**Held τ** (`TauPrior.Held()`, D03's second mode, for the completion contract). τ is held at a declared value, under the g-form unit information 1/N by default, which is BARS's convention and the family's limit a → ∞. The state then has no τ coordinate and the sweep no τ block. It keeps the practical intent of R's Jeffreys mode: the coefficients are essentially unshrunk (by the factor N/(N + 1)), and structures are selected by the unit-information Occam factor, about ½·ln(N + 1) per term, the BIC penalty. It is coherent: the collapsed marginal with no τ block.
- *R's Jeffreys mode is not ported;* it is a defect fixture (F-O6). It holds τ at 0, so the coefficient prior is flat, and its moves and its σ² draw imply different coefficient priors, σ^{−2k}|G|^{3/4} (|G| the Gram's determinant) for the moves and σ^{−(k−1)} for the σ² draw (confirmed by quadrature), so it has no single target. The root cause is general: an improper coefficient prior's constant changes with the dimension (Bartlett's paradox), so any such mode picks an arbitrary constant per term.
- *Not adopted:* a flat prior with a declared per-term constant, coherent but an information-criterion penalty in disguise.

**Propriety and the exact-fit rule** (D03). The structure space is finite, so the posterior is proper exactly when every admitted structure's marginal is finite.
- *Under the Jeffreys σ² prior* and Zellner–Siow τ, integrating a and σ² leaves $\int p(\tau)\,(\tau/(1+\tau))^{k/2}\,(S(\tau)/2)^{-N/2}\,d\tau$, with $S(\tau)=\lVert r\rVert^2+\tau Q/(1+\tau)$. It is finite whenever ‖r‖² > 0. When a structure fits y exactly, the integrand behaves as τ^{(k−N−1)/2} near zero, so it diverges for k < N. On F-O1's 40-point fixture with k = 3 and an exact fit, the integral grew by 18 per unit of ln(1/τ_lower) as its lower limit τ_lower moved down, the predicted (N − k − 1)/2, and converged once g₂ = 10⁻³. The Liang convention behaves the same in g.
- *A proper σ² prior* (g₂ > 0) replaces S by 2g₂ + S > 0, so every marginal is finite, and the default needs no rule.
- *The rule, under the Jeffreys option only:* F excludes a structure whose residual is zero to working precision, ‖r‖ ≤ 10³·ε·‖y − level‖ from the explicit residual, refused as `BassRefusal.ExactFit`. Deterministic simulators with piecewise-linear structure produce such fits: two opposite hinges at one knot span x exactly. A Jeffreys chain on fixture E's exact response, over a larger design, refused 246 exact-fit candidates by reason and never held one, while the scaled default visits those fits; among fixture E's 299 states, Jeffreys refuses exactly the oracle's 27.
- *The response's own refusals* (`BassResolved`): a zero response under the zero centre; a constant one under level centring or the flat intercept, which the intercept absorbs whole; and a constant one under the zero centre with Jeffreys, which the intercept alone fits exactly.

**The tempered forms** (D05). The default path is likelihood-tempered, π_b ∝ prior(structure, τ)·p_b(y | structure, τ): the count prior, the term measure and τ's prior stay untempered at every rung, so the prior factors cancel from swaps, and the hot rungs approach the actual prior. The target gives the likelihood in two forms (`TemperedForm`), each with its slope in b, for path sampling:
- *integrated,* the default: p_b, the coefficients and σ² integrated against the tempered likelihood, through the block's $S=b(\lVert r\rVert^2+\tfrac{\tau}{\tau+b}Q)$ and determinant term. Its slope comes from $\partial S/\partial b=\lVert r\rVert^2+(\tau/(\tau+b))^2Q$ and the determinant term's slope (`GPriorMarginal.LogMarginalSlope`). At b = 0 under the defaults its log likelihood is 0, so the rung is the prior; its slope there is infinite in expectation under Zellner–Siow, since it carries −(k/2)·E[1/τ].
- *linear,* selectable: b·ln p₁, the collapsed marginal powered as a whole, whose slope is ln p₁ itself.

`Defines(b)` says where each form exists.
- *The linear form* is defined on all of [0, 1]. Its b = 0 member is the prior over (structure, τ), which stays proper under the improper options, so they keep an exact reference rung and an SMC start, though no evidence.
- *The integrated form* is defined on [0, 1] under the proper defaults, for b > 0 under Jeffreys, and where g₁ + (Nb − 1)/2 > 0 under the flat intercept, whose α integral spends one observation. Outside that, an evaluation throws.
- *Which to use.* Under the proper defaults the choice is efficiency alone, the ground on which BARS kept the integrated form: its linear form's edge, about 1.2×, fell short of the pre-registered threshold.
- *The powered path* stays selectable (D05). Its b → 0 limit is its reference measure, which is not normalizable here, so it has no exact reference rung; and since BASS's knots sit at data values, the counting measure over structures concentrates at the cap.
- *The ladder* that uses these forms is S8's, not built yet.

**Checked.** `BassTargetTests`, and `BassTargetParityTests` against the posterior mode's fixtures. Each state is loaded into a design from its emitted terms, and the marginal is integrated over u by the test's own Gauss–Legendre rule.
- *Fixture A at M_max = 2* (299 states) under each of six options (level centring, flat, zero, Jeffreys, flat with Jeffreys, held τ): ‖r‖² to 5·10⁻¹⁶ and Q to 3·10⁻¹⁵ of ‖y − level‖², the structure prior to 2·10⁻¹⁴, the marginal to 7·10⁻¹⁴.
- *Fixture C* (1,949 states): the mixed order prior to 9·10⁻¹⁴.
- *Fixture B's rungs,* b = 0, 0.1, 0.5 and 1 on both forms: log Z_b to 10⁻⁹, and the slope's mean wherever it is finite.
- *R's own statistics* under its zero centre at τ = 0.05 (`bass_statistics_oracle.R`): qf, yᵀy and ‖r‖² to 2·10⁻¹⁶ relative and S/2 to 4·10⁻¹⁶, the parity where R is right.
- *Invariances:* offsets of 0, 10 and −250 move no marginal under level centring or Liang (10⁻⁹), while the zero centre's Bayes factors move by more than a nat; and the scale invariance above.
- *τ's density on u* integrates to one, with E[τ] = 1/N and E[ln τ] = ψ(a) − ln(aN) at a = ½ and 2; each form's slope matches central differences.

### BASS: the coefficient prior's two axes, and the ridge

`BassCoefficientPrior(Form, Scale)` (D03; its two axes decided on 2026-10-29) holds a covariance form and τ's prior, so every form takes every member of one scale family.
- *The forms* (`CoefficientForm`):
  - `GForm`, a | σ², τ ~ N(μ₀, σ²G_c⁻¹/τ), the default ([above](#bass-the-collapsed-target));
  - `RidgeForm`, the normalized ridge, β | σ², τ ~ N(0, σ²I/τ) on the centred columns scaled to unit root mean square.
  - `BassCoefficientPrior.G(scale)` and `Ridge(scale)` build them.
- *The scale* is the shared family `TauPrior` (`Maths.Regression.GPrior`, beside `Sigma2Prior`). Its members are stated on the form's unit-information precision τ₀ = 1/n₀: n₀ = N under the g-form, 1 under the ridge (`BassCoefficientPrior.UnitObservations`; `TauPrior.For(n₀)` → `ResolvedTauPrior`).

| Member | Law | Proper | Note |
|---|---|---|---|
| `UnitInformationGamma(a)` | τ ~ Gamma(a, rate a/τ₀), mean τ₀ | yes | a = ½ is Zellner–Siow, the default |
| `HyperG(a)` | on g = 1/τ, p(g) = ((a − 2)/2)(1 + g)^(−a/2) | for a > 2 | Liang, Paulo, Molina, Clyde and Berger 2008; stated on g, not on τ₀ |
| `HyperGOverN(a)` | the same law on τ₀/τ | for a > 2 | Liang et al.'s hyper-g/n; on the ridge, where τ₀ = 1, it is hyper-g |
| `Held(τ)` | τ held | yes | null is τ₀ under the g-form; a declared value is required under the ridge |

- *Each member* gives its density on u = ln τ with the Jacobian, an exact draw where it is proper, and its propriety.
- *The shape, as built.* The decision named the scale as τ sampled under a declared prior, or held. The held value is a member of the family, so the scale axis is the family itself, with no second `Held` inside a sampled wrapper.
- *Improper hyper-g* (1 < a ≤ 2) is accepted where τ's posterior stays proper.
  - There every structure keeps a column inside the prior, so τ's likelihood vanishes as τ^(k/2) at 0, and the posterior on u converges for a > 2 − k.
  - It is refused under the flat intercept, whose intercept-only structure has k = 0, and wherever a proper prior is needed: a b = 0 rung, SMC from the prior and the evidence.
  - Such a run reports no evidence, and a chain starts at ln τ₀.
  - At a ≤ 1 even k = 1 leaves τ no posterior, so the constructor refuses it.

**The normalized ridge.** With D = diag(G_c)/N, the ridge on unit-RMS columns is the generalized ridge τD on the raw centred columns.
- *Why unit RMS* (D03). τ = 1 is then one observation's worth for every term, whatever its knot. On raw columns a hinge's squared norm varies as about n(1 − t)³/3 with its knot t, so one τ would shrink terms by where their knots sit.
- *The marginal* (`CollapsedGPrior.CentredRidgeMarginalAt`, `CollapsedRidgePrior`). At inverse temperature b the precision is bG_c + τD. With G̃ = D^(−½)G_cD^(−½) = VΛVᵀ, whose diagonal is N, w = VᵀD^(−½)Z_cᵀỹ and qⱼ = wⱼ²/λⱼ:

$$\tfrac12\ln\frac{|\tau D|}{|bG_c+\tau D|}=\tfrac12\sum_j\ln\frac{\tau}{b\lambda_j+\tau},\qquad S=b\Big(\lVert r\rVert^2+\sum_j q_j\,\frac{\tau}{b\lambda_j+\tau}\Big).$$

  - σ² is integrated against its declared prior, as under the g-form.
  - S is a sum of non-negative terms: the exponent's b²·Σⱼwⱼ²/(bλⱼ + τ), rewritten without cancellation.
  - The terms read the columns only through G̃ and w, so the marginal does not depend on the columns' scale.
  - On orthogonal columns every λⱼ is N, and the ridge is the g-form at τ/N.
- *The intercept,* a unit-RMS column already, takes the same prior under level centring and the zero centre: α | σ² ~ N(level, σ²/τ), one more direction with λ₀ = N, counted in k. Under Liang's flat intercept it stays flat. Both forms then give the intercept unit information at τ₀.
- *Where the forms differ* (read from the formulas, not measured).
  - The g-form's Occam factor is ½ln(1 + b/τ) per term whatever the columns' correlation, and its marginal depends on the column space alone.
  - The ridge's is ½ln(1 + bλⱼ/τ) per eigen-direction. A term nearly collinear with the current ones adds a small λ, costs little, and is shrunk hard. So the ridge charges a term by the new direction it brings, and its marginal depends on the terms themselves, not only on the space they span.
  - Which regime favours which was F-R1's pre-registered comparison. KHAOS found its fixed ridge best on noise-free simulators, and its g-prior better on most noisy ones. F-R1 ran on 2026-10-31 and read no regime: at S10's lengths every untempered fit failed its convergence gate (predictive R̂ above 1.05, median 2.14 noise-free), so no benchmark was scored and both forms stay options.
- *The fast path* (chosen by measured cost). A design under the ridge keeps its current structure's spectrum, decomposed lazily once per structure (cyclic Jacobi on G̃/N, `DenseEigen`).
  - A birth, death or change is evaluated through that spectrum's eigenvectors (`RidgeSpectrum.Birth`, `Death`, `Change`). The parent's precision is Schur-reduced by the removed column, then bordered by the new one. That costs O(M²) to start and O(M) per (τ, b), with no factor of its own and no refactorization when τ moves.
  - Where admission's strict factorization sits within a factor of two of s_min, the candidate is decomposed fresh.
  - Measured per proposal at N = 500: 0.7–17 µs bordered for M = 5–80, at or below the design's own border products (4–19 µs), against 8 µs–5.6 ms fresh. On S10's Friedman, 3,500 sweeps took 0.31 s against the g-form's 0.25 s.
- *What else the ridge changes:*
  - the coefficient draws, β | σ², τ ~ N((bG_c + τD)⁻¹bZ_cᵀỹ, σ²(bG_c + τD)⁻¹), mapped to the uncentred basis as the g-form's are;
  - σ²'s conditional, with the ridge's S;
  - the collapsed pointwise content, a Student-t per structural draw;
  - the evidence's record: `BassEvidencePriors` carries the form, the declared scale and the resolved one with its τ₀;
  - the prior sampler and the reference rung, where τ is drawn from its prior.
- *Admission is unchanged.* D06's gates read the design's strict factorization under both forms, so F is the same, and priors are compared on one model space.
- *A held ridge's hot end* auto-picks min(τ/N, 1/N), the b at which the held τ halves a unit-RMS column's shrinkage τ/(τ + bN). Under the g-form it stays min(τ, 1/N).
- *The default is unchanged* and draws bit-identically: the g-form with Zellner–Siow. Every configuration's fingerprint moved with the new serialization.

**Checked** (`TauPriorTests`, `BassCoefficientPriorTests`, `BassRidgeTests`, `BassRidgeDrawTests`, and the oracle classes under the ridge; [the chip's report](../../../issues/regression/bass/ridge-scale-prior-brief.md#reports)).
- *The family:*
  - each density normalizes on u and matches its law on τ (or on g′) with the Jacobian, to 10⁻¹²;
  - E[τ] = τ₀ and E[ln τ] = ψ(a) − ln(a·n₀), and hyper-g's Beta moments, to 10⁻¹¹;
  - draws within 5 SE;
  - hyper-g's Bayes factor against Liang et al.'s eq. (17), to 10⁻⁹;
  - R's τ-prior oracle takes both hyper-g members.
- *τ's posterior on a fixed structure* against quadrature under both forms, within 5 batch SE: the Gamma at a = ½ and 2, hyper-g, hyper-g/n and improper hyper-g at 1.5.
- *The marginal:*
  - against a dense factorization under every centre and σ² option, at M = 0, 1 and 3 and b ∈ {1, 0.3, 0.02}, to 10⁻¹¹;
  - under rescaled columns (×10⁻³ to ×10⁴), to 10⁻¹¹;
  - against the g-form at τ/N on orthogonal columns, to 10⁻¹².
- *Bordered against fresh* over 600-step randomized move sequences, close knots at s_min = 10⁻¹⁰ included: within 10⁻⁸, worst seen 4·10⁻¹⁰.
- *R's enumeration oracle* (`form=ridge`, an independent route through `eigen()`): fixture A's 299 marginals to 7·10⁻¹⁴ under seven options, and fixture B's rungs' log Z_b to 10⁻⁸.
- *Against the exact posterior:* the chain on fixtures A and C; each rung of a reference ladder; the self-sized ladder's evidence and SMC's normalizer on B and C, within 4–5 SE plus small absolute floors (0.003–0.01).
- *The rest:*
  - the draws against a joint-Gibbs reference on a fixed structure, within 5 SE;
  - the pointwise contents against dense references and leave-one-out refits;
  - resume, bit-identical;
  - end-to-end fits under the ridge, running every view.

### BASS: the sweep

**The rule** (D04). The coefficients a and σ² are integrated out of every block up to their own draw, and drawn immediately before the first block that conditions on them. Every block then samples a conditional of the same joint posterior, and no block reads a statistic from a stale draw. R's τ update is the pattern this rules out: it reads ‖Rβ‖² (R the Gram's Cholesky factor) from a coefficient draw made before σ² was redrawn (the derivations review, §2).

**The chain** (`BassChain`). The state, `BassState`, is the structure (its terms in the design's column order, with their fit) and u = ln τ. Each sweep is a `BlockScan` over the `ReversibleJumpChain`:
1. one structure move on D01's schedule (`BassStructureMoves`): a birth, death or change on p(y | structure, τ), with a and σ² integrated out;
2. in the g mode, τ's slice step on u (`BassTauBlock`) at the chain's own temperature: width 1, at most 32 stepping-out steps, as BAPS's component smoothing does, each evaluation O(1) from the structure's cached fit. The density is already on u, so no Jacobian is added.

No block conditions on a, σ² or λ. Against the derivations review's sweep, λ's conjugate block leaves the chain, since D02 integrates λ out, and returns as a draw at retention.
- *τ under tempering* (the stage-4 plan's choice 2). τ's step is a chain block, so it can run unchanged on a tempering replica between its moves and the swaps, through the per-replica blocks `ParallelTempering` gained for it (`5ce7453`); S8 seats it there. The alternative, τ as a Metropolis move inside the palette, would have changed D04's block.
- *The later slots* (F-X2). `BassSweepSlots` writes this order once, for the chain, the tempering replicas and the SMC kernel, in D04's general form. The slots after τ belong to an observation model with latent variables, and Gaussian BASS leaves them empty ([the extension seams](#bass-the-extension-seams)).

**The design follows the chain.** A proposal leaves its candidate pending in the chain's design. When the chain next shows a state holding that candidate's structure, the chain accepted it, and the design commits it; a state's structure identity is its terms' array. A τ update shares the structure and needs nothing. A restored chain takes its design back from its snapshot (S7), and a replica handed a state by a swap adopts its giver's (S8): a whole load differs from the incrementally updated design in the last bits and, near D06's angle floor, could refuse what the giver admitted. Only a structure no design holds, a reference draw, is loaded whole, by the admission that drew it.
- *A known approximation* (S6's review): an accepted state keeps its candidate's statistics from the update formulas, while the design refits explicitly, about 10⁻¹² apart.

**Composition at retention.** A retained draw is composed from the chain's state. Each part is an exact draw given (structure, τ), so no argument about the scan order is needed:
- σ² from its collapsed conditional, IG(g₁ + ν/2, g₂ + S/2) with ν = N (N − 1 under the flat intercept), which under the default is IG(c + N/2, c·s² + S/2);
- the coefficients given σ², from the block's conditional on the accepted design (`FromCentredLeastSquares`);
- λ | M ~ Gamma(h₁ + M, rate h₂ + 1).

Composition draws from its own stream, beside the moves' and τ's, so the trajectory does not depend on how often draws are retained.

**The draw** (`BassDraw`) is every draw's generative state: the structure, τ, σ², the coefficients and λ.
- *The structure* (`BassStructure`) is sorted ids into a `BassTermTable`, an unordered collection, so a retained draw never aliases mutable state.
- *The coefficients* are uncentred, in R's basis: the intercept, then one per id, so f = a₀ + Σⱼ aⱼBⱼ. Each column's mean is recomputed exactly as the design formed it (a case column's mean, or āh̄ for a direct term) to move the centred intercept to a₀, at the cost of one column evaluation per term per retained draw. Prediction (stage 5) and the pointwise density (S9) then read the basis directly, with no centring to undo.
- *Starts.* `BassStructureMoves` refuses a start outside F, the target's exact-fit refusal included. Starting each chain from an exact prior draw on F, rather than at R's M = 0 (the stage-4 plan's choice 3), is the run's (S7).

**Checked.**
- *Under a constant likelihood* (`BassSweepTests`): the laws under the move schedules above, on both schedules, with R's thirds asserted corrected.
- *For a fixed structure:* E[τ] against quadrature of its marginal; and the composed τ, σ², σ⁴, coefficients and a₁² against a joint-Gibbs reference that draws a, σ² and τ in turn on the explicit design [1, Z_c] for 400,000 sweeps, all within five errors.
- *The design follows the chain* through 4,000 sweeps and a swapped-in state, its statistics matching a fresh load to 10⁻¹⁰. A retained draw's uncentred coefficients match the same draw in the centred basis to 10⁻¹². On a direct fixture (n_x = 25, n_r = 7), a case-only, an output-only and a mixed term reproduce the centred fit at all 175 observations to 10⁻¹².
- *Against the exact posterior* (`BassChainParityTests`; four chains of 40,000 sweeps from the intercept alone):
  - fixture A at M_max = 3: the count law for M = 1–3 at 0.165, 0.762 and 0.073 against 0.176, 0.751 and 0.073 (± 0.009); E[τ] 0.01537 against 0.01548 (± 0.00007); all 48 terms' inclusions within a third of their allowance;
  - fixture C at M_max = 2: the count law, E[τ] 0.04996 against 0.04999, and all 90 inclusions;
  - fixed g, the flat intercept and the Jeffreys σ² prior on fixture A at M_max = 2, each against its own exact posterior.

### BASS: the local-move options (F-M2)

F-M2's Phase 1 (decided 2026-10-29; built by a chip, [its report](../../../issues/regression/bass/local-move-options-brief.md#reports)) adds kernels beside D01's, each an option, off by default and then bit-identical. Phase 0 and the mode-exit pricing ([the design draft](../../../issues/regression/bass/local-move-mixing-brief.md)) found why the default mixes poorly: the chains hold local modes 12–42 nats apart whose terms are each load-bearing, and births are bounded by the Occam cost, about −½ ln N per term. Which option becomes a default is decided by F-M2's pre-registered pass (Phase 2):
- eight arms: the baseline, the local steps, the fixed-M family, and both; the prior-cancelling schedule on the baseline and on both; tempered with the split rule, on the baseline and on both;
- four fixtures (S10's three, and Friedman at n = 2,000), and five master seeds × four chains;
- the rule, fixed before any run: an arm's change becomes the default when its geometric-mean ratio of minimum ESS per second, over the fixtures, is at least 1.5 with no fixture below 0.8.

**The pass's verdict** (2026-10-29; [its report](../../../issues/regression/bass/local-move-pass-brief.md#reports), rescored independently): no arm clears, so every option stays an option and the defaults stay D01's.
- *The scores:* the fixed-M family 1.031; both levers 0.998; prior-cancelling with both 0.984; the local steps 0.929; prior-cancelling 0.864, with the mixed fit at 0.793, under the floor; tempered 0.078 and, with both levers, 0.067. The verdict holds with any one fixture left out, and with CPU time in place of wall time.
- *What the metric read.* At S10's lengths (5,000 + 2,000 sweeps, four chains) the chains hold their own bands of M, as S10 and Phase 0 found.
  - So every untempered fit's minimum bulk ESS sits at the estimator's floor for split chains that disagree, about their count over two: a median of 4.96 over 120 fits, 101 of them below 6. M's R̂ has a median of 1.96–2.45 per arm.
  - The untempered ratios therefore mostly compare times. No untempered option lifts the chains off that floor at these lengths, and the pass cannot see what any of them does within a mode.
- *Tempering* lifts the minimum ESS 1.8–3.8-fold, holds M's R̂ to at most 1.55 and covers σ in all 40 of its fits. Untempered, Friedman's σ is covered in 0–2 of 5 fits per arm at n = 500. It costs 24–43 times the time, completes no round trip, and loses 13–15-fold in ESS per second.
- *Prior-cancelling* births accept 3.5–6.2 times as often as the default's and its deaths less often. It misses σ on both Friedman fixtures in all 20 of its fits there.

**Local steps** (`BassSampling.LocalSteps`, `LocalSteps(Probability = 0.5, WindowFraction = 0.02)`).
- *The kernel.* Inside a change, with probability ρ, a hinge's knot steps j′ = j + δ on its variable's distinct values, with its sign kept. |δ| is uniform on {1, …, w}, with w = max(1, ⌈ω·n_v⌉) per variable (`BassResolved.StepWindows`). A categorical factor instead toggles one level, uniform among the n_v. Otherwise, with probability 1 − ρ, the factor is redrawn from the prior, as D01's change does.
- *The ratio.* The step is symmetric and the knot prior uniform, so the acceptance is Δℓ. A toggle's adds the subset measure's ln π_v(D′) − ln π_v(D) through the prior: under R's measure ln C(n_v, |D|) − ln C(n_v, |D′|), under the uniform one 0.
- *Null moves:* a step off the knots, and a toggle to an empty or full subset.
- *Trade-off.* The coefficients are integrated out, so a nearby knot keeps most of the fit.
  - Steps accept 0.83 of admitted proposals on Friedman and the mixed fit, against the redraw's 0.06 and 0.15 with the options off, and 0.39–0.40 on the functional fit, where w is 1–2.
  - They refine within a mode. The pricing's drift endpoints are 14–53 distinct-value steps away on the functional fit, but 137–1,237 on Friedman, a long random walk.

**The fixed-M family** (`BassSampling.FactorMoves`, `FactorMoves(Change = 1, PruneOrGrow = 1, VariableChange = 1)`): a factor's prune and grow, and a factor's variable change. These moves are beyond R.
- *Selection.* A term uniformly, then a part uniformly among the fit's P parts.
  - A prune takes one of the part's K factors uniformly.
  - A grow takes a variable uniformly among the part's unused ones, with its sign and knot, or its subset, from the prior.
  - A variable change takes one of the K factors, then v′ uniformly among the unused variables, with v′'s prior.
- *Null moves and refusals.* A prune from an empty part, a grow into a full one (at its cap, or with every variable used), and a variable change with no unused variable are null moves. A prune that would empty the whole term is refused, since the empty term is the intercept, outside F.
- *The ratio is Δℓ alone.* The prior's ratio cancels the proposals'.
  - For a prune at K ≥ 2, by C(p, K)/C(p, K − 1) = (p − K + 1)/K.
  - For a part's K = 1 → 0, allowed while the term keeps another part, because each part's order is uniform on {0, …, cap} and C(p, 1) = p.
  - For the variable change, the prior's n_v/n_v′ cancels the proposals' n_v′/n_v.
- *The schedule.* The family shares the change's slot, by the weights' shares (`BassResolved.FactorShares`, thirds by default, confirmed by the user on 2026-10-29), so the births' and deaths' weights stay the schedule's own. The prune and the grow reverse each other at equal weights, and the variable change reverses itself.
  - The alternative, a slot of its own in the schedule, would lower the births' share and confound the family's reading with F-M1's.
- *Trade-off.* The pricing found near-free two-step exits where a worse mode differs from a better one by a factor's variable or a redundant factor: the mixed fit's chain 2, and the functional fit's chain 3.
  - It found none on Friedman's worse modes, whose cheapest exit costs 9.8–19.2 nats.
  - In the smoke reading the family accepts at most 0.03 of admitted proposals. On Friedman 69–71% of its prunes are refused, since every main effect's prune would leave the intercept.

**Null moves** (with the options on only) propose the state itself with a log ratio of −∞. The chain attempts and rejects them, the state and the design untouched, and draws its uniform as for any proposal. This keeps each kernel's selection symmetric at its edges.

**Not built**, by the decisions:
- several structure moves per sweep untempered, which is close to `Thin` = k;
- data-informed births, which lose to the Occam cost and not to inert inputs (×1.07 pooled);
- the pair death, since no pair of terms was jointly redundant in any priced state.

F-X4's proposal-kernel contract is built with KHAOS's Bernoulli-inclusion proposal instead (the extensions brief): `IBassProposalKernel`, whose hinge member is D01's proposal, its draws unchanged, and whose Legendre member is KHAOS's (2026-11-07; [Polynomial chaos](../PolynomialChaos.md#the-proposal-f-x4-and-khaoss-bernoulli-inclusion-e4)).

**The tallies** (`BassStructureMoves.Kinds`; `BassDiagnostics.MoveKinds`, a `BassMoveTally`). For each `BassMoveKind` (birth, death, redraw, local step, toggle, prune, grow, variable change), the attempts, acceptances, null moves and refusals, each chain's adding up to its chain's own counts.

**Checked** (`BassKernelEnumeration`, `BassLocalStepTests`, `BassFactorMoveTests`, `BassMoveOptionsTests`, `BassBitIdentityTests`).
- *Detailed balance by enumeration,* at fixed τ = 1/N.
  - Every proposal is enumerated from the moves' stated selection rules. The production code contributes only the schedule's weights, the target and each move's ratio.
  - π(S)K(S, S′) = π(S′)K(S′, S), and πK = π, hold per kernel and over the whole structure kernel with births and deaths. The worst gap is 5.7·10⁻¹⁴ in logs.
  - The fixtures: the 2,263-state continuous one, a categorical one (2,918 states, under both subset measures) and a functional one (593 states). Both schedules, at b = 1 and 0.3.
  - Negative controls fail it: a step ratio 0.05 off, a toggle without the measure's ratio, and a part-emptying prune with an extra ln 2.
- *The draws against the statement:* every outcome within five standard errors of its stated probability.
- *The constant-likelihood laws* on both schedules, with each option and with all on.
- *Visit frequencies:* the real chain against the exact posterior. With all options on, M = 1–3 reads 0.8145, 0.1696 and 0.0159, against 0.8153, 0.1688 and 0.0160.
- *Bit-identity:* eight fits' hashes, recorded before any production edit, hold with the options off: S10's three fixtures on both schedules, Friedman tempered self-sizing, and BARS tempered.
- *The design follows the chain* through every new move, on every replica.

### BASS: the birth mixture and the history birth (E4)

Built 2026-11-08 by the proposal-kernels brief's slice 1 ([the brief](../../../issues/regression/bass/proposal-kernels-brief.md)), as an engine kernel for hinge and Legendre terms alike. It changes no default: a fit opts in.

**The mixture** (`BassSampling.Births = new BirthMixture { History = new HistoryBirth(...) }`; `BassBirthMixture`).

$$q(t\mid S)=(1-\omega_H)\,q_{\text{base}}(t\mid S)+\omega_H\,q_H(t\mid S).$$

- q_base is the fit's own birth, D01's or KHAOS's, unchanged. Its share 1 − ω_H stays positive, so every structure it reaches stays reachable, and the mixture is bounded below by it.
- One uniform picks the component, then the component draws. ln q is the log-sum-exp of the components' exact ln q, and a death's reverse is the same sum under the counts without the dying term.
- Deaths stay uniform. Changes, KHAOS's mutations and the fixed-M family stay the fit's own.
- The mixture lives with each chain's moves (`BassStructureMoves`), over the model's kernel (F-X4). Without it the moves call the model's kernel as before, and every pinned hash holds.

**The history birth** (Nott and Kohn 2005, *Biometrika* 92).
- ĥ_u is the share of the recorded states in which some term uses variable u, per part (`BassInclusionHistory`). It is Nott and Kohn's running mean of the inclusion indicators, read over BASS's terms.
- The weights are max(ĥ_u, δ), with δ = 0.01 by default: their truncation, whose upper bound is moot for a weight.
- *Hinge births.* The weights replace z_u ∝ w₂ + (terms using u) at every order, K = 1 included, where D01 draws its variable uniformly. Wallenius' law gives the set's probability.
- *Legendre births.* The weights take the place of (aⱼ^α + ε) in `make_weights`, whose passes bring η to sum q₀. The inclusion law and the degree and composition laws are KHAOS's.
- *The order law* is the base kernel's: the history touches variables only, as E4 states it.
- Before the first refresh there are no weights, and q_H = q_base.

**D10's policy, as built.**
- *Chains.* The run's burn-in records the state after every sweep and refreshes the weights every `Refresh` sweeps, from the past only. At burn-in's end the history refreshes once more, from every record, and freezes.
  - Production then runs a fixed kernel. Nott and Kohn's condition (3), each step's balance, and (4), a uniform minorization through the base share, hold at every step. Their (5), diminishing change, is not needed.
- *A ladder.* Every rung reads its b = 1 rung's history, which records the b = 1 state and freezes with the ladder. A frozen kernel leaves each rung's tempered law invariant.
  - The brief planned a history per rung. One per chain was decided in the slice: a self-sizing ladder rebuilds rungs as it resizes, and a shared history survives that.
- *SMC.* `BassSmcKernel.Prepare` sets the weights from each stage's weighted population, and they stay fixed through the stage's mutations. That is the adaptive SMC of Beskos, Jasra, Kantas and Thiery (2016), as the adaptive count of sweeps is.
- *The flat histogram.* Each chain learns through the learning phase and freezes when production starts. Its pilot ladder does not learn and runs the base birth.
- *Checkpoints* carry the counts, the records, ĥ and the phase. Every resume is bit for bit: a chain, a ladder, SMC, and the flat histogram across learning and production.
- *The readout.* Each chain's frozen ĥ and its record count (`BassChainAccount.Inclusion`, `BassInclusionReading`): Nott and Kohn's estimate of each variable's marginal inclusion probability.

**Checks** (`HistoryBirthTests`).
- *The kernel's own law,* with a deliberately skewed ĥ that leaves the main part's first variable at its floor. On three continuous inputs, the categorical and functional fixtures and KHAOS's fixture:
  - Σ_t q(t | S) plus the empty draw is 1 to 10⁻¹³, for the mixture and for the history component alone;
  - each ln q against the statement written out (the base law with each part's variable law replaced by Wallenius' over the weights, or `make_weights` over them) to 9·10⁻¹⁶;
  - 200,000 draws against q by χ², z from −1.39 to 0.29.
- *Detailed balance, exactly.* The sweep's flows against the enumerated law, each birth's forward probability from the statement and the ratios from the production mixture: a worst gap of 5.7·10⁻¹⁴ in logs over seven cases, a learned ĥ among them. A mixture whose ratio reads another share fails it.
- *The chain's law* after a learned freeze, 400,000 sweeps against the enumeration: worst |z| 2.07 on a hinge fit and 2.21 on a Legendre fit.
- *Geweke,* with proper priors and a skewed ĥ, 200,000 marginal draws against 300,000 successive steps: worst |z| 2.32 (hinge, three inputs at order 2) and 1.45 (Legendre).
- *The reach.*
  - A ladder's rungs share the b = 1 history through a self-sizing ladder's resizes, and each rung's count law and E[ln τ] meet the enumeration's tempered law, hinge and Legendre.
  - SMC's stage weights are the population's weighted share, and its ln p(y) meets the enumeration's: −72.2912 ± 0.0137 against −72.2910 under the raw ridge.
  - The flat histogram matches the enumeration.
  - Resume is bit for bit in every sampler.

**The measurement** (pre-registered in the brief; `ProposalKernelMeasurementTests`, the table [`proposal-kernels-history.csv`](../../../issues/regression/bass/proposal-kernels-history.csv)).
- *The design.* Four master seeds of four chains, 5,000 burn-in and 5,000 kept sweeps, each fit on one thread, with ω_H = ¼, δ 0.01 and a refresh every 100 sweeps. The reading is bulk ESS per CPU second, its mean over the seeds, the history birth's over the fit's own. Tiered compilation was off, so no fit carries the process's compilation.

| Fixture | ESS per CPU second, history over base | Beyond the seeds' spread | Births accepted | CPU per fit |
| --- | --- | --- | --- | --- |
| Friedman, p 10, the g-form | M 0.73, ln π 0.91, σ² 0.76, fitted values 0.50–0.96 | M lower (z −4.2), f₂ lower (z −2.2) | 1.86% → 2.03% | +7% |
| Friedman, p 50 (45 inert) | M 1.19, ln π 0.59, σ² 2.84, fitted values 0.92–1.50 | f₀ higher (z 2.1) | 2.08% → 1.49% | +9% |
| Ishigami, the raw ridge | M 1.38, ln π 1.41, σ² 1.27, fitted values 1.48–1.51 | none (\|z\| ≤ 0.95) | 1.25% → 1.49% | +7% |
| Ishigami, the modified g-prior | M 0.88, ln π 0.80, σ² 1.00, fitted values 0.70–2.20 | none | 13.8% → 13.9% | +4% |

- *Reading.*
  - No fixture gains across D10's groups.
  - On Friedman the structural group fails D10's verdict under both kernels: R̂ of M is 1.3–3.0. The history birth does not move chains between Friedman's representations, as the local-move mixing brief's reframe expects: there the problem is moving between representations, not the births.
  - On p 10 it slows M's mixing, and on p 50 it helps one fitted value, one comparison of 27.
  - τ's ESS under the g-form falls 13× (p 10) and 25× (p 50), beside a base spread wide enough that neither is beyond it (z −1.7, −1.0).
  - On Ishigami under the raw ridge every ratio is 1.3–1.5, within the seeds' spread.
- *The decision it changes:* no default candidate. The history birth stays opt-in.
- *The cost,* on the cycle clock, indicative. The kernel's work is 1.6–1.8× per birth (1.1 to 1.8 µs at the cost report's base, 6.3 to 10.7 µs on Ishigami) and 2.0–3.2× per death's ratio. The design's work is unchanged, and a whole fit takes 4–9% more CPU.

**Not built** (each recorded with its decision in the brief):
- Nott and Kohn's best linear predictor, the inclusion of u predicted from the current γ₋ᵤ by the history's covariance;
- adaptation that continues past burn-in, under their condition (5);
- history weights in the swaps and variable changes.

### BASS: the residual birth (E4)

Built 2026-11-09 by the proposal-kernels brief's slice 2 ([the brief](../../../issues/regression/bass/proposal-kernels-brief.md)), a component of [the birth mixture](#bass-the-birth-mixture-and-the-history-birth-e4) for hinge and Legendre terms alike. Since 2026-11-09 it is the default for KHAOS's raw-ridge fits of a scalar response, at ω_R = ¼: the user's decision on the five-fixture measurement below. A config that declares no mixture takes it there (`BassResolved.Births`, logged as an auto-pick), and `new BirthMixture()` turns it off. Elsewhere a fit opts in with `BirthMixture { Residual = new ResidualBirth(share) }`, beside the history birth or alone.

$$q(t\mid S)=(1-\omega_H-\omega_R)\,q_{\text{base}}+\omega_H\,q_H+\omega_R\,q_R .$$

**The component** (`BassResidualBirth`).
- *The residual.* r(S) is the raw response's least-squares residual on [1, B_S], the structure's raw columns in a canonical order. It is a function of the terms alone: a birth's forward probability and the reverse of the death that undoes it read the same number.
  - A tiny ridge, 10⁻¹² of the Gram's largest diagonal, keeps a collinear structure's solve defined.
  - An exact fit's residual is rounding: at or below 10⁻²⁴ of the response's own spread it weighs every candidate alike.
- *The parent* p is drawn uniformly over P(S): the intercept, and the structure's terms that can take one more continuous factor. That means a hinge within the continuous cap, or a Legendre factor within q_max and d_max.
- *The new factor* f is drawn ∝ exp(½·b·Δℓ(f)), with Δℓ(f) = −(n/2)·ln(1 − ΔQ(f)/‖r‖²), the single column's profile-likelihood gain.
  - ΔQ(f) = (c̃ᵀr)²/‖c̃‖², c̃ the centred column B_p ∘ φ_f. A hinge's scale cancels, so its shape alone, (s(x − t))₊^q, is scored.
  - The exponent ½ is the locally balanced √ weighting (Zanella 2020). b is the move's inverse temperature: a ladder sets each rung's every scan, and SMC each stage's. At b near 0 the component is near uniform over candidates, as the flat tempered law is.
  - *Hinge candidates:* every variable p leaves, both signs, every knot of the prior's set. They are scored at once per variable by cumulative sums of B·r·xʲ, B·xʲ and B²·xʲ along its sort order, any degree by binomial expansion.
  - *Legendre candidates:* every variable p leaves, every degree within d_max − d_p.
  - Categorical and functional factors are the base birth's to propose.
- *The forward probability, exactly:* q_R(t | S) = Σ_{p ∈ sub(t) ∩ P(S)} (1/|P(S)|)·w_p(t ∖ p)/W_p. Here sub(t) are the terms that t less one continuous factor leaves, the intercept for a one-factor term.
- *The death's reverse* reads the reduced structure's residual through F-X4's birth context (`BassBirthContext`: the structure a birth starts from and b), which the moves build only for the mixture. The hinge and Legendre kernels take the context's default, the counts alone, so their draws are unchanged.
- *Caches,* all pure, so they change no draw:
  - per instance, the residual by structure and each ln W_p and ln q_R by structure, parent or term, and b, with the last full table kept for draws;
  - per model and shared, the raw columns.

**The reach.**

| Path | Under the residual birth | Gate or reason |
|---|---|---|
| A chain | Runs | The chain's law against the enumeration; Geweke |
| The ladder | Runs, each rung at its own b | Each rung's law against the enumeration, hinge and Legendre |
| SMC | Runs, each stage at its b | ln p(y) against the enumeration under the raw ridge and the modified g-prior, flags clear |
| The flat histogram | Runs at b = 1 | Against the enumeration |
| Resume | Bit for bit | A chain, a ladder, both components together, and SMC across stages |
| GBASS's noise, the random-intercept block | Runs | The score reads the raw response's residual. That is an approximate score under weights or whitening, with exact probabilities |
| A direct functional response; several outputs (E34) | Refused at `BassModel.Create` | The residual is a matrix, and no score here reads one |
| Data with no continuous predictor | Refused at `BassModel.Create` | There is no continuous factor to add |

**Checks** (`ResidualBirthTests`).
- *Primitives.* Every candidate's score from the cumulative-sum pass against its dense column, over twelve structures per fixture at b 1 and 0.3: a worst relative gap of 1.9·10⁻¹¹. The fixtures cover ties at knots, hinge degrees 1–3, inputs shifted by 10⁶, columns that vanish on the cases, Legendre and categorical fits.
  - The residual against modified Gram–Schmidt to 2.6·10⁻¹² of ‖r‖, the same bit for bit in any term order.
  - An exact fit weighs every candidate alike.
- *The kernel's own law,* at b 1 and 0.4. Σ_t q_R(t | S) = 1 and the mixture's Σ_t q with its empty draw = 1, to 10⁻¹⁴. ln q_R against the statement written out (dense columns, the statement's parents) to 1.1·10⁻¹¹. 100,000 draws against q_R by χ², z from −1.83 to 0.93.
- *Detailed balance, exactly,* against the enumerated law. The births' forward probabilities come from the statement and the ratios from the production mixture. The worst gap is 1.3·10⁻¹¹ in logs over eight cases: hinge at b 1 and 0.4, with the history birth beside it, categorical, three inputs at order 2, quadratic hinges, and Legendre with and without the history birth. A death that reads the full structure's residual fails it.
- *The chain's law,* 300,000 sweeps against the enumeration: worst |z| 1.80 (hinge) and 1.30 (Legendre).
- *Geweke,* with proper priors at ω_R = ½, the redrawn response moving the residual every step: worst |z| 0.56 (hinge) and 2.05 (Legendre).
- *The reach:* the ladder's rung laws, SMC's ln p(y) (−72.2983 ± 0.0126 against −72.2910 under the raw ridge, −69.5983 ± 0.0143 against −69.5869 under the modified g-prior), the flat histogram, and resume, as in the table.

**The measurements** (pre-registered in the brief; ω_R = ¼, the exponent ½, tiered compilation off).

*Mixing per second* ([`proposal-kernels-residual.csv`](../../../issues/regression/bass/proposal-kernels-residual.csv)): four seeds of four chains, 5,000 + 5,000 sweeps, the residual birth over the fit's own.

| Fixture | ESS per sweep, residual over base | ESS per CPU second | CPU per fit | Births accepted |
| --- | --- | --- | --- | --- |
| Friedman, p 10 | 0.47–1.57, none beyond the spread | 0.03–0.07 | 0.62 → 11.9 s | 1.86% → 1.74% |
| Friedman, p 50 | 0.08–2.25, none beyond the spread | 0.00–0.02 | 0.60 → 52.7 s | 2.08% → 1.36% |
| Ishigami, the raw ridge | 1.83–2.23 for M, σ² and the fitted values (z 1.9–2.4); ln π 1.34 | 1.86–2.36 (z 2.0–2.5); ln π 1.24 | 1.64 → 1.71 s | 1.25% → 1.38% |
| Ishigami, the modified g-prior | 0.77–2.49; f₁ 2.49 (z 2.2), f₂ 1.55 (z 2.9) | 0.63–2.07; M 0.63 (z −2.0) | 3.61 → 4.46 s | 13.8% → 14.5% |

- On Ishigami under the raw ridge the residual chains also agree: R̂ of M at most 1.015, against up to 1.37 for the base.
- On Friedman the residual birth leaves per-sweep mixing where it was, as the local-move mixing brief's reframe expects at a mode. Its cost is the scoring: O(n·p) per parent, with a logarithm and an exponential for every knot.
- *The cost row,* on the cycle clock, indicative. A warm birth's kernel work is 9.3 µs against 3.6 on Ishigami, and 0.97 ms against 0.83 µs at the cost report's hinge base (n 1,000, M 20). A warm death's ratio, once its probability is cached, is 3.1 against 2.1 µs and 5.4 against 0.3 µs.

*The SMC trap* ([`proposal-kernels-trap.csv`](../../../issues/regression/bass/proposal-kernels-trap.csv)): the raw-ridge Ishigami fixture, N 500 at 5 sweeps a stage, 16 replicates each, against slice 5's table ([its addendum](../../../issues/regression/bass/khaos-brief.md#slice-5-addendum-smcs-flags-chip-2026-11-08)).
- *The engine's schedule:* 0 of 16 reach the posterior. Each still jumps over the collapse from b below 0.03, where the residual weights are nearly flat.
- *The 1.25 cap:* 16 of 16 reach it, at M 10.09–10.18 and σ² 0.222–0.227 against the long chains' 10.15 and 0.2245. The cap alone reached 11 of 16. Births accept 8.4%, against about 2%.
- *The evidence still spreads:* ln Ẑ −287.0 ± 6.4 over 21 nats, against per-replicate errors of 0.4–0.8, so the disagreement flag reads it.

*MCMC's climb from empty* on the same fixture, at b = 1, eight seeds each: to the long chain's band, the 5% quantile of its ln π, in a median 356 sweeps (139–615) against 4,384 (2,241–7,557). That is 12 times fewer.

**The decisions they change,** each the user's:
- *On KHAOS under the raw ridge,* the residual birth was first reported as a candidate for that class: ESS per CPU second about 2× on five of six functionals at the same cost, and the climb from empty 12× shorter.
  - Then measured on five fixtures by a rule fixed before the run (the proposal-kernels brief's addendum, 2026-11-09): Ishigami at n 120 and 500, Ishigami with seven inert inputs, Sobol's g-function and the borehole, eight seeds each.
  - Four gains (median ratios 1.49–8.65), no loss, and the borehole neutral: its log posterior, which mixes poorly under both births, read 0.50 while every other functional gained. The climbs were 8–39× shorter, at 26–100% more CPU per fit.
  - By that rule it went to the user as a recommended default for the class, and the user made it the default on 2026-11-09 (`0738183`). A config that declares no mixture takes it on KHAOS's raw-ridge fits of a scalar response; an empty mixture turns it off; a declared one runs as declared.
  - The fit's own birth stays the default for KHAOS's other priors, the modified g-prior and the normalized ridge, and for hinge fits. The fixture behind the earlier SMC measurements (`BassSmcMutationQualificationTests.Khaos`) now declares the fit's own birth, so they reproduce as they ran.
- *For SMC's trap,* the pair of the 1.25 cap and the residual birth reaches the posterior in every replicate, and the evidence stays dispersed. The remedy table's third row, a reference nearer the posterior, is what the evidence would need.
- *On hinge fits* it stays opt-in: no per-sweep gain on Friedman, at 19–88× the CPU.

**Not built** (each recorded with its decision in the brief):
- a cheaper hinge score (a knot grid, or a residual updated rather than refitted), which matters only where a hinge problem gains per sweep;
- residual-informed knots in the change moves;
- a sharpness other than ½;
- a functional or several-output residual.

### BASS: the run

**The fit** (`Bass.Fit(config, data, run)` → `BassFit`; stage 4's S7). One `BassChain` per chain on the ensemble shell (`ChainEnsemble`), or with tempering one `BassLadder`, the chains in parallel within `BassRunOptions.MaxDegreeOfParallelism`. A ladder's rungs step in order unless `BassRunOptions.RungWorkers` asks for more: −1 takes what the chains leave of the bound (`BassRunOptions.Split`: min(chains, P) chains at once, ⌊P / that⌋ rung workers each, P the processor count when unbounded), and n at most n. The draws are the same bit for bit at every setting ([Samplers](../Samplers.md#tempering), rungs at once; [the report](../../../issues/regression/bass/rung-parallelism-report.md)). Burn-in counts sweeps, and each retained draw follows `Thin` more. The fit is immutable: the config and its resolution, the input transform, the merged term table, each chain's `DrawStore<BassDraw>`, the diagnostics, the evidence, the residual check and, when asked for, the checkpoint.
- *Starts are exact prior draws on F* (`BassPriorSampler`; the stage-4 plan's choice 3), not R's M = 0. M is drawn from its count law with its orderings, p(M) ∝ Γ(h₁ + M)/(M!(h₂ + 1)^M) truncated at M_max (λ integrated out, so no Poisson sampler is needed), then M terms from D01's measure and τ from its prior; a collection is kept only in F. Overdispersed starts are what make R̂ informative. The sampler refuses after 10⁶ attempts, which happens when F holds little of the prior's mass (a count prior pushed toward a cap on few points, for one).
- *Streams.* Each chain's seed fans out into four, by the seed tree: the start's draws, the moves, τ's block and a ladder's own. The shell's readout seed drives composition. So the trajectory depends neither on how often draws are retained nor on how many attempts the start took, and a fit is identical under every parallelism bound.
- *Terms.* Each chain interns its draws' terms in its own table; the fit merges them in chain order and remaps each draw's ids, so ids never depend on the schedule.

**The diagnostics** (`BassDiagnostics`), over the sampling phase only:
- rank-normalized R̂ and bulk and tail ESS for M, τ in the g mode, each retained draw's log posterior (the collapsed log π(structure, u | y) at b = 1, the chain's own density; D10), σ², and each draw's fitted value at the monitored points (`BassSampling.MonitoredPoints`: declared cases, with an output location for a direct response, or by default up to three cases at the quartiles of the data's order, at the middle location; decided 2026-10-28);
- acceptance per move type, every proposal by its kind (`MoveKinds`: attempts, acceptances, null moves and refusals; F-M2's options), refusals by reason, empty births, the support outcome of each candidate term by its orders, the mass at the cap, and each chain's own account;
- with a positive `RHatTarget`, the run stops early on a batch boundary once its gating groups meet it (`StoppedEarly`, `SamplesUsed`);
- *the verdict* (D10, built 2026-10-29 by N14; `BassDiagnostics.Verdict`), read on every run, against 1.01 when no target is set (see [Samplers](../Samplers.md)):
  - *Three groups:* structural (M, τ), level (the log posterior) and predictive (σ², the fitted values), each reported with its failing quantities.
  - *The scope:* `BassSampling.StoppingScope` names the groups a target gates, `All` by default or `Predictive` alone. A predictive-only stop can end a run whose chains sit at different posterior levels, so it is a declared choice.
  - *Beside it, the mode-separation readings* (`ModeSeparation`): the gap between the chains' mean log posteriors, and the largest distance between two chains' mean fitted values. Nothing stops on them.
  - *On S10's fits* every group fails. The log-posterior gap is 42.8 nats on Friedman, 21.8 on the mixed fit and 30.9 on the functional fit; tempered, Friedman's falls to 27.0.

**Checkpoints** (`Bass.Resume(checkpoint, config, data, run)`, which takes the fit's own `EnsembleCheckpoint`). Each chain saves its state, its streams, its counters, tallies, term table, draws and start, and its design as it is, not the structure to reload. A strict reload differs from the incrementally updated design in the last bits after deaths and changes (about 10⁻¹⁶ in ‖r‖²), which would break bit-identity, so the snapshot carries an exact copy of the factor (`BorderedCholesky.Copy`). A resumed run is bit-identical to the uninterrupted one; a config, seed or data change is refused by the fingerprint.

**Checked** (chip A's [report](../../../issues/regression/bass/stage4-run-brief.md#reports)).
- The prior sampler against the enumerated F: the count law, the order law, χ² over every state (357 on 339), τ, and the acceptance rate against F's exact share (0.6540 against 0.6550 ± 0.0010); a mixed and a direct functional term space, and Jeffreys on exact-fit data.
- Each chain rebuilt by hand from the seed tree reproduces the fit's draws exactly; the fit is the same under parallelism bounds 1, 2, 3 and none; thirteen options run end to end with every draw in F.
- Against the oracle's exact posterior (`BassFitParityTests`, four chains of 40,000 draws): on fixture A at M_max = 3 the count law 0.1717, 0.7498, 0.0785 against 0.1759, 0.7508, 0.0733, E[τ] within its error and all 48 inclusions; on fixture C all 90 inclusions.
- Resume bit for bit under eight settings, twice over, and the refusals.

### BASS: tempering and the evidence

**The ladder** (`TemperingSpec` in `BassSampling`, null by default: untempered, D05; S8). Each chain becomes a `ParallelTempering` ladder over replicas of D04's sweep, the b = 1 replica the sampler, whose draws are retained.
- *Path and form.* The likelihood path, prior · p_b, by default, with D05's integrated form p_b (b tempers the likelihood inside the a, σ² integral) or the linear form b·ln p₁; the powered path selectable. Swaps accept on ℓ_b alone, since the prior cancels. DEO's even-odd swaps follow every sweep, from the first scan. R's `start.temper`, which delays them by 1,000 iterations or half the burn-in, whichever is smaller, is not carried (the user's ruling, 2026-11-01): it changes no target, and the self-sizing ladder reads its swap rejections during burn-in, which a late start would starve. `BassResolved` refuses a rung where the form is undefined (`BassTarget.Defines`): b = 0 under Jeffreys or the flat intercept, and under the flat intercept b ≤ (1 − 2c)/N.
- *Replicas.* Each runs its own structure moves and its own τ block at its rung. A state a swap hands over brings its design: the receiver adopts the giver's, its factor and fit copied and its columns shared (a design never writes into a column once evaluated), so a structure one replica admitted is neither refactored, re-admitted nor re-evaluated by another; a re-admission near D06's angle floor could refuse it. A reference draw, which no design holds, is loaded whole by the admission that drew it.
- *A declared ladder* (`InverseTemperatures`, R's `temp.ladder`, 2026-11-02): the rungs themselves, in place of `Levels` and `HotEnd`, strictly descending from exactly b = 1 and ending above 0, or at exactly 0 with the reference rung. A list declares its spacing, so it runs as given, as R's does (the user's ruling, 2026-11-02). `Adaptive` null, its default, reads false for a list and true for a built ladder (`TemperingSpec.Adapts`), and `Adaptive = true` re-spaces the interior from the list. Its count is fixed, so it takes no `Split`.
- *The reference rung* (opt-in): the hot replica sits at b = 0 and draws exactly from the prior on F every scan (`BassPriorSampler`), which makes the ladder's low end normalized.

**A self-sizing ladder** (decided 2026-10-28, chip E). A fixed count cannot serve: for a power posterior the barrier grows as Λ ≈ √(d/2π)·ln N (a Bernstein–von Mises reading: ℓ fluctuates as χ²_d/(2b)). BASS's swaps read the collapsed ℓ_b, so its effective d is smaller, but the count must still follow the data. On the functional fixture (N = 175), S8's default 8 levels rejected 0.96–1.00 of swaps per pair. So `Levels` null, the default, sizes the ladder during burn-in by Syed et al.'s N ≈ 2Λ (`TemperingLadder.Size`):
- it grows at once to `RecommendedLevels(Λ̂)` when that exceeds the count, since a saturated ladder reads Λ̂ at most its pairs; it shrinks only when the reading asks for two fewer (a one-rung deadband), since early rounds over-read; and it stays within [2 (3 with the rung), `MaxLevels` = 64], a compute budget;
- the rungs are placed on the cumulative barrier (`Place`, `Equalize`'s monotone cubic inverted at the new count), each new rung forked from its nearest rung on the barrier with its own copy of the term array and its source's design;
- the count is settled one round before burn-in ends, so the last round, half the burn-in, equalizes it; the ladder is frozen after burn-in (D10);
- an explicit `Levels` fixes the count, and then the run is S8's bit for bit.

The rounds read the barrier high while the replicas relax: on the functional fixture Λ̂ read 5–7.5 in burn-in against 3.8 in sampling, so the ladder settled at 11–16 rungs where about 9 would do. That is the safe direction; it costs compute.
- *The alternatives.* A fixed count is cheaper and fails as N grows. A grow-only rule is simpler and keeps the early over-reads. SMC's adaptive temperatures (stage 8) size themselves by construction and are the other route to the evidence.

**The split rule** (`TemperingSpec.Split`, `SplitRule(Threshold = 0.9)`, an option for a self-sizing ladder; F-M2's lever 4b, decided 2026-10-29).
- *The fault it addresses.* Λ̂ = Σr counts a pair at most 1, however steep the band between its rungs, so the sizing under-reads a steep band. On Friedman the band lies at b ≈ 0.16–0.32. There Phase 0's profile at Δb = 0.01 read Σr 4.7–5.2 over [0.10, 0.40], against 4.0–4.3 on the self-sized ladders.
- *The rule.* A round in which a pair rejects past the threshold splits each steep pair at its barrier midpoint (`TemperingLadder.Split`, on `Place`'s monotone cubic), steepest first, within `MaxLevels`. Every current rung and its replica is kept. Once no pair is steep, the ladder sizes as before, on the finer readings. The run's last round adapts and never resizes, so a split ladder is equalized before it freezes.
- *The alternatives* (Phase 0's M5):
  - an uncensored increment, h(r) = √(2/π)·Φ⁻¹((1 + r)/2), over-read the finer ladder's Σr by 0.6–1.8;
  - thermodynamic length reads a rung's sd low when the rung holds one regime of a two-regime band.
  - A split measures the band rather than extrapolating an under-read barrier.
- *Read once* (one seed, the chip's smoke reading). On Friedman the rule settles at 19–21 rungs, against 20–23 without it, and raises the ideal round-trip rate 6.5-fold, to 0.0091 per scan against 0.0014. Neither ladder completes a round trip in 2,000 scans, so the local moves' efficiency binds. A synthetic band, a step of 6 in Λ, resolves in 4 rounds.

**The hot end** (`HotEnd` null, the default, auto-picked and logged): 1/N. At b = E[τ] = 1/N the shrinkage τ/(τ + b) is one half, so the tempered likelihood weighs one observation, as the unit-information prior does; below it the replica is prior-dominated. Under the flat intercept's integrated form it is 2/N, a full 1/N above the pole (ν = Nb − 1 = 1); under fixed g the smaller of the held τ and that. With the reference rung it is the adaptation's start and moves with it.

**The evidence** (`BassFit.Evidence`). Reported when asked for (`ReportEvidence`), with the reference rung, under proper priors; null under the flat intercept or Jeffreys σ², whose constants leave it undefined.
- *Stepping-stone* gives `LogEvidence`, each chain's over its own frozen ladder. *Path sampling* runs beside it on the linear form, whose slope is ℓ₁; on the integrated form it is left out, since its slope at b = 0 carries −(k/2)·E[1/τ], infinite under Zellner–Siow, and its discretization diagnostic misses E_b[∂²ℓ_b/∂b²].
- *The error is the larger of two* (`LadderEvidence.Pool`, shared with BARS): the within-chain error √Σs_c²/C and the between-chain error, the standard error of the chains' mean. Each chain adapts its own ladder, so the chains are independent replicates, and their disagreement catches rungs that do not communicate, which the within-chain errors cannot see. In a calibration of four chains spreading with sd 1 while each reports 0.1, the within-chain error covered 16% of the time and the reported one 85% (a t on 3 degrees of freedom gives 86%).
- *Flags, never a null:* the largest pair rejection, `Saturated` above 0.9 (a "look at the chains" mark: stepping-stone weights' relative variance grows as e^{σ²/2} − 1, about 220 at a rejection of 0.9), and `AtMaxLevels` when the budget binds.
- *`BassEvidencePriors`* records what the value depends on beyond the data: the coefficient prior as declared and resolved, the count prior with its cap, the term measure and the caps, the support rule and D06's gates (which define F, where the prior is normalized), the centre and its level, and σ²'s prior. Two evidences compare as a Bayes factor only where these agree.

**The readout** (`BassTemperingReadout`): each chain's pair rejections and acceptances, its frozen ladder and its count, the chains' mean barrier and ideal round-trip rate (Syed et al.'s Theorem 1), round trips and annealed restarts per scan, and the budget's flag.

**Checked** (chip D's and chip E's reports).
- Each rung's law against enumeration over both forms and paths, the reference rung, fixed g and the improper options; each b = 0 rung's τ against the prior's closed form.
- The swap's exact balance against the oracle's laws on fixture A (spread at most 3.2·10⁻¹⁴ over F); R's swap density, a defect fixture, misses its own rung laws by 1.13–1.88.
- The evidence against the oracle's exact ln Z: fixture B −28.554 ± 0.070 (linear) and −28.609 ± 0.086 (integrated) against −28.589, fixture C −30.882 ± 0.075 against −30.903, self-sized.
- The functional fixture, four chains of 4,000 burn-in and 4,000 draws: self-sized, ln Z 192.59 ± 0.52 (integrated) and 193.45 ± 0.41 (linear, path sampling 193.42 ± 0.53), against 193.40 ± 0.37 and 193.03 ± 0.30 at fixed 32 and 64 rungs; fixed 8 read 189.68 ± 1.51, the between-chain error showing what the within-chain 0.43 hid.
- Resume bit for bit across resizes and after burn-in; the clone invariant; a binding budget flagged and still reported.

### BASS: SMC

Stage 8's C4 (2026-10-29; [the chip's report](../../../issues/regression/bass/stage8-fitting-smc-brief.md#reports)) fits BASS by SMC on D05's sequence. It is an option beside MCMC: `BassSampling.Smc`, a `SmcSpec`, with null meaning MCMC, the default.

**The sequence** (D05):
- it starts from exact prior draws on F (`BassPriorSampler`) and runs to b = 1 on the likelihood path;
- the integrated form is the default (`Form`). The linear form is selectable, and it is the only one under D03's improper options, where the normalizer is no evidence;
- each particle's incremental weight is ℓ_{b′} − ℓ_b;
- stages are placed by conditional ESS (`TargetConditionalEss` 0.99), each step from b > 0 bounded at r·b when `StepCap` declares r ([the step cap](#bass-smc), below; null by default);
- the population is resampled multinomially when its ESS falls below half (`ResampleBelow` 0.5);
- the mutation is `MutationSteps` 5 of D04's sweep at the stage's b. `AdaptiveMutation` (null by default) instead sets each stage's count from the previous stage's fraction of particles whose structure changed, within 5 per stage taken; measured against the fixed count on 2026-11-08, it stays an option ([Samplers.md](Samplers.md), "SMC's sweeps per stage"). The kernel itself reads nothing from the population.

**A particle** is a `BassState` with its own moves and design. A fork copies the design as the ladder's hand-over does, the factor copied and the columns shared.

**Replicates** (`Replicates` 4, of `Particles` 1,000 each) are independent runs, the SMC analogue of chains, run in parallel.
- *Refused with `Smc`:* `Tempering` and an R̂ target.
- *Not read:* `Chains`, `BurnIn`, `Samples`, `Thin` and the stop's settings, which the resolution logs as not read.
- *The stage budget:* `MaxStages`, 1,000 per replicate. A replicate that spends it before b = 1 is flagged (`BassSmcReadout.ReachedPosterior`), and its posterior views refuse. `Bass.Resume` continues it, bit for bit.

**The draws keep their weights** (`BassFit.DrawWeights`, each replicate's population at weight 1/R in all), and every view reads them. `PointwiseLogDensity.Weights` puts them into WAIC and PSIS-LOO (the pointwise section below).
- *Shared ancestry* makes errors computed as though the draws were independent read optimistic. The readout reports each replicate's families beside its ESS, and the honest error is the spread between replicates.

**The readout** (`BassDiagnostics.Smc`).
- *Per replicate:* its stage reports (ESS, conditional ESS, resampling, families, acceptance, the normalizer's increments and genealogy variance), and its final ESS and families.
- *Across replicates* (`Spread`): the spread of each monitored quantity's posterior mean, in posterior-sd units.
- *No chain R̂:* it is unavailable under SMC, and the verdict reads undetermined in every group. Split R̂ means nothing for a population, which has no order to split.

**The evidence** (`BassEvidence.OfSmc`), under D03's proper defaults.
- Each replicate's normalizer, with its genealogy error (Lee and Whiteley 2018).
- Pooled by `LadderEvidence.Pool`'s rule, the larger of the within-run and between-run errors.
- Stepping-stone over the traced stages (`TraceStages`) equals the normalizer, and path sampling joins on the linear form.

**The flags** (2026-11-08, on KHAOS's trap: [Polynomial chaos](../PolynomialChaos.md#the-engines-reach-slice-5)). The fit raises two flags on itself from what the run already holds. Neither changes a draw, weight or normalizer: `BassBitIdentityTests.Smc_KeepsTheBaseline` pins an SMC fit's draws, weights and normalizers as they stood before the flags. The coordinating session chose flags over a remedy, the user unsure and the call reversible.
- *A jump over a collapse* (`BassSmcReplicate.JumpedOverACollapse`, counted in `BassSmcReadout.JumpedReplicates`). The last step to b = 1 left from b below ½ (`JumpBelow`), over a population of at most 5% of N distinct structures (`CollapsedShare`).
  - Such a population reweights alike, so the conditional ESS sees no spread and allows the jump, and ln Ẑ is read across a region the population never visited.
  - The run counts the structures once, before the step that reaches b = 1, and a checkpoint carries the count.
- *Disagreeing replicates* (`BassSmcReadout.Overdispersed`, mirrored in `BassEvidence.Overdispersed`). The replicates' ln Ẑ between-to-within ratio is above 3 (`BassEvidence.OverdispersionThreshold`).
  - The evidence's error is then their spread, the pooling rule's larger error. The genealogy errors are not to be read: each sees only the spread its own population holds.
  - The flag is null with one replicate, or where a genealogy error is NaN. The readout's `Flags` say so in words, rather than reading clear.
- *Why both:* replicates trapped alike agree with each other, so only the jump reads them. Replicates that settle apart need no jump to be wrong, and their disagreement reads them.
  - On KHAOS's trap fixture every group of four replicates, a default fit's, was flagged one way or the other.
  - On Friedman's band (below) the replicates' 6.5–8.4 nats of disagreement against genealogy errors of at most 0.83 read as disagreement.
- *The thresholds* were fixed from the trap's measurement before any healthy fit was read.
  - Trapped replicates jumped from b 0.003–0.35 over 1–9 structures of 500. Those that reached the posterior left from b 0.98 and up.
  - If each genealogy error is right, the ratio's 3 is passed with probability below 0.3% at any R ≥ 2.
- *Measured clear:* GBASS's t(5) and skewed GH evidence fixtures, R's enumeration fixtures, KHAOS's exact-law fixture and the stage's own fixture. Between-to-within read 0.19–1.44 there, and every last step left from b 0.948 or above.
- *The collapse share alone does not separate.* Healthy small spaces held 2.1–8.8% of N distinct before their last step, against the trapped runs' 0.2–1.8%, so the jump carries the separation in every reading. A small space whose data barely inform could leave from b < ½ over so few structures, and read as jumped.
- *What neither sees:* a wrong region that every replicate reaches alike, by a schedule that never jumps. Only a reference outside the run reads it.

**The step cap** (`StepCap`, opt-in, built 2026-11-09 for KHAOS's trap, `0ffc0e0`). With r set, each step from b > 0 takes b′ = min(the rule's b′, r·b, 1); the first step, from b = 0, stays the rule's alone.
- *Why.* On the trap the rule jumps to b = 1 from a population collapsed at near-empty structures. A capped schedule cannot jump: its last step leaves from b ≥ 1/r, so the jump flag never fires on it.
- *What it costs.* At least ln(1/b₁)/ln r stages, b₁ the first step's b, and the fit's draws differ from the uncapped fit's wherever it binds. `BassSmcReplicate.CappedSteps` counts the steps it bound, read from the stage reports, so a resumed fit counts the steps before its checkpoint too.
- *What it reaches* (the proposal-kernels brief's slice 2; [the residual birth](#bass-the-residual-birth-e4)).
  - On the trap fixture at r = 1.25, the cap alone took 11 of 16 replicates to the posterior, and the cap with E4's residual birth 16 of 16.
  - Their evidence still spread over 21 nats against per-replicate errors under 1, which the disagreement flag reads. The cap settles the posterior's reach, not the evidence.
- *Checked* (`SmcStepCapTests`).
  - Every step from b > 0 is within r·b, and the first is the uncapped fit's. The readout's count matches a recount. At r = 1.05 a hinge fit took 146 and 158 stages against 14 and 15 uncapped.
  - ln p(y) on KHAOS's exact-law fixture at r = 1.25: −72.3035 ± 0.0128 (raw ridge) and −72.2910 ± 0.0172 (with the residual birth) against −72.2910, and −69.5909 ± 0.0172 against −69.5869 (the modified g-prior with the residual birth). The flags were clear, and 19–27 of 43–49 steps were bound.
  - Resume is bit for bit through 25 resumes. Without the cap, SMC's pinned baseline holds.
  - A production fit on the trap fixture (qualification): r = 1.25 with the residual birth, N 500 at 5 sweeps, four replicates. All four reach the posterior, at σ² 0.2245–0.2273 and M 10.12–10.14, none jumping. Their ln Ẑ (−286.4 to −291.1) disagree at 3.66 times their errors, and the flag says so.

**Checked:**
- the kernel's law against enumeration at three b, on both forms and under fixed g;
- forks independent of their parents, and parallel runs equal to serial ones;
- the normalizer against the oracle's exact ln Z: fixture B −28.597 ± 0.027 on the integrated form and −28.593 ± 0.032 on the linear, against −28.589; fixture C −30.888 ± 0.019, against −30.903;
- fixture A's posterior: the count law 0.176, 0.750, 0.074 against 0.176, 0.751, 0.073, E[τ], and all 48 inclusions;
- the weighted WAIC and PSIS-LOO against an equal-weight resample;
- resume bit for bit;
- MCMC fits bit-identical with `Smc` null.

**The qualification** (D10's stage 8 question; exploration). Four replicates of 1,000 particles, with S10's seeds, on data of n ≈ 250, S10's size, 1,000 and 4,000 cases. The particle count was not varied: the chip read its brief's N as the data size, so D10's question of how the particle count matters is still unread.
- *Mixed and functional: SMC finds the better representation.* At S10's sizes every replicate holds the structure the best MCMC chains hold, where one mixed chain and two functional chains sat 22–29 nats below. The replicates agree within 1.3–2.4 nats.
- *Friedman: it does not.*
  - The replicates' log posteriors span −897 to −915, against the best chain's −892.4.
  - They mostly hold three x1·x2 terms where the best chain holds two.
  - They disagree by 13–66 nats, more at larger N.
  - σ's interval misses the truth at N = 250–1,000, lying above it as untempered MCMC's does, and covers it at 4,000.
- *Lineages collapse.* Friedman's replicates end descended from one or two initial particles, and the others thin as n grows.
  - So the comparison between modes happens among the few lineages alive at low b, not near b = 1, where the modes' mass ratios would decide it.
  - The genealogy error saturates at √ln 2 = 0.83 once one family is left. The replicates' spread, 1–21 times larger, is the error to read.
- *Prediction agrees better than structure:* σ²'s spread is 0.04–0.4 posterior sd, while M's and the log posterior's reach 0.6–1.0 on Friedman.
- *The cost:* about 125–300 stages, and 0.5–5 minutes a fit.
- *No default changes:* MCMC stays the default.

**Friedman's band, measured** (2026-10-29; [the measurement's report](../../../issues/regression/bass/friedman-band-brief.md#reports); exploration).
- *A first-order transition on the tempering path.* At fixed b inside the band (0.23–0.28), ℓ₁ is bimodal within each chain.
  - The two phases sit 120–190 nats apart even at a fixed M, so they differ in what the terms are, not in how many.
  - A second split at b ≈ 0.12–0.15 divides M classes (4 and 5 against 6).
  - Λ̂ is 2.0–3.0 over the band and 8.6–10.2 over the whole ladder.
- *SMC does not cross it reliably at any particle count read.* At 1,000, 4,000 and 16,000 particles the replicates' ln Ẑ spreads by 8.4, 7.1 and 6.5 nats, and still rises with N. Lineages collapse inside the band at every N, to 1–10 families by b = 0.32. One 16,000-particle replicate found a basin, four x1·x2 terms with no x1 or x2 main effect, with the largest ln Ẑ, which no MCMC chain visits.
- *Mutation buys more than particles:* 20 sweeps per stage at 1,000 particles keep as many lineages as 16,000 particles at 5 sweeps, at a tenth of the time, but the replicates still disagree by 14 nats.
- *Not the support's doing:* the angle floor refuses nothing, and halving the support minimum leaves the chains apart.
- *So Friedman's structural posterior is not computed reliably by any sampler here,* and its evidence is unresolved. Predictions agree better than structure, and D10's verdict reports the disagreement. *(2026-10-30)* The flat-histogram sampler has since computed it, on three master seeds by 2026-10-31: a1, with ln Z(1) ≈ −841 ([below](#bass-the-flat-histogram-sampler-n17-phase-2-2026-10-30)).

**A path from a pilot reference** (2026-10-29, [its report](../../../issues/regression/bass/friedman-reference-brief.md#reports); exploration, tests only). SMC ran on the geometric path from q = ε·p₀ + (1 − ε)·e to the posterior, with e the pilots' structures (M4's chains, S10's tempered fit and from-prior SMC) and ε 0.05 or 0.5.
- *Agreement, on the pilots' evidence.* At 4,000 particles every arm's replicates agree: sd 0.26–0.58 nats, pooled −886.4 to −887.5. That is, within 0.8 nats, the evidence of e's own structures, −887.15 by quadrature with no move.
- *A lower bound it misses.*
  - An SMC normalizer is unbiased, so by Markov's inequality P(Ẑ ≥ e^t·p(y)) ≤ e^(−t) for each replicate.
  - The largest of R2's twelve from-prior replicates, −855.22, then puts ln p(y) above −864.6 at 99.9% (t = ln 12,000).
  - The arms lie 22–23 nats below that bound, and the shortfall lies inside basin A, not between basins.
- *Why.* e stands 95–175 nats above ε·p₀ on its own structures, so leaving them costs (1 − b)·h at b. The population holds e's structures through b = 0.75 and spreads only in the last stages. No transition is crossed on this path, and 20–174 families survive to b = 1, against 1–5 from the prior.
- *So agreement from such a start measures how well the reference covers the posterior, not ln p(y).* Every arm's 95% interval for σ misses the true σ = 1, which R2's 4,000- and 16,000-particle runs covered.
- *The structure prior's normalization.* `BassTarget`'s structure prior is D02's count weight w, normalized neither by C, the count law's normalizer, nor on F. An SMC normalizer from a normalized reference other than the prior on F, as here, estimates ∫w·p(u)·L, so it subtracts ln C + ln Z_F. An estimator that normalizes against the prior itself cancels both: SMC from `BassPriorSampler`'s draws, as in C4's runs; a ladder's ratio to its reference rung; and the flat-histogram ratio estimator.

### BASS: the flat-histogram sampler (N17 Phase 2, 2026-10-30)

`BassSampling.FlatHistogram` fits BASS by the flat-histogram engine (`docs/Maths/Samplers.md`, "Flat-histogram sampling"). It learns multicanonical weights over ℓ₁, the collapsed log marginal at b = 1, then reweighs the production chains' draws to the posterior. It is null by default, so every other configuration draws as before; every fingerprint moved with the new field, and no draw did. Built by a chip ([its brief and report](../../../issues/regression/bass/flat-histogram-phase2-brief.md#reports)) and verified by the main session.

**The configuration** (`BassFlatHistogramSpec`):
- `Run`, the engine's `FlatHistogramSpec`: chains, the learning's and production's lengths, the start (flat or WHAM), and the window, b_lo, K and the gain where they are not read from the pilot.
- `Pilot`, a `BassFlatHistogramPilot`, for a WHAM start: its `Tempering`, by default the linear form with the reference rung, self-sizing; `Chains` (4); and `BurnIn` and `TracedScans`, both required.
- `MaxIterations`, a per-call budget, as `SmcSpec.MaxStages` is for SMC. A fit short of its lengths is flagged incomplete, reports no evidence, refuses the posterior views, and its checkpoint continues it. A complete fit carries no checkpoint, since the engine refuses a complete one.
- *An iteration is `Sampling.Thin` sweeps* of every chain, MCMC's meaning. Production keeps one draw per iteration, so `Thin` sets memory as well as mixing.
- *Validation refuses:*
  - the spec beside `Sampling.Smc`, beside a top-level `Sampling.Tempering`, or with an R̂ target;
  - a pilot missing where the engine needs one, or present where nothing reads it;
  - a pilot on the integrated form or the powered path, since WHAM needs rung laws prior·e^{b·ℓ₁};
  - an improper τ prior wherever b_lo = 0, since the prior over (structure, τ) is then improper.
  - Jeffreys σ² and the flat intercept run from b_lo = 0, since the prior over (structure, τ) stays proper (D05), but report no evidence.

**The chain** (`BassMulticanonicalChain`). It is `BassChain` on `MulticanonicalTarget<BassState>` over the collapsed target, so D04's sweep runs unchanged on π_w:
- a structure move accepts on D01's proposal ratios and the inner prior, with ln w(ℓ₁) in the likelihood's place;
- τ's slice step targets p(u)·w(ℓ₁(S, u));
- `Reweigh` re-reads ln w from the current evaluation's cached ℓ₁, with no new marginal.
- A production draw is composed at b = 1, as an MCMC draw is, so it keeps its full state (D04). Its functionals (M, τ, the log posterior from the inner prior and ℓ₁, σ², and the fitted values) fall in D10's three groups.

**The fit.**
- `DrawWeights` are the b = 1 weights e^{ℓ₁}/w(ℓ₁), self-normalized over every chain's draws, in the stores' pooled order. Every view reads them as it reads an SMC fit's: the summary, the pointwise contents, prediction, sensitivity, effects, the training fit and the residual check.
- `Traces()` refuses, as for SMC: the draws' order is π_w's, not the posterior's. `BassDiagnostics.FlatHistogram` (`BassFlatHistogramReadout`) takes the traces' place, with:
  - the weights as resolved (window, K, b_lo), t₀ and the gain's start;
  - the learning's visits and weight trace;
  - production's visits and round trips, pooled and per chain;
  - the b = 1 Kish ESS, pooled and per chain;
  - the production draws, for reweighting to any b;
  - the evidence curve, and the pilot's readout, with the nodes a warm start tied.
- *D10's verdict* reads π_w's functionals over production, the chains' own target. The mode separation reads the chains' reweighted means at b = 1.

**The evidence** (`BassFlatHistogramEvidence`). It is the ratio estimator's ln Z(b) on 101 points from b_lo to 1.
- *`BassFit.Evidence`* holds ln Z(1) only from b_lo = 0, under proper σ², intercept and τ priors. From b_lo > 0, the readout holds the ratio ln Z(1) − ln Z(b_lo), flagged, and `BassFit.Evidence` is null, so a ratio never reads as an evidence. The alternative was to report the ratio in `BassFit.Evidence` with a flag, as the brief drafted it; it was not taken, since every other sampler's `LogEvidence` is ln p(y).
- *No normalization is added.* `BassTarget`'s structure prior is D02's unnormalized count weight, but the ratio estimator normalizes against the prior itself, so it gives ln p(y) under the prior on F directly.
- *`BassEvidence.SteppingStone` is nullable,* since a flat-histogram evidence has none; `LogEvidence` and `MonteCarloError` read its own estimate.

**The pilot** (`BassFlatHistogramFit.Pilot`). Its ladders run on the linear form, each from an exact draw of the prior on F, burned as a tempered fit burns them, then frozen and traced.
- Self-sizing ladders each freeze their own rungs, so each ladder is read alone by `TemperingPilot.FromLadders`. The rungs then pool by b: a rung several ladders share holds their draws, with inefficiency Σn/Σ(n/g).
- *Where SAMC's fixed point refuses the pilot's ĝ,* the warm start ties the nodes under ten pilot draws (2026-10-31; `docs/Maths/Samplers.md`, "The warm start"). A start still refused, or a pilot too thin for two nodes to hold ten draws, throws, naming the window, K, b_lo and the nodes the pilot left empty, with the warm start's reason. The readout says whether a start was tied (`BassFlatHistogramPilotReadout.TiedNodes`, 2026-10-31: the nodes tied, empty where ĝ was accepted, the same through a cut and resumed fit). A refusal names the segment where flat occupancy leaves the hull of ĝ's atoms' hats, and by how much.

**Checked** (`BassFlatHistogramTests`; [the report](../../../issues/regression/bass/flat-histogram-phase2-brief.md#reports)).
- *The chain on π_w at fixed weights,* on fixture A's enumerable F, against π_w by enumeration, with τ held and sampled, within 5 batch SE plus small slacks. Every 97th sweep the chain's density equals the inner prior plus ln w(ℓ₁).
- *The fit on fixture A,* under τ held and sampled, both coefficient forms, and b_lo = 0 and 0.3:
  - the learned weights against SAMC's limit on the exact density of states, 0.03–0.10 RMS;
  - the reweighted posterior (M's law, the inclusions, the heaviest structures, E[τ]) within 5 SE of the chains' spread;
  - ln Z(b) on five b against the enumeration within 5 errors.
- *Fixtures B and C:* ln Ẑ(1) −28.529 ± 0.047 and −30.863 ± 0.054, against the oracle's exact −28.589 and −30.903.
- *End to end,* a scalar fit and a direct functional fit: every view on the weighted draws, each draw's weight rebuilt from its own state, and `Traces()` refused.
- *Resume* across learning and production, bit for bit; the improper options with no evidence; a chain that never tunnels, flagged.
- *The warm start:* WHAM's start meets SAMC's limit to 0.097 RMS over eight pilots. The pilot is worth about 8,000 flat iterations (median RMS over 32 chain seeds: the bare start 0.081, flat at 8,000 0.089). From the default gain start the warm learning first gives some back (0.091 at 2,000), then falls more slowly than flat's: 0.069 and 0.057 at 8,000 and 16,000, against 0.089, 0.063 and 0.043 at 8,000, 16,000 and 32,000, the pilot's residual error not worked off. So it leads at equal iterations, on 30 of 32 seeds at 8,000 and 21 at 16,000. One seed had read half the iterations beating flat's full length; that holds on 17 of 32 seeds at 8,000 and 7 at 32,000 (2026-11-01; `ComparativeSweep`), as on the toys.
- *The readout and the window's top* (2026-10-31): the gap pilots' fits, accepted (7803) and tied (7802), are pinned by hash before the top became adjustable, and hold with the options off; 7802's tied nodes, [0, 7, 10, 11, 13, 14, 15, 16], equal the warm start's own list, through a cut and resumed fit. The fixed point's refusals name their segments: gap pilot 7801's at segment 13, by 0.074 of the mass, and 7802's at segment 16, by 0.127. On fixture A with τ held, a top three spacings low rises to F's highest ℓ₁, −15.80, its learned ln w 0.102 from θ* and ln Z(b) within the enumeration's errors. A top margin there adds nodes past F's highest ℓ₁, which no state reaches, and the WHAM start is refused even tied.

**Friedman, qualified** (2026-10-30; exploration, one run per seed; [the report](../../../issues/regression/bass/flat-histogram-phase2-brief.md#reports)). What would count was written into the harness before the run:
- every production chain tunnels between the window's ends;
- the chains' ln Ẑ_c(1) agree within about a nat;
- the pooled ln Ẑ(1) reaches −864.6, the 99.9% lower bound that R2's from-prior normalizers set ("BASS: SMC", "A path from a pilot reference").

*On master seed 8017, from the WHAM start, all three hold,* the first sampler to meet them.
- *Tunnelling.* Every chain made 4–17 round trips across a window of 447–773 nats, which the pilots' ladders, at 16–21 rungs, never crossed. That is fewer than the 20 sized for.
- *Agreement.* The chains' ln Ẑ_c(1) spread 0.34 nats (SD) from b_lo = 0, and 1.09 with the gain counted from 0. The ratio from b_lo = 0.10 spread 0.14.
- *The evidence.* ln Ẑ(1) = −841.4 ± 0.2 (−841.0 ± 0.5 with the gain from 0), 23 nats above the bound and 45 above the reference arms.
- *The posterior it computes* is a1, 4 × x1·x2 + 2 × x3 + x4 + x5, at 66–68%, with a1 + x6·x8 at 11%; four x1·x2 terms carry 95–98%. BASS fits sin(πx1x2) with four x1·x2 products, not with main effects in x1 and x2. A, B, G, S and T, the basins every MCMC chain and the reference path held, carry 1.3% or less together. σ's 95% interval covers the true 1. Only one of R2's 16,000-particle replicates had found a1.
- *Cost.* 13–15 minutes an arm, 45M sweeps a chain, against 9–11 minutes for R2's 16,000-particle SMC, whose replicates disagreed by 6.5 nats.

*Through SMC to b = 0.10* (a reading added after the run). The ratio ln Z(1) − ln Z(0.10) agrees between the two windows (−715.25 ± 0.07 against −715.5). But SMC from the prior on the linear form reads ln Z(0.10) about 5.5 nats low, so that arm's −846.6 is low by as much. Below b = 0.05 the flat histogram, SMC and exact importance sampling from the prior agree to a tenth of a nat, which pins the ratio estimator's denominator; above it importance sampling falls first and SMC next, as estimators missing rare mass do.

*On both forms* (2026-10-31; a reading, written before its run; the loose ends' chip, [its report](../../../issues/regression/bass/flat-histogram-loose-ends-brief.md#reports)). SMC ran from the prior with a stage forced at b = 0.03, 0.05, 0.07 and 0.10 (4 replicates × 4,000 particles, C4's defaults), at 5 sweeps a stage and at 20. Its references are importance sampling from 200,000 exact prior draws and seed 8017's flat-histogram run, rerun bit for bit. On the integrated form, D05's SMC path, that run's 500,000 draws are reweighted by the form's tempered marginal, Ẑ(b) = Σ p_b(S, τ)/w(ℓ₁) over Σ 1/w(ℓ₁); every draw's rebuilt ℓ₁ equalled its level, and on the linear form the same code reproduces the curve to 6 digits.

| b | form | importance sampling | flat histogram | SMC, 5 sweeps | SMC, 20 sweeps | families, 5; 20 sweeps |
| --- | --- | --- | --- | --- | --- | --- |
| 0.07 | linear | −98.84 ± 0.50 | −97.61 ± 0.01 | −97.89 ± 0.12 | −97.73 ± 0.03 | 478–524; 881–1,316 |
| 0.10 | linear | −136.52 ± 0.85 | −125.93 ± 0.06 | −131.51 ± 0.22 | −129.83 ± 0.99 | 188–217; 69–639 |
| 0.07 | integrated | −113.19 ± 0.10 | −112.95 ± 0.01 | −112.95 ± 0.04 | −112.94 ± 0.01 | 1,359–1,428; 2,081–2,149 |
| 0.10 | integrated | −153.66 ± 0.33 | −150.17 ± 0.02 | −150.59 ± 0.14 | −150.58 ± 0.03 | 442–576; 823–855 |

- *The 5.5 nats are the linear form's:* 5.6 at b = 0.10. Twenty sweeps a stage narrow the loss to 3.9 nats without closing it. One replicate came within a nat of the flat histogram as its families fell to 69, while the other three stayed 4.2–5.3 low: rare high-likelihood mass that a population reaches by chance, not mixing at a fixed b that more sweeps repair.
- *On the integrated form SMC loses about 0.4 nats by b = 0.10,* keeping 442–855 families, and twenty sweeps tighten its replicates without moving that. At 0.07 and below it agrees with the flat histogram within 0.01. Whether the 0.4 nats is SMC's or the seed's is not settled: the reference is one run reweighted (Kish ESS 33,452 at 0.10, 2,991 at 0.03).
- At b = 0.05 and below every estimator agrees on both forms. Importance sampling reads low first above it, SMC next, as before.
- The integrated form's Z(b) is the lower of the two between b = 0 and 1: it tempers the likelihood inside the a, σ² integral, E[L^b], which Jensen's inequality puts below the linear form's (E L)^b.

*A second seed, before the warm start was made robust* (the main session's, 2026-10-30).
- On master seed 8217, SAMC's fixed point refused the pilot's ĝ: no node was empty, but one stood 125 from flat. So the WHAM warm start is not robust across seeds on Friedman.
- Rerun from a flat start over that pilot's window, at twice the sweeps, the learning made no round trip in 100,000 iterations, and production none: flatness 0.98, b = 1 Kish ESS 8 of 500,000, no draw above the window. By the first criterion it reads nothing.
- What it does show agrees with seed 8017 without confirming it. Its ln Z(b) curve matches to within a nat up to b = 0.8, and every chain found a1. But its chains disagree on a1's share (46–81%), and its ln Ẑ(1), −843.6 ± 0.5, rests on a handful of draws.

*Replicated* (2026-10-31; exploration, one run per seed; the warm start's chip, [its report](../../../issues/regression/bass/flat-histogram-warmstart-brief.md#reports), verified by the main session). Seed 8217's refusal was an ill-posed limit, not a failure of Newton's: flat occupancy lay outside the hull of the pilot's atoms' hats ("The warm start" in `docs/Maths/Samplers.md`). With the thin nodes tied, seed 8217 ran by its test as committed before Phase 2's run, and a third master seed, 8417, declared before any run, by the same criteria at the same settings. Both meet every criterion: every chain tunnels; the chains' ln Ẑ_c(1) SD is within a nat, with a1 passing the χ² across chains; ln Ẑ(1) reaches the bound; and against seed 8017, ln Ẑ(1) lies within a nat plus both runs' errors and a1's share within 5 points.

| | seed 8017 (`Thin` 200) | seed 8217 (`Thin` 400) | seed 8417 (`Thin` 400) |
| --- | --- | --- | --- |
| the warm start | accepted as it stood | nodes 0 and 1 tied | accepted as it stood |
| production round trips per chain | 9–17 | 14–17 | 6–21 |
| the chains' ln Ẑ_c(1) SD | 0.34 | 0.15 | 0.63 |
| pooled ln Ẑ(1) | −841.36 ± 0.17 | −841.55 ± 0.07 | −841.14 ± 0.32 |
| a1; a1 + x6·x8 | 65.9%; 10.6% | 66.6%; 9.8% | 68.4%; 10.1% |
| b = 1 Kish ESS of 500,000 | 34,256 | 10,937 | 7,684 |
| time, pilot to production | 13 min | 39 min | 36 min |

- The ln Ẑ(b) curves agree within 0.4 nats at every b read, σ's 95% interval covers 1, and D10's verdict on π_w's functionals is converged, on each seed. Four x1·x2 terms carry 96–97%; A, B, G, S and T together 1.3% or less.
- *Twice the sweeps bought more crossings, not 20 a chain:* one chain of the eight reached the 20 the lengths were sized for.
- *The b = 1 Kish ESS fell* to 32% and 22% of seed 8017's. A reading: these pilots' windows top out higher (−773.6 and −771.7, against −778.3), so fewer production draws fall above the window, where the extension holds the posterior's own law. The top was then read as two options on seed 8017 (below).
- *Reported, not ruled:* on seed 8417 the unlabelled a1 + x6·x8's share disagrees across chains (13.8 to 6.7%, χ² p < 0.001); the criteria read labelled basins only.
- *A flat start at t₀ = 100·K* (a reading, one run, on seed 8217's window) learned the window in the warm start's budget: ln Ẑ(1) = −841.38 ± 0.27, a1 at 67.1%. Its gain was still at t₀/t = 0.13 when learning stopped, and production's flatness, 1.36, was two to four times the warm starts'. At t₀ = 10·K it had learned nothing.

*So* the flat histogram computes Friedman's structural posterior on three master seeds, within 0.41 nats and 2.5 points of one another, each run's error resting on four chains. It stays an option; no default changed.

*The window's top, read on seed 8017* (2026-10-31; exploration, one run each at the qualification's lengths, written before the runs; the loose ends' chip, [its report](../../../issues/regression/bass/flat-histogram-loose-ends-brief.md#reports)). The replications' pilots topped out higher and lost b = 1 draws; seed 8017's pilot never found the posterior's basin, and its learning spent 29% of its visits above the window, where SAMC does not update. The engine's two options for the top (`docs/Maths/Samplers.md`, "The run"), against today's window:

| | today's window | `ExtendTopUntil` 50,000 | `WindowTopMargin` 3 |
| --- | --- | --- | --- |
| the window's top; K | −778.29, fixed; 204 | −751.65 by iteration 4,000; 211 | −767.73, fixed; 206 |
| the start | accepted as it stood | the same | tied at both ends (its top nodes held no pilot draw) |
| learning: visits above the window | 29.2% | 0.1% | 1.1% |
| production: round trips, all chains (fewest) | 47 (9) | 3 (0) | 24 (4) |
| production: draws above the window | 6.4% | 0.0% | 0.7% |
| b = 1 Kish ESS of 500,000 | 34,256 | 7,462 | 14,543 |
| ln Ẑ(1) (chains' SD) | −841.36 ± 0.17 (0.34) | −841.45 ± 0.19 (0.39) | −841.11 ± 0.27 (0.54) |
| a1; a1 + x6·x8 (χ² p across chains) | 65.9%; 10.6% (0.197) | 69.6%; 8.5% (0.518) | 65.5%; 10.9% (< 0.001) |
| time, pilot to production (a shared machine) | 15 min | 16 min | 15 min |

- *Raising the top empties the learning's inert visits and costs production's b = 1 draws,* to 22% and 42% of today's ESS: the trade the replications foresaw.
- *The extension as built climbed 26.6 nats within 4,000 iterations,* into a tail of rare states. Production then rarely reached the top bins and made 3 round trips, two chains none, failing the qualification's first criterion, though ln Ẑ(1) agrees within 0.1 nat and a1 within 4 points.
- *The margin tunnels,* at 42% of today's ESS and within 0.25 nats, but a1 + x6·x8 parts across its chains (6.6–14.8%), as on seed 8417.
- So neither option buys a reading today's window lacks on this seed. A gentler extension rule is the open variant. No default changed.

### BASS: the pointwise density and the residual check

**The pointwise density** (`BassPointwise.LogDensity(model, terms, draws, units, content)`, `BassFit.PointwiseLogDensity(data, units, content)`; S9) is OQ17's draws × units matrix, with the chains' lengths for PSIS's r_eff. A unit is an observation of a scalar response, or a case of a direct one.
- *Retained parameters* (the default): the Normal density at the draw's coefficients and σ², every latent conditioned on.
- *The collapsed predictive:* the draw's structure and τ kept, a and σ² integrated against their posterior given the observations outside the unit, a multivariate t (`LeaveOutLogDensity`), computed once per distinct structure's design. The units' quadratic forms in the Gram are τ-free, so they are taken once per structure and serve every draw's τ (`LeaveOutQuadratics`), bit for bit; a draw at a τ already seen (every draw of a structure under fixed g) copies its row. A unit must leave an observation outside it, beyond the one the flat intercept spends. S10's leave-one-out check (the qualification section) reads which content's k̂ to trust.
- *Under a noise law* the retained content is the law's density with each observation's latent integrated out, reported `Integrated` with observation scope: its prior is independent of every other observation's given the globals, so the term is in sample, WAIC reads it and PSIS-LOO has an lppd to subtract from. Until 2026-11-06 it was reported `IntegratedGivenOthers`, which made `LeavesUnitsOut` true for every noise-law fit, so WAIC refused them and p_loo was NaN.

**The candidate density** (`BassFit.CandidateLogDensity(data, cases, grids, strata, unseenLevels)`; E16's seam for the [conformal layer](../ConformalPrediction.md#bass-the-candidate-density)) is log p(y | x_t, draw) for every retained draw and every value y of each new case's grid, in the pointwise density's rows, with its chain lengths, weights and latent report (`BassPointwise.Latents`, shared). The surface at the cases is prediction's, and the density is the training terms': N(y | f, σ²), or the noise law's with each latent integrated. A candidate is one observation of a scalar response. A weighted fit, a draw whose working response moved, and a direct functional fit are refused.

**The residual check** (`ResidualDependence`; F11, S9 and S10), over the posterior-mean residuals of a direct functional fit. A model whose noise is independent across locations (D07) reads about its expectation; a misfit smooth across the locations, as a deterministic simulator's is, reads high, and then the intervals come out too narrow.
- *Along an order* (`ResidualDependence.LagOne`): each case's r₁ = Σₖ(eₖ − ē)(eₖ₊₁ − ē)/Σₖ(eₖ − ē)², about its own mean, so a case's offset does not read as dependence; about −1/m under independence, sd about 1/√m.
- *Over output coordinates* (`ResidualDependence.MoransI`): Moran's I per case over each location's K nearest neighbours (binary weights, ties to the lower location), −1/(m − 1) in expectation under independence.
- *The index is declared on the call* (`ResidualIndex`, decided 2026-10-28): along an output variable, along an explicit order, or nearest neighbours in chosen output variables' scaled coordinates. A fit with one output variable reads along it by default, the fit's `ResidualDependence` field; with several, no index is decided until one is declared (`BassFit.ResidualDependenceFor`), since an output variable may be a covariate known at each location and not a coordinate.
- *The posterior mean goes by term:* each term's coefficients summed over the draws, then its column once, O(S·M + T·N) for T distinct terms, in place of every draw's fit, O(S·M·N).

**Checked.** The retained content against the Normal density at the explicit design, and the collapsed content against explicit refits that leave each unit out (10⁻⁸), over 54 combinations of intercept, σ² prior, g mode and data kind. The residual check hand-worked (r₁ and Moran's I), against a direct computation under equal and given weights, each declared index against its direct reading, a case-only fit of a deterministic sine reading exactly the sine's r₁ (0.9458), and independent noise reading its expectation: r₁ −0.089 ± 0.044 against −0.062, and Moran's I on a 5 × 5 output grid −0.041 ± 0.019 against −0.042, where a smooth misfit reads 0.43.

### BASS: qualification

S10 (2026-10-28) ran BASS's defaults against known truths. Each run sits behind `BASS_QUALIFICATION` = 1 in `tests/regression/BassQualification*Tests.cs`, and the [checklist's report](../../../issues/regression/bass/port-fitting-checklist.md#reports) holds what each read. R's posteriors are exploration only (the parity rule).

**Which content's k̂ to trust** (D04's leave-one-out check). On fixture A's design (N = 27), where F enumerates, the exact leave-one-out predictive is logΣ_S w_S − logΣ_S w_S/p(yᵢ | y₋ᵢ, S), the prior held at the full data's as PSIS holds it, so PSIS on both contents is scored unit by unit with no Monte Carlo error in the reference. Three responses (fixture A's surface, an outlier, a corner spike), four chains of 1,000 and 2,000, pre-registered verdicts as BARS's (`LooBenchmark`). Under fixed g, where the collapsed ratios clump exactly:
- *the retained content's k̂ is no false negative:* none of the 80 units it passes misses by |z| > 3 (the largest |z| 2.5), one unit flagged;
- *the collapsed content's k̂ is pessimistic:* it flags 72 of 81 units, all 72 within |z| ≤ 3; its estimates are right where its k̂ says they are not.

Under sampled g, the default, each draw at its own τ:
- *the retained content's k̂ is no false negative:* none of its 79 passed units misses by |z| > 3 (the largest 2.2); it flags two, both within |z| ≤ 3;
- *the collapsed content's k̂ is pessimistic, less so:* it flags 21 of 81 units, all 21 within |z| ≤ 3, and none of its 60 passed units misses (the largest |z| 1.6). τ's spread breaks the exact clumping.

So the retained content's k̂ is the diagnostic to read under either g, as BARS found. Under sampled g the reference integrates u = ln τ for each of F's 464,094 structures, in parallel (about 40 minutes). A total's error is NaN whenever a unit is flagged, as loo reports it.

**Recovery** under the defaults (interactions to 3, the count prior, the proper σ² and level centring, sampled g), four chains of 5,000 burn-in and 2,000 draws:

| surface | N | RMSE | 95% coverage of f | σ's 95% interval | mean M | largest R̂ |
| --- | --- | --- | --- | --- | --- | --- |
| Friedman, p = 10 (5 inert), σ = 1 | 500 | 0.335 | 91.6% | [1.013, 1.161] | 9.6 | 2.35 (M) |
| Friedman, tempered (self-sized, 20–23 rungs) | 500 | 0.282 | 97.2% | [0.991, 1.139] | 10.2 | 1.45 (M) |
| mixed categorical, σ = 0.5 | 300 | 0.107 | 98.0% | [0.472, 0.558] | 5.8 | 1.54 |
| direct functional, σ = 0.2 | 1,200 | 0.021 | 99.0% | [0.188, 0.203] | 6.1 | 2.55 |

- *The chains disagree in M.* Each untempered Friedman chain held a narrow band (10–12, 8–10, 10–12, 9–11 over 2,000 draws): births and deaths accepted 1.5% of attempts and changes 5%, about twenty changes of dimension per chain. R BASS 1.3.1 on its own Friedman sample behaves alike (chains at 10–12, 8–11, 9–10, 9–10; σ's intervals missing 1 in three of four), so this is the method under its defaults, not the port.
- *Friedman's σ reads high untempered,* its interval missing 1, the pre-registered criterion failed: the posterior absorbs the fit's bias (1 + 0.33² gives σ ≈ 1.05) and the chains do not cross between structures. Tempering brings R̂ of M from 2.35 to 1.45 and covers σ, but completes no round trip. Two things bound it. Every chain's ladder holds one steep band, b ≈ 0.16–0.32, whose pairs reject 0.94–1.000 even split to Δb ≈ 0.04: the tempered posterior changes sharply there, and Λ̂ reads such a band censored, each pair counting at most 1, so the self-sizing places too few rungs in it. And the cold rungs' local moves hold each chain in its own band of M, as untempered (chip E saw this alone on its fixture, round trips 0.0001–0.0008 per scan against an ideal of a few hundredths).

**Count against column mass** (D06's fixture). On 200 points with an edge kink of 24 points, which the count admits and the mass (about 12) refuses, at noise sd 0.3 and 0.02: the posterior sd at the edges is the same under both rules (0.06 and 0.006), the RMSEs within about 0.01, and column mass uses about one more term. Neither rule shows the edge instability here; the count default stands.

**Cost** (µs per sweep, one thread, untempered, composition included; a rich low-noise surface, the cap read back as the mean M):

| n \ M (p = 2) | ≈ 5 | ≈ 20 | ≈ 45–50 |
| --- | --- | --- | --- |
| 200 | 68 | 129–154 | — |
| 1,000 | 89 | 190 | 376 |
| 5,000 | 186 | 465 | 1,184 |

The time grows about linearly in N·M at large n, over a fixed 50–60 µs per sweep (τ's slice, composition, bookkeeping). With p = 10 the realized M stayed below its cap at n ≤ 1,000 (3 and 11 terms), and at n = 5,000 a sweep cost somewhat more at the same M (544–659 against 465 at M ≈ 20). R BASS took about 630 µs per sweep on Friedman at n = 500, M ≈ 10, where this baseline puts C# near 100. A tempered run costs about its rung count times this.

Most of that N·M was the composition's: each retained draw re-evaluated every term's column for its mean (F-M2's Phase 0, M1). Since `9ed600a` the chain's moves keep each term's mean, and every draw is bit-identical. Rerun back to back on one machine (p = 2), a sweep at n = 5,000 went from 117, 358 and 694 µs to 69, 88 and 148 (M ≈ 5, 20, 50), and at n = 1,000, M = 20, from 207 to 88; at n = 200 the gain is 1.2–2×. The rerun's earlier times sit below the table's (694 against 1,184), so compare within the rerun. In a Release build (F-M2's pass, 2026-10-29) the defaults took about 7–9 µs per chain sweep at n ≤ 500 and 25 at n = 2,000, the whole fit over its sweeps. The builds behind the table and the rerun were not recorded, so compare only within one build.

**Cost, qualified** (the cost half of performance qualification, 2026-10-31; a chip, [its brief and report](../../../issues/regression/bass/performance-cost-brief.md#reports), verified by the main session). Release, one thread, `Thin` = 1, each point read at chain 0's burned state by M1's timing method. The grid moves one axis at a time about n = 1,000, p = 10, M_max = 20 and K = 2, on a rich surface whose chains reach their cap, with S10's surface along p beside it. It supersedes the table above for today's costs.

| axis | µs per sweep |
| --- | --- |
| n (M_max 20) | 200: 7 (M 7); 1,000: 22; 5,000: 57; 10,000: 217; 50,000: 1,066 |
| M_max at n = 10,000 | 5: 82; 20: 217; 50: 255; 100: 741; 200: 1,388 |
| p and K | flat: 17–25 over p = 2–200 (past p = 50 at a lower M), 20–22 over K = 1–5 |
| L (n = 4,096) | 67–118 over L = 2–1,024, not monotone |
| direct n_x × m | 32 at 100 × 10; 877 at 500 × 1,000 |

- *The cost model at the cap:* a sweep is about 5–8 µs (τ, a composition, a death, the target) plus half a change plus the commits. At the cap the chain proposes no births, so half its moves are deaths and half changes. A change is a column, n(ō + 4) element passes at 2–2.5 ns, plus its cross products, n(M + 3) multiply-adds at 0.3 ns in cache and 0.9 ns beyond: 54–61% of a sweep at n ≈ 1,000 and 79–96% from n = 10,000. Below the cap, births join it, at 16–28% of a sweep.
- *The composition no longer dominates.* M1 read it at 34–94% of a sweep. Since `9ed600a` it costs about 2 µs a draw at M = 20 and 44 at M = 200; with τ it is 18% of the base sweep, 36% at n = 200, and under 3% from n = 10,000.
- *A direct functional fit's commit* rebuilds the residual explicitly (`Refit`, O(N·M) per accepted move, N = n_x·m): 30–56% of a sweep from N = 20,000, 19.7 ms a commit at 500 × 1,000.
- *What runs at the machine's rate.* The production dot product runs at the vector rate, 0.28–0.46 ns per multiply-add, and so do the cross products in cache; beyond the caches they run at the DRAM rate. Columns run at 3–7× the machine's rate, each kernel making its own pass. Categorical membership costs 9–65 ns per case from L = 2 to 1,024, a binary search of the sorted levels, against about 2 ns for a hinge factor. Commits and factor work run near the scalar rate, and the rebuild cadence costs at most 0.12% of a sweep.
- *Memory:* 9.6 KB allocated a sweep at the base, growing as the change's column, 8n bytes; GC pauses take at most 6.5% of a sweep anywhere. A retained draw holds about 120 B plus 13 B a term, and a fit keeps no O(n) state. The peak working set reached 171 MiB (n = 50,000; M = 200 at n = 10,000) over an 85 MiB host.
- *The outputs, per draw:* prediction 0.06–0.13 µs per case at p ≤ 100; Sobol 2.9–39 µs (M from 5 to about 48); basis Sobol superlinear in the rank, 20 µs at r = 1 and 1,459 at r = 30; effects 4.5–9.6 ns per cell.
- *The other backends, at the base:* a tempered rung-sweep costs 14–17 µs, and a swap round 11–24% of a scan. An SMC particle's mutation step costs about twice the untempered step, and a fork more than a mutation. The flat histogram's own work is under 1% of an iteration. Calibration's σ² proposal (N15) is 97–99.8% construction, the Jacobi thin SVD; the capacitance Cholesky's count is 7–8× below a single Jacobi sweep.
- *Parallel and R, as context.* Four chains in parallel ran 2.0–3.4× faster than serially, their draws identical. R BASS 1.3.1, run as written (`truncdist` pinned in `r/`'s renv, `91155ea`) under R's reference BLAS, took 691 µs a sweep on Friedman at n = 500, and 825 and 525 at two grid points (n = 5,000; direct 100 × 50) that the C# read at 57 and 35. R's figure covers its whole run, so it is context, not a comparison.
- *How far to trust it.* No other test host ran during a timed block, but other processes kept 1.6–6.9 cores busy, and two full passes 25 minutes apart differ by a median 1.29× (0.70–2.60). Read absolute times to about ±30% and the parts' shares to about ±5 points. The timed parts leave a rest of −18% to +25% of a sweep, positive at small sweeps, where fixed work outside the timed parts shows. The grid reads chains at their cap by design; the order axis held mostly one-way terms (mean held order 1.05 at the base); the sensitivity moment work alone is not isolated.
- *The efficiency pass's levers, by share of the base sweep:* prefix sums for a main-effect hinge's cross products, which would also skip the proposal's column until acceptance, 35% (16% for the cross products alone; 43–62% with the column at n ≥ 10,000 and M ≥ 50); the column kernels, 25%; the composition and τ, 18%; the direct commit's explicit residual, 4% (first in large direct fits, at 56%); the allocation, `CentredBorder`'s re-centring and the rebuild cadence, a few percent or less. The first two overlap on the column, so their shares do not add. A vectorized dot product is worth nothing. Nothing was built in that report, and no default changed; the pass that built them is the next bullet.
- *The efficiency pass, built (2026-11-04; [brief and report](../../../issues/regression/bass/efficiency-pass-brief.md#reports)).* The five levers, in the recorded order. **Bit for bit** (against verbatim copies of the old forms over seeded sweeps with far ranges, specials and NaN, and with hardware intrinsics off): a term's hinge factors in one pass each in four lanes, with the column's support count, centring and whitening in lanes (lever 2); a categorical subset's membership from a bitset over its level codes, 1.3–2 ns per case against the binary search's 9–65 (lever 3); the commit's explicit residual in column blocks of four lanes (lever 4); and the terms of τ's marginal that neither τ nor the structure moves, `g₁ln g₂`, `lnΓ(g₁ + ν/2)`, `lnΓ(g₁)` and the count prior's collection weight by M, kept once (lever 5). **To a tolerance** (lever 1): a one-way hinge's Gram entries against the held columns are looked up in per-(column, variable, sign) tables, each built by one anchored scan over the cases in the variable's order, within 4·K·ε + 8ε·log₂ n of the columns' norms (93 ε at most, observed), in a plain chain of a scalar fit with at most 16 continuous predictors and 500 cases or more. The tables are a pure cache, so a restored chain computes the numbers a live one does; ladders, SMC, direct fits, a prior apart from the Gram and terms of order two keep the dots. Over 160,000 sweeps at the base the structure sequences with tables on and off were identical, and the two untempered Friedman pins were re-baselined (`BASS_TABLES=off` restores the old bits). Read on the thread's cycle clock under load (indicative; no quiet window came), the base point's sweep fell to about half, L = 1,024 to 0.53 and direct 500 × 1,000 to 0.59 of their times, the commit to 0.24–0.30 and τ's block to 0.40–0.50. The wall-clock pairs await a quiet machine.

### BASS: decomposition, components and prediction

Stage 5 (2026-10-29, D07) adds R's `bassPCA`, `bassBasis` and `predict`. The code is under `spline/bass/output/` and `spline/bass/posterior/`.

**Decomposition** (`BassDecomposition.FromPca(y, BassPcaSpec)`, `FromBasis(K, W, offsets, scales, …)`):
- *PCA.* The responses are centred by their column means y_m and scaled by y_s: 1 by default, or the column sds with `ScaleByColumnSd`.
  - The [thin SVD](../LinearAlgebra.md#thin-svd-thinsvd) Y_cᵀ = U D Vᵀ gives the basis K = U_r D_r and the coefficient responses W = V_rᵀ, R's convention.
  - Under D03's scale-free defaults a component's fit does not depend on W's scale. The convention only fixes the units each component's σ_k is read in.
- *The rank.*
  - It is the smallest r whose cumulative variance share reaches `VarianceShare`, capped at the numerical rank. The comparison is ≥, where R's strict one fails at 100%.
  - Every rank is judged against one tolerance, max(n, m)·ε·d₁, and a share of 1 is the numerical rank itself.
  - A constant output column is centred, with y_s = 1.
  - For comparisons, each component's sign makes its largest |K| entry positive. Ties in the singular values are compared through K Kᵀ and K W.
- *A supplied basis.* K is m × r and may be nonorthogonal. It comes with prepared coefficient responses W, offsets and scales.
- *The representation* is a `ReconstructionMap`; only `LinearReconstruction` is built.
  - It reconstructs y = y_m + y_s ∘ (K w + e).
  - Its `ReconstructionDomain` says where it can be evaluated: a PCA basis exists only at the training locations.
  - It keeps a `BasisUncertainty` slot, null for now.

**Truncation residuals** (D07). A predictive draw needs a new case's error outside the basis.
- *The default, `LeaveOneCaseOut`.*
  - Case i's residual is taken against the rank-r basis built without it, centred by the other cases' mean, with the rank capped at the other cases' numerical rank.
  - Every quantity lies in the span of the training cases, so the n refits run in case coordinates from the one SVD: O(n⁴) for the refits, then O(n·m·ρ) to map back.
  - Under the mean centre, the left-out case's coordinates are n/(n − 1)·cᵢ and the others' cⱼ + cᵢ/(n − 1).
- *`InSample`,* R's, measures each case against a basis it helped build. It understates the error, and at full retention it is zero.
- *The diagnostic.* The ratio of the two RMS values is the basis-stability diagnostic: near 1 when the basis is settled (1.10 at n = 40 on O-O4's fixture), large as r nears n − 1 (2.81 at n = 5).
- *A supplied basis* takes in-sample residuals from its training responses, or residuals given with it, since it has no construction to refit.

**Components** (`BassComponents.Fit(decomposition, predictors, config, …)` → `BassBasisFit`):
- *Each coefficient response* is fitted by `Bass.Fit` as a `ScalarResponse` on the shared predictors.
  - All components share one config by default; a list of configs is the option.
  - Component k takes child k of the master seed.
  - Components, chains and replicas share one parallelism bound P: max(1, ⌊P/C⌋) components run at once, each with min(P, C) chains.
- *Joint draws* (`BassJointDraws`).
  - Independent fits' equal-weight draws pair by index.
  - Under the R̂ stop the fits end at different lengths, so each is thinned evenly to the shortest, keeping draw ⌊(2j + 1)L/(2L′)⌋ within each chain. The components must share their chain count.
  - Weighted draws (SMC, stage 8) are resampled within each component, multinomially by default, down to the smallest count. They are then shuffled, so that the ancestors' sorted order does not pair across components.

**Prediction** (`BassFit.Predict(data, cases, options)` → `BassPrediction`; `BassBasisFit.Predict` → `BassBasisPrediction`):
- *Latent draws* evaluate each draw's composition coefficients (D04) at the cases.
  - That is f(x), or f(x, t) for a direct functional fit, at the training output locations or at new ones (`OutputLocations`).
  - Arrays are draw × case × location, in the requested draw order, with singleton axes kept.
  - At the training inputs they equal the stage-4 fitted values, bit for bit.
- *Predictive draws* add the noise.
  - For a direct fit: N(0, σ²) independently for each case, location and draw.
  - For a basis fit: y = y_m + y_s ∘ (K(w(x) + ε) + e). Here ε_k ~ N(0, σ_k²) for each case and draw, and e is a whole training residual vector drawn uniformly for each case and draw, which keeps its within-output correlation.
  - `BassPredictOptions.BasisParts`, a `BasisPredictiveParts`, takes ε alone (`Noise`, R's defaults) or e alone (`Truncation`), R's `nugget` and `trunc.error`, and the covariance follows (2026-11-02). It shows which source drives an interval's width. R's pair of neither is `Latent`. D07's full draw stays the default, and a part left out draws nothing from the stream. `BassBasisPrediction.Parts` records the choice, and a direct fit, with no truncation residuals, refuses any other.
  - R reuses one ε across cases and indexes its residuals by the row count (O-O1).
- *Inputs* pass through the fit's transform, with columns matched by name.
  - Values outside the training range extrapolate the hinges linearly and are flagged per case; new output coordinates outside the range are flagged per location.
  - An unseen level is refused, or placed outside every subset (`UnseenLevel`).
- *Intervals* are quantiles of the predictive mixture.
  - For a direct fit, the Gaussian mixture's exact quantile.
  - For a basis fit, the weighted quantile of its draws (`DrawSample`, type 7 under equal weights), since its mixture runs over draws × residuals.
  - A basis fit's covariance in physical units is diag(y_s)·C·diag(y_s).
- *The mean surface* (`MeanSurface`, off by default) averages each draw's conditional mean E[f | S, τ, y].
  - At the training inputs under the mean centre that mean is ȳ + (ŷ_OLS − ȳ)/(1 + τ).
  - The average is a Rao–Blackwellized estimate of the posterior mean, never a draw.
- *A call takes its own `Seed`,* and every operation takes the draw weights.

**The density for basis fits** (OQ17; `BassBasisFit.PointwiseLogDensity(units, space)`) is by case.
- *By default, in coefficient space:* Σ_k log N(w_ik | w_k(xᵢ), σ_k²). It is exact for comparing fits that share a basis.
- *`BasisDensitySpace.Physical`* is the declared approximation for comparing across bases.
  - It adds a diagonal Gaussian for the truncation residual, whose mean is the residuals' mean and whose variance is their uniform-draw variance, and subtracts Σ log y_s.
  - It is refused where that variance is zero.

**Checks**, against O-O1's and O-O4's corrected expectations:
- the reconstruction reproduces Y to 2·10⁻¹⁶, where R's misses by 3.09;
- the mixture intervals are ±2.77 and ±2.64, against R's ±3.91 and ±5.20;
- the component noise is uncorrelated across cases (correlation 0.006, R's 1);
- both residual sets match the oracle to 2·10⁻¹⁴ at n = 2, 5, 10 and 40, and at full retention.

### BASS: sensitivity

Stage 6 (2026-10-29, D08) adds R's `sobol` for direct fits and `sobolBasis` for basis fits, under `spline/bass/sensitivity/`.

The same moment service carries E8's derivative moments and C = E[∇f∇fᵀ], the active-subspace matrix, on its own page: [Active subspaces](../ActiveSubspaces.md).

**What is decomposed.** For each posterior draw, the variance of the latent mean surface under a declared product measure over the inputs, with observation noise outside it. That gives each index a posterior distribution. The surface is a declared input (`DecomposedSurface`); only the mean is built.

**Measures** (`SensitivitySpec.Measures`), by predictor name and in each input's own units, as R declares them.
- *Continuous inputs:*
  - `UniformMeasure`, by default over the box: the declared domain, else the training range;
  - `MixtureMeasure` of `NormalComponent` and `StudentComponent`, truncated to its interval. A Student component takes a location-scale `Scale`, with a constructor from a marginal sd for ν > 2. The weights apply before truncation by default (R's convention), or `ConditionalOnTruncation`;
  - `EmpiricalMeasure`: the training values, equally weighted.
- *Categorical inputs:* `CategoricalMeasure`, level probabilities by name, uniform over the training levels by default. R ignores supplied probabilities (O-O2).
- *A support beyond the training range* is allowed, and flagged; an undeclared measure over a declared domain wider than the data is one.
- *Which draws* (`SensitivitySpec.Draws`, R's `mcmc.use`, 2026-11-02): indices into the pooled draws, a basis fit's joint draws, as `BassPredictOptions.Draws` takes them: in the order wanted, repeats allowed, each keeping its weight. Effects evaluate only those draws, so a selection fits a `CellBudget` the whole run would exceed.
- *Functional responses.*
  - `FunctionalForm.Integrate` treats the output variables as more inputs. Their `OutputMeasures` are uniform over each range by default, or `EmpiricalMeasure` for the training grid, taken as the product of each variable's marginal.
  - `Condition` decomposes at declared output coordinates, by default the training locations. The coefficients become aⱼ·hⱼ(t), and each coordinate has its own decomposition.
  - `ConditionOn` names the output variables held (2026-11-02; R's `func.var`, which holds one; several may be named). The others stay inputs under `OutputMeasures`, as R keeps them. hⱼ then carries the held factors with the output side's normalizer, one constant over the whole part, and a term left with no factor on an input folds into the level. A coordinate holds one value per held variable, in the order `CoordinateVariables` gives on `SobolAnalysis` and `EffectAnalysis`, by default their distinct training values. A basis fit refuses it, since a basis has no output coordinates.
  - Both forms report main and total indices; R's functional path omits the totals.

**Moments** (`BassMoments`):
- *A term is κ·Π_v φ_v,* a product of one-variable factors: hinge powers and subset indicators.
  - Each factor is a `PiecewisePolynomial`, split at its knots and integrated piece by piece against the measure's interval moments, about the centre of each piece clipped to the support.
  - Centred products are built directly, so a covariance is integrated as one polynomial, never as E[φφ′] − E[φ]E[φ′].
- *Under the product measure,* a term's component on a subset u of inputs is Π_{v∈u}(φ_v − Eφ_v)·Π_{v∉u}Eφ_v. Two terms' components therefore covary only on subsets of their shared variables, with Cov_u = κκ′·Π_{v∈u}Cov(φ_v, φ′_v)·Π_{v∉u}Eφ_v·Eφ′_v.
  - The variance, each subset's main variance and each input's total are sums of aⱼaₗ times these cached pair shares. An input's total sums the subsets that hold it.
  - Draws repeat their structures, so factor, term and pair moments are cached across draws, and a draw costs one pass over its pairs.
- *Categorical moments* are sums of the supplied probabilities, intersections included. Under (0.9, 0.05, 0.05) the subsets {a} and {a, b} and their intersection have moments 0.9, 0.95 and 0.9, where R gives ⅓, ⅔ and ⅓ (O-O2).
- *Every hinge degree* is covered.
  - Uniform, normal and empirical measures give interval moments with no quadrature error. Student components use Gauss–Legendre at every ν (Distributions; see [Probability and Information](../ProbabilityAndInformation.md)).
  - An analysis reports its largest quadrature error estimate (`MomentError`).
  - R's Student closed form needs a 2F1 evaluator, covers degree 1 only and is NaN at ν = 1.

**Order, budget and zero variance:**
- *Main indices* are reported for every subset up to `MaxOrder`. By default that is the posterior's largest term, counted in variables over all parts, so the mains sum to 1.
- *Only subsets inside some term are enumerated;* every other subset's main is exactly 0. More than `SubsetBudget` subsets (10,000) refuses the call before anything is computed.
- *Zero variance.* A draw whose total variance is zero (at most the rounding tolerance times its moment magnitudes, such as an intercept-only draw) has status `ZeroVariance` and NaN indices. It is left out of the means rather than counted as 0.
- *Rounding.* A negative contribution above −tolerance × total is clamped, and below it flagged.
- *Posterior summaries.* The indices' posterior means are weighted over the defined draws, and stage 9's `SobolSummary` gives their quantiles through `DrawSample`.

**Effect functions** (`BassFit.Effects`):
- *What they are:* the main and two-way effects E[f | x_v] − f₀ and E[f | x_v, x_w] − f_v − f_w − f₀, under the declared measure at the fit's degree.
- *Grids:* `GridPoints` (100) values over the measure's support, or the levels of a categorical input. By default every main effect, and every pair a posterior term holds (R's `combs`).
- *A cell budget* of 10⁷.
- *R's effects* use the U(0, 1) moments at degree 1, whatever the measure (O-O2).

**Effect functions for basis fits** (stage 9; `BassBasisFit.Effects(spec, effects)` → `BasisEffectAnalysis`):
- *An effect is linear in the surface,* so it has no cross-component terms.
  - At ℓ, the effect on u is Σ_k g_ℓk·e_k,u, with the mean y_m(ℓ) + g_ℓᵀE[w].
  - In the integrated form, the input effects are weighted by ḡ; the location's main effect is y_m(ℓ) + g_ℓᵀE[w] − f₀; and each one-input effect with the location is Σ_k (g_ℓk − ḡ_k)·e_k,v.
- *Defaults and refusals:*
  - the default set is every main effect, and each pair some component's terms hold;
  - output coordinates and output measures are refused, as basis sensitivity refuses them;
  - mixed hinge degrees are allowed, since no cross terms arise.
- *Checks:* one component reproduces its direct effects scaled by y_s(ℓ)·K_ℓ1, to 10⁻¹³; tensor quadrature of the reconstructed output agrees, to 1.2·10⁻¹⁴.

**Basis sensitivity** (O7; `BassBasisFit.Sensitivity(spec)` → `BasisSobolAnalysis`), R's `sobolBasis`:
- *The surface.* At output location ℓ, the reconstructed latent output is y_m(ℓ) + Σ_k g_ℓk·w_k(x), with g_ℓk = y_s(ℓ)·K_ℓk and w_k(x) = Σ_j a_kj B_kj(x) under a joint draw's paired component draws. The truncation residual and the noise lie outside the decomposition.
- *Independent fits still covary under the input measure.* Their components share inputs, so Var_x at ℓ is Σ_{k,k′} g_ℓk g_ℓk′·Cov(w_k, w_k′), with Cov(w_k, w_k′) = Σ_{j,l} a_kj a_k′l·Cov(B_kj, B_k′l).
  - Each reported subset u of the inputs has its r × r matrix C_u of component covariances, from `BassMoments.Pair` across components as within one, and each input has its total matrix.
  - These matrices are free of ℓ, so each draw keeps them once (`ComponentCovariance`, `ComponentMainCovariance`, `ComponentTotalCovariance`, `ComponentMeans`), and every view is a weighting of them.
  - In the tests the cross terms carry more than 10⁻³ of some locations' variance.
- *Two views,* each a `SobolDecomposition` in output units.
  - `AtLocation(ℓ)` (`PerLocation()` for all) is R's form. The variance is gᵀC g, and the mains and totals are weighted alike.
  - `OtherVariance` and `OtherIndex` (2026-11-02; R's `other`) hold each draw's variance outside the reported subsets, the total less their main variances: what `MaxOrder` leaves out, 0 when every order is reported, and a difference within the rounding tolerance reads 0. R's also holds what `nind` drops; the port has a budget that refuses instead (the outputs brief), so here it is the higher orders alone.
  - `Integrated` adds the location as one more input, `"location"`, with the training locations equally weighted. With ḡ the locations' mean of g and Σ_G their covariance:
    - a subset u has main variance ḡᵀC_uḡ, and u with the location has tr(Σ_G C_u);
    - the location alone carries the variance over ℓ of the location means y_m(ℓ) + g_ℓᵀE[w], so a strongly varying mean curve dominates it;
    - an input's total is the locations' mean of its totals.
  - The law of total variance links the two views: the integrated variance is the locations' mean variance plus the variance of the location means.
- *Units.* Contributions are in output units, y_s applied. R's `sobolBasis` reports them in decomposition units, a factor of y_s² too small: up to 0.96 of the physical variance on O-O2's fixture, whose physical values the port matches to 10⁻¹³. The indices are dimensionless, and R's agree with them.
- *Order, budget and status* follow the direct path.
  - The default order is the largest term over all components, plus the location in the integrated form.
  - `SubsetBudget` bounds the subsets of the inputs.
  - The location means' part is a sum of squared deviations, so it is zero when its sd is at most the tolerance times the means' magnitude, not that magnitude squared. A 10⁻³ variation on offsets of 10⁶ stays defined.
- *Refused:* `Coordinates` (a basis has none), `OutputMeasures` (the locations weigh equally), and components fitted at different hinge degrees.
- *Checks:*
  - a one-component fit reproduces the direct decomposition scaled by y_s(ℓ)²·K_ℓ1², to 10⁻¹³;
  - tensor quadrature of the reconstructed output under a truncated normal mixture agrees to 10⁻¹⁰ of the variance, in both views;
  - the cross terms satisfy Cauchy–Schwarz, and the subsets' matrices sum to the whole;
  - the identities and the law of total variance hold on every draw;
  - thinned and resampled pairings match each component's direct decomposition;
  - mutation runs: dropping the cross terms, y_s or the pairing each fails the tests.

**Checks** (direct):
- O-O3's 156 moment cases under 9 measures, to 6·10⁻¹⁵ relative;
- term moments at degrees 1–3 against a tensor quadrature of the design's columns, to 10⁻¹²;
- the empirical measure against the columns' averages, to 10⁻¹⁴;
- the Sobol identities on every draw, to 10⁻¹²;
- the variance, mains, totals and pair variances within 5 standard errors of pick-freeze Monte Carlo (Saltelli 2010 for the mains, Jansen for the totals) at N = 200,000;
- the law of total variance between the two functional forms, to 10⁻¹¹;
- effects at degree 2 against quadrature, to 10⁻¹¹.

### BASS: calibration

Stage 7 (2026-10-29) builds D09's calibration under `spline/bass/calibration/`. It samples the cut posterior p(e | D_sim)·p(θ | y, e) of the calibration inputs θ, the measurement variances and the discrepancy, given emulator draws e that the observations never update. R's `calibrate.bassBasis` is the reference for the API and the θ kernel; its defects (C-O1) are corrected, none ported.

**The emulator** (`ICalibrationEmulator`, C-X2).
- *What it serves:* retained draws (weighted, when they are), and predictions in physical units at θ in physical coordinates.
- *Each draw's noise:* `ComponentNoise(y_s ∘ K, σ²_e,k)` for a basis fit, `IndependentNoise(σ²_e)`, on the diagonal, otherwise.
- *D07's truncation residuals,* in physical units.
- *Adapters:* `BassBasisEmulator` and `BassFitEmulator` serve BASS's fits. Each predictor is a calibration parameter, matched by name, or a fixed input carried by the observations.
- A non-BASS emulator implements the same contract.

**Parameters.** `ContinuousParameter(Name, Lower, Upper)` is sampled in [0, 1] and reported in physical coordinates through its recorded transform. `CategoricalParameter(Name, Levels)`.

**The measurement model** (`GaussianMeasurement`, C-X1). Variance groups partition the outputs, and each declares its prior:
- `InverseGammaVariance(ν, ŝ)`, R's IG(ν/2, (ν/2)·ŝ²), with ν > 0;
- or `HalfCauchyDeviation(ŝ)`, a half-Cauchy on σ_g (Gelman 2006), which is proper and scale-equivariant.
- *R's −log(1 + σ²) kernel is a defect fixture.* It is improper, and its scale, fixed at 1, makes the posterior depend on y's units. Between data in units 10⁻³ and 10³, σ's posterior median moves from 0.2691 to 0.2386, where the half-Cauchy gives 0.2362 in both.

**The covariance** (C1), declared and validated at construction.
- *Functional* (R's type 1): diag(σ²_{g(j)}) + the emulator's noise + the truncation term.
  - The truncation term is `TruncationTerm.ResidualSet` by default: (1/n)·Σᵢ (y_s ∘ eᵢ)(y_s ∘ eᵢ)ᵀ, the covariance of D07's uniform draw over its residual vectors. It is low rank, and compressed once to its numerical rank.
  - `Diagonal` is R's option: the centred n − 1 variances.
  - Solves and log-determinants go through [`DiagonalPlusLowRank`](../LinearAlgebra.md#diagonal-plus-low-rank-diagonalpluslowrank).
  - The discrepancy enters the mean: δ = μ_δ + D·v, v ~ N(0, s_v²·I), with s_v a required declaration. v is drawn from its conditional on the likelihood path, with precision I/s_v² + b·DᵀΣ⁻¹D.
- *Multivariate* (R's type 2): Σ_trunc + diag(σ_g)·R·diag(σ_g) + Σ_δ + the emulator's noise, on the dense factor.
  - R is an m × m correlation over the outputs, as R's `meas.error.cor` is.
  - Σ_δ and the mean shift are fixed.
  - The truncation term defaults to the residual set's second moment about zero, as in the functional form.

**Categorical parameters** (C2).
- *Below the cap:* the discrete conditional is enumerated whole, up to `EnumerationCap` = 1,024 combinations, with each term's categorical part evaluated once per combination.
- *Above it:* single-site Metropolis steps. On the fixture they accept 0.15, with an integrated autocorrelation of 12.4 scans.
- *The alternative:* a per-site Gibbs draw, which needs no rejection.

**The conditional sampler** (`CalibrationSampler`): a block scan for θ and the rest, given one emulator draw.
- *θ:* R's adaptive random walk in [0, 1], with its 23% target, its scale schedule, its 300-scan switch and its every-100-scans independence decorrelation step.
  - It starts at 0.1² per coordinate rather than R's 10⁻⁶, which hardly moves.
  - It freezes at the end of burn-in, under D10's adaptation policy.
- *Each variance group:* a random-walk step on ln σ²_g with its Jacobian untempered, adapted towards 0.44 acceptance in burn-in.
- *The categorical block, and the discrepancy's Gaussian draw.*
- *Tempering* is an option with no default ladder. It runs through `ParallelTempering` on D05's likelihood path:
  - ℓ_b = b·ℓ₁, with even-odd adjacent swaps whose ratio is (b₁ − b₂)(ℓ(x₂) − ℓ(x₁));
  - a b = 0 rung is refused unless every prior is proper;
  - R, by contrast, tempers its variance Jacobians and swaps off its own target (C-O1: swap gaps of 234 and 140).

**The nested reference sampler** (`NestedCalibrationSampler`, C3).
- *The outer draws:* S emulator draws, every retained draw by default, or by weight with replacement when the draws are weighted.
- *Each inner run:* an ensemble on p(θ | y, e_s), run until N14's verdict reads converged (R̂ 1.01 and ESS 100 per chain by default), within `MaxSamples`. Each inner run's verdict is reported.
- *Pooling:* `Pooled` weights each converged run's draws 1/(S′·n_s), so every emulator draw counts once. It leaves out not-converged and undetermined runs unless asked.
- *The inner-run error.* On a fixture whose exact mixture E[θ] is known, the pooled mean's RMS error is 0.0115 at 25 kept draws per chain, with its spread 15% low. From 100 on it is Monte Carlo noise: 0.0015, 0.0023 and 0.0009 at 100, 400 and 1,600.

**Recovery.** The simulator is η(t; θ) = θ₁ sin(2π(t + θ₂/4)) + θ₂t/2 on 20 locations, run on a 6 × 6 design, with θ* = (1.3, 0.4), noise sd 0.05, and five seeds. The 95% intervals cover θ*:
- 5/5 and 5/5 on a PCA basis fit;
- 4/5 and 5/5 on a direct fit, the less accurate emulator, whose θ₁ means are 1.41–1.47;
- 5/5 with a declared discrepancy, and with two variance groups under each prior.
- At a declared gate of R̂ 1.05, ESS 50 per chain and 1,500 draws, 3 to 6 of the 6 inner runs converged per seed. The rest were reported and not pooled. How often inner runs converge under the default gate is not yet measured.

**Checked** (54 calibration tests, 9 parity tests):
- the adapters' predictions equal the fits' own, bit for bit;
- the low-rank covariance against dense (10⁻¹¹), C-O2's references at its conditionings, and the units (c-scaled outputs move ℓ by −m ln c);
- the target against dense Gaussians with every prior written out (10⁻¹⁰);
- v's conditional against the joint Gaussian at b = 1, 0.35 and 0;
- ln σ² against quadrature;
- the categorical enumeration against brute force;
- the chain's laws against quadrature: untempered, multivariate, and on a four-rung ladder;
- each of C-O1's defects asserted with its correction.

**A cost to watch.** `DiagonalPlusLowRank` factors by a Jacobi thin SVD of the m × k low-rank factor on every σ² proposal. That took 1.4 ms at m = 20 with 39 columns before the residual set's compression, and 33 µs after. Where sizes make it binding, the Cholesky factor of the capacitance matrix, O(m·k² + k³), is the cheaper route (N15).

**Cut SMC** (stage 8's C5, 2026-10-29; [the chip's report](../../../issues/regression/bass/stage8-cut-smc-brief.md#reports)). It is an option beside the nested reference: `CalibrationConfig.Smc`, a `CutSmcSpec`, with null keeping C3's sampler. It follows Mathews, Gopalan, Gattiker, Smith and Francom (arXiv 2406.07555).
- *The first leg* tempers from an exact prior draw to the first emulator draw's conditional μ_1, on D05's likelihood path.
- *The tour* orders the S outer draws, nearest neighbour then 2-opt.
  - The distance d(s, s′) is the weighted spread of ℓ(θ; e_s) − ℓ(θ; e_s′) over the first leg's particles: the quantity a bridge's conditional ESS falls with.
  - Any order gives a valid estimator; the tour changes only the cost. `Order.Random` is the comparison.
- *The bridges* (`CalibrationBridge`): log π + ℓ_{s−1} + t·Δ, with Δ = ℓ_s − ℓ_{s−1}, a target on the likelihood path with t placed by conditional ESS.
  - Each bridge's first step is bisected on cached increments; the core's tempering then takes it to t = 1.
  - One `SmcRun` spans the tour, so its families count the lineages over every bridge.
- *The kernel* (`CutSmcKernel`) is the conditional sampler's block scan, aimed at the bridge:
  - the θ walk and the variance steps are Metropolis on the bridged density, with scales set per stage from the population;
  - the categorical block enumerates the bridged conditional;
  - the discrepancy draws from its bridged Gaussian, with precision I/s_v² + Σ_k w_k DᵀΣ_k⁻¹D.
  - R's independence decorrelation step is left out, and every block keeps its stage invariant.
- *The cut posterior* pools each visit's population at weight 1/S, whatever its evidence.
  - A draw that misfits, with evidence 5 nats below the others', still pools at 1/3, where a joint posterior would weigh it at 0.36%.
  - Weighted emulator draws are resampled by weight into the outer draws, as C3's are, so a draw resampled k times carries k/S.
- *Checked:*
  - the bridge's density, each block's invariance at t = 0.3 and 0.7 by quadrature, and the bridged Gaussian and enumeration against closed forms;
  - the pooled posterior against the exact mixture (E[θ] 0.5359 against 0.5340, 1.5 SE), and each visit's evidence against quadrature;
  - against C3 on stage 7's fixtures: every pooled moment within 2 SE;
  - the cut kept, S = 1, and a collapsing bridge (45 stages, 196 families left, its population on its conditional).
- *Compared* (exploration, on stage 7's recovery fixture):
  - Against a random order, the tour halves the bridge stages (53 against 118 at S = 20) and keeps more lineages (341 against 241).
  - At S = 6, C3 under the default gate reaches cut SMC's error on E[θ] for less processor time. SMC's first leg is most of its cost, about 70 of its 85–94 stages.
  - Each further draw costs about 2.8 bridge stages, so SMC's advantage would come at larger S, which is untested.
  - θ*'s coverage matches C3's.
- *Two readings for stage 7:*
  - every C3 inner run converged under the default gate (R̂ 1.01, ESS 100 per chain), in at most 2,750 draws a chain;
  - on the direct fit, C3 at 1,000 fixed draws a chain sat 2.8 SE from its own 8,000-draw value: a short-run bias that its R̂ and ESS gate does not see.
- *The cost* is the covariance factor's. A bridge sweep builds about four factors, at 164 µs each, against 1.6 µs for a prediction (N15).

### BASS: presentation

Stage 9 (2026-10-29) builds the data R's `plot`, `print` and `summary` methods show, as values in `Maths.Regression`, under `spline/bass/presentation/`. Rendering and text formatting belong above maths (the layer order). No consumer above maths takes BASS yet, and viz v2's contract shapes are open, so neither is built.

| R's method and purpose | The data | The view |
|---|---|---|
| `plot.bass`: traces of M and σ² | `BassTraces`: one `BassChainTrace` per chain, with `Iterations` (each draw's sweep after burn-in), `Count`, `Sigma2`, `Tau`, `Lambda` | `BassFit.Traces()` |
| `plot.bass`: fitted against observed, intervals, residuals | `BassTrainingFit`: `Observed`, `Mean`, `Quantiles` (probability × case × location), `Residuals`, `NoiseSd`, and the method behind its `Intervals` | `BassFit.TrainingFit(data, probabilities)`, through `Predict` at the training cases |
| `plot.bassBasis`: traces per component; fit and residuals in physical units | the same records: a `BassTraces` per component, and a `BassTrainingFit` whose observed values are reconstructed exactly from the decomposition | `BassBasisFit.Traces()`, `TrainingFit(probabilities, seed)` |
| `print.bass`, `print.bassBasis` | `BassMetadata`: `BassDataShape`, the fit's `Config`, `Resolved` (with its auto-picks), `Identity`, `Diagnostics` and `Evidence`, and `BassRunLengths`. `BassBasisMetadata` adds the PCA account, the truncation summary, the joint draws and each component's metadata | `Metadata()` on both fits |
| `summary.bass`, `summary.bassBasis` | `BassSummary`, `BassBasisSummary`: the metadata and `BassPosteriorSummary` (M, σ², τ when sampled, λ) | `Summary(probabilities)` on both fits |
| `plot.bassSob`: boxplots of S and T; functional panels | `SobolSummary` (each subset's main index and each input's total, the index and variance outside the reported subsets (`OtherIndex`, `OtherVariance`), the variances, the draw counts), `SobolAnalysisSummary` (one per decomposition), `BasisSobolSummary` (per location and integrated) | `BassFit.SensitivitySummary(data, spec, probabilities)`, `BassBasisFit.SensitivitySummary(spec, probabilities)` |

- *Every scalar summary* is a `PosteriorSummary(Mean, Sd, Min, Max, Probabilities, Quantiles, Draws)`, taken under the weights through `DrawSample`.
  - `Min` and `Max` are R's boxplot whiskers at `range = 0`.
  - The default probabilities are 0.025, 0.25, 0.5, 0.75 and 0.975.
- *One weighted quantile.* `DrawSample.Quantile` (see [Samplers](../Samplers.md)) serves every summary and a basis fit's predictive quantiles. Stage 5's step-CDF form (R's type 1) is gone, and a direct fit keeps the Gaussian mixture's exact quantile.
  - It is Akinshin's weighted type 7 (2023, §6; stage 8, `d0a0fd8`). It is exactly type 7 under equal weights at any scale, and continuous in the weights, a weight going to 0 included, so a negligible SMC particle moves no quantile.
  - Under unequal weights it is not anchored at the extremes: Q(0) and Q(1) lie inside them once the end draws' weights fall below 1/n*. No rule can be both anchored and continuous, since the maximum over the positive weights jumps when the top draw's weight reaches 0.
  - The rule it replaced, the midpoint rule, jumped between equal and near-equal weights, and let a weight of 10⁻⁹ act as a full node. An atom at 100 among atoms at 0 and 1 moved Q(0.99) from 0.99 to 97.03.
- *An SMC fit* has no traces, since its draws are a population with no order; its stage readings are `Diagnostics.Smc`. Its run lengths report its `SmcSpec`, with the sweep counts null.
- *Defaults* (the stage 9 chip's report; confirmed by the user on 2026-10-29):
  - indices are summarized over the defined draws, variances over every draw;
  - `NoiseSd` is E[σ], as R's code computes it;
  - the trace axis counts sweeps after burn-in, a ladder's scan counting as one;
  - a basis fit's training-fit intervals come from one predictive draw per joint draw, and its observed values from the decomposition, since the fit does not store Y.
- *Intervals* are the predictive mixture's (D07). On O-O1's construction the training fit gives ±2.771 and ±2.644, against R's per-draw endpoints of ±3.91 and ±5.20.
- *What an adapter will need.* Every record is arrays and values with explicit axes and names: trace panels from `Iterations`; scatter and interval segments from `Observed`, `Mean` and two `Quantiles` rows; residual bins with N(0, `NoiseSd`²); boxes from `Min`, the quantiles and `Max`; effect bands as a reduction of each cell's draws. `BassConfig.ToJson()` stands in for R's printed call. A basis fit holds no output coordinates, so the caller supplies its grid.

### BASS: the extension seams

Each of the port's segments listed room for extensions. These are the contracts for them, built on 2026-10-31 ([the report](../../../issues/regression/bass/extensions-preparation-brief.md#reports), with every shape's alternatives). Today's code runs through them, and no extension fills them yet. Each was a refactor that changes no draw: the eight bit-identity hashes (six BASS fits, Friedman tempered, the tempered BARS fit) held at every seam.

| Seam | Contract | Today | Filled by |
|---|---|---|---|
| Term types (F-X3) | `IDesignKernel<TPoints>` | `DesignKernels`: the continuous hinge product, the categorical subsets, the functional hinge product | ABMRS's tensor B-splines; a directional hinge |
| Term derivatives (O-X4) | `BassBasis.CaseDerivative`, `OutputDerivative`; `BassMoments.FactorDerivative` | Hinges differentiate and a categorical part refuses; nothing reads them yet | Active subspaces |
| Observation model (F-X1) | `IBassObservationModel.Working` → `WorkingResponse`, and `Law` | `GaussianObservation`: zero offset, unit weights; `VarianceMeanMixture`: Student t ([below](../VarianceMeanMixtures.md)) | GBMARS |
| Sweep slots (F-X2) | `BassSweepSlots`, `IBassObservationModel.Blocks(slot)` | The slots after τ are empty | GBMARS |
| Output-side inner product (F-X6) | `IInnerProduct`, on `WorkingResponse.Output` | `IdentityInnerProduct`, Q = I | Residuals correlated along the output index |
| Face-splitting contraction (N-X1) | `KhatriRao.FaceSplittingGram` and its kin | Only tests call it | ABMRS |
| Wide-band factor (N-X2) | `WideBandCholesky`, `SpdFactors.ForBand` | Nothing selects it; BASS keeps the dense factor | ABMRS |
| GIG limits (N-X3) | `GigParameters`, the GIG forms of `ConjugateUpdates` | The two updates keep their formulas, as the forms' limits | GBMARS |

The extensions themselves are in [the extensions brief](../../../issues/regression/bass/extensions-brief.md). E34's several outputs on one basis, Σ integrated into a K × K block with the design's K-output side beside its scalar statistics, are on their own page, [Several outputs](../MultiOutput.md). F-X4, the proposal-kernel contract, is not among these: it is built with KHAOS's Bernoulli-inclusion proposal, not before it (decided 2026-10-29), and was so on 2026-11-07 ([Polynomial chaos](../PolynomialChaos.md)). The shared pieces are documented where they live: [the contractions and the inner product](../LinearAlgebra.md#khatrirao-contractions-khatrirao), [the wide-band factor](../LinearAlgebra.md#wide-band-cholesky-widebandcholesky) and [the GIG limits](../Samplers.md#block-scan-and-conjugate-primitives-oq18-2026-10-24).

**Term types (F-X3).** A term is a product of parts. Each part's type sits behind `IDesignKernel<TPoints>` on its side of the Khatri–Rao column, `CaseInputs` or `OutputInputs`, and supplies:
- its column evaluation (`Multiply`, applied in order to a column of ones);
- its derivative in one variable (`MultiplyDerivative`, for a `Differentiable` part);
- its normalizing constant κ (`Scale`) and its separable `Factors`, each a one-variable `DesignFactor` with its value, derivative and moments.

`BassBasis` composes columns from `DesignKernels` and `BassMoments` composes moments from them, so sensitivity serves every type by one rule: a term is κ·Πᵥ φᵥ(xᵥ).
- *Per part and generic in the points type, with `BassTerm` closed.* A new part kind adds a typed list to `BassTerm` and a member to `DesignKernels`, and the contract stays as it is. The alternatives were one non-generic interface with a side enum and both inputs, and `BassTerm` holding polymorphic parts, which touches every reader of `Continuous`, `Categorical` and `Functional`.
- *Support stays with the composed column* (`BassBasis.Support`). The admission gate reads the support of the product, hinge times indicator, so a per-kernel count would need the joint column anyway.
- *Moments are separable factors.* `HingeDesignFactor` and `SubsetDesignFactor` are value records, so `BassMoments` caches a factor's moments, and a pair's, by the factor itself. Kernel-level `Mean` and `Cross` would move the ANOVA bookkeeping into every kernel.
- *`DesignKernels` is a fixed set.* A registration hook is the extension's to add.
- *How the members join.* A tensor B-spline part is a case-side kernel whose column is the product of one-variable B-splines, whose derivative differentiates one factor, and whose factors are piecewise polynomials that `PiecewisePolynomial` integrates as it does a hinge. A directional hinge, (wᵀx − t)₊, is differentiable in every variable of w and is no product of one-variable factors, so it declares none and supplies its moments through a further interface that the moment service tests for. A test implements this kernel against today's contract, which shows the signatures carry it.

**Term derivatives (O-X4).** `BassBasis.CaseDerivative` and `OutputDerivative` write a term's derivative in one variable over the points, and `BassMoments.FactorDerivative` sits beside `FactorValue` for one factor. Nothing in sensitivity reads them yet; the active-subspace work writes its gradient integrals from them.
- *At a kink, the mean of the one-sided limits:* s/2 for a degree-one hinge at its knot, 0 from degree two up. That is what a central difference returns, and it matters because knots are training inputs, so cases sit on them. The alternatives were the right derivative and the Clarke interval.
- *A variable belongs to one part,* so a derivative is the differentiable part's derivative times the other parts' values, with no general product rule.
- *A categorical part refuses* (`NotSupportedException`). Zero would claim a derivative in a level index.

**The observation model (F-X1).** `IBassObservationModel.Working(response)` hands the Gaussian core a `WorkingResponse`: `Values`, y − offset (the response itself, uncopied, when there is no offset), `Weights`, and `Output` (below). `BassDesign` and `BassResolved.Resolve` read it, so an offset is the shifted data in the resolved record, the target and the design, bit for bit.
- *The core holds weights* since 2026-11-02 ([the weighted core](#bass-the-weighted-core)); until then it refused them.
- *A runtime argument to the model:* `BassModel.Create(config, data, observation, form)`. A `BassConfig` field would change the serialized config and the fingerprint. The implication is that the seam stops at the model: `Bass.Fit`, `Bass.Resume` and `Bass.Fingerprint` take no observation model, so a run through them is Gaussian, and the extension that fills the seam adds it to the run and to the fingerprint. E1 did so on 2026-11-05, as a declarative `BassConfig.Noise` that `BassModel.Create` reads ([Student t noise](../VarianceMeanMixtures.md)).
- *What reads the offset.* The pointwise densities, the training fit and the residual-dependence check read each draw's working response ([a design whose working response moves](#bass-a-design-whose-working-response-moves)). Prediction adds no offset back, which is the extension's.
- *How GBMARS fills it:* y − β√w·v and weights 1/(w·v) at the current latent state, on the weighted core.

**The sweep's slots (F-X2).** `BassSweepSlot` is D04's general order: `Structure`, `Tau`, `Coefficients`, `Latent`, `Scale`. `BassSweepSlots.Blocks` is the one place it is written. A chain starts at the structure's slot; a tempering replica or an SMC particle, whose moves the engine runs itself, starts at τ's. The observation model's `Blocks(slot)` fills the later three. Gaussian BASS returns none, so its scan has the two blocks of [the sweep](#bass-the-sweep), and an empty slot draws nothing and consumes no random number.
- *The latent state lives with the observation model,* whose blocks are the existing `IBlock<ReversibleJumpChain<BassState>>`. A latent field on `BassState` would change every snapshot, so it waits for an extension that needs it.
- *How GBMARS fills it:* a block for a, one for the latent scales, one for σ² and the skew.

**The output-side inner product (F-X6).** `IInnerProduct` is Q on the output index, for a residual with covariance σ²(I ⊗ R), Q = R⁻¹. `WorkingResponse.Output` carries it, and `BassDesign` routes ‖ỹ‖², the residual's sum of squares, the output-side Gram entries of a border and the cross products through it, so a declared identity takes today's arithmetic exactly. A scalar response refuses a Q.
- *The centred forms stay at Q = I,* and `BassDesign` refuses another Q with the reason. Under Q the intercept's orthogonal complement is the Q-centred one, so that centring, the level and σ²'s prior scale are the extension's.
- *Q rides the working response,* with the offset and the weights, as the noise's part of what the core fits. The alternatives were a `BassConfig` field, which changes the serialized form, and a `BassData` member. So like the offset it reaches the model, not a run.
- *How correlated residuals fill it:* a tridiagonal `IInnerProduct` for the Markov residual, with the Q-centred forms.

**What each extension still writes.**

| Extension | Beyond filling its contracts |
|---|---|
| GBMARS | Built for Student t, the asymmetric Laplace, the Normal-Wald, the skewed t, the generalized hyperbolic and the generalized beta prime with the horseshoe, on chains, tempered ladders and SMC with the evidence, on a direct functional response, and on BARS, skewed and quantile BARS among them, with the noise's scale interwoven under the z parameterization, opt-in, and a move along β·v, on by default under the skewed horseshoe alone at the user's "Default there" of 2026-11-08 (2026-11-05 to 08; [their page](../VarianceMeanMixtures.md)): the latent record on the working response with σ², β and γ on it, the coefficients', latent, scale and noise-shape blocks, the GIG sampler, the modified half-normal, log Bessel K and the generalized hyperbolic, the noise's law in pointwise, prediction and the training fit, the noise on the config, the run and the fingerprint, and retained draws that keep their latents alone. Still to write: the flat histogram under a noise, refused, with two routes on the page |
| KHAOS | Slices 1, 2, 4 and 5 built (2026-11-07 and 08; [its page](../PolynomialChaos.md)): the Legendre part on `BassTerm` and among the case-side kernels, its term space and measure, F-X4 with KHAOS's proposal, the Poisson-binomial, declared input domains, the prior-precision builder with KHAOS's raw ridge and modified g-prior, g₀² by τ's slice step under a Gamma on τ, the ladders, SMC, the flat histogram, the pointwise and candidate densities and GBASS's noise gated on Legendre terms; Sobol and the Jacobi derivative measure from the coefficients under the uniform law, the other measures' moments gated, and the gradient readout E8 takes. SMC's trap on a raw-ridge fixture is flagged, the coordinating session's call. Still to write: ridge blocks (slice 3) |
| ABMRS | A `BassTerm` slot or its own term type; the knot moves; band output from the face-splitting accumulation, with a `DecomposeBanded` on the wide-band factor; the factor's use by its target |
| Correlated residuals | The Q-centred forms; the Q-aware level and σ² prior scale; the tridiagonal `IInnerProduct` and its θ updates |
| Active subspaces | The gradient integrals from `FactorDerivative`; the moments of a non-separable kernel |

**Checked.** The shared pieces' checks are with them. Here:
- *The term types* (`BassDesignKernelTests`). The factors with the scale reproduce the column to 10⁻¹³ relative at degrees 1–3 on the scalar and functional fixtures' terms, and to 10⁻¹² on every term of a direct functional fit. Derivatives match central differences: to 10⁻⁶ at interior points, and to 5·10⁻⁵ at the training inputs, where a degree-one derivative jumps, on both sides of the column. The directional hinge matches to six digits.
- *The observation model* (`BassObservationModelTests`). A design on y − offset equals a design on the shifted data exactly, scalar and direct: the births' proposals, the sums of squares, the coefficients and the residual. So do the resolved record's level, mean square and σ² prior, and two 150-step chains, state by state.
- *The slots* (`BassSweepSlotTests`). A model that fills every slot with nothing walks Gaussian's chain, streams and block counters word for word over 100 steps. A filled model runs its blocks after τ in slot order, at every rung of a three-rung ladder. The SMC kernel's reading rests on the shared assembly and its Gaussian parity tests, not on a filled model.
- *The inner product in the design* (`BassOutputProductTests`). A declared identity equals the default direct design exactly, through four births and a rebuild.
- *Not checked,* since the core refuses it: a fit under a Q other than the identity. Weights are checked [below](#bass-the-weighted-core).

**The shapes are open to review.** They went to the user on 2026-10-31 ([the fitting checklist's report](../../../issues/regression/bass/port-fitting-checklist.md#reports)), and no ruling is recorded. Only the tests and today's members lean on them, so revising one costs little.

### BASS: the weighted core

Built 2026-11-02 ([the report](../../../issues/regression/bass/extensions-checklist.md#reports)), ahead of the consumers that need it: E1, E20, E22 through E20, E31's multiplier and E32. The model is y ~ N(α + Z_cβ, σ²W⁻¹), with W from `WorkingResponse.Weights`. Unit weights take today's arithmetic: the hashes hold. Nothing fills the seam yet, and `Bass.Fit` stays Gaussian.

**The rows are whitened.** `BassDesign` holds each column as W^½(z − z̄_w), centred under W, and the response as W^½(y − ȳ_w). Every inner product already in the design then reads the weighted form, so the Gram G_c = Z_cᵀWZ_c, the cross products, the residual ‖r‖²_W, bordering, rebuilds and the admission gate carry the weights with no branch of their own. The intercept's Gram entry is Σwᵢ, its column W-orthogonal to the others. A draw uncentres at z̄_w.
- *The alternative* was an inner product carrying W through each dot, on the `IInnerProduct` pattern. It touches every dot, including the unit path's, for the same numbers.
- *The implication:* `Residual` reads W^½r, not r.
- *The gate reads W,* as this upgrade asks: two columns that differ only where the weights nearly vanish are one column to the fit. So F depends on W, and a latent move could put the current structure below s_min. A latent block treats that as off the support ([the options, and BARS's other reading](#bass-a-design-whose-working-response-moves)).
- *A direct functional design holds its stacked columns,* since per-row weights do not contract through the Khatri–Rao pair. That is O(N) per inner product and N·M doubles, against the unit path's contraction. Weights separable as case × location would keep the contraction. They are not built, since the seam's contract is per row.
- *The support threshold* reads the raw parts: rows, not weights, so a latent scale does not move it. Replication weights count on a scalar response (below).

**What a weight states** (`WeightMeaning`). The two readings share the statistics and differ in what the likelihood counts.

| | `Precision`, the default | `Replication` |
|---|---|---|
| Likelihood | Πᵢ N(yᵢ \| fᵢ, σ²/wᵢ) | Πᵢ N(yᵢ \| fᵢ, σ²)^wᵢ |
| σ²'s ν | Nb | (Σwᵢ)b |
| Normalizer | +(b/2)Σ ln wᵢ | None |
| Support | Rows | Σwᵢ over a scalar response's nonzero rows |
| τ₀ under the weighted prior, the auto hot end | 1/N | 1/Σwᵢ |
| Who hands it | A latent scale, Pólya–Gamma, a mixture | Aggregated rows, resampling (E18) |

Only `Replication` makes integer weights equal duplicated rows when σ² is integrated out. Under `Precision`, given σ², the two differ by −(b/2)[(Σwᵢ − N)ln(2πσ²) + Σ ln wᵢ], free of the structure and τ, so their Bayes factors agree; σ²'s integral separates them, since it counts a different ν. Precision alone would have left that check true only given σ².

**Which columns the coefficient prior reads** (`WeightedPriorGram`).

| | `Unweighted`, the default (D03's robust row) | `Weighted` |
|---|---|---|
| Prior precision | τ·XᵀX on the columns as they stand; the ridge's D from them | τ·XᵀWX; the ridge's D from them, Σwᵢ for N |
| Level, s² | Unweighted | Under W |
| A latent scale's conditional | Closed form: the prior does not read W | Gains \|ZᵀWZ\|^½ and a quadratic, so needs a Metropolis correction |
| At unit weights | Gaussian BASS | Gaussian BASS |
| Duplicated rows | No: the prior stays on the N rows | Yes, under `Replication` |
| Cost of each (τ, b) | O(p) after O(p³) per structure | O(1) (g-form), O(M) (ridge), as unit |

A third reading, `Frozen`, is a prior on a fixed W₀ ([priors that stop reading a moving response](#bass-priors-that-stop-reading-a-moving-response)).

So "integer weights equal duplicated rows" holds for the likelihood's statistics under every reading, and for the marginal only on `Replication` with `Weighted`. D03's default cannot satisfy it, by construction, since a duplicated row's own g-prior reads its copies.

`Weighted` needs nothing new: the unit block's closed forms hold on the weighted statistics, with `GPriorWeights` carrying the count, the mass and Σ ln wᵢ. `Unweighted` is the pencil.
- *The pencil* (`PriorPencil`, `CollapsedPencilPrior`). With prior precision τD and likelihood Gram A, D = LLᵀ and L⁻¹AL⁻ᵀ = VΛVᵀ give bA + τD = LV(bΛ + τ)VᵀLᵀ. So the determinant term is ½Σⱼ ln(τ/(bλⱼ + τ)) and S = b(‖r‖²_W + Σⱼ qⱼτ/(bλⱼ + τ)), the ridge's `RidgeTerms` over other eigenvalues, with their slopes in b. At unit weights every λⱼ = 1, the g-form's τ/(τ + b).
- *The coordinates* stay (α, β) with columns centred under W, the design's own. The prior reads α on the columns' own centring, α_u = α + δᵀβ with δ = z̄ − z̄_w, so D = D_α(1, δ)(1, δ)ᵀ + diag(0, D_β), and under weights α and β are drawn jointly. A flat intercept has no prior to couple and stands apart, with its integral (2πσ²/(Σwᵢ·b))^½.
- *The prior's Gram clears the gate too,* since it is a precision.
- *Why not `FromWeighted`,* BARS's weighted block. It is the same prior, but built at one τ and b = 1 with no slope, so BASS's τ step and tempered swaps could not read it. The pencil gives any (τ, b), and BARS's weighted tempering could read it too, on a design that holds both Grams.
- *Each admitted proposal costs its own O(p³),* with no bordering on the parent, since dropping a column of a full D re-triangulates L. A latent move already costs O(NM²) per sweep, so this is the same order for the models that want this prior. The alternatives are a bordered pencil, a Cholesky of bA + τD per (τ, b) (a third the flops per evaluation, but τ's slice step needs several and each loses the slopes), and D04's τ given a. Each is a lever for E1 once it measures its own runs.

**The marginal given σ².** ln p_b(y | structure, τ, σ²) = −(ν/2)ln(2πσ²) + T − S/(2σ²), from the same S and determinant term T as the collapsed form, which is its integral against σ²'s prior (`GPriorMarginal.LogMarginalGiven`, `LogMarginalGivenSlope`). It is finite on an exact fit and exactly 0 at b = 0, where the flat intercept's pole has already refused. So F has no exact-fit refusal, the values are complete (`HasCompleteValues`), and a ladder's evidence reads them.

It is selected two ways, and no default changes:
- *`Sigma2Prior.Held(v)`, a config value:* σ² is a point mass, as probit's identification or a known noise states it. It is the shape τ's `TauPrior.Held` already has, and a `Bass.Fit` can run it. `Value` is left out of the JSON when null, so no other kind's form or fingerprint moves. It lifts the refusals of a constant or zero response, which need σ²'s posterior.
- *`BassTarget.Given(σ²)`, at runtime:* the target given a σ² that moves, handed to the chain by `ReplaceTarget`, the pattern of BARS's robust sweep. This is the form for E1's σ² block and E32's outside σ². `Sigma2` stays the prior such a block draws against.
- *The alternative* was σ² on `BassState`, which changes every snapshot, the decision the seams left to E1. A target per draw costs one small allocation.
- *The chain's own target.* `BassTauBlock` and the composition read the target the chain holds now (`ReversibleJumpChain.Target`), a lone chain's, a ladder's b = 1 replica's or a run's, so a sampled τ sits beside a σ² that moves.

The composition draws β given σ² (`GPriorConditional.SampleGiven`), and the collapsed pointwise content is a Gaussian predictive (`LeaveOutLogDensityGiven`).

**What reads the design.** `Compose`, the conditional means and `CentredDesign` centre under W (`WorkingResponse.ColumnMean`). The pointwise densities read the working response, N(y | f, σ²/w) under `Precision` and w copies of N(y | f, σ²) under `Replication`, whose collapsed predictive refuses, since a row left out is all its copies at once. Prediction's noise reads σ² as a unit-weight observation's. `BassResolved.Resolve` reads the weights once, and the level and s² stay as read, or as a [frozen reference](#bass-priors-that-stop-reading-a-moving-response) states them.

**Checked** (`BassWeightedCoreTests`, `BassBitIdentityTests`).
- *Unit weights.* The eight hashes hold. Declared unit weights under all four readings equal Gaussian to 10⁻¹⁰, through the marginal and its slope.
- *Integer weights and duplicated rows.* The design's statistics, proposals and rebuilds agree to 10⁻¹⁰ under all four readings. On `Replication` with `Weighted`, ℓ_b and its slope agree at each of two forms, three intercepts and four σ² priors. A chain on the weights walks a chain on the duplicated rows state by state for 400 sweeps on one seed, draws included, on scattered inputs. On a grid some terms span the same columns, so a move's ratio is an exact tie, and the sign of its rounding decides whether the chain draws a uniform: the streams then part with no difference in the target.
- *Dense algebra.* S, the determinant term, the conditional mean and the fit's variances, for each form, intercept and prior, against matrices built from the columns. Per-row weights on a direct design. The retained density row by row, and the collapsed predictive against the joint Gaussian by conditioning.
- *σ²'s integral.* The marginal given σ², integrated over σ²'s prior by a trapezoid on ln σ², reproduces the collapsed marginal to 10⁻⁹, under a scaled, a declared and the Jeffreys prior, at b = 1 and 0.4. Its slope matches a finite difference.
- *Enumeration.* All 2,263 structures of F, for a weighted target under either prior and for targets given σ²: detailed balance and stationarity to 10⁻¹³, the log posterior against dense algebra to 10⁻⁹, and a chain's visits against π. With σ² drawn apart by a block and the target replaced each sweep, a chain's visits match the collapsed posterior.
- *The rest.* The gate refuses columns collinear under W that it admits unweighted. Snapshots and hand-overs continue bit for bit, a ladder's replicas included. A held σ²'s run retains it and reports an evidence.
- *Cost, a rough reading under load.* The unit path's draws are bit-identical, so block k does the same work in either binary, and blocks pair across clean builds of the commits before and after. Eight alternating rounds of nine blocks of 20,000 sweeps, pinned to the P-core threads, with each block's load recorded (`BassWeightsCostTests`). No block met the statistical brief's quiet rule: one other test host throughout, and other processes kept 1.6–14 cores busy, median 3.
  - *The unit path,* wall-clock after/before over 72 pairs: Friedman at M = 9–10, a median of 0.99 (IQR 0.94–1.05), about 10 µs per sweep. Direct functional at M = 5, 1.07 (IQR 0.98–1.26), about 7 µs, which one burst of load inflates: 0.99 on the pairs where others kept at most 4 cores busy. A difference under about 5% cannot be read here.
  - *A slowdown this found and removed.* An earlier pass read Friedman at +6 to +14%. The support count had taken a nullable lookup and a floating-point add per nonzero element, on every proposal. The unit path has its own loop again.
  - *The readings, on Friedman:* the prior on the weighted columns about 13 µs at M = 9–12, the unit path's order; the pencil about 44 µs at M = 8–10, so 3–4×, growing as M³. These ran under more load (others' median 7–9 cores), so read them against each other, not against the unit path's figures.
- *Not checked.* Weights with a Q other than the identity, which the core still refuses. A run through `Bass.Fit` under weights is checked since E1's Student t ([below](../VarianceMeanMixtures.md)).

### BASS: priors that stop reading a moving response

Built 2026-11-02 ([the options](../../../issues/regression/bass/moving-response-priors-options.md), [the report](../../../issues/regression/bass/extensions-checklist.md#reports)), steps 1 and 2 of the note's order, as a family of selectable forms on one seam rather than one pick. D03 stays Gaussian BASS's default, and no default changes. Which route a consumer takes is the user's call, per engine.

**The problem.** A moving working response is not data. y − β√w·v is a function of β, w and v, so a prior that reads it is conditioned on parameters. Two things are coherent: read something that does not move, or carry the conditional into every block. An uncorrected re-read does neither, since its conditionals share no joint.

**The constraint.** Noise with a location-scale group, y → ky + m, covers Gaussian, t, Laplace, asymmetric Laplace and Normal-Wald. Under that group a prior can have at most two of proper, reading no data, and equivariant, because a proper measure cannot be invariant under a non-compact group. Binary responses have no such group, so a declared prior costs them nothing. Each route gives one up:

| Route | Selected by | Proper | Reads | Equivariant |
|---|---|---|---|---|
| Invariant | `GPriorCentre.Flat`, `Sigma2Prior.Jeffreys` | No | Nothing | Yes |
| Frozen read | `BassPriorReference` under `Mean`, `ScaledToResponse` | Yes | Observed y, once | Yes |
| Declared | `Zero` or a stated level, `InverseGamma` or `Held`, κ | Yes | Nothing | No |
| Hierarchical | `Sigma2Prior.Hierarchical`, `HalfCauchy` | Yes | s², for the hyperprior | Yes |

The conditional on the latents and the pilot read wait for the latent blocks, and the fractional rung for the run ([the note](../../../issues/regression/bass/moving-response-priors-options.md#what-to-build-in-order)).

**The reference** (`BassPriorReference`, on `WorkingResponse`) holds what the priors read, apart from the working response, so the target stays one target whatever moves later. The observation model names the rule, since the read has to mean something for its noise:

| Rule | Level | s² | For |
|---|---|---|---|
| `Moments` | ȳ | Mean square about it | D03's read |
| `About` | Stated | Mean square about it | The zero centre, at 0 |
| `Robust` | Median | (MAD/Φ⁻¹(¾))², or the mean square about the median when over half the values tie | Heavy tails |
| `Quantile` | The ⌈qN⌉-th smallest | (Mean check loss)² | The asymmetric Laplace's level and scale, by maximum likelihood |
| `FrozenGram` | Mean under W₀ | Mean square under W₀ | A prior on a frozen W₀ |
| `Link` | None | None | κ alone, for a response with no group |

- *What it leaves null the config must not read:* level centring needs a level, and a σ² prior that reads s² needs a mean square. Each is refused otherwise. The zero centre refuses a level other than 0.
- *D03's refusals read it.* A working response that happens to be constant is no reason to refuse once the reference has spread. A reference with none is refused, unless σ² is held. The Jeffreys prior still refuses a constant working response, since its posterior does not exist there. The exact-fit floor reads s² too.
- *The alternatives.* Standardizing y at entry gives the same numbers with every output to map back. Re-reading the moving values has no joint.
- *The implication:* the data are used twice, and an evidence stays empirical-Bayes in what is read, as `BassEvidencePriors` already records for D03's default.

**The level's energy.** Under a frozen level c, ȳ_w moves away from it, and the gap enters as the zero centre's energy does at c = 0: S = b(‖r‖² + τ/(τ + b)·(Q_c + Σwᵢ(ȳ_w − c)²)), with τ/(bΣwᵢ + τ) for the ridge's intercept. With no reference the centre is ȳ_w and no arithmetic is added. The pencil carried this term already.

**The coefficient prior's constant κ.** The covariance is σ²κ(τG)⁻¹, so every form reads τ/κ. For a σ² held at 1, Sabanés Bové and Held's (2011) κ = v(h(0))/h′(0)² anchors the prior at an imaginary sample, so it reads no data: 4 for the logit, π/2 for the probit, π²/4 for the cauchit, 1 for the Poisson log.
- *Why on the reference, not on τ's prior:* the link supplies it, so the observation model states it. Reading τ/κ reaches a held τ and hyper-g as it reaches the unit-information Gamma, where changing τ₀ would reach only the last. τ₀ keeps its meaning, and a draw's τ stays the state's.
- *The implication:* the automatic hot end is 1/(κN), or the held τ's halving point over κ. The flat intercept's pole stays at 2/N.

**The frozen Gram** (`WeightedPriorGram.Frozen`) puts the prior on the information at a fixed point (Li and Clyde 2018), free of the moving W: columns centred under W₀, D_α = Σw₀ and D_β = Z_c₀ᵀW₀Z_c₀, or under the ridge D_α = 1 and D_β = diag(G₀)/Σw₀. It takes the same pencil as `Unweighted`.
- *W₀ reads in information per observation,* so τ₀ stays 1/N, and under the g-form W₀ = (1/κ)·1 is κ on the unweighted columns. That is why κ and W₀ are one idea: κ is the information at β = 0. The alternative, τ₀ = 1/Σw₀, would normalize W₀'s units away and lose the identity.
- *The ridge's normalization strips W₀'s units,* as under `Weighted`, so there only κ carries them.
- *With no likelihood weights* the working response reads unit weights, at the whitened path's cost.

**The hierarchical prior** (`Sigma2Prior.Hierarchical(a, h)`). σ² | r ~ IG(a, r), r ~ Gamma(h, rate h/(a·s²)), whose mean is the scaled prior's a·s² at every h. This one family connects the other routes:
- *h → ∞ is the frozen read:* r pins to a·s².
- *h → 0 is the invariant prior, for any a:* IG(a, r) mixed over 1/r gives 1/σ².
- *a = h = ½ is σ ~ half-Cauchy(s)* (`HalfCauchy()`; Wand et al. 2011).

So s² sets the hyperprior's scale alone, and a poor read hurts less than under the plug-in. The choices made:
- *r is a coordinate of the state, as τ is* (`BassState.LogRate`). The target stays one immutable target, so swaps, SMC's particles, the reference rung and the flat histogram carry r with no change of their own. The alternative, a target per replica through `ReplaceTarget`, would make every swap re-evaluate under the receiver's target.
- *σ² stays integrated, and the marginal is given r,* as it is given τ. A slice step on v = ln r (`BassSigma2RateBlock`) runs last in the sweep's scale slot, at τ's cost per evaluation: O(1) under the g-form, O(M) under the ridge, O(p) on a pencil.
- *The alternative is Wand's conjugate two-block form,* σ² | · then r | σ² ~ Gamma(a + h, rate + 1/σ²). It is exact only on the integrated form's rungs, and it draws σ² every sweep. The slice step is exact on every path and form, so it is the one built.
- *With σ² drawn apart* (E1's w), the marginal does not read r, so the state's prior carries ln IG(σ²; a, r), and the same step reads it. The conjugate draw is exact there, should the slice's cost ever matter.
- *A chain starts at r's mean* unless told otherwise; a run starts from its prior draw of r. r is not monitored: R̂ reads σ², which depends on it.
- *Below the doubles.* Under a tiny h much of r's prior lies below the doubles, about half at h = 10⁻³ and 93% at 10⁻⁴, and all of it counts, since the prior is normalized over all of r. So nothing is cut, and v stands for r:
  - `LogDensity` reads any v, and `SampleLogRate` draws v through `Gamma.SampleLog`.
  - `PriorAt` reads r no lower than the smallest normal double r_min. S, far above either, cannot tell r from r_min, so the integrated marginal is exact once `LogRemainder` adds a·(v − ln r_min). At b = 0 it reads no r.
  - A draw retains r at r_min or above (`RateAt`), and the pointwise densities read it exactly, since rᵃ cancels in them.
  - A chain starts no lower than r_min (`BassState.StartLogRate`), since from far below a step's reach it could not climb within burn-in. A prior rung and SMC's particles keep the draw itself.
  - *The alternative* was to draw again on an underflow and add ln P(r in the doubles) to the evidence. It is quieter while a is not small, since the draws it drops weigh nothing. But once a + h is small r's own posterior reaches below the doubles, and that mass would be lost: about 0.25 nats at a = h = 10⁻³ (derived here, not run). Drawing again with no correction is what this replaced, whose ln Z read high by 0.6 at h = 10⁻³ and 2.6 at 10⁻⁴.
  - *The implication:* of order h/a of the prior's draws land where the likelihood is, so a tiny h costs an evidence precision, and the run reports it. At a = ½, h = 10⁻⁴, against the enumeration's −23.29, a ladder of 8 levels read −23.75 ± 0.36 on 3,000 samples and −23.25 ± 0.13 on 30,000, with its largest pair's rejection at 0.998 or above, so `Saturated`. SMC read −23.47 ± 0.14 on 4,000 particles and −23.30 ± 0.07 on 16,000.

**What was checked** (`BassPriorReferenceTests`).
- *Unasked, nothing changes:* a reference that restates today's read walks today's chain, draw for draw, under the g-form's three centres and the ridge's level and zero centres.
- *Dense algebra.* S, the determinant term and the conditional mean for each form and Gram, and the collapsed predictive against the joint normal, each to 10⁻⁹; the slope in b against a finite difference.
- *The trilemma on the targets.* Under y → ky + m, across snapshots, the moments, robust, quantile, frozen-Gram, invariant and hierarchical routes keep their log posterior differences to 10⁻⁸, and a declared prior's move. Chains on y and on ky + m visit the same structures, sweep by sweep, through a moving response. A moving response here is a model per snapshot, each segment starting where the last ended, and [a chain moved in place](#bass-a-design-whose-working-response-moves) passes the same test.
- *The family against outside readings.* At a = h = ½, against `HalfCauchy`'s own density, integrated on σ, to 10⁻⁸. Across structures both limits hold: at a = 2 the gap between two structures is 1.51143 at h = 10⁴ against the plug-in's 1.51160, and 1.21895 at h = 10⁻⁵ against Jeffreys'.
- *Enumeration.* Detailed balance on every structure of F for a frozen level and a frozen Gram, against dense algebra. Chains' visits against the posterior, including with r integrated by quadrature, with σ² integrated and drawn apart. The rate's own conditional against quadrature, on both paths.
- *The evidence,* against the enumeration's −16.294 with r integrated by quadrature: −16.296 on a ladder, −16.281 on a linear ladder, −16.302 on SMC.
- *The evidence under a tiny h,* each against its own enumeration, with r integrated whole by a closed form below the quadrature, reading no inverse gamma: `ATinyHyperShape_ReportsTheStatedPriorsEvidence`, whose check allows the run's reported error up to a ceiling, and a ladder's flag must read the thin first stone at a = ½. The readings are in [the report](../../../issues/regression/bass/extensions-checklist.md#reports). `BelowTheDoubles_TheTargetReadsTheRate` holds the marginal against its closed form in logs, tempered and given σ².

### BASS: a design whose working response moves

Built 2026-11-02 ([the report](../../../issues/regression/bass/extensions-checklist.md#reports)), ahead of the consumers that need it: E1, E20, E23, E32. Nothing fills it yet, since a latent block is the consumer's and `Bass.Fit` stays Gaussian. No default changes, and Gaussian BASS's draws hold their hashes.

**What moves.** A working response is immutable, so a move is a new one: `WorkingResponse.Moved(offset, weights)` is the model's reading at another latent state, its meaning, prior's Gram, inner product and reference kept. A chain takes it between sweeps, two ways:
- *`BassStructureMoves.ProposeFollow(state, working)`* returns the candidate for a block to weigh. Like a structure move's candidate it is pending until the chain shows it, and a refused one is scored off the support. A Metropolis step reads this.
- *`BassStructureMoves.Follow(chain, working)`,* and `BassChain.Follow`, moves the chain, for a draw from the conditional or a response given from outside. It returns the refusal with everything as it was.

A block reaches its replica's moves by `BassStructureMoves.Of(chain)`, so the same block serves a lone chain, each rung of a ladder and an SMC particle.

**What follows.**

| | An offset alone (probit's z, E23, E32) | Weights too (a latent scale) |
|---|---|---|
| ỹ, ȳ_w, ‖ỹ‖²_W | Read again, O(N) | Read again, O(N) |
| Columns | Kept | Whitened again, O(NM) |
| Gram, factor, gate | Kept: F does not read the values | Formed and factored strictly, O(NM²), and gated |
| b_c, β̂, r, ‖r‖², Q | O(NM), and a solve in O(M²) | The same, on the new factor |
| A prior's terms | Cross products only, O(M²) on a spectrum and O(p²) on a pencil | Decomposed, O(M³) or O(p³) |
| A draw's column means | Kept | Read again |

The weights' path runs the arithmetic of a design built there and loaded, to the bit. The offset's keeps the factor's incremental updates, so it differs from that design in rounding.

**Where the statistics live.** A state's fit names the working response it was fitted on (`BassStructureFit.Working`), and the target reads ȳ_w, the weighting and the exact-fit floor from it. So the target stays one target, holding what the priors read, frozen.
- *Why on the state:* σ²'s rate's reason. A swap hands over a state with its evaluation, an SMC stage holds one target for its particles, and the flat histogram wraps one. A state that carries what it read serves each of them with no change of theirs.
- *The alternative* was a target per working response through `ReplaceTarget`, the pattern of BARS's robust sweep and `BassTarget.Given`. It needs a target per replica, which a swap would have to exchange, and it has no seat on SMC.
- *The implication:* a state is its structure on its working response, and a design shows a state only when it holds both. Replicas that share a term array on different working responses hand over their designs as they do for different structures.
- *The floor is no prior.* It is the round-off of the current squares, so under Jeffreys it reads them, with s² standing in under a frozen read that states one.

**Decided: F under a moving W.** The gate reads the weighted Gram, so a latent draw can put the current structure below s_min, where its posterior is 0. Two options, and this is a design call for the user's review:

| | The gate reads W (built) | The gate reads the prior's columns alone |
|---|---|---|
| F | Moves with W | Fixed |
| The joint | π(S, z) ∝ π(S)·p(z)·1[S ∈ F(W(z))]·p(y \| S, τ, z): coherent, but F reads a parameter, so integrating z out gives the stated noise law only up to the mass the gate removes | The stated law, exactly |
| A latent block | A draw off the support is refused and the old state kept. A conjugate draw is then an independence proposal from the conditional truncated to F, which accepts every draw inside F, so the step is exact | Unchanged |
| A singular G | Never reached: the factor of G stays well conditioned | Admitted. ‖r‖²_W, Q and β̂ come from that factor, so the design needs a path with no factor of G: S = b(‖ỹ‖²_W − Σⱼ wⱼ²·b/(bλⱼ + τ)) on the pencil alone, proposals with no bordering, and the draw from the pencil |
| Reach | Every reading of the prior | `Unweighted` and `Frozen` only. Under `Weighted` the prior's columns are the weighted ones, so a singular G is a singular prior and the gate must read W |
| Cost per move | None beyond the Gram the move forms anyway | The same order, once the path exists |

- *BARS went the other way.* After its review follow-up its robust block admits on G = ZᵀZ alone, since a joint gate is exact only with a weight update that rejects draws refusing the current knots, and it makes the knot support depend on latent outlier weights ([the admission's notes](#bars-observation-models)). That is this table's left column, rejected. BARS could do so cheaply: it factors P_a = Zᵀ(W + τI)Z, which is positive definite whenever G is, and never holds the likelihood's Gram alone. BASS's design does hold it.
- *Why built this way here:* it needs no second path through the design, it covers all three readings, and it can be exact. The right column adds to this one, with no seam changed, so taking it later loses nothing.
- *What it costs, read on the two fixtures and a small N.* On a Gaussian chain's structures at s_min = 10⁻¹⁰, with scales drawn from their conditionals given the least-squares residuals, Friedman (N = 500, M = 11) and direct functional (N = 1,200, M = 5) refused none of 3,600 draws each over t(3), t(1) and a Laplace's. The smallest angle was 9.6·10⁻⁴, on a Laplace's scales, nearly seven orders above the floor. At N = 30 with a support floor of one row, where a column can rest on a single row, t(3) and t(1) refused none of 1,200 each, with angles down to 1.2·10⁻⁶, and a Laplace's refused 16 of 1,200, or 1.3%. So the truncation is out of reach at the fixtures' sizes and measurable once few rows carry a column. A refused draw costs a block one step.
- *What this reading leaves out.* Its residuals are least-squares ones at a Gaussian chain's structures, not a latent chain's own. Noise with a Laplace's latent scale in a chain that draws a can hold residuals smaller still, so E1 should read its own refusals; a block's accept counter reads them.
- *An option refused:* gating at structure moves alone and never at a latent move. F would then depend on the chain's history, and no joint has that support.

**Decided: what a model that moves must state.** With no reference, level centring and a σ² prior that reads s² would read the snapshot, conditioning the prior on parameters. `BassModel.PriorsReadResponse` says so, and the first move throws, naming the routes. The routes that read nothing are wider than the invariant one: a flat or zero intercept with Jeffreys, a declared inverse gamma or a held σ², so E32's and a mixture's zero centre qualifies.
- *When it is refused:* at the move, not at `BassModel.Create`. A model that never moves is unaffected, and the seam gains no member.
- *The alternatives.* A declaration on `IBassObservationModel` refused at `Create` would catch it one sweep sooner, for a member every implementer states. Freezing `BassPriorReference.Moments(y)` at the first move would pick a route for the user, and defaults are the user's. A warning leaves conditionals that share no joint.
- *`Weighted` under a moving W is allowed:* the prior moves with W by its own statement. A block that weighs the marginal is exact on that joint, as the enumeration below checks; a conjugate draw needs the correction of [the options' step 3](../../../issues/regression/bass/moving-response-priors-options.md#what-to-build-in-order).
- *What stays fixed.* The model's meaning, Gram and reference. A model with no weights takes no weights, since a design's layout depends on holding them, so a model whose weights move declares unit weights at the start. Weights that count rows do not move, since the support threshold and τ₀ read their sum at resolution, and nothing that hands them moves them.

**Decided: which target a composition reads.** A moved response asks nothing of the target, since the composition reads the design at the state. What was left is a σ² drawn apart, and a ladder's and a run's composition now read the target their sampling chain holds, as a lone chain's did.
- *So a replica's target has to be its state's.* A swap hands over a state with an evaluation that read its giver's σ². A replica that kept its own target would then weigh its next moves against a σ² drawn for another state, so it takes the handed one (`BassLadder`), and the swap moves the pair.
- *SMC and the flat histogram keep the model's,* since the first refuses another and the second wraps it: a σ² drawn apart has no seat on either yet.
- *Two entries read the model's target for a moment:* a resize's new replica, in burn-in, until its block runs, and a state from the reference rung, which carries no σ² of its own.

**What reads it.** A draw names the working response it was composed on (`BassDraw.Working`, null for the response itself), so these views read each draw's own.
- *The pointwise densities.* The retained content is N(y − offset | f, σ²/w) at the draw's own offset and weights. The collapsed content takes a block per structure and working response. Where the draws differ, `Latents` reports the offset and weights as observation-level latents conditioned on. Integrating them out is the noise's own density, which arrives with its observation model.
- *The training fit.* `Offset` is the draws' mean offset, and `Residuals` read y − offset − f̄. `Observed` stays y, so the name stays true, and `Mean` and `Quantiles` stay the surface's, as prediction's are. `NoiseSd` is null under weights, whose noise varies by row, as a basis fit's does by location.
- *The alternatives for it.* Reporting y − offset as `Observed` would make the name false. Adding the offset to the mean would read well for E23 and E32, whose offset is shared, and badly for a latent scale, whose offset has read yᵢ itself, so the fit would sit near the data it conditioned on. The field gives either view.
- *The residual check and the conditional means* read them too.
- *Prediction adds no offset back.* At a new input no offset exists; a consumer with a g(x) or an f_lo(x) adds its own.
- *A draw read on other data* reads that model's working response, since its own offset and weights say nothing there.
- *The cost:* a retained draw shares its working response with the chain, so N to 3N doubles per distinct snapshot retained, about 24 MB at N = 1,000 over 1,000 draws under a latent scale. A draw `with { Working = null }` drops it, and the views then read the model's. The collapsed content and the conditional means take a design per structure and working response, so under a latent scale about one per draw.

**Checked** (`BassMovingResponseTests`; `BassPriorReferenceTests.AMovingResponse_WalksTheSameStructures_OnKYPlusM`).
- *In place against built.* After each of three moves a design holds what a design built there holds, on both forms and all readings, scalar and direct with the Khatri–Rao pair: ȳ_w, the weighting, ‖r‖², Q, β̂, r, the angle and the conditional mean to 10⁻⁹, S and the marginal at three (τ, b) against a rebuilt model's target, and the proposals after a move. The weights' path to the bit.
- *In place against rebuilt.* A chain moved in place walks a chain rebuilt per snapshot for 320 sweeps, fifteen routes and forms: the same structures, τ, r and marginal, and each snapshot's draw to 10⁻⁸.
- *Equivariance.* Chains on y and on ky + m visit the same structures, through a moved chain as through a model per snapshot.
- *At each snapshot.* Detailed balance on all 2,263 structures of the rebuilt model to 10⁻¹³, the posterior against dense algebra to 10⁻⁹, and every structure through the moved path against that enumeration to 10⁻⁹, for level centring, the frozen and weighted Grams, and the flat and zero routes.
- *The joint.* A latent of three values, each its own working response, with a prior that evens their posterior. A chain's visits over (structure, z) against the enumerated joint, by Metropolis and by Gibbs on the latent, on seven readings and a ladder.
- *Off the support.* A move refused exactly where the enumeration says F leaves the structure, 720 refusals among 6,543 moves. In the joint, 3.6% of the mass sits on structures a move would leave F. A refused move leaves a chain bit for bit as if it was never tried.
- *What reads it.* Each view against a rebuilt model's, by snapshot.
- *Restart.* A chain's snapshot after a move, its columns re-evaluated under the moved weights, and a ladder's, with each replica's working response, continue bit for bit.
- *A ladder under a σ² drawn apart.* 25,000 swaps, the visits against the collapsed posterior, each replica's target the one its evaluation read, and each composed σ² the one held.
- *Whether the checks can fail.* 28 deliberate breaks, one at a time in an exported tree: the target reading the model's weighting, mean or floor; stale columns, responses or prior terms on a commit; a sync that ignores the working response; a view reading the model's. All are caught. Two survived the first pass, each a gap now closed: the conditional's draw reads Σw where its mean does not, and no test held a state no design holds.
- *Cost.* A rough reading under load. The unit path's draws are bit-identical, so blocks pair across clean builds of the commit before and this one: eight alternating rounds of nine blocks of 200,000 sweeps, pinned to the P-core threads, M identical in all 144 pairs (`BassWeightsCostTests`). Others kept about 2 cores busy, and 5 pairs met the quiet rule.
  - *The unit path,* after/before. Friedman at M = 9–10: wall 1.005 (IQR 0.93–1.13), CPU 1.017 (IQR 0.93–1.11). Direct functional at M = 5: wall 1.009 (IQR 0.95–1.09), CPU 1.028 (IQR 0.93–1.11). On minima, wall 0.974 and 1.025, CPU 0.978 and 1.018. The signs flip by fixture, so no slowdown is resolved, and a change of a few percent could hide in it.
  - *A move,* CPU µs per sweep with one move each sweep, on Friedman at N = 500 and M = 8–9, each line's walk its own, so read them against each other: unit weights 9.4 with none and 19.5 with an offset (21.9 under the ridge); the prior on the unweighted columns 45 with none, 62 with an offset and 188 with weights too; the prior on the weighted columns 7.8 with none and 85 with weights too. So an offset costs about 10–16 µs and weights about 80–140 µs. A latent scale's draw of N values is its own, on top.
- *Not checked.* A move on SMC or the flat histogram: the route exists and reads the same fit, but neither draws its start's latents, so neither would be exact until the latent blocks arrive. For the same reason a ladder's reference rung sits at the model's working response, which is checked only for its mechanics. A run's composition under a target a block handed it, since `Bass.Fit` takes no observation model. A σ² drawn apart across a restart, which a chain's snapshot does not hold.

### BASS: noise beyond the Gaussian (E1)

GBASS's noise on this core, a normal variance-mean mixture whose working response and weights move with its latents, has its own page: [Variance-mean mixture noise](../VarianceMeanMixtures.md), with Student t (robust BASS) and the asymmetric Laplace (quantile regression) built so far. It fills the observation-model and sweep-slot seams above, on the weighted core and the moving working response below.

### BASS: random-effect intercepts (E21)

Built 2026-11-08 as KHAOS's slice 3 ([the report](../../../issues/regression/bass/khaos-brief.md#reports)), as a capability of the core: any grouping, any coefficient form, hinge and Legendre terms alike. KHAOS is one consumer. Its extensions, a drawn s, several groupings, slopes, the block under GBASS's noises and its conformal reading, are on [RandomEffects.md](../RandomEffects.md).

**The one-group model.** `BassConfig.RandomIntercepts = [new BassRandomIntercept("site", s)]`, on data that hold the named grouping, `new BassData(predictors, response, [site])`. The grouping is a `CategoricalPredictor` beside the predictors, and no term reads it. Named crossed and nested collections extend this form in [RandomEffects.md](../RandomEffects.md#named-crossed-and-nested-intercepts-slice-2).

$$y_i=\alpha+\sum_j\beta_jB_j(x_i)+u_{g(i)}+\varepsilon_i,\qquad u\mid\sigma^2\sim N(0,\ \sigma^2s^2\Pi)\ \text{on the }L\text{ observed levels}.$$

- *s* is the effects' standard deviation in units of σ, and it is held.
- *The convention* (`RandomInterceptConvention.ObservationCentred`, the only one built).
  - Π projects onto Σ_g n_g u_g = 0. Zu then has mean zero over the cases, and the intercept alone carries the level, as level centring reads it.
  - The unweighted sum-to-zero is not built. On unbalanced levels Zu keeps a nonzero mean, which is the overlap with the intercept that the convention exists to remove.
  - The qualification shows the difference: the fitted level sits at the surface's mean plus Σ n_g u*_g / n, not plus the levels' unweighted mean.
- *The levels* are the ones the training cases take. A declared level no case takes has no effect, and prediction treats it as unseen.

**Integrated by whitening** (`RandomInterceptBlock`).
- *The covariance.* With Z the level indicators, D = ZᵀZ and H an orthonormal basis of n's complement, y has covariance σ²V, with V = I + s²(ZH)(ZH)ᵀ.
  - V1 = 1, so V^(−½) keeps the intercept's column, every column's mean and ȳ.
  - With HᵀDH = EΛEᵀ, V^(−½)x = x + Z·M·Zᵀx, where M = HE·diag(((1 + s²λ)^(−½) − 1)/λ)·EᵀHᵀ. That costs O(n + L²) a column, and ln|V| = Σ ln(1 + s²λ).
- *What is whitened.* The model whitens the working response, as an offset y − V^(−½)y, and the design whitens each column as it evaluates it. The support check reads the raw column.
- *The price* −½ln|V| rides on the weighting's Σ ln w: `WorkingResponse`'s `logJacobian`, which the model sets to −ln|V|. Every marginal then carries it at its b: each rung, each particle, and the evidence.
- *σ²'s scaled prior* (`Sigma2Prior.ScaledToResponse`) reads the whitened response's mean square. That is σ²'s own scale there, the noise within the levels.
- *The coefficient forms* read the whitened design as they read any design. The g-form's and the modified g-prior's Gram becomes ΨᵀV⁻¹Ψ, and the ridges take their columns whitened.
- *Why whitening, not block columns.* The brief placed the block on slice 2's prior-precision builder, `BassPriorPrecision`. As columns of a ridge block there, the effects would be L more columns in every structure's factor and residual sum of squares, through every move, rung and particle.
  - Whitening leaves the factor, the samplers and their caches as they are, and its marginal is exact.
  - Its cost is a held s: a drawn s re-whitens every column, O(nM) per value.
  - This was decided in slice 3, beyond the brief.
- *The group effects* are drawn last in each composed draw, from their exact conditional given α, β and σ²:
  - γ ~ N((HᵀDH + I/s²)⁻¹HᵀZᵀr, σ²(HᵀDH + I/s²)⁻¹) and u = Hγ, with r = y − α − Xβ on the raw columns. The pair is a joint draw.
  - `BassDraw.Groups` holds each named grouping's scale and u over its observed levels.
  - A fit without the block draws as it did, and its hashes hold.

**Readouts.**
- *Prediction.*
  - A case whose grouping column names an observed level reads f + u_g.
  - A case without the column, a population case, reads f, and its predictive draw adds its own effect σs·z.
  - An unseen level is refused (`UnseenLevel.Refuse`) or placed at the population (`OutsideEverySubset`).
  - The mean surface adds E[u | S, τ, y]: the effects' conditional mean at the coefficients' conditional means, since it is linear in them and free of σ².
- *The training fit* is the prediction at the training cases, f + u_g. It puts back no offset, since the whitening is not one.
- *Sensitivity and the gradient* read f alone, since u is flat in x.
- *The pointwise density.* The retained content is N(yᵢ | fᵢ + u_g, σ²), and it reports the group effects as conditioned on, at the group's scope.

**The engine's reach.**

| Path | Under the block | Gate or reason |
|---|---|---|
| A chain; the ladder, its reference rung included | Runs | Each rung's law against the enumeration, on a hinge fit (the g-form, τ drawn) and a Legendre fit (the modified g-prior, g₀² drawn) |
| SMC and its flags | Runs | ln p(y) against the enumeration, hinge and Legendre, with the flags clear |
| The ladder's evidence | Runs | Against the enumeration, hinge and Legendre |
| The flat histogram | Runs | Against the enumeration, Legendre |
| Resume | Runs | Bit for bit, the group effects included |
| The retained pointwise content | Runs | Against the dense Normal at f + u_g |
| The collapsed pointwise content | Refused | Leaving a unit out would leave it out of the block too, whose whitening couples every case of its level |
| The candidate density, and so conformal prediction | Refused | A new case's group effect, and the exchangeability across groups that the conformal layer reads, need a reading of their own |
| Noise beyond the Gaussian (E1) | Refused | The block whitens under unit weights, and a noise's latent weights move every sweep |
| A direct functional response; a decomposition's components | Refused | The block serves a scalar response, and a decomposition's components share no grouping |

**Checks** (`RandomInterceptTests`, `RandomInterceptQualificationTests`).
- *Dense algebra.*
  - The whitening: W·V·W = I to 3·10⁻¹⁵, W·1 = 1, and ln|V| against a dense Cholesky.
  - The effects' conditional: the mean to 10⁻¹⁵; the covariance and Σ n_g u_g = 0 over 200,000 draws.
  - The target's ln p(y | S, τ) against the dense multivariate t, absolutely, to 10⁻¹³. This holds for the raw ridge and the modified g-prior, with and without the block.
  - The mean surface against the dense joint posterior mean, to 2·10⁻¹⁴.
- *Recovery from outside the model.*
  - The generator: twelve unbalanced levels (4 to 37 cases), and effects from a centred exponential, sorted onto the levels by size, which the block's Normal prior does not hold. There are four datasets and s = 2 against a generator spread of 0.8/0.3.
  - The hinge fit's 90% intervals cover 44 of 48 effects and the Legendre fit's 43, with standardized errors' mean squares of 1.03 and 1.04.
  - The level sits within 2 posterior sd of the convention's target. Its squared gaps sum to 0.002, against 0.35 from the unweighted mean's.
  - σ comes out at 0.29 to 0.33 against the generator's 0.3. Without the block it comes out at 0.38 to 1.0.
  - At held-out cases of the observed levels, the RMSE against the true mean is 0.11 to 0.15. Without the block it is 0.39 to 1.0.

**Not built** (follow-ups):
- a drawn s, a variance component, which re-whitens per value;
- more than one grouping, crossed or nested;
- random slopes;
- the collapsed pointwise content, the candidate density and conformal prediction;
- noise, weights and a functional response under the block.

