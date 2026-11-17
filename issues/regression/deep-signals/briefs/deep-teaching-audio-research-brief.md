# Deep Teaching principles for an unsupervised audio model

Date: 2026-11-17. Status: **research proposal**. This brief records the research direction and candidate mathematical questions; it commissions no implementation, probes, or Lean slices. The proposed results below have not been proved. Their novelty relative to the wider literature remains to be assessed.

## Purpose

Develop a native audio foundation model from human conversation, without human annotations and with efficient use of recordings. Investigate whether principles suggested by Voroninski's inverse-problem work and Helm.ai's public disclosures can connect three operations: learning a signal prior, recovering useful structure through that prior, and teaching a model from the recovered structure.

The research hypothesis is that appropriately chosen observable views and structural constraints can supply reliable training targets, including distributions when a point target is not identifiable. A useful result must explain which information the target retains, how errors enter it, and why learning from it helps the intended prediction task.

This develops the [approach brief](approach-brief.md). It leaves the generator architecture open. Speech-production physics can inform constraints and observation operators; a neural generator can supply the expressive prior. A full articulatory simulator is not a prerequisite established by this research.

## What the public evidence supports

### Voroninski and collaborators

The relevant local sources are in `D:/aghado01/codex-scientiae/supellex-private/gauntlet/voroninski`. Versioned directories contain the TeX sources. The following public links identify the results, including papers written during the Helm.ai period.

| Source | Established contribution | Use here and boundary |
| --- | --- | --- |
| [PhaseLift, 1109.4499](https://arxiv.org/abs/1109.4499) | Exact recovery modulo phase and noise stability for specified random intensity measurements via convex lifting. | Missing information in one observation can be recovered from a suitable collection. This does not establish the same result for a chosen audio window family. |
| [ShapeFit, 1506.01437](https://arxiv.org/abs/1506.01437) | Recovery of global locations from corrupted pairwise directions under stated geometric and graph assumptions. | Redundancy and consistency can support robust inference. Its location-recovery theorem is not an audio theorem. |
| [Global Guarantees, 1705.07576](https://arxiv.org/abs/1705.07576) | Landscape guarantees for inversion through specified random expansive generators; descent away from neighborhoods of two distinguished latent points. | A generator supplies a restricted signal class and coordinates for optimization. A generic claim of no troublesome stationary points would overstate the result. |
| [Provably Convergent Scheme, 1812.04176](https://arxiv.org/abs/1812.04176) | Algorithmic recovery to the noise level under random expansive-network assumptions. | Connects favorable geometry to an actual solver. Trained audio generators require their own hypotheses. |
| [Generative phase retrieval, 1807.04261](https://arxiv.org/abs/1807.04261), [2008.10579](https://arxiv.org/abs/2008.10579) | Landscape and algorithmic results with measurement complexity scaling with latent dimension, with depth and logarithmic factors. | Structural priors can change recovery complexity. The generator is supplied to the recovery problem. |
| [Deep Denoising / Rate-Optimal Denoising, 1805.08855](https://arxiv.org/abs/1805.08855) | Denoising guarantees for specified generator-based optimization and encoder–generator procedures. | Connects structural inference to feedforward estimation. It supplies no general guarantee for learning the generator from unlabeled speech. |

**Measurement efficiency and learning efficiency are separate claims.** A bound on measurements per reconstructed signal, given a generator, does not bound the number of recordings required to learn that generator. Likewise, population identifiability does not by itself give a stable estimator or a finite-sample learning rate.

### Public Helm.ai mechanisms

The four US grants examined belong to three families. They disclose concrete procedures, rather than a complete reproducible specification of the original Deep Teaching system. Each audio connection in this table is a research inference.

| Disclosure | Procedure described | Possible audio connection |
| --- | --- | --- |
| Camera-invariant 3D property formulation: [US12579725B1](https://patents.justia.com/patent/12579725), continuation [US12694601B1](https://patents.justia.com/patent/12694601) | Predict normalized geometric quantities and combine them with calibration. Improve training labels with tracks and physical measurements. Surrounding frames, including future frames, can denoise labels for a single-image predictor. | Give an offline teacher more observations and computation than the streaming student; explicitly handle known recording transformations. |
| Label image synthesis using generative AI: [US12205197B1](https://patents.google.com/patent/US12205197B1/en), claims 1 and 10 | Condition synthesis on segmentation labels and scene prompts, then use the image and conditioning labels as training data. | Generate variations with known target correspondence. The source of the initial structure and its validity on real audio remain separate questions. |
| Road surface detection and prediction: [US12643557B1](https://patents.google.com/patent/US12643557B1/en) | Train image-based road-condition prediction using physically sensed traction information. | Use an additional observable, possibly available only during training, to supply targets for a restricted-input predictor. |

Helm's [2020 announcement](https://helm.ai/post/deep-teaching) describes training without human annotation or simulation. Voroninski's [December 2025 account](https://helm.ai/post/factored-embodied-ai) describes perception learned from unlabeled video, a geometric interface to policy, and subsequent simulation in that representation space. These are company descriptions. The later simulation and synthesis components should not be substituted for a disclosed recipe for the original unsupervised method, and the company's performance claims are not consequences of the recovery theorems.

## Mathematical setup and information contracts

Let a waveform be x, a candidate generator be G_theta(z), and the known view operators be M_j. A frozen-generator recovery step has the schematic form

\[
\widehat z\in\operatorname*{argmin}_z\sum_j
\ell_j(M_j(G_\theta(z)),y_j),\qquad y_j=M_j(x)+e_j.
\]

Learning theta is an additional statistical problem. The target of that learning might be a waveform law, a conditional future law, or a useful function of either. Latent coordinates themselves are generally nonunique under changes of parameterization; every identification claim must name its equivalence relation.

The views may be deliberately constructed from a fully observed recording. Magnitude-only observations are one choice of learning objective, not a restriction imposed on every stage of the project. A negative result about separate population power laws does not imply impossibility when raw waveforms or paired temporal observations are retained.

Every candidate theorem must specify:

- **Training information:** raw records, paired view tuples, or separately labeled marginal laws; exact population information versus finite samples.
- **Deployment information:** the causal input U available to the student, distinguished from the richer offline observation V.
- **Target:** an observable Y, a recoverable functional, a conditional law, or a latent with a defined interpretation.
- **Model class:** source family, generator restrictions, stationarity or temporal dependence, positivity and signal-amplitude bounds.
- **Allowed transformations:** which changes preserve the target, change it by a known rule, or remove useful information. Speaker identity, prosody and loudness are not universally disposable nuisance variables.
- **Sampling unit:** independent recordings or an explicitly justified dependence model. Multiple windows from one recording do not automatically supply independent examples.

An offline teacher T(V) may produce a point estimate, distribution, or set. Under squared loss, an unrestricted population student trained on T learns E[T|U]. It predicts E[Y|U] only when the appropriate conditional-bias condition holds. A teacher using future observations cannot transfer their realized unpredictable content to a causal student. These conditional-expectation facts are standard tools; the research concerns when the proposed inference procedure satisfies the needed conditions.

The direct predictor of a held-back observable is an essential comparator: its target is already available without human labels. A structural teacher must justify its additional assumptions through reduced estimation error, useful constraints, or a demonstrated learning-efficiency advantage.

## Existing mathematical boundary

These summaries were checked against the current briefs and declaration names. Their verification status comes from the linked proof reports; this documentation pass did not rerun Lean.

| Existing result | Consequence for new work |
| --- | --- |
| [Measured stability](phaseless-audio/04-measured-stability-brief.md): quotient-error transfer from a positive measured lower gain; explicit cut, scaling and blind-window obstructions. | Prove a concrete gain or a target-specific substitute on a meaningful class. Repeating a residual triangle inequality is not a new slice. |
| [Paired model identification](phaseless-audio/05-paired-model-identification-brief.md): paired powers identify laws on a common compact sign quotient when they separate its points. | Establish separation for the chosen class, or determine the smaller target that remains identifiable. |
| [Unpaired model identification](phaseless-audio/06-unpaired-model-identification-brief.md): finite-dictionary support and marginal-incidence criteria, with physical Fourier witnesses. | Use the criteria to design measurements and quantify ambiguity on a specified source class. |
| [Context-tree memory](context-tree-memory-identifiability-brief.md): under stationary sign symmetry, positive block masses and a known memory bound, separate exact rectangular DC-power laws identify the binary source through depth three. At every greater depth there is a stationary collision, surviving full rectangular powers jointly within each width through D+1. One asymmetric taper distinguishes that constructed pair. | A pair-separating taper does not prove recovery of the whole class. A constructive completion theorem is the next identification question. |
| [Conformal audio](conformal-audio/README.md): observable coverage, ambiguity, transport, recorded risk and alarms have separate contracts. | Calibration can assess a specified observable target. It cannot resolve a measurement collision or supply unproved synthetic-to-real conditional-law agreement. |

The [phaseless sequence](phaseless-audio/README.md) locates the distinct contributions of Jaganathan–Eldar–Hassibi, Grohs–Rathmair and Tachella–Chen–Davies. In particular, model-set identification, probability-law identification and identification of a generator parameterization must remain separate.

## Proposed Lean questions, in dependency order

### 1. Complete the first unresolved measurement boundary

**Question:** For every stationary, sign-symmetric binary source of memory at most four with positive length-five word masses, do the separate rectangular DC-power laws at widths one through five, augmented by the two four-site DC-power laws below, determine the transition kernel and its minimal context tree?

Use signs X_i in {-1,1}, zero-based positions in a length-five block, and unnormalized powers. Define

\[
Y_w=\left(\sum_{i=0}^{w-1}X_i\right)^2,\qquad
Z_S=\left(\sum_{i\in S}X_i\right)^2,
\qquad A=\{0,1,2,4\},\ B=\{0,1,3,4\},\ C=\{0,2,3,4\}.
\]

The proposed augmentation observes the laws of Z_A and Z_B separately. It needs no correspondence between their individual samples. The proposed sharpness statement is that **one extra four-site indicator taper on the same five sites cannot identify the entire class**, even if the original rectangular observations retain all frequency powers jointly within each width. This minimality claim is only within that taper class. Arbitrary real weights, extra windows or paired measurements are different experiments.

There is a concrete algebraic route to test before commissioning the proof:

1. Existing power-mean second differences recover r_1 through r_4, where r_k=E[X_0 X_k]. Width-four extreme power recovers q=E[X_0 X_1 X_2 X_3]; stationarity also supplies the shifted four-way moment.
2. Write q_S=E[product of X_i over i in S]. Each selected four-site law gives
   \[
   q_S=8\Pr(Z_S=16)-1-\sum_{i<j\in S}r_{j-i}.
   \]
3. The width-five extreme event gives
   \[
   16\Pr(Y_5=25)=1+4r_1+3r_2+2r_3+r_4+2q+q_A+q_B+q_C.
   \]
   Thus the two new laws appear to recover the two remaining independent degrees of freedom. Sign symmetry removes odd moments; finite Walsh inversion would recover the length-five law, then positivity permits division to obtain transitions.
4. For necessity, examine the two-parameter family
   \[
   p_{a,b}(x)=2^{-5}\{1+a\chi_A(x)+b\chi_B(x)-(a+b)\chi_C(x)\},
   \qquad |a|+|b|+|a+b|<1,
   \]
   where chi_S is the product over S. Each character contains both endpoints, so proper contiguous marginals are fair and both length-four endpoint marginals agree. Check that this defines a positive stationary order-at-most-four source. Cyclic-shift invariance of the full rectangular Fourier-power vector is the proposed reason its law cancels the zero-sum quartic perturbation. Any one four-site indicator taper appears to retain only one linear combination of a and b.

All sufficiency and necessity steps above are **candidate derivations**, not a verification report. Connect them to actual `framePower` observations and actual stationary sources. An acceptable result is a proved completion/minimality theorem or a counterexample with the corrected boundary. Counting moments alone is insufficient.

**Why this matters:** it would identify an exact, modest change in an observation objective that makes previously invisible temporal structure learnable. It extends the existing depth cutoff to a constructive design decision. General depth D is a subsequent question; two tapers at depth four imply no uniform result for larger memory.

### 2. Determine the predictive cost of discarded structure

**Question:** What unavoidable next-symbol prediction error is imposed by training information that cannot distinguish the stationary sources above? Which predictive quantities can still be learned without identifying the full source?

For two indistinguishable sources with a shared history law pi and next-positive probabilities q_0(h), q_1(h), define Brier excess

\[
\mathcal E_i(a)=\sum_h\pi(h)(a(h)-q_i(h))^2.
\]

The standard two-point argument gives max_i E_i(a) at least one quarter of the squared L2(pi) distance between q_0 and q_1; the midpoint predictor attains the bound for that known two-source decision problem. The substantive task is to derive the exact separation from the stationary Fourier collision, characterize which predictive functionals escape it, and show how an admissible measurement augmentation changes the bound. Do not commission the abstract two-point inequality alone.

The learner in this statement receives only the specified population summaries, plus the declared deployment history. If it may inspect full ordered training recordings, the indistinguishability premise must be re-established or withdrawn.

**Why this matters:** it connects the observation objective to an actual prediction loss. A representation can be adequate for a target even when waveform or source identification fails; the target must be named.

### 3. Quantify stability on the completed class

**Question:** How accurately must the observable law summaries be estimated to control transition predictions or teacher targets?

Start with the class and measurement family settled in question 1. Specify the norm on histogram or moment error and assume a positive lower bound eta on history mass when dividing to form transition probabilities. Derive explicit conditioning constants, attainable bad cases, and their dependence on eta. Separate recovery of the block law from division into conditional laws. Determine whether a predictive target has better conditioning than full-law recovery.

The later neural-prior version must establish conditioning on the chosen generator range or target, and quantify model mismatch. A hypothesis asserting the desired inverse bound is not itself that result. Existing `MeasuredStability.residual_error` already handles transfer once a lower gain is available.

**Why this matters:** exact population identifiability can coexist with a teacher too unstable to be useful. This result gives tolerances for an eventual estimator or probe and a reason to compare measurement designs.

### 4. Transfer an offline teacher to a causal student

**Question:** When does inference with richer observations yield the correct conditional prediction from the student's restricted input?

Specify Y, U and V, with U measurable from V where a nested-information argument is used. An exact teacher T=E[Y|V] has E[T|U]=E[Y|U]. For an approximate teacher T=E[Y|V]+e, its unrestricted population squared-loss student differs from the Bayes predictor by E[e|U]. Approximation and finite-sample student errors are additional terms.

These identities are standard scaffolding. The research target is to connect the measurement and model errors from question 3 to the relevant conditional bias for a specified temporal source family. Include an explicit failure case: internally consistent inverse solutions may select a latent or source branch that biases causal predictions. If a target is not determined by observations, define a distribution or set target and its loss before asserting a transfer guarantee.

**Why this matters:** this is the precise mathematical connection suggested by the patent's richer training observations. An accurate offline reconstruction and a useful causal teaching target are distinct properties.

### 5. Establish a learning-efficiency claim

**Question:** Under what conditions does the structural teacher or representation reduce the number of independent recordings needed for a specified prediction accuracy?

Choose the same observable target and deployment input for the teacher-assisted and direct held-back-prediction estimators. Specify hypothesis classes, teacher training data, estimation procedures and dependence assumptions. Account for learning the prior itself, including any separate data it consumes. The desired bound separates target ambiguity, teacher bias, estimation variance and student approximation error; it must say where an improvement occurs and exhibit a regime in which it does not.

This follows the preceding questions rather than inheriting a sample-complexity claim from compressive recovery. A result about a finite source class can establish a mechanism worth generalizing. No result about arbitrary speech networks is presumed.

## Scope, sequencing and acceptance

Begin by reviewing the candidate depth-four completion theorem and its minimality class. Questions 2 and 3 then measure what the added observations buy. Question 4 connects a reliable teacher to deployment inputs. Question 5 addresses the ultimate data-efficiency claim. No conversation recording is needed to resolve these first mathematical questions.

For each later slice brief, require a concrete source class, actual observation operators, exact quantifiers, a useful positive statement and an obstruction or sharpness case where applicable. Standard factorization, conditional-expectation, linear-algebra and risk identities are supporting lemmas. The main deliverable must settle a project decision that those identities alone leave open.

Use the [Lean harness conventions](../../../../lean/README.md), scoped imports and the existing modules. Record hypotheses remaining as premises, declaration paths, axiom checks, and build results in the owning slice's report. Distinguish established literature, project-proved statements, candidate extensions and empirical assumptions. A formal proof establishes the stated result, not its historical novelty or its empirical applicability to speech.

Generator selection and engine implementation follow the resulting contracts. Borrowing a voice-network architecture and borrowing pretrained weights are different research choices; the latter must account for any transcript or other annotation supervision used in pretraining. Synthetic labels are exact relative to the chosen generator, while their relevance to real audio remains a model-validation question.

## Preparation record

- Reviewed the public paper and patent material in the preceding research discussion; source links above preserve that evidence.
- Checked current reports for the phaseless slices and context-tree result, and located the corresponding Lean declarations. Reused established results instead of assigning them again. No new proof or numerical claim is reported here.
- Read the existing approach and conformal sequence. This brief supplements their research direction without changing their reported outcomes.
- Investigation anchors are retained under `artifacts/doc-dive/deep-teaching-research/`; authored findings identify this coordinating session.

### Doc changes owed

After a proposed theorem is verified, update the approach's measurement and teacher claims with its precise source class, observation regime and stability assumptions. Keep this research brief's conjectures visibly separate from subsequent proof reports.
