# KITT, Formalized: The Math and Engineering of Coupling a Helm.ai-Style Factored Driving Stack with a Dessa-Like Full-Duplex Voice Model

The cleanest way to get a KITT-like car is to treat the voice model as a **second, read-mostly consumer of the driving stack's factored semantic-geometric state**. In practice that means adding a synchronized, rate-aligned *ego stream* to a Moshi/DuplexSLA-style multi-stream duplex transformer, and a narrow, typed *intent stream* that the planner treats as a bounded soft preference behind a safety filter. Do not build a single end-to-end "drive-and-talk" model. The factored interface is what makes three things possible: authority separation, testable faithfulness of explanations, and semantic-space simulation for conversation as well as driving.

## TL;DR

- **Architecture:** Helm.ai's split of perception and policy over a "semantic geometry" interface is exactly the seam a conversational model should attach to. Tokenize that interface (object slots, predicted trajectories, planned maneuver, criticality, planner attributions) as an extra conditioning stream in a full-duplex model, and emit a rate-limited intent channel, as DuplexSLA does with its ≤10-token-per-160 ms action channel. The planner consumes intents only as bounded cost terms inside a hard-constraint set enforced by a control-barrier-function (CBF) safety filter. Intents may shrink the admissible set but never enlarge it.
- **The hard science problems are timing, faithfulness and alignment, not plumbing.** (1) Turn-taking must be conditioned on scene criticality and driver workload. This is formalizable as optimal stopping or constrained RL, and trainable with DuplexPO-style event-window GRPO. (2) Verbal explanations must be *read out* from what drives the planner. A May 2026 study (arXiv 2605.17268, accepted at the CVPR 2026 DriveX workshop) found that for Alpamayo-R1 "overall reasoning fidelity is only 42.5%," with chain-of-causation reasoning "matching scene reality less than half the time," which shows that generated rationales cannot be trusted by default. (3) The speech and driving representations need geometric regression tests (principal angles, Procrustes, CKA) whenever either model is updated.
- **Feasible now on one H200:** LoRA-fine-tune an open Moshi-family model with a synthetic ego/workload stream. Train it so that speech onsets are suppressed in high-criticality windows and hazard alerts barge in. Evaluate with context-conditioned timing metrics, Full-Duplex-Bench regression, and causal patching of the ego embedding. Dessa itself cannot be used: Block has released no weights or API. Helm.ai's architecture is known only from press releases and four narrow patents, so treat its performance claims as unverified.

---

## Key Findings

1. **What is verifiable about Helm.ai is the architecture's shape, not its performance.**
   - Helm.ai describes a "Factored Embodied AI" framework (Dec 11–12, 2025). A perception layer and "Geometric Reasoning Engine" extract 3D structure and semantic segmentation. A policy is trained on that "semantic geometry" in "semantic space" simulation. A world model projects "ghost trails" to predict agent intent.\[1\]\[2\]\[3\]
   - The 1,000-hour real-world fine-tuning figure and "zero-shot" urban steering come from Helm.ai's Dec 11, 2025 Business Wire release, which describes "a continuous 20-minute drive without steering disengagement" in Torrance. Neither has been independently benchmarked.
   - Deep Teaching is described only as "proprietary unsupervised learning," with a vendor-stated link to "applied mathematics and compressive sensing."\[4\]\[5\]
   - Helm.ai's four granted US patents cover label-to-image generative synthesis, camera-invariant monocular 3D estimation, and traction-supervised road-surface prediction.\[6\] None discloses the Deep Teaching method.
2. **Dessa is a capability signal, not a component.**
   - Block (Oct 23, 2026) says Dessa is full-duplex, trained from scratch on human conversational audio with no text or pretrained model, trained on a single GB200 NVL72, and runs locally on a MacBook Pro M4 Max or better.\[7\]\[8\]\[9\]
   - No weights, API, parameter count, or paper have been released.\[10\]\[11\] All performance claims are vendor claims.
   - Engineering work today must use the open Moshi lineage (Moshi, PersonaPlex, Lychee-FD, DuplexSLA).
3. **The open duplex literature already contains every building block except the ego stream.**
   - Moshi has parallel user/model streams plus a time-aligned text "inner monologue," with 160 ms theoretical and 200 ms practical latency.\[12\]
   - DuplexSLA co-decodes assistant audio with a structured action stream on a shared 160 ms chunk clock.\[13\]
   - DuplexPO trains *timing* separately from *content* with RL.\[14\]
   - GPT-Live delegates hard reasoning to a background model while the conversation continues.\[15\]
4. **Human-factors evidence favors a "passenger-like" assistant.**
   - In Drews, Pasupathi & Strayer's 2008 simulator study (J. Exp. Psych.: Applied 14(4)), "the number of driving errors was highest in the cell phone condition," and a University of Utah release reported that cell-phone drivers "drift out of their lanes and missed exits more frequently" than drivers talking with a passenger.
   - Passengers referenced traffic more often, and both partners modulated the complexity of their speech with driving demand.\[16\]\[17\]
   - An ego-aware voice model can reproduce that modulation. A bolted-on assistant structurally cannot.
5. **Regulation and safety engineering push toward strict partitioning.**
   - EU Advanced Driver Distraction Warning (ADDW) has been mandatory for all newly registered vehicles since July 7, 2026, and must work without biometric identification.\[18\]\[19\]
   - Automotive SoCs such as DRIVE Thor advertise multi-domain isolation for running driving and infotainment workloads concurrently.\[20\]\[21\]
   - The voice model should be a QM-class (non-safety-rated) consumer with a validated, one-way data path from the ego stack and an arbitrated return path.

---

## Details

### 1. Ground truth: what each side actually is

| Component | Verified (papers / filings / open weights) | Vendor claim only |
|---|---|---|
| Helm.ai perception | Patents on camera-invariant 3D estimation (US 12,579,725; US 12,694,601), label-to-image synthesis (US 12,205,197), and traction-supervised road-surface prediction (US 12,643,557).\[6\] Founders' peer-reviewed math work: PhaseLift (CPAM 2013), deep generative priors (COLT 2018), phase retrieval under generative priors (NeurIPS 2018)\[22\]\[23\] | "Geometric Reasoning Engine," Deep Teaching as "unsupervised learning" from non-driving data |
| Helm.ai policy | None public | Planner trained on semantic geometry in semantic simulation; 1,000 h of real-world data; zero-shot steering (Torrance, Redwood City demos); "ISO 26262-certifiable" path\[24\]\[25\] |
| Helm.ai world model | None public | Intent prediction via "ghost trails"; generative adversarial scenarios in semantic space\[3\] |
| Dessa | None public | Speech-only training; full-duplex; backchannels and interruptions learned from data; local real-time on M4 Max; "hundreds of times fewer words"\[7\]\[8\]\[26\] |
| Open duplex models | Moshi (arXiv 2410.00037, CC-BY weights), PersonaPlex (Moshi-based, NVIDIA Open Model License), Lychee-FD (ACL 2026, open weights), DuplexSLA (arXiv 2605.20755; at time of a September 2026 citation, model not yet released)\[10\]\[12\]\[27\]\[28\]\[29\] | — |

The founders' academic background is real and relevant. The intellectual lineage of compressed sensing and generative priors suggests Deep Teaching *might* exploit low-dimensional structure priors for label-free learning.\[30\]\[31\]\[32\] That link is an inference from the founder's publication record. Helm.ai has not technically disclosed it.

### 2. Architecture overview

```
 Cameras ──► [Helm-style Perception + Geometric Reasoning] ──► Semantic-geometric state S_k (10–20 Hz)
                                                                   │            │
                                    ┌──────────────────────────────┘            │
                                    ▼                                           ▼
                     [World model: agent intent, ghost trails]        [Policy / Planner] ─► [CBF safety filter] ─► actuators
                                    │                                     ▲   │ (attributions, active constraints,
                                    ▼                                     │   │  accepted/rejected intents)
                ┌───────── Ego tokenizer (typed, versioned interface) ◄───┼───┘
                │                         E_t  @ 12.5 Hz                  │
                ▼                                                         │
 Mic array ► AEC/zone/speaker attribution ► U_t ─► [Full-duplex speech model] ─► A_t (speech) ─► cabin speakers
                                                   │   streams: U, A, M (monologue), E (ego), I (intent)
                                                   ├─► I_t (intent tokens, ≤ few/frame) ─► [Intent validator/arbiter, ASIL-rated] ─► planner prefs
                                                   └─► delegation calls ─► [Slow deliberator: trip context, search, reasoning]
```

Design invariants:
- **Read-rich, write-narrow:** E is high-bandwidth and one-way into the voice model. I is low-bandwidth, typed, and arbitrated.
- **The ego stack owns all safety-critical speech.** Takeover requests and collision warnings come from deterministic HMI paths. The generative model must yield the floor immediately when they fire.
- **Every utterance about the vehicle's behavior should be traceable** to an E-stream field such as an attribution, an active constraint, or an intent-acceptance status.

### 3. Component math

#### 3.1 Ego state as a synchronized stream in a multi-stream duplex model

**Streams.** Use a frame index $t$ with $\Delta_f = 80$ ms (Mimi at 12.5 Hz). Define:
- $U_t = (u_t^{1..Q})$: user audio codebook tokens (Q residual codebooks).
- $A_t = (a_t^{1..Q})$: assistant audio tokens.
- $M_t$: inner-monologue text token(s).
- $E_t$: ego tokens (conditioning only).
- $I_t$: intent tokens (output).

**Factorization.** Moshi-style, with a temporal transformer across frames and a depth transformer within a frame:

$$
p(A, M, I \mid U, E) = \prod_{t} p_\theta\big(M_t, I_t, A_t \mid A_{<t}, M_{<t}, I_{<t}, U_{\le t}, E_{\le t}\big)
$$

Within a frame, decode in the order $M_t \rightarrow I_t \rightarrow a_t^{1} \rightarrow \dots \rightarrow a_t^{Q}$ (semantic before acoustic, with Moshi's acoustic delay of 1–2 frames on the fine codebooks):

$$
p(M_t, I_t, A_t\mid\cdot)=p(M_t\mid h_t)\,p(I_t\mid M_t,h_t)\prod_{q=1}^{Q}p(a_t^{q}\mid a_t^{<q},I_t,M_t,h_t),
$$

where $h_t$ is the temporal transformer's context vector. $U$ and $E$ are never predicted at inference. During training you may add an auxiliary world-model loss $\mathcal{L}_{\text{aux}} = -\log p(E_{t+k}\mid h_t)$ to force the backbone to model ego dynamics (predicting a criticality spike before it happens). Keep its weight small to avoid diluting conversational competence. Lychee-FD's finding that gradient conflict between acoustic and semantic objectives degrades knowledge is a warning that heads competing for one backbone interfere.\[29\]\[33\]

**Precedent.** DuplexSLA already runs a three-channel structure on one 160 ms clock: two 80 ms causal user features, four 40 ms assistant tokens with a text anchor, and up to 10 action tokens per chunk.\[13\] Adding E is a fourth, input-only channel. NemotronLabs VoiceChat explicitly prefers parallel specialized streams over serializing heterogeneous tokens into one channel, for latency reasons.\[27\] This argues for E as its own embedding path rather than text inserted into the monologue.

**Rate alignment.** Perception and planning emit snapshots $S_k$ at times $t^{e}_k$ (10–20 Hz) with pipeline latency $\lambda$. Speech frame $t$ may consume only snapshots available by its start:

$$
\kappa(t) = \max\{k : t^{e}_k + \lambda \le t\,\Delta_f\},\qquad E_t = \mathrm{Enc}\big(S_{\kappa(t)}, \{S_k\}_{t^{e}_k + \lambda \in ((t-1)\Delta_f,\, t\Delta_f]}\big).
$$

- At 10 Hz versus 12.5 Hz, the hyperperiod is 400 ms (4 ego frames, 5 speech frames). Some speech frames receive no new snapshot (zero-order hold) and some receive one.
- The information age is bounded: $\text{age}(E_t)\in[\lambda, \lambda + T_e]$, with $T_e$ the perception period. At 20 Hz some frames receive two snapshots, so pool them with attention or take the latest.
- Train with randomized $\lambda$ and jitter so the model is robust to scheduler variance on the vehicle.
- Include $\text{age}$ as an explicit feature. A model that "knows" its information is 180 ms stale can hedge ("looks like that car is merging").

**Tokenizing the semantic-geometric state.** Represent $S_k$ as a set of $K$ object slots plus global tokens. For object $j$ in the ego frame:

$$
o_j = \big[\text{class}_j,\; x_j, y_j,\; v_j,\; \psi_j,\; \{\hat\tau_j^{(m)}, \pi_j^{(m)}\}_{m=1}^{M},\; \text{TTC}_j,\; \sigma_j,\; \alpha_j\big]
$$

Here $\hat\tau^{(m)}$ are the world model's predicted trajectory modes with probabilities $\pi^{(m)}$, $\sigma_j$ is uncertainty, and $\alpha_j$ is the planner's attribution to object $j$ (see §3.6b). Global tokens carry the planned maneuver (a meta-action such as yield, lane change, or stop, as in Alpamayo 2), ego speed, the route segment, scene criticality $c_t \in [0,1]$, driver workload $\omega_t$, and intent acknowledgments.\[34\]\[35\]

Embed each slot as $e_j = \mathrm{MLP}(o_j) + \mathrm{PE}_{\text{Fourier}}(x_j, y_j)$. Use continuous Fourier features rather than discretizing coordinates into vocabulary tokens, because text-ified coordinates (the GPT-Driver style) waste sequence length at 12.5 Hz.

**Injection: summed embeddings vs. cross-attention.**

- *Summed (Moshi/PersonaPlex style):* $h^{0}_t = \sum_{q}\mathrm{Emb}(u_t^q) + \sum_q \mathrm{Emb}(a_{t-1}^q) + \mathrm{Emb}(M_{t-1}) + W_g\, g_t + W_s\,\mathrm{Pool}(\{e_j\})$. This costs $O(1)$ per frame and keeps sequence length constant, but it compresses a set into one vector. It is suitable for low-dimensional context ($c_t$, $\omega_t$, maneuver, age). PersonaPlex sums 17 per-codebook embeddings into one 4096-d vector per frame, so this is native to the architecture.\[36\]
- *Gated cross-attention (Flamingo style):* in selected layers $\ell$, $h \leftarrow h + \tanh(\beta_\ell)\cdot \mathrm{CrossAttn}(h, \{e_j\}_{\text{window}})$ with $\beta_\ell$ initialized to 0, so the pretrained duplex behavior is preserved at step 0. This costs $O(K)$ per layer per frame and is needed for *grounded reference* ("the cyclist on your right").
- **Recommendation:** use both. Sum the global context into the input embedding, and add cross-attention to object slots in the upper third of layers. DyaPlex reads a frozen PersonaPlex's per-layer hidden states into a separate motion tower.\[36\] That shows the reverse direction (reading the speech model's residual stream for another modality) is also tractable.

#### 3.2 Intent and preference passing under ego-stack authority

**Intent vocabulary.** $I_t \in \{\varnothing\} \cup \{(\text{type}, \text{args}, \text{conf})\}$ with types such as `SPEED_PREF(Δv)`, `LANE_PREF(L|R)`, `STOP_REQUEST(region)`, `ROUTE_CHANGE(poi)`, `STYLE(comfort|assertive)`, and `FOLLOW_DIST(Δ)`. Natural-language preferences of exactly this type are what Li Auto advertises for its VLA ("drive slower," "change to the right lane," remembered speed preferences) and what Wayve's LINGO-2 demonstrated as closed-loop language-adapted driving.\[37\]\[38\]\[39\]\[40\]\[41\]

**Arbitration as bounded cost shaping over a hard-constraint set.** Let $\mathcal{T}$ be candidate trajectories, $g_i(\tau)\le 0$ the hard constraints (collision with predicted agents under uncertainty, road boundary, traffic rules, dynamics), and $\mathcal{T}_{\text{safe}}=\{\tau: g_i(\tau)\le 0\ \forall i\}$. Then:

$$
\tau^{*} = \arg\min_{\tau\in\mathcal{T}_{\text{safe}}(I)} \; J_{\text{drive}}(\tau) + \sum_{i} w_i\, \phi_i(\tau; I_i),\qquad 0 \le w_i \le \bar w_i\,\mathrm{conf}_i\,\rho(c_t)
$$

- $\phi_i$ is a preference penalty, e.g. $\phi_{\text{speed}}=\sum_k \max(0, v_k - (v_{\text{ref}}-\Delta v))^2$.
- $\rho(c_t)$ down-weights preferences as criticality rises.
- $\bar w_i$ is a hard cap, so no preference can dominate the driving objective.

**Monotonicity (authority) rule:** $\mathcal{T}_{\text{safe}}(I) \subseteq \mathcal{T}_{\text{safe}}(\varnothing)$ for all $I$. Intents may add constraints (a lower speed cap, a larger following gap) but never relax them.

**Safety filter.** Downstream of the planner, a CBF-QP enforces forward invariance of the safe set $\mathcal{C}=\{x: h(x)\ge 0\}$:

$$
u^{\text{safe}} = \arg\min_{u\in\mathcal{U}} \|u - u^{\text{nom}}\|^2 \quad \text{s.t.}\quad L_f h(x) + L_g h(x)\,u + \gamma\big(h(x)\big) \ge 0
$$

In discrete time, the constraint is $h(x_{k+1}) \ge (1-\gamma)h(x_k)$ with $\gamma\in(0,1]$.

**Proposition (conditional authority guarantee).** Suppose $h$ is a valid CBF on $\mathcal{C}$, the QP is feasible, and intents enter only through $u^{\text{nom}}$ (via $J$) or through constraint-tightening. Then $x_0\in\mathcal{C}$ implies $x_k\in\mathcal{C}$ for every intent sequence. The conversational model can therefore change *which* safe behavior is chosen, but not *whether* behavior is safe.

The guarantee is only as good as $h$. With learned perception, $h$ depends on perceived geometry, so perception errors break it. That is exactly the SOTIF concern, and it applies to the driving stack regardless of the voice model.

Literature anchors:
- *Semantically parameterized MPC-CBF* maps language only to validated, bounded CBF decay parameters and explicitly excludes language-generated control commands.\[42\] Under the monotonicity rule, only the conservative direction ($\gamma' \le \gamma$) should be allowed.
- *LA-RL* pairs an LLM decision-maker with MPC + discrete CBFs (with a slack mechanism) for highway driving.\[43\]
- *Language-conditioned latent safety filters* (arXiv 2608.00315) learn a single filter conditioned on a language embedding.\[44\]
- *From Words to Safety* (arXiv 2511.05889) grounds instructions into constraints enforced by an MPC safety filter.\[44\]

**Making it work with an end-to-end learned policy (the Helm case).** A learned policy has no explicit $J$ to add to. There are two options.

- **(a) Preference-conditioned policy.** Train $\pi(\tau \mid S, I)$ in semantic simulation with randomized preferences, using rewards $r = r_{\text{drive}} - \sum_i w_i\phi_i$. This is a natural fit for semantic-space simulation, because preference randomization costs nothing when there is no rendering.
- **(b) Sample-and-rescore.** Draw $N$ candidates from the policy and choose
$$\tau^* = \arg\min_{n:\,\tau_n\in\mathcal{T}_{\text{safe}}} \big[-\log\pi(\tau_n\mid S) + \textstyle\sum_i w_i\phi_i(\tau_n;I_i)\big].$$
This leaves the certified policy untouched and is easier to argue in a safety case.

Either way, the planner returns an **acknowledgment token** (accepted / deferred / rejected + reason code, e.g. `REJECT: BLIND_SPOT_OCCUPIED`) into E. That closes the loop so the voice model can say "Not yet, there's a car in your right blind spot," grounded in the planner's own reason.

#### 3.3 Context-conditioned turn-taking

**Action space per frame:** $a_t \in \{\textsf{silent}, \textsf{backchannel}, \textsf{continue}, \textsf{start\_turn}, \textsf{yield}, \textsf{alert}\}$. ASPIRin projects raw token RL onto a binary speak/silent action to protect semantics. DuplexPO keeps the original action space and restricts updates to dynamics-critical windows.\[45\]\[46\]

**State:** conversation state $z_t$, criticality $c_t$, driver workload $\omega_t$ (from DMS gaze, maneuver phase, and speech-derived load), hazard urgency $\eta_t$, and a pending message with value $V_0$ and urgency time-constant $T_m$.

**Objective (constrained):**

$$
\max_{\pi}\ \mathbb{E}\Big[\sum_t r_{\text{conv}}(a_t,z_t) + r_{\text{info}}(a_t)\Big]\quad\text{s.t.}\quad \mathbb{E}\Big[\sum_t \mathbb{1}\{a_t\in\mathcal{A}_{\text{demand}}\}\,L(\omega_t,c_t)\Big]\le d,
$$

- $\mathcal{A}_{\text{demand}}=\{\textsf{start\_turn},\textsf{continue}\}$ (actions that impose listening load or expect a reply).
- $L$ is a load cost, e.g. $L=\omega_t c_t$.
- A Lagrangian relaxation gives a per-frame penalty $\lambda L(\omega_t,c_t)$, with $\lambda$ adapted by dual ascent.

**Hold-off as optimal stopping.** A non-urgent message loses value as $V(t)=V_0 e^{-(t-t_0)/T_m}$. Deliver at the first $t$ with

$$
V(t) - \lambda L(\omega_t,c_t) \ \ge\ \mathbb{E}\big[\max_{t'>t}\big(V(t') - \lambda L(\omega_{t'},c_{t'})\big)\big].
$$

Because the world model predicts $c_{t'}$ (the merge ends in about 4 s), the assistant can *anticipate* the low-load window, which is what a good passenger does.

**Hazard barge-in** is a separate regime. If $\eta_t > \eta^*$, then $\textsf{alert}$ preempts regardless of floor state. For safety-relevant alerts, the alert content should come from the deterministic HMI path. The generative model only yields and, optionally, adds short follow-up speech.

**Content modulation.** Drews et al. (2008) found that passengers referenced traffic and both partners modulated speech complexity with driving difficulty.\[16\]\[47\] Make utterance length and lexical complexity a function of $c_t$: a token budget $B_t = B_{\max}(1-c_t)$ for the current turn, enforced softly via reward.

**RL recipe (DuplexPO-style).** Sample event windows centered on turn transitions, backchannels, barge-ins, *and* criticality onsets in the synthetic ego stream. For $G$ rollouts per window, compute

$$
R = R_{\text{turn}}+R_{\text{bc}}+R_{\text{yield}} \underbrace{-\ \beta\, c_t\,\omega_t\,\mathbb{1}\{\textsf{start\_turn}\}\ -\ \beta_h\,\mathbb{1}\{\eta_t>\eta^*\wedge a_t\neq\textsf{yield}\}}_{R_{\text{ctx}}},\qquad
\hat A_i=\frac{R_i-\mathrm{mean}(R_{1..G})}{\mathrm{std}(R_{1..G})},
$$

and optimize the clipped GRPO objective with a KL penalty to the SFT reference. DuplexPO reports that this decoupling improves backchannels, turn-taking and barge-in handling without degrading reasoning. Its authors also note that factorized rewards "may miss subtle pragmatic factors,"\[14\]\[45\] which matters even more in a car.

**Interaction with L3 takeover.** Takeover requests (TORs) and mode-confusion risks argue that the voice model should never be the channel that *announces* a TOR.\[48\] It should, however, be conditioned on TOR state so it falls silent, and it can later *explain* the TOR from E-stream reason codes.

#### 3.4 Dual-rate cognition

**Two coupled processes:**
- A fast process $\mathcal{F}$ on the 80 ms clock (the duplex model).
- A slow deliberator $\mathcal{D}$ that runs asynchronously on events (a user question, a route deviation, a trip-context update) with latency of 0.5–10 s.

$\mathcal{D}$ produces results $r_k$ at times $t_k$, each tagged with the ego snapshot ID it reasoned over, $s_{k}$. The fast model conditions on the latest result:

$$
p(A_t\mid\cdot) = p_\theta\big(A_t \mid h_t, E_t, r_{k(t)}\big),\quad k(t)=\max\{k: t_k\le t\}.
$$

**Staleness gate.** Before verbalizing $r_k$, check $d(S_{s_k}, S_{\kappa(t)}) < \epsilon$ on the relevant fields (e.g., "the exit is still ahead"). If the check fails, re-query or reframe. GPT-Live implements the delegation pattern commercially: the live model "delegates to our latest frontier model behind the scenes and brings the result back into the conversation when it's ready," while continuing to talk.\[15\] DuplexSLA emits planning text and tool calls on its action channel without halting speech.\[49\] A recent frontend-backend architecture (arXiv 2609.19334) uses a delegation token from a duplex model to hand transcripts to a backend LLM.\[50\]

**Mapping to driving dual systems.** DriveVLM and Li Auto's earlier "end-to-end + VLM" dual system use the same fast/slow split for *driving*.\[51\]\[52\] In the KITT design, the slow conversational deliberator must never sit in the driving loop. It reads E, never writes to the planner except through I, and runs at lowest scheduling priority.

#### 3.5 Self-supervised alignment between speech and driving state

**Peri-event analysis** (directly from single-unit practice):
- Define driving events $\{e_j\}$: hard braking ($a_x<-0.4g$), TTC < 1.5 s, merge onset, lane change.
- Take any speech feature $x(t)$: voice-activity probability, energy, f0, a laughter/gasp detector score, or a projection of Mimi latents.
- Compute the peri-event time histogram $\bar x(\tau)=\frac1J\sum_j x(e_j+\tau)$ for $\tau\in[-T,T]$.
- Build a null by jittering event times within the same drive by ±U(5–30 s). This preserves route and time-of-day confounds while breaking time-locking.
- Use wavelet coherence between $c(t)$ and speech activity to localize the coupling in time and frequency. Expect silence during merges (suppressed speech onset) and paralinguistic spikes after hard braking.

These are free labels for workload and arousal.

**Contrastive alignment** (InfoNCE with lag and hard negatives):

$$
\mathcal{L}_{\text{NCE}} = -\mathbb{E}\left[\log\frac{\exp(\mathrm{sim}(f(x_{t-W:t}), g(S_{t+\delta}))/T)}{\sum_{n\in\mathcal{N}}\exp(\mathrm{sim}(f(x_{t-W:t}), g(S_n))/T)}\right],\qquad I(x;S)\ \ge\ \log|\mathcal{N}| - \mathcal{L}_{\text{NCE}}.
$$

Draw negatives from the *same drive* at distant times, so the model cannot shortcut on speaker identity, cabin acoustics or route. Sweep $\delta$ to estimate the lag structure (speech often *follows* events by 0.3–1.5 s). Be honest about effect size: most conversation is not about driving, so the mutual information is small and concentrated around events. Event-window sampling, the same trick DuplexPO uses for timing, is the right sampler.

**JEPA-style predictive objective.** Predict future ego latents from speech plus past ego, with a stop-gradient EMA target:

$$
\mathcal{L}_{\text{JEPA}} = \big\| P_\phi\big(f(x_{\le t}), z^{e}_{\le t}\big) - \mathrm{sg}\big(\bar z^{e}_{t+\Delta}\big)\big\|_2^2.
$$

The improvement over an ego-only predictor measures how much speech carries predictive information about the scene, for instance a passenger saying "watch out" before perception has confirmed the hazard.

**Semantic-space simulation for dialogue.** This is the strongest argument for the factored design. Because the voice model consumes the *same E tokens* in simulation and on the road, the sim-to-real gap for the conversational grounding problem lies only in the distribution of E, not in pixels. Audio realism comes from real conversational pretraining, which is what Dessa and Moshi supply.

Pipeline:
1. Generate scenario timelines in semantic simulation.
2. Generate time-stamped dialogue scripts conditioned on those timelines: questions, explanations, hold-offs, alerts.
3. Render the scripts as two-channel audio with TTS, or with duplex self-play.
4. Train.

This mirrors Helm's own claim that semantic-space training "works immediately in the real world" for driving.\[1\] That claim is unverified for driving, and the conversational analog is untested.

#### 3.6 Representation geometry

**(a) Interface alignment and drift monitoring.** Fix a probe set of $n$ scenes. Let $X\in\mathbb{R}^{n\times d}$ be the voice model's layer-$\ell$ activations for those scenes under model version $v$, and $Y$ the same under version $v'$ (or the ego encoder's embeddings before and after a perception update).

- *Principal angles.* With orthonormal bases $Q_X, Q_Y\in \mathrm{St}(d,k)$ of the top-$k$ PCA subspaces, take the SVD $Q_X^\top Q_Y = P\,\mathrm{diag}(\cos\theta_1,\dots,\cos\theta_k)\,R^\top$. The Grassmann geodesic distance is $d_G=\|\theta\|_2$ and the projection distance is $\|\sin\theta\|_2$.
- *Orthogonal Procrustes.* $R^*=\arg\min_{R\in O(d)}\|XR-Y\|_F = UV^\top$, where $X^\top Y=U\Sigma V^\top$. The residual $\|XR^*-Y\|_F/\|Y\|_F$ measures non-rotational change.
- *CCA / SVCCA / PWCCA.* CCA maximizes $\mathrm{corr}(Xa, Yb)$. SVCCA applies CCA after SVD truncation. PWCCA weights canonical correlations by how much of the representation each direction accounts for.
- *Linear CKA.* For centered $X$ and $Y$: $\mathrm{CKA}(X,Y)=\dfrac{\|Y^\top X\|_F^2}{\|X^\top X\|_F\,\|Y^\top Y\|_F}$. It is invariant to orthogonal transforms and isotropic scaling.

**Regression-test contract.** Run these checks when either model updates:
1. Schema check: typed fields and units.
2. Ego-embedding CKA ≥ threshold on the probe set.
3. Principal angles of the voice model's *E-reading subspace* below a threshold. Use the row space of the cross-attention $W_K$ and $W_V$ restricted to E, or the top singular vectors of $\partial h_\ell/\partial E$.
4. Behavioral grounded-QA accuracy and explanation-faithfulness scores (b) unchanged within confidence intervals.

Geometric drift without behavioral drift is a warning. Behavioral drift without geometric drift points to the decoder.

**(b) Faithfulness of verbal explanations.** Let explanation $\varepsilon$ cite a set of objects $C$ and a reason $\rho$ ("slowing for the pedestrian on the right").

1. *Planner-side causal relevance.* For each cited object $j\in C$, $\Delta_j = D\big(\pi(\cdot\mid S),\ \pi(\cdot\mid \mathrm{do}(S\setminus o_j))\big)$ must exceed a threshold, and non-cited objects should have smaller $\Delta$. In a factored stack this intervention is trivial: delete the slot in semantic space. That is a major advantage over pixel-space VLAs.
2. *Readout consistency.* Under the same intervention, the voice model's explanation must change accordingly: the explanation logit for $\rho$ drops when $o_j$ is removed. Activation patching swaps the E-slot of object $j$ from the counterfactual run and measures the change in explanation logits. *Directional ablation* projects out a direction $v$ from $h_\ell$ ($h\leftarrow h - vv^\top h$) and measures the same.
3. *Subspace overlap.* Let $\mathcal{S}_\pi$ be the planner's decision-relevant subspace (top singular vectors of $\partial\tau/\partial z$ over scenes) and $\mathcal{S}_{\text{read}}$ the subspace through which the voice model reads E. Small principal angles between their images in a shared coordinate system (the E-embedding space) indicate the explainer is looking where the planner looks. This is necessary but not sufficient, so pair it with 1–2.

**Design rule: explanation by readout, not rationalization.** Pass the planner's attribution vector $\alpha$ (ablation deltas $\Delta_j$, which are cheap in semantic space) and active-constraint IDs in E. Then require $C \subseteq \mathrm{top}\text{-}k(\alpha)$ as a checkable runtime invariant.

The empirical motivation is strong:
- A CVPR 2026 DriveX workshop study of 300 Alpamayo-R1-10B inferences across 100 scenarios reported 42.5% overall fidelity, 94 missed pedestrians, 53.3% low reasoning-action consistency, 97.7% perturbation fragility and 14.3% silent failures. It concluded that chain-of-causation reasoning "cannot be relied upon for safety assurance."\[53\]
- A counterfactual reasoning-action coupling study (arXiv 2610.30953) concludes that "plausible, visually responsive reasoning is not by itself a faithful explanation."\[54\]
- Mitigations exist: rule-grounded traces from classical planners (Neuro-Symbolic Drive) and faithfulness critics for RL post-training (Pinocchio, which reports improving faithfulness by 4% and 18% over alignment and trajectory-error baselines).\[55\]\[56\]\[57\]

**(c) Off-manifold detection for driver state and scenes.**
- *Scenes:* compute the Mahalanobis distance $d_M(z)=\sqrt{(z-\mu)^\top\Sigma^{-1}(z-\mu)}$ in the ego-embedding space, or use kNN distance or the log-density of a normalizing flow.
- *Driver state:* voice is idiosyncratic, so use per-driver baselines $(\mu_d,\Sigma_d)$ estimated over the first $N$ drives, plus a population prior with shrinkage $\Sigma_d \leftarrow (1-s)\hat\Sigma_d + s\Sigma_{\text{pop}}$.
- *Fusion:* fuse voice with DMS telemetry by summing calibrated log-likelihood ratios.
- *Calibration:* calibrate the alarm threshold conformally. With scores on $n$ in-distribution calibration windows, set $q=$ the $\lceil (1-\alpha)(n+1)\rceil$-th smallest score, which guarantees a false-alarm rate ≤ α under exchangeability.
- *Regulatory limits:* the EU mandates are defined differently. ADDW is gaze-based (a warning after 6 s of distraction at ≥ 20 km/h). Industry summaries tie DDAW to the Karolinska Sleepiness Scale (KSS) level 8.\[18\]\[58\] Voice paralinguistics can therefore only *supplement* compliant systems. Voice-based impairment detection is not a validated regulatory measure.

**(d) Equivariance.**
- *Pitch.* A pitch shift by factor $c$ is a translation by $\log c$ on a log-frequency axis (CQT), so convolution along that axis is translation-equivariant. Use it for speaker-normalized paralinguistic features, so that distress detection is not confounded by voice register.
- *Scene geometry.* The ego tokenizer should canonicalize to the ego frame, which gives invariance to global SE(2) pose. Note that *reflection* is not a symmetry of driving, because drive-on-right versus drive-on-left breaks it. Mirror augmentation, flipping the scene and swapping "left" and "right" in the dialogue, is valid only together with a handedness token. Otherwise it teaches wrong traffic semantics.
- *Equivariant networks.* For internal E encoders, SE(2)-equivariant message passing over object slots satisfies $f(g\cdot S)=\rho(g)f(S)$. Language about direction is ego-relative, so feed the voice model canonicalized (invariant) features plus explicit ego-relative direction tokens.

### 4. Engineering requirements

**Compute and scheduling.**
- NVIDIA positions DRIVE Thor as a centralized computer able to run driving, cabin AI and infotainment with "multi-compute domain isolation," using MIG-derived GPU partitioning so time-critical processes are not interrupted. NVIDIA materials cite both 1,000 and 2,000 TFLOPS depending on configuration and precision,\[20\]\[21\]\[59\] so treat the headline figure as precision-dependent marketing.
- Allocate a fixed MIG slice to the duplex model.
- The real-time requirement is a real-time factor < 1: one 80 ms frame (temporal transformer plus depth transformer plus codec decode) must compute in well under 80 ms, with headroom for the cross-attention over E.
- The deliberator runs preemptibly at lowest priority.

**Latency budget** (targets, not measured):

| Path | Budget |
|---|---|
| Perception → E available | λ ≈ 50–150 ms (stack-dependent) |
| Speech frame compute | < 60 ms of the 80 ms frame |
| User end-of-turn → assistant onset | ~200 ms (Moshi's practical figure); DuplexPO reports 0.24 s on Full-Duplex-Bench v3\[12\]\[60\] |
| Hazard onset → generative model yields floor | ≤ 1 frame after E flag (80–160 ms) |
| Safety alert | Deterministic HMI path, independent of the voice model |

**Functional safety partitioning.**
- The voice model is QM (non-safety-rated).
- Freedom from interference is required in space (memory and MIG isolation), time (bounded GPU and memory-bandwidth contention, a known weak point with shared memory controllers) and communication (one-way E publication, validated I).
- The intent validator/arbiter (schema, range and plausibility checks, rate limiting, source-zone authorization) belongs on the safety side.
- SOTIF hazards specific to voice:
  - overtrust caused by fluent but unfaithful explanations;
  - mode confusion ("I asked it to slow down, so it is handling it");
  - distraction induced by badly timed speech.
  Each needs a triggering-condition analysis.

**Interface versioning.** Treat E as a versioned API: semantic-versioned schema, units, coordinate frames, and a probe-set regression suite (§3.6a). Any perception update that shifts E-embeddings beyond tolerance triggers re-validation of the voice adapter, exactly as a sensor-calibration change would.

**Cabin acoustics.**
- Use multi-zone microphone arrays with beamforming for per-seat streams.
- Acoustic echo cancellation needs reference signals for *both* the assistant's own voice and media playback.
- Speaker diarization and zone attribution map directly to authorization: only the driver zone, or a configured authorized zone, may issue driving intents.
- Road noise varies with speed and surface. AISHELL-5 (arXiv 2505.23036) provides "over 100 hours of multi-channel speech data recorded in an electric vehicle across more than 60 real driving scenarios," captured by far-field microphones on each car door, plus 40 hours of real-world environmental noise recordings for development. Self-listening matters too: recent duplex work ("What Did I Just Say?") targets the model tracking its own output.\[61\]

**Security.** Auditory prompt injection is demonstrated, not hypothetical:
- AudioHijack (Chen et al., arXiv 2604.14604) reports "average success rates of 79%-96% on unseen user contexts" across 13 large audio-language models, and showed that "commercial voice agents from Mistral AI and Microsoft Azure can be induced to execute unauthorized actions."
- Concurrent injections can overlap ongoing user speech.\[62\]

Mitigations:
- Anything played through the car's own speakers (radio, media, the assistant itself) is known via the AEC reference and is tagged non-authoritative.
- Intents require zone-localized origin plus semantic plausibility against E.
- No audio-derived intent bypasses the arbiter.
- Intents are rate-limited.

**Privacy.**
- Process on device by default. Dessa's local real-time claim and Moshi's on-device lineage support this.\[7\]
- ADDW must work without biometric identification.\[18\] Voiceprint-based authorization would be biometric processing under GDPR and needs a separate legal basis, so prefer zone-based authorization.

### 5. Data

| Corpus | What it pairs | Fit |
|---|---|---|
| CIAIR in-car corpus | Synchronized 16-mic audio, 3 cameras, vehicle signals; about 60-min dialogue drives by each of 200 subjects\[63\] | Closest to "speech + vehicle state" but prompted dialogue, older |
| AMIE (Intel) | ~30 h, 30 passengers, Wizard-of-Oz in-cabin commands (incl. GoSlower intents), multimodal\[64\] | Intent grounding, not natural duplex conversation |
| AISHELL-5 | 100+ h in-car multi-channel natural conversation + noise\[65\] | Acoustics and diarization; no ego state |
| Commentary driving (2204.09109), NARRATE (2608.14767) | Think-aloud driver narration time-aligned with sensors\[66\]\[67\] | Explanation grounding and peri-event studies |
| BDD-X, DriveLM, nuScenes-QA, Talk2Car | Text explanations, QA, commands over driving scenes | Grounded language; text, not speech |
| NVIDIA PhysicalAI + CoC labels | Large multi-sensor driving dataset with chain-of-causation reasoning\[68\] | Reasoning supervision; text |
| Fisher / two-channel conversational corpora | Natural full-duplex conversation | Base duplex dynamics (used by DuplexPO)\[69\] |

**Gap:** no public corpus pairs *natural full-duplex conversation* with *full ego state*. The practical route is synthesis: overlay real ego logs or semantic-simulation timelines on real two-channel conversation and edit the timing (§6), plus small instrumented collections for validation.

### 6. A feasible single-H200 experiment

**Question:** can a duplex model learn *context-conditioned* turn-taking from an aligned ego stream, and is the effect causal and mechanistically localized?

**Base model.** Moshi (7B, CC-BY weights) or PersonaPlex, LoRA-fine-tuned. Weights in bf16 take about 15 GB, which fits comfortably on a 141 GB H200 along with activations for multi-minute contexts at 12.5 Hz. A 5-minute clip is about 3,750 frames. A smaller from-scratch dGSLM-style model is a fallback for full ablations.

**Data construction.**
1. *Ego timelines.* Generate them from a lightweight 2D semantic simulator (highway merges, intersections, sudden braking). Per frame, compute $c_t$ (from TTC and maneuver phase), $\omega_t$, maneuver one-hot, and a hazard flag $\eta_t$.
2. *Conversations.* Take 200–500 h of two-channel conversation and assign one channel as the "assistant."
3. *Timing edits conditioned on the timeline:*
   - (i) Assistant turn onsets that fall in high-$c_t$ windows are delayed to the next low-$c_t$ window by time-stretching silence, producing a hold-off target.
   - (ii) At hazard onsets, a short alert utterance (TTS, in the assistant voice) is inserted within 1–2 frames, and any ongoing assistant speech is truncated.
   - (iii) Under high $c_t$, the length of assistant utterances is reduced by selecting shorter turns from the pool.
4. *Controls:* identical audio with temporally shuffled E (aligned versus misaligned), so a model can only succeed by using alignment, not the marginals.

**Training.**
- Stage 1: SFT with the E stream injected as a summed embedding $W_g g_t$, a new zero-initialized projection, so the step-0 model equals base Moshi.
- Stage 2 (optional): DuplexPO-style event-window GRPO with $R_{\text{ctx}}$ (§3.3). The public DuplexPO reference implementation of event construction and event-level GRPO is a starting point.\[70\]

**Metrics.**

| Metric | Definition | Success criterion |
|---|---|---|
| Onset hazard ratio | Rate of assistant turn onsets in high-$c$ windows ÷ rate in low-$c$ windows | ≪ 1 for aligned E; ≈ base for shuffled E |
| Alert latency | Hazard onset → assistant yield + alert onset | ≤ 160 ms median |
| Hold-off cost | Delay of non-urgent content vs. the next low-$c$ window | Close to the oracle optimal-stopping delay |
| Complexity modulation | Syllables/word and turn length vs. $c_t$ (slope) | Negative slope, mirroring Drews et al. |
| Base dynamics regression | Full-Duplex-Bench v1/v1.5 turn-taking, backchannel, interruption | No significant degradation |
| Content regression | Spoken QA (e.g., Llama Questions / WebQ) | No significant degradation |

**Mechanistic analysis** (using your Grassmann/Stiefel tooling):
1. Linear probes for $c_t$ across layers.
2. Principal angles between the $c_t$-probe subspace and the subspace that separates speak/silent decisions.
3. Activation patching: swap $W_g g_t$ between high-$c$ and low-$c$ runs and measure the change in onset probability.
4. Directional ablation of the $c$-direction, which should restore base timing.
5. CKA between base and fine-tuned residual streams on non-driving audio, to quantify collateral drift.

**Experiment B (no training).** Run peri-event and wavelet-coherence analysis of speech activity and paralinguistics around braking, merge and TTC events on any instrumented corpus you can access (CIAIR-style, NARRATE-style, or self-collected). This establishes the natural coupling that the synthetic edits assume.

The compute estimates are engineering judgment, not measured: LoRA SFT over ~300 h of audio should take on the order of days on one H200.

---

## Open Problems

1. **Faithfulness under learned policies.** Semantic-space ablation makes planner-side attribution cheap, but attribution methods for transformer planners still disagree. "Top-k(α) contains cited objects" is necessary, not sufficient.
2. **Authority guarantees inherit perception error.** The CBF proposition holds in the *perceived* state. No voice-layer design fixes a missed pedestrian. The voice layer should instead never *increase* the trust the driver places in perception.
3. **Speech-only cognition vs. grounding.** Dessa-style speech-only training yields natural prosody and timing. Whether such a model can be grounded in structured E tokens without text supervision is unknown. Moshi's text monologue and DuplexSLA's text action channel suggest text remains the practical bridge.
4. **Evaluation of context-conditioned timing.** No benchmark exists. Full-Duplex-Bench measures context-free dynamics, and DuplexSLA-Bench measures tool timing.
5. **Workload estimation from voice** is plausible from peri-event effects but unvalidated against regulatory measures (gaze-based ADDW, KSS-anchored DDAW).
6. **Helm-specific unknowns.** Deep Teaching's mechanism, the planner's internal representation, and whether its interface exposes uncertainty and attributions are all unpublished.
7. **Interference on shared SoCs.** MIG partitions compute, but memory-bandwidth and thermal coupling between a 7B voice model and the driving stack needs empirical worst-case analysis.

## Recommendations (prioritized)

1. **Commit to the factored coupling.** Use E as a read-only, versioned, rate-aligned stream, I as a typed, bounded, arbitrated preference channel with the monotonicity rule, and deterministic safety HMI. Reject end-to-end "drive-and-talk" models for anything beyond research.
2. **Run Experiment B first** (peri-event and wavelet coherence). It is cheap, plays to your neurophysiology strengths, and validates the core assumption that conversation is time-locked to driving demand.
3. **Then run the single-H200 experiment** (§6) with shuffled-E controls and mechanistic patching. A clean causal result that "ego-stream alignment causally gates onset" would be a publishable contribution, since no benchmark or paper covers context-conditioned duplex timing.
4. **Build the faithfulness harness alongside the model.** Semantic-slot ablation, readout consistency and subspace overlap should become regression tests from day one.
5. **Track Dessa and Helm.ai disclosures.** Re-evaluate if Block releases weights or a paper, or if Helm.ai publishes its policy interface. Until then, treat both companies' claims as vendor claims.

## Reading list (prioritized)

1. Moshi: speech-text foundation model for real-time dialogue — https://arxiv.org/abs/2410.00037 (multi-stream, Mimi, inner monologue)
2. DuplexSLA: synchronized speech, language and action — https://arxiv.org/abs/2605.20755 (action channel on a shared clock; template for E/I streams)
3. DuplexPO: decoupling conversational dynamics via RL — https://arxiv.org/abs/2607.07148 (event-window GRPO for timing)
4. Is VLA Reasoning Faithful? (Alpamayo-R1) — https://arxiv.org/abs/2605.17268 ; Reasoning–Action Coupling via Counterfactuals — https://arxiv.org/html/2610.30953
5. Helm.ai, Factored Embodied AI (vendor) — https://helm.ai/post/factored-embodied-ai ; Helm.ai Driver release (vendor) — https://helm.ai/post/helmai-driver
6. Lychee-FD: hierarchical acoustic-semantic modeling — https://arxiv.org/abs/2607.06540
7. ASPIRin: action-space projection RL for duplex timing — https://arxiv.org/abs/2604.10065
8. OpenAI, Introducing GPT-Live (delegation) — https://openai.com/index/introducing-gpt-live/
9. NemotronLabs VoiceChat (parallel streams, tool calling) — https://arxiv.org/abs/2609.21967 ; Frontend–backend tool calls in duplex models — https://arxiv.org/html/2609.19334v1
10. DyaPlex (reading PersonaPlex hidden states into another modality) — https://arxiv.org/abs/2606.03874
11. Semantically parameterized MPC-CBF with bounded LLM parameterization — https://link.springer.com/article/10.1007/s44443-026-01234-z ; LA-RL — https://arxiv.org/abs/2512.05686 ; Language-conditioned latent safety filters — https://arxiv.org/abs/2608.00315 ; From Words to Safety — https://arxiv.org/abs/2511.05889
12. Do VLAs Mean What They Say? (Pinocchio) — https://arxiv.org/abs/2607.04681 ; Neuro-Symbolic Drive — https://arxiv.org/html/2606.23938
13. Drews, Pasupathi & Strayer (2008), Passenger and cell phone conversations in simulated driving, J. Exp. Psych.: Applied 14(4) — https://doi.org/10.1037/a0013119 \[17\]
14. Wayve LINGO-2 — https://wayve.ai/thinking/lingo-2-driving-with-language/ ; SimLingo — https://arxiv.org/abs/2503.09594 \[38\]\[71\] ; NVIDIA Alpamayo 2 — https://huggingface.co/blog/nvidia/nvidia-alpamayo-2
15. AudioHijack: auditory prompt injection — https://arxiv.org/abs/2604.14604 ; AudioJailbreak — https://arxiv.org/abs/2505.14103
16. AISHELL-5 in-car speech — https://arxiv.org/abs/2505.23036 ; AMIE in-cabin intents — https://arxiv.org/abs/1904.10500 ; NARRATE — https://arxiv.org/html/2608.14767
17. Foundations: PhaseLift — https://arxiv.org/abs/1109.4499 ; Hand & Voroninski, deep generative priors — https://arxiv.org/abs/1705.07576 ; CKA (Kornblith et al.) — https://arxiv.org/abs/1905.00414 ; SVCCA — https://arxiv.org/abs/1706.05806 ; CBF theory (Ames et al.) — https://arxiv.org/abs/1903.11199 ; CPC/InfoNCE — https://arxiv.org/abs/1807.03748 ; I-JEPA — https://arxiv.org/abs/2301.08243 ; dGSLM — https://arxiv.org/abs/2203.16502

## Caveats

- **Vendor claims:**
  - Helm.ai: 1,000 hours, zero-shot steering, the Deep Teaching mechanism, ISO 26262-certifiability.
  - Block/Dessa: speech-only training, realism, local real-time performance.
  - OpenAI GPT-Live: benchmark gains.
  - NVIDIA: Thor TFLOPS.
  None has been independently verified.
- **Conflicting figures:** NVIDIA describes Alpamayo 2 Super as both 32B (press release) and 34B (developer blog / model card) parameters.\[34\]\[68\]\[72\]
- **Secondary sourcing:** some DuplexPO benchmark figures (0.24 s latency, 100% yield on Full-Duplex-Bench v3) come from secondary summaries of the paper.\[60\]
- **Unverified models:** Diffuplex was not verified in this research, and PersonaPlex details come from citing papers rather than its own paper.
- **Original proposals:** the proposition in §3.2, the rate-alignment scheme, the turn-taking objective, the faithfulness criteria and the experiment design are engineering formalizations proposed here, not results from the literature. Their parameters (thresholds, budgets, compute times) are starting points to be tuned empirically.
- **Regulatory summaries:** the descriptions of ADDW/DDAW come from industry summaries of Regulation (EU) 2019/2144 and Delegated Regulation (EU) 2023/2590.\[73\] Consult the legal texts before relying on specific thresholds.

## Sources

1. [Helm.ai releases new architectural framework for autonomous vehicles - The Robot Report](https://www.therobotreport.com/helm-ai-releases-new-architectural-framework-for-autonomous-vehicles/)
2. [Helm.ai Driver Achieves Vision-Only Urban Autonomy, Unlocking Scalability from Level 2+ through Level 4](https://www.businesswire.com/news/home/20260225868470/en/Helm.ai-Driver-Achieves-Vision-Only-Urban-Autonomy-Unlocking-Scalability-from-Level-2-through-Level-4)
3. [Factored Embodied AI: Breaking the Data Wall with Geometric Reasoning and Semantic Simulation](https://helm.ai/post/factored-embodied-ai)
4. [Helm.ai Driver Achieves Vision-Only Urban Autonomy, Unlocking Scalability from Level 2+ through Level 4](https://finance.yahoo.com/news/helm-ai-driver-achieves-vision-140000710.html)
5. [13459: SAE TOMORROW TODAY: Unsupervised Learning Makes Scalable Autonomy a Reality - Podcast](https://saemobilus.sae.org/media/sae-tomorrow-today-unsupervised-learning-makes-scalable-autonomy-a-reality-13459)
6. [Patents Assigned to Helm.ai, Inc. - Justia Patents Search](https://patents.justia.com/assignee/helm-ai-inc)
7. [Block introduces Dessa, a speech-to-speech AI research preview · Digg](https://digg.com/ai/qq1uwed6)
8. [Jack Dorsey’s Block launches Dessa, an AI voice model trained on human conversation - Indian Television Dot Com](https://indiantelevision.com/mam/jack-dorseys-block-launches-dessa-an-ai-voice-model-trained-on-human-conversation/)
9. [Block Dessa: full-duplex speech-to-speech model trained only on human conversation](https://www.datastudios.org/post/block-dessa-full-duplex-speech-model-trained-only-on-human-conversation)
10. [Block launches Dessa, a full-duplex speech model trained only on audio](https://www.voiceaispace.com/articles/block-launches-dessa-a-full-duplex-speech-model-trained-only-on-audio)
11. [Dessa Full-Duplex Speech Model: 0 Text Tokens \[2026\]](https://shattered.io/block-dessa-full-duplex-speech-model-no-text-2026/)
12. [Moshi: a speech-text foundation model for real-time dialogue](https://arxiv.org/abs/2410.00037)
13. [DuplexSLA: A Full-Duplex Spoken Language Model with Synchronized Speech, Language, and Action](https://arxiv.org/html/2605.20755v1)
14. [Decoupling Conversational Dynamics in Full-Duplex Spoken Models through Reinforcement Learning](https://www.alphaxiv.org/abs/2607.07148)
15. [Introducing GPT-Live](https://openai.com/index/introducing-gpt-live/)
16. [\[PDF\] Passenger and cell phone conversations in simulated driving.](https://www.semanticscholar.org/paper/Passenger-and-cell-phone-conversations-in-simulated-Drews-Pasupathi/def7c6f86c476a5f0f1dbf19cd715231389a4438)
17. [Passenger and Cell Phone Conversa... preview & related info](https://www.mendeley.com/catalogue/59430960-cb14-3955-be56-dc7ad471a2ec/)
18. [ADDW Systems Explained: EU's Driver Distraction Mandate - Captain Compliance](https://captaincompliance.com/education/addw-systems-explained-eus-driver-distraction-mandate/)
19. [EU mandates driver distraction warning tech in cars](https://www.traffictechnologytoday.com/news/enforcement/eu-mandates-driver-distraction-warning-tech-in-cars.html)
20. [NVIDIA Unveils DRIVE Thor — Centralized Car Computer Unifying Cluster, Infotainment, Automated Driving, and Parking in a Single, Cost-Saving System](https://nvidianews.nvidia.com/news/nvidia-unveils-drive-thor-centralized-car-computer-unifying-cluster-infotainment-automated-driving-and-parking-in-a-single-cost-saving-system)
21. [NVIDIA DRIVE Thor Strikes AI Performance Balance, Uniting AV and Cockpit on a Single Computer](https://blogs.nvidia.com/blog/drive-thor/)
22. [Vladislav Voroninski - Google Scholar](https://scholar.google.com/citations?user=P8NRbbYAAAAJ&hl=en)
23. <https://api.crossref.org/works?query.bibliographic=PhaseLift+exact+and+stable+signal+recovery+Candes+Strohmer+Voroninski&rows=2>
24. [Helm.ai Driver achieves vision-only urban autonomy - Robotics 24/7](https://www.robotics247.com/article/helm.ai-driver-achieves-vision-only-urban-autonomy/technologies)
25. [Helm.ai unveils low-data approach to autonomous driving](https://www.automotiveworld.com/news/helm-ai-unveils-low-data-approach-to-autonomous-driving/)
26. [Block Unveils Dessa, a Full-Duplex Speech Model Trained on Conversation](https://www.unite.ai/block-unveils-dessa-a-full-duplex-speech-model-trained-on-conversation/)
27. [NemotronLabs VoiceChat: An Open Full-duplex Speech-to-Speech Model with Tool Calling Capabilities](https://arxiv.org/html/2609.21967)
28. [HIT-TMG/Lychee-FD · Hugging Face](https://huggingface.co/HIT-TMG/Lychee-FD)
29. [Yunxin Li on X: "🚀ACL 2026 Outstanding Paper 🏆 Lychee-FD, a fully open-source native end-to-end full-duplex spoken large language model. Built on our language-centric paradigm: split text & speech into dual core carriers • 📝 Text: logical reasoning, knowledge, multi-modal fusion backbone • 🎙️ … / X](https://x.com/LyxTg/status/2077645299234046382)
30. [Global Guarantees for Enforcing Deep Generative Priors by Empirical Risk](https://arxiv.org/abs/1705.07576)
31. [Deep Teaching from Helm.ai designed to ease autonomous vehicle training](https://www.therobotreport.com/deep-teaching-helm-ai-applies-level-2-4-autonomous-vehicles/)
32. [Phase Retrieval Under a Generative Prior](https://arxiv.org/abs/1807.04261)
33. [Lychee Hierarchical Acoustic-Semantic Modeling: Modality](https://arxiv.org/pdf/2607.06540)
34. [NVIDIA Launches Alpamayo 2 Super Open Reasoning Model for Robotaxis](https://www.taiwannews.com.tw/news/6373726)
35. [Taking Alpamayo to New Heights with Driving Foundation Models and Closed-Loop Training](https://huggingface.co/blog/nvidia/nvidia-alpamayo-2)
36. [DyaPlex: Full-Duplex Speech-Motion Model for Dyadic Interaction](https://arxiv.org/pdf/2606.03874)
37. [Wayve launches multimodal driving model Lingo-2](https://www.autonomousvehicleinternational.com/news/ai-sensor-fusion/wayve-launches-multimodal-driving-model-lingo-2.html)
38. [Wayve on X: "At Wayve, we recognize the power of language in robot interaction 🚙💬 especially for building trust in self-driving technology. This drives our ongoing exploration of breakthrough uses of #LLMs in autonomous driving systems. Last September, we debuted LINGO-1, the first" / X](https://x.com/wayve_ai/status/1780567719655469092)
39. [Wayve debuts multimodal driving model Lingo-2](https://icavcluster.co.uk/wayve-debuts-multimodal-driving-model-lingo-2/)
40. [LI Auto i8 Launches: Is VLA the Future of Driving? - 42HOW](https://en.42how.com/2025/07/30/li-auto-i8-launches-is-vla-the-future-of-driving/)
41. [LINGO-2: Driving with Natural Language - Wayve](https://wayve.ai/thinking/lingo-2-driving-with-language/)
42. [Semantically parameterized MPC-CBF control for safety-critical obstacle avoidance in autonomous vehicles](https://link.springer.com/article/10.1007/s44443-026-01234-z)
43. [LA-RL: Language Action-guided Reinforcement Learning with Safety Guarantees for Autonomous Highway Driving](https://arxiv.org/html/2512.05686)
44. [Towards General Language-Conditioned Latent Safety Filters](https://arxiv.org/pdf/2608.00315)
45. [Decoupling Conversational Dynamics in Full-Duplex Spoken Models through Reinforcement Learning](https://arxiv.org/html/2607.07148v1)
46. [Decoupling Conversational Dynamics in Full-Duplex Spoken Models through Reinforcement Learning](https://arxiv.org/pdf/2607.07148)
47. [Drivers Distracted More by Cell Phones Than by Passengers](https://archive.unews.utah.edu/news_releases/drivers-distracted-more-by-cell-phones-than-by-passengers/)
48. [Driver-Net: Multi-Camera Fusion for Assessing Driver Take-Over Readiness in Automated Vehicles](https://arxiv.org/pdf/2507.04139)
49. [GitHub - hyzhang24/DuplexSLA: DuplexSLA: A Full-Duplex Spoken Language Model with Synchronized Speech, Language, and Action · GitHub](https://github.com/hyzhang24/DuplexSLA)
50. [A FRONTEND-BACKEND ARCHITECTURE FOR TOOL CALLS IN FULL-DUPLEX SPEECH MODELS](https://arxiv.org/html/2609.19334v1)
51. [Li Auto unveils next-gen autonomous driving architecture MindVLA](https://autonews.gasgoo.com/articles/news/li-auto-unveils-next-gen-autonomous-driving-architecture-mindvla-70036319)
52. [Li Auto: Advancing competitiveness in AI, future of mobility](https://subsites.chinadaily.com.cn/wic/2026-02/26/c_1162494.htm)
53. [Is VLA Reasoning Faithful? Probing Safety of Chain-of-Causation in Autonomous Driving Models](https://arxiv.org/html/2605.17268v1)
54. [Measuring the Reasoning–Action Coupling in DrivingVision–Language–Action Models through Counterfactuals](https://arxiv.org/html/2610.30953)
55. [Neuro-Symbolic Drive: Rule-Grounded Faithful Reasoning for Driving VLAs](https://arxiv.org/html/2606.23938)
56. [Do Vision-Language-Action Models Mean What They Say? On the Role of Faithfulness in Embodied Reasoning](https://arxiv.org/pdf/2607.04681)
57. [\[2607.04681\] Do Vision-Language-Action Models Mean What They Say? On the Role of Faithfulness in Embodied Reasoning](https://arxiv.org/abs/2607.04681)
58. [A Brief Analysis of EU DDAW and ADDW Regulations](https://www.ulvisions.com/blog/eu-ddaw-and-addw-regulations.html)
59. [Nvidia's new automotive superchip covers all SAE levels, meets FuSa standards](https://www.eenewseurope.com/en/nvidias-new-automotive-superchip-covers-all-sae-levels-meets-fusa-standards/)
60. [Decoupling Dynamics in Full-Duplex Models](https://www.emergentmind.com/papers/2607.07148)
61. [What Did I Just Say? Self-Listening for Full-Duplex Speech Models](https://arxiv.org/pdf/2609.05592)
62. [Piggybacking on Perception: Stealthy Concurrent Audio Prompt Injections against Multimodal LLM Agents](https://arxiv.org/html/2607.28165)
63. [(PDF) Multimedia data collection of in-car speech communication](https://www.academia.edu/23594999/Multimedia_data_collection_of_in_car_speech_communication)
64. [Natural Language Interactions in Autonomous Vehicles: Intent Detection and Slot Filling from Passenger Utterances](https://arxiv.org/pdf/1904.10500)
65. [AISHELL-5: The First Open-Source In-Car Multi-Channel Multi-Speaker Speech Dataset for Automatic Speech Diarization and Recognition](https://arxiv.org/pdf/2505.23036)
66. [\[2204.09109\] From Spoken Thoughts to Automated Driving Commentary: Predicting and Explaining Intelligent Vehicles’ Actions](https://ar5iv.labs.arxiv.org/html/2204.09109)
67. [NARRATE: A Multimodal Real-World Australian Driving Dataset for Human-Centred Explanations in Automated Driving](https://arxiv.org/html/2608.14767)
68. [NVIDIA Alpamayo: Open AI for Level 4 Autonomous Vehicles](https://www.nvidia.com/en-eu/solutions/autonomous-vehicles/alpamayo/)
69. [\[Literature Review\] Decoupling Conversational Dynamics in Full-Duplex Spoken Models through Reinforcement Learning](https://www.themoonlight.io/en/review/decoupling-conversational-dynamics-in-full-duplex-spoken-models-through-reinforcement-learning)
70. [GitHub - LiYuxin44/DuplexPO\_official · GitHub](https://github.com/LiYuxin44/DuplexPO_official)
71. [arxiv.org](https://arxiv.org/pdf/2503.09594)
72. [Generate Trajectories, Reasoning Traces, and Auto-Labels with NVIDIA Alpamayo 2 Super](https://developer.nvidia.com/blog/generate-trajectories-reasoning-traces-and-auto-labels-with-nvidia-alpamayo-2-super/)
73. [EU ADDW mandate: What drivers need to know](https://www.adt.media/vehicle-connectivity/what-the-eus-new-addw-mandate-means-for-drivers/2718653)
