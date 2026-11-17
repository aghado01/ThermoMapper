# Full-Duplex Speech-to-Speech Models at the "Dessa Level": Requirements, Lineage, and the Most Promising Data-Efficiency Ideas (as of Oct 28, 2026)

Getting to a Dessa-class model depends mainly on three things: the right representation of speech (low-rate, semantically rich tokens), a very large amount of two-channel conversational audio, and rack-scale compute. Physical priors or wavelet front ends matter much less. Block's claim of learning language from speech alone with "hundreds of times fewer words" goes against the best published evidence, and right now it can't be checked. Your single-H200 plan can realistically reproduce dGSLM-to-Moshi-style conversational dynamics. It cannot reproduce semantic competence learned from audio alone.

## TL;DR

- **Dessa is still a vendor claim, not a verified result.** Block (Oct 23, 2026) says it trained a full-duplex speech-to-speech model from scratch on conversational audio only, on one GB200 NVL72 rack, and that it runs faster than real time on an M4 Max MacBook Pro. Block has published no parameter count, tokenizer, data size, benchmarks or paper.\[1\] The published literature (Cuervo & Marxer 2024; Stanford's SODA 2026; Apple's 16B continuous-diffusion SLM 2026) consistently finds that speech-only models learn syntax and semantics far more slowly than text LMs.\[2\]\[3\]\[4\] If Dessa holds up, it most likely depends on a much better tokenizer plus very large conversational data (by my estimate, on the order of 10⁵–10⁶+ hours) and roughly 10²³ FLOPs of compute.
- **The "learn more from less" techniques with real evidence are mostly about representation and data, not architecture.** In rough order of evidence: (1) semantic, low-frame-rate tokens (syllable-level ~4–8 Hz; HuBERT-25 Hz beats 50 Hz; keep acoustic detail out of the language-model stream); (2) synthetic, story-like spoken data; (3) warm-starting from text LMs, which is the single largest lever but rules out "speech-only"; (4) preference optimization or RL after pretraining; (5) compute-optimal allocation. Learnable or wavelet/constant-Q front ends and DDSP/source-filter priors have no published evidence of improving language learning in SLMs. Their plausible payoff is in codec and vocoder quality, pitch/prosody handling, and efficient streaming.
- **For your project:** freeze Mimi, model both channels dGSLM/Moshi-style on Open Yap 1K (plus Fisher if licensed), evaluate with dGSLM turn-taking statistics and Full-Duplex-Bench, and add DuplexPO-style RL for timing. On one H200, expect natural timing, backchannels, laughter and prosody but semantically shallow content, which is also what dGSLM found. The most original contribution you could make is a rigorous study of pitch-equivariant (log-frequency) or streaming CWT front ends as a codec or prosody component, compared against Mimi at matched bitrate and latency.

## Key Findings

### 1. What is actually known about Dessa

**What Block disclosed (all vendor claims):**
- Dessa was announced Oct 23, 2026 as a research preview. It is audio in, audio out, and full-duplex: it "listens and talks at the same time," and it learned "when to cut in and when to backchannel" from real conversation.\[5\]
- Training was "from scratch on human conversation… No text, no pretrained models." Block says "we figured out how to train Dessa on only one NVIDIA GB200 NVL72" (72 Blackwell GPUs plus 36 Grace CPUs in one NVLink domain).\[1\]\[6\]
- Block says it "learned language from hundreds of times fewer words" than LLMs. Its justification is that a spoken word takes thousands of times more bytes than the written word. No corpus size or comparison baseline is given.\[1\]
- It runs faster than real time locally on a MacBook Pro M4 Max or better. Jack Dorsey's announcement on X called it "the most authentic full-duplex speech model ever built." Thurairatnam's X post adds "consumer GPUs" and calls it "a GPT3 moment for speech."
- The demos are live recordings with two synchronized speaker tracks, subtitles added afterwards. They show laughter mid-sentence, getting louder when asked to speak up, finishing a user's sentence when the user blanks on a movie title, gasps, and "teary gratitude." Block says the model chooses its own voice and personality. One clip pairs Dessa with OpenAI's GPT-Live-1.\[1\]
- Access is limited to early access for a small group of AI safety researchers plus a waitlist. No weights, API, pricing, parameter count, latency figures, benchmarks or technical report have been published. As of Oct 27, 2026, Block's page was still the only official write-up.\[1\]\[7\]

**Team lineage (verifiable).** The lab was founded by former members of Dessa, the Toronto startup Square acquired in 2020. It says it has built generative speech models for more than 8 years and reports directly to Jack Dorsey.\[8\] Hashiam Kadhim and Rayhane Mama built RealTalk (May 2019, with Joseph Palermo), a text-to-speech clone of Joe Rogan's voice trained on 10 hours of audio from his show (Dessa's RealTalk page: "For our final dataset, we used 10 hours of audio from Joe Rogan's show"). It generated breaths and "um"s, and the team withheld code and data out of misuse concerns. Rayhane Mama is also widely known for a popular open-source Tacotron-2 implementation. The team is therefore strong in high-fidelity generative speech and vocoding. That background makes a custom codec or tokenizer a plausible "secret sauce," but this is my inference, not something Block disclosed.

**What I could not find:** a paper, talk, model card, tokenizer description or data-provenance statement from the team as of Oct 28, 2026. Secondary write-ups (TechCrunch, Unite.AI, DataStudios, eesel, voiceaispace) all trace back to Block's page and the X posts.

### 2. Is "speech-only, hundreds of times fewer words" plausible?

**Back-of-envelope (my estimates, not Block's):**
- Conversational speech runs about 196 words per minute for both speakers combined over elapsed time, per Yuan, Liberman & Cieri (Interspeech 2006) on 2,438 Switchboard calls ("we find an overall average rate of 196 words per minute"), so roughly 12k words per conversation-hour.
- Frontier text LMs train on more than 10¹³ tokens (Meta says Llama 3 "is pretrained on over 15T tokens"). "Hundreds of times fewer" therefore implies roughly 10¹⁰–10¹¹ words, which is about 1–10 million hours of conversation.
- That is far beyond any public two-channel corpus. The largest open one, DuplexChat, has 282,634 h of English recovered from podcasts by source separation.\[9\]
- Compute: one GB200 NVL72 running at about 40% utilization for about 3 months gives very roughly 10²³–5×10²³ BF16 FLOPs. That is about an order of magnitude more than the largest published cold-start audio-LM sweep (SODA, 1.3×10²² FLOPs).
- Faster-than-real-time inference on an M4 Max is consistent with a roughly 1–8B-parameter model at a low frame rate, quantized. Kyutai's 7B Moshi already ships int4/int8 MLX builds for Apple silicon.\[10\]

**What the literature says:**
- **Cuervo & Marxer (EMNLP 2024, arXiv 2404.00685).** They trained more than 50 HuBERT-token SLMs (20M–823M parameters, up to 82B tokens, about 160k hours including synthetic data). Linguistic performance scales with compute exponents of 0.021 (sBLIMP), 0.025 (topic cloze) and 0.017 (story cloze), compared with 0.066 / 0.039 / 0.046 for text LMs (Pythia). For syntax, "SLMs require 10^3.14 ΔC to get the same ΔQ" as text.\[3\] This is the source of the "up to three orders of magnitude" claim. Note that it is a ratio of compute-scaling exponents, not "3× more data," as some later papers misquote it.\[11\]\[12\] Their best 823M model reached sBLIMP 61.3 and StoryCloze 56.7, against Pythia-6.9B at 80.0/76.2 and humans at about 79.9 on story cloze.\[3\]
- **SODA (Stanford, Feb 2026, arXiv 2602.16687).** Models were trained from scratch on about 95% speech (Yodas plus Emilia, about 300k h) as Mimi tokens at 100 tokens/s, from 135M to 4B parameters on 500B tokens. sBLIMP stayed near chance in the IsoFLOP sweep. The final 4B base model scored sBLIMP 52.4. The authors "do not observe strong emergent skills in speech or audio capabilities from pre-training alone." Optimal data grows faster than model size (D* ∝ C^0.579, N* ∝ C^0.367), at roughly 13 tokens per parameter at 10²⁰ FLOPs.\[2\]
- **Apple (2026, arXiv 2604.24416).** A 16B continuous-diffusion speech-only model trained on "tens of millions of hours of conversational speech" still found "long-form linguistic coherence remains elusive." The authors suggest further scaling "is impractical unless a new speech representation or modeling paradigm emerge." They compare current SLMs to the proficiency of "a three- to four-year-old child."\[4\]
- **BabySLM (Interspeech 2023, arXiv 2306.01506):** the speech LM STELA trained on 1,024 h of child-centered audio scored at chance (49.5% lexical, 50.3% syntactic), while the text model BabyBERTa trained on 5M words scored 70.4% syntactic.

**My assessment.** Dessa's demos mostly show conversational pragmatics: timing, backchannels, affect and short coherent replies. Learning those from audio alone is well established, and dGSLM showed it in 2022 with 2,000 h of Fisher.\[13\] The extraordinary part of the claim is "learned language." That needs evidence on knowledge, multi-turn coherence and factuality, and Block has not provided any.\[1\] If the claim is true, the likely ingredients are:
1. a tokenizer far better aligned to linguistic units (low frame rate, semantic, possibly variable-rate or syllabic);
2. far more conversational audio than any public corpus, possibly separated from podcasts or other media, as DuplexChat does;
3. about 10²³ FLOPs;
4. a training recipe that keeps acoustic detail from crowding out semantics (separate heads or streams, as Lychee-FD and SODA's token-type ablation suggest).

Treat "no text" carefully. Text-derived tools (ASR filtering, TTS-generated data) could still be in the data pipeline, and Block has not addressed this.

### 3. Landscape and lineage of full-duplex models

**Model-level vs system-level duplex.** In system-level (cascaded) duplex, a VAD or turn detector gates an ASR→LLM→TTS pipeline, or a small state predictor (e.g., Freeze-Omni's frozen LLM with chunk-wise state prediction, FireRedChat, SoulX-DuPlug, DuplexCascade) decides when to speak. In model-level duplex, one network consumes the user stream and emits its own stream every frame, so timing decisions are learned. OpenAI explicitly describes GPT-Live as removing "the turn detector from the audio path."\[14\]

| Family | Representative models | Mechanism | Notes |
|---|---|---|---|
| Textless dual-channel | dGSLM (TACL 2023) | HuBERT + k-means units at 50 Hz per channel; dual-tower Transformer with cross-attention; HiFi-GAN | 2,000 h Fisher, no text; natural turn-taking, laughter; weak semantics |
| Parallel multi-stream + text scaffold | Moshi (Kyutai 2024), PersonaPlex (NVIDIA 2026) | Mimi tokens for user and model streams + "inner monologue" text stream; RQ-Transformer (temporal 7B + small depth transformer over codebooks); acoustic delay | Helium 7B text LM backbone; 160 ms theoretical / ~200 ms practical latency |
| Next-token-pair / dual-channel AR | NTPP (2025) | Predicts both speakers' tokens jointly each step | Evaluated with dGSLM turn-taking stats |
| Time-chunked interleaving | SyncLLM (Meta/UW 2024), OmniFlatten (2024) | Interleave/flatten both speakers into a single sequence in fixed time chunks | Reuses a text LLM directly |
| Listen-while-speaking fusion | LSLM (2024) | Fuses listening channel into embeddings; predicts only assistant tokens | Fewer IPUs/gaps in dGSLM-style eval |
| Codec-free / continuous | SALMONN-omni (2024–25) | No codec injection into the LLM | Standalone duplex speech LLM |
| Decoupled control vs content (2026) | DuplexPO, Diffuplex, DiffuPlex (rolling masked diffusion), Lychee-FD, Listen-Write-Speak, DuplexSLA, HiPLEX | Separate "when to speak" from "what to say"; diffusion for content; separate acoustic/semantic heads; visible text channel; synchronized actions | Reported trend: decoupling resolves the "intelligence–dynamics trade-off" |
| Commercial | OpenAI GPT-Live-1 / mini (Jul 8, 2026), Sesame CSM, Google | Full-duplex voice model; GPT-Live delegates hard reasoning to GPT-5.5 in background | Closed; few technical details |

**Key 2026 papers:**
- **DuplexPO** (arXiv 2607.07148, NTU/NTU-Taiwan/HKUST/NVIDIA) argues that the intelligence–dynamics trade-off "is not fundamental."\[15\] It applies GRPO-style RL with timing-only, factorized rewards on "dynamics-critical windows" (turn transitions, backchannels, barge-ins) drawn from Fisher and Seamless. It reports that supervised fine-tuning alone on natural dynamics *regressed* the dynamics metrics.\[16\] Code is public.\[17\]
- **Diffuplex** (arXiv 2610.17780, Oct 13, 2026) runs an autoregressive interaction policy that tracks the user while a discrete-diffusion backbone (DiffusionGemma) generates and refines the response in parallel.\[18\] It reports first-audio latency on par with other full-duplex models while "exceeding existing open-source systems by more than 20 points in spoken intelligence."\[19\]
- **DiffuPlex** (Univ. of Seoul, arXiv 2610.12214) uses rolling masked diffusion to draft future frames, revising only the unplayed part. It reports 1.46×/1.59× wall-clock speedups, with every backbone call finishing inside the 80 ms Mimi frame.\[20\]
- **Lychee-FD** (ACL 2026, arXiv 2607.06540) attributes "modality interference" to "inherent gradient conflicts between acoustic and semantic modeling when the two modalities are forced to share a deep parameter space."\[21\] It separates the deep layers into acoustic and semantic heads. Weights and code are open.\[22\] This is directly relevant to a speech-only design.
- **PersonaPlex** (arXiv 2602.06053, Jan 2026) is Moshi plus role and voice prompts. It was fine-tuned on 7,303 Fisher conversations (1,217 h) plus synthetic dialogues: 39,322 assistant (410 h) and 105,410 customer-service (1,840 h), with LLM-written scripts and TTS voices. Reported training took **6 h on 8×A100** (per The Decoder), which shows how cheap behavioral fine-tuning is once a pretrained duplex backbone exists.\[23\]\[24\]\[25\]
- **NemotronLabs VoiceChat** (arXiv 2609.21967) is NVIDIA's open full-duplex model with tool calling.
- **F-Actor** (arXiv 2601.11329) covers controllable backchannel and interruption behavior, trained on the 2,164-h synthetic Behavior-SD corpus. The authors deliberately restrict it to a closed speaker set to prevent cloning.\[26\]

### 4. Codecs and tokenizers: the actual foundation

- **SoundStream / EnCodec / DAC** are RVQ-GAN acoustic codecs. Their many codebooks at 50–75 Hz carry everything (phase, timbre, content), which makes them hard for an LM to predict.\[27\]
- **Mimi (Kyutai)** works at 24 kHz with 12.5 Hz frames, 8 codebooks × 2048 entries, and 1.1 kbps. It is fully streaming, with 80 ms latency (one frame). It has Transformers in both encoder and decoder and uses adversarial-only training plus feature matching.\[28\]\[29\] Its first codebook is distilled from WavLM, which makes it "semantic," a SpeechTokenizer-style split. Kyutai chose 12.5 Hz to get closer to text's "~3–4 Hz" rate and reduce the number of AR steps.\[29\]
- **Syllable-level tokens:**
  - SyllableLM (arXiv 2410.04029) produces units at 5–8.33 Hz and as low as 60 bps. It cut resynthesis WER from 37% (SD-HuBERT units) to 7%, and its 90M-parameter LM matches or beats SOTA SpeechLMs.\[30\]
  - Sylber (arXiv 2410.07168) averages 4.27 tokens/s, 6–7× fewer than HuBERT.\[31\]
  - A 2025 scaling study (arXiv 2509.26634) found Sylber-token SLMs match or beat HuBERT-token SLMs with 5× fewer tokens, >2× less training time and 5× fewer FLOPs.\[32\]\[33\]
  - Sylber 2.0 (arXiv 2601.22306, Jan 2026) is multilingual at about 5 Hz and can resynthesize at 24 kHz.\[34\]
- **Counter-evidence:** Cuervo & Marxer found that unigram (BPE-like) compression of HuBERT tokens gave *worse* downstream scaling, and StoryCloze stopped scaling with compute.\[3\] So coarser is not automatically better. *Linguistically aligned* coarsening, as with syllables, seems to be what matters.
- **Low-bitrate single-codebook codecs:** WavTokenizer, BigCodec, TAAE (FSQ at 50 Hz, 0.8 kbps), SNAC (multi-scale), and 2026 work like ReLMCodec ("predictable speech tokens").\[35\] These improve reconstruction per bit, but predictability by an LM, not fidelity, is the bottleneck for semantics.
- **Wavelets / constant-Q:**
  - Direct evidence is thin. "Wavelet as Tokenizer" (arXiv 2606.02631) is explicitly preliminary and "do[es] not establish a production-quality tokenizer."\[36\]
  - Related ideas: SpecTokenizer (a compressed-spectrum-domain streaming codec), Vocos (a Fourier-domain vocoder), and SNAC (multi-scale temporal resolution, conceptually close to a dyadic decomposition).
  - I found no paper showing that a CQT or CWT front end improves SLM linguistic learning.

### 5. "Learning more from less": what moves the needle

| Technique | Evidence | Strength | Fits "speech-only"? |
|---|---|---|---|
| Text warm-start (TWIST) | TWIST trained on 10% of data matched or beat cold-init trained on 100%; Slam (Qwen-0.5B init) hit sBLIMP 58.86 / tSC 82.04 on **one A5000 in 24 h** | Very strong | No |
| Semantic-only LM tokens | SODA at 3×10²⁰ FLOPs: semantic tokens sBLIMP 58.6 / sWUGGY 72.1 vs 50.9 / 59.0 when acoustic tokens added | Strong | Yes |
| Lower, linguistically aligned frame rate | TWIST: 25 Hz > 50 Hz on sWUGGY/sBLIMP; syllabic tokens 5× fewer training tokens | Strong | Yes |
| Synthetic story-like speech | sTinyStories beats audiobooks "across all model scales" (Cuervo & Marxer); real+synthetic best overall (Maimon et al. 2025) | Strong | Partly (TTS needs text) |
| Preference optimization / RL | Slam DPO: sSC 55.59→58.04; Align-SLM (RLAIF); DuplexPO for timing | Moderate–strong | Yes (with audio-native rewards) |
| Compute allocation | Interleaved SLMs: favor model size (7B on 4.2B tokens beat Cuervo's 463M/72B suggestion at 2e20 FLOPs); textless SODA: data grows faster than model size | Moderate | Yes |
| Decoupling acoustic and semantic parameters | Lychee-FD gradient-conflict analysis; DuplexPO/Diffuplex control vs content split | Moderate | Yes |
| In-domain SSL encoder | dGSLM: HuBERT trained on Fisher gave better phonetic ABX than LibriSpeech HuBERT | Moderate | Yes |
| SSL objectives (HuBERT, wav2vec 2.0, data2vec, BEST-RQ) | Basis of all semantic tokens; data2vec/JEPA-style latent prediction is the closest relative of "predict abstractions, not samples" | Foundational | Yes |
| Child-acquisition framing (BabySLM, BabyLM) | Benchmarks exist; speech still lags text badly at child-scale input | Diagnostic, not yet a fix | Yes |
| Learnable / physical front ends (SincNet, LEAF, CQT/CWT, DDSP, source-filter) | Gains in classification and synthesis efficiency; no SLM-semantics evidence found | Weak / untested | Yes |

How to read the table: every proven "less data" win either (a) imports text knowledge or (b) reduces the information rate the LM must model per linguistic unit. Physical priors fit under (b) only if they let you discard or factor out nuisance variation (pitch, speaker, channel) *before* the LM. That is where your ideas belong.

### 6. Turn-taking, evaluation, and latency

- **Corpus statistics (dGSLM):** inter-pausal units (IPUs), pauses, gaps and overlaps, measured as counts and cumulative durations per minute, compared against held-out Fisher.\[37\] NTPP and others reuse this protocol, and it is the cheapest first eval for your project.\[38\] Open Yap 1K metadata already includes per-conversation turns per minute (2.7–19.8) and turn-taking gaps (320–860 ms).\[39\]
- **Full-Duplex-Bench (Lin et al., NTU/CMU):**
  - v1 (ASRU 2025, arXiv 2503.04721): pause handling, backchanneling, smooth turn-taking, user interruption. Metrics include takeover rate, response latency and backchannel frequency.
  - v1.5 (ICASSP 2026, arXiv 2507.23159): overlap handling.
  - v2 (arXiv 2510.07838): multi-turn, with an automated spoken-LLM examiner across Daily, Correction, Entity Tracking and Safety tasks. Its finding: systems "get confused when people talk at the same time, struggle to handle corrections… lose track of who or what is being talked about."\[40\] A real-time WebRTC orchestrator was released Feb 23, 2026.\[41\]
  - v3 (arXiv 2604.04847): tool use under disfluency.
- **Others:** Talking Turns, Game-Time (temporal dynamics), MTR-DuplexBench (multi-round), FD-Bench, TurnBench, EchoChain, the ICASSP 2026 HumDial challenge, and paralinguistics benchmarks (ParaS2S, SpeechParaling-Bench).
- **Latency definitions:** Moshi's floor is 80 ms (one Mimi frame) plus 80 ms acoustic delay, which gives the 160 ms theoretical figure and about 200 ms measured on an L4. Full-duplex benchmarks report response onset after user end-of-turn, and yield time after barge-in. Report both, plus real-time factor.

### 7. Engineering requirements and budgets

| System | Params | Audio data | Compute |
|---|---|---|---|
| dGSLM (2022) | Dual-tower Transformer (small by today's standards) | 2,000 h Fisher, two-channel | Academic scale |
| Moshi (2024) | 7B temporal + depth transformer; Mimi | 7M h pretraining (Whisper-transcribed), 2,000 h Fisher (8 kHz→24 kHz via AudioSR), 170 h Kyutai conversations, 20k h synthetic | 127 DGX nodes = 1,016 H100s |
| PersonaPlex (2026) | 7B (from Moshi) | 1,217 h Fisher + ~2,250 h synthetic | 6 h on 8×A100 (fine-tune) |
| SODA (2026) | 135M–4B | ~300k h, 500B tokens | up to 1.3×10²² FLOPs |
| Apple continuous SLM (2026) | 16B | tens of millions of hours | not given here |
| Dessa (claimed) | undisclosed | undisclosed conversational audio | 1× GB200 NVL72 |

**Two-channel data:**
- **Open Yap 1K:** 1,000 h, 1,602 conversations, 239 speakers, 48 kHz/16-bit, one track per speaker, Deepgram Nova-3 word timings. Speakers talk with people they already know, so expect more interruptions and backchannels.\[42\]\[43\] The full corpus is free for commercial and research use under a DUA; the 8.9-h HF sample is CC BY 4.0.\[44\]
- **Fisher:** about 2,000 h, 8 kHz telephone, LDC license.
- **otoSpeech 280h:** per your notes; I did not verify it.
- **DuplexChat:** 282,634 h English plus 132,723 h Japanese, separated from podcasts with DialogueSidon.\[9\] Check the license and provenance carefully, since separation artifacts exist and podcast rights are unclear.
- **Behavior-SD:** 2,164 h synthetic.
- **Seamless Interaction / Naturalistic (Meta):** used in DuplexPO.

**Echo and self-voice.** Model-level duplex models receive their own output as a separate stream, so they "know" what they said. In open-mic deployment, though, the user channel also contains the model's playback, so you need acoustic echo cancellation (e.g., WebRTC AEC3) or headphones, or you need to train with echo-augmented user channels. Lychee-FD calls "unconstrained, real-world open-mic scenarios" an open frontier.\[22\]

**On-device.** Moshi provides int4/int8/bf16 builds in PyTorch, MLX and Rust.\[10\] Faster than real time at 12.5 Hz means each frame step (temporal plus depth transformer plus Mimi decode) must finish in under 80 ms. DiffuPlex shows that compute budget is the binding constraint.

**What one H200 can do (my estimates):**
- Token counts: Open Yap 1K with Mimi gives 3.6M s × 12.5 Hz × 8 codebooks × 2 channels ≈ 720M tokens, or about 90M if you model only the semantic codebook per channel.
- A 100–300M-parameter dual-stream model at about 6ND FLOPs trains for several epochs in hours to a day.
- A week on one H200 at about 40% utilization is roughly 2×10²⁰ FLOPs. That is SODA's IsoFLOP regime, where cold-start speech semantics stayed near chance.
- So expect dGSLM-class behavior. Semantics would need either a text warm-start (Slam shows one GPU-day suffices) or fine-tuning an open duplex backbone (Moshi or PersonaPlex weights).

### 8. Safety

- Block frames Dessa's purpose partly as safety. It says current AI is "nearly deaf" to emotion and gives the example of an assistant that carries on normally while the user is in pain.\[1\] Distress detection from prosody is the stated benefit.
- The risks are the mirror image: a convincingly human, emotionally reactive model that runs locally invites impersonation and voice phishing. Block describes no watermarking, cloning restrictions or release criteria, only restricted access for safety researchers.\[1\]
- A speech-only model also has no intermediate transcript to log, filter or audit before it speaks.\[1\]
- The team's 2019 RealTalk project was itself a deliberate deepfake-awareness exercise.\[45\]\[46\]
- For your work, follow F-Actor's practice: train or serve with a closed speaker set and don't release voice-prompting cloning ability.\[26\]

## Recommendations

**Project roadmap (single H200):**
1. **Baseline: dGSLM-style replication on Mimi.** Freeze Mimi and precompute tokens for Open Yap 1K (plus Fisher, upsampled, if you have LDC access). Model only codebook-1 (semantic) for both channels with a dual-stream Transformer, either summed embeddings as in Moshi or cross-attention as in dGSLM. Evaluate IPU/pause/gap/overlap distributions against held-out data.
2. **Add acoustics.** Add the remaining 7 codebooks through a small depth transformer with an acoustic delay of 1–2 frames, following Moshi's RQ-Transformer. Measure the real-time factor on the H200 and later with MLX on a Mac.
3. **Separate content and acoustics in the network,** following Lychee-FD's separate deep heads, and run an ablation. This is cheap and directly tests the strongest "speech-only" hypothesis.
4. **Timing RL.** Port DuplexPO's public reward and windowing code to tune backchannels, yields and turn onsets. Evaluate on Full-Duplex-Bench v1/v1.5.
5. **Optional semantic track.** Run two arms and report the gap honestly: (a) cold-start speech-only pretraining on monologue audio, e.g. LibriLight-scale, then dialogue fine-tuning; (b) a Slam-style Qwen-0.5B warm-start. This gap is your own small-scale measurement of "the Dessa question."

**Most promising open directions for an independent researcher:**
- **Pitch-equivariant tokenization (best fit for your background).** On a log-frequency (CQT/CWT) axis, a pitch shift becomes a translation. Convolution along that axis can then factor F0 and speaker out of the content stream before the LM, which is exactly the "reduce information rate per linguistic unit" lever the evidence supports. Test it as a codec encoder or a semantic-distillation target and compare against Mimi at matched bitrate and latency. Report the LM's sBLIMP/sWUGGY/BabySLM scores and turn-taking statistics, not just codec MUSHRA scores. Evidence today: none for SLMs, so this is an open question.
- **Streaming "scrolling spectrogram" with a fixed analysis lag.** CWT's frequency-dependent time support conflicts with a fixed causal latency budget, since low frequencies need long windows. Frame the design as an explicit lag-versus-resolution trade-off against Mimi's 80 ms frame. Complex coefficients and phase derivatives (instantaneous frequency) could make a better vocoder-side target, in the spirit of Vocos, than for the LM.
- **Source-filter / DDSP decoders** as the waveform generator conditioned on low-rate tokens. These give interpretable F0 and loudness controls (e.g., "speak up") and need fewer decoder parameters, which matters for on-device inference. They are unlikely to affect linguistic learning.
- **Variable-rate / syllabic duplex tokens.** No one has yet published a full-duplex model on syllable-rate tokens. The open problem is reconciling variable-rate units with fixed-clock duplex timing.
- **"Deep Teaching" analogue** (unsupervised learning plus physical priors). The closest validated speech versions are SSL with in-domain data plus physically structured decoders. Pitch the idea as "priors reduce nuisance entropy," not as "priors teach language."

## Prioritized Reading List

1. dGSLM: Generative Spoken Dialogue Language Modeling — https://arxiv.org/abs/2203.16502
2. Moshi technical report — https://arxiv.org/abs/2410.00037 ; code and Mimi — https://github.com/kyutai-labs/moshi ; model card — https://huggingface.co/kyutai/moshiko-pytorch-bf16
3. Cuervo & Marxer, Scaling Properties of Speech Language Models — https://arxiv.org/abs/2404.00685
4. SODA, Scaling Open Discrete Audio Foundation Models — https://arxiv.org/abs/2602.16687
5. Apple, Scaling Properties of Continuous Diffusion SLMs — https://arxiv.org/abs/2604.24416
6. Slamming (one GPU, one day) — https://arxiv.org/abs/2502.15814 ; code — https://github.com/slp-rl/slamkit
7. Scaling Analysis of Interleaved Speech-Text LMs — https://arxiv.org/abs/2504.02398
8. TWIST — https://arxiv.org/abs/2305.13009
9. SyllableLM — https://arxiv.org/abs/2410.04029 ; Sylber — https://arxiv.org/abs/2410.07168 ; Sylber 2.0 — https://arxiv.org/abs/2601.22306 ; Syllabic scaling — https://arxiv.org/abs/2509.26634
10. PersonaPlex — https://arxiv.org/abs/2602.06053
11. DuplexPO — https://arxiv.org/abs/2607.07148 ; code — https://github.com/LiYuxin44/DuplexPO_official
12. Lychee-FD — https://arxiv.org/abs/2607.06540
13. Diffuplex (discrete diffusion) — https://arxiv.org/abs/2610.17780 ; DiffuPlex (rolling masked diffusion) — https://arxiv.org/abs/2610.12214
14. NTPP — https://arxiv.org/abs/2506.00975 ; SyncLLM — https://arxiv.org/abs/2409.15594 ; Freeze-Omni — https://arxiv.org/abs/2411.00774 ; SALMONN-omni — https://arxiv.org/abs/2411.18138
15. Full-Duplex-Bench v1/v1.5/v2/v3 — https://arxiv.org/abs/2503.04721 , https://arxiv.org/abs/2507.23159 , https://arxiv.org/abs/2510.07838 , https://arxiv.org/abs/2604.04847 ; repo — https://github.com/DanielLin94144/Full-Duplex-Bench
16. Survey of Full-Duplex Spoken Dialogue Systems (2026) — https://arxiv.org/abs/2606.19453 ; curated list — https://github.com/Ruiqi-Yan/Awesome-Full-Duplex-SDM
17. DuplexChat corpus — https://arxiv.org/abs/2607.04941 ; Open Yap 1K — https://huggingface.co/datasets/TheAgenticDataCompany/open-yap-1k
18. BabySLM — https://github.com/MarvinLvn/BabySLM
19. SSL and codec foundations: HuBERT https://arxiv.org/abs/2106.07447 ; wav2vec 2.0 https://arxiv.org/abs/2006.11477 ; data2vec https://arxiv.org/abs/2202.03555 ; BEST-RQ https://arxiv.org/abs/2202.01855 ; SoundStream https://arxiv.org/abs/2107.03312 ; EnCodec https://arxiv.org/abs/2210.13438 ; DAC https://arxiv.org/abs/2306.06546 ; SpeechTokenizer https://arxiv.org/abs/2308.16692 ; BigCodec https://arxiv.org/abs/2409.05377 ; SNAC https://arxiv.org/abs/2410.14411
20. Priors and front ends: DDSP https://arxiv.org/abs/2001.04643 ; LEAF https://arxiv.org/abs/2101.08596 ; SincNet https://arxiv.org/abs/1808.00158 ; Vocos https://arxiv.org/abs/2306.00814 ; Wavelet as Tokenizer https://arxiv.org/abs/2606.02631

## Caveats

- Everything about Dessa's internals is unknown. My data and compute figures are order-of-magnitude inferences from Block's phrasing, not disclosures.
- Several 2026 papers cited here are recent arXiv preprints without independent replication. Diffuplex's "+20 points" and DuplexPO's "breaks the trade-off" are the authors' own claims.
- PersonaPlex's 8×A100 / 6 h figure comes from secondary press (The Decoder), and the synthetic-hour totals differ slightly between secondary sources.
- I did not verify OmniFlatten, LSLM, MiniCPM-o, Sesame CSM or Google's duplex systems in primary sources for this report, and they are only summarized here. I did not verify otoSpeech-full-duplex-280h at all.
- The wavelet/CQT and physical-prior directions are hypotheses with no direct supporting evidence for spoken language modeling. Present them as such.

## Sources

1. [Block Dessa: full-duplex speech-to-speech model trained only on human conversation](https://www.datastudios.org/post/block-dessa-full-duplex-speech-model-trained-only-on-human-conversation)
2. <https://arxiv.org/pdf/2602.16687>
3. <https://arxiv.org/pdf/2404.00685v2>
4. <https://arxiv.org/pdf/2604.24416>
5. [Block Unveils Dessa, a Full-Duplex Speech Model Trained on Conversation](https://www.unite.ai/block-unveils-dessa-a-full-duplex-speech-model-trained-on-conversation/)
6. [Block Dessa: The Speech Model Trained Without Text](https://cellcog.ai/blog/block-dessa/)
7. [Block Dessa explained: Block's speech-only AI voice model](https://www.eesel.ai/blog/block-dessa)
8. [Block launches Dessa, a full-duplex speech model trained only on audio](https://www.voiceaispace.com/articles/block-launches-dessa-a-full-duplex-speech-model-trained-only-on-audio)
9. [DuplexChat: Constructing Speaker-Separated Full-Duplex Dialogue Speech at Scale forSpoken Dialogue Language Modeling](https://arxiv.org/html/2607.04941v1)
10. [Philipp Schmid on X: "Kyutais Voice Model Moshi is here! Moshi is a real-time speech-to-speech-text foundation model achieving 160-200ms end-to-end latency and can be run on-device. 🤯 @kyutai\_labs releases the first set of open Weights, Code, and Technical Reports under a permissive cc-by-4.0 license. TL;DR: 🧠 7B parameter Temporal Transformer 🤗 Available under cc-by-4.0 on @huggingface ⏱️ latency of 160ms - 200ms on an L4 GPU 🎙️ Moshi uses Mimi, a state-of-the-art streaming neural audio codec 🚀 Achieves 12.5 Hz representation with just 1.1 kbps bandwidth 🗣️ Supports two audio streams: one for Moshi, one for the user 📝 Predicts text tokens of its own speech for improved accuracy 🌟 Outperforms existing codecs like SpeechTokenizer and SemantiCodec 🖥️ Available in Python (PyTorch and MLX) and Rust 🤖 Includes fine-tuned models for male (Moshiko) and female (Moshika) voices 📊 Various quantization options (int4, int8, bf16) for different platforms Model: https://t.co/5uRoEAkwPS Github: https://t.co/8mLPCkMz83 Paper: https://t.co/VscBNqt7bF Demo: https://t.co/CHLTgPgDjR" / X](https://x.com/_philschmid/status/1836431576768413701)
11. [Scaling Analysis of Interleaved Speech-Text Language Models](https://arxiv.org/html/2504.02398)
12. [Slamming: Training a Speech Language Model on One GPU in a Day](https://arxiv.org/pdf/2502.15814)
13. [Generative Spoken Dialogue Language Modeling - ACL Anthology](https://aclanthology.org/2023.tacl-1.15/)
14. [How we built a realtime system for responsive voice AI in six months](https://openai.com/index/continuous-voice-interaction-with-gpt-live/)
15. [Decoupling Dynamics in Full-Duplex Models](https://www.emergentmind.com/papers/2607.07148)
16. [Decoupling Conversational Dynamics in Full-Duplex Spoken Models through Reinforcement Learning](https://arxiv.org/html/2607.07148v1)
17. [GitHub - LiYuxin44/DuplexPO\_official · GitHub](https://github.com/LiYuxin44/DuplexPO_official)
18. [\[2610.17780\] Full-duplex Spoken Language Modeling with Discrete Diffusion](https://arxiv.org/abs/2610.17780)
19. [Full-duplex Spoken Language Modeling with Discrete Diffusion](https://arxiv.org/html/2610.17780v1)
20. [DiffuPlex: Accelerating Full-Duplex Spoken Dialog Models via Rolling Masked Diffusion](https://www.alphaxiv.org/abs/2610.12214)
21. [\[2607.06540\] Hierarchical Acoustic-Semantic Modeling: Modality Separation and Semantic Coherence for Full-Duplex SLMs](https://arxiv.org/abs/2607.06540)
22. [Hierarchical Acoustic-Semantic Modeling: Modality Separation and Semantic Coherence for Full-Duplex SLMs](https://arxiv.org/html/2607.06540v3)
23. [PersonaPlex-7B-v1 NVIDIA: Human-Like Speech Model](https://collabnix.com/nvidia-personaplex-7b-v1-the-full-duplex-speech-model-that-talks-like-a-human/)
24. [Nvidia releases open model PersonaPlex, a voice AI that listens and talks at the same time](https://the-decoder.com/nvidia-open-sources-personaplex-a-voice-ai-that-listens-and-talks-at-the-same-time/)
25. [NVIDIA PersonaPlex, What It Is and How to Try the Full-Duplex Voice AI](https://pasqualepillitteri.it/en/news/26826/nvidia-personaplex-explained-full-duplex-voice-ai)
26. [F-Actor: Controllable Conversational Behaviour in Full-Duplex Models](https://arxiv.org/pdf/2601.11329)
27. [Vec-Tok Speech: speech vectorization and tokenization for neural speech generation](https://arxiv.org/pdf/2310.07246)
28. [Moshi: Unified Speech-Text Dialogue Model](https://www.emergentmind.com/topics/moshi-a-speech-text-foundation-model)
29. [GitHub - kyutai-labs/moshi: Moshi is a speech-text foundation model and full-duplex spoken dialogue framework. It uses Mimi, a state-of-the-art streaming neural audio codec. · GitHub](https://github.com/kyutai-labs/moshi)
30. [SyllableLM: Learning Coarse Semantic Units for Speech Language Models](https://arxiv.org/html/2410.04029.pdf)
31. [sylber: syllabic embedding representation](https://arxiv.org/pdf/2410.07168)
32. [Scaling Spoken Language Models with Syllabic Speech Tokenization](https://arxiv.org/html/2509.26634)
33. [Scaling Spoken Language Models with Syllabic Speech Tokenization](https://arxiv.org/pdf/2509.26634)
34. [arxiv.org](https://arxiv.org/pdf/2601.22306)
35. [AffectCodec: Emotion-Preserving Neural Speech Codec for Expressive Speech Modeling](https://arxiv.org/pdf/2605.11098)
36. [Wavelet as Tokenizer: Preliminary Results on a Shared Wavelet Token Schema for Natural Signals](https://arxiv.org/html/2606.02631v1)
37. [Generative Spoken Dialogue Language Modeling](https://direct.mit.edu/tacl/article/doi/10.1162/tacl_a_00545/115240/Generative-Spoken-Dialogue-Language-Modeling)
38. [NTPP: Generative Speech Language Modeling for Dual-Channel Spoken Dialogue via Next-Token-Pair Prediction](https://arxiv.org/pdf/2506.00975)
39. [24892400294f/open-yap-1k · Datasets at Hugging Face](https://huggingface.co/datasets/24892400294f/open-yap-1k)
40. [\[2510.07838\] Full-Duplex-Bench-v2: A Multi-Turn Evaluation Framework for Duplex Dialogue Systems with an Automated Examiner](https://arxiv.org/abs/2510.07838)
41. [GitHub - DanielLin94144/Full-Duplex-Bench: A Benchmark for Evaluating Turn-Taking and Overlap Handling in Full-Duplex Spoken Dialogue Models · GitHub](https://github.com/DanielLin94144/Full-Duplex-Bench)
42. [GitHub - The-Agentic-Data-Co/open-yap-1k: The world's largest open dataset of natural two-speaker English conversation: 1000 hours, commercially usable. · GitHub](https://github.com/The-Agentic-Data-Co/open-yap-1k)
43. [Open Yap 1K - The Agentic Data Company](https://theagenticdatacompany.com/open-yap-1k)
44. [Open Yap 1K: 1,000 hours of full-duplex natural conversation, free for commercial use](https://huggingface.co/blog/TheAgenticDataCompany/open-yap-1k)
45. [RealTalk: We Recreated Joe Rogan's Voice Using Artificial Intelligence](https://ghostarchive.org/varchive/DWK_iYBl8cA)
46. [Dessa - RealTalk: Lifelike Synthesized Speech](https://www.dessa.com/realtalk)
