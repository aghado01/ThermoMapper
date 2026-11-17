# How to Train a Model Like Dessa
Dessa is a full‑duplex speech‑to‑speech model trained from scratch on raw human conversation audio, with no text or pretrained language models, and it learns timing, emotion, and turn‑taking directly from sound 

### Core Training Approach
Dessa’s uniqueness lies in its speech‑only foundation:

- Data: Hundreds of hours of real human conversation audio, with no transcripts or text tokens 
- No pretrained models: Trained from scratch, unlike most LLMs that start from text‑based pretraining 
- Full‑duplex: Learns to listen and speak simultaneously, capturing pauses, laughter, interruptions, and emotional prosody 
- Efficiency: Acquires fluent language from “hundreds of times fewer words” than text LLMs 

## Technical Steps (Conceptual)
While Block’s exact architecture is proprietary, the general pipeline for a similar model would involve:

### Data Collection
- Gather large, diverse, high‑quality human conversation recordings (clean and noisy environments, multiple speakers, varied topics and emotions) 
- Ensure balanced representation of turn‑taking, interruptions, and emotional cues.

### Preprocessing
- Convert audio to waveform or feature representations (e.g., mel‑spectrograms, raw PCM).
- Normalize, trim silence, and segment into utterance‑level chunks.
- Label turn boundaries and optional metadata (speaker ID, emotion tags) if needed for supervision.

### Model Architecture
- Use a speech‑to‑speech architecture: encoder for incoming audio, decoder for outgoing audio.
- Incorporate full‑duplex design: parallel listening and speaking streams, possibly with attention over both inputs and outputs 
- Include components for prosody, timing, and emotional modulation.

### Training

Train end‑to‑end on raw audio pairs (input: previous utterance + user input; output: model’s response).
- Use large‑scale GPU clusters (Block used a single NVIDIA GB200 NVL72) 
- Optimize for both linguistic fluency and natural conversational behavior (timing, backchanneling, emotional expression).

### Evaluation

Test on unsteered live conversations, measuring fluency, latency, emotional responsiveness, and turn‑taking accuracy 

### Practical Considerations
- Scale: You need massive, high‑quality audio datasets; text‑only data won’t capture prosody and timing.
- Compute: Training such a model requires significant GPU resources; Block’s approach shows it can be done on a single high‑end rack, but still far beyond consumer‑grade hardware 


### How to Start Your Own Experiment
If you want to explore similar ideas:

- Begin with speech‑only datasets (e.g., LibriSpeech, Common Voice) and build a speech‑to‑speech model.
- Use audio feature encoders (e.g., wav2vec, HuBERT) and autoregressive decoders for speech.
- Focus on full‑duplex behavior by training on paired input/output audio streams.
- Start small: train on short, curated conversational clips before scaling to large datasets.

In short, training a model like Dessa requires massive, high‑quality speech data, a speech‑to‑speech architecture, and a focus on natural conversational dynamics—not just language understanding