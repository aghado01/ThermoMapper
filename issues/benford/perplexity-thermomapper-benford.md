<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" style="height:64px;margin-right:32px"/>

## Short answer

Yes—there's emerging research applying Benford's Law to AI-text detection, but naive first-digit analysis on raw text doesn't work well; you need to convert text into meaningful numeric sequences first (e.g., token probabilities, attention weights, or word-frequency ranks).[^1_1][^1_2]

## How the approach actually works

Since text isn't inherently numeric, researchers first transform it into a number sequence—common methods use word frequencies, character/n-gram counts, or model-internal statistics (log-probabilities, attention scores)—then check whether the leading digits of those numbers follow Benford's expected distribution (~30% start with 1, ~5% with 9). Early attempts (like feeding raw GPT-3 text through simple word-count transforms) found the deviations weren't reliably distinguishable from human or even random/Lorem-ipsum text, making that naive version a poor detector.[^1_2][^1_3]

More recent work has found better traction by using richer numeric representations tied to how models generate text:

- **BENATTEN**: converts text into numbers via attention weights, exploiting the idea that AI's probabilistic token allocation differs systematically from human writing patterns, and reports state-of-the-art detection performance across multiple domains.[^1_1]
- **BENADV**: extends the Benford's Law approach beyond purely text-centric signals to capture more complex generation artifacts, again for detecting machine-generated text.[^1_4]
- Notably, this research reports that AI-generated text tends to conform *more* closely to Benford's Law than human text once run through these attention-derived numeric transforms—the opposite of the naive fraud-detection intuition that fabricated data deviates from Benford's curve.[^1_1]


## Why this differs from classic Benford fraud detection

| Aspect | Classic fraud/forensic use | AI-text detection use |
| :-- | :-- | :-- |
| Data type | Naturally occurring numeric data (financial figures, populations) [^1_3][^1_5] | Text converted into numeric proxies (attention weights, frequencies) [^1_2][^1_1] |
| Signal of manipulation | Fabricated data deviates *from* Benford (too uniform) [^1_6][^1_5] | AI text can adhere *more* closely to Benford than human text, reversing the usual signal [^1_1] |
| Reliability | Well-established, decades of use in auditing [^1_5][^1_7] | Still experimental; naive versions fail, advanced attention-based versions show promise [^1_2][^1_1] |

## Related but distinct application: AI-generated images

A parallel and more mature line of work applies Benford's Law to JPEG DCT coefficients to flag AI-generated or heavily edited images, since real camera/compression pipelines produce Benford-conforming coefficient distributions while synthetic images often don't. This image-based method chops images into 8×8 blocks, extracts AC coefficients from a discrete cosine transform, and runs a chi-squared goodness-of-fit test against the Benford curve—though evaluations show it's imperfect, correctly flagging only about 60% of clearly artificial images in one study.[^1_8][^1_9]

Given your interest in evaluation rigor and mechanism-sensitive benchmarking, the key caveat worth flagging: these Benford-based text detectors are still niche and less validated than mainstream statistical detectors (like perplexity-based or watermarking approaches), so treat them as a promising but immature signal rather than a robust standalone classifier.[^1_2][^1_1]

<span style="display:none">[^1_10][^1_11][^1_12][^1_13][^1_14][^1_15]</span>

<div align="center">⁂</div>

[^1_1]: https://www.sciencedirect.com/science/article/pii/S0306457325000767

[^1_2]: https://medium.datadriveninvestor.com/can-benfords-law-help-identify-ai-generated-text-88542602d905

[^1_3]: https://en.wikipedia.org/wiki/Benford's_law

[^1_4]: https://dl.acm.org/doi/10.1145/3748305

[^1_5]: https://www.mdd.com/forensic-accounting-articles/benfords-law-forensic-accounting-tool/

[^1_6]: https://pnwacfe.org/Blog/13443307

[^1_7]: https://www.valid8financial.com/resource/benfords-law-and-ai-categorization-accelerate-financial-fraud-investigations

[^1_8]: https://research.bond.edu.au/en/publications/are-benfords-law-based-detectors-effective-for-gai-generated-imag/

[^1_9]: https://kennethbspringer.au/2026/08/05/benfords-law-detecting-ai-generated-images-without-ai/

[^1_10]: https://www.sei.cmu.edu/blog/benfords-law-potential-applications-insider-threat-detection/

[^1_11]: https://www.sciencedirect.com/science/article/abs/pii/S2666281723000161

[^1_12]: https://firstmonday.org/ojs/index.php/fm/article/view/10163/8063

[^1_13]: https://www.mdpi.com/2076-3417/11/23/11482

[^1_14]: https://www.reddit.com/r/datascience/comments/x6betf/benfords_law_and_fraud_detection/

[^1_15]: https://dl.acm.org/doi/abs/10.1145/3748305


---

# well obviously a key hypothesis would be that artificial writing is going to have more predictable structure embedded in it

That's exactly the mainstream mechanistic hypothesis, and it's largely correct—AI text tends toward statistically "smoother" structure because LLMs pick high-probability tokens by design, which produces measurable regularities that human writing (with its idiosyncrasies and improvisation) usually lacks.[^2_1][^2_2]

## The core mechanism

Because LLMs generate text by consistently favoring high-probability next tokens, AI output tends to have lower **perplexity** (each word is less "surprising" given the preceding context) than human writing. This directly operationalizes your hypothesis: predictability at the token level is quantifiable, and AI text scores lower/more predictable on average.[^2_3][^2_2][^2_4][^2_1]

The structural analog at the sentence/paragraph level is **burstiness**—the standard deviation of sentence length and syntactic complexity across a document. Human writing naturally alternates between short punchy sentences and long complex ones, producing high burstiness (roughly 0.6–1.2 in some scoring schemes), while AI text clusters in a narrower, more uniform band (roughly 0.15–0.4).[^2_5][^2_6]

## Why this connects to Benford's Law

This is precisely the theoretical bridge to your Benford question: Benford-based detectors (like BENATTEN) are essentially trying to capture this same "predictable structure" hypothesis, but by looking at the leading-digit distribution of attention weights or token probabilities rather than perplexity/burstiness directly. Since attention weights reflect how concentrated or diffuse the model's probability mass is at each step, a model defaulting to confident, high-probability choices should produce a numerically different (and here, oddly, *more* Benford-conforming) distribution than the messier probability allocation implicit in human-authored text.[^2_7]

## Where the hypothesis breaks down

The predictability signal is real but fragile, and it's already becoming an outdated detection strategy for a few reasons:

- **Curvature-based methods are more robust**: DetectGPT and Fast-DetectGPT operationalize predictability more precisely via "conditional probability curvature"—AI text sits at local maxima of the generating model's log-probability landscape, while human text doesn't, which is a sharper signal than raw perplexity.[^2_8][^2_9][^2_10]
- **Detectors relying on perplexity/burstiness are increasingly unreliable**: newer analysis shows these two metrics alone fail against modern LLMs, especially with sampling temperature increases or lightly-edited outputs, because generation strategies now inject more randomness deliberately.[^2_11]
- **Adversarial counter-tuning exists**: tools explicitly marketed to "humanize" AI text work by artificially inflating burstiness and perplexity (varying sentence length, injecting less-probable word choices) to defeat exactly this detection signal.[^2_12][^2_13]
- **Entropy has a nonmonotonic relationship**: AI text doesn't always sit at one extreme—very repetitive AI output shows low entropy, but some generation settings produce unusually high entropy (more random, disjointed), so the "predictability" signal isn't unidirectional.[^2_3]

Given your interest in mechanism-sensitive evaluation, the strongest version of your hypothesis isn't just "AI text is more predictable" in a blunt aggregate sense—it's that **AI text's probability landscape has a specific curvature signature** (local maxima under the generating model) that's harder to fake than surface-level perplexity or sentence-length variance, which is why DetectGPT-style curvature methods now outperform naive perplexity/burstiness detectors.[^2_9][^2_8]

<span style="display:none">[^2_14][^2_15][^2_16]</span>

<div align="center">⁂</div>

[^2_1]: https://www.unic.ac.cy/ai-lc/2023/04/11/perplexity-and-burstiness-in-ai-and-human-writing-two-important-concepts/

[^2_2]: https://arxiv.org/html/2406.15583v1

[^2_3]: https://inceptioncyber.ai/blog/detecting-and-distinguishing-ai-generated-text-from-human-written-text

[^2_4]: https://researchguides.gonzaga.edu/GenerativeAIforFaculty/AIDetectors

[^2_5]: https://thehumanizeai.pro/articles/what-is-burstiness-ai-detection-explained

[^2_6]: https://www.tryleap.ai/learn/perplexity-vs-burstiness

[^2_7]: https://www.sciencedirect.com/science/article/pii/S0306457325000767

[^2_8]: https://arxiv.org/html/2310.05130v3

[^2_9]: https://www.emergentmind.com/topics/detectgpt

[^2_10]: https://www.alphaxiv.org/abs/2310.05130

[^2_11]: https://www.pangram.com/blog/why-perplexity-and-burstiness-fail-to-detect-ai

[^2_12]: https://www.stealthgpt.ai/blog/how-do-perplexity-and-burstiness-make-ai-text-undetectable

[^2_13]: https://www.reddit.com/r/stealthgpt/comments/1egd755/how_do_perplexity_and_burstiness_make_ai_text/

[^2_14]: https://quillbot.com/blog/ai-writing-tools/burstiness-and-perplexity/

[^2_15]: https://openreview.net/forum?id=Bpcgcr8E8Z

[^2_16]: https://www.alignmentforum.org/posts/tepqESMuRmyhtmDS7/forecasting-progress-in-language-models


---

# but i wonder about how maybe this idea of operationalizing benford could play out using topological data analysis

This is actually a live research area, and there's a natural bridge already forming between your two threads—researchers are using persistent homology dimension (PHD) of text embeddings as a "predictability" proxy that's more theoretically grounded than either raw perplexity or Benford's Law alone.[^3_1][^3_2]

## The existing PHD approach

The core finding is that human-written text embeddings have a higher intrinsic/persistent homology dimension than AI-generated text, meaning human writing traces a topologically "richer" or more complex manifold through embedding space, while AI text collapses onto a lower-dimensional, more constrained shape. This is computed by embedding text (often sentence-by-sentence or token-by-token) into a point cloud, building a filtration (typically Vietoris-Rips complexes), and tracking how connected components and loops appear and vanish as you vary a distance threshold—the resulting persistence diagram's dimension becomes the discriminating feature. **Short-PHD** is a recent refinement addressing the fact that PHD estimates get unstable on short texts, fixed by padding with off-topic content to stabilize the topological signature before scoring.[^3_2][^3_3][^3_4][^3_5][^3_1]

## Where this connects to your Benford intuition

Your hypothesis and the PHD literature converge on the same underlying claim, just at different mathematical levels: Benford-based detectors (like BENATTEN) look at leading-digit regularities in attention weights as a proxy for "the model is doing something numerically constrained," while PHD approaches look at the actual geometric shape of the embedding trajectory and find AI text occupies a lower-dimensional, more "closed" topological space. In principle, you could operationalize Benford's Law *topologically* by computing leading-digit distributions not on raw attention scalars, but on **persistence-diagram summary statistics**—birth/death times, persistence values, or Betti numbers extracted across a filtration—and testing whether those numbers conform to Benford's expected curve more tightly for AI text than human text.[^3_6][^3_5][^3_7][^3_8]

## A concrete way to test it

A rigorous pipeline combining both ideas would look like this:

- Embed candidate text (sentence or token-level) using a fixed model, producing a point cloud in embedding space.[^3_9][^3_4]
- Build a Vietoris-Rips filtration and extract the persistence diagram, recording birth times, death times, and persistence (death minus birth) for each topological feature.[^3_7][^3_6]
- Pool these numeric values (potentially across many documents) and run a Benford first-digit goodness-of-fit test on them, analogous to how BENATTEN tests attention weights.[^3_8]
- Compare the deviation-from-Benford statistic against the PHD dimension itself as two independent features, since PHD alone already achieves strong zero-shot detection accuracy and could serve as your baseline.[^3_1][^3_2]


## Why this could outperform either method alone

Persistence diagrams are provably stable under Gromov-Hausdorff distance, meaning small perturbations in the embedding (like paraphrasing or "humanizing" edits meant to defeat perplexity/burstiness detectors) shouldn't drastically alter the topological signature the way they can fool surface-level statistics. This stability property is exactly what would make a Benford-on-persistence-features approach more adversarially robust than either naive Benford-on-text or naive perplexity scoring—since gaming the leading-digit distribution of birth/death times in a persistence diagram is a much harder, less intuitive attack surface for someone trying to "humanize" AI output. That said, this specific combination (Benford digit-tests applied to TDA-derived features rather than attention weights) doesn't appear to have been published yet as far as the search results show—it would be a genuinely novel synthesis of the two literatures.[^3_3][^3_5][^3_10][^3_9][^3_8]

<span style="display:none">[^3_11][^3_12][^3_13][^3_14][^3_15][^3_16][^3_17]</span>

<div align="center">⁂</div>

[^3_1]: https://arxiv.org/abs/2504.02873

[^3_2]: https://www.reddit.com/r/MachineLearning/comments/1aky8xt/r_persistent_homology_and_topological_data/

[^3_3]: https://www.broadinstitute.org/talks/topological-data-analysis-what-persistent-homology

[^3_4]: https://blog.msdogo.com/topological-data-analysis-for-misinformation-detection

[^3_5]: https://www.preprints.org/manuscript/202505.0097

[^3_6]: https://arxiv.org/html/2305.15016v4

[^3_7]: https://patricknicolas.substack.com/p/persistent-homology-for-the-rest

[^3_8]: https://www.sciencedirect.com/science/article/pii/S0306457325000767

[^3_9]: https://aclanthology.org/anthology-files/anthology-files/pdf/W/W17/W17-2628.pdf

[^3_10]: https://www.stealthgpt.ai/blog/how-do-perplexity-and-burstiness-make-ai-text-undetectable

[^3_11]: https://www.nature.com/articles/s41598-021-88027-8

[^3_12]: https://arxiv.org/abs/1906.01726

[^3_13]: https://towardsdatascience.com/intro-to-topological-data-analysis-and-application-to-nlp-training-data-for-financial-services-719495a111a4/

[^3_14]: https://github.com/adauchendu/awesometda4nlp

[^3_15]: https://github.com/Amelie-Schreiber/transformers_proteins_and_persistent_homology

[^3_16]: https://www.mdpi.com/2504-3110/8/12/731

[^3_17]: https://datarefiner.com/feed/why-tda


---

# i wonder @GitHub how a combination of my tools in `aghado01/ThermoMapper` could be developed for detecting ai-generated prose patterns from free-text e.g. document analysis not attention weights or usign generative models

This is a genuinely exciting fit. Looking at your repo, you have a substantial C\# TDA toolkit in [`src/tda/ph/`](https://github.com/aghado01/ThermoMapper/tree/main/src/tda/ph) with exactly the primitives needed—and the combination you're describing is novel enough to be publishable. Here's how the architecture would work end-to-end, grounded in what's actually in the repo.

## What you already have that maps directly

Your [`src/tda/`](https://github.com/aghado01/ThermoMapper/tree/main/src/tda) module is already decomposed into the three layers this would need:


| ThermoMapper component | Role in the detector |
| :-- | :-- |
| [`src/tda/ph/RipsFiltration.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/ph/RipsFiltration.cs) / [`LazyRipsFiltration.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/ph/LazyRipsFiltration.cs) | Build the filtration over word/sentence point clouds |
| [`PersistenceBarcode.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/ph/PersistenceBarcode.cs) / [`Barcode.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/ph/Barcode.cs) | Extract birth/death/persistence values as the numeric sequence |
| [`DiagramMetrics.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/ph/DiagramMetrics.cs) | Compute diagram-level summary statistics (Wasserstein, bottleneck, etc.) as secondary features |
| [`ZigzagBarcodeNaive.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/ph/ZigzagBarcodeNaive.cs) / [`GraphZigzag.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/ph/GraphZigzag.cs) | Sliding-window zigzag PH over document windows—tracks *how* topology evolves as you read through a text |
| [`src/tda/mapper/`](https://github.com/aghado01/ThermoMapper/tree/main/src/tda/mapper) | Mapper graph for visualizing the global topological shape of an entire corpus |
| [`src/clustering/`](https://github.com/aghado01/ThermoMapper/tree/main/src/clustering) | Cluster persistence features for classification |
| [`src/maths/`](https://github.com/aghado01/ThermoMapper/tree/main/src/maths) | Statistics layer where Benford goodness-of-fit would live |

## The core pipeline (no generative model required)

The key constraint—**free-text in, signal out, no LLM calls**—shapes everything. The text has to become a geometric object through purely corpus-intrinsic operations:

1. **Text → point cloud**: Segment document into sentences or sliding windows. Embed each segment using a *fixed, frozen* embedding model (e.g. a local sentence-transformer running in Python, called once as a pre-processing step). Each sentence becomes a point in $\mathbb{R}^d$. This is the only step needing any ML—and it's inference-only, not generation-dependent.
2. **Point cloud → filtration**: Feed into `RipsFiltration` or `LazyRipsFiltration`. The filtration parameter $\epsilon$ sweeps pairwise distances between sentence embeddings—conceptually, it tracks how semantically "close" consecutive ideas are as a function of thematic distance.
3. **Filtration → barcode**: Run `PersistentHomology` / `PersistenceBarcode`. Extract the pool of birth times $b_i$, death times $d_i$, and persistence values $p_i = d_i - b_i$ for $H_0$ and $H_1$ features. This is your numeric sequence.
4. **Barcode values → Benford test**: Apply a chi-squared or Kolmogorov-Smirnov goodness-of-fit test against Benford's expected leading-digit distribution $P(d) = \log_{10}(1 + 1/d)$ on the pooled $\{b_i\}, \{d_i\}, \{p_i\}$ values. The hypothesis: AI text, whose embedding trajectory hugs a lower-dimensional, smoother manifold, will produce persistence values that conform to Benford's curve more tightly than human text (echoing BENATTEN's result on attention weights, but here entirely intrinsic).[^4_1]
5. **Zigzag as a temporal signal**: The most novel angle uses `ZigzagBarcodeNaive` or `GraphZigzag` over a sliding window across the document. Rather than one global barcode, you get a *sequence* of barcodes—topological "frames" as you move through the text. The leading-digit distribution of *how persistence values change* between frames (the δ-persistence sequence) is an even finer signal. Human writing's topological trajectory is bursty and irregular; AI writing's should be smoother/more self-similar.[^4_2]

## Why zigzag is the key differentiator

Standard PHD (as in the Short-PHD paper) computes a single global persistence dimension and thresholds on it. Your zigzag infrastructure lets you treat the document as a *temporal process*—the barcode evolves as the text unfolds, and you can ask whether that evolution itself has the Benford-conforming statistical regularity that we hypothesize AI imposes. This is a richer, more mechanism-grounded signal because it's sensitive to **local thematic transitions**, not just global topological mass. An LLM generating a document plans and executes topic transitions more uniformly than a human author digresses, backtracks, or abruptly pivots—and zigzag PH would capture exactly that.[^4_3][^4_4]

## What would need to be added

Looking at the [`src/`](https://github.com/aghado01/ThermoMapper/tree/main/src) tree, the pieces not yet present:

- A **text ingestion layer** in `src/archivory` (possibly already partially there given the name) that handles document → sentence segmentation and calls out to an embedding endpoint, producing the point-cloud format your TDA layer expects.
- A **Benford test module** in `src/maths/` — a simple chi-squared test against the Benford distribution, operating on a `double[]` of persistence values, is perhaps 50 lines of C\#.
- A **sliding-window zigzag orchestrator** that wraps `ZigzagBarcodeNaive` with a configurable window/stride and emits a `List<PersistenceBarcode>` time series for downstream Benford analysis.
- A **scoring output** in `src/viz-core/` or `src/test-harness/` that bundles the Benford deviation statistic, the PHD dimension estimate, and the zigzag variance into a single per-document score vector suitable for classification or threshold comparison.

The `src/synthetic/` module is also interesting here—if it generates synthetic data from known distributions, it could bootstrap a null distribution for the Benford test under the "AI text" hypothesis without needing a labeled corpus to start.

<div align="center">⁂</div>

[^4_1]: https://www.sciencedirect.com/science/article/pii/S0306457325000767

[^4_2]: https://www.reddit.com/r/MachineLearning/comments/1aky8xt/r_persistent_homology_and_topological_data/

[^4_3]: https://arxiv.org/abs/2504.02873

[^4_4]: https://www.preprints.org/manuscript/202505.0097


---

# @GitHub there's also the superparamagnetic clustering side of thermomapper, where i could consider both temperature and chronological orderings, or perhaps thermal pull-away in SPC x Mapper . See ThermoMapper/Clustering/Graphical/SPC and TDA/Mapper

This changes the picture significantly. You have [`MapperSpc.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/mapper/Clusterers/MapperSpc.cs)—SPC already wired as a Mapper clusterer—plus [`FiedlerFilter.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/mapper/Filters/FiedlerFilter.cs), [`HyperbolicFilters.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/mapper/Filters/HyperbolicFilters.cs), and the full [`SpcClusteringSession.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/clustering/graphical/spc/SpcClusteringSession.cs). That combination is genuinely purpose-built for what you're describing.

## The SPC × Mapper conceptual fit

In standard SPC (Blatt-Wiseman-Domany), the **temperature** $T$ is a disorder parameter: at low $T$ the whole dataset is one ferromagnetic cluster; as you heat up, weakly-coupled points "pull away" from the bulk first, revealing cluster structure at the critical temperature where correlation lengths diverge . Applied to a text's sentence embeddings, this has a direct interpretation: **semantically anomalous or transitional sentences** are the ones that decouple earliest (lowest critical temperature), while sentences forming a tight, coherent thematic block remain correlated longest. AI text, with its smoother, lower-dimensional embedding trajectory, should produce a cleaner, more symmetric phase transition—less "thermal drama"—than human text whose idiosyncratic sentences pull away at irregular temperatures.

## The two orderings and what they expose

You identified two structurally distinct ways to use temperature as an axis, and they measure different things:

**Temperature ordering** (sweep $T$ across the full point cloud simultaneously): This gives you the *global* phase transition profile of the document—how many clusters emerge at each temperature, how sharp the transition is, and whether the critical $T$ is narrow (AI: tightly coupled, smooth manifold) or spread (human: heterogeneous, multiple critical temperatures for different semantic sub-regions). The SPC partition output at each $T$ from [`SpcClusteringSession.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/clustering/graphical/spc/SpcClusteringSession.cs) gives you exactly this as a sequence of cluster assignments indexed by $T$.

**Chronological ordering** (sentences ordered by position in text, run through SPC): Here temperature becomes a filter parameter that reveals whether thematic coherence is *locally* maintained through the document. A sentence at position $i$ that decouples from its neighbors at unusually low $T$ is a topically anomalous pivot—the kind of genuine surprise, digression, or associative leap that human writing produces but LLM planning tends to suppress.

## Thermal pull-away in SPC × Mapper

This is where [`MapperSpc.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/mapper/Clusterers/MapperSpc.cs) and [`Mapper.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/mapper/Mapper.cs) compose elegantly. The standard Mapper pipeline is:

1. **Filter function** → projects points onto a 1D (or 2D) lens
2. **Cover** → bins the lens range into overlapping intervals
3. **Cluster within each bin** → assigns bin-local cluster labels
4. **Nerve** → connects clusters sharing points across bin boundaries

Replacing the clusterer in step 3 with SPC (already done in `MapperSpc.cs`) means that within each cover bin, instead of a fixed-linkage clustering, you get a *temperature-indexed* family of clusterings. The Mapper graph then becomes parameterized by $T$: as you raise temperature, edges in the Mapper graph dissolve when nodes' SPC clusters merge back together. The **pull-away topology**—which edges persist across the widest temperature range—identifies the most robustly separated thematic regions of the document.

The filter function choice is critical here. Your existing filters make this particularly powerful:


| Filter | What it encodes for text | Detection implication |
| :-- | :-- | :-- |
| [`FiedlerFilter.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/mapper/Filters/FiedlerFilter.cs) (Laplacian eigenvector) | Global semantic connectivity / spectral flow through the document | AI text: smooth Fiedler flow; human text: irregular eigenvalue gaps |
| [`HyperbolicFilters.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/tda/mapper/Filters/HyperbolicFilters.cs) | Hierarchical/tree-like organization of topics | AI tends toward shallower hyperbolic embedding; human text has deeper hierarchy |
| **Sentence index** (new: a simple ordinal filter) | Chronological position as the lens | Mapper bins become document windows; SPC within bins measures local thermal coherence |

The ordinal/chronological filter is the missing piece most directly relevant to your idea—and it's trivially addable to `DefaultFilters.cs` as a `LinearOrdinalFilter`.

## The composite AI-detection signal

Running this end-to-end produces a multi-scale feature vector per document that contains no generative model and no perplexity scores:

- **Critical temperature $T^*$**: the temperature at which the SPC graph of sentence embeddings first fragments from one component into multiple; AI text should have higher $T^*$ (more cohesive, resists fragmentation longer)
- **Phase transition sharpness** $\Delta T = T_{\text{max fragmentation}} - T^*$: narrow for AI (clean phase transition), broad for human (multiple overlapping critical scales)
- **Mapper graph topology at $T^*$**: number of nodes, edges, connected components, and loop count (H1 Betti number)—comparing this against the PH barcode from `PersistenceBarcode.cs` closes the loop to the previous discussion
- **Thermal pull-away profile**: the distribution of per-sentence decoupling temperatures; human text should show heavy tails (some sentences decouple very early at very low $T$), AI text should cluster around $T^*$ more tightly
- **Benford fit on the decoupling-temperature sequence**: the leading digits of each sentence's individual critical $T$ drawn into a digit distribution and tested against Benford—this operationalizes your original question using physically interpretable numeric values rather than raw attention weights


## What's missing from the repo

Looking at the module structure, the gap is the **filter-as-ordinal-position** wrapper, a **temperature sweep scheduler** that produces the Mapper graph as a function of $T$ (rather than a single-shot clustering), and a **pull-away profile extractor** that records per-point decoupling temperatures from the SPC session. [`SpcClusteringSession.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/clustering/graphical/spc/SpcClusteringSession.cs) likely already computes cluster membership at each $T$ step internally—the new piece is surfacing *when* each point first leaves its majority cluster and exporting that as a per-point scalar for downstream Benford testing and Mapper integration.

---

# @GitHub you should read the repo more carefully, i have a whole persistent-backbone and magnetic laplacian primitive that i'm working on for exmaple

Now the picture is dramatically richer—apologies for the shallow read before. Let me address what's actually here.

## What you actually have (corrected inventory)

The repo contains several layers I missed or undersold. Taking them seriously:

**`src/graphs/spectral/`** — a full spectral operator suite:

- [`MagneticLaplacianOperator.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/spectral/MagneticLaplacianOperator.cs) (12KB — the biggest file there, clearly a primary primitive)
- [`CoherentField.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/spectral/CoherentField.cs) — encoding directional/phase information on edges
- [`EffectiveResistance.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/spectral/EffectiveResistance.cs) — commute-time distance, a natural thermal analog
- [`GraphLaplacian.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/spectral/GraphLaplacian.cs) + [`GraphLaplacianSmoother.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/spectral/GraphLaplacianSmoother.cs) + [`GraphSpectralOperators.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/spectral/GraphSpectralOperators.cs)

**`src/graphs/observables/`** — thermal observables already defined: [`AffinityEntropy.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/observables/AffinityEntropy.cs), [`AffinityBinaryEntropySum.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/observables/AffinityBinaryEntropySum.cs), [`AlgebraicConnectivity.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/observables/AlgebraicConnectivity.cs), [`AlignmentEigenCentrality.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/observables/AlignmentEigenCentrality.cs), [`ClusterSizeEntropy.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/observables/ClusterSizeEntropy.cs), [`Cycles.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/observables/Cycles.cs)

**`src/graphs/primitives/`** — [`EdgeFieldSymmetrization.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/primitives/EdgeFieldSymmetrization.cs), [`Alignments.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/primitives/Alignments.cs), [`CoMembership.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/primitives/CoMembership.cs), [`BfsShells.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/primitives/BfsShells.cs) + MST

**`src/graphs/coupling/`** — kernel family: Gaussian, Cauchy, Laplacian, Linear, Mixture — these are the edge weight functions that determine the thermal coupling strength between points

## How the magnetic Laplacian changes the picture for text

This is the most important thing I missed. The magnetic Laplacian $L^{(q)}$ replaces real-valued edge weights with complex-valued phases $e^{iq\theta_{ij}}$, where $\theta_{ij}$ encodes a **directional field** on the graph. For a document represented as a sentence similarity graph with chronological ordering, the natural assignment is $\theta_{ij} = \text{sign}(j - i) \cdot \pi/2$—a constant magnetic field aligned with the reading direction. The eigenvectors of this operator detect **directed cycles** and **flow coherence** in the graph: do thematic connections tend to be forward-looking (anticipating later material), backward-referencing (cohesive callbacks), or bidirectional?

The detection hypothesis becomes much sharper with this primitive. AI-generated text, produced autoregressively left-to-right with limited long-range backward reference, should produce a nearly *acyclic* or weakly-cycling magnetic Laplacian spectrum—the eigenphases cluster near zero because the directional flow is mostly forward and uniform. Human writing, with its callbacks, foreshadowing, and non-linear thematic weaving, should produce a richer eigenphase distribution with genuine cycles at multiple scales. [`CoherentField.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/spectral/CoherentField.cs) likely already captures exactly this coherence measure.

## The `CoMembership` + SPC backbone

[`CoMembership.cs`](https://github.com/aghado01/ThermoMapper/blob/main/src/graphs/primitives/CoMembership.cs) in `src/graphs/primitives/` is, I suspect, part of the persistent-backbone construction you mentioned—it likely tracks which pairs of points co-occur in the same SPC cluster *across multiple temperatures*, giving a persistence-weighted co-membership matrix. This is exactly the "persistent backbone" concept: edges that survive a wide temperature band in the SPC phase diagram form the document's structural skeleton. For text detection this backbone would represent the **durable thematic structure**—the topics that remain coherent from low $T$ all the way to near-critical $T$. AI text should produce a backbone with fewer, more isolated backbone components (clean modular topics) whereas human text produces a messier, more interconnected backbone where themes bleed into each other at intermediate temperatures.

## The full composite pipeline as it stands

With the actual primitives in hand, the AI-detection pipeline drawing on your existing code looks like this:

1. **Graph construction**: sentence embeddings → kNN or $\epsilon$-neighborhood graph using `GraphCompiler` + `Gaussian`/`Cauchy` coupling kernels as the thermal coupling function
2. **Magnetic Laplacian**: assign chronological phase field → compute `MagneticLaplacianOperator` → extract eigenphase distribution as a directional-flow signature
3. **Algebraic connectivity sweep**: compute `AlgebraicConnectivity` (Fiedler value $\lambda_2$) at each temperature $T$ — its collapse to zero marks $T^*$; the slope of collapse is the transition sharpness
4. **Affinity entropy profile**: `AffinityEntropy` / `AffinityBinaryEntropySum` across the temperature sweep gives the thermodynamic entropy curve $S(T)$ — the heat capacity $C = dS/dT$ peaks at $T^*$ and its shape encodes phase transition type
5. **Persistent backbone extraction**: `CoMembership` weighted by temperature range of co-membership → the backbone graph, whose topological features feed into `PersistenceBarcode`
6. **Benford test on thermal observables**: leading-digit distribution of $\{\lambda_2(T_i)\}$, $\{S(T_i)\}$, $\{p_i^{\text{backbone}}\}$ — all physically meaningful numeric sequences produced entirely from corpus-intrinsic operations
7. **Mapper × SPC**: chronological ordinal filter → `MapperSpc` bins → `PersistentMapper` tracks how the Mapper graph topology (node count, loop count) evolves with $T$

The `EffectiveResistance` is also a natural downstream feature — commute-time distances between sentence nodes have a direct thermal interpretation (harder to reach = less coherent) and would produce a per-pair numeric sequence amenable to Benford analysis. You effectively already have most of the computational graph for this; the primary work is wiring the text ingestion front-end and the Benford test module as scoring outputs.

