using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>What a weight states about its observation.</summary>
public enum WeightMeaning
{
    /// <summary>yᵢ ~ N(fᵢ, σ²/wᵢ): one observation at a known relative precision, what a latent scale or Pólya–Gamma hands the
    /// core. σ²'s shape counts the rows, and the likelihood carries ½Σ ln wᵢ.</summary>
    Precision,

    /// <summary>yᵢ counted wᵢ times: the likelihood is Πᵢ N(yᵢ | fᵢ, σ²)^wᵢ, so σ²'s shape counts Σwᵢ and integer weights are
    /// duplicated rows.</summary>
    Replication,
}

/// <summary>Which columns the coefficient prior reads under weights (D03's robust row).</summary>
public enum WeightedPriorGram
{
    /// <summary>The columns as they stand, N rows, free of the weights: the g-form's τG_c on the columns centred over N rows, and
    /// the ridge's D from them. The level and s² read the working response unweighted, so a latent scale's conditional stays
    /// in closed form, at a decomposition per structure (<see cref="PriorPencil"/>). D03's default under robust noise.</summary>
    Unweighted,

    /// <summary>The weighted columns: τZ_cᵀWZ_c and its D, with the level and s² weighted. The closed forms of the unweighted core
    /// hold on the weighted statistics, at O(1) per (τ, b), and integer replication weights are duplicated rows. A latent scale's
    /// conditional gains |Z_cᵀWZ_c|^½, since the prior moves with W.</summary>
    Weighted,

    /// <summary>The columns under a frozen W₀ (<see cref="BassPriorReference.GramWeights"/>), the information at a fixed point, free of
    /// the moving W: τZ_c₀ᵀW₀Z_c₀ with the columns centred under W₀, and the ridge's D from them. W₀ reads in information per
    /// observation: under the g-form W₀ = (1/κ)·1 is <see cref="BassPriorReference.CoefficientConstant"/> κ on the unweighted columns.
    /// At a decomposition per structure, as <see cref="Unweighted"/>'s.</summary>
    Frozen,
}

/// <summary>
/// What the priors read of the response, frozen apart from it, for a working response that moves (y − β√w·v, the latent z of a
/// probit, y − ρf_lo). A moving working response is not data: a prior that re-reads it is conditioned on parameters, with
/// conditionals that share no joint. So the observation model states here what the priors read, from the observed response, and the
/// target stays one target whatever moves later. What it leaves null the config must not read. The model names the rule, since the
/// read has to mean something for its noise:
/// <list type="bullet">
/// <item><see cref="Moments(ReadOnlySpan{double}, double)"/>: the mean and the mean square about it, D03's read.</item>
/// <item><see cref="About"/>: the mean square about a stated level, as the zero centre reads 0.</item>
/// <item><see cref="Robust"/>: the median and MAD², for heavy tails.</item>
/// <item><see cref="Quantile"/>: the sample q-quantile, for the asymmetric Laplace, whose f is the q-quantile.</item>
/// <item><see cref="FrozenGram"/>: moments under a frozen W₀, with the prior's Gram on it.</item>
/// <item><see cref="Link"/>: only the coefficient prior's constant, for a binary or count response, which has no location-scale
/// group for a read to serve.</item>
/// </list>
/// Each read is equivariant under y → ky + m, so the frozen route keeps BASS's equivariance. It uses the data twice, so an evidence
/// stays empirical-Bayes in what is read. The alternatives and their costs are in docs/Maths/BayesianRegression.md.
/// </summary>
public sealed class BassPriorReference
{
    private readonly double[]? _gramWeights;

    /// <summary>A reference that states <paramref name="level"/>, <paramref name="meanSquare"/> (non-negative),
    /// <paramref name="coefficientConstant"/> (positive) and <paramref name="gramWeights"/> (positive, one per observation; none leaves
    /// the Gram to <see cref="WeightedPriorGram"/>'s other readings), each finite. The weights are copied.</summary>
    public BassPriorReference(double? level = null, double? meanSquare = null, double coefficientConstant = 1.0, ReadOnlySpan<double> gramWeights = default)
    {
        if (level is double l && !double.IsFinite(l)) throw new ArgumentOutOfRangeException(nameof(level), "The level must be finite.");
        if (meanSquare is double s && !(s >= 0.0 && double.IsFinite(s))) throw new ArgumentOutOfRangeException(nameof(meanSquare), "The mean square must be non-negative and finite.");
        if (!(coefficientConstant > 0.0 && double.IsFinite(coefficientConstant))) throw new ArgumentOutOfRangeException(nameof(coefficientConstant), "The constant must be positive and finite.");
        if (!gramWeights.IsEmpty)
        {
            _gramWeights = gramWeights.ToArray();
            foreach (double w in _gramWeights)
                if (!(w > 0.0 && double.IsFinite(w))) throw new ArgumentException("A frozen weight is positive and finite.", nameof(gramWeights));
        }
        Level = level;
        MeanSquare = meanSquare;
        CoefficientConstant = coefficientConstant;
    }

    /// <summary>The mean of <paramref name="values"/> and their mean square about it, D03's read.</summary>
    public static BassPriorReference Moments(ReadOnlySpan<double> values, double coefficientConstant = 1.0)
    {
        double level = PairwiseSum.Of(values) / values.Length;
        return About(values, level, coefficientConstant);
    }

    /// <summary>The mean square of <paramref name="values"/> about <paramref name="level"/>, stated, as the zero centre reads 0.</summary>
    public static BassPriorReference About(ReadOnlySpan<double> values, double level, double coefficientConstant = 1.0)
    {
        var squares = new double[values.Length];
        for (int i = 0; i < values.Length; i++) squares[i] = (values[i] - level) * (values[i] - level);
        return new BassPriorReference(level, PairwiseSum.Of(squares) / values.Length, coefficientConstant);
    }

    /// <summary>
    /// The median of <paramref name="values"/> and (MAD/Φ⁻¹(¾))², the normal-consistent scale, for heavy tails: an outlier a latent
    /// scale discounts cannot move either. When over half the values tie, the MAD is 0 and says nothing of the spread, so the mean
    /// square about the median stands in.
    /// </summary>
    public static BassPriorReference Robust(ReadOnlySpan<double> values, double coefficientConstant = 1.0)
    {
        double median = Middle(values.ToArray());
        var deviations = new double[values.Length];
        for (int i = 0; i < values.Length; i++) deviations[i] = Math.Abs(values[i] - median);
        double scale = Middle(deviations) / Normal.StandardQuantile(0.75);
        return scale > 0.0 ? new BassPriorReference(median, scale * scale, coefficientConstant) : About(values, median, coefficientConstant);
    }

    /// <summary>
    /// The sample <paramref name="q"/>-quantile of <paramref name="values"/>, the ⌈qN⌉-th smallest, which minimizes the check loss, and
    /// the squared mean check loss: the asymmetric Laplace's maximum-likelihood level and scale, whose f is the q-quantile and whose
    /// ȳ is not its level.
    /// </summary>
    public static BassPriorReference Quantile(ReadOnlySpan<double> values, double q, double coefficientConstant = 1.0)
    {
        if (!(q > 0.0 && q < 1.0)) throw new ArgumentOutOfRangeException(nameof(q), "q lies in (0, 1).");
        double[] sorted = values.ToArray();
        Array.Sort(sorted);
        double level = sorted[Math.Max((int)Math.Ceiling(q * sorted.Length) - 1, 0)];
        var losses = new double[values.Length];
        for (int i = 0; i < values.Length; i++)
        {
            double e = values[i] - level;
            losses[i] = e < 0.0 ? (q - 1.0) * e : q * e;
        }
        double scale = PairwiseSum.Of(losses) / values.Length;
        return new BassPriorReference(level, scale * scale, coefficientConstant);
    }

    /// <summary>The mean of <paramref name="values"/> under W₀ = <paramref name="gramWeights"/> and their mean square under W₀, with the
    /// prior's Gram on W₀ (<see cref="WeightedPriorGram.Frozen"/>).</summary>
    public static BassPriorReference FrozenGram(ReadOnlySpan<double> values, ReadOnlySpan<double> gramWeights, double coefficientConstant = 1.0)
    {
        if (gramWeights.Length != values.Length) throw new ArgumentException($"The frozen weights hold {gramWeights.Length} values for {values.Length} observations.", nameof(gramWeights));
        double mass = PairwiseSum.Of(gramWeights), level = PairwiseSum.Dot(gramWeights, values) / mass;
        var squares = new double[values.Length];
        for (int i = 0; i < values.Length; i++) squares[i] = (values[i] - level) * (values[i] - level);
        return new BassPriorReference(level, PairwiseSum.Dot(gramWeights, squares) / mass, coefficientConstant, gramWeights);
    }

    /// <summary>
    /// Sabanés Bové and Held's (2011) constant alone, their c, κ = v(h(0))/h′(0)², from a link's response variance <paramref name="variance"/>
    /// and slope <paramref name="slope"/> at 0: 4 for the logit, π/2 for the probit, π²/4 for the cauchit, 1 for the Poisson log. It
    /// anchors the prior at an imaginary sample, so it reads no data.
    /// </summary>
    public static BassPriorReference Link(double variance, double slope) => new(coefficientConstant: variance / (slope * slope));

    /// <summary>The level centring reads; null states none.</summary>
    public double? Level { get; }

    /// <summary>s², which the scaled and hierarchical σ² priors read and D03's refusal and the exact-fit floor gate on; null states none.</summary>
    public double? MeanSquare { get; }

    /// <summary>κ in the coefficient prior's covariance σ²κ(τG)⁻¹ (Sabanés Bové and Held 2011's c), whatever the form: the unit-information
    /// prior for a noise whose information per observation is 1/κ. 1 is the Gaussian's.</summary>
    public double CoefficientConstant { get; }

    /// <summary>W₀, the prior's frozen Gram weights; empty when it holds none.</summary>
    public ReadOnlySpan<double> GramWeights => _gramWeights;

    /// <summary>W₀'s own array, not to be written; null when it holds none.</summary>
    internal double[]? GramWeightArray => _gramWeights;

    private static double Middle(double[] values)
    {
        Array.Sort(values);
        int half = values.Length / 2;
        return values.Length % 2 == 1 ? values[half] : 0.5 * (values[half - 1] + values[half]);
    }
}

/// <summary>
/// What an observation model hands the Gaussian core (F-X1, F-X6): the working response y − offset, the observation
/// weights W, the precisions the core's normal equations carry (ZᵀWZ, ZᵀWỹ), and, for a direct functional response, the
/// inner product Q on the output index that a residual correlated along it declares (covariance σ²(I ⊗ R), Q = R⁻¹).
/// Gaussian BASS hands the response itself, unit weights and Q = I, so its marginal and draws are the core's own. A
/// variance-mean mixture (generalized-hyperbolic noise) hands y − β√w·v with weights 1/(w·v) at the current latent state,
/// where a scale mixture hands weights alone.
/// <para>The core fits the working response under its weights (<see cref="BassDesign"/>), read as <see cref="Meaning"/> says,
/// with the coefficient prior on <see cref="PriorGram"/>'s columns, and the priors read <see cref="Reference"/> when the model states
/// one, the working response as it stands otherwise. A Q other than the identity is carried by the contract and
/// refused by today's core, which holds the independent Gram: the centring under Q and the marginal's Σᵢ yᵢᵀQyᵢ and log|R|
/// terms arrive with the extension that supplies them.</para>
/// </summary>
public sealed class WorkingResponse
{
    private readonly double[]? _shifted;
    private readonly double[]? _offset, _weights;
    private readonly double _responseMean;
    private readonly double[]? _outputMeans;   // several outputs' means; null for one response (E34)
    private double _responseSquares = double.NaN;   // ‖y − ȳ‖²_W, when first read; two readers compute the one value

    /// <summary>The working response of <paramref name="response"/> less <paramref name="offset"/> (one value per
    /// observation, case by case and location by location; none is zero), with <paramref name="weights"/> (positive, one per
    /// observation; none is unit), read as <paramref name="meaning"/> says with the coefficient prior on
    /// <paramref name="priorGram"/>'s columns, <paramref name="output"/>, the inner product on the output index (none is the
    /// identity; a direct functional response's only), and <paramref name="reference"/>, what the priors read in place of this
    /// working response (none reads it as it stands). <see cref="WeightedPriorGram.Frozen"/> takes its W₀ from the reference, and alone
    /// does; with no weights of its own it reads unit weights. <paramref name="latents"/> are the observation model's latents at this
    /// reading (one value per observation), none for a model with none. <paramref name="logJacobian"/> is the log-Jacobian of a linear map
    /// the offset applied to the response, which the likelihood carries beside the weights' Σ ln wᵢ: a random-intercept block's −ln|V| for
    /// its whitening (<see cref="RandomInterceptBlock.Whitened"/>), 0 otherwise. The arrays are copied.</summary>
    public WorkingResponse(BassResponse response, ReadOnlySpan<double> offset = default, ReadOnlySpan<double> weights = default,
                           IInnerProduct? output = null, WeightMeaning meaning = WeightMeaning.Precision,
                           WeightedPriorGram priorGram = WeightedPriorGram.Unweighted, BassPriorReference? reference = null,
                           ObservationLatents? latents = null, double logJacobian = 0.0)
    {
        if (!double.IsFinite(logJacobian)) throw new ArgumentOutOfRangeException(nameof(logJacobian));
        ArgumentNullException.ThrowIfNull(response);
        if (meaning is not (WeightMeaning.Precision or WeightMeaning.Replication)) throw new ArgumentOutOfRangeException(nameof(meaning));
        if (priorGram is not (WeightedPriorGram.Unweighted or WeightedPriorGram.Weighted or WeightedPriorGram.Frozen)) throw new ArgumentOutOfRangeException(nameof(priorGram));
        int n = response.Observations;
        if ((priorGram == WeightedPriorGram.Frozen) != reference?.GramWeightArray is not null)
            throw new ArgumentException("The frozen Gram and the reference's W₀ come together.", nameof(priorGram));
        if (reference?.GramWeightArray is { } frozen && frozen.Length != n)
            throw new ArgumentException($"The reference's W₀ holds {frozen.Length} values for {n} observations.", nameof(reference));
        if (latents is not null && latents.Values.Length != n)
            throw new ArgumentException($"The latents hold {latents.Values.Length} values for {n} observations.", nameof(latents));
        Response = response;
        Meaning = meaning;
        PriorGram = priorGram;
        Reference = reference;
        Latents = latents;
        Weighting = GPriorWeights.Unit(n);
        Output = output ?? IdentityInnerProduct.Instance;
        if (!Output.IsIdentity)
        {
            if (response is not FunctionalResponse) throw new ArgumentException("A scalar response has no output index.", nameof(output));
            if (Output.Dimension != response.Locations)
                throw new ArgumentException($"The inner product acts on {Output.Dimension} output locations, the response has {response.Locations}.", nameof(output));
        }
        (_offset, _shifted) = Shifted(response, offset);
        if (!weights.IsEmpty)
        {
            if (weights.Length != n) throw new ArgumentException($"The weights hold {weights.Length} values for {n} observations.", nameof(weights));
            _weights = weights.ToArray();
            var logs = new double[n];
            for (int i = 0; i < n; i++)
            {
                double w = _weights[i];
                if (!(w > 0.0 && double.IsFinite(w))) throw new ArgumentException("A weight is positive and finite.", nameof(weights));
                logs[i] = Math.Log(w);
            }
            double mass = PairwiseSum.Of(_weights);
            Weighting = meaning == WeightMeaning.Replication ? new GPriorWeights(mass, mass, 0.0) : new GPriorWeights(n, mass, PairwiseSum.Of(logs));
        }
        else if (priorGram == WeightedPriorGram.Frozen)
        {
            _weights = new double[n];
            Array.Fill(_weights, 1.0);
        }
        if (logJacobian != 0.0) Weighting = Weighting with { SumLogWeights = Weighting.SumLogWeights + logJacobian };
        if (response is MultiOutputResponse outputs)
        {
            // Several outputs (E34) read the response itself in slice 1: an offset, weights, latents, a frozen reference and a log-Jacobian
            // arrive with the noise and the random intercepts (slice 4).
            if (!offset.IsEmpty || !weights.IsEmpty || latents is not null || reference is not null || logJacobian != 0.0 || priorGram != WeightedPriorGram.Unweighted)
                throw new NotSupportedException("Several outputs read the response itself: offsets, weights, latents and a frozen reference arrive with E34's noise (slice 4).");
            _outputMeans = new double[outputs.Outputs];
            ReadOnlySpan<double> values = outputs.Values;
            int k = outputs.Outputs;
            var column = new double[n];
            for (int l = 0; l < k; l++)
            {
                for (int i = 0; i < n; i++) column[i] = values[i * k + l];
                _outputMeans[l] = PairwiseSum.Of(column) / n;
            }
            _responseMean = double.NaN;
            return;
        }
        _responseMean = MeanOf(_weights, Values, Weighting, n);
    }

    // `source` at another offset with `latents`: its weights shared, every other part as it stands.
    private WorkingResponse(WorkingResponse source, ReadOnlySpan<double> offset, ObservationLatents? latents)
    {
        if (latents is not null && latents.Values.Length != source.Response.Observations)
            throw new ArgumentException($"The latents hold {latents.Values.Length} values for {source.Response.Observations} observations.", nameof(latents));
        Response = source.Response;
        Meaning = source.Meaning;
        PriorGram = source.PriorGram;
        Reference = source.Reference;
        Latents = latents;
        Weighting = source.Weighting;
        Output = source.Output;
        Whitening = source.Whitening;
        _weights = source._weights;
        if (Response is MultiOutputResponse)
            throw new NotSupportedException("Several outputs read the response itself: offsets and latents arrive with E34's noise (slice 4).");
        (_offset, _shifted) = Shifted(Response, offset);
        _responseMean = MeanOf(_weights, Values, Weighting, Response.Observations);
    }

    private static (double[]? Offset, double[]? Shifted) Shifted(BassResponse response, ReadOnlySpan<double> offset)
    {
        if (offset.IsEmpty) return (null, null);
        int n = response.Observations;
        if (offset.Length != n) throw new ArgumentException($"The offset has {offset.Length} values for {n} observations.", nameof(offset));
        double[] copied = offset.ToArray();
        var shifted = new double[n];
        ReadOnlySpan<double> y = response.Values;
        for (int i = 0; i < n; i++)
        {
            if (!double.IsFinite(copied[i])) throw new ArgumentException("The offset has a non-finite value.", nameof(offset));
            shifted[i] = y[i] - copied[i];
        }
        return (copied, shifted);
    }

    private static double MeanOf(double[]? weights, ReadOnlySpan<double> values, GPriorWeights weighting, int n)
        => weights is null ? PairwiseSum.Of(values) / n : PairwiseSum.Dot(weights, values) / weighting.Mass;

    /// <summary>
    /// This reading of the response at another latent state: <paramref name="offset"/> (none is zero), <paramref name="weights"/>
    /// (none keeps these, shared rather than copied) and the latents there, <paramref name="latents"/> (none carries none), with the
    /// meaning, the prior's Gram, the inner product and the reference as they stand, since those are the model's. A chain moves to it by
    /// <see cref="BassStructureMoves.Follow"/>. A design's layout depends on whether it holds weights, so a model whose weights move declares
    /// them from the start, unit weights if need be.
    /// </summary>
    public WorkingResponse Moved(ReadOnlySpan<double> offset, ReadOnlySpan<double> weights = default, ObservationLatents? latents = null)
    {
        if (weights.IsEmpty) return new WorkingResponse(this, offset, latents);
        if (_weights is null) throw new ArgumentException("This working response holds no weights, and a design's layout depends on whether it does; declare them from the start.", nameof(weights));
        return new WorkingResponse(Response, offset, weights, Output, Meaning, PriorGram, Reference, latents);
    }

    /// <summary>The response the working response derives from.</summary>
    public BassResponse Response { get; }

    /// <summary>
    /// The random-effect block whose whitening V^(−½) this reading's offset applies, its ln|V| in the log-Jacobian
    /// (<see cref="RandomInterceptBlock.Whitened"/>); null where none does. The design whitens each column by the block of the reading it
    /// is on, so a reading at another scale moves every column as moved weights do (E21).
    /// </summary>
    public RandomInterceptBlock? Whitening { get; internal init; }

    /// <summary>y − offset, case by case and location by location; the response itself when there is no offset.</summary>
    public ReadOnlySpan<double> Values => _shifted ?? Response.Values;

    /// <summary>The offset; empty when there is none, which is the zero offset.</summary>
    public ReadOnlySpan<double> Offset => _offset;

    /// <summary>The weights W; empty when there are none, which is unit weights.</summary>
    public ReadOnlySpan<double> Weights => _weights;

    /// <summary>Whether the weights are unit.</summary>
    public bool IsUnweighted => _weights is null;

    /// <summary>What the weights state; unit weights state nothing.</summary>
    public WeightMeaning Meaning { get; }

    /// <summary>Which columns the coefficient prior reads under weights; unit weights make the first two one.</summary>
    public WeightedPriorGram PriorGram { get; }

    /// <summary>What the priors read in place of this working response; null reads it as it stands.</summary>
    public BassPriorReference? Reference { get; }

    /// <summary>The observation model's latents at this reading, the state's (E1's decision 1); null for a model with none.</summary>
    public ObservationLatents? Latents { get; }

    /// <summary>The count, the mass and Σ ln wᵢ the marginal reads beyond the weighted statistics.</summary>
    public GPriorWeights Weighting { get; }

    /// <summary>Whether a coefficient prior apart from the likelihood's Gram is needed: weights, with the prior free of them.</summary>
    public bool PriorStandsApart => _weights is not null && PriorGram != WeightedPriorGram.Weighted;

    /// <summary>ȳ, the working response's mean under W over all N observations: Σwᵢyᵢ/Σwᵢ, the mean at unit weights. Several outputs
    /// have one each (<see cref="OutputMeans"/>), and refuse this.</summary>
    public double ResponseMean => _outputMeans is null ? _responseMean
        : throw new InvalidOperationException("Several outputs have a mean each: read OutputMeans.");

    /// <summary>K, the outputs (<see cref="BassResponse.Outputs"/>).</summary>
    public int Outputs => Response.Outputs;

    /// <summary>Each output's mean over the N observations, for several outputs (<see cref="MultiOutputResponse"/>); empty otherwise.</summary>
    public ReadOnlySpan<double> OutputMeans => _outputMeans;

    /// <summary>Whether this is the response itself under unit weights: no offset, no weights.</summary>
    public bool IsResponseItself => _offset is null && _weights is null;

    /// <summary>‖y − ȳ‖²_W, the working response's squares about its mean under W, which the exact-fit floor reads once a chain has
    /// moved; taken when first read.</summary>
    internal double ResponseSquares
    {
        get
        {
            if (!double.IsNaN(_responseSquares)) return _responseSquares;
            ReadOnlySpan<double> y = Values;
            var squares = new double[y.Length];
            for (int i = 0; i < y.Length; i++) squares[i] = (y[i] - ResponseMean) * (y[i] - ResponseMean);
            return _responseSquares = _weights is null ? PairwiseSum.Of(squares) : PairwiseSum.Dot(_weights, squares);
        }
    }

    /// <summary>The weights' own array, not to be written; null at unit weights.</summary>
    internal double[]? WeightArray => _weights;

    /// <summary>A column's mean under W, the design's centring, from its case side <paramref name="a"/> and its output side
    /// <paramref name="h"/> for a direct response: āh̄ at unit weights.</summary>
    internal double ColumnMean(double[] a, double[]? h)
    {
        if (_weights is null)
        {
            double mean = PairwiseSum.Of(a) / a.Length;
            return h is null ? mean : mean * (PairwiseSum.Of(h) / h.Length);
        }
        if (h is null) return PairwiseSum.Dot(_weights, a) / Weighting.Mass;
        var byCase = new double[a.Length];
        for (int i = 0; i < a.Length; i++) byCase[i] = PairwiseSum.Dot(_weights.AsSpan(i * h.Length, h.Length), h);
        return PairwiseSum.Dot(a, byCase) / Weighting.Mass;
    }

    /// <summary>Q, the inner product on the output index; the identity for independent noise.</summary>
    public IInnerProduct Output { get; }
}

/// <summary>
/// The seam between the noise and the Gaussian core (F-X1, F-X2): the observation model supplies the core its
/// <see cref="WorkingResponse"/> and, for a noise with latent variables, the blocks that update them in the sweep's slots
/// (<see cref="BassSweepSlot"/>). Gaussian BASS is the one member (<see cref="GaussianObservation"/>); generalized-hyperbolic
/// noise joins by supplying an offset and weights at its latent state, the blocks of its latent scales and its law
/// (<see cref="VarianceMeanMixture"/>).
/// </summary>
public interface IBassObservationModel
{
    /// <summary>The working response and weights the core starts on for <paramref name="response"/>; a block moves a chain from there
    /// (<see cref="BassStructureMoves.Follow"/>).</summary>
    WorkingResponse Working(BassResponse response);

    /// <summary>The blocks that fill <paramref name="slot"/> of the sweep, in the order they run; none leaves the slot empty,
    /// which draws nothing and consumes no random number. Asked of <see cref="BassSweepSlot.Coefficients"/>,
    /// <see cref="BassSweepSlot.Latent"/> and <see cref="BassSweepSlot.Scale"/> only: the structure's slot and τ's are BASS's.</summary>
    IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot);

    /// <summary>The noise's law, each observation's latent integrated out, which the pointwise density, prediction and the training fit
    /// read (<see cref="IBassNoiseLaw"/>); null, the default, is the Gaussian N(0, σ²/w) those readers hold themselves.</summary>
    IBassNoiseLaw? Law => null;

    /// <summary>Whether the model rebuilds its working response from its latents (<see cref="Rebuild"/>), so its retained draws keep the
    /// latents alone; false, the default, keeps the working response on the draw.</summary>
    bool Rebuilds => false;

    /// <summary>The working response at <paramref name="latents"/>, rebuilt from <paramref name="model"/>, the model's own: what a retained
    /// draw that keeps its latents alone reads, equal to the one its chain held. Null, the default, for a model that cannot, whose draws
    /// keep their working response.</summary>
    WorkingResponse? Rebuild(WorkingResponse model, ObservationLatents latents) => null;
}

/// <summary>Independent Gaussian noise: the working response is the response, with zero offset and unit weights, so the
/// core's marginal and draws are BASS's own, and every latent slot is empty.</summary>
public sealed class GaussianObservation : IBassObservationModel
{
    /// <summary>The one Gaussian observation model.</summary>
    public static GaussianObservation Instance { get; } = new();

    private GaussianObservation() { }

    /// <inheritdoc/>
    public WorkingResponse Working(BassResponse response) => new(response);

    /// <inheritdoc/>
    public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => [];
}

/// <summary>
/// The slots of the sweep (F-X2), in D04's general order, the robust and generalized-hyperbolic one: the structure given the
/// latents, with the coefficients integrated out; τ; the coefficients a, drawn every sweep just before the latents that read
/// the residuals; the latent scales; then σ² when it is not integrated, and the skew. Gaussian BASS fills the first two and
/// draws a and σ² at retention by composition (D04), so its other slots are empty.
/// </summary>
public enum BassSweepSlot
{
    /// <summary>The reversible-jump structure moves.</summary>
    Structure,

    /// <summary>τ's slice step, when τ is sampled.</summary>
    Tau,

    /// <summary>The coefficients a, drawn each sweep when a latent block reads their residuals.</summary>
    Coefficients,

    /// <summary>The per-observation latent scales.</summary>
    Latent,

    /// <summary>σ² when it is not integrated out, the skew, and σ²'s rate under a hierarchical prior.</summary>
    Scale,
}

/// <summary>The sweep assembled from its slots, in <see cref="BassSweepSlot"/> order: the one place the order is written, for
/// the chain, the tempering replicas and the SMC kernel.</summary>
public static class BassSweepSlots
{
    /// <summary>The slots in the order the sweep runs them.</summary>
    public static IReadOnlyList<BassSweepSlot> Order { get; } =
        [BassSweepSlot.Structure, BassSweepSlot.Tau, BassSweepSlot.Coefficients, BassSweepSlot.Latent, BassSweepSlot.Scale];

    /// <summary>
    /// The blocks of a sweep of <paramref name="model"/> from slot <paramref name="first"/> on, each reading the chain's own target: the
    /// structure's move step, τ's slice step when τ is sampled and a drawn group scale's step (<see cref="BassGroupScaleBlock"/>), then
    /// the observation model's blocks for the slots after τ, and last
    /// in the scale slot σ²'s rate's step under a hierarchical prior, so it reads a σ² the model has just drawn. A chain starts at
    /// the structure slot; a tempering replica or an SMC particle, whose moves the engine runs itself, at τ's.
    /// </summary>
    public static List<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassModel model, BassSweepSlot first = BassSweepSlot.Structure)
    {
        ArgumentNullException.ThrowIfNull(model);
        var blocks = new List<IBlock<ReversibleJumpChain<BassState>>>();
        foreach (BassSweepSlot slot in Order)
        {
            if (slot < first) continue;
            switch (slot)
            {
                case BassSweepSlot.Structure:
                    blocks.Add(new MoveStep<BassState>());
                    break;
                case BassSweepSlot.Tau:
                    if (model.Target.SamplesTau) blocks.Add(new BassTauBlock());
                    if (model.ScaleDrawn) blocks.Add(new BassGroupScaleBlock(model));   // a drawn group scale's step on the marginal (E21)
                    break;
                case BassSweepSlot.Scale:
                    blocks.AddRange(model.Observation.Blocks(slot));
                    if (model.Target.Sigma2Rate is not null) blocks.Add(new BassSigma2RateBlock());
                    break;
                default:
                    blocks.AddRange(model.Observation.Blocks(slot));
                    break;
            }
        }
        return blocks;
    }
}
