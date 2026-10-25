using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Baps;

/// <summary>Posterior summary of a hierarchical P-spline fit: the population (Bayes) curve, each replicate's
/// shrunken curve, and the variance-component posteriors (posterior means of σ, τ_u and τ_0). Every summary is read
/// from <see cref="PosteriorDraws"/>, every retained draw per chain in draw order, pooled over the chains.</summary>
/// <remarks>
/// <see cref="PopulationCoefficients"/> and <see cref="ReplicateCoefficients"/> average the drawn curves.
/// <see cref="RaoBlackwellizedPopulationCoefficients"/> and <see cref="RaoBlackwellizedReplicateCoefficients"/> average
/// instead the conditional means E[β | ·] the sweep computed. Both estimate the same posterior means. The Rao–Blackwellized
/// ones drop the conditional draw's own noise, which usually lowers the Monte Carlo error; for a scan of more than two
/// blocks that is typical rather than guaranteed (Liu, Wong &amp; Kong 1994 prove it for two).
/// <para><see cref="NoiseVarianceConvergence"/>, <see cref="DeviationVarianceConvergence"/> and
/// <see cref="PopulationSmoothingVarianceConvergence"/> are the ensemble shell's rank-normalized diagnostics (split R̂,
/// bulk and tail ESS) on log σ², log τ_u² and log τ_0², the functionals an early stop waits on.
/// <see cref="RHatPopulation"/> is the classic R̂ of each population coefficient across chains, a consensus map, and
/// <see cref="RHatPopulationMax"/> its least favourable reading. <see cref="ChainSeeds"/> are the resolved per-chain
/// seeds, and <see cref="Checkpoint"/> the run's restart state when it was requested.</para>
/// <para>With anisotropic smoothing each draw's τ_0² is the geometric mean of its components' τ_0j², so the τ_0
/// summaries read the overall smoothing. Each τ_0j²'s posterior is <see cref="ComponentSmoothingVarianceInterval"/>, and
/// <see cref="ComponentSmoothingVarianceConvergence"/> holds the shell's diagnostics on each log τ_0j², which also gate
/// an early stop.</para>
/// </remarks>
public sealed record HierarchicalResult(
    double[] PopulationCoefficients,
    double[][] ReplicateCoefficients,
    double NoiseSd,
    double DeviationSd,
    double PopulationSmoothingSd,
    int Draws,
    IReadOnlyList<IReadOnlyList<HierarchicalDraw>> PosteriorDraws,
    double[] RaoBlackwellizedPopulationCoefficients,
    double[][] RaoBlackwellizedReplicateCoefficients,
    ConvergenceDiagnostic NoiseVarianceConvergence,
    ConvergenceDiagnostic DeviationVarianceConvergence,
    ConvergenceDiagnostic PopulationSmoothingVarianceConvergence,
    RHatDiagnostic[] RHatPopulation,
    RHatDiagnostic RHatPopulationMax,
    int[] ChainSeeds,
    EnsembleCheckpoint? Checkpoint,
    IReadOnlyList<ConvergenceDiagnostic>? ComponentSmoothingVarianceConvergence = null)
{
    /// <summary>σ²'s posterior median and 95% interval, read from the retained draws.</summary>
    public PosteriorInterval NoiseVarianceInterval => Interval(d => d.NoiseVariance);

    /// <summary>τ_u²'s posterior median and 95% interval, read from the retained draws.</summary>
    public PosteriorInterval DeviationVarianceInterval => Interval(d => d.DeviationVariance);

    /// <summary>τ_0²'s posterior median and 95% interval, read from the retained draws.</summary>
    public PosteriorInterval PopulationSmoothingVarianceInterval => Interval(d => d.PopulationSmoothingVariance);

    /// <summary>With anisotropic smoothing, component j's population smoothing variance τ_0j²: its posterior median and
    /// 95% interval.</summary>
    public PosteriorInterval ComponentSmoothingVarianceInterval(int component)
    {
        if (PosteriorDraws.Count == 0 || PosteriorDraws[0].Count == 0
            || PosteriorDraws[0][0].ComponentSmoothingVariance is not double[] first)
            throw new InvalidOperationException("Only a model with anisotropic smoothing has component smoothing variances.");
        if ((uint)component >= (uint)first.Length) throw new ArgumentOutOfRangeException(nameof(component));
        return Interval(d => d.ComponentSmoothingVariance![component]);
    }

    /// <summary>The pointwise band on the population curve at the rows of <paramref name="design"/>: the mean is the
    /// design times <see cref="PopulationCoefficients"/>; the sd and 2.5% and 97.5% quantiles are those of the drawn
    /// population curves.</summary>
    public PointwiseBand PopulationBand(double[,] design)
    {
        var draws = new List<double[]>(Draws);
        foreach (IReadOnlyList<HierarchicalDraw> chain in PosteriorDraws)
            foreach (HierarchicalDraw d in chain) draws.Add(d.Population);
        return PointwiseBand.Of(design, PopulationCoefficients, draws);
    }

    /// <summary>The pointwise band on one replicate's curve, as <see cref="PopulationBand"/>. It needs the replicate
    /// draws, so a run without them refuses.</summary>
    public PointwiseBand ReplicateBand(int replicate, double[,] design)
    {
        if (replicate < 0 || replicate >= ReplicateCoefficients.Length) throw new ArgumentOutOfRangeException(nameof(replicate));
        var draws = new List<double[]>(Draws);
        foreach (IReadOnlyList<HierarchicalDraw> chain in PosteriorDraws)
            foreach (HierarchicalDraw d in chain)
                draws.Add(d.Replicates?[replicate]
                          ?? throw new InvalidOperationException("The run did not retain replicate draws (retainReplicates: false)."));
        return PointwiseBand.Of(design, ReplicateCoefficients[replicate], draws);
    }

    private PosteriorInterval Interval(Func<HierarchicalDraw, double> value)
    {
        var sample = new double[Draws];
        int s = 0;
        foreach (IReadOnlyList<HierarchicalDraw> chain in PosteriorDraws)
            foreach (HierarchicalDraw d in chain) sample[s++] = value(d);
        return PosteriorInterval.Of(sample);
    }
}

/// <summary>
/// One retained draw of the hierarchical P-spline: the population curve β_0, each replicate's curve β_r, the noise
/// variance σ², the deviation variance τ_u² and the population smoothing variance τ_0², all from the end of the sweep.
/// It also keeps the conditional means the sweep computed on the way: E[β_0 | ·] and each E[β_r | ·], which
/// Rao–Blackwellized estimates average. <see cref="Replicates"/> and <see cref="ReplicateConditionalMeans"/> are null
/// when the run did not retain replicate draws. With anisotropic smoothing, <see cref="ComponentSmoothingVariance"/> holds
/// each component's τ_0j², and <see cref="PopulationSmoothingVariance"/> is their geometric mean; otherwise it is null.
/// </summary>
public sealed record HierarchicalDraw(
    double[] Population, double[][]? Replicates, double NoiseVariance, double DeviationVariance,
    double PopulationSmoothingVariance, double[] PopulationConditionalMean, double[][]? ReplicateConditionalMeans,
    double[]? ComponentSmoothingVariance = null);

/// <summary>How each chain of a <see cref="HierarchicalPSpline"/> run starts.</summary>
public enum HierarchicalStart
{
    /// <summary>
    /// Overdispersed, the default: each variance at the responses' pooled variance v times exp(2z), z ~ N(0, 1), so the
    /// chains' variances start orders of magnitude apart, and the population curve's coefficients at independent
    /// N(0, v) draws. R̂ compares chains that started dispersed relative to the posterior (Gelman &amp; Rubin 1992).
    /// </summary>
    Overdispersed,

    /// <summary>Every curve at 0 and every variance at 1, the sampler's original start. Every chain starts alike, which a
    /// check of the first sweep against its known conditionals needs.</summary>
    Fixed,
}

/// <summary>
/// The sweep a <see cref="HierarchicalPSpline"/> run uses. All target the same posterior under a given prior; they differ
/// in how fast they mix. The first four serve the inverse-gamma prior on τ_u², the last the half-t prior on τ_u.
/// </summary>
public enum HierarchicalScheme
{
    /// <summary>The original centred Gibbs sweep: β_r, β_0, σ², τ_u², τ_0². τ_u² mixes slowly when it is small, since it
    /// is drawn from deviations β_r − β_0 that it holds small (the funnel of a centred parameterization).</summary>
    Centred,

    /// <summary>β_0 drawn with every β_r integrated out, then β_r given β_0: (β_0, β_r) as one joint block, which removes
    /// the coupling between the population curve and the replicates. The marginal precision is formed through the
    /// eigendecomposition of ZᵀZ, which stays stable as τ_u² → 0.</summary>
    CollapsedPopulation,

    /// <summary>The centred sweep plus a second draw of τ_u in the non-centred parameterization: after the centred τ_u²
    /// draw, the deviations are standardized, ũ_r = (β_r − β_0)/τ_u, and τ_u is drawn given ũ, where the likelihood is
    /// Gaussian in τ_u (an ancillarity–sufficiency interweaving of the scale, Yu &amp; Meng 2011). The draw is a Metropolis
    /// step with that Gaussian as its proposal, so it accepts on the prior ratio alone.</summary>
    Interweaved,

    /// <summary>The collapsed population block and the scale interweaving in one sweep: β_0 with the β_r integrated out,
    /// β_r given β_0, the variances, then the non-centred τ_u draw. The two cure different couplings, the collapse the
    /// population curve's and the interweaving the deviation scale's, so together they serve both.</summary>
    CollapsedInterweaved,

    /// <summary>The half-t prior's conjugate sweep by parameter expansion (Gelman 2006): β_r = β_0 + ξη_r with
    /// η_r ~ N(0, σ_η² I), ξ ~ N(0, A²) and σ_η² ~ IG(ν/2, ν/2), so τ_u = |ξ|σ_η has the half-t prior and every block is
    /// conjugate. The redundant scale ξ is itself a mixing device (Liu &amp; Wu 1999).</summary>
    ParameterExpanded,
}

/// <summary>
/// A half-t prior on the replicate deviation sd τ_u, with <see cref="Degrees"/> ν and <see cref="Scale"/> A; null scale
/// takes the responses' pooled sd. ν = 1 is the half-Cauchy Gelman (2006) recommends: weakly informative, without the
/// inverse-gamma's sensitivity to its own hyperparameters when τ_u is small or the replicates few. It samples with
/// <see cref="HierarchicalScheme.ParameterExpanded"/>.
/// </summary>
public sealed record HalfTDeviationPrior(double Degrees = 1.0, double? Scale = null);

/// <summary>
/// Hierarchical P-spline for multiple related curves (Behseta–Kass–Wallstrom 2005, the penalized form): replicates
/// <c>y_r = Z β_r + ε</c> share a smooth population mean <c>β_0</c> (penalized, <c>τ_0² P⁻</c>) with replicate
/// deviations <c>β_r ~ N(β_0, τ_u² I)</c>, so each curve borrows strength from the others — shrunk toward the
/// population where its own data is weak. A fully conjugate banded Gibbs sampler, run as a <see cref="BlockScan{TState}"/>
/// of five blocks: <c>β_r | ·</c> for every replicate and <c>β_0 | ·</c> are banded Gaussian draws
/// (<see cref="GaussianConditional"/> on a <see cref="BandCholesky"/> factor), and the three variance components
/// <c>(σ², τ_u², τ_0²)</c> are inverse-gamma (<see cref="ConjugateUpdates.NormalVariance"/>). Each sweep draws β_r,
/// then β_0, then the variances, every draw conditioning on the current values of the others.
/// <c>Q_r = ZᵀZ/σ² + I/τ_u²</c> is shared across replicates so it factors once per sweep; the identity deviation
/// term keeps <c>β_0</c>'s posterior proper. (Unstructured deviations are the clean first cut; smooth deviations
/// would need a sum-to-zero constraint.) As a client of the ensemble shell (<see cref="ChainEnsemble"/>) its chains run
/// in parallel from one seed tree and overdispersed starts, with rank-normalized R̂ and ESS on the log variances, a
/// consensus map of classic R̂ on the population curve, an optional early stop, and checkpoints.
/// <para>
/// With anisotropic smoothing (an <see cref="IAnisotropicPenalty"/>, such as a surface's x- and y-penalties) the
/// population curve's prior precision is <c>Σⱼ Pⱼ/τ_0j²</c>, one smoothing variance per component, and their logs are
/// updated as <see cref="BapsScheme.AnisotropicGibbs"/> updates its smoothing variances: a Newton-proposal Metropolis
/// step on their full conditional, then a random-walk step on each (Bach &amp; Klein 2025). Every scheme accepts it.
/// </para>
/// </summary>
public sealed class HierarchicalPSpline
{
    private readonly double[,] _z;
    private readonly double[][] _responses;
    private readonly IBandPenalty _penalty;
    private readonly BandedDesign _design;
    private readonly int _n;
    private readonly int _nu;
    private readonly int _r;          // the penalty's null-space dimension
    private readonly int _penaltyBw;  // the penalty's half-bandwidth
    private readonly int _gramBw;
    private readonly double _a0;
    private readonly double _b0;
    private readonly double[,] _gram;           // ZᵀZ — constant
    private readonly double[][] _zty;           // each Zᵀy_r — constant
    private readonly double _responseVariance;  // the pooled variance of every response, the overdispersed start's scale
    private readonly double[,] _gramDense;      // ZᵀZ, dense
    private readonly double[] _gramValues;      // its eigenvalues λ_k
    private readonly double[][] _gramVectors;   // its eigenvectors v_k, as rows
    private readonly double[] _sumZtyInEigenbasis;   // v_kᵀ Σ_r Zᵀy_r
    private readonly double[,] _penaltyDense;   // P, dense
    private readonly IAnisotropicPenalty? _anisotropic;   // set when each penalty component has its own τ_0j²
    private readonly double[][,]? _componentDense;       // each component Pⱼ, dense, for the collapsed population block
    private readonly ISmoothingVariancePrior _componentPrior;   // the prior on each τ_0j² with anisotropic smoothing
    private readonly SmoothingVarianceUpdate _componentUpdate;   // and the update of their logs
    private readonly HalfTDeviationPrior? _halfT;   // null: the inverse-gamma prior on τ_u²
    private readonly double _halfTScale;        // the half-t prior's resolved scale A

    /// <param name="design">The banded design Z shared by every replicate (n × ν).</param>
    /// <param name="responses">The replicates' responses, each of length n.</param>
    /// <param name="penalty">The roughness penalty on the population curve: a difference penalty, or any banded
    /// penalty such as the tensor penalty of a surface.</param>
    /// <param name="priorShape">The inverse-gamma prior's shape on each variance (σ², τ_0², and τ_u² unless
    /// <paramref name="deviationPrior"/> is set).</param>
    /// <param name="priorScale">The inverse-gamma prior's scale on each variance, as <paramref name="priorShape"/>.</param>
    /// <param name="deviationPrior">A half-t prior on τ_u in place of the inverse-gamma prior on τ_u²; null (the
    /// default) keeps the inverse-gamma prior.</param>
    /// <param name="anisotropicSmoothing">When true, each component of the penalty, which must be an
    /// <see cref="IAnisotropicPenalty"/>, gets its own population smoothing variance τ_0j², each under the inverse-gamma
    /// prior; false (the default) smooths every component with one τ_0².</param>
    /// <param name="componentPrior">With anisotropic smoothing, the prior on each τ_0j²; null takes the inverse gamma
    /// IG(<paramref name="priorShape"/>, <paramref name="priorScale"/>).</param>
    /// <param name="componentUpdate">With anisotropic smoothing, how the ln τ_0j² are updated
    /// (<see cref="SmoothingVarianceUpdate"/>).</param>
    public HierarchicalPSpline(double[,] design, double[][] responses, IBandPenalty penalty,
                               double priorShape = 1e-3, double priorScale = 1e-3, HalfTDeviationPrior? deviationPrior = null,
                               bool anisotropicSmoothing = false, ISmoothingVariancePrior? componentPrior = null,
                               SmoothingVarianceUpdate componentUpdate = SmoothingVarianceUpdate.Newton)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(responses);
        ArgumentNullException.ThrowIfNull(penalty);
        if (responses.Length == 0) throw new ArgumentException("Need at least one replicate.", nameof(responses));
        _z = design;
        _responses = responses;
        _penalty = penalty;
        _design = new BandedDesign(design);
        _n = _design.Rows;
        _nu = _design.Dimension;
        _r = penalty.Nullity;
        _penaltyBw = penalty.Bandwidth;
        _gramBw = _design.Bandwidth;
        foreach (double[] y in responses)
            if (y.Length != _n) throw new ArgumentException("Every replicate must match the design row count.", nameof(responses));
        _a0 = priorShape;
        _b0 = priorScale;
        _componentPrior = componentPrior ?? new InverseGammaVariancePrior(priorShape, priorScale);
        if (componentUpdate is not (SmoothingVarianceUpdate.Newton or SmoothingVarianceUpdate.Slice))
            throw new ArgumentOutOfRangeException(nameof(componentUpdate));
        _componentUpdate = componentUpdate;

        _gram = new double[_gramBw + 1, _nu];
        _design.Accumulate(null, _gram, null, null);
        _zty = new double[responses.Length][];
        double sum = 0.0, sumSq = 0.0;
        for (int rr = 0; rr < responses.Length; rr++)
        {
            _zty[rr] = new double[_nu];
            _design.AccumulateRhs(null, responses[rr], _zty[rr]);
            foreach (double v in responses[rr]) { sum += v; sumSq += v * v; }
        }
        double count = (double)responses.Length * _n, mean = sum / count, variance = sumSq / count - mean * mean;
        _responseVariance = variance > 0.0 && double.IsFinite(variance) ? variance : 1.0;

        _gramDense = new double[_nu, _nu];
        for (int d = 0; d <= _gramBw; d++)
            for (int j = 0; j + d < _nu; j++) { _gramDense[j + d, j] = _gram[d, j]; _gramDense[j, j + d] = _gram[d, j]; }
        EigenResult eigen = Eigen.DecomposeSymmetric(_gramDense);
        _gramValues = eigen.Eigenvalues;
        _gramVectors = eigen.Eigenvectors;
        _sumZtyInEigenbasis = new double[_nu];
        for (int k = 0; k < _nu; k++)
            foreach (double[] zty in _zty)
                for (int j = 0; j < _nu; j++) _sumZtyInEigenbasis[k] += _gramVectors[k][j] * zty[j];
        var penaltyBand = new double[_penaltyBw + 1, _nu];
        penalty.AccumulateInto(penaltyBand, _nu, 1.0);
        _penaltyDense = Symmetric(penaltyBand, _penaltyBw, _nu);
        if (anisotropicSmoothing)
        {
            _anisotropic = penalty as IAnisotropicPenalty ?? throw new ArgumentException(
                $"Anisotropic smoothing gives each penalty component its own variance; {penalty.GetType().Name} is not an IAnisotropicPenalty.",
                nameof(penalty));
            _componentDense = new double[_anisotropic.Components][,];
            for (int c = 0; c < _componentDense.Length; c++)
            {
                var lambdas = new double[_anisotropic.Components];
                lambdas[c] = 1.0;
                var band = new double[_penaltyBw + 1, _nu];
                _anisotropic.AccumulateComponents(band, _nu, lambdas);
                _componentDense[c] = Symmetric(band, _penaltyBw, _nu);
            }
        }

        if (deviationPrior is not null)
        {
            if (!(deviationPrior.Degrees > 0.0 && double.IsFinite(deviationPrior.Degrees)))
                throw new ArgumentOutOfRangeException(nameof(deviationPrior), "The half-t prior's degrees of freedom must be positive.");
            if (deviationPrior.Scale is double a && !(a > 0.0 && double.IsFinite(a)))
                throw new ArgumentOutOfRangeException(nameof(deviationPrior), "The half-t prior's scale must be positive.");
        }
        _halfT = deviationPrior;
        _halfTScale = deviationPrior?.Scale ?? Math.Sqrt(_responseVariance);
    }

    /// <param name="chains">Independent chains; R̂ and the ESS need at least two.</param>
    /// <param name="masterSeed">The seed every chain's streams fan out from.</param>
    /// <param name="burn">Burn-in sweeps per chain.</param>
    /// <param name="samples">Retained draws per chain; an early stop may keep fewer.</param>
    /// <param name="retainReplicates">When true (the default), every draw keeps each replicate's curve and conditional
    /// mean, R × ν numbers twice over. When false, draws keep the population curve and the variances only, and the
    /// replicate summaries accumulate as the chains run; use it when R × ν × draws outgrows memory.</param>
    /// <param name="start">How each chain starts (<see cref="HierarchicalStart"/>); overdispersed by default.</param>
    /// <param name="scheme">The sweep (<see cref="HierarchicalScheme"/>); null takes the prior's own: centred for the
    /// inverse-gamma prior, parameter-expanded for the half-t prior.</param>
    /// <param name="rHatTarget">When &gt; 0, stop after the first round in which the log variances and the population
    /// coefficients meet this R̂ target and the log variances' bulk and tail ESS reach <paramref name="minEssPerChain"/>
    /// per chain (<see cref="ChainEnsemble"/>'s rule). It needs at least two chains.</param>
    /// <param name="batchSize">Draws per round between stop checks; ≤ 0 is a single round.</param>
    /// <param name="minEssPerChain">The early stop's ESS gate, per chain.</param>
    /// <param name="checkpoint">When true, the result carries the run's restart state.</param>
    /// <param name="resumeFrom">Continue that checkpoint's chains for <paramref name="samples"/> more draws, exactly as the
    /// uninterrupted run would have; the result covers every draw since burn-in. It needs the same data, settings,
    /// master seed and chain count.</param>
    public HierarchicalResult Run(int chains = 4, int masterSeed = 1, int burn = 1000, int samples = 2000,
                                  bool retainReplicates = true, HierarchicalStart start = HierarchicalStart.Overdispersed,
                                  HierarchicalScheme? scheme = null, double rHatTarget = 0.0, int batchSize = 0,
                                  double minEssPerChain = ChainEnsemble.MinEssPerChain, bool checkpoint = false,
                                  EnsembleCheckpoint? resumeFrom = null)
    {
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (samples < 1 || burn < 0) throw new ArgumentOutOfRangeException(nameof(samples));
        if (start is not (HierarchicalStart.Overdispersed or HierarchicalStart.Fixed)) throw new ArgumentOutOfRangeException(nameof(start));
        HierarchicalScheme sweep = scheme ?? (_halfT is null ? HierarchicalScheme.Centred : HierarchicalScheme.ParameterExpanded);
        if (sweep is not (HierarchicalScheme.Centred or HierarchicalScheme.CollapsedPopulation or HierarchicalScheme.Interweaved
                          or HierarchicalScheme.CollapsedInterweaved or HierarchicalScheme.ParameterExpanded))
            throw new ArgumentOutOfRangeException(nameof(scheme));
        if ((sweep == HierarchicalScheme.ParameterExpanded) != (_halfT is not null))
            throw new ArgumentException(_halfT is null
                ? "The parameter-expanded sweep samples the half-t prior; construct the model with a HalfTDeviationPrior."
                : "The half-t prior samples with the parameter-expanded sweep.", nameof(scheme));
        int reps = _responses.Length;

        (EnsembleRun summary, IReadOnlyList<IChainRun<SweepState>> runs) = new ChainEnsemble().Run<SweepState>(
            new HierarchicalEnsemble(this, retainReplicates, start, sweep), chains, masterSeed, burn, samples, rHatTarget,
            batchSize, minEssPerChain, checkpoint, resumeFrom);

        // Every summary is read from the retained draws, chain by chain in draw order; without retained replicate draws,
        // the replicate sums were accumulated as each chain ran.
        int total = chains * summary.SamplesUsed;
        var posteriorDraws = new IReadOnlyList<HierarchicalDraw>[chains];
        var beta0Sum = new double[_nu];
        var mean0Sum = new double[_nu];
        var betaSum = new double[reps][];
        var meanSum = new double[reps][];
        for (int rr = 0; rr < reps; rr++) { betaSum[rr] = new double[_nu]; meanSum[rr] = new double[_nu]; }
        double sigSum = 0.0, tauUSum = 0.0, tau0Sum = 0.0;
        for (int c = 0; c < chains; c++)
        {
            var run = (HierarchicalChainRun)runs[c];
            posteriorDraws[c] = run.Draws.ToArray();
            for (int rr = 0; rr < reps; rr++)
                for (int j = 0; j < _nu; j++)
                {
                    betaSum[rr][j] += run.ReplicateSums[rr][j];
                    meanSum[rr][j] += run.ReplicateMeanSums[rr][j];
                }
            foreach (HierarchicalDraw d in run.Draws)
            {
                for (int j = 0; j < _nu; j++) beta0Sum[j] += d.Population[j];
                for (int j = 0; j < _nu; j++) mean0Sum[j] += d.PopulationConditionalMean[j];
                if (d.Replicates is double[][] replicates)
                    for (int rr = 0; rr < reps; rr++)
                        for (int j = 0; j < _nu; j++) betaSum[rr][j] += replicates[rr][j];
                if (d.ReplicateConditionalMeans is double[][] means)
                    for (int rr = 0; rr < reps; rr++)
                        for (int j = 0; j < _nu; j++) meanSum[rr][j] += means[rr][j];
                sigSum += Math.Sqrt(d.NoiseVariance);
                tauUSum += Math.Sqrt(d.DeviationVariance);
                tau0Sum += Math.Sqrt(d.PopulationSmoothingVariance);
            }
        }

        for (int j = 0; j < _nu; j++) beta0Sum[j] /= total;
        for (int j = 0; j < _nu; j++) mean0Sum[j] /= total;
        var repMeans = new double[reps][];
        var repConditionalMeans = new double[reps][];
        for (int rr = 0; rr < reps; rr++)
        {
            repMeans[rr] = new double[_nu];
            repConditionalMeans[rr] = new double[_nu];
            for (int j = 0; j < _nu; j++) repMeans[rr][j] = betaSum[rr][j] / total;
            for (int j = 0; j < _nu; j++) repConditionalMeans[rr][j] = meanSum[rr][j] / total;
        }
        int components = Components;
        RHatDiagnostic[] rHatPopulation = summary.FunctionalRHat[(3 + components)..];
        return new HierarchicalResult(beta0Sum, repMeans, sigSum / total, tauUSum / total, tau0Sum / total, total,
                                      posteriorDraws, mean0Sum, repConditionalMeans,
                                      summary.Monitored[0], summary.Monitored[1], summary.Monitored[2], rHatPopulation,
                                      RHatDiagnostic.Worst(rHatPopulation), summary.ChainSeeds, summary.Checkpoint,
                                      _anisotropic is null ? null : summary.Monitored[3..]);
    }

    /// <summary>
    /// log p(y_unit | draw) for every retained draw of a result this model produced and every unit (OQ17). Observations are
    /// numbered replicate by replicate, o = r·n + i. A replicate's curve β_r is a group latent: a unit that holds all n of
    /// replicate r's observations integrates β_r out, and its term for them is N(y_r; Zβ_0, σ²I + τ_u²ZZᵀ), the question
    /// of a new curve; a unit that holds only some of them conditions on β_r, the question of a new point on a seen curve,
    /// which needs the retained replicate draws. β_0, σ² and τ_u² are global and always conditioned on. With
    /// ZᵀZ = V diag(λ) Vᵀ the marginal's log-determinant is n·ln σ² + Σ ln(1 + τ_u²λ_k/σ²) and its quadratic
    /// (‖e‖² − Σ (v_kᵀZᵀe)²/(σ²/τ_u² + λ_k))/σ² for e = y_r − Zβ_0.
    /// </summary>
    /// <param name="result">A run of this model, whose retained draws are scored.</param>
    /// <param name="units">The partition of the R·n observations; null takes each observation as its own unit, and
    /// <see cref="UnitPartition.ConsecutiveBlocks"/>(R, n) takes each curve as one.</param>
    public PointwiseLogDensity PointwiseLogDensity(HierarchicalResult result, UnitPartition? units = null)
    {
        ArgumentNullException.ThrowIfNull(result);
        int reps = _responses.Length, total = reps * _n;
        UnitPartition partition = units ?? UnitPartition.EachObservation(total);
        if (partition.Observations != total)
            throw new ArgumentException("The partition must cover the R·n observations, replicate by replicate.", nameof(units));

        // Per unit, its runs of observations from one replicate, and whether each run is the whole replicate.
        var runs = new List<(int Replicate, int[] Points, bool Whole)>[partition.Units];
        var tally = new TreatmentTally();
        bool anyConditioned = false;
        for (int u = 0; u < partition.Units; u++)
        {
            runs[u] = new List<(int, int[], bool)>();
            int[] members = partition.Members[u];
            for (int k = 0; k < members.Length;)
            {
                int r = members[k] / _n, end = k;
                while (end < members.Length && members[end] / _n == r) end++;
                var points = new int[end - k];
                for (int j = k; j < end; j++) points[j - k] = members[j] % _n;
                bool whole = points.Length == _n;
                runs[u].Add((r, points, whole));
                if (whole) tally.Integrated(); else { tally.Conditioned(); anyConditioned = true; }
                k = end;
            }
        }

        int draws = 0;
        foreach (IReadOnlyList<HierarchicalDraw> chain in result.PosteriorDraws) draws += chain.Count;
        var values = new double[draws, partition.Units];
        var marginal = new double[reps];
        var marginalDone = new bool[reps];
        var zte = new double[_nu];
        int row = 0;
        foreach (IReadOnlyList<HierarchicalDraw> chain in result.PosteriorDraws)
            foreach (HierarchicalDraw d in chain)
            {
                if (anyConditioned && d.Replicates is null)
                    throw new InvalidOperationException("Units that split a replicate condition on its curve, which needs retained replicate draws (retainReplicates: true).");
                double[] fit0 = DesignTimes(d.Population);
                Array.Clear(marginalDone);
                double s2 = d.NoiseVariance, ratio = s2 / d.DeviationVariance;
                for (int u = 0; u < partition.Units; u++)
                {
                    double sum = 0.0;
                    foreach ((int r, int[] points, bool whole) in runs[u])
                    {
                        double[] y = _responses[r];
                        if (whole)
                        {
                            if (!marginalDone[r])
                            {
                                double ee = 0.0;
                                Array.Clear(zte);
                                for (int i = 0; i < _n; i++)
                                {
                                    double e = y[i] - fit0[i];
                                    ee += e * e;
                                    for (int j = 0; j < _nu; j++) zte[j] += _z[i, j] * e;
                                }
                                double logDet = _n * Math.Log(s2), reduction = 0.0;
                                for (int k = 0; k < _nu; k++)
                                {
                                    double lambda = Math.Max(0.0, _gramValues[k]), proj = 0.0;
                                    for (int j = 0; j < _nu; j++) proj += _gramVectors[k][j] * zte[j];
                                    logDet += Math.Log(1.0 + lambda / ratio);
                                    reduction += proj * proj / (ratio + lambda);
                                }
                                marginal[r] = -0.5 * (_n * Math.Log(2.0 * Math.PI) + logDet + (ee - reduction) / s2);
                                marginalDone[r] = true;
                            }
                            sum += marginal[r];
                        }
                        else
                        {
                            double[] beta = d.Replicates![r];
                            foreach (int i in points)
                            {
                                double f = 0.0;
                                for (int j = 0; j < _nu; j++) f += _z[i, j] * beta[j];
                                double e = y[i] - f;
                                sum += -0.5 * Math.Log(2.0 * Math.PI * s2) - 0.5 * e * e / s2;
                            }
                        }
                    }
                    values[row, u] = sum;
                }
                row++;
            }

        var latents = new[]
        {
            new LatentReport("replicate curves β_r", LatentScope.Group, tally.Result),
            new LatentReport("population curve β_0", LatentScope.Global, LatentTreatment.Conditioned),
            new LatentReport("noise variance σ²", LatentScope.Global, LatentTreatment.Conditioned),
            new LatentReport("deviation variance τ_u²", LatentScope.Global, LatentTreatment.Conditioned),
        };
        return new PointwiseLogDensity(values, partition, latents);
    }

    /// <summary>
    /// The ensemble-shell client: it starts and restores chains, fingerprints the run for checkpoints, and names the
    /// functionals. The three log variances lead, then with anisotropic smoothing each log τ_0j², all monitored (split R̂,
    /// bulk and tail ESS, the early stop's ESS gate); the population curve's coefficients follow with classic R̂, a
    /// consensus map an early stop also waits on.
    /// </summary>
    private sealed class HierarchicalEnsemble(HierarchicalPSpline owner, bool retainReplicates, HierarchicalStart start,
                                              HierarchicalScheme scheme)
        : IEnsembleModel<SweepState>
    {
        public int FunctionalDim => 3 + owner.Components + owner._nu;

        public int MonitoredDim => 3 + owner.Components;

        // One child seed per consumer: the start's draws and the sweep's stream. The model has no readout randomness.
        public IChainRun<SweepState> StartChain(int chainSeed, int readoutSeed)
        {
            int[] seeds = SeedTree.Derive(chainSeed, 2);
            int reps = owner._responses.Length;
            SweepState state = start == HierarchicalStart.Fixed
                ? SweepState.Start(reps, owner._nu, owner.Components)
                : SweepState.Overdispersed(reps, owner._nu, owner.Components, owner._responseVariance, new Xoshiro256PlusPlus(seeds[0]));
            var scan = new BlockScan<SweepState>(owner.Blocks(scheme), state, new Xoshiro256PlusPlus(seeds[1]), x => x.Copy());
            return new HierarchicalChainRun(scan, retainReplicates, new DrawStore<HierarchicalDraw>(), Zeros(reps, owner._nu),
                                            Zeros(reps, owner._nu));
        }

        /// <summary>The design, the responses, the penalty order, the priors, the start, the scheme and whether replicate
        /// draws are retained.</summary>
        public string Fingerprint
        {
            get
            {
                var f = new RunFingerprint().Add("HierarchicalPSpline").Add(owner._n).Add(owner._nu).Add(owner._r)
                                            .Add(owner._penaltyBw).Add(owner._penalty.GetType().FullName!).Add(owner._responses.Length);
                foreach (double v in owner._z) f.Add(v);
                foreach (double[] y in owner._responses) f.Add(y);
                f.Add(owner._a0).Add(owner._b0).Add(owner._halfT is null ? 0 : 1).Add(owner._halfT?.Degrees ?? 0.0)
                 .Add(owner._halfTScale);
                return f.Add(retainReplicates ? 1 : 0).Add((int)start).Add((int)scheme).Add(owner.Components)
                        .Add(owner._componentPrior.ToString()!).Add((int)owner._componentUpdate).ToString();
            }
        }

        public IChainRun<SweepState> RestoreChain(IChainState state) => state is HierarchicalChainState s
            ? new HierarchicalChainRun(BlockScan<SweepState>.Restore(owner.Blocks(scheme), s.Scan, x => x.Copy()), retainReplicates,
                                       s.Draws.Copy(), Clone(s.ReplicateSums), Clone(s.ReplicateMeanSums))
            : throw new ArgumentException("Not a hierarchical P-spline chain's restart state.", nameof(state));
    }

    /// <summary>One chain's restart state: the scan's snapshot (its state, stream and counters), its retained draws, and
    /// the replicate sums a run without retained replicate draws accumulates.</summary>
    private sealed record HierarchicalChainState(
        BlockScanSnapshot<SweepState> Scan, DrawStore<HierarchicalDraw> Draws, double[][] ReplicateSums,
        double[][] ReplicateMeanSums) : IChainState;

    /// <summary>
    /// One chain: its block scan, its retained draws and its restart state. Burn-in steps the scan without recording;
    /// each sampling step records the draw, so the chain's stream is the same however its steps are split into rounds or
    /// runs.
    /// </summary>
    private sealed class HierarchicalChainRun(BlockScan<SweepState> scan, bool retainReplicates, DrawStore<HierarchicalDraw> draws,
                                              double[][] replicateSums, double[][] replicateMeanSums)
        : IResumableChainRun<SweepState>
    {
        public DrawStore<HierarchicalDraw> Draws => draws;           // every post-burn draw, in order
        public double[][] ReplicateSums => replicateSums;       // Σ β_r, accumulated only without retained replicate draws
        public double[][] ReplicateMeanSums => replicateMeanSums;

        public IChain<SweepState> Chain => scan;

        public void Burn(int steps)
        {
            for (int i = 0; i < steps; i++) scan.Step();
        }

        public void Accumulate(in SweepState state, Span<double> functionals)
        {
            functionals[0] = Math.Log(state.Sigma2);
            functionals[1] = Math.Log(state.TauU2);
            functionals[2] = Math.Log(state.Tau02);
            int components = state.LogTau02.Length;
            for (int c = 0; c < components; c++) functionals[3 + c] = state.LogTau02[c];
            for (int j = 0; j < state.Population.Length; j++) functionals[3 + components + j] = state.Population[j];

            int reps = state.Replicates.Length;
            double[][]? replicates = null, conditionalMeans = null;
            if (retainReplicates)
            {
                replicates = new double[reps][];
                for (int rr = 0; rr < reps; rr++) replicates[rr] = (double[])state.Replicates[rr].Clone();
                conditionalMeans = (double[][])state.ReplicateMeans.Clone();   // each sweep's means are fresh arrays
            }
            else
            {
                for (int rr = 0; rr < reps; rr++)
                    for (int j = 0; j < state.Population.Length; j++)
                    {
                        replicateSums[rr][j] += state.Replicates[rr][j];
                        replicateMeanSums[rr][j] += state.ReplicateMeans[rr][j];
                    }
            }
            double[]? componentTau02 = null;
            if (components > 0)
            {
                componentTau02 = new double[components];
                for (int c = 0; c < components; c++) componentTau02[c] = Math.Exp(state.LogTau02[c]);
            }
            draws.Add(new HierarchicalDraw((double[])state.Population.Clone(), replicates, state.Sigma2, state.TauU2,
                                           state.Tau02, state.PopulationMean, conditionalMeans, componentTau02));
        }

        public IChainState SaveState() => new HierarchicalChainState(scan.Snapshot(), draws.Copy(),
                                                                     Clone(replicateSums), Clone(replicateMeanSums));
    }

    // The number of smoothing components with their own variance: 0 for one τ_0².
    private int Components => _anisotropic?.Components ?? 0;

    // P/τ_0², or with anisotropic smoothing Σⱼ Pⱼ/τ_0j², added into band storage; weights is the caller's scratch.
    private void AccumulatePopulationPenalty(double[,] band, SweepState state, double[] weights)
    {
        if (_anisotropic is null)
        {
            _penalty.AccumulateInto(band, _nu, 1.0 / state.Tau02);
            return;
        }
        for (int c = 0; c < weights.Length; c++) weights[c] = Math.Exp(-state.LogTau02[c]);
        _anisotropic.AccumulateComponents(band, _nu, weights);
    }

    // A symmetric matrix from its lower-band storage.
    private static double[,] Symmetric(double[,] band, int bandwidth, int nu)
    {
        var m = new double[nu, nu];
        for (int d = 0; d <= bandwidth; d++)
            for (int j = 0; j + d < nu; j++) { m[j + d, j] = band[d, j]; m[j, j + d] = band[d, j]; }
        return m;
    }

    private static double[][] Zeros(int rows, int columns)
    {
        var m = new double[rows][];
        for (int i = 0; i < rows; i++) m[i] = new double[columns];
        return m;
    }

    private static double[][] Clone(double[][] m)
    {
        var copy = new double[m.Length][];
        for (int i = 0; i < m.Length; i++) copy[i] = (double[])m[i].Clone();
        return copy;
    }

    // The scheme's blocks in sweep order, each with its own scratch, so every chain builds its own.
    private IBlock<SweepState>[] Blocks(HierarchicalScheme scheme) => scheme switch
    {
        HierarchicalScheme.CollapsedPopulation => new IBlock<SweepState>[]
        {
            new CollapsedPopulationBlock(this), new ReplicateBlock(this, _gram, _zty), new NoiseBlock(this),
            new DeviationBlock(this), new SmoothingBlock(this),
        },
        HierarchicalScheme.Interweaved => new IBlock<SweepState>[]
        {
            new ReplicateBlock(this, _gram, _zty), new PopulationBlock(this), new NoiseBlock(this), new DeviationBlock(this),
            new ScaleInterweaveBlock(this), new SmoothingBlock(this),
        },
        HierarchicalScheme.CollapsedInterweaved => new IBlock<SweepState>[]
        {
            new CollapsedPopulationBlock(this), new ReplicateBlock(this, _gram, _zty), new NoiseBlock(this),
            new DeviationBlock(this), new ScaleInterweaveBlock(this), new SmoothingBlock(this),
        },
        HierarchicalScheme.ParameterExpanded => new IBlock<SweepState>[]
        {
            new ExpandedDeviationBlock(this), new ExpandedPopulationBlock(this), new ExpansionScaleBlock(this),
            new ExpansionVarianceBlock(this), new NoiseBlock(this), new SmoothingBlock(this),
        },
        _ => new IBlock<SweepState>[]
        {
            new ReplicateBlock(this, _gram, _zty), new PopulationBlock(this), new NoiseBlock(this), new DeviationBlock(this),
            new SmoothingBlock(this),
        },
    };

    // Z·β, the curve a coefficient vector draws at the design rows.
    private double[] DesignTimes(double[] beta)
    {
        var fit = new double[_n];
        for (int i = 0; i < _n; i++)
        {
            double f = 0.0;
            for (int j = 0; j < _nu; j++) f += _z[i, j] * beta[j];
            fit[i] = f;
        }
        return fit;
    }

    /// <summary>
    /// The sweep's state: each replicate's curve β_r, the population curve β_0, the three variances, and the conditional
    /// means the last sweep drew around (fresh arrays each sweep, so a retained draw can keep them). The parameter-expanded
    /// sweep also carries η_r, ξ and σ_η², with β_r = β_0 + ξη_r and τ_u² = ξ²σ_η² kept in step. With anisotropic
    /// smoothing it carries each ln τ_0j², and τ_0² is their geometric mean.
    /// </summary>
    private sealed class SweepState
    {
        public double[][] Replicates = null!;
        public double[] Population = null!;
        public double Sigma2, TauU2, Tau02;
        public double[] LogTau02 = null!;
        public double[][] ReplicateMeans = null!;
        public double[] PopulationMean = null!;
        public double[][] Eta = null!;
        public double Xi = 1.0, EtaVariance = 1.0;

        // Overdispersed (HierarchicalStart.Overdispersed): each variance at v·exp(2z) and the population coefficients at
        // N(0, v) draws. The replicate curves are drawn first in every sweep, so their start is never read.
        public static SweepState Overdispersed(int reps, int nu, int components, double v, Xoshiro256PlusPlus rng)
        {
            SweepState state = Start(reps, nu, components);
            state.Sigma2 = v * Math.Exp(2.0 * Normal.SampleStandard(rng));
            state.TauU2 = v * Math.Exp(2.0 * Normal.SampleStandard(rng));
            if (components == 0) state.Tau02 = v * Math.Exp(2.0 * Normal.SampleStandard(rng));
            else
            {
                for (int c = 0; c < components; c++) state.LogTau02[c] = Math.Log(v) + 2.0 * Normal.SampleStandard(rng);
                state.Tau02 = GeometricMean(state.LogTau02);
            }
            state.EtaVariance = state.TauU2;   // ξ = 1, so τ_u² = σ_η²
            double sd = Math.Sqrt(v);
            for (int j = 0; j < nu; j++) state.Population[j] = sd * Normal.SampleStandard(rng);
            return state;
        }

        // Every curve at 0 and every variance at 1 (HierarchicalStart.Fixed).
        public static SweepState Start(int reps, int nu, int components)
        {
            var state = new SweepState
            {
                Replicates = new double[reps][], Population = new double[nu], Sigma2 = 1.0, TauU2 = 1.0, Tau02 = 1.0,
                LogTau02 = new double[components],
                ReplicateMeans = new double[reps][], PopulationMean = new double[nu], Eta = new double[reps][],
            };
            for (int rr = 0; rr < reps; rr++)
            {
                state.Replicates[rr] = new double[nu];
                state.ReplicateMeans[rr] = new double[nu];
                state.Eta[rr] = new double[nu];
            }
            return state;
        }

        public SweepState Copy()
        {
            var copy = new SweepState
            {
                Replicates = new double[Replicates.Length][], Population = (double[])Population.Clone(),
                Sigma2 = Sigma2, TauU2 = TauU2, Tau02 = Tau02, LogTau02 = (double[])LogTau02.Clone(),
                ReplicateMeans = new double[ReplicateMeans.Length][], PopulationMean = (double[])PopulationMean.Clone(),
                Eta = new double[Eta.Length][], Xi = Xi, EtaVariance = EtaVariance,
            };
            for (int rr = 0; rr < Replicates.Length; rr++)
            {
                copy.Replicates[rr] = (double[])Replicates[rr].Clone();
                copy.ReplicateMeans[rr] = (double[])ReplicateMeans[rr].Clone();
                copy.Eta[rr] = (double[])Eta[rr].Clone();
            }
            return copy;
        }

        // exp of the mean log: the geometric mean of the component variances.
        public static double GeometricMean(double[] logs)
        {
            double sum = 0.0;
            foreach (double l in logs) sum += l;
            return Math.Exp(sum / logs.Length);
        }
    }

    // β_r | · ~ N(Q_r⁻¹(Zᵀy_r/σ² + β_0/τ_u²), Q_r⁻¹) with Q_r = ZᵀZ/σ² + I/τ_u², shared by all replicates, so it is
    // factored once per sweep.
    private sealed class ReplicateBlock(HierarchicalPSpline owner, double[,] gram, double[][] zty) : IBlock<SweepState>
    {
        private readonly double[,] _band = new double[owner._gramBw + 1, owner._nu];
        private readonly BandCholesky _factor = new(owner._nu, owner._gramBw, BandFactorization.Ldlt);
        private readonly double[] _rhs = new double[owner._nu];
        private readonly GaussianConditional _gaussian = new(owner._nu);

        public string Name => "replicates";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            int nu = owner._nu;
            for (int d = 0; d <= owner._gramBw; d++)
                for (int j = 0; j < nu; j++) _band[d, j] = gram[d, j] / state.Sigma2;
            for (int j = 0; j < nu; j++) _band[0, j] += 1.0 / state.TauU2;
            _factor.DecomposeBanded(_band);
            if (_factor.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                throw new InvalidOperationException($"ZᵀZ/σ² + I/τ_u² is not numerically positive definite ({_factor.Status}).");
            for (int rr = 0; rr < zty.Length; rr++)
            {
                for (int j = 0; j < nu; j++) _rhs[j] = zty[rr][j] / state.Sigma2 + state.Population[j] / state.TauU2;
                state.ReplicateMeans[rr] = _gaussian.Draw(_factor, _rhs, rng, state.Replicates[rr]);
            }
            return true;
        }
    }

    // β_0 | · ~ N(Q_0⁻¹ Σβ_r/τ_u², Q_0⁻¹) with Q_0 = (R/τ_u²) I + P/τ_0², or Σⱼ Pⱼ/τ_0j² with anisotropic smoothing: the
    // identity term keeps it proper despite P's null space.
    private sealed class PopulationBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        private readonly double[] _weights = new double[owner.Components];
        private readonly double[,] _band = new double[owner._penaltyBw + 1, owner._nu];
        private readonly BandCholesky _factor = new(owner._nu, owner._penaltyBw, BandFactorization.Ldlt);
        private readonly double[] _rhs = new double[owner._nu];
        private readonly GaussianConditional _gaussian = new(owner._nu);

        public string Name => "population";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            int nu = owner._nu, reps = state.Replicates.Length;
            Array.Clear(_band);
            for (int j = 0; j < nu; j++) _band[0, j] = reps / state.TauU2;
            owner.AccumulatePopulationPenalty(_band, state, _weights);
            _factor.DecomposeBanded(_band);
            if (_factor.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                throw new InvalidOperationException($"(R/τ_u²)I + P/τ_0² is not numerically positive definite ({_factor.Status}).");
            Array.Clear(_rhs);
            for (int rr = 0; rr < reps; rr++)
                for (int j = 0; j < nu; j++) _rhs[j] += state.Replicates[rr][j];
            for (int j = 0; j < nu; j++) _rhs[j] /= state.TauU2;
            state.PopulationMean = _gaussian.Draw(_factor, _rhs, rng, state.Population);
            return true;
        }
    }

    // σ² | · from every replicate's residual sum of squares, R·n terms.
    private sealed class NoiseBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        public string Name => "noise variance";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            double sumResid = 0.0;
            for (int rr = 0; rr < state.Replicates.Length; rr++) sumResid += owner.ResidualSS(owner._responses[rr], state.Replicates[rr]);
            state.Sigma2 = ConjugateUpdates.NormalVariance(owner._a0, owner._b0, state.Replicates.Length * owner._n, sumResid)
                                           .Sample(rng);
            return true;
        }
    }

    // τ_u² | · from the deviations of the replicates from the population curve just drawn, R·ν terms.
    private sealed class DeviationBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        public string Name => "deviation variance";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            double sumDev = 0.0;
            for (int rr = 0; rr < state.Replicates.Length; rr++)
                for (int j = 0; j < owner._nu; j++) { double d = state.Replicates[rr][j] - state.Population[j]; sumDev += d * d; }
            state.TauU2 = ConjugateUpdates.NormalVariance(owner._a0, owner._b0, state.Replicates.Length * owner._nu, sumDev)
                                          .Sample(rng);
            return true;
        }
    }

    // τ_0² | · from the population curve's roughness, ν − r terms (the penalty's rank). With anisotropic smoothing, the
    // components' ln τ_0j² | β_0 by ComponentSmoothing's Metropolis update instead, and τ_0² is their geometric mean.
    private sealed class SmoothingBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        private readonly ComponentSmoothing? _components =
            owner._anisotropic is null ? null : new ComponentSmoothing(owner._anisotropic, owner._componentPrior, owner._componentUpdate);

        public string Name => "population smoothing variance";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            if (_components is not null)
            {
                bool accepted = _components.Update(state.LogTau02, state.Population, rng);
                state.Tau02 = SweepState.GeometricMean(state.LogTau02);
                return accepted;
            }
            state.Tau02 = ConjugateUpdates.NormalVariance(owner._a0, owner._b0, owner._nu - owner._r,
                                                          owner._penalty.Roughness(state.Population)).Sample(rng);
            return true;
        }
    }

    // β_0 | σ², τ_u², τ_0², y with every β_r integrated out (HierarchicalScheme.CollapsedPopulation). Each y_r is then
    // N(Zβ_0, σ²I + τ_u²ZZᵀ), so β_0's precision is P/τ_0² + R·Zᵀ(σ²I + τ_u²ZZᵀ)⁻¹Z and its linear term
    // Σ_r Zᵀ(σ²I + τ_u²ZZᵀ)⁻¹y_r. With ZᵀZ = V diag(λ) Vᵀ these are R·V diag(λ/(σ² + τ_u²λ)) Vᵀ and
    // V diag(1/(σ² + τ_u²λ)) Vᵀ Σ_r Zᵀy_r, with no cancellation as τ_u² → 0. The precision is dense (ν × ν).
    private sealed class CollapsedPopulationBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        private readonly double[,] _precision = new double[owner._nu, owner._nu];
        private readonly DenseCholesky _factor = new(owner._nu);
        private readonly double[] _rhs = new double[owner._nu];
        private readonly GaussianConditional _gaussian = new(owner._nu);

        public string Name => "population (collapsed)";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            int nu = owner._nu, reps = state.Replicates.Length;
            if (owner._componentDense is null)
                for (int i = 0; i < nu; i++)
                    for (int j = 0; j < nu; j++) _precision[i, j] = owner._penaltyDense[i, j] / state.Tau02;
            else
            {
                Array.Clear(_precision);
                for (int c = 0; c < owner._componentDense.Length; c++)
                {
                    double weight = Math.Exp(-state.LogTau02[c]);
                    double[,] component = owner._componentDense[c];
                    for (int i = 0; i < nu; i++)
                        for (int j = 0; j < nu; j++) _precision[i, j] += component[i, j] * weight;
                }
            }
            Array.Clear(_rhs);
            for (int k = 0; k < nu; k++)
            {
                double lambda = Math.Max(0.0, owner._gramValues[k]), denominator = state.Sigma2 + state.TauU2 * lambda;
                double weight = reps * lambda / denominator, coefficient = owner._sumZtyInEigenbasis[k] / denominator;
                double[] v = owner._gramVectors[k];
                for (int i = 0; i < nu; i++)
                {
                    _rhs[i] += coefficient * v[i];
                    for (int j = 0; j < nu; j++) _precision[i, j] += weight * v[i] * v[j];
                }
            }
            _factor.Decompose(_precision);
            if (_factor.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                throw new InvalidOperationException($"β_0's collapsed precision is not numerically positive definite ({_factor.Status}).");
            state.PopulationMean = _gaussian.Draw(_factor, _rhs, rng, state.Population);
            return true;
        }
    }

    // The non-centred τ_u draw of the interweaving (HierarchicalScheme.Interweaved). With ũ_r = (β_r − β_0)/τ_u fixed,
    // y_r − Zβ_0 = τ_u·Zũ_r + ε is a regression on τ_u, Gaussian with mean B/A and variance σ²/A, A = Σ‖Zũ_r‖² and
    // B = Σ (y_r − Zβ_0)ᵀZũ_r. That Gaussian proposes τ_u′ over the whole line, and the step accepts on the prior ratio:
    // IG(a, b) on τ_u² is ∝ |τ_u|^(−2a−1) exp(−b/τ_u²) in τ_u. The sign of τ_u′ is unidentified with ũ's, so τ_u′ < 0
    // flips the deviations and leaves τ_u² = τ_u′²; the posterior is symmetric under the flip, so the step leaves it
    // invariant.
    private sealed class ScaleInterweaveBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        private readonly double[][] _standardized = Zeros(owner._responses.Length, owner._nu);

        public string Name => "deviation scale (non-centred)";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            int nu = owner._nu, n = owner._n;
            double tau = Math.Sqrt(state.TauU2);
            if (!(tau > 0.0)) return false;
            double[] fit0 = owner.DesignTimes(state.Population);
            double a = 0.0, b = 0.0;
            for (int rr = 0; rr < state.Replicates.Length; rr++)
            {
                double[] u = _standardized[rr];
                for (int j = 0; j < nu; j++) u[j] = (state.Replicates[rr][j] - state.Population[j]) / tau;
                double[] w = owner.DesignTimes(u);
                double[] y = owner._responses[rr];
                for (int i = 0; i < n; i++) { a += w[i] * w[i]; b += (y[i] - fit0[i]) * w[i]; }
            }
            if (!(a > 0.0)) return false;
            double proposal = b / a + Math.Sqrt(state.Sigma2 / a) * Normal.SampleStandard(rng);
            if (proposal == 0.0) return false;
            double logRatio = LogPrior(Math.Abs(proposal)) - LogPrior(tau);
            if (!(Math.Log(1.0 - rng.NextDouble()) < logRatio)) return false;
            for (int rr = 0; rr < state.Replicates.Length; rr++)
                for (int j = 0; j < nu; j++) state.Replicates[rr][j] = state.Population[j] + proposal * _standardized[rr][j];
            state.TauU2 = proposal * proposal;
            return true;
        }

        private double LogPrior(double tau) => -(2.0 * owner._a0 + 1.0) * Math.Log(tau) - owner._b0 / (tau * tau);
    }

    // η_r | ξ, σ_η², β_0, σ², y (HierarchicalScheme.ParameterExpanded): y_r − Zβ_0 = ξZη_r + ε with η_r ~ N(0, σ_η² I), so
    // the precision ξ²ZᵀZ/σ² + I/σ_η² is shared by the replicates and banded, and the linear term is ξZᵀ(y_r − Zβ_0)/σ².
    private sealed class ExpandedDeviationBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        private readonly double[,] _band = new double[owner._gramBw + 1, owner._nu];
        private readonly BandCholesky _factor = new(owner._nu, owner._gramBw, BandFactorization.Ldlt);
        private readonly double[] _rhs = new double[owner._nu];
        private readonly double[] _gramBeta0 = new double[owner._nu];
        private readonly GaussianConditional _gaussian = new(owner._nu);

        public string Name => "deviations (expanded)";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            int nu = owner._nu;
            double xi = state.Xi;
            for (int d = 0; d <= owner._gramBw; d++)
                for (int j = 0; j < nu; j++) _band[d, j] = xi * xi * owner._gram[d, j] / state.Sigma2;
            for (int j = 0; j < nu; j++) _band[0, j] += 1.0 / state.EtaVariance;
            _factor.DecomposeBanded(_band);
            if (_factor.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                throw new InvalidOperationException($"ξ²ZᵀZ/σ² + I/σ_η² is not numerically positive definite ({_factor.Status}).");
            for (int i = 0; i < nu; i++)
            {
                double g = 0.0;
                for (int j = 0; j < nu; j++) g += owner._gramDense[i, j] * state.Population[j];
                _gramBeta0[i] = g;
            }
            for (int rr = 0; rr < owner._zty.Length; rr++)
            {
                for (int j = 0; j < nu; j++) _rhs[j] = xi * (owner._zty[rr][j] - _gramBeta0[j]) / state.Sigma2;
                double[] mean = _gaussian.Draw(_factor, _rhs, rng, state.Eta[rr]);
                var replicateMean = new double[nu];
                for (int j = 0; j < nu; j++)
                {
                    state.Replicates[rr][j] = state.Population[j] + xi * state.Eta[rr][j];
                    replicateMean[j] = state.Population[j] + xi * mean[j];
                }
                state.ReplicateMeans[rr] = replicateMean;
            }
            return true;
        }
    }

    // β_0 | η, ξ, σ², τ_0², y (HierarchicalScheme.ParameterExpanded): y_r − ξZη_r = Zβ_0 + ε, so the precision is
    // R·ZᵀZ/σ² + P/τ_0², banded, and the linear term Σ_r (Zᵀy_r − ξZᵀZη_r)/σ².
    private sealed class ExpandedPopulationBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        private readonly double[] _weights = new double[owner.Components];
        private readonly double[,] _band = new double[Math.Max(owner._gramBw, owner._penaltyBw) + 1, owner._nu];
        private readonly BandCholesky _factor = new(owner._nu, Math.Max(owner._gramBw, owner._penaltyBw), BandFactorization.Ldlt);
        private readonly double[] _rhs = new double[owner._nu];
        private readonly GaussianConditional _gaussian = new(owner._nu);

        public string Name => "population (expanded)";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            int nu = owner._nu, reps = state.Replicates.Length;
            Array.Clear(_band);
            for (int d = 0; d <= owner._gramBw; d++)
                for (int j = 0; j < nu; j++) _band[d, j] = reps * owner._gram[d, j] / state.Sigma2;
            owner.AccumulatePopulationPenalty(_band, state, _weights);
            _factor.DecomposeBanded(_band);
            if (_factor.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                throw new InvalidOperationException($"R·ZᵀZ/σ² + P/τ_0² is not numerically positive definite ({_factor.Status}).");
            for (int i = 0; i < nu; i++)
            {
                double s = 0.0;
                for (int rr = 0; rr < reps; rr++)
                {
                    double g = 0.0;
                    for (int j = 0; j < nu; j++) g += owner._gramDense[i, j] * state.Eta[rr][j];
                    s += owner._zty[rr][i] - state.Xi * g;
                }
                _rhs[i] = s / state.Sigma2;
            }
            state.PopulationMean = _gaussian.Draw(_factor, _rhs, rng, state.Population);
            for (int rr = 0; rr < reps; rr++)
                for (int j = 0; j < nu; j++) state.Replicates[rr][j] = state.Population[j] + state.Xi * state.Eta[rr][j];
            return true;
        }
    }

    // ξ | η, β_0, σ², y (HierarchicalScheme.ParameterExpanded): y_r − Zβ_0 = ξ·Zη_r + ε with ξ ~ N(0, A²), a scalar
    // regression with precision Σ‖Zη_r‖²/σ² + 1/A².
    private sealed class ExpansionScaleBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        public string Name => "expansion scale ξ";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            int nu = owner._nu, n = owner._n;
            double[] fit0 = owner.DesignTimes(state.Population);
            double precision = 1.0 / (owner._halfTScale * owner._halfTScale), linear = 0.0;
            for (int rr = 0; rr < state.Replicates.Length; rr++)
            {
                double[] w = owner.DesignTimes(state.Eta[rr]);
                double[] y = owner._responses[rr];
                for (int i = 0; i < n; i++) { precision += w[i] * w[i] / state.Sigma2; linear += (y[i] - fit0[i]) * w[i] / state.Sigma2; }
            }
            state.Xi = linear / precision + Normal.SampleStandard(rng) / Math.Sqrt(precision);
            for (int rr = 0; rr < state.Replicates.Length; rr++)
                for (int j = 0; j < nu; j++) state.Replicates[rr][j] = state.Population[j] + state.Xi * state.Eta[rr][j];
            state.TauU2 = state.Xi * state.Xi * state.EtaVariance;
            return true;
        }
    }

    // σ_η² | η (HierarchicalScheme.ParameterExpanded): IG(ν/2, ν/2) prior, R·ν terms.
    private sealed class ExpansionVarianceBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        public string Name => "expansion variance σ_η²";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            double ss = 0.0;
            foreach (double[] eta in state.Eta)
                foreach (double e in eta) ss += e * e;
            double degrees = owner._halfT!.Degrees;
            state.EtaVariance = ConjugateUpdates.NormalVariance(0.5 * degrees, 0.5 * degrees, state.Eta.Length * owner._nu, ss)
                                                .Sample(rng);
            state.TauU2 = state.Xi * state.Xi * state.EtaVariance;
            return true;
        }
    }

    private double ResidualSS(double[] y, double[] beta)
    {
        double ss = 0.0;
        for (int i = 0; i < _n; i++)
        {
            double fit = 0.0;
            for (int j = 0; j < _nu; j++) fit += _z[i, j] * beta[j];
            double d = y[i] - fit;
            ss += d * d;
        }
        return ss;
    }
}
