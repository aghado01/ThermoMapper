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
    EnsembleCheckpoint? Checkpoint)
{
    /// <summary>σ²'s posterior median and 95% interval, read from the retained draws.</summary>
    public PosteriorInterval NoiseVarianceInterval => Interval(d => d.NoiseVariance);

    /// <summary>τ_u²'s posterior median and 95% interval, read from the retained draws.</summary>
    public PosteriorInterval DeviationVarianceInterval => Interval(d => d.DeviationVariance);

    /// <summary>τ_0²'s posterior median and 95% interval, read from the retained draws.</summary>
    public PosteriorInterval PopulationSmoothingVarianceInterval => Interval(d => d.PopulationSmoothingVariance);

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
/// when the run did not retain replicate draws.
/// </summary>
public sealed record HierarchicalDraw(
    double[] Population, double[][]? Replicates, double NoiseVariance, double DeviationVariance,
    double PopulationSmoothingVariance, double[] PopulationConditionalMean, double[][]? ReplicateConditionalMeans);

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
/// </summary>
public sealed class HierarchicalPSpline
{
    private readonly double[,] _z;
    private readonly double[][] _responses;
    private readonly DifferencePenalty _penalty;
    private readonly BandedDesign _design;
    private readonly int _n;
    private readonly int _nu;
    private readonly int _r;          // penalty null-space dim
    private readonly int _gramBw;
    private readonly double _a0;
    private readonly double _b0;
    private readonly double[,] _gram;           // ZᵀZ — constant
    private readonly double[][] _zty;           // each Zᵀy_r — constant
    private readonly double _responseVariance;  // the pooled variance of every response, the overdispersed start's scale

    public HierarchicalPSpline(double[,] design, double[][] responses, DifferencePenalty penalty,
                               double priorShape = 1e-3, double priorScale = 1e-3)
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
        _r = penalty.Order;
        _gramBw = _design.Bandwidth;
        foreach (double[] y in responses)
            if (y.Length != _n) throw new ArgumentException("Every replicate must match the design row count.", nameof(responses));
        _a0 = priorShape;
        _b0 = priorScale;

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
    }

    /// <param name="chains">Independent chains; R̂ and the ESS need at least two.</param>
    /// <param name="masterSeed">The seed every chain's streams fan out from.</param>
    /// <param name="burn">Burn-in sweeps per chain.</param>
    /// <param name="samples">Retained draws per chain; an early stop may keep fewer.</param>
    /// <param name="retainReplicates">When true (the default), every draw keeps each replicate's curve and conditional
    /// mean, R × ν numbers twice over. When false, draws keep the population curve and the variances only, and the
    /// replicate summaries accumulate as the chains run; use it when R × ν × draws outgrows memory.</param>
    /// <param name="start">How each chain starts (<see cref="HierarchicalStart"/>); overdispersed by default.</param>
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
                                  double rHatTarget = 0.0, int batchSize = 0,
                                  double minEssPerChain = ChainEnsemble.MinEssPerChain, bool checkpoint = false,
                                  EnsembleCheckpoint? resumeFrom = null)
    {
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (samples < 1 || burn < 0) throw new ArgumentOutOfRangeException(nameof(samples));
        if (start is not (HierarchicalStart.Overdispersed or HierarchicalStart.Fixed)) throw new ArgumentOutOfRangeException(nameof(start));
        int reps = _responses.Length;

        (EnsembleRun summary, IReadOnlyList<IChainRun<SweepState>> runs) = new ChainEnsemble().Run<SweepState>(
            new HierarchicalEnsemble(this, retainReplicates, start), chains, masterSeed, burn, samples, rHatTarget, batchSize,
            minEssPerChain, checkpoint, resumeFrom);

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
        RHatDiagnostic[] rHatPopulation = summary.FunctionalRHat[3..];
        return new HierarchicalResult(beta0Sum, repMeans, sigSum / total, tauUSum / total, tau0Sum / total, total,
                                      posteriorDraws, mean0Sum, repConditionalMeans,
                                      summary.Monitored[0], summary.Monitored[1], summary.Monitored[2], rHatPopulation,
                                      RHatDiagnostic.Worst(rHatPopulation), summary.ChainSeeds, summary.Checkpoint);
    }

    /// <summary>
    /// The ensemble-shell client: it starts and restores chains, fingerprints the run for checkpoints, and names the
    /// functionals. The three log variances lead and are monitored (split R̂, bulk and tail ESS, the early stop's ESS
    /// gate); the population curve's coefficients follow with classic R̂, a consensus map an early stop also waits on.
    /// </summary>
    private sealed class HierarchicalEnsemble(HierarchicalPSpline owner, bool retainReplicates, HierarchicalStart start)
        : IEnsembleModel<SweepState>
    {
        public int FunctionalDim => 3 + owner._nu;

        public int MonitoredDim => 3;

        // One child seed per consumer: the start's draws and the sweep's stream. The model has no readout randomness.
        public IChainRun<SweepState> StartChain(int chainSeed, int readoutSeed)
        {
            int[] seeds = SeedTree.Derive(chainSeed, 2);
            int reps = owner._responses.Length;
            SweepState state = start == HierarchicalStart.Fixed
                ? SweepState.Start(reps, owner._nu)
                : SweepState.Overdispersed(reps, owner._nu, owner._responseVariance, new Xoshiro256PlusPlus(seeds[0]));
            var scan = new BlockScan<SweepState>(owner.Blocks(), state, new Xoshiro256PlusPlus(seeds[1]), x => x.Copy());
            return new HierarchicalChainRun(scan, retainReplicates, new List<HierarchicalDraw>(), Zeros(reps, owner._nu),
                                            Zeros(reps, owner._nu));
        }

        /// <summary>The design, the responses, the penalty order, the hyperparameters, the start and whether replicate
        /// draws are retained.</summary>
        public string Fingerprint
        {
            get
            {
                var f = new RunFingerprint().Add("HierarchicalPSpline").Add(owner._n).Add(owner._nu).Add(owner._r)
                                            .Add(owner._responses.Length);
                foreach (double v in owner._z) f.Add(v);
                foreach (double[] y in owner._responses) f.Add(y);
                return f.Add(owner._a0).Add(owner._b0).Add(retainReplicates ? 1 : 0).Add((int)start).ToString();
            }
        }

        public IChainRun<SweepState> RestoreChain(IChainState state) => state is HierarchicalChainState s
            ? new HierarchicalChainRun(BlockScan<SweepState>.Restore(owner.Blocks(), s.Scan, x => x.Copy()), retainReplicates,
                                       new List<HierarchicalDraw>(s.Draws), Clone(s.ReplicateSums), Clone(s.ReplicateMeanSums))
            : throw new ArgumentException("Not a hierarchical P-spline chain's restart state.", nameof(state));
    }

    /// <summary>One chain's restart state: the scan's snapshot (its state, stream and counters), its retained draws, and
    /// the replicate sums a run without retained replicate draws accumulates.</summary>
    private sealed record HierarchicalChainState(
        BlockScanSnapshot<SweepState> Scan, List<HierarchicalDraw> Draws, double[][] ReplicateSums,
        double[][] ReplicateMeanSums) : IChainState;

    /// <summary>
    /// One chain: its block scan, its retained draws and its restart state. Burn-in steps the scan without recording;
    /// each sampling step records the draw, so the chain's stream is the same however its steps are split into rounds or
    /// runs.
    /// </summary>
    private sealed class HierarchicalChainRun(BlockScan<SweepState> scan, bool retainReplicates, List<HierarchicalDraw> draws,
                                              double[][] replicateSums, double[][] replicateMeanSums)
        : IResumableChainRun<SweepState>
    {
        public List<HierarchicalDraw> Draws => draws;           // every post-burn draw, in order
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
            for (int j = 0; j < state.Population.Length; j++) functionals[3 + j] = state.Population[j];

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
            draws.Add(new HierarchicalDraw((double[])state.Population.Clone(), replicates, state.Sigma2, state.TauU2,
                                           state.Tau02, state.PopulationMean, conditionalMeans));
        }

        public IChainState SaveState() => new HierarchicalChainState(scan.Snapshot(), new List<HierarchicalDraw>(draws),
                                                                     Clone(replicateSums), Clone(replicateMeanSums));
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

    // The five blocks in sweep order, each with its own scratch, so every chain builds its own.
    private IBlock<SweepState>[] Blocks() => new IBlock<SweepState>[]
    {
        new ReplicateBlock(this, _gram, _zty), new PopulationBlock(this), new NoiseBlock(this), new DeviationBlock(this),
        new SmoothingBlock(this),
    };

    /// <summary>
    /// The sweep's state: each replicate's curve β_r, the population curve β_0, the three variances, and the conditional
    /// means the last sweep drew around (fresh arrays each sweep, so a retained draw can keep them).
    /// </summary>
    private sealed class SweepState
    {
        public double[][] Replicates = null!;
        public double[] Population = null!;
        public double Sigma2, TauU2, Tau02;
        public double[][] ReplicateMeans = null!;
        public double[] PopulationMean = null!;

        // Overdispersed (HierarchicalStart.Overdispersed): each variance at v·exp(2z) and the population coefficients at
        // N(0, v) draws. The replicate curves are drawn first in every sweep, so their start is never read.
        public static SweepState Overdispersed(int reps, int nu, double v, Xoshiro256PlusPlus rng)
        {
            SweepState state = Start(reps, nu);
            state.Sigma2 = v * Math.Exp(2.0 * Normal.SampleStandard(rng));
            state.TauU2 = v * Math.Exp(2.0 * Normal.SampleStandard(rng));
            state.Tau02 = v * Math.Exp(2.0 * Normal.SampleStandard(rng));
            double sd = Math.Sqrt(v);
            for (int j = 0; j < nu; j++) state.Population[j] = sd * Normal.SampleStandard(rng);
            return state;
        }

        // Every curve at 0 and every variance at 1 (HierarchicalStart.Fixed).
        public static SweepState Start(int reps, int nu)
        {
            var state = new SweepState
            {
                Replicates = new double[reps][], Population = new double[nu], Sigma2 = 1.0, TauU2 = 1.0, Tau02 = 1.0,
                ReplicateMeans = new double[reps][], PopulationMean = new double[nu],
            };
            for (int rr = 0; rr < reps; rr++) { state.Replicates[rr] = new double[nu]; state.ReplicateMeans[rr] = new double[nu]; }
            return state;
        }

        public SweepState Copy()
        {
            var copy = new SweepState
            {
                Replicates = new double[Replicates.Length][], Population = (double[])Population.Clone(),
                Sigma2 = Sigma2, TauU2 = TauU2, Tau02 = Tau02,
                ReplicateMeans = new double[ReplicateMeans.Length][], PopulationMean = (double[])PopulationMean.Clone(),
            };
            for (int rr = 0; rr < Replicates.Length; rr++)
            {
                copy.Replicates[rr] = (double[])Replicates[rr].Clone();
                copy.ReplicateMeans[rr] = (double[])ReplicateMeans[rr].Clone();
            }
            return copy;
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

    // β_0 | · ~ N(Q_0⁻¹ Σβ_r/τ_u², Q_0⁻¹) with Q_0 = (R/τ_u²) I + P/τ_0²: the identity term keeps it proper despite P's
    // null space.
    private sealed class PopulationBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        private readonly double[,] _band = new double[owner._r + 1, owner._nu];
        private readonly BandCholesky _factor = new(owner._nu, owner._r, BandFactorization.Ldlt);
        private readonly double[] _rhs = new double[owner._nu];
        private readonly GaussianConditional _gaussian = new(owner._nu);

        public string Name => "population";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            int nu = owner._nu, reps = state.Replicates.Length;
            Array.Clear(_band);
            for (int j = 0; j < nu; j++) _band[0, j] = reps / state.TauU2;
            owner._penalty.AccumulateInto(_band, nu, 1.0 / state.Tau02);
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

    // τ_0² | · from the population curve's roughness, ν − r terms (the penalty's rank).
    private sealed class SmoothingBlock(HierarchicalPSpline owner) : IBlock<SweepState>
    {
        public string Name => "population smoothing variance";

        public bool Update(SweepState state, Xoshiro256PlusPlus rng)
        {
            state.Tau02 = ConjugateUpdates.NormalVariance(owner._a0, owner._b0, owner._nu - owner._r,
                                                          owner._penalty.Roughness(state.Population)).Sample(rng);
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
