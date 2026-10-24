using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Baps;

/// <summary>The scheme a <see cref="Baps"/> run samples with: the model, and the sampler embedded in the method.</summary>
public enum BapsScheme
{
    /// <summary>
    /// Conjugate Gibbs over the variance components, a block Gibbs scan: <c>β | σ²,τ²</c> from its banded Gaussian full
    /// conditional, <c>σ² | β</c> and <c>τ² | β</c> from inverse-gamma full conditionals (λ = σ²/τ² derived). Full Bayesian.
    /// </summary>
    Gibbs,

    /// <summary>
    /// Metropolis on <c>log λ</c> against the profiled REML marginal evidence (β integrated out) under a uniform prior
    /// on the effective degrees of freedom (<see cref="PenalizedSpline.EdfSlope"/>). Each retained λ then draws σ² and β
    /// exactly from their conditionals, and the coefficient mean is Rao–Blackwellized on β̂(λ). Marginalizes the spline.
    /// The evidence levels off as λ → ∞, so a flat prior on <c>log λ</c> would be improper; the edf prior decays like
    /// 1/λ and makes the target proper. A collapsed Metropolis sampler.
    /// </summary>
    Collapsed,

    /// <summary>
    /// Locally-adaptive Gibbs (the "A" in BAPS): as <see cref="Gibbs"/> but each difference carries its own mean-1
    /// multiplier <c>η_i ~ Gamma(φ, φ)</c>, so the smoothing <c>λ_i = λ_global·η_i</c> varies across the domain
    /// — heavy where the fit is flat, light where it must track a feature. Every conditional stays conjugate
    /// (the prior normalizer factorizes as <c>(Πη_i)|DDᵀ|</c>, giving each <c>η_i</c> a Gamma full conditional). A block
    /// Gibbs scan with a fourth block, the multipliers.
    /// </summary>
    LocallyAdaptiveGibbs
}

/// <summary>Posterior summary of a BAPS run: the Bayes fit coefficients and the smoothing / noise posteriors.</summary>
/// <remarks>
/// <paramref name="LambdaMedian"/>/<paramref name="LambdaLo"/>/<paramref name="LambdaHi"/> are the 50%, 2.5% and 97.5%
/// quantiles of the pooled draws of the global smoothing (in the locally adaptive scheme, λ_global). A median, not a mean,
/// because λ can be heavy-tailed enough that its posterior mean does not exist. <paramref name="LocalSmoothing"/> is
/// the posterior-mean local smoothing field λ_i (length ν−r) — non-null only for <see cref="BapsScheme.LocallyAdaptiveGibbs"/>; its dips
/// mark where the fit was allowed to track features. <paramref name="RHatLogLambda"/> is classic R̂ on log λ across chains
/// with its <see cref="DiagnosticState"/> (unavailable for a single chain). <paramref name="PosteriorDraws"/> keeps
/// every retained draw (<see cref="BapsDraw"/>), per chain in draw order; every summary above is read from it.
/// <paramref name="LogLambdaConvergence"/> and <paramref name="LogSigma2Convergence"/> are the ensemble shell's
/// rank-normalized diagnostics (split R̂, bulk and tail ESS) on log λ and log σ², the functionals an early stop waits
/// on. <paramref name="AcceptanceRate"/> is the sampling phase's: 1 for the Gibbs schemes, the Metropolis rate on log λ
/// for the collapsed scheme. <paramref name="ChainSeeds"/> are the resolved per-chain seeds, and
/// <paramref name="Checkpoint"/> the run's restart state when it was requested.
/// </remarks>
public sealed record BapsResult(
    double[] Coefficients,
    double LambdaMedian,
    double LambdaLo,
    double LambdaHi,
    double SigmaMean,
    RHatDiagnostic RHatLogLambda,
    int Draws,
    IReadOnlyList<IReadOnlyList<BapsDraw>> PosteriorDraws,
    ConvergenceDiagnostic LogLambdaConvergence,
    ConvergenceDiagnostic LogSigma2Convergence,
    double AcceptanceRate,
    int[] ChainSeeds,
    EnsembleCheckpoint? Checkpoint,
    double[]? LocalSmoothing = null)
{
    /// <summary>σ²'s posterior median and 95% interval, read from the retained draws.</summary>
    public PosteriorInterval Sigma2Interval => Interval(d => d.Sigma2);

    /// <summary>τ²'s posterior median and 95% interval, read from the retained draws (τ² = σ²/λ in every scheme).</summary>
    public PosteriorInterval Tau2Interval => Interval(d => d.Tau2);

    /// <summary>
    /// The pointwise band on the fitted curve at the rows of <paramref name="design"/> (any rows of the coefficients'
    /// basis), read from the retained draws. <see cref="PointwiseBand.Mean"/> is the design times
    /// <see cref="Coefficients"/>, the affine readout, Rao–Blackwellized in the collapsed scheme. The sd and the
    /// 2.5% and 97.5% quantiles are those of the draws' fits Zβ, so they carry the coefficients' uncertainty as well as
    /// λ's.
    /// </summary>
    public PointwiseBand FitBand(double[,] design)
    {
        var draws = new List<double[]>(Draws);
        foreach (IReadOnlyList<BapsDraw> chain in PosteriorDraws)
            foreach (BapsDraw d in chain) draws.Add(d.Coefficients);
        return PointwiseBand.Of(design, Coefficients, draws);
    }

    private PosteriorInterval Interval(Func<BapsDraw, double> value)
    {
        var sample = new double[Draws];
        int s = 0;
        foreach (IReadOnlyList<BapsDraw> chain in PosteriorDraws)
            foreach (BapsDraw d in chain) sample[s++] = value(d);
        return PosteriorInterval.Of(sample);
    }
}

/// <summary>
/// One retained BAPS draw: the coefficients, the global smoothing λ, the noise variance σ², the penalty variance τ²
/// (λ = σ²/τ²), and the locally adaptive scheme's local multipliers η.
/// <list type="bullet">
///   <item>In the Gibbs schemes all of them come from the end of the sweep.</item>
///   <item>The collapsed scheme samples λ, then draws σ² and β exactly from their conditionals given λ, and
///   sets τ² = σ²/λ. <see cref="ConditionalMean"/> holds β̂(λ) = E[β | λ, y] beside the draw.</item>
///   <item><see cref="Eta"/> is set in the locally adaptive scheme only, where the local smoothing is λᵢ = λ·ηᵢ.</item>
/// </list>
/// </summary>
public sealed record BapsDraw(double[] Coefficients, double[]? ConditionalMean, double Lambda, double Sigma2, double Tau2,
                              double[]? Eta);

/// <summary>
/// BAPS (Bayesian adaptive P-splines) — the measure-side dual of free-knot BARS. A rich fixed B-spline basis with a
/// difference penalty whose strength λ is inferred (Lang &amp; Brezger 2004), sampled either by a block Gibbs scan over
/// the variance components or by collapsed Metropolis against the REML marginal evidence (<see cref="BapsScheme"/>).
/// Both samplers are built on the shared blocks (<see cref="BlockScan{TState}"/>, <see cref="GaussianConditional"/>,
/// <see cref="ConjugateUpdates"/>, <see cref="RandomWalkMetropolis"/>), with the model's blocks here. Every inner solve
/// is banded (<c>A = ZᵀZ + λP</c> via <see cref="BandCholesky"/>'s root-free LDLᵀ, whose factor also supplies the
/// Gaussian β-draw, the block update of a Gaussian Markov random field, Rue 2001). As a client of the ensemble shell
/// (<see cref="ChainEnsemble"/>) its chains run in
/// parallel from one seed tree, with rank-normalized R̂ and ESS on log λ and log σ², an optional early stop, and
/// checkpoints. Every draw is retained, and the Bayes fit, the λ posterior (median + 95% interval), σ and the bands are
/// read from the draws.
/// </summary>
public sealed class Baps
{
    private readonly double[,] _z;
    private readonly double[] _y;
    private readonly DifferencePenalty _penalty;
    private readonly BandedDesign _design;
    private readonly PenalizedSpline _model;
    private readonly int _n;
    private readonly int _nu;
    private readonly int _bw;
    private readonly int _r;
    private readonly double[,] _gram0;   // ZᵀZ band — constant across draws
    private readonly double[] _zty;      // Zᵀy — constant
    private readonly double _priorShape; // inverse-gamma hyperprior on each variance component
    private readonly double _priorScale;
    private readonly double _mhLogStep;  // Metropolis log-λ proposal sd (the collapsed scheme)
    private readonly double _phi;        // adaptivity: η_i ~ Gamma(φ, φ); φ→∞ recovers global smoothing

    public Baps(double[,] design, double[] y, DifferencePenalty penalty,
                       double priorShape = 1e-3, double priorScale = 1e-3, double mhLogStep = 0.6, double phi = 1.0)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(penalty);
        _z = design;
        _y = y;
        _penalty = penalty;
        _design = new BandedDesign(design);
        _model = new PenalizedSpline(design, y, penalty);
        _n = _design.Rows;
        _nu = _design.Dimension;
        _r = penalty.Order;
        if (y.Length != _n)
            throw new ArgumentException("Response length must match the design row count.", nameof(y));
        _bw = Math.Max(_design.Bandwidth, penalty.Order);
        _gram0 = new double[_bw + 1, _nu];
        _zty = new double[_nu];
        _design.Accumulate(null, _gram0, y, _zty);   // ZᵀZ and Zᵀy, once
        _priorShape = priorShape;
        _priorScale = priorScale;
        _mhLogStep = mhLogStep;
        _phi = phi;
    }

    /// <param name="scheme">The sampling scheme (<see cref="BapsScheme"/>).</param>
    /// <param name="chains">Independent chains; R̂ and the ESS need at least two.</param>
    /// <param name="masterSeed">The seed every chain's streams fan out from.</param>
    /// <param name="burn">Burn-in sweeps per chain.</param>
    /// <param name="samples">Retained draws per chain; an early stop may keep fewer.</param>
    /// <param name="rHatTarget">When &gt; 0, stop after the first round in which log λ and log σ² meet this R̂ target
    /// and their bulk and tail ESS reach <paramref name="minEssPerChain"/> per chain (<see cref="ChainEnsemble"/>'s
    /// rule). It needs at least two chains.</param>
    /// <param name="batchSize">Draws per round between stop checks; ≤ 0 is a single round.</param>
    /// <param name="minEssPerChain">The early stop's ESS gate, per chain.</param>
    /// <param name="checkpoint">When true, the result carries the run's restart state.</param>
    /// <param name="resumeFrom">Continue that checkpoint's chains for <paramref name="samples"/> more draws, exactly as the
    /// uninterrupted run would have; the result covers every draw since burn-in. It needs the same data, settings, scheme,
    /// master seed and chain count.</param>
    public BapsResult Run(BapsScheme scheme, int chains = 4, int masterSeed = 1, int burn = 500, int samples = 1500,
                          double rHatTarget = 0.0, int batchSize = 0, double minEssPerChain = ChainEnsemble.MinEssPerChain,
                          bool checkpoint = false, EnsembleCheckpoint? resumeFrom = null)
    {
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (samples < 2 || burn < 0) throw new ArgumentOutOfRangeException(nameof(samples));
        if (scheme is not (BapsScheme.Gibbs or BapsScheme.Collapsed or BapsScheme.LocallyAdaptiveGibbs))
            throw new ArgumentOutOfRangeException(nameof(scheme));

        // The shell gives each chain a kernel stream and a readout stream. Only the collapsed scheme's conditional draws
        // use the readout stream, so the λ chain does not depend on them.
        (EnsembleRun summary, IReadOnlyList<IChainRun<BapsState>> runs) = new ChainEnsemble().Run<BapsState>(
            new BapsEnsemble(this, scheme), chains, masterSeed, burn, samples, rHatTarget, batchSize, minEssPerChain,
            checkpoint, resumeFrom);

        // Every summary is read from the retained draws: each chain's sums in draw order, then over chains, and the pooled
        // λ sample for its quantiles. R̂ on log λ is the shell's classic reading of the same draws; the rank-normalized
        // diagnostics come from its traces of log λ and log σ².
        double total = chains * (double)summary.SamplesUsed;
        var coef = new double[_nu];
        double sigSum = 0.0;
        var allLam = new double[(int)total];
        int nLoc = _nu - _r;
        double[]? localField = scheme == BapsScheme.LocallyAdaptiveGibbs && nLoc > 0 ? new double[nLoc] : null;
        var posteriorDraws = new IReadOnlyList<BapsDraw>[chains];
        int idx = 0;
        for (int c = 0; c < chains; c++)
        {
            DrawStore<BapsDraw> draws = ((BapsChainRun)runs[c]).Draws;
            posteriorDraws[c] = draws.ToArray();
            var coefSum = new double[_nu];
            double sigmaSum = 0.0;
            double[]? localSum = localField is null ? null : new double[nLoc];
            foreach (BapsDraw d in draws)
            {
                double[] coefficients = d.ConditionalMean ?? d.Coefficients;   // Rao–Blackwellized where the scheme has it
                for (int p = 0; p < _nu; p++) coefSum[p] += coefficients[p];
                allLam[idx++] = d.Lambda;
                sigmaSum += Math.Sqrt(d.Sigma2);
                if (localSum is not null)
                    for (int k = 0; k < nLoc; k++) localSum[k] += d.Lambda * d.Eta![k];
            }
            for (int p = 0; p < _nu; p++) coef[p] += coefSum[p];
            sigSum += sigmaSum;
            if (localSum is not null)
                for (int k = 0; k < nLoc; k++) localField![k] += localSum[k];
        }
        for (int p = 0; p < _nu; p++) coef[p] /= total;
        if (localField is not null)
            for (int k = 0; k < nLoc; k++) localField[k] /= total;

        Array.Sort(allLam);

        return new BapsResult(
            coef, PosteriorInterval.Quantile(allLam, 0.5), PosteriorInterval.Quantile(allLam, 0.025),
            PosteriorInterval.Quantile(allLam, 0.975),
            sigSum / total, summary.FunctionalRHat[0], (int)total, posteriorDraws,
            summary.Monitored[0], summary.Monitored[1], summary.AcceptanceRate, summary.ChainSeeds, summary.Checkpoint,
            localField);
    }

    // log p(u | y) up to a constant on u = ln λ: the profiled evidence plus the uniform-edf prior density |d edf/du|.
    private double MarginalLogTarget(double lambda) => _model.RemlLogEvidence(lambda) + Math.Log(_model.EdfSlope(lambda));

    /// <summary>
    /// The ensemble-shell client for one scheme: it starts and restores chains, fingerprints the run for checkpoints, and
    /// names the functionals the shell monitors. Those are log λ and log σ², which get split R̂ and bulk and tail ESS and
    /// gate an early stop.
    /// </summary>
    private sealed class BapsEnsemble(Baps owner, BapsScheme scheme) : IEnsembleModel<BapsState>
    {
        public int FunctionalDim => 2;

        public int MonitoredDim => 2;

        // The kernel stream drives the scan; the readout stream serves the collapsed scheme's conditional draws.
        public IChainRun<BapsState> StartChain(int chainSeed, int readoutSeed) => new BapsChainRun(owner, scheme,
            new BlockScan<BapsState>(owner.Blocks(scheme), BapsState.Start(owner, scheme), new Xoshiro256PlusPlus(chainSeed),
                                     state => state.Copy()),
            new Xoshiro256PlusPlus(readoutSeed), new DrawStore<BapsDraw>());

        /// <summary>The design, the responses, the penalty order, the hyperparameters and the scheme.</summary>
        public string Fingerprint
        {
            get
            {
                var f = new RunFingerprint().Add("BAPS").Add((int)scheme).Add(owner._n).Add(owner._nu).Add(owner._r);
                foreach (double v in owner._z) f.Add(v);
                return f.Add(owner._y).Add(owner._priorShape).Add(owner._priorScale).Add(owner._mhLogStep)
                        .Add(owner._phi).ToString();
            }
        }

        public IChainRun<BapsState> RestoreChain(IChainState state) => state is BapsChainState s
            ? new BapsChainRun(owner, scheme, BlockScan<BapsState>.Restore(owner.Blocks(scheme), s.Scan, x => x.Copy()),
                               Xoshiro256PlusPlus.Restore(s.Readout.Item1, s.Readout.Item2, s.Readout.Item3, s.Readout.Item4),
                               s.Draws.Copy())
            : throw new ArgumentException("Not a BAPS chain's restart state.", nameof(state));
    }

    /// <summary>One BAPS chain's restart state: the scan's snapshot (the state its next sweep conditions on, the kernel
    /// stream and the blocks' counters), the readout stream's words and the retained draws.</summary>
    private sealed record BapsChainState(
        BlockScanSnapshot<BapsState> Scan, (ulong, ulong, ulong, ulong) Readout, DrawStore<BapsDraw> Draws) : IChainState;

    /// <summary>
    /// One BAPS chain: its scheme's block scan, its retained draws and its restart state. Burn-in steps the scan without
    /// recording; each sampling step records the draw, so the chain's streams are the same however its steps are split
    /// into rounds or runs.
    /// </summary>
    private sealed class BapsChainRun(Baps owner, BapsScheme scheme, BlockScan<BapsState> scan, Xoshiro256PlusPlus readout,
                                      DrawStore<BapsDraw> draws) : IResumableChainRun<BapsState>
    {
        private readonly GaussianConditional _gaussian = new(owner._nu);   // the collapsed scheme's conditional β draw

        public DrawStore<BapsDraw> Draws => draws;           // every post-burn draw, in order

        public IChain<BapsState> Chain => scan;

        public void Burn(int steps)
        {
            for (int i = 0; i < steps; i++) scan.Step();
        }

        public void Accumulate(in BapsState state, Span<double> functionals)
        {
            BapsDraw draw = Record(state);
            functionals[0] = Math.Log(draw.Lambda);
            functionals[1] = Math.Log(draw.Sigma2);
            draws.Add(draw);
        }

        public IChainState SaveState() => new BapsChainState(scan.Snapshot(), readout.SaveState(), draws.Copy());

        // The draw at the end of the sweep. The collapsed scheme adds its exact conjugate draws given λ here, on the readout
        // stream: σ² | λ, y ~ IG((n − r)/2, RSS_pen/2) under the Jeffreys prior the evidence assumes, then
        // β | σ², λ, y ~ N(β̂, σ²A⁻¹), with β̂(λ) kept for the Rao–Blackwellized coefficient mean.
        private BapsDraw Record(BapsState state)
        {
            if (scheme != BapsScheme.Collapsed)
                return new BapsDraw((double[])state.Beta.Clone(), null, state.Sigma2 / state.Tau2, state.Sigma2, state.Tau2,
                                    scheme == BapsScheme.LocallyAdaptiveGibbs ? (double[])state.Eta.Clone() : null);

            double lambda = Math.Exp(state.LogLambda);
            (double[] betaHat, double rssPen, BandCholesky factor) = owner._model.SolveFactored(lambda);
            double sigma2 = ConjugateUpdates.NormalVariance(0.0, 0.0, owner._n - owner._r, rssPen).Sample(readout);
            var beta = new double[owner._nu];
            _gaussian.Perturb(factor, betaHat, readout, beta, Math.Sqrt(sigma2));   // v ~ N(0, A⁻¹), scaled by σ
            return new BapsDraw(beta, betaHat, lambda, sigma2, sigma2 / lambda, null);
        }
    }

    /// <summary>
    /// A BAPS chain's state: the coefficients β and the variances σ² and τ² (the Gibbs schemes), the local multipliers η
    /// and the squared differences they weigh (the locally adaptive scheme), and log λ with its log target (the collapsed
    /// scheme).
    /// </summary>
    private sealed class BapsState
    {
        public double[] Beta = null!, Eta = null!, SquaredDifferences = null!;
        public double Sigma2 = 1.0, Tau2 = 1.0, LogLambda, Evidence;

        // β at 0, σ² = τ² = 1, η at 1 (unadapted, the global smoothing), and log λ = 0.
        public static BapsState Start(Baps owner, BapsScheme scheme)
        {
            int nLoc = owner._nu - owner._r;
            var state = new BapsState { Beta = new double[owner._nu], Eta = new double[nLoc], SquaredDifferences = new double[nLoc] };
            Array.Fill(state.Eta, 1.0);
            if (scheme == BapsScheme.Collapsed) state.Evidence = owner.MarginalLogTarget(1.0);
            return state;
        }

        public BapsState Copy() => new()
        {
            Beta = (double[])Beta.Clone(), Eta = (double[])Eta.Clone(), SquaredDifferences = (double[])SquaredDifferences.Clone(),
            Sigma2 = Sigma2, Tau2 = Tau2, LogLambda = LogLambda, Evidence = Evidence,
        };
    }

    // The scheme's blocks in sweep order: the block Gibbs scan of the Gibbs schemes, or the collapsed Metropolis step on
    // log λ. Each chain builds its own, since they hold scratch.
    private IBlock<BapsState>[] Blocks(BapsScheme scheme) => scheme switch
    {
        BapsScheme.Collapsed => new IBlock<BapsState>[] { new LogLambdaBlock(this) },
        BapsScheme.LocallyAdaptiveGibbs => new IBlock<BapsState>[]
        {
            new CoefficientBlock(this, adaptive: true), new NoiseBlock(this), new PenaltyBlock(this, adaptive: true), new LocalBlock(this),
        },
        _ => new IBlock<BapsState>[] { new CoefficientBlock(this, adaptive: false), new NoiseBlock(this), new PenaltyBlock(this, adaptive: false) },
    };

    // β | σ², τ², y ~ N(β̂, σ²A⁻¹) with A = ZᵀZ + λP and λ = σ²/τ², or, in the locally adaptive scheme,
    // A = ZᵀZ + Dᵀdiag(λ_i)D with λ_i = λ·η_i; β̂ = A⁻¹Zᵀy.
    private sealed class CoefficientBlock(Baps owner, bool adaptive) : IBlock<BapsState>
    {
        private readonly double[,] _band = new double[owner._bw + 1, owner._nu];
        private readonly BandCholesky _factor = new(owner._nu, owner._bw, BandFactorization.Ldlt);
        private readonly double[] _weights = new double[owner._nu - owner._r];   // λ_i = (σ²/τ²)·η_i
        private readonly GaussianConditional _gaussian = new(owner._nu);

        public string Name => "coefficients";

        public bool Update(BapsState state, Xoshiro256PlusPlus rng)
        {
            int nu = owner._nu;
            double lambda = state.Sigma2 / state.Tau2;
            Array.Copy(owner._gram0, _band, owner._gram0.Length);
            if (adaptive)
            {
                for (int i = 0; i < _weights.Length; i++) _weights[i] = lambda * state.Eta[i];
                owner._penalty.AccumulateInto(_band, nu, _weights);
                _factor.DecomposeBanded(_band);
                if (_factor.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                    throw new InvalidOperationException($"ZᵀZ + Dᵀdiag(λᵢ)D is not numerically positive definite at global λ = {lambda} ({_factor.Status}).");
            }
            else
            {
                owner._penalty.AccumulateInto(_band, nu, lambda);
                _factor.DecomposeBanded(_band);
                if (_factor.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                    throw new InvalidOperationException($"ZᵀZ + λP is not numerically positive definite at λ = {lambda} ({_factor.Status}).");
            }
            _gaussian.Draw(_factor, owner._zty, rng, state.Beta, Math.Sqrt(state.Sigma2));
            return true;
        }
    }

    // σ² | β from the residual sum of squares, n terms.
    private sealed class NoiseBlock(Baps owner) : IBlock<BapsState>
    {
        public string Name => "noise variance";

        public bool Update(BapsState state, Xoshiro256PlusPlus rng)
        {
            state.Sigma2 = ConjugateUpdates.NormalVariance(owner._priorShape, owner._priorScale, owner._n, owner.ResidualSS(state.Beta))
                                           .Sample(rng);
            return true;
        }
    }

    // τ² | β from the roughness, ν − r terms; in the locally adaptive scheme from the η-weighted roughness Σ η_i (Δβ)_i².
    private sealed class PenaltyBlock(Baps owner, bool adaptive) : IBlock<BapsState>
    {
        public string Name => "penalty variance";

        public bool Update(BapsState state, Xoshiro256PlusPlus rng)
        {
            if (!adaptive)
            {
                state.Tau2 = ConjugateUpdates.NormalVariance(owner._priorShape, owner._priorScale, owner._nu - owner._r,
                                                             owner._penalty.Roughness(state.Beta)).Sample(rng);
                return true;
            }
            owner._penalty.SquaredDifferencesInto(state.Beta, state.SquaredDifferences);   // (Δ^rβ)_i²
            double weightedRough = 0.0;
            for (int i = 0; i < state.Eta.Length; i++) weightedRough += state.Eta[i] * state.SquaredDifferences[i];
            state.Tau2 = ConjugateUpdates.NormalVariance(owner._priorShape, owner._priorScale, state.Eta.Length, weightedRough)
                                         .Sample(rng);
            return true;
        }
    }

    // η_i | β, τ² ~ Gamma(φ + ½, φ + (Δβ)_i²/(2τ²)), the normal-precision update of one term under the Gamma(φ, φ) prior
    // (the +½ is the prior normalizer's Jacobian).
    private sealed class LocalBlock(Baps owner) : IBlock<BapsState>
    {
        public string Name => "local multipliers";

        public bool Update(BapsState state, Xoshiro256PlusPlus rng)
        {
            for (int i = 0; i < state.Eta.Length; i++)
                state.Eta[i] = ConjugateUpdates.NormalPrecision(owner._phi, owner._phi, 1.0, state.SquaredDifferences[i] / state.Tau2)
                                               .Sample(rng);
            return true;
        }
    }

    // Random-walk Metropolis on log λ against the REML evidence under the uniform-edf prior, β and σ² integrated out.
    private sealed class LogLambdaBlock(Baps owner) : IBlock<BapsState>
    {
        private readonly Func<double, double> _target = u => owner.MarginalLogTarget(Math.Exp(u));

        public string Name => "log λ (collapsed)";

        public bool Update(BapsState state, Xoshiro256PlusPlus rng)
            => RandomWalkMetropolis.Step(ref state.LogLambda, ref state.Evidence, owner._mhLogStep, _target, rng);
    }

    private double ResidualSS(double[] beta)
    {
        double ss = 0.0;
        for (int i = 0; i < _n; i++)
        {
            double fit = 0.0;
            for (int j = 0; j < _nu; j++) fit += _z[i, j] * beta[j];
            double d = _y[i] - fit;
            ss += d * d;
        }
        return ss;
    }
}
