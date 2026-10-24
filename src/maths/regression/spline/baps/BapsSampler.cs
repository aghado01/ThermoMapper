using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;
using Maths.Samplers.Ensemble;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Baps;

/// <summary>How the smoothing parameter λ is sampled in <see cref="BapsSampler"/>.</summary>
public enum BapsLambdaUpdate
{
    /// <summary>
    /// Conjugate Gibbs over the variance components: <c>β | σ²,τ²</c> from its banded Gaussian full conditional,
    /// <c>σ² | β</c> and <c>τ² | β</c> from inverse-gamma full conditionals (λ = σ²/τ² derived). Full Bayesian.
    /// </summary>
    Gibbs,

    /// <summary>
    /// Metropolis on <c>log λ</c> against the profiled REML marginal evidence (β integrated out) under a uniform prior
    /// on the effective degrees of freedom (<see cref="PenalizedSpline.EdfSlope"/>). Each retained λ then draws σ² and β
    /// exactly from their conditionals, and the coefficient mean is Rao–Blackwellized on β̂(λ). Marginalizes the spline.
    /// The evidence levels off as λ → ∞, so a flat prior on <c>log λ</c> would be improper; the edf prior decays like
    /// 1/λ and makes the target proper.
    /// </summary>
    MarginalEvidence,

    /// <summary>
    /// Locally-adaptive Gibbs (the "A" in BAPS): as <see cref="Gibbs"/> but each difference carries its own mean-1
    /// multiplier <c>η_i ~ Gamma(φ, φ)</c>, so the smoothing <c>λ_i = λ_global·η_i</c> varies across the domain
    /// — heavy where the fit is flat, light where it must track a feature. Every conditional stays conjugate
    /// (the prior normalizer factorizes as <c>(Πη_i)|DDᵀ|</c>, giving each <c>η_i</c> a Gamma full conditional).
    /// </summary>
    AdaptiveGibbs
}

/// <summary>Posterior summary of a BAPS run: the Bayes fit coefficients and the smoothing / noise posteriors.</summary>
/// <remarks>
/// <paramref name="LambdaMedian"/>/<paramref name="LambdaLo"/>/<paramref name="LambdaHi"/> are the 50%, 2.5% and 97.5%
/// quantiles of the pooled draws of the global smoothing (in the adaptive mode, λ_global). A median, not a mean,
/// because λ can be heavy-tailed enough that its posterior mean does not exist. <paramref name="LocalSmoothing"/> is
/// the posterior-mean local smoothing field λ_i (length ν−r) — non-null only for <see cref="BapsLambdaUpdate.AdaptiveGibbs"/>; its dips
/// mark where the fit was allowed to track features. <paramref name="RHatLogLambda"/> is classic R̂ on log λ across chains
/// with its <see cref="DiagnosticState"/> (unavailable for a single chain). <paramref name="PosteriorDraws"/> keeps
/// every retained draw (<see cref="BapsDraw"/>), per chain in draw order; every summary above is read from it.
/// <paramref name="LogLambdaConvergence"/> and <paramref name="LogSigma2Convergence"/> are the ensemble shell's
/// rank-normalized diagnostics (split R̂, bulk and tail ESS) on log λ and log σ², the functionals an early stop waits
/// on. <paramref name="AcceptanceRate"/> is the sampling phase's: 1 for the Gibbs modes, the Metropolis rate on log λ
/// for the marginal-evidence mode. <paramref name="ChainSeeds"/> are the resolved per-chain seeds, and
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

    /// <summary>τ²'s posterior median and 95% interval, read from the retained draws (τ² = σ²/λ in every mode).</summary>
    public PosteriorInterval Tau2Interval => Interval(d => d.Tau2);

    /// <summary>
    /// The pointwise band on the fitted curve at the rows of <paramref name="design"/> (any rows of the coefficients'
    /// basis), read from the retained draws. <see cref="PointwiseBand.Mean"/> is the design times
    /// <see cref="Coefficients"/>, the affine readout, Rao–Blackwellized in the marginal-evidence mode. The sd and the
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
/// (λ = σ²/τ²), and the adaptive mode's local multipliers η.
/// <list type="bullet">
///   <item>In the Gibbs modes all of them come from the end of the sweep.</item>
///   <item>The marginal-evidence mode samples λ, then draws σ² and β exactly from their conditionals given λ, and
///   sets τ² = σ²/λ. <see cref="ConditionalMean"/> holds β̂(λ) = E[β | λ, y] beside the draw.</item>
///   <item><see cref="Eta"/> is set in the adaptive mode only, where the local smoothing is λᵢ = λ·ηᵢ.</item>
/// </list>
/// </summary>
public sealed record BapsDraw(double[] Coefficients, double[]? ConditionalMean, double Lambda, double Sigma2, double Tau2,
                              double[]? Eta);

/// <summary>
/// The BAPS (Bayesian adaptive P-spline) sampler — the measure-side dual of free-knot BARS. A rich fixed B-spline
/// basis with a difference penalty whose strength λ is inferred, sampled either by conjugate Gibbs over the
/// variance components or by Metropolis against the REML marginal evidence (<see cref="BapsLambdaUpdate"/>). Every
/// inner solve is banded (<c>A = ZᵀZ + λP</c> via <see cref="BandCholesky"/>'s root-free LDLᵀ, whose factor also
/// supplies the Gaussian β-draw). As a client of the ensemble shell (<see cref="ChainEnsemble"/>) its chains run in
/// parallel from one seed tree, with rank-normalized R̂ and ESS on log λ and log σ², an optional early stop, and
/// checkpoints. Every draw is retained, and the Bayes fit, the λ posterior (median + 95% interval), σ and the bands are
/// read from the draws.
/// </summary>
public sealed class BapsSampler
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
    private readonly double _mhLogStep;  // Metropolis log-λ proposal sd (marginal-evidence mode)
    private readonly double _phi;        // adaptivity: η_i ~ Gamma(φ, φ); φ→∞ recovers global smoothing

    public BapsSampler(double[,] design, double[] y, DifferencePenalty penalty,
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

    /// <param name="mode">How λ is sampled.</param>
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
    /// uninterrupted run would have; the result covers every draw since burn-in. It needs the same data, settings, mode,
    /// master seed and chain count.</param>
    public BapsResult Run(BapsLambdaUpdate mode, int chains = 4, int masterSeed = 1, int burn = 500, int samples = 1500,
                          double rHatTarget = 0.0, int batchSize = 0, double minEssPerChain = ChainEnsemble.MinEssPerChain,
                          bool checkpoint = false, EnsembleCheckpoint? resumeFrom = null)
    {
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (samples < 2 || burn < 0) throw new ArgumentOutOfRangeException(nameof(samples));
        if (mode is not (BapsLambdaUpdate.Gibbs or BapsLambdaUpdate.MarginalEvidence or BapsLambdaUpdate.AdaptiveGibbs))
            throw new ArgumentOutOfRangeException(nameof(mode));

        // The shell gives each chain a kernel stream and a readout stream. Only the marginal-evidence mode's conditional
        // draws use the readout stream, so the λ chain does not depend on them.
        (EnsembleRun summary, IReadOnlyList<IChainRun<BapsDraw>> runs) = new ChainEnsemble().Run<BapsDraw>(
            new BapsEnsemble(this, mode), chains, masterSeed, burn, samples, rHatTarget, batchSize, minEssPerChain,
            checkpoint, resumeFrom);

        // Every summary is read from the retained draws: each chain's sums in draw order, then over chains, and the pooled
        // λ sample for its quantiles. R̂ on log λ is the shell's classic reading of the same draws; the rank-normalized
        // diagnostics come from its traces of log λ and log σ².
        double total = chains * (double)summary.SamplesUsed;
        var coef = new double[_nu];
        double sigSum = 0.0;
        var allLam = new double[(int)total];
        int nLoc = _nu - _r;
        double[]? localField = mode == BapsLambdaUpdate.AdaptiveGibbs && nLoc > 0 ? new double[nLoc] : null;
        var posteriorDraws = new IReadOnlyList<BapsDraw>[chains];
        int idx = 0;
        for (int c = 0; c < chains; c++)
        {
            List<BapsDraw> draws = ((BapsChainRun)runs[c]).Draws;
            posteriorDraws[c] = draws.ToArray();
            var coefSum = new double[_nu];
            double sigmaSum = 0.0;
            double[]? localSum = localField is null ? null : new double[nLoc];
            foreach (BapsDraw d in draws)
            {
                double[] coefficients = d.ConditionalMean ?? d.Coefficients;   // Rao–Blackwellized where the mode has it
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
    /// The ensemble-shell client for one mode: it starts and restores chains, fingerprints the run for checkpoints, and
    /// names the functionals the shell monitors. Those are log λ and log σ², which get split R̂ and bulk and tail ESS and
    /// gate an early stop.
    /// </summary>
    private sealed class BapsEnsemble : IEnsembleModel<BapsDraw>
    {
        private readonly BapsSampler _owner;
        private readonly BapsLambdaUpdate _mode;

        public BapsEnsemble(BapsSampler owner, BapsLambdaUpdate mode)
        {
            _owner = owner;
            _mode = mode;
        }

        public int FunctionalDim => 2;

        public int MonitoredDim => 2;

        public IChainRun<BapsDraw> StartChain(int chainSeed, int readoutSeed) =>
            new BapsChainRun(_owner, _mode, new Xoshiro256PlusPlus(chainSeed), new Xoshiro256PlusPlus(readoutSeed));

        /// <summary>The design, the responses, the penalty order, the hyperparameters and the mode.</summary>
        public string Fingerprint
        {
            get
            {
                var f = new RunFingerprint().Add("BAPS").Add((int)_mode).Add(_owner._n).Add(_owner._nu).Add(_owner._r);
                foreach (double v in _owner._z) f.Add(v);
                return f.Add(_owner._y).Add(_owner._priorShape).Add(_owner._priorScale).Add(_owner._mhLogStep)
                        .Add(_owner._phi).ToString();
            }
        }

        public IChainRun<BapsDraw> RestoreChain(IChainState state) => state is BapsChainState s
            ? new BapsChainRun(_owner, _mode, s)
            : throw new ArgumentException("Not a BAPS chain's restart state.", nameof(state));
    }

    /// <summary>One BAPS chain's restart state: both streams' words, the state its next sweep conditions on (σ², τ² and η
    /// for the Gibbs modes; log λ and its target value for the marginal-evidence mode), its counters and its retained
    /// draws.</summary>
    private sealed record BapsChainState(
        (ulong, ulong, ulong, ulong) Rng, (ulong, ulong, ulong, ulong) Readout, double Sigma2, double Tau2, double[] Eta,
        double LogLambda, double Evidence, long Accepted, long Attempts, List<BapsDraw> Draws) : IChainState;

    /// <summary>
    /// One BAPS chain: its mode's sweep, its retained draws and its restart state. Burn-in advances without recording;
    /// each sampling step advances once and records the draw, so the chain's streams are the same however its steps
    /// are split into rounds or runs.
    /// </summary>
    private sealed class BapsChainRun : IChain<BapsDraw>, IResumableChainRun<BapsDraw>
    {
        private readonly BapsSampler _o;
        private readonly BapsLambdaUpdate _mode;
        private readonly Xoshiro256PlusPlus _rng;       // the sweep's stream
        private readonly Xoshiro256PlusPlus _readout;   // the marginal-evidence mode's conditional draws
        private readonly double[,] _band;
        private readonly BandCholesky _chol;
        private readonly double[] _z, _innov, _beta, _sqDiff, _weights;
        private readonly double[] _eta;                 // adaptive mode's local multipliers
        private double _sigma2 = 1.0, _tau2 = 1.0;      // the Gibbs modes' variance components
        private double _logLam, _ev;                    // the marginal-evidence chain's state and its log target

        public readonly List<BapsDraw> Draws;           // every post-burn draw, in order

        public BapsChainRun(BapsSampler owner, BapsLambdaUpdate mode, Xoshiro256PlusPlus rng, Xoshiro256PlusPlus readout)
            : this(owner, mode, rng, readout, new List<BapsDraw>())
        {
            Array.Fill(_eta, 1.0);   // start unadapted (= global)
            if (mode == BapsLambdaUpdate.MarginalEvidence) _ev = owner.MarginalLogTarget(1.0);   // log λ = 0
        }

        public BapsChainRun(BapsSampler owner, BapsLambdaUpdate mode, BapsChainState s)
            : this(owner, mode, Restore(s.Rng), Restore(s.Readout), new List<BapsDraw>(s.Draws))
        {
            Array.Copy(s.Eta, _eta, _eta.Length);
            _sigma2 = s.Sigma2;
            _tau2 = s.Tau2;
            _logLam = s.LogLambda;
            _ev = s.Evidence;
            Accepted = s.Accepted;
            Attempts = s.Attempts;
        }

        private BapsChainRun(BapsSampler owner, BapsLambdaUpdate mode, Xoshiro256PlusPlus rng, Xoshiro256PlusPlus readout,
                             List<BapsDraw> draws)
        {
            _o = owner;
            _mode = mode;
            _rng = rng;
            _readout = readout;
            Draws = draws;
            int nu = owner._nu, nLoc = owner._nu - owner._r;
            _band = new double[owner._bw + 1, nu];
            _chol = new BandCholesky(nu, owner._bw, BandFactorization.Ldlt);
            _z = new double[nu];
            _innov = new double[nu];
            _beta = new double[nu];
            _sqDiff = new double[nLoc];
            _weights = new double[nLoc];   // λ_i = (σ²/τ²)·η_i
            _eta = new double[nLoc];
        }

        private static Xoshiro256PlusPlus Restore((ulong, ulong, ulong, ulong) w) =>
            Xoshiro256PlusPlus.Restore(w.Item1, w.Item2, w.Item3, w.Item4);

        // IChain
        public long Accepted { get; private set; }
        public long Attempts { get; private set; }

        public BapsDraw Step()
        {
            Advance();
            return Record();
        }

        // IChainRun
        public IChain<BapsDraw> Chain => this;

        public void Burn(int steps)
        {
            for (int i = 0; i < steps; i++) Advance();
        }

        public void Accumulate(in BapsDraw draw, Span<double> functionals)
        {
            functionals[0] = Math.Log(draw.Lambda);
            functionals[1] = Math.Log(draw.Sigma2);
            Draws.Add(draw);
        }

        public IChainState SaveState() => new BapsChainState(
            _rng.SaveState(), _readout.SaveState(), _sigma2, _tau2, (double[])_eta.Clone(), _logLam, _ev,
            Accepted, Attempts, new List<BapsDraw>(Draws));

        private void Advance()
        {
            switch (_mode)
            {
                case BapsLambdaUpdate.Gibbs: AdvanceGibbs(); break;
                case BapsLambdaUpdate.MarginalEvidence: AdvanceMarginal(); break;
                default: AdvanceAdaptiveGibbs(); break;
            }
            Attempts++;
        }

        // The draw at the end of the sweep. The marginal-evidence mode adds its exact conjugate draws given λ here, on the
        // readout stream: σ² | λ, y ~ IG((n − r)/2, RSS_pen/2) under the Jeffreys prior the evidence assumes, then
        // β | σ², λ, y ~ N(β̂, σ²A⁻¹), with β̂(λ) kept for the Rao–Blackwellized coefficient mean.
        private BapsDraw Record()
        {
            if (_mode != BapsLambdaUpdate.MarginalEvidence)
                return new BapsDraw((double[])_beta.Clone(), null, _sigma2 / _tau2, _sigma2, _tau2,
                                    _mode == BapsLambdaUpdate.AdaptiveGibbs ? (double[])_eta.Clone() : null);

            int nu = _o._nu;
            double lambda = Math.Exp(_logLam);
            (double[] betaHat, double rssPen, BandCholesky factor) = _o._model.SolveFactored(lambda);
            double sigma2 = InverseGamma.ShapeScale(0.5 * (_o._n - _o._r), 0.5 * rssPen).Sample(_readout);
            for (int p = 0; p < nu; p++) _z[p] = Normal.SampleStandard(_readout);
            factor.SampleInnovation(_z, _innov);                 // v ~ N(0, A⁻¹)
            double sigma = Math.Sqrt(sigma2);
            var beta = new double[nu];
            for (int p = 0; p < nu; p++) beta[p] = betaHat[p] + sigma * _innov[p];
            return new BapsDraw(beta, betaHat, lambda, sigma2, sigma2 / lambda, null);
        }

        // ── Gibbs: β | σ²,τ²  then  σ² | β  then  τ² | β ───────────────────────────
        private void AdvanceGibbs()
        {
            int nu = _o._nu;
            double lambda = _sigma2 / _tau2;
            Array.Copy(_o._gram0, _band, _o._gram0.Length);       // A = ZᵀZ + λP
            _o._penalty.AccumulateInto(_band, nu, lambda);
            _chol.DecomposeBanded(_band);
            if (_chol.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                throw new InvalidOperationException($"ZᵀZ + λP is not numerically positive definite at λ = {lambda} ({_chol.Status}).");
            double[] betaHat = _chol.Solve(_o._zty);             // posterior mean A⁻¹Zᵀy

            for (int p = 0; p < nu; p++) _z[p] = Normal.SampleStandard(_rng);
            _chol.SampleInnovation(_z, _innov);                  // v ~ N(0, A⁻¹)
            double sigma = Math.Sqrt(_sigma2);
            for (int p = 0; p < nu; p++) _beta[p] = betaHat[p] + sigma * _innov[p];

            _sigma2 = InverseGamma.ShapeScale(_o._priorShape + 0.5 * _o._n, _o._priorScale + 0.5 * _o.ResidualSS(_beta)).Sample(_rng);
            _tau2 = InverseGamma.ShapeScale(_o._priorShape + 0.5 * (nu - _o._r), _o._priorScale + 0.5 * _o._penalty.Roughness(_beta)).Sample(_rng);
            Accepted++;
        }

        // ── Marginal evidence: Metropolis on log λ against REML evidence × uniform-edf prior ──
        private void AdvanceMarginal()
        {
            double prop = _logLam + _o._mhLogStep * Normal.SampleStandard(_rng);
            double evProp = _o.MarginalLogTarget(Math.Exp(prop));
            if (Math.Log(1.0 - _rng.NextDouble()) < evProp - _ev)
            {
                _logLam = prop;
                _ev = evProp;
                Accepted++;
            }
        }

        // ── Adaptive Gibbs: as Gibbs, plus per-difference multipliers η_i ~ Gamma(φ, φ) ──
        private void AdvanceAdaptiveGibbs()
        {
            int nu = _o._nu, nLoc = _eta.Length;
            double lamGlobal = _sigma2 / _tau2;
            for (int i = 0; i < nLoc; i++) _weights[i] = lamGlobal * _eta[i];
            Array.Copy(_o._gram0, _band, _o._gram0.Length);       // A = ZᵀZ + Dᵀdiag(λ_i)D
            _o._penalty.AccumulateInto(_band, nu, _weights);
            _chol.DecomposeBanded(_band);
            if (_chol.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                throw new InvalidOperationException($"ZᵀZ + Dᵀdiag(λᵢ)D is not numerically positive definite at global λ = {lamGlobal} ({_chol.Status}).");
            double[] betaHat = _chol.Solve(_o._zty);

            for (int p = 0; p < nu; p++) _z[p] = Normal.SampleStandard(_rng);
            _chol.SampleInnovation(_z, _innov);
            double sigma = Math.Sqrt(_sigma2);
            for (int p = 0; p < nu; p++) _beta[p] = betaHat[p] + sigma * _innov[p];

            _o._penalty.SquaredDifferencesInto(_beta, _sqDiff);   // (Δ^rβ)_i²

            _sigma2 = InverseGamma.ShapeScale(_o._priorShape + 0.5 * _o._n, _o._priorScale + 0.5 * _o.ResidualSS(_beta)).Sample(_rng);

            double weightedRough = 0.0;                           // Σ η_i (Δβ)_i²
            for (int i = 0; i < nLoc; i++) weightedRough += _eta[i] * _sqDiff[i];
            _tau2 = InverseGamma.ShapeScale(_o._priorShape + 0.5 * nLoc, _o._priorScale + 0.5 * weightedRough).Sample(_rng);

            // η_i | β, τ² ~ Gamma(φ + ½, φ + (Δβ)_i²/(2τ²)) — the +½ is the prior-normalizer Jacobian.
            for (int i = 0; i < nLoc; i++)
                _eta[i] = Gamma.ShapeRate(_o._phi + 0.5, _o._phi + 0.5 * _sqDiff[i] / _tau2).Sample(_rng);
            Accepted++;
        }
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
