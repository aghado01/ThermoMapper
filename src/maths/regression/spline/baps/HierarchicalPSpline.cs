using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;
using Maths.Samplers.Blocks;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Baps;

/// <summary>Posterior summary of a hierarchical P-spline fit: the population (Bayes) curve, each replicate's
/// shrunken curve, and the variance-component posteriors (posterior means of σ, τ_u and τ_0). Every summary is read
/// from <see cref="PosteriorDraws"/>, every retained draw per chain in draw order (one chain for now).</summary>
/// <remarks>
/// <see cref="PopulationCoefficients"/> and <see cref="ReplicateCoefficients"/> average the drawn curves.
/// <see cref="RaoBlackwellizedPopulationCoefficients"/> and <see cref="RaoBlackwellizedReplicateCoefficients"/> average
/// instead the conditional means E[β | ·] the sweep computed. Both estimate the same posterior means. The Rao–Blackwellized
/// ones drop the conditional draw's own noise, which usually lowers the Monte Carlo error; for a scan of more than two
/// blocks that is typical rather than guaranteed (Liu, Wong &amp; Kong 1994 prove it for two).
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
    double[][] RaoBlackwellizedReplicateCoefficients)
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
/// would need a sum-to-zero constraint.)
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
    }

    /// <param name="burn">Burn-in sweeps.</param>
    /// <param name="samples">Retained draws.</param>
    /// <param name="seed">The chain's seed.</param>
    /// <param name="retainReplicates">When true (the default), every draw keeps each replicate's curve and conditional
    /// mean, R × ν numbers twice over. When false, draws keep the population curve and the variances only, and the
    /// replicate summaries accumulate as the chain runs; use it when R × ν × draws outgrows memory.</param>
    public HierarchicalResult Run(int burn = 1000, int samples = 2000, int seed = 1, bool retainReplicates = true)
    {
        if (samples < 1 || burn < 0) throw new ArgumentOutOfRangeException(nameof(samples));
        int reps = _responses.Length;
        var scan = new BlockScan<SweepState>(Blocks(), SweepState.Start(reps, _nu), new Xoshiro256PlusPlus(seed),
                                             state => state.Copy());

        var draws = new List<HierarchicalDraw>(samples);
        var betaSum = new double[reps][];                  // the replicate sums, accumulated from the draws when retained
        var meanSum = new double[reps][];                  // and the sums of their conditional means
        for (int rr = 0; rr < reps; rr++) { betaSum[rr] = new double[_nu]; meanSum[rr] = new double[_nu]; }

        int steps = burn + samples;
        for (int it = 0; it < steps; it++)
        {
            SweepState state = scan.Step();
            if (it < burn) continue;

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
                    for (int j = 0; j < _nu; j++)
                    {
                        betaSum[rr][j] += state.Replicates[rr][j];
                        meanSum[rr][j] += state.ReplicateMeans[rr][j];
                    }
            }
            draws.Add(new HierarchicalDraw((double[])state.Population.Clone(), replicates, state.Sigma2, state.TauU2,
                                           state.Tau02, state.PopulationMean, conditionalMeans));
        }

        // Every summary is read from the retained draws, each sum in draw order as the running sums used to run; without
        // retained replicate draws, the replicate sums were accumulated as the chain ran.
        var beta0Sum = new double[_nu];
        var mean0Sum = new double[_nu];
        double sigSum = 0.0, tauUSum = 0.0, tau0Sum = 0.0;
        foreach (HierarchicalDraw d in draws)
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

        for (int j = 0; j < _nu; j++) beta0Sum[j] /= samples;
        for (int j = 0; j < _nu; j++) mean0Sum[j] /= samples;
        var repMeans = new double[reps][];
        var repConditionalMeans = new double[reps][];
        for (int rr = 0; rr < reps; rr++)
        {
            repMeans[rr] = new double[_nu];
            repConditionalMeans[rr] = new double[_nu];
            for (int j = 0; j < _nu; j++) repMeans[rr][j] = betaSum[rr][j] / samples;
            for (int j = 0; j < _nu; j++) repConditionalMeans[rr][j] = meanSum[rr][j] / samples;
        }
        return new HierarchicalResult(beta0Sum, repMeans, sigSum / samples, tauUSum / samples, tau0Sum / samples, samples,
                                      new IReadOnlyList<HierarchicalDraw>[] { draws.ToArray() }, mean0Sum, repConditionalMeans);
    }

    // The five blocks in sweep order, over this fit's constants: ZᵀZ and each Zᵀy_r.
    private IBlock<SweepState>[] Blocks()
    {
        int reps = _responses.Length;
        var gram = new double[_gramBw + 1, _nu];
        _design.Accumulate(null, gram, null, null);
        var zty = new double[reps][];
        for (int rr = 0; rr < reps; rr++)
        {
            zty[rr] = new double[_nu];
            _design.AccumulateRhs(null, _responses[rr], zty[rr]);
        }
        return new IBlock<SweepState>[]
        {
            new ReplicateBlock(this, gram, zty), new PopulationBlock(this), new NoiseBlock(this), new DeviationBlock(this),
            new SmoothingBlock(this),
        };
    }

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

        // Every curve at 0 and every variance at 1.
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
