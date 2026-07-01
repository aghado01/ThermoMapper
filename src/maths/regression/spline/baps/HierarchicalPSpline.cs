using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Baps;

/// <summary>Posterior summary of a hierarchical P-spline fit: the population (Bayes) curve, each replicate's
/// shrunken curve, and the variance-component posteriors.</summary>
public sealed record HierarchicalResult(
    double[] PopulationCoefficients,
    double[][] ReplicateCoefficients,
    double NoiseSd,
    double DeviationSd,
    double PopulationSmoothingSd,
    int Draws);

/// <summary>
/// Hierarchical P-spline for multiple related curves (Behseta–Kass–Wallstrom 2005, the penalized form): replicates
/// <c>y_r = Z β_r + ε</c> share a smooth population mean <c>β_0</c> (penalized, <c>τ_0² P⁻</c>) with replicate
/// deviations <c>β_r ~ N(β_0, τ_u² I)</c>, so each curve borrows strength from the others — shrunk toward the
/// population where its own data is weak. A fully conjugate banded Gibbs on the BAPS machinery: <c>β_r | ·</c> and
/// <c>β_0 | ·</c> are banded Gaussian draws (<see cref="BandCholesky"/> factor + <see cref="BandCholesky.SampleInnovation"/>),
/// the three variance components <c>(σ², τ_u², τ_0²)</c> inverse-gamma (reciprocal of <see cref="Gamma"/>).
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

    public HierarchicalResult Run(int burn = 1000, int samples = 2000, int seed = 1)
    {
        if (samples < 1 || burn < 0) throw new ArgumentOutOfRangeException(nameof(samples));
        var rng = new Xoshiro256PlusPlus(seed);
        int reps = _responses.Length;

        var gram = new double[_gramBw + 1, _nu];          // ZᵀZ — constant
        _design.Accumulate(null, gram, null, null);
        var zty = new double[reps][];                      // Zᵀy_r — constant
        for (int rr = 0; rr < reps; rr++)
        {
            zty[rr] = new double[_nu];
            _design.AccumulateRhs(null, _responses[rr], zty[rr]);
        }

        var beta = new double[reps][];
        for (int rr = 0; rr < reps; rr++) beta[rr] = new double[_nu];
        var beta0 = new double[_nu];
        double sigma2 = 1.0, tauU2 = 1.0, tau02 = 1.0;

        var beta0Sum = new double[_nu];
        var betaSum = new double[reps][];
        for (int rr = 0; rr < reps; rr++) betaSum[rr] = new double[_nu];
        double sigSum = 0.0, tauUSum = 0.0, tau0Sum = 0.0;

        var qrBand = new double[_gramBw + 1, _nu];
        var q0Band = new double[_r + 1, _nu];
        var cholR = new BandCholesky(_nu, _gramBw, BandFactorization.Ldlt);
        var chol0 = new BandCholesky(_nu, _r, BandFactorization.Ldlt);
        var rhs = new double[_nu];
        var z = new double[_nu];
        var innov = new double[_nu];

        int steps = burn + samples;
        for (int it = 0; it < steps; it++)
        {
            // Q_r = ZᵀZ/σ² + I/τ_u² — shared by all replicates, factor once.
            for (int d = 0; d <= _gramBw; d++)
                for (int j = 0; j < _nu; j++) qrBand[d, j] = gram[d, j] / sigma2;
            for (int j = 0; j < _nu; j++) qrBand[0, j] += 1.0 / tauU2;
            cholR.DecomposeBanded(qrBand);

            double sumResid = 0.0, sumDev = 0.0;
            for (int rr = 0; rr < reps; rr++)
            {
                for (int j = 0; j < _nu; j++) rhs[j] = zty[rr][j] / sigma2 + beta0[j] / tauU2;
                double[] mr = cholR.Solve(rhs);            // E[β_r | ·] = Q_r⁻¹ rhs
                for (int j = 0; j < _nu; j++) z[j] = StandardNormal(rng);
                cholR.SampleInnovation(z, innov);          // ~ N(0, Q_r⁻¹)
                double[] br = beta[rr];
                for (int j = 0; j < _nu; j++) br[j] = mr[j] + innov[j];

                sumResid += ResidualSS(_responses[rr], br);
                for (int j = 0; j < _nu; j++) { double d = br[j] - beta0[j]; sumDev += d * d; }
            }

            // Q_0 = (R/τ_u²) I + (1/τ_0²) P — identity term keeps it proper despite P's null space.
            Array.Clear(q0Band);
            for (int j = 0; j < _nu; j++) q0Band[0, j] = reps / tauU2;
            _penalty.AccumulateInto(q0Band, _nu, 1.0 / tau02);
            chol0.DecomposeBanded(q0Band);
            for (int j = 0; j < _nu; j++) rhs[j] = 0.0;
            for (int rr = 0; rr < reps; rr++)
                for (int j = 0; j < _nu; j++) rhs[j] += beta[rr][j];
            for (int j = 0; j < _nu; j++) rhs[j] /= tauU2;
            double[] m0 = chol0.Solve(rhs);
            for (int j = 0; j < _nu; j++) z[j] = StandardNormal(rng);
            chol0.SampleInnovation(z, innov);
            for (int j = 0; j < _nu; j++) beta0[j] = m0[j] + innov[j];

            sigma2 = InvGamma(rng, _a0 + 0.5 * reps * _n, _b0 + 0.5 * sumResid);
            tauU2 = InvGamma(rng, _a0 + 0.5 * reps * _nu, _b0 + 0.5 * sumDev);
            tau02 = InvGamma(rng, _a0 + 0.5 * (_nu - _r), _b0 + 0.5 * _penalty.Roughness(beta0));

            if (it >= burn)
            {
                for (int j = 0; j < _nu; j++) beta0Sum[j] += beta0[j];
                for (int rr = 0; rr < reps; rr++)
                    for (int j = 0; j < _nu; j++) betaSum[rr][j] += beta[rr][j];
                sigSum += Math.Sqrt(sigma2);
                tauUSum += Math.Sqrt(tauU2);
                tau0Sum += Math.Sqrt(tau02);
            }
        }

        for (int j = 0; j < _nu; j++) beta0Sum[j] /= samples;
        var repMeans = new double[reps][];
        for (int rr = 0; rr < reps; rr++)
        {
            repMeans[rr] = new double[_nu];
            for (int j = 0; j < _nu; j++) repMeans[rr][j] = betaSum[rr][j] / samples;
        }
        return new HierarchicalResult(beta0Sum, repMeans, sigSum / samples, tauUSum / samples, tau0Sum / samples, samples);
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

    private static double InvGamma(Xoshiro256PlusPlus rng, double shape, double scale)
        => 1.0 / Gamma.Sample(rng, shape, 1.0 / scale);

    private static double StandardNormal(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
