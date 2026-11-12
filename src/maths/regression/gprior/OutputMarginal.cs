using System;
using Maths.Distributions;
using Maths.LinAlg;

namespace Maths.Regression.GPrior;

/// <summary>
/// The collapsed g-prior block's marginal for K outputs on one design (E34), at one (τ, b): Y (N × K) = 1μᵀ + Z_cA + E with each row of
/// E ∼ N(0, Σ), A | Σ ∼ MN(0, (τG_c)⁻¹, Σ) and the intercept placed by its convention, with A and the intercept integrated against the
/// likelihood raised to b. What is left is the K × K scale S, the determinant term and the count ν that Σ's posterior
/// IW(ν₀ + ν, Ψ₀ + S) adds, the matrix-t's statistics. Every piece is the scalar block's (<see cref="GPriorMarginal"/>) carried to K
/// outputs: S = b(R + τ/(τ + b)·(F + E)), from the least-squares fit's K × K residual form R and fitted form F in the centred basis and
/// the levels' energy E (the zero centre's Nȳȳᵀ; none under level centring), so a small τ on a close fit keeps its digits; the
/// determinant term is K times the scalar one; and ν = Nb, or Nb − 1 under the flat intercept, whose integral spends an observation.
/// At K = 1 it is the scalar block's marginal under IG(ν₀/2, Ψ₀/2), to rounding. A value: the τ updates evaluate it many times.
/// </summary>
public sealed class OutputMarginal
{
    private static readonly double LogPi = Math.Log(Math.PI), LogTwoPi = Math.Log(2.0 * Math.PI);

    private readonly double[,] _s;
    private readonly double[,] _slope;
    private readonly double _determinantSlope;
    private readonly int _spent;   // observations a flat intercept's integral spends: 0 or 1
    private readonly double _count;   // the observations ν counts: N under precision weights, as the scalar block's

    private OutputMarginal(double[,] s, double[,] slope, double logDeterminantTerm, double determinantSlope,
                           int rows, double tau, double inverseTemperature, int spent, double? count = null)
    {
        _s = s;
        _slope = slope;
        _determinantSlope = determinantSlope;
        LogDeterminantTerm = logDeterminantTerm;
        Rows = rows;
        _count = count ?? rows;
        Tau = tau;
        InverseTemperature = inverseTemperature;
        _spent = spent;
    }

    /// <summary>
    /// The marginal from the centred least-squares fit's statistics: <paramref name="centredColumns"/> centred columns beside the
    /// intercept, the K × K residual form R = Σᵢ rᵢrᵢᵀ and fitted form F = B_cᵀG_c⁻¹B_c, and the outputs' means ȳ over
    /// <paramref name="rows"/> observations, with the intercept placed by <paramref name="intercept"/>: level centring, inside the g-prior
    /// and centred on ȳ (or on <paramref name="levels"/>, frozen apart from the response, whose gap's energy then joins F); the zero
    /// centre, inside and centred at 0; or Liang's flat intercept, apart and integrated. The median centre is the robust path's and is
    /// refused, as the scalar block refuses it. <paramref name="weights"/> carries the reading's row weights as the scalar block takes them
    /// (<see cref="GPriorWeights"/>; E34's slice 4): their mass in the intercept's direction, their Σ ln wᵢ (a random-effect block's −ln|V|
    /// among them) once per output, and the count ν reads; the forms and means are then the weighted reading's. Null is unit weights.
    /// </summary>
    public static OutputMarginal CentredAt(int centredColumns, double[,] residualForm, double[,] centredFittedForm, int rows, double[] responseMeans,
                                           GPriorCentre intercept, double tau, double inverseTemperature = 1.0, double[]? levels = null,
                                           GPriorWeights? weights = null)
    {
        ArgumentNullException.ThrowIfNull(residualForm);
        ArgumentNullException.ThrowIfNull(centredFittedForm);
        ArgumentNullException.ThrowIfNull(responseMeans);
        int k = responseMeans.Length;
        if (k < 1) throw new ArgumentException("There is at least one output.", nameof(responseMeans));
        if (residualForm.GetLength(0) != k || residualForm.GetLength(1) != k || centredFittedForm.GetLength(0) != k || centredFittedForm.GetLength(1) != k)
            throw new ArgumentException($"The forms are {k} × {k}, one row and column per output.", nameof(residualForm));
        if (centredColumns < 0) throw new ArgumentOutOfRangeException(nameof(centredColumns));
        if (rows < 1) throw new ArgumentOutOfRangeException(nameof(rows));
        if (levels is not null)
        {
            if (intercept != GPriorCentre.Mean) throw new ArgumentException("A frozen level is level centring's.", nameof(levels));
            if (levels.Length != k) throw new ArgumentException("One level per output is required.", nameof(levels));
        }
        if (!(tau > 0.0) || double.IsPositiveInfinity(tau)) throw new ArgumentOutOfRangeException(nameof(tau), "τ must be positive and finite.");
        double b = inverseTemperature;
        if (!(b >= 0.0) || double.IsPositiveInfinity(b))
            throw new ArgumentOutOfRangeException(nameof(inverseTemperature), "The inverse temperature must be non-negative and finite.");

        var gap = new double[k];   // ȳ less the intercept's centre: 0 under level centring and the flat intercept, ȳ under the zero centre
        for (int a = 0; a < k; a++)
        {
            if (!double.IsFinite(responseMeans[a])) throw new ArgumentOutOfRangeException(nameof(responseMeans));
            gap[a] = intercept switch
            {
                GPriorCentre.Mean => levels is null ? 0.0 : responseMeans[a] - levels[a],
                GPriorCentre.Flat => 0.0,
                GPriorCentre.Zero => responseMeans[a],
                GPriorCentre.Median => throw new ArgumentException(
                    "The median centre is the robust path's, on a basis that carries the constant; a design with its intercept apart takes the mean, "
                    + "the zero centre or the flat intercept.", nameof(intercept)),
                _ => throw new ArgumentOutOfRangeException(nameof(intercept)),
            };
        }
        GPriorWeights w = weights ?? GPriorWeights.Unit(rows);
        double mass = w.Mass, logWeights = k * 0.5 * w.SumLogWeights;   // each output's likelihood carries ½Σ ln wᵢ
        double shrink = tau / (tau + b);
        var s = new double[k, k];
        var slope = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int d = 0; d <= a; d++)
            {
                double r = residualForm[a, d], f = centredFittedForm[a, d];
                if (r != residualForm[d, a] || f != centredFittedForm[d, a]) throw new ArgumentException("The forms are symmetric.", nameof(residualForm));
                s[a, d] = s[d, a] = b * (r + shrink * (f + mass * gap[a] * gap[d]));
                slope[a, d] = slope[d, a] = r + shrink * shrink * (f + mass * gap[a] * gap[d]);
            }
        for (int a = 0; a < k; a++)
            if (!(residualForm[a, a] >= 0.0) || !(centredFittedForm[a, a] >= 0.0)) throw new ArgumentOutOfRangeException(nameof(residualForm), "A form's diagonal is not negative.");

        if (intercept != GPriorCentre.Flat)
            return new OutputMarginal(s, slope, k * (0.5 * (centredColumns + 1) * Math.Log(shrink)) + b * logWeights,
                                      -0.5 * k * (centredColumns + 1) / (tau + b) + logWeights, rows, tau, b, 0, w.Count);
        // The flat intercept's integral (2π)^(K/2)|Σ|^(½)(Nb)^(−K/2): ½ln 2π moves into the count's 2π term, −½ln(Nb) stays here.
        if (!(b > 0.0))
            throw new ArgumentOutOfRangeException(nameof(inverseTemperature),
                "b = 0 is out of reach on the flat intercept's integrated form, whose (2π|Σ|/(Nb))^(K/2) diverges there.");
        return new OutputMarginal(s, slope, k * (0.5 * centredColumns * Math.Log(shrink) - 0.5 * Math.Log(mass * b)) + b * logWeights,
                                  -0.5 * k * (centredColumns / (tau + b) + 1.0 / b) + logWeights, rows, tau, b, 1, w.Count);
    }

    /// <summary>K, the outputs.</summary>
    public int Outputs => _s.GetLength(0);

    /// <summary>Observations N.</summary>
    public int Rows { get; }

    /// <summary>Prior precision τ.</summary>
    public double Tau { get; }

    /// <summary>Likelihood inverse temperature b (1 is untempered).</summary>
    public double InverseTemperature { get; }

    /// <summary>ν, the count Σ's posterior adds: Nb, or Nb − 1 under the flat intercept.</summary>
    public double DegreesOfFreedom => _count * InverseTemperature - _spent;

    /// <summary>The K × K scale S = b(R + τ/(τ + b)·(F + E)), copied.</summary>
    public double[,] S => (double[,])_s.Clone();

    /// <summary>The coefficients' integral: K·(k/2)·ln(τ/(τ + b)) with k the g-prior's columns, the intercept's among them inside it; under
    /// the flat intercept K·(M/2)·ln(τ/(τ + b)) − (K/2)·ln(Nb).</summary>
    public double LogDeterminantTerm { get; }

    // The pencil's form (slice 2), with the reading's weights taken as CentredAt takes them (E34's slice 4): the flat intercept's mass,
    // (Kb/2)Σ ln wᵢ and the count ν reads.
    internal static OutputMarginal FromPencil(double[,] residual, PriorPencil[] pencils, int rows, double tau, double b, GPriorWeights? weights = null)
    {
        GPriorWeights w = weights ?? GPriorWeights.Unit(rows);
        int k = pencils.Length, p = pencils[0].Eigenvalues.Length;
        bool flat = !pencils[0].InterceptInside;
        if (!(tau > 0 && double.IsFinite(tau)) || !(b >= 0 && double.IsFinite(b)) || (flat && b == 0))
            throw new ArgumentOutOfRangeException(nameof(b), "The coefficient precision and temperature must define the block.");
        var s = new double[k, k];
        var slope = new double[k, k];
        double logWeights = k * 0.5 * w.SumLogWeights;
        double determinant = (flat ? -0.5 * k * Math.Log(w.Mass * b) : 0) + b * logWeights;
        double determinantSlope = (flat ? -0.5 * k / b : 0) + logWeights;
        for (int a = 0; a < k; a++)
            for (int d = 0; d < k; d++)
            {
                s[a, d] = b * residual[a, d];
                slope[a, d] = residual[a, d];
            }
        for (int j = 0; j < p; j++)
        {
            double lambda = pencils[0].Eigenvalues[j];
            double shrink = tau / (tau + b * lambda);
            determinant += 0.5 * k * Math.Log(shrink);
            determinantSlope -= 0.5 * k * lambda / (tau + b * lambda);
            for (int a = 0; a < k; a++)
                for (int d = 0; d <= a; d++)
                {
                    double energy = pencils[a].TransformedCross[j] * pencils[d].TransformedCross[j] / lambda;
                    s[a, d] += b * shrink * energy;
                    slope[a, d] += shrink * shrink * energy;
                    s[d, a] = s[a, d];
                    slope[d, a] = slope[a, d];
                }
        }
        return new OutputMarginal(s, slope, determinant, determinantSlope, rows, tau, b, flat ? 1 : 0, w.Count);
    }

    /// <summary>
    /// Σ | Y ∼ IW(ν₀ + ν, Ψ₀ + S) under <paramref name="prior"/> IW(ν₀, Ψ₀), or IW(ν, S) under the improper |Σ|^(−(K + 1)/2) (null).
    /// Refused where it is improper, ν₀ + ν ≤ K − 1: under the improper prior at small b, and under the flat intercept's hottest rungs.
    /// </summary>
    public InverseWishart SigmaPosterior(InverseWishart? prior)
    {
        int k = Outputs;
        CheckPrior(prior);
        double nu = (prior?.DegreesOfFreedom ?? 0.0) + DegreesOfFreedom;
        if (!(nu > k - 1)) throw new InvalidOperationException($"Σ's posterior IW(ν₀ + ν = {nu:R}, ·) is improper on {k} × {k} matrices.");
        var scale = (double[,])_s.Clone();
        if (prior is not null)
        {
            double[,] psi = prior.Scale;
            for (int a = 0; a < k; a++)
                for (int d = 0; d < k; d++) scale[a, d] += psi[a, d];
        }
        return InverseWishart.Of(nu, scale);
    }

    /// <summary>
    /// The complete ln p_b(Y | Z, τ) with A, the intercept and Σ integrated, ν = <see cref="DegreesOfFreedom"/>:
    /// <c>−(νK/2) ln π + LogDeterminantTerm + (ν₀/2) ln|Ψ₀| − ln Γ_K(ν₀/2) + ln Γ_K((ν₀ + ν)/2) − ((ν₀ + ν)/2) ln|Ψ₀ + S|</c> under
    /// IW(ν₀, Ψ₀), exactly 0 at b = 0; under the improper |Σ|^(−(K + 1)/2) (null), up to its constant,
    /// <c>−(νK/2) ln 2π + LogDeterminantTerm + (νK/2) ln 2 + ln Γ_K(ν/2) − (ν/2) ln|S|</c>, which needs ν &gt; K − 1 and S positive definite.
    /// The flat intercept's is up to its prior's constant too. Refused where Σ's posterior is improper.
    /// </summary>
    public double LogMarginal(InverseWishart? prior)
    {
        int k = Outputs;
        CheckPrior(prior);
        double nu = DegreesOfFreedom;
        if (prior is null)
        {
            if (!(nu > k - 1)) throw new InvalidOperationException($"The improper Σ prior needs ν = {nu:R} > {k - 1}.");
            double logDetS = LogDeterminant(_s);
            if (double.IsNaN(logDetS)) throw new InvalidOperationException("The improper Σ prior needs S positive definite.");
            return -0.5 * nu * k * LogTwoPi + LogDeterminantTerm + 0.5 * nu * k * Math.Log(2.0) + SpecialFunctions.LogMultivariateGamma(k, 0.5 * nu) - 0.5 * nu * logDetS;
        }
        double nu0 = prior.DegreesOfFreedom, posterior = nu0 + nu;
        if (!(posterior > k - 1)) throw new InvalidOperationException($"Σ's posterior IW(ν₀ + ν = {posterior:R}, ·) is improper on {k} × {k} matrices.");
        if (InverseTemperature == 0.0) return 0.0;   // no likelihood, and a proper prior integrates to one; the terms below cancel only to rounding
        double[,] psi = prior.Scale;
        for (int a = 0; a < k; a++)
            for (int d = 0; d < k; d++) psi[a, d] += _s[a, d];
        double logDetTotal = LogDeterminant(psi);
        if (double.IsNaN(logDetTotal)) throw new InvalidOperationException("Ψ₀ + S is not positive definite.");
        return -0.5 * nu * k * LogPi + LogDeterminantTerm + 0.5 * nu0 * prior.LogDeterminantScale - SpecialFunctions.LogMultivariateGamma(k, 0.5 * nu0)
               + SpecialFunctions.LogMultivariateGamma(k, 0.5 * posterior) - 0.5 * posterior * logDetTotal;
    }

    /// <summary>The derivative of the integrated log likelihood with respect to b, for path diagnostics.</summary>
    public double LogLikelihoodSlope(InverseWishart? prior)
    {
        CheckPrior(prior);
        int k = Outputs;
        double posterior = (prior?.DegreesOfFreedom ?? 0) + DegreesOfFreedom;
        if (!(posterior > k - 1)) throw new InvalidOperationException("Σ's posterior is improper.");
        double[,] total = prior?.Scale ?? new double[k, k];
        for (int a = 0; a < k; a++)
            for (int d = 0; d < k; d++) total[a, d] += _s[a, d];
        var factor = new DenseCholesky(k);
        factor.Decompose(total);
        if (factor.FailedAt >= 0) throw new InvalidOperationException("Σ's posterior scale is not positive definite.");
        double digamma = 0, trace = 0;
        var column = new double[k];
        for (int a = 0; a < k; a++)
        {
            digamma += SpecialFunctions.Digamma(0.5 * (posterior - a));
            for (int d = 0; d < k; d++) column[d] = _slope[d, a];
            factor.SolveInPlace(column);
            trace += column[a];
        }
        return _determinantSlope + 0.5 * _count * (digamma - k * LogPi - factor.LogDet) - 0.5 * posterior * trace;
    }

    // ln|A| by Cholesky, its lower triangle read; NaN when A is not positive definite.
    private static double LogDeterminant(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        double logDet = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = a[i, j];
                for (int r = 0; r < j; r++) sum -= l[i, r] * l[j, r];
                if (i != j)
                {
                    l[i, j] = sum / l[j, j];
                    continue;
                }
                if (!(sum > 0.0) || double.IsPositiveInfinity(sum)) return double.NaN;
                l[i, i] = Math.Sqrt(sum);
                logDet += Math.Log(sum);
            }
        return logDet;
    }

    private void CheckPrior(InverseWishart? prior)
    {
        if (prior is not null && prior.Dimension != Outputs)
            throw new ArgumentException($"The prior is on {prior.Dimension} × {prior.Dimension} matrices, the marginal on {Outputs} outputs.", nameof(prior));
    }
}
