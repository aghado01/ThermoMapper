using System;
using Maths.LinAlg;

namespace Maths.Regression.GPrior;

/// <summary>
/// A centred design whose coefficient prior stands apart from its likelihood's Gram, decomposed: under weights, D03's robust
/// row, a prior precision τD that reads the unweighted columns against a likelihood Gram A that reads the weighted ones. With
/// D = LLᵀ and L⁻¹AL⁻ᵀ = VΛVᵀ, the posterior precision bA + τD is LV(bΛ + τ)VᵀLᵀ, so the coefficient integral and the shrinkage
/// are <see cref="RidgeTerms"/> over the λⱼ at O(p) per (τ, b), as the normalized ridge's are over its own. Unit weights make A
/// and D one matrix, every λⱼ = 1, which is the g-form's closed form. Each structure costs its own O(p³) decomposition: no
/// candidate is bordered on its parent, as <see cref="RidgeSpectrum"/>'s are, since a change of D re-triangulates L.
/// Immutable: a design's current structure shares one with every state that holds it.
/// </summary>
public sealed class PriorPencil : RidgeColumns
{
    private const double JacobiTolerance = 1e-15;

    private readonly int _columns;
    private readonly double[] _lambda;      // eigenvalues of L⁻¹AL⁻ᵀ, descending
    private readonly double[][] _vectors;   // _vectors[j][i]: component i of eigenvector j
    private readonly double[,] _inverseRoot;   // L⁻¹, lower triangular
    private readonly double[] _w;           // VᵀL⁻¹c
    private readonly double[] _q;           // wⱼ²/λⱼ
    private readonly double _logPriorDeterminant;   // ln|D|
    private readonly double[,] _likelihood, _prior;   // A and D, kept for the factor's exact condition

    private PriorPencil(int columns, bool interceptInside, double[] lambda, double[][] vectors, double[,] inverseRoot, double[] w, double logPriorDeterminant,
                        double[,] likelihood, double[,] prior)
    {
        _likelihood = likelihood;
        _prior = prior;
        _columns = columns;
        InterceptInside = interceptInside;
        _lambda = lambda;
        _vectors = vectors;
        _inverseRoot = inverseRoot;
        _w = w;
        _logPriorDeterminant = logPriorDeterminant;
        _q = new double[lambda.Length];
        for (int j = 0; j < lambda.Length; j++) _q[j] = w[j] * w[j] / lambda[j];
    }

    /// <summary>
    /// The pencil of <paramref name="likelihood"/> (A, positive definite) and <paramref name="prior"/> (D), with
    /// <paramref name="cross"/> the fit's cross products; null when D or the pencil is not positive definite, which puts the
    /// structure outside F. The order is <paramref name="centredColumns"/>, one more when
    /// <paramref name="interceptInside"/> puts the intercept first, inside the prior.
    /// </summary>
    public static PriorPencil? Of(double[,] likelihood, double[,] prior, ReadOnlySpan<double> cross, int centredColumns, bool interceptInside)
    {
        ArgumentNullException.ThrowIfNull(likelihood);
        ArgumentNullException.ThrowIfNull(prior);
        int p = centredColumns + (interceptInside ? 1 : 0);
        if (likelihood.GetLength(0) != p || likelihood.GetLength(1) != p || prior.GetLength(0) != p || prior.GetLength(1) != p)
            throw new ArgumentException("The two matrices must have the pencil's order.", nameof(likelihood));
        if (cross.Length != p) throw new ArgumentException("One cross product per direction is required.", nameof(cross));
        if (p == 0) return new PriorPencil(0, false, [], [], new double[0, 0], [], 0.0, likelihood, prior);
        if (InverseRoot(prior, out double logDeterminant) is not { } t) return null;

        // L⁻¹AL⁻ᵀ, scaled to a unit largest diagonal for the eigensolver's absolute tolerance.
        var scratch = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
            {
                double sum = 0.0;
                for (int k = 0; k <= i; k++) sum += t[i, k] * likelihood[k, j];
                scratch[i, j] = sum;
            }
        var unit = new double[p * p];
        double top = 0.0;
        for (int i = 0; i < p; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = 0.0;
                for (int k = 0; k <= j; k++) sum += scratch[i, k] * t[j, k];
                unit[i * p + j] = unit[j * p + i] = sum;
                if (i == j) top = Math.Max(top, sum);
            }
        if (!(top > 0.0) || double.IsPositiveInfinity(top)) return null;
        for (int k = 0; k < unit.Length; k++) unit[k] /= top;
        EigenResult eigen = DenseEigen.DecomposeSymmetric(unit, p, tol: JacobiTolerance);
        var lambda = new double[p];
        var whitened = new double[p];
        for (int i = 0; i < p; i++)
        {
            double sum = 0.0;
            for (int k = 0; k <= i; k++) sum += t[i, k] * cross[k];
            whitened[i] = sum;
        }
        var w = new double[p];
        for (int j = 0; j < p; j++)
        {
            lambda[j] = top * eigen.Eigenvalues[j];
            if (!(lambda[j] > 0.0)) return null;
            double[] v = eigen.Eigenvectors[j];
            double sum = 0.0;
            for (int i = 0; i < p; i++) sum += v[i] * whitened[i];
            w[j] = sum;
        }
        return new PriorPencil(centredColumns, interceptInside, lambda, eigen.Eigenvectors, t, w, logDeterminant, likelihood, prior);
    }

    /// <summary>This pencil with the cross products <paramref name="cross"/> in place of its own: the same A and D under another
    /// response, in O(p²).</summary>
    public PriorPencil WithCross(ReadOnlySpan<double> cross)
    {
        int p = _lambda.Length;
        if (cross.Length != p) throw new ArgumentException("One cross product per direction is required.", nameof(cross));
        Span<double> whitened = p <= 128 ? stackalloc double[p] : new double[p];
        for (int i = 0; i < p; i++)
        {
            double sum = 0.0;
            for (int k = 0; k <= i; k++) sum += _inverseRoot[i, k] * cross[k];
            whitened[i] = sum;
        }
        var w = new double[p];
        for (int j = 0; j < p; j++)
        {
            double[] v = _vectors[j];
            double sum = 0.0;
            for (int i = 0; i < p; i++) sum += v[i] * whitened[i];
            w[j] = sum;
        }
        return new PriorPencil(_columns, InterceptInside, _lambda, _vectors, _inverseRoot, w, _logPriorDeterminant, _likelihood, _prior);
    }

    // L⁻¹ for D = LLᵀ, and ln|D|; null when D is not positive definite.
    private static double[,]? InverseRoot(double[,] d, out double logDeterminant)
    {
        int p = d.GetLength(0);
        var l = new double[p, p];
        logDeterminant = 0.0;
        for (int i = 0; i < p; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = d[i, j];
                for (int k = 0; k < j; k++) sum -= l[i, k] * l[j, k];
                if (i != j)
                {
                    l[i, j] = sum / l[j, j];
                    continue;
                }
                if (!(sum > 0.0) || double.IsPositiveInfinity(sum)) return null;
                l[i, i] = Math.Sqrt(sum);
                logDeterminant += Math.Log(sum);
            }
        var t = new double[p, p];
        for (int i = 0; i < p; i++)
        {
            t[i, i] = 1.0 / l[i, i];
            for (int j = 0; j < i; j++)
            {
                double sum = 0.0;
                for (int k = j; k < i; k++) sum -= l[i, k] * t[k, j];
                t[i, j] = sum / l[i, i];
            }
        }
        return t;
    }

    /// <inheritdoc/>
    public override int Count => _columns;

    /// <summary>Whether the intercept is the first direction, inside the prior; a flat intercept stands apart.</summary>
    public bool InterceptInside { get; }

    /// <summary>The pencil's eigenvalues λⱼ, descending: each column's information under the likelihood against the prior's.</summary>
    public ReadOnlySpan<double> Eigenvalues => _lambda;

    private protected override RidgeTerms Evaluate(double tau, double b) => SpectralTerms(_lambda, _q, tau, b);

    // Eigenbasis coordinates to the pencil's own, L⁻ᵀV·x.
    private void Raw(ReadOnlySpan<double> eigenCoordinates, Span<double> result)
    {
        int p = _lambda.Length;
        Span<double> x = p <= 128 ? stackalloc double[p] : new double[p];
        for (int i = 0; i < p; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < p; j++) sum += _vectors[j][i] * eigenCoordinates[j];
            x[i] = sum;
        }
        for (int i = 0; i < p; i++)
        {
            double sum = 0.0;
            for (int k = i; k < p; k++) sum += _inverseRoot[k, i] * x[k];
            result[i] = sum;
        }
    }

    /// <summary>The precision bA + τD, as a factor read through the eigenvectors.</summary>
    internal ISpdFactor Factor(double tau, double b) => new SpectralFactor(this, tau, b);

    /// <summary>The posterior mean about the prior's centre, b(bA + τD)⁻¹c, into <paramref name="mean"/>.</summary>
    internal void PosteriorMean(double tau, double b, Span<double> mean)
    {
        int p = _lambda.Length;
        Span<double> x = p <= 128 ? stackalloc double[p] : new double[p];
        for (int j = 0; j < p; j++) x[j] = b * _w[j] / (b * _lambda[j] + tau);
        Raw(x, mean);
    }

    private sealed class SpectralFactor : ISpdFactor
    {
        private readonly PriorPencil _pencil;
        private readonly double[] _precision;   // bλⱼ + τ
        private readonly double _tau, _b;

        public SpectralFactor(PriorPencil pencil, double tau, double b)
        {
            _pencil = pencil;
            _tau = tau;
            _b = b;
            int p = pencil._lambda.Length;
            _precision = new double[p];
            double logDet = pencil._logPriorDeterminant, min = double.PositiveInfinity, max = double.NegativeInfinity;
            for (int j = 0; j < p; j++)
            {
                _precision[j] = b * pencil._lambda[j] + tau;
                logDet += Math.Log(_precision[j]);
                min = Math.Min(min, _precision[j]);
                max = Math.Max(max, _precision[j]);
            }
            LogDet = logDet;
            Status = new FactorStatus(true, 0, min, max, p == 0 ? 1.0 : min / max);   // the whitened precision's eigenvalues stand in for pivots
        }

        public int Dimension => _pencil._lambda.Length;

        public FactorStatus Status { get; }

        public double LogDet { get; }

        public double[] Solve(ReadOnlySpan<double> b)
        {
            double[] x = b.ToArray();
            SolveInPlace(x);
            return x;
        }

        public void SolveInPlace(Span<double> x)
        {
            int p = Dimension;
            if (x.Length != p) throw new ArgumentException("The right-hand side must have the factor's order.", nameof(x));
            Span<double> whitened = p <= 128 ? stackalloc double[p] : new double[p];
            for (int i = 0; i < p; i++)
            {
                double sum = 0.0;
                for (int k = 0; k <= i; k++) sum += _pencil._inverseRoot[i, k] * x[k];
                whitened[i] = sum;
            }
            Span<double> z = p <= 128 ? stackalloc double[p] : new double[p];
            for (int j = 0; j < p; j++)
            {
                double sum = 0.0;
                double[] v = _pencil._vectors[j];
                for (int i = 0; i < p; i++) sum += v[i] * whitened[i];
                z[j] = sum / _precision[j];
            }
            _pencil.Raw(z, x);
        }

        public void SampleInnovation(ReadOnlySpan<double> z, Span<double> result)
        {
            int p = Dimension;
            if (z.Length != p || result.Length != p) throw new ArgumentException("The draw must have the factor's order.");
            Span<double> scaled = p <= 128 ? stackalloc double[p] : new double[p];
            for (int j = 0; j < p; j++) scaled[j] = z[j] / Math.Sqrt(_precision[j]);
            _pencil.Raw(scaled, result);
        }

        public double EquilibratedCondition()
        {
            // Exact: E = S·A·S with S = diag(A)^(−½), and E⁻¹ = S⁻¹A⁻¹S⁻¹, A⁻¹ read by its columns' solves.
            int p = Dimension;
            if (p == 0) return 1.0;
            var inverse = new double[p, p];
            var column = new double[p];
            for (int k = 0; k < p; k++)
            {
                Array.Clear(column);
                column[k] = 1.0;
                SolveInPlace(column);
                for (int i = 0; i < p; i++) inverse[i, k] = column[i];
            }
            double normE = 0.0, normInverse = 0.0;
            for (int k = 0; k < p; k++)
            {
                double sumE = 0.0, sumInverse = 0.0;
                for (int i = 0; i < p; i++)
                {
                    double scale = Math.Sqrt(Precision(i, i) * Precision(k, k));
                    sumE += Math.Abs(Precision(i, k)) / scale;
                    sumInverse += Math.Abs(inverse[i, k]) * scale;
                }
                normE = Math.Max(normE, sumE);
                normInverse = Math.Max(normInverse, sumInverse);
            }
            return normE * normInverse;
        }

        private double Precision(int i, int k) => _b * _pencil._likelihood[i, k] + _tau * _pencil._prior[i, k];
    }
}

/// <summary>
/// The collapsed block under a coefficient prior apart from the likelihood's Gram (<see cref="PriorPencil"/>), in the coordinates
/// the design holds, (α, β_c) with the columns centred under W: any (τ, b) in O(p). With the intercept inside the prior the
/// coefficients are jointly drawn, since the prior couples α to the columns through the centring; a flat intercept stands apart,
/// α ~ N(ȳ_w, σ²/(Σwᵢ·b)) given σ².
/// </summary>
public sealed class CollapsedPencilPrior
{
    private readonly PriorPencil _pencil;
    private readonly double _rss, _responseMean, _level;
    private readonly int _rows;
    private readonly GPriorWeights _weights;

    private CollapsedPencilPrior(PriorPencil pencil, double rss, int rows, GPriorWeights weights, double responseMean, double level, GPriorCentre intercept)
    {
        _pencil = pencil;
        _rss = rss;
        _rows = rows;
        _weights = weights;
        _responseMean = responseMean;
        _level = level;
        Intercept = intercept;
    }

    /// <summary>
    /// The block for <paramref name="pencil"/>, whose weighted fit leaves <paramref name="residualSumOfSquares"/> over
    /// <paramref name="rows"/> observations, with ȳ_w = <paramref name="responseMean"/>, the prior centred on
    /// <paramref name="level"/> (read under the intercept inside the prior) and <paramref name="intercept"/> its convention.
    /// </summary>
    public static CollapsedPencilPrior FromCentred(PriorPencil pencil, double residualSumOfSquares, int rows, GPriorWeights weights, double responseMean,
                                                   double level, GPriorCentre intercept)
    {
        ArgumentNullException.ThrowIfNull(pencil);
        if (!(residualSumOfSquares >= 0.0)) throw new ArgumentOutOfRangeException(nameof(residualSumOfSquares));
        return new CollapsedPencilPrior(pencil, residualSumOfSquares, rows, weights, responseMean, level, intercept);
    }

    /// <summary>The intercept's convention.</summary>
    public GPriorCentre Intercept { get; }

    /// <summary>Coefficients: the intercept, then the centred columns.</summary>
    public int Columns => _pencil.Count + 1;

    /// <summary>The marginal at (<paramref name="tau"/>, <paramref name="inverseTemperature"/>)
    /// (<see cref="CollapsedGPrior.CentredPencilMarginalAt"/>).</summary>
    public GPriorMarginal MarginalAt(double tau, double inverseTemperature = 1.0)
        => CollapsedGPrior.CentredPencilMarginalAt(_pencil, _rss, _rows, _weights, Intercept, tau, inverseTemperature);

    /// <summary>
    /// The conditional posterior at (<paramref name="tau"/>, <paramref name="inverseTemperature"/>): the coefficients inside the
    /// prior ~ N(centre + b·A_P⁻¹c, σ²A_P⁻¹) with A_P = bA + τD, a factor read through the pencil.
    /// </summary>
    public GPriorConditional At(double tau, double inverseTemperature = 1.0)
    {
        GPriorMarginal marginal = MarginalAt(tau, inverseTemperature);
        double b = inverseTemperature;
        var mean = new double[Columns];
        if (_pencil.InterceptInside)
        {
            _pencil.PosteriorMean(tau, b, mean);
            mean[0] += _level;
            return new GPriorConditional(_pencil.Factor(tau, b), mean, marginal, 1.0, Columns);
        }
        mean[0] = _responseMean;
        _pencil.PosteriorMean(tau, b, mean.AsSpan(1));
        return new GPriorConditional(_pencil.Factor(tau, b), mean, marginal, 1.0, Columns, 1.0 / b, _weights.Mass);
    }
}
