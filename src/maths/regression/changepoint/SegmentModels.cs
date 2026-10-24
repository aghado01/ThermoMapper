using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Changepoint;

/// <summary>
/// The per-segment model behind the exact change-point DP (Hutter 2005, "any noise / level prior"): for a segment
/// covering data points a..b it supplies the single-segment evidence (the level integrated out) and the posterior
/// level mean/variance (for the regression curve). The DP itself is model-agnostic — it only reads these tables —
/// so swapping a Gaussian likelihood for a heavy-tailed (outlier-robust) one is a model swap, not a DP change.
/// </summary>
public interface ISegmentModel
{
    /// <summary>Number of data points n.</summary>
    int Count { get; }

    /// <summary>
    /// Fill the upper triangle (1 ≤ a ≤ b ≤ n): <paramref name="logEvidence"/>[a,b] = log D(a,b) (level integrated
    /// out), <paramref name="levelMean"/>[a,b] = E[μ | segment], <paramref name="levelVar"/>[a,b] = Var[μ | segment].
    /// Tables are (n+1)×(n+1).
    /// </summary>
    void Fill(double[,] logEvidence, double[,] levelMean, double[,] levelVar);

    /// <summary>
    /// A draw of segment a..b's level (1-based, inclusive) from its posterior given that it is a segment, consistent with
    /// <see cref="Fill"/>'s evidence and moments. Callable without <see cref="Fill"/>.
    /// </summary>
    double SampleLevel(int a, int b, Xoshiro256PlusPlus rng);
}

/// <summary>
/// Gaussian noise + Gaussian segment-level prior — the conjugate, closed-form model. The single-segment evidence
/// is the marginal <c>N(y_{a..b}; ν1, σ²I + ρ²11ᵀ)</c>, and the posterior level is the standard conjugate
/// Gaussian; both are O(1) per segment from prefix sums.
/// </summary>
public sealed class GaussianSegmentModel : ISegmentModel
{
    private readonly double[] _y;
    private readonly double _sigma2;
    private readonly double _nu;
    private readonly double _rho2;

    public GaussianSegmentModel(double[] y, double sigma2, double nu, double rho2)
    {
        ArgumentNullException.ThrowIfNull(y);
        if (!(sigma2 > 0.0)) throw new ArgumentOutOfRangeException(nameof(sigma2));
        if (!(rho2 > 0.0)) throw new ArgumentOutOfRangeException(nameof(rho2));
        _y = y; _sigma2 = sigma2; _nu = nu; _rho2 = rho2;
    }

    public int Count => _y.Length;

    public void Fill(double[,] logEvidence, double[,] levelMean, double[,] levelVar)
    {
        int n = _y.Length;
        var rPre = new double[n + 1];
        var qPre = new double[n + 1];
        var yPre = new double[n + 1];
        for (int i = 1; i <= n; i++)
        {
            double ri = _y[i - 1] - _nu;
            rPre[i] = rPre[i - 1] + ri;
            qPre[i] = qPre[i - 1] + ri * ri;
            yPre[i] = yPre[i - 1] + _y[i - 1];
        }

        double logSigma2 = Math.Log(_sigma2);
        double half2Pi = 0.5 * Math.Log(2.0 * Math.PI);
        double invRho2 = 1.0 / _rho2;

        for (int a = 1; a <= n; a++)
            for (int b = a; b <= n; b++)
            {
                int len = b - a + 1;
                double sr = rPre[b] - rPre[a - 1];
                double sq = qPre[b] - qPre[a - 1];
                double t = _sigma2 + len * _rho2;
                logEvidence[a, b] = -len * (half2Pi + 0.5 * logSigma2)
                                    - 0.5 * (Math.Log(t) - logSigma2)
                                    - 0.5 / _sigma2 * (sq - _rho2 / t * sr * sr);
                double v = 1.0 / (invRho2 + len / _sigma2);
                levelVar[a, b] = v;
                levelMean[a, b] = (_nu * invRho2 + (yPre[b] - yPre[a - 1]) / _sigma2) * v;
            }
    }

    /// <summary>The conjugate draw: N((ν/ρ² + Σy/σ²)·v, v) with v = 1/(1/ρ² + L/σ²).</summary>
    public double SampleLevel(int a, int b, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (a < 1 || b < a || b > _y.Length) throw new ArgumentOutOfRangeException(nameof(b));
        double sum = 0.0;
        for (int i = a; i <= b; i++) sum += _y[i - 1];
        double v = 1.0 / (1.0 / _rho2 + (b - a + 1) / _sigma2);
        return (_nu / _rho2 + sum / _sigma2) * v + Math.Sqrt(v) * Normal.SampleStandard(rng);
    }
}

/// <summary>
/// Cauchy (heavy-tailed) noise + Gaussian level prior — the outlier-robust model (Hutter: "Cauchy can handle
/// outliers"). A single gross outlier contributes only boundedly to a Cauchy likelihood, so it neither inflates
/// the segment level nor spawns a spurious one-point segment. No conjugacy, so the single-segment evidence and
/// posterior level moments are computed by 1-D quadrature over the level μ on a shared uniform grid (<c>scale</c> is
/// the Cauchy noise width γ, the robust analogue of σ). The grid's range and spacing are chosen so that every
/// segment's integral is resolved and nothing material is cut off:
/// <list type="bullet">
///   <item><b>Range</b> <c>[min(ν, y_min) − 6ρ, max(ν, y_max) + 6ρ]</c>. Beyond it the prior and every Cauchy factor
///   both decrease, so each segment loses at most the prior's 6ρ tail ratio, about 2·10⁻⁹ of its mass per side.</item>
///   <item><b>Spacing</b> <c>h = ½/√C</c> with <c>C = 1/ρ² + 2n/γ²</c>, a bound on the log integrand's curvature (a
///   Cauchy factor's log bends by at most 2/γ² either way), so the narrowest segment posterior, about γ/√(2n) wide,
///   spans at least two nodes. <see cref="GridPoints"/> is the resulting count, never below
///   <see cref="MinAutomaticGridPoints"/>.</item>
///   <item><b>Windowing.</b> Per segment the log integrand is read at every 16th node; with its curvature bounded by
///   <c>1/ρ² + 2L/γ²</c> for a segment of L points, that bounds the maximum between them, and only the intervals that
///   can come within 50 nats of it are summed node by node. What is skipped is below e⁻⁵⁰ per node.</item>
/// </list>
/// The per-point log-likelihoods are prefix-summed at every node, so each node of a segment costs O(1); the table holds
/// <c>G·(n+1)</c> numbers, at most <see cref="MaxTableEntries"/>. <c>gridPoints</c> overrides the count on the same
/// range; <see cref="Resolution"/> then reports how coarse it is.
/// </summary>
public sealed class CauchySegmentModel : ISegmentModel
{
    /// <summary>The fewest nodes the automatic grid uses (the former fixed default).</summary>
    public const int MinAutomaticGridPoints = 512;

    /// <summary>The largest prefix table, in numbers (nodes × (n+1)), the model builds before refusing the grid.</summary>
    public const long MaxTableEntries = 1L << 25;

    private const int Stride = 16;          // the window's coarse read: every Stride-th node
    private const double Cutoff = 50.0;     // nodes this many nats below the segment's maximum are skipped

    private readonly double[] _y;
    private readonly double _scale;
    private readonly double _nu;
    private readonly double _rho2;
    private readonly double _lo;
    private readonly double _step;
    private readonly int _grid;
    private readonly double _curvature;
    private readonly Lazy<(double[] Cum, double[] LogPrior)> _table;

    /// <param name="gridPoints">Null (the default) sizes the grid from the data, the scale and the prior (see the
    /// summary); a value fixes the node count on the same range, at least 8.</param>
    public CauchySegmentModel(double[] y, double scale, double nu, double rho2, int? gridPoints = null)
    {
        ArgumentNullException.ThrowIfNull(y);
        if (!(scale > 0.0)) throw new ArgumentOutOfRangeException(nameof(scale));
        if (!(rho2 > 0.0)) throw new ArgumentOutOfRangeException(nameof(rho2));
        if (gridPoints is int fixedCount && fixedCount < 8) throw new ArgumentOutOfRangeException(nameof(gridPoints));
        _y = y; _scale = scale; _nu = nu; _rho2 = rho2;

        double rho = Math.Sqrt(rho2), yMin = nu, yMax = nu;
        foreach (double v in y) { yMin = Math.Min(yMin, v); yMax = Math.Max(yMax, v); }
        _lo = yMin - 6.0 * rho;
        double hi = yMax + 6.0 * rho;
        _curvature = 1.0 / rho2 + 2.0 * Math.Max(1, y.Length) / (scale * scale);
        double nodes = gridPoints is int count
            ? count
            : Math.Max(MinAutomaticGridPoints, Math.Ceiling((hi - _lo) * 2.0 * Math.Sqrt(_curvature)) + 1.0);
        if (!(nodes * (y.Length + 1) <= MaxTableEntries))
            throw new ArgumentOutOfRangeException(gridPoints is null ? nameof(scale) : nameof(gridPoints),
                $"The quadrature needs {nodes:G3} nodes × {y.Length + 1} prefix entries; widen the scale or pass fewer gridPoints.");
        _grid = (int)nodes;
        _step = (hi - _lo) / (_grid - 1);
        _table = new Lazy<(double[], double[])>(BuildTable);
    }

    public int Count => _y.Length;

    /// <summary>The quadrature's node count.</summary>
    public int GridPoints => _grid;

    /// <summary>The node spacing against the log integrand's curvature scale, <c>h·√C</c>: at most ½ on the automatic
    /// grid. Errors grow quickly above about 1, where the narrowest posterior falls between nodes.</summary>
    public double Resolution => _step * Math.Sqrt(_curvature);

    public void Fill(double[,] logEvidence, double[,] levelMean, double[,] levelVar)
    {
        int n = _y.Length;
        (double[] cum, double[] logPrior) = _table.Value;
        var coarse = new double[(_grid - 2) / Stride + 2];
        var ranges = new List<(int Start, int End)>();
        for (int a = 1; a <= n; a++)
            for (int b = a; b <= n; b++)
            {
                Window(cum, logPrior, a, b, coarse, ranges);
                int lower = (a - 1) * _grid, upper = b * _grid;
                double top = double.NegativeInfinity;
                foreach ((int start, int end) in ranges)
                    for (int h = start; h < end; h++) top = Math.Max(top, cum[upper + h] - cum[lower + h] + logPrior[h]);

                // Σ e, Σ (μ−ν) e and Σ (μ−ν)² e with e = exp(l − top); centring at ν keeps the variance free of cancellation.
                double s0 = 0.0, s1 = 0.0, s2 = 0.0;
                foreach ((int start, int end) in ranges)
                    for (int h = start; h < end; h++)
                    {
                        double e = Math.Exp(cum[upper + h] - cum[lower + h] + logPrior[h] - top);
                        double d = _lo + h * _step - _nu;
                        s0 += e; s1 += d * e; s2 += d * d * e;
                    }
                logEvidence[a, b] = top + Math.Log(s0 * _step);            // log ∫ exp(·) dμ
                double offset = s1 / s0;
                levelMean[a, b] = _nu + offset;
                levelVar[a, b] = Math.Max(0.0, s2 / s0 - offset * offset);
            }
    }

    /// <summary>
    /// A draw from the quadrature's posterior over the level: a node, in proportion to its weight in the segment's
    /// integral. Expectations over the draws are the quadrature's own, so they agree with <see cref="Fill"/>'s moments and
    /// are as accurate as its rule.
    /// </summary>
    public double SampleLevel(int a, int b, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (a < 1 || b < a || b > _y.Length) throw new ArgumentOutOfRangeException(nameof(b));
        (double[] cum, double[] logPrior) = _table.Value;
        var ranges = new List<(int Start, int End)>();
        Window(cum, logPrior, a, b, new double[(_grid - 2) / Stride + 2], ranges);
        int lower = (a - 1) * _grid, upper = b * _grid;
        double top = double.NegativeInfinity;
        foreach ((int start, int end) in ranges)
            for (int h = start; h < end; h++) top = Math.Max(top, cum[upper + h] - cum[lower + h] + logPrior[h]);
        double total = 0.0;
        foreach ((int start, int end) in ranges)
            for (int h = start; h < end; h++) total += Math.Exp(cum[upper + h] - cum[lower + h] + logPrior[h] - top);

        double target = rng.NextDouble() * total, acc = 0.0;
        int node = ranges[0].Start;
        foreach ((int start, int end) in ranges)
            for (int h = start; h < end; h++)
            {
                acc += Math.Exp(cum[upper + h] - cum[lower + h] + logPrior[h] - top);
                node = h;
                if (acc >= target) return _lo + node * _step;
            }
        return _lo + node * _step;
    }

    // cum[i·G + h] = Σ_{j≤i} log Cauchy(y_j | μ_h), node-contiguous so a segment's nodes are read in order, and the log
    // prior at each node.
    private (double[] Cum, double[] LogPrior) BuildTable()
    {
        int n = _y.Length, g = _grid;
        double logScalePi = Math.Log(Math.PI * _scale);
        double logPriorNorm = -0.5 * Math.Log(2.0 * Math.PI * _rho2);
        var cum = new double[(n + 1) * g];
        var logPrior = new double[g];
        for (int h = 0; h < g; h++)
        {
            double m = _lo + h * _step, d = m - _nu;
            logPrior[h] = logPriorNorm - 0.5 * d * d / _rho2;
            double c = 0.0;
            for (int i = 1; i <= n; i++)
            {
                double z = (_y[i - 1] - m) / _scale;
                c += -logScalePi - Math.Log(1.0 + z * z);
                cum[i * g + h] = c;
            }
        }
        return (cum, logPrior);
    }

    // The node ranges, [Start, End), that carry segment [a..b]'s integral. The log integrand is read at every Stride-th
    // node (and the last); between two reads it exceeds the larger by at most C·H²/8, with C = 1/ρ² + 2L/γ² its curvature
    // bound and H the read spacing. So an interval can hold a node within Cutoff of the maximum only if its larger read is
    // within Cutoff + C·H²/8 of the largest read, and only those intervals are kept.
    private void Window(double[] cum, double[] logPrior, int a, int b, double[] coarse, List<(int Start, int End)> ranges)
    {
        int last = _grid - 1, reads = (last + Stride - 1) / Stride + 1;
        int lower = (a - 1) * _grid, upper = b * _grid;
        double best = double.NegativeInfinity;
        for (int k = 0; k < reads; k++)
        {
            int h = Math.Min(k * Stride, last);
            coarse[k] = cum[upper + h] - cum[lower + h] + logPrior[h];
            best = Math.Max(best, coarse[k]);
        }
        double spacing = Stride * _step;
        double slack = (1.0 / _rho2 + 2.0 * (b - a + 1) / (_scale * _scale)) * spacing * spacing / 8.0;
        double threshold = best - Cutoff - slack;

        ranges.Clear();
        for (int k = 0; k + 1 < reads; k++)
        {
            if (Math.Max(coarse[k], coarse[k + 1]) < threshold) continue;
            int start = k * Stride, end = k + 2 == reads ? last + 1 : (k + 1) * Stride;
            if (ranges.Count > 0 && ranges[^1].End == start) ranges[^1] = (ranges[^1].Start, end);
            else ranges.Add((start, end));
        }
    }
}
