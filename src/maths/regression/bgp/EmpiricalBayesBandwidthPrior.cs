using System;

namespace Maths.Regression.Bgp;

/// <summary>
/// The empirical-Bayes prior on the GP kernel bandwidth t (Tang, Wu, Cheng &amp; Dunson 2025, §4.3, eqs. 12–14) —
/// the paper's central contribution. <c>p(t) ∝ t^{−a₀} exp(−b₀ / v̂_n(t))</c> on the support <c>(γ₁ T_n², 1]</c>
/// and 0 elsewhere, where v̂_n is the kernel affinity (<see cref="GpRegression.KernelAffinity"/>) and T_n is the
/// averaged k-NN distance. Both statistics carry the intrinsic dimension d implicitly — v̂_n(t) ∼ t^{d/2} and
/// T_n² ∼ n^{−2/d} — so the prior concentrates t near the dimension-adaptive scale n^{−2/(2s+d)} <em>without ever
/// estimating d</em> (their Prop. 4.4). It is an observable over the data that induces a density over the
/// bandwidth field; the normalizer is dropped because the downstream Metropolis–Hastings sampler never needs it.
/// </summary>
public sealed class EmpiricalBayesBandwidthPrior
{
    private readonly GpRegression _gp;
    private readonly double _a0;
    private readonly double _b0;
    private readonly double _lowerBound;   // γ₁ T_n²

    /// <summary>The averaged k-NN distance T_n (scales as n^{−1/d}).</summary>
    public double Tn { get; }

    /// <summary>Lower edge of the support, γ₁ T_n²; the prior is 0 at or below it.</summary>
    public double LowerBound => _lowerBound;

    /// <summary>The neighbor order k actually used.</summary>
    public int K { get; }

    /// <param name="x">n × D predictors — the same data wrapped by <paramref name="gp"/>.</param>
    /// <param name="gp">The GP core supplying the kernel-affinity statistic v̂_n(t).</param>
    /// <param name="a0">Tail exponent a₀ &gt; 0 (arbitrary hyperparameter).</param>
    /// <param name="b0">Affinity weight b₀ &gt; 0 (arbitrary hyperparameter).</param>
    /// <param name="gamma1">Support-floor constant γ₁ (paper uses 1/4).</param>
    /// <param name="gamma2">Neighbor-order constant γ₂ (paper uses 1/4); k = ⌈γ₂ (ln n)²⌉, or 2 when n &lt; 200.</param>
    /// <param name="k">Explicit neighbor order, overriding the γ₂ rule (e.g. for scaling diagnostics).</param>
    public EmpiricalBayesBandwidthPrior(double[,] x, GpRegression gp,
        double a0 = 1.0, double b0 = 1.0, double gamma1 = 0.25, double gamma2 = 0.25, int? k = null)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(gp);
        if (!(a0 > 0.0)) throw new ArgumentOutOfRangeException(nameof(a0));
        if (!(b0 > 0.0)) throw new ArgumentOutOfRangeException(nameof(b0));
        if (!(gamma1 > 0.0)) throw new ArgumentOutOfRangeException(nameof(gamma1));
        int n = x.GetLength(0), dim = x.GetLength(1);
        if (gp.Count != n) throw new ArgumentException("gp must wrap the same data as x.", nameof(gp));
        _gp = gp; _a0 = a0; _b0 = b0;

        K = k ?? (n < 200 ? 2 : (int)Math.Ceiling(gamma2 * Math.Log(n) * Math.Log(n)));
        K = Math.Clamp(K, 2, n);

        // T_n = (1/n) Σ_i R̂_k(X_i); R̂_k(x) = distance to the k-th nearest neighbor (x is its own 1st NN, dist 0).
        var d2 = new double[n];
        double sum = 0.0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                double s = 0.0;
                for (int c = 0; c < dim; c++) { double diff = x[i, c] - x[j, c]; s += diff * diff; }
                d2[j] = s;
            }
            Array.Sort(d2);
            sum += Math.Sqrt(d2[K - 1]);   // k-th smallest (index 0 = self = 0)
        }
        Tn = sum / n;
        _lowerBound = gamma1 * Tn * Tn;
    }

    /// <summary>
    /// log p(t) up to an additive constant (the normalizer Ẑ_n is dropped — only ratios are used). Returns −∞
    /// outside the support (γ₁ T_n², 1].
    /// </summary>
    public double LogDensity(double t)
    {
        if (t <= _lowerBound || t > 1.0) return double.NegativeInfinity;
        return -_a0 * Math.Log(t) - _b0 / _gp.KernelAffinity(t);
    }
}
