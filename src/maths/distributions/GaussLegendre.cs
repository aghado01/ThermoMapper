using System;
using System.Collections.Concurrent;

namespace Maths.Distributions;

/// <summary>
/// The n-point Gauss–Legendre rule on [−1, 1], exact for polynomials of degree 2n − 1: nodes by Newton's method on the
/// Legendre polynomial Pₙ from Tricomi's starting values, weights 2/((1 − x²)Pₙ′(x)²), mirrored so the rule is exactly
/// symmetric. Each rule is computed once and shared; <see cref="StudentT"/>'s interval moments integrate with it.
/// </summary>
public sealed class GaussLegendre
{
    private static readonly ConcurrentDictionary<int, GaussLegendre> Rules = new();
    private readonly double[] _nodes, _weights;

    private GaussLegendre(int n)
    {
        _nodes = new double[n];
        _weights = new double[n];
        for (int i = 0; i < (n + 1) / 2; i++)
        {
            double x = Math.Cos(Math.PI * (i + 0.75) / (n + 0.5));
            double slope = 0.0;
            for (int iteration = 0; iteration < 100; iteration++)
            {
                (double p, double dp) = Legendre(n, x);
                slope = dp;
                double step = p / dp;
                x -= step;
                if (Math.Abs(step) <= 1e-16 * Math.Max(1.0, Math.Abs(x))) break;
            }
            slope = Legendre(n, x).Slope;
            double w = 2.0 / ((1.0 - x * x) * slope * slope);
            _nodes[i] = -x;
            _nodes[n - 1 - i] = x;
            _weights[i] = _weights[n - 1 - i] = w;
        }
        if (n % 2 == 1) _nodes[n / 2] = 0.0;
    }

    /// <summary>The <paramref name="n"/>-point rule.</summary>
    public static GaussLegendre Of(int n)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), "A rule needs a node.");
        return Rules.GetOrAdd(n, static k => new GaussLegendre(k));
    }

    /// <summary>The nodes n.</summary>
    public int Count => _nodes.Length;

    /// <summary>The nodes on [−1, 1], ascending.</summary>
    public ReadOnlySpan<double> Nodes => _nodes;

    /// <summary>The weights, which sum to 2.</summary>
    public ReadOnlySpan<double> Weights => _weights;

    // Pₙ(x) and Pₙ′(x) by the three-term recurrence.
    private static (double Value, double Slope) Legendre(int n, double x)
    {
        double p0 = 1.0, p1 = x;
        for (int k = 2; k <= n; k++)
        {
            double p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
            p0 = p1;
            p1 = p2;
        }
        if (n == 0) return (1.0, 0.0);
        return (p1, n * (x * p1 - p0) / (x * x - 1.0));
    }
}
