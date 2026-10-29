using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>The order a cut SMC run visits its emulator draws in: the tour (the default), or a random order, the
/// tour's comparison.</summary>
public enum CutSmcOrder
{
    /// <summary>Nearest neighbour then 2-opt on the spread of the log-likelihood differences
    /// (<see cref="CutSmcOrdering.Tour"/>).</summary>
    Tour,

    /// <summary>The first leg's draw, then the rest in a uniformly random order.</summary>
    Random,
}

/// <summary>
/// The order of a cut SMC run's emulator draws (D09, stage 8; Mathews et al.'s travelling-salesman ordering), so that
/// successive conditionals stay close. The distance d(s, s′) is the weighted spread (standard deviation) of
/// ℓ(θ_i; e_s) − ℓ(θ_i; e_s′) over the first leg's final particles: the quantity a bridge's conditional ESS falls with, so
/// a short tour is a cheap one. It costs one S × N evaluation, kept as a matrix, then S²·N arithmetic. The tour is an
/// open path from the first leg's draw, by nearest neighbour then 2-opt. Every μ_s is targeted exactly whatever the
/// order, so letting the order depend on the population changes the run's efficiency only.
/// <para>The recorded fallback, not built: the L₂ distance between the draws' predictions on a fixed 64-point
/// space-filling design in the unit θ box, cheaper but measured where the posterior has no mass.</para>
/// </summary>
public static class CutSmcOrdering
{
    /// <summary>
    /// The S × N matrix ℓ(θ_i; e_s) of each target's untempered log-likelihood at each particle, every pair evaluated
    /// on a fresh copy of the particle so the population's caches stay its own.
    /// </summary>
    public static double[,] LogLikelihoods(IReadOnlyList<CalibrationTarget> targets, IReadOnlyList<CalibrationState> particles, int maxDegreeOfParallelism = 1)
    {
        ArgumentNullException.ThrowIfNull(targets);
        ArgumentNullException.ThrowIfNull(particles);
        var l = new double[targets.Count, particles.Count];
        Parallel.For(0, particles.Count, new ParallelOptions { MaxDegreeOfParallelism = maxDegreeOfParallelism }, i =>
        {
            CalibrationState p = particles[i];
            var fresh = new CalibrationState(p.Unit.ToArray(), p.Levels.ToArray(), p.Variances.ToArray(), p.Coefficients.ToArray());
            for (int s = 0; s < targets.Count; s++) l[s, i] = targets[s].EvaluateState(fresh).LogLikelihoodValue;
        });
        return l;
    }

    /// <summary>
    /// d(s, s′) = √Σ_i W_i·(D_i − D̄)², D_i = ℓ_{s,i} − ℓ_{s′,i}, D̄ = Σ_i W_i·D_i, for every pair of rows of
    /// <paramref name="logLikelihoods"/> under the normalized <paramref name="weights"/>. Each row is centred first, so the
    /// differences' spread is read without the cancellation of the rows' large common level.
    /// </summary>
    public static double[,] Distances(double[,] logLikelihoods, ReadOnlySpan<double> weights)
    {
        ArgumentNullException.ThrowIfNull(logLikelihoods);
        int s = logLikelihoods.GetLength(0), n = logLikelihoods.GetLength(1);
        if (weights.Length != n) throw new ArgumentException("One weight per particle.", nameof(weights));
        var centred = new double[s, n];
        for (int r = 0; r < s; r++)
        {
            double mean = 0.0;
            for (int i = 0; i < n; i++)
            {
                if (!double.IsFinite(logLikelihoods[r, i])) throw new ArgumentException("The log-likelihoods must be finite.", nameof(logLikelihoods));
                mean += weights[i] * logLikelihoods[r, i];
            }
            for (int i = 0; i < n; i++) centred[r, i] = logLikelihoods[r, i] - mean;
        }
        var d = new double[s, s];
        for (int a = 0; a < s; a++)
            for (int b = 0; b < a; b++)
            {
                double sum = 0.0;
                for (int i = 0; i < n; i++)
                {
                    double e = centred[a, i] - centred[b, i];
                    sum += weights[i] * e * e;
                }
                d[a, b] = d[b, a] = Math.Sqrt(sum);
            }
        return d;
    }

    /// <summary>The tour over the rows of <paramref name="distances"/> from <paramref name="start"/>: nearest neighbour,
    /// then 2-opt with the start held.</summary>
    public static int[] Tour(double[,] distances, int start) => TwoOpt(distances, NearestNeighbour(distances, start));

    /// <summary>The open path from <paramref name="start"/> that moves each time to the nearest unvisited node, the
    /// lowest index among ties.</summary>
    public static int[] NearestNeighbour(double[,] distances, int start)
    {
        int s = Check(distances);
        if ((uint)start >= (uint)s) throw new ArgumentOutOfRangeException(nameof(start));
        var path = new int[s];
        var visited = new bool[s];
        path[0] = start;
        visited[start] = true;
        for (int k = 1; k < s; k++)
        {
            int from = path[k - 1], next = -1;
            for (int c = 0; c < s; c++)
                if (!visited[c] && (next < 0 || distances[from, c] < distances[from, next])) next = c;
            path[k] = next;
            visited[next] = true;
        }
        return path;
    }

    /// <summary>
    /// 2-opt on an open path whose first node is held: reverse a segment p[i..j], i ≥ 1, whenever that shortens the path,
    /// until no reversal does by more than rounding. It never lengthens the path.
    /// </summary>
    public static int[] TwoOpt(double[,] distances, int[] path)
    {
        int s = Check(distances);
        ArgumentNullException.ThrowIfNull(path);
        if (path.Length != s) throw new ArgumentException("The path must visit every node.", nameof(path));
        int[] p = (int[])path.Clone();
        double tolerance = 1e-12 * (1.0 + Length(distances, p));
        bool improved = true;
        while (improved)
        {
            improved = false;
            for (int i = 1; i < s - 1; i++)
                for (int j = i + 1; j < s; j++)
                {
                    double before = distances[p[i - 1], p[i]], after = distances[p[i - 1], p[j]];
                    if (j + 1 < s)
                    {
                        before += distances[p[j], p[j + 1]];
                        after += distances[p[i], p[j + 1]];
                    }
                    if (before - after > tolerance)
                    {
                        Array.Reverse(p, i, j - i + 1);
                        improved = true;
                    }
                }
        }
        return p;
    }

    /// <summary>The length of the open path <paramref name="path"/>: the sum of its consecutive distances.</summary>
    public static double Length(double[,] distances, ReadOnlySpan<int> path)
    {
        double length = 0.0;
        for (int k = 1; k < path.Length; k++) length += distances[path[k - 1], path[k]];
        return length;
    }

    /// <summary><paramref name="start"/>, then the other nodes of 0 … <paramref name="count"/> − 1 in a uniformly random
    /// order.</summary>
    public static int[] Random(int count, int start, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if ((uint)start >= (uint)count) throw new ArgumentOutOfRangeException(nameof(start));
        var path = new int[count];
        path[0] = start;
        for (int k = 1, c = 0; k < count; c++)
            if (c != start) path[k++] = c;
        for (int k = count - 1; k > 1; k--)   // Fisher–Yates over positions 1 … count − 1
        {
            int j = 1 + rng.NextInt(k);
            (path[k], path[j]) = (path[j], path[k]);
        }
        return path;
    }

    private static int Check(double[,] distances)
    {
        ArgumentNullException.ThrowIfNull(distances);
        int s = distances.GetLength(0);
        if (s == 0 || distances.GetLength(1) != s) throw new ArgumentException("The distances must be a square matrix over at least one node.", nameof(distances));
        return s;
    }
}
