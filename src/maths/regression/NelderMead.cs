using System;

namespace Maths.Regression;

/// <summary>
/// Nelder–Mead simplex minimization in a box, with the standard coefficients (reflection 1, expansion 2, contraction ½,
/// shrink ½; Lagarias, Reeds, Wright &amp; Wright 1998) and every trial point clamped into the box. Derivative-free, for the
/// regression engines' small smooth searches, such as the anisotropic REML's two log smoothings. It stops when the
/// simplex's values agree within <c>tolerance·(1 + |f_best|)</c> and every vertex lies within <c>tolerance</c> of the best
/// in every coordinate, or when the evaluation budget runs out. +∞ (or NaN) marks an infeasible point; the start must be
/// feasible.
/// </summary>
internal static class NelderMead
{
    /// <summary>The best vertex, its value, the evaluations spent, and whether the convergence test passed.</summary>
    internal readonly record struct Result(double[] Point, double Value, int Evaluations, bool Converged);

    public static Result Minimize(Func<double[], double> f, double[] start, double[] step, double[] lower, double[] upper,
                                  double tolerance, int maxEvaluations)
    {
        ArgumentNullException.ThrowIfNull(f);
        ArgumentNullException.ThrowIfNull(start);
        int d = start.Length;
        if (d == 0 || step.Length != d || lower.Length != d || upper.Length != d)
            throw new ArgumentException("The start, steps and bounds need one entry per coordinate.", nameof(start));
        for (int k = 0; k < d; k++)
            if (!(lower[k] <= start[k] && start[k] <= upper[k]) || !(step[k] > 0.0))
                throw new ArgumentException("The start must lie in the box, and every step must be positive.", nameof(start));
        if (!(tolerance > 0.0)) throw new ArgumentOutOfRangeException(nameof(tolerance));
        if (maxEvaluations < d + 1) throw new ArgumentOutOfRangeException(nameof(maxEvaluations));

        int evaluations = 0;
        double Evaluate(double[] p)
        {
            evaluations++;
            double v = f(p);
            return double.IsNaN(v) ? double.PositiveInfinity : v;
        }
        void Clamp(double[] p)
        {
            for (int k = 0; k < d; k++) p[k] = Math.Clamp(p[k], lower[k], upper[k]);
        }

        var x = new double[d + 1][];
        var fx = new double[d + 1];
        x[0] = (double[])start.Clone();
        fx[0] = Evaluate(x[0]);
        if (double.IsPositiveInfinity(fx[0])) throw new ArgumentException("The start must be feasible.", nameof(start));
        for (int i = 0; i < d; i++)
        {
            var v = (double[])start.Clone();
            v[i] = start[i] + step[i] <= upper[i] ? start[i] + step[i] : start[i] - step[i];
            Clamp(v);
            x[i + 1] = v;
            fx[i + 1] = Evaluate(v);
        }

        var centroid = new double[d];
        while (true)
        {
            Sort(x, fx);
            if (Converged(x, fx, tolerance)) return new Result((double[])x[0].Clone(), fx[0], evaluations, true);
            if (evaluations >= maxEvaluations) return new Result((double[])x[0].Clone(), fx[0], evaluations, false);

            Array.Clear(centroid);
            for (int i = 0; i < d; i++)
                for (int k = 0; k < d; k++) centroid[k] += x[i][k] / d;
            double[] worst = x[d];

            double[] reflected = Along(centroid, worst, -1.0, Clamp);
            double fr = Evaluate(reflected);
            if (fr < fx[0])
            {
                double[] expanded = Along(centroid, worst, -2.0, Clamp);
                double fe = Evaluate(expanded);
                if (fe < fr) { x[d] = expanded; fx[d] = fe; }
                else { x[d] = reflected; fx[d] = fr; }
                continue;
            }
            if (fr < fx[d - 1]) { x[d] = reflected; fx[d] = fr; continue; }

            bool outside = fr < fx[d];
            double[] contracted = outside ? Along(centroid, reflected, 0.5, Clamp) : Along(centroid, worst, 0.5, Clamp);
            double fc = Evaluate(contracted);
            if (outside ? fc <= fr : fc < fx[d]) { x[d] = contracted; fx[d] = fc; continue; }

            for (int i = 1; i <= d; i++)   // shrink towards the best vertex
            {
                x[i] = Along(x[0], x[i], 0.5, Clamp);
                fx[i] = Evaluate(x[i]);
            }
        }
    }

    // c + t·(p − c), clamped into the box: t = −1 reflects p through c, −2 expands, ½ contracts.
    private static double[] Along(double[] c, double[] p, double t, Action<double[]> clamp)
    {
        var q = new double[c.Length];
        for (int k = 0; k < c.Length; k++) q[k] = c[k] + t * (p[k] - c[k]);
        clamp(q);
        return q;
    }

    private static void Sort(double[][] x, double[] fx)
    {
        for (int i = 1; i < fx.Length; i++)
        {
            double f = fx[i];
            double[] p = x[i];
            int j = i - 1;
            while (j >= 0 && fx[j] > f) { fx[j + 1] = fx[j]; x[j + 1] = x[j]; j--; }
            fx[j + 1] = f;
            x[j + 1] = p;
        }
    }

    private static bool Converged(double[][] x, double[] fx, double tolerance)
    {
        double best = fx[0];
        for (int i = 1; i < fx.Length; i++)
        {
            if (!(Math.Abs(fx[i] - best) <= tolerance * (1.0 + Math.Abs(best)))) return false;
            for (int k = 0; k < x[0].Length; k++)
                if (!(Math.Abs(x[i][k] - x[0][k]) <= tolerance)) return false;
        }
        return true;
    }
}
