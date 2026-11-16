using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace SplineSpaceProbe;

/// <summary>
/// The probe's own spline arithmetic, independent of the engine: clamped knot vectors, Cox–de Boor evaluation, the
/// supports of (4.133) and Theorem 4.63, the greedy ordered matching (S7d), one-sided Jacobi singular values, Jacobi
/// eigenvalues, a bipartite matching by augmenting paths, and a Cholesky least-squares solve.
/// </summary>
public static class SplineTools
{
    public static string F(double v, string f = "G6") => v.ToString(f, CultureInfo.InvariantCulture);

    /// <summary>The clamped knot vector of order m on [0, 1]: m zeros, the interior knots, m ones.</summary>
    public static double[] Clamped(double[] interior, int order)
    {
        var t = new double[interior.Length + 2 * order];
        for (int i = 0; i < order; i++) { t[i] = 0.0; t[t.Length - 1 - i] = 1.0; }
        Array.Copy(interior, 0, t, order, interior.Length);
        return t;
    }

    /// <summary>All n = len(t) − m B-splines of order m at x, Cox–de Boor, right-closed at the last knot.</summary>
    public static double[] Basis(double[] t, int order, double x)
    {
        int n = t.Length - order;
        // Order 1: the indicator of [t_j, t_{j+1}), with the last nonempty interval closed at its right end.
        var b = new double[t.Length - 1];
        int last = -1;
        for (int j = 0; j < t.Length - 1; j++) if (t[j] < t[j + 1]) last = j;
        for (int j = 0; j < t.Length - 1; j++)
            b[j] = t[j] < t[j + 1] && ((x >= t[j] && x < t[j + 1]) || (j == last && x == t[j + 1])) ? 1.0 : 0.0;
        for (int r = 2; r <= order; r++)
        {
            var next = new double[t.Length - r];
            for (int j = 0; j < next.Length; j++)
            {
                double v = 0.0;
                if (t[j + r - 1] > t[j]) v += (x - t[j]) / (t[j + r - 1] - t[j]) * b[j];
                if (t[j + r] > t[j + 1]) v += (t[j + r] - x) / (t[j + r] - t[j + 1]) * b[j + 1];
                next[j] = v;
            }
            b = next;
        }
        return b[..n];
    }

    public static double[,] Design(double[] t, int order, double[] x)
    {
        int n = t.Length - order;
        var z = new double[x.Length, n];
        for (int i = 0; i < x.Length; i++)
        {
            double[] row = Basis(t, order, x[i]);
            for (int j = 0; j < n; j++) z[i, j] = row[j];
        }
        return z;
    }

    /// <summary>Support σ_j of B_j (4.133, Theorem 4.63): (t_j, t_{j+m}), closed on the left when t_j = t_{j+m−1},
    /// and closed at 1 for the last B-spline.</summary>
    public readonly record struct Support(double Lo, double Hi, bool LeftClosed, bool RightClosed)
    {
        public bool Contains(double x) => (LeftClosed ? x >= Lo : x > Lo) && (RightClosed ? x <= Hi : x < Hi);
    }

    public static Support[] Supports(double[] t, int order)
    {
        int n = t.Length - order;
        var s = new Support[n];
        for (int j = 0; j < n; j++)
            s[j] = new Support(t[j], t[j + order], t[j] == t[j + order - 1], j == n - 1);
        return s;
    }

    /// <summary>S7d's greedy ordered matching: for each support in turn, the first unused site after the previous pick
    /// that lies in it. True when every support is served.</summary>
    public static bool Matches(double[] sortedSites, Support[] supports)
    {
        int next = 0;
        foreach (Support s in supports)
        {
            while (next < sortedSites.Length && !s.Contains(sortedSites[next])) next++;
            if (next == sortedSites.Length) return false;
            next++;
        }
        return true;
    }

    /// <summary>The singular values of an N × n matrix by one-sided Jacobi on its columns, descending.</summary>
    public static double[] SingularValues(double[,] a)
    {
        int rows = a.GetLength(0), cols = a.GetLength(1);
        var u = (double[,])a.Clone();
        for (int sweep = 0; sweep < 80; sweep++)
        {
            double off = 0.0;
            for (int p = 0; p < cols - 1; p++)
                for (int q = p + 1; q < cols; q++)
                {
                    double alpha = 0, beta = 0, gamma = 0;
                    for (int i = 0; i < rows; i++)
                    {
                        alpha += u[i, p] * u[i, p];
                        beta += u[i, q] * u[i, q];
                        gamma += u[i, p] * u[i, q];
                    }
                    if (gamma == 0.0 || Math.Abs(gamma) <= 1e-17 * Math.Sqrt(alpha * beta)) continue;
                    off = Math.Max(off, Math.Abs(gamma) / Math.Sqrt(alpha * beta));
                    double zeta = (beta - alpha) / (2 * gamma);
                    double tt = Math.Sign(zeta) / (Math.Abs(zeta) + Math.Sqrt(1 + zeta * zeta));
                    if (zeta == 0) tt = 1.0;
                    double c = 1 / Math.Sqrt(1 + tt * tt), s = c * tt;
                    for (int i = 0; i < rows; i++)
                    {
                        double up = u[i, p], uq = u[i, q];
                        u[i, p] = c * up - s * uq;
                        u[i, q] = s * up + c * uq;
                    }
                }
            if (off < 1e-15) break;
        }
        var sv = new double[cols];
        for (int j = 0; j < cols; j++)
        {
            double sum = 0;
            for (int i = 0; i < rows; i++) sum += u[i, j] * u[i, j];
            sv[j] = Math.Sqrt(sum);
        }
        if (rows < cols)   // only `rows` of them can be nonzero
        {
            Array.Sort(sv);
            Array.Reverse(sv);
            for (int j = rows; j < cols; j++) sv[j] = 0.0;
        }
        return sv.OrderByDescending(v => v).ToArray();
    }

    /// <summary>The eigenvalues of a symmetric matrix by cyclic Jacobi, ascending.</summary>
    public static double[] Eigenvalues(double[,] s)
    {
        int n = s.GetLength(0);
        var a = (double[,])s.Clone();
        for (int sweep = 0; sweep < 100; sweep++)
        {
            double off = 0, diag = 0;
            for (int i = 0; i < n; i++)
            {
                diag += a[i, i] * a[i, i];
                for (int j = i + 1; j < n; j++) off += a[i, j] * a[i, j];
            }
            if (off <= 1e-32 * diag) break;
            for (int p = 0; p < n - 1; p++)
                for (int q = p + 1; q < n; q++)
                {
                    if (a[p, q] == 0) continue;
                    double theta = (a[q, q] - a[p, p]) / (2 * a[p, q]);
                    double t = Math.Sign(theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1));
                    if (theta == 0) t = 1;
                    double c = 1 / Math.Sqrt(t * t + 1), sn = t * c;
                    for (int k = 0; k < n; k++)
                    {
                        double akp = a[k, p], akq = a[k, q];
                        a[k, p] = c * akp - sn * akq;
                        a[k, q] = sn * akp + c * akq;
                    }
                    for (int k = 0; k < n; k++)
                    {
                        double apk = a[p, k], aqk = a[q, k];
                        a[p, k] = c * apk - sn * aqk;
                        a[q, k] = sn * apk + c * aqk;
                    }
                }
        }
        var e = new double[n];
        for (int i = 0; i < n; i++) e[i] = a[i, i];
        Array.Sort(e);
        return e;
    }

    /// <summary>The 2-norm condition number of the equilibrated Gram D^{-1/2} ZᵀZ D^{-1/2}.</summary>
    public static double EquilibratedCondition(double[,] z)
    {
        int rows = z.GetLength(0), cols = z.GetLength(1);
        var g = new double[cols, cols];
        for (int i = 0; i < rows; i++)
            for (int p = 0; p < cols; p++)
            {
                double zp = z[i, p];
                if (zp == 0) continue;
                for (int q = p; q < cols; q++) g[p, q] += zp * z[i, q];
            }
        for (int p = 0; p < cols; p++) for (int q = 0; q < p; q++) g[p, q] = g[q, p];
        var d = new double[cols];
        for (int p = 0; p < cols; p++) d[p] = g[p, p] > 0 ? 1 / Math.Sqrt(g[p, p]) : 0;
        for (int p = 0; p < cols; p++) for (int q = 0; q < cols; q++) g[p, q] *= d[p] * d[q];
        double[] e = Eigenvalues(g);
        return e[0] <= 0 ? double.PositiveInfinity : e[^1] / e[0];
    }

    /// <summary>A maximum bipartite matching of columns to rows (Kuhn's augmenting paths); returns its size.</summary>
    public static int MaximumMatching(bool[,] adjacent)
    {
        int rows = adjacent.GetLength(0), cols = adjacent.GetLength(1);
        var rowOf = Enumerable.Repeat(-1, cols).ToArray();
        var colOf = Enumerable.Repeat(-1, rows).ToArray();
        int size = 0;
        for (int c = 0; c < cols; c++)
        {
            var seen = new bool[rows];
            if (Augment(c)) size++;

            bool Augment(int col)
            {
                for (int r = 0; r < rows; r++)
                {
                    if (!adjacent[r, col] || seen[r]) continue;
                    seen[r] = true;
                    if (colOf[r] < 0 || Augment(colOf[r]))
                    {
                        colOf[r] = col;
                        rowOf[col] = r;
                        return true;
                    }
                }
                return false;
            }
        }
        return size;
    }

    /// <summary>Least squares by Cholesky on the normal equations; returns the coefficients.</summary>
    public static double[] LeastSquares(double[,] z, double[] y)
    {
        int rows = z.GetLength(0), cols = z.GetLength(1);
        var g = new double[cols, cols];
        var r = new double[cols];
        for (int i = 0; i < rows; i++)
            for (int p = 0; p < cols; p++)
            {
                double zp = z[i, p];
                if (zp == 0) continue;
                r[p] += zp * y[i];
                for (int q = p; q < cols; q++) g[p, q] += zp * z[i, q];
            }
        for (int p = 0; p < cols; p++) for (int q = 0; q < p; q++) g[p, q] = g[q, p];
        var l = new double[cols, cols];
        for (int i = 0; i < cols; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = g[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                if (i == j) l[i, i] = Math.Sqrt(Math.Max(s, 1e-300));
                else l[i, j] = s / l[j, j];
            }
        var w = new double[cols];
        for (int i = 0; i < cols; i++)
        {
            double s = r[i];
            for (int k = 0; k < i; k++) s -= l[i, k] * w[k];
            w[i] = s / l[i, i];
        }
        var b = new double[cols];
        for (int i = cols - 1; i >= 0; i--)
        {
            double s = w[i];
            for (int k = i + 1; k < cols; k++) s -= l[k, i] * b[k];
            b[i] = s / l[i, i];
        }
        return b;
    }

    public static double[] Multiply(double[,] z, double[] b)
    {
        int rows = z.GetLength(0), cols = z.GetLength(1);
        var y = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            double s = 0;
            for (int j = 0; j < cols; j++) s += z[i, j] * b[j];
            y[i] = s;
        }
        return y;
    }

    public static double Rms(double[] a, double[] b) => Math.Sqrt(a.Zip(b, (u, v) => (u - v) * (u - v)).Average());

    public static double Slope(IReadOnlyList<double> x, IReadOnlyList<double> y)
    {
        double mx = x.Average(), my = y.Average(), sxy = 0, sxx = 0;
        for (int i = 0; i < x.Count; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) * (x[i] - mx); }
        return sxy / sxx;
    }
}
