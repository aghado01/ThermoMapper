using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using TDA.Ph;
using Xunit;

namespace Maths.Oracle.Tests;

public sealed class TdaParityTests
{
    /// <summary>
    /// Our full-Rips PH — a complete distance graph → <see cref="RipsFiltration"/> →
    /// <see cref="PersistentHomology"/> — must reproduce the persistence diagram Ripser computes
    /// (via <c>TDAstats::calculate_homology</c>, the exact engine <c>TDAkit::diagRips</c> wraps).
    /// Finite H0/H1 bars match within tolerance; C# additionally carries the one essential H0 bar
    /// Ripser omits. Cross-checks the boundary-matrix reduction against the gold-standard library.
    /// </summary>
    [Fact]
    public void FullRips_PH_Matches_Ripser_H0_H1()
    {
        if (!ROracle.IsAvailable) return;   // opt-in: no-op without the portable R toolchain

        // Small n on purpose: a complete graph has C(n,3) triangles, and the dense full-Rips reduction
        // is superlinear — 12 points (β0=1, β1=1) is plenty to exercise H0+H1 and stays fast.
        const int n = 12;
        double[][] pts = NoisyCircle(n, seed: 7);
        string csv = WriteCsv(pts);
        try
        {
            // C#: full Rips (maxDim = 2 builds the triangle fillers) → PH.
            Barcode cs = PersistentHomology.Compute(FullRips.Build(pts, maxDimension: 2), 2);
            List<(double b, double d)> csH0 = FiniteBars(cs, 0);
            List<(double b, double d)> csH1 = FiniteBars(cs, 1);
            int csEssentialH0 = cs.Bars.Count(b => b.Dimension == 0 && b.IsInfinite);

            // R: Ripser full Rips (threshold past the diameter so nothing is cut).
            double thr = Diameter(pts) * 1.5;
            JsonElement r = ROracle.Run("oracles/tda_oracle.R", csv, "1",
                thr.ToString("R", CultureInfo.InvariantCulture));
            (List<(double b, double d)> rH0, List<(double b, double d)> rH1) = ParseFiniteBars(r);

            Assert.Equal(1, csEssentialH0);            // the one component that lives forever
            AssertBarsMatch(csH0, rH0, tol: 1e-4, dim: 0);
            AssertBarsMatch(csH1, rH1, tol: 1e-4, dim: 1);
        }
        finally { File.Delete(csv); }
    }

    private static List<(double b, double d)> FiniteBars(Barcode bc, int dim)
    {
        // Ripser omits diagonal zero-persistence intervals; the explicit reducer materializes them.
        const double diagonalTol = 1e-10;
        var list = new List<(double, double)>();
        foreach (Bar bar in bc.Bars)
            if (bar.Dimension == dim && !bar.IsInfinite && bar.Death - bar.Birth > diagonalTol)
                list.Add((bar.Birth, bar.Death));
        return list;
    }

    private static (List<(double b, double d)> h0, List<(double b, double d)> h1) ParseFiniteBars(JsonElement r)
    {
        int[] dim = r.GetProperty("dimension").EnumerateArray().Select(x => (int)x.GetDouble()).ToArray();
        double[] birth = r.GetProperty("birth").EnumerateArray().Select(x => x.GetDouble()).ToArray();
        double[] death = r.GetProperty("death").EnumerateArray().Select(x => x.GetDouble()).ToArray();
        var h0 = new List<(double, double)>();
        var h1 = new List<(double, double)>();
        for (int i = 0; i < dim.Length; i++)
        {
            if (death[i] < 0) continue;                 // essential (sentinel −1)
            if (dim[i] == 0) h0.Add((birth[i], death[i]));
            else if (dim[i] == 1) h1.Add((birth[i], death[i]));
        }
        return (h0, h1);
    }

    private static void AssertBarsMatch(List<(double b, double d)> a, List<(double b, double d)> bR, double tol, int dim)
    {
        a.Sort(CmpBar);
        bR.Sort(CmpBar);
        Assert.True(a.Count == bR.Count, $"H{dim} finite-bar count: C#={a.Count} R={bR.Count}");
        for (int i = 0; i < a.Count; i++)
        {
            Assert.True(Math.Abs(a[i].b - bR[i].b) < tol, $"H{dim}[{i}] birth: C#={a[i].b} R={bR[i].b}");
            Assert.True(Math.Abs(a[i].d - bR[i].d) < tol, $"H{dim}[{i}] death: C#={a[i].d} R={bR[i].d}");
        }
    }

    private static int CmpBar((double b, double d) x, (double b, double d) y)
    {
        int c = x.b.CompareTo(y.b);
        return c != 0 ? c : x.d.CompareTo(y.d);
    }

    private static double[][] NoisyCircle(int n, int seed)
    {
        var rng = new Random(seed);
        var pts = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double t = 2.0 * Math.PI * i / n;
            pts[i] = new[]
            {
                Math.Cos(t) + 0.02 * (rng.NextDouble() - 0.5),
                Math.Sin(t) + 0.02 * (rng.NextDouble() - 0.5),
            };
        }
        return pts;
    }

    private static double Dist(double[] a, double[] b)
    {
        double s = 0.0;
        for (int k = 0; k < a.Length; k++) { double d = a[k] - b[k]; s += d * d; }
        return Math.Sqrt(s);
    }

    private static double Diameter(double[][] p)
    {
        double m = 0.0;
        for (int i = 0; i < p.Length; i++)
            for (int j = i + 1; j < p.Length; j++)
            {
                double d = Dist(p[i], p[j]);
                if (d > m) m = d;
            }
        return m;
    }

    private static string WriteCsv(double[][] data)
    {
        string path = Path.Combine(Path.GetTempPath(), $"cloud_{Guid.NewGuid():N}.csv");
        var sb = new StringBuilder();
        foreach (double[] row in data)
            sb.AppendLine(string.Join(",", row.Select(v => v.ToString("R", CultureInfo.InvariantCulture))));
        File.WriteAllText(path, sb.ToString());
        return path;
    }
}
