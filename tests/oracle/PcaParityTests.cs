using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using Maths.Geometry;
using Maths.LinAlg;
using Xunit;

namespace Maths.Oracle.Tests;

public sealed class PcaParityTests
{
    /// <summary>
    /// The same fixture through Maths.LinAlg.Pca and base R's prcomp must agree: identical mean,
    /// the same principal subspace (sign-/rotation-agnostic, measured by Grassmann distance), and
    /// matching eigenvalues once R's (n−1) denominator is rescaled to the C# MLE (n).
    /// </summary>
    [Fact]
    public void Pca_Subspace_And_Eigenvalues_Match_BaseR_prcomp()
    {
        // Opt-in oracle: when the portable R toolchain isn't present this parity check is a no-op
        // (it never reds the suite); it runs for real wherever $PORTABLE_ROOT/rlang exists.
        if (!ROracle.IsAvailable) return;

        const int n = 40, d = 4, k = 2;
        double[][] data = BuildData(n, d, seed: 17);
        string csv = WriteCsv(data);

        try
        {
            PcaResult cs = Pca.Compute(data, numComponents: k, center: true, whiten: false);
            JsonElement r = ROracle.Run("oracles/pca_oracle.R", csv, k.ToString(CultureInfo.InvariantCulture));

            // Mean — direct.
            double[] rMean = ToArray(r.GetProperty("mean"));
            for (int j = 0; j < d; j++) Assert.Equal(cs.Mean[j], rMean[j], 9);

            // Subspace — sign-agnostic + rotation-invariant via Grassmann distance.
            double[][] rComps = ToMatrix(r.GetProperty("components"));   // k x d
            var grass = new GrassmannManifold(d, k);
            double dGr = grass.Distance(PackFrame(cs.Components, d, k), PackFrame(rComps, d, k));
            Assert.InRange(dGr, 0.0, 1e-6);

            // Eigenvalues — rescale R's (n−1) denominator to the C# MLE (n).
            double[] rEig = ToArray(r.GetProperty("eigenvalues"));
            for (int c = 0; c < k; c++)
                Assert.Equal(cs.Eigenvalues[c], rEig[c] * (n - 1.0) / n, 6);
        }
        finally
        {
            File.Delete(csv);
        }
    }

    private static double[][] BuildData(int n, int d, int seed)
    {
        var rng = new Random(seed);
        var data = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double z1 = NextGaussian(rng), z2 = NextGaussian(rng), z3 = NextGaussian(rng), z4 = NextGaussian(rng);
            // Correlated with a well-separated spectrum → the top-2 subspace is unambiguous
            // and non-axis-aligned (so the eigensolve is genuinely exercised, not just axes).
            data[i] = new[]
            {
                5.0 * z1 + 1.0 * z2,
                3.0 * z2 - 0.5 * z1,
                1.5 * z3,
                0.5 * z4,
            };
        }
        return data;
    }

    private static double NextGaussian(Random rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static string WriteCsv(double[][] data)
    {
        string path = Path.Combine(Path.GetTempPath(), $"fixture_{Guid.NewGuid():N}.csv");
        var sb = new StringBuilder();
        foreach (double[] row in data)
            sb.AppendLine(string.Join(",", row.Select(v => v.ToString("R", CultureInfo.InvariantCulture))));
        File.WriteAllText(path, sb.ToString());
        return path;
    }

    private static double[] ToArray(JsonElement e) => e.EnumerateArray().Select(x => x.GetDouble()).ToArray();
    private static double[][] ToMatrix(JsonElement e) => e.EnumerateArray().Select(ToArray).ToArray();

    // PCA components (k vectors of length d) → column-major d×k frame for GrassmannManifold.
    private static double[] PackFrame(double[][] comps, int d, int k)
    {
        var frame = new double[d * k];
        for (int c = 0; c < k; c++)
            for (int row = 0; row < d; row++)
                frame[c * d + row] = comps[c][row];
        return frame;
    }
}
