using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using Maths.Geometry;
using Maths.Geometry.Estimators.Intrinsic;
using Maths.Geometry.Solver;
using Maths.LinAlg;
using Xunit;

namespace Maths.Oracle.Tests;

public sealed class GrassmannMedianParityTests
{
    [Fact]
    public void GrassmannMedian_Matches_Riemann_IntrinsicMedian()
    {
        if (!ROracle.IsAvailable) return;

        const int d = 5, k = 2;
        double[][] frames = BuildClusteredFrames(d, k);
        string csv = WriteFramesCsv(frames);

        try
        {
            var grass = new GrassmannManifold(d, k);
            double[] weights = Enumerable.Repeat(1.0 / frames.Length, frames.Length).ToArray();
            double[] csMedian = InitialFrame(frames, d, k);
            GeometricMedian.Compute(
                grass,
                frames,
                weights,
                csMedian,
                IrlsOptions.Default with
                {
                    MaxIterations = 200,
                    Tolerance = 1e-8,
                    HybridMode = HybridMode.WeiszfeldOnly,
                });

            JsonElement r = ROracle.Run(
                "oracles/mom_oracle.R",
                csv,
                d.ToString(CultureInfo.InvariantCulture),
                k.ToString(CultureInfo.InvariantCulture),
                "200",
                "1e-5");
            double[] rMedian = MatrixRowsToColumnMajor(r.GetProperty("median"), d, k);

            double distance = grass.Distance(csMedian, rMedian);
            Assert.InRange(distance, 0.0, 2e-3);
        }
        finally
        {
            File.Delete(csv);
        }
    }

    private static double[][] BuildClusteredFrames(int d, int k)
    {
        double[][] offsets =
        [
            [ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00 ],
            [ 0.05, -0.02,  0.03, -0.04,  0.02,  0.01 ],
            [-0.04,  0.03, -0.02,  0.03, -0.01,  0.02 ],
            [ 0.03,  0.02, -0.04,  0.01,  0.04, -0.02 ],
            [-0.02, -0.04,  0.02,  0.04, -0.03,  0.03 ],
            [ 0.08, -0.05,  0.05, -0.08,  0.03,  0.04 ],
            [-0.07,  0.04, -0.06,  0.06, -0.04, -0.03 ],
        ];

        var frames = new double[offsets.Length][];
        for (int i = 0; i < offsets.Length; i++)
            frames[i] = FrameFromOffsets(d, k, offsets[i]);
        return frames;
    }

    private static double[] FrameFromOffsets(int d, int k, double[] offsets)
    {
        var frame = new double[d * k];
        frame[0] = 1.0;
        frame[d + 1] = 1.0;

        int idx = 0;
        for (int col = 0; col < k; col++)
            for (int row = k; row < d; row++)
                frame[col * d + row] = offsets[idx++];

        MatrixOps.Orthonormalize(frame, d, k);
        return frame;
    }

    private static double[] InitialFrame(double[][] frames, int d, int k)
    {
        var init = new double[d * k];
        foreach (double[] frame in frames)
            for (int i = 0; i < init.Length; i++)
                init[i] += frame[i] / frames.Length;
        MatrixOps.Orthonormalize(init, d, k);
        return init;
    }

    private static string WriteFramesCsv(double[][] frames)
    {
        string path = Path.Combine(Path.GetTempPath(), $"grassmann_{Guid.NewGuid():N}.csv");
        var sb = new StringBuilder();
        foreach (double[] frame in frames)
            sb.AppendLine(string.Join(",", frame.Select(v => v.ToString("R", CultureInfo.InvariantCulture))));
        File.WriteAllText(path, sb.ToString());
        return path;
    }

    private static double[] MatrixRowsToColumnMajor(JsonElement matrix, int rows, int cols)
    {
        double[][] byRow = matrix.EnumerateArray()
            .Select(row => row.EnumerateArray().Select(x => x.GetDouble()).ToArray())
            .ToArray();

        var packed = new double[rows * cols];
        for (int row = 0; row < rows; row++)
            for (int col = 0; col < cols; col++)
                packed[col * rows + row] = byRow[row][col];
        return packed;
    }
}
