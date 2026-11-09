using System;
using System.Collections.Generic;

namespace Maths.Regression.Bgp;

public sealed partial class BgpSampler
{
    /// <summary>
    /// BGP's candidate density (E16's seam, the conformal layer's slice 5): per retained draw (t, η), the predictive of a new observation at
    /// each case given the training data, log p(y | x_t, Y, t, η), on each case's grid (<see cref="GpRegression.PredictiveLogDensities"/>),
    /// in <see cref="PointwiseLogDensity"/>'s rows with its chain lengths and latent report. Its content is the collapsed one, BGP's only:
    /// f, the constant mean and s² are integrated against their posterior given the training data, never drawn.
    /// <para>The split form takes it (<see cref="Conformal.SplitConformal"/>, CB-2): a fit on the training split read at the calibration
    /// cases and the test cases. The full step refuses it, as it refuses the training density: a draw (t, η) leaves the units dependent
    /// through f, so each term leaves its unit out rather than reading it in sample (CB-4).</para>
    /// <para>At the run's own test points (<see cref="BgpResult.TestInputs"/>), each draw's fit's variance there was kept
    /// (<see cref="BgpDraw.FitVariance"/>), so the seam reads it and the mean at O(n) per case, the same numbers bit for bit. Elsewhere
    /// draws that repeat their predecessor's (t, η) reuse its row, and each accepted state costs one O(n³) refit.</para>
    /// </summary>
    /// <param name="result">A run of this sampler, whose retained draws are read.</param>
    /// <param name="xCases">The new cases' inputs (cases × D), in the original units.</param>
    /// <param name="grids">Each case's candidate values, ascending.</param>
    /// <param name="strata">Each case's stratum in the calibration's partition; null puts every case in stratum 0.</param>
    /// <param name="cdf">Adds each draw's log CDF at each candidate, for the conformal predictive distribution.</param>
    public CandidateLogDensity CandidateLogDensity(BgpResult result, double[,] xCases, IReadOnlyList<double[]> grids, int[]? strata = null,
                                                   bool cdf = false)
    {
        ArgumentNullException.ThrowIfNull(result);
        ArgumentNullException.ThrowIfNull(xCases);
        if (xCases.GetLength(1) != _gp.Dimension) throw new ArgumentException("The cases' inputs must have the training dimension.", nameof(xCases));
        int[] offsets = CandidateSeam.Offsets(grids, xCases.GetLength(0));
        double[,] mapped = _gp.Transform.MapInputs(xCases);
        int[]? kept = KeptRows(result, xCases);
        var variance = kept is null ? null : new double[kept.Length];
        int columns = offsets[^1];
        var values = new double[result.Draws, columns];
        double[,]? logCdfs = cdf ? new double[result.Draws, columns] : null;
        var rowValues = new double[columns];
        double[]? rowCdfs = cdf ? new double[columns] : null;
        GpFit? previous = null;
        int row = 0;
        foreach (IReadOnlyList<BgpDraw> chain in result.BandwidthDraws)
            foreach (BgpDraw d in chain)
            {
                if (previous is null || d.Fit.Bandwidth != previous.Bandwidth || d.Fit.NoiseRatio != previous.NoiseRatio)
                {
                    if (kept is null) _gp.PredictiveInFrame(d.Fit, mapped, grids, rowValues, rowCdfs);
                    else
                    {
                        for (int t = 0; t < kept.Length; t++) variance![t] = d.FitVariance![kept[t]];
                        _gp.Predictive(d.Fit, mapped, variance!, grids, rowValues, rowCdfs);
                    }
                    previous = d.Fit;
                }
                for (int c = 0; c < columns; c++)
                {
                    values[row, c] = rowValues[c];
                    if (rowCdfs is not null) logCdfs![row, c] = rowCdfs[c];
                }
                row++;
            }
        return new CandidateLogDensity(values, grids, Latents(), strata, ChainLengths.Of(result.BandwidthDraws)) { LogCdf = logCdfs };
    }

    // Each case's row among the run's test points, where every draw kept its fit's variance; null if a case is not among them or a
    // draw kept none. Rows match on their values exactly.
    private static int[]? KeptRows(BgpResult result, double[,] xCases)
    {
        if (result.TestInputs is not { } inputs || inputs.GetLength(1) != xCases.GetLength(1)) return null;
        foreach (IReadOnlyList<BgpDraw> chain in result.BandwidthDraws)
            foreach (BgpDraw d in chain)
                if (d.FitVariance is null || d.FitVariance.Length != inputs.GetLength(0)) return null;
        var rows = new Dictionary<Row, int>();
        for (int i = inputs.GetLength(0) - 1; i >= 0; i--) rows[new Row(inputs, i)] = i;
        var kept = new int[xCases.GetLength(0)];
        for (int t = 0; t < kept.Length; t++)
            if (!rows.TryGetValue(new Row(xCases, t), out kept[t])) return null;
        return kept;
    }

    // A matrix row as a key, its values compared exactly (0 and −0 alike).
    private readonly struct Row(double[,] matrix, int index) : IEquatable<Row>
    {
        public bool Equals(Row other)
        {
            for (int c = 0; c < matrix.GetLength(1); c++)
                if (!(matrix[index, c] == other.Matrix[other.Index, c])) return false;
            return true;
        }

        public override bool Equals(object? obj) => obj is Row other && Equals(other);

        public override int GetHashCode()
        {
            var hash = new HashCode();
            for (int c = 0; c < matrix.GetLength(1); c++) hash.Add(matrix[index, c] == 0.0 ? 0.0 : matrix[index, c]);
            return hash.ToHashCode();
        }

        private double[,] Matrix => matrix;
        private int Index => index;
    }

    // The collapsed content's latents: the bandwidth and the noise ratio drawn and conditioned on; f, the mean and s² integrated against
    // their posterior given the observations outside the unit.
    private List<LatentReport> Latents()
    {
        var latents = new List<LatentReport> { new("bandwidth t", LatentScope.Global, LatentTreatment.Conditioned) };
        if (_unknownScale) latents.Add(new("noise ratio η", LatentScope.Global, LatentTreatment.Conditioned));
        latents.Add(new("regression function f", LatentScope.Global, LatentTreatment.IntegratedGivenOthers));
        latents.Add(new("constant mean m", LatentScope.Global, LatentTreatment.IntegratedGivenOthers));
        if (_unknownScale) latents.Add(new("signal variance s²", LatentScope.Global, LatentTreatment.IntegratedGivenOthers));
        return latents;
    }
}
