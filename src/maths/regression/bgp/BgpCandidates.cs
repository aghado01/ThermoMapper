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
    /// through f, so each term leaves its unit out rather than reading it in sample (CB-4). Draws that repeat their predecessor's (t, η)
    /// reuse its row, so a run costs one O(n³) refit per accepted state.</para>
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
                    _gp.PredictiveInFrame(d.Fit, mapped, grids, rowValues, rowCdfs);
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
