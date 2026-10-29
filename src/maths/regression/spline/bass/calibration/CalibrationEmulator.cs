using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// The emulator-draw contract the cut posterior reads (D09, C-X2): retained draws from p(e | D_sim), each one's
/// prediction at θ and its noise covariance, all in physical units. Calibration and stage 8's cut SMC are written against
/// it, not against BASS's term arrays, so a non-BASS emulator (BGP, ABMRS) can serve them. The observations do not update
/// it: that is the cut.
/// <list type="bullet">
/// <item>θ arrives in physical coordinates: the continuous parameters' values in declaration order among the continuous
/// ones, and each categorical parameter's level as an index into its declared levels.</item>
/// <item>The outputs are the observations' order; the emulator owns any fixed inputs the observations carry.</item>
/// <item>Implementations are safe to call from several chains at once.</item>
/// </list>
/// </summary>
public interface ICalibrationEmulator
{
    /// <summary>The calibration parameters the emulator's inputs take, in declaration order.</summary>
    IReadOnlyList<CalibrationParameter> Parameters { get; }

    /// <summary>The retained draws S.</summary>
    int DrawCount { get; }

    /// <summary>The draws' weights, one per draw; null when they weigh equally (MCMC's).</summary>
    IReadOnlyList<double>? DrawWeights { get; }

    /// <summary>The outputs m, one per observation.</summary>
    int Outputs { get; }

    /// <summary>The draw's noise covariance in physical units (<see cref="EmulatorNoise"/>).</summary>
    EmulatorNoise Noise(int draw);

    /// <summary>
    /// The truncation residual vectors the prediction adds (D07), in physical units, n × m, one case per row: y_s ∘ eᵢ; a
    /// fresh copy. Null for an emulator without a basis.
    /// </summary>
    double[,]? TruncationResiduals();

    /// <summary>Draw <paramref name="draw"/>'s latent prediction at θ (<paramref name="continuous"/>,
    /// <paramref name="levels"/>), in physical units, into <paramref name="output"/> (m values).</summary>
    void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output);

    /// <summary>
    /// Draw <paramref name="draw"/>'s prediction at the continuous values with each of <paramref name="combinations"/>'
    /// levels, into the rows of <paramref name="output"/> (one combination per row, m columns): the categorical
    /// conditional's enumeration. The default predicts each combination in turn; an emulator that can evaluate the
    /// continuous part once (Francom et al. 2019) overrides it.
    /// </summary>
    void PredictLevels(int draw, ReadOnlySpan<double> continuous, IReadOnlyList<int[]> combinations, double[,] output)
    {
        ArgumentNullException.ThrowIfNull(combinations);
        ArgumentNullException.ThrowIfNull(output);
        int m = Outputs;
        if (output.GetLength(0) != combinations.Count || output.GetLength(1) != m)
            throw new ArgumentException($"The output needs one row per combination and {m} columns.", nameof(output));
        var row = new double[m];
        for (int c = 0; c < combinations.Count; c++)
        {
            Predict(draw, continuous, combinations[c], row);
            for (int j = 0; j < m; j++) output[c, j] = row[j];
        }
    }
}

/// <summary>An emulator draw's noise covariance over the outputs, in physical units (D09).</summary>
public abstract record EmulatorNoise
{
    private protected EmulatorNoise() { }
}

/// <summary>
/// A basis emulator's component noise, (y_s ∘ K)·diag(σ²_{e,k})·(y_s ∘ K)ᵀ: a low-rank term over the basis in physical
/// units. <paramref name="Factor"/> (m × r) is shared by every draw and read only.
/// </summary>
public sealed record ComponentNoise(double[,] Factor, double[] Variances) : EmulatorNoise;

/// <summary>An emulator without a basis: its noise variance σ²_e on every output, σ²_e·I (D09's amendment; a rank-1 term
/// would assert one error shared by every output).</summary>
public sealed record IndependentNoise(double Variance) : EmulatorNoise;
