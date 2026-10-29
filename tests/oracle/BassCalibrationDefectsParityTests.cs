using System;
using System.Linq;
using System.Text.Json;
using Maths.Regression.Spline.Bass.Calibration;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// R's calibration defects (C-O1, <c>bass_calibration_defects_oracle.R</c>) against D09's corrections (stage 7). The
/// oracle runs <c>calibrate.bassBasis</c> as written for one iteration on a mock emulator whose two components are
/// linear in θ; each test asserts R's defect as the fixture records it and the corrected behaviour on this side, on the
/// same emulator rebuilt from the fixture's inputs: the branch chosen by shape and the units of the covariance. Opt-in:
/// it runs where the R toolchain and the BASS checkout are present.
/// </summary>
public sealed class BassCalibrationDefectsParityTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationDefectsParityTests(ITestOutputHelper output) => _out = output;

    private static bool Available => BassTargetParityTests.Available;

    internal static readonly Lazy<JsonElement> Oracle = new(() =>
        BassOutputTruncationParityTests.Run("oracles/bass_calibration_defects_oracle.R", BassCalibrationCovarianceParityTests.Fixture()));

    private static double[] Vector(JsonElement e) => BassCalibrationCovarianceParityTests.Vector(e);
    private static double[][] Rows(JsonElement e) => BassCalibrationCovarianceParityTests.Rows(e);

    // The fixture's mock emulator in physical coordinates: w_k(θ) = c_k0 + c_k1 θ₁ + c_k2 θ₂, its 50 draws of component
    // variances, and the residuals y_s ∘ e, the trunc.error columns.
    internal static OracleEmulator Emulator(JsonElement inputs)
    {
        (double[] offset, double[] scale, double[,] basis) = BassCalibrationCovarianceParityTests.Pieces();
        double[][] coef = Rows(inputs.GetProperty("coef"));
        double[][] trunc = Rows(inputs.GetProperty("trunc_error"));   // m × n
        int m = offset.Length, n = trunc[0].Length;
        var residuals = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++) residuals[i, l] = scale[l] * trunc[l][i];
        double[][] bounds = Rows(inputs.GetProperty("bounds"));   // 2 × p: lower, upper
        return new OracleEmulator(
            [new ContinuousParameter("theta1", bounds[0][0], bounds[1][0]), new ContinuousParameter("theta2", bounds[0][1], bounds[1][1])],
            offset, scale, basis,
            (theta, _) => [.. coef.Select(c => c[0] + c[1] * theta[0] + c[2] * theta[1])],
            Rows(inputs.GetProperty("emu_s2")), residuals);
    }

    internal static MeasurementModel Groups(JsonElement inputs, VariancePrior first, VariancePrior second)
    {
        double[] groups = Vector(inputs.GetProperty("groups"));
        return new MeasurementModel(new GaussianMeasurement([
            new VarianceGroup([.. Enumerable.Range(0, groups.Length).Where(j => groups[j] == 1)], first),
            new VarianceGroup([.. Enumerable.Range(0, groups.Length).Where(j => groups[j] == 2)], second)]), groups.Length);
    }

    private static MeasurementModel RsGroups(JsonElement inputs) => Groups(inputs, new InverseGammaVariance(10, 0.1), new InverseGammaVariance(10, 0.1));

    [Fact]
    public void BranchByShape_TheDeclaredFormulationHolds()
    {
        if (!Available) return;
        JsonElement o = Oracle.Value, branch = o.GetProperty("branch_by_shape"), inputs = o.GetProperty("inputs");
        // R: type 1 with a square discrepancy runs its multivariate branch, type 2 with a two-column one its functional branch.
        Assert.Equal("multivariate", branch.GetProperty("type1_square").GetProperty("r_branch").GetString());
        Assert.Equal("functional", branch.GetProperty("type2_narrow").GetProperty("r_branch").GetString());

        OracleEmulator emulator = Emulator(inputs);
        int m = emulator.Outputs;
        double[][] identity = [.. Enumerable.Range(0, m).Select(i => Enumerable.Range(0, m).Select(j => i == j ? 1.0 : 0.0).ToArray())];
        double[][] narrow = Rows(inputs.GetProperty("D_narrow"));
        // Corrected: a functional formulation with a square basis samples m discrepancy coefficients ...
        CalibrationCovariance square = CalibrationCovariance.Create(new FunctionalFormulation { Discrepancy = new FunctionalDiscrepancy(identity, 1.0) },
                                                                    RsGroups(inputs), emulator);
        Assert.Equal(m, square.DiscrepancyBasis!.GetLength(1));
        // ... and a multivariate one given a two-column discrepancy covariance is refused, never run as functional.
        string refusal = Assert.Throws<ArgumentException>(() =>
            CalibrationCovariance.Create(new MultivariateFormulation(identity, narrow), RsGroups(inputs), emulator)).Message;
        Assert.Contains($"hold {m} values", refusal);
        _out.WriteLine($"R: type 1 square → multivariate, type 2 narrow → functional; here: functional with {m} coefficients, and \"{refusal}\"");
    }

    [Fact]
    public void Units_TheCovarianceIsInPhysicalUnits()
    {
        if (!Available) return;
        JsonElement o = Oracle.Value, inputs = o.GetProperty("inputs");
        OracleEmulator emulator = Emulator(inputs);
        CalibrationCovariance covariance = CalibrationCovariance.Create(new FunctionalFormulation { Truncation = TruncationTerm.Diagonal },
                                                                        RsGroups(inputs), emulator);
        JsonElement[] rungs = [.. o.GetProperty("unit_mismatch").EnumerateArray()];
        JsonElement[] states = [.. o.GetProperty("run").EnumerateArray()];
        double worstPhysical = 0.0, closestR = double.PositiveInfinity;
        for (int t = 0; t < rungs.Length; t++)
        {
            double rLdet = rungs[t].GetProperty("r_ldet").GetDouble();
            Assert.Equal(rungs[t].GetProperty("dense_r_units").GetDouble(), rLdet, 1e-10);   // R's is the dense form in R's units
            var noise = new ComponentNoise(((ComponentNoise)emulator.Noise(0)).Factor, Vector(states[t].GetProperty("emu_vars")));
            double ldet = covariance.Factor(Vector(states[t].GetProperty("s2")), noise).LogDet;
            worstPhysical = Math.Max(worstPhysical, Math.Abs(ldet - rungs[t].GetProperty("dense_physical").GetDouble()));
            closestR = Math.Min(closestR, Math.Abs(ldet - rLdet));
        }
        _out.WriteLine($"{rungs.Length} rungs: this side vs the dense physical log-determinant ≤ {worstPhysical:E1}; vs R's ≥ {closestR:F3}");
        Assert.InRange(worstPhysical, 0, 1e-10);
        Assert.True(closestR > 0.1);
    }
}
