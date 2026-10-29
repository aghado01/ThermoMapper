using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass.Calibration;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// An emulator in physical units for calibration's parity tests: y_m + y_s ∘ (K w(θ)) with w given, component noise over
/// y_s ∘ K at each draw's variances, and declared truncation residuals.
/// </summary>
internal sealed class OracleEmulator : ICalibrationEmulator
{
    private readonly double[] _offset, _scale;
    private readonly double[,] _basis, _factor;
    private readonly Func<double[], int[], double[]> _w;
    private readonly double[][] _variances;
    private readonly double[,]? _residuals;
    private readonly List<double[]> _calls = [];

    public OracleEmulator(IReadOnlyList<CalibrationParameter> parameters, double[] offset, double[] scale, double[,] basis,
                          Func<double[], int[], double[]> w, double[][] variances, double[,]? residuals, double[,]? factor = null)
    {
        Parameters = parameters;
        _offset = offset;
        _scale = scale;
        _basis = basis;
        _w = w;
        _variances = variances;
        _residuals = residuals;
        int m = offset.Length, r = basis.GetLength(1);
        _factor = factor ?? new double[m, r];
        if (factor is null)
            for (int l = 0; l < m; l++)
                for (int k = 0; k < r; k++) _factor[l, k] = scale[l] * basis[l, k];
    }

    public IReadOnlyList<CalibrationParameter> Parameters { get; }
    public int DrawCount => _variances.Length;
    public IReadOnlyList<double>? DrawWeights => null;
    public int Outputs => _offset.Length;

    /// <summary>The continuous inputs of every prediction, in call order.</summary>
    public IReadOnlyList<double[]> Calls { get { lock (_calls) return [.. _calls]; } }

    public EmulatorNoise Noise(int draw) => new ComponentNoise(_factor, (double[])_variances[draw].Clone());
    public double[,]? TruncationResiduals() => (double[,]?)_residuals?.Clone();

    public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output)
    {
        lock (_calls) _calls.Add(continuous.ToArray());
        double[] w = _w(continuous.ToArray(), levels.ToArray());
        for (int l = 0; l < Outputs; l++)
        {
            double s = 0.0;
            for (int k = 0; k < w.Length; k++) s += _basis[l, k] * w[k];
            output[l] = _offset[l] + _scale[l] * s;
        }
    }
}

/// <summary>
/// Calibration's covariance formulations (stage 7, C1; D09) against the C-O2 references (<c>bass_covariance_oracle.R</c>):
/// the functional formulation's diagonal option at C-O2's well-conditioned case, and its low-rank solver at the moderate
/// and extreme conditionings, against R's dense Cholesky solve and log-determinant, with the structured backward error
/// where no dense reference certifies; the multivariate formulation against its dense reference; and units, the
/// physical log-likelihood on C-O2's θ grid against R's physical one, whose posterior the standardized one shares, with
/// R's units, for contrast, far off. Opt-in: it runs where the R toolchain and the BASS checkout are present.
/// </summary>
public sealed class BassCalibrationCovarianceParityTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationCovarianceParityTests(ITestOutputHelper output) => _out = output;

    private static bool Available => BassTargetParityTests.Available;

    // The stage-1 checks' fixture: y_m = 5 + t, y_s = 0.5 + 2t, K = (sin πt, cos πt) at 20 locations.
    internal static double[,] Fixture()
    {
        var b = new double[20, 4];
        for (int l = 0; l < 20; l++)
        {
            double t = l / 19.0;
            b[l, 0] = 5 + t;
            b[l, 1] = 0.5 + 2 * t;
            b[l, 2] = Math.Sin(Math.PI * t);
            b[l, 3] = Math.Cos(Math.PI * t);
        }
        return b;
    }

    internal static (double[] Offset, double[] Scale, double[,] Basis) Pieces()
    {
        double[,] f = Fixture();
        int m = f.GetLength(0);
        var basis = new double[m, 2];
        for (int l = 0; l < m; l++) (basis[l, 0], basis[l, 1]) = (f[l, 2], f[l, 3]);
        return ([.. Enumerable.Range(0, m).Select(l => f[l, 0])], [.. Enumerable.Range(0, m).Select(l => f[l, 1])], basis);
    }

    private static readonly Lazy<JsonElement> Oracle = new(() => BassOutputTruncationParityTests.Run("oracles/bass_covariance_oracle.R", Fixture()));

    internal static double[] Vector(JsonElement e) => [.. e.EnumerateArray().Select(v => v.GetDouble())];
    internal static double[][] Rows(JsonElement e) => [.. e.EnumerateArray().Select(Vector)];

    private static double RelativeError(double[] a, double[] reference)
        => a.Zip(reference, (x, y) => Math.Abs(x - y)).Max() / reference.Max(Math.Abs);

    // te (m × 10) as the emulator's residual rows in physical units, y_s ∘ te[, i].
    private static double[,] PhysicalResiduals(double[][] te, double[] scale, bool centredForRsCovariance = false)
    {
        int m = te.Length, n = te[0].Length;
        var r = new double[n, m];
        for (int l = 0; l < m; l++)
        {
            double mean = centredForRsCovariance ? te[l].Average() : 0.0;
            // (1/n)·Σ rᵢrᵢᵀ over rᵢ = √(n/(n − 1))·y_s ∘ (eᵢ − ē) is R's cov(): the residual set fed R's estimator.
            double factor = centredForRsCovariance ? Math.Sqrt(n / (n - 1.0)) : 1.0;
            for (int i = 0; i < n; i++) r[i, l] = factor * scale[l] * (te[l][i] - mean);
        }
        return r;
    }

    private static MeasurementModel Halves(int m, double[] groups)
    {
        int[] first = [.. Enumerable.Range(0, m).Where(j => groups[j] == 1)], second = [.. Enumerable.Range(0, m).Where(j => groups[j] == 2)];
        return new MeasurementModel(new GaussianMeasurement([new VarianceGroup(first, new InverseGammaVariance(10, 0.1)),
                                                             new VarianceGroup(second, new InverseGammaVariance(10, 0.1))]), m);
    }

    [Fact]
    public void Functional_MatchesTheDenseReferences_AtC02sConditionings()
    {
        if (!Available) return;
        JsonElement o = Oracle.Value, inputs = o.GetProperty("inputs");
        (double[] offset, double[] scale, double[,] basis) = Pieces();
        int m = offset.Length;
        double[][] te = Rows(inputs.GetProperty("te"));
        double[] v = Vector(inputs.GetProperty("v")), s2 = Vector(inputs.GetProperty("s2_g"));
        MeasurementModel measurement = Halves(m, Vector(inputs.GetProperty("groups")));
        JsonElement[] cases = [.. o.GetProperty("functional").EnumerateArray()];

        // The well case through the formulation: the diagonal option (R's y_s² ∘ var(te)), the groups and the components.
        var emulator = new OracleEmulator([new ContinuousParameter("theta", 0, 1)], offset, scale, basis, (_, _) => [0.0, 0.0],
                                          [[0.3, 0.7]], PhysicalResiduals(te, scale));
        CovarianceFactor well = CalibrationCovariance.Create(new FunctionalFormulation { Truncation = TruncationTerm.Diagonal }, measurement, emulator)
                                                     .Factor(s2, emulator.Noise(0));
        JsonElement reference = cases[0].GetProperty("reference");
        double ldetGap = Math.Abs(well.LogDet - reference.GetProperty("ldet").GetDouble());
        double solveGap = RelativeError(well.Solve(v), Vector(reference.GetProperty("solve")));
        _out.WriteLine($"well, through the formulation: ldet {ldetGap:E1}, solve {solveGap:E1}");
        Assert.InRange(ldetGap, 0, 1e-12);
        Assert.InRange(solveGap, 0, 1e-12);

        // Each case's d and e on the solver the formulation uses, against the dense reference where it certifies and by
        // the structured backward error everywhere.
        var u = new double[m, 2];
        for (int l = 0; l < m; l++) (u[l, 0], u[l, 1]) = (scale[l] * basis[l, 0], scale[l] * basis[l, 1]);
        double[] tolerance = [1e-12, 1e-8, double.PositiveInfinity];
        for (int c = 0; c < 3; c++)
        {
            JsonElement inputCase = inputs.GetProperty("cases").GetProperty(cases[c].GetProperty("case").GetString()!);
            double[] d = Vector(inputCase.GetProperty("d")), e = Vector(inputCase.GetProperty("e"));
            var sigma = new DiagonalPlusLowRank(d, u, e);
            double[] x = sigma.Solve(v);
            reference = cases[c].GetProperty("reference");
            double ldet = Math.Abs(sigma.LogDet - reference.GetProperty("ldet").GetDouble());
            double solve = RelativeError(x, Vector(reference.GetProperty("solve")));
            double backward = 0.0;
            for (int i = 0; i < m; i++)
            {
                double ux = 0.0, uxAbs = 0.0;
                for (int k = 0; k < 2; k++)
                {
                    double ck = 0.0, ckAbs = 0.0;
                    for (int j = 0; j < m; j++) { ck += u[j, k] * x[j]; ckAbs += Math.Abs(u[j, k]) * Math.Abs(x[j]); }
                    ux += u[i, k] * e[k] * ck;
                    uxAbs += Math.Abs(u[i, k]) * e[k] * ckAbs;
                }
                backward = Math.Max(backward, Math.Abs(d[i] * x[i] + ux - v[i]) / (Math.Abs(d[i] * x[i]) + uxAbs + Math.Abs(v[i])));
            }
            _out.WriteLine($"{cases[c].GetProperty("case").GetString(),-8} (cond {cases[c].GetProperty("condition").GetDouble():E0}): ldet {ldet:E1}, " +
                           $"solve {solve:E1}, backward {backward:E1}; R's low rank {cases[c].GetProperty("low_rank").GetProperty("solve_rel").GetDouble():E1}");
            Assert.InRange(ldet, 0, tolerance[c]);
            Assert.InRange(solve, 0, tolerance[c]);
            Assert.InRange(backward, 0, 1e-14);
        }
    }

    [Fact]
    public void Multivariate_MatchesItsDenseReference()
    {
        if (!Available) return;
        JsonElement o = Oracle.Value, inputs = o.GetProperty("inputs");
        (double[] offset, double[] scale, double[,] basis) = Pieces();
        int m = offset.Length;
        double[][] te = Rows(inputs.GetProperty("te"));
        var emulator = new OracleEmulator([new ContinuousParameter("theta", 0, 1)], offset, scale, basis, (_, _) => [0.0, 0.0],
                                          [[0.3, 0.7]], PhysicalResiduals(te, scale, centredForRsCovariance: true));
        var formulation = new MultivariateFormulation(Rows(inputs.GetProperty("R_groups")), Rows(inputs.GetProperty("Sigma_delta")));
        CovarianceFactor factor = CalibrationCovariance.Create(formulation, Halves(m, Vector(inputs.GetProperty("groups"))), emulator)
                                                       .Factor(Vector(inputs.GetProperty("s2_g")), emulator.Noise(0));
        JsonElement mv = o.GetProperty("multivariate");
        double ldet = Math.Abs(factor.LogDet - mv.GetProperty("ldet").GetDouble());
        double solve = RelativeError(factor.Solve(Vector(inputs.GetProperty("v"))), Vector(mv.GetProperty("solve")));
        _out.WriteLine($"multivariate (cond {mv.GetProperty("condition").GetDouble():E0}): ldet {ldet:E1}, solve {solve:E1}");
        Assert.InRange(ldet, 0, 1e-10);
        Assert.InRange(solve, 0, 1e-10);
    }

    [Fact]
    public void Units_ThePhysicalLikelihoodIsRsPhysicalOne()
    {
        if (!Available) return;
        JsonElement o = Oracle.Value, inputs = o.GetProperty("inputs"), units = o.GetProperty("units");
        (double[] offset, double[] scale, double[,] basis) = Pieces();
        int m = offset.Length;
        double[][] te = Rows(inputs.GetProperty("te"));
        var emulator = new OracleEmulator([new ContinuousParameter("theta", 10, 20)], offset, scale, basis,
                                          (theta, _) => [Math.Sin(theta[0] / 3), Math.Cos(theta[0] / 4)],
                                          [Vector(inputs.GetProperty("e_units"))], PhysicalResiduals(te, scale));
        CovarianceFactor factor = CalibrationCovariance.Create(new FunctionalFormulation { Truncation = TruncationTerm.Diagonal },
                                                               Halves(m, Vector(inputs.GetProperty("groups"))), emulator)
                                                       .Factor(Vector(inputs.GetProperty("s2_g")), emulator.Noise(0));
        double[] y = Vector(units.GetProperty("y")), grid = Vector(units.GetProperty("theta"));
        double[] rPhysical = Vector(units.GetProperty("loglik_physical"));
        var ll = new double[grid.Length];
        var pred = new double[m];
        for (int g = 0; g < grid.Length; g++)
        {
            emulator.Predict(0, [grid[g]], [], pred);
            double[] r = [.. Enumerable.Range(0, m).Select(j => y[j] - pred[j])];
            ll[g] = -0.5 * (factor.LogDet + factor.InverseQuadraticForm(r));
        }
        double gap = ll.Zip(rPhysical, (a, b) => Math.Abs(a - b)).Max();
        double max = ll.Max();
        double[] w = [.. ll.Select(l => Math.Exp(l - max))];
        double total = w.Sum();
        double[] posterior = [.. w.Select(x => x / total)];
        double vsPhysical = posterior.Zip(Vector(units.GetProperty("posterior_physical")), (a, b) => Math.Abs(a - b)).Max();
        double vsStandardized = posterior.Zip(Vector(units.GetProperty("posterior_standardized")), (a, b) => Math.Abs(a - b)).Max();
        double tvRUnits = 0.5 * posterior.Zip(Vector(units.GetProperty("posterior_r_units")), (a, b) => Math.Abs(a - b)).Sum();
        _out.WriteLine($"log-likelihood vs R's physical {gap:E1}; posterior vs physical {vsPhysical:E1}, vs standardized {vsStandardized:E1}; " +
                       $"vs R's units TV {tvRUnits:F2}");
        Assert.InRange(gap, 0, 1e-9);
        Assert.InRange(vsPhysical, 0, 1e-12);
        Assert.InRange(vsStandardized, 0, 1e-12);
        Assert.True(tvRUnits > 0.5);
    }
}
