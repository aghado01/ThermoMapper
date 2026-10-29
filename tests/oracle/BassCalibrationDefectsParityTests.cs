using System;
using System.Linq;
using System.Text.Json;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// R's calibration defects (C-O1, <c>bass_calibration_defects_oracle.R</c>) against D09's corrections (stage 7). The
/// oracle runs <c>calibrate.bassBasis</c> as written for one iteration on a mock emulator whose two components are
/// linear in θ; each test asserts R's defect as the fixture records it and the corrected behaviour on this side, on the
/// same emulator rebuilt from the fixture's inputs: the branch chosen by shape, the units of the covariance, the stale
/// likelihood after a discrepancy draw, the coordinates θ is reported in and the all-groups <c>s2.df</c> switch. Opt-in:
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

    // R's run as a config: the functional formulation with R's diagonal truncation term, its two-column discrepancy at
    // R's fixed scale (s_v = 1 reproduces R), and both groups IG(5, 0.05).
    internal static CalibrationModel Model(JsonElement inputs, OracleEmulator emulator)
    {
        MeasurementModel halves = RsGroups(inputs);
        var config = new CalibrationConfig
        {
            Parameters = emulator.Parameters,
            Measurement = new GaussianMeasurement([new VarianceGroup([.. halves.OutputsOf(0)], new InverseGammaVariance(10, 0.1)),
                                                   new VarianceGroup([.. halves.OutputsOf(1)], new InverseGammaVariance(10, 0.1))]),
            Covariance = new FunctionalFormulation
            {
                Truncation = TruncationTerm.Diagonal,
                Discrepancy = new FunctionalDiscrepancy(Rows(inputs.GetProperty("D_narrow")), 1.0),
            },
        };
        return new CalibrationModel(config, Vector(inputs.GetProperty("y")), emulator);
    }

    /// <summary>A θ move that always leaves the box, so the chain rejects it, as R's rung 5 rejected its θ move.</summary>
    private sealed class Rejected : IRjMove<CalibrationState>
    {
        public string Key => "theta";
        public string ReverseKey => "theta";
        public double Weight(CalibrationState state) => 1.0;
        public Proposal<CalibrationState>? Propose(CalibrationState current, Xoshiro256PlusPlus rng)
            => new Proposal<CalibrationState>(current.WithUnit([.. current.Unit.ToArray().Select(u => u + 2)]), 0.0, 0.0);
    }

    [Fact]
    public void StaleLikelihood_TheCachedLikelihoodIsAlwaysTheStatesOwn()
    {
        if (!Available) return;
        JsonElement o = Oracle.Value, stale = o.GetProperty("stale_likelihood"), inputs = o.GetProperty("inputs");
        bool[] rejected = [.. stale.GetProperty("theta_rejected").EnumerateArray().Select(v => v.GetBoolean())];
        double[] gap = Vector(stale.GetProperty("cached_minus_recomputed"));
        // R: the rung whose θ move was rejected kept the likelihood from before its discrepancy draw.
        Assert.Contains(true, rejected);
        for (int t = 0; t < gap.Length; t++)
            if (rejected[t]) Assert.True(Math.Abs(gap[t]) > 1e-3);
            else Assert.InRange(Math.Abs(gap[t]), 0, 1e-10);

        // Corrected: R's iteration on each rung's state (the discrepancy draw, then a rejected θ move) leaves the chain's
        // cached evaluation equal to a fresh one of its state.
        OracleEmulator emulator = Emulator(inputs);
        CalibrationModel model = Model(inputs, emulator);
        JsonElement[] states = [.. o.GetProperty("run").EnumerateArray()];
        double worst = 0.0;
        for (int t = 0; t < states.Length; t++)
        {
            CalibrationTarget target = model.Target(states[t].GetProperty("ii").GetInt32() - 1);
            var start = new CalibrationState(Vector(states[t].GetProperty("theta_start")), [], Vector(states[t].GetProperty("s2")), [0.0, 0.0]);
            var chain = new ReversibleJumpChain<CalibrationState>([new Rejected()], target, start, new Xoshiro256PlusPlus(t), Temperature.Posterior);
            var discrepancy = new DiscrepancyBlock(target);
            var rng = new Xoshiro256PlusPlus(100 + t);
            for (int i = 0; i < 3; i++)
            {
                discrepancy.Update(chain, rng);
                chain.Step();
                CalibrationState held = chain.Current;
                var copy = new CalibrationState(held.Unit.ToArray(), held.Levels.ToArray(), held.Variances.ToArray(), held.Coefficients.ToArray());
                double cached = ((CalibrationEvaluation)chain.CurrentEvaluation).LogLikelihoodValue;
                worst = Math.Max(worst, Math.Abs(cached - target.EvaluateState(copy).LogLikelihoodValue));
                Assert.Equal(Temperature.Posterior.LogDensity(target.EvaluateState(copy)), chain.CurrentLogDensity, 1e-12);
            }
            Assert.Equal(0, chain.Accepted);
        }
        _out.WriteLine($"R: cached minus recomputed {string.Join(" ", gap.Select(g => g.ToString("F3")))}; this side, cached minus fresh ≤ {worst:E1}");
        Assert.InRange(worst, 0, 1e-12);
    }

    [Fact]
    public void Coordinates_ThetaIsReportedAsTheEmulatorSawIt()
    {
        if (!Available) return;
        JsonElement o = Oracle.Value, coordinates = o.GetProperty("coordinates"), inputs = o.GetProperty("inputs");
        double[][] rRange = Rows(coordinates.GetProperty("r_returned_range"));      // 2 × p: R's returned min and max
        double[][] seen = Rows(coordinates.GetProperty("emulator_input_range"));
        // R: returned θ in [0, 1] while the emulator saw physical values; its start is not where it first predicted.
        Assert.All(rRange.SelectMany(r => r), v => Assert.InRange(v, 0.0, 1.0));
        Assert.True(seen[0][0] > 10 && seen[1][0] < 20);
        Assert.True(coordinates.GetProperty("r_start_vs_first_call").GetDouble() > 1);
        Assert.InRange(coordinates.GetProperty("corrected_start_vs_first_call").GetDouble(), 0, 1e-12);

        OracleEmulator emulator = Emulator(inputs);
        CalibrationModel model = Model(inputs, emulator);
        CalibrationTarget target = model.Target(0);
        foreach (JsonElement state in o.GetProperty("run").EnumerateArray())
        {
            var s = new CalibrationState(Vector(state.GetProperty("theta_start")), [], Vector(state.GetProperty("s2")), [0.0, 0.0]);
            _ = target.EvaluateState(s);
            double[] call = emulator.Calls[^1];
            double[] reported = model.Physical(s);
            Assert.Equal(call, reported);
            for (int k = 0; k < 2; k++)
            {
                var p = (ContinuousParameter)model.Config.Parameters[k];
                Assert.Equal(p.Lower + s.Unit[k] * (p.Upper - p.Lower), reported[k], 1e-12);
                Assert.InRange(reported[k], p.Lower, p.Upper);
            }
        }
        _out.WriteLine($"R returns θ₁ in [{rRange[0][0]:F2}, {rRange[1][0]:F2}] while its emulator saw [{seen[0][0]:F2}, {seen[1][0]:F2}]; " +
                       "here θ is reported as the emulator saw it, through the recorded transform");
    }

    [Fact]
    public void S2DfSwitch_EachGroupKeepsItsOwnProperPrior()
    {
        if (!Available) return;
        JsonElement o = Oracle.Value, sw = o.GetProperty("s2_df_switch"), inputs = o.GetProperty("inputs");
        double[] sigma2 = Vector(sw.GetProperty("sigma2"));
        double[] r = Vector(sw.GetProperty("r_group2")), corrected = Vector(sw.GetProperty("corrected_group2"));
        // R: with s2.df = (0, 10) the second group, whose IG is proper, takes −log(1 + σ²) too; with (10, 10) its IG.
        for (int i = 0; i < sigma2.Length; i++) Assert.Equal(-Math.Log(1 + sigma2[i]), r[i], 1e-15);
        Assert.Equal(corrected, Vector(sw.GetProperty("r_group2_when_all_positive")));

        // Corrected, per group: the first group (s2.df = 0) cannot declare an inverse gamma and declares the half-Cauchy;
        // the second keeps IG(5, 0.05), whose log density is R's corrected kernel up to its normalizing constant.
        Assert.Contains("half-Cauchy", Assert.Throws<ArgumentException>(() => Groups(inputs, new InverseGammaVariance(0, 0.1), new InverseGammaVariance(10, 0.1))).Message);
        MeasurementModel groups = Groups(inputs, new HalfCauchyDeviation(0.1), new InverseGammaVariance(10, 0.1));
        Assert.True(groups.IsProper);
        double constant = groups.LogPrior(1, sigma2[0]) - corrected[0];
        for (int i = 1; i < sigma2.Length; i++) Assert.Equal(corrected[i], groups.LogPrior(1, sigma2[i]) - constant, 1e-10);
        _out.WriteLine($"R's group 2 at σ² = {string.Join(", ", sigma2)}: {string.Join(" ", r.Select(v => v.ToString("F3")))}; its IG and this side's: " +
                       $"{string.Join(" ", corrected.Select(v => v.ToString("F2")))}");
    }

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
