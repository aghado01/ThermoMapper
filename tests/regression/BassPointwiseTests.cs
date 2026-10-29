using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's pointwise log density (stage 4, S9; OQ17; D04), on draws from <see cref="BassChain"/> under every intercept
/// convention, σ² prior and g mode, on scalar, categorical and direct functional data. The retained content against the
/// Normal density at an explicit uncentred design; the collapsed content against explicit refits that leave each unit
/// out (the normal–inverse-gamma posterior from the rest's rows under the full data's prior, then the unit's
/// multivariate Student t), at the same structure and τ; each unit partition; the empty structure; the latent reports;
/// WAIC and PSIS-LOO reading the result with r_eff from the chains; and the refusals.
/// </summary>
public sealed class BassPointwiseTests
{
    private readonly ITestOutputHelper _out;
    public BassPointwiseTests(ITestOutputHelper output) => _out = output;

    public enum Kind { Scalar, Categorical, Functional }

    // Scalar: x₀ and x₁ on a five-value grid (n = 25), y = 1 + 2x₀ + 3(x₀ − ½)₊(½ − x₁)₊ + 0.3ε.
    // Categorical: x₀ on 24 values and z on three levels, y = 1 + 2x₀ + 1.5·1(z = b) − 0.8·1(z = c)·x₀ + 0.3ε.
    // Functional: x₀ on four values and x₁ on three (n = 12), t on six locations (N = 72),
    // y = 1 + 2x₀ + 1.5·t·x₁ + 2(x₀ − ½)₊(t − ½)₊ + 0.2ε.
    internal static BassData Data(Kind kind, int seed = 5)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        switch (kind)
        {
            case Kind.Scalar:
            {
                var x0 = new double[25];
                var x1 = new double[25];
                var y = new double[25];
                for (int i = 0; i < 25; i++)
                {
                    x0[i] = (i / 5) / 4.0;
                    x1[i] = (i % 5) / 4.0;
                    y[i] = 1 + 2 * x0[i] + 3 * Math.Max(x0[i] - 0.5, 0) * Math.Max(0.5 - x1[i], 0) + 0.3 * Normal.SampleStandard(rng);
                }
                return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
            }
            case Kind.Categorical:
            {
                var x0 = new double[24];
                var z = new int[24];
                var y = new double[24];
                for (int i = 0; i < 24; i++)
                {
                    x0[i] = (i / 3) / 7.0;
                    z[i] = i % 3;
                    y[i] = 1 + 2 * x0[i] + (z[i] == 1 ? 1.5 : 0) - (z[i] == 2 ? 0.8 * x0[i] : 0) + 0.3 * Normal.SampleStandard(rng);
                }
                return new BassData([new ContinuousPredictor("x0", x0), new CategoricalPredictor("z", z, ["a", "b", "c"])], new ScalarResponse(y));
            }
            default:
            {
                const int n = 12, m = 6;
                double[] t = [.. Enumerable.Range(0, m).Select(r => r / (m - 1.0))];
                var x0 = new double[n];
                var x1 = new double[n];
                var y = new double[n, m];
                for (int i = 0; i < n; i++)
                {
                    x0[i] = (i % 4) / 3.0;
                    x1[i] = (i / 4) / 2.0;
                    for (int r = 0; r < m; r++)
                        y[i, r] = 1 + 2 * x0[i] + 1.5 * t[r] * x1[i] + 2 * Math.Max(x0[i] - 0.5, 0) * Math.Max(t[r] - 0.5, 0) + 0.2 * Normal.SampleStandard(rng);
                }
                return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)],
                                    new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
            }
        }
    }

    internal static BassConfig Config(GPriorCentre intercept = GPriorCentre.Mean, Sigma2Prior? sigma2 = null, bool fixedG = false, int maxBasis = 6)
        => new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = maxBasis },
            Coefficients = fixedG ? new FixedG() : new SampledG(),
            Sigma2 = sigma2 ?? Sigma2Prior.ScaledToResponse(),
            Intercept = intercept,
        };

    /// <summary>Draws per chain from <see cref="BassChain"/>: <paramref name="kept"/> per chain, one every
    /// <paramref name="thin"/> sweeps after <paramref name="burn"/>, composed into one table.</summary>
    internal static (BassModel Model, BassTermTable Table, List<List<BassDraw>> Draws) Run(BassConfig config, BassData data, int seed,
                                                                                           int chains = 2, int burn = 100, int kept = 20, int thin = 5)
    {
        var model = BassModel.Create(config, data);
        var table = new BassTermTable();
        int[] seeds = SeedTree.Derive(seed, 3 * chains);
        var draws = new List<List<BassDraw>>();
        for (int c = 0; c < chains; c++)
        {
            var chain = new BassChain(model, [], Math.Log(1.0 / model.Resolved.Observations),
                                      new Xoshiro256PlusPlus(seeds[3 * c]), new Xoshiro256PlusPlus(seeds[3 * c + 1]));
            var composer = new Xoshiro256PlusPlus(seeds[3 * c + 2]);
            for (int i = 0; i < burn; i++) chain.Step();
            var list = new List<BassDraw>();
            for (int s = 0; s < kept; s++)
            {
                for (int i = 0; i < thin; i++) chain.Step();
                list.Add(chain.Compose(composer, table));
            }
            draws.Add(list);
        }
        return (model, table, draws);
    }

    // A term's uncentred column B at all N observations, case by case and location by location: a(xᵢ)·h(t_r).
    internal static double[] Column(BassModel model, BassTerm term)
    {
        int n = model.Data.Response.Cases, m = model.Data.Response.Locations;
        var a = new double[n];
        model.Basis.CaseColumn(term, model.Inputs.Cases, a);
        var h = new double[m];
        if (model.Inputs.Outputs is { } outputs) model.Basis.OutputColumn(term, outputs, h);
        else h[0] = 1.0;
        var b = new double[n * m];
        for (int i = 0; i < n; i++)
            for (int r = 0; r < m; r++) b[i * m + r] = a[i] * h[r];
        return b;
    }

    // The draw's Normal log term at each observation, from its explicit uncentred design.
    internal static double[] DenseRetained(BassModel model, BassTermTable table, BassDraw draw)
    {
        double[] y = model.Data.Response.Values.ToArray();
        double[][] columns = [.. draw.Structure.TermIds.Select(id => Column(model, table[id]))];
        var log = new double[y.Length];
        for (int k = 0; k < y.Length; k++)
        {
            double f = draw.Coefficients[0];
            for (int q = 0; q < columns.Length; q++) f += draw.Coefficients[q + 1] * columns[q][k];
            log[k] = Normal.MeanVariance(f, draw.Sigma2).LogPdf(y[k]);
        }
        return log;
    }

    // The explicit centred design X = [1, Z_c], each column centred over all N rows directly.
    private static double[,] CentredDesign(BassModel model, BassTermTable table, int[] ids)
    {
        int n = model.Data.Response.Observations;
        var x = new double[n, ids.Length + 1];
        for (int k = 0; k < n; k++) x[k, 0] = 1.0;
        for (int j = 0; j < ids.Length; j++)
        {
            double[] b = Column(model, table[ids[j]]);
            double mean = b.Average();
            for (int k = 0; k < n; k++) x[k, j + 1] = b[k] - mean;
        }
        return x;
    }

    /// <summary>
    /// ln p(y_U | y_rest) by an explicit refit on the rest's rows of X = [1, Z_c], under the full data's prior: precision
    /// Λ₀ = τXᵀX (its intercept row and column zero under the flat intercept, which leaves k₀ = p − 1 proper columns),
    /// centre μ₀ = (ȳ, 0) under level centring and 0 under the zero centre, and σ² ~ IG(a, b) (a = b = 0 for Jeffreys).
    /// The posterior is Λₙ = X_rᵀX_r + Λ₀, μₙ = Λₙ⁻¹(X_rᵀy_r + Λ₀μ₀), aₙ = a + (n_r + k₀ − p)/2,
    /// bₙ = b + ½(‖y_r − X_rμₙ‖² + (μₙ − μ₀)ᵀΛ₀(μₙ − μ₀)), and the unit's predictive the multivariate t
    /// lnΓ(aₙ + m/2) − lnΓ(aₙ) − (m/2)ln 2π − ½ln|Σ_U| + aₙ ln bₙ − (aₙ + m/2) ln(bₙ + ½eᵀΣ_U⁻¹e), with
    /// Σ_U = I + X_UΛₙ⁻¹X_Uᵀ and e = y_U − X_Uμₙ.
    /// </summary>
    private static double Refit(double[,] x, double[] y, int[] unit, double tau, GPriorCentre intercept, InverseGamma? prior)
    {
        int n = y.Length, p = x.GetLength(1), m = unit.Length;
        var inUnit = new bool[n];
        foreach (int i in unit) inUnit[i] = true;
        var lambda0 = new double[p, p];
        for (int a = 0; a < p; a++)
            for (int b = 0; b < p; b++)
            {
                double g = 0;
                for (int k = 0; k < n; k++) g += x[k, a] * x[k, b];
                lambda0[a, b] = intercept == GPriorCentre.Flat && (a == 0 || b == 0) ? 0.0 : tau * g;
            }
        var mu0 = new double[p];
        if (intercept == GPriorCentre.Mean) mu0[0] = y.Average();
        int k0 = intercept == GPriorCentre.Flat ? p - 1 : p;

        var lambdaN = (double[,])lambda0.Clone();
        var rhs = new double[p];
        for (int a = 0; a < p; a++)
            for (int b = 0; b < p; b++) rhs[a] += lambda0[a, b] * mu0[b];
        for (int k = 0; k < n; k++)
        {
            if (inUnit[k]) continue;
            for (int a = 0; a < p; a++)
            {
                rhs[a] += x[k, a] * y[k];
                for (int b = 0; b < p; b++) lambdaN[a, b] += x[k, a] * x[k, b];
            }
        }
        var factor = new DenseCholesky(p);
        factor.Decompose(lambdaN);
        Assert.True(factor.Status.PositiveDefinite);
        double[] muN = factor.Solve(rhs);

        double shapeN = (prior?.Shape ?? 0.0) + 0.5 * (n - m + k0 - p), scaleN = prior?.Scale ?? 0.0;
        double resid = 0, penalty = 0;
        for (int k = 0; k < n; k++)
        {
            if (inUnit[k]) continue;
            double e = y[k];
            for (int a = 0; a < p; a++) e -= x[k, a] * muN[a];
            resid += e * e;
        }
        for (int a = 0; a < p; a++)
            for (int b = 0; b < p; b++) penalty += (muN[a] - mu0[a]) * lambda0[a, b] * (muN[b] - mu0[b]);
        scaleN += 0.5 * (resid + penalty);

        var sigmaU = new double[m, m];
        var eU = new double[m];
        var solved = new double[m][];
        for (int u = 0; u < m; u++)
        {
            var row = new double[p];
            for (int a = 0; a < p; a++) row[a] = x[unit[u], a];
            solved[u] = factor.Solve(row);
            eU[u] = y[unit[u]];
            for (int a = 0; a < p; a++) eU[u] -= row[a] * muN[a];
        }
        for (int u = 0; u < m; u++)
            for (int v = 0; v < m; v++)
            {
                double q = 0;
                for (int a = 0; a < p; a++) q += x[unit[u], a] * solved[v][a];
                sigmaU[u, v] = (u == v ? 1.0 : 0.0) + q;
            }
        var unitFactor = new DenseCholesky(m);
        unitFactor.Decompose(sigmaU);
        double[] w = unitFactor.Solve(eU);
        double quad = 0;
        for (int u = 0; u < m; u++) quad += eU[u] * w[u];
        return SpecialFunctions.LogGamma(shapeN + 0.5 * m) - SpecialFunctions.LogGamma(shapeN) - 0.5 * m * Math.Log(2 * Math.PI)
               - 0.5 * unitFactor.LogDet + shapeN * Math.Log(scaleN) - (shapeN + 0.5 * m) * Math.Log(scaleN + 0.5 * quad);
    }

    private static void Close(double expected, double actual, double tolerance, string what)
        => Assert.True(Math.Abs(expected - actual) <= tolerance * Math.Max(1.0, Math.Abs(expected)), $"{what}: {actual:R} against {expected:R}");

    // Checks both contents of every draw against the dense references, on the partition given.
    private void CheckBoth(BassModel model, BassTermTable table, IReadOnlyList<IReadOnlyList<BassDraw>> draws, UnitPartition? units, string label)
    {
        UnitPartition partition = units ?? BassPointwise.DefaultUnits(model);
        PointwiseLogDensity retained = BassPointwise.LogDensity(model, table, draws, units);
        PointwiseLogDensity collapsed = BassPointwise.LogDensity(model, table, draws, units, PointwiseContent.CollapsedPredictive);
        BassDraw[] pooled = [.. draws.SelectMany(c => c)];
        Assert.Equal(pooled.Length, retained.Draws);
        Assert.Equal(partition.Units, retained.Values.GetLength(1));
        Assert.Equal(draws.Select(c => c.Count), retained.ChainLengths!);
        Assert.Equal(draws.Select(c => c.Count), collapsed.ChainLengths!);
        if (units is not null) Assert.Same(units, retained.Units);

        double[] y = model.Data.Response.Values.ToArray();
        int structures = pooled.Select(d => string.Join(",", d.Structure.TermIds)).Distinct().Count();
        for (int s = 0; s < pooled.Length; s++)
        {
            BassDraw d = pooled[s];
            double[] log = DenseRetained(model, table, d);
            double[,] x = CentredDesign(model, table, d.Structure.TermIds);
            for (int u = 0; u < partition.Units; u++)
            {
                int[] members = partition.Members[u];
                Close(members.Sum(i => log[i]), retained.Values[s, u], 1e-10, $"{label} retained, draw {s}, unit {u}");
                double refit = Refit(x, y, members, d.Tau, model.Target.Intercept, model.Target.Sigma2);
                Close(refit, collapsed.Values[s, u], 1e-8, $"{label} collapsed, draw {s} (M = {d.Structure.TermIds.Length}, τ = {d.Tau:G3}), unit {u}");
            }
        }
        _out.WriteLine($"{label}: {pooled.Length} draws over {structures} structures, {partition.Units} units");
    }

    public static IEnumerable<object[]> Conventions()
    {
        Sigma2Prior[] priors = [Sigma2Prior.ScaledToResponse(), Sigma2Prior.Jeffreys, Sigma2Prior.InverseGamma(2.0, 0.5)];
        foreach (Kind kind in Enum.GetValues<Kind>())
            foreach (GPriorCentre intercept in new[] { GPriorCentre.Mean, GPriorCentre.Zero, GPriorCentre.Flat })
                for (int p = 0; p < priors.Length; p++)
                    foreach (bool fixedG in new[] { false, true })
                        yield return [kind, intercept, p, fixedG];
    }

    [Theory]
    [MemberData(nameof(Conventions))]
    public void BothContents_MatchDenseReferencesAndRefits(Kind kind, GPriorCentre intercept, int prior, bool fixedG)
    {
        Sigma2Prior sigma2 = prior switch { 0 => Sigma2Prior.ScaledToResponse(), 1 => Sigma2Prior.Jeffreys, _ => Sigma2Prior.InverseGamma(2.0, 0.5) };
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = Run(Config(intercept, sigma2, fixedG), Data(kind), seed: 11, kept: 8);
        Assert.Contains(draws.SelectMany(c => c), d => d.Structure.TermIds.Length > 0);
        CheckBoth(model, table, draws, null, $"{kind}, {intercept}, {sigma2}, {(fixedG ? "fixed g" : "g")}");
    }

    [Fact]
    public void Partitions_MatchDenseReferencesAndSumOverTheirObservations()
    {
        // Scalar: each observation (the default), and consecutive pairs with a last unit of one.
        (BassModel scalar, BassTermTable st, List<List<BassDraw>> sd) = Run(Config(), Data(Kind.Scalar), seed: 21, kept: 6);
        Assert.Equal(25, BassPointwise.DefaultUnits(scalar).Units);
        CheckBoth(scalar, st, sd, null, "scalar, observations");
        CheckBoth(scalar, st, sd, new UnitPartition([.. Enumerable.Range(0, 25).Select(i => i / 2)]), "scalar, pairs");

        // Functional: each case (the default), each observation, and each location across the cases (not contiguous).
        (BassModel functional, BassTermTable ft, List<List<BassDraw>> fd) = Run(Config(), Data(Kind.Functional), seed: 22, kept: 6);
        UnitPartition cases = BassPointwise.DefaultUnits(functional);
        Assert.Equal(12, cases.Units);
        Assert.All(cases.Members, members => Assert.Equal(6, members.Length));
        Assert.Equal([6, 7, 8, 9, 10, 11], cases.Members[1]);
        CheckBoth(functional, ft, fd, null, "functional, cases");
        CheckBoth(functional, ft, fd, UnitPartition.EachObservation(72), "functional, observations");
        var locations = new UnitPartition([.. Enumerable.Range(0, 72).Select(k => k % 6)]);
        CheckBoth(functional, ft, fd, locations, "functional, locations");

        // The retained content factorizes over observations: a case's term is the sum of its observations' terms.
        PointwiseLogDensity byCase = BassPointwise.LogDensity(functional, ft, fd);
        PointwiseLogDensity byObservation = BassPointwise.LogDensity(functional, ft, fd, UnitPartition.EachObservation(72));
        for (int s = 0; s < byCase.Draws; s++)
            for (int i = 0; i < 12; i++)
                Close(Enumerable.Range(6 * i, 6).Sum(k => byObservation.Values[s, k]), byCase.Values[s, i], 1e-12, $"draw {s}, case {i}");
    }

    [Fact]
    public void EmptyStructure_MatchesUnderBothContents()
    {
        // The intercept-only draw, composed at M = 0 under each intercept convention, on scalar and functional data.
        foreach (Kind kind in new[] { Kind.Scalar, Kind.Functional })
            foreach (GPriorCentre intercept in new[] { GPriorCentre.Mean, GPriorCentre.Zero, GPriorCentre.Flat })
            {
                var model = BassModel.Create(Config(intercept), Data(kind));
                var moves = new BassStructureMoves(model, [], Math.Log(0.05));
                var table = new BassTermTable();
                var rng = new Xoshiro256PlusPlus(3);
                List<List<BassDraw>> draws = [[moves.Compose(moves.Start, rng, table), moves.Compose(moves.Start, rng, table)]];
                Assert.All(draws[0], d => Assert.Empty(d.Structure.TermIds));
                CheckBoth(model, table, draws, null, $"{kind}, {intercept}, M = 0");
            }
    }

    [Fact]
    public void ExactFits_RefusedUnderJeffreysAndMatchedUnderTheScaledDefault()
    {
        // y = 1 + 2(x₀ − ½)₊ exactly, which the hinge at x₀ = ½ spans. Under Jeffreys that structure lies outside F (D03's
        // propriety rule), so the collapsed content refuses its draw; under the scaled default it is proper, and both
        // contents match their references there: the deterministic-simulator regime.
        var x0 = new double[25];
        var x1 = new double[25];
        var y = new double[25];
        for (int i = 0; i < 25; i++)
        {
            x0[i] = (i / 5) / 4.0;
            x1[i] = (i % 5) / 4.0;
            y[i] = 1 + 2 * Math.Max(x0[i] - 0.5, 0);
        }
        var data = new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
        var hinge = new BassTerm([new HingeFactor(0, 1, 2)]);
        var jeffreys = BassModel.Create(Config(sigma2: Sigma2Prior.Jeffreys), data);
        var table = new BassTermTable();
        table.Intern(hinge);
        var exact = new BassDraw(new BassStructure([0]), 0.05, 1e-6, [1.0, 1.0], 1.0);
        Assert.Throws<ArgumentException>("draws", () => BassPointwise.LogDensity(jeffreys, table, [[exact]], content: PointwiseContent.CollapsedPredictive));
        BassPointwise.LogDensity(jeffreys, table, [[exact]]);   // the retained content reads the draw as given

        var scaled = BassModel.Create(Config(), data);
        var moves = new BassStructureMoves(scaled, [hinge], Math.Log(0.05));
        var composed = new BassTermTable();
        var rng = new Xoshiro256PlusPlus(4);
        CheckBoth(scaled, composed, [[moves.Compose(moves.Start, rng, composed), moves.Compose(moves.Start, rng, composed)]], null, "exact fit, scaled σ²");
    }

    [Fact]
    public void DrawStores_WithUnequalChains()
    {
        // The fit's own shape, a list of DrawStore per chain, here of 5, 4 and 3 draws.
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = Run(Config(), Data(Kind.Categorical), seed: 33, chains: 3, kept: 5);
        List<DrawStore<BassDraw>> stores = [.. draws.Select((chain, k) =>
        {
            var store = new DrawStore<BassDraw>();
            foreach (BassDraw d in chain.Take(5 - k)) store.Add(d);
            return store;
        })];
        CheckBoth(model, table, stores, null, "categorical, unequal chains");
        PointwiseLogDensity collapsed = BassPointwise.LogDensity(model, table, stores, content: PointwiseContent.CollapsedPredictive);
        Assert.Equal([5, 4, 3], collapsed.ChainLengths!);
        Assert.True(double.IsFinite(PredictiveChecks.PsisLoo(collapsed).Elpd));
    }

    [Fact]
    public void AFitsDraws_FeedBothEntryPoints()
    {
        // The run's own output: its merged table and per-chain stores, on the model rebuilt from its config and data.
        BassData data = Data(Kind.Functional);
        BassConfig config = Config() with { Sampling = new BassSampling { Chains = 2, BurnIn = 100, Samples = 6, MasterSeed = 3 } };
        BassFit fit = Bass.Fit(config, data);
        var model = BassModel.Create(config, data);
        CheckBoth(model, fit.Terms, fit.Draws, null, "a fit's draws");
        ResidualDependence check = BassResidualDependence.Of(model, fit.Terms, fit.Draws, weights: fit.DrawWeights)!;
        Assert.Equal(12, check.PerCase.Length);
    }

    [Fact]
    public void Latents_SayWhatWasConditionedAndWhatIntegrated()
    {
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = Run(Config(), Data(Kind.Scalar), seed: 31, kept: 3);
        PointwiseLogDensity retained = BassPointwise.LogDensity(model, table, draws);
        Assert.Equal(["structure", "prior precision τ", "coefficients a", "noise variance σ²"], retained.Latents.Select(l => l.Name));
        Assert.All(retained.Latents, l => Assert.Equal(LatentScope.Global, l.Scope));
        Assert.All(retained.Latents, l => Assert.Equal(LatentTreatment.Conditioned, l.Treatment));
        Assert.False(retained.LeavesUnitsOut);

        PointwiseLogDensity collapsed = BassPointwise.LogDensity(model, table, draws, content: PointwiseContent.CollapsedPredictive);
        Assert.Equal([LatentTreatment.Conditioned, LatentTreatment.Conditioned, LatentTreatment.IntegratedGivenOthers, LatentTreatment.IntegratedGivenOthers],
                     collapsed.Latents.Select(l => l.Treatment));
        Assert.True(collapsed.LeavesUnitsOut);

        // Fixed g holds τ, which is then no latent.
        (BassModel fixedModel, BassTermTable ft, List<List<BassDraw>> fd) = Run(Config(fixedG: true), Data(Kind.Scalar), seed: 32, kept: 3);
        Assert.Equal(["structure", "coefficients a", "noise variance σ²"], BassPointwise.LogDensity(fixedModel, ft, fd).Latents.Select(l => l.Name));
        Assert.All(fd.SelectMany(c => c), d => Assert.Equal(1.0 / 25, d.Tau, 15));
    }

    [Fact]
    public void PredictiveChecks_ReadBothContentsWithTheChains()
    {
        foreach (Kind kind in new[] { Kind.Scalar, Kind.Functional })
        {
            (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = Run(Config(), Data(kind), seed: 41, chains: 3, burn: 300, kept: 150, thin: 2);
            PointwiseLogDensity retained = BassPointwise.LogDensity(model, table, draws);
            PointwiseLogDensity collapsed = BassPointwise.LogDensity(model, table, draws, content: PointwiseContent.CollapsedPredictive);
            Assert.Equal([150, 150, 150], retained.ChainLengths!);

            PredictiveCheck waic = PredictiveChecks.Waic(retained);
            PredictiveCheck looRetained = PredictiveChecks.PsisLoo(retained);
            PredictiveCheck looCollapsed = PredictiveChecks.PsisLoo(collapsed);
            Assert.Throws<ArgumentException>(() => PredictiveChecks.Waic(collapsed));
            foreach (PredictiveCheck check in new[] { waic, looRetained, looCollapsed })
            {
                Assert.True(double.IsFinite(check.Elpd));
                Assert.Equal(retained.Values.GetLength(1), check.PerUnit.Length);
            }
            // r_eff comes from the chains: estimated per unit, not the independent draws' 1 everywhere.
            Assert.Contains(looRetained.RelativeEfficiency!, r => r != 1.0);
            Assert.Contains(looCollapsed.RelativeEfficiency!, r => r != 1.0);
            Assert.True(double.IsNaN(looCollapsed.EffectiveParameters));

            double threshold = PredictiveChecks.ParetoKThreshold(retained.Draws);
            int flaggedRetained = looRetained.ParetoK!.Count(k => k > threshold);
            int flaggedCollapsed = looCollapsed.ParetoK!.Count(k => k > threshold);
            _out.WriteLine($"{kind}: WAIC {waic.Elpd:F3}; PSIS-LOO retained {looRetained.Elpd:F3} ± {looRetained.MonteCarloError:F3} "
                           + $"({flaggedRetained} units above k̂ {threshold:F2}), collapsed {looCollapsed.Elpd:F3} "
                           + $"± {looCollapsed.JointMonteCarloError:F3} ({flaggedCollapsed} above)");
        }
    }

    [Fact]
    public void Refusals()
    {
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = Run(Config(), Data(Kind.Scalar), seed: 51, kept: 3);
        BassDraw good = draws[0][^1];

        Assert.Throws<ArgumentNullException>("model", () => BassPointwise.LogDensity(null!, table, draws));
        Assert.Throws<ArgumentNullException>("terms", () => BassPointwise.LogDensity(model, null!, draws));
        Assert.Throws<ArgumentNullException>("draws", () => BassPointwise.LogDensity(model, table, null!));
        Assert.Throws<ArgumentNullException>("model", () => BassPointwise.DefaultUnits(null!));
        Assert.Throws<ArgumentOutOfRangeException>("content", () => BassPointwise.LogDensity(model, table, draws, content: (PointwiseContent)7));
        Assert.Throws<ArgumentException>("units", () => BassPointwise.LogDensity(model, table, draws, UnitPartition.EachObservation(24)));
        Assert.Throws<ArgumentException>("draws", () => BassPointwise.LogDensity(model, table, []));
        Assert.Throws<ArgumentException>("draws", () => BassPointwise.LogDensity(model, table, [draws[0], []]));
        Assert.Throws<ArgumentException>("draws", () => BassPointwise.LogDensity(model, table, [[null!]]));
        Assert.Throws<ArgumentException>("draws", () => BassPointwise.LogDensity(model, table, [draws[0], null!]));

        void Refused(BassDraw bad, PointwiseContent content = PointwiseContent.RetainedParameters)
            => Assert.Throws<ArgumentException>("draws", () => BassPointwise.LogDensity(model, table, [[good, bad]], content: content));
        int[] two = table.Count >= 2 ? [0, 1] : throw new InvalidOperationException("The fixture needs two terms.");
        Refused(good with { Structure = default });
        Refused(good with { Structure = new BassStructure([table.Count]), Coefficients = [1.0, 1.0] });
        Refused(good with { Structure = new BassStructure([1, 0]), Coefficients = [1.0, 1.0, 1.0] });
        Refused(good with { Structure = new BassStructure([0, 0]), Coefficients = [1.0, 1.0, 1.0] });
        Refused(good with { Structure = new BassStructure(two), Coefficients = [1.0, 1.0] });
        Refused(good with { Coefficients = [.. good.Coefficients.Select((a, j) => j == 0 ? double.NaN : a)] });
        Refused(good with { Coefficients = [.. good.Coefficients.Select((a, j) => j == good.Coefficients.Length - 1 ? double.PositiveInfinity : a)] });
        Refused(good with { Sigma2 = 0.0 });
        Refused(good with { Sigma2 = double.NaN });
        Refused(good with { Sigma2 = double.PositiveInfinity });
        Refused(good with { Tau = 0.0 }, PointwiseContent.CollapsedPredictive);
        Refused(good with { Tau = double.PositiveInfinity }, PointwiseContent.CollapsedPredictive);
        BassPointwise.LogDensity(model, table, [[good with { Tau = double.NaN }]]);   // the retained content reads no τ

        // A table from other inputs: a term the model's basis does not have.
        var foreign = new BassTermTable();
        foreign.Intern(new BassTerm([new HingeFactor(5, 1, 0)]));
        BassDraw alien = good with { Structure = new BassStructure([0]), Coefficients = [1.0, 1.0] };
        Assert.Throws<ArgumentException>("terms", () => BassPointwise.LogDensity(model, foreign, [[alien]]));

        // A structure outside F: a term below the support threshold, which the collapsed content's design refuses.
        var sparse = new BassTermTable();
        sparse.Intern(new BassTerm([new HingeFactor(0, 1, model.Inputs.ContinuousKnots(0).Length - 1)]));
        BassDraw unsupported = good with { Structure = new BassStructure([0]), Coefficients = [1.0, 1.0] };
        BassPointwise.LogDensity(model, sparse, [[unsupported]]);   // the retained content needs no design
        Assert.Throws<ArgumentException>("draws", () => BassPointwise.LogDensity(model, sparse, [[unsupported]], content: PointwiseContent.CollapsedPredictive));

        // A unit that leaves no observation outside it to condition on, under any σ² prior: all of them, or under the flat
        // intercept, whose integral spends one, all but one. Under Jeffreys that is where σ²'s shape given the rest reaches 0.
        var all = new UnitPartition(new int[25]);
        var nearlyAll = new UnitPartition([.. Enumerable.Range(0, 25).Select(i => i == 0 ? 1 : 0)]);
        Assert.Throws<ArgumentException>("units", () => BassPointwise.LogDensity(model, table, draws, all, PointwiseContent.CollapsedPredictive));
        BassPointwise.LogDensity(model, table, draws, nearlyAll, PointwiseContent.CollapsedPredictive);
        BassPointwise.LogDensity(model, table, draws, all);   // the retained content conditions on everything
        foreach (Sigma2Prior prior in new[] { Sigma2Prior.Jeffreys, Sigma2Prior.ScaledToResponse() })
        {
            (BassModel mean, BassTermTable mt, List<List<BassDraw>> md) = Run(Config(sigma2: prior), Data(Kind.Scalar), seed: 52, kept: 2);
            Assert.Throws<ArgumentException>("units", () => BassPointwise.LogDensity(mean, mt, md, all, PointwiseContent.CollapsedPredictive));
            BassPointwise.LogDensity(mean, mt, md, nearlyAll, PointwiseContent.CollapsedPredictive);
            (BassModel flat, BassTermTable lt, List<List<BassDraw>> ld) = Run(Config(GPriorCentre.Flat, prior), Data(Kind.Scalar), seed: 53, kept: 2);
            Assert.Throws<ArgumentException>("units", () => BassPointwise.LogDensity(flat, lt, ld, nearlyAll, PointwiseContent.CollapsedPredictive));
            BassPointwise.LogDensity(flat, lt, ld, new UnitPartition([.. Enumerable.Range(0, 25).Select(i => i < 2 ? 1 : 0)]),
                                     PointwiseContent.CollapsedPredictive);
        }
    }
}
