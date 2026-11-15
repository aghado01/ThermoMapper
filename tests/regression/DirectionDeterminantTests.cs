using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E36's slice 4: the determinant direction prior, π(Θ | M) = Π π(θ_m)·det(ΘᵀΘ)/c_M, a labelled option beside the paper's independent
/// directions (decision 2), with c_M = p!/((p − M)!p^M) stated and its tilt on M reported (<see cref="DirectionDeterminant"/>).
/// </summary>
public sealed class DirectionDeterminantTests(ITestOutputHelper output)
{
    // A directional fit on p inputs uniform on [0, 1] with A = a, under the slice 2 fixture's proper priors.
    private static BassModel Model(int p, int a, DirectionPrior directions = DirectionPrior.Independent, int maxBasis = 8)
    {
        var rng = new Xoshiro256PlusPlus(3702 + p);
        const int n = 60;
        ContinuousPredictor[] x = [.. Enumerable.Range(0, p).Select(j => new ContinuousPredictor($"x{j}", [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())], new InputDomain(0.0, 1.0)))];
        double[] y = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        BassConfig config = DirectionalChainTests.Config(maxBasis);
        config = config with { Structure = config.Structure with { MaxInteraction = a, Directional = new DirectionalBasis { Directions = directions } } };
        return BassModel.Create(config, new BassData([.. x], new ScalarResponse(y)));
    }

    private static BassConfig Determinant(int maxBasis = 4, BassCountPrior? count = null)
    {
        BassConfig config = DirectionalChainTests.Config(maxBasis);
        return config with
        {
            Structure = config.Structure with { Directional = new DirectionalBasis { Directions = DirectionPrior.Determinant } },
            Count = count ?? config.Count,
        };
    }

    private static double[,] Gram(IReadOnlyList<ProjectionDirection> directions, int p)
    {
        int m = directions.Count;
        var g = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                g[i, j] = Enumerable.Range(0, p).Sum(v => directions[i].WeightOn(v) * directions[j].WeightOn(v));
        return g;
    }

    // The determinant by Gaussian elimination with partial pivoting, apart from the engine's Cholesky.
    private static double HandDeterminant(double[,] source)
    {
        int m = source.GetLength(0);
        var a = (double[,])source.Clone();
        double det = 1.0;
        for (int c = 0; c < m; c++)
        {
            int pivot = c;
            for (int r = c + 1; r < m; r++)
                if (Math.Abs(a[r, c]) > Math.Abs(a[pivot, c])) pivot = r;
            if (a[pivot, c] == 0.0) return 0.0;
            if (pivot != c)
            {
                for (int k = 0; k < m; k++) (a[c, k], a[pivot, k]) = (a[pivot, k], a[c, k]);
                det = -det;
            }
            det *= a[c, c];
            for (int r = c + 1; r < m; r++)
            {
                double f = a[r, c] / a[c, c];
                for (int k = c; k < m; k++) a[r, k] -= f * a[c, k];
            }
        }
        return det;
    }

    // sin² of θ's angle to the span of `others`: its squared residual after modified Gram–Schmidt against them.
    private static double SineSquared(ProjectionDirection theta, IReadOnlyList<ProjectionDirection> others, int p)
    {
        var basis = new List<double[]>();
        foreach (ProjectionDirection d in others)
        {
            double[] v = [.. Enumerable.Range(0, p).Select(d.WeightOn)];
            foreach (double[] q in basis)
            {
                double dot = Enumerable.Range(0, p).Sum(k => v[k] * q[k]);
                for (int k = 0; k < p; k++) v[k] -= dot * q[k];
            }
            double norm = Math.Sqrt(v.Sum(t => t * t));
            if (norm > 1e-12) basis.Add([.. v.Select(t => t / norm)]);
        }
        double[] r = [.. Enumerable.Range(0, p).Select(theta.WeightOn)];
        foreach (double[] q in basis)
        {
            double dot = Enumerable.Range(0, p).Sum(k => r[k] * q[k]);
            for (int k = 0; k < p; k++) r[k] -= dot * q[k];
        }
        return r.Sum(t => t * t);
    }

    // ---- c_M ------------------------------------------------------------------------------------------------------------------------

    [Fact]
    public void TheNormalizer_IsTheProductForm_AndMonteCarloAgrees_ForDenseAndSparseDirections()
    {
        double worstForm = 0.0;
        for (int p = 1; p <= 12; p++)
        {
            double product = 1.0;
            for (int m = 0; m <= p; m++)
            {
                double c = Math.Exp(DirectionDeterminant.LogNormalizer(p, m));
                worstForm = Math.Max(worstForm, Math.Abs(c - product) / product);
                product *= (double)(p - m) / p;
            }
            Assert.Equal(double.NegativeInfinity, DirectionDeterminant.LogNormalizer(p, p + 1));
        }
        output.WriteLine($"c_M against Π_(j < M) (p − j)/p, p = 1 … 12, M = 0 … p: worst relative {worstForm:E2}");
        Assert.True(worstForm <= 1e-13, $"{worstForm:E2}");

        // E[det ΘᵀΘ] by Monte Carlo: dense directions drawn here (a normalized Gaussian over every coordinate), sparse ones from the engine's
        // term measure (a uniform on 1 … A, 𝒥 uniform, θ_𝒥 uniform on its subsphere). Both laws give E[θθᵀ] = I/p, so c_M is the one form.
        const int draws = 40_000;
        var rng = new Xoshiro256PlusPlus(3701);
        double worstZ = 0.0;
        void Read(string law, int p, int m, Func<ProjectionDirection> draw)
        {
            double sum = 0.0, squares = 0.0;
            for (int k = 0; k < draws; k++)
            {
                double det = HandDeterminant(Gram([.. Enumerable.Range(0, m).Select(_ => draw())], p));
                sum += det;
                squares += det * det;
            }
            double mean = sum / draws, sd = Math.Sqrt((squares - draws * mean * mean) / (draws - 1));
            double c = Math.Exp(DirectionDeterminant.LogNormalizer(p, m)), z = (mean - c) / (sd / Math.Sqrt(draws));
            worstZ = Math.Max(worstZ, Math.Abs(z));
            output.WriteLine($"{law,-22} p {p,2}, M {m,2}: c_M {c:F6}, Monte Carlo {mean:F6} ± {sd / Math.Sqrt(draws):F6}, z {z,6:F2}");
        }
        ProjectionDirection Dense(int p)
        {
            double[] g = [.. Enumerable.Range(0, p).Select(_ => Normal.SampleStandard(rng))];
            double norm = Math.Sqrt(g.Sum(t => t * t));
            return new ProjectionDirection(Enumerable.Range(0, p), g.Select(t => t / norm));
        }
        foreach ((int p, int m) in new[] { (3, 2), (3, 3), (5, 3), (6, 6), (10, 5) }) Read("dense", p, m, () => Dense(p));
        foreach ((int p, int a, int m) in new[] { (6, 3, 2), (6, 3, 3), (6, 3, 6), (10, 3, 3), (10, 3, 5), (4, 2, 4), (6, 1, 3) })
        {
            BassModel model = Model(p, a);
            Read($"engine's, A = {a}", p, m, () => model.Prior.DrawTerm(rng)!.Directional!.Direction);
        }
        output.WriteLine($"worst |z| {worstZ:F2}");
        Assert.True(worstZ < 4.0, $"{worstZ:F2}");
    }

    // ---- the determinant ------------------------------------------------------------------------------------------------------------

    [Fact]
    public void TheDeterminant_IsTheGramsByHand_ABirthsFactorIsTheSineSquared_AndADependenceHasNone()
    {
        BassModel model = Model(6, 3);
        var rng = new Xoshiro256PlusPlus(3703);
        double worstDet = 0.0, worstSine = 0.0;
        for (int k = 0; k < 400; k++)
        {
            int m = 1 + k % 6;
            ProjectionDirection[] d = [.. Enumerable.Range(0, m).Select(_ => model.Prior.DrawTerm(rng)!.Directional!.Direction)];
            double engine = DirectionDeterminant.LogDeterminant(d), hand = HandDeterminant(Gram(d, 6));
            if (hand < 1e-10) continue;   // a draw near a dependence; its own cases follow
            worstDet = Math.Max(worstDet, Math.Abs(engine - Math.Log(hand)));
            // A birth's factor det(Θ′ᵀΘ′)/det(ΘᵀΘ), the new direction's sin² to the others' span.
            double factor = Math.Exp(engine - DirectionDeterminant.LogDeterminant(d[..^1]));
            worstSine = Math.Max(worstSine, Math.Abs(factor - SineSquared(d[^1], d[..^1], 6)) / Math.Max(1e-3, factor));
        }
        output.WriteLine($"ln det against elimination: worst |Δ| {worstDet:E2}; a birth's factor against sin² by Gram–Schmidt: worst {worstSine:E2}");
        Assert.True(worstDet <= 1e-10, $"{worstDet:E2}");
        Assert.True(worstSine <= 1e-10, $"{worstSine:E2}");

        // No mass on a dependence: a sign pair, a repeat (its pivot rounding within 10⁻¹³ of the diagonal), three directions in a 2-plane, M past p.
        ProjectionDirection e1 = new([1], [1.0]), minusE1 = new([1], [-1.0]), tilted = new([0, 1], [0.6, 0.8]);
        Assert.Equal(double.NegativeInfinity, DirectionDeterminant.LogDeterminant([e1, minusE1]));
        Assert.Equal(double.NegativeInfinity, DirectionDeterminant.LogDeterminant([tilted, tilted]));
        Assert.Equal(double.NegativeInfinity, DirectionDeterminant.LogDeterminant([tilted, new([0, 1], [0.8, -0.6]), new([0, 1], [Math.Sqrt(0.5), Math.Sqrt(0.5)])]));
        BassTerm Term(ProjectionDirection d) => new(directional: new DirectionalHinge(d, 0.0));
        Assert.Equal(double.NegativeInfinity, DirectionDeterminant.LogCoupling([Term(e1), Term(new([0], [1.0])), Term(tilted)], 2));
        // Directions on disjoint coordinates do not repel: det 1, and the coupling is −ln c_M alone.
        ProjectionDirection e0 = new([0], [1.0]), across = new([2, 3], [0.6, -0.8]);
        Assert.Equal(0.0, DirectionDeterminant.LogDeterminant([e0, e1, across]));
        Assert.Equal(-DirectionDeterminant.LogNormalizer(6, 3), DirectionDeterminant.LogCoupling([Term(e0), Term(e1), Term(across)], 6));
        Assert.Equal(0.0, DirectionDeterminant.LogCoupling([], 6));
    }

    // ---- the option off -------------------------------------------------------------------------------------------------------------

    [Fact]
    public void Independently_TheTermsPriorIsTheSumAsBefore_AndTheConfigCarriesTheOption()
    {
        BassModel model = Model(6, 3);
        var rng = new Xoshiro256PlusPlus(3704);
        BassTerm[] terms = [.. Enumerable.Range(0, 5).Select(_ => model.Prior.DrawTerm(rng)!)];
        double sum = 0.0;
        foreach (BassTerm t in terms) sum += model.Prior.LogPrior(t);
        Assert.Equal(DirectionPrior.Independent, model.Prior.Directions);
        Assert.Equal(BitConverter.DoubleToInt64Bits(sum), BitConverter.DoubleToInt64Bits(model.Prior.LogTermsPrior(terms)));
        Assert.Equal(0.0, model.Prior.LogCoupling(terms));
        Assert.Equal(new DirectionalBasis(), new DirectionalBasis { Directions = DirectionPrior.Independent });

        // Under the option the terms' prior gains the coupling; the config round-trips with it and refuses a value it does not know.
        BassModel repulsive = Model(6, 3, DirectionPrior.Determinant);
        Assert.Equal(sum + DirectionDeterminant.LogCoupling(terms, 6), repulsive.Prior.LogTermsPrior(terms));
        BassConfig config = Determinant();
        Assert.Equal(DirectionPrior.Determinant, BassConfig.FromJson(config.ToJson()).Structure.Directional!.Directions);
        Assert.Throws<ArgumentException>(() => (config with { Structure = config.Structure with { Directional = new DirectionalBasis { Directions = (DirectionPrior)7 } } }).Validate());
    }

    // ---- each move's ratio ----------------------------------------------------------------------------------------------------------

    [Fact]
    public void EachMovesRatio_UnderTheDeterminant_MatchesADirectEvaluation()
    {
        BassData data = DirectionalChainTests.Data();
        BassModel model = BassModel.Create(Determinant(), data);
        BassStructurePrior prior = model.Prior;
        DirectionalTermSpace space = prior.Directional!;
        IBassProposalKernel kernel = model.Kernel;
        BassCountPrior count = model.Config.Count;
        const int p = 4;
        double Target(BassTerm[] s)
        {
            BassEvaluation e = DirectionalChainTests.Evaluate(model, s);
            return e.LogPrior + e.LogMarginal;
        }
        ProjectionDirection[] Of(IEnumerable<BassTerm> s) => [.. s.Select(t => t.Directional!.Direction)];

        ProjectionDirection[] directions = DirectionalChainTests.Directions;
        BassTerm[] s = [DirectionalChainTests.Within(space, directions[2], 0.7), DirectionalChainTests.Within(space, directions[1], 0.4)];
        var counts = new BassTermCounts(prior);
        foreach (BassTerm t in s) counts.Add(t);
        double worst = 0.0;
        int dependent = 0;
        foreach (ProjectionDirection d in directions)
        {
            // Slice 2's birth ratio with the determinant's: the new direction's sin² to the others' span, over c_(M+1)/c_M = (p − M)/p.
            BassTerm born = DirectionalChainTests.Within(space, d, 0.55);
            BassTerm[] grown = [.. s, born];
            int a = d.Count, m = s.Length;
            double engine = Target(grown) - Target(s) + kernel.LogBirthRatio(born, counts);
            double sine = SineSquared(d, Of(s), p);
            if (sine < 1e-12)
            {
                // A repeat of a current direction ((0.6, 0.8) on x₀, x₁ and −e₂ here): no prior mass.
                Assert.Equal(double.NegativeInfinity, engine);
                dependent++;
                output.WriteLine($"birth of {born}: a dependence, engine {engine}");
                continue;
            }
            double direct = DirectionalChainTests.Dense(model, data, grown) - DirectionalChainTests.Dense(model, data, s)
                          + Math.Log((m + count.Shape) / ((m + 1.0) * (count.Rate + 1.0)))
                          - Math.Log(3.0) - Math.Log(new[] { 1, 4, 6, 4 }[a]) - DirectionalChainTests.HandLogOrderAndVariables(counts, [.. d.Variables], 0)
                          + Math.Log(sine) - Math.Log((p - m) / (double)p);
            BassTermCounts after = counts.Clone();
            after.Add(born);
            double death = Target(s) - Target(grown) + kernel.LogDeathRatio(born, after);
            worst = Math.Max(worst, Math.Max(Math.Abs(engine - direct), Math.Abs(death + direct)));
            output.WriteLine($"birth of {born}: engine {engine:F12}, direct {direct:F12}; death {death:F12}; sin² {sine:F6}");
        }

        // The change: slice 2's collapsed likelihood ratio with the determinant's ratio, det of the changed directions over the current.
        var rng = new Xoshiro256PlusPlus(3705);
        foreach (int index in new[] { 0, 1 })
            for (int k = 0; k < 4; k++)
            {
                BassTermMove move = kernel.ProposeChange(rng, s[index], counts);
                BassTerm moved = move.Term!;
                BassTerm[] changed = [.. s];
                changed[index] = moved;
                double engine = Target(changed) - Target(s) + move.LogProposalRatio;
                double direct = DirectionalChainTests.Dense(model, data, changed) - DirectionalChainTests.Dense(model, data, s)
                              + Math.Log(HandDeterminant(Gram(Of(changed), p))) - Math.Log(HandDeterminant(Gram(Of(s), p)));
                worst = Math.Max(worst, Math.Abs(engine - direct));
                output.WriteLine($"change of {s[index]} to {moved}: engine {engine:F12}, direct {direct:F12}");
            }
        output.WriteLine($"worst |Δ| over the moves {worst:E2}; births on a dependence {dependent}");
        Assert.Equal(2, dependent);
        Assert.True(worst <= 1e-9, $"{worst:E2}");
    }

    // ---- the tilt on M --------------------------------------------------------------------------------------------------------------

    [Fact]
    public void TheTiltOnM_IsTheCountLawCutAtP_AsResolvedAndAsThePriorSamplerDraws()
    {
        // A count law with λ near 4, so M = p = 4 binds: p(M) ∝ Γ(h₁ + M)/(M!(h₂ + 1)^M) at h₁ = 40, h₂ = 10, declared cap 8.
        var count = new BassCountPrior(40.0, 10.0);
        BassData data = DirectionalChainTests.Data(80);
        BassModel model = BassModel.Create(Determinant(maxBasis: 8, count), data);
        const int p = 4, cap = 8;
        double[] law = [.. Enumerable.Range(0, cap + 1).Select(m => Math.Exp(SpecialFunctions.LogGamma(40.0 + m) - SpecialFunctions.LogGamma(m + 1.0) - m * Math.Log(11.0)))];
        double total = law.Sum(), past = law.Skip(p + 1).Sum() / total;

        // The resolution: M's cap at p, c_M stated, and the count law's mass past p reported.
        Assert.Equal(p, model.MaxBasis);
        Assert.Contains(model.Resolved.AutoPicks, a => a.Setting == "MaxBasis" && a.Reason.Contains("no mass past M = p"));
        BassAutoPick directions = Assert.Single(model.Resolved.AutoPicks, a => a.Setting == "Directions");
        output.WriteLine($"{directions.Setting}: {directions.Value}, {directions.Reason}");
        Assert.Contains($"{past:G4}", directions.Reason);
        foreach (double c in new[] { 0.75, 0.375, 0.09375 }) Assert.Contains($"{c:G4}", directions.Reason);

        // The prior sampler's M against the count law cut at p, the tilt's definition: F holds every drawn collection on this fixture, so the
        // law on M is p(M)·1{M ≤ p} renormalized, the directions' det(ΘᵀΘ)/c_M integrating to one at each M.
        var sampler = new BassPriorSampler(model);
        var rng = new Xoshiro256PlusPlus(3706);
        const int draws = 40_000;
        var seen = new int[p + 1];
        for (int k = 0; k < draws; k++) seen[sampler.DrawState(rng).State.Count]++;
        double kept = law.Take(p + 1).Sum(), chi = 0.0;
        for (int m = 0; m <= p; m++)
        {
            double expected = draws * law[m] / kept;
            chi += (seen[m] - expected) * (seen[m] - expected) / expected;
            output.WriteLine($"M = {m}: drawn {seen[m],6}, the count law cut at p {expected,9:F1}; uncut {draws * law[m] / total,9:F1}; c_M {Math.Exp(DirectionDeterminant.LogNormalizer(p, m)):F5}");
        }
        output.WriteLine($"χ² {chi:F2} on 4 degrees of freedom; the count law's mass past p {past:F4}; determinant rejections {sampler.DeterminantRejections} "
                         + $"of {sampler.Attempts} attempts; refusals {sampler.Refusals.Total}");
        Assert.Equal(0, sampler.Refusals.Total);
        Assert.True(sampler.DeterminantRejections > 0);
        Assert.True(chi < 20.0, $"χ² {chi:F2}");   // 4 degrees of freedom: p ≈ 5·10⁻⁴
    }

    // ---- Geweke's test with the prior alone -----------------------------------------------------------------------------------------

    [Theory]
    [InlineData(3707)]
    [InlineData(3717)]
    public void Geweke_UnderTheDeterminant_ThePriorChainLeavesThePriorInvariant(int seed)
    {
        BassData data = DirectionalChainTests.Data(80);
        BassModel model = BassModel.Create(Determinant(maxBasis: 8, new BassCountPrior(40.0, 10.0)), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        string[] statistics = [.. DirectionalChainTests.Statistics, "det(ΘᵀΘ)", "M = 4"];
        double[] Functions(BassState state)
        {
            double[] f = DirectionalChainTests.Functions(state, space);
            double det = Math.Exp(DirectionDeterminant.LogDeterminant([.. state.Terms.Select(t => t.Directional!.Direction)]));
            return [.. f, det, state.Count == 4 ? 1.0 : 0.0];
        }
        var rng = new Xoshiro256PlusPlus(seed);
        var sampler = new BassPriorSampler(model);
        const int draws = 40_000, burn = 5_000, steps = 200_000, batches = 40;
        var direct = new double[draws][];
        for (int k = 0; k < draws; k++) direct[k] = Functions(sampler.DrawState(rng).State);
        var chain = new BassChain(model, [], Math.Log(model.Target.HeldTau!.Value), new Xoshiro256PlusPlus(seed + 1), new Xoshiro256PlusPlus(seed + 2));
        chain.Moves.SetTemperature(new Temperature(0.0, TemperingPath.Likelihood));
        for (int t = 0; t < burn; t++) chain.Step();
        var successive = new double[steps][];
        for (int t = 0; t < steps; t++) successive[t] = Functions(chain.Step());
        double worst = 0.0;
        for (int g = 0; g < statistics.Length; g++)
        {
            double[] a = [.. direct.Select(r => r[g])], b = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1);
            double meanB = b.Average(), se = Math.Sqrt(varA / a.Length + Math.Pow(DirectionalChainTests.BatchError(b, batches), 2));
            double z = (meanA - meanB) / se;
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{statistics[g],-16} prior {meanA,9:F5}  chain {meanB,9:F5}  z {z,6:F2}");
        }
        output.WriteLine($"worst |z| {worst:F2}; moves {string.Join(", ", chain.Structure.Kinds.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Kind} {e.Accepted}/{e.Attempts}"))}; "
                         + $"prior sampler: acceptance {sampler.AcceptanceRate:F3}, determinant rejections {sampler.DeterminantRejections}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- a fit --------------------------------------------------------------------------------------------------------------------

    // Under the chain and SMC (whose particles start from the prior sampler), every draw's directions are independent and number p at most.
    // The option's regime is ridge blocks (K = 4), one of which follows a single index. With one-column hinges the single index needs several
    // hinges along nearly one direction, which the determinant taxes (the Grassmann discussion's reading), so that case is read at a looser bound.
    [Theory]
    [InlineData("chain", 4)]
    [InlineData("smc", 4)]
    [InlineData("chain", 1)]
    public void AFit_UnderTheDeterminant_DrawsIndependentDirections_AndFollowsTheFunction(string sampler, int columns)
    {
        ProjectionPursuitReference reference = ProjectionPursuitReferences.SingleIndex();
        var rng = new Xoshiro256PlusPlus(3710);
        const int n = 300;
        var x = new double[4][];
        for (int j = 0; j < 4; j++) x[j] = new double[n];
        var y = new double[n];
        var truth = new double[n];
        var point = new double[4];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < 4; j++) point[j] = x[j][i] = rng.NextDouble();
            truth[i] = reference.Function.Value(point);
            y[i] = truth[i] + 0.1 * Normal.SampleStandard(rng);
        }
        var data = new BassData([.. Enumerable.Range(0, 4).Select(j => new ContinuousPredictor($"x{j}", x[j], new InputDomain(0.0, 1.0)))], new ScalarResponse(y));
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 3, Directional = new DirectionalBasis { Directions = DirectionPrior.Determinant, Columns = columns } },
            Sampling = new BassSampling { Chains = 2, MasterSeed = 3711, BurnIn = 4000, Samples = 400, Thin = 2 },
        };
        if (sampler == "smc")
            config = config with { Sampling = config.Sampling with { Smc = new SmcSpec { Particles = 200, Replicates = 2, MutationSteps = 12, TargetConditionalEss = 0.9 } } };
        BassFit fit = Bass.Fit(config, data);
        Assert.Equal(4, BassModel.Create(config, data).MaxBasis);   // the default cap, min(1000, N − 1), cut at p
        int most = 0;
        foreach (BassDraw draw in fit.Draws.SelectMany(c => c))
        {
            // A ridge block's directions are its first column's; a hinge's its own.
            ProjectionDirection[] d = [.. draw.Structure.TermIds.Select(id => fit.Terms[id]).Where(t => t.Ridge is not { Column: > 1 })
                                                           .Select(t => t.Directional?.Direction ?? t.Ridge!.Direction)];
            most = Math.Max(most, d.Length);
            Assert.True(d.Length <= 4);
            Assert.True(double.IsFinite(DirectionDeterminant.LogDeterminant(d)));
        }
        double[,,] latent = fit.Predict(data, PredictionCases.Of(data)).Latent;
        double rmse = Math.Sqrt(Enumerable.Range(0, n).Average(i => Enumerable.Range(0, latent.GetLength(0)).Average(d => latent[d, i, 0]) is var f ? (f - truth[i]) * (f - truth[i]) : 0.0));
        double spread = Math.Sqrt(truth.Average(t => (t - truth.Average()) * (t - truth.Average())));
        output.WriteLine($"{sampler}, K = {columns}: surface against the function: root mean square {rmse:F4} of a spread {spread:F3}; most directions in a draw {most}");
        Assert.True(rmse <= (columns == 1 ? 0.10 : 0.05) * spread, $"{rmse:F4}");
    }
}
