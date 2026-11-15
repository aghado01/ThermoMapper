using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Maths.Rng;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E36's slice 3 gates on the ridge block as one structural unit (the BPPR brief): a ridge basis (<see cref="DirectionalBasis.Columns"/> = K)
/// holds each structure term (θ, t₀) as its block of K columns (<see cref="RidgeExpansion"/>), so the chain, the prior and the proposals are
/// slice 2's and the fit reads columns:
/// <list type="bullet">
/// <item>the design holds each term as its block, in order, and reads the terms back; a draw's columns load as they are;</item>
/// <item>every block proposal, birth, death and change, is the proposed structure's strict load, bit for bit, before and after acceptance;</item>
/// <item>the count prior reads structure terms, not columns;</item>
/// <item>each move's ratio against a direct evaluation: the dense marginal over every block's columns by the paper's formulas, the count
/// law's and the birth's π/q by hand;</item>
/// <item>Geweke's test with the prior alone, the chain at b = 0 against the prior sampler;</item>
/// <item>the config's and the resolution's refusals;</item>
/// <item>a fit end to end: whole blocks in every draw, the surface against the function, and the direction through E8.</item>
/// </list>
/// </summary>
public sealed class DirectionalRidgeTests(ITestOutputHelper output)
{
    private const int K = 4;

    private static readonly ProjectionDirection[] Directions =
    [
        new([1], [1.0]),
        new([0, 1], [0.6, 0.8]),
        new([1, 3], [-0.8, 0.6]),
        new([0, 2, 3], [0.48, -0.6, 0.64]),
    ];

    private static BassConfig Config(int maxBasis = 6, int columns = K)
    {
        BassConfig c = DirectionalChainTests.Config(maxBasis);
        return c with { Structure = c.Structure with { Directional = new DirectionalBasis { Columns = columns } } };
    }

    private static BassTerm Unit(ProjectionDirection d, double knot) => new(directional: new DirectionalHinge(d, knot));

    // A structure term with its knot at fraction `at` of its range from L.
    private static BassTerm Within(DirectionalTermSpace space, ProjectionDirection d, double at)
    {
        KnotRange r = space.KnotRangeOf(d)!.Value;
        return Unit(d, r.Lower + at * r.Width);
    }

    // Every block's columns by the paper's formulas at the projected training cases, z by the compensated standardization, the knots
    // the expansion's (checked against R's quantiles by hand in RidgeColumnTests).
    private static double[][] HandColumns(BassModel model, BassData data, IEnumerable<BassTerm> units)
    {
        double[][] raw = DirectionalHingeTests.Raw(data);
        (double Mean, double Deviation)[] m = DirectionalHingeTests.Moments(data);
        var columns = new List<double[]>();
        foreach (BassTerm unit in units)
        {
            DirectionalHinge h = unit.Directional!;
            double[] knots = model.Prior.Directional!.RidgeKnots(h, K)!;
            double[] u = [.. Enumerable.Range(0, raw[0].Length).Select(i =>
                Enumerable.Range(0, h.Direction.Count).Sum(k => h.Direction.Weights[k] * (raw[h.Direction.Variables[k]][i] - m[h.Direction.Variables[k]].Mean) / m[h.Direction.Variables[k]].Deviation))];
            for (int l = 1; l <= K; l++) columns.Add([.. u.Select(x => RidgeColumnTests.Paper(knots, l, x).Value)]);
        }
        return [.. columns];
    }

    private static double Dense(BassModel model, BassData data, IEnumerable<BassTerm> units)
        => DirectionalChainTests.DenseLogMarginal(((ScalarResponse)data.Response).Values.ToArray(), HandColumns(model, data, units), model.Target.HeldTau!.Value, 3.0, 2.0);

    private static BassEvaluation Evaluate(BassModel model, BassTerm[] units)
    {
        var chain = new BassChain(model, units, Math.Log(model.Target.HeldTau!.Value), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        return model.Target.Evaluate(chain.Current);
    }

    // ---- the design ----------------------------------------------------------------------------------------------------------------

    [Fact]
    public void TheDesign_HoldsEachTermAsItsBlock_AndReadsTheTermsBack()
    {
        BassData data = DirectionalChainTests.Data();
        BassModel model = BassModel.Create(Config(), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        Assert.NotNull(model.Expansion);
        BassTerm[] units = [Within(space, Directions[1], 0.3), Within(space, Directions[3], 0.6)];
        BassDesign design = model.NewDesign();
        Assert.Null(design.Load(units));
        Assert.True(design.Expands);
        Assert.Equal(units, design.Terms);
        Assert.Equal(2, design.Units);
        Assert.Equal(2 * K, design.Count);
        for (int b = 0; b < units.Length; b++)
        {
            DirectionalHinge h = units[b].Directional!;
            double[] knots = space.RidgeKnots(h, K)!;
            for (int l = 1; l <= K; l++) Assert.Equal(new BassTerm(ridge: new RidgeColumn(h.Direction, knots, l)), design.ColumnTerms[b * K + l - 1]);
        }

        // A draw's terms are its columns: they load as they are, to the same design.
        BassDesign columns = model.NewDesign();
        Assert.Null(columns.Load(design.ColumnTerms));
        Assert.Equal(units, columns.Terms);
        Assert.Equal(design.ResidualSquares, columns.ResidualSquares);
        Assert.Equal(design.FittedSquares, columns.FittedSquares);

        // The count prior reads the structure's terms, not its columns.
        BassEvaluation e = Evaluate(model, units);
        double expected = model.Prior.LogCollectionWeight(2) + units.Sum(model.Prior.LogPrior);
        output.WriteLine($"log prior {e.LogPrior:R} against M = 2's {expected:R}; M = 8's would be {model.Prior.LogCollectionWeight(8) + units.Sum(model.Prior.LogPrior):R}");
        Assert.Equal(expected, e.LogPrior, 12);
    }

    [Fact]
    public void BlockProposals_AreTheProposedStructuresStrictLoads_BitForBit()
    {
        BassData data = DirectionalChainTests.Data();
        BassModel model = BassModel.Create(Config(), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        BassTerm[] units = [Within(space, Directions[1], 0.3), Within(space, Directions[3], 0.6)];
        BassTerm born = Within(space, Directions[2], 0.5), other = Within(space, Directions[0], 0.4);
        BassDesign design = model.NewDesign();
        Assert.Null(design.Load(units));

        void Same(DesignProposal proposal, BassTerm[] proposed, bool death = false)
        {
            BassDesign fresh = model.NewDesign();
            Assert.Null(fresh.Load(proposed));
            Assert.True(proposal.Admitted);
            Assert.Equal(fresh.Count, proposal.Columns);
            Assert.Equal(fresh.ResidualSquares, proposal.ResidualSquares);
            Assert.Equal(fresh.FittedSquares, proposal.FittedSquares);
            Assert.Equal(fresh.MinAngle, proposal.MinAngle);
            // The newest angle: the least over the new block's columns, the last K; none for a death.
            Assert.Equal(death ? double.NaN : Enumerable.Range(fresh.Count - K, K).Min(fresh.Angle), proposal.NewestAngle);
            design.Accept();
            Assert.Equal(proposed, design.Terms);
            Assert.Equal(fresh.ColumnTerms, design.ColumnTerms);
            Assert.Equal(fresh.ResidualSquares, design.ResidualSquares);
            Assert.Equal(fresh.FittedSquares, design.FittedSquares);
        }

        Same(design.ProposeBirth(born), [.. units, born]);
        Same(design.ProposeChange(0, other), [units[1], born, other]);
        Same(design.ProposeDeath(1), [units[1], other], death: true);
        Assert.Equal(BassRefusal.Numerical, design.ProposeBirth(other).Refusal);   // a repeated block's columns repeat
    }

    // ---- each move's ratio ---------------------------------------------------------------------------------------------------------

    [Fact]
    public void EachMovesRatio_MatchesADirectEvaluation()
    {
        BassData data = DirectionalChainTests.Data();
        BassModel model = BassModel.Create(Config(), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        IBassProposalKernel kernel = model.Kernel;
        BassCountPrior count = model.Config.Count;
        double Target(BassTerm[] s)
        {
            BassEvaluation e = Evaluate(model, s);
            return e.LogPrior + e.LogMarginal;
        }

        // The collapsed law given frozen blocks against the dense marginal over their columns.
        BassTerm[] s = [Within(space, Directions[1], 0.7), Within(space, Directions[0], 0.4)];
        double worstMarginal = 0.0;
        for (int m = 0; m <= s.Length; m++)
        {
            double dense = Dense(model, data, s[..m]), engine = Evaluate(model, s[..m]).LogMarginal;
            worstMarginal = Math.Max(worstMarginal, Math.Abs(engine - dense) / Math.Max(1.0, Math.Abs(dense)));
        }
        output.WriteLine($"the collapsed law against the dense marginal over the blocks' columns: worst relative {worstMarginal:E2}");
        Assert.True(worstMarginal <= 1e-11, $"{worstMarginal:E2}");

        var counts = new BassTermCounts(model.Prior);
        foreach (BassTerm t in s) counts.Add(t);
        double worst = 0.0;
        int[] choose = [1, 4, 6, 4];
        foreach (ProjectionDirection d in Directions.Skip(2))
        {
            BassTerm born = Within(space, d, 0.55);
            BassTerm[] grown = [.. s, born];
            int a = d.Count, m = s.Length;
            double engine = Target(grown) - Target(s) + kernel.LogBirthRatio(born, counts);
            double lq = Math.Log(5.0 + counts.Order(BassPart.Directional, a)) - Math.Log(15.0 + m);
            if (a == 1) lq -= Math.Log(4.0);
            else
            {
                double[] z = [.. Enumerable.Range(0, 4).Select(u => 5.0 + counts.Appearances(BassPart.Directional, u))];
                lq += Wallenius.LogProbability(z, [.. d.Variables]);
            }
            double direct = Dense(model, data, grown) - Dense(model, data, s) + Math.Log((m + count.Shape) / ((m + 1.0) * (count.Rate + 1.0)))
                          - Math.Log(3.0) - Math.Log(choose[a]) - lq;
            BassTermCounts after = counts.Clone();
            after.Add(born);
            double death = Target(s) - Target(grown) + kernel.LogDeathRatio(born, after);
            worst = Math.Max(worst, Math.Max(Math.Abs(engine - direct), Math.Abs(death + direct)));
            output.WriteLine($"birth of {born}: engine {engine:F12}, direct {direct:F12}; death {death:F12}");
        }
        var rng = new Xoshiro256PlusPlus(3651);
        foreach (int index in new[] { 0, 1 })
            for (int k = 0; k < 4; k++)
            {
                BassTermMove move = kernel.ProposeChange(rng, s[index], counts);
                if (move.Term is not { } moved) continue;
                BassTerm[] changed = [.. s];
                changed[index] = moved;
                double engine = Target(changed) - Target(s) + move.LogProposalRatio, direct = Dense(model, data, changed) - Dense(model, data, s);
                worst = Math.Max(worst, Math.Abs(engine - direct));
                output.WriteLine($"change of {s[index]} to {moved}: engine {engine:F12}, direct {direct:F12}");
            }
        output.WriteLine($"worst |Δ| over the moves {worst:E2}");
        Assert.True(worst <= 1e-9, $"{worst:E2}");
    }

    // ---- Geweke's test with the prior alone ----------------------------------------------------------------------------------------

    private static readonly string[] Statistics = ["M", "a = 1", "a = 2", "a = 3", "uses x0", "uses x1", "uses x2", "uses x3", "θ₁ > 0", "knot position"];

    private static double[] Functions(BassState state, DirectionalTermSpace space)
    {
        var f = new double[Statistics.Length];
        f[0] = state.Count;
        foreach (BassTerm t in state.Terms)
        {
            DirectionalHinge h = t.Directional!;
            f[h.Direction.Count]++;
            foreach (int v in h.Direction.Variables) f[4 + v]++;
            if (h.Direction.Weights[0] > 0.0) f[8]++;
            KnotRange r = space.KnotRangeOf(h.Direction)!.Value;
            f[9] += (r.Upper - h.Knot) / r.Width;
        }
        return f;
    }

    private static double BatchError(double[] x, int batches)
    {
        int size = x.Length / batches;
        double[] means = [.. Enumerable.Range(0, batches).Select(b => x.Skip(b * size).Take(size).Average())];
        double grand = means.Average();
        return Math.Sqrt(means.Sum(v => (v - grand) * (v - grand)) / (batches - 1) / batches);
    }

    [Fact]
    public void Geweke_ThePriorChain_LeavesThePriorInvariant()
    {
        BassData data = DirectionalChainTests.Data(80);
        BassModel model = BassModel.Create(Config(maxBasis: 3), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        var rng = new Xoshiro256PlusPlus(3652);
        var sampler = new BassPriorSampler(model);
        const int draws = 30_000, burn = 3_000, steps = 120_000, batches = 40;
        var direct = new double[draws][];
        for (int k = 0; k < draws; k++) direct[k] = Functions(sampler.DrawState(rng).State, space);
        var chain = new BassChain(model, [], Math.Log(model.Target.HeldTau!.Value), new Xoshiro256PlusPlus(3653), new Xoshiro256PlusPlus(3654));
        chain.Moves.SetTemperature(new Temperature(0.0, TemperingPath.Likelihood));
        for (int t = 0; t < burn; t++) chain.Step();
        var successive = new double[steps][];
        for (int t = 0; t < steps; t++) successive[t] = Functions(chain.Step(), space);
        double worst = 0.0;
        for (int g = 0; g < Statistics.Length; g++)
        {
            double[] a = [.. direct.Select(r => r[g])], b = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1);
            double meanB = b.Average(), se = Math.Sqrt(varA / a.Length + Math.Pow(BatchError(b, batches), 2));
            double z = (meanA - meanB) / se;
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{Statistics[g],-14} prior {meanA,9:F5}  chain {meanB,9:F5}  z {z,6:F2}");
        }
        output.WriteLine($"worst |z| {worst:F2}; moves {string.Join(", ", chain.Structure.Kinds.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Kind} {e.Accepted}/{e.Attempts}"))}; "
                         + $"prior sampler acceptance {sampler.AcceptanceRate:F3}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");

        // The order counts are structure terms': each birth's and change's candidate counted once, at its active count, none at order 0.
        BassSupportTally support = chain.Structure.Support;
        long counted = Enumerable.Range(1, 3).Sum(a => support.Attempts(BassPart.Directional, a));
        BassKindCount births = chain.Structure.Kinds[BassMoveKind.Birth], changes = chain.Structure.Kinds[BassMoveKind.DirectionChange];
        Assert.Equal(0, support.Attempts(BassPart.Directional, 0));
        Assert.Equal(births.Attempts - births.Null + changes.Attempts - changes.Null, counted);
    }

    // ---- coincident knots and admission all or none -------------------------------------------------------------------------------

    [Fact]
    public void CoincidentKnots_PutABlockOutsideTheSupport_AndARefusedBlockLeavesTheDesignAsItWas()
    {
        // x₀ takes 0.5 on every third case, so on +e₀ the knots at the quartiles of the cases above a low t₀ fall in that tie.
        const int n = 120;
        var rng = new Xoshiro256PlusPlus(3656);
        var x = new double[4][];
        for (int j = 0; j < 4; j++) x[j] = [.. Enumerable.Range(0, n).Select(i => j == 0 && i % 3 == 0 ? 0.5 : rng.NextDouble())];
        var y = new double[n];
        for (int i = 0; i < n; i++) y[i] = Math.Exp(0.6 * x[0][i] + 0.8 * x[1][i]) + 0.05 * Normal.SampleStandard(rng);
        var data = new BassData([.. Enumerable.Range(0, 4).Select(j => new ContinuousPredictor($"x{j}", x[j], new InputDomain(0.0, 1.0)))], new ScalarResponse(y));
        BassModel model = BassModel.Create(Config(), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        var along = new ProjectionDirection([0], [1.0]);
        BassTerm tied = Unit(along, space.KnotRangeOf(along)!.Value.Lower);
        double[] knots = space.RidgeKnots(tied.Directional!, K)!;
        output.WriteLine($"knots on +e₀ from its range's lower end: {string.Join(", ", knots.Select(t => t.ToString("G6")))}");
        Assert.Contains(Enumerable.Range(1, K - 1), j => knots[j] == knots[j + 1]);   // two interior knots coincide

        // Two coincident interior knots make two columns equal, which D06's gate refuses: the block lies outside F, all of it.
        Assert.Equal(BassRefusal.Numerical, model.NewDesign().Load([tied]));
        // Coincident upper knots leave d_K undefined: no block, outside the support, as a range with no width is.
        BassModel tiedTop = BassModel.Create(Config(), DirectionalChainTests.Data(tied: true));
        var down = new ProjectionDirection([3], [-1.0]);
        BassTerm topTie = Unit(down, tiedTop.Prior.Directional!.KnotRangeOf(down)!.Value.Lower);
        Assert.Null(tiedTop.Expansion!.Expand(topTie));
        Assert.Equal(BassRefusal.Support, tiedTop.NewDesign().Load([topTie]));

        // A birth or a change whose block is refused adds none of it: the design's terms, columns and fit stay as they were.
        BassTerm[] units = [Within(space, Directions[1], 0.3), Within(space, Directions[3], 0.6)];
        BassDesign design = model.NewDesign();
        Assert.Null(design.Load(units));
        BassTerm[] columns = [.. design.ColumnTerms];
        double rss = design.ResidualSquares, fitted = design.FittedSquares;
        Assert.Equal(BassRefusal.Numerical, design.ProposeBirth(tied).Refusal);
        Assert.Equal(BassRefusal.Numerical, design.ProposeChange(0, tied).Refusal);
        Assert.Equal(units, design.Terms);
        Assert.Equal(columns, design.ColumnTerms);
        Assert.Equal(rss, design.ResidualSquares);
        Assert.Equal(fitted, design.FittedSquares);
        Assert.Throws<InvalidOperationException>(() => design.Accept());   // nothing admitted is pending
    }

    // ---- the degree, and K = 1 ------------------------------------------------------------------------------------------------------

    [Fact]
    public void AtOneColumn_TheBasisIsSlice2sHinge_AtTheFitsDegree()
    {
        // K = 1 is no expansion: each term is its own column, slice 2's directional hinge (θᵀz − t)₊^q/c at the fit's degree q.
        BassData data = DirectionalChainTests.Data();
        BassConfig hinge = DirectionalChainTests.Config(degree: 2);
        Assert.Equal(new DirectionalBasis(), new DirectionalBasis { Columns = 1 });
        BassModel model = BassModel.Create(hinge, data);
        Assert.Null(model.Expansion);
        Assert.Equal(1, model.Config.Structure.ColumnsPerTerm);
        DirectionalTermSpace space = model.Prior.Directional!;
        BassTerm[] units = [Within(space, Directions[1], 0.3), Within(space, Directions[3], 0.6)];
        BassDesign design = model.NewDesign();
        Assert.Null(design.Load(units));
        Assert.False(design.Expands);
        Assert.Equal(units, design.ColumnTerms);
        Assert.Equal(2, model.Basis.Degree);
        double worst = 0.0;
        foreach (BassTerm unit in units)
        {
            var column = new double[data.Response.Cases];
            model.Basis.CaseColumn(unit, model.Inputs.Cases, column);
            DirectionalHinge h = unit.Directional!;
            (double[] hand, _) = DirectionalHingeTests.Hand(data, DirectionalHingeTests.Raw(data), h.Direction, h.Knot, 2, HingeNormalization.Box);
            for (int i = 0; i < column.Length; i++) worst = Math.Max(worst, Math.Abs(column[i] - hand[i]) / Math.Max(1.0, Math.Abs(hand[i])));
        }
        output.WriteLine($"K = 1 at degree 2: worst |column − slice 1's hand hinge| {worst:E2}");
        Assert.True(worst <= 1e-13, $"{worst:E2}");
    }

    // ---- the refusals --------------------------------------------------------------------------------------------------------------

    [Fact]
    public void TheConfigAndTheResolution_RefuseWhatARidgeBasisDoesNotRead()
    {
        BassConfig good = Config();
        good.Validate();
        BassStructureSpec s = good.Structure;
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Directional = new DirectionalBasis { Columns = 0 } } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Degree = 2 } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { Support = new SupportRule(SupportKind.ColumnMass) } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Structure = s with { HingeNormalization = HingeNormalization.TrainingRange } }).Validate());
        Assert.Throws<ArgumentException>(() => (good with { Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Held(1.0)) }).Validate());
        (Config(columns: 1) with { Structure = s with { Directional = new DirectionalBasis(), Degree = 2 } }).Validate();   // the hinge alone keeps q

        BassData data = DirectionalChainTests.Data();
        BassModel model = BassModel.Create(good with { Structure = s with { MaxBasis = null } }, data);
        Assert.Equal((120 - 1) / K, model.MaxBasis);   // N − 1 columns at K a term
        BassModel few = BassModel.Create(good with { Structure = s with { Directional = new DirectionalBasis { Columns = K, ActiveCases = 2.0 } } }, data);
        Assert.Equal(1.0 - (K + 1.0) / 120.0, few.Resolved.KnotLevel);   // BayesPPR's floor: K + 1 cases above U
        Assert.Contains(few.Resolved.AutoPicks, p => p.Setting == "Directional" && p.Reason.Contains("K + 1"));
    }

    // ---- a fit end to end ----------------------------------------------------------------------------------------------------------

    // Under the chain, tempering (whose swaps hand a design over by its snapshot, which holds columns) and SMC (whose particles fork theirs).
    [Theory]
    [InlineData("chain")]
    [InlineData("tempering")]
    [InlineData("smc")]
    public void AFit_DrawsWholeBlocks_FollowsTheFunction_AndE8FindsItsDirection(string sampler)
    {
        ProjectionPursuitReference reference = ProjectionPursuitReferences.SingleIndex();
        var rng = new Xoshiro256PlusPlus(3655);
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
            Structure = new BassStructureSpec { MaxInteraction = 3, MaxBasis = 10, Directional = new DirectionalBasis { Columns = K } },
            Sampling = new BassSampling { Chains = 2, MasterSeed = 3655, BurnIn = 4000, Samples = 400, Thin = 2 },
        };
        config = sampler switch
        {
            "tempering" => config with { Sampling = config.Sampling with { Tempering = new TemperingSpec { InverseTemperatures = [1.0, 0.5, 0.25, 0.1] } } },
            "smc" => config with { Sampling = config.Sampling with { Smc = new SmcSpec { Particles = 200, Replicates = 2, MutationSteps = 12, TargetConditionalEss = 0.9 } } },
            _ => config,
        };
        BassFit fit = Bass.Fit(config, data);
        foreach (BassDraw draw in fit.Draws.SelectMany(c => c))
        {
            BassTerm[] terms = [.. draw.Structure.TermIds.Select(id => fit.Terms[id])];
            Assert.Equal(0, terms.Length % K);
            Assert.All(terms, t => Assert.NotNull(t.Ridge));
            Assert.Equal(terms.Length / K, terms.Count(t => t.Ridge!.Column == 1));
        }
        double[,,] latent = fit.Predict(data, PredictionCases.Of(data)).Latent;
        double rmse = Math.Sqrt(Enumerable.Range(0, n).Average(i => Enumerable.Range(0, latent.GetLength(0)).Average(d => latent[d, i, 0]) is var f ? (f - truth[i]) * (f - truth[i]) : 0.0));
        double spread = Math.Sqrt(truth.Average(t => (t - truth.Average()) * (t - truth.Average())));
        ActiveSubspaceAnalysis analysis = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Points = 1024 });
        double[,] vectors = analysis.Posterior.MeanCVectors;
        double[,] leading = new double[4, 1], target = new double[4, 1];
        for (int j = 0; j < 4; j++) { leading[j, 0] = vectors[j, 0]; target[j, 0] = reference.Directions[0][j]; }
        double angle = SubspacePosterior.PrincipalAngles(leading, target)[0];
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        double[] weights = fit.DrawWeights is { } w ? [.. w] : [.. Enumerable.Repeat(1.0 / pooled.Length, pooled.Length)];
        double meanM = pooled.Select((d, i) => weights[i] * d.Structure.TermIds.Length / K).Sum() / weights.Sum();
        // M's readouts count ridges, not columns: the posterior summary, and a chain's traces draw by draw (an SMC fit has none).
        Assert.Equal(meanM, fit.Summary().Posterior.Count.Mean, 10);
        if (sampler != "smc")
        {
            BassTraces traces = fit.Traces();
            for (int c = 0; c < fit.Draws.Count; c++)
                for (int i = 0; i < fit.Draws[c].Count; i++) Assert.Equal(fit.Draws[c][i].Structure.TermIds.Length / K, traces.Chains[c].Count[i]);
            Assert.True(traces.Chains.SelectMany(t => t.Count).Max() <= fit.Config.Structure.MaxBasis);
        }
        // An inclusion count: the candidates' tally by order counts each ridge a birth or a change proposes once, at its a (the sweep's
        // candidate adds the structure term), so it equals those moves' attempts with a term, not K times them.
        BassKindCount births = fit.Diagnostics.MoveKinds[BassMoveKind.Birth], changes = fit.Diagnostics.MoveKinds[BassMoveKind.DirectionChange];
        long tallied = fit.Diagnostics.Support.Entries.Where(e => e.Part == BassPart.Directional).Sum(e => e.Attempts);
        Assert.True(tallied > 0);
        Assert.Equal(births.Attempts - births.Null + changes.Attempts - changes.Null, tallied);
        Assert.All(fit.Diagnostics.Support.Entries.Where(e => e.Attempts > 0), e => Assert.InRange(e.Order, 1, 3));
        output.WriteLine($"{sampler}: surface against the function: root mean square {rmse:F4} of a spread {spread:F3}; E8's angle to the direction {angle:F4} rad; mean M {meanM:F2} blocks");
        Assert.True(rmse <= 0.05 * spread, $"{rmse:F4}");
        Assert.True(angle <= 0.05, $"{angle:F4}");
    }
}
