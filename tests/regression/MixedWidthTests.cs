using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MixedWidthTests(ITestOutputHelper output)
{
    internal static BassConfig Config(int max = 6) => CategoricalRidgeTests.Config(columns: 4, maxBasis: max) with
    {
        Structure = CategoricalRidgeTests.Config(columns: 4, maxBasis: max).Structure with
        { Directional = new DirectionalBasis { Columns = 4, MixedWidths = true } },
    };

    private static BassTerm At(BassModel model, ProjectionDirection d, double at, int width = 1, params CategoricalFactor[] cells)
    {
        KnotRange r = model.Prior.Directional!.KnotRangeOf(d, cells)!.Value;
        return new BassTerm(categorical: cells, directional: new DirectionalHinge(d, r.Lower + at * r.Width, width));
    }

    private static readonly ProjectionDirection Pair = new([0, 1], [0.6, 0.8]);
    private static readonly ProjectionDirection Axis = new([2], [1.0]);

    [Fact]
    public void RejectedMixedRemoval_IsAtomic_AndPurePriorsRejectMixedMarkers()
    {
        BassData data = CategoricalRidgeTests.Data();
        BassModel model = BassModel.Create(Config(), data);
        BassTerm held = At(model, Pair, .4), absent = At(model, Pair, .4, 4);
        var counts = new BassTermCounts(model.Prior); counts.Add(held);
        Assert.Throws<ArgumentException>(() => counts.Remove(absent));
        Assert.Equal(1, counts.Terms);
        Assert.Equal(1, counts.Order(BassPart.Directional, 2));
        Assert.Equal(1, counts.Appearances(BassPart.Directional, 0));
        Assert.Equal(1, counts.ForKind(1).Terms);
        Assert.Equal(0, counts.ForKind(0).Terms);
        Assert.True(double.IsFinite(model.Kernel.LogDeathRatio(held, counts)));
        counts.Remove(held);
        Assert.Equal(0, counts.Terms);
        Assert.Equal(0, counts.Order(BassPart.Directional, 2));
        BassModel pure = BassModel.Create(CategoricalRidgeTests.Config(columns: 4), data);
        Assert.Throws<ArgumentException>(() => pure.Prior.LogPrior(absent));
    }

    [Fact]
    public void Representation_ReadsWidthsPerTerm_AndRestoresMixedColumns()
    {
        BassModel model = BassModel.Create(Config(), CategoricalRidgeTests.Data());
        BassTerm[] units = [At(model, Pair, 0.4, 4), At(model, Axis, 0.7), new(categorical: [new CategoricalFactor(0, [0])])];
        BassDesign design = model.NewDesign();
        Assert.Null(design.Load(units));
        Assert.Equal(6, design.ColumnTerms.Count);
        Assert.Equal(units, design.Terms);
        BassTerm[] columns = [.. design.ColumnTerms];
        Assert.Equal(units, model.Expansion!.Units(columns));
        BassDesign restored = model.NewDesign();
        Assert.Null(restored.Load(columns));
        Assert.Equal(units, restored.Terms);
        Assert.NotEqual(units[0], new BassTerm(directional: new DirectionalHinge(Pair, units[0].Directional!.Knot)));
        Assert.Equal(3, new HashSet<BassTerm>(units).Count);
        var table = new BassTermTable();
        var structure = new BassStructure([.. columns.Select(table.Intern)]);
        Assert.Equal((1, 1, 1), structure.DirectionalTermCounts(table));
        Assert.Equal(3, structure.TermCount(table));
        Assert.Equal((0, 0, 0), new BassStructure([]).DirectionalTermCounts(table));
        int ordinary = table.Intern(new BassTerm(continuous: [new HingeFactor(0, 1, 0)]));
        Assert.Throws<ArgumentException>(() => new BassStructure([ordinary]).DirectionalTermCounts(table));
    }

    [Fact]
    public void StructuralAndColumnCaps_AreIndependent_IncludingAWidthUpMove()
    {
        BassData data = DirectionalChainTests.Data(12);
        BassConfig config = Config(max: 11) with { Structure = Config(11).Structure with { MaxInteraction = 1, MaxVarianceInflation = null, AngleFloor = 1e-12 } };
        BassModel model = BassModel.Create(config, data);
        Assert.Equal(11, model.MaxBasis);
        Assert.Equal(11, model.Admission.MaxBasis);
        var units = new List<BassTerm>();
        for (int v = 0; v < 3; v++)
            for (int k = 0; k < 3; k++) units.Add(At(model, new ProjectionDirection([v], [1.0]), 0.55 + 0.12 * k));
        BassDesign design = model.NewDesign();
        // Nine one-column terms satisfy both caps even though 4 * 9 exceeds N - 1.
        Assert.Null(design.Load(units));
        Assert.Equal(9, design.Terms.Count);
        BassTerm first = units[0];
        var up = new BassTerm(directional: new DirectionalHinge(first.Directional!.Direction, first.Directional.Knot, 4));
        Assert.Equal(BassRefusal.Cap, design.ProposeChange(0, up).Refusal);
        BassModel capped = BassModel.Create(Config(2), CategoricalRidgeTests.Data());
        Assert.Equal(BassRefusal.Cap, capped.NewDesign().Load([At(capped, Pair, 0.4), At(capped, Axis, 0.7), new(categorical: [new CategoricalFactor(0, [0])])]));
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void ConditionedKindsAndIntegratedShares_MatchHandLaws(bool categories)
    {
        BassData data = CategoricalRidgeTests.Data();
        if (!categories) data = CategoricalRidgeTests.Continuous(data);
        BassModel model = BassModel.Create(Config(), data);
        BassStructurePrior prior = model.Prior;
        int r = categories ? 3 : 2;
        Assert.Equal(r, prior.AvailableKinds);
        var units = new List<BassTerm> { At(model, Pair, 0.4, 4), At(model, Axis, 0.7), At(model, Pair, 0.6, 4) };
        if (categories) units.Add(new(categorical: [new CategoricalFactor(0, [0])]));
        // Counts (2,1[,1]): (r-1)! * 2! / (M+r-1)!, computed without gamma functions.
        double expected = 2.0;
        for (int j = r; j < units.Count + r; j++) expected /= j;
        Assert.Equal(Math.Log(expected), prior.LogKindShares(units), 12);
        BassTerm hinge = units[1];
        double hand = -Math.Log(2.0) - Math.Log(3.0) - Math.Log(2.0)
                    + prior.Directional!.LogKnot(hinge.Directional!, []);
        if (categories) hand -= Math.Log(2.0); // categorical order is uniform on 0..1.
        Assert.Equal(hand, prior.LogPrior(hinge), 12);
        if (categories) Assert.Equal(-Math.Log(2.0) - Math.Log(3.0), prior.LogPrior(units[^1]), 12);
        var seen = new int[r];
        var rng = new Xoshiro256PlusPlus(3811);
        for (int j = 0; j < 12_000; j++)
        {
            BassTerm term = prior.DrawTerm(rng)!;
            Assert.NotNull(term);
            seen[BassStructurePrior.KindOf(term)]++;
            if (term.Directional is null) Assert.NotEmpty(term.Categorical);
        }
        foreach (int count in seen) Assert.InRange(count, 12_000 / r - 200, 12_000 / r + 200);
        output.WriteLine($"available kinds {r}, conditional-law draws {string.Join('/', seen)}; log shares {prior.LogKindShares(units):R}");
    }

    [Fact]
    public void BirthDeathAndWidthChange_UseTheActualForwardAndReverseProbabilities()
    {
        BassModel model = BassModel.Create(Config(), CategoricalRidgeTests.Data());
        BassStructurePrior prior = model.Prior;
        BassTerm[] terms = [At(model, Pair, 0.4, 4), At(model, Axis, 0.7), new(categorical: [new CategoricalFactor(0, [0])])];
        var counts = new BassTermCounts(prior);
        foreach (BassTerm t in terms) counts.Add(t);
        double w = model.Config.Structure.OrderWeight;
        foreach (int width in new[] { 1, 4 })
        {
            BassTerm born = At(model, new ProjectionDirection([0], [-1.0]), 0.6, width, new CategoricalFactor(0, [1]));
            // One own-kind term: smooth a=2, hinge a=1; both have categorical order zero. Variable at a=1 is uniform.
            int ownOrderCount = width == 1 ? 1 : 0;
            double logQ = -Math.Log(3.0) + Math.Log(w + ownOrderCount) - Math.Log(2.0 * w + 1.0)
                - Math.Log(3.0) + prior.Directional!.LogDirection(1) + prior.Directional.LogKnot(born.Directional!, born.Categorical)
                + Math.Log(w) - Math.Log(2.0 * w + 1.0) - Math.Log(2.0) - Math.Log(3.0);
            Assert.Equal(logQ, model.Kernel.LogBirth(born, counts), 12);
            Assert.Equal(-Math.Log(4.0) - logQ, model.Kernel.LogBirthRatio(born, counts), 12);
            var after = counts.Clone(); after.Add(born);
            Assert.Equal(logQ + Math.Log(4.0), model.Kernel.LogDeathRatio(born, after), 12);
        }
        var rng = new Xoshiro256PlusPlus(3812);
        int flips = 0, signs = 0, changes = 0;
        for (int j = 0; j < 2000; j++)
        {
            BassTerm old = terms[j % 2];
            BassTermMove move = model.Kernel.ProposeChange(rng, old, counts);
            Assert.NotNull(move.Term);
            BassTerm moved = move.Term!;
            if (move.Kind == BassMoveKind.WidthChange)
            {
                flips++;
                Assert.Equal(0.0, move.LogProposalRatio);
                Assert.Equal(old.Directional!.Direction, moved.Directional!.Direction);
                Assert.Equal(old.Directional.Knot, moved.Directional.Knot);
                BassTerm[] next = [.. terms]; next[j % 2] = moved;
                // From (1,1,1) to (2,0,1) or (0,2,1): integrated shares ratio 2.
                Assert.Equal(Math.Log(2.0), prior.LogTermsPrior(next) - prior.LogTermsPrior(terms), 11);
            }
            else
            {
                changes++;
                Assert.Equal(old.Directional!.Columns, moved.Directional!.Columns);
                Assert.Equal(0.0, prior.LogPrior(moved) - prior.LogPrior(old) + move.LogProposalRatio, 11);
                if (old.Directional.Direction.Count == 1 && moved.Directional.Direction.Weights[0] < 0.0) signs++;
            }
        }
        Assert.InRange(flips, 900, 1100);
        Assert.InRange(signs, 180, 320);
        output.WriteLine($"width flips {flips}, direction changes {changes}, negative a=1 changes {signs}");
        // Categorical-only birth: kind 1/3, order k=1 surely, subset size and set as declared; no phantom direction factor.
        BassTerm categorical = new(categorical: [new CategoricalFactor(0, [1])]);
        Assert.Equal(-Math.Log(3.0) - Math.Log(2.0) - Math.Log(3.0), model.Kernel.LogBirth(categorical, counts), 12);
        // A two-coordinate birth uses the own-kind Wallenius law; the other kinds' appearances do not enter.
        foreach (int width in new[] { 1, 4 })
        {
            BassTerm born = At(model, new ProjectionDirection([1, 2], [.6, -.8]), .6, width);
            double wv = model.Config.Structure.VariableWeight;
            double[] variableWeights = width == 4 ? [wv+1, wv+1, wv] : [wv, wv, wv+1];
            double logQ = -Math.Log(3.0) + Math.Log(w + (width == 4 ? 1 : 0)) - Math.Log(2*w + 1)
                + Math.Log(DirectionalChainTests.WalleniusByQuadrature(variableWeights, [1,2]))
                + prior.Directional!.LogDirection(2) + prior.Directional.LogKnot(born.Directional!, [])
                + Math.Log(w+1) - Math.Log(2*w+1);
            Assert.Equal(logQ, model.Kernel.LogBirth(born, counts), 10);
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void MixedColumns_KeepMarginalsAndPredictions_WhenOneColumnIsRescaled(bool ridge)
    {
        BassConfig config = Config();
        if (ridge) config = config with { Coefficients = BassCoefficientPrior.Ridge(TauPrior.Held(1.0)) };
        BassModel model = BassModel.Create(config, CategoricalRidgeTests.Data());
        BassTerm[] terms = [At(model, Pair, .4, 4), At(model, Axis, .7), new(categorical: [new CategoricalFactor(0, [0])])];
        BassDesign design = model.NewDesign(); Assert.Null(design.Load(terms));
        double[][] columns = BassRidgeTests.CentredColumns(model, design.ColumnTerms);
        double[] y = model.Data.Response.Values.ToArray(); double mean = y.Average();
        double total = y.Sum(v => (v-mean)*(v-mean)); int n = y.Length, m = columns.Length;
        GPriorConditional Conditional(double scale, double tau)
        {
            var g = new double[m,m]; var c = new double[m];
            for (int j = 0; j < m; j++)
            {
                double sj = j == 2 ? scale : 1;
                for (int k = 0; k < m; k++) g[j,k] = sj*(k == 2 ? scale : 1)*columns[j].Zip(columns[k], (a,b) => a*b).Sum();
                c[j] = sj*columns[j].Zip(y, (a,b) => a*(b-mean)).Sum();
            }
            var factor = new DenseCholesky(m); factor.Decompose(g); double[] ols = factor.Solve(c);
            double q = ols.Zip(c, (a,b) => a*b).Sum();
            if (ridge) return CollapsedRidgePrior.FromCentred(RidgeSpectrum.Of(g,c,n), total-q, mean, model.Config.Intercept).At(tau);
            return CollapsedGPrior.FromCentredLeastSquares(factor, ols, total-q, q, n, mean, model.Config.Intercept).At(tau);
        }
        foreach (double tau in new[] {.03,1.0,30.0})
        {
            GPriorConditional unit = Conditional(1,tau);
            Assert.Equal(model.Target.Marginal(BassStructureFit.Of(design,0),tau).LogMarginal(model.Target.Sigma2), unit.LogMarginal(model.Target.Sigma2), 7);
            foreach (double scale in new[] {.001,-7.0,1000.0})
            {
                GPriorConditional moved = Conditional(scale,tau);
                Assert.Equal(unit.LogMarginal(model.Target.Sigma2), moved.LogMarginal(model.Target.Sigma2), 7);
                for (int i = 0; i < n; i++)
                {
                    double a = unit.Mean[0], b = moved.Mean[0];
                    for (int j = 0; j < m; j++)
                    {
                        a += unit.Mean[j+1]*columns[j][i];
                        b += moved.Mean[j+1]*(j == 2 ? scale : 1)*columns[j][i];
                    }
                    Assert.True(Math.Abs(a-b) < 1e-8, $"prediction change {a-b:E3}");
                }
            }
        }
        if (ridge)
        {
            var marginal = model.Target.Marginal(BassStructureFit.Of(design,0), 1.0);
            var dense = BassRidgeTests.Dense(model, design.ColumnTerms, 1.0, 1.0);
            Assert.Equal(dense.S, marginal.S, 8);
        }
        else Assert.Equal(DirectionalChainTests.DenseLogMarginal(y, [.. design.ColumnTerms.Select(t => { var c = new double[n]; model.Basis.CaseColumn(t, model.Inputs.Cases,c); return c; })],
            model.Target.HeldTau!.Value,3,2), model.Target.Marginal(BassStructureFit.Of(design,0),model.Target.HeldTau.Value).LogMarginal(model.Target.Sigma2), 7);
        Assert.Throws<ArgumentException>(() => (config with { Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Held(1.0)) }).Validate());
    }

    [Fact]
    public void Determinant_CountsBothDirectionalWidths_AndExcludesCategoricalOnly()
    {
        BassConfig config = Config() with { Structure = Config().Structure with { Directional = new DirectionalBasis { Columns = 4, MixedWidths = true, Directions = DirectionPrior.Determinant } } };
        BassModel model = BassModel.Create(config, CategoricalRidgeTests.Data());
        BassTerm smooth = At(model, Pair, 0.4, 4), hinge = At(model, Axis, 0.7);
        BassTerm categorical = new(categorical: [new CategoricalFactor(0, [0])]);
        Assert.Equal(Math.Log(1.5), model.Prior.LogCoupling([smooth, hinge, categorical]), 12);
        Assert.Equal(double.NegativeInfinity, model.Prior.LogCoupling([smooth, At(model, Pair, 0.6)]));
    }

    [Fact]
    public void StructureMoves_IncludeBoundaryMenus_AllTermDeaths_AndSymmetricWidthFlips()
    {
        BassModel model = BassModel.Create(Config(3), CategoricalRidgeTests.Data());
        BassTerm[] pool = [At(model, Pair, .4, 4), At(model, Axis, .7)];
        var rng = new Xoshiro256PlusPlus(3819);
        foreach (int count in new[] {0,1,2})
        {
            var moves = new BassStructureMoves(model, pool.Take(count).ToArray(), Math.Log(model.Target.HeldTau!.Value));
            BassState before = moves.Start;
            var counts = new BassTermCounts(model.Prior); foreach (BassTerm term in before.Terms) counts.Add(term);
            var birth = moves.Moves[0]; var death = moves.Moves[1];
            var proposed = birth.Propose(before, rng)!.Value;
            BassState after = proposed.Candidate;
            BassTerm born = after.Terms[^1];
            double forward = Math.Log(birth.Weight(before)/moves.Moves.Sum(m => m.Weight(before))) + model.Kernel.LogBirth(born, counts);
            double reverse = Math.Log(death.Weight(after)/moves.Moves.Sum(m => m.Weight(after))) - Math.Log(after.Count);
            double selection = Math.Log(death.Weight(after)/moves.Moves.Sum(m => m.Weight(after)))
                - Math.Log(birth.Weight(before)/moves.Moves.Sum(m => m.Weight(before)));
            Assert.Equal(reverse-forward, proposed.LogProposalRatio + selection, 11);
            Assert.Equal(0.0, proposed.LogJacobian);
            output.WriteLine($"M {count}->{after.Count}: log forward {forward:F8}, log reverse {reverse:F8}, menu ratio {selection:F8}");
        }
        var changes = new BassStructureMoves(model, pool, Math.Log(model.Target.HeldTau!.Value));
        bool sawWidth = false;
        for (int j = 0; j < 40; j++)
        {
            var move = changes.Moves[2].Propose(changes.Start, rng)!.Value;
            BassTerm moved = move.Candidate.Terms[^1];
            BassTerm? old = pool.SingleOrDefault(t => t.Directional!.Direction.Equals(moved.Directional!.Direction) && t.Directional.Knot == moved.Directional.Knot);
            if (old is null || old.Directional!.Columns == moved.Directional!.Columns) continue;
            sawWidth = true;
            Assert.Equal(0.0, move.LogProposalRatio);
            Assert.Equal(0.0, move.LogJacobian);
            Assert.Equal(changes.Moves[2].Weight(changes.Start), changes.Moves[2].Weight(move.Candidate));
            // Selection is uniform among M terms, then width with probability 1/2 in both states.
            Assert.Equal(changes.Start.Count, move.Candidate.Count);
        }
        Assert.True(sawWidth);
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void PriorOnlyChain_ReproducesPerKindCounts(bool categories, bool determinant)
    {
        BassData data = CategoricalRidgeTests.Data(80);
        if (!categories) data = CategoricalRidgeTests.Continuous(data);
        BassConfig config = Config(4);
        if (determinant) config = config with { Structure = config.Structure with { Directional = config.Structure.Directional! with { Directions = DirectionPrior.Determinant } } };
        BassModel model = BassModel.Create(config, data);
        int[] seeds = SeedTree.Derive(3813, 3);
        var rng = new Xoshiro256PlusPlus(seeds[0]);
        var sampler = new BassPriorSampler(model);
        double[] Stats(BassState state)
        {
            var result = new double[6]; result[0] = state.Count;
            foreach (BassTerm t in state.Terms)
            {
                result[1 + BassStructurePrior.KindOf(t)]++;
                if (t.Directional is { } h) { result[4] += h.Direction.Count; result[5] += h.Direction.Weights[0] > 0 ? 1 : 0; }
            }
            return result;
        }
        const int draws = 16_000, steps = 100_000;
        double[][] direct = [.. Enumerable.Range(0, draws).Select(_ => Stats(sampler.DrawState(rng).State))];
        var chain = new BassChain(model, [], Math.Log(model.Target.HeldTau!.Value), new Xoshiro256PlusPlus(seeds[1]), new Xoshiro256PlusPlus(seeds[2]));
        chain.Moves.SetTemperature(new Temperature(0.0, TemperingPath.Likelihood));
        for (int j = 0; j < 3000; j++) chain.Step();
        double[][] successive = [.. Enumerable.Range(0, steps).Select(_ => Stats(chain.Step()))];
        for (int k = 0; k < 6; k++)
        {
            double[] a = [.. direct.Select(s => s[k])], b = [.. successive.Select(s => s[k])];
            double ma = a.Average(), mb = b.Average();
            double variance = a.Sum(x => (x - ma) * (x - ma)) / (a.Length - 1);
            double se = Math.Sqrt(variance / draws + Math.Pow(DirectionalChainTests.BatchError(b, 40), 2));
            double z = se > 0.0 ? (ma - mb) / se : 0.0;
            output.WriteLine($"stat {k}: independent {ma:F5}, chain {mb:F5}, SE {se:F5}, z {z:F2}");
            Assert.True(Math.Abs(z) < 4.5, $"stat {k}: z={z:F2}");
        }
        output.WriteLine($"prior acceptance {sampler.AcceptanceRate:F4}; determinant rejections {sampler.DeterminantRejections}; width attempts {chain.Structure.Kinds.Entries.Single(e => e.Kind == BassMoveKind.WidthChange).Attempts}");
    }

    [Fact]
    public void MixedSamplerResolutionAndCheckpoint_KeepWidthsAndContinueExactly()
    {
        BassConfig config = Config(4) with { Sampling = Config().Sampling with { BurnIn = 120, Samples = 20, MasterSeed = 3814 } };
        BassData data = CategoricalRidgeTests.Data(80);
        Assert.NotNull(BassModel.Create(config, data).Config.Sampling.Tempering);
        BassConfig plain = config with { Sampling = config.Sampling with { Sampler = BassSamplerChoice.PlainChains } };
        Assert.NotNull(BassModel.Create(plain, data).SamplerCaution);
        BassFit plainFit = Bass.Fit(plain, data);
        Assert.Null(plainFit.Diagnostics.Tempering);
        Assert.NotNull(plainFit.Diagnostics.SamplerCaution);
        BassConfig smc = config with { Sampling = config.Sampling with { Smc = new SmcSpec { Particles = 40, Replicates = 1, MutationSteps = 2 } } };
        Assert.Same(smc, BassModel.Create(smc, data).Config);
        Assert.True(Bass.Fit(smc, data).Diagnostics.Smc!.ReachedPosterior);
        BassFit whole = Bass.Fit(config with { Sampling = config.Sampling with { Samples = 40 } }, data);
        BassFit half = Bass.Fit(config, data, new BassRunOptions { Checkpoint = true, MaxDegreeOfParallelism = 1 });
        BassFit resumed = Bass.Resume(half.Checkpoint!, config, data);
        BassRunTests.AssertSameFit(whole, resumed);
        Assert.True(BassConfig.FromJson(config.ToJson()).Structure.Directional!.MixedWidths);
        BassConfig legacy = config with { Structure = config.Structure with { Directional = new DirectionalBasis { Columns = 4 } } };
        Assert.DoesNotContain("MixedWidths", legacy.ToJson());
        Assert.Equal(legacy.ToJson(), BassConfig.FromJson(legacy.ToJson()).ToJson());
    }

    [Theory]
    [InlineData("plain")]
    [InlineData("pt")]
    [InlineData("smc")]
    public void RetainedKindCountsAndWidthTally_ReadActualMixedFits(string sampler)
    {
        BassConfig config = Config(4) with { Sampling = Config().Sampling with
            { Chains = 2, BurnIn = 120, Samples = 80, MasterSeed = 3818 } };
        if (sampler == "plain") config = config with { Sampling = config.Sampling with { Sampler = BassSamplerChoice.PlainChains } };
        if (sampler == "smc") config = config with { Sampling = config.Sampling with
            { Chains = 1, Smc = new SmcSpec { Particles = 32, Replicates = 2, MutationSteps = 2 } } };
        BassData data = CategoricalRidgeTests.Data(80);
        BassModel model = BassModel.Create(config, data);
        BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        double[]? originalWeights = fit.DrawWeights?.ToArray();
        int tableSize = fit.Terms.Count, pooledIndex = 0;
        double weightedKinds = 0, weightedM = 0;
        foreach (IReadOnlyList<BassDraw> start in fit.Draws)
            foreach (BassDraw draw in start)
            {
                var actual = draw.Structure.DirectionalTermCounts(fit.Terms);
                IReadOnlyList<BassTerm> units = model.Expansion!.Units([.. draw.Structure.TermIds.Select(id => fit.Terms[id])]);
                Assert.Equal(units.Count(t => t.Directional?.Columns == 4), actual.Blocks);
                Assert.Equal(units.Count(t => t.Directional?.Columns == 1), actual.Hinges);
                Assert.Equal(units.Count(t => t.Directional is null), actual.CategoricalOnly);
                int total = actual.Blocks + actual.Hinges + actual.CategoricalOnly;
                Assert.Equal(draw.Structure.TermCount(fit.Terms), total);
                Assert.Equal(4 * actual.Blocks + actual.Hinges + actual.CategoricalOnly, draw.Structure.TermIds.Length);
                double weight = fit.DrawWeights?[pooledIndex++] ?? 1.0;
                weightedKinds += weight * total;
                weightedM += weight * units.Count;
            }
        Assert.Equal(weightedM, weightedKinds);
        Assert.Equal(tableSize, fit.Terms.Count);
        if (originalWeights is not null) Assert.Equal(originalWeights, fit.DrawWeights);
        BassKindCount width = fit.Diagnostics.MoveKinds[BassMoveKind.WidthChange];
        Assert.True(width.Attempts > 0);
        BassKindCount[] accounts = fit.Diagnostics.Smc is { } smc
            ? [.. smc.Replicates.Select(r => r.Account.Kinds[BassMoveKind.WidthChange])]
            : [.. fit.Diagnostics.Chains.Select(c => c.MoveKinds[BassMoveKind.WidthChange])];
        Assert.Equal(2, accounts.Length);
        Assert.Equal(width.Attempts, accounts.Sum(c => c.Attempts));
        Assert.Equal(width.Accepted, accounts.Sum(c => c.Accepted));
        Assert.Equal(width.Null, accounts.Sum(c => c.Null));
        Assert.Equal(width.Refused, accounts.Sum(c => c.Refused));
        output.WriteLine($"{sampler}: width attempts {width.Attempts}, accepted {width.Accepted}; weighted M {weightedM:R}");
    }
}
