using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.BassTemperingTests;

namespace Maths.Regression.Tests;

/// <summary>
/// E21's ridge blocks beyond KHAOS's first parts (<c>issues/regression/bass/ridge-blocks-brief.md</c>), slice by slice.
/// <list type="bullet">
/// <item><b>Slice 0, the whitening on the working response:</b> a held-s block's fits keep their draws bit for bit, pinned by hash on a
/// hinge chain, a Legendre chain and a tempered hinge fit with its reference rung; the hinge prefix tables stand aside under a block, whose
/// columns they would read unwhitened; a design follows a reading at another scale as it follows moved weights, to the fit of a design
/// built there.</item>
/// <item><b>Slice 1, s drawn:</b> the half-Cauchy on s, declared apart and serialized only when declared, its reading's priors frozen on the
/// data; the scale's step keeps its exact law in ln s on a frozen structure; each draw carries its s and reads at it; each rung of a ladder
/// with its reference rung samples its law over F, τ and ln s, and SMC meets its evidence, against an enumeration that integrates ln s; a
/// tempered fit resumes bit for bit.</item>
/// </list>
/// </summary>
public sealed class RandomEffectsTests(ITestOutputHelper output)
{
    // A fit's hash: the bit-identity tests' (structures, τ, σ², coefficients, the move counts, the ladder), with each draw's group effects.
    private static string Hash(BassFit fit)
    {
        var hash = new RunFingerprint().Add(BassBitIdentityTests.Hash(fit));
        foreach (DrawStore<BassDraw> chain in fit.Draws)
            foreach (BassDraw d in chain) hash.Add(d.GroupEffects ?? []);
        return hash.ToString();
    }

    private static (BassConfig Config, BassData Data) Pinned(string fit)
    {
        var sampling = new BassSampling { Chains = 2, MasterSeed = 2108, BurnIn = 300, Samples = 200 };
        if (fit == "tempered hinge")
        {
            var spec = new TemperingSpec { Levels = 3, HotEnd = 0.2, ReferenceRung = true };
            (BassModel model, BassConfig config) = RandomInterceptTests.Fixture("hinge", spec);
            return (config with { Sampling = sampling with { Tempering = spec } }, model.Data);
        }
        (BassModel m, BassConfig c) = RandomInterceptTests.Fixture(fit);
        return (c with { Sampling = sampling }, m.Data);
    }

    // Slice 0's gate: moving the whitening onto the working response leaves a held-s block's fits as they were.
    [Theory]
    [Trait("Tier", "Pinned")]
    [InlineData("hinge", "b9b95aee57ef459c")]
    [InlineData("legendre", "a1f4528886afa5f7")]
    [InlineData("tempered hinge", "05122582e8f56f75")]
    public void AHeldScaleFit_KeepsItsDraws(string fit, string expected)
    {
        (BassConfig config, BassData data) = Pinned(fit);
        string hash = Hash(Bass.Fit(config, data));
        output.WriteLine($"{fit}: {hash}");
        Assert.Equal(expected, hash);
    }

    // The design follows a reading at another scale (WorkingResponse.Whitening) as it follows moved weights: every column whitened afresh and
    // the Gram factored strictly, so its fit is a design's built on that reading and loaded with the same terms. Back on its first reading
    // it holds its first fit.
    [Fact]
    public void TheDesign_FollowsTheWhitening_ToAnotherScale()
    {
        (BassModel model, _) = RandomInterceptTests.Fixture("hinge");
        WorkingResponse At(double s) => RandomInterceptBlock.Of(model.Data.Grouping!, new BassRandomIntercept(s)).Whitened(new WorkingResponse(model.Data.Response));
        BassDesign moving = model.NewDesign();
        var rng = new Xoshiro256PlusPlus(seed: 41);
        for (int k = 0; k < 200 && moving.Count < 4; k++)
        {
            var counts = new BassTermCounts(model.Prior);
            foreach (BassTerm t in moving.Terms) counts.Add(t);
            if (model.Proposal.DrawBirth(rng, counts) is { } term && moving.ProposeBirth(term).Admitted) moving.Accept();
        }
        Assert.True(moving.Count >= 2, $"the design holds {moving.Count} terms");
        double rss0 = moving.ResidualSquares, fit0 = moving.FittedSquares;
        double[] beta0 = moving.Coefficients.ToArray();

        WorkingResponse other = At(0.4);
        Assert.True(moving.ProposeFollow(other).Admitted);
        moving.Accept();
        Assert.Same(other.Whitening, moving.Working.Whitening);
        BassDesign fresh = model.NewDesign(other);
        foreach (BassTerm t in moving.Terms)
        {
            Assert.True(fresh.ProposeBirth(t).Admitted);
            fresh.Accept();
        }
        double Rel(double a, double b) => Math.Abs(a - b) / Math.Max(Math.Abs(a), Math.Abs(b));
        double worst = Math.Max(Rel(moving.ResidualSquares, fresh.ResidualSquares), Rel(moving.FittedSquares, fresh.FittedSquares));
        double scale = Math.Max(1.0, fresh.Coefficients.ToArray().Max(Math.Abs));
        for (int j = 0; j < fresh.Count; j++) worst = Math.Max(worst, Math.Abs(moving.Coefficients[j] - fresh.Coefficients[j]) / scale);
        Assert.NotEqual(rss0, moving.ResidualSquares);   // the scale moved the fit
        output.WriteLine($"{moving.Count} terms; at s = 0.4 against a design built there, the worst relative gap {worst:E2}; ‖r‖² {rss0:G6} at 1.3, {moving.ResidualSquares:G6} at 0.4");
        Assert.True(worst <= 1e-12, $"{worst:E2}");

        Assert.True(moving.ProposeFollow(model.Working).Admitted);
        moving.Accept();
        double back = Math.Max(Rel(moving.ResidualSquares, rss0), Rel(moving.FittedSquares, fit0));
        for (int j = 0; j < beta0.Length; j++) back = Math.Max(back, Math.Abs(moving.Coefficients[j] - beta0[j]) / scale);
        Assert.True(back <= 1e-12, $"back at s = 1.3: {back:E2}");
        Assert.Throws<ArgumentException>(() => moving.ProposeFollow(new WorkingResponse(model.Data.Response)));   // no whitening: another layout
    }

    // Under a block the design holds whitened columns, V^(−½)xⱼ, and a birth's border is their inner product with the new column's
    // whitened image. The hinge prefix tables sum a held column against the raw hinge, which would give xⱼᵀV^(−½)h for xⱼᵀV⁻¹h, so they
    // stand aside: with the tables asked for, every proposal is the dots' own.
    [Fact]
    public void UnderABlock_TheTablesStandAside()
    {
        var point = BassCost.Base with { N = 600, P = 6, Cap = 15 };
        BassData plain = BassCost.Data(point);
        var data = new BassData(plain.Predictors, plain.Response, RandomInterceptTests.Grouping(point.N, 8, seed: 11));
        BassModel model = BassModel.Create(BassCost.Config(point) with { RandomIntercept = new BassRandomIntercept(1.3) }, data);
        BassDesign off = model.NewDesign(), on = model.NewDesign();
        off.TablesOn = false;
        on.TablesOn = true;
        var rng = new Xoshiro256PlusPlus(seed: 5151);
        int proposals = 0;
        double worst = 0.0;
        for (int step = 0; step < 1500; step++)
        {
            var counts = new BassTermCounts(model.Prior);
            foreach (BassTerm t in off.Terms) counts.Add(t);
            DesignProposal a, b;
            if ((rng.NextDouble() < 0.6 || off.Count == 0) && off.Count < model.MaxBasis)
            {
                if (model.Proposal.DrawBirth(rng, counts) is not { } t) continue;
                a = off.ProposeBirth(t);
                b = on.ProposeBirth(t);
            }
            else if (off.Count > 0)
            {
                int i = (int)(rng.NextUInt64() % (ulong)off.Count);
                a = off.ProposeDeath(i);
                b = on.ProposeDeath(i);
            }
            else continue;
            proposals++;
            Assert.Equal(a.Refusal, b.Refusal);
            if (!a.Admitted) continue;
            worst = Math.Max(worst, Math.Abs(a.ResidualSquares - b.ResidualSquares) / a.ResidualSquares);
            if (rng.NextDouble() < 0.5)
            {
                off.Accept();
                on.Accept();
            }
        }
        output.WriteLine($"{proposals} proposals; the tables read {on.TableLookups} entries; the worst relative ‖r‖² between tables and dots {worst:E2}");
        Assert.Equal(0L, on.TableLookups);
        Assert.Equal(0.0, worst);
    }

    // ---- Slice 1: s drawn ----------------------------------------------------------------------------------------------------------------

    // KHAOS's fixtures with s drawn from C⁺(0, 1), starting at its held 1.3.
    private static (BassModel Model, BassConfig Config) Drawn(string fit, TemperingSpec? spec = null)
    {
        (BassModel held, BassConfig config) = RandomInterceptTests.Fixture(fit, spec);
        BassConfig drawn = config with { RandomIntercept = config.RandomIntercept! with { ScalePrior = new BassGroupScalePrior() } };
        return (BassModel.Create(drawn, held.Data), drawn);
    }

    [Fact]
    public void TheScalePrior_IsDeclaredApart_ItsPriorsFrozenOnTheData_AndRefusedWhereNotBuilt()
    {
        (BassModel held, BassConfig config) = RandomInterceptTests.Fixture("hinge");
        (BassModel model, BassConfig drawn) = Drawn("hinge");
        Assert.DoesNotContain("ScalePrior", config.ToJson());
        Assert.Contains("\"ScalePrior\"", drawn.ToJson());
        Assert.Equal(drawn, BassConfig.FromJson(drawn.ToJson()));
        Assert.NotEqual(Bass.Fingerprint(config, held.Data), Bass.Fingerprint(drawn, held.Data));
        // The reading carries the priors' frozen reference, the response's own moments, so a chain may follow it to another s.
        Assert.True(model.ScaleDrawn && !held.ScaleDrawn);
        Assert.NotNull(model.Working.Reference);
        Assert.False(model.PriorsReadResponse);
        Assert.Equal(model.Data.Response.Values.ToArray().Average(), model.Working.Reference!.Level!.Value, 12);
        Assert.Equal(["tau", "group scale"], BassSweepSlots.Blocks(model, BassSweepSlot.Tau).Select(b => b.Name));
        Assert.Contains(model.Resolved.AutoPicks, p => p.Setting == "RandomIntercept" && p.Value.Contains("s drawn"));
        Assert.Contains("RandomIntercept.ScalePrior.Scale", Assert.ThrowsAny<ArgumentException>(
            () => (drawn with { RandomIntercept = drawn.RandomIntercept! with { ScalePrior = new BassGroupScalePrior(0.0) } }).Validate()).Message);
        BassConfig flat = drawn with
        {
            Sampling = drawn.Sampling with
            {
                FlatHistogram = new BassFlatHistogramSpec { Run = new FlatHistogramSpec { LearningIterations = 3_000, ProductionDraws = 500, Nodes = 20 } },
            },
        };
        Assert.Contains("RandomIntercept.ScalePrior", Assert.ThrowsAny<ArgumentException>(() => flat.Validate()).Message);

        // The half-Cauchy: its density's mass and its draws' median, A.
        var prior = new BassGroupScalePrior(0.7);
        double mass = 0.0;
        for (int k = 0; k < 400_000; k++)
        {
            double t = -25.0 + 50.0 * (k + 0.5) / 400_000;   // in ln s, its Jacobian eᵗ
            mass += Math.Exp(prior.LogDensity(Math.Exp(t)) + t) * 50.0 / 400_000;
        }
        var rng = new Xoshiro256PlusPlus(17);
        int below = Enumerable.Range(0, 200_000).Count(_ => prior.Sample(rng) < 0.7);
        output.WriteLine($"C⁺(0, 0.7): mass {mass:F8}; draws below its median {below / 200_000.0:F4}");
        Assert.Equal(1.0, mass, 6);
        Assert.InRange(below / 200_000.0, 0.4955, 0.5045);
    }

    // The step's statement on a frozen structure and τ: its stationary law in u = ln s is exactly ∝ π_b(state at eᵘ)·p(eᵘ)·eᵘ, the marginal
    // with a, σ² and u integrated read by a refit at each point of a grid, which its chain meets in mean and mean square (batch means).
    [Theory]
    [InlineData("hinge")]
    [InlineData("legendre")]
    public void TheScalesStep_KeepsItsLaw_OnAFrozenStructure(string fit)
    {
        (BassModel model, _) = Drawn(fit);
        var start = new BassChain(model, [], 0.0, new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        for (int k = 0; k < 300; k++) start.Step();
        var moves = new BassStructureMoves(model, start.Current.Terms, start.Current.LogTau);
        var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(3),
                                                       new Temperature(1.0, TemperingPath.Likelihood));
        var block = new BassGroupScaleBlock(model);
        RandomInterceptBlock at0 = BassGroupScaleBlock.Reading(chain).Whitening!;
        BassGroupScalePrior prior = model.RandomIntercept!.Spec.ScalePrior!;
        double LogLaw(double u)
        {
            double s = Math.Exp(u);
            BassState candidate = moves.ProposeFollow(chain.Current, at0.At(s).Whitened(model.BaseReading));
            return chain.Target.Evaluate(candidate).LogLikelihood(1.0) + prior.LogDensity(s) + u;
        }
        (double Mean, double Square) Moments(int radius, double step)
        {
            double[] grid = [.. Enumerable.Range(0, (int)(2 * radius / step) + 1).Select(k => -radius + step * k)];
            double[] logs = [.. grid.Select(LogLaw)];
            double top = logs.Max();
            Assert.True(Math.Exp(logs[0] - top) < 1e-8 && Math.Exp(logs[^1] - top) < 1e-8, "the grid holds the law's mass");
            double[] mass = [.. logs.Select((l, k) => Math.Exp(l - top) * (k == 0 || k == logs.Length - 1 ? 0.5 : 1.0))];
            double z = mass.Sum();
            return (mass.Select((w, k) => w * grid[k]).Sum() / z, mass.Select((w, k) => w * grid[k] * grid[k]).Sum() / z);
        }
        var coarse = Moments(24, 0.05);
        (double mean, double square) = Moments(28, 0.025);
        output.WriteLine($"grid refinement: ΔE[ln s] {mean - coarse.Mean:E2}, ΔE[(ln s)²] {square - coarse.Square:E2}");
        Assert.InRange(Math.Abs(mean - coarse.Mean), 0, 1e-7);
        Assert.InRange(Math.Abs(square - coarse.Square), 0, 1e-7);

        const int steps = 60_000, batches = 50;
        var us = new double[steps];
        int accepted = 0;
        var rng = new Xoshiro256PlusPlus(20261109);
        for (int k = 0; k < steps; k++)
        {
            if (block.Update(chain, rng)) accepted++;
            us[k] = Math.Log(BassGroupScaleBlock.Reading(chain).Whitening!.Spec.Scale);
        }
        Assert.Equal(chain.Current.Terms.Select(t => t.ToString()), start.Current.Terms.Select(t => t.ToString()));
        double Z(Func<double, double> g, double expected)
        {
            int size = steps / batches;
            double[] means = [.. Enumerable.Range(0, batches).Select(b => us.Skip(b * size).Take(size).Average(g))];
            (double m, double error) = Summary(means);
            return (m - expected) / error;
        }
        double zMean = Z(u => u, mean), zSquare = Z(u => u * u, square);
        output.WriteLine($"{fit}: {start.Current.Count} terms; E[ln s] {us.Average():F4} against {mean:F4} (z {zMean:F2}), E[(ln s)²] z {zSquare:F2}; " +
                         $"sd {Math.Sqrt(square - mean * mean):F3}; accepted {accepted} of {steps}");
        Assert.True(Math.Abs(zMean) < 4.0 && Math.Abs(zSquare) < 4.0, $"z {zMean:F2}, {zSquare:F2}");
        Assert.InRange(accepted / (double)steps, 0.15, 0.85);
    }

    [Fact]
    public void ADrawnScalesDraws_CarryTheirScale_AndReadAtIt()
    {
        (BassModel model, BassConfig config) = Drawn("hinge");
        config = config with { Sampling = new BassSampling { Chains = 2, MasterSeed = 7, BurnIn = 300, Samples = 300 } };
        BassFit fit = Bass.Fit(config, model.Data);
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        Assert.All(draws, d => Assert.NotNull(d.GroupScale));
        int distinct = draws.Select(d => d.GroupScale).Distinct().Count();
        MoveStat step = fit.Diagnostics.Chains[0].Blocks!.Single(k => k.Key == "group scale");
        output.WriteLine($"{distinct} distinct s among {draws.Length} draws; the step accepted {step.Accepted} of {step.Attempts}");
        Assert.True(distinct > draws.Length / 10, "s moves");
        var options = new BassPredictOptions { Predictive = true, Seed = 9 };
        BassPrediction prediction = fit.Predict(model.Data, PredictionCases.Of(model.Data.Predictors), options);
        double secondMoment = 0.0;
        for (int d = 0; d < prediction.DrawCount; d++)
        {
            BassDraw draw = draws[prediction.Draws[d]];
            double variance = draw.Sigma2 * (1 + Math.Pow(draw.GroupScale!.Value, 2));
            for (int i = 0; i < model.Data.Response.Cases; i++)
                secondMoment += Math.Pow(prediction.Predictive![d, i, 0] - prediction.Latent[d, i, 0], 2) / variance;
        }
        Assert.InRange(secondMoment / (model.Data.Response.Cases * prediction.DrawCount), 0.94, 1.06);
        foreach (BassDraw d in draws.Take(40))
        {
            // The draw keeps its s, and its reading is rebuilt at it; its effects keep the convention there.
            Assert.Null(d.Working);
            RandomInterceptBlock at = d.WorkingOn(model).Whitening!;
            Assert.Equal(d.GroupScale!.Value, at.Spec.Scale);
            double centred = 0.0, spread = 0.0;
            for (int g = 0; g < at.Levels; g++)
            {
                centred += at.Counts[g] * d.GroupEffects![g];
                spread += at.Counts[g] * Math.Abs(d.GroupEffects[g]);
            }
            Assert.True(Math.Abs(centred) <= 1e-12 * spread, $"Σ n_g u_g = {centred:E2}");
        }
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void ADrawnScaleFit_ResumesBitForBit(bool tempered)
    {
        var spec = tempered ? new TemperingSpec { Levels = 3, HotEnd = 0.2, ReferenceRung = true } : null;
        (BassModel model, BassConfig config) = Drawn("hinge", spec);
        BassConfig whole = config with
        {
            Sampling = new BassSampling { Chains = 2, MasterSeed = 31, BurnIn = 100, Samples = 120, Thin = 2, Tempering = spec, ReportEvidence = tempered },
        };
        BassConfig half = whole with { Sampling = whole.Sampling with { Samples = 60 } };
        BassFit full = Bass.Fit(whole, model.Data);
        BassFit first = Bass.Fit(half, model.Data, new BassRunOptions { Checkpoint = true });
        BassFit resumed = Bass.Resume(first.Checkpoint!, half, new BassData(model.Data.Predictors, new ScalarResponse(model.Data.Response.Values.ToArray()), model.Data.Grouping));
        BassRunTests.AssertSameFit(full, resumed);
        for (int c = 0; c < full.Draws.Count; c++)
            for (int d = 0; d < full.Draws[c].Count; d++) Assert.Equal(full.Draws[c][d].GroupScale, resumed.Draws[c][d].GroupScale);
        Assert.True(full.Draws.SelectMany(c => c).Select(d => d.GroupScale).Distinct().Count() > 1, "s moves");
    }

    [Fact]
    public void PoweredRungs_ReadTheLogScalePriorInTheirJointDensity()
    {
        (BassModel model, _) = Drawn("hinge");
        var moves = new BassStructureMoves(model, [], 0.0);
        var temperature = new Temperature(0.4, TemperingPath.Powered);
        var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(7), temperature);
        RandomInterceptBlock at = BassGroupScaleBlock.Reading(chain).Whitening!;
        var update = new BassGroupScaleBlock(model);
        double t = Math.Log(2.0), h = update.Scale(at);
        var proposed = update.Propose(chain, at, h, t)!.Value;
        var next = proposed.Candidate.Fit.Working!.Whitening!;
        double reverse = update.Scale(next);
        // For C+(0,1), p(ln s) ∝ s/(1+s²).
        double logPriorRatio = t + Math.Log((1 + at.Spec.Scale * at.Spec.Scale) / (1 + next.Spec.Scale * next.Spec.Scale));
        ITargetEvaluation old = model.Target.Evaluate(chain.Current);
        Assert.Equal(logPriorRatio, proposed.Evaluation.LogPrior - old.LogPrior, 12);
        double powered = 0.4 * (logPriorRatio + proposed.Evaluation.LogLikelihood(1.0) - old.LogLikelihood(1.0));
        Assert.Equal(powered, temperature.LogDensity(proposed.Evaluation) - chain.CurrentLogDensity, 12);
        Assert.Equal(powered + Math.Log(h / reverse) + t * t * (1 / (h * h) - 1 / (reverse * reverse)) / 2, proposed.LogAcceptance, 12);
    }

    // Each rung's law over F, τ and v = ln s: per structure the target on a grid of (ln τ, v) by the trapezoid, with s's half-Cauchy and its
    // Jacobian eᵛ, a design per s on the reading whitened there; ln Z against the prior normalized on F.
    private sealed record ScaleLaw(double LogZ, double[] Count, double[] Structures, double LogTau, double LogScale);

    // F is checked at each s. A structure refused at the starting scale may enter it elsewhere.
    private static List<(BassTerm[] Terms, double LogPrior)> SmallStructures(BassModel model)
    {
        Assert.Equal(2, model.MaxBasis);
        List<BassTerm> terms = BassTermSpace.Terms(model);
        var structures = new List<(BassTerm[], double)>();
        void Add(params BassTerm[] set) => structures.Add((set, model.Prior.LogCollectionWeight(set.Length) + set.Sum(model.Prior.LogPrior)));
        Add();
        for (int i = 0; i < terms.Count; i++)
        {
            Add(terms[i]);
            for (int j = i + 1; j < terms.Count; j++) Add(terms[i], terms[j]);
        }
        return structures;
    }

    private static ScaleLaw[] ExactWithScale(BassModel model, double[] rungs, int radius = 22, double step = 0.1)
    {
        BassTarget target = model.Target;
        List<(BassTerm[] Terms, double LogPrior)> states = SmallStructures(model);
        RandomInterceptBlock block = model.RandomIntercept!;
        BassGroupScalePrior prior = block.Spec.ScalePrior!;
        double[] taus = target.SamplesTau ? [.. Enumerable.Range(0, (int)(68 / step) + 1).Select(i => -60 + step * i)] : [Math.Log(target.HeldTau!.Value)];
        double[] tauPrior = [.. taus.Select(u => target.SamplesTau ? target.LogTauDensity(u) : 0.0)];
        double[] scales = [.. Enumerable.Range(0, (int)(2 * radius / step) + 1).Select(i => -radius + step * i + Math.Log(prior.Scale))];
        double logCell = (target.SamplesTau ? Math.Log(step) : 0.0) + Math.Log(step);
        int cells = states.Count * scales.Length;
        var priorMass = new double[cells];
        Array.Fill(priorMass, double.NegativeInfinity);
        var logMass = new double[rungs.Length][];
        var tauMean = new double[rungs.Length][];
        for (int r = 0; r < rungs.Length; r++) (logMass[r], tauMean[r]) = (new double[cells], new double[cells]);
        var values = new double[taus.Length];
        for (int v = 0; v < scales.Length; v++)
        {
            double s = Math.Exp(scales[v]), logScalePrior = prior.LogDensity(s) + scales[v];
            if (v == 0 || v == scales.Length - 1) logScalePrior -= Math.Log(2.0);
            BassDesign design = model.NewDesign(block.At(s).Whitened(model.BaseReading));
            for (int k = 0; k < states.Count; k++)
            {
                BassTerm[] terms = states[k].Terms;
                int cell = k * scales.Length + v;
                if (design.Load(terms) is not null)
                {
                    foreach (double[] mass in logMass) mass[cell] = double.NegativeInfinity;
                    continue;
                }
                BassStructureFit fit = BassStructureFit.Of(design, terms.Sum(model.Prior.LogPrior));
                if (target.Refusal(fit) is null) priorMass[cell] = states[k].LogPrior + logScalePrior + Math.Log(step);
                BassEvaluation[] e = [.. taus.Select(u => target.SamplesTau ? target.Evaluate(fit, u) : target.Evaluate(fit))];
                for (int r = 0; r < rungs.Length; r++)
                {
                    for (int j = 0; j < taus.Length; j++)
                        values[j] = e[j].OffSupport ? double.NegativeInfinity : states[k].LogPrior + tauPrior[j]
                            + e[j].LogLikelihood(rungs[r]) + logScalePrior + logCell
                            - (target.SamplesTau && (j == 0 || j == taus.Length - 1) ? Math.Log(2.0) : 0.0);
                    if (values.All(double.IsNegativeInfinity))
                    {
                        logMass[r][k * scales.Length + v] = double.NegativeInfinity;
                        continue;
                    }
                    double m = LogSumExp(values), u = 0.0;
                    for (int j = 0; j < taus.Length; j++) u += Math.Exp(values[j] - m) * taus[j];
                    (logMass[r][k * scales.Length + v], tauMean[r][k * scales.Length + v]) = (m, u);
                }
            }
        }
        double logPriorMass = LogSumExp(priorMass);
        var laws = new ScaleLaw[rungs.Length];
        for (int r = 0; r < rungs.Length; r++)
        {
            double z = LogSumExp(logMass[r]), logTau = 0.0, logScale = 0.0;
            var count = new double[model.MaxBasis + 1];
            var structures = new double[states.Count];
            for (int k = 0; k < states.Count; k++)
                for (int v = 0; v < scales.Length; v++)
                {
                    double w = Math.Exp(logMass[r][k * scales.Length + v] - z);
                    count[states[k].Terms.Length] += w;
                    structures[k] += w;
                    logTau += w * tauMean[r][k * scales.Length + v];
                    logScale += w * scales[v];
                }
            laws[r] = new ScaleLaw(z - logPriorMass, count, structures, logTau, logScale);
        }
        return laws;
    }

    // Each rung's batch means of the count's indicators, of ln τ and of ln s, as BassTemperingTests.RungBatches reads the first two.
    private static (double[][,] Count, double[][,] Structures, double[][] LogTau, double[][] LogScale) RungBatches(BassModel model, TemperingSpec spec, double[] rungs, int ladders,
                                                                                         int batches, int burn, int scans)
    {
        int levels = rungs.Length + 1, per = scans / batches;   // last row: an ordinary posterior chain
        var count = new double[levels][,];
        var indices = SmallStructures(model).Select((x, k) => (Key: BassPriorSamplerTests.Key(x.Terms), Index: k)).ToDictionary(x => x.Key, x => x.Index);
        var structures = new double[levels][,];
        var logTau = new double[levels][];
        var logScale = new double[levels][];
        for (int r = 0; r < levels; r++)
        {
            count[r] = new double[model.MaxBasis + 1, ladders * batches];
            structures[r] = new double[indices.Count, ladders * batches];
            logTau[r] = new double[ladders * batches];
            logScale[r] = new double[ladders * batches];
        }
        var sampler = new BassPriorSampler(model);
        for (int l = 0; l < ladders; l++)
        {
            BassState start = sampler.DrawState(new Xoshiro256PlusPlus(100 + l)).State;
            var ladder = new BassLadder(model, spec, rungs, start.Terms, start.LogTau, seed: 200 + l);
            for (int t = 0; t < burn; t++) ladder.Step();
            // Hand the cold state to an ordinary chain after tempered burn-in. Its production samples use no swaps.
            ReversibleJumpChain<BassState> cold = ladder.Tempering.Replicas[0];
            var single = new BassChain(model, [], cold.Current.LogTau, new Xoshiro256PlusPlus(300 + l), new Xoshiro256PlusPlus(400 + l));
            single.Structure.Adopt(cold.Current, BassStructureMoves.Of(cold).DesignFor(cold.Current));
            single.Moves.SetState(cold.Current, model.Target.Evaluate(cold.Current));
            const int singleLength = 10;
            for (int t = 0; t < singleLength * burn; t++) single.Step();
            void Observe(ReversibleJumpChain<BassState> replica, int r, int batch, int size)
            {
                count[r][replica.Current.Count, batch] += 1.0 / size;
                structures[r][indices[BassPriorSamplerTests.Key(replica.Current.Terms)], batch] += 1.0 / size;
                logTau[r][batch] += replica.Current.LogTau / size;
                logScale[r][batch] += Math.Log(replica.Current.GroupScale!.Value) / size;
            }
            for (int t = 0; t < scans; t++)
            {
                ladder.Step();
                for (int r = 0; r < rungs.Length; r++)
                    Observe(ladder.Tempering.Replicas[r], r, l * batches + t / per, per);
            }
            for (int t = 0; t < singleLength * scans; t++)
            {
                single.Step();
                Observe(single.Moves, levels - 1, l * batches + t / (singleLength * per), singleLength * per);
            }
        }
        return (count, structures, logTau, logScale);
    }

    [Theory]
    [InlineData("hinge")]
    [InlineData("legendre")]
    public void Ladder_EachRungSamplesItsLaw_WithTheScaleDrawn(string fit)
    {
        double[] rungs = [1.0, 0.4, 0.0];
        var spec = new TemperingSpec { Levels = rungs.Length, HotEnd = 0.4, ReferenceRung = true };
        (BassModel model, _) = Drawn(fit, spec);
        const int ladders = 4, batches = 20;
        var (count, structures, logTau, logScale) = RungBatches(model, spec, rungs, ladders, batches, burn: 1_000, scans: 20_000);
        ScaleLaw[] exact = ExactWithScale(model, rungs);
        ScaleLaw[] coarse = ExactWithScale(model, rungs, radius: 16, step: 0.2);
        for (int r = 0; r < rungs.Length; r++)
        {
            double delta = new[] { exact[r].LogZ - coarse[r].LogZ, exact[r].LogScale - coarse[r].LogScale, exact[r].LogTau - coarse[r].LogTau }
                .Concat(exact[r].Structures.Zip(coarse[r].Structures, (a, b) => a - b)).Max(Math.Abs);
            output.WriteLine($"b = {rungs[r]}: max quadrature change {delta:E2} in ln Z, the two log moments and structure probabilities");
            Assert.InRange(delta, 0, 1e-5);
        }
        // The unit half-Cauchy is symmetric in log coordinates. Check that conditioning on F in these fixtures
        // preserves its zero log mean to quadrature accuracy, and that the reference rung integrates to one.
        Assert.InRange(Math.Abs(exact[^1].LogScale), 0, 1e-7);
        Assert.InRange(Math.Abs(exact[^1].LogZ), 0, 1e-6);
        exact = [.. exact, exact[0]];
        rungs = [.. rungs, 1.0];
        for (int r = 0; r < rungs.Length; r++)
        {
            for (int k = 0; k < exact[r].Structures.Length; k++)
            {
                (double mean, double error) = Summary([.. Enumerable.Range(0, ladders * batches).Select(i => structures[r][k, i])]);
                Assert.True(Math.Abs(mean - exact[r].Structures[k]) <= 5 * error + 0.004, $"{fit}, row {r}, structure {k}: {mean:F4} against {exact[r].Structures[k]:F4}");
            }
            for (int m = 0; m <= model.MaxBasis; m++)
            {
                (double mean, double error) = Summary([.. Enumerable.Range(0, ladders * batches).Select(k => count[r][m, k])]);
                Assert.True(Math.Abs(mean - exact[r].Count[m]) <= 5 * error + 0.004, $"{fit}, b = {rungs[r]}, M = {m}: {mean:F4} against {exact[r].Count[m]:F4}");
            }
            (double tau, double tauError) = Summary(logTau[r]);
            (double scale, double scaleError) = Summary(logScale[r]);
            output.WriteLine($"{fit}, b = {rungs[r]:G3}: E[ln τ] {tau:F4} ± {tauError:F4} against {exact[r].LogTau:F4}; E[ln s] {scale:F4} ± {scaleError:F4} against {exact[r].LogScale:F4}");
            Assert.True(Math.Abs(tau - exact[r].LogTau) <= 5 * tauError + 0.02, $"b = {rungs[r]}: E[ln τ] {tau:F4} against {exact[r].LogTau:F4}");
            Assert.True(Math.Abs(scale - exact[r].LogScale) <= 5 * scaleError + 0.02, $"b = {rungs[r]}: E[ln s] {scale:F4} against {exact[r].LogScale:F4}");
        }
    }

    [Theory]
    [InlineData("hinge")]
    [InlineData("legendre")]
    public void Smc_MeetsTheEnumeratedEvidence_WithTheScaleDrawn(string fit)
    {
        (BassModel model, BassConfig config) = Drawn(fit);
        config = config with
        {
            Sampling = config.Sampling with { Chains = 1, MasterSeed = 3, ReportEvidence = true, Smc = new SmcSpec { Particles = 2000, Replicates = 4, MutationSteps = 4 }, Tempering = null },
        };
        BassFit result = Bass.Fit(config, model.Data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        double exact = ExactWithScale(model, [1.0])[0].LogZ;
        BassEvidence evidence = result.Evidence!;
        BassSmcReadout readout = Assert.IsType<BassSmcReadout>(result.Diagnostics.Smc);
        output.WriteLine($"{fit}: ln p(y) {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against the enumeration's {exact:F4}; flags {readout.Flags.Count}");
        Assert.True(Math.Abs(evidence.LogEvidence - exact) <= 4 * evidence.MonteCarloError + 0.01, $"{evidence.LogEvidence:F4} against {exact:F4}");
        Assert.Empty(readout.Flags);
    }

    [Theory]
    [InlineData(20261109)]
    [InlineData(20261110)]
    public void Geweke_TheDrawnScaleSweepPreservesTheJointLaw(int seed)
    {
        if (!BassQualification.Enabled) return;
        const int n = 18, levels = 4, pairs = 24_000;
        BassConfig config = KhaosGPriorTests.Config(scale: TauPrior.Gamma(2, 2), cases: n, sigma2: Sigma2Prior.InverseGamma(3, 2)) with
        {
            RandomIntercept = new BassRandomIntercept(1.0, ScalePrior: new BassGroupScalePrior()),
        };
        BassData template = RandomInterceptTests.Grouped(KhaosBasisTests.UnitData(n, 3, 77), levels);
        BassModel prior = BassModel.Create(config, template);
        var rng = new Xoshiro256PlusPlus(seed);
        var sampler = new BassPriorSampler(prior);
        double[] counts = prior.RandomIntercept!.Counts.ToArray();
        double norm = counts.Sum(x => x * x);
        string[] names = ["ln s", "ln RMS(u)", "M", "ln τ", "ln σ²"];
        double[][] diff = [.. names.Select(_ => new double[pairs])];
        double[] Values(double scale, double[] u, int m, double logTau, double sigma2)
            => [Math.Log(scale), 0.5 * Math.Log(u.Sum(x => x * x) / levels), m, logTau, Math.Log(sigma2)];

        // Independent joint draws, paired with their successors: simulate y from that draw, then update at fixed y.
        // Each paired difference has expectation zero, and pairs are independent. Use log RMS: the Cauchy scale
        // gives the raw group spread no finite prior mean.
        for (int j = 0; j < pairs; j++)
        {
            BassState state = sampler.DrawPrior(rng).State;
            RandomInterceptBlock block = state.Fit.Working!.Whitening!;
            double scale = block.Spec.Scale;
            var start = new BassChain(prior, [], state.LogTau, rng, rng);
            start.Structure.Adopt(state, sampler.LastDesign());
            start.Moves.SetState(state, prior.Target.Evaluate(state));
            BassSweepFit fit = start.Structure.DrawSweepFit(state, rng, prior.Target, power: 0);
            double[] z = [.. counts.Select(_ => Normal.SampleStandard(rng))];
            double along = counts.Zip(z, (a, b) => a * b).Sum() / norm;
            double sd = Math.Sqrt(fit.Sigma2);
            double[] u = [.. z.Select((v, g) => sd * scale * (v - along * counts[g]))];
            double[] y = [.. fit.Fitted.Select((f, i) => f + u[block.GroupOf(i)] + sd * Normal.SampleStandard(rng))];
            double[] before = Values(scale, u, state.Count, state.LogTau, fit.Sigma2);
            var atScale = config with { RandomIntercept = config.RandomIntercept! with { Scale = scale } };
            BassModel model = BassModel.Create(atScale, new BassData(template.Predictors, new ScalarResponse(y), template.Grouping));
            var chain = new BassChain(model, state.Terms, state.LogTau, rng, rng);
            for (int t = 0; t < 4; t++) chain.Step();
            BassDraw draw = chain.Compose(rng, new BassTermTable());
            double[] after = Values(draw.GroupScale!.Value, draw.GroupEffects!, chain.Current.Count, chain.Current.LogTau, draw.Sigma2);
            for (int k = 0; k < names.Length; k++) diff[k][j] = after[k] - before[k];
        }
        double worst = 0.0;
        for (int k = 0; k < names.Length; k++)
        {
            (double mean, double se) = Summary(diff[k]);
            output.WriteLine($"{names[k]}: paired difference {mean:G5} ± {se:G5}, z = {mean / se:F2}");
            worst = Math.Max(worst, Math.Abs(mean / se));
        }
        Assert.True(worst < 4.0, $"Geweke worst |z| {worst:F2}");
    }

    [Theory]
    [InlineData("hinge")]
    [InlineData("legendre")]
    public void Cost_PerSweepWithScaleHeldAndDrawn(string basis)
    {
        if (!BassQualification.Enabled) return;
        (BassModel held, _) = RandomInterceptTests.Fixture(basis);
        (BassModel drawn, _) = Drawn(basis);
        const int sweeps = 20_000;
        double Time(BassModel model, int seed)
        {
            var chain = new BassChain(model, [], 0.0, new Xoshiro256PlusPlus(seed), new Xoshiro256PlusPlus(seed + 1));
            for (int i = 0; i < 1_000; i++) chain.Step();
            long start = Stopwatch.GetTimestamp();
            for (int i = 0; i < sweeps; i++) chain.Step();
            return Stopwatch.GetElapsedTime(start).TotalMicroseconds / sweeps;
        }
        var heldTimes = new double[3];
        var drawnTimes = new double[3];
        for (int k = 0; k < heldTimes.Length; k++)
        {
            if (k % 2 == 0) heldTimes[k] = Time(held, 500 + k);
            drawnTimes[k] = Time(drawn, 500 + k);
            if (k % 2 == 1) heldTimes[k] = Time(held, 500 + k);
        }
        double baseline = heldTimes.Order().ElementAt(1), sampled = drawnTimes.Order().ElementAt(1);
        output.WriteLine($"{basis}, n={held.Data.Response.Cases}, L={held.RandomIntercept!.Levels}, cap={held.MaxBasis}: " +
                         $"median of 3 × {sweeps} sweeps: held {baseline:F2} µs, drawn {sampled:F2} µs, ratio {sampled / baseline:F2}");
    }
}
