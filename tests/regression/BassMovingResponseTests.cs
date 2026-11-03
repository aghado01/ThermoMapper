using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A chain whose working response moves between sweeps (<see cref="BassStructureMoves.Follow"/>).
/// <list type="bullet">
/// <item><b>In place against rebuilt.</b> A design moved in place holds what a design built there holds, through its statistics, its
/// prior's terms, the target's reading of them and the proposals after. A chain moved in place walks a chain rebuilt per snapshot,
/// state by state, draws included.</item>
/// <item><b>At each snapshot.</b> Every structure's log posterior through the moved path against an enumeration, itself checked for
/// detailed balance and against dense algebra.</item>
/// <item><b>The joint.</b> With a small latent whose values hand different working responses, a chain's visits over (structure, latent)
/// against the joint by enumeration, where F moves with the weights too, on a chain and on a ladder.</item>
/// <item><b>Off the support.</b> A move that leaves F is refused, the chain continuing bit for bit as if it was never tried.</item>
/// <item><b>What reads it.</b> The pointwise densities, the conditional means, the training fit and the residual check read each draw's
/// working response. A snapshot continues bit for bit.</item>
/// <item><b>What refuses it.</b> A model whose priors read the working response, and a reading that is not the model's.</item>
/// </list>
/// </summary>
public sealed class BassMovingResponseTests
{
    private readonly ITestOutputHelper _out;
    public BassMovingResponseTests(ITestOutputHelper output) => _out = output;

    private static void Close(double expected, double actual, double relative, string what) => BassWeightedCoreTests.Close(expected, actual, relative, what);

    private static double[] Offset(int n, int snapshot) => BassPriorReferenceTests.Offset(n, snapshot);

    private static double[] Unit(int n) => [.. Enumerable.Repeat(1.0, n)];

    private static BassTerm Hinge(int variable, int sign, int knot) => new(continuous: [new HingeFactor(variable, sign, knot)]);

    private static CoefficientForm Form(bool ridge) => ridge ? new RidgeForm() : new GForm();

    // An observation model whose blocks a test sets once the model exists, since they read its working response.
    private sealed class Latent(Func<BassResponse, WorkingResponse> working, BassSweepSlot slot = BassSweepSlot.Latent) : IBassObservationModel
    {
        public IBlock<ReversibleJumpChain<BassState>>? Block { get; set; }

        public WorkingResponse Working(BassResponse response) => working(response);

        public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot asked) => asked == slot && Block is { } block ? [block] : [];
    }

    // ---- A design moved in place ----

    private static void Same(BassDesign a, BassDesign b, string step)
    {
        Assert.Equal(b.Count, a.Count);
        Assert.Equal(b.Weighting.Count, a.Weighting.Count);
        Close(b.Weighting.Mass, a.Weighting.Mass, 1e-13, $"{step}: Σw");
        Close(b.Weighting.SumLogWeights, a.Weighting.SumLogWeights, 1e-13, $"{step}: Σ ln w");
        Close(b.ResponseMean, a.ResponseMean, 1e-13, $"{step}: ȳ");
        Close(b.PriorLevel, a.PriorLevel, 1e-13, $"{step}: the prior's level");
        Close(b.ResponseSquares, a.ResponseSquares, 1e-12, $"{step}: ‖ỹ‖²");
        Close(b.ResidualSquares, a.ResidualSquares, 1e-10, $"{step}: ‖r‖²");
        Close(b.FittedSquares, a.FittedSquares, 1e-10, $"{step}: Q");
        if (a.Count > 0) Close(b.MinAngle, a.MinAngle, 1e-9, $"{step}: the smallest angle");
        for (int j = 0; j < a.Count; j++) Close(b.Coefficients[j], a.Coefficients[j], 1e-9, $"{step}: β̂_{j}");
        for (int i = 0; i < b.Residual.Length; i++) Close(b.Residual[i], a.Residual[i], 1e-10, $"{step}: r_{i}");
    }

    private static void Proposed(DesignProposal p, DesignProposal q, string step)
    {
        Assert.Equal(q.Refusal, p.Refusal);
        if (q.Refusal is not null) return;
        Assert.Equal(q.Columns, p.Columns);
        Close(q.ResidualSquares, p.ResidualSquares, 1e-10, $"{step}: the proposal's ‖r‖²");
        Close(q.FittedSquares, p.FittedSquares, 1e-10, $"{step}: the proposal's Q");
        Close(q.MinAngle, p.MinAngle, 1e-9, $"{step}: the proposal's angle");
    }

    // Both targets on both fits at several (τ, b): what the model's target reads of a moved fit is what a rebuilt model's reads of its own.
    private static void SameMarginal(BassTarget target, BassDesign moved, BassTarget rebuilt, BassDesign built, string step)
    {
        BassStructureFit a = BassStructureFit.Of(moved, 0.0), b = BassStructureFit.Of(built, 0.0);
        foreach ((double tau, double t) in new[] { (0.2, 1.0), (0.05, 0.6), (3.0, 0.1) })
        {
            GPriorMarginal x = target.Marginal(a, tau, t), y = rebuilt.Marginal(b, tau, t);
            Close(y.S, x.S, 1e-9, $"{step}: S at ({tau}, {t})");
            Close(y.LogDeterminantTerm, x.LogDeterminantTerm, 1e-9, $"{step}: the determinant term at ({tau}, {t})");
            Assert.Equal(y.DegreesOfFreedom, x.DegreesOfFreedom, 1e-12);
            Close(y.LogMarginal(rebuilt.Sigma2), x.LogMarginal(target.Sigma2), 1e-9, $"{step}: the marginal at ({tau}, {t})");
        }
        // The conditional's draw reads what its mean does not: the intercept's variance reads Σw.
        GPriorConditional c = target.CoefficientBlock(moved)(0.2, 1.0), d = rebuilt.CoefficientBlock(built)(0.2, 1.0);
        Assert.Equal(d.Columns, c.Columns);
        var u = new double[c.Columns];
        var v = new double[d.Columns];
        Close(d.Sample(rebuilt.Sigma2, new Xoshiro256PlusPlus(3), v), c.Sample(target.Sigma2, new Xoshiro256PlusPlus(3), u), 1e-9, $"{step}: σ²'s draw");
        for (int j = 0; j < c.Columns; j++)
        {
            Close(d.Mean[j], c.Mean[j], 1e-9, $"{step}: the conditional mean's {j}");
            Close(v[j], u[j], 1e-8, $"{step}: the conditional's draw of {j}");
        }
    }

    [Theory]
    [InlineData(false, "unit", false)]
    [InlineData(true, "unit", false)]
    [InlineData(false, "apart", false)]
    [InlineData(true, "apart", false)]
    [InlineData(false, "apart", true)]
    [InlineData(true, "apart", true)]
    [InlineData(false, "weighted", false)]
    [InlineData(true, "weighted", false)]
    [InlineData(false, "weighted", true)]
    [InlineData(true, "weighted", true)]
    [InlineData(false, "frozen", false)]
    [InlineData(false, "frozen", true)]
    [InlineData(true, "frozen", true)]
    public void AMovedDesign_IsADesignBuiltThere(bool ridge, string reading, bool reweights)
    {
        BassData data = BassSweepTests.Data(large: true);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        WeightedPriorGram gram = reading switch { "weighted" => WeightedPriorGram.Weighted, "frozen" => WeightedPriorGram.Frozen, _ => WeightedPriorGram.Unweighted };
        BassPriorReference reference = reading == "frozen" ? BassPriorReference.FrozenGram(y, BassWeightedCoreTests.Weights(n, seed: 9)) : BassPriorReference.Moments(y);
        double[]? start = reading == "unit" ? null : BassWeightedCoreTests.Weights(n, seed: 21);
        BassConfig config = BassWeightedCoreTests.Config(Form(ridge), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2, maxBasis: 6);
        BassModel Reading(double[] offset, double[]? w) => BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, offset, w, priorGram: gram, reference: reference)));
        BassModel model = Reading(Offset(n, 0), start);
        BassDesign a = model.NewDesign();
        foreach (BassTerm term in BassWeightedCoreTests.Terms[..3])
        {
            Assert.True(a.ProposeBirth(term).Admitted);
            a.Accept();
        }
        BassTerm other = Hinge(0, -1, 3);
        for (int t = 1; t <= 3; t++)
        {
            double[] offset = Offset(n, t);
            double[]? w = reweights ? BassWeightedCoreTests.Weights(n, seed: 30 + t) : start;
            WorkingResponse working = model.Working.Moved(offset, reweights ? w : default);
            WorkingResponse before = a.Working;
            double rss = a.ResidualSquares;
            DesignProposal p = a.ProposeFollow(working);
            Assert.True(p.Admitted);
            Assert.Same(before, a.Working);   // nothing moves until it is accepted
            Assert.Equal(rss, a.ResidualSquares);
            a.Accept();
            Assert.Same(working, a.Working);

            BassModel rebuilt = Reading(offset, w);
            BassDesign b = rebuilt.NewDesign();
            Assert.Null(b.Load(a.Terms));
            Same(a, b, $"snapshot {t}");
            Close(b.ResidualSquares, p.ResidualSquares, 1e-10, "the proposal's ‖r‖²");
            Close(b.FittedSquares, p.FittedSquares, 1e-10, "the proposal's Q");
            SameMarginal(model.Target, a, rebuilt.Target, b, $"snapshot {t}");
            BassDesign built = model.NewDesign(working);
            Assert.Null(built.Load(a.Terms));
            Same(a, built, $"snapshot {t}, on the same working response");
            if (reweights)
            {
                // The weights' path is a design built there and loaded, to the bit.
                Assert.Equal(built.ResidualSquares, a.ResidualSquares);
                Assert.Equal(built.FittedSquares, a.FittedSquares);
                Assert.True(built.Residual.SequenceEqual(a.Residual));
                Assert.True(built.Coefficients.SequenceEqual(a.Coefficients));
                Assert.Equal(built.MinAngle, a.MinAngle);
            }

            // The moves after it read it.
            Proposed(a.ProposeBirth(BassWeightedCoreTests.Terms[3]), b.ProposeBirth(BassWeightedCoreTests.Terms[3]), $"snapshot {t}: birth");
            a.Accept();
            b.Accept();
            Same(a, b, $"snapshot {t}: birth");
            SameMarginal(model.Target, a, rebuilt.Target, b, $"snapshot {t}: birth");
            Proposed(a.ProposeChange(0, other), b.ProposeChange(0, other), $"snapshot {t}: change");
            BassTerm old = a.Terms[0];
            a.Accept();
            b.Accept();
            Same(a, b, $"snapshot {t}: change");
            Proposed(a.ProposeChange(a.Count - 1, old), b.ProposeChange(b.Count - 1, old), $"snapshot {t}: change back");
            a.Accept();
            b.Accept();
            int newest = a.Terms.ToList().IndexOf(BassWeightedCoreTests.Terms[3]);
            Proposed(a.ProposeDeath(newest), b.ProposeDeath(newest), $"snapshot {t}: death");
            a.Accept();
            b.Accept();
            Same(a, b, $"snapshot {t}: death");
            SameMarginal(model.Target, a, rebuilt.Target, b, $"snapshot {t}: death");
        }
    }

    // A direct response at unit weights keeps its Khatri–Rao pair, and under weights its stacked columns.
    [Theory]
    [InlineData(false, false)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void AMovedFunctionalDesign_IsADesignBuiltThere(bool weighted, bool reweights)
    {
        BassData data = BassKernelEnumeration.FunctionalData();
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        BassPriorReference reference = BassPriorReference.Moments(y);
        double[]? start = weighted ? BassWeightedCoreTests.Weights(n, seed: 21) : null;
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2) with
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 5, Support = new SupportRule(SupportKind.Count, 1.0, 1.0) },
        };
        BassModel Reading(double[] offset, double[]? w) => BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, offset, w, reference: reference)));
        BassModel model = Reading(Offset(n, 0), start);
        BassTerm[] terms =
        [
            new(continuous: [new HingeFactor(0, 1, 0)]),
            new(continuous: [new HingeFactor(0, 1, 1)], functional: [new HingeFactor(0, 1, 1)]),
            new(functional: [new HingeFactor(1, 1, 0)]),
        ];
        BassDesign a = model.NewDesign();
        foreach (BassTerm term in terms)
        {
            Assert.True(a.ProposeBirth(term).Admitted, $"{term}");
            a.Accept();
        }
        BassTerm next = new(continuous: [new HingeFactor(0, -1, 2)], functional: [new HingeFactor(1, -1, 1)]);
        for (int t = 1; t <= 3; t++)
        {
            double[] offset = Offset(n, t);
            double[]? w = reweights ? BassWeightedCoreTests.Weights(n, seed: 30 + t) : start;
            Assert.True(a.ProposeFollow(model.Working.Moved(offset, reweights ? w : default)).Admitted);
            a.Accept();
            BassModel rebuilt = Reading(offset, w);
            BassDesign b = rebuilt.NewDesign();
            Assert.Null(b.Load(a.Terms));
            Same(a, b, $"snapshot {t}");
            SameMarginal(model.Target, a, rebuilt.Target, b, $"snapshot {t}");
            Proposed(a.ProposeBirth(next), b.ProposeBirth(next), $"snapshot {t}: birth");
            a.Accept();
            b.Accept();
            Same(a, b, $"snapshot {t}: birth");
            Proposed(a.ProposeDeath(a.Count - 1), b.ProposeDeath(b.Count - 1), $"snapshot {t}: death");
            a.Accept();
            b.Accept();
            Same(a, b, $"snapshot {t}: death");
        }
    }

    // ---- A chain moved in place against a chain rebuilt per snapshot ----

    // State by state, with each snapshot's draws. Under a change of weights the design's arithmetic is a rebuilt design's own; under an
    // offset alone the factor carries its incremental updates, so the two differ only in rounding.
    [Theory]
    [InlineData("moments", false, true)]
    [InlineData("moments", true, true)]
    [InlineData("moments", false, false)]
    [InlineData("moments", true, false)]
    [InlineData("robust", false, true)]
    [InlineData("frozen", false, true)]
    [InlineData("frozen", true, true)]
    [InlineData("weighted", false, true)]
    [InlineData("weighted", true, true)]
    [InlineData("invariant", false, true)]
    [InlineData("invariant", true, false)]
    [InlineData("hierarchical", false, true)]
    [InlineData("hierarchical", true, true)]
    [InlineData("zero", false, true)]
    [InlineData("held", false, true)]
    public void AChainMovedInPlace_WalksAChainRebuiltPerSnapshot(string route, bool ridge, bool reweights)
    {
        BassData data = BassWeightedCoreTests.Scattered(60, 11);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        (BassConfig config, Func<double[], BassPriorReference?> rule, WeightedPriorGram gram, _) = BassPriorReferenceTests.Route(route, ridge, n);
        config = config with { Coefficients = new BassCoefficientPrior(config.Coefficients.Form, TauPrior.UnitInformationGamma()) };
        BassPriorReference? reference = rule(y);
        double[] Weights(int t) => BassWeightedCoreTests.Weights(n, seed: 30 + (reweights ? t : 0));
        BassModel Reading(int t) => BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, Offset(n, t), Weights(t), priorGram: gram, reference: reference)));
        BassModel model = Reading(0);
        var moveA = new Xoshiro256PlusPlus(5);
        var blockA = new Xoshiro256PlusPlus(6);
        var moveB = new Xoshiro256PlusPlus(5);
        var blockB = new Xoshiro256PlusPlus(6);
        var inPlace = new BassChain(model, [], Math.Log(0.05), moveA, blockA);
        IReadOnlyList<BassTerm> terms = [];
        double logTau = Math.Log(0.05), logRate = double.NaN;
        int changes = 0;
        BassTerm[] last = [];
        var table = new BassTermTable();
        var builtTable = new BassTermTable();
        for (int t = 0; t < 8; t++)
        {
            if (t > 0) Assert.Null(inPlace.Follow(model.Working.Moved(Offset(n, t), reweights ? Weights(t) : default)));
            BassModel rebuilt = Reading(t);
            var chain = new BassChain(rebuilt, terms, logTau, moveB, blockB, logRate: logRate);
            for (int s = 0; s < 40; s++)
            {
                BassState x = inPlace.Step(), z = chain.Step();
                Assert.Equal(z.Terms, x.Terms);
                Close(z.LogTau, x.LogTau, 1e-9, $"snapshot {t}, sweep {s}: ln τ");
                if (model.Target.Sigma2Rate is not null) Close(z.LogRate, x.LogRate, 1e-9, $"snapshot {t}, sweep {s}: ln r");
                Close(z.Fit.ResidualSquares, x.Fit.ResidualSquares, 1e-9, $"snapshot {t}, sweep {s}: ‖r‖²");
                Close(rebuilt.Target.Evaluate(z).LogMarginal, model.Target.Evaluate(x).LogMarginal, 1e-9, $"snapshot {t}, sweep {s}: the marginal");
                if (!x.Terms.SequenceEqual(last)) changes++;
                last = [.. x.Terms];
            }
            BassDraw p = inPlace.Compose(new Xoshiro256PlusPlus(70 + t), table), q = chain.Compose(new Xoshiro256PlusPlus(70 + t), builtTable);
            Close(q.Sigma2, p.Sigma2, 1e-8, $"snapshot {t}: σ²");
            Assert.Equal(q.Coefficients.Length, p.Coefficients.Length);
            for (int j = 0; j < p.Coefficients.Length; j++) Close(q.Coefficients[j], p.Coefficients[j], 1e-8, $"snapshot {t}: coefficient {j}");
            terms = chain.Current.Terms;
            logTau = chain.Current.LogTau;
            logRate = chain.Current.LogRate;
        }
        Assert.True(changes > 15, $"The chain changed structure {changes} times.");
    }

    // ---- Enumeration at each snapshot ----

    // Every structure of F at the snapshot through the moved path, against an enumeration checked for detailed balance and, on level
    // centring, against dense algebra.
    [Theory]
    [InlineData("moments", false)]
    [InlineData("moments", true)]
    [InlineData("frozen", true)]
    [InlineData("weighted", true)]
    [InlineData("invariant", true)]
    [InlineData("zero", true)]
    public void Enumeration_AtEachSnapshot(string route, bool reweights)
    {
        const double tau = 0.1;
        double logTau = Math.Log(tau);
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        double[] w0 = BassWeightedCoreTests.Weights(n, seed: 22);
        WeightedPriorGram gram = route switch { "frozen" => WeightedPriorGram.Frozen, "weighted" => WeightedPriorGram.Weighted, _ => WeightedPriorGram.Unweighted };
        BassPriorReference? reference = route switch
        {
            "frozen" => BassPriorReference.FrozenGram(y, w0),
            "moments" or "weighted" => BassPriorReference.Moments(y),
            _ => null,
        };
        BassConfig config = BassWeightedCoreTests.EnumerationConfig(route switch { "invariant" => Sigma2Prior.Jeffreys, "zero" => Sigma2Prior.InverseGamma(2.0, 0.3), _ => Sigma2Prior.ScaledToResponse(0.5) });
        if (route is "invariant" or "zero") config = config with { Intercept = route == "zero" ? GPriorCentre.Zero : GPriorCentre.Flat };
        double[] Weights(int t) => BassWeightedCoreTests.Weights(n, seed: 21 + (reweights ? t : 0));
        BassModel Reading(int t) => BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, Offset(n, t), Weights(t), priorGram: gram, reference: reference)));
        BassModel model = Reading(0);
        var moves = new BassStructureMoves(model, [], logTau);
        for (int t = 1; t <= 2; t++)
        {
            double[] offset = Offset(n, t), w = Weights(t);
            BassModel rebuilt = Reading(t);
            BassWeightedCoreTests.CheckKernel(rebuilt, $"{route}, snapshot {t}", _out);
            BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(rebuilt);
            if (reference is not null)
                BassWeightedCoreTests.DensePosterior(space, [.. y.Zip(offset, (a, b) => a - b)], w, gram == WeightedPriorGram.Frozen ? w0 : gram == WeightedPriorGram.Weighted ? w : Unit(n),
                                                     reference.Level, tau, null, rebuilt.Target.Sigma2, w.Sum(Math.Log));
            WorkingResponse working = model.Working.Moved(offset, reweights ? w : default);
            int compared = 0;
            foreach (BassKernelEnumeration.State s in space.States)
            {
                if (BassMixing.Loaded(model, s.Terms, logTau, out _) is not { } at) continue;
                BassState moved = moves.ProposeFollow(at, working);
                BassEvaluation evaluation = model.Target.Evaluate(moved);
                Assert.Null(evaluation.Refusal);
                double logPrior = model.Prior.LogCollectionWeight(s.Terms.Length) + s.Terms.Sum(model.Prior.LogPrior);
                Close(s.LogPi, logPrior + evaluation.LogLikelihood(1.0), 1e-9, $"snapshot {t}, {string.Join(", ", s.Terms.Select(v => v.ToString()))}");
                Close(evaluation.LogPrior, model.Prior.LogCollectionWeight(s.Terms.Length) + s.Terms.Sum(model.Prior.LogPrior), 1e-12, "the prior");
                compared++;
            }
            Assert.Equal(space.States.Count, compared);
            _out.WriteLine($"{route}, snapshot {t}: {compared} structures through the moved path");
        }
    }

    // ---- The joint over the structure and a latent ----

    // A latent z with a working response per value, moved by a uniform independence proposal, or by a Gibbs step over its values: the
    // block F's off-support rule reads, since a value that puts the structure outside F has posterior 0.
    private sealed class LatentBlock(WorkingResponse[] workings, double[] prior, bool gibbs) : IBlock<ReversibleJumpChain<BassState>>
    {
        public string Name => "latent";

        public bool Update(ReversibleJumpChain<BassState> chain, Xoshiro256PlusPlus rng)
        {
            BassStructureMoves moves = BassStructureMoves.Of(chain);
            if (gibbs)
            {
                var logs = new double[workings.Length];
                for (int z = 0; z < logs.Length; z++)
                    logs[z] = Math.Log(prior[z]) + chain.Temperature.LogDensity(chain.Target.Evaluate(moves.ProposeFollow(chain.Current, workings[z])));
                double top = logs.Max(), total = logs.Sum(v => Math.Exp(v - top)), u = rng.NextDouble() * total, sum = 0.0;
                int drawn = logs.Length - 1;
                for (int z = 0; z < logs.Length; z++)
                {
                    sum += Math.Exp(logs[z] - top);
                    if (u < sum)
                    {
                        drawn = z;
                        break;
                    }
                }
                return BassStructureMoves.Follow(chain, workings[drawn]) is null;
            }
            int current = Array.IndexOf(workings, chain.Current.Fit.Working), proposed = rng.NextInt(workings.Length);
            BassState candidate = moves.ProposeFollow(chain.Current, workings[proposed]);
            ITargetEvaluation evaluation = chain.Target.Evaluate(candidate);
            double ratio = Math.Log(prior[proposed]) - Math.Log(prior[current]) + chain.Temperature.LogDensity(evaluation) - chain.CurrentLogDensity;
            if (Math.Log(rng.NextDouble()) >= ratio) return false;
            chain.SetState(candidate, evaluation);
            return true;
        }
    }

    // The joint over (z, structure), π ∝ p(z)·π(S)·p(y | S, τ, z) on F(z), by an enumeration per value. Its values are comparable across
    // z: what they leave out, under an improper prior, is free of z. A z's prior is what makes its posterior about even, so each value
    // holds mass and a target that misreads one shows: at unit prior they run from 0.95 down to 0.001.
    private static (double[] Prior, double[] Joint, int[] Starts, BassKernelEnumeration.Space[] Spaces) Joint(BassModel[] rebuilt)
    {
        BassKernelEnumeration.Space[] spaces = [.. rebuilt.Select(m => BassKernelEnumeration.Enumerate(m))];
        double[] evidence = [.. spaces.Select(space =>
        {
            double largest = space.States.Max(x => x.LogPi);
            return largest + Math.Log(space.States.Sum(x => Math.Exp(x.LogPi - largest)));
        })];
        double lowest = evidence.Min();
        double[] flat = [.. evidence.Select(v => Math.Exp(lowest - v))];
        double sum = flat.Sum();
        double[] prior = [.. flat.Select(v => v / sum)];
        var starts = new int[spaces.Length];
        var logs = new List<double>();
        for (int z = 0; z < spaces.Length; z++)
        {
            starts[z] = logs.Count;
            logs.AddRange(spaces[z].States.Select(x => Math.Log(prior[z]) + x.LogPi));
        }
        double top = logs.Max();
        double[] p = [.. logs.Select(v => Math.Exp(v - top))];
        double total = p.Sum();
        return (prior, [.. p.Select(v => v / total)], starts, spaces);
    }

    // The joint's mass on structures in F at their z and outside it at another: what a move would put off the support.
    private static double MassLeavingF(double[] joint, int[] starts, BassKernelEnumeration.Space[] spaces)
    {
        double mass = 0.0;
        for (int z = 0; z < spaces.Length; z++)
            for (int i = 0; i < spaces[z].States.Count; i++)
                if (Enumerable.Range(0, spaces.Length).Any(other => spaces[other].Find(spaces[z].States[i].Terms) < 0)) mass += joint[starts[z] + i];
        return mass;
    }

    private void CheckJointVisits(double[] joint, int[] starts, BassKernelEnumeration.Space[] spaces, WorkingResponse[] workings, string label,
                                  int chains, int sweeps, Func<int, Func<BassState>> start)
    {
        const int burn = 2_000;
        var visits = new double[joint.Length];
        for (int c = 0; c < chains; c++)
        {
            Func<BassState> step = start(c);
            for (int s = 0; s < burn; s++) step();
            for (int s = 0; s < sweeps; s++)
            {
                BassState state = step();
                int z = Array.IndexOf(workings, state.Fit.Working);
                Assert.True(z >= 0, "The state names a working response the block never handed it.");
                int found = spaces[z].Find(state.Terms);
                Assert.True(found >= 0, $"The chain holds a structure outside F at z = {z}.");
                visits[starts[z] + found] += 1.0;
            }
        }
        double total = (double)chains * sweeps;
        for (int z = 0; z < spaces.Length; z++)
        {
            int count = spaces[z].States.Count;
            double pz = 0.0, vz = 0.0;
            for (int i = 0; i < count; i++)
            {
                pz += joint[starts[z] + i];
                vz += visits[starts[z] + i] / total;
            }
            _out.WriteLine($"{label}: z = {z}: {count} structures, π {pz:F4}, visits {vz:F4}");
            Assert.True(Math.Abs(pz - vz) <= 0.03, $"{label}: z = {z}: π {pz:F4}, visits {vz:F4}");
        }
        int cap = spaces.Max(s => s.States.Max(x => x.Terms.Length));
        for (int m = 0; m <= cap; m++)
        {
            double pm = 0.0, vm = 0.0;
            for (int z = 0; z < spaces.Length; z++)
                for (int i = 0; i < spaces[z].States.Count; i++)
                    if (spaces[z].States[i].Terms.Length == m)
                    {
                        pm += joint[starts[z] + i];
                        vm += visits[starts[z] + i] / total;
                    }
            Assert.True(Math.Abs(pm - vm) <= 0.03, $"{label}: M = {m}: π {pm:F4}, visits {vm:F4}");
        }
        foreach (int i in Enumerable.Range(0, joint.Length).OrderByDescending(i => joint[i]).Take(10))
            Assert.True(Math.Abs(joint[i] - visits[i] / total) <= 0.02, $"{label}: state {i}: π {joint[i]:F4}, visits {visits[i] / total:F4}");
    }

    // The readings a joint runs on: the starting one, then one per latent value.
    private static (BassConfig Config, Func<int, double[]> Offsets, Func<int, double[]?> Weights, WeightedPriorGram Gram, BassPriorReference? Reference) JointReading(string reading, double[] y, bool ridge)
    {
        int n = y.Length;
        double[] Uneven(int z)   // one row in three loud and two quiet, their product 1
        {
            var w = new double[n];
            for (int i = 0; i < n; i++) w[i] = z == 0 ? 1.0 : (i + z) % 3 == 0 ? 4.0 : 0.5;
            return w;
        }
        double[] Offsets(int z) => [.. Offset(n, z).Select(o => 0.6 * o)];
        double[] Weights(int z) => BassWeightedCoreTests.Weights(n, seed: 41 + z);
        BassConfig Enumerated(Sigma2Prior sigma2, GPriorCentre intercept = GPriorCentre.Mean) => BassWeightedCoreTests.EnumerationConfig(sigma2) with
        {
            Intercept = intercept,
            Coefficients = new BassCoefficientPrior(Form(ridge), TauPrior.Held(0.1)),
        };
        BassConfig scaled = Enumerated(Sigma2Prior.ScaledToResponse(0.5));
        return reading switch
        {
            "offset" => (scaled, Offsets, _ => null, WeightedPriorGram.Unweighted, BassPriorReference.Moments(y)),
            "apart" => (scaled, Offsets, Weights, WeightedPriorGram.Unweighted, BassPriorReference.Moments(y)),
            "weighted" => (scaled, Offsets, Weights, WeightedPriorGram.Weighted, BassPriorReference.Moments(y)),
            "frozen" => (scaled, Offsets, Weights, WeightedPriorGram.Frozen, BassPriorReference.FrozenGram(y, BassWeightedCoreTests.Weights(n, seed: 22))),
            "flat" => (Enumerated(Sigma2Prior.InverseGamma(2.0, 0.3), GPriorCentre.Flat), Offsets, Weights, WeightedPriorGram.Unweighted, null),
            "zero" => (Enumerated(Sigma2Prior.InverseGamma(2.0, 0.3), GPriorCentre.Zero), Offsets, Weights, WeightedPriorGram.Unweighted, null),
            // F moves with the weights: the gate high, and the weights uneven. Far more uneven and the values of z become modes apart,
            // which a chain crosses too rarely to read here (at 100 and 0.1, each chain's z = 1 reads 0.28 to 0.35 over 5,000,000 sweeps).
            // The floor sits on no angle this grid gives: at 0.6 some are exact ties, which rounding decides by the columns' order.
            "gate" => (scaled with { Structure = scaled.Structure with { AngleFloor = 0.61 } }, Offsets, Uneven, WeightedPriorGram.Unweighted, BassPriorReference.Moments(y)),
            _ => throw new ArgumentOutOfRangeException(nameof(reading)),
        };
    }

    [Theory]
    [InlineData("offset", false, false)]
    [InlineData("offset", false, true)]
    [InlineData("apart", false, false)]
    [InlineData("apart", false, true)]
    [InlineData("apart", true, false)]
    [InlineData("weighted", false, true)]
    [InlineData("weighted", true, true)]
    [InlineData("frozen", false, false)]
    [InlineData("flat", false, false)]
    [InlineData("zero", false, true)]
    [InlineData("gate", false, false)]
    [InlineData("gate", false, true)]
    public void TheJoint_OverTheStructureAndALatent(string reading, bool ridge, bool gibbs)
    {
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        (BassConfig config, Func<int, double[]> offsets, Func<int, double[]?> weights, WeightedPriorGram gram, BassPriorReference? reference) = JointReading(reading, y, ridge);
        Func<BassResponse, WorkingResponse> At(int z) => r => new WorkingResponse(r, offsets(z), weights(z), priorGram: gram, reference: reference);
        var observation = new Latent(At(0));
        BassModel model = BassModel.Create(config, data, observation);
        WorkingResponse[] workings = [model.Working, .. Enumerable.Range(1, 2).Select(z => model.Working.Moved(offsets(z), weights(z)))];
        BassModel[] rebuilt = [.. Enumerable.Range(0, workings.Length).Select(z => BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(At(z))))];
        (double[] prior, double[] joint, int[] starts, BassKernelEnumeration.Space[] spaces) = Joint(rebuilt);
        observation.Block = new LatentBlock(workings, prior, gibbs);
        double leaving = MassLeavingF(joint, starts, spaces);
        _out.WriteLine($"{reading}: F holds {string.Join(", ", spaces.Select(x => x.States.Count))} structures; mass on those a move would leave F {leaving:F4}");
        if (reading == "gate") Assert.True(leaving > 0.03, $"The gate case should put mass where F moves, not {leaving:F4}.");
        CheckJointVisits(joint, starts, spaces, workings, $"{reading}, {(ridge ? "ridge" : "g")}, {(gibbs ? "Gibbs" : "independence")}", 4, 30_000, c =>
        {
            var chain = new BassChain(model, [], Math.Log(0.1), new Xoshiro256PlusPlus(300 + c), new Xoshiro256PlusPlus(310 + c));
            return chain.Step;
        });
    }

    // On a ladder each replica holds its own latent, and a swap hands it over with the state: the cold replica's visits are the joint's.
    [Theory]
    [InlineData("apart")]
    [InlineData("gate")]
    public void TheJoint_OnALadder(string reading)
    {
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        (BassConfig config, Func<int, double[]> offsets, Func<int, double[]?> weights, WeightedPriorGram gram, BassPriorReference? reference) = JointReading(reading, y, false);
        Func<BassResponse, WorkingResponse> At(int z) => r => new WorkingResponse(r, offsets(z), weights(z), priorGram: gram, reference: reference);
        var observation = new Latent(At(0));
        BassModel model = BassModel.Create(config, data, observation);
        WorkingResponse[] workings = [model.Working, .. Enumerable.Range(1, 2).Select(z => model.Working.Moved(offsets(z), weights(z)))];
        BassModel[] rebuilt = [.. Enumerable.Range(0, workings.Length).Select(z => BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(At(z))))];
        (double[] prior, double[] joint, int[] starts, BassKernelEnumeration.Space[] spaces) = Joint(rebuilt);
        observation.Block = new LatentBlock(workings, prior, false);
        long swaps = 0;
        CheckJointVisits(joint, starts, spaces, workings, $"a ladder, {reading}", 3, 20_000, c =>
        {
            var ladder = new BassLadder(model, new TemperingSpec { Levels = 3, HotEnd = 0.2 }, [], Math.Log(0.1), 500 + c);
            return () =>
            {
                ladder.Step();
                swaps = ladder.Tempering.PairStats.Sum(p => p.Accepts);
                return ladder.Current;
            };
        });
        Assert.True(swaps > 100, $"The ladder swapped {swaps} times.");
    }

    // ---- Off the support ----

    // Exactly: a structure of F at one value of z moved to another is refused where the enumeration leaves it out there, and otherwise
    // reads the enumeration's log posterior.
    [Fact]
    public void AMove_IsRefusedExactlyWhereFLeavesTheStructure()
    {
        const double tau = 0.1;
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        (BassConfig config, Func<int, double[]> offsets, Func<int, double[]?> weights, WeightedPriorGram gram, BassPriorReference? reference) = JointReading("gate", y, false);
        Func<BassResponse, WorkingResponse> At(int z) => r => new WorkingResponse(r, offsets(z), weights(z), priorGram: gram, reference: reference);
        BassModel model = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(At(0)));
        WorkingResponse[] workings = [model.Working, .. Enumerable.Range(1, 2).Select(z => model.Working.Moved(offsets(z), weights(z)))];
        BassKernelEnumeration.Space[] spaces = [.. Enumerable.Range(0, 3).Select(z => BassKernelEnumeration.Enumerate(BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(At(z)))))];
        int refused = 0, admitted = 0;
        for (int z = 0; z < 3; z++)
            foreach (BassKernelEnumeration.State s in spaces[z].States)
                for (int other = 0; other < 3; other++)
                {
                    BassDesign design = model.NewDesign(workings[z]);
                    Assert.Null(design.Load(s.Terms));
                    DesignProposal proposal = design.ProposeFollow(workings[other]);
                    int found = spaces[other].Find(s.Terms);
                    Assert.Equal(found < 0, proposal.Refusal is not null);
                    if (found < 0)
                    {
                        refused++;
                        Assert.Equal(BassRefusal.Numerical, proposal.Refusal);
                        continue;
                    }
                    admitted++;
                    double logTermPrior = s.Terms.Sum(model.Prior.LogPrior);
                    BassEvaluation evaluation = model.Target.Evaluate(new BassStructureFit(proposal.Columns, logTermPrior, proposal.ResidualSquares, proposal.FittedSquares, null, proposal.Ridge, workings[other]));
                    Close(spaces[other].States[found].LogPi, evaluation.LogPrior + evaluation.LogLikelihood(1.0), 1e-9, $"z = {z} to {other}");
                    Assert.Equal(tau, evaluation.Tau);
                }
        _out.WriteLine($"{admitted} moves admitted, {refused} refused");
        Assert.True(refused > 100, $"F should move with the weights: {refused} refused.");
    }

    private static (BassData Data, double[] Silencing, BassTerm First, BassTerm Second) NearlyOneColumn()
    {
        BassData data = BassSweepTests.Data(large: true);
        ReadOnlySpan<double> x0 = ((ContinuousPredictor)data.Predictors[0]).Values.Span;
        var w = new double[x0.Length];
        for (int i = 0; i < w.Length; i++) w[i] = x0[i] < 0.5 ? 1e-22 : 1.0;   // the two hinges below differ by a constant from ½ up
        return (data, w, Hinge(0, 1, 1), Hinge(0, 1, 2));
    }

    // A move that puts the structure below s_min is refused, whatever the prior's Gram, and the chain is as if it was never tried.
    [Theory]
    [InlineData(WeightedPriorGram.Unweighted)]
    [InlineData(WeightedPriorGram.Weighted)]
    [InlineData(WeightedPriorGram.Frozen)]
    public void AMoveOutsideF_IsRefused_AndTheChainContinuesAsIfNeverTried(WeightedPriorGram gram)
    {
        (BassData data, double[] silencing, BassTerm first, BassTerm second) = NearlyOneColumn();
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        BassPriorReference reference = gram == WeightedPriorGram.Frozen ? BassPriorReference.FrozenGram(y, BassWeightedCoreTests.Weights(n, seed: 9)) : BassPriorReference.Moments(y);
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2) with
        {
            Coefficients = new BassCoefficientPrior(new GForm(), TauPrior.UnitInformationGamma()),
        };
        BassModel model = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: Unit(n), priorGram: gram, reference: reference)));
        WorkingResponse silenced = model.Working.Moved(default, silencing);
        BassChain Chain() => new(model, [first, second], Math.Log(0.05), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        BassChain tried = Chain(), plain = Chain();

        BassState before = tried.Current;
        BassDesign design = tried.Structure.DesignAt(before);
        double rss = design.ResidualSquares;
        Assert.Equal(BassRefusal.Numerical, tried.Follow(silenced));
        Assert.Same(before, tried.Current);
        Assert.Same(model.Working, design.Working);
        Assert.Equal(rss, design.ResidualSquares);
        BassState refused = tried.Structure.ProposeFollow(before, silenced);
        Assert.Equal(BassRefusal.Numerical, refused.Fit.Refusal);
        BassEvaluation evaluation = model.Target.Evaluate(refused);
        Assert.True(evaluation.OffSupport);
        Assert.True(double.IsNegativeInfinity(Temperature.Posterior.LogDensity(evaluation)));

        int attempts = 0;
        for (int s = 0; s < 400; s++)
        {
            BassState x = tried.Step(), z = plain.Step();
            Assert.Equal(z.Terms, x.Terms);
            Assert.Equal(z.LogTau, x.LogTau);
            if (x.Terms.Contains(first) && x.Terms.Contains(second))
            {
                Assert.Equal(BassRefusal.Numerical, tried.Follow(silenced));
                attempts++;
            }
        }
        Assert.True(attempts > 0, "The chain should hold both columns for a while.");
    }

    [Fact]
    public void AVarianceInflationCap_ReadsTheMoveToo()
    {
        (BassData data, double[] silencing, BassTerm first, BassTerm second) = NearlyOneColumn();
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        double[] partly = [.. silencing.Select(v => v == 1.0 ? 1.0 : 1e-4)];
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2);
        BassModel Capped(double cap) => BassModel.Create(config with { Structure = config.Structure with { MaxVarianceInflation = cap } }, data,
                                                         new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: Unit(n), reference: BassPriorReference.Moments(y))));
        BassModel loose = Capped(1e12);
        var chain = new BassChain(loose, [first, second], Math.Log(0.2), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        double unit = chain.Structure.DesignAt(chain.Current).MinAngle;
        Assert.Null(chain.Follow(loose.Working.Moved(default, partly)));
        double moved = chain.Structure.DesignAt(chain.Current).MinAngle;
        Assert.True(moved < unit, $"The weights should close the angle: {unit} to {moved}.");
        BassModel capped = Capped(1.0 / Math.Sqrt(unit * moved));   // the cap clears the start and not the move
        var tight = new BassChain(capped, [first, second], Math.Log(0.2), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        BassState before = tight.Current;
        Assert.Equal(BassRefusal.VarianceInflation, tight.Follow(capped.Working.Moved(default, partly)));
        Assert.Same(before, tight.Current);
    }

    // The exact-fit floor reads the working response a state is on, as a model built there does, and a fit that is exact there is refused.
    [Theory]
    [InlineData("invariant", false)]
    [InlineData("invariant", true)]
    [InlineData("zero", false)]
    [InlineData("zero", true)]
    [InlineData("reference", false)]
    [InlineData("reference", true)]
    public void TheExactFitFloor_FollowsTheWorkingResponse(string route, bool weighted)
    {
        const double multiple = 1e3, epsilon = 2.220446049250313e-16;
        BassData data = BassSweepTests.Data(large: true);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        BassPriorReference? reference = route == "reference" ? BassPriorReference.Moments(y) : null;
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), route == "zero" ? GPriorCentre.Zero : route == "invariant" ? GPriorCentre.Flat : GPriorCentre.Mean, Sigma2Prior.Jeffreys, 0.2);
        double[]? w = weighted ? BassWeightedCoreTests.Weights(n, seed: 21) : null;
        BassModel Reading(double[] offset, double[]? weights) => BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, offset, weights, reference: reference)));
        BassModel model = Reading(Offset(n, 0), w);
        BassTerm term = Hinge(0, 1, 1);
        double[] offset = Offset(n, 2);
        double[]? moved = weighted ? BassWeightedCoreTests.Weights(n, seed: 31) : null;
        WorkingResponse working = model.Working.Moved(offset, moved);
        BassModel rebuilt = Reading(offset, moved);
        double floor = multiple * epsilon * (multiple * epsilon) * rebuilt.Resolved.Weighting.Mass * rebuilt.Resolved.WeightedMeanSquare;
        Assert.True(floor > 0.0);
        foreach (double factor in new[] { 0.999, 1.001 })
        {
            BassRefusal? expected = factor < 1.0 ? BassRefusal.ExactFit : null;
            Assert.Equal(expected, rebuilt.Target.Refusal(new BassStructureFit(0, 0.0, factor * floor, 0.0)));
            Assert.Equal(expected, model.Target.Refusal(new BassStructureFit(0, 0.0, factor * floor, 0.0, Working: working)));
        }

        // A working response this structure fits exactly: y − offset = 1 + 2z.
        var chain = new BassChain(model, [term], Math.Log(0.2), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        double[] column = new double[n];
        model.Basis.CaseColumn(term, model.Basis.Inputs.Cases, column);
        double[] exactly = [.. y.Select((v, i) => v - 1.0 - 2.0 * column[i])];
        BassState before = chain.Current;
        Assert.Equal(BassRefusal.ExactFit, chain.Follow(model.Working.Moved(exactly, moved)));
        Assert.Same(before, chain.Current);
        Assert.Null(chain.Follow(working));
    }

    // ---- What refuses it ----

    public static TheoryData<GPriorCentre, string, bool, bool> Reads => new()
    {
        { GPriorCentre.Mean, "scaled", false, true },
        { GPriorCentre.Flat, "scaled", false, true },
        { GPriorCentre.Zero, "scaled", false, true },
        { GPriorCentre.Mean, "jeffreys", false, true },
        { GPriorCentre.Mean, "declared", false, true },
        { GPriorCentre.Mean, "held", false, true },
        { GPriorCentre.Mean, "hierarchical", false, true },
        { GPriorCentre.Flat, "hierarchical", false, true },
        { GPriorCentre.Flat, "jeffreys", false, false },
        { GPriorCentre.Flat, "held", false, false },
        { GPriorCentre.Flat, "declared", false, false },
        { GPriorCentre.Zero, "jeffreys", false, false },
        { GPriorCentre.Zero, "declared", false, false },
        { GPriorCentre.Zero, "held", false, false },
        { GPriorCentre.Mean, "scaled", true, false },
        { GPriorCentre.Mean, "hierarchical", true, false },
        { GPriorCentre.Mean, "jeffreys", true, false },
        { GPriorCentre.Flat, "scaled", true, false },
    };

    [Theory]
    [MemberData(nameof(Reads))]
    public void AModelWhosePriorsReadTheResponse_RefusesToMove(GPriorCentre intercept, string sigma2, bool stated, bool refused)
    {
        BassData data = BassSweepTests.Data(large: true);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        Sigma2Prior prior = sigma2 switch
        {
            "scaled" => Sigma2Prior.ScaledToResponse(0.5),
            "jeffreys" => Sigma2Prior.Jeffreys,
            "declared" => Sigma2Prior.InverseGamma(2.0, 0.3),
            "held" => Sigma2Prior.Held(0.3),
            _ => Sigma2Prior.HalfCauchy(),
        };
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), intercept, prior, 0.2);
        BassPriorReference? reference = stated ? BassPriorReference.Moments(y) : null;
        BassModel model = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, Offset(n, 0), reference: reference)));
        Assert.Equal(refused, model.PriorsReadResponse);
        var chain = new BassChain(model, [Hinge(0, 1, 1)], Math.Log(0.2), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        WorkingResponse working = model.Working.Moved(Offset(n, 1));
        if (refused)
        {
            var e = Assert.Throws<InvalidOperationException>(() => chain.Follow(working));
            Assert.Contains("BassPriorReference", e.Message);
            return;
        }
        Assert.Null(chain.Follow(working));
        Assert.Same(working, chain.Current.Fit.Working);
    }

    // A chain's working response stays its model's reading: the same response, what a weight states, the prior's Gram and the reference.
    [Fact]
    public void AMove_KeepsTheModelsReading()
    {
        BassData data = BassSweepTests.Data(large: true), other = BassSweepTests.Data(large: true, seed: 6);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        BassPriorReference reference = BassPriorReference.Moments(y);
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2);
        double[] w = BassWeightedCoreTests.Weights(n, seed: 21);
        BassModel weighted = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: w, reference: reference)));
        BassDesign design = weighted.NewDesign();
        Assert.Throws<ArgumentException>(() => design.ProposeFollow(new WorkingResponse(other.Response, weights: w, reference: reference)));
        Assert.Throws<ArgumentException>(() => design.ProposeFollow(new WorkingResponse(data.Response, weights: w, reference: BassPriorReference.Moments(y))));
        Assert.Throws<ArgumentException>(() => design.ProposeFollow(new WorkingResponse(data.Response, weights: w)));
        Assert.Throws<ArgumentException>(() => design.ProposeFollow(new WorkingResponse(data.Response, weights: w, meaning: WeightMeaning.Replication, reference: reference)));
        Assert.Throws<ArgumentException>(() => design.ProposeFollow(new WorkingResponse(data.Response, weights: w, priorGram: WeightedPriorGram.Weighted, reference: reference)));
        Assert.Throws<ArgumentException>(() => design.ProposeFollow(new WorkingResponse(data.Response, reference: reference)));
        Assert.Null(design.ProposeFollow(new WorkingResponse(data.Response, Offset(n, 1), w, reference: reference)).Refusal);

        // A model that holds no weights takes an offset, not weights.
        BassModel unit = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, reference: reference)));
        Assert.Throws<ArgumentException>(() => unit.Working.Moved(Offset(n, 1), w));
        Assert.Throws<ArgumentException>(() => unit.NewDesign().ProposeFollow(new WorkingResponse(data.Response, weights: w, reference: reference)));
        Assert.Null(unit.NewDesign().ProposeFollow(unit.Working.Moved(Offset(n, 1))).Refusal);

        // Weights that count rows stay; an offset beside them moves.
        BassModel counted = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: Unit(n), meaning: WeightMeaning.Replication, reference: reference)));
        Assert.Throws<NotSupportedException>(() => counted.NewDesign().ProposeFollow(counted.Working.Moved(default, w)));
        Assert.Null(counted.NewDesign().ProposeFollow(counted.Working.Moved(Offset(n, 1))).Refusal);

        // Moved keeps the reading, and shares the weights it keeps.
        WorkingResponse moved = weighted.Working.Moved(Offset(n, 1));
        Assert.Same(reference, moved.Reference);
        Assert.Equal(weighted.Working.Meaning, moved.Meaning);
        Assert.Equal(weighted.Working.PriorGram, moved.PriorGram);
        Assert.Same(weighted.Working.Output, moved.Output);
        Assert.Equal(weighted.Working.Weighting, moved.Weighting);
        Assert.True(weighted.Working.Weights == moved.Weights);
        Assert.Throws<ArgumentException>(() => weighted.Working.Moved(new double[n - 1]));
        Assert.Throws<ArgumentException>(() => weighted.Working.Moved(default, new double[n]));
        Assert.Equal(new WorkingResponse(data.Response, Offset(n, 1), w, reference: reference).ResponseMean, moved.ResponseMean);
    }

    // ---- What reads it ----

    // Draws composed across snapshots of one chain, each read on its own working response: a rebuilt model's density of that snapshot's draws.
    [Theory]
    [InlineData(PointwiseContent.RetainedParameters, false, true)]
    [InlineData(PointwiseContent.RetainedParameters, true, true)]
    [InlineData(PointwiseContent.CollapsedPredictive, false, true)]
    [InlineData(PointwiseContent.CollapsedPredictive, true, false)]
    public void ThePointwiseDensities_ReadEachDrawsWorkingResponse(PointwiseContent content, bool ridge, bool reweights)
    {
        BassData data = BassWeightedCoreTests.Scattered(40, 11);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        BassPriorReference reference = BassPriorReference.Moments(y);
        BassConfig config = BassWeightedCoreTests.Config(Form(ridge), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2);
        double[] Weights(int t) => BassWeightedCoreTests.Weights(n, seed: 30 + (reweights ? t : 0));
        BassModel Reading(int t) => BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, Offset(n, t), Weights(t), reference: reference)));
        BassModel model = Reading(0);
        var chain = new BassChain(model, [], Math.Log(0.2), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        var table = new BassTermTable();
        var rng = new Xoshiro256PlusPlus(7);
        var draws = new List<BassDraw>();
        var bySnapshot = new List<List<BassDraw>>();
        for (int t = 0; t < 3; t++)
        {
            if (t > 0) Assert.Null(chain.Follow(model.Working.Moved(Offset(n, t), reweights ? Weights(t) : default)));
            var these = new List<BassDraw>();
            for (int k = 0; k < 6; k++)
            {
                for (int s = 0; s < 10; s++) chain.Step();
                BassDraw draw = chain.Compose(rng, table);
                Assert.Same(chain.Current.Fit.Working, draw.Working);
                these.Add(draw);
                draws.Add(draw);
            }
            bySnapshot.Add(these);
        }
        PointwiseLogDensity all = BassPointwise.LogDensity(model, table, [draws], content: content);
        Assert.Contains(all.Latents, l => l.Name.Contains("working response") && l.Scope == LatentScope.Observation && l.Treatment == LatentTreatment.Conditioned);
        int row = 0;
        for (int t = 0; t < 3; t++)
        {
            List<BassDraw> bare = [.. bySnapshot[t].Select(d => d with { Working = null })];
            PointwiseLogDensity alone = BassPointwise.LogDensity(Reading(t), table, [bare], content: content);
            Assert.DoesNotContain(alone.Latents, l => l.Name.Contains("working response"));
            for (int k = 0; k < alone.Draws; k++, row++)
                for (int u = 0; u < n; u++) Close(alone.Values[k, u], all.Values[row, u], 1e-9, $"snapshot {t}, draw {k}, unit {u}");
        }
        // A model that never moves reports no latent, whatever it holds.
        var fixedChain = new BassChain(model, [], Math.Log(0.2), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        for (int s = 0; s < 20; s++) fixedChain.Step();
        BassDraw kept = fixedChain.Compose(new Xoshiro256PlusPlus(8), table);
        Assert.Same(model.Working, kept.Working);
        Assert.DoesNotContain(BassPointwise.LogDensity(model, table, [[kept]], content: content).Latents, l => l.Name.Contains("working response"));
    }

    // A Gaussian chain's draws name no working response, so its fits and hashes hold what they held.
    [Fact]
    public void AGaussianDraw_NamesNoWorkingResponse()
    {
        BassModel model = BassModel.Create(BassSweepTests.Config(maxBasis: 3), BassSweepTests.Data());
        var chain = new BassChain(model, [], Math.Log(0.05), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        for (int s = 0; s < 20; s++) chain.Step();
        Assert.Null(chain.Compose(new Xoshiro256PlusPlus(8), new BassTermTable()).Working);
        Assert.True(model.Working.IsResponseItself);
    }

    private static (List<BassDraw> Draws, BassTermTable Table, double[] MeanOffset, bool Weighted) Composed(BassData data, BassConfig config, bool weighted, bool moving, int snapshots = 3)
    {
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        BassPriorReference reference = BassPriorReference.Moments(y);
        double[]? Weights(int t) => weighted ? BassWeightedCoreTests.Weights(n, seed: 30 + t) : null;
        BassModel model = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, Offset(n, 0), Weights(0), reference: reference)));
        var chain = new BassChain(model, [], Math.Log(0.2), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        var table = new BassTermTable();
        var rng = new Xoshiro256PlusPlus(7);
        var draws = new List<BassDraw>();
        var mean = new double[n];
        for (int t = 0; t < snapshots; t++)
        {
            if (t > 0 && moving) Assert.Null(chain.Follow(model.Working.Moved(Offset(n, t), Weights(t))));
            double[] offset = Offset(n, moving ? t : 0);
            for (int k = 0; k < 4; k++)
            {
                for (int s = 0; s < 10; s++) chain.Step();
                draws.Add(chain.Compose(rng, table));
                for (int i = 0; i < n; i++) mean[i] += offset[i] / (4.0 * snapshots);
            }
        }
        return (draws, table, mean, weighted);
    }

    // The residuals read y − offset − f̄, with the offset the draws' mean; the surface's mean is the draws'. Observed stays y, and under
    // weights no one reference normal fits the residuals.
    [Theory]
    [InlineData(false, false)]
    [InlineData(false, true)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void TheTrainingFit_ReadsWhatTheSurfaceFits(bool weighted, bool moving)
    {
        (BassData data, BassFit fit) = BassPredictTests.FitOf(BassPointwiseTests.Kind.Scalar);
        (List<BassDraw> draws, BassTermTable table, double[] mean, _) = Composed(data, fit.Config, weighted, moving);
        var store = new DrawStore<BassDraw>();
        foreach (BassDraw draw in draws) store.Add(draw);
        BassFit composed = fit with { Terms = table, Draws = [store], DrawWeights = null };
        BassTrainingFit t = composed.TrainingFit(data);
        BassPrediction prediction = composed.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false });
        double[,] surface = prediction.Mean();
        ReadOnlySpan<double> y = data.Response.Values;
        Assert.NotNull(t.Offset);
        for (int i = 0; i < t.Cases; i++)
        {
            Assert.Equal(y[i], t.Observed[i, 0]);
            Assert.Equal(surface[i, 0], t.Mean[i, 0]);
            Close(mean[i], t.Offset![i, 0], 1e-12, $"case {i}: the offset");
            Close(y[i] - mean[i] - surface[i, 0], t.Residuals[i, 0], 1e-12, $"case {i}: the residual");
        }
        Assert.Equal(weighted, t.NoiseSd is null);

        // Gaussian BASS's training fit holds no offset, and its residuals are y less the mean.
        BassTrainingFit plain = fit.TrainingFit(data);
        Assert.Null(plain.Offset);
        Assert.NotNull(plain.NoiseSd);
        for (int i = 0; i < plain.Cases; i++) Assert.Equal(plain.Observed[i, 0] - plain.Mean[i, 0], plain.Residuals[i, 0]);
    }

    // The residual check reads the same residuals as the training fit.
    [Fact]
    public void TheResidualCheck_ReadsTheWorkingResponse()
    {
        BassData data = BassKernelEnumeration.FunctionalData();
        double[] y = data.Response.Values.ToArray();
        int n = y.Length, cases = data.Response.Cases, locations = data.Response.Locations;
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2) with
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 5, Support = new SupportRule(SupportKind.Count, 1.0, 1.0) },
        };
        (List<BassDraw> draws, BassTermTable table, double[] offset, _) = Composed(data, config, false, true);
        BassModel model = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, reference: BassPriorReference.Moments(y))));
        ResidualIndex index = new ResidualIndex.AlongVariable(0);
        ResidualDependence read = BassResidualDependence.Of(model, table, [draws], index)!;
        var fit = new BassDrawFit(model, table);
        var mean = new double[n];
        var f = new double[n];
        foreach (BassDraw draw in draws)
        {
            fit.Evaluate(draw, f);
            for (int i = 0; i < n; i++) mean[i] += f[i] / draws.Count;
        }
        var residuals = new double[cases, locations];
        for (int i = 0; i < cases; i++)
            for (int r = 0; r < locations; r++) residuals[i, r] = y[i * locations + r] - offset[i * locations + r] - mean[i * locations + r];
        ResidualDependence hand = ResidualDependence.LagOne(residuals, ((FunctionalResponse)data.Response).OutputVariables[0].Values.Span);
        for (int i = 0; i < cases; i++)
            if (!double.IsNaN(hand.PerCase[i])) Close(hand.PerCase[i], read.PerCase[i], 1e-10, $"case {i}");

        // Bare draws read the model's own offset, and draws of another response read none of theirs.
        double[] fixedOffset = Offset(n, 3);
        BassModel holding = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, fixedOffset, reference: BassPriorReference.Moments(y))));
        List<BassDraw> bare = [.. draws.Select(d => d with { Working = null })];
        for (int i = 0; i < cases; i++)
            for (int r = 0; r < locations; r++) residuals[i, r] = y[i * locations + r] - fixedOffset[i * locations + r] - mean[i * locations + r];
        ResidualDependence onHolding = BassResidualDependence.Of(holding, table, [bare], index)!;
        hand = ResidualDependence.LagOne(residuals, ((FunctionalResponse)data.Response).OutputVariables[0].Values.Span);
        for (int i = 0; i < cases; i++)
            if (!double.IsNaN(hand.PerCase[i])) Close(hand.PerCase[i], onHolding.PerCase[i], 1e-10, $"case {i}, a model's own offset");
        BassData elsewhere = BassKernelEnumeration.FunctionalData(seed: 23);
        BassModel plain = BassModel.Create(config, elsewhere);
        Assert.NotNull(BassResidualDependence.Of(plain, table, [draws], index));
    }

    [Fact]
    public void TheConditionalMeans_ReadEachDrawsWorkingResponse()
    {
        BassData data = BassWeightedCoreTests.Scattered(40, 11);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        BassPriorReference reference = BassPriorReference.Moments(y);
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2);
        double[] Weights(int t) => BassWeightedCoreTests.Weights(n, seed: 30 + t);
        BassModel Reading(int t) => BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, Offset(n, t), Weights(t), reference: reference)));
        BassModel model = Reading(0);
        var chain = new BassChain(model, [BassWeightedCoreTests.Terms[0], BassWeightedCoreTests.Terms[1]], Math.Log(0.2), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        var table = new BassTermTable();
        Assert.Null(chain.Follow(model.Working.Moved(Offset(n, 2), Weights(2))));
        BassDraw draw = chain.Compose(new Xoshiro256PlusPlus(8), table);
        double[] moved = new BassConditionalMeans(model, table).Coefficients(draw);
        double[] rebuilt = new BassConditionalMeans(Reading(2), table).Coefficients(draw with { Working = null });
        for (int j = 0; j < moved.Length; j++) Close(rebuilt[j], moved[j], 1e-10, $"coefficient {j}");
    }

    // ---- Restart ----

    // A chain's snapshot taken after a move, its design's columns re-evaluated under the moved weights, continues bit for bit.
    [Theory]
    [InlineData("moments", false, true)]
    [InlineData("moments", true, true)]
    [InlineData("moments", false, false)]
    [InlineData("frozen", false, true)]
    [InlineData("invariant", true, true)]
    public void ASnapshotAfterAMove_ContinuesBitForBit(string route, bool ridge, bool reweights)
    {
        BassData data = BassWeightedCoreTests.Scattered(60, 11);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        (BassConfig config, Func<double[], BassPriorReference?> rule, WeightedPriorGram gram, _) = BassPriorReferenceTests.Route(route, ridge, n);
        BassPriorReference? reference = rule(y);
        double[] w = BassWeightedCoreTests.Weights(n, seed: 31);
        BassModel model = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, Offset(n, 0), w, priorGram: gram, reference: reference)));
        var chain = new BassChain(model, [], Math.Log(0.2), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        for (int s = 0; s < 60; s++) chain.Step();
        WorkingResponse working = model.Working.Moved(Offset(n, 1), reweights ? BassWeightedCoreTests.Weights(n, seed: 32) : default);
        Assert.Null(chain.Follow(working));
        for (int s = 0; s < 3; s++) chain.Step();
        BassChain restored = BassChain.Restore(model, chain.Snapshot());
        Assert.Same(working, restored.Structure.DesignAt(restored.Current).Working);
        for (int s = 0; s < 200; s++)
        {
            BassState a = chain.Step(), b = restored.Step();
            Assert.Equal(a.Terms, b.Terms);
            Assert.Equal(a.LogTau, b.LogTau);
            Assert.Equal(a.Fit.ResidualSquares, b.Fit.ResidualSquares);
            Assert.Same(working, b.Fit.Working);
        }
        Assert.Equal(chain.Compose(new Xoshiro256PlusPlus(8), new BassTermTable()).Coefficients, restored.Compose(new Xoshiro256PlusPlus(8), new BassTermTable()).Coefficients);
    }

    // A ladder's replicas each hold their own working response, and a restart continues each, swaps and all.
    [Fact]
    public void ALaddersSnapshot_ContinuesEachReplicasWorkingResponse()
    {
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        double[] prior = [0.34, 0.33, 0.33];
        (BassConfig config, Func<int, double[]> offsets, Func<int, double[]?> weights, WeightedPriorGram gram, BassPriorReference? reference) = JointReading("apart", y, false);
        var observation = new Latent(r => new WorkingResponse(r, offsets(0), weights(0), priorGram: gram, reference: reference));
        BassModel model = BassModel.Create(config, data, observation);
        WorkingResponse[] workings = [model.Working, .. Enumerable.Range(1, prior.Length - 1).Select(z => model.Working.Moved(offsets(z), weights(z)))];
        observation.Block = new LatentBlock(workings, prior, false);
        var spec = new TemperingSpec { Levels = 3, HotEnd = 0.2 };
        var ladder = new BassLadder(model, spec, [], Math.Log(0.1), 500);
        for (int s = 0; s < 300; s++) ladder.Step();
        BassLadder restored = BassLadder.Restore(model, spec, ladder.Snapshot());
        var held = new HashSet<int>();
        for (int s = 0; s < 300; s++)
        {
            ladder.Step();
            restored.Step();
            for (int r = 0; r < 3; r++)
            {
                BassState a = ladder.Tempering.Replicas[r].Current, b = restored.Tempering.Replicas[r].Current;
                Assert.Equal(a.Terms, b.Terms);
                Assert.Same(a.Fit.Working, b.Fit.Working);
                Assert.Equal(a.Fit.ResidualSquares, b.Fit.ResidualSquares);
                held.Add(Array.IndexOf(workings, a.Fit.Working));
            }
        }
        Assert.Equal(prior.Length, held.Count);
    }

    // A state no design holds is loaded on the working response its fit names, whichever the design was on: what a state from the
    // reference rung reaches, on a replica that has moved.
    [Fact]
    public void AStateNoDesignHolds_IsLoadedOnItsOwnWorkingResponse()
    {
        BassData data = BassSweepTests.Data(large: true);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2);
        BassModel model = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, Offset(n, 0), BassWeightedCoreTests.Weights(n, seed: 21), reference: BassPriorReference.Moments(y))));
        WorkingResponse here = model.Working.Moved(Offset(n, 1), BassWeightedCoreTests.Weights(n, seed: 31));
        WorkingResponse there = model.Working.Moved(Offset(n, 2), BassWeightedCoreTests.Weights(n, seed: 32));
        BassTerm[] terms = BassWeightedCoreTests.Terms[..3];
        var moved = new BassStructureMoves(model, terms, Math.Log(0.2));
        var other = new BassStructureMoves(model, terms, Math.Log(0.2));
        BassState stateHere = moved.ProposeFollow(moved.Start, here), stateThere = other.ProposeFollow(other.Start, there);
        Assert.Same(here, moved.DesignAt(stateHere).Working);
        BassDesign reference = other.DesignAt(stateThere);
        Assert.Same(there, reference.Working);
        Assert.False(moved.Holds(stateThere));
        BassDesign shown = moved.DesignAt(stateThere);
        Assert.Same(there, shown.Working);
        Same(shown, reference, "the state shown");
        // Back to the model's own, as a state from the reference rung would be.
        BassState start = new BassStructureMoves(model, terms, Math.Log(0.2)).Start;
        Assert.Same(model.Working, start.Fit.Working);
        Assert.Same(model.Working, moved.DesignAt(start).Working);
    }

    // A ladder with a reference rung and a latent block: each replica's design stays on the working response its state names, through
    // swaps and states loaded from the prior. Mechanics only: nothing draws the latent at b = 0, so this ladder is not exact yet.
    [Fact]
    public void AReferenceRung_KeepsEachDesignOnItsStatesWorkingResponse()
    {
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        double[] prior = [0.34, 0.33, 0.33];
        (BassConfig config, Func<int, double[]> offsets, Func<int, double[]?> weights, WeightedPriorGram gram, BassPriorReference? reference) = JointReading("apart", y, false);
        var observation = new Latent(r => new WorkingResponse(r, offsets(0), weights(0), priorGram: gram, reference: reference));
        BassModel model = BassModel.Create(config, data, observation);
        WorkingResponse[] workings = [model.Working, .. Enumerable.Range(1, prior.Length - 1).Select(z => model.Working.Moved(offsets(z), weights(z)))];
        observation.Block = new LatentBlock(workings, prior, false);
        var ladder = new BassLadder(model, new TemperingSpec { Levels = 4, HotEnd = 0.2, ReferenceRung = true }, [], Math.Log(0.1), 500);
        var seen = new HashSet<WorkingResponse>();
        for (int s = 0; s < 400; s++)
        {
            ladder.Step();
            for (int r = 0; r < ladder.Tempering.Replicas.Count; r++)
            {
                BassState state = ladder.Tempering.Replicas[r].Current;
                Assert.Same(state.Fit.Working, ladder.Structures[r].DesignAt(state).Working);
                seen.Add(state.Fit.Working!);
            }
        }
        Assert.Equal(workings.Length, seen.Count);
    }

    // ---- A ladder's composition ----

    // σ² from its tempered conditional at the replica's rung, the replica handed its target given that σ².
    private sealed class TemperedSigma2Block(Func<BassTarget> target) : IBlock<ReversibleJumpChain<BassState>>
    {
        public string Name => "sigma2";

        public bool Update(ReversibleJumpChain<BassState> chain, Xoshiro256PlusPlus rng)
        {
            BassTarget t = target();
            double sigma2 = t.Marginal(chain.Current.Fit, t.HeldTau!.Value, chain.Temperature.InverseTemperature).Sigma2Posterior(t.Sigma2).Sample(rng);
            chain.ReplaceTarget(t.Given(sigma2));
            return true;
        }
    }

    // σ² drawn apart at every rung travels with the state through a swap, so the cold replica's visits are the collapsed posterior's, and
    // its composition reads the σ² it holds.
    [Fact]
    public void ALadder_ComposesUnderTheTargetItHolds()
    {
        BassData data = BassSweepTests.Data();
        BassConfig config = BassWeightedCoreTests.EnumerationConfig(Sigma2Prior.InverseGamma(2.0, 0.3));
        BassModel model = null!;
        var observation = new Latent(r => new WorkingResponse(r), BassSweepSlot.Scale) { Block = new TemperedSigma2Block(() => model.Target) };
        model = BassModel.Create(config, data, observation);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(BassModel.Create(config, data));
        double[] pi = space.Posterior();
        const int chains = 3, burn = 2_000, sweeps = 20_000;
        var visits = new double[pi.Length];
        long swaps = 0;
        var table = new BassTermTable();
        for (int c = 0; c < chains; c++)
        {
            var ladder = new BassLadder(model, new TemperingSpec { Levels = 3, HotEnd = 0.2 }, [], Math.Log(0.1), 700 + c);
            for (int s = 0; s < burn; s++) ladder.Step();
            for (int s = 0; s < sweeps; s++)
            {
                ladder.Step();
                visits[space.Find(ladder.Current.Terms)] += 1.0;
                if (s % 500 != 0) continue;
                var held = (BassTarget)ladder.Moves.Target;
                Assert.Equal(held.GivenSigma2, ladder.Compose(new Xoshiro256PlusPlus(s), table).Sigma2);
                foreach (ReversibleJumpChain<BassState> replica in ladder.Tempering.Replicas)
                    Assert.Same(replica.Target, ((BassEvaluation)replica.CurrentEvaluation).Target);
            }
            swaps += ladder.Tempering.PairStats.Sum(p => p.Accepts);
        }
        double total = (double)chains * sweeps;
        foreach (int i in Enumerable.Range(0, pi.Length).OrderByDescending(i => pi[i]).Take(10))
            Assert.True(Math.Abs(pi[i] - visits[i] / total) <= 0.02, $"state {i}: π {pi[i]:F4}, visits {visits[i] / total:F4}");
        Assert.True(swaps > 100, $"The ladders swapped {swaps} times.");
        _out.WriteLine($"largest π {pi.Max():F4}, its visits {visits[Array.IndexOf(pi, pi.Max())] / total:F4}; {swaps} swaps");
    }

    // ---- How often F's gate reads a latent draw ----

    // A reading, not a check: how often a draw of the latent scales, conjugate given the residuals, puts a posterior structure outside F.
    // t's scales are Gamma, a Laplace's inverse Gaussian, whose weights grow as a residual shrinks. The gate should bite most where few
    // rows carry each column, so a small N with the support floor low sits beside the two fixtures. Run it behind BASS_FOLLOW_GATE = 1.
    [Theory]
    [InlineData(BassMixing.Friedman)]
    [InlineData(BassMixing.Functional)]
    [InlineData("small N")]
    public void HowOften_ALatentDrawLeavesF(string name)
    {
        if (Environment.GetEnvironmentVariable("BASS_FOLLOW_GATE") != "1") return;
        BassData data;
        BassConfig config;
        if (name == "small N")
        {
            data = BassWeightedCoreTests.Scattered(30, 11);
            config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2, maxBasis: 8) with
            {
                Coefficients = new BassCoefficientPrior(new GForm(), TauPrior.UnitInformationGamma()),
            };
        }
        else
        {
            BassMixing.Fixture f = BassMixing.Get(name, null);
            data = f.Data;
            config = f.Config;
        }
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        var reference = BassPriorReference.Moments(y);
        var observation = new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: Unit(n), reference: reference));
        BassModel model = BassModel.Create(config, data, observation);
        var chain = new BassChain(model, [], Math.Log(1.0 / n), new Xoshiro256PlusPlus(11), new Xoshiro256PlusPlus(12));
        for (int s = 0; s < 3000; s++) chain.Step();
        var rng = new Xoshiro256PlusPlus(13);
        (string Label, Func<double, double, double> Draw)[] laws =
        [
            ("t(3)", (r2, s2) => Maths.Distributions.Gamma.ShapeRate(2.0, 0.5 * (3.0 + r2 / s2)).Sample(rng)),
            ("t(1)", (r2, s2) => Maths.Distributions.Gamma.ShapeRate(1.0, 0.5 * (1.0 + r2 / s2)).Sample(rng)),
            ("Laplace", (r2, s2) => InverseGaussian(Math.Sqrt(2.0 * s2 / Math.Max(r2, 1e-300)), 2.0, rng)),
        ];
        var counts = new long[laws.Length];
        var refused = new long[laws.Length];
        var angles = new List<double>[laws.Length];
        for (int l = 0; l < laws.Length; l++) angles[l] = [];
        var rows = new double[n];
        const int states = 60, drawsPerState = 20;
        for (int k = 0; k < states; k++)
        {
            for (int s = 0; s < 100; s++) chain.Step();
            BassDesign design = chain.Structure.DesignAt(chain.Current);
            double sigma2 = design.ResidualSquares / Math.Max(n - design.Count - 1, 1);
            ReadOnlySpan<double> r = design.Residual;
            for (int l = 0; l < laws.Length; l++)
                for (int d = 0; d < drawsPerState; d++)
                {
                    double max = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        rows[i] = Math.Max(laws[l].Draw(r[i] * r[i], sigma2), 1e-300);
                        max = Math.Max(max, rows[i]);
                    }
                    if (!double.IsFinite(max)) continue;
                    counts[l]++;
                    DesignProposal proposal = design.ProposeFollow(model.Working.Moved(default, rows));
                    if (proposal.Refusal is not null) refused[l]++;
                    else if (design.Count > 0) angles[l].Add(proposal.MinAngle);
                }
        }
        for (int l = 0; l < laws.Length; l++)
        {
            double[] sorted = [.. angles[l].Order()];
            double At(double q) => sorted.Length == 0 ? double.NaN : sorted[(int)(q * (sorted.Length - 1))];
            _out.WriteLine(string.Create(CultureInfo.InvariantCulture,
                $"GATE\t{name}\tN={n}\t{laws[l].Label}\tM={chain.Current.Count}\t{refused[l]}/{counts[l]}\tsmallest angle: min {At(0):E2}, 1% {At(0.01):E2}, median {At(0.5):E2}"));
        }
    }

    // A draw of the inverse Gaussian (Michael, Schucany and Haas 1976).
    private static double InverseGaussian(double mean, double shape, Xoshiro256PlusPlus rng)
    {
        double normal = BassQualification.Gaussian(rng), v = normal * normal;
        double x = mean + mean * mean * v / (2.0 * shape) - mean / (2.0 * shape) * Math.Sqrt(4.0 * mean * shape * v + mean * mean * v * v);
        return rng.NextDouble() <= mean / (mean + x) ? x : mean * mean / x;
    }
}
