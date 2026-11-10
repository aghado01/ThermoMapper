using System;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

public sealed partial class MultiOutputLawTests
{
    [Fact]
    public void ProposalFitsAgainstDenseRefits()
    {
        var rng = new Xoshiro256PlusPlus(20261123);
        var y = new double[GewekeN, 2];
        for (int i = 0; i < GewekeN; i++)
        {
            y[i, 0] = 400 + 300 * Math.Max((i * 7 % GewekeN) / 11.0 - 0.4, 0) + Math.Sin(i);
            y[i, 1] = -200 + 100 * Math.Max(i / 11.0 - 0.2, 0) + Math.Cos(i);
        }
        var model = BassModel.Create(GewekeConfig, GewekeData(y));
        var prior = new BassPriorSampler(model);
        var iw = new ReferenceSigmaPrior(6, GewekePsi);
        var refit = model.NewDesign();
        double[] nodes = [-25, -15, -7, 0];
        int compared = 0, refused = 0;
        double worst = 0, minAngle = 1;
        for (int i = 0; i < 16_000; i++)
        {
            var start = prior.DrawState(rng).State;
            var moves = new BassStructureMoves(model, start.Terms, start.LogTau);
            foreach (var move in moves.Moves)
            {
                if (move.Weight(moves.Start) == 0 || move.Propose(moves.Start, rng) is not { } p) continue;
                var s = p.Candidate;
                var refusal = refit.Load(s.Terms);
                Assert.Equal(refusal is not null, model.Target.Evaluate(s).OffSupport);
                if (refusal is not null) { refused++; continue; }
                var xc = MultiOutputCoreTests.Columns(model, refit);
                minAngle = Math.Min(minAngle, refit.MinAngle);
                foreach (double u in nodes)
                {
                    double actual = model.Target.Evaluate(s.WithLogTau(u)).LogMarginal;
                    double expected = MatrixTReference.LogMarginal(MatrixTReference.Collapse(
                        MultiOutputCase.GForm(y, xc, GPriorCentre.Zero, Math.Exp(u), 1.0)), iw);
                    worst = Math.Max(worst, Math.Abs(actual - expected));
                }
                compared++;
            }
        }
        output.WriteLine($"{compared} candidate fits x {nodes.Length} scales, {refused} matching refusals, min angle {minAngle:G5}, max absolute marginal error {worst:G5} nats");
        Assert.True(compared > 10_000 && worst < 1e-6);
    }
}
