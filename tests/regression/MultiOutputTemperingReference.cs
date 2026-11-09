using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;

namespace Maths.Regression.Tests;

// Dense matrix-normal integration, then inverse-Wishart integration, then quadrature over ln τ.
// The structural weights come from the enumerated prior on F; normalizing them is essential for evidence.
internal static class MultiOutputTemperingReference
{
    internal const int Cases = 8;
    internal static readonly double[,] Psi = { { 1.2, -0.3 }, { -0.3, 0.7 } };

    internal static BassData Data()
    {
        double[] x = [.. Enumerable.Range(0, Cases).Select(i => i / (Cases - 1.0))];
        var y = new double[Cases, 2];
        for (int i = 0; i < Cases; i++)
        {
            y[i, 0] = 0.4 + Math.Sin(2.5 * x[i]) + 0.12 * Math.Cos(7 * i);
            y[i, 1] = -0.3 + 0.8 * Math.Abs(x[i] - 0.4) + 0.15 * Math.Sin(5 * i);
        }
        return new BassData([new ContinuousPredictor("x", x)], new MultiOutputResponse(y, ["a", "b"]));
    }

    internal static BassConfig Config(GPriorCentre centre = GPriorCentre.Zero, double nu = 6.0) => new()
    {
        Structure = new BassStructureSpec { MaxBasis = 1, MaxInteraction = 1 },
        Intercept = centre,
        OutputCovariance = OutputCovariancePrior.InverseWishart(nu, Psi),
    };

    internal static BassTemperingTests.RungLaw[] Exact(BassModel model, double nu, double[] rungs, double panelWidth = 0.5)
    {
        var states = BassTermSpace.F(model);
        var response = (MultiOutputResponse)model.Data.Response;
        var y = new double[response.Cases, response.Outputs];
        for (int i = 0; i < response.Cases; i++)
            for (int l = 0; l < response.Outputs; l++) y[i, l] = response.Values[i * response.Outputs + l];
        var prior = new ReferenceSigmaPrior(nu, Psi);
        var nodes = BassTargetTests.Nodes(-70, 8, panelWidth).ToArray();
        var marginal = rungs.Select(_ => new double[states.Count]).ToArray();
        var means = rungs.Select(_ => new double[states.Count]).ToArray();
        var tauMeans = rungs.Select(_ => new double[states.Count]).ToArray();
        var design = model.NewDesign();
        double logRate = Math.Log(Cases / 2.0);
        var values = new double[nodes.Length];
        for (int s = 0; s < states.Count; s++)
        {
            if (design.Load(states[s].Terms) is not null) throw new InvalidOperationException("An enumerated structure was refused.");
            double[,] xc = MultiOutputCoreTests.Columns(model, design);
            for (int r = 0; r < rungs.Length; r++)
            {
                for (int j = 0; j < nodes.Length; j++)
                {
                    var (u, weight) = nodes[j];
                    var c = MultiOutputCase.GForm(y, xc, model.Target.Intercept, Math.Exp(u), rungs[r]);
                    double likelihood = MatrixTReference.LogMarginal(MatrixTReference.Collapse(c), prior);
                    double logPrior = 0.5 * (logRate + u) - SpecialFunctions.LogGamma(0.5) - Cases / 2.0 * Math.Exp(u);
                    values[j] = Math.Log(weight) + logPrior + likelihood;
                }
                double z = BassTemperingTests.LogSumExp(values);
                marginal[r][s] = states[s].LogPrior + z;
                means[r][s] = nodes.Select((p, j) => p.U * Math.Exp(values[j] - z)).Sum();
                tauMeans[r][s] = nodes.Select((p, j) => Math.Exp(p.U + values[j] - z)).Sum();
            }
        }
        double priorMass = BassTemperingTests.LogSumExp(states.Select(s => s.LogPrior));
        return Enumerable.Range(0, rungs.Length).Select(r =>
        {
            double z = BassTemperingTests.LogSumExp(marginal[r]);
            var count = new double[model.MaxBasis + 1];
            double mean = 0, tau = 0;
            for (int s = 0; s < states.Count; s++)
            {
                double weight = Math.Exp(marginal[r][s] - z);
                count[states[s].Terms.Length] += weight;
                mean += weight * means[r][s];
                tau += weight * tauMeans[r][s];
            }
            return new BassTemperingTests.RungLaw(z - priorMass, count, tau, mean);
        }).ToArray();
    }
}
