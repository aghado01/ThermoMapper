using System;
using Maths.LinAlg;

namespace Maths.Regression.GPrior;

/// <summary>A shared prior pencil with K response cross products, for a matrix-normal coefficient block.</summary>
public sealed class OutputPencil
{
    private readonly PriorPencil[] _outputs;
    private readonly double[] _levels;

    private OutputPencil(PriorPencil[] outputs, double[] levels)
    {
        _outputs = outputs;
        _levels = levels;
    }

    // `mass` is the intercept's likelihood Gram entry: N at unit weights, Σwᵢ under the reading's weights (E34's slice 4).
    // `levels`, when given, are the intercept's levels frozen apart from the reading (a noise's, E34's slice 4); else the centre's rule.
    internal static OutputPencil Of(PriorPencil pencil, double[][] crosses, ReadOnlySpan<double> means, GPriorCentre centre, double mass, double[]? levels = null)
    {
        int k = crosses.Length, offset = pencil.InterceptInside ? 1 : 0;
        var outputs = new PriorPencil[k];
        var levelsOf = new double[k];
        for (int l = 0; l < k; l++)
        {
            var cross = new double[pencil.Count + offset];
            if (offset > 0)
            {
                levelsOf[l] = centre == GPriorCentre.Zero ? 0 : levels?[l] ?? means[l];
                cross[0] = mass * (means[l] - levelsOf[l]);
            }
            crosses[l].CopyTo(cross, offset);
            outputs[l] = pencil.WithCross(cross);
        }
        return new OutputPencil(outputs, levelsOf);
    }

    internal OutputMarginal At(double[,] residual, int rows, double tau, double b, GPriorWeights? weights = null)
        => OutputMarginal.FromPencil(residual, _outputs, rows, tau, b, weights);

    internal (ISpdFactor Factor, double[,] Mean) Conditional(double tau, double b)
    {
        int p = _outputs[0].Eigenvalues.Length, k = _outputs.Length;
        var mean = new double[p, k];
        var column = new double[p];
        for (int l = 0; l < k; l++)
        {
            _outputs[l].PosteriorMean(tau, b, column);
            for (int j = 0; j < p; j++) mean[j, l] = column[j];
            if (_outputs[l].InterceptInside) mean[0, l] += _levels[l];
        }
        return (_outputs[0].Factor(tau, b), mean);
    }
}
