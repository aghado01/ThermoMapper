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

    internal static OutputPencil Of(PriorPencil pencil, double[][] crosses, ReadOnlySpan<double> means, GPriorCentre centre, int rows)
    {
        int k = crosses.Length, offset = pencil.InterceptInside ? 1 : 0;
        var outputs = new PriorPencil[k];
        var levels = new double[k];
        for (int l = 0; l < k; l++)
        {
            var cross = new double[pencil.Count + offset];
            if (offset > 0)
            {
                levels[l] = centre == GPriorCentre.Zero ? 0 : means[l];
                cross[0] = rows * (means[l] - levels[l]);
            }
            crosses[l].CopyTo(cross, offset);
            outputs[l] = pencil.WithCross(cross);
        }
        return new OutputPencil(outputs, levels);
    }

    internal OutputMarginal At(double[,] residual, int rows, double tau, double b)
        => OutputMarginal.FromPencil(residual, _outputs, rows, tau, b);

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
