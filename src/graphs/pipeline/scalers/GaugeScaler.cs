using System;
using Graphs.Primitives;

namespace Graphs.Pipeline.Scalers;

/// <summary>
/// Stage 5, last: multiplies the finished couplings by one global factor. With
/// <see cref="GaugeKind.InverseMeanDegree"/> the factor is 1/K-hat, K-hat the mean degree of the
/// graph the inner scaler returned, so it follows every threshold and the floor never depends on
/// it; the result is the inner graph with temperature rescaled, and nothing else.
/// </summary>
public sealed class GaugeScaler : IEdgeScaler
{
    private readonly IEdgeScaler _inner;
    private readonly GaugeKind _gauge;

    public GaugeScaler(IEdgeScaler inner, GaugeKind gauge)
    {
        _inner = inner ?? throw new ArgumentNullException(nameof(inner));
        _gauge = gauge;
    }

    public ScalerResult Scale(NeighborSelection refined, int n)
    {
        ScalerResult inner = _inner.Scale(refined, n);
        if (_gauge == GaugeKind.None) return inner;

        CsrGraph graph = inner.Graph;
        if (graph.NodeCount == 0 || graph.Targets.Length == 0) return inner;

        double factor = graph.NodeCount / (double)graph.Targets.Length;
        var weights = new double[graph.Weights.Length];
        for (int i = 0; i < weights.Length; i++) weights[i] = graph.Weights[i] * factor;

        // In log space the factor is a shift, ln g, and the gauge stays a pure rescaling of temperature.
        double[]? logWeights = null;
        if (graph.LogWeights is not null)
        {
            double shift = Math.Log(factor);
            logWeights = new double[graph.LogWeights.Length];
            for (int i = 0; i < logWeights.Length; i++) logWeights[i] = graph.LogWeights[i] + shift;
        }

        return inner with
        {
            Graph = new CsrGraph
            {
                Targets = graph.Targets,
                RowPointers = graph.RowPointers,
                NodeCount = graph.NodeCount,
                Weights = weights,
                LogWeights = logWeights,
            },
            GaugeFactor = factor,
        };
    }
}
