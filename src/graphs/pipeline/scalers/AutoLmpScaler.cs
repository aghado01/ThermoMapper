using System;
using System.Collections.Generic;
using Graphs;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Primitives;
using Graphs.Proximity;

namespace Graphs.Pipeline.Scalers;

/// <summary>
/// Stage 5 — automatic Local Mutual Proximity wrapping.
/// Evaluates the inner scaler first, inspects the tentative graph for
/// near-zero edge collapse, and applies LMP only when the graph appears
/// numerically fragile.
/// </summary>
internal sealed class AutoLmpScaler : IEdgeScaler
{
    private readonly IEdgeScaler _innerScaler;
    private readonly DiagnosticsLog _log;
    private readonly double _nearZeroEdgeRatioThreshold;
    private readonly Func<NeighborSelection, int, IReadOnlySet<(int Lo, int Hi)>?>? _protectedProvider;

    public AutoLmpScaler(
        IEdgeScaler innerScaler,
        DiagnosticsLog log,
        double nearZeroEdgeRatioThreshold = 0.05,
        Func<NeighborSelection, int, IReadOnlySet<(int Lo, int Hi)>?>? protectedEdgeProvider = null)
    {
        _innerScaler = innerScaler ?? throw new ArgumentNullException(nameof(innerScaler));
        _log = log ?? throw new ArgumentNullException(nameof(log));
        if (nearZeroEdgeRatioThreshold < 0.0 || nearZeroEdgeRatioThreshold > 1.0)
            throw new ArgumentOutOfRangeException(nameof(nearZeroEdgeRatioThreshold));
        _nearZeroEdgeRatioThreshold = nearZeroEdgeRatioThreshold;
        _protectedProvider = protectedEdgeProvider;
    }

    public ScalerResult Scale(NeighborSelection refined, int n)
    {
        ScalerResult tentative = _innerScaler.Scale(refined, n);
        double ratio = ComputeNearZeroRatio(tentative.Graph);

        if (ratio > _nearZeroEdgeRatioThreshold)
        {
            _log.Warning(
                "Scaling",
                $"Auto-LMP enabled because tentative global scaling produced a near-zero edge ratio " +
                $"{ratio:P1}, which exceeds the threshold {_nearZeroEdgeRatioThreshold:P1}.");

            IReadOnlySet<(int Lo, int Hi)>? protectedEdges = _protectedProvider?.Invoke(refined, n);
            CsrGraph rescaled = LocalMutualProximity.ApplyLocalScaling(
                tentative.Graph,
                weightsAreCouplings: true,
                protectedEdges: protectedEdges);
            return tentative with { Graph = rescaled };
        }

        if (ratio > _nearZeroEdgeRatioThreshold * 0.5)
        {
            _log.Warning(
                "Scaling",
                $"Tentative global scaling produced {ratio:P1} near-zero edges; graph may still be numerically fragile.");
        }
        else
        {
            _log.Info(
                "Scaling",
                $"Tentative near-zero edge ratio is {ratio:F4}; LMP is not required.");
        }

        return tentative;
    }

    private static double ComputeNearZeroRatio(CsrGraph graph)
    {
        const double NearZero = 1e-8;
        int count = 0;
        double[] weights = graph.Weights;
        for (int i = 0; i < weights.Length; i++)
        {
            if (weights[i] < NearZero)
            {
                count++;
            }
        }

        return weights.Length == 0 ? 0.0 : (double)count / weights.Length;
    }
}
