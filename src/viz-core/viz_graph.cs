using System;

namespace Viz;

public enum VizKernel
{
    Gaussian,
    Cauchy,
    Laplacian,
    Linear,
}

public abstract record ProximitySpec(NeighborRule Kind);

public sealed record KnnSpec(int K) : ProximitySpec(NeighborRule.Knn);
public sealed record MutualKnnSpec(int K) : ProximitySpec(NeighborRule.MutualKnn);
public sealed record EpsilonBallSpec(double Epsilon) : ProximitySpec(NeighborRule.EpsilonBall);
public sealed record MstAugmentedSpec(int K) : ProximitySpec(NeighborRule.MstAugmented);

public enum NeighborRule
{
    Knn,
    MutualKnn,
    EpsilonBall,
    MstAugmented,
}
