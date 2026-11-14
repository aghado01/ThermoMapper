namespace Clustering.Geometric.EMeans;

/// <summary>How the initial centroids are chosen when the caller supplies none.</summary>
public enum EMeansInit
{
    /// <summary><c>k</c> distinct training points drawn uniformly (the paper's default).</summary>
    RandomPoints,

    /// <summary>k-means++ seeding (<see cref="Clustering.Geometric.KMeans.KMeansPlusPlusSeeding"/>): one pass over
    /// the data per centroid before training starts.</summary>
    KMeansPlusPlus,
}

/// <summary>How the batch size <c>B</c> is resolved.</summary>
public enum BatchStrategy
{
    /// <summary>The batch law of <see cref="EMeansScalingLaw"/>; the manifest records the value and whether the
    /// problem lies outside the range the law was fitted on.</summary>
    Auto,

    /// <summary><see cref="BatchConfig.Size"/> as given.</summary>
    Explicit,
}

/// <summary>Batch size of one EMA update.</summary>
public sealed record BatchConfig
{
    public BatchStrategy Strategy { get; init; } = BatchStrategy.Auto;

    /// <summary>Points per batch under <see cref="BatchStrategy.Explicit"/>, between 1 and <c>n</c>.</summary>
    public int? Size { get; init; }
}

/// <summary>How the peak step <c>α₀ = 1 − γ₀</c> is resolved.</summary>
public enum StepStrategy
{
    /// <summary>The step law of <see cref="EMeansScalingLaw"/> at the resolved batch size.</summary>
    Auto,

    /// <summary><see cref="StepConfig.Peak"/> as given.</summary>
    Explicit,
}

/// <summary>Shape of the step across the run's updates.</summary>
public enum StepSchedule
{
    /// <summary>Cosine anneal from the peak on the first update towards <see cref="StepConfig.FloorFraction"/> of
    /// it on the last.</summary>
    Cosine,

    /// <summary>The peak step on every update.</summary>
    Constant,
}

/// <summary>
/// The EMA step <c>α = 1 − γ</c>: the fraction of the way the stored statistics move towards a batch's statistics.
/// </summary>
public sealed record StepConfig
{
    public StepStrategy Strategy { get; init; } = StepStrategy.Auto;

    /// <summary>Peak step <c>α₀</c> under <see cref="StepStrategy.Explicit"/>, in <c>(0, 1]</c>.</summary>
    public double? Peak { get; init; }

    public StepSchedule Schedule { get; init; } = StepSchedule.Cosine;

    /// <summary>The cosine schedule's floor <c>ε</c> as a fraction of the peak, in <c>[0, 1]</c>.</summary>
    public double FloorFraction { get; init; } = 0.1;
}

/// <summary>
/// Dead-cluster revival: a cluster whose count has fallen below <see cref="Threshold"/> times the mean count takes
/// half of a live donor, the donor's centroid splitting into two children either side of it.
/// </summary>
public sealed record RevivalConfig
{
    public bool Enabled { get; init; } = true;

    /// <summary>Updates between revival checks.</summary>
    public int Period { get; init; } = 5;

    /// <summary>Fraction <c>τ</c> of the mean count below which a cluster counts as dead, in <c>[0, 1)</c>.</summary>
    public double Threshold { get; init; } = 0.1;

    /// <summary>Half-distance <c>δ</c> between the two children of a split, as a fraction of the donor's
    /// root-mean-square radius.</summary>
    public double Perturbation { get; init; } = 1.0 / 1024.0;
}

/// <summary>Which centroid set a run returns as <see cref="EMeansResult.Centroids"/>.</summary>
public enum EMeansFinalization
{
    /// <summary>Each cluster's mean over the points assigned to it during the last epoch (the paper's default,
    /// aimed at the training error).</summary>
    LastEpochMean,

    /// <summary>The training EMA's own centroids <c>S/N</c>.</summary>
    Ema,

    /// <summary>The slow EMA's centroids (aimed at the error on points outside the training set).</summary>
    SlowEma,
}

/// <summary>
/// All settings of an e-means run (Bessa et al. 2026, arXiv:2611.09501): mini-batch k-means whose centroids are the
/// ratio of exponentially decayed sums and counts. Pure declarative data; initial centroids supplied by the caller
/// ride the call to <see cref="EMeans.Cluster"/>.
/// </summary>
public sealed record EMeansConfig
{
    /// <summary>Number of clusters <c>k</c>, between 1 and <c>n</c>.</summary>
    public required int K { get; init; }

    /// <summary>Epoch budget <c>T</c>, one epoch being a pass over all points in a fresh random order. The step
    /// schedule spans the whole budget, so this sets the method's course as well as its cost.</summary>
    public required int Epochs { get; init; }

    /// <summary>Master seed; the initialization, shuffle and revival streams are derived from it by
    /// <see cref="Maths.Rng.SeedTree"/>.</summary>
    public int Seed { get; init; } = 42;

    public EMeansInit Init { get; init; } = EMeansInit.RandomPoints;

    /// <summary>Count each cluster starts with, its initial centroid entering as a virtual point of this mass.</summary>
    public double InitialCount { get; init; } = 1e-10;

    public BatchConfig Batch { get; init; } = new();

    public StepConfig Step { get; init; } = new();

    public RevivalConfig Revival { get; init; } = new();

    public EMeansFinalization Finalization { get; init; } = EMeansFinalization.LastEpochMean;

    /// <summary>The slow EMA's constant step as a multiple of the schedule's floor step <c>ε·α₀</c>.</summary>
    public double SlowStepMultiplier { get; init; } = 1.0;
}
