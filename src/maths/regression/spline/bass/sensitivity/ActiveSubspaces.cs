using System;
using System.Collections.Generic;
using Maths.Regression.Subspaces;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>How C per draw is computed (E8, decision 1).</summary>
public enum ActiveSubspaceRoute
{
    /// <summary>The exact route where it holds (a product law over terms that have factors), the sampled route otherwise.</summary>
    Auto,

    /// <summary>The moment service's derivative moments (<see cref="BassMoments.GradientOuterProduct"/>): exact.</summary>
    Exact,

    /// <summary>The gradient readout averaged over points drawn from the law, shared by every draw, with its Monte Carlo error.</summary>
    Sampled,
}

/// <summary>
/// The input metric C is read in (E8, decision 2). For metric coordinates z with box coordinates x = Lz + c, C_z = LᵀC_xL, so the
/// subspace moves under any non-orthogonal change of units; every readout states its metric.
/// </summary>
public enum InputMetric
{
    /// <summary>Each continuous input scaled to its box, the declared domain or else the training range, as the fit's transform does.</summary>
    Box,

    /// <summary>The inputs' own units.</summary>
    Raw,

    /// <summary>Whitened by the law's covariance, affine-equivariant: under a product law, each input scaled by its standard deviation,
    /// so C_z is the generalized problem C v = λ Σ⁻¹ v in z's coordinates. Under a product uniform law it is the box scaling up to a
    /// constant.</summary>
    Whitened,
}

/// <summary>
/// What E8's readout reads (<see cref="BassFit.ActiveSubspace"/>): the law, the metric, the subspace's dimension k, the route and the
/// draws. A declarative record, JSON-serializable.
/// </summary>
public sealed record ActiveSubspaceSpec
{
    /// <summary>The predictors' measures by name, in their own units, as <see cref="SensitivitySpec.Measures"/> takes them; a predictor
    /// left out takes the uniform law on its box, or uniform level probabilities.</summary>
    public IReadOnlyDictionary<string, InputMeasure>? Measures { get; init; }

    /// <summary>The metric C and the subspaces are read in.</summary>
    public InputMetric Metric { get; init; } = InputMetric.Box;

    /// <summary>k, the subspaces' dimension, from 1 to the number of continuous predictors.</summary>
    public int Dimension { get; init; } = 1;

    /// <summary>The route.</summary>
    public ActiveSubspaceRoute Route { get; init; } = ActiveSubspaceRoute.Auto;

    /// <summary>The sampled route's points N, drawn once from the law and shared by every draw.</summary>
    public int Points { get; init; } = 4096;

    /// <summary>The seed of the sampled route's points.</summary>
    public int Seed { get; init; } = 1;

    /// <summary>The draws, as indices into the fit's draws pooled in chain order, repeats allowed; every draw when null.</summary>
    public IReadOnlyList<int>? Draws { get; init; }
}

/// <summary>
/// Decision 6's per-draw consistency checks, reported and not enforced, for each draw and continuous predictor in the readout's
/// order. Poincaré on a uniform input on [a, b]: total(v) ≤ ((b − a)²/π²)·C_vv in box units (NaN where the input's law is not
/// uniform). On a Legendre fit under the product uniform law on the boxes, SG-4(b)'s total(v) ≤ ½·E[x_v(1 − x_v)(∂_v f)²] as well.
/// A ratio above one beyond rounding, or beyond three Monte Carlo errors on the sampled route, is a violation, which is a defect.
/// </summary>
/// <param name="Total">total(v), the input's unnormalized total variance, per draw (the moment route's).</param>
/// <param name="Bound">The Poincaré bound per draw.</param>
/// <param name="Ratio">Total over Bound.</param>
/// <param name="JacobiRatio">Total over ½·the Jacobi measure, for a Legendre fit under the product uniform law; null otherwise.</param>
/// <param name="Violations">The ratios above one beyond their tolerance.</param>
public sealed record PoincareCheck(double[,] Total, double[,] Bound, double[,] Ratio, double[,]? JacobiRatio, int Violations);

/// <summary>
/// E8's readout (<see cref="BassFit.ActiveSubspace"/>): the active subspaces of a fit's draws under a law, in a declared metric.
/// </summary>
/// <param name="Inputs">The continuous predictors, in C's order.</param>
/// <param name="Excluded">The categorical predictors, which take no part in C (decision 3): their effect is Sobol's.</param>
/// <param name="Route">The route that ran.</param>
/// <param name="Metric">The metric C and the subspaces are in.</param>
/// <param name="Map">L, p × p: box coordinates are L times metric coordinates (plus a shift), so C here is LᵀC_boxL.</param>
/// <param name="Measure">The law, in the fit's transformed units.</param>
/// <param name="Points">The sampled route's points N; 0 on the exact route.</param>
/// <param name="Posterior">The draws' C, subspaces and summaries, in the metric.</param>
/// <param name="MonteCarloError">The sampled route's standard error of each entry of each draw's C, in the metric; null on the exact route.</param>
/// <param name="Poincare">The per-draw consistency checks.</param>
/// <param name="DrawIndex">Each draw's index among the fit's pooled draws.</param>
public sealed record ActiveSubspaceAnalysis(
    IReadOnlyList<string> Inputs, IReadOnlyList<string> Excluded, ActiveSubspaceRoute Route, InputMetric Metric, double[,] Map,
    ProductMeasure Measure, int Points, SubspacePosterior Posterior, double[][,]? MonteCarloError, PoincareCheck Poincare,
    IReadOnlyList<int> DrawIndex);
