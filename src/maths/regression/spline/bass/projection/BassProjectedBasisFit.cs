using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// Rotate-and-refit for a basis fit (E8, slice 8c; decisions E1 to E3, the recommended options): the components of a
/// <see cref="BassBasisFit"/> refitted on the active variables of one <see cref="BassProjection"/>, the frame read from the basis fit's
/// own E[C] over its reconstruction under the output metric (E1), so every component keeps one input space and the reconstruction's
/// readouts stay whole.
/// <list type="bullet">
/// <item><b>The decomposition is reused</b> (E2): for the refit on the full data, K depends on the outputs alone and the refit changes
/// only the inputs, so the native fit's decomposition is the refit's, the same object.</item>
/// <item><b>The component mode is carried over</b>: components fitted apart are refitted apart, a joint fit jointly (E34's), each under
/// its native config with the basis fit's master seed; a hinge config that leaves the normalization unset takes the training range
/// (slice 7d).</item>
/// <item><b>Prediction</b> takes cases in the original units and projects them; <b>the gradient</b> of the reconstruction
/// y_ℓ = o_ℓ + s_ℓ·Σ_k K_ℓk w_k is ∂y_ℓ/∂x_i = s_ℓ Σ_k K_ℓk Σⱼ (a_ij/width_i)·∂w_k/∂u_j per joint draw.</item>
/// <item><b>E8's readout</b> runs on the refit in u at a pushforward joint sample and maps to the original inputs as
/// <see cref="BassProjectedFit"/>'s does.</item>
/// </list>
/// The frame is a plug-in, as the PCA basis is.
/// </summary>
public sealed class BassProjectedBasisFit
{
    private BassProjectedBasisFit(BassProjection projection, BassBasisFit fit)
    {
        Projection = projection;
        Fit = fit;
    }

    /// <summary>The frame and the map.</summary>
    public BassProjection Projection { get; }

    /// <summary>The refit, in the active variables, on the native fit's decomposition.</summary>
    public BassBasisFit Fit { get; }

    /// <summary>
    /// The refit of <paramref name="native"/>'s components on <paramref name="projection"/>'s active variables, its decomposition reused,
    /// its component mode and configs carried over (a config per component, or the joint fit's), with the basis fit's master seed.
    /// </summary>
    public static BassProjectedBasisFit Of(BassProjection projection, BassBasisFit native, BassRunOptions? run = null)
    {
        ArgumentNullException.ThrowIfNull(projection);
        ArgumentNullException.ThrowIfNull(native);
        IReadOnlyList<BassPredictor> projected = projection.ProjectPredictors(native.ComponentData(0).Predictors);
        BassBasisFit fit = native.Joint is { } joint
            ? BassComponents.FitJoint(native.Decomposition, projected,
                                      BassProjectedFit.TrainingRange(joint.Config) with { Sampling = joint.Config.Sampling with { MasterSeed = native.MasterSeed } }, run)
            : BassComponents.Fit(native.Decomposition, projected, [.. native.Components.Select(c => BassProjectedFit.TrainingRange(c.Config))], native.MasterSeed, run);
        return new BassProjectedBasisFit(projection, fit);
    }

    /// <summary>The refit's posterior at <paramref name="cases"/>, given in the original inputs' units: the reconstruction.</summary>
    public BassBasisPrediction Predict(PredictionCases cases, BassPredictOptions? options = null) => Fit.Predict(Projection.ProjectCases(cases), options);

    /// <summary>
    /// The reconstruction's partial derivatives at <paramref name="cases"/> in the original continuous inputs, in their own units, per
    /// joint draw (all of them when <paramref name="draws"/> is null): draws × cases × output locations per input, by the chain rule
    /// through the components' gradients in u (<see cref="BassFit.Gradient"/>) and the reconstruction's map.
    /// <see cref="BassGradient.CaseOutsideRange"/> is any component's: a case outside some active variable's training range.
    /// </summary>
    public BassGradient Gradient(PredictionCases cases, IReadOnlyList<int>? draws = null)
    {
        ArgumentNullException.ThrowIfNull(cases);
        int[] joint = draws is null ? [.. Enumerable.Range(0, Fit.JointDraws.Count)] : [.. draws];
        foreach (int b in joint)
            if ((uint)b >= (uint)Fit.JointDraws.Count) throw new ArgumentException($"Joint draw {b} is not among the fit's {Fit.JointDraws.Count}.", nameof(draws));
        PredictionCases projected = Projection.ProjectCases(cases);
        int r = Fit.Components.Count;
        var g = new BassGradient[r];
        for (int k = 0; k < r; k++)
        {
            int kk = k;
            g[k] = Fit.Components[k].Gradient(Fit.ComponentData(k), projected, [.. joint.Select(b => Fit.JointDraws.Draw(kk, b))]);
        }
        LinearReconstruction map = Fit.Decomposition.Map;
        double[,] a = Projection.Map;
        IReadOnlyList<double> widths = Projection.Widths;
        int p = Projection.Inputs.Count, d = Projection.Dimension, m = map.Locations;
        int[] active = [.. Projection.Names.Select(name => IndexOf(g[0].Inputs, name))];
        int n = g[0].Partials[active[0]].GetLength(1);
        var partials = new double[p][,,];
        for (int i = 0; i < p; i++)
        {
            var values = new double[joint.Length, n, m];
            for (int j = 0; j < d; j++)
            {
                double factor = a[i, j] / widths[i];
                if (factor == 0.0) continue;
                for (int k = 0; k < r; k++)
                {
                    double[,,] dw = g[k].Partials[active[j]];
                    for (int l = 0; l < m; l++)
                    {
                        double coefficient = factor * map.Scales[l] * map.Basis(l, k);
                        if (coefficient == 0.0) continue;
                        for (int b = 0; b < joint.Length; b++)
                            for (int c = 0; c < n; c++) values[b, c, l] += coefficient * dw[b, c, 0];
                    }
                }
            }
            partials[i] = values;
        }
        bool[] outside = new bool[n];
        foreach (BassGradient gk in g)
            for (int c = 0; c < n; c++) outside[c] |= gk.CaseOutsideRange[c];
        return new BassGradient(Projection.Inputs, partials, null, Array.AsReadOnly(outside));
    }

    /// <summary>
    /// E8's readout of the refit (<see cref="ProjectedActiveSubspace"/>), as <see cref="BassProjectedFit.ActiveSubspace"/> reads it: the
    /// basis fit's C = E[JᵀW_yJ] over the reconstruction in the active variables, on the sampled route at a pushforward joint sample,
    /// mapped to the original inputs in the analysis's metric; the same rules on the spec.
    /// </summary>
    public ProjectedActiveSubspace ActiveSubspace(ActiveSubspaceSpec? spec = null)
    {
        spec ??= new ActiveSubspaceSpec();
        JointSample sample = BassProjectedFit.Pushforward(Projection, spec);
        ActiveSubspaceAnalysis active = Fit.ActiveSubspace(spec with { Sample = sample, Route = ActiveSubspaceRoute.Sampled });
        return BassProjectedFit.ToInputs(Projection, active, spec);
    }

    /// <summary>
    /// This refit beside <paramref name="native"/>, the basis fit <paramref name="analysis"/> read the frame from (slice 8c,
    /// <see cref="BassProjectionComparison"/>): each fit's size (basis functions per joint draw, summed over the components), PSIS-LOO on
    /// the coefficients, the held-out error of the reconstruction, the refit's subspace against the frame and the native draws' spread
    /// about it. A basis fit has no evidence and no held-out density here (C1). Cross-fitting a basis fit would decompose each fold's
    /// training outputs again (E2), which needs the decomposition's own spec; it is refused.
    /// </summary>
    public BassProjectionComparison Compare(BassBasisFit native, ActiveSubspaceAnalysis analysis, ProjectionComparisonSpec? spec = null)
        => BassProjectionComparisons.Of(this, native, analysis, spec ?? new ProjectionComparisonSpec());

    private static int IndexOf(IReadOnlyList<string> names, string name)
    {
        for (int i = 0; i < names.Count; i++)
            if (names[i] == name) return i;
        throw new InvalidOperationException($"The refit's gradient names no input '{name}'.");
    }
}
