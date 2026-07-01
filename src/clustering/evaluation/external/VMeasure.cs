using System;

namespace Clustering.Evaluation.External;

/// <summary>
/// V-measure: harmonic mean of <see cref="Homogeneity"/> and
/// <see cref="Completeness"/>.
/// <code>V = 2·h·c / (h + c)</code>
/// </summary>
/// <remarks>
/// <para><b>Range.</b> <c>[0, 1]</c>; <c>1.0</c> requires both perfect
/// homogeneity and perfect completeness simultaneously.
/// <b>Higher is better.</b></para>
///
/// <para><b>Symmetric.</b> Unlike its two constituents, V-measure is
/// symmetric in its arguments — homogeneity and completeness swap
/// under argument swap, and the harmonic mean is invariant to that
/// swap.</para>
///
/// <para><b>When to prefer V-measure.</b> Use this when you want a
/// single scalar that punishes both over-clustering (low homogeneity)
/// and under-clustering (low completeness). For asymmetric questions —
/// "does my clustering preserve class identity but it's OK if it
/// splits classes across clusters" — use <see cref="Homogeneity"/>
/// alone.</para>
///
/// <para><b>Edge cases.</b> Returns 0.0 when either component is 0
/// (the harmonic mean of 0 and anything is 0). Returns 0.0 for empty
/// input.</para>
/// </remarks>
public sealed class VMeasure : IExternalClusterEvaluator
{
    public string Name => "VMeasure";

    private static readonly Homogeneity HomogeneityEvaluator = new();
    private static readonly Completeness CompletenessEvaluator = new();

    public double Evaluate(int[] predictedLabels, int[] referenceLabels)
    {
        double h = HomogeneityEvaluator.Evaluate(predictedLabels, referenceLabels);
        double c = CompletenessEvaluator.Evaluate(predictedLabels, referenceLabels);

        double denom = h + c;
        if (denom <= 0.0) return 0.0;

        return Math.Clamp(2.0 * h * c / denom, 0.0, 1.0);
    }
}
