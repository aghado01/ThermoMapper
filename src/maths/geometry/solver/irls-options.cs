// ============================================================================
// Optimization/IrlsOptions.cs
// ============================================================================
namespace Maths.Geometry.Solver
{
    public enum HybridMode
    {
        /// <summary>
        /// Use Weiszfeld steps normally; switch to projected subgradient when the
        /// current iterate falls within SubgradientThreshold of a data point.
        /// </summary>
        Hybrid,

        /// <summary>Force Weiszfeld steps only (useful for testing / profiling).</summary>
        WeiszfeldOnly,

        /// <summary>Force projected subgradient only (useful for testing).</summary>
        SubgradientOnly
    }

    public enum SingularityPolicy
    {
        /// <summary>
        /// Floor the distance to Epsilon before computing the IRLS weight.
        /// Smoothly regularises the singularity; iterate drifts away naturally.
        /// </summary>
        Regularise,

        /// <summary>
        /// Check the Weiszfeld optimality condition at exact coincidence.
        /// If the gradient norm is within the coincident point's weight, the
        /// iterate IS the median; return immediately.
        /// </summary>
        OptimalityCheck
    }

    public enum ConvergenceCriterion
    {
        /// <summary>Converged when ||x_{k+1} - x_k|| &lt; Tolerance.</summary>
        Absolute,

        /// <summary>Converged when ||x_{k+1} - x_k|| / ||x_k|| &lt; Tolerance.</summary>
        RelativeToNorm
    }

    public readonly struct IrlsOptions
    {
        public int MaxIterations { get; init; }
        public double Tolerance { get; init; }

        /// <summary>
        /// Distance floor used by the Regularise singularity policy and the
        /// coincidence check in OptimalityCheck.
        /// </summary>
        public double Epsilon { get; init; }

        public HybridMode HybridMode { get; init; }

        /// <summary>
        /// Distance below which the solver switches to subgradient steps
        /// (Hybrid mode only).
        /// </summary>
        public double SubgradientThreshold { get; init; }

        /// <summary>Initial step size η₀ for the decaying subgradient schedule η_k = η₀/√k.</summary>
        public double Eta0 { get; init; }

        public SingularityPolicy SingularityPolicy { get; init; }
        public ConvergenceCriterion ConvergenceCriterion { get; init; }

        public static IrlsOptions Default => new()
        {
            MaxIterations = 200,
            Tolerance = 1e-8,
            Epsilon = 1e-10,
            HybridMode = HybridMode.Hybrid,
            SubgradientThreshold = 1e-4,
            Eta0 = 1.0,
            SingularityPolicy = SingularityPolicy.OptimalityCheck,
            ConvergenceCriterion = ConvergenceCriterion.Absolute,
        };
    }
}
