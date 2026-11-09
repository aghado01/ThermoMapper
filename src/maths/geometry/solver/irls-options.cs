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
        /// <summary>Converged when the shift is below Tolerance, read in the <see cref="IrlsScale"/> unit λ.</summary>
        Absolute,

        /// <summary>Converged when the shift is below Tolerance·(λ + ||x_{k+1}||), with ||·|| the ambient norm.</summary>
        RelativeToNorm
    }

    /// <summary>
    /// The unit λ in which the distance-valued constants (Tolerance, Epsilon, SubgradientThreshold, Eta0) are read.
    /// </summary>
    public enum IrlsScale
    {
        /// <summary>
        /// λ = σ, the spread: the weighted median distance from the current iterate to the data, recomputed each
        /// iteration (their weighted mean when half the weight sits on the iterate), so the median and mean are
        /// equivariant under rescaling the data. The radii are read from the Log map, the shift is the norm of the
        /// tangent step, flat steps accumulate the differences xᵢ − p, the subgradient step is also divided by Σw, and
        /// the solver stops as well when a step no longer moves the iterate.
        /// </summary>
        Spread,

        /// <summary>
        /// λ = 1: the constants are in the data's own units, the solver's semantics before the spread unit.
        /// </summary>
        Data
    }

    public readonly struct IrlsOptions
    {
        public int MaxIterations { get; init; }

        /// <summary>Stopping radius for the shift, in units of λ (see <see cref="Scale"/>).</summary>
        public double Tolerance { get; init; }

        /// <summary>
        /// Coincidence radius, in units of λ: the distance floor used by the Regularise singularity
        /// policy and the coincidence check in OptimalityCheck.
        /// </summary>
        public double Epsilon { get; init; }

        public HybridMode HybridMode { get; init; }

        /// <summary>
        /// Distance below which the solver switches to subgradient steps, in units of λ
        /// (Hybrid mode only).
        /// </summary>
        public double SubgradientThreshold { get; init; }

        /// <summary>Initial step size η₀ for the decaying subgradient schedule η_k = λη₀/√(k+1).</summary>
        public double Eta0 { get; init; }

        public SingularityPolicy SingularityPolicy { get; init; }
        public ConvergenceCriterion ConvergenceCriterion { get; init; }

        /// <summary>The unit of the distance-valued constants; <see cref="IrlsScale.Spread"/> by default.</summary>
        public IrlsScale Scale { get; init; }

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
            Scale = IrlsScale.Spread,
        };
    }
}
