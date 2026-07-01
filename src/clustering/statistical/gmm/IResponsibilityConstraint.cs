namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// Blends an external supervisor's N×K confidence matrix sᵢₖ into the EM
    /// E-step output: r̂ᵢₖ ← (1 − λ) · r̂ᵢₖ + λ · sᵢₖ, each row renormalised.
    /// See docs/gmm.md for the semi-supervised path.
    /// </summary>
    public interface IResponsibilityConstraint
    {
        /// <summary>
        /// Blends <paramref name="responsibilities"/> in-place with the external
        /// signal for the current EM iteration. Each row must sum to 1 on return.
        /// </summary>
        /// <param name="responsibilities">N×K E-step output, modified in-place.</param>
        /// <param name="n">Number of data points.</param>
        /// <param name="k">Number of components.</param>
        /// <param name="iteration">Zero-based iteration index within the current Fit call.</param>
        /// <param name="maxIterations">Total iteration budget for this Fit call.</param>
        void Apply(double[,] responsibilities, int n, int k, int iteration, int maxIterations);
    }
}
