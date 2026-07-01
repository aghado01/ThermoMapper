namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// One row of a <see cref="BicKSweep"/> result: the K value, its BIC score
    /// (BIC = −2·logL + p·ln(N)), the converged log-likelihood, iteration count,
    /// convergence flag, and the fitted model.
    /// </summary>
    public sealed record KSweepResult(
        int K,
        double Bic,
        double LogLikelihood,
        int NumIterations,
        bool IsConverged,
        GaussianMixtureModel Model);
}
