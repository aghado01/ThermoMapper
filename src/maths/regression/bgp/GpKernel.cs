using System;

namespace Maths.Regression.Bgp;

/// <summary>
/// A radial GP covariance kernel as a function of the squared <em>ambient</em> distance ‖x−x'‖² and a bandwidth t:
/// k_t(x,x') = h(‖x−x'‖²/t). The kernel is computed from the high-D Euclidean coordinates directly — it carries no
/// knowledge of any low-dimensional structure; adaptivity to the intrinsic dimension comes entirely from the prior
/// on t (Tang, Wu, Cheng &amp; Dunson 2025). The theory holds for a class of h (their Assumption A.3); the
/// squared-exponential is the representative example.
/// </summary>
public interface IGpKernel
{
    /// <summary>k_t value from the squared ambient distance <paramref name="sqDistance"/> = ‖x−x'‖² and bandwidth t.</summary>
    double Evaluate(double sqDistance, double t);
}

/// <summary>Squared-exponential (Gaussian RBF) covariance k_t(x,x') = exp(−‖x−x'‖²/(2t)) (TWCD2025 eq. 2).</summary>
public sealed class SquaredExponentialKernel : IGpKernel
{
    public double Evaluate(double sqDistance, double t) => Math.Exp(-sqDistance / (2.0 * t));
}
