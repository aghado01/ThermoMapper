namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// A Legendre fit's sensitivity read from its coefficients (KHAOS slice 4; <see cref="BassSensitivity.FromCoefficients"/>), under
/// the product uniform law on each input's box: its declared domain, or its training range when it declares none. On that law the
/// ψ_α are orthonormal with mean 0 for α ≠ 0, so for f = b₀ + Σ_α b_α ψ_α (Lean SG-4(b), <c>BassIntegrals.lean</c>):
/// <list type="bullet">
/// <item><b><see cref="Sobol"/>,</b> the variance decomposition, the moment route's type: Var f = Σ_α b_α², a subset u's main
/// variance the sum over the terms whose variables are exactly u, and an input's total variance Σ_{α_v ≥ 1} b_α², which is
/// E[(f − E_v f)²] itself (<c>tensorTotal_eq</c>).</item>
/// <item><b><see cref="Jacobi"/>,</b> the Jacobi-weighted derivative measure E[x_v(1 − x_v)(∂_v f)²] = Σ_α b_α² α_v(α_v + 1)
/// (<c>integral_jacobi_tensorPartial</c>), in the box's units; under the affine map to raw units it reads
/// E[(z − a)(b − z)(∂_z f)²], the same number. Each total is at most half its Jacobi measure, with equality where every term
/// that holds the input has degree 1 in it (<c>tensorTotal_le</c>, <c>tensorTotal_eq_half</c>).</item>
/// </list>
/// <para>These are exact only under that law and on this basis. Under another measure (<see cref="MixtureMeasure"/>,
/// <see cref="EmpiricalMeasure"/>, a uniform on another interval) the ψ_α are not orthonormal, and the decomposition is the moment
/// route's (<see cref="BassSensitivity.Of"/>), exact through each measure's interval moments; the derivative measure under any
/// measure, a dependent one included, is the gradient readout's (<see cref="BassFit.Gradient"/>) averaged over the measure.</para>
/// </summary>
/// <param name="Sobol">The decomposition, its measure the product uniform law on the boxes.</param>
/// <param name="Jacobi">Draws × inputs: E[x_v(1 − x_v)(∂_v f)²] for each input v in <see cref="SobolAnalysis.Measure"/>'s numbering.</param>
public sealed record PolynomialChaosSensitivity(SobolAnalysis Sobol, double[,] Jacobi);
