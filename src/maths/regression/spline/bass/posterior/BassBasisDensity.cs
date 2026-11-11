using System;
using System.Collections.Generic;
using Maths.LinAlg;

namespace Maths.Regression.Spline.Bass;

/// <summary>Where a basis fit's pointwise density is taken (OQ17, D07).</summary>
public enum BasisDensitySpace
{
    /// <summary>
    /// The default: coefficient space, Σ_k log N(w_ik | w_k(xᵢ), σ_k²) over a unit's coefficients. For a fixed basis it
    /// differs from the density in y-space by a constant and a residual term the same for every fit, so it is exact
    /// for comparing fits that share a basis.
    /// </summary>
    Coefficients,

    /// <summary>
    /// A declared approximation for comparing across bases or engines: y-space in physical units, each case
    /// y_i ~ N(y_m + y_s ∘ (K w(xᵢ) + ē), diag(y_s)(K diag(σ²) Kᵀ + diag(v))diag(y_s)), the truncation residual modelled as
    /// a Gaussian with the declared residual set's mean ē and diagonal variance v (calibration's functional form).
    /// </summary>
    ResponseWithDiagonalTruncation,
}

/// <summary>
/// A basis fit's pointwise log density (OQ17; D07), over its joint draws, on the shared <see cref="PointwiseLogDensity"/>
/// seam the <see cref="PredictiveChecks"/> read. The rows follow the joint draws, which keep the components' chains when
/// they were paired by index, so PSIS reads each unit's relative efficiency from them.
/// <list type="bullet">
/// <item><b>Coefficient space</b> (<see cref="BasisDensitySpace.Coefficients"/>, the default): the observations are the
/// n·r coefficient responses w_ik, case by case and, within a case, component by component; by default each case is a
/// unit. The retained parameters' content: every structure, coefficient and σ_k² is conditioned on.</item>
/// <item><b>y-space under the diagonal Gaussian</b>: the observations are the n·m training responses, case by case and
/// location by location, read back from W and the in-sample residuals; a unit must hold whole cases, since the
/// locations of a case are correlated through K. The component noise and the truncation residual of a case are
/// integrated; the density is refused where a location's truncation variance is zero, as at a constant column.</item>
/// <item><b>Under the joint option</b> (E34's decision 8, <see cref="BassBasisFit.Joint"/>) a case's component noises share Σ:
/// in coefficient space a case's density is N(wᵢ | w(xᵢ), Σ), so a unit holds whole cases there too, and in y-space the
/// component term is KΣKᵀ.</item>
/// </list>
/// </summary>
internal static class BassBasisDensity
{
    public static PointwiseLogDensity LogDensity(BassBasisFit fit, UnitPartition? units, BasisDensitySpace space)
    {
        if (space is not (BasisDensitySpace.Coefficients or BasisDensitySpace.ResponseWithDiagonalTruncation))
            throw new ArgumentOutOfRangeException(nameof(space));
        BassDecomposition decomposition = fit.Decomposition;
        BassJointDraws joint = fit.JointDraws;
        int n = decomposition.Cases, r = decomposition.Rank, m = decomposition.Locations, s = joint.Count;
        bool coefficients = space == BasisDensitySpace.Coefficients;
        int width = coefficients ? r : m;
        UnitPartition partition = units ?? UnitPartition.ConsecutiveBlocks(n, width);
        if (partition.Observations != n * width)
            throw new ArgumentException($"The partition must cover the {n * width} observations, case by case and {(coefficients ? "component by component" : "location by location")}.",
                                        nameof(units));

        // Each component's fit at the training cases under each joint draw, and its σ².
        var fitted = new double[s, n, r];
        var sigma2 = new double[s, r];
        var f = new double[n];
        bool samplesTau = false;
        for (int k = 0; k < r; k++)
        {
            BassFit component = fit.Components[k];
            BassModel model = fit.Model(k);
            samplesTau |= model.Target.SamplesTau;
            (List<BassDraw> pooled, _) = BassDrawFit.Pool(component.Draws, component.Terms, false, "fit");
            var drawFit = new BassDrawFit(model, component.Terms);
            for (int d = 0; d < s; d++)
            {
                BassDraw draw = pooled[joint.Draw(k, d)];
                sigma2[d, k] = draw.Sigma2;
                drawFit.Evaluate(draw, f);
                for (int i = 0; i < n; i++) fitted[d, i, k] = f[i];
            }
        }

        // Under the joint option the components' noises share Σ (E34's decision 8): each joint draw's factor L, LLᵀ = Σ.
        double[][,]? roots = null;
        if (fit.Joint is not null)
        {
            roots = new double[s][,];
            for (int d = 0; d < s; d++) roots[d] = fit.NoiseRoot(d)!;
        }

        double[,] values = coefficients ? Coefficients(decomposition, fitted, sigma2, roots, partition)
                                        : Response(decomposition, fitted, sigma2, roots, partition);
        var latents = new List<LatentReport>
        {
            new(roots is null ? "component structures" : "the components' shared structure", LatentScope.Global, LatentTreatment.Conditioned),
        };
        if (samplesTau) latents.Add(new LatentReport(roots is null ? "component prior precisions τ_k" : "the shared prior precision τ", LatentScope.Global, LatentTreatment.Conditioned));
        latents.Add(new LatentReport("component coefficients", LatentScope.Global, LatentTreatment.Conditioned));
        latents.Add(new LatentReport(roots is null ? "component noise variances σ_k²" : "the components' noise covariance Σ", LatentScope.Global, LatentTreatment.Conditioned));
        if (!coefficients)
        {
            latents.Add(new LatentReport("component noise ε", LatentScope.Group, LatentTreatment.Integrated));
            latents.Add(new LatentReport("truncation residual e (diagonal Gaussian)", LatentScope.Group, LatentTreatment.Integrated));
        }
        return new PointwiseLogDensity(values, partition, latents, joint.ChainLengths);
    }

    // Σ over each unit's (case, component) observations of log N(w_ik | w_k(xᵢ), σ_k²); under the joint option, over each unit's
    // cases of log N(wᵢ | w(xᵢ), Σ), so a unit holds whole cases.
    private static double[,] Coefficients(BassDecomposition decomposition, double[,,] fitted, double[,] sigma2, double[][,]? roots,
                                          UnitPartition partition)
    {
        if (roots is not null) return JointCoefficients(decomposition, fitted, roots, partition);
        int s = fitted.GetLength(0), n = fitted.GetLength(1), r = fitted.GetLength(2);
        var values = new double[s, partition.Units];
        var log = new double[n * r];
        for (int d = 0; d < s; d++)
        {
            for (int i = 0; i < n; i++)
                for (int k = 0; k < r; k++)
                {
                    double v = sigma2[d, k], e = decomposition.Coefficient(k, i) - fitted[d, i, k];
                    log[i * r + k] = -0.5 * Math.Log(2.0 * Math.PI * v) - 0.5 * e * e / v;
                }
            for (int u = 0; u < partition.Units; u++)
            {
                double sum = 0.0;
                foreach (int o in partition.Members[u]) sum += log[o];
                values[d, u] = sum;
            }
        }
        return values;
    }

    // The joint option's coefficient-space density: log N(wᵢ | w(xᵢ), Σ) = −(r/2) ln 2π − Σ_k ln L_kk − ½‖L⁻¹(wᵢ − w(xᵢ))‖², per case.
    private static double[,] JointCoefficients(BassDecomposition decomposition, double[,,] fitted, double[][,] roots, UnitPartition partition)
    {
        int s = fitted.GetLength(0), n = fitted.GetLength(1), r = fitted.GetLength(2);
        int[] caseUnit = new int[n];
        for (int i = 0; i < n; i++)
        {
            caseUnit[i] = partition.UnitOf[i * r];
            for (int k = 1; k < r; k++)
                if (partition.UnitOf[i * r + k] != caseUnit[i])
                    throw new ArgumentException($"Case {i}'s components fall in more than one unit; under the joint option their noises share Σ, so a unit holds whole cases.", "units");
        }
        var values = new double[s, partition.Units];
        var u = new double[r];
        for (int d = 0; d < s; d++)
        {
            double[,] root = roots[d];
            double logDet = 0.0;
            for (int k = 0; k < r; k++) logDet += Math.Log(root[k, k]);
            double constant = -0.5 * r * Math.Log(2.0 * Math.PI) - logDet;
            for (int i = 0; i < n; i++)
            {
                double q = 0.0;
                for (int k = 0; k < r; k++)   // u = L⁻¹e by forward substitution
                {
                    double e = decomposition.Coefficient(k, i) - fitted[d, i, k];
                    for (int j = 0; j < k; j++) e -= root[k, j] * u[j];
                    u[k] = e / root[k, k];
                    q += u[k] * u[k];
                }
                values[d, caseUnit[i]] += constant - 0.5 * q;
            }
        }
        return values;
    }

    // Each case's y-space log density under the diagonal Gaussian, summed over a unit's cases: the standardized
    // covariance Σ = diag(v) + K diag(σ²) Kᵀ held as a diagonal plus a low-rank term (calibration's form), so no m × m
    // matrix is formed; the physical density subtracts Σ log y_s. Under the joint option the low-rank term is (KL)(KL)ᵀ, LLᵀ the
    // components' noise covariance.
    private static double[,] Response(BassDecomposition decomposition, double[,,] fitted, double[,] sigma2, double[][,]? roots, UnitPartition partition)
    {
        TruncationResiduals truncation = decomposition.Truncation;
        double[,] inSample = truncation.InSample
            ?? throw new ArgumentException("The y-space density needs the training responses, which supplied residuals do not carry.", "space");
        LinearReconstruction map = decomposition.Map;
        double[,] k = map.BasisArray;
        int s = fitted.GetLength(0), n = fitted.GetLength(1), r = fitted.GetLength(2), m = map.Locations;

        int[] caseUnit = new int[n];
        for (int i = 0; i < n; i++)
        {
            caseUnit[i] = partition.UnitOf[i * m];
            for (int l = 1; l < m; l++)
                if (partition.UnitOf[i * m + l] != caseUnit[i])
                    throw new ArgumentException($"Case {i}'s locations fall in more than one unit; under K they are correlated, so a unit holds whole cases.", "units");
        }

        double[] mean = truncation.Mean();
        double[,] covariance = truncation.Covariance();
        var v = new double[m];
        double logScales = 0.0;
        for (int l = 0; l < m; l++)
        {
            v[l] = covariance[l, l];
            if (!(v[l] > 0.0))
                throw new ArgumentException($"The truncation variance at location {l} is zero, so the diagonal Gaussian has no density there.", "space");
            logScales += Math.Log(map.Scales[l]);
        }
        // Each case's standardized response, Kwᵢ + eᵢ (the in-sample residual reconstructs it).
        var yc = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++)
            {
                double sum = inSample[i, l];
                for (int q = 0; q < r; q++) sum += k[l, q] * decomposition.Coefficient(q, i);
                yc[i, l] = sum;
            }

        var values = new double[s, partition.Units];
        var x = new double[m];
        var weights = new double[r];
        double constant = -0.5 * m * Math.Log(2.0 * Math.PI) - logScales;
        double[,]? kl = roots is null ? null : new double[m, r];
        for (int d = 0; d < s; d++)
        {
            DiagonalPlusLowRank sigma;
            if (roots is null)
            {
                for (int q = 0; q < r; q++) weights[q] = sigma2[d, q];
                sigma = new DiagonalPlusLowRank(v, k, weights);
            }
            else
            {
                double[,] root = roots[d];
                for (int l = 0; l < m; l++)
                    for (int q = 0; q < r; q++)
                    {
                        double a = 0.0;
                        for (int j = q; j < r; j++) a += k[l, j] * root[j, q];
                        kl![l, q] = a;
                    }
                Array.Fill(weights, 1.0);
                sigma = new DiagonalPlusLowRank(v, kl!, weights);
            }
            for (int i = 0; i < n; i++)
            {
                for (int l = 0; l < m; l++)
                {
                    double reconstructed = mean[l];
                    for (int q = 0; q < r; q++) reconstructed += k[l, q] * fitted[d, i, q];
                    x[l] = yc[i, l] - reconstructed;
                }
                values[d, caseUnit[i]] += constant - 0.5 * (sigma.LogDet + sigma.InverseQuadraticForm(x));
            }
        }
        return values;
    }
}
