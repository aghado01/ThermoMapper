using System;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Subspaces;
using Xunit;

namespace Maths.Regression.Tests;

public sealed class RelativeEigenMigrationTests
{
    /// <summary>
    /// The relative rule follows any rescaling (spectral-solver brief, B1): c·A's eigenvalues divided by c match A's from 10⁻²⁴⁰ to
    /// 10²⁴⁰, which the retired unit-scale division used to supply. A prior pencil whose independently whitened likelihood is A reads
    /// its determinant, shrinkage and evidence from A's spectrum.
    /// </summary>
    [Fact]
    public void RelativeSolve_FollowsScale_AndThePencilReadsTheWhitenedSpectrum()
    {
        double[,] a = { { 2, 0.2, 0.3 }, { 0.2, 3, -0.4 }, { 0.3, -0.4, 4 } };
        EigenResult spectrum = DenseEigen.DecomposeSymmetric(a, tol: 1e-14);
        foreach (int exponent in new[] { -240, -120, 0, 120, 240 })
        {
            double scale = Math.Pow(10, exponent);
            var scaled = new double[3,3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) scaled[i,j] = scale * a[i,j];
            EigenResult now = DenseEigen.DecomposeSymmetric(scaled, tol: 1e-14);
            for (int j = 0; j < 3; j++) Assert.InRange(Math.Abs(now.Eigenvalues[j] / (scale * spectrum.Eigenvalues[j]) - 1), 0, 2e-13);
        }

        // A nontrivial diagonal prior, whose independently whitened likelihood is a.
        double[] d = [2, 3, 4], cross = [0.3, -0.2, 0.5];
        var prior = new double[3,3];
        var likelihood = new double[3,3];
        var rawCross = new double[3];
        for (int i = 0; i < 3; i++)
        {
            prior[i,i] = d[i];
            rawCross[i] = Math.Sqrt(d[i]) * cross[i];
            for (int j = 0; j < 3; j++)
            {
                likelihood[i,j] = Math.Sqrt(d[i]) * a[i,j] * Math.Sqrt(d[j]);
            }
        }
        PriorPencil pencil = Assert.IsType<PriorPencil>(PriorPencil.Of(likelihood, prior, rawCross, 2, true));
        foreach (double tau in new[] { 1e-8, 0.1, 1, 100, 1e8 })
            foreach (double b in new[] { 0.0, 0.4, 1.0 })
            {
                double logDet = 0, shrink = 0;
                for (int j = 0; j < 3; j++)
                {
                    double w = 0;
                    for (int i = 0; i < 3; i++) w += spectrum.Eigenvectors[j][i] * cross[i];
                    double lambda = spectrum.Eigenvalues[j];
                    logDet -= 0.5 * double.LogP1(b * lambda / tau);
                    shrink += w*w/lambda * (tau/(tau+b*lambda));
                }
                RidgeTerms terms = pencil.At(tau,b);
                Assert.InRange(Math.Abs(terms.LogDeterminant-logDet),0,2e-13);
                Assert.InRange(Math.Abs(terms.Shrinkage-shrink),0,2e-13);
                GPriorMarginal marginal = CollapsedGPrior.CentredPencilMarginalAt(pencil, 5, 10, GPriorWeights.Unit(10), GPriorCentre.Zero, tau, b);
                // The conditional log evidence from the spectrum's determinant and shrinkage, sigma² = 2.
                double oldEvidence = -0.5*b*10*Math.Log(4*Math.PI) + logDet - 0.25*b*(5+shrink);
                Assert.InRange(Math.Abs(marginal.LogMarginalGiven(2)-oldEvidence),0,2e-13);
            }
    }

    [Fact]
    public void PriorPencil_WeakDirectionControlsShrinkageAndEvidence()
    {
        double[,] a = { { 1, 0, 0 }, { 0, 1e-20, 0.5e-20 }, { 0, 0.5e-20, 1e-20 } };
        double[,] identity = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        PriorPencil pencil = Assert.IsType<PriorPencil>(PriorPencil.Of(a, identity, [0, 5e-11, -5e-11], 2, true));
        Assert.InRange(Math.Abs(pencil.Eigenvalues[2]/0.5e-20-1),0,1e-14);
        foreach (double tau in new[] { 1e-22, 1e-20, 1.0 })
        {
            RidgeTerms terms = pencil.At(tau);
            double expected = tau/(tau+0.5e-20);
            Assert.InRange(Math.Abs(terms.Shrinkage/expected-1),0,2e-14);
            double logDet = -0.5*(double.LogP1(1/tau)+double.LogP1(1.5e-20/tau)+double.LogP1(0.5e-20/tau));
            Assert.InRange(Math.Abs(terms.LogDeterminant-logDet),0,2e-14);
        }
    }
}
