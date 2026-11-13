using System;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Subspaces;
using Xunit;

namespace Maths.Regression.Tests;

public sealed class RelativeEigenMigrationTests
{
    [Fact]
    public void FirstCallers_AgreeWithUnitScaleOnWellConditionedFixtures()
    {
        double[,] a = { { 2, 0.2, 0.3 }, { 0.2, 3, -0.4 }, { 0.3, -0.4, 4 } };
        foreach (int exponent in new[] { -240, -120, 0, 120, 240 })
        {
            double scale = Math.Pow(10, exponent);
            var scaled = new double[3,3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) scaled[i,j] = scale * a[i,j];
            EigenResult old = ScaledEigen.Of(scaled);
            EigenResult now = DenseEigen.DecomposeSymmetric(scaled, tol: 1e-14, stoppingRule: JacobiStoppingRule.ElementWiseRelative);
            for (int j = 0; j < 3; j++) Assert.InRange(Math.Abs(old.Eigenvalues[j] / now.Eigenvalues[j] - 1), 0, 2e-13);
        }

        // A nontrivial diagonal prior, whose independently whitened likelihood is a.
        double[] d = [2, 3, 4], cross = [0.3, -0.2, 0.5];
        var prior = new double[3,3];
        var likelihood = new double[3,3];
        var rawCross = new double[3];
        var flat = new double[9];
        for (int i = 0; i < 3; i++)
        {
            prior[i,i] = d[i];
            rawCross[i] = Math.Sqrt(d[i]) * cross[i];
            for (int j = 0; j < 3; j++)
            {
                likelihood[i,j] = Math.Sqrt(d[i]) * a[i,j] * Math.Sqrt(d[j]);
                flat[j*3+i] = a[i,j];
            }
        }
        PriorPencil pencil = Assert.IsType<PriorPencil>(PriorPencil.Of(likelihood, prior, rawCross, 2, true));
        EigenResult legacy = PriorPencil.DecomposeAtUnitScale(flat, 3, 4);
        foreach (double tau in new[] { 1e-8, 0.1, 1, 100, 1e8 })
            foreach (double b in new[] { 0.0, 0.4, 1.0 })
            {
                double logDet = 0, shrink = 0;
                for (int j = 0; j < 3; j++)
                {
                    double w = 0;
                    for (int i = 0; i < 3; i++) w += legacy.Eigenvectors[j][i] * cross[i];
                    double lambda = legacy.Eigenvalues[j];
                    logDet -= 0.5 * double.LogP1(b * lambda / tau);
                    shrink += w*w/lambda * (tau/(tau+b*lambda));
                }
                RidgeTerms terms = pencil.At(tau,b);
                Assert.InRange(Math.Abs(terms.LogDeterminant-logDet),0,2e-13);
                Assert.InRange(Math.Abs(terms.Shrinkage-shrink),0,2e-13);
                GPriorMarginal marginal = CollapsedGPrior.CentredPencilMarginalAt(pencil, 5, 10, GPriorWeights.Unit(10), GPriorCentre.Zero, tau, b);
                // Same conditional log evidence from the former pencil's determinant and shrinkage, sigma² = 2.
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
