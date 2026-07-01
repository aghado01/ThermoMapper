using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The skewed-Poisson β refinement (WLS2008 §5.4): for small Poisson means the conditional posterior of the
/// log-rate is skewed, so the Laplace/normal approximation (centered at the MLE) misplaces its mean. The
/// independence-MH correction in <c>PoissonModel.SampleCoefficients</c> recovers the true posterior mean.
/// Tested on an intercept-only model where the 1-D posterior mean is available by quadrature.
/// </summary>
public sealed class SkewedPoissonTests
{
    private readonly ITestOutputHelper _out;
    public SkewedPoissonTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void SkewedPoisson_BetaDraw_RecoversTruePosteriorMean()
    {
        const int n = 12;
        var design = new double[n, 1];
        for (int i = 0; i < n; i++) design[i, 0] = 1.0;        // intercept-only ⇒ β is the log-rate
        var y = new double[] { 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 };   // small counts, S = 3 → skewed
        double s = 0.0;
        foreach (double v in y) s += v;

        // True posterior over β: π(β) ∝ exp(Sβ − n·e^β − β²/200) — the same target SampleCoefficients uses
        // (Poisson log-likelihood + the weak ridge prior). Mean by fine quadrature.
        double LogPost(double b) => s * b - n * Math.Exp(b) - 0.5 * b * b / 100.0;
        const double lo = -8.0, hi = 3.0;
        const int g = 200_000;
        double dh = (hi - lo) / g;
        double maxlp = double.NegativeInfinity;
        for (int i = 0; i <= g; i++) maxlp = Math.Max(maxlp, LogPost(lo + i * dh));
        double z = 0.0, weighted = 0.0;
        for (int i = 0; i <= g; i++)
        {
            double b = lo + i * dh;
            double p = Math.Exp(LogPost(b) - maxlp);
            z += p;
            weighted += b * p;
        }
        double trueMean = weighted / z;
        double mle = Math.Log(s / n);                          // MLE = posterior mode

        // Full independence-MH (suspect shortcut disabled via +∞ threshold) ⇒ draws from the true posterior.
        var model = new PoissonModel();
        var rng = new Xoshiro256PlusPlus(seed: 7);
        const int draws = 20_000;
        double sum = 0.0;
        for (int d = 0; d < draws; d++)
            sum += model.SampleCoefficients(design, y, rng, mhIterations: 30, suspectThreshold: double.PositiveInfinity)[0];
        double mhMean = sum / draws;

        _out.WriteLine($"[skew-poisson] trueMean={trueMean:F4} MLE(β̂)={mle:F4} MH mean={mhMean:F4} " +
                       $"(|MH−true|={Math.Abs(mhMean - trueMean):F4}, |MLE−true|={Math.Abs(mle - trueMean):F4})");

        // The MH-corrected mean is closer to the true (skewed) posterior mean than the MLE/normal-approx center,
        // and recovers it.
        Assert.True(Math.Abs(mhMean - trueMean) < Math.Abs(mle - trueMean),
            $"MH mean {mhMean:F4} should beat the MLE {mle:F4} as an estimate of the true mean {trueMean:F4}");
        Assert.True(Math.Abs(mhMean - trueMean) < 0.03,
            $"MH mean {mhMean:F4} should recover the true posterior mean {trueMean:F4}");
    }
}
