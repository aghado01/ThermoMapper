using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// S10's support-rule fixture (D06, owed since S4): the count rule, BASS's default, against column mass, at the data's
/// edges. A hinge whose k nonzero points sit at an edge has count k but mass about k/2 (its values rise from 0 to its
/// normalized maximum of one), so at one threshold b the mass rule asks an edge hinge for about 2b points, while
/// spread-out hinges pass both. Count admits the more local edge terms, which Francom et al. 2018 call unstable; mass
/// refuses them, and then an edge feature narrower than 2b points is out of reach.
/// <list type="bullet">
///   <item><b>Data:</b> n = 200 on (i + ½)/n, noise sd 0.3 and 0.02, the threshold R's min(20, 0.1·n) = 20. A smooth surface,
///   sin 2πx; and an edge kink, sin 2πx + 8(x − 0.88)₊, whose 24 points meet the count but whose mass, about 12, does
///   not.</item>
///   <item><b>Fits:</b> four chains of 2,000 burn-in and 2,000 draws under each rule.</item>
///   <item><b>Read:</b> on the edge regions (x below 0.1 and above 0.9) and the interior, the posterior sd of f
///   averaged over the training points, the posterior mean's RMSE against the truth, and the 95% intervals' coverage of
///   it; with the mean M.</item>
/// </list>
/// Expectations, recorded before the run: on the smooth surface, mass reads a lower edge sd at a like RMSE; on the kink,
/// count reads the lower upper-edge RMSE. The first run, at sd 0.3, separated nothing: the posterior kept three or four
/// broad terms, which neither threshold binds. The sd 0.02 arm was added after it, for the low-noise regime of computer
/// models where Francom et al. place the edge instability. A qualification run (<see cref="BassQualification.Enabled"/>).
/// </summary>
public sealed class BassQualificationColumnMassTests
{
    private readonly ITestOutputHelper _out;
    public BassQualificationColumnMassTests(ITestOutputHelper output) => _out = output;

    private const int N = 200;

    [Fact]
    public void CountAgainstColumnMass_AtTheEdges()
    {
        if (!BassQualification.Enabled) return;   // heavy
        double[] x = [.. Enumerable.Range(0, N).Select(i => (i + 0.5) / N)];
        foreach (double noise in new[] { 0.3, 0.02 })
        foreach ((string surface, Func<double, double> f) in new (string, Func<double, double>)[]
                 {
                     ("smooth", t => Math.Sin(2 * Math.PI * t)),
                     ("edge kink", t => Math.Sin(2 * Math.PI * t) + 8 * Math.Max(t - 0.88, 0)),
                 })
        {
            string name = $"{surface}, sd {noise}";
            var rng = new Xoshiro256PlusPlus(seed: 41);
            double[] truth = [.. x.Select(f)];
            double[] y = [.. truth.Select(v => v + noise * BassQualification.Gaussian(rng))];
            var data = new BassData([new ContinuousPredictor("x", x)], new ScalarResponse(y));
            foreach (SupportKind kind in new[] { SupportKind.Count, SupportKind.ColumnMass })
            {
                var config = new BassConfig
                {
                    Structure = new BassStructureSpec { Support = new SupportRule(kind) },
                    Sampling = new BassSampling { Chains = 4, MasterSeed = 43, BurnIn = 2000, Samples = 2000 },
                };
                BassFit fit = Bass.Fit(config, data);
                double[][] fitted = BassQualification.Fitted(BassModel.Create(config, data), fit);
                var at = Enumerable.Range(0, N).Select(i => BassQualification.At(fitted, i)).ToArray();
                string Region(Func<double, bool> inside)
                {
                    int[] points = [.. Enumerable.Range(0, N).Where(i => inside(x[i]))];
                    double sd = points.Average(i => at[i].Sd);
                    double rmse = Math.Sqrt(points.Average(i => (at[i].Mean - truth[i]) * (at[i].Mean - truth[i])));
                    double covered = points.Count(i => at[i].Lower <= truth[i] && truth[i] <= at[i].Upper) / (double)points.Length;
                    Assert.True(double.IsFinite(sd) && double.IsFinite(rmse));
                    return $"sd {sd:F3}, RMSE {rmse:F3}, covered {covered:P0}";
                }
                double count = fit.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length);
                _out.WriteLine($"{name}, {kind}: mean M {count:F2}; lower edge {Region(t => t < 0.1)}; upper edge {Region(t => t > 0.9)}; " +
                               $"interior {Region(t => t >= 0.1 && t <= 0.9)}");
            }
        }
    }
}
