using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Declared monitored points (<see cref="BassSampling.MonitoredPoints"/>, decided 2026-10-28): validated in the config,
/// resolved against the data (a direct response's missing location takes the middle one; null takes the defaults,
/// logged), carried in the fingerprint, and read by the run at each point's own case and location.
/// </summary>
public sealed class BassMonitoredPointsTests
{
    [Fact]
    public void Config_ValidatesAndSerializesThePoints()
    {
        BassConfig config = BassRunTests.Config() with
        {
            Sampling = BassRunTests.Config().Sampling with { MonitoredPoints = [new BassMonitoredPoint(3, 1), new BassMonitoredPoint(7, null)] },
        };
        BassConfig back = BassConfig.FromJson(config.ToJson());
        Assert.Equal(config.Sampling.MonitoredPoints, back.Sampling.MonitoredPoints);
        BassData data = BassRunTests.FunctionalData();
        Assert.NotEqual(Bass.Fingerprint(BassRunTests.Config(), data), Bass.Fingerprint(config, data));

        BassConfig With(params BassMonitoredPoint[] points)
            => config with { Sampling = config.Sampling with { MonitoredPoints = points } };
        Assert.Throws<ArgumentException>(() => With().Validate());
        Assert.Throws<ArgumentException>(() => With(new BassMonitoredPoint(-1, null)).Validate());
        Assert.Throws<ArgumentException>(() => With(new BassMonitoredPoint(2, -1)).Validate());
        Assert.Throws<ArgumentException>(() => With(new BassMonitoredPoint(2, 1), new BassMonitoredPoint(2, 1)).Validate());
    }

    [Fact]
    public void Resolve_ChecksThePointsAgainstTheData()
    {
        BassData functional = BassRunTests.FunctionalData();   // 25 cases, 7 locations
        BassData scalar = BassSweepTests.Data(large: true);    // 75 cases
        BassConfig With(params BassMonitoredPoint[] points)
            => BassRunTests.Config() with { Sampling = BassRunTests.Config().Sampling with { MonitoredPoints = points } };

        BassResolved r = BassResolved.Resolve(With(new BassMonitoredPoint(4, 0), new BassMonitoredPoint(9, null)), functional);
        Assert.Equal([new BassMonitoredPoint(4, 0), new BassMonitoredPoint(9, 3)], r.MonitoredPoints);
        Assert.DoesNotContain(r.AutoPicks, p => p.Setting == "MonitoredPoints");
        Assert.Throws<ArgumentException>(() => BassResolved.Resolve(With(new BassMonitoredPoint(25, 0)), functional));
        Assert.Throws<ArgumentException>(() => BassResolved.Resolve(With(new BassMonitoredPoint(4, 7)), functional));
        Assert.Throws<ArgumentException>(() => BassResolved.Resolve(With(new BassMonitoredPoint(4, 3), new BassMonitoredPoint(4, null)), functional));   // null fills in 3
        Assert.Throws<ArgumentException>(() => BassResolved.Resolve(With(new BassMonitoredPoint(4, 0)), scalar));
        Assert.Throws<ArgumentException>(() => BassResolved.Resolve(With(new BassMonitoredPoint(75, null)), scalar));

        BassResolved defaults = BassResolved.Resolve(BassRunTests.Config(), scalar);
        Assert.Equal(BassDiagnostics.DefaultPoints(75, null), defaults.MonitoredPoints);
        Assert.Contains(defaults.AutoPicks, p => p.Setting == "MonitoredPoints");
    }

    [Fact]
    public void Run_MonitorsEachPointAtItsOwnLocation()
    {
        BassMonitoredPoint[] points = [new(2, 0), new(13, 6), new(20, 2)];
        BassConfig config = BassRunTests.Config(chains: 2, samples: 200);
        config = config with { Sampling = config.Sampling with { MonitoredPoints = points } };
        BassData data = BassRunTests.FunctionalData();
        BassFit fit = Bass.Fit(config, data);
        Assert.Equal(points, fit.Diagnostics.Fitted.Select(f => f.Point));

        var model = BassModel.Create(config, data);
        BassInputs inputs = model.Inputs;
        var cache = new Dictionary<int, (double[] A, double[] H)>();
        double At(int id, BassMonitoredPoint point)
        {
            if (!cache.TryGetValue(id, out var columns))
            {
                columns = (new double[inputs.Cases.Cases], new double[inputs.Outputs!.Locations]);
                model.Basis.CaseColumn(fit.Terms[id], inputs.Cases, columns.A);
                model.Basis.OutputColumn(fit.Terms[id], inputs.Outputs, columns.H);
                cache[id] = columns;
            }
            return columns.A[point.Case] * columns.H[point.Location!.Value];
        }
        for (int p = 0; p < points.Length; p++)
        {
            double[][] traces = [.. fit.Draws.Select(chain => chain.Select(d => d.Structure.TermIds
                .Select((id, k) => d.Coefficients[k + 1] * At(id, points[p])).Aggregate(d.Coefficients[0], (f, t) => f + t)).ToArray())];
            ConvergenceDiagnostic reference = ChainDiagnostics.Diagnose(traces), actual = fit.Diagnostics.Fitted[p].Diagnostic;
            Assert.Equal(reference.RHat.Value, actual.RHat.Value, 9);
            Assert.Equal(reference.BulkEss.Value, actual.BulkEss.Value, 6);
        }
    }
}
