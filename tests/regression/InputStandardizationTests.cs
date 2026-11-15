using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The transform's <see cref="Standardization"/> beside each continuous input's box (E36's decision 3): the training mean and sample
/// deviation against a compensated two-pass computation, and <see cref="BassInputs.Standardizing"/>, the map z = a + b·u from the
/// scaled value, against (x − mean)/deviation at the training cases, under a declared domain too, where the box is not the range.
/// </summary>
public sealed class InputStandardizationTests(ITestOutputHelper output)
{
    // x on [−2, 5] skewed toward −2, y on [10, 13] inside a declared domain [9, 14], z a categorical predictor.
    private static BassData Data()
    {
        const int n = 50;
        var x = new double[n];
        var y = new double[n];
        var z = new int[n];
        var r = new double[n];
        for (int i = 0; i < n; i++)
        {
            double u = i * 0.4142135623730950 % 1.0, v = i * 0.7320508075688772 % 1.0;
            x[i] = -2.0 + 7.0 * u * u;
            y[i] = 10.0 + 3.0 * v;
            z[i] = i % 2;
            r[i] = x[i] + y[i];
        }
        return new BassData([new ContinuousPredictor("x", x), new CategoricalPredictor("z", z, ["a", "b"]), new ContinuousPredictor("y", y, new InputDomain(9.0, 14.0))],
                            new ScalarResponse(r));
    }

    private static double Compensated(double[] values)
    {
        double sum = 0.0, carry = 0.0;
        foreach (double v in values)
        {
            double next = sum + v;
            carry += Math.Abs(sum) >= Math.Abs(v) ? sum - next + v : v - next + sum;
            sum = next;
        }
        return sum + carry;
    }

    [Fact]
    public void TheTransform_KeepsTheTrainingMeanAndSampleDeviation_BesideTheBox()
    {
        BassData data = Data();
        BassInputTransform transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        double worst = 0.0;
        foreach ((int predictor, int variable) in new[] { (0, 0), (2, 1) })
        {
            double[] raw = ((ContinuousPredictor)data.Predictors[predictor]).Values.ToArray();
            double mean = Compensated(raw) / raw.Length;
            double deviation = Math.Sqrt(Compensated([.. raw.Select(v => (v - mean) * (v - mean))]) / (raw.Length - 1));
            var scale = (RangeScale)transform.Predictors[predictor];
            Standardization s = scale.Standardization!;
            Assert.Equal(mean, s.Mean, 14);
            Assert.True(Math.Abs(s.Deviation / deviation - 1.0) <= 1e-14, $"{scale.Name}: {s.Deviation:R} against {deviation:R}");
            Assert.Equal((raw[3] - s.Mean) / s.Deviation, s.Standardize(raw[3]));
            (double offset, double slope) = inputs.Standardizing(variable)!.Value;
            Assert.Equal((scale.Max - scale.Min) / s.Deviation, slope);
            ReadOnlySpan<double> u = inputs.Cases.Continuous(variable);
            for (int i = 0; i < u.Length; i++) worst = Math.Max(worst, Math.Abs(offset + slope * u[i] - (raw[i] - mean) / deviation));
        }
        output.WriteLine($"worst |a + b·u − (x − mean)/sd| over the training cases: {worst:E2}");
        Assert.True(worst <= 1e-14, $"{worst:E2}");

        // The declared domain moves the box and so the map's offset and slope, not the standardization.
        var y = (RangeScale)transform.Predictors[2];
        Assert.Equal(9.0, y.Min);
        Assert.Equal((9.0 - y.Standardization!.Mean) / y.Standardization.Deviation, inputs.Standardizing(1)!.Value.Offset, 14);
    }

    [Fact]
    public void AScaleBuiltWithoutData_CarriesNone_AndAStandardizationRefusesWhatIsNotOne()
    {
        Assert.Null(new RangeScale("x", 0.0, 1.0).Standardization);
        Assert.Equal(new Standardization(1.0, 2.0), new RangeScale("x", 0.0, 1.0, null, new Standardization(1.0, 2.0)).Standardization);
        Assert.Throws<ArgumentOutOfRangeException>(() => new Standardization(0.0, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Standardization(0.0, -1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Standardization(double.NaN, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Standardization(0.0, double.PositiveInfinity));
    }
}
