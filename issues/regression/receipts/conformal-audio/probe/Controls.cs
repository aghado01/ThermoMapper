using Maths.Regression.Conformal;

namespace ConformalAudio;

internal static class Controls
{
    internal static double[] Dft(double[] signal) => Enumerable.Range(0, signal.Length).Select(k => {
        double re = 0, im = 0;
        for (int n = 0; n < signal.Length; n++) { double phase = 2 * Math.PI * k * n / signal.Length; re += signal[n] * Math.Cos(phase); im -= signal[n] * Math.Sin(phase); }
        return re * re + im * im;
    }).ToArray();
    internal static void Run()
    {
        double alpha = App.D(App.Settings, "Alpha"), tolerance = App.D(App.Settings, "FourierTolerance");
        var cases = new List<object>(); int frameMisses = 0, groupMisses = 0;
        for (int calibration = 0; calibration <= 1; calibration++)
            for (int test = 0; test <= 1; test++)
            {
                SplitCalibration frames = SplitConformal.Calibrate(Enumerable.Repeat(-(double)calibration, 20).ToArray(), alpha: alpha);
                SplitCalibration groups = SplitConformal.Calibrate([-(double)calibration], alpha: alpha);
                bool fa = frames.Accepts(0, -test), ga = groups.Accepts(0, -test);
                if (!fa) frameMisses++; if (!ga) groupMisses++;
                cases.Add(new { Calibration = calibration, Test = test, FrameAccepted = fa, GroupAccepted = ga });
            }
        App.Require(frameMisses == 1 && groupMisses == 0, "Duplicated-frame rank enumeration.");
        int k18 = SplitConformal.Calibrate(new double[18], alpha: alpha).Level(0), k19 = SplitConformal.Calibrate(new double[19], alpha: alpha).Level(0);
        App.Require(k18 == 0 && k19 == 1, "Group-count boundary.");
        double[] plus = Dft([1, 0, 1]), minus = Dft([1, 0, -1]);
        App.Require(Math.Abs(plus[0] - 4) < tolerance && minus[0] == 0 && Math.Abs(plus[1] - 1) < tolerance && Math.Abs(minus[1] - 3) < tolerance, "DFT witness.");
        foreach (int start in new[] { 0, 1 })
        {
            double[] a = Dft(new[] { 1.0, 0, 1 }[start..(start + 2)]), b = Dft(new[] { 1.0, 0, -1 }[start..(start + 2)]);
            App.Require(a.Zip(b, (x, y) => Math.Abs(x - y)).Max() < tolerance, "Adjacent short views.");
        }
        double coverage = 0, size = 0;
        foreach (int sign in new[] { -1, 1 }) foreach (bool both in new[] { false, true })
        {
            double probability = 0.5 * (both ? 0.9 : 0.1); int target = sign == 1 ? 4 : 0;
            coverage += probability * (both || target == 0 ? 1 : 0); size += probability * (both ? 2 : 1);
        }
        App.Require(Math.Abs(coverage - 0.95) < 1e-12 && Math.Abs(size - 1.9) < 1e-12, "Alphabet sharpness.");
        App.Write("controls.json", new { Passed = true, RankCases = cases, FrameMissProbability = frameMisses / 4.0,
            GroupMissProbability = groupMisses / 4.0, LevelAt18 = k18, LevelAt19 = k19, PlusDft = plus, MinusDft = minus,
            RandomizedCoverage = coverage, RandomizedCardinality = size, DistinguishingOracleCoverage = 1, DistinguishingOracleCardinality = 1 });
    }
}
