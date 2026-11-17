using System.Text.Json;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;
using LatentRecoveryProbe;

namespace GenerativePrior;

internal sealed record AudioSettings(int MasterSeed, int Repetitions, int Secants, int[] RealRows,
    int[] Windows, int SingleWindow, int[] KeptBins, double Rate, double LowHz, double HighHz,
    double Floor, double RankCutoff, double WhiteningTolerance, double ScaleTolerance, double DecisionEpsilon);

internal static class AudioGeometry
{
    private sealed record Secant(double[] A, double[] B, double[] Unit, double Norm);
    private sealed record Whitening(double[][] Rows, double[][] ProjectorBasis, int Rank, double Condition,
        double OrthogonalityError, double EigenResidual);
    private sealed record Reading(string Manifold, int Repetition, int Seed, string Map, int Rows, int Rank,
        int LinearAtomRank, bool RankIsJacobian,
        double Scale, double Epsilon, double Condition, double WhiteningError, double EigenResidual,
        double ScaleCheckError, double[] MeasuredInnerProducts, double[] OriginalInnerProducts);

    private static double Dot(double[] a, double[] b) => PairwiseSum.Dot(a, b);
    private static double Norm(double[] x) => Math.Sqrt(Dot(x, x));
    private static double[] Apply(double[][] rows, double[] x) => rows.Select(row => Dot(row, x)).ToArray();

    private static double[] Frequencies(int count, AudioSettings s) => Enumerable.Range(0, count)
        .Select(i => count == 1 ? Math.Sqrt(s.LowHz * s.HighHz) :
            s.LowHz * Math.Pow(s.HighHz / s.LowHz, i / (count - 1.0))).ToArray();

    private static double[][] Atoms(int length, int window, double[] frequencies, double rate)
    {
        double[] taper = Windows.Hann(window);
        int start = (length - window) / 2;
        var rows = new List<double[]>();
        foreach (double frequency in frequencies)
        {
            var real = new double[length];
            var imaginary = new double[length];
            for (int j = 0; j < window; j++)
            {
                double phase = 2 * Math.PI * frequency * (start + j) / rate;
                real[start + j] = taper[j] * Math.Cos(phase);
                imaginary[start + j] = -taper[j] * Math.Sin(phase);
            }
            rows.Add(real);
            rows.Add(imaginary);
        }
        return rows.ToArray();
    }

    private static double[][] Combine(double[,] coefficients, double[][] rows)
    {
        var result = new double[coefficients.GetLength(0)][];
        for (int i = 0; i < result.Length; i++)
        {
            result[i] = new double[rows[0].Length];
            for (int j = 0; j < rows.Length; j++)
                for (int t = 0; t < result[i].Length; t++) result[i][t] += coefficients[i, j] * rows[j][t];
        }
        return result;
    }

    private static Whitening Whiten(double[][] rows, AudioSettings s)
    {
        int m = rows.Length;
        var gram = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j <= i; j++) gram[j, i] = gram[i, j] = Dot(rows[i], rows[j]);
        EigenResult eigen = DenseEigen.DecomposeSymmetric(gram);
        if (!eigen.Convergence.Converged) throw new ArithmeticException("Row-Gram eigensolver did not converge.");
        double maximum = eigen.Eigenvalues.Max();
        int[] kept = Enumerable.Range(0, m).Where(i => eigen.Eigenvalues[i] > s.RankCutoff * maximum).ToArray();
        double[][] white = Combine(MatrixFunctions.FactorRows(eigen, kept, value => 1 / Math.Sqrt(value)), rows);
        double[][] rawBasis = Combine(MatrixFunctions.FactorRows(eigen, kept, _ => 1), rows);
        var basis = new List<double[]>();
        foreach (double[] raw in rawBasis)
        {
            double[] v = (double[])raw.Clone();
            for (int pass = 0; pass < 2; pass++)
                foreach (double[] q in basis)
                {
                    double projection = Dot(v, q);
                    for (int j = 0; j < v.Length; j++) v[j] -= projection * q[j];
                }
            double norm = Norm(v);
            if (!(norm > 0)) throw new ArithmeticException("Retained row-space basis lost rank.");
            for (int j = 0; j < v.Length; j++) v[j] /= norm;
            basis.Add(v);
        }
        double error = 0, residual = 0;
        for (int i = 0; i < white.Length; i++)
            for (int j = 0; j <= i; j++) error = Math.Max(error, Math.Abs(Dot(white[i], white[j]) - (i == j ? 1 : 0)));
        foreach (int v in kept)
            for (int i = 0; i < m; i++)
            {
                double product = 0;
                for (int j = 0; j < m; j++) product += gram[i, j] * eigen.Eigenvectors[v][j];
                residual = Math.Max(residual, Math.Abs(product - eigen.Eigenvalues[v] * eigen.Eigenvectors[v][i]));
            }
        double condition = kept.Length == 0 ? 0 : maximum / kept.Min(i => eigen.Eigenvalues[i]);
        return new(white, basis.ToArray(), kept.Length, condition, error, residual);
    }

    private static double Objective(double[] a, double[] c, double scale)
    {
        double maximum = 0;
        for (int i = 0; i < a.Length; i++) maximum = Math.Max(maximum, Math.Abs(scale * a[i] - c[i]));
        return maximum;
    }

    private static bool Interval(double[] a, double[] c, double error, out double lower, out double upper)
    {
        lower = 0; upper = double.PositiveInfinity;
        for (int i = 0; i < a.Length; i++)
        {
            if (a[i] == 0)
            {
                if (Math.Abs(c[i]) > error) return false;
                continue;
            }
            double lo = (c[i] - error) / a[i], hi = (c[i] + error) / a[i];
            lower = Math.Max(lower, Math.Min(lo, hi));
            upper = Math.Min(upper, Math.Max(lo, hi));
            if (upper < lower) return false;
        }
        return true;
    }

    private static (double Scale, double Error, double CheckError) Fit(double[] a, double[] c, bool check)
    {
        double lo = 0, hi = c.Max(Math.Abs);
        for (int i = 0; i < 80; i++)
        {
            double mid = 0.5 * (lo + hi);
            if (Interval(a, c, mid, out _, out _)) hi = mid;
            else lo = mid;
        }
        if (!Interval(a, c, hi, out double lower, out double upper))
            throw new ArithmeticException("Scale interval was lost.");
        double scale = double.IsFinite(upper) ? 0.5 * (lower + upper) : Math.Max(1, lower);
        double error = Objective(a, c, scale), reference = error;
        if (check)
        {
            reference = Objective(a, c, 0);
            var lines = a.Zip(c, (x, y) => (Slope: x, Intercept: -y)).SelectMany(
                line => new[] { line, (Slope: -line.Slope, Intercept: -line.Intercept) }).ToArray();
            for (int i = 0; i < lines.Length; i++)
                for (int j = 0; j < i; j++)
                {
                    double denominator = lines[i].Slope - lines[j].Slope;
                    if (denominator == 0) continue;
                    double candidate = (lines[j].Intercept - lines[i].Intercept) / denominator;
                    if (candidate >= 0 && double.IsFinite(candidate))
                        reference = Math.Min(reference, Objective(a, c, candidate));
                }
        }
        return (scale, error, Math.Abs(error - reference));
    }

    private static double[] LogPower(double[][] rows, double[] signal, double floor)
    {
        double[] linear = Apply(rows, signal);
        return Enumerable.Range(0, linear.Length / 2)
            .Select(j => Math.Log(linear[2 * j] * linear[2 * j] + linear[2 * j + 1] * linear[2 * j + 1] + floor)).ToArray();
    }

    private static double[][] LogJacobian(double[][] rows, double[] signal, double floor)
    {
        double[] linear = Apply(rows, signal);
        var result = new double[rows.Length / 2][];
        for (int i = 0; i < result.Length; i++)
        {
            result[i] = new double[signal.Length];
            double denominator = linear[2 * i] * linear[2 * i] + linear[2 * i + 1] * linear[2 * i + 1] + floor;
            for (int j = 0; j < signal.Length; j++)
                result[i][j] = 2 * (linear[2 * i] * rows[2 * i][j] +
                    linear[2 * i + 1] * rows[2 * i + 1][j]) / denominator;
        }
        return result;
    }

    public static int Run(AudioSettings s, string output)
    {
        int generatorSeed = SeedTree.Derive(1600914390, 30)[0];
        var generatorSettings = new Settings(80, 250, 600, 150, 30, 1600914390, 1,
            Math.ScaleB(1, -52), 50000, 0.001, 1e-10);
        var generator = new Network(generatorSettings, SeedTree.Derive(generatorSeed, 5));
        var synth = new ResonanceSynthesizer(3);
        var readings = new List<Reading>();
        using var secantFile = new StreamWriter(Path.Combine(output, "secants.jsonl"));
        int[] seeds = SeedTree.Derive(s.MasterSeed, s.Repetitions);
        for (int rep = 0; rep < seeds.Length; rep++)
            foreach (string manifold in new[] { "resonance", "g1" })
            {
                int[] branches = SeedTree.Derive(seeds[rep], 4);
                var inputRng = new Xoshiro256PlusPlus(branches[manifold == "g1" ? 0 : 1]);
                var secants = new List<Secant>();
                for (int i = 0; i < s.Secants; i++)
                {
                    double[] MakeCode() => manifold == "resonance" ? synth.Draw(inputRng) :
                        Enumerable.Range(0, 80).Select(_ => Normal.SampleStandard(inputRng)).ToArray();
                    double[] Signal(double[] code)
                    {
                        if (manifold == "resonance") return synth.Waveform(code);
                        var state = new Forward(generatorSettings);
                        generator.Predict(code, state);
                        return state.Output;
                    }
                    double[] ca = MakeCode(), cb = MakeCode(), a = Signal(ca), b = Signal(cb);
                    double[] delta = a.Zip(b, (x, y) => x - y).ToArray();
                    double norm = Norm(delta);
                    if (!(norm > 0)) throw new ArithmeticException("Zero secant.");
                    secants.Add(new(a, b, delta.Select(v => v / norm).ToArray(), norm));
                    secantFile.WriteLine(JsonSerializer.Serialize(new { Manifold = manifold, Repetition = rep + 1,
                        Seed = seeds[rep], Secant = i, CodeA = ca, CodeB = cb, Norm = norm }));
                }
                int length = secants[0].A.Length;
                var pairs = Enumerable.Range(0, s.Secants).Select(i => (i, i)).Concat(
                    Enumerable.Range(0, s.Secants / 2).Select(i => (2 * i, 2 * i + 1))).ToArray();
                double[] original = pairs.Select(pair => pair.Item1 == pair.Item2 ? 1.0 :
                    Math.Clamp(Dot(secants[pair.Item1].Unit, secants[pair.Item2].Unit), -1, 1)).ToArray();
                var whiteningCache = new Dictionary<double[][], Whitening>();
                Whitening GetWhitening(double[][] map)
                {
                    if (!whiteningCache.TryGetValue(map, out Whitening? result))
                    {
                        result = Whiten(map, s);
                        whiteningCache.Add(map, result);
                    }
                    return result;
                }
                void Measure(string name, double[][] map, bool whiten, bool nonlinear = false)
                {
                    Whitening metadata = GetWhitening(map);
                    Whitening? white = whiten ? metadata : null;
                    double[][] rows = white?.Rows ?? map;
                    double[][] transformed = secants.Select(secant => nonlinear ?
                        LogPower(map, secant.A, s.Floor).Zip(LogPower(map, secant.B, s.Floor),
                            (x, y) => (x - y) / secant.Norm).ToArray() : Apply(rows, secant.Unit)).ToArray();
                    double[] inner = pairs.Select(pair => Dot(transformed[pair.Item1], transformed[pair.Item2])).ToArray();
                    double whiteningError = white?.OrthogonalityError ?? 0;
                    if (white is not null)
                    {
                        double[][] projected = secants.Select(secant => Apply(white.ProjectorBasis, secant.Unit)).ToArray();
                        for (int i = 0; i < pairs.Length; i++)
                            whiteningError = Math.Max(whiteningError, Math.Abs(inner[i] -
                                Dot(projected[pairs[i].Item1], projected[pairs[i].Item2])));
                    }
                    var fit = Fit(inner, original, rep == 0);
                    if (whiteningError > s.WhiteningTolerance || fit.CheckError > s.ScaleTolerance)
                        throw new ArithmeticException("Whitening/projector or scale control failed.");
                    Whitening spectral = nonlinear ? Whiten(LogJacobian(map, secants[0].A, s.Floor), s) : metadata;
                    int rank = spectral.Rank;
                    readings.Add(new(manifold, rep + 1, seeds[rep], name, nonlinear ? map.Length / 2 : map.Length,
                        rank, metadata.Rank, nonlinear, fit.Scale, fit.Error, spectral.Condition, whiteningError,
                        spectral.EigenResidual, fit.CheckError, inner, original));
                }
                foreach (int m in s.RealRows)
                {
                    Measure("gaussian", new Matrix(m, length, branches[2]).Rows, false);
                    double[][] one = Atoms(length, s.SingleWindow, Frequencies(m / 2, s), s.Rate);
                    Measure("single", one, false);
                    Measure("single-white", one, true);
                    double[][] multi = s.Windows.SelectMany(window =>
                        Atoms(length, window, Frequencies(m / (2 * s.Windows.Length), s), s.Rate)).ToArray();
                    Measure("multi", multi, false);
                    Measure("multi-white", multi, true);
                }
                int[] order = Enumerable.Range(0, synth.L).ToArray();
                var subsetRng = new Xoshiro256PlusPlus(branches[3]);
                for (int i = order.Length - 1; i > 0; i--)
                {
                    int j = subsetRng.NextInt(i + 1);
                    (order[i], order[j]) = (order[j], order[i]);
                }
                foreach (int bins in s.KeptBins)
                {
                    double[] frequencies = order[..bins].Order().Select(i => synth.Frequencies[i]).ToArray();
                    double[][] kept = Atoms(length, length, frequencies, s.Rate);
                    Measure("kept-complex", kept, false);
                    Measure("kept-complex-white", kept, true);
                    Measure("kept-log-power", kept, false, true);
                }
                Console.WriteLine($"G4 {manifold}, repetition {rep + 1}/{s.Repetitions} completed.");
            }
        File.WriteAllLines(Path.Combine(output, "readings.jsonl"), readings.Select(row => JsonSerializer.Serialize(row)));
        var groups = readings.GroupBy(r => new { r.Manifold, r.Map, r.Rows }).Select(g => {
            double[] eps = g.Select(r => r.Epsilon).Order().ToArray();
            return new { g.Key.Manifold, g.Key.Map, g.Key.Rows, RankMinimum = g.Min(r => r.Rank), RankMaximum = g.Max(r => r.Rank),
                MedianEpsilon = eps[eps.Length / 2], MinimumEpsilon = eps[0], MaximumEpsilon = eps[^1],
                MinimumScale = g.Min(r => r.Scale), MaximumScale = g.Max(r => r.Scale),
                SampledNearIsometry = eps[eps.Length / 2] <= s.DecisionEpsilon };
        }).ToArray();
        LabProgram.Write(Path.Combine(output, "summary.json"), new { LabProgram.Revision, LabProgram.SettingsSha256,
            NumericalControlsPassed = true, Groups = groups, MaximumWhiteningError = readings.Max(r => r.WhiteningError),
            MaximumScaleCheckError = readings.Max(r => r.ScaleCheckError) });
        return 0;
    }
}
