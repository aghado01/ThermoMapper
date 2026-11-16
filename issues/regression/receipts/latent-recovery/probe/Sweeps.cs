using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace LatentRecoveryProbe;

/// <summary>Energy-normalized analysis windows, Σ g² = 1 (the Lean chip's S2 normalization).</summary>
public static class Windows
{
    public static double[] Rectangular(int w) => Enumerable.Repeat(1.0 / Math.Sqrt(w), w).ToArray();

    /// <summary>The periodic Hann window, 0.5·(1 − cos(2πj/w)), scaled to unit energy.</summary>
    public static double[] Hann(int w)
    {
        var g = new double[w];
        double energy = 0.0;
        for (int j = 0; j < w; j++)
        {
            g[j] = 0.5 * (1.0 - Math.Cos(2.0 * Math.PI * j / w));
            energy += g[j] * g[j];
        }
        double s = 1.0 / Math.Sqrt(energy);
        for (int j = 0; j < w; j++) g[j] *= s;
        return g;
    }

    public static double[] Of(string name, int w) => name switch
    {
        "rect" => Rectangular(w),
        "hann" => Hann(w),
        _ => throw new ArgumentException($"Unknown window '{name}'."),
    };
}

/// <summary>Small statistics helpers shared by the sweeps and the pilot.</summary>
public static class Stats
{
    public static double Mean(IReadOnlyList<double> v) => v.Average();

    public static double Sd(IReadOnlyList<double> v)
    {
        double m = v.Average();
        return Math.Sqrt(v.Sum(x => (x - m) * (x - m)) / (v.Count - 1));
    }

    /// <summary>Least-squares slope of y on x.</summary>
    public static double Slope(IReadOnlyList<double> x, IReadOnlyList<double> y)
    {
        double mx = x.Average(), my = y.Average(), sxy = 0, sxx = 0;
        for (int i = 0; i < x.Count; i++)
        {
            sxy += (x[i] - mx) * (y[i] - my);
            sxx += (x[i] - mx) * (x[i] - mx);
        }
        return sxy / sxx;
    }

    public static string F(double v, string format = "G6") => v.ToString(format, CultureInfo.InvariantCulture);
}

/// <summary>
/// B1: the covariance of log powers. Part one draws jointly circular complex Gaussian pairs (a, b), E|a|² = E|b|² = 1 and
/// E[a·conj(b)] = √κ·e^{iφ}, and compares Cov(ln|a|², ln|b|²) with Li₂(κ). Part two computes short-time Fourier cells of
/// white noise and compares the correlation of their log powers with Li₂(κ)/Li₂(1), κ the squared modulus of the two
/// cells' atoms' Gram entry.
/// </summary>
public static class LogPowerCovarianceSweep
{
    public sealed record Settings(int MasterSeed, double[] Kappas, double Phase, int PairsPerKappa, int PairBatches,
                                  int[] WindowLengths, int SignalLength, int Centre, int Replicates, int CellBatches, int Threads);

    public static Settings Default => new(
        MasterSeed: 1116001,
        Kappas: [0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0],
        Phase: 1.0,
        PairsPerKappa: 4_000_000,
        PairBatches: 200,
        WindowLengths: [64, 128, 256],
        SignalLength: 1024,
        Centre: 512,
        Replicates: 200_000,
        CellBatches: 100,
        Threads: 6);

    // Seed layout from the master: [0, K) the pairs' main seeds, [K, 2K) their rerun seeds, 2K the complex noise, 2K+1
    // the real noise, 2K+2 and 2K+3 their rerun seeds.
    public static int[] Seeds(Settings s) => SeedTree.Derive(s.MasterSeed, 2 * s.Kappas.Length + 4);

    public static (double Estimate, double Se) PairCovariance(double kappa, double phase, int seed, int n, int batches)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double rho = Math.Sqrt(kappa), cr = rho * Math.Cos(phase), ci = rho * Math.Sin(phase), s = Math.Sqrt(1.0 - kappa);
        double h = Math.Sqrt(0.5);
        int per = n / batches;
        var covs = new double[batches];
        for (int b = 0; b < batches; b++)
        {
            double sa = 0, sb = 0, sab = 0;
            for (int i = 0; i < per; i++)
            {
                double x1 = h * Normal.SampleStandard(rng), y1 = h * Normal.SampleStandard(rng);
                double x2 = h * Normal.SampleStandard(rng), y2 = h * Normal.SampleStandard(rng);
                double br = cr * x1 - ci * y1 + s * x2, bi = cr * y1 + ci * x1 + s * y2;
                double la = Math.Log(x1 * x1 + y1 * y1), lb = Math.Log(br * br + bi * bi);
                sa += la;
                sb += lb;
                sab += la * lb;
            }
            covs[b] = (sab - sa * sb / per) / (per - 1);
        }
        return (Stats.Mean(covs), Stats.Sd(covs) / Math.Sqrt(batches));
    }

    public static void RunPairs(Settings s, string outDir, int[]? only = null, bool rerun = false)
    {
        int[] seeds = Seeds(s);
        int k = s.Kappas.Length;
        var rows = new string[k];
        var indices = only ?? Enumerable.Range(0, k).ToArray();
        Parallel.ForEach(indices, new ParallelOptions { MaxDegreeOfParallelism = s.Threads }, i =>
        {
            int seed = rerun ? seeds[k + i] : seeds[i];
            (double est, double se) = PairCovariance(s.Kappas[i], s.Phase, seed, s.PairsPerKappa, s.PairBatches);
            double pred = Dilogarithm.Li2(s.Kappas[i]);
            double linear = s.Kappas[i] * Dilogarithm.ZetaTwo;
            rows[i] = string.Join(",", Stats.F(s.Kappas[i]), seed, Stats.F(est, "F6"), Stats.F(se, "E3"), Stats.F(pred, "F6"),
                                  Stats.F((est - pred) / se, "F3"), Stats.F(linear, "F6"), Stats.F((est - linear) / se, "F1"));
        });
        string file = Path.Combine(outDir, rerun ? "b1-pairs-rerun.csv" : "b1-pairs.csv");
        File.WriteAllLines(file, ["kappa,seed,estimate,se,li2,z,linear_kappa_zeta2,z_linear", .. indices.Select(i => rows[i])]);
        Console.WriteLine(File.ReadAllText(file));
    }

    /// <summary>A short-time Fourier cell: the window of length W centred at sample Centre, at Frequency cycles/sample.</summary>
    public sealed record Cell(int Centre, double Frequency, int W)
    {
        public int Start => Centre - W / 2;

        /// <summary>conj(a_n) over the support, a_n = g(n − Start)·e^{2πi f (n − Centre)}: X = Σ x_n·conj(a_n).</summary>
        public Complex[] Analysis(double[] g)
        {
            var c = new Complex[W];
            for (int j = 0; j < W; j++)
            {
                double phase = -2.0 * Math.PI * Frequency * (Start + j - Centre);
                c[j] = g[j] * new Complex(Math.Cos(phase), Math.Sin(phase));
            }
            return c;
        }
    }

    public sealed record Pair(string Noise, string Family, Cell Reference, Cell Other);

    public static List<Pair> Pairs(Settings s)
    {
        var pairs = new List<Pair>();
        int c0 = s.Centre;
        foreach (string noise in new[] { "complex", "real" })
        {
            foreach (int w in s.WindowLengths)
            {
                var references = new List<(double F, string Family)> { (0.25, "mid") };
                if (noise == "real") references.Add((1.0 / w, "near-dc"));
                foreach ((double f0, string family) in references)
                {
                    var r = new Cell(c0, f0, w);
                    int[] dts = [0, w / 8, w / 4, 3 * w / 8, w / 2, 3 * w / 4, w];
                    double[] dfs = family == "near-dc" ? [-1.0 / w, 0, 0.5 / w, 1.0 / w, 2.0 / w] : [0, 0.5 / w, 1.0 / w, 2.0 / w];
                    foreach (int dt in dts)
                        foreach (double df in dfs)
                        {
                            if (dt == 0 && df == 0) continue;
                            pairs.Add(new Pair(noise, family, r, new Cell(c0 + dt, f0 + df, w)));
                        }
                }
            }
            (int, int)[] cross = [(64, 128), (128, 256), (64, 256)];
            foreach ((int w1, int w2) in cross)
            {
                var r = new Cell(c0, 0.25, w1);
                foreach (int dt in new[] { 0, w1 / 2 })
                    foreach (double df in new[] { 0.0, 1.0 / w2 })
                        pairs.Add(new Pair(noise, "cross-w", r, new Cell(c0 + dt, 0.25 + df, w2)));
            }
        }
        return pairs;
    }

    /// <summary>The Gram entry Σ conj(a1)·a2 and the pseudo entry Σ conj(a1)·conj(a2) of two cells.</summary>
    public static (Complex Gram, Complex Pseudo) Entries(Cell c1, Cell c2)
    {
        Complex[] a1 = c1.Analysis(Windows.Hann(c1.W)), a2 = c2.Analysis(Windows.Hann(c2.W));
        Complex gram = 0, pseudo = 0;
        for (int n = Math.Max(c1.Start, c2.Start); n < Math.Min(c1.Start + c1.W, c2.Start + c2.W); n++)
        {
            Complex x1 = a1[n - c1.Start], x2 = a2[n - c2.Start];   // x = conj(a)
            gram += x1 * Complex.Conjugate(x2);                      // conj(a1)·a2
            pseudo += x1 * x2;                                       // conj(a1)·conj(a2)
        }
        return (gram, pseudo);
    }

    public static void RunCells(Settings s, string outDir, string noise, bool rerun = false)
    {
        int[] seeds = Seeds(s);
        int k = s.Kappas.Length;
        int seed = noise == "complex" ? seeds[2 * k + (rerun ? 2 : 0)] : seeds[2 * k + 1 + (rerun ? 2 : 0)];
        List<Pair> pairs = Pairs(s).Where(p => p.Noise == noise).ToList();
        var cells = pairs.SelectMany(p => new[] { p.Reference, p.Other }).Distinct().ToList();
        var index = cells.Select((c, i) => (c, i)).ToDictionary(t => t.c, t => t.i);
        var analysis = cells.Select(c => c.Analysis(Windows.Hann(c.W))).ToArray();
        int np = pairs.Count;
        int[] ri = pairs.Select(p => index[p.Reference]).ToArray(), oi = pairs.Select(p => index[p.Other]).ToArray();
        int batches = s.CellBatches, per = s.Replicates / batches;
        int[] batchSeeds = SeedTree.Derive(seed, batches);
        var sums = new double[batches, np, 5];
        Parallel.For(0, batches, new ParallelOptions { MaxDegreeOfParallelism = s.Threads }, b =>
        {
            var rng = new Xoshiro256PlusPlus(batchSeeds[b]);
            var re = new double[s.SignalLength];
            var im = new double[s.SignalLength];
            var logs = new double[cells.Count];
            double h = noise == "complex" ? Math.Sqrt(0.5) : 1.0;
            for (int r = 0; r < per; r++)
            {
                for (int n = 0; n < s.SignalLength; n++)
                {
                    re[n] = h * Normal.SampleStandard(rng);
                    im[n] = noise == "complex" ? h * Normal.SampleStandard(rng) : 0.0;
                }
                for (int c = 0; c < cells.Count; c++)
                {
                    Complex[] a = analysis[c];
                    int start = cells[c].Start;
                    double xr = 0, xi = 0;
                    for (int j = 0; j < a.Length; j++)
                    {
                        double ar = a[j].Real, ai = a[j].Imaginary, vr = re[start + j], vi = im[start + j];
                        xr += vr * ar - vi * ai;
                        xi += vr * ai + vi * ar;
                    }
                    logs[c] = Math.Log(xr * xr + xi * xi);
                }
                for (int p = 0; p < np; p++)
                {
                    double l1 = logs[ri[p]], l2 = logs[oi[p]];
                    sums[b, p, 0] += l1;
                    sums[b, p, 1] += l2;
                    sums[b, p, 2] += l1 * l1;
                    sums[b, p, 3] += l2 * l2;
                    sums[b, p, 4] += l1 * l2;
                }
            }
        });
        static double Corr(double n, double s1, double s2, double s11, double s22, double s12)
        {
            double c12 = s12 - s1 * s2 / n, c11 = s11 - s1 * s1 / n, c22 = s22 - s2 * s2 / n;
            return c12 / Math.Sqrt(c11 * c22);
        }
        var lines = new List<string> { "noise,family,w_ref,f_ref,dt,df,w_other,kappa,pseudo_ref_other,pseudo_ref_self,predicted,estimate,se,z,seed" };
        for (int p = 0; p < np; p++)
        {
            var batchCorr = new double[batches];
            double t1 = 0, t2 = 0, t11 = 0, t22 = 0, t12 = 0;
            for (int b = 0; b < batches; b++)
            {
                batchCorr[b] = Corr(per, sums[b, p, 0], sums[b, p, 1], sums[b, p, 2], sums[b, p, 3], sums[b, p, 4]);
                t1 += sums[b, p, 0]; t2 += sums[b, p, 1]; t11 += sums[b, p, 2]; t22 += sums[b, p, 3]; t12 += sums[b, p, 4];
            }
            double est = Corr((double)per * batches, t1, t2, t11, t22, t12);
            double se = Stats.Sd(batchCorr) / Math.Sqrt(batches);
            Pair pr = pairs[p];
            (Complex gram, Complex pseudo) = Entries(pr.Reference, pr.Other);
            (_, Complex self) = Entries(pr.Reference, pr.Reference);
            double kappa = Math.Min(1.0, gram.Magnitude * gram.Magnitude);
            double pred = Dilogarithm.Li2(kappa) / Dilogarithm.ZetaTwo;
            lines.Add(string.Join(",", pr.Noise, pr.Family, pr.Reference.W, Stats.F(pr.Reference.Frequency, "F6"), pr.Other.Centre - pr.Reference.Centre,
                                  Stats.F(pr.Other.Frequency - pr.Reference.Frequency, "F6"), pr.Other.W, Stats.F(kappa, "F6"), Stats.F(pseudo.Magnitude, "F6"),
                                  Stats.F(self.Magnitude, "F6"), Stats.F(pred, "F6"), Stats.F(est, "F6"), Stats.F(se, "E3"), Stats.F((est - pred) / se, "F3"), seed));
        }
        string file = Path.Combine(outDir, $"b1-cells-{noise}{(rerun ? "-rerun" : "")}.csv");
        File.WriteAllLines(file, lines);
        Console.WriteLine($"{file}: {np} pairs over {cells.Count} cells, {per * batches} replicates");
    }
}

/// <summary>
/// B2: peak power against window length for a tone, an impulse, white noise, a decaying resonance and a linear chirp, under
/// energy-normalized windows (S2: Σ g² = 1, power |Σ g_j·x_j·e^{−2πi f j}|²).
/// </summary>
public static class WindowScalingSweep
{
    public sealed record Settings(int MasterSeed, int[] WindowLengths, string[] WindowNames, double ToneFrequency, int NoiseFrames,
                                  double[] DecayTimes, double[] ChirpRates, int Threads);

    // Even lengths on a quarter-octave grid from 4 to 4096.
    public static int[] Grid() => Enumerable.Range(0, 41).Select(j => 2 * (int)Math.Round(2.0 * Math.Pow(2.0, j / 4.0))).Distinct().ToArray();

    public static Settings Default => new(
        MasterSeed: 1116002,
        WindowLengths: Grid(),
        WindowNames: ["rect", "hann"],
        ToneFrequency: 0.1234567,
        NoiseFrames: 4000,
        DecayTimes: [64, 128, 256],
        ChirpRates: [1e-4, 2.5e-4, 6.25e-4],
        Threads: 6);

    private static int NextPow2(int n)
    {
        int p = 1;
        while (p < n) p <<= 1;
        return p;
    }

    /// <summary>The largest |FFT|² of the windowed segment zero-padded to at least pad·w points.</summary>
    public static double GridPeak(Complex[] windowed, int pad)
    {
        int n = NextPow2(pad * windowed.Length);
        var buf = new Complex[n];
        Array.Copy(windowed, buf, windowed.Length);
        Fft.Forward(buf);
        double best = 0;
        foreach (Complex c in buf) best = Math.Max(best, c.Real * c.Real + c.Imaginary * c.Imaginary);
        return best;
    }

    public static void Run(Settings s, string outDir, bool rerunNoise = false, int[]? noiseOnly = null)
    {
        int[] seeds = SeedTree.Derive(s.MasterSeed, 4 * s.WindowLengths.Length);   // [noise complex | noise real | reruns]
        var lines = new List<string> { "signal,window,param,w,power,formula,grid_peak,extra" };
        object gate = new();
        void Add(string line) { lock (gate) lines.Add(line); }

        if (!rerunNoise)
        {
            foreach (string wn in s.WindowNames)
                foreach (int w in s.WindowLengths)
                {
                    double[] g = Windows.Of(wn, w);
                    // Tone, A = 1: the S2 formula |Σ g|², the DFT at f0 summed directly, and the padded grid's peak.
                    var seg = new Complex[w];
                    Complex atF0 = 0;
                    for (int j = 0; j < w; j++)
                    {
                        Complex x = Complex.FromPolarCoordinates(1.0, 2 * Math.PI * s.ToneFrequency * j);
                        seg[j] = g[j] * x;
                        atF0 += seg[j] * Complex.FromPolarCoordinates(1.0, -2 * Math.PI * s.ToneFrequency * j);
                    }
                    double sumG = g.Sum();
                    Add(string.Join(",", "tone", wn, "", w, Stats.F(atF0.Magnitude * atF0.Magnitude, "R"), Stats.F(sumG * sumG, "R"), Stats.F(GridPeak(seg, 16), "R"), ""));
                    // Impulse, a = 1, best placement: max over positions of g², against max g².
                    double best = 0;
                    for (int j = 0; j < w; j++) best = Math.Max(best, g[j] * g[j]);
                    Add(string.Join(",", "impulse", wn, "", w, Stats.F(best, "R"), Stats.F(best, "R"), "", ""));
                }
        }

        // White noise: the mean power over frames and padded bins, against σ² = 1; and the mean of each frame's peak bin.
        var noiseJobs = (from wn in s.WindowNames from kind in new[] { "complex", "real" } from wi in Enumerable.Range(0, s.WindowLengths.Length) select (wn, kind, wi))
            .Where(t => noiseOnly is null || noiseOnly.Contains(t.wi)).ToList();
        Parallel.ForEach(noiseJobs, new ParallelOptions { MaxDegreeOfParallelism = s.Threads }, job =>
        {
            (string wn, string kind, int wi) = job;
            int w = s.WindowLengths[wi];
            int seedIndex = (kind == "complex" ? 0 : 1) * s.WindowLengths.Length + wi + (rerunNoise ? 2 * s.WindowLengths.Length : 0);
            int seed = SeedTree.Derive(seeds[seedIndex], 2)[wn == "rect" ? 0 : 1];
            var rng = new Xoshiro256PlusPlus(seed);
            double[] g = Windows.Of(wn, w);
            int n = NextPow2(w);
            var frameMeans = new double[s.NoiseFrames];
            double peakSum = 0;
            double h = kind == "complex" ? Math.Sqrt(0.5) : 1.0;
            var buf = new Complex[n];
            for (int f = 0; f < s.NoiseFrames; f++)
            {
                Array.Clear(buf);
                for (int j = 0; j < w; j++)
                    buf[j] = g[j] * new Complex(h * Normal.SampleStandard(rng), kind == "complex" ? h * Normal.SampleStandard(rng) : 0.0);
                Fft.Forward(buf);
                double sum = 0, peak = 0;
                foreach (Complex c in buf)
                {
                    double p = c.Real * c.Real + c.Imaginary * c.Imaginary;
                    sum += p;
                    peak = Math.Max(peak, p);
                }
                frameMeans[f] = sum / n;
                peakSum += peak;
            }
            double mean = Stats.Mean(frameMeans), se = Stats.Sd(frameMeans) / Math.Sqrt(s.NoiseFrames);
            Add(string.Join(",", "noise-" + kind, wn, seed, w, Stats.F(mean, "R"), "1", Stats.F(peakSum / s.NoiseFrames, "R"), Stats.F(se, "R")));
        });

        if (!rerunNoise)
        {
            // Decaying resonance x_n = e^{−n/τ}·e^{2πi f0 n}, n ≥ 0: g·envelope is non-negative, so the peak over f is at f0 and
            // equals (Σ_j g_j·env_{s+j})², maximized here over the window's start s.
            var resonanceJobs = (from wn in s.WindowNames from tau in s.DecayTimes from w in s.WindowLengths select (wn, tau, w)).ToList();
            Parallel.ForEach(resonanceJobs, new ParallelOptions { MaxDegreeOfParallelism = s.Threads }, job =>
            {
                (string wn, double tau, int w) = job;
                double[] g = Windows.Of(wn, w);
                double best = 0;
                int bestStart = 0;
                int last = (int)Math.Ceiling(4 * tau);
                for (int st = -w + 1; st <= last; st++)
                {
                    double acc = 0;
                    for (int j = Math.Max(0, -st); j < w; j++) acc += g[j] * Math.Exp(-(st + j) / tau);
                    if (acc * acc > best) { best = acc * acc; bestStart = st; }
                }
                // Check the peak over f on the padded grid at the best start (f0 = 0.1234567).
                var seg = new Complex[w];
                for (int j = 0; j < w; j++)
                {
                    int nn = bestStart + j;
                    seg[j] = nn < 0 ? 0 : g[j] * Complex.FromPolarCoordinates(Math.Exp(-nn / tau), 2 * Math.PI * s.ToneFrequency * nn);
                }
                double rectFormula = 0;
                if (wn == "rect")
                {
                    double a = Math.Exp(-1 / tau);
                    double sum = (1 - Math.Pow(a, w)) / (1 - a);
                    rectFormula = sum * sum / w;
                }
                Add(string.Join(",", "resonance", wn, Stats.F(tau), w, Stats.F(best, "R"), wn == "rect" ? Stats.F(rectFormula, "R") : "",
                                Stats.F(GridPeak(seg, 16), "R"), bestStart));
            });

            // Linear chirp x_j = e^{iπ c j²}, window centred at j = 0, lengths kept below 0.8/c so the instantaneous frequency
            // stays inside the band; the peak over f on a padded grid.
            var chirpJobs = (from wn in s.WindowNames from c in s.ChirpRates from w in s.WindowLengths where w <= 0.8 / c select (wn, c, w)).ToList();
            Parallel.ForEach(chirpJobs, new ParallelOptions { MaxDegreeOfParallelism = s.Threads }, job =>
            {
                (string wn, double c, int w) = job;
                double[] g = Windows.Of(wn, w);
                var seg = new Complex[w];
                for (int j = 0; j < w; j++)
                {
                    double t = j - w / 2;
                    seg[j] = g[j] * Complex.FromPolarCoordinates(1.0, Math.PI * c * t * t);
                }
                Add(string.Join(",", "chirp", wn, Stats.F(c), w, Stats.F(GridPeak(seg, 16), "R"), "", "", ""));
            });
        }

        string file = Path.Combine(outDir, rerunNoise ? "b2-rerun.csv" : "b2.csv");
        File.WriteAllLines(file, lines);
        Console.WriteLine($"{file}: {lines.Count - 1} rows");
    }
}
