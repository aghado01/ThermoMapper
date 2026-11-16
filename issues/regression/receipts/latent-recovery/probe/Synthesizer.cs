using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;

namespace LatentRecoveryProbe;

/// <summary>
/// The pilot's synthesizer, outside the estimator: a code θ (natural logs of physical quantities, on a box) maps to a
/// waveform on one frame, the frame to a magnitude spectrum at fixed analysis frequencies through a unit-energy periodic
/// Hann window, and the spectrum to the feature y = ln(|X(f)|² + ε). Observations add real white Gaussian noise of
/// variance σ²ₙ per sample to the waveform before the transform, so y's noise is the log of a noisy power, not Gaussian.
/// </summary>
public abstract class Synthesizer
{
    protected Synthesizer(int frame, double sampleRate, double floor, double[] frequencies, string[] names, double[] lower, double[] upper)
    {
        Frame = frame;
        SampleRate = sampleRate;
        Floor = floor;
        Frequencies = frequencies;
        Names = names;
        Lower = lower;
        Upper = upper;
        Window = Windows.Hann(frame);
        (_cos, _sin) = Table(frequencies);
    }

    public int Frame { get; }
    public double SampleRate { get; }
    public double Floor { get; }
    public double[] Frequencies { get; }
    public string[] Names { get; }
    public double[] Lower { get; }
    public double[] Upper { get; }
    public double[] Window { get; }
    public int K => Names.Length;
    public int L => Frequencies.Length;

    private readonly double[,] _cos, _sin;

    /// <summary>g_n·cos(2πfn/fs) and −g_n·sin(2πfn/fs) for each analysis frequency.</summary>
    public (double[,] Cos, double[,] Sin) Table(IReadOnlyList<double> frequencies)
    {
        var c = new double[frequencies.Count, Frame];
        var s = new double[frequencies.Count, Frame];
        for (int l = 0; l < frequencies.Count; l++)
            for (int n = 0; n < Frame; n++)
            {
                double phase = 2.0 * Math.PI * frequencies[l] * n / SampleRate;
                c[l, n] = Window[n] * Math.Cos(phase);
                s[l, n] = -Window[n] * Math.Sin(phase);
            }
        return (c, s);
    }

    /// <summary>The clean waveform of code θ.</summary>
    public abstract double[] Waveform(ReadOnlySpan<double> theta);

    /// <summary>|X(f)|² at the table's frequencies.</summary>
    public static void Power(double[] x, double[,] cos, double[,] sin, Span<double> power)
    {
        int frame = x.Length;
        for (int l = 0; l < power.Length; l++)
        {
            double re = 0, im = 0;
            for (int n = 0; n < frame; n++)
            {
                re += cos[l, n] * x[n];
                im += sin[l, n] * x[n];
            }
            power[l] = re * re + im * im;
        }
    }

    /// <summary>The feature ln(|X|² + ε) of a waveform at every analysis frequency.</summary>
    public double[] Feature(double[] x)
    {
        var p = new double[L];
        Power(x, _cos, _sin, p);
        for (int l = 0; l < L; l++) p[l] = Math.Log(p[l] + Floor);
        return p;
    }

    /// <summary>The noiseless feature of code θ at every analysis frequency.</summary>
    public double[] Clean(ReadOnlySpan<double> theta) => Feature(Waveform(theta));

    /// <summary>The noisy feature: σ²ₙ-variance real white noise added to the waveform.</summary>
    public double[] Noisy(ReadOnlySpan<double> theta, double noiseVariance, Xoshiro256PlusPlus rng)
    {
        double[] x = Waveform(theta);
        double sd = Math.Sqrt(noiseVariance);
        for (int n = 0; n < x.Length; n++) x[n] += sd * Normal.SampleStandard(rng);
        return Feature(x);
    }

    /// <summary>A code drawn uniformly from the box.</summary>
    public double[] Draw(Xoshiro256PlusPlus rng)
    {
        var t = new double[K];
        for (int j = 0; j < K; j++) t[j] = Lower[j] + rng.NextDouble() * (Upper[j] - Lower[j]);
        return t;
    }

    /// <summary>A Latin hypercube of n codes on the box.</summary>
    public double[][] LatinHypercube(int n, Xoshiro256PlusPlus rng)
    {
        var design = new double[n][];
        for (int i = 0; i < n; i++) design[i] = new double[K];
        for (int j = 0; j < K; j++)
        {
            int[] perm = Enumerable.Range(0, n).ToArray();
            for (int i = n - 1; i > 0; i--)
            {
                int r = rng.NextInt(i + 1);
                (perm[i], perm[r]) = (perm[r], perm[i]);
            }
            for (int i = 0; i < n; i++)
                design[i][j] = Lower[j] + (perm[i] + rng.NextDouble()) / n * (Upper[j] - Lower[j]);
        }
        return design;
    }

    /// <summary>The clean signal's energy Σ x².</summary>
    public double Energy(ReadOnlySpan<double> theta) => Waveform(theta).Sum(v => v * v);
}

/// <summary>
/// Stage 1: one or two impulse-excited decaying resonances, x(t) = Σ_r G_r·e^{−πB_r t}·sin(2πF_r t) for t ≥ 0, the onset a
/// quarter of the way into a 512-sample frame at 8 kHz, analysed at 64 log-spaced frequencies from 150 Hz to 3.6 kHz. The
/// code's coordinates, in order, are ln F₁, ln B₁, ln G₁, ln F₂, ln B₂; the first k are free, the rest fixed (B₁ = 120 Hz,
/// G₁ = 1, the second resonance absent below k = 4, B₂ = 160 Hz at k = 4, G₂ = 0.7 at k ≥ 4).
/// </summary>
public sealed class ResonanceSynthesizer : Synthesizer
{
    public const int FrameLength = 512, Onset = 128, Bins = 64;
    public const double Rate = 8000.0, FloorValue = 1e-4, LowHz = 150.0, HighHz = 3600.0;

    private static readonly string[] AllNames = ["lnF1", "lnB1", "lnG1", "lnF2", "lnB2"];
    private static readonly double[] AllLower = [Math.Log(300), Math.Log(60), Math.Log(0.5), Math.Log(1400), Math.Log(80)];
    private static readonly double[] AllUpper = [Math.Log(1200), Math.Log(240), Math.Log(2.0), Math.Log(3000), Math.Log(320)];

    public static double[] AnalysisFrequencies() => [.. Enumerable.Range(0, Bins).Select(l => LowHz * Math.Pow(HighHz / LowHz, l / (Bins - 1.0)))];

    public ResonanceSynthesizer(int k)
        : base(FrameLength, Rate, FloorValue, AnalysisFrequencies(), AllNames[..k], AllLower[..k], AllUpper[..k])
    {
        if (k is < 1 or > 5) throw new ArgumentOutOfRangeException(nameof(k));
    }

    public override double[] Waveform(ReadOnlySpan<double> theta)
    {
        int k = K;
        double f1 = Math.Exp(theta[0]);
        double b1 = k >= 2 ? Math.Exp(theta[1]) : 120.0;
        double g1 = k >= 3 ? Math.Exp(theta[2]) : 1.0;
        bool second = k >= 4;
        double f2 = second ? Math.Exp(theta[3]) : 0.0;
        double b2 = k >= 5 ? Math.Exp(theta[4]) : 160.0;
        const double g2 = 0.7;
        var x = new double[Frame];
        for (int n = Onset; n < Frame; n++)
        {
            double t = (n - Onset) / SampleRate;
            double v = g1 * Math.Exp(-Math.PI * b1 * t) * Math.Sin(2 * Math.PI * f1 * t);
            if (second) v += g2 * Math.Exp(-Math.PI * b2 * t) * Math.Sin(2 * Math.PI * f2 * t);
            x[n] = v;
        }
        return x;
    }
}

/// <summary>
/// Stage 2: a harmonic source of fundamental f₀ through one resonance, steady over a 1024-sample frame at 8 kHz:
/// x_n = Σ_h A_h·sin(2π h f₀ n/fs), A_h = R(h f₀), R(f) = 1/√(1 + (2(f − F)/B)²), harmonics up to 3.8 kHz, analysed at 128
/// linearly spaced frequencies from 100 Hz to 2 kHz. The code's coordinates are ln f₀, ln F, ln B; the first k are free
/// (B = 200 Hz when fixed).
/// </summary>
public sealed class HarmonicSynthesizer : Synthesizer
{
    public const int FrameLength = 1024, Bins = 128;
    public const double Rate = 8000.0, FloorValue = 1e-2, LowHz = 100.0, HighHz = 2000.0, TopHarmonicHz = 3800.0;

    private static readonly string[] AllNames = ["lnf0", "lnF", "lnB"];
    private static readonly double[] AllLower = [Math.Log(100), Math.Log(400), Math.Log(100)];
    private static readonly double[] AllUpper = [Math.Log(400), Math.Log(1600), Math.Log(400)];

    public static double[] AnalysisFrequencies() => [.. Enumerable.Range(0, Bins).Select(l => LowHz + (HighHz - LowHz) * l / (Bins - 1.0))];

    public HarmonicSynthesizer(int k)
        : base(FrameLength, Rate, FloorValue, AnalysisFrequencies(), AllNames[..k], AllLower[..k], AllUpper[..k])
    {
        if (k is < 2 or > 3) throw new ArgumentOutOfRangeException(nameof(k));
    }

    public override double[] Waveform(ReadOnlySpan<double> theta)
    {
        double f0 = Math.Exp(theta[0]), centre = Math.Exp(theta[1]), width = K >= 3 ? Math.Exp(theta[2]) : 200.0;
        int top = (int)Math.Floor(TopHarmonicHz / f0);
        var x = new double[Frame];
        for (int h = 1; h <= top; h++)
        {
            double f = h * f0, u = 2 * (f - centre) / width, a = 1.0 / Math.Sqrt(1 + u * u);
            // sin(w(n+1)) = 2cos(w)·sin(wn) − sin(w(n−1)), from sin 0 and sin w.
            double w = 2 * Math.PI * f / SampleRate, c2 = 2 * Math.Cos(w), previous = 0.0, current = Math.Sin(w);
            for (int n = 1; n < Frame; n++)
            {
                x[n] += a * current;
                double next = c2 * current - previous;
                previous = current;
                current = next;
            }
        }
        return x;
    }
}
