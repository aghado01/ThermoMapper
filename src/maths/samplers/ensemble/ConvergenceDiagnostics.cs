using System;
using System.Globalization;
using System.Numerics;
using Maths.Distributions;
using Maths.LinAlg;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// One functional's effective sample size with its <see cref="DiagnosticState"/>: the estimate when
/// <see cref="DiagnosticState.Estimated"/>, and NaN otherwise. A constant functional has no estimate and nothing to
/// sample; frozen and unavailable readings fail every requirement.
/// </summary>
public readonly record struct EssDiagnostic(DiagnosticState State, double Value) : IFormattable
{
    public static EssDiagnostic Unavailable => new(DiagnosticState.Unavailable, double.NaN);
    public static EssDiagnostic Constant => new(DiagnosticState.Constant, double.NaN);
    public static EssDiagnostic Frozen => new(DiagnosticState.Frozen, double.NaN);
    public static EssDiagnostic Estimate(double value) => new(DiagnosticState.Estimated, value);

    /// <summary>Whether this reading satisfies a minimum: an estimate at or above it, or a constant functional (exempt).</summary>
    public bool MeetsMinimum(double minimum)
        => State == DiagnosticState.Constant || (State == DiagnosticState.Estimated && Value >= minimum);

    /// <summary>The least favourable reading of a pair: constant, then estimates by value (smaller is worse), then
    /// frozen, then unavailable.</summary>
    public static EssDiagnostic Worst(EssDiagnostic a, EssDiagnostic b)
    {
        int ra = Rank(a.State), rb = Rank(b.State);
        if (ra != rb) return ra > rb ? a : b;
        return a.State == DiagnosticState.Estimated && a.Value < b.Value ? a : b;
    }

    private static int Rank(DiagnosticState state) => state switch
    {
        DiagnosticState.Constant => 0,
        DiagnosticState.Estimated => 1,
        DiagnosticState.Frozen => 2,
        _ => 3,
    };

    public override string ToString() => ToString(null, CultureInfo.InvariantCulture);

    public string ToString(string? format, IFormatProvider? formatProvider) => State switch
    {
        DiagnosticState.Estimated => Value.ToString(format, formatProvider),
        DiagnosticState.Constant => "constant",
        DiagnosticState.Frozen => "frozen",
        _ => "unavailable",
    };
}

/// <summary>
/// One functional's rank-normalized convergence diagnostics (Vehtari, Gelman, Simpson, Carpenter and Bürkner 2021):
/// split R̂ on the rank-normalized draws (<see cref="BulkRHat"/>) and on their folded distances from the median
/// (<see cref="TailRHat"/>), <see cref="RHat"/> the less favourable of the two; and the effective sample size of the
/// bulk (<see cref="BulkEss"/>) and of the 5% and 95% quantiles (<see cref="TailEss"/>, the smaller).
/// </summary>
public readonly record struct ConvergenceDiagnostic(
    RHatDiagnostic RHat, RHatDiagnostic BulkRHat, RHatDiagnostic TailRHat, EssDiagnostic BulkEss, EssDiagnostic TailEss)
{
    /// <summary>The less favourable of the bulk and tail effective sample sizes.</summary>
    public EssDiagnostic Ess => EssDiagnostic.Worst(BulkEss, TailEss);

    internal static ConvergenceDiagnostic All(DiagnosticState state) => state switch
    {
        DiagnosticState.Constant => new(RHatDiagnostic.Constant, RHatDiagnostic.Constant, RHatDiagnostic.Constant,
                                        EssDiagnostic.Constant, EssDiagnostic.Constant),
        DiagnosticState.Frozen => new(RHatDiagnostic.Frozen, RHatDiagnostic.Frozen, RHatDiagnostic.Frozen,
                                      EssDiagnostic.Frozen, EssDiagnostic.Frozen),
        _ => new(RHatDiagnostic.Unavailable, RHatDiagnostic.Unavailable, RHatDiagnostic.Unavailable,
                 EssDiagnostic.Unavailable, EssDiagnostic.Unavailable),
    };
}

public static partial class ChainDiagnostics
{
    /// <summary>
    /// Rank-normalized split R̂ and bulk and tail ESS for one scalar functional, each row of <paramref name="chains"/>
    /// one chain's draws, all of the same length (Vehtari et al. 2021). Each chain is split into halves, so one chain
    /// suffices. The failure policy is R̂'s: every draw the same is <see cref="DiagnosticState.Constant"/>; every
    /// split chain constant but not all at one value is <see cref="DiagnosticState.Frozen"/>; fewer than four draws
    /// per chain, or a non-finite draw, is <see cref="DiagnosticState.Unavailable"/>.
    /// <para>The effective sample size follows Stan: within-chain autocovariances by FFT, combined across chains,
    /// Geyer's initial positive and then monotone sequence of paired autocorrelations, and the bound
    /// τ̂ ≥ 1/log₁₀ S, so antithetic chains may exceed S but not S·log₁₀ S. The cost is O(S log S).</para>
    /// </summary>
    public static ConvergenceDiagnostic Diagnose(double[][] chains)
    {
        ArgumentNullException.ThrowIfNull(chains);
        int c = chains.Length;
        if (c < 1) return ConvergenceDiagnostic.All(DiagnosticState.Unavailable);
        int n = chains[0].Length;
        foreach (double[] chain in chains)
        {
            ArgumentNullException.ThrowIfNull(chain);
            if (chain.Length != n) throw new ArgumentException("Every chain must hold the same number of draws.", nameof(chains));
            foreach (double v in chain)
                if (!double.IsFinite(v)) return ConvergenceDiagnostic.All(DiagnosticState.Unavailable);
        }
        if (n < 4) return ConvergenceDiagnostic.All(DiagnosticState.Unavailable);

        double[][] split = Split(chains);
        switch (Spread(split))
        {
            case DiagnosticState.Constant: return ConvergenceDiagnostic.All(DiagnosticState.Constant);
            case DiagnosticState.Frozen: return ConvergenceDiagnostic.All(DiagnosticState.Frozen);
        }

        double[][] z = RankNormalize(split);
        RHatDiagnostic bulkRHat = ClassicRHat(z);
        EssDiagnostic bulkEss = StanEss(z);

        double median = PooledQuantile(split, 0.5);
        double[][] folded = Map(split, v => Math.Abs(v - median));
        RHatDiagnostic tailRHat = Spread(folded) switch
        {
            DiagnosticState.Constant => RHatDiagnostic.Constant,
            DiagnosticState.Frozen => RHatDiagnostic.Frozen,
            _ => ClassicRHat(RankNormalize(folded)),
        };

        double q05 = PooledQuantile(split, 0.05), q95 = PooledQuantile(split, 0.95);
        EssDiagnostic tailEss = EssDiagnostic.Worst(IndicatorEss(split, q05), IndicatorEss(split, q95));

        ReadOnlySpan<RHatDiagnostic> both = [bulkRHat, tailRHat];
        return new ConvergenceDiagnostic(RHatDiagnostic.Worst(both), bulkRHat, tailRHat, bulkEss, tailEss);
    }

    // Each chain's first and last ⌊n/2⌋ draws (an odd middle draw is dropped).
    private static double[][] Split(double[][] chains)
    {
        int half = chains[0].Length / 2;
        var split = new double[2 * chains.Length][];
        for (int j = 0; j < chains.Length; j++)
        {
            split[2 * j] = chains[j][..half];
            split[2 * j + 1] = chains[j][^half..];
        }
        return split;
    }

    // Constant: one value everywhere. Frozen: each chain constant, the chains not all equal. Otherwise estimated.
    private static DiagnosticState Spread(double[][] chains)
    {
        double first = chains[0][0];
        bool allEqual = true, eachConstant = true;
        foreach (double[] chain in chains)
            for (int i = 0; i < chain.Length; i++)
            {
                if (chain[i] != first) allEqual = false;
                if (chain[i] != chain[0]) { eachConstant = false; break; }
            }
        if (allEqual) return DiagnosticState.Constant;
        return eachConstant ? DiagnosticState.Frozen : DiagnosticState.Estimated;
    }

    // z = Φ⁻¹((r − 3/8)/(S + 1/4)) with r the pooled rank, ties averaged (Blom's offsets, as Vehtari et al. use).
    private static double[][] RankNormalize(double[][] chains)
    {
        int m = chains.Length, len = chains[0].Length, s = m * len;
        var values = new double[s];
        var index = new int[s];
        for (int j = 0, k = 0; j < m; j++)
            for (int i = 0; i < len; i++, k++) { values[k] = chains[j][i]; index[k] = k; }
        Array.Sort(values, index);

        var ranks = new double[s];
        for (int a = 0; a < s;)
        {
            int b = a;
            while (b + 1 < s && values[b + 1] == values[a]) b++;
            double rank = 0.5 * (a + b) + 1.0;   // the average of the 1-based ranks a+1 … b+1
            for (int k = a; k <= b; k++) ranks[index[k]] = rank;
            a = b + 1;
        }

        var z = new double[m][];
        for (int j = 0, k = 0; j < m; j++)
        {
            z[j] = new double[len];
            for (int i = 0; i < len; i++, k++) z[j][i] = Normal.StandardQuantile((ranks[k] - 0.375) / (s + 0.25));
        }
        return z;
    }

    private static RHatDiagnostic ClassicRHat(double[][] chains)
    {
        int m = chains.Length;
        var means = new double[m];
        var variances = new double[m];
        for (int j = 0; j < m; j++)
        {
            RunningMoments mo = RunningMoments.Of(chains[j]);
            means[j] = mo.Mean;
            variances[j] = mo.Variance;
        }
        return RHat(means, variances, chains[0].Length);
    }

    // The ESS of the indicator x ≤ q; an indicator that is the same everywhere has nothing to estimate.
    private static EssDiagnostic IndicatorEss(double[][] chains, double q)
    {
        double[][] ind = Map(chains, v => v <= q ? 1.0 : 0.0);
        return Spread(ind) switch
        {
            DiagnosticState.Constant => EssDiagnostic.Constant,
            DiagnosticState.Frozen => EssDiagnostic.Frozen,
            _ => StanEss(ind),
        };
    }

    // Stan's effective sample size of M chains of N draws (the posterior package's ess_rfun).
    private static EssDiagnostic StanEss(double[][] chains)
    {
        int m = chains.Length, len = chains[0].Length;
        var acov = new double[m][];
        var chainMean = new double[m];
        for (int j = 0; j < m; j++) (acov[j], chainMean[j]) = Autocovariance(chains[j]);

        double meanVar = 0.0;
        for (int j = 0; j < m; j++) meanVar += acov[j][0];
        meanVar = meanVar / m * len / (len - 1.0);
        double varPlus = meanVar * (len - 1.0) / len;
        if (m > 1)
        {
            double grand = 0.0;
            foreach (double cm in chainMean) grand += cm;
            grand /= m;
            double between = 0.0;
            foreach (double cm in chainMean) between += (cm - grand) * (cm - grand);
            varPlus += between / (m - 1);
        }

        double Rho(int t)
        {
            double a = 0.0;
            for (int j = 0; j < m; j++) a += acov[j][t];
            return 1.0 - (meanVar - a / m) / varPlus;
        }

        var rho = new double[len];
        rho[0] = 1.0;
        double even = 1.0, odd = Rho(1);
        rho[1] = odd;
        int t = 0;
        while (t < len - 5 && !double.IsNaN(even + odd) && even + odd > 0.0)
        {
            t += 2;
            even = Rho(t);
            odd = Rho(t + 1);
            if (even + odd >= 0.0)
            {
                rho[t] = even;
                rho[t + 1] = odd;
            }
        }
        int maxT = t;
        if (even > 0.0) rho[maxT] = even;   // the improved estimate's last term

        for (int u = 2; u <= maxT - 2; u += 2)   // Geyer's initial monotone sequence on the pairs
            if (rho[u] + rho[u + 1] > rho[u - 2] + rho[u - 1])
            {
                rho[u] = 0.5 * (rho[u - 2] + rho[u - 1]);
                rho[u + 1] = rho[u];
            }

        double tau = -1.0 + rho[maxT];
        for (int u = 0; u < maxT; u++) tau += 2.0 * rho[u];
        double total = m * (double)len;
        tau = Math.Max(tau, 1.0 / Math.Log10(total));
        return EssDiagnostic.Estimate(total / tau);
    }

    // The biased autocovariance (1/N) Σ (x_i − x̄)(x_{i+t} − x̄) at every lag, by FFT on a zero-padded copy.
    private static (double[] Acov, double Mean) Autocovariance(double[] x)
    {
        int len = x.Length;
        double mean = 0.0;
        foreach (double v in x) mean += v;
        mean /= len;
        int size = 1;
        while (size < 2 * len) size <<= 1;
        var data = new Complex[size];
        for (int i = 0; i < len; i++) data[i] = new Complex(x[i] - mean, 0.0);
        Fft.Forward(data);
        for (int i = 0; i < size; i++) data[i] = new Complex(data[i].Real * data[i].Real + data[i].Imaginary * data[i].Imaginary, 0.0);
        Fft.Inverse(data);
        var acov = new double[len];
        for (int i = 0; i < len; i++) acov[i] = data[i].Real / len;
        return (acov, mean);
    }

    // Quantile type 7 (linear interpolation between order statistics) of all draws together.
    private static double PooledQuantile(double[][] chains, double p)
    {
        int len = chains[0].Length;
        var all = new double[chains.Length * len];
        for (int j = 0; j < chains.Length; j++) Array.Copy(chains[j], 0, all, j * len, len);
        Array.Sort(all);
        double h = (all.Length - 1) * p;
        int lo = (int)Math.Floor(h);
        return lo + 1 < all.Length ? all[lo] + (h - lo) * (all[lo + 1] - all[lo]) : all[lo];
    }

    private static double[][] Map(double[][] chains, Func<double, double> f)
    {
        var mapped = new double[chains.Length][];
        for (int j = 0; j < chains.Length; j++)
        {
            mapped[j] = new double[chains[j].Length];
            for (int i = 0; i < chains[j].Length; i++) mapped[j][i] = f(chains[j][i]);
        }
        return mapped;
    }
}
