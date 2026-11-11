using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Bgp;

public sealed partial class BgpSampler
{
    /// <summary>
    /// BGP's retained-parameters content (E16, the conformal layer's slice 11): for each retained draw (t, η, s²), the regression
    /// function g = m + f drawn from its posterior given the data, jointly at the training inputs and at the inputs of
    /// <paramref name="xCases"/>, so every unit's term reads the draw in sample, N(y | g(x), σ²), as BASS's terms do. The full step takes
    /// it, where the collapsed content (<see cref="PointwiseLogDensity(BgpResult, UnitPartition?)"/>) leaves each unit out.
    /// <list type="bullet">
    /// <item><b>The draw.</b> Given (t, η) and s², g at the points is Gaussian with the kriging mean and s² times the universal-kriging
    /// covariance (<see cref="GpRegression"/>, the constant mean's flat prior included), and s² | t, η, Y is the run's own exact draw, so
    /// (t, η, s², g) is a draw from the joint posterior. A draw repeating its predecessor's (t, η) shares its mean and factor, O((n + T)³)
    /// once per accepted state, and takes its own standard normals.</item>
    /// <item><b>The points.</b> The training inputs and the cases' inputs, each distinct mapped input once: a case at a training input reads
    /// that input's g, so its candidate term is the training term, bit for bit, and two units at one input share one value.</item>
    /// <item><b>The factor.</b> The posterior covariance is near singular where inputs are close next to the bandwidth, so it is factored by
    /// Cholesky with diagonal pivoting, stopped where every remaining pivot is below 10⁻¹² of the largest variance (Harbrecht, Peters and
    /// Schneider 2012): each point's variance falls short by at most that share, and a draw takes one normal per column kept.</item>
    /// <item><b>Pairing.</b> The training density and the candidates must come from one drawn function: a second call draws g afresh, and
    /// its candidates would not be drawn jointly with the first call's training terms.</item>
    /// </list>
    /// The normals come from a stream per chain derived from <paramref name="seed"/>, so the draws repeat for one seed and one set of cases.
    /// </summary>
    /// <param name="result">A run of this sampler, whose retained draws are read.</param>
    /// <param name="xCases">The cases' inputs (cases × D) in the original units, where the candidate density will be read; null for none.</param>
    /// <param name="seed">The draws' seed.</param>
    public BgpDrawnFunction DrawFunction(BgpResult result, double[,]? xCases = null, int seed = 0)
    {
        ArgumentNullException.ThrowIfNull(result);
        int d = _gp.Dimension, n = _gp.Count;
        if (xCases is not null && xCases.GetLength(1) != d) throw new ArgumentException("The cases' inputs must have the training dimension.", nameof(xCases));
        double[,] training = _gp.MappedInputs;
        double[,] cases = xCases is null ? new double[0, d] : _gp.Transform.MapInputs(xCases);
        int count = cases.GetLength(0);

        // Each distinct mapped input once: the training inputs first, then the cases'.
        var index = new Dictionary<Row, int>();
        var points = new List<(double[,] From, int Row)>();
        int At(double[,] from, int row)
        {
            var key = new Row(from, row);
            if (index.TryGetValue(key, out int at)) return at;
            index.Add(key, points.Count);
            points.Add((from, row));
            return points.Count - 1;
        }
        var trainingAt = new int[n];
        for (int i = 0; i < n; i++) trainingAt[i] = At(training, i);
        var caseAt = new int[count];
        for (int t = 0; t < count; t++) caseAt[t] = At(cases, t);
        int p = points.Count;
        var mapped = new double[p, d];
        for (int a = 0; a < p; a++)
            for (int c = 0; c < d; c++) mapped[a, c] = points[a].From[points[a].Row, c];

        int total = 0;
        foreach (IReadOnlyList<BgpDraw> chain in result.BandwidthDraws) total += chain.Count;
        var values = new double[total, p];
        var noise = new double[total];
        var factor = new double[p * p];
        var z = new double[p];
        double[] mean = [];
        int rank = 0;
        GpFit? previous = null;
        int draw = 0, chainIndex = 0;
        foreach (IReadOnlyList<BgpDraw> chain in result.BandwidthDraws)
        {
            Xoshiro256PlusPlus rng = Stream(seed, chainIndex++);
            foreach (BgpDraw dr in chain)
            {
                if (previous is null || dr.Fit.Bandwidth != previous.Bandwidth || dr.Fit.NoiseRatio != previous.NoiseRatio)
                {
                    (mean, double[,] covariance) = _gp.JointPosteriorInFrame(dr.Fit, mapped);
                    rank = PivotedFactor(covariance, factor);
                    previous = dr.Fit;
                }
                double scale = Math.Sqrt(dr.SignalVariance);
                for (int k = 0; k < rank; k++) z[k] = Normal.SampleStandard(rng);
                for (int a = 0; a < p; a++)
                {
                    double s = 0.0;
                    int row = a * p;
                    for (int k = 0; k < rank; k++) s += factor[row + k] * z[k];
                    values[draw, a] = mean[a] + scale * s;
                }
                noise[draw] = dr.NoiseVariance;
                draw++;
            }
        }
        return new BgpDrawnFunction(values, trainingAt, caseAt, noise, [.. _gp.Responses], ChainLengths.Of(result.BandwidthDraws), _unknownScale);
    }

    // A pivot below this share of the largest variance ends the factor: the rest of the covariance is left out.
    private const double PivotTolerance = 1e-12;

    // Cholesky with diagonal pivoting of a positive semidefinite covariance given by its lower triangle (Harbrecht, Peters and Schneider
    // 2012): each step takes the largest remaining diagonal, and the factor stops where every remaining one is at most PivotTolerance of
    // the largest variance. Row a of the factor (p × p, row-major, in the points' own order) holds a's loadings on the columns kept, so
    // the covariance is LLᵀ up to the remaining diagonals, each below the tolerance. Returns the columns kept. An unpivoted factor with a
    // floor fails here: the posterior covariance of a smooth kernel at close inputs has a spectrum that falls below rounding, and floored
    // pivots let its rounding grow without bound.
    internal static int PivotedFactor(double[,] lower, double[] factor)
    {
        int p = lower.GetLength(0);
        Array.Clear(factor, 0, p * p);
        var remaining = new double[p];
        var order = new int[p];
        double largest = 0.0;
        for (int a = 0; a < p; a++)
        {
            remaining[a] = lower[a, a];
            order[a] = a;
            largest = Math.Max(largest, remaining[a]);
        }
        double tolerance = PivotTolerance * largest;
        int k = 0;
        for (; k < p; k++)
        {
            int best = k;
            for (int i = k + 1; i < p; i++) if (remaining[order[i]] > remaining[order[best]]) best = i;
            double pivot = remaining[order[best]];
            if (!(pivot > tolerance)) break;
            (order[k], order[best]) = (order[best], order[k]);
            int at = order[k];
            double root = Math.Sqrt(pivot);
            factor[at * p + k] = root;
            for (int i = k + 1; i < p; i++)
            {
                int q = order[i];
                double s = q >= at ? lower[q, at] : lower[at, q];
                for (int c = 0; c < k; c++) s -= factor[q * p + c] * factor[at * p + c];
                double v = s / root;
                factor[q * p + k] = v;
                remaining[q] -= v * v;
            }
        }
        return k;
    }

    // A stream per chain from the seed: SplitMix64 over the seed and the chain's index gives xoshiro's four words.
    private static Xoshiro256PlusPlus Stream(int seed, int chain)
    {
        ulong x = unchecked((ulong)(uint)seed * 0x9E3779B97F4A7C15UL ^ ((ulong)(uint)chain + 1UL) * 0xD1B54A32D192ED03UL);
        ulong Next()
        {
            x = unchecked(x + 0x9E3779B97F4A7C15UL);
            ulong v = x;
            v = unchecked((v ^ (v >> 30)) * 0xBF58476D1CE4E5B9UL);
            v = unchecked((v ^ (v >> 27)) * 0x94D049BB133111EBUL);
            return v ^ (v >> 31);
        }
        ulong s0 = Next(), s1 = Next(), s2 = Next(), s3 = Next();
        return Xoshiro256PlusPlus.Restore(s0, s1, s2, s3 | 1UL);
    }
}

/// <summary>
/// BGP's retained-parameters content (<see cref="BgpSampler.DrawFunction"/>): per retained draw, the regression function g = m + f at
/// each training input and each case's input, drawn jointly from its posterior, and the noise variance σ² = s²η. Each unit's term reads
/// the draw in sample, N(y | g(x), σ²), so the training density (<see cref="PointwiseLogDensity"/>) and the candidate density
/// (<see cref="CandidateLogDensity"/>) serve the full step and the split form alike; WAIC reads the training density too. The two must
/// come from this one object. Rows are draws, pooled over chains in chain order.
/// </summary>
public sealed class BgpDrawnFunction
{
    private readonly double[,] _values;   // draws × distinct points: g
    private readonly int[] _trainingAt, _caseAt;
    private readonly double[] _noise, _response;
    private readonly IReadOnlyList<int> _chains;
    private readonly bool _unknownScale;

    internal BgpDrawnFunction(double[,] values, int[] trainingAt, int[] caseAt, double[] noise, double[] response,
                              IReadOnlyList<int> chains, bool unknownScale)
    {
        _values = values;
        _trainingAt = trainingAt;
        _caseAt = caseAt;
        _noise = noise;
        _response = response;
        _chains = chains;
        _unknownScale = unknownScale;
    }

    /// <summary>The retained draws.</summary>
    public int Draws => _noise.Length;

    /// <summary>The cases g was drawn at, in the order given.</summary>
    public int Cases => _caseAt.Length;

    /// <summary>The training observations.</summary>
    public int Observations => _trainingAt.Length;

    /// <summary>The distinct inputs g was drawn at, the training inputs' and the cases', each once.</summary>
    public int Points => _values.GetLength(1);

    /// <summary>g at training observation <paramref name="observation"/> under draw <paramref name="draw"/>.</summary>
    public double AtObservation(int draw, int observation) => _values[draw, _trainingAt[observation]];

    /// <summary>g at case <paramref name="case"/> under draw <paramref name="draw"/>.</summary>
    public double AtCase(int draw, int @case) => _values[draw, _caseAt[@case]];

    /// <summary>σ² = s²η under draw <paramref name="draw"/> (σ² itself for a known σ²).</summary>
    public double NoiseVariance(int draw) => _noise[draw];

    /// <summary>
    /// The training density in the retained-parameters content: per draw and unit, Σ over the unit's observations of log N(yᵢ | g(xᵢ), σ²),
    /// with the latents reported conditioned on, none left out (<see cref="Maths.Regression.PointwiseLogDensity.LeavesUnitsOut"/> false).
    /// With <paramref name="cdf"/>, each observation's log Φ((yᵢ − g(xᵢ))/σ), units of one observation each.
    /// </summary>
    /// <param name="units">The partition; null gives each observation its own unit.</param>
    /// <param name="cdf">Adds each unit's log CDF, for the full step's predictive distribution.</param>
    public PointwiseLogDensity PointwiseLogDensity(UnitPartition? units = null, bool cdf = false)
    {
        int n = Observations;
        UnitPartition partition = units ?? UnitPartition.EachObservation(n);
        if (partition.Observations != n) throw new ArgumentException("The partition must cover the fitted observations.", nameof(units));
        if (cdf) CandidateSeam.OneObservationEach(partition, nameof(units));
        int s = Draws;
        var values = new double[s, partition.Units];
        double[,]? logCdfs = cdf ? new double[s, partition.Units] : null;
        for (int d = 0; d < s; d++)
        {
            (double logNorm, double variance, double sd) = Law(d);
            for (int u = 0; u < partition.Units; u++)
            {
                double sum = 0.0;
                foreach (int i in partition.Members[u])
                {
                    double e = _response[i] - _values[d, _trainingAt[i]];
                    sum += Term(logNorm, variance, e);
                    if (logCdfs is not null) logCdfs[d, u] = Normal.StandardLogCdf(e / sd);
                }
                values[d, u] = sum;
            }
        }
        return new PointwiseLogDensity(values, partition, Latents(), _chains) { LogCdf = logCdfs };
    }

    /// <summary>
    /// The candidate density: per draw, log N(y | g(x_t), σ²) at each value of each case's grid, with the training density's latent report,
    /// so the full step takes it beside <see cref="PointwiseLogDensity"/> from this object; at a training input with that observation's value
    /// it is the training term, bit for bit.
    /// </summary>
    /// <param name="grids">Each read case's candidate values, ascending.</param>
    /// <param name="strata">Each read case's stratum; null puts every case in stratum 0.</param>
    /// <param name="cdf">Adds each draw's log CDF at each candidate.</param>
    /// <param name="cases">The drawn cases to read, in this order (one grid each); null reads every case in order.</param>
    public CandidateLogDensity CandidateLogDensity(IReadOnlyList<double[]> grids, int[]? strata = null, bool cdf = false, int[]? cases = null)
    {
        ArgumentNullException.ThrowIfNull(grids);
        int[] read = cases ?? [.. System.Linq.Enumerable.Range(0, Cases)];
        foreach (int t in read)
            if ((uint)t >= (uint)Cases) throw new ArgumentOutOfRangeException(nameof(cases), $"Case {t} was not drawn: {Cases} cases were.");
        int[] offsets = CandidateSeam.Offsets(grids, read.Length);
        int s = Draws, columns = offsets[^1];
        var values = new double[s, columns];
        double[,]? logCdfs = cdf ? new double[s, columns] : null;
        for (int d = 0; d < s; d++)
        {
            (double logNorm, double variance, double sd) = Law(d);
            int at = 0;
            for (int c = 0; c < read.Length; c++)
            {
                double g = _values[d, _caseAt[read[c]]];
                foreach (double y in grids[c])
                {
                    double e = y - g;
                    values[d, at] = Term(logNorm, variance, e);
                    if (logCdfs is not null) logCdfs[d, at] = Normal.StandardLogCdf(e / sd);
                    at++;
                }
            }
        }
        return new CandidateLogDensity(values, grids, Latents(), strata, _chains) { LogCdf = logCdfs };
    }

    // One draw's Normal: −½ ln(2πσ²), σ² and σ, which the training and candidate terms share.
    private (double LogNorm, double Variance, double Sd) Law(int draw)
    {
        double v = _noise[draw];
        return (-0.5 * Math.Log(2.0 * Math.PI * v), v, Math.Sqrt(v));
    }

    private static double Term(double logNorm, double variance, double e) => logNorm - 0.5 * e * e / variance;

    // The retained content's latents: everything drawn and conditioned on, nothing integrated.
    private List<LatentReport> Latents()
    {
        var latents = new List<LatentReport> { new("bandwidth t", LatentScope.Global, LatentTreatment.Conditioned) };
        if (_unknownScale) latents.Add(new("noise ratio η", LatentScope.Global, LatentTreatment.Conditioned));
        latents.Add(new("regression function f", LatentScope.Global, LatentTreatment.Conditioned));
        latents.Add(new("constant mean m", LatentScope.Global, LatentTreatment.Conditioned));
        if (_unknownScale) latents.Add(new("signal variance s²", LatentScope.Global, LatentTreatment.Conditioned));
        return latents;
    }
}
