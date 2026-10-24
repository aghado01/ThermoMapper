using System;
using System.Collections.Generic;
using Maths.Rng;

namespace Maths.Regression.Changepoint;

/// <summary>
/// One exact draw from the changepoint posterior: a segmentation and a level for each segment.
/// <see cref="SegmentEnds"/> holds each segment's last index, 0-based and ascending, the last being n − 1, so a break
/// follows every entry but the last. <see cref="Levels"/>[s] is segment s's level.
/// </summary>
public sealed record ChangepointDraw(int[] SegmentEnds, double[] Levels)
{
    /// <summary>The draw's piecewise-constant curve: each point's segment level.</summary>
    public double[] Curve()
    {
        var curve = new double[SegmentEnds[^1] + 1];
        int start = 0;
        for (int s = 0; s < SegmentEnds.Length; s++)
        {
            for (int i = start; i <= SegmentEnds[s]; i++) curve[i] = Levels[s];
            start = SegmentEnds[s] + 1;
        }
        return curve;
    }
}

/// <summary>
/// The exact changepoint posterior as the dynamic program leaves it: the segment model's tables (each segment's log
/// evidence and level moments), the forward and backward tables, the segment-count posterior, and each segment's
/// posterior probability. Indices are 0-based and segments inclusive.
/// <para>Draws are exact and independent, by forward filtering and backward sampling: k from P(k | y), then the
/// segments from the right end, the start of the s-th given its end b in proportion to F_{s−1}(a)·D(a+1, b), then each
/// level from the segment model (<see cref="ISegmentModel.SampleLevel"/>). Given k the prior is uniform over
/// segmentations, so this is the posterior over segmentations exactly.</para>
/// The tables hold 4(n+1)² + 2(K+1)(n+1) numbers, and the posterior keeps the segment model for the level draws.
/// </summary>
public sealed class ChangepointPosterior
{
    private readonly ISegmentModel _model;
    private readonly double[,] _logD, _mean, _var, _fwd, _bwd, _segment;
    private readonly double[] _countPosterior;

    internal ChangepointPosterior(ISegmentModel model, double[,] logD, double[,] mean, double[,] var, double[,] fwd,
                                  double[,] bwd, double[,] segment, double[] countPosterior, double logEvidence)
    {
        _model = model;
        _logD = logD; _mean = mean; _var = var; _fwd = fwd; _bwd = bwd; _segment = segment;
        _countPosterior = countPosterior;
        LogEvidence = logEvidence;
    }

    /// <summary>The number of data points n.</summary>
    public int Count => _model.Count;

    /// <summary>The largest segment count K the prior allows.</summary>
    public int MaxSegments => _countPosterior.Length;

    /// <summary>log P(y), as <see cref="ChangepointResult.LogEvidence"/>.</summary>
    public double LogEvidence { get; }

    /// <summary>log D(first..last): the segment's evidence with its level integrated out.</summary>
    public double SegmentLogEvidence(int first, int last)
    {
        CheckSegment(first, last);
        return _logD[first + 1, last + 1];
    }

    /// <summary>The posterior mean and variance of the segment's level, given that it is a segment.</summary>
    public (double Mean, double Variance) SegmentLevel(int first, int last)
    {
        CheckSegment(first, last);
        return (_mean[first + 1, last + 1], _var[first + 1, last + 1]);
    }

    /// <summary>P(first..last is a segment | y).</summary>
    public double SegmentProbability(int first, int last)
    {
        CheckSegment(first, last);
        return _segment[first + 1, last + 1];
    }

    /// <summary>The forward table F_s(length): log Σ, over segmentations of the first <paramref name="length"/> points
    /// into <paramref name="segments"/> segments, of Π D(segment).</summary>
    public double LogPrefixEvidence(int segments, int length)
    {
        if (segments < 0 || segments > MaxSegments) throw new ArgumentOutOfRangeException(nameof(segments));
        if (length < 0 || length > Count) throw new ArgumentOutOfRangeException(nameof(length));
        return _fwd[segments, length];
    }

    /// <summary>The backward table B_s(start): the same sum over the points from index <paramref name="start"/> on.</summary>
    public double LogSuffixEvidence(int segments, int start)
    {
        if (segments < 0 || segments > MaxSegments) throw new ArgumentOutOfRangeException(nameof(segments));
        if (start < 0 || start > Count) throw new ArgumentOutOfRangeException(nameof(start));
        return _bwd[segments, start];
    }

    /// <summary>
    /// log p(y_unit | draw) for each of <paramref name="draws"/> and each unit (OQ17). A draw's segmentation is global, so
    /// every unit conditions on it. A segment's level is a group latent: a unit that holds all of the segment's
    /// observations integrates the level out, and its term for them is the segment's evidence D(a, b) from the tables
    /// (the level under its prior); a unit that holds only some of them conditions on the drawn level, each observation's
    /// term the segment model's likelihood there.
    /// </summary>
    /// <param name="units">The partition of the n observations; null takes each observation as its own unit.</param>
    public PointwiseLogDensity PointwiseLogDensity(IReadOnlyList<ChangepointDraw> draws, UnitPartition? units = null)
    {
        ArgumentNullException.ThrowIfNull(draws);
        int n = Count;
        UnitPartition partition = units ?? UnitPartition.EachObservation(n);
        if (partition.Observations != n) throw new ArgumentException("The partition must cover the n observations.", nameof(units));
        var values = new double[draws.Count, partition.Units];
        var segmentOf = new int[n];
        var tally = new TreatmentTally();
        for (int d = 0; d < draws.Count; d++)
        {
            ChangepointDraw draw = draws[d];
            int start = 0;
            for (int s = 0; s < draw.SegmentEnds.Length; s++)
            {
                for (int i = start; i <= draw.SegmentEnds[s]; i++) segmentOf[i] = s;
                start = draw.SegmentEnds[s] + 1;
            }
            for (int u = 0; u < partition.Units; u++)
            {
                int[] members = partition.Members[u];
                double sum = 0.0;
                // Members ascend and segments are contiguous, so each segment's members in the unit form one run.
                for (int k = 0; k < members.Length;)
                {
                    int s = segmentOf[members[k]], run = k;
                    while (run < members.Length && segmentOf[members[run]] == s) run++;
                    int first = s == 0 ? 0 : draw.SegmentEnds[s - 1] + 1, last = draw.SegmentEnds[s];
                    if (run - k == last - first + 1)
                    {
                        sum += _logD[first + 1, last + 1];
                        tally.Integrated();
                    }
                    else
                    {
                        for (int j = k; j < run; j++) sum += _model.LogLikelihood(members[j], draw.Levels[s]);
                        tally.Conditioned();
                    }
                    k = run;
                }
                values[d, u] = sum;
            }
        }
        var latents = new[]
        {
            new LatentReport("segmentation", LatentScope.Global, LatentTreatment.Conditioned),
            new LatentReport("segment levels", LatentScope.Group, tally.Result),
        };
        return new PointwiseLogDensity(values, partition, latents);
    }

    /// <summary><paramref name="draws"/> exact, independent draws on one stream seeded by <paramref name="seed"/>.</summary>
    public IReadOnlyList<ChangepointDraw> Sample(int draws, int seed)
    {
        if (draws < 0) throw new ArgumentOutOfRangeException(nameof(draws));
        var rng = new Xoshiro256PlusPlus(seed);
        var result = new ChangepointDraw[draws];
        for (int d = 0; d < draws; d++) result[d] = Sample(rng);
        return result;
    }

    /// <summary>One exact draw: the segment count, the segmentation given it, then the levels.</summary>
    public ChangepointDraw Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        int k = Pick(_countPosterior, rng) + 1;
        var ends = new int[k];
        var levels = new double[k];

        // The s-th segment ends at b (1-based); its start a+1 takes a in [s−1, b−1] with weight F_{s−1}(a)·D(a+1, b)/F_s(b).
        int b = Count;
        for (int s = k; s > 1; s--)
        {
            double target = rng.NextDouble(), acc = 0.0, total = _fwd[s, b];
            int start = b - 1;
            for (int a = s - 1; a < b; a++)
            {
                double w = Math.Exp(_fwd[s - 1, a] + _logD[a + 1, b] - total);
                if (!(w > 0.0)) continue;
                acc += w;
                start = a;
                if (acc >= target) break;
            }
            ends[s - 1] = b - 1;
            levels[s - 1] = _model.SampleLevel(start + 1, b, rng);
            b = start;
        }
        ends[0] = b - 1;
        levels[0] = _model.SampleLevel(1, b, rng);
        return new ChangepointDraw(ends, levels);
    }

    // An index drawn in proportion to the probabilities; rounding that leaves the sum short of the draw falls to the last
    // index with positive probability.
    private static int Pick(double[] probabilities, Xoshiro256PlusPlus rng)
    {
        double target = rng.NextDouble(), acc = 0.0;
        int pick = -1;
        for (int i = 0; i < probabilities.Length; i++)
        {
            if (!(probabilities[i] > 0.0)) continue;
            acc += probabilities[i];
            pick = i;
            if (acc >= target) break;
        }
        return pick;
    }

    private void CheckSegment(int first, int last)
    {
        if (first < 0 || first >= Count) throw new ArgumentOutOfRangeException(nameof(first));
        if (last < first || last >= Count) throw new ArgumentOutOfRangeException(nameof(last));
    }
}
