using System;
using Maths.Samplers.FlatHistogram;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>
/// The sampler suite's comparative assertions, each with the seeds its test uses, rerun across seeds by
/// <see cref="SeedSweep"/> (step 3.5). An assertion appears as written in its test, as Lower &lt; Upper. Set
/// <c>SEED_SWEEP</c> to the number of seeds per axis; these tests have no simulated data, so only the chain axis runs.
/// </summary>
public sealed class ComparativeSweep
{
    private readonly ITestOutputHelper _out;
    public ComparativeSweep(ITestOutputHelper output) => _out = output;

    internal static readonly SweptClaim[] Claims =
    [
        new("ParallelTemperingTests.Tempering_VisitsBothModes", DataSeed: null, ChainSeeds: [1, 7], (_, c) =>
        {
            var r = ParallelTemperingTests.SingleAgainstTempered(c[0], c[1]);
            return [new("single chain's right-mode fraction < 0.05", r.SingleFraction, 0.05),
                    new("|cold replica's right-mode fraction − 0.5| < 0.15", Math.Abs(r.TemperedFraction - 0.5), 0.15)];
        }),
        new("SwapScheduleTests.AdaptLadder_EqualizesRejection", DataSeed: null, ChainSeeds: [3, 5], (_, c) =>
        {
            var r = SwapScheduleTests.AdaptAndMeasure(c[0], c[1]);
            return [new("rejection max/min: adapted < uniform ladder", SwapScheduleTests.Spread(r.After), SwapScheduleTests.Spread(r.Before))];
        }),
        new("SamcLearningTests.DoubleWell_LearnsSamcsLimit_CloserTheLongerItLearns", DataSeed: null, ChainSeeds: [21], (_, c) =>
        {
            double e1 = SamcLearningTests.WellLearningError(4_000, c[0]), e2 = SamcLearningTests.WellLearningError(16_000, c[0]),
                   e3 = SamcLearningTests.WellLearningError(64_000, c[0]);
            return [new("RMS to θ*: 16k iterations < 4k", e2, e1), new("RMS to θ*: 64k iterations < 16k", e3, e2)];
        }),
        new("WhamTests.DoubleWell_LearningFromTheWarmStart_LeadsTheFlatStart_AtEqualIterations", DataSeed: null, ChainSeeds: [57, 61], (_, c) =>
        {
            var r = WhamTests.WarmAgainstFlat(WhamTests.WellError, WhamTests.WellPilot(c[0]), 32_000, c[1]);
            return [new("RMS to θ*: warm start's 16k < flat start's 16k", r.WarmHalf, r.FlatCurve[1])];
        }),
        new("WhamTests.Potts_LearningFromTheWarmStart_LeadsTheFlatStart_AtEqualIterations", DataSeed: null, ChainSeeds: [59, 63], (_, c) =>
        {
            var r = WhamTests.WarmAgainstFlat(WhamTests.PottsError, WhamTests.PottsPilot(c[0]), 100_000, c[1]);
            return [new("RMS to θ*: warm start's 50k < flat start's 50k", r.WarmHalf, r.FlatCurve[1])];
        }),
        // A reading beside the claims, no test's assertion: the same half-length comparison along the learning's length.
        new("WhamTests.WarmAgainstFlat_AlongTheLength_DoubleWell", DataSeed: null, ChainSeeds: [57, 61], (_, c) =>
            AlongTheLength(WhamTests.WellError, WhamTests.WellPilot(c[0]), 32_000, c[1])),
        new("WhamTests.WarmAgainstFlat_AlongTheLength_Potts", DataSeed: null, ChainSeeds: [59, 63], (_, c) =>
            AlongTheLength(WhamTests.PottsError, WhamTests.PottsPilot(c[0]), 100_000, c[1])),
    ];

    // The warm start's error at each length against the flat start's at twice it, from the bare start (against L/4) to 2L.
    private static Comparison[] AlongTheLength(Func<long, TemperingPilot?, int, double> error, TemperingPilot pilot, long length, int seed)
    {
        var near = WhamTests.WarmAgainstFlat(error, pilot, length, seed);
        var far = WhamTests.WarmAgainstFlat(error, pilot, 4 * length, seed);
        return [new($"warm start's 0 < flat start's {length / 4}", near.WarmCurve[0], near.FlatCurve[0]),
                new($"warm start's {length / 4} < flat start's {length / 2}", near.WarmCurve[1], near.FlatCurve[1]),
                new($"warm start's {length / 2} < flat start's {length}", near.WarmCurve[2], near.FlatCurve[2]),
                new($"warm start's {length} < flat start's {2 * length}", far.WarmCurve[1], far.FlatCurve[1]),
                new($"warm start's {2 * length} < flat start's {4 * length}", far.WarmCurve[2], far.FlatCurve[2])];
    }

    [Fact]
    public void ComparativeAssertions_AcrossSeeds()
    {
        if (SeedSweep.Seeds is not int seeds) return;   // heavy; set SEED_SWEEP to run
        SeedSweep.Run("Maths.Samplers.Tests", Claims, seeds, _out);
    }
}
