using System;
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
    ];

    [Fact]
    public void ComparativeAssertions_AcrossSeeds()
    {
        if (SeedSweep.Seeds is not int seeds) return;   // heavy; set SEED_SWEEP to run
        SeedSweep.Run("Maths.Samplers.Tests", Claims, seeds, _out);
    }
}
