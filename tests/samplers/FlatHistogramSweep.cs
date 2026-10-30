using System;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>
/// The flat-histogram engine's comparative assertions, each with the seeds its test uses, rerun across seeds by
/// <see cref="SeedSweep"/> beside <see cref="ComparativeSweep"/>'s, which the main session may fold them into. Set
/// <c>SEED_SWEEP</c> to the number of seeds per axis; only the chain axis runs.
/// </summary>
public sealed class FlatHistogramSweep
{
    private readonly ITestOutputHelper _out;
    public FlatHistogramSweep(ITestOutputHelper output) => _out = output;

    internal static readonly SweptClaim[] Claims =
    [
        new("SamcLearningTests.DoubleWell_LearnsSamcsLimit_CloserTheLongerItLearns", DataSeed: null, ChainSeeds: [21], (_, c) =>
        {
            double e1 = SamcLearningTests.WellLearningError(4_000, c[0]), e2 = SamcLearningTests.WellLearningError(16_000, c[0]),
                   e3 = SamcLearningTests.WellLearningError(64_000, c[0]);
            return [new("RMS to θ*: 16k iterations < 4k", e2, e1), new("RMS to θ*: 64k iterations < 16k", e3, e2)];
        }),
    ];

    [Fact]
    public void ComparativeAssertions_AcrossSeeds()
    {
        if (SeedSweep.Seeds is not int seeds) return;   // heavy; set SEED_SWEEP to run
        SeedSweep.Run("Maths.Samplers.Tests.FlatHistogram", Claims, seeds, _out);
    }
}
