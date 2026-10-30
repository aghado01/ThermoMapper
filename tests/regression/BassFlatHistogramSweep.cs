using Maths.Samplers.Tests;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's flat-histogram comparative assertions (N17 Phase 2), each with the seeds its test uses, rerun across seeds by
/// <see cref="SeedSweep"/> beside <see cref="ComparativeSweep"/>'s, which the main session may fold them into. Set
/// <c>SEED_SWEEP</c> to the number of seeds per axis; only the chain axis runs.
/// </summary>
public sealed class BassFlatHistogramSweep
{
    private readonly ITestOutputHelper _out;
    public BassFlatHistogramSweep(ITestOutputHelper output) => _out = output;

    internal static readonly SweptClaim[] Claims =
    [
        new("BassFlatHistogramTests.Pilot_LearningFromTheWarmStart_ReachesTheFlatStartsAccuracy_InHalfTheIterations", DataSeed: null,
            ChainSeeds: [7701, 7703], (_, c) =>
            {
                var r = BassFlatHistogramTests.WarmAgainstFlat(c[0], c[1], 8_000);
                return [new("RMS to θ*: warm start's 4k < flat start's 8k", r.WarmHalf, r.Flat)];
            }),
    ];

    [Fact]
    public void ComparativeAssertions_AcrossSeeds()
    {
        if (SeedSweep.Seeds is not int seeds) return;   // heavy; set SEED_SWEEP to run
        SeedSweep.Run("Maths.Regression.Tests.BassFlatHistogram", Claims, seeds, _out);
    }
}
