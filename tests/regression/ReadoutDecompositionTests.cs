using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The decompositions a BASS sensitivity readout runs, counted at the LinAlg kernels (spectral-solver brief, C). Equal outputs also
/// pass while duplicate work remains, so these count the decompositions the reuse removed. Symmetric solves are counted whichever
/// kernel the build dispatches them to.
/// </summary>
public sealed class ReadoutDecompositionTests
{
    /// <summary>C2: on the sampled route with a declared W_y and output patterns asked, W_y (m × m) is decomposed once, though its
    /// validation, the output side's symmetric root and the sampled route's factor each read it (three solves before).</summary>
    [Fact]
    public void DirectReadout_DecomposesTheOutputMetricOnce()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int m = data.Response.Locations;
        var spec = new ActiveSubspaceSpec
        {
            OutputMetric = new OutputMetric(FunctionalActiveSubspaceTests.RandomMetric(m, 13)), OutputPatterns = 1,
            Route = ActiveSubspaceRoute.Sampled, Points = 64, Draws = [0, 1, 2],
        };
        IReadOnlyList<Decomposition> work = Run(() => fit.ActiveSubspace(data, spec));

        Assert.Equal(1, work.Count(d => Symmetric(d) && d.Rows == m));
    }

    /// <summary>C3: on the exact route the output side's term factor is decomposed once per structure. Three copies of one draw add
    /// two symmetric solves when patterns are asked, W_y's symmetric root and the one structure's term Gram (four before).</summary>
    [Fact]
    public void ExactReadout_DecomposesEachStructuresTermGramOnce()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        var three = new ActiveSubspaceSpec { Draws = [0, 0, 0] };
        int without = Run(() => fit.ActiveSubspace(data, three)).Count(Symmetric);
        int with = Run(() => fit.ActiveSubspace(data, three with { OutputPatterns = 1 })).Count(Symmetric);

        Assert.Equal(2, with - without);
    }

    /// <summary>C1 and C2: adequacy reads the direct readout's output posterior and W_y's root as the readout built them, so it adds
    /// only the span's SVD to the readout's decompositions (a second posterior and two metric solves before).</summary>
    [Fact]
    public void Adequacy_AddsOnlyTheSpansDecomposition()
    {
        (BassData data, BassFit direct, BassPredictor[] predictors, double[,] y, _, _) = RankAndOutputSideTests.AdequacyFixture.Value;
        BassBasisFit basis = BassComponents.Fit(BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 2 }), predictors, BassOutputComponentsTests.Config());
        var spec = new ActiveSubspaceSpec
        {
            OutputMetric = new OutputMetric(FunctionalActiveSubspaceTests.RandomMetric(data.Response.Locations, 5)), OutputPatterns = 1,
            Draws = [0, 3, 6],
        };
        IReadOnlyList<Decomposition> read = Run(() => direct.ActiveSubspace(data, spec));
        IReadOnlyList<Decomposition> adequacy = Run(() => basis.Adequacy(direct, data, spec));

        Assert.Equal(read.Count(Symmetric), adequacy.Count(Symmetric));
        Assert.Equal(read.Count(d => !Symmetric(d)) + 1, adequacy.Count(d => !Symmetric(d)));
    }

    private static bool Symmetric(Decomposition d) => d.Kernel != DecompositionKernel.ThinSvd;

    private static IReadOnlyList<Decomposition> Run(Action readout)
    {
        using DecompositionLog log = DecompositionLog.Begin();
        readout();
        return log.Entries;
    }
}
