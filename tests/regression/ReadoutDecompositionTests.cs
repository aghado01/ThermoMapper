using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Maths.Rng;
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

    /// <summary>C4: a sampled draw's C is the Gram the route's running sums accumulated, kept by <see cref="DrawSubspace.FromFactor"/>
    /// rather than formed again from the factor. The subspace still comes from the factor's thin SVD, and with fewer factor rows than
    /// the dimension the draw falls back to C's own eigenspace.</summary>
    [Fact]
    public void FromFactor_KeepsTheAccumulatedGram()
    {
        var factor = new double[,] { { 1.0, 2.0, 0.5 }, { 0.3, -1.0, 2.0 }, { 0.7, 0.1, -0.4 }, { -1.2, 0.4, 0.9 } };
        double[,] c = Gram(factor);
        DrawSubspace draw = default!;
        IReadOnlyList<Decomposition> work = Run(() => draw = DrawSubspace.FromFactor(factor, c, 2));
        DrawSubspace fromC = DrawSubspace.FromC(c, 2);

        Assert.Same(c, draw.C);
        Assert.Equal([new Decomposition(DecompositionKernel.ThinSvd, 4, 3)], work);
        for (int i = 0; i < 3; i++) Assert.InRange(Math.Abs(draw.Eigenvalues[i] - fromC.Eigenvalues[i]), 0.0, 1e-13 * fromC.Eigenvalues[0]);

        var thin = new double[,] { { 1.0, 2.0, 0.5 } };
        double[,] ct = Gram(thin);
        DrawSubspace fallback = DrawSubspace.FromFactor(thin, ct, 2), direct = DrawSubspace.FromC(ct, 2);
        Assert.Same(ct, fallback.C);
        Assert.Equal(direct.Eigenvalues, fallback.Eigenvalues);
        Assert.Equal(direct.Basis, fallback.Basis);
        Assert.Throws<ArgumentException>(() => DrawSubspace.FromFactor(factor, new double[2, 2], 2));
    }

    /// <summary>C5: a projection of a whitened analysis reads the inverse root the readout built from the covariance's eigensystem, and
    /// decomposes nothing symmetric (one solve of L before). An analysis built by hand from the same parts, which carries no inverse,
    /// derives it from its map, and both project alike.</summary>
    [Fact]
    public void WhitenedProjection_ReadsTheReadoutsInverseRoot()
    {
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        double[][] x = JointSampleActiveSubspaceTests.Triangle(new Xoshiro256PlusPlus(67), 3000);
        var sample = new JointSample(new PredictionCases([new ContinuousColumn("x0", x[0]), new ContinuousColumn("x1", x[1])]));
        ActiveSubspaceAnalysis white = f.ActiveSubspace(df, new ActiveSubspaceSpec { Draws = [0, 5, 10], Sample = sample, Dimension = 2, Metric = InputMetric.Whitened });
        BassProjection projection = default!;
        IReadOnlyList<Decomposition> work = Run(() => projection = BassProjection.Of(f.Transform, white, 2));

        Assert.Equal(0, work.Count(Symmetric));
        double[,] l = white.Map, inverse = white.MapInverse;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                Assert.InRange(Math.Abs(l[i, 0] * inverse[0, j] + l[i, 1] * inverse[1, j] - (i == j ? 1.0 : 0.0)), 0.0, 1e-13);

        var byHand = new ActiveSubspaceAnalysis(white.Inputs, white.Excluded, white.Route, white.Metric, white.Map, white.Measure, white.Points,
                                                white.Posterior, white.MonteCarloError, white.Poincare, white.DrawIndex, white.Ranks, white.Outputs,
                                                white.MatrixPoincare, white.MatrixPoincareStatus, white.Boxes);
        BassProjection derived = default!;
        Assert.Equal(1, Run(() => derived = BassProjection.Of(f.Transform, byHand, 2)).Count(Symmetric));
        double[,] a = projection.Map, b = derived.Map;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) Assert.InRange(Math.Abs(a[i, j] - b[i, j]), 0.0, 1e-12 * Math.Max(1.0, Math.Abs(a[i, j])));
        Assert.Equal(1, Run(() => _ = (white with { Map = (double[,])white.Map.Clone() }).MapInverse).Count(Symmetric));
    }

    private static double[,] Gram(double[,] factor)
    {
        int n = factor.GetLength(0), p = factor.GetLength(1);
        var c = new double[p, p];
        for (int r = 0; r < n; r++)
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) c[i, j] += factor[r, i] * factor[r, j];
        return c;
    }

    private static bool Symmetric(Decomposition d) => d.Kernel != DecompositionKernel.ThinSvd;

    private static IReadOnlyList<Decomposition> Run(Action readout)
    {
        using DecompositionLog log = DecompositionLog.Begin();
        readout();
        return log.Entries;
    }
}
