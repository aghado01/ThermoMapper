using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>Independent dense Gaussian references for E21's crossed and nested intercept blocks.</summary>
public sealed class RandomEffectsWhiteningTests(ITestOutputHelper output)
{
    internal sealed record Fixture(CategoricalPredictor[] Groups, CategoricalPredictor?[] Parents, double[] Scales)
    {
        public RandomEffectGeometry Geometry() => new([.. Groups.Select((g, k) => Parents[k] is { } p
            ? RandomEffectGrouping.WithinParent(g, p) : RandomEffectGrouping.ObservationCentred(g))]);
    }

    private static CategoricalPredictor Column(string name, int[] codes, int? declared = null)
        => new(name, codes, [.. Enumerable.Range(0, declared ?? codes.Max() + 1).Select(g => $"L{g}")]);

    internal static Fixture Example(string name)
    {
        if (name == "crossed")
            return new Fixture(
                [Column("school", [0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2]),
                 Column("reader", [0,1,2,0,1,0,1,2,2,3,3,0,0,1,2,3,3])],
                [null, null], [0.7, 1.4]);
        if (name == "nested")
        {
            var parent = Column("parent", [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1]);
            var child = Column("child", [0,0,1,1,1,2,2,2,3,4,4,4,4,5,5,5,5,5]);
            return new Fixture([parent, child], [null, parent], [0.9, 0.6]);
        }
        if (name == "alias")
            return new Fixture([Column("first", [0,0,0,1,1,2,2,2,2], declared: 4),
                                Column("second", [0,0,0,1,1,2,2,2,2])], [null, null], [0.4, 1.7]);
        if (name == "mixed alias")
        {
            Fixture nested = Example("nested");
            var crossed = Column("reader", nested.Groups[1].Codes.ToArray());
            return new Fixture([.. nested.Groups, crossed], [.. nested.Parents, null], [.. nested.Scales, 1.7]);
        }
        var grouping = Column("one", [0,0,0,0,1,1,2,2,2,3,3,3,3,3]);
        return new Fixture([grouping], [null], [1.3]);
    }

    internal sealed record Dense(double[,] Indicators, double[,] Prior, double[,] Covariance, int[] Offsets);

    // The dense oracle constructs each level-space projector directly from the count constraint. It does not use
    // the production contrast basis, its thin factorization, or the retained observation-space rank.
    internal static Dense Reference(Fixture fixture)
    {
        int n = fixture.Groups[0].Length;
        int[][] observed = [.. fixture.Groups.Select(g => g.Codes.ToArray().Distinct().Order().ToArray())];
        int[] offsets = new int[fixture.Groups.Length + 1];
        for (int k = 0; k < observed.Length; k++) offsets[k + 1] = offsets[k] + observed[k].Length;
        var indicators = new double[n, offsets[^1]];
        var prior = new double[offsets[^1], offsets[^1]];
        for (int k = 0; k < fixture.Groups.Length; k++)
        {
            int[] codes = fixture.Groups[k].Codes.ToArray();
            int[] parent = fixture.Parents[k]?.Codes.ToArray() ?? new int[n];
            double[] counts = [.. observed[k].Select(g => (double)codes.Count(c => c == g))];
            int[] parentOfLevel = [.. observed[k].Select(g => parent[Array.FindIndex(codes, c => c == g)])];
            for (int i = 0; i < n; i++) indicators[i, offsets[k] + Array.IndexOf(observed[k], codes[i])] = 1;
            for (int g = 0; g < observed[k].Length; g++)
                for (int h = 0; h < observed[k].Length; h++)
                {
                    double denominator = Enumerable.Range(0, observed[k].Length)
                        .Where(j => parentOfLevel[j] == parentOfLevel[g]).Sum(j => counts[j] * counts[j]);
                    double projector = (g == h ? 1 : 0) - (parentOfLevel[g] == parentOfLevel[h] ? counts[g] * counts[h] / denominator : 0);
                    prior[offsets[k] + g, offsets[k] + h] = fixture.Scales[k] * fixture.Scales[k] * projector;
                }
        }
        double[,] covariance = Multiply(Multiply(indicators, prior), Transpose(indicators));
        for (int i = 0; i < n; i++) covariance[i, i] += 1;
        return new Dense(indicators, prior, covariance, offsets);
    }

    [Theory]
    [InlineData("one")]
    [InlineData("crossed")]
    [InlineData("nested")]
    [InlineData("alias")]
    [InlineData("mixed alias")]
    public void JointWhitening_AgreesWithTheFullDenseCovariance(string name)
    {
        Fixture fixture = Example(name);
        RandomEffectGeometry geometry = fixture.Geometry();
        RandomEffectWhitening whitening = geometry.At(fixture.Scales);
        Dense dense = Reference(fixture);
        double[,] actual = Matrix(geometry.Cases, whitening.Whiten);
        EigenResult reference = Eigen.DecomposeSymmetric(dense.Covariance, tol: 1e-15);
        var expected = new double[geometry.Cases, geometry.Cases];
        for (int i = 0; i < geometry.Cases; i++)
            for (int j = 0; j < geometry.Cases; j++)
                for (int k = 0; k < geometry.Cases; k++)
                    expected[i, j] += reference.Eigenvectors[k][i] * reference.Eigenvectors[k][j] / Math.Sqrt(reference.Eigenvalues[k]);
        double whiteGap = Difference(actual, expected);
        double inverseGap = Difference(Multiply(Multiply(actual, dense.Covariance), actual), Identity(geometry.Cases));
        double symmetry = Difference(actual, Transpose(actual));
        var factor = new DenseCholesky(geometry.Cases);
        factor.Decompose(dense.Covariance);
        double[] ones = [.. Enumerable.Repeat(1.0, geometry.Cases)];
        whitening.Whiten(ones);
        output.WriteLine($"{name}: n={geometry.Cases}, q={geometry.Dimension}, r={geometry.Rank}; full inverse square root gap {whiteGap:E3}, " +
                         $"WVW-I {inverseGap:E3}, symmetry {symmetry:E3}, log determinant gap {whitening.LogDeterminant - factor.LogDet:E3}");
        Assert.InRange(whiteGap, 0, 1e-12);
        Assert.InRange(inverseGap, 0, 1e-12);
        Assert.InRange(symmetry, 0, 1e-12);
        Assert.InRange(Math.Abs(whitening.LogDeterminant - factor.LogDet), 0, 1e-12);
        Assert.All(ones, x => Assert.InRange(Math.Abs(x - 1), 0, 1e-12));
        if (name == "alias")
        {
            Assert.Equal(2, geometry.Rank);
            Assert.Equal(4, geometry.Dimension);
            Assert.Equal(3, geometry.Aliases.Count);
            Assert.Equal(-1, geometry.Groups[0].ObservedOf(3));
        }
        if (name == "mixed alias") Assert.Equal(6, geometry.Aliases.Count);
    }

    [Theory]
    [InlineData("crossed")]
    [InlineData("nested")]
    [InlineData("alias")]
    [InlineData("mixed alias")]
    public void JointEffects_MatchTheDenseConditionalCovariance(string name)
    {
        Fixture fixture = Example(name);
        RandomEffectGeometry geometry = fixture.Geometry();
        RandomEffectWhitening whitening = geometry.At(fixture.Scales);
        Dense dense = Reference(fixture);
        double[,] inverse = Inverse(dense.Covariance);
        double[,] gain = Multiply(Multiply(dense.Prior, Transpose(dense.Indicators)), inverse);
        const double sigma2 = 0.7;
        double[] residual = [.. Enumerable.Range(0, geometry.Cases).Select(i => Math.Sin(0.4 * i) + 0.1 * i)];
        double[] mean = Product(gain, residual);
        double[,] reduction = Multiply(Multiply(gain, dense.Indicators), dense.Prior);
        var covariance = new double[mean.Length, mean.Length];
        for (int g = 0; g < mean.Length; g++)
            for (int h = 0; h < mean.Length; h++) covariance[g, h] = sigma2 * (dense.Prior[g, h] - reduction[g, h]);
        double[] actualMean = [.. whitening.MeanEffects(residual).SelectMany(x => x)];
        double gap = actualMean.Zip(mean, (a, b) => Math.Abs(a - b)).Max();
        Assert.InRange(gap, 0, 1e-12);
        const int draws = 50_000;
        var first = new double[mean.Length];
        var second = new double[mean.Length, mean.Length];
        int seed = name == "mixed alias" ? 210335 : SeedTree.Derive(20261130, 3)[Array.IndexOf(new[] { "crossed", "nested", "alias" }, name)];
        var rng = new Xoshiro256PlusPlus(seed);
        double worstConstraint = 0;
        for (int j = 0; j < draws; j++)
        {
            double[][] effects = whitening.DrawEffects(residual, sigma2, rng);
            double[] u = [.. effects.SelectMany(x => x)];
            worstConstraint = Math.Max(worstConstraint, Constraint(fixture, effects));
            for (int g = 0; g < u.Length; g++)
            {
                first[g] += u[g] / draws;
                for (int h = 0; h < u.Length; h++) second[g, h] += (u[g] - mean[g]) * (u[h] - mean[h]) / draws;
            }
        }
        double worstMean = 0, worstCovariance = 0;
        for (int g = 0; g < mean.Length; g++)
        {
            worstMean = Math.Max(worstMean, Math.Abs(first[g] - mean[g]) / Math.Sqrt(covariance[g, g] / draws));
            for (int h = 0; h < mean.Length; h++)
            {
                double se = Math.Sqrt((covariance[g, g] * covariance[h, h] + covariance[g, h] * covariance[g, h]) / draws);
                worstCovariance = Math.Max(worstCovariance, Math.Abs(second[g, h] - covariance[g, h]) / se);
            }
        }
        output.WriteLine($"{name}: exact mean gap {gap:E3}; {draws} joint draws, worst mean z {worstMean:F3}, covariance z {worstCovariance:F3}, " +
                         $"largest centred contribution {worstConstraint:E3}");
        Assert.InRange(worstMean, 0, 6);
        Assert.InRange(worstCovariance, 0, 6);
        Assert.InRange(worstConstraint, 0, 1e-12);
    }

    [Theory]
    [InlineData(0.05)]
    [InlineData(0.4)]
    [InlineData(1.3)]
    [InlineData(4.0)]
    public void OneGrouping_MatchesTheExistingClosedForm(double scale)
    {
        Fixture fixture = Example("one");
        RandomEffectWhitening general = fixture.Geometry().At([scale]);
        RandomInterceptBlock closed = RandomInterceptBlock.Of(fixture.Groups[0], new BassRandomIntercept(fixture.Groups[0].Name, scale));
        double[,] a = Matrix(fixture.Groups[0].Length, general.Whiten), b = Matrix(fixture.Groups[0].Length, closed.Whiten);
        double gap = Difference(a, b);
        output.WriteLine($"s={scale}: general versus closed form {gap:E3}, log determinant gap {general.LogDeterminant - closed.LogDeterminant:E3}");
        Assert.InRange(gap, 0, 1e-12);
        Assert.InRange(Math.Abs(general.LogDeterminant - closed.LogDeterminant), 0, 1e-12);
    }

    [Theory]
    [InlineData(1e-6, 0.6)]
    [InlineData(1e6, 0.6)]
    [InlineData(0.9, 1e6)]
    public void UnequalCrossedScales_KeepTheSmallComponent(double firstScale, double secondScale)
    {
        int[] first = [0,0,0,1,1,1,1], second = [0,1,1,0,1,0,1];
        var geometry = new RandomEffectGeometry([RandomEffectGrouping.ObservationCentred(Column("first", first)),
            RandomEffectGrouping.ObservationCentred(Column("second", second))]);
        var at = geometry.At([firstScale, secondScale]);
        // Independent two-contrast reference. Form the small eigenvalue from det(C)/λ_max,
        // using the triangular factor's positive product instead of subtracting two large numbers.
        double[] ColumnOf(int[] codes)
        {
            double zero = codes.Count(g => g == 0), one = codes.Length - zero, norm = Math.Sqrt(zero*zero + one*one);
            return [.. codes.Select(g => g == 0 ? one/norm : -zero/norm)];
        }
        double[] a = ColumnOf(first), b = ColumnOf(second);
        double aa = a.Sum(v => v*v), ab = a.Zip(b, (x,y) => x*y).Sum(), bb = b.Sum(v => v*v);
        double normA = Math.Sqrt(aa), cross = ab/normA, normB = Math.Sqrt(bb-cross*cross);
        double[] q0 = [.. a.Select(v => v/normA)], q1 = [.. b.Select((v,i) => (v-cross*q0[i])/normB)];
        double c00 = firstScale*firstScale*aa + secondScale*secondScale*cross*cross;
        double c01 = secondScale*secondScale*cross*normB, c11 = secondScale*secondScale*normB*normB;
        double large = 0.5*(c00+c11+Math.Sqrt((c00-c11)*(c00-c11)+4*c01*c01));
        double small = Math.Pow(firstScale*secondScale*normA*normB, 2)/large;
        double angle = 0.5*Math.Atan2(2*c01,c00-c11), cosine = Math.Cos(angle), sine = Math.Sin(angle);
        double[] e0 = [.. q0.Select((v,i) => cosine*v+sine*q1[i])], e1 = [.. q0.Select((v,i) => -sine*v+cosine*q1[i])];
        double[] residual = [.. Enumerable.Range(0, first.Length).Select(i => Math.Sin(0.4*i)+0.1*i)];
        double dot0 = e0.Zip(residual, (x,y) => x*y).Sum(), dot1 = e1.Zip(residual, (x,y) => x*y).Sum();
        double[] expected = [.. residual.Select((v,i) => v + e0[i]*dot0*(1/Math.Sqrt(1+large)-1) + e1[i]*dot1*(1/Math.Sqrt(1+small)-1))];
        double[] actual = (double[])residual.Clone();
        at.Whiten(actual);
        double gap = actual.Zip(expected, (x,y) => Math.Abs(x-y)).Max();
        double p00 = aa+1/(firstScale*firstScale), p11 = bb+1/(secondScale*secondScale), det = p00*p11-ab*ab;
        double ar = a.Zip(residual, (x,y) => x*y).Sum(), br = b.Zip(residual, (x,y) => x*y).Sum();
        double[] mean = [(p11*ar-ab*br)/det, (p00*br-ab*ar)/det];
        double[][] effects = at.MeanEffects(residual);
        for (int i = 0; i < residual.Length; i++)
        {
            gap = Math.Max(gap, Math.Abs(effects[0][first[i]] - a[i]*mean[0]));
            gap = Math.Max(gap, Math.Abs(effects[1][second[i]] - b[i]*mean[1]));
        }
        const double sigma2 = 0.7;
        double[,] covariance = {{sigma2*p11/det, -sigma2*ab/det}, {-sigma2*ab/det, sigma2*p00/det}};
        var actualCovariance = new double[2,2];
        for (int k = 0; k < 2; k++)
        {
            double[] e = new double[2]; e[k] = 1;
            double[][] u = geometry.Effects(at.ConditionalNoise(e));
            double[] gamma = [u[0][0]/a[0], u[1][0]/b[0]];
            for (int g = 0; g < 2; g++)
                for (int h = 0; h < 2; h++) actualCovariance[g,h] += sigma2*gamma[g]*gamma[h];
        }
        double covarianceGap = Difference(covariance, actualCovariance);
        const int draws = 50_000;
        var rng = new Xoshiro256PlusPlus(210336);
        double[] sampleMean = new double[2]; var sampleCovariance = new double[2,2];
        for (int j = 0; j < draws; j++)
        {
            double[][] u = at.DrawEffects(residual, sigma2, rng);
            double[] gamma = [u[0][0]/a[0], u[1][0]/b[0]];
            for (int g = 0; g < 2; g++)
            {
                sampleMean[g] += gamma[g]/draws;
                for (int h = 0; h < 2; h++) sampleCovariance[g,h] += (gamma[g]-mean[g])*(gamma[h]-mean[h])/draws;
            }
        }
        double worstZ = 0;
        for (int g = 0; g < 2; g++)
        {
            worstZ = Math.Max(worstZ, Math.Abs(sampleMean[g]-mean[g])/Math.Sqrt(covariance[g,g]/draws));
            for (int h = 0; h < 2; h++)
            {
                double se = Math.Sqrt((covariance[g,g]*covariance[h,h]+covariance[g,h]*covariance[g,h])/draws);
                worstZ = Math.Max(worstZ, Math.Abs(sampleCovariance[g,h]-covariance[g,h])/se);
            }
        }
        output.WriteLine($"scales {firstScale:G3}, {secondScale:G3}: whitening/mean gap {gap:E3}, covariance gap {covarianceGap:E3}; {draws} draws, max mean/covariance z {worstZ:F3}");
        Assert.InRange(gap, 0, 1e-9);
        Assert.InRange(covarianceGap, 0, 1e-9);
        Assert.InRange(worstZ, 0, 6);
    }

    [Fact]
    public void CrossedWhitenings_NeedTheJointCovariance_AndNestedSpansAreOrthogonal()
    {
        foreach (string name in new[] { "crossed", "nested" })
        {
            Fixture fixture = Example(name);
            RandomEffectGeometry geometry = fixture.Geometry();
            var first = new RandomEffectGeometry([geometry.Groups[0]]).At([fixture.Scales[0]]);
            var second = new RandomEffectGeometry([geometry.Groups[1]]).At([fixture.Scales[1]]);
            double[,] sequential = Multiply(Matrix(geometry.Cases, first.Whiten), Matrix(geometry.Cases, second.Whiten));
            double gap = Difference(Multiply(Multiply(sequential, Reference(fixture).Covariance), Transpose(sequential)), Identity(geometry.Cases));
            output.WriteLine($"{name}: sequential whitening covariance gap {gap:E3}");
            if (name == "nested") Assert.InRange(gap, 0, 1e-12);
            else Assert.True(gap > 1e-3);
        }
    }

    [Fact]
    public void WithinParent_RefusesCrossingLevels_AndCanHaveNoContrast()
    {
        var parent = Column("parent", [0,0,1,1]);
        var crossed = Column("child", [0,1,0,1]);
        Assert.Contains("more than one", Assert.Throws<ArgumentException>(() => RandomEffectGrouping.WithinParent(crossed, parent)).Message);
        var singleChild = Column("child", [0,0,1,1]);
        RandomEffectGrouping group = RandomEffectGrouping.WithinParent(singleChild, parent);
        var geometry = new RandomEffectGeometry([group]);
        RandomEffectWhitening whitening = geometry.At([2.0]);
        Assert.Equal(0, geometry.Dimension);
        Assert.Equal(0, geometry.Rank);
        double[] x = [1,2,3,4], original = [.. x];
        whitening.Whiten(x);
        Assert.Equal(original, x);
        Assert.Equal(0.0, whitening.LogDeterminant);
        Assert.All(whitening.DrawEffects(x, 0.7, new Xoshiro256PlusPlus(1))[0], u => Assert.Equal(0, u));
    }

    private static double Constraint(Fixture fixture, double[][] effects)
    {
        double worst = 0;
        for (int k = 0; k < effects.Length; k++)
        {
            int[] codes = fixture.Groups[k].Codes.ToArray();
            int[] parent = fixture.Parents[k]?.Codes.ToArray() ?? new int[codes.Length];
            int[] observed = [.. codes.Distinct().Order()];
            foreach (int p in parent.Distinct())
            {
                double sum = 0;
                for (int i = 0; i < codes.Length; i++)
                    if (parent[i] == p) sum += effects[k][Array.IndexOf(observed, codes[i])];
                worst = Math.Max(worst, Math.Abs(sum));
            }
        }
        return worst;
    }

    private delegate void Transform(Span<double> x);
    private static double[,] Matrix(int n, Transform transform)
    {
        var matrix = new double[n, n];
        for (int j = 0; j < n; j++)
        {
            var column = new double[n];
            column[j] = 1;
            transform(column);
            for (int i = 0; i < n; i++) matrix[i, j] = column[i];
        }
        return matrix;
    }

    private static double[,] Identity(int n)
    {
        var matrix = new double[n, n];
        for (int i = 0; i < n; i++) matrix[i, i] = 1;
        return matrix;
    }

    private static double[,] Transpose(double[,] a)
    {
        var result = new double[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) result[j, i] = a[i, j];
        return result;
    }

    private static double[,] Multiply(double[,] a, double[,] b)
    {
        var result = new double[a.GetLength(0), b.GetLength(1)];
        for (int i = 0; i < result.GetLength(0); i++)
            for (int j = 0; j < result.GetLength(1); j++)
                for (int k = 0; k < a.GetLength(1); k++) result[i, j] += a[i, k] * b[k, j];
        return result;
    }

    private static double[] Product(double[,] matrix, double[] x)
    {
        var result = new double[matrix.GetLength(0)];
        for (int i = 0; i < result.Length; i++)
            for (int j = 0; j < x.Length; j++) result[i] += matrix[i, j] * x[j];
        return result;
    }

    private static double[,] Inverse(double[,] matrix)
    {
        var factor = new DenseCholesky(matrix.GetLength(0));
        factor.Decompose(matrix);
        Assert.Equal(-1, factor.FailedAt);
        var result = new double[matrix.GetLength(0), matrix.GetLength(1)];
        for (int j = 0; j < result.GetLength(1); j++)
        {
            var basis = new double[result.GetLength(0)];
            basis[j] = 1;
            double[] solved = factor.Solve(basis);
            for (int i = 0; i < solved.Length; i++) result[i, j] = solved[i];
        }
        return result;
    }

    private static double Difference(double[,] a, double[,] b)
    {
        double worst = 0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) worst = Math.Max(worst, Math.Abs(a[i, j] - b[i, j]));
        return worst;
    }
}
