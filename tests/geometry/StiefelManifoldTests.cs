using System;
using System.Linq;
using Maths.Geometry;
using Maths.Geometry.Estimators.Intrinsic;
using Maths.LinAlg;
using Xunit;

namespace Maths.Geometry.Tests;

/// <summary>
/// <see cref="StiefelManifold"/> under the canonical metric, against references from outside its own code: the sphere at r = 1, the
/// orthogonal group at r = n, a turn of the frame within its span, <see cref="GrassmannManifold"/>'s distance as a lower bound met by
/// horizontal geodesics, and the first-variation identity grad ½d(·, Z)² = −log(Z) in the canonical metric. A seeded sweep over shapes and
/// lengths checks the round trip, and IRLS finds the centre of a symmetric cloud.
/// </summary>
public sealed class StiefelManifoldTests
{
    [Theory]
    [InlineData(1e-9)]
    [InlineData(0.4)]
    [InlineData(2.5)]
    public void Distance_SeesATurnOfTheFrame(double phi)
    {
        int n = 6, r = 3;
        var st = new StiefelManifold(n, r);
        double[] y = MatrixOps.RandomOrthonormal(n, r, 31);
        // Z = YQ, Q a turn by φ in the plane of Y's first two columns: ‖log Q‖_F = √2·φ, and the canonical metric halves its square.
        double[] z = (double[])y.Clone();
        for (int i = 0; i < n; i++)
        {
            z[i] = Math.Cos(phi) * y[i] + Math.Sin(phi) * y[i + n];
            z[i + n] = -Math.Sin(phi) * y[i] + Math.Cos(phi) * y[i + n];
        }

        Assert.Equal(phi, st.Distance(y, z), 1e-15);
        Assert.InRange(new GrassmannManifold(n, r).Distance(y, z), 0.0, 1e-15);
        Assert.InRange(st.Distance(y, y), 0.0, 1e-15);
    }

    [Theory]
    [InlineData(1e-10)]
    [InlineData(0.3)]
    [InlineData(2.0)]
    [InlineData(3.0)]
    public void FrameWidthOne_IsTheSphere(double theta)
    {
        int n = 5;
        var st = new StiefelManifold(n, 1);
        double[] basis = MatrixOps.RandomOrthonormal(n, 2, 41);
        double[] y = basis[..n], u = basis[n..];
        double[] z = y.Select((yi, i) => Math.Cos(theta) * yi + Math.Sin(theta) * u[i]).ToArray();

        Assert.Equal(theta, st.Distance(y, z), 1e-15 * Math.Max(1.0, 1.0 / (Math.PI - theta)));
        var log = new double[n];
        st.LogMap(y, z, log);
        for (int i = 0; i < n; i++) Assert.Equal(theta * u[i], log[i], 1e-14);
    }

    [Fact]
    public void FullWidth_IsTheOrthogonalGroup()
    {
        int n = 4;
        var st = new StiefelManifold(n, n);
        double[] y = MatrixOps.RandomOrthonormal(n, n, 51);
        double[] planes = MatrixOps.RandomOrthonormal(n, n, 52);
        double[] theta = { 2.2, -0.7 };
        var l = new double[n, n];
        for (int j = 0; j < 2; j++)
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++)
                    l[a, b] += theta[j] * (planes[2 * j * n + a] * planes[(2 * j + 1) * n + b]
                                           - planes[(2 * j + 1) * n + a] * planes[2 * j * n + b]);
        double[,] q = MatrixFunctions.SkewExponential(l);
        double[] z = Times(y, n, n, q);
        double[] expected = Times(y, n, n, l);

        // The bi-invariant distance ‖L‖_F, halved in square by the canonical metric: √(θ₁² + θ₂²).
        Assert.Equal(Math.Sqrt(theta[0] * theta[0] + theta[1] * theta[1]), st.Distance(y, z), 1e-13);
        var log = new double[n * n];
        st.LogMap(y, z, log);
        for (int i = 0; i < n * n; i++) Assert.Equal(expected[i], log[i], 1e-13);
    }

    [Fact]
    public void Sweep_RoundTripLengthAndSymmetry()
    {
        (int N, int R)[] shapes = { (3, 1), (7, 1), (4, 2), (5, 2), (3, 3), (5, 4), (6, 3), (8, 4) };
        double[] lengths = { 1e-10, 1e-3, 0.5, 1.5, 2.5 };
        var rng = new Random(20261113);
        double worstRoundTrip = 0.0, worstLength = 0.0, worstSymmetry = 0.0, worstOrth = 0.0;
        foreach (var (n, r) in shapes)
        {
            var st = new StiefelManifold(n, r);
            foreach (double length in lengths)
                for (int trial = 0; trial < 4; trial++)
                {
                    double[] y = MatrixOps.RandomOrthonormal(n, r, rng.Next(1, int.MaxValue));
                    double[] delta = Tangent(st, y, n, r, rng, length);
                    var z = new double[n * r];
                    st.ExpMap(y, delta, z);
                    worstOrth = Math.Max(worstOrth, OrthonormalityDefect(z, n, r));

                    var log = new double[n * r];
                    st.LogMap(y, z, log);
                    for (int i = 0; i < n * r; i++) worstRoundTrip = Math.Max(worstRoundTrip, Math.Abs(log[i] - delta[i]));
                    double d = st.Distance(y, z);
                    worstLength = Math.Max(worstLength, Math.Abs(d - length));
                    worstSymmetry = Math.Max(worstSymmetry, Math.Abs(d - st.Distance(z, y)));
                }
        }
        Assert.InRange(worstRoundTrip, 0.0, 1e-12);
        Assert.InRange(worstLength, 0.0, 1e-12);
        Assert.InRange(worstSymmetry, 0.0, 1e-12);
        Assert.InRange(worstOrth, 0.0, 1e-14);
    }

    [Fact]
    public void Exp_StartsAtYAlongTheTangent()
    {
        int n = 6, r = 2;
        var st = new StiefelManifold(n, r);
        double[] y = MatrixOps.RandomOrthonormal(n, r, 61);
        double[] delta = Tangent(st, y, n, r, new Random(62), 0.8);

        var z = new double[n * r];
        st.ExpMap(y, new double[n * r], z);
        Assert.Equal(y, z);

        double h = 1e-5;
        var plus = new double[n * r];
        var minus = new double[n * r];
        st.ExpMap(y, delta.Select(x => h * x).ToArray(), plus);
        st.ExpMap(y, delta.Select(x => -h * x).ToArray(), minus);
        for (int i = 0; i < n * r; i++) Assert.Equal(delta[i], (plus[i] - minus[i]) / (2 * h), 1e-9);
    }

    [Fact]
    public void Distance_IsInvariantUnderRotationsOfSpaceAndFrame()
    {
        int n = 6, r = 3;
        var st = new StiefelManifold(n, r);
        double[] y = MatrixOps.RandomOrthonormal(n, r, 71);
        var z = new double[n * r];
        st.ExpMap(y, Tangent(st, y, n, r, new Random(72), 1.2), z);
        double d = st.Distance(y, z);

        double[] u = MatrixOps.RandomOrthonormal(n, n, 73);
        Assert.Equal(d, st.Distance(LeftTimes(u, n, y, r), LeftTimes(u, n, z, r)), 1e-13);

        double[] frame = MatrixOps.RandomOrthonormal(r, r, 74);
        var rot = new double[r, r];
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++) rot[i, j] = frame[i + j * r];
        Assert.Equal(d, st.Distance(Times(y, n, r, rot), Times(z, n, r, rot)), 1e-13);
    }

    [Fact]
    public void Grassmann_BoundsBelow_AndHorizontalGeodesicsMeetIt()
    {
        int n = 7, r = 3;
        var st = new StiefelManifold(n, r);
        var gr = new GrassmannManifold(n, r);
        var rng = new Random(81);
        for (int trial = 0; trial < 20; trial++)
        {
            double[] y = MatrixOps.RandomOrthonormal(n, r, rng.Next(1, int.MaxValue));
            var z = new double[n * r];
            st.ExpMap(y, Tangent(st, y, n, r, rng, 0.2 + 2.0 * rng.NextDouble()), z);
            // The projection to the Grassmannian is a Riemannian submersion: it shortens no curve.
            Assert.True(st.Distance(y, z) >= gr.Distance(y, z) - 1e-13);

            // A horizontal tangent (YᵀΔ = 0) of norm below π/2 projects to a minimizing Grassmann geodesic of the same length.
            double[] horizontal = Tangent(st, y, n, r, rng, 1.0);
            double[,] yth = AtB(y, horizontal, n, r);
            for (int j = 0; j < r; j++)
                for (int l = 0; l < r; l++)
                    for (int i = 0; i < n; i++) horizontal[i + j * n] -= y[i + l * n] * yth[l, j];
            double length = 1.4 * rng.NextDouble() + 0.01;
            double scale = length / st.Norm(y, horizontal);
            for (int i = 0; i < n * r; i++) horizontal[i] *= scale;
            st.ExpMap(y, horizontal, z);
            Assert.Equal(length, st.Distance(y, z), 1e-12);
            Assert.Equal(length, gr.Distance(y, z), 1e-12);
        }
    }

    [Fact]
    public void Gradient_OfHalfSquaredDistance_IsMinusTheLog()
    {
        int n = 6, r = 3;
        var st = new StiefelManifold(n, r);
        var rng = new Random(91);
        for (int trial = 0; trial < 5; trial++)
        {
            double[] y = MatrixOps.RandomOrthonormal(n, r, rng.Next(1, int.MaxValue));
            var z = new double[n * r];
            st.ExpMap(y, Tangent(st, y, n, r, rng, 0.5 + rng.NextDouble()), z);
            double[] xi = Tangent(st, y, n, r, rng, 1.0);
            var log = new double[n * r];
            st.LogMap(y, z, log);

            double h = 1e-5;
            var plus = new double[n * r];
            var minus = new double[n * r];
            st.ExpMap(y, xi.Select(x => h * x).ToArray(), plus);
            st.ExpMap(y, xi.Select(x => -h * x).ToArray(), minus);
            double dPlus = st.Distance(plus, z), dMinus = st.Distance(minus, z);
            double derivative = (0.5 * dPlus * dPlus - 0.5 * dMinus * dMinus) / (2 * h);
            Assert.Equal(-CanonicalInner(y, log, xi, n, r), derivative, 1e-8);
        }
    }

    [Fact]
    public void Log_RefusesFramesOnTheCutLocus()
    {
        int n = 4;
        double[] y = MatrixOps.RandomOrthonormal(n, 2, 101);
        var sphere = new StiefelManifold(n, 1);
        Assert.Throws<InvalidOperationException>(() => sphere.Distance(y[..n], y[..n].Select(x => -x).ToArray()));

        // A reflection of one column within the span: det(YᵀZ) = −1.
        var st = new StiefelManifold(n, 2);
        double[] z = (double[])y.Clone();
        for (int i = 0; i < n; i++) z[i + n] = -y[i + n];
        Assert.Throws<InvalidOperationException>(() => st.Distance(y, z));
    }

    [Fact]
    public void Irls_FindsTheCentreOfASymmetricCloud()
    {
        int n = 5, r = 2;
        var st = new StiefelManifold(n, r);
        double[] centre = MatrixOps.RandomOrthonormal(n, r, 111);
        var rng = new Random(112);
        var cloud = new double[8][];
        for (int i = 0; i < 4; i++)
        {
            double[] delta = Tangent(st, centre, n, r, rng, 0.1 + 0.3 * rng.NextDouble());
            cloud[2 * i] = new double[n * r];
            cloud[2 * i + 1] = new double[n * r];
            st.ExpMap(centre, delta, cloud[2 * i]);
            st.ExpMap(centre, delta.Select(x => -x).ToArray(), cloud[2 * i + 1]);
        }
        double[] weights = Enumerable.Repeat(1.0, cloud.Length).ToArray();
        var start = new double[n * r];
        st.ExpMap(centre, Tangent(st, centre, n, r, rng, 0.05), start);

        // By symmetry the logs from the centre cancel, so it is the Karcher mean and the geometric median.
        double[] mean = (double[])start.Clone();
        GeometricMean.Compute(st, cloud, weights, mean);
        Assert.InRange(st.Distance(mean, centre), 0.0, 1e-9);

        double[] median = (double[])start.Clone();
        GeometricMedian.Compute(st, cloud, weights, median);
        Assert.InRange(st.Distance(median, centre), 0.0, 1e-8);
    }

    // A random tangent at Y (D − Y·sym(YᵀD)) of the given canonical length.
    private static double[] Tangent(StiefelManifold st, double[] y, int n, int r, Random rng, double length)
    {
        var d = new double[n * r];
        for (int i = 0; i < d.Length; i++) d[i] = 2.0 * rng.NextDouble() - 1.0;
        double[,] ytd = AtB(y, d, n, r);
        for (int j = 0; j < r; j++)
            for (int l = 0; l < r; l++)
            {
                double s = 0.5 * (ytd[l, j] + ytd[j, l]);
                for (int i = 0; i < n; i++) d[i + j * n] -= y[i + l * n] * s;
            }
        double scale = length / st.Norm(y, d);
        for (int i = 0; i < d.Length; i++) d[i] *= scale;
        return d;
    }

    // ⟨a, b⟩ = tr(aᵀb) − ½·tr((Yᵀa)ᵀ(Yᵀb)).
    private static double CanonicalInner(double[] y, double[] a, double[] b, int n, int r)
    {
        double s = 0.0;
        for (int i = 0; i < n * r; i++) s += a[i] * b[i];
        double[,] ya = AtB(y, a, n, r), yb = AtB(y, b, n, r);
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++) s -= 0.5 * ya[i, j] * yb[i, j];
        return s;
    }

    private static double[,] AtB(double[] a, double[] b, int n, int r)
    {
        var c = new double[r, r];
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++)
            {
                double s = 0.0;
                for (int l = 0; l < n; l++) s += a[l + i * n] * b[l + j * n];
                c[i, j] = s;
            }
        return c;
    }

    // Column-major n×k block times a k×k matrix.
    private static double[] Times(double[] y, int n, int k, double[,] q)
    {
        var z = new double[n * k];
        for (int j = 0; j < k; j++)
            for (int l = 0; l < k; l++)
                for (int i = 0; i < n; i++) z[i + j * n] += y[i + l * n] * q[l, j];
        return z;
    }

    // Column-major n×n U times column-major n×r Y.
    private static double[] LeftTimes(double[] u, int n, double[] y, int r)
    {
        var z = new double[n * r];
        for (int j = 0; j < r; j++)
            for (int l = 0; l < n; l++)
                for (int i = 0; i < n; i++) z[i + j * n] += u[i + l * n] * y[l + j * n];
        return z;
    }

    private static double OrthonormalityDefect(double[] z, int n, int r)
    {
        double[,] g = AtB(z, z, n, r);
        double worst = 0.0;
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++) worst = Math.Max(worst, Math.Abs(g[i, j] - (i == j ? 1.0 : 0.0)));
        return worst;
    }
}
