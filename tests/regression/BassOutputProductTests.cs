using System;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The declared output-side inner product on a direct fit (F-X6): Q = I is the design's own arithmetic, bit for bit, whether
/// it is declared or left out; and the core, which holds the independent Gram, refuses a Q it has no centring for, with the
/// reason, while the working response checks that the Q it carries fits the response.
/// </summary>
public sealed class BassOutputProductTests
{
    private sealed class Diagonal(double[] d) : IInnerProduct
    {
        public bool IsIdentity => false;
        public int Dimension => d.Length;
        public double LogDeterminant => d.Sum(Math.Log);

        public double Dot(ReadOnlySpan<double> u, ReadOnlySpan<double> v)
        {
            double sum = 0.0;
            for (int i = 0; i < u.Length; i++) sum += u[i] * d[i % d.Length] * v[i];
            return sum;
        }
    }

    [Fact]
    public void ADeclaredIdentity_OnADirectFit_IsTheDesignsOwnArithmetic()
    {
        BassData data = SensitivityFixture.Functional();
        BassModel model = BassModel.Create(SensitivityFixture.Config(), data);
        var declared = new WorkingResponse(data.Response, output: IdentityInnerProduct.Instance);
        Assert.True(declared.Output.IsIdentity);
        Assert.Same(IdentityInnerProduct.Instance, new WorkingResponse(data.Response).Output);

        BassDesign a = model.NewDesign();
        var b = new BassDesign(model.Basis, data, model.Admission, working: declared);
        BassTerm[] terms =
        [
            new([new HingeFactor(0, 1, 2)]),
            new(functional: [new HingeFactor(0, -1, 3)]),
            new([new HingeFactor(1, -1, 4)], functional: [new HingeFactor(0, 1, 1)]),
            new([new HingeFactor(0, 1, 2), new HingeFactor(1, 1, 1)], functional: [new HingeFactor(0, -1, 5)]),
        ];
        Assert.Equal(a.ResponseSquares, b.ResponseSquares);
        foreach (BassTerm term in terms)
        {
            DesignProposal pa = a.ProposeBirth(term), pb = b.ProposeBirth(term);
            Assert.Equal(pa, pb);
            Assert.True(pa.Admitted);
            a.Accept();
            b.Accept();
        }
        Assert.Equal(a.ResidualSquares, b.ResidualSquares);
        Assert.Equal(a.FittedSquares, b.FittedSquares);
        Assert.True(a.Coefficients.SequenceEqual(b.Coefficients));
        Assert.True(a.Residual.SequenceEqual(b.Residual));
        a.Rebuild();
        b.Rebuild();
        Assert.Equal(a.ResidualSquares, b.ResidualSquares);
    }

    [Fact]
    public void ACorrelatedQ_IsCarriedByTheContract_AndRefusedByTheCore()
    {
        BassData functional = SensitivityFixture.Functional();
        int locations = functional.Response.Locations;
        var q = new Diagonal([.. Enumerable.Repeat(2.0, locations)]);
        var working = new WorkingResponse(functional.Response, output: q);
        Assert.Same(q, working.Output);

        BassModel model = BassModel.Create(SensitivityFixture.Config(), functional);
        var e = Assert.Throws<NotSupportedException>(() => new BassDesign(model.Basis, functional, model.Admission, working: working));
        Assert.Contains("F-X6", e.Message);

        // Q acts on the response's output index: a scalar response has none, and the length must match.
        Assert.Throws<ArgumentException>(() => new WorkingResponse(SensitivityFixture.Scalar().Response, output: new Diagonal([1.0])));
        Assert.Throws<ArgumentException>(() => new WorkingResponse(functional.Response, output: new Diagonal([.. Enumerable.Repeat(2.0, locations + 1)])));
    }
}
