using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The dense product against its definition, the scalar loop of fused multiply-adds in inner order, bit for bit: chosen
/// shapes around the tile (4 rows, 8 columns) and the chunk (256 inner indices), a seeded sweep of shapes and magnitudes with
/// far ranges (subnormal and huge entries, signed zeros, exact zeros), equal columns giving equal columns, and a packed operand
/// reused and refilled.
/// </summary>
public sealed class DenseProductTests
{
    // The definition: c_ij = fma(a_i,K−1, b_K−1,j, … fma(a_i0, b_0j, +0)).
    private static double[] Reference(double[] a, int rows, int inner, double[] b, int columns)
    {
        var c = new double[rows * columns];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < inner; k++) sum = Math.FusedMultiplyAdd(a[i * inner + k], b[k * columns + j], sum);
                c[i * columns + j] = sum;
            }
        return c;
    }

    private static double[] Random(Xoshiro256PlusPlus rng, int length, bool far)
    {
        var x = new double[length];
        for (int i = 0; i < length; i++)
        {
            double u = rng.NextDouble();
            x[i] = !far ? 2 * rng.NextDouble() - 1
                 : u < 0.1 ? 0.0
                 : u < 0.15 ? -0.0
                 : u < 0.25 ? double.Epsilon * (1 + 1000 * rng.NextDouble())          // subnormal
                 : (rng.NextDouble() < 0.5 ? -1 : 1) * Math.Pow(2, 1000 * rng.NextDouble() - 500);
        }
        return x;
    }

    private static void AssertBitwise(double[] expected, double[] actual)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (int i = 0; i < expected.Length; i++)
            Assert.True(BitConverter.DoubleToInt64Bits(expected[i]) == BitConverter.DoubleToInt64Bits(actual[i]),
                        $"entry {i}: {expected[i]:R} expected, {actual[i]:R}");
    }

    [Theory]
    [InlineData(1, 1, 1)]
    [InlineData(4, 1, 8)]
    [InlineData(3, 5, 7)]
    [InlineData(5, 255, 9)]
    [InlineData(4, 256, 8)]
    [InlineData(9, 257, 17)]
    [InlineData(8, 600, 16)]
    [InlineData(13, 513, 3)]
    [InlineData(2, 0, 5)]
    public void Product_IsTheFusedLoop_BitForBit(int rows, int inner, int columns)
    {
        var rng = new Xoshiro256PlusPlus(seed: 100 * rows + 10 * columns + inner);
        double[] a = Random(rng, rows * inner, far: false), b = Random(rng, inner * columns, far: false);
        var c = new double[rows * columns];
        DenseProduct.Multiply(a, rows, DenseProduct.Pack(b, inner, columns), c);
        AssertBitwise(Reference(a, rows, inner, b, columns), c);
    }

    [Fact]
    public void SeededSweep_FarRanges_BitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 2026_11_06);
        for (int trial = 0; trial < 60; trial++)
        {
            int rows = 1 + (int)(rng.NextDouble() * 11), columns = 1 + (int)(rng.NextDouble() * 21);
            int inner = (int)(rng.NextDouble() * 700);
            bool far = trial % 2 == 1;
            double[] a = Random(rng, rows * inner, far), b = Random(rng, inner * columns, far);
            var c = new double[rows * columns];
            DenseProduct.Multiply(a, rows, DenseProduct.Pack(b, inner, columns), c);
            AssertBitwise(Reference(a, rows, inner, b, columns), c);
        }
    }

    [Fact]
    public void EqualColumns_GiveEqualColumns_AndFusedDotIsAnEntry()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7);
        int rows = 7, inner = 301, columns = 13;
        double[] a = Random(rng, rows * inner, far: false), b = Random(rng, inner * columns, far: false);
        for (int k = 0; k < inner; k++) { b[k * columns + 11] = b[k * columns + 2]; b[k * columns + 12] = b[k * columns + 2]; }
        var c = new double[rows * columns];
        DenseProduct.Multiply(a, rows, DenseProduct.Pack(b, inner, columns), c);
        var column = new double[inner];
        for (int k = 0; k < inner; k++) column[k] = b[k * columns + 2];
        for (int i = 0; i < rows; i++)
        {
            Assert.Equal(BitConverter.DoubleToInt64Bits(c[i * columns + 2]), BitConverter.DoubleToInt64Bits(c[i * columns + 11]));
            Assert.Equal(BitConverter.DoubleToInt64Bits(c[i * columns + 2]), BitConverter.DoubleToInt64Bits(c[i * columns + 12]));
            Assert.Equal(BitConverter.DoubleToInt64Bits(c[i * columns + 2]), BitConverter.DoubleToInt64Bits(DenseProduct.FusedDot(a.AsSpan(i * inner, inner), column)));
        }
    }

    [Fact]
    public void PackedOperand_ServesManyLeftOperands_AndRefills()
    {
        var rng = new Xoshiro256PlusPlus(seed: 11);
        int inner = 270, columns = 10;
        var packed = new PackedColumns(inner, columns);
        for (int round = 0; round < 3; round++)
        {
            double[] b = Random(rng, inner * columns, far: round == 2);
            packed.Fill(b);
            foreach (int rows in new[] { 1, 4, 6 })
            {
                double[] a = Random(rng, rows * inner, far: false);
                var c = new double[rows * columns];
                DenseProduct.Multiply(a, rows, packed, c);
                AssertBitwise(Reference(a, rows, inner, b, columns), c);
            }
        }
    }

    [Fact]
    public void TwoDimensionalOverload_MatchesAndRefusesMismatch()
    {
        var a = new double[,] { { 1, 2, 3 }, { 4, 5, 6 } };
        var b = new double[,] { { 1, 0 }, { 0, 1 }, { 1, 1 } };
        double[,] c = DenseProduct.Multiply(a, b);
        Assert.Equal(new double[,] { { 4, 5 }, { 10, 11 } }, c);
        Assert.Throws<ArgumentException>(() => DenseProduct.Multiply(a, new double[2, 2]));
        Assert.Equal(new double[2, 4], DenseProduct.Multiply(new double[2, 0], new double[0, 4]));
    }
}
