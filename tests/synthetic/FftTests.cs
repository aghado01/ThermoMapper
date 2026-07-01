using System;
using System.Numerics;
using Maths.LinAlg;
using Xunit;

namespace Synthetic.Tests;

public sealed class FftTests
{
    [Fact]
    public void Forward_ThenInverse_RecoversInput()
    {
        var rng = new Random(1);
        var data = new Complex[8];
        var original = new Complex[8];
        for (int i = 0; i < 8; i++)
        {
            data[i] = new Complex(rng.NextDouble(), rng.NextDouble());
            original[i] = data[i];
        }

        Fft.Forward(data);
        Fft.Inverse(data);

        for (int i = 0; i < 8; i++)
        {
            Assert.Equal(original[i].Real, data[i].Real, 10);
            Assert.Equal(original[i].Imaginary, data[i].Imaginary, 10);
        }
    }

    [Fact]
    public void Forward_UnitImpulse_GivesFlatSpectrum()
    {
        // DFT of a unit impulse at index 0 is all-ones.
        var data = new Complex[8];
        data[0] = Complex.One;

        Fft.Forward(data);

        for (int i = 0; i < 8; i++)
        {
            Assert.Equal(1.0, data[i].Real, 10);
            Assert.Equal(0.0, data[i].Imaginary, 10);
        }
    }

    [Fact]
    public void Forward_PureCosine_ConcentratesInTwoBins()
    {
        // cos(2π·k0·n/N) has spectral mass N/2 only at bins k0 and N−k0.
        const int n = 16, k0 = 3;
        var data = new Complex[n];
        for (int i = 0; i < n; i++)
            data[i] = new Complex(Math.Cos(2.0 * Math.PI * k0 * i / n), 0.0);

        Fft.Forward(data);

        for (int i = 0; i < n; i++)
        {
            double mag = data[i].Magnitude;
            if (i == k0 || i == n - k0)
                Assert.Equal(n / 2.0, mag, 8);
            else
                Assert.True(mag < 1e-8, $"bin {i} should be ~0 but was {mag}");
        }
    }

    [Fact]
    public void Forward3D_ThenInverse3D_RecoversInput()
    {
        const int n = 4, nn = n * n * n;
        var rng = new Random(7);
        var data = new Complex[nn];
        var original = new Complex[nn];
        for (int i = 0; i < nn; i++)
        {
            data[i] = new Complex(rng.NextDouble(), 0.0);
            original[i] = data[i];
        }

        Fft.Forward3D(data, n, n, n);
        Fft.Inverse3D(data, n, n, n);

        for (int i = 0; i < nn; i++)
        {
            Assert.Equal(original[i].Real, data[i].Real, 10);
            Assert.Equal(original[i].Imaginary, data[i].Imaginary, 10);
        }
    }

    [Fact]
    public void ForwardND_ThenInverseND_RecoversInput_4D()
    {
        int[] dims = { 4, 2, 4, 2 }; // mixed power-of-two axis lengths
        int total = 4 * 2 * 4 * 2;
        var rng = new Random(13);
        var data = new Complex[total];
        var original = new Complex[total];
        for (int i = 0; i < total; i++)
        {
            data[i] = new Complex(rng.NextDouble(), rng.NextDouble());
            original[i] = data[i];
        }

        Fft.ForwardND(data, dims);
        Fft.InverseND(data, dims);

        for (int i = 0; i < total; i++)
        {
            Assert.Equal(original[i].Real, data[i].Real, 10);
            Assert.Equal(original[i].Imaginary, data[i].Imaginary, 10);
        }
    }

    [Fact]
    public void Forward_NonPowerOfTwoLength_Throws()
    {
        var data = new Complex[6];
        Assert.Throws<ArgumentException>(() => Fft.Forward(data));
    }
}
