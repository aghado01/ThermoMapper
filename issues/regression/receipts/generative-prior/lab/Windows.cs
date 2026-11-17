namespace LatentRecoveryProbe;

// Periodic unit-energy Hann, using the pilot receipt's arithmetic.
internal static class Windows
{
    internal static double[] Hann(int w)
    {
        var g = new double[w];
        double energy = 0.0;
        for (int j = 0; j < w; j++)
        {
            g[j] = 0.5 * (1.0 - Math.Cos(2.0 * Math.PI * j / w));
            energy += g[j] * g[j];
        }
        double s = 1.0 / Math.Sqrt(energy);
        for (int j = 0; j < w; j++) g[j] *= s;
        return g;
    }
}
