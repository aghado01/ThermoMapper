namespace GenerativePrior;

// Scalar row sums and the forward derivative along each coordinate.
// Only the sampled weights, code and observed data are shared with the descent path.
internal static class Reference
{
    private static double Dot(double[] a, double[] b)
    {
        double sum = 0;
        for (int i = 0; i < a.Length; i++) sum += a[i] * b[i];
        return sum;
    }

    private static double[] Apply(Matrix matrix, double[] x) =>
        matrix.Rows.Select(row => Dot(row, x)).ToArray();

    private static double Discrepancy(double actual, double expected) =>
        Math.Abs(actual - expected) / Math.Max(1, Math.Abs(expected));

    public static double Check(Network net, double[] code, double[] y, double objective,
                               double[]? gradient = null)
    {
        double[] h = Apply(net.W1, code).Select(v => Math.Max(0, v)).ToArray();
        double[] g = Apply(net.W2, h).Select(v => Math.Max(0, v)).ToArray();
        double[] r = Apply(net.A, g).Zip(y, (a, b) => a - b).ToArray();
        double discrepancy = Discrepancy(objective, 0.5 * Dot(r, r));
        if (gradient is null) return discrepancy;

        for (int j = 0; j < code.Length; j++)
        {
            double[] dh = net.W1.Rows.Select((row, i) => h[i] > 0 ? row[j] : 0).ToArray();
            double[] dg = Apply(net.W2, dh);
            for (int i = 0; i < g.Length; i++)
                if (g[i] == 0) dg[i] = 0;
            double component = Dot(r, Apply(net.A, dg));
            discrepancy = Math.Max(discrepancy, Discrepancy(gradient[j], component));
        }
        return discrepancy;
    }
}
