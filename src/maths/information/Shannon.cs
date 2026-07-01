using System;

namespace Maths.Information;

public static class Shannon
{
    public static double EntropyNats(ReadOnlySpan<int> counts)
    {
        long total = 0;
        for (int i = 0; i < counts.Length; i++)
        {
            if (counts[i] > 0)
                total += counts[i];
        }

        if (total <= 0)
            return 0.0;

        double invTotal = 1.0 / total;
        double h = 0.0;

        for (int i = 0; i < counts.Length; i++)
        {
            if (counts[i] <= 0)
                continue;

            double p = counts[i] * invTotal;
            h -= p * Math.Log(p);
        }

        return h;
    }

    public static double EntropyNats(ReadOnlySpan<double> probabilities)
    {
        double sum = 0.0;
        for (int i = 0; i < probabilities.Length; i++)
        {
            if (probabilities[i] > 0.0)
                sum += probabilities[i];
        }

        if (sum <= 0.0)
            return 0.0;

        double invSum = 1.0 / sum;
        double h = 0.0;

        for (int i = 0; i < probabilities.Length; i++)
        {
            if (probabilities[i] <= 0.0)
                continue;

            double p = probabilities[i] * invSum;
            h -= p * Math.Log(p);
        }

        return h;
    }
}
