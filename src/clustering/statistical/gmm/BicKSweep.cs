using System;

namespace Clustering.Statistical.GMM
{
    /// <summary>
    /// K-selection by BIC: fits a GMM for each K in a range and returns one
    /// <see cref="KSweepResult"/> per K, with the BIC score for ranking. Lower is
    /// better. See docs/gmm.md for K-selection guidance and the planned
    /// algorithm-agnostic criteria (AIC, silhouette).
    /// </summary>
    public static class BicKSweep
    {
        /// <summary>
        /// Fits a GMM for each K ∈ [<paramref name="kMin"/>, <paramref name="kMax"/>]
        /// with <paramref name="restarts"/> random restarts per K, keeping the
        /// highest-LL fit. Pass <paramref name="randomSeed"/> for reproducibility.
        /// </summary>
        public static KSweepResult[] Run(
            double[][] data,
            int dimension,
            int kMin = 1,
            int kMax = 10,
            int maxIterations = 100,
            double tolerance = 1e-4,
            int restarts = 1,
            int? randomSeed = null)
        {
            if (kMin < 1) throw new ArgumentOutOfRangeException(nameof(kMin), "Must be ≥ 1.");
            if (kMax < kMin) throw new ArgumentOutOfRangeException(nameof(kMax), "Must be ≥ kMin.");
            if (restarts < 1) throw new ArgumentOutOfRangeException(nameof(restarts), "Must be ≥ 1.");

            int n = data.Length;
            double logN = Math.Log(n);

            // Free parameters for a K-component full-covariance GMM in D dimensions:
            // K·D means + K·D(D+1)/2 covariances + (K-1) weights.
            static int FreeParams(int k, int d)
                => k * d
                 + k * (d * (d + 1) / 2)
                 + (k - 1);

            var results = new KSweepResult[kMax - kMin + 1];

            for (int ki = 0; ki < results.Length; ki++)
            {
                int k = kMin + ki;
                int p = FreeParams(k, dimension);

                GaussianMixtureModel? bestModel = null;
                double bestLogL = double.NegativeInfinity;

                int baseSeed = randomSeed ?? Environment.TickCount;
                for (int restart = 0; restart < restarts; restart++)
                {
                    var model = new GaussianMixtureModel(k, dimension);
                    model.RandomInitialize(data, new Random(baseSeed + restart * 7919));
                    model.Fit(data, maxIterations, tolerance);

                    if (model.FinalLogLikelihood > bestLogL)
                    {
                        bestLogL = model.FinalLogLikelihood;
                        bestModel = model;
                    }
                }

                double bic = -2.0 * bestLogL + p * logN;

                results[ki] = new KSweepResult(
                    K: k,
                    Bic: bic,
                    LogLikelihood: bestLogL,
                    NumIterations: bestModel!.NumIterations,
                    IsConverged: bestModel.IsConverged,
                    Model: bestModel);
            }

            return results;
        }

        /// <summary>Returns the result with the lowest BIC.</summary>
        public static KSweepResult BestByBic(KSweepResult[] results)
        {
            KSweepResult best = results[0];
            for (int i = 1; i < results.Length; i++)
                if (results[i].Bic < best.Bic) best = results[i];
            return best;
        }

    }
}
