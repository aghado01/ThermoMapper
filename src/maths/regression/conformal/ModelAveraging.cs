using System;
using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>
/// The model-averaged full step's sets (<see cref="ModelAveraging.FullStep"/>), with each candidate's model weights: for case t,
/// <paramref name="ModelWeights"/>[t][g, k] is q_k(y_g), the posterior probability of model k given the n units and the candidate.
/// </summary>
public sealed record AveragedPrediction(ConformalPrediction Prediction, IReadOnlyList<double[,]> ModelWeights);

/// <summary>
/// Conformal Bayesian model averaging (Bhagwat, Kong and Jiang, CBMA, arXiv 2511.16924; the conformal layer's slice 17): the full step's
/// score averaged over K models, σᵢ = Σ_k q_k(y)·σᵢ^{M_k}, each model's in-sample predictive density under its add-one-in weights weighed
/// by q_k(y) ∝ p(M_k | Z₁:ₙ)·p_{M_k}(y | x, Z₁:ₙ), the model's posterior probability given the n units times its predictive density at the
/// candidate. That is p(M_k | Z₁:ₙ₊₁), so σᵢ is the posterior predictive density under the hierarchical model that draws a model, then its
/// parameters (their lemma): symmetric in the n + 1 units, and CB-1's rule on it covers at least 1 − α (CB-1, their Theorem 1). Under the
/// right model among them the weights concentrate on it and the set approaches its own (their Theorem 2); under none, on the one nearest
/// in Kullback–Leibler divergence.
/// <para>Each model brings its full-step seams over the same units, strata and grids, and its log posterior probability given the n units,
/// ln p(M_k | Z₁:ₙ) = ln p(M_k) + ln p(Z₁:ₙ | M_k) up to a constant: its log evidence. Exact evidence keeps the score exactly symmetric;
/// an estimated one is a further error beside the draws' own, as the full step's Monte Carlo error is.</para>
/// </summary>
public static class ModelAveraging
{
    /// <summary>
    /// CB-1's sets on the averaged score at level α for the models' seams (<paramref name="models"/>) and their log posterior
    /// probabilities given the n units (<paramref name="logModelWeights"/>, normalized here): π(y) = (B + T + 1)/m on the averaged scores,
    /// accepted above k/m, k = ⌊α·m⌋.
    /// </summary>
    public static AveragedPrediction FullStep(IReadOnlyList<(PointwiseLogDensity Training, CandidateLogDensity Candidates)> models,
                                              IReadOnlyList<double> logModelWeights, double alpha = 0.1)
    {
        ArgumentNullException.ThrowIfNull(models);
        ArgumentNullException.ThrowIfNull(logModelWeights);
        SplitConformal.CheckAlpha(alpha);
        int k = models.Count;
        if (k == 0) throw new ArgumentException("At least one model is required.", nameof(models));
        if (logModelWeights.Count != k) throw new ArgumentException($"One log weight per model is required: {logModelWeights.Count} for {k} models.", nameof(logModelWeights));
        double top = double.NegativeInfinity;
        foreach (double w in logModelWeights)
        {
            if (double.IsNaN(w) || double.IsPositiveInfinity(w)) throw new ArgumentException("Each log model weight is a number below +∞.", nameof(logModelWeights));
            top = Math.Max(top, w);
        }
        if (double.IsNegativeInfinity(top)) throw new ArgumentException("Some model needs positive weight.", nameof(logModelWeights));
        CandidateLogDensity first = models[0].Candidates;
        int units = models[0].Training.Units.Units;
        for (int m = 1; m < k; m++)
        {
            CandidateLogDensity c = models[m].Candidates;
            if (models[m].Training.Units.Units != units || c.Cases != first.Cases || c.Values.GetLength(1) != first.Values.GetLength(1))
                throw new ArgumentException($"Model {m}'s seams cover other units or candidates than model 0's.", nameof(models));
            for (int t = 0; t < first.Cases; t++)
                if (c.StratumOf[t] != first.StratumOf[t] || !ConformalBayes.SameSequence(c.Grids[t], first.Grids[t]))
                    throw new ArgumentException($"Model {m}'s case {t} has another stratum or grid than model 0's.", nameof(models));
        }
        // Each model's draws and their normalized log weights, for the candidate's predictive density given the n units.
        var rows = new int[k][];
        var drawWeights = new double[]?[k];
        for (int m = 0; m < k; m++) (rows[m], drawWeights[m]) = PredictiveChecks.Rows(models[m].Training);

        var sets = new ConformalSet[first.Cases];
        var weights = new double[first.Cases][,];
        var logQ = new double[k];
        var predictive = new double[1];
        for (int t = 0; t < first.Cases; t++)
        {
            double[] grid = first.Grids[t];
            int offset = first.Offset(t), g = grid.Length;
            int[]? members = null;
            int[] below = new int[g], tied = new int[g];
            var p = new double[g];
            var own = new double[g];
            var q = new double[g, k];
            int m1 = 0, level = 0;
            for (int j = 0; j < g; j++)
            {
                int column = offset + j;
                var scores = new ConformalScores[k];
                double qTop = double.NegativeInfinity;
                for (int m = 0; m < k; m++)
                {
                    scores[m] = ConformalBayes.Scores(models[m].Training, models[m].Candidates, column);
                    if (members is null) members = scores[m].Units;
                    else if (!SameUnits(members, scores[m].Units))
                        throw new ArgumentException($"Model {m}'s stratum holds other units than model 0's.", nameof(models));
                    SplitConformal.LogPredictive(models[m].Candidates.Values, rows[m], drawWeights[m], column, 1, predictive);
                    logQ[m] = logModelWeights[m] + predictive[0];
                    qTop = Math.Max(qTop, logQ[m]);
                }
                // q_k ∝ p(M_k | Z₁:ₙ)·p_{M_k}(y | x, Z₁:ₙ); a candidate no model supports takes the models' prior weights.
                double qTotal = 0.0;
                for (int m = 0; m < k; m++)
                {
                    logQ[m] = double.IsNegativeInfinity(qTop) ? logModelWeights[m] - top : logQ[m] - qTop;
                    qTotal += Math.Exp(logQ[m]);
                }
                for (int m = 0; m < k; m++)
                {
                    logQ[m] -= Math.Log(qTotal);
                    q[j, m] = Math.Exp(logQ[m]);
                }
                int n = members!.Length;
                m1 = n + 1;
                level = (int)Math.Floor(alpha * m1);
                double candidate = Mixture(scores, logQ, -1);
                own[j] = candidate;
                for (int i = 0; i < n; i++)
                {
                    double unit = Mixture(scores, logQ, i);
                    if (unit < candidate) below[j]++;
                    else if (unit == candidate) tied[j]++;
                }
                p[j] = (below[j] + tied[j] + 1) / (double)m1;
            }
            weights[t] = q;
            sets[t] = level == 0
                ? new ConformalSet(t, first.StratumOf[t], m1, level, grid, null, null, null, null, null)
                : new ConformalSet(t, first.StratumOf[t], m1, level, grid, p, below, tied, own, null);
        }
        return new AveragedPrediction(new ConformalPrediction(sets, alpha), weights);
    }

    // ln Σ_k q_k σ^{M_k}: unit i's averaged score, or the candidate's at i = −1, by log-sum-exp over the models.
    private static double Mixture(ConformalScores[] scores, double[] logQ, int i)
    {
        double top = double.NegativeInfinity;
        for (int m = 0; m < scores.Length; m++) top = Math.Max(top, logQ[m] + (i < 0 ? scores[m].Candidate : scores[m].Scores[i]));
        if (double.IsNegativeInfinity(top)) return top;
        double sum = 0.0;
        for (int m = 0; m < scores.Length; m++) sum += Math.Exp(logQ[m] + (i < 0 ? scores[m].Candidate : scores[m].Scores[i]) - top);
        return top + Math.Log(sum);
    }

    private static bool SameUnits(int[] a, int[] b)
    {
        if (a.Length != b.Length) return false;
        for (int i = 0; i < a.Length; i++) if (a[i] != b[i]) return false;
        return true;
    }
}
