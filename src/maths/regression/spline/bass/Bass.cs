using System;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS, Bayesian adaptive spline surfaces (Francom and Sansó 2020): the entry points. The fit and its resumption join
/// in stage 4's S7; the run's fingerprint is here from S1, since the contracts define it.
/// </summary>
public static class Bass
{
    /// <summary>
    /// The run's fingerprint, which a checkpoint must match to resume: FNV-1a (<see cref="RunFingerprint"/>) over the
    /// config's canonical JSON (<see cref="BassConfig.ToJson"/>) and the data's exact bits, predictors in order with
    /// their kinds and names, then the response and its output variables. The config is all values, so the
    /// fingerprint covers every setting, not only the components' types.
    /// </summary>
    public static string Fingerprint(BassConfig config, BassData data)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(data);
        RunFingerprint f = new RunFingerprint().Add("BASS").Add(config.ToJson()).Add(data.Predictors.Count);
        foreach (BassPredictor p in data.Predictors)
            switch (p)
            {
                case ContinuousPredictor c:
                    f.Add("continuous").Add(c.Name).Add(c.Values.Span);
                    break;
                case CategoricalPredictor k:
                    f.Add("categorical").Add(k.Name).Add(k.Levels.Count);
                    foreach (string level in k.Levels) f.Add(level);
                    foreach (int code in k.Codes.Span) f.Add(code);
                    break;
            }
        BassResponse response = data.Response;
        f.Add(response is FunctionalResponse ? "functional" : "scalar").Add(response.Cases).Add(response.Locations).Add(response.Values);
        if (response is FunctionalResponse functional)
        {
            f.Add(functional.OutputVariables.Count);
            foreach (ContinuousPredictor v in functional.OutputVariables) f.Add(v.Name).Add(v.Values.Span);
        }
        return f.ToString();
    }
}
