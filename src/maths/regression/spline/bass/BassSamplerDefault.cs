using System;
using Maths.Regression.Spline.Bars;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// The sampler a fit runs when its config declares none (<see cref="BassSampling.Sampler"/>; the user's decision of 2026-11-15, E36): plain
/// chains on every basis but one that admits ridge blocks (<see cref="DirectionalBasis.Columns"/> ≥ 2), which takes self-sizing parallel
/// tempering. E36's pure-basis readings found a ridge basis's plain chains holding M where they started (R̂(M) 5 to 8 over four starts) where
/// tempering's chains and SMC's replicates agreed.
/// <list type="bullet">
/// <item><b>The ladder.</b> Self-sizing (<see cref="TemperingSpec.Levels"/> null) on the likelihood path in the form the model runs (the
/// integrated form unless its caller names one), with the exact reference rung where the declared priors and that form admit it, as the
/// engine's own validation and resolution read them. Otherwise a positive hot end, auto-picked by the resolution, and the reference rung's
/// refusal is recorded as the reason.</item>
/// <item><b>What wins.</b> A declared sampler (<see cref="BassSampling.Tempering"/>, <see cref="BassSampling.Smc"/>,
/// <see cref="BassSampling.FlatHistogram"/>) is run as declared. Declared plain chains (<see cref="BassSamplerChoice.PlainChains"/>) run, and on
/// a ridge basis carry a caution. A resumed checkpoint runs the sampler that produced it, tempered or not, whatever the rule now says.</item>
/// <item><b>Every other fit</b> is untouched: its config passes through as it came, bit for bit, and nothing is logged.</item>
/// </list>
/// The rule reads the configured basis, so a fit whose draws hold no block is resolved as one that holds some.
/// </summary>
internal static class BassSamplerDefault
{
    private const string Caution =
        "plain chains on a ridge basis can hold M where they started, and f's intervals can then be too narrow (E36's pure-basis readings: "
        + "R̂(M) 5 to 8 over four starts); tempering or SMC is the default for this basis";

    /// <summary>Whether <paramref name="config"/>'s basis admits ridge blocks.</summary>
    public static bool AdmitsRidgeBlocks(BassConfig config) => config.Structure.Directional is { Columns: >= 2 };

    /// <summary>
    /// The config the run takes, its log row (null where nothing was resolved) and the caution a plain chain on a ridge basis carries.
    /// <paramref name="form"/> is the tempered form the model's caller names, null for the config's own; <paramref name="resumedTempered"/>
    /// is a resumed checkpoint's own sampler: true for tempered chains, false for plain ones, null for a fresh run.
    /// </summary>
    public static (BassConfig Config, BassAutoPick? Pick, string? Caution) Resolve(BassConfig config, BassData data, WorkingResponse? working,
                                                                                   TemperedForm? form = null, bool? resumedTempered = null)
    {
        ArgumentNullException.ThrowIfNull(config);
        BassSampling sampling = config.Sampling;
        if (!AdmitsRidgeBlocks(config) || sampling.Tempering is not null || sampling.Smc is not null || sampling.FlatHistogram is not null)
            return (config, null, null);
        if (sampling.Sampler == BassSamplerChoice.PlainChains)
            return (config, new BassAutoPick("Sampler", "plain chains", "as declared (Sampling.Sampler = PlainChains); " + Caution), Caution);
        if (resumedTempered == false)
            return (config, new BassAutoPick("Sampler", "plain chains", "as the resumed checkpoint ran; " + Caution), Caution);

        // A ridge basis, sampler automatic: self-sizing parallel tempering, with the reference rung where it is admitted.
        TemperedForm runs = form ?? TemperedForm.Integrated;
        BassConfig referenced = With(config, new TemperingSpec { ReferenceRung = true, Form = runs });
        string? refusal = Refusal(referenced, data, working);
        if (refusal is null)
            return (referenced, new BassAutoPick("Sampler", "parallel tempering, self-sizing, exact reference rung",
                                                 "automatic for a ridge basis (K ≥ 2): its plain chains can hold M where they started"), null);
        BassConfig hot = With(config, new TemperingSpec { Form = runs });
        string? hotRefusal = Refusal(hot, data, working);
        if (hotRefusal is null)
            return (hot, new BassAutoPick("Sampler", "parallel tempering, self-sizing, auto-picked hot end",
                                          $"automatic for a ridge basis (K ≥ 2); no reference rung, since {refusal}"), null);
        // Neither ladder is admitted: the config's own refusal surfaces from the resolution as it would have, so run it as declared.
        return (config, new BassAutoPick("Sampler", "plain chains", $"no ladder is admitted ({hotRefusal}); " + Caution), Caution);
    }

    private static BassConfig With(BassConfig config, TemperingSpec tempering) => config with { Sampling = config.Sampling with { Tempering = tempering } };

    // Why the engine refuses `candidate` on `data`, by its own validation and resolution; null where it admits it.
    private static string? Refusal(BassConfig candidate, BassData data, WorkingResponse? working)
    {
        try
        {
            BassResolved.Resolve(candidate, data, working);
            return null;
        }
        catch (ArgumentException e)
        {
            string suffix = $" (Parameter '{e.ParamName}')";
            return e.Message.EndsWith(suffix, StringComparison.Ordinal) ? e.Message[..^suffix.Length] : e.Message;
        }
    }
}
