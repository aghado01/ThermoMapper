using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.ExceptionServices;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A basis fit (R's <c>bassPCA</c> and <c>bassBasis</c>; P02, P03): the decomposition, the shared predictors and one
/// <see cref="BassFit"/> per component, each fitted to its coefficient response as a scalar response, with the joint
/// draws that pair them (<see cref="BassJointDraws"/>). Immutable; the components' models are built once, on first use.
/// </summary>
public sealed class BassBasisFit
{
    private readonly BassPredictor[] _predictors;
    private readonly BassFit[] _components;
    private readonly Lazy<BassModel>[] _models;

    internal BassBasisFit(BassDecomposition decomposition, BassPredictor[] predictors, BassFit[] components, BassJointDraws joint, int masterSeed)
    {
        Decomposition = decomposition;
        _predictors = predictors;
        _components = components;
        JointDraws = joint;
        MasterSeed = masterSeed;
        _models = new Lazy<BassModel>[components.Length];
        for (int k = 0; k < components.Length; k++)
        {
            int component = k;
            _models[k] = new Lazy<BassModel>(() => BassModel.Create(_components[component].Config, ComponentData(component)));
        }
    }

    /// <summary>The decomposition.</summary>
    public BassDecomposition Decomposition { get; }

    /// <summary>The predictors every component shares, in the caller's order.</summary>
    public IReadOnlyList<BassPredictor> Predictors => _predictors;

    /// <summary>Each component's fit, in component order.</summary>
    public IReadOnlyList<BassFit> Components => _components;

    /// <summary>The draws paired across the components.</summary>
    public BassJointDraws JointDraws { get; }

    /// <summary>The master seed the components' seeds and the pairing's stream derive from.</summary>
    public int MasterSeed { get; }

    /// <summary>Component <paramref name="component"/>'s data: the shared predictors and its coefficient response.</summary>
    public BassData ComponentData(int component) => new(_predictors, new ScalarResponse(Decomposition.CoefficientResponse(component)));

    internal BassModel Model(int component) => _models[component].Value;

    /// <summary>
    /// Prediction at <paramref name="cases"/> (D07, P08; <see cref="BassBasisPrediction"/>), over the joint draws:
    /// y = y_m + y_s ∘ (K(w(x) + ε) + e). The basis exists only at its training locations, so new output locations are
    /// refused.
    /// </summary>
    public BassBasisPrediction Predict(PredictionCases cases, BassPredictOptions? options = null)
        => BassBasisPredictor.Predict(this, cases, options ?? new BassPredictOptions());

    /// <summary>
    /// The draws × units log density over the joint draws (OQ17; <see cref="BassBasisDensity"/>): in coefficient space by
    /// default, each case its own unit; or in physical units under the declared diagonal Gaussian truncation model.
    /// </summary>
    public PointwiseLogDensity PointwiseLogDensity(UnitPartition? units = null, BasisDensitySpace space = BasisDensitySpace.Coefficients)
        => BassBasisDensity.LogDensity(this, units, space);

    /// <summary>
    /// The Sobol decomposition of the reconstructed output over the joint draws (D08, O7; <see cref="BassBasisSensitivity.Of"/>)
    /// under <paramref name="spec"/>'s measure: per output location and integrated over the locations, in output units, the
    /// cross-component covariances kept. The components share this fit's predictors, so no data is passed.
    /// </summary>
    public BasisSobolAnalysis Sensitivity(SensitivitySpec? spec = null) => BassBasisSensitivity.Of(this, spec);
}

/// <summary>
/// Component execution (D07; R's <c>bassBasis</c>): each coefficient response fitted as a <see cref="ScalarResponse"/> on
/// the shared predictors through <see cref="Bass.Fit"/>, one config for every component by default.
/// <list type="bullet">
/// <item><b>Seeds.</b> The master seed fans out through the seed tree: component k runs at child k, and the joint draws'
/// resampling stream is child r, so a component's draws do not depend on how many others there are beyond it.</item>
/// <item><b>One bound on parallelism.</b> <see cref="BassRunOptions.MaxDegreeOfParallelism"/> bounds components and
/// chains together: with P threads and C chains per fit, max(1, ⌊P/C⌋) components run at once, each with min(P, C)
/// chains at once; unbounded, both are left to the machine. Results do not depend on it.</item>
/// <item><b>Progress</b> comes per component, from each one's own fit.</item>
/// </list>
/// </summary>
public static class BassComponents
{
    /// <summary>Fits every component of <paramref name="decomposition"/> on <paramref name="predictors"/> under
    /// <paramref name="config"/>, the master seed its <see cref="BassSampling.MasterSeed"/>.</summary>
    public static BassBasisFit Fit(BassDecomposition decomposition, IReadOnlyList<BassPredictor> predictors, BassConfig config,
                                   BassRunOptions? run = null)
    {
        ArgumentNullException.ThrowIfNull(decomposition);
        ArgumentNullException.ThrowIfNull(config);
        return Fit(decomposition, predictors, Enumerable.Repeat(config, decomposition.Rank).ToArray(), config.Sampling.MasterSeed, run);
    }

    /// <summary>
    /// Fits component k under <paramref name="configs"/>[k], its seed derived from <paramref name="masterSeed"/> by the
    /// seed tree, which replaces the config's own; the seeds belong to the basis fit, not to the components.
    /// </summary>
    public static BassBasisFit Fit(BassDecomposition decomposition, IReadOnlyList<BassPredictor> predictors, IReadOnlyList<BassConfig> configs,
                                   int masterSeed, BassRunOptions? run = null)
    {
        ArgumentNullException.ThrowIfNull(decomposition);
        ArgumentNullException.ThrowIfNull(predictors);
        ArgumentNullException.ThrowIfNull(configs);
        int r = decomposition.Rank;
        if (configs.Count != r) throw new ArgumentException($"Give one config per component ({r}).", nameof(configs));
        BassRunOptions options = run ?? new BassRunOptions();
        if (options.MaxDegreeOfParallelism is 0 or < -1)
            throw new ArgumentOutOfRangeException(nameof(run), "MaxDegreeOfParallelism is −1, no bound, or at least 1.");
        int[] seeds = SeedTree.Derive(masterSeed, r);
        var seeded = new BassConfig[r];
        for (int k = 0; k < r; k++)
        {
            BassConfig config = configs[k] ?? throw new ArgumentException($"Config {k} is null.", nameof(configs));
            seeded[k] = config with { Sampling = config.Sampling with { MasterSeed = seeds[k] } };
        }
        BassPredictor[] shared = Checked(predictors, decomposition.Cases);

        int bound = options.MaxDegreeOfParallelism, chains = seeded.Max(c => c.Sampling.Chains);
        int atOnce = bound == -1 ? -1 : Math.Max(1, bound / chains);
        BassRunOptions perFit = bound == -1 ? options : options with { MaxDegreeOfParallelism = Math.Min(bound, chains) };
        var fits = new BassFit[r];
        var parallel = new ParallelOptions { MaxDegreeOfParallelism = atOnce, CancellationToken = options.Cancellation };
        try
        {
            Parallel.For(0, r, parallel, k =>
                fits[k] = Bass.Fit(seeded[k], new BassData(shared, new ScalarResponse(decomposition.CoefficientResponse(k))), perFit));
        }
        catch (AggregateException e) when (Shared(e) is { } failure)
        {
            ExceptionDispatchInfo.Throw(failure);
            throw;
        }
        return Assemble(decomposition, shared, fits, masterSeed);
    }

    /// <summary>
    /// The basis fit of components already fitted, <paramref name="components"/>[k] to component k's data (refused when
    /// its fingerprint says otherwise), with the joint draws' stream derived from <paramref name="masterSeed"/>.
    /// </summary>
    public static BassBasisFit Assemble(BassDecomposition decomposition, IReadOnlyList<BassPredictor> predictors, IReadOnlyList<BassFit> components,
                                        int masterSeed)
    {
        ArgumentNullException.ThrowIfNull(decomposition);
        ArgumentNullException.ThrowIfNull(components);
        int r = decomposition.Rank;
        if (components.Count != r) throw new ArgumentException($"Give one fit per component ({r}).", nameof(components));
        BassPredictor[] shared = Checked(predictors, decomposition.Cases);
        var fits = new BassFit[r];
        for (int k = 0; k < r; k++)
        {
            BassFit fit = components[k] ?? throw new ArgumentException($"Component {k}'s fit is null.", nameof(components));
            string fingerprint = Bass.Fingerprint(fit.Config, new BassData(shared, new ScalarResponse(decomposition.CoefficientResponse(k))));
            if (fingerprint != fit.Identity.Fingerprint)
                throw new ArgumentException($"Component {k}'s fit was run on other data or another config than its component's.", nameof(components));
            fits[k] = fit;
        }
        BassJointDraws joint = BassJointDraws.Pair(fits, SeedTree.Derive(masterSeed, r + 1)[r]);
        return new BassBasisFit(decomposition, shared, fits, joint, masterSeed);
    }

    private static BassPredictor[] Checked(IReadOnlyList<BassPredictor> predictors, int cases)
    {
        ArgumentNullException.ThrowIfNull(predictors);
        var shared = predictors.ToArray();
        foreach (BassPredictor? p in shared)
        {
            if (p is null) throw new ArgumentException("A predictor is null.", nameof(predictors));
            if (p.Length != cases) throw new ArgumentException($"Predictor '{p.Name}' has {p.Length} values for the decomposition's {cases} cases.", nameof(predictors));
        }
        return shared;
    }

    // The components' one kind of failure, which the parallel fan-out wraps: the first, when every inner exception has its type.
    private static Exception? Shared(AggregateException e)
    {
        IReadOnlyList<Exception> inner = e.Flatten().InnerExceptions;
        return inner.Count > 0 && inner.All(x => x.GetType() == inner[0].GetType()) ? inner[0] : null;
    }
}
