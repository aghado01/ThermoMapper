using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.ExceptionServices;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A basis fit (R's <c>bassPCA</c> and <c>bassBasis</c>; P02, P03): the decomposition, the shared predictors and one
/// <see cref="BassFit"/> per component, each fitted to its coefficient response as a scalar response, with the joint
/// draws that pair them (<see cref="BassJointDraws"/>). Immutable; the components' models are built once, on first use.
/// <para>
/// Under the joint option (E34's decision 8, <see cref="BassComponents.FitJoint"/>) the scores are fitted together as several
/// outputs on one structure, <see cref="Joint"/>, and each component is a view of it: the joint fit's draws with component k's
/// coefficients and Σ_kk as σ². The readers of the components' surfaces read the views as they read components fitted apart; the
/// readers of their noise read Σ whole, its covariances across the components included.
/// </para>
/// </summary>
public sealed class BassBasisFit
{
    private readonly BassPredictor[] _predictors;
    private readonly BassFit[] _components;
    private readonly Lazy<BassModel>[] _models;
    private readonly Lazy<BassDraw[]>? _jointPooled;

    internal BassBasisFit(BassDecomposition decomposition, BassPredictor[] predictors, BassFit[] components, BassJointDraws joint, int masterSeed,
                          BassFit? jointFit = null, BassModel[]? models = null)
    {
        Decomposition = decomposition;
        _predictors = predictors;
        _components = components;
        JointDraws = joint;
        MasterSeed = masterSeed;
        Joint = jointFit;
        _models = new Lazy<BassModel>[components.Length];
        for (int k = 0; k < components.Length; k++)
        {
            int component = k;
            _models[k] = models is null
                ? new Lazy<BassModel>(() => BassModel.Create(_components[component].Config, ComponentData(component)))
                : new Lazy<BassModel>(() => models[component]);
        }
        if (jointFit is not null) _jointPooled = new Lazy<BassDraw[]>(() => [.. jointFit.Draws.SelectMany(chain => chain)]);
    }

    /// <summary>
    /// The fit of every component's coefficient response together, as several outputs on one structure with their noise covariance
    /// Σ integrated (E34's decision 8; <see cref="BassComponents.FitJoint"/>), whose views <see cref="Components"/> holds; null for
    /// components fitted apart. Its own readers (<see cref="BassMultiOutput"/>) read it on <see cref="JointData"/>, and
    /// <see cref="Bass.Resume"/> continues it.
    /// </summary>
    public BassFit? Joint { get; }

    /// <summary>The shared predictors and every component's coefficient response as one response of several outputs, the components
    /// in order and named w0, w1, …: the data <see cref="Joint"/> was run on.</summary>
    public BassData JointData() => BassComponents.JointData(Decomposition, _predictors);

    /// <summary>Joint draw <paramref name="jointDraw"/>'s noise covariance across the components, Σ (r × r), under the joint option;
    /// null for components fitted apart, whose noises are independent.</summary>
    internal double[,]? NoiseCovariance(int jointDraw)
        => _jointPooled is null ? null : _jointPooled.Value[JointDraws.Draw(0, jointDraw)].Sigma;

    /// <summary>The lower Cholesky factor of <see cref="NoiseCovariance"/>, L with LLᵀ = Σ; null for components fitted apart.</summary>
    internal double[,]? NoiseRoot(int jointDraw) => NoiseCovariance(jointDraw) is { } sigma ? BassComponents.Root(sigma) : null;

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

    /// <summary>
    /// The active subspaces of the reconstructed output over the joint draws (E8's slice 3; <see cref="BassFit.ActiveSubspace"/>):
    /// C = E[JᵀW_yJ] with J the Jacobian of y = y_m + y_s∘(Kw), so C = Σ_kk′ (GᵀW_yG)_kk′ E[∇w_k∇w_k′ᵀ] with G = diag(y_s)K, under the
    /// spec's output metric W_y, the identity by default. On the sampled route, product law or joint sample, each component's gradient
    /// readout at the shared points, with W_y's root taken in component space, GᵀW_yG (E8's slice 4).
    /// </summary>
    public ActiveSubspaceAnalysis ActiveSubspace(ActiveSubspaceSpec? spec = null) => BassActiveSubspaces.Of(this, spec ?? new ActiveSubspaceSpec());

    /// <summary>
    /// The co-active subspaces of two outputs of the reconstruction (E8's slice 3), f = a_fᵀy and g = a_gᵀy for the spec's output functionals
    /// <see cref="CoActiveSpec.Output"/> and <see cref="CoActiveSpec.OtherOutput"/> over the training output locations, over the joint
    /// draws, on either route.
    /// </summary>
    public CoActiveAnalysis CoActiveSubspace(CoActiveSpec spec) => BassActiveSubspaces.CoActive(this, spec);

    /// <summary>
    /// Basis adequacy (E8's slice 5): beside <paramref name="direct"/>, a direct-mode fit of the same cases, the share of its
    /// G = E[J W_x Jᵀ] outside this fit's span, G = diag(y_s)K, and the containment defect of G's leading patterns in it, per draw of
    /// the direct fit and for E[G], under the spec's law, input metric and output metric. The spec's <see cref="ActiveSubspaceSpec.OutputPatterns"/>
    /// gives the patterns read, 1 when 0. Did the truncation drop structure that moves with the inputs? This fit's own G lies inside its
    /// span by construction, so the direct fit answers.
    /// </summary>
    public Maths.Regression.Subspaces.BasisAdequacy Adequacy(BassFit direct, BassData directData, ActiveSubspaceSpec? spec = null)
        => BassActiveSubspaces.Adequacy(this, direct, directData, spec ?? new ActiveSubspaceSpec());

    /// <summary>
    /// The main and two-way effect functions of the reconstructed output over the joint draws (O8 for basis fits;
    /// <see cref="BassBasisEffects.Of"/>) under <paramref name="spec"/>'s measure: per output location and integrated over
    /// the locations, in output units.
    /// </summary>
    public BasisEffectAnalysis Effects(SensitivitySpec? spec = null, EffectSpec? effects = null) => BassBasisEffects.Of(this, spec, effects);

    /// <summary>Each component's traces, in component order (stage 9, P09; <see cref="BassTraces"/>), R's per-component
    /// trace panels.</summary>
    public IReadOnlyList<BassTraces> Traces() => [.. _components.Select(c => c.Traces())];

    /// <summary>
    /// Fitted against observed at the training cases in physical units (stage 9, P11; <see cref="BassTrainingFit"/>): the
    /// posterior mean of the latent output and the quantiles at <paramref name="probabilities"/> (each in (0, 1);
    /// <see cref="PosteriorSummary.DefaultProbabilities"/> when null) of the predictive draws <see cref="Predict"/> makes from
    /// <paramref name="seed"/>, with the residuals.
    /// </summary>
    public BassTrainingFit TrainingFit(IReadOnlyList<double>? probabilities = null, int seed = 0) => BassTrainingFit.Of(this, probabilities, seed);

    /// <summary>The basis fit's metadata (stage 9, P13; <see cref="BassBasisMetadata"/>): the data's shape, the
    /// decomposition, the joint draws and each component's metadata.</summary>
    public BassBasisMetadata Metadata() => BassBasisMetadata.Of(this);

    /// <summary>The metadata with each component's posterior summaries at <paramref name="probabilities"/> (stage 9, P15;
    /// <see cref="BassBasisSummary"/>).</summary>
    public BassBasisSummary Summary(IReadOnlyList<double>? probabilities = null) => BassBasisSummary.Of(this, probabilities);

    /// <summary>The posterior summary of <see cref="Sensitivity"/>'s indices and variances, per location and integrated, at
    /// <paramref name="probabilities"/> (stage 9, P11; <see cref="BasisSobolSummary"/>).</summary>
    public BasisSobolSummary SensitivitySummary(SensitivitySpec? spec = null, IReadOnlyList<double>? probabilities = null)
        => BasisSobolSummary.Of(Sensitivity(spec), probabilities);
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
            if (config.HasRandomEffects)
                throw new ArgumentException($"Config {k} declares random intercepts; a decomposition's components share predictors and no grouping, and a block across components is not built (E21).", nameof(configs));
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

    /// <summary>
    /// Fits every component of <paramref name="decomposition"/> together (E34's decision 8): the r coefficient responses as r outputs
    /// of one fit on <paramref name="predictors"/> under <paramref name="config"/>, one structure for all of them with their noise
    /// covariance Σ integrated (<see cref="MultiOutputResponse"/>), the truncation residuals as for components fitted apart. Under
    /// decision 2's Σ prior the structure's and τ's posterior depend on the scores' span alone, not on the basis that spans it; the
    /// trade is one structure for every component in place of a structure each. The joint fit runs at the seed tree's child 0, a
    /// component fitted apart's first, so a rank-1 decomposition fits as <see cref="Fit(BassDecomposition, IReadOnlyList{BassPredictor}, BassConfig, BassRunOptions?)"/> does.
    /// </summary>
    public static BassBasisFit FitJoint(BassDecomposition decomposition, IReadOnlyList<BassPredictor> predictors, BassConfig config,
                                        BassRunOptions? run = null)
    {
        ArgumentNullException.ThrowIfNull(decomposition);
        ArgumentNullException.ThrowIfNull(config);
        if (config.HasRandomEffects)
            throw new ArgumentException("The config declares random intercepts; a decomposition's components share predictors and no grouping (E21), " +
                                        "and several outputs' block arrives with E34's slice 4.", nameof(config));
        int r = decomposition.Rank;
        BassPredictor[] shared = Checked(predictors, decomposition.Cases);
        BassConfig seeded = config with { Sampling = config.Sampling with { MasterSeed = SeedTree.Derive(config.Sampling.MasterSeed, r)[0] } };
        BassFit joint = Bass.Fit(seeded, JointData(decomposition, shared), run);
        return AssembleJoint(decomposition, shared, joint, config.Sampling.MasterSeed);
    }

    /// <summary>
    /// The basis fit of a joint fit already run, <paramref name="joint"/> on <paramref name="decomposition"/>'s scores as several
    /// outputs (refused when its fingerprint says otherwise), as <see cref="FitJoint"/> assembles it: one view per component, and the
    /// joint draws each one of <paramref name="joint"/>'s draws, the pairing's stream derived from <paramref name="masterSeed"/>.
    /// </summary>
    public static BassBasisFit AssembleJoint(BassDecomposition decomposition, IReadOnlyList<BassPredictor> predictors, BassFit joint, int masterSeed)
    {
        ArgumentNullException.ThrowIfNull(decomposition);
        ArgumentNullException.ThrowIfNull(joint);
        int r = decomposition.Rank;
        BassPredictor[] shared = Checked(predictors, decomposition.Cases);
        BassData data = JointData(decomposition, shared);
        if (Bass.Fingerprint(joint.Config, data) != joint.Identity.Fingerprint)
            throw new ArgumentException("The joint fit was run on other data or another config than the decomposition's scores.", nameof(joint));
        var views = new BassFit[r];
        var models = new BassModel[r];
        for (int k = 0; k < r; k++) (views[k], models[k]) = View(joint, k, ComponentDataOf(decomposition, shared, k));
        BassJointDraws pairs = BassJointDraws.Shared(joint, r, SeedTree.Derive(masterSeed, r + 1)[r]);
        return new BassBasisFit(decomposition, shared, views, pairs, masterSeed, joint, models);
    }

    /// <summary>The shared predictors and the r coefficient responses as one response of r outputs, named w0, w1, … in component
    /// order: the data a joint fit runs on.</summary>
    internal static BassData JointData(BassDecomposition decomposition, BassPredictor[] predictors)
    {
        int r = decomposition.Rank, n = decomposition.Cases;
        var y = new double[n, r];
        for (int k = 0; k < r; k++)
        {
            double[] w = decomposition.CoefficientResponse(k);
            for (int i = 0; i < n; i++) y[i, k] = w[i];
        }
        return new BassData(predictors, new MultiOutputResponse(y, [.. Enumerable.Range(0, r).Select(k => $"w{k}")]));
    }

    private static BassData ComponentDataOf(BassDecomposition decomposition, BassPredictor[] predictors, int component)
        => new(predictors, new ScalarResponse(decomposition.CoefficientResponse(component)));

    // Component k of a joint fit as a fit of its own response: the joint draws with column k of the coefficients and Σ_kk as σ², on a
    // scalar config whose model reads what the components' readers read (the basis, the inputs, the coefficients' conditional means,
    // which do not depend on Σ). Σ's prior is the joint fit's, so the view's config drops it, and its births are the ones that ran. The
    // evidence, the checkpoint and the residual check are the joint fit's, so the view holds none.
    private static (BassFit View, BassModel Model) View(BassFit joint, int k, BassData data)
    {
        BassSampling sampling = joint.Config.Sampling;
        BirthMixture? births = sampling.Births ?? (joint.Resolved.Births is null ? new BirthMixture() : null);
        BassConfig config = joint.Config with { OutputCovariance = null, Sampling = sampling with { Births = births } };
        BassModel model = BassModel.Create(config, data);
        var stores = new DrawStore<BassDraw>[joint.Draws.Count];
        for (int c = 0; c < stores.Length; c++)
        {
            var store = new DrawStore<BassDraw>();
            foreach (BassDraw draw in joint.Draws[c]) store.Add(ViewDraw(draw, k));
            stores[c] = store;
        }
        BassFit view = joint with
        {
            Config = config, Resolved = model.Resolved, Identity = joint.Identity with { Fingerprint = Bass.Fingerprint(config, data) },
            Draws = stores, Evidence = null, ResidualDependence = null, Checkpoint = null,
        };
        return (view, model);
    }

    private static BassDraw ViewDraw(BassDraw draw, int k)
    {
        double[,] a = draw.OutputCoefficients ?? throw new ArgumentException("A joint fit's draw holds several outputs' coefficients.");
        double[,] sigma = draw.Sigma ?? throw new ArgumentException("A joint fit's draw holds Σ.");
        var coefficients = new double[a.GetLength(0)];
        for (int j = 0; j < coefficients.Length; j++) coefficients[j] = a[j, k];
        return new BassDraw(draw.Structure, draw.Tau, sigma[k, k], coefficients, draw.Lambda);
    }

    // The lower Cholesky factor of a drawn Σ, positive definite by construction.
    internal static double[,] Root(double[,] sigma)
    {
        int r = sigma.GetLength(0);
        var l = new double[r, r];
        for (int i = 0; i < r; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = sigma[i, j];
                for (int q = 0; q < j; q++) sum -= l[i, q] * l[j, q];
                if (i == j)
                {
                    if (!(sum > 0.0)) throw new InvalidOperationException("A draw of Σ is not positive definite in doubles.");
                    l[i, i] = Math.Sqrt(sum);
                }
                else l[i, j] = sum / l[j, j];
            }
        return l;
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
