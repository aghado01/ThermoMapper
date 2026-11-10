# Random effects on the BASS core (E21's ridge blocks)

> Random effects as ridge blocks on BASS's collapsed core: group intercepts and slopes, crossed and nested groupings, their scales held or drawn, under the Gaussian likelihood and GBASS's noises, and their reading by the conformal layer. Built slice by slice from [the brief](../../issues/regression/bass/ridge-blocks-brief.md). KHAOS's first parts, the prior-precision builder and random intercepts on one grouping with s held, are on [BayesianRegression.md](./BayesianRegression.md#bass-random-effect-intercepts-e21), which this page extends and does not repeat.

## How a block is carried: the whitening on the working response (slice 0)

Built 2026-11-08. A random-effect block is integrated out of every structure's marginal by whitening the rows: with V = I + s²(ZH)(ZH)ᵀ, the response and every column are taken through V^(−½), and the log density pays −½ln|V| ([KHAOS's form](./BayesianRegression.md#bass-random-effect-intercepts-e21)). V depends on s, and under GBASS's noises on the weights, both of which a fit may move with its state. So the whitening is carried by the reading it belongs to.
- *The working response carries the block* (`WorkingResponse.Whitening`): the block whose V^(−½) its offset applies, its ln|V| in the log-Jacobian beside the weights' Σ ln wᵢ. `RandomEffectsBlock.Whitened` sets it, and a reading moved to another offset keeps it.
- *The design whitens each column by the reading it is on,* not by a block of its own. Whether a design's columns are whitened is its layout, fixed when it is built, as whether it holds weights is.
- *A reading at another scale moves the design as moved weights do* (`BassDesign.ProposeFollow`): every column is evaluated and whitened afresh, and the Gram is factored strictly, by the arithmetic of a design built on that reading and loaded with the same terms. This is what a drawn s (slice 1) and the block under the noises (slice 4) ride.
- *What it changed:* nothing a fit draws. A held-s block's fits keep their draws bit for bit (`RandomEffectsTests`: a hinge chain, a Legendre chain and a tempered hinge fit with its reference rung, pinned by hash before the change), and so does every fit without a block.

**The hinge prefix tables stand aside under a block** (a defect fixed with the change). The tables give a one-way hinge's border by summing each held column against the raw hinge. Under a block the held columns are whitened, V^(−½)xⱼ, so the tables read xⱼᵀV^(−½)h where the border is xⱼᵀV⁻¹h. They switch on at 500 cases and up to 16 continuous predictors, which KHAOS's gates (n 246) never reached. Building them under a block was already refused (`WarmTables`); reading them was not, and on HEAD a proposal the dots admitted came back refused as numerical. A design with whitened columns now reads every border by the dots: on 600 cases, 1,500 proposals agree with the dots' to the bit, the tables read 0 entries.

**Checked** (`RandomEffectsTests`). The three pins. A design loaded at s = 1.3, followed to a reading at s = 0.4, against a design built at 0.4 with the same terms: the residual and fitted squares and the coefficients to the bit (the gate holds them to 10⁻¹²); followed back to 1.3, its first fit within 10⁻¹²; a reading without a whitening refused as another layout. The tables under a block, as above.

## Drawing the scale of one Gaussian intercept grouping (slice 1)

`BassRandomIntercept.Scale` is the initial group standard deviation in units of the residual standard deviation σ. Supply `ScalePrior` to draw s under a half-Cauchy of scale A. The prior defaults to A = 1. Omitting `ScalePrior` holds s at the specified value. Each declaration names its data grouping.

```csharp
RandomIntercepts = [new BassRandomIntercept("site", 1.3, ScalePrior: new BassGroupScalePrior(1.0))]
```

The Gaussian sweep moves the structure, updates τ if sampled, then proposes a change in v = ln s with the coefficients, effects and noise variance integrated as in the structure's marginal. The proposed reading starts from the design's raw observations, whitens them at the new s, and carries the design's frozen prior reference. In particular a response-scaled σ² prior reads the original observations' mean square throughout the fit. A restored replica retains its checkpointed response and reference objects, including after a swap.

The joint target includes ln p(s) + v. The scale move uses the rung's density and the reverse-to-forward proposal density ratio. For the proposal width, the eigenvalues λ of the group count matrix on the observation-centred subspace give the Gaussian information I = 2 Σ (s²λ/(1+s²λ))². With r = (s/A)² the width is 2.4 / √max(1, I + 4r/(1+r)²). The implementation shares the grouping's eigensystem across scales. Each proposal re-whitens the current structure and refactors its Gram.

Reference rungs and SMC particles draw s from its proper prior along with the structure and τ, then admit the structure under that whitening. Each record in `BassDraw.Groups` carries its grouping name, scale and conditional effects. Known groups use their draw's effect. A new group's Gaussian predictive variance given the draw is σ²(1+s²). The flat-histogram configuration refuses a drawn s because its weights and WHAM reconstruction have not been qualified on this moving marginal.

## Named crossed and nested intercepts (slice 2)

`BassData.Groupings` is a collection of named categorical columns beside the predictors. `BassConfig.RandomIntercepts` names the columns it reads, in the order of the scale coordinates. Data columns, declarations and retained records resolve by name. A grouping used only as a parent need not have its own intercept block. Unknown names, duplicate names and unused data groupings are refused.

```csharp
var data = new BassData(predictors, response, [site, classroom]);
var config = new BassConfig
{
    RandomIntercepts =
    [
        new("site", 1.0, ScalePrior: new BassGroupScalePrior()),
        new("classroom", 0.8, RandomInterceptConvention.WithinParent,
            ScalePrior: new BassGroupScalePrior(), ParentGrouping: "site"),
    ],
};
```

Each grouping declares one convention:

- `ObservationCentred`, the default, enforces Σ_g n_g u_g = 0. For crossed groupings each contribution has mean zero over the training cases, and the intercept carries the level.
- `WithinParent` enforces Σ_{g in p} n_g u_g = 0 separately inside every observed parent level. Each child level must belong to one parent. The parent contribution carries each parent's mean. A parent with one observed child has zero child contrast there.

Let H_k be an orthonormal basis of the level-space constraints, Z_k the indicators, and S the diagonal matrix repeating each s_k for its contrast coordinates. The joint row covariance is V = I + Z̃ S² Z̃ᵀ, with Z̃ = [Z₁H₁ … Z_KH_K]. One thin SVD of the scale-free design gives Q = U_r and R = diag(σ_r)V_rᵀ. At each scale vector a reduced SVD of R·S supplies the r-dimensional covariance's eigenvectors and squared singular values. Forming RS²Rᵀ explicitly would square the condition number and lose a small component next to a much larger one; the reduced factor avoids that loss. Whitening costs O(nr + r²) per column and O(r³ + r²q) per scale vector; no n × n covariance is factored in production.

Separate group whitenings cannot be composed for crossed groupings. Under the within-parent convention the parent and child spans are orthogonal, but the same joint machinery serves both shapes. A single observation-centred block delegates to the existing closed form, retaining its arithmetic and random draw order.

An exact alias between levels of crossed groupings is reported in the resolution log, including a crossed grouping that aliases a within-parent child. It adds no constraint. The conditional simulation retains every original contrast coordinate, including latent directions lost from the row rank. With R·S = U·diag(d)·Vᵀ, the mean uses S·V·diag(d/(1+d²))·UᵀQᵀr; the standardized covariance square root is I + V·diag((1+d²)^(−½)−1)·Vᵀ. This retains the full latent complement without a subtraction of large prior effects. Unequal proper scales determine how the effects share the signal and its uncertainty.

The sweep updates each drawn ln s_k separately after τ. Its information is ½ tr[(V⁻¹ ∂V/∂ln s_k)²], with the same prior curvature and Hastings correction as the single-group step. The joint target sums each drawn scale's log density and log-coordinate Jacobian. `BassState.GroupScales` reads the fitted whitening's scale vector. `BassDraw.Groups` contains named `BassGroupDraw` records with `Scale` and `Effects`; effects follow the observed levels' declaration order. Prediction adds every known effect, and each missing grouping contributes its draw's population spread. These blocks currently take a scalar Gaussian response; moving noise weights, slopes and the conformal candidate density belong to later slices.

## Migrating a saved one-group configuration

Replace the singular object and supply its data grouping's name:

```json
{"RandomIntercept":{"Scale":1.3,"Convention":"ObservationCentred","ScalePrior":{"Scale":1}}}
```

becomes

```json
{"RandomIntercepts":[{"Grouping":"site","Scale":1.3,"Convention":"ObservationCentred","ScalePrior":{"Scale":1}}]}
```

Pass `[site]` to `BassData` in place of `site`. Read `draw.Groups[0].Effects` and `.Scale` in place of `GroupEffects` and `GroupScale`, resolving by `Grouping` when order is not known. The old singular JSON key is rejected.

Grouped configuration JSON changes the run fingerprint. Converting a saved config alone therefore cannot resume an old grouped checkpoint under the existing fingerprint guard; start a new grouped run. Fits without groupings retain their previous fingerprints. New checkpoints carry the named representation.
