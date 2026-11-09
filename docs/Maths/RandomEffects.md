# Random effects on the BASS core (E21's ridge blocks)

> Random effects as ridge blocks on BASS's collapsed core: group intercepts and slopes, crossed and nested groupings, their scales held or drawn, under the Gaussian likelihood and GBASS's noises, and their reading by the conformal layer. Built slice by slice from [the brief](../../issues/regression/bass/ridge-blocks-brief.md). KHAOS's first parts, the prior-precision builder and random intercepts on one grouping with s held, are on [BayesianRegression.md](./BayesianRegression.md#bass-random-effect-intercepts-e21), which this page extends and does not repeat.

## How a block is carried: the whitening on the working response (slice 0)

Built 2026-11-08. A random-effect block is integrated out of every structure's marginal by whitening the rows: with V = I + s²(ZH)(ZH)ᵀ, the response and every column are taken through V^(−½), and the log density pays −½ln|V| ([KHAOS's form](./BayesianRegression.md#bass-random-effect-intercepts-e21)). V depends on s, and under GBASS's noises on the weights, both of which a fit may move with its state. So the whitening is carried by the reading it belongs to.
- *The working response carries the block* (`WorkingResponse.Whitening`): the block whose V^(−½) its offset applies, its ln|V| in the log-Jacobian beside the weights' Σ ln wᵢ. `RandomInterceptBlock.Whitened` sets it, and a reading moved to another offset keeps it.
- *The design whitens each column by the reading it is on,* not by a block of its own. Whether a design's columns are whitened is its layout, fixed when it is built, as whether it holds weights is.
- *A reading at another scale moves the design as moved weights do* (`BassDesign.ProposeFollow`): every column is evaluated and whitened afresh, and the Gram is factored strictly, by the arithmetic of a design built on that reading and loaded with the same terms. This is what a drawn s (slice 1) and the block under the noises (slice 4) ride.
- *What it changed:* nothing a fit draws. A held-s block's fits keep their draws bit for bit (`RandomEffectsTests`: a hinge chain, a Legendre chain and a tempered hinge fit with its reference rung, pinned by hash before the change), and so does every fit without a block.

**The hinge prefix tables stand aside under a block** (a defect fixed with the change). The tables give a one-way hinge's border by summing each held column against the raw hinge. Under a block the held columns are whitened, V^(−½)xⱼ, so the tables read xⱼᵀV^(−½)h where the border is xⱼᵀV⁻¹h. They switch on at 500 cases and up to 16 continuous predictors, which KHAOS's gates (n 246) never reached. Building them under a block was already refused (`WarmTables`); reading them was not, and on HEAD a proposal the dots admitted came back refused as numerical. A design with whitened columns now reads every border by the dots: on 600 cases, 1,500 proposals agree with the dots' to the bit, the tables read 0 entries.

**Checked** (`RandomEffectsTests`). The three pins. A design loaded at s = 1.3, followed to a reading at s = 0.4, against a design built at 0.4 with the same terms: the residual and fitted squares and the coefficients to the bit (the gate holds them to 10⁻¹²); followed back to 1.3, its first fit within 10⁻¹²; a reading without a whitening refused as another layout. The tables under a block, as above.

## Drawing the scale of one Gaussian intercept grouping (slice 1)

`BassRandomIntercept.Scale` is the initial group standard deviation in units of the residual standard deviation σ. Supply `ScalePrior` to draw s under a half-Cauchy of scale A. The prior defaults to A = 1. Omitting `ScalePrior` holds s at the specified value and retains the held-scale configuration's serialization and fingerprint.

```csharp
RandomIntercept = new BassRandomIntercept(1.3, ScalePrior: new BassGroupScalePrior(1.0))
```

The Gaussian sweep moves the structure, updates τ if sampled, then proposes a change in v = ln s with the coefficients, effects and noise variance integrated as in the structure's marginal. The proposed reading starts from the design's raw observations, whitens them at the new s, and carries the design's frozen prior reference. In particular a response-scaled σ² prior reads the original observations' mean square throughout the fit. A restored replica retains its checkpointed response and reference objects, including after a swap.

The joint target includes ln p(s) + v. The scale move uses the rung's density and the reverse-to-forward proposal density ratio. For the proposal width, the eigenvalues λ of the group count matrix on the observation-centred subspace give the Gaussian information I = 2 Σ (s²λ/(1+s²λ))². With r = (s/A)² the width is 2.4 / √max(1, I + 4r/(1+r)²). The implementation shares the grouping's eigensystem across scales. Each proposal re-whitens the current structure and refactors its Gram.

Reference rungs and SMC particles draw s from its proper prior along with the structure and τ, then admit the structure under that whitening. `BassDraw.GroupScale` carries s alongside the conditional effects. Known groups use their draw's effect. A new group's Gaussian predictive variance given the draw is σ²(1+s²). The flat-histogram configuration refuses a drawn s because its weights and WHAM reconstruction have not been qualified on this moving marginal.
