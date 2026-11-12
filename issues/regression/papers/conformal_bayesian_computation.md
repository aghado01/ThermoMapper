# Conformal Bayesian Computation

Edwin Fong  Affiliation: University of Oxford  Affiliation: The Alan Turing Institute  Email: [edwin.fong@stats.ox.ac.uk](mailto:)    Chris Holmes  Affiliation: University of Oxford  Affiliation: The Alan Turing Institute  Email: [cholmes@stats.ox.ac.uk](mailto:)

###### Abstract

We develop scalable methods for producing conformal Bayesian predictive intervals with finite sample calibration guarantees. Bayesian posterior predictive distributions, p⁡(y∣x)p(y\\mid x), characterize subjective beliefs on outcomes of interest, yy, conditional on predictors, xx. Bayesian prediction is well-calibrated when the model is true, but the predictive intervals may exhibit poor empirical coverage when the model is misspecified, under the so called ℳ{\\cal{M}}-open perspective. In contrast, conformal inference provides finite sample frequentist guarantees on predictive confidence intervals without the requirement of model fidelity. Using ‘add-one-in’ importance sampling, we show that conformal Bayesian predictive intervals are efficiently obtained from re-weighted posterior samples of model parameters. Our approach contrasts with existing conformal methods that require expensive refitting of models or data-splitting to achieve computational efficiency. We demonstrate the utility on a range of examples including extensions to partially exchangeable settings such as hierarchical models.

## 1 Introduction

We consider Bayesian prediction using training data Z1:n={Yi,Xi}i=1:nZ\_{1:n}=\\{Y\_{i},X\_{i}\\}\_{i=1:n} for an outcome of interest YiY\_{i} and covariates Xi∈ℝdX\_{i}\\in\\mathbb{R}^{d}. Given a model likelihood fθ​(y∣x)f\_{\\theta}(y\\mid x) and prior on parameters, π⁡(θ)\\pi(\\theta) for θ∈ℝp\\theta\\in\\mathbb{R}^{p}, the posterior predictive distribution for the response at a new Xn+1=xn+1X\_{n+1}=x\_{n+1} takes on the form

|  |  |  |  |
| --- | --- | --- | --- |
|  | p(y∣xn+1,Z1:n)=∫fθ(y∣xn+1)π(θ∣Z1:n)dθ,p(y\\mid x\_{n+1},Z\_{1:n})=\\int f\_{\\theta}(y\\mid x\_{n+1})\\,\\pi(\\theta\\mid Z\_{1:n})\\,d\\theta\\,, |  | (1) |

where π(θ∣Z1:n)\\pi(\\theta\\mid Z\_{1:n}) is the Bayesian posterior. Asymptotically exact samples from the posterior can be obtained through Markov chain Monte Carlo (MCMC) and the above density can be computed through Monte Carlo (MC), or by direct sampling from an approximate model. Given a Bayesian predictive distribution, one can then construct the highest density 100×(1−α)%100\\times(1-\\alpha)\\% posterior predictive credible intervals, which are the shortest intervals to contain (1−α)(1-\\alpha) of the predictive probability. Alternatively, the central 100×(1−α)%100\\times(1-\\alpha)\\% credible interval can be computed using the α/2{\\alpha}/{2} and 1−α/21-{\\alpha}/{2} quantiles. Posterior predictive distributions condition on the observed Z1:nZ\_{1:n} and represent subjective and coherent beliefs. However, it is well known that model misspecification can lead Bayesian intervals to be poorly calibrated in the frequentist sense ([Dawid, 1982](#bib.bib10); [Fraser et al., 2011](#bib.bib15)), that is the long run proportion of the observed data lying in the (1−α)(1-\\alpha) Bayes predictive interval is not necessarily equal to (1−α)(1-\\alpha). This has consequences for the robustness of such approaches and trust in using Bayesian models to aid decisions.

Alternatively, one can seek intervals around a point prediction from the model, y^=μ^​(x)\\widehat{y}=\\widehat{\\mu}(x), that have the correct frequentist coverage of (1−α)(1-\\alpha) . This is precisely what is offered by the conformal prediction framework of [Vovk et al., 2005](#bib.bib37), which allows the construction of prediction bands with finite sample validity without assumptions on the generative model beyond exchangeability of the data. Formally, for Zi={Yi,Xi}i=1:nZ\_{i}=\\{Y\_{i},X\_{i}\\}\_{i=1:n}, Zi​∼iid​ℙZ\_{i}\\overset{\\mathrm{iid}}{\\sim}\\mathbb{P} and miscoverage level α\\alpha, conformal inference allows us to construct a confidence set Cα​(Xn+1)C\_{\\alpha}(X\_{n+1}) from Z1:nZ\_{1:n} and Xn+1X\_{n+1} such that

|  |  |  |  |
| --- | --- | --- | --- |
|  | ℙ⁡(Yn+1∈Cα​(Xn+1))≥1−α\\mathbb{P}(Y\_{n+1}\\in C\_{\\alpha}(X\_{n+1}))\\geq 1-\\alpha |  | (2) |

noting that ℙ\\mathbb{P} is over Z1:n+1Z\_{1:n+1}. In this paper we develop computationally efficient conformal inference methods for Bayesian models including extensions to hierarchical settings. A general theme of our work is that, somewhat counter-intuitively, Bayesian models are well suited for the conformal method.

Conformal inference for calibrating Bayesian models was previously suggested in [Melluish et al., 2001](#bib.bib27), [Vovk et al., 2005](#bib.bib37), [Wasserman, 2011](#bib.bib38) and [Burnaev and Vovk, 2014](#bib.bib5), where it is referred to as “de-Bayesing”, “frequentizing” and “conformalizing”, but only in the context of conjugate models. Here, we present a scalable MC method for conformal Bayes, implementing full conformal Bayesian prediction using an ‘add-one-in’ importance sampling algorithm. The automated method can construct conformal predictive intervals from any Bayesian model given only samples of model parameter values from the posterior θ∼π(θ∣Z1:n)\\theta\\sim\\pi(\\theta\\mid Z\_{1:n}), up to MC error. Such samples are readily available in most Bayesian analyses from probabilistic programming languages such as Stan ([Carpenter et al., 2017](#bib.bib7)) and PyMC3 ([Salvatier et al., 2016](#bib.bib30)). We also extend conformal inference to partially exchangeable settings which utilize the important class of Bayesian hierarchical models, and note the connection to Mondrian conformal prediction ([Vovk et al., 2005](#bib.bib37), Chapter 4.5). Previously, the extension of conformal prediction to random effects was introduced in [Dunn et al., 2020](#bib.bib12) in a non-Bayesian setting, with a focus on prediction in new groups, as well as within-group predictions without covariates. We will see that the Bayesian hierarchical model allows for a natural sharing of information between groups for within-group predictions with covariates. We discuss the motivation behind using the Bayesian posterior predictive density as the conformity measure for both the Bayesian and the frequentist, and demonstrate the benefits in a number of examples.

### 1.1 Background

The conformal inference framework was first introduced by [Gammerman et al., 1998](#bib.bib16), followed by the thorough book of [Vovk et al., 2005](#bib.bib37). Full conformal prediction is computationally expensive, requiring the whole model to be retrained at each test covariate xn+1x\_{n+1} and for each value in a reference grid of potential outcomes, e.g. y∈ℝy\\in\\mathbb{R} for regression. This makes the task computationally infeasible beyond a few special cases where we can shortcut the evaluation along the outcome reference grid, e.g. ridge regression ([Vovk et al., 2005](#bib.bib37); [Burnaev and Vovk, 2014](#bib.bib5)) and lasso ([Lei, 2019](#bib.bib22)). Shrinking the search grid is possible, but still requires many refittings of the model ([Chen et al., 2016](#bib.bib9)). The split conformal prediction method ([Lei et al., 2018](#bib.bib23)) is a useful alternative method which only requires a single model fit, but increases variability by dividing the data into a training and test set that includes randomness in the choice of the split, and has a tendency for wider intervals. Methods based on cross-validation such as cross-conformal prediction ([Vovk, 2015](#bib.bib36)) and the jacknife+ ([Barber et al., 2021](#bib.bib1)) lie in between the split and full conformal method in terms of computation. A detailed discussion of computational costs of various conformal methods are provided in [Barber et al., 2021](#bib.bib1). A review of recent advances in conformal prediction is given in [Zeni et al., 2020](#bib.bib40), and interesting extensions have been developed by works such as [Tibshirani and Foygel, 2019](#bib.bib33); [Romano et al., 2019](#bib.bib29); [Candès et al., 2021](#bib.bib6).

## 2 Conformal Bayes

### 2.1 Full Conformal Prediction

We begin by summarizing the full conformal prediction algorithm discussed in [Vovk et al., 2005](#bib.bib37); [Lei et al., 2018](#bib.bib23). Firstly, a conformity (goodness-of-fit) measure,

|  |  |  |
| --- | --- | --- |
|  | σi:=σ(Z1:n+1;Zi),\\sigma\_{i}:=\\sigma(Z\_{1:n+1};Z\_{i}), |  |

takes as input a set of data points Z1:n+1Z\_{1:n+1}, and computes how similar the data point ZiZ\_{i} is for i=1,…,n+1i=1,\\ldots,n+1. A typical conformity measure for regression would be the negative squared error arising from a point prediction −{yi−μ^​(xi)}2-\\left\\{y\_{i}-\\widehat{\\mu}(x\_{i})\\right\\}^{2}, where μ^​(x)\\widehat{\\mu}(x) is the point predictor fit to the augmented dataset Z1:n+1Z\_{1:n+1}, assumed to be symmetric with respect to the permutation of the input dataset. The key property of any conformity measure is that it is exchangeable in the first argument, i.e. the conformity measure for ZiZ\_{i} is invariant to the permutation of Z1:n+1Z\_{1:n+1}. Under the assumption that Z1:n+1Z\_{1:n+1} is exchangeable, we then have that σ1:n+1\\sigma\_{1:n+1} is also exchangeable, and its rank is uniform among {1,…,n+1}\\{1,\\ldots,n+1\\} (assuming continuous σ1:n+1\\sigma\_{1:n+1}). From this, we have that the rank of σn+1\\sigma\_{n+1} is a valid pp-value. If we now consider a plug-in value Yn+1=yY\_{n+1}=y (where Xn+1X\_{n+1} is known), we can denote the rank of σn+1\\sigma\_{n+1} among σ1:n+1\\sigma\_{1:n+1} as

|  |  |  |
| --- | --- | --- |
|  | π⁡(y)=1n+1​∑i=1n+1𝟙​(σi≤σn+1).\\pi(y)=\\frac{1}{n+1}\\sum\_{i=1}^{n+1}\\mathbbm{1}\\left(\\sigma\_{i}\\leq\\sigma\_{n+1}\\right). |  |

For miscoverage level α\\alpha, the full conformal predictive set,

|  |  |  |  |
| --- | --- | --- | --- |
|  | Cα​(Xn+1)={y∈ℝ:π⁡(y)>α},C\_{\\alpha}(X\_{n+1})=\\{y\\in\\mathbb{R}:\\pi(y)>\\alpha\\}, |  | (3) |

satisfies the desired frequentist coverage as in ([2](#S1.E2 "In 1 Introduction ‣ Conformal Bayesian Computation")). Intuitively, we are reporting the values of yy which conform better than the fraction α\\alpha of observed conformity scores in the augmented dataset. A formal proof can be found in [Vovk et al., 2005](#bib.bib37). For continuous σ1:n+1\\sigma\_{1:n+1}, we also have from [Lei et al., 2018](#bib.bib23) that the conformal predictive set does not significantly over-cover.

In practice, beyond a few exceptions, the function π⁡(y)\\pi(y) must be computed on a fine grid y∈𝒴gridy\\in\\mathcal{Y}\_{\\text{grid}}, for example of size 100, in which case the model must be retrained 100 times to the augmented dataset to compute σ1:n+1\\sigma\_{1:n+1}, with plug-in values for yn+1y\_{n+1} on the grid. This is illustrated in Algorithm [1](#algorithm1 "In 2.1 Full Conformal Prediction ‣ 2 Conformal Bayes ‣ Conformal Bayesian Computation") below. We note here that the grid method only provides approximate coverage, as yy may be selected even if it lies between two grid points that are not selected. This is formalized in [Chen et al., 2018](#bib.bib8), but we do not discuss this further. In the Appendix, we provide an empirical comparison of the grid effects. This is also valid for binary classification where we now have a finite 𝒴grid={0,1}\\mathcal{Y}\_{\\text{grid}}=\\{0,1\\}, and so the grid method for full conformal prediction is exact and feasible.

Observed data is Z1:n,Xn+1Z\_{1:n},X\_{n+1}; Specify miscoverge level α\\alpha

for *each y∈𝒴gridy\\in\\mathcal{Y}\_{\\textnormal{grid}}* do

Fit model to augmented dataset {Z1,…,Zn,{y,Xn+1}}\\{Z\_{1},\\ldots,Z\_{n},\\{y,X\_{n+1}\\}\\}

Compute σ1:n\\sigma\_{1:n} and σn+1\\sigma\_{n+1}

Store the rank, π⁡(y)\\pi(y) , of σn+1\\sigma\_{n+1} among σ1:n+1\\sigma\_{1:n+1}

end for

Return the set Cα​(Xn+1)={y∈𝒴grid:π⁡(y)>α}C\_{\\alpha}(X\_{n+1})=\\{y\\in\\mathcal{Y}\_{\\text{grid}}:\\pi(y)>\\alpha\\}.

Algorithm 1 Full Conformal Prediction

### 2.2 Conformal Bayes and Add-One-In Importance Sampling

In a Bayesian model, a natural suggestion for the conformity score, as noted in [Vovk et al., 2005](#bib.bib37); [Wasserman, 2011](#bib.bib38), is the posterior predictive density ([1](#S1.E1 "In 1 Introduction ‣ Conformal Bayesian Computation")), that is

|  |  |  |
| --- | --- | --- |
|  | σ(Z1:n+1;Zi)=p(Yi∣Xi,Z1:n+1).\\sigma(Z\_{1:n+1};Z\_{i})=p(Y\_{i}\\mid X\_{i},Z\_{1:n+1}). |  |

This is a valid conformity score, as we have π(θ∣Z1:n+1)∝π(θ)∏i=1n+1fθ(Yi∣Xi),\\pi(\\theta\\mid Z\_{1:n+1})\\propto\\pi(\\theta)\\prod\_{i=1}^{n+1}f\_{\\theta}(Y\_{i}\\mid X\_{i}), and so σ\\sigma is indeed invariant to the permutation of Z1:n+1Z\_{1:n+1}. We denote this method as conformal Bayes (CB), and we will see shortly that the exchangeability structure of Bayesian models is key to constructing conformity scores in the partial exchangeability scenario.

Beyond conjugate models, we are usually able to obtain (asymptotically exact) posterior samples θ(1:T)∼π(θ∣Z1:n)\\theta^{(1:T)}\\sim\\pi(\\theta\\mid Z\_{1:n}), e.g through MCMC, where TT is a large integer. Such samples are typically available as standard output from Bayesian model fitting. The posterior predictive can then be computed up to Monte Carlo error through

|  |  |  |
| --- | --- | --- |
|  | p^(Yn+1∣Xn+1,Z1:n)=1T∑t=1Tfθ(t)(Yn+1∣Xn+1).\\widehat{p}(Y\_{n+1}\\mid X\_{n+1},Z\_{1:n})=\\frac{1}{T}\\sum\_{t=1}^{T}f\_{\\theta^{(t)}}(Y\_{n+1}\\mid X\_{n+1}). |  |

The key insight is that refitting the Bayesian model with {Z1,…,Zn,{y,Xn+1}}\\{Z\_{1},\\ldots,Z\_{n},\\{y,X\_{n+1}\\}\\} is well approximated through importance sampling (IS), as only {y,Xn+1}\\{y,X\_{n+1}\\} changes between refits. This leads immediately to an IS approach to full conformal Bayes, where we just need to compute ‘add-one-in’ (AOI) predictive densities. Here AOI refers to the inclusion of {Yn+1,Xn+1}\\{Y\_{n+1},X\_{n+1}\\} into the training set, named in relation to ‘leave-one-out’ (LOO) cross-validation. Specifically, for Yn+1=yY\_{n+1}=y and θ(1:T)∼π(θ∣Z1:n)\\theta^{(1:T)}\\sim\\pi(\\theta\\mid Z\_{1:n}), we can compute

|  |  |  |  |
| --- | --- | --- | --- |
|  | p^(Yi∣Xi,Z1:n+1)=∑t=1Tw~(t)fθ(t)(Yi∣Xi)\\widehat{p}(Y\_{i}\\mid X\_{i},Z\_{1:n+1})=\\sum\_{t=1}^{T}\\widetilde{w}^{(t)}f\_{\\theta^{(t)}}(Y\_{i}\\mid X\_{i}) |  | (4) |

where w~(t)\\widetilde{w}^{(t)} are our self-normalized importance weights of the form

|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
|  | w(t)\\displaystyle w^{(t)} | =fθ(t)​(y∣Xn+1),w~(t)\\displaystyle=f\_{\\theta^{(t)}}(y\\mid X\_{n+1}),\\quad\\widetilde{w}^{(t)} | =w(t)∑t′=1Tw(t′)⋅\\displaystyle=\\frac{w^{(t)}}{\\sum\_{t^{\\prime}=1}^{T}w^{(t^{\\prime})}}\\cdot |  | (5) |

We see that the unnormalized importance weights have the intuitive form of the predictive likelihood at the reference point {y,Xn+1}\\{y,X\_{n+1}\\} given the model parameters θ(t)\\theta^{(t)}.

The use of AOI importance sampling has similarities to the computation of Bayesian leave-one-out cross-validation (LOOCV) predictive densities ([Vehtari et al., 2017](#bib.bib34)), which is also used in accounting for model misspecification. An interesting aspect of AOI in comparison with LOO is that AOI predictive densities are less vulnerable to importance weight instability for the following reasons:

- •

  In LOOCV, the target π⁡(θ∣Z−i)\\pi(\\theta\\mid Z\_{-i}) generally has thinner tails than the proposal π(θ∣Z1:n)\\pi(\\theta\\mid Z\_{1:n}), leading to importance weight instability. In contrast, AOI uses the posterior π(θ∣Z1:n)\\pi(\\theta\\mid Z\_{1:n}) as a proposal for the thinner-tailed π(θ∣Z1:n+1)\\pi(\\theta\\mid Z\_{1:n+1}). For LOOCV the importance weights are proportional to 1/fθ​(y∣x)1/f\_{\\theta}(y\\mid x), in contrast to the typically bounded fθ​(y∣x)f\_{\\theta}(y\\mid x) for AOI.
- •

  For AOI, we are predicting ZiZ\_{i} given Z1:n+1Z\_{1:n+1} which is always in-sample unlike in LOOCV where the datum is out-of-sample, so we can expect greater stability with AOI.
- •

  The IS weight stability is governed by Yn+1=yY\_{n+1}=y, which is not random as we select it for the grid. For sufficiently large α\\alpha, we will not need to compute the AOI predictive density for extreme values of yy.

We provide some IS weight diagnostics in the experiments and find that they are stable. In difficult settings such as very high-dimensions, one can make use of the recommendations of [Vehtari et al., 2015](#bib.bib35) for assessing and Pareto-smoothing the importance weights if necessary.

### 2.3 Computational complexity

Given the posterior samples, we must compute the likelihood for each θ(t)\\theta^{(t)} at Z1:nZ\_{1:n}, as well at {y,Xn+1}\\{y,X\_{n+1}\\} for y∈𝒴gridy\\in\\mathcal{Y}\_{\\text{grid}}. The additional computation required for CB for each Xn+1X\_{n+1} is thus T×(n+ngrid)T\\times(n+n\_{\\text{grid}}) likelihood evaluations, which is relatively cheap. This is then followed by the dot product of an (n+1)×T(n+1)\\times T matrix with a TT vector for each yy, which is 𝒪⁡(n​T)\\mathcal{O}(nT), so the overall complexity is 𝒪⁡(ngrid​T​n)\\mathcal{O}(n\_{\\text{grid}}Tn). The values ngridn\_{\\text{grid}} and TT are constants, though we may want to increase TT with the dimensionality of the model to reduce importance sampling variance. The large matrices involved in computing the AOI predictives suggests we can take advantage of GPU computation, and machine learning packages such as JAX ([Bradbury et al., 2018](#bib.bib4)) are highly suitable for this application.

### 2.4 Motivation

Much has been written on the contrasting foundations and interpretation of Bayes versus frequentist measures of uncertainty ([Little, 2006](#bib.bib25); [Shafer and Vovk, 2008](#bib.bib32); [Bernardo and Smith, 2009](#bib.bib3); [Wasserman, 2011](#bib.bib38)), and we provide a summary in the Appendix. Here we motivate CB predictive intervals from both a Bayesian and frequentist perspective.

The pragmatic Bayesian, aware of the potential for model misspecification in either the prior or likelihood, may be interested in conformal inference as a countermeasure. CB predictive intervals with guaranteed frequentist coverage can be provided as a supplement to the usual Bayesian predictive intervals. The difference between the Bayesian and conformal interval may also serve as an informal diagnostic for model evaluation (e.g. [Gelman et al., 2013](#bib.bib17)). Posterior samples through MCMC or direct sampling are typically available, and so CB through automated AOI carries little overhead.

The frequentist may also wish to use a Bayesian model as a tool for constructing predictive confidence intervals. Firstly, the likelihood can take into account skewness, heteroscedasticity unlike the usual residual conformity score. Secondly, features such as sparsity, support, and regularization can be included through priors, while CB ensures correct coverage. Finally, a subtle issue that arises in full conformal prediction is that we lose validity if hyperparameter selection is not symmetric with respect to Zn+1Z\_{n+1}, e.g. if we estimate the lasso penalty λ\\lambda using only Z1:nZ\_{1:n} before computing the full conformal intervals with said λ(Z1:n)\\lambda(Z\_{1:n}). For CB, a prior on hyperparameters induces weighting of the hyperparameter values by implicit cross-validation for each refit ([Gneiting and Raftery, 2007](#bib.bib19); [Fong and Holmes, 2020](#bib.bib14)). We highlight here that this issue does not affect the split conformal method.

## 3 Partial Exchangeability and Hierarchical Models

A setting of particular interest is for grouped data, which corresponds to a weakening of exchangeability often denoted as partial exchangeability ([Bernardo and Smith, 2009](#bib.bib3), Chapter 4.6). Assume that we observe data from JJ groups, each of size njn\_{j}, where again Zi,j={Yi,j,Xi,j}Z\_{i,j}=\\{Y\_{i,j},X\_{i,j}\\}. We denote the full dataset as Z={Zi,j:i=1,…,nj,j=1,…,J}Z=\\{Z\_{i,j}:i=1,\\ldots,n\_{j},j=1,\\ldots,J\\}. We may not expect the entire sequence ZZ to be exchangeable, instead only that data points are exchangeable within groups. Formally, this means that

|  |  |  |  |
| --- | --- | --- | --- |
|  | p(Z1:n1,1,…,Z1:nJ,J)=p(Zπ1​(1):π1​(n1),1,…,ZπJ​(1):πJ​(nJ),J)p(Z\_{1:n\_{1},1},\\ldots,Z\_{1:n\_{J},J})=p(Z\_{\\pi\_{1}(1):\\pi\_{1}(n\_{1}),1},\\ldots,Z\_{\\pi\_{J}(1):\\pi\_{J}(n\_{J}),J}) |  | (6) |

for any permutations πj\\pi\_{j} of 1,…,nj1,\\ldots,n\_{j}, for j=1,…,Jj=1,\\ldots,J. Alternatively, we can enforce the usual definition of exchangeability but only consider permutations π\\pi of 1,…,n1,\\ldots,n such that the groupings are preserved. A simple example of this partial exchangeability is if Zi,j​∼iid​PjZ\_{i,j}\\overset{\\mathrm{iid}}{\\sim}P\_{j} for i=1,…,nj,j=1,…,Ji=1,\\ldots,n\_{j},\\,j=1,\\ldots,J, where PjP\_{j} can now be distinct.

Partial exchangeability is useful in multilevel modelling, e.g. where Z1:nj,jZ\_{1:n\_{j},j} records exam results on students within school jj, for schools j=1,…,Jj=1,\\ldots,J. Students may be deemed exchangeable within schools, but not between schools. Further examples may be found in [Gelman and Hill, 2006](#bib.bib18).

### 3.1 Group Conformal Prediction

Given a new Xnj+1,jX\_{n\_{j}+1,j} belonging to group jj for j∈{1,…,J}j\\in\\{1,\\ldots,J\\}, we seek to construct a (1−αj)(1-\\alpha\_{j}) confidence interval for Ynj+1,jY\_{n\_{j}+1,j}. We define a within-group conformity score as

|  |  |  |
| --- | --- | --- |
|  | σi,j:=σZ−j(Z1:nj+1,j;Zi,j)\\sigma\_{i,j}:=\\sigma\_{Z\_{-j}}(Z\_{1:n\_{j}+1,j};Z\_{i,j}) |  |

for i=1,…,nj+1i=1,\\ldots,n\_{j}+1. We denote Z−jZ\_{-j} as the dataset without group jj, and the subscript indicates the dependence of the conformity score on this, which we motivate in the next subsection. For each Z−jZ\_{-j}, we require the score to be invariant with respect to the permutation of Z1:nj+1,jZ\_{1:n\_{j}+1,j}. For Znj+1,j={y,Xnj+1,j}Z\_{n\_{j}+1,j}=\\{y,X\_{n\_{j}+1,j}\\}, the conformal predictive set is then defined

|  |  |  |  |
| --- | --- | --- | --- |
|  | πj(y)=1nj+1∑i=1nj+1𝟙(σi,j≤σnj+1,j),Cαj(Xnj+1,j)={y∈ℝ:πj(y)>αj)\\displaystyle\\pi\_{j}(y)=\\frac{1}{n\_{j}+1}\\sum\_{i=1}^{n\_{j}+1}\\mathbbm{1}\\left(\\sigma\_{i,j}\\leq\\sigma\_{n\_{j}+1,j}\\right),\\quad C\_{\\alpha\_{j}}\\left(X\_{n\_{j}+1,j}\\right)=\\{y\\in\\mathbb{R}:\\pi\_{j}(y)>\\alpha\_{j}) |  | (7) |

In other words, we rank the conformity scores σ1:nj+1,j\\sigma\_{1:n\_{j}+1,j} within the group jj, and compute the conformal interval as usual with Algorithm [1](#algorithm1 "In 2.1 Full Conformal Prediction ‣ 2 Conformal Bayes ‣ Conformal Bayesian Computation"). The interval is valid from the following.

###### Proposition 1.

Assume that {Z,Znj+1,j}\\{Z,Z\_{n\_{j}+1,j}\\} is partially exchangeable as in ([6](#S3.E6 "In 3 Partial Exchangeability and Hierarchical Models ‣ Conformal Bayesian Computation")), and the conformity measure σi,j\\sigma\_{i,j} for group jj is invariant to the permutation of Z1:nj+1,jZ\_{1:n\_{j}+1,j}. We then have

|  |  |  |
| --- | --- | --- |
|  | ℙ⁡(Ynj+1,j∈Cαj​(Xnj+1,j))≥1−αj\\mathbb{P}\\left(Y\_{n\_{j}+1,j}\\in C\_{\\alpha\_{j}}\\left(X\_{n\_{j}+1,j}\\right)\\right)\\geq 1-\\alpha\_{j} |  |

where Cαj​(Xnj+1,j)C\_{\\alpha\_{j}}\\left(X\_{n\_{j}+1,j}\\right) is defined in ([7](#S3.E7 "In 3.1 Group Conformal Prediction ‣ 3 Partial Exchangeability and Hierarchical Models ‣ Conformal Bayesian Computation")), and ℙ\\mathbb{P} is over {Z,Znj+1,j}\\{Z,Z\_{n\_{j}+1,j}\\}.

###### Proof.

Conditional on Z−jZ\_{-j}, the observations Z1:nj+1,jZ\_{1:n\_{j}+1,j} are still exchangeable, and thus so are σ1:nj+1,j\\sigma\_{1:n\_{j}+1,j} from the invariance of the conformity measure. The usual conformal guarantee then holds:

|  |  |  |
| --- | --- | --- |
|  | ℙ⁡(Ynj+1,j∈Cαj​(Xnj+1,j)∣Z−j)≥1−αj.\\mathbb{P}\\left(Y\_{n\_{j}+1,j}\\in C\_{\\alpha\_{j}}\\left(X\_{n\_{j}+1,j}\\right)\\mid Z\_{-j}\\right)\\geq 1-\\alpha\_{j}. |  |

Taking the expectation with respect to Z−jZ\_{-j} gives us the result. ∎

It is interesting to note that the above group conformal predictor coincides with the attribute-conditional Mondrian conformal predictor of [Vovk et al., 2005](#bib.bib37), with the group allocations as the taxonomy. Validity under the relaxed Mondrian-exchangeability of [Vovk et al., 2005](#bib.bib37) is key for us here.

### 3.2 Conformal Hierarchical Bayes

Under this setting, a hierarchical Bayesian model can be defined of the form

|  |  | ∼iidfθj,τ(⋅∣Xi,j)i=1,…,nj,\\displaystyle\\overset{\\mathrm{iid}}{\\sim}f\_{\\theta\_{j},\\tau}(\\cdot\\mid X\_{i,j})\\quad i=1,\\ldots,n\_{j},\\quad | j=1,…,J\\displaystyle j=1,\\ldots,J |  |
| --- | --- | --- | --- | --- |
|  | \[θj∣ϕ\]\\displaystyle\[\\theta\_{j}\\mid\\phi\] | ∼iidπ(⋅∣ϕ)\\displaystyle\\overset{\\mathrm{iid}}{\\sim}\\pi(\\cdot\\mid\\phi)\\quad | j=1,…,J\\displaystyle j=1,\\ldots,J |  |
|  | ϕ∼π⁡(ϕ),\\displaystyle\\phi\\sim\\pi(\\phi), | τ∼π⁡(τ).\\displaystyle\\tau\\sim\\pi(\\tau). |  |  |

Here τ\\tau is a common parameter across groups (e.g. a common standard deviation for the residuals under homoscedastic errors). The desired partial exchangeability structure is clearly preserved in the Bayesian model ([Bernardo, 1996](#bib.bib2)). De Finetti representation theorems are also available for partially exchangeable sequences (when defined in a slightly different manner to the above), which motivate the specification of hierarchical Bayesian models ([Bernardo and Smith, 2009](#bib.bib3), Chapter 4.6).

The posterior predictive is once again a natural choice for the conformity measure. Denoting Z¯y\\bar{Z}\_{y} as the entire dataset augmented with Znj+1,j={y,Xnj+1,j}Z\_{n\_{j}+1,j}=\\{y,X\_{n\_{j}+1,j}\\}, we have

|  |  |  |  |
| --- | --- | --- | --- |
|  | σi,j=p⁡(Yi,j∣Xi,j,Z¯y)=∫fθj,τ​(Yi,j∣Xi,j)​π​(θj,τ∣Z¯y)​d​θj​𝑑τ\\sigma\_{i,j}=p(Y\_{i,j}\\mid X\_{i,j},\\bar{Z}\_{y})=\\int f\_{\\theta\_{j},\\tau}(Y\_{i,j}\\mid X\_{i,j})\\,\\pi(\\theta\_{j},\\tau\\mid\\bar{Z}\_{y})\\,d\\theta\_{j}\\,d\\tau |  | (8) |

for i=1,…,nj+1i=1,\\ldots,n\_{j}+1. The within-group permutation invariance follows as the likelihood is exchangeable within groups, and thus so is the posterior and resulting posterior predictive. Practically, this structure allows for independent coefficients θj\\theta\_{j} for each group, but partial pooling through π⁡(θ∣ϕ)\\pi(\\theta\\mid\\phi) allows information to be shared between groups. A fully pooled model, whilst still valid, is usually too simple and predicts poorly, whereas a no-pooling conformity score ignores information sharing between groups. More details on hierarchical models can be found in [Gelman et al., 2013](#bib.bib17). We point out that we can select a separate coverage level αj\\alpha\_{j} for each group, which will be useful when group sizes njn\_{j} vary - we provide a demonstration of this in the Appendix. Computation of σi,j\\sigma\_{i,j} is again straightforward, where MCMC now returns \[θ1:J(1:T),ϕ(1:T),τ(1:T)\]∼π(θ1:J,ϕ,τ∣Z)\[\\theta\_{1:J}^{(1:T)},\\phi^{(1:T)},\\tau^{(1:T)}\]\\sim\\pi(\\theta\_{1:J},\\phi,\\tau\\mid Z). We can then estimate ([8](#S3.E8 "In 3.2 Conformal Hierarchical Bayes ‣ 3 Partial Exchangeability and Hierarchical Models ‣ Conformal Bayesian Computation")) using AOI importance sampling as in ([4](#S2.E4 "In 2.2 Conformal Bayes and Add-One-In Importance Sampling ‣ 2 Conformal Bayes ‣ Conformal Bayesian Computation")) and ([5](#S2.E5 "In 2.2 Conformal Bayes and Add-One-In Importance Sampling ‣ 2 Conformal Bayes ‣ Conformal Bayesian Computation")) using the marginal samples {θj(1:T),τ(1:T)}∼π(θj,τ∣Z)\\{\\theta\_{j}^{(1:T)},\\tau^{(1:T)}\\}\\sim\\pi(\\theta\_{j},\\tau\\mid Z) and weights w(t)=fθj(t),τ(t)​(y∣Xnj+1,j)w^{(t)}=f\_{\\theta\_{j}^{(t)},\\tau^{(t)}}(y\\mid X\_{n\_{j}+1,j}).

In the above, we consider predictive intervals within groups with covariates, extending the within-group approach of [Dunn et al., 2020](#bib.bib12). Predictive intervals for new groups are possible with the Bayesian model, but a conformal predictor would require additional stronger assumptions of exchangeability to ensure validity. The relevant assumptions and methods based on pooling and subsampling for new group predictions are discussed in [Dunn et al., 2020](#bib.bib12), but would require rerunning MCMC in our case. We leave this for future work, noting that utilizing the Bayesian predictive density directly here seems nontrivial due to different group sizes.

## 4 Experiments

We run and time all examples on an Azure NC6 Virtual Machine, which has 6 Intel Xeon E5-2690 v3 vCPUs and a one-half Tesla K80 GPU card. We use PyMC3 ([Salvatier et al., 2016](#bib.bib30)) for MCMC and sklearn ([Pedregosa et al., 2011](#bib.bib28)) for the regular conformal predictor; both are run on the CPU. Computation of the CB and Bayes intervals is implemented in JAX ([Bradbury et al., 2018](#bib.bib4)), and run on the GPU. The code is available online11 1 [https://github.com/edfong/conformal\_bayes](https://github.com/edfong/conformal_bayes) and further examples are provided in the Appendix.

### 4.1 Sparse Regression

We first demonstrate our method under a sparse linear regression model on the diabetes dataset ([Efron et al., 2004](#bib.bib13)) considered by [Lei, 2019](#bib.bib22). The dataset is available in sklearn, and consists of n=442n=442 subjects, where the response variable is a continuous diabetes progression and the d=10d=10 covariates consist of patient readings such as blood serum measurements. We standardize all covariates and the response to have mean 0 and standard deviation 1.

The Bayesian model we consider is

|  | fθ​(y∣x)=𝒩⁡(y∣θT​x+θ0,τ2)\\displaystyle f\_{\\theta}(y\\mid x)=\\mathcal{N}(y\\mid\\theta^{\\mathrm{\\scriptscriptstyle T}}x+\\theta\_{0},\\tau^{2}) |  | (9) |
| --- | --- | --- | --- |
|  | π⁡(θj)=Laplace​(0,b),π⁡(θ0)∝1,π⁡(b)=Gamma​(1,1)\\displaystyle\\pi(\\theta\_{j})=\\text{Laplace}(0,b),\\quad\\pi(\\theta\_{0})\\propto 1,\\quad\\pi(b)=\\text{Gamma}(1,1) | π⁡(τ)=𝒩+​(0,c)\\displaystyle\\pi(\\tau)=\\mathcal{N}^{+}(0,c) |  |

for j=1,…,dj=1,\\ldots,d, and where bb is the scale parameter and 𝒩+\\mathcal{N}^{+} is the half-normal distribution. Note that a hyperprior on bb has removed the need for cross-validation that is required for lasso. We consider two values of cc for the hyperprior on τ\\tau, which correspond to a well-specified (c=1c=1) and poorly-specified (c=0.02)(c=0.02) prior; in the latter case our posterior on τ\\tau will be heavily weighted towards a small value. This model is well-specified for the diabetes dataset ([Jansen, 2013](#bib.bib21), Chapter 4.5) under a reasonable prior (c=1c=1). We compute the central (1−α)(1-\\alpha) credible interval from the Bayesian posterior predictive CDF estimated using Monte Carlo and the same grid as for CB.

To check coverage, we repeatedly divide into a training and test dataset for 5050 repeats, with 30%30\\% of the dataset in the test split. We evaluate the conformal prediction set on a grid of size ngrid=100n\_{\\text{grid}}=100 between \[ymin−2,ymax+2\]\[y\_{\\text{min}}-2,y\_{\\text{max}}+2\], where ymin,ymaxy\_{\\text{min}},y\_{\\text{max}} is computed from each training dataset. The average coverage, length and run-times (excluding MCMC) with standard errors are given in Table [1](#S4.T1 "Table 1 ‣ 4.1 Sparse Regression ‣ 4 Experiments ‣ Conformal Bayesian Computation") for α=0.2\\alpha=0.2. MCMC induced an average overhead of 21.921.9s for a=1a=1 and 26.826.8s for c=0.02c=0.02 for the Bayes and CB interval, where we simulate T=8000T=8000 posterior samples. The CB intervals are only slightly slower than the Bayes intervals, and still a small fraction of the time required for MCMC, and is thus an efficient post-processing step. For c=1c=1, the Bayesian intervals have coverage close to (1−α)(1-\\alpha) with the smallest expected length, with CB slightly wider and more conservative. However, when the prior is misspecified with c=0.02c=0.02, the Bayes intervals severely undercover, whilst the CB coverage and length remain unchanged from the c=1c=1 case.

As baselines, we compare to the split and full conformal method using the non-Bayesian lasso as the predictor, with the usual residual as the nonconformity score. For the split method, we fit lasso with cross-validation on the subset of size ntrain/2n\_{\\text{train}}/2 to obtain the lasso penalty λ\\lambda. For the full conformal method, we use the grid method for fair timing, as other estimators beyond lasso would not have the shortcut of [Lei, 2019](#bib.bib22). As setting a default λ=1\\lambda=1 gives poor average lengths, we estimate λ=0.004\\lambda=0.004 on cross-validation on one of the training sets, and use this value over the 50 repeats. However, we must emphasize again that this is somewhat misleading, as discussed in Section [2.4](#S2.SS4 "2.4 Motivation ‣ 2 Conformal Bayes ‣ Conformal Bayesian Computation"). A fairer approach would involve fitting lasso with CV for each of the 100 grid values and 133 test values, but this is infeasible as each fit requires around 80ms, resulting in a total run-time of 17 minutes. On the other hand, the AOI scheme of CB is equivalent to refitting bb for each grid/test value. In terms of performance, the split method has wider intervals than CB/full, but performs well given the low computational costs. The full conformal method performs as well as CB, but is comparable in time as MCMC + CB, whilst not refitting λ\\lambda. We note that the value of cc does not affect the split/full method.

Table 1: Diabetes; Coverage values not within 3 standard errors (in brackets) of the target coverage (1−α)=0.8(1-\\alpha)=0.8 are in red.

|  |  | Bayes | CB | Split | Full (λ=0.004\\lambda=0.004) |
| --- | --- | --- | --- | --- | --- |
| Coverage | c=1c=1 | 0.806 (0.005) | 0.808 (0.006) | 0.816 (0.006) | 0.808 (0.006) |
|  | c=0.02c=0.02 | 0.563 (0.006) | 0.809 (0.006) | / | / |
| Length | c=1c=1 | 1.84 (0.01) | 1.87 (0.01) | 1.95 (0.02) | 1.86 (0.01) |
|  | c=0.02c=0.02 | 1.14 (0.00) | 1.87 (0.01) | / | / |
| Run-time | c=1c=1 | 0.488 (0.107) | 0.702 (0.019) | 0.065 (0.001) | 11.529 (0.232) |
| (secs) | c=0.02c=0.02 | 0.373 (0.002) | 0.668 (0.003) | / | / |

#### 4.1.1 Importance weights

For the diabetes dataset, we look at the effective sample size (ESS) of the self-normalized importance weights ([5](#S2.E5 "In 2.2 Conformal Bayes and Add-One-In Importance Sampling ‣ 2 Conformal Bayes ‣ Conformal Bayesian Computation")), which can be computed as ESS=1/∑t=1T{w(t)}2\\text{ESS}=1/\\sum\_{t=1}^{T}\\{w^{(t)}\\}^{2} for each xn+1x\_{n+1} and yy. The ESS as a function of yy for a single xn+1x\_{n+1} is shown in Figure [1](#S4.F1 "Figure 1 ‣ 4.1.1 Importance weights ‣ 4.1 Sparse Regression ‣ 4 Experiments ‣ Conformal Bayesian Computation") for the two cases c=1,0.02c=1,0.02, with the CB conformal bands given for α=0.2,0.5\\alpha=0.2,0.5. We have scaled the ESS plots by ESSMCMC/T\\text{ESS}\_{\\text{MCMC}}/T, where T=8000T=8000 is the number of posterior samples and ESSMCMC\\text{ESS}\_{\\text{MCMC}} is the minimum ESS out of all posterior parameters return by PyMC3. We observe the ESS is well behaved and stable across the range of yy values. In both cases, the ESS for α=0.2\\alpha=0.2 is sufficiently large for a reliable estimate of the conformity scores. However, for c=0.02c=0.02, the ESS decays more quickly with yy as the Bayes predictive intervals are too narrow, which the CB corrects for. Other values of xn+1x\_{n+1} produce similar behaviour.

Figure 1: Effective sample sizes of IS weights with CB conformal bands for diabetes dataset with (left) c=1c=1 and (right) c=0.02c=0.02.

### 4.2 Sparse Classification

In this section, we analyze the Wisconsin breast cancer ([Wolberg and Mangasarian, 1990](#bib.bib39)), again available in sklearn. The dataset is of size 569569, where the binary response variable corresponds to a malignant or benign tumour. The 3030 covariates consist of measurements of cell nuclei. Again, we standardize all covariates to have mean 0 and standard deviation 1.

We consider the logistic likelihood fθ​(y=1∣x)=\[1+exp⁡{−(θT​x+θ0)}\]−1,f\_{\\theta}(y=1\\mid x)=\[1+\\exp\\left\\{-\\left(\\theta^{\\mathrm{\\scriptscriptstyle T}}x+\\theta\_{0}\\right)\\right\\}\]^{-1}, with the same priors for θ,θ0\\theta,\\theta\_{0} as in ([9](#S4.E9 "In 4.1 Sparse Regression ‣ 4 Experiments ‣ Conformal Bayesian Computation")). The Bayesian predictive set is the smallest set from {0},{1},{0,1}\\{0\\},\\{1\\},\\{0,1\\} that contains at least (1−α)(1-\\alpha) of the posterior predictive probability. The conformal baselines are as above but with L1L\_{1}-penalized logistic regression, and for the full conformal method we have λ=1\\lambda=1. We again have 50 repeats with 70-30 train-test split, and set α=0.2\\alpha=0.2. The grid method is now exact, and the size of the CB intervals can take on the values {0,1,2}\\{0,1,2\\}. The results are provided in Table [2](#S4.T2 "Table 2 ‣ 4.2 Sparse Classification ‣ 4 Experiments ‣ Conformal Bayesian Computation"), where MCMC required an average of 45.445.4s to produce T=8000T=8000 samples. We see that even with reasonable priors, Bayes can over-cover substantially, which CB corrects in roughly the same amount of time as it takes to compute the usual Bayes interval. However, we point out that CB may produce empty prediction sets, whereas Bayes cannot, and we investigate this in the Appendix.

Table 2: Breast Cancer; Coverage values not within 3 standard errors (in brackets) of the target coverage (1−α)=0.8(1-\\alpha)=0.8 are in red. “Size” denotes the average number of elements in the conformal prediction set, averaged over the test points and repetitions.

|  | Bayes | CB | Split | Full |
| --- | --- | --- | --- | --- |
| Coverage | 0.990 (0.001) | 0.812 (0.005) | 0.814 (0.006) | 0.811 (0.005) |
| Size | 1.06 (0.00) | 0.81 (0.00) | 0.82 (0.01) | 0.81 (0.00) |
| Run-time (secs) | 0.364 (0.007) | 0.665 (0.012) | 0.079 (0.002) | 1.008 (0.016) |

### 4.3 Hierarchical Model

We now demonstrate Bayesian conformal inference using a hierarchical Bayesian model for multilevel data. We stick to the varying intercept and varying slope model ([Gelman et al., 2013](#bib.bib17)), that is for j=1,…,Jj=1,\\ldots,J:

|  | fθj,τ​(yi,j)=𝒩⁡(yi,j∣θjT​Xi,j+θ0,j,τ2)\\displaystyle f\_{\\theta\_{j},\\tau}(y\_{i,j})=\\mathcal{N}(y\_{i,j}\\mid\\theta^{\\mathrm{\\scriptscriptstyle T}}\_{j}X\_{i,j}+\\theta\_{0,j},\\tau^{2}) |  | (10) |
| --- | --- | --- | --- |
|  | π⁡(θj)=𝒩⁡(ϕ,s2),π⁡(θ0,j)=𝒩⁡(ϕ0,s02)\\displaystyle\\pi(\\theta\_{j})=\\mathcal{N}(\\phi,s^{2}),\\quad\\pi(\\theta\_{0,j})=\\mathcal{N}(\\phi\_{0},s\_{0}^{2}) |  |  |

with hyperpriors 𝒩⁡(0,1)\\mathcal{N}(0,1) on the location parameters ϕ,ϕ0\\phi,\\phi\_{0} and Exp​(1)\\text{Exp}(1) on the standard deviations s,s0,τs,s\_{0},\\tau. We now apply this to a simulated example, and an application to the radon dataset of [Gelman and Hill, 2006](#bib.bib18) is given in the Appendix.

We consider two simulation scenarios, with J=5J=5 groups and nj=10n\_{j}=10 elements per group:

1. 1.

   Well-specified: We generate group slopes θj​∼iid​𝒩​(0,1)\\theta\_{j}\\overset{\\mathrm{iid}}{\\sim}\\mathcal{N}(0,1) for j=1,…,Jj=1,\\ldots,J. For each jj, we generate Xi,j∼𝒩⁡(0,1)X\_{i,j}\\sim\\mathcal{N}(0,1) and Yi,j∼𝒩⁡(θj​Xi,j,1)Y\_{i,j}\\sim\\mathcal{N}\\left(\\theta\_{j}X\_{i,j},1\\right).
2. 2.

   Misspecified: We generate group slopes and variances θj​∼iid​𝒩​(0,1),τj​∼iid​Exp​(1)\\theta\_{j}\\overset{\\mathrm{iid}}{\\sim}\\mathcal{N}(0,1),\\tau\_{j}\\overset{\\mathrm{iid}}{\\sim}\\text{Exp}(1) for j=1,…,Jj=1,\\ldots,J. For each jj, we generate Xi,j∼𝒩⁡(0,1)X\_{i,j}\\sim\\mathcal{N}(0,1) and Yi,j∼𝒩⁡(θj​Xi,j,τj2)Y\_{i,j}\\sim\\mathcal{N}\\left(\\theta\_{j}X\_{i,j},\\tau\_{j}^{2}\\right).

The first scenario has homoscedastic noise between groups as assumed in the model ([10](#S4.E10 "In 4.3 Hierarchical Model ‣ 4 Experiments ‣ Conformal Bayesian Computation")) whereas the second scenario is heteroscedastic between groups. To evaluate coverage, we only draw θ1:J,τ1:J\\theta\_{1:J},\\tau\_{1:J} once (and not per repeat), giving us the values

|  |  |  |
| --- | --- | --- |
|  | θ1:J=\[1.33,−0.77,−0.32,−0.99,−1.07\],τ1:J=\[1.24,2.30,0.76,0.28,1.11\].\\theta\_{1:J}=\[1.33,-0.77,-0.32,-0.99,-1.07\],\\quad\\tau\_{1:J}=\[1.24,2.30,0.76,0.28,1.11\]. |  |

For each of the 5050 repeats, we draw nj=10n\_{j}=10 training and test data points from each group using the above θ1:J\\theta\_{1:J} (and τ1:J\\tau\_{1:J} for scenario 2), and report test coverage and lengths within each group. We use a grid of size 100 between \[−10,10\]\[-10,10\]. The group-wise average lengths and coverage are given in Table [3](#S4.T3 "Table 3 ‣ 4.3 Hierarchical Model ‣ 4 Experiments ‣ Conformal Bayesian Computation") again with α=0.2\\alpha=0.2. Again run-times are given post-MCMC, where MCMC required an average of 90.190.1s and 78.478.4s for scenarios 1 and 2 respectively to generate T=8000T=8000 samples. The Bayes interval is again the central (1−α)(1-\\alpha) credible interval. The CB and Bayes methods have comparable run-times, likely due to the small nn. As a reference, fitting a linear mixed-effects model in statsmodels ([Seabold and Perktold, 2010](#bib.bib31)) to the dataset takes around 200ms, so the full conformal method, which requires refitting for each of the 100 grid value and 50 test values, would take a total of 17 minutes. For scenario 1, both Bayes and CB provide close to (1−α)(1-\\alpha) coverage, with the Bayes lengths being smaller. This is unsurprising, as the Bayesian model is well-specified. In scenario 2, the Bayes intervals noticeably over/under-cover depending on the value of τ1:J\\tau\_{1:J} in relation to the Bayes posterior mean τ¯≈1.3\\bar{\\tau}\\approx 1.3. CB is robust to this, adapting its interval lengths accordingly (in particular for Groups 2 and 4) and providing within-group validity.

Table 3: Simulated grouped dataset; Coverage values not within 3 standard errors (in brackets) of the target coverage (1−α)=0.8(1-\\alpha)=0.8 are in red.

|  |  | Scenario 1 | Scenario 2 |  |  |
| --- | --- | --- | --- | --- | --- |
|  | Group | Bayes | CB | Bayes | CB |
| Coverage | 1 | 0.808 (0.020) | 0.794 (0.022) | 0.826 (0.020) | 0.786 (0.025) |
|  | 2 | 0.800 (0.019) | 0.812 (0.024) | 0.522 (0.027) | 0.812 (0.024) |
|  | 3 | 0.824 (0.017) | 0.824 (0.022) | 0.974 (0.008) | 0.824 (0.020) |
|  | 4 | 0.786 (0.017) | 0.798 (0.022) | 1.000 (0.000) | 0.836 (0.021) |
|  | 5 | 0.772 (0.019) | 0.810 (0.020) | 0.826 (0.022) | 0.796 (0.022) |
|  | Overall | 0.798 (0.009) | 0.808 (0.009) | 0.830 (0.010) | 0.811 (0.009) |
| Length | 1 | 2.80 (0.05) | 3.19 (0.13) | 3.65 (0.08) | 4.01 (0.17) |
|  | 2 | 2.76 (0.05) | 3.21 (0.15) | 3.61 (0.08) | 7.27 (0.33) |
|  | 3 | 2.75 (0.04) | 3.07 (0.13) | 3.59 (0.08) | 2.28 (0.09) |
|  | 4 | 2.75 (0.05) | 3.05 (0.12) | 3.57 (0.08) | 1.23 (0.04) |
|  | 5 | 2.78 (0.05) | 3.14 (0.11) | 3.61 (0.08) | 3.47 (0.12) |
|  | Overall | 2.77 (0.04) | 3.13 (0.06) | 3.61(0.08) | 3.65 (0.09) |
| Run-time (secs) | Overall | 0.222 (0.002) | 0.381 (0.009) | 0.221 (0.002) | 0.375 (0.002) |

## 5 Discussion

In this work, we have introduced the AOI importance sampling scheme for conformal Bayesian computation, which allow us to construct frequentist-valid predictive intervals from a baseline Bayesian model using the output of an MCMC sampler. This extends naturally to the partially exchangeable setting and hierarchical Bayesian models.

Under model misspecification, or the ℳ\\mathcal{M}-open scenario ([Bernardo and Smith, 2009](#bib.bib3)), CB can produce calibrated intervals from the Bayesian model. In the partially exchangeable case, CB can remain valid within groups. We find that even under reasonable priors, Bayesian predictives can over-cover, and CB can help reduce the length of intervals to get closer to nominal coverage. Diagnosing Bayesian miscalibration is in general non-trivial, but CB automatically corrects for this. When posterior samples of model parameters are available, AOI importance sampling is only a minor increase in computation, and interestingly is much faster than the split method which would require another run of MCMC. For the frequentist, CB intervals enjoy the tightness of the full conformal method, for a single expensive fit with MCMC followed by a cheap refitting process. We are also free to incorporate prior information, and use more complex likelihoods or priors, as well as automatically fitting hyperparameters.

There are however limitations to our approach, dictated by the realities of MCMC and IS. Firstly, the intervals are approximate up to MC error and reliant on representative MC samples not disrupting exchangeability of the conformity scores. The stability of AOI importance sampling also depends on the posterior predictive being a good proposal, which may break down if the addition of the new datum {y,Xn+1}\\{y,X\_{n+1}\\} has very high leverage on the posterior.

If only approximate posterior samples are available, e.g. through variational Bayes (VB), then an AOI scheme may still be feasible, where one includes an additional correction term in the IS weights for the VB approximation, e.g. in [Magnusson et al., 2019](#bib.bib26). However, this remains to be investigated. Combining this with the Pareto-smoothed IS method of [Vehtari et al., 2015](#bib.bib35) may lead to additional scalability with dimensionality. In our experience, CB intervals tend to be a single connected interval, which may allow for computational shortcuts in adapting the search grid. It would also be interesting to pursue the theoretical connections between the Bayesian and CB intervals, in a similar light to [Burnaev and Vovk, 2014](#bib.bib5).

## Acknowledgements

Fong is funded by The Alan Turing Institute Doctoral Studentship, under the EPSRC grant EP/N510129/1. Holmes is supported by The Alan Turing Institute, the Health Data Research, U.K., the Li Ka Shing Foundation, the Medical Research Council, and the U.K. Engineering and Physical Sciences Research Council through the Bayes4Health programme grant.

## References

- Barber et al., (2021)  Barber, R. F., Candes, E. J., Ramdas, A., Tibshirani, R. J., et al. (2021).  Predictive inference with the jackknife+.  Annals of Statistics, 49(1):486–507.
- Bernardo, (1996)  Bernardo, J. M. (1996).  The concept of exchangeability and its applications.  Far East Journal of Mathematical Sciences, 4:111–122.
- Bernardo and Smith, (2009)  Bernardo, J. M. and Smith, A. F. (2009).  Bayesian theory, volume 405.  John Wiley & Sons.
- Bradbury et al., (2018)  Bradbury, J., Frostig, R., Hawkins, P., Johnson, M. J., Leary, C., Maclaurin, D., and Wanderman-Milne, S. (2018).  JAX: composable transformations of Python+NumPy programs.
- Burnaev and Vovk, (2014)  Burnaev, E. and Vovk, V. (2014).  Efficiency of conformalized ridge regression.  In Conference on Learning Theory, pages 605–622.
- Candès et al., (2021)  Candès, E. J., Lei, L., and Ren, Z. (2021).  Conformalized survival analysis.  arXiv preprint arXiv:2103.09763.
- Carpenter et al., (2017)  Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., Brubaker, M. A., Guo, J., Li, P., and Riddell, A. (2017).  Stan: a probabilistic programming language.  Grantee Submission, 76(1):1–32.
- Chen et al., (2018)  Chen, W., Chun, K.-J., and Barber, R. F. (2018).  Discretized conformal prediction for efficient distribution-free inference.  Stat, 7(1):e173.
- Chen et al., (2016)  Chen, W., Wang, Z., Ha, W., and Barber, R. F. (2016).  Trimmed conformal prediction for high-dimensional models.  arXiv preprint arXiv:1611.09933.
- Dawid, (1982)  Dawid, A. P. (1982).  The well-calibrated Bayesian.  Journal of the American Statistical Association, 77(379):605–610.
- Dua and Graff, (2017)  Dua, D. and Graff, C. (2017).  UCI machine learning repository.
- Dunn et al., (2020)  Dunn, R., Wasserman, L., and Ramdas, A. (2020).  Distribution-free prediction sets with random effects.  arXiv preprint arXiv:1809.07441v2.
- Efron et al., (2004)  Efron, B., Hastie, T., Johnstone, I., Tibshirani, R., et al. (2004).  Least angle regression.  Annals of statistics, 32(2):407–499.
- Fong and Holmes, (2020)  Fong, E. and Holmes, C. (2020).  On the marginal likelihood and cross-validation.  Biometrika, 107(2):489–496.
- Fraser et al., (2011)  Fraser, D. A. et al. (2011).  Is Bayes posterior just quick and dirty confidence?  Statistical Science, 26(3):299–316.
- Gammerman et al., (1998)  Gammerman, A., Vovk, V., and Vapnik, V. (1998).  Learning by transduction.  In Proceedings of the Fourteenth conference on Uncertainty in artificial intelligence, pages 148–155.
- Gelman et al., (2013)  Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., and Rubin, D. B. (2013).  Bayesian data analysis.  CRC press.
- Gelman and Hill, (2006)  Gelman, A. and Hill, J. (2006).  Data analysis using regression and multilevel/hierarchical models.  Cambridge university press.
- Gneiting and Raftery, (2007)  Gneiting, T. and Raftery, A. E. (2007).  Strictly proper scoring rules, prediction, and estimation.  Journal of the American Statistical Association, 102(477):359–378.
- Harrison Jr and Rubinfeld, (1978)  Harrison Jr, D. and Rubinfeld, D. L. (1978).  Hedonic housing prices and the demand for clean air.  Journal of environmental economics and management, 5(1):81–102.
- Jansen, (2013)  Jansen, L. (2013).  Robust Bayesian inference under model misspecification.  PhD thesis, Master’s thesis, Leiden University.
- Lei, (2019)  Lei, J. (2019).  Fast exact conformalization of the lasso using piecewise linear homotopy.  Biometrika, 106(4):749–764.
- Lei et al., (2018)  Lei, J., G’Sell, M., Rinaldo, A., Tibshirani, R. J., and Wasserman, L. (2018).  Distribution-free predictive inference for regression.  Journal of the American Statistical Association, 113(523):1094–1111.
- Little et al., (2008)  Little, M., McSharry, P., Hunter, E., Spielman, J., and Ramig, L. (2008).  Suitability of dysphonia measurements for telemonitoring of parkinson’s disease.  Nature Precedings, pages 1–1.
- Little, (2006)  Little, R. J. (2006).  Calibrated Bayes: a Bayes/frequentist roadmap.  The American Statistician, 60(3):213–223.
- Magnusson et al., (2019)  Magnusson, M., Andersen, M., Jonasson, J., and Vehtari, A. (2019).  Bayesian leave-one-out cross-validation for large data.  In International Conference on Machine Learning, pages 4244–4253. PMLR.
- Melluish et al., (2001)  Melluish, T., Saunders, C., Nouretdinov, I., and Vovk, V. (2001).  Comparing the Bayes and typicalness frameworks.  In European Conference on Machine Learning, pages 360–371. Springer.
- Pedregosa et al., (2011)  Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., Blondel, M., Prettenhofer, P., Weiss, R., Dubourg, V., Vanderplas, J., Passos, A., Cournapeau, D., Brucher, M., Perrot, M., and Duchesnay, E. (2011).  Scikit-learn: Machine learning in Python.  Journal of Machine Learning Research, 12:2825–2830.
- Romano et al., (2019)  Romano, Y., Patterson, E., and Candès, E. J. (2019).  Conformalized quantile regression.  arXiv preprint arXiv:1905.03222.
- Salvatier et al., (2016)  Salvatier, J., Wiecki, T. V., and Fonnesbeck, C. (2016).  Probabilistic programming in Python using PyMC3.  PeerJ Computer Science, 2:e55.
- Seabold and Perktold, (2010)  Seabold, S. and Perktold, J. (2010).  Statsmodels: Econometric and statistical modeling with python.  In Proceedings of the 9th Python in Science Conference, volume 57, page 61. Austin, TX.
- Shafer and Vovk, (2008)  Shafer, G. and Vovk, V. (2008).  A tutorial on conformal prediction.  Journal of Machine Learning Research, 9(Mar):371–421.
- Tibshirani and Foygel, (2019)  Tibshirani, R. and Foygel, R. (2019).  Conformal prediction under covariate shift.  Advances in neural information processing systems.
- Vehtari et al., (2017)  Vehtari, A., Gelman, A., and Gabry, J. (2017).  Practical Bayesian model evaluation using leave-one-out cross-validation and waic.  Statistics and computing, 27(5):1413–1432.
- Vehtari et al., (2015)  Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2015).  Pareto smoothed importance sampling.  arXiv preprint arXiv:1507.02646.
- Vovk, (2015)  Vovk, V. (2015).  Cross-conformal predictors.  Annals of Mathematics and Artificial Intelligence, 74(1):9–28.
- Vovk et al., (2005)  Vovk, V., Gammerman, A., and Shafer, G. (2005).  Algorithmic learning in a random world.  Springer Science & Business Media.
- Wasserman, (2011)  Wasserman, L. (2011).  Frasian inference.  Statistical Science, pages 322–325.
- Wolberg and Mangasarian, (1990)  Wolberg, W. H. and Mangasarian, O. L. (1990).  Multisurface method of pattern separation for medical diagnosis applied to breast cytology.  Proceedings of the national academy of sciences, 87(23):9193–9196.
- Zeni et al., (2020)  Zeni, G., Fontana, M., and Vantini, S. (2020).  Conformal prediction: a unified review of theory and new challenges.  arXiv preprint arXiv:2005.07972.

## Appendix A Bayesian and Frequentist Intervals

Bayesian predictive intervals are conditioned on the specific observed sequence Z1:nZ\_{1:n} and make statements on the next value \[Yn+1∣Xn+1\]\[Y\_{n+1}\\mid X\_{n+1}\]. In contrast, a conformal (frequentist) interval relates to the properties of intervals returned by the algorithm if run repeatedly across different data sets of size nn. Subjective Bayesian statements on predictions are non-refutable, and are in this sense unscientific, but are optimal according to decision theoretic foundations. Meanwhile, frequentist statements are in principle verifiable and hence refutable.

To us, in the hands of an expert analyst with careful prior elicitation, the Bayesian conditional argument is the more persuasive for posterior and predictive uncertainty. The Bayesian predictive provides statements of uncertainty conditional on what has been observed, and so decisions pertain to each specific dataset. However, to make such strong statements, the Bayesian must usually make the strict assumption of the model being well-specified. If we wish to ensure that the predictive coverage of reported intervals is calibrated on average across repeats under weaker assumptions, then the conformal intervals are much more suitable. More details contrasting probability and confidence can be found in [Shafer and Vovk, 2008](#bib.bib32).

At the end of the day, the Bayes and frequentist answer different questions, and the common confusion arises when treating them as answering the same. As long as we are aware they are addressing different needs, we believe both solutions are informative and useful, and indeed that is our recommendation in this paper.

## Appendix B Derivation of IS weights for Hierarchical Models

We provide a quick derivation for the importance weights in Section [3.1](#S3.SS1 "3.1 Group Conformal Prediction ‣ 3 Partial Exchangeability and Hierarchical Models ‣ Conformal Bayesian Computation") to estimate

|  |  |  |
| --- | --- | --- |
|  | p⁡(Yi,j∣Xi,j,Z¯y)=∫fθj,τ​(Yi,j∣Xi,j)​π​(θj,τ∣Z¯y)​d​θj​𝑑τp(Y\_{i,j}\\mid X\_{i,j},\\bar{Z}\_{y})=\\int f\_{\\theta\_{j},\\tau}(Y\_{i,j}\\mid X\_{i,j})\\,\\pi(\\theta\_{j},\\tau\\mid\\bar{Z}\_{y})\\,d\\theta\_{j}\\,d\\tau |  |

from posterior samples \[θ1:J(1:T),τ(1:T),ϕ(1:T)\]∼π(θ1:J,τ,ϕ∣Z)\[\\theta\_{1:J}^{(1:T)},\\tau^{(1:T)},\\phi^{(1:T)}\]\\sim\\pi(\\theta\_{1:J},\\tau,\\phi\\mid Z). We can write the above as

|  | p⁡(Yi,j∣Xi,j,Z¯y)=\\displaystyle p(Y\_{i,j}\\mid X\_{i,j},\\bar{Z}\_{y})= | ∫fθj,τ(Yi,j∣Xi,j)π(θ1:J,ϕ,τ∣Z¯y)dθ1:Jdϕdτ\\displaystyle\\int f\_{\\theta\_{j},\\tau}(Y\_{i,j}\\mid X\_{i,j})\\,\\pi(\\theta\_{1:J},\\phi,\\tau\\mid\\bar{Z}\_{y})\\,d\\theta\_{1:J}\\,d\\phi\\,d\\tau |  |
| --- | --- | --- | --- |
|  |  | =∫fθj,τ(Yi,j∣Xi,j)π(θ1:J,ϕ,τ∣Z¯y)π(θ1:J,ϕ,τ∣Z)π(θ1:J,ϕ,τ∣Z)dθ1:Jdϕdτ\\displaystyle=\\int f\_{\\theta\_{j},\\tau}(Y\_{i,j}\\mid X\_{i,j})\\,\\frac{\\pi(\\theta\_{1:J},\\phi,\\tau\\mid\\bar{Z}\_{y})}{\\pi(\\theta\_{1:J},\\phi,\\tau\\mid Z)}\\pi(\\theta\_{1:J},\\phi,\\tau\\mid Z)\\,d\\theta\_{1:J}d\\phi\\,d\\tau |  |

where

|  |  |  |
| --- | --- | --- |
|  | π(θ1:J,ϕ,τ∣Z¯y)π(θ1:J,ϕ,τ∣Z)∝fθj,τ(y∣Xnj+1,j).\\frac{\\pi(\\theta\_{1:J},\\phi,\\tau\\mid\\bar{Z}\_{y})}{\\pi(\\theta\_{1:J},\\phi,\\tau\\mid Z)}\\propto f\_{\\theta\_{j},\\tau}(y\\mid X\_{n\_{j}+1,j}). |  |

As the importance weight only depends on θj,τ\\theta\_{j},\\tau, we only require the marginal posterior samples \[θj(1:T),τ(1:T)\]\[\\theta\_{j}^{(1:T)},\\tau^{(1:T)}\] and we have that

|  | p^​(Yi,j∣Xi,j,Z¯y)=∑t=1Tw(t)​fθj(t),τ(t)​(Yi,j∣Xi,j)\\displaystyle\\hat{p}(Y\_{i,j}\\mid X\_{i,j},\\bar{Z}\_{y})=\\sum\_{t=1}^{T}w^{(t)}f\_{\\theta\_{j}^{(t)},\\tau^{(t)}}(Y\_{i,j}\\mid X\_{i,j}) |  |
| --- | --- | --- |
|  | w(t)=fθj(t),τ(t)(y∣Xnj+1,j),w~(t)=w(t)∑t′=1Tw(t′)⋅\\displaystyle w^{(t)}=f\_{\\theta\_{j}^{(t)},\\tau^{(t)}}(y\\mid X\_{n\_{j}+1,j}),\\quad\\widetilde{w}^{(t)}=\\frac{w^{(t)}}{\\sum\_{t^{\\prime}=1}^{T}w^{(t^{\\prime})}}\\cdot |  |

## Appendix C Datasets, Licenses and Societal Impact

We demonstrate our examples on 5 datasets, namely the the diabetes dataset ([Efron et al., 2004](#bib.bib13)), the Boston housing dataset ([Harrison Jr and Rubinfeld, 1978](#bib.bib20)), the Wisconsin Breast cancer dataset ([Wolberg and Mangasarian, 1990](#bib.bib39)), the Parkinson’s dataset ([Little et al., 2008](#bib.bib24)) and the Radon dataset ([Gelman and Hill, 2006](#bib.bib18)). The first 3 datasets are available in sklearn, the Parkinson’s dataset can be found on the UCI machine learning repository22 2 [https://archive.ics.uci.edu/ml/datasets/Parkinson’s](https://archive.ics.uci.edu/ml/datasets/Parkinson's) ([Dua and Graff, 2017](#bib.bib11)), and the Radon dataset is available on Andrew Gelman’s website33 3 [http://www.stat.columbia.edu/~gelman/arm/examples/radon/](http://www.stat.columbia.edu/~gelman/arm/examples/radon/). Details on data acquisition is provided in the relevant references. We verified that the datasets do not contain personally identifiable information or offensive content by manual checking. The package sklearn is distributed under the 3-Clause BSD license. JAX and PyMC3 are both distributed under the Apache License, V2.

The conformal method relies on the weak assumption of exchangeability. In terms of negative societal impacts, it may be tempting to apply the method blindly to real world problems without challenging this assumption as it seems quite weak. Applications where calibration is very important but data is not exchangeable would then be at risk.

## Appendix D Additional Experiments

### D.1 Experimental Details

For all experiments, we repeat train-test splits or simulations 50 times, where the 70-30 train-test splits are random. For each repeat, we compute the average coverage and lengths for the test set. Means and standard errors are then computed from the 50 test set average coverages/lengths. For all MCMC examples, we generate T=8000T=8000 samples, with 40004000 tune steps for sparse regression/classification and 80008000 for the hierarchical example.

### D.2 Sparse Regression

#### D.2.1 Diabetes

We repeat analysis on the diabetes dataset, but this time with priors

|  | fθ​(y∣x)\\displaystyle f\_{\\theta}(y\\mid x) | =𝒩⁡(y∣θT​x+θ0,τ2)\\displaystyle=\\mathcal{N}(y\\mid\\theta^{\\mathrm{\\scriptscriptstyle T}}x+\\theta\_{0},\\tau^{2}) |  | (11) |
| --- | --- | --- | --- | --- |
|  | π⁡(θj)=Normal​(0,d),π⁡(θ0)\\displaystyle\\pi(\\theta\_{j})=\\text{Normal}(0,d),\\quad\\pi(\\theta\_{0}) | =Normal​(0,d),π⁡(τ)=𝒩+​(0,1)\\displaystyle=\\text{Normal}(0,d),\\quad\\pi(\\tau)=\\mathcal{N}^{+}(0,1) |  |  |

where we have different values d=5,0.001d=5,0.001 which corresponds to weak and strong regularization towards 0. For the baselines, we instead use ridge regression, with and without cross-validation for split/full as before. We emphasize that it is not exactly a fair comparison for the d=0.001d=0.001 case, as the baselines are tuning the parameter λ\\lambda, whereas CB is subject to the misspecified prior. We still include them as baselines however, but highlight that they are not affected by the value of dd.

MCMC required 30.8s and 13.2s for d=5d=5 and d=0.001d=0.001 respectively. The effect on coverage of setting d=0.001d=0.001 is not as detrimental as before as seen in Table [4](#A4.T4 "Table 4 ‣ D.2.1 Diabetes ‣ D.2 Sparse Regression ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation"), as the posterior on τ\\tau compensates by increasing in value; the posterior mean is τ¯=1\\bar{\\tau}=1 for d=0.001d=0.001 versus τ¯=0.71\\bar{\\tau}=0.71 for d=5d=5.

Table 4: Diabetes; Coverage values not within 3 standard errors (in brackets) of the target coverage (1−α)=0.8(1-\\alpha)=0.8 are in red.

|  |  | Bayes | CB | Split | Full (λ=0.004\\lambda=0.004) |
| --- | --- | --- | --- | --- | --- |
| Coverage | d=5d=5 | 0.805 (0.005) | 0.809 (0.005) | 0.816 (0.006) | 0.809 (0.005) |
|  | d=0.001d=0.001 | 0.779 (0.006) | 0.809 (0.006) | / | / |
| Length | d=5d=5 | 1.85 (0.01) | 1.86 (0.01) | 1.94 (0.02) | 1.86 (0.01) |
|  | d=0.001d=0.001 | 2.56 (0.01) | 2.60 (0.01) | / | / |
| Run-time | d=1d=1 | 0.417 (0.002) | 0.677 (0.003) | 0.024 (0.000) | 8.409 (0.007) |
| (secs) | d=0.02d=0.02 | 0.540 (0.116) | 0.692 (0.008) | / | / |

#### D.2.2 Boston Housing

The Boston housing dataset ([Harrison Jr and Rubinfeld, 1978](#bib.bib20)) is of size n=506n=506, consisting of d=13d=13 predictors relating to housing such as demographic and air quality, with the response as the median value of owner-occupied homes. We use the same Bayesian model as in ([9](#S4.E9 "In 4.1 Sparse Regression ‣ 4 Experiments ‣ Conformal Bayesian Computation")), again considering c=1,0.02c=1,0.02. For c=1c=1, the model is already misspecified for the Boston housing dataset as the errors are non-normal and have heavy tails ([Jansen, 2013](#bib.bib21)). All experimental settings are the same as in Section [4.1](#S4.SS1 "4.1 Sparse Regression ‣ 4 Experiments ‣ Conformal Bayesian Computation").   

MCMC required an average of 22.822.8s and 24.424.4s for c=1,0.02c=1,0.02 to produce T=8000T=8000 posterior samples. Again, in Table [5](#A4.T5 "Table 5 ‣ D.2.2 Boston Housing ‣ D.2 Sparse Regression ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation") we see similar behaviour to the diabetes dataset case, but we note that even under c=1c=1, the Bayesian model over-covers. This is likely due to the presence of heavy tails in the residuals, leading to more conservative Bayesian predictive intervals. Here for c=1c=1, CB attains very close to nominal coverage and has a noticeably smaller average length. For c=0.02c=0.02, CB is not affected much but the Bayes interval under-covers.

Table 5: Boston; Coverage values not within 3 standard errors (in brackets) of the target coverage (1−α)=0.8(1-\\alpha)=0.8 are in red.

|  |  | Bayes | CB | Split | Full (λ=0.004\\lambda=0.004) |
| --- | --- | --- | --- | --- | --- |
| Coverage | c=1c=1 | 0.860 (0.004) | 0.800 (0.005) | 0.805 (0.006) | 0.799 (0.005) |
|  | c=0.02c=0.02 | 0.728 (0.005) | 0.799 (0.005) | / | / |
| Length | c=1c=1 | 1.35 (0.01) | 1.12 (0.01) | 1.22 (0.02) | 1.12 (0.01) |
|  | c=0.02c=0.02 | 0.96 (0.00) | 1.13 (0.01) | / | / |
| Run-time | c=1c=1 | 0.414 (0.003) | 0.746 (0.011) | 0.061 (0.000) | 12.448 (0.042) |
| (secs) | c=0.02c=0.02 | 0.406 (0.003) | 0.744 (0.003) | / | / |

### D.3 Sparse Classification

#### D.3.1 Parkinson’s Disease

We provide an additional demonstration on the Parkinson’s dataset ([Little et al., 2008](#bib.bib24)), which consists of n=195n=195 voice recordings (after removing missing data) of patients with or without Parkinson’s disease encoded in the binary response. The covariates consist of d=22d=22 different voice recording properties.   

The experimental setup is identical to Section [4.2](#S4.SS2 "4.2 Sparse Classification ‣ 4 Experiments ‣ Conformal Bayesian Computation"), and MCMC required 29.2s to produce T=8000T=8000 samples. Again, in Table [6](#A4.T6 "Table 6 ‣ D.3.1 Parkinson’s Disease ‣ D.3 Sparse Classification ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation") we see that Bayes over-covers even for reasonable priors, and CB produces tighter intervals that are closer to nominal coverage.

Table 6: Parkinson’s; Coverage values not within 3 standard errors (in brackets) of the target coverage (1−α)=0.8(1-\\alpha)=0.8 are in red. “Size” denotes the average number of elements in the conformal prediction set, averaged over the test points and repetitions.

|  | Bayes | CB | Split | Full |
| --- | --- | --- | --- | --- |
| Coverage | 0.955 (0.004) | 0.815 (0.008) | 0.842 (0.010) | 0.816 (0.008) |
| Size | 1.31 (0.01) | 0.93 (0.01) | 1.05 (0.02) | 0.95 (0.01) |
| Times | 0.203 (0.003) | 0.379 (0.008) | 0.478 (0.057) | 0.168 (0.003) |

#### D.3.2 Uninformative Predictions

As the Bayesian model returns p:=p(y=1∣x,Z1:n)p:=p(y=1\\mid x,Z\_{1:n}), we compute (1−α)(1-\\alpha) predictive sets by returning the smallest set of {0},{1},{0,1}\\{0\\},\\{1\\},\\{0,1\\} such that it contains at least (1−α)(1-\\alpha) of the predictive probability mass. In other words, we return:

|  |  |  |  |
| --- | --- | --- | --- |
|  | {{0}if ​(1−p)≥(1−α){1}if ​p≥(1−α){0,1}if ​max⁡{(1−p),p}≤(1−α)\\begin{cases}\\{0\\}\\quad&\\text{if }(1-p)\\geq(1-\\alpha)\\\\ \\{1\\}\\quad&\\text{if }p\\geq(1-\\alpha)\\\\ \\{0,1\\}\\quad&\\text{if }\\max\\{(1-p),p\\}\\leq(1-\\alpha)\\\\ \\end{cases} |  | (12) |

As this process is quite conservative, it is unsuprising that Bayes overcovers. The set {0,1}\\{0,1\\} is clearly uninformative, as it is always correct. On the other hand, the conformal sets can take on the empty set {}\\{\\} as well, which we know to be incorrect. As discussed in [Melluish et al., 2001](#bib.bib27), empty set predictions correspond to being unable to make a prediction at the desired confidence level. [Shafer and Vovk, 2008](#bib.bib32) discusses the notion of confidence and credibility, which correspond to the greatest (1−α)(1-\\alpha) such that the conformal set is of size 1 and the greatest α\\alpha such that the conformal set is empty respectively.

We can decompose the informative and uninformative predictive sets and look at the misclassification rate, which is the error percentage within single element predictions. In comparison, we can look at the percentage of uninformative predictives (either both elements or empty). This is shown in Tables [7](#A4.T7 "Table 7 ‣ D.3.2 Uninformative Predictions ‣ D.3 Sparse Classification ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation"), [8](#A4.T8 "Table 8 ‣ D.3.2 Uninformative Predictions ‣ D.3 Sparse Classification ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation"), where the target coverage is (1−α)=0.8(1-\\alpha)=0.8 as before. For the breast cancer dataset, CB has a very small misclassification rate, but almost 19% of all prediction sets are empty and 0% are both, so the coverage is attained by making either single correct predictions or empty ones. Bayes on the other consists of more misclassifications but fewer uninformative predictions, but the attained coverage is a much higher value of 0.99. For the Parkinson’s dataset, CB makes very few uninformative predictions, but has a relatively high misclassification rate. Bayes on the other hand is very conservative, with 31% uninformative predictions, hence the high average length and over-coverage. It is interesting to note the two sorts of behaviours attained by CB, which likely depends on the Bayesian model that was used to construct the CB intervals.

In Figure [2](#A4.F2 "Figure 2 ‣ D.3.2 Uninformative Predictions ‣ D.3 Sparse Classification ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation"), we see the distributions of pi:=p⁡(yi∣xi,Z)p\_{i}:=p(y\_{i}\\mid x\_{i},Z) of the Bayesian model with the corresponding CB interval length. We see that for CB intervals of length 1, the values of pip\_{i} tend to be heavily skewed towards 00 or 11, which corresponds to the Bayesian model being strongly predictive. For empty CB intervals, in both cases the probability mass is distributed away from 00 and 11; for the breast cancer dataset it is evenly distributed on (0,1)(0,1) whereas for Parkinson’s it is concentrated around 0.50.5. When given a CB interval of length 0, it may be more informative to actually return the value pip\_{i}, which is the corresponding Bayesian prediction.

Table 7: Misclassification rates

| Dataset | Bayes | CB |
| --- | --- | --- |
| Breast Cancer | 0.011 (0.001) | 0.002 (0.000) |
| Parkinson’s | 0.064 (0.006) | 0.124 (0.004) |

Table 8: Uninformative Rates

|  | Both | Empty |  |  |
| --- | --- | --- | --- | --- |
| Dataset | Bayes | CB | Bayes | CB |
| Breast Cancer | 0.059 (0.002) | 0.000 (0.000) | 0 | 0.186 (0.005) |
| Parkinson’s | 0.312 (0.009) | 0.003 (0.001) | 0 | 0.070 (0.006) |

Figure 2: Distribution of p⁡(yi∣xi,Z)p(y\_{i}\\mid x\_{i},Z) on test data for CB intervals of length 11 or 00 for breast cancer (left) and Parkinson’s (right).

### D.4 Hierarchical

#### D.4.1 Radon

Using the same model as in ([10](#S4.E10 "In 4.3 Hierarchical Model ‣ 4 Experiments ‣ Conformal Bayesian Computation")), we analyze the radon dataset44 4 We base this example on the PyMC3 notebook here: [https://docs.pymc.io/notebooks/multilevel\_modeling.html/](https://docs.pymc.io/notebooks/multilevel_modeling.html), introduced in [Gelman and Hill, 2006](#bib.bib18). The dataset consists of 919 home radon levels in Minnesota, where the covariate is the location of measurement, with x=0x=0 corresponding to basement and x=1x=1 to the first floor. The groups are the 85 counties in which the homes are located, and vary significantly in group size. Around half of the counties contain nj≤4n\_{j}\\leq 4 measurements, with the smallest county containing one value and the largest containing 119.

As many of the group sizes are quite small, we do not repeat train-test splits and evaluate coverage. Instead, we compare the CB and Bayes intervals on the entire dataset for different floor values xx and counties, and discuss the effects of njn\_{j} on the choice of αj\\alpha\_{j}. As each x∈{0,1}x\\in\\{0,1\\}, we specify xtestx\_{\\text{test}} as all possible group indicators and predictors, resulting in 85×2=17085\\times 2=170 test values. For the predictive intervals, we use a grid of size 100 between \[-6,6\]. MCMC for the radon example required around 156s, and computing the 170 predictive intervals took 0.65s and 2.69s for Bayes and CB respectively, where we have excluded the first run compilation time for JAX.

As we need αj≥1/(nj+1)\\alpha\_{j}\\geq{1}/({n\_{j}+1}) to get intervals that are not the entire real line, we set αj=1.1/(nj+1)\\alpha\_{j}={1.1}/({n\_{j}+1}) (for numerical reasons) and compare the Bayes and CB intervals. The average CB length is 2.66 compared to 2.17 for the Bayes intervals, noting that we are averaging over all possibilities instead of the distribution of xtestx\_{\\text{test}}. In Figure [3](#A4.F3 "Figure 3 ‣ D.4.1 Radon ‣ D.4 Hierarchical ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation"), we plot πj​(y)\\pi\_{j}(y) for the two value of x∈{0,1}x\\in\\{0,1\\} for two groups. For the group size nj=4n\_{j}=4, we see that πj​(y)≥0.2\\pi\_{j}(y)\\geq 0.2, so any α<0.2\\alpha<0.2 would return us the real line as the confidence set. For nj=52n\_{j}=52, the ranks are much smoother, giving us more resolution in the confidence sets with respect to α\\alpha. In Figure [4](#A4.F4 "Figure 4 ‣ D.4.1 Radon ‣ D.4 Hierarchical ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation"), we show the rank plots for nj=1n\_{j}=1, which only contain the ranks {0.5,1}\\{0.5,1\\}. Interestingly, for county 41, x=1x=1 returns the empty set for α≥0.5\\alpha\\geq 0.5 and the real line for α<0.5\\alpha<0.5, which is a consequence of the small group size. CB is able to return non empty sets for county 49 with α≥0.5\\alpha\\geq 0.5. All CB sets appear to be connected.

Figure 3: Plot of rank πj​(y)\\pi\_{j}(y) for x∈{0,1}x\\in\\{0,1\\} with nj=4n\_{j}=4 (left) and nj=52n\_{j}=52 (right).

Figure 4: Plot of rank πj​(y)\\pi\_{j}(y) for x∈{0,1}x\\in\\{0,1\\} with nj=1n\_{j}=1 for two groups.

As a reference, fitting a linear mixed-effects model in statsmodels ([Seabold and Perktold, 2010](#bib.bib31)) to the whole dataset takes around 600ms, so the full conformal method, which would require refitting for each of the 100 grid value and 170 test values, would require 170 minutes in total. On the other hand, CB requires much less time and has similar group structure.

### D.5 MCMC Times

In Table [9](#A4.T9 "Table 9 ‣ D.5 MCMC Times ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation"), we report the average times (and standard errors) for running MCMC on the NC6 virtual machine. We point out that the times for the hierarchical methods are longer as we needed to increase the tuning steps and acceptance probability to prevent divergences in the chains.

Table 9: Run-time in seconds for MCMC

| Dataset | MCMC |
| --- | --- |
| Diabetes (c=1)(c=1) | 21.868 (0.135) |
| Diabetes (c=0.02)(c=0.02) | 26.790 (0.365) |
| Diabetes (d=5)(d=5) | 30.825 (0.214) |
| Diabetes (d=0.001)(d=0.001) | 13.166 (0.072) |
| Boston (c=1)(c=1) | 22.827 (0.036) |
| Boston (c=0.02)(c=0.02) | 24.362 (0.429) |
| Breast Cancer | 45.418 (0.804) |
| Parkinson’s | 29.239 (0.302) |
| Scenario 1 | 90.109 (1.605) |
| Scenario 2 | 78.403 (1.544) |
| Radon | 156.087 (0.000) |

### D.6 Grid effects

To quantify the grid effects, we also compute the coverage by directly evaluating π⁡(Yn+1)\\pi(Y\_{n+1}) for each test point and checking if it satisfies condition ([3](#S2.E3 "In 2.1 Full Conformal Prediction ‣ 2 Conformal Bayes ‣ Conformal Bayesian Computation")). Of course in practice this is not possible as we do not observe Yn+1Y\_{n+1}.

For the grid conformal method, we compute the y∈𝒴gridy\\in\\mathcal{Y}\_{\\text{grid}} that is nearest to Yn+1Y\_{n+1}, and report 00 or 11 if this grid value is in the conformal prediction set. Note that this implementation of the grid method can both under and over cover. Denote δ\\delta as the resolution of the grid, and the smallest grid value in the conformal prediction set as aa. If a−δ<Yn+1<a−δ/2a-\\delta<Y\_{n+1}<a-\\delta/2, we may incorrectly reject Yn+1Y\_{n+1} if it is truly in the set and a−δa-\\delta is not. Similarly, if a−δ/2<Yn+1<aa-\\delta/2<Y\_{n+1}<a we can incorrectly accept if Yn+1Y\_{n+1} is not actually in the set but aa is. Note that the estimated average length is also affected by this.

We compare the grid and exact method in Tables [10](#A4.T10 "Table 10 ‣ D.6 Grid effects ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation"), [11](#A4.T11 "Table 11 ‣ D.6 Grid effects ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation"), [12](#A4.T12 "Table 12 ‣ D.6 Grid effects ‣ Appendix D Additional Experiments ‣ Conformal Bayesian Computation"). The largest discrepancy in average coverage is only 0.0080.008, which is quite negligible. However, we expect this discrepancy to increase as |𝒴grid||\\mathcal{Y}\_{\\text{grid}}| decreases.

Table 10: Diabetes; Grid versus exact coverage, with target (1−α)=0.8(1-\\alpha)=0.8

|  |  | CB Grid | CB Exact |
| --- | --- | --- | --- |
| Coverage | c=1c=1 | 0.808 (0.006) | 0.810 (0.005) |
|  | c=0.02c=0.02 | 0.809 (0.006) | 0.810 (0.006) |

Table 11: Boston; Grid versus exact coverage, with target (1−α)=0.8(1-\\alpha)=0.8

|  |  | CB Grid | CB Exact |
| --- | --- | --- | --- |
| Coverage | c=1c=1 | 0.800 (0.005) | 0.800 (0.005) |
|  | c=0.02c=0.02 | 0.799 (0.005) | 0.799 (0.005) |

Table 12: Simulated grouped dataset; Grid versus exact coverage, with target (1−α)=0.8(1-\\alpha)=0.8

|  |  | Scenario 1 | Scenario 2 |  |  |
| --- | --- | --- | --- | --- | --- |
|  | Group | CB Grid | CB Exact | CB Grid | CB Exact |
| Coverage | 1 | 0.794 (0.022) | 0.786 (0.023) | 0.786 (0.025) | 0.790 (0.024) |
|  | 2 | 0.812 (0.024) | 0.816 (0.023) | 0.812 (0.024) | 0.818 (0.023) |
|  | 3 | 0.824 (0.022) | 0.820 (0.022) | 0.824 (0.020) | 0.824 (0.020) |
|  | 4 | 0.798 (0.022) | 0.796 (0.021) | 0.836 (0.021) | 0.838 (0.022) |
|  | 5 | 0.810 (0.020) | 0.812 (0.019) | 0.796 (0.022) | 0.792 (0.022) |