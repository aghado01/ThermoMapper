Junhui He, Department of Mathematical Sciences, Tsinghua University

<!-- README.md is generated from README.Rmd. Please edit that file -->

# ABMRS

<!-- badges: start -->

<!-- badges: end -->

The goal of ABMRS is to implement the adaptive Bayesian multivariate
regression spline algorithm of the article ***Adaptive Bayesian
Multivariate Spline Knot Inference with Prior Specifications on Model
Complexity***. In spline regression, the number and location of knots
influence the performance and interpretability significantly. We propose
a fully Bayesian approach for knot inference in multivariate spline
regression. The knot inference is of interest in many problems, such as
change point detection. We can estimate the knot number and location
simultaneously and accurately by the proposed method in univariate or
multivariate splines.

We specify a prior on the knot number to take into account the
complexity of the model space and derive an analytic formula in the
normal model. In the non-normal cases, we utilize the extended Bayesian
information criterion (EBIC) to approximate the posterior density. The
samples are simulated in the space with differing dimensions via
reversible jump Markov chain Monte Carlo (RJMCMC). Experiments
demonstrate the splendid capability of the algorithm, especially in
function fitting with jumping discontinuity.

More examples refer to
[exampleABMRS](https://github.com/junhuihe2000/exampleABMRS).

## Installation

You can install the development version of ABMRS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("junhuihe2000/ABMRS")
```

## Example

This is a basic example which shows you how to make knot inference:

1.  Generate data

``` r
library(ABMRS)
library(splines)
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.3.3
library(ggpubr)

f = function(x) {
  knots = c(0.2,0.2,0.5,0.7)
  beta = c(0,-1,1,0,1,0)
  B = bs(x,knots=knots,degree=1,intercept=TRUE,Boundary.knots=c(0,1))
  y = B%*%beta
  return(y)
}

set.seed(1234)
m = 200
x_train = runif(m,0,1)

y_train = f(x_train)
y_obs = y_train + rnorm(m,0,0.4)
p = ggplot() + geom_point(aes(x_train,y_obs),size=0.6,col="red") + geom_line(aes(x_train,y_train)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank())

annotate_figure(p, top="Linear splines with four knots")
```

<img src="man/figures/README-data-1.png" width="50%" style="display: block; margin: auto;" />

2.  Knot inference

``` r
set.seed(1234)

m = 500; sd = 0.4
x = runif(m)
y = f(x) + rnorm(m,0,sd)
# ABMRS
ebars = ebars(x,y,gamma=1,c=0.3,n=1000,degree=1,intercept=TRUE)
ebars$rjmcmc(burns=5000,steps=5000)
knots = ebars$knots()
nums = sapply(knots, function(xi) {return(length(xi))})
points = unlist(knots)
p1 = ggplot(mapping=aes(x=nums)) + geom_bar(aes(y=after_stat(prop))) + ylab("") + xlab("") + geom_vline(aes(xintercept=mean(nums)), color="red", linetype="dashed") + scale_x_continuous(breaks=seq(min(nums),max(nums)),labels=seq(min(nums),max(nums)))
p2 = ggplot(mapping=aes(x=points)) + geom_density(adjust=1) + ylab("") + xlab("")

annotate_figure(ggarrange(p1,p2,ncol=2,nrow=1,align="hv"),top="The posterior distribution of knots by ABMRS")
```

<img src="man/figures/README-inference-1.png" width="100%" />
