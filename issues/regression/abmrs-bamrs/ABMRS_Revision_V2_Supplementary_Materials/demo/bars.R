#BIC approximation of marginal likelihood
library(splines)

# sample (k,xi) through reversible jump mcmc
Birth <- function(k,c=0.4,lam=3){
  # compute b_k
  # k: the number of splines' knots
  return(c*min(1,dpois(k+1,lam)/dpois(k,lam)))
}
Death <- function(k,c=0.4,lam=3){
  # compute d_k
  return(c*min(1,dpois(k-1,lam)/dpois(k,lam)))
}
Relocate <- function(k,c=0.4,lam=3){
  # compute r_k
  return(1-Birth(k)-Death(k))
}

# predict y values in x through BARS estimation
bars_fit <- function(x,conf){
  # conf: return value of BARS function
  # value: \hat{f}(x)
  n = length(x)
  t = (x - conf$range[1])/(conf$range[2] - conf$range[1])
  B = bs(t,knots=conf$knots,Boundary.knots=c(0,1),intercept=conf$intercept)
  y = B%*%conf$beta
  return(y)
}

#Bayesian adaptive regression splines
BARS <- function(x,y,c=0.3,mu=50,burns_in=1000,effective_length=1000,lam=1,intercept=TRUE){
  # return knots' configuration including the number and location of knots
  # and corresponding beta_mle
  # x: d*1 vector, y: d*1 vector, d=length(x)=length(y)
  # mu: hyper-parameter for beta distribution
  # lam: location parameter for the prior Poisson distribution
  # value: conf, a list of n objects, each object is a list(num,knot,beta)

  init.time = Sys.time()

  x.min = min(x)
  x.max = max(x)
  t = (x - x.min)/(x.max - x.min)

  y = matrix(y, ncol=1)
  nsamp = length(y)

  #initialize
  k = 1;xi = runif(1)
  # cat("xi.init: ",xi,"\n")

  n = 0
  # B_hat: length(x)*(length(knots)+3) basis matrix on x with knots
  B_hat = bs(t, knots = xi, Boundary.knots = c(0,1), intercept = intercept)
  BTB = t(B_hat)%*%B_hat
  beta_hat = chol2inv(chol(BTB+0.01*sum(diag(BTB))/nrow(BTB)
                                *diag(nrow(BTB))))%*%t(B_hat)%*%y
  rss_sd = norm(y-B_hat%*%beta_hat,type='f')

  while(n<burns_in+effective_length){
    # construct mcmc to sample (k,xi) from the posterior distribution
    # p(k,xi|y) based on birth and death chain

    sche = runif(1,0,1)
    birth = Birth(k,c,lam)
    death = Death(k,c,lam)
    i = sample(1:k,1)

    if(k==1) {
      death = 0
    }

    # birth scheme
    if(sche<birth){
      xi_can = (rbeta(1,xi[i]*mu,(1-xi[i])*mu)-0.5)*(1-1e-5)+0.5
      xi_next = sort(c(xi,xi_can))
      k_next = k+1
      # compute proposal density ratio
      mix = 0
      for(j in 1:k){
        mix = mix + dbeta(xi_can, xi[j]*mu, (1-xi[j])*mu)
      }
      q_ratio = (Death(k_next,c,lam)/(k+1))/(Birth(k,c,lam)/k*mix)
    }

    # death scheme
    else if(sche<birth+death){
      xi_del = xi[i]
      xi_next = xi[-i]
      k_next = k - 1
      # compute proposal density ratio
      mix = 0
      for(j in 1:k_next){
        mix = mix + dbeta(xi_del,xi_next[j]*mu,(1-xi_next[j])*mu)
      }
      q_ratio = (Birth(k_next,c,lam)/k_next*mix)/(Death(k,c,lam)/k)
    }

    # relocation scheme
    else {
      xi_chan = xi[i]
      xi_can = (rbeta(1,xi_chan*mu,(1-xi_chan)*mu)-0.5)*(1-1e-5)+0.5
      xi_next = sort(c(xi[-i],xi_can))
      k_next = k
      # compute proposal density ratio
      q_ratio = dbeta(xi_chan,xi_can*mu,(1-xi_can)*mu)/
        dbeta(xi_can,xi_chan*mu,(1-xi_chan)*mu)
    }

    #compute accept probability to decide whether to move
    B_hat_next = bs(t, knots = xi_next, Boundary.knots = c(0,1), intercept = intercept)
    BTB_next = t(B_hat_next)%*%B_hat_next
    beta_hat_next = chol2inv(chol(BTB_next+0.01*sum(diag(BTB_next))/(nrow(BTB_next))
                                  *diag(nrow(BTB_next))))%*%t(B_hat_next)%*%y

    #margin_like_ratio = n^((k-k')/2)*(RSS_k/RSS_k')^(n/2)
    rss_sd_next = norm(y-B_hat_next%*%beta_hat_next,type='f')
    like_ratio = nsamp^((k-k_next)/2)*(rss_sd/(rss_sd_next+1e-7))^nsamp
    #acc_ratio = min(1,A) with A = q_ratio*margin_like_ratio*lambda^(k'-k)
    acc_ratio = q_ratio*like_ratio*lam^(k_next-k)
    acc = runif(1,0,1)

    #decide whether to move or not
    if(acc<acc_ratio){
      n = n + 1
      k = k_next;xi = xi_next
      B_hat = B_hat_next;beta_hat = beta_hat_next;rss_sd = rss_sd_next
    }
  }

  result = list(knots=xi,beta=beta_hat,range=c(x.min, x.max),intercept=intercept)
  return(result)
}
