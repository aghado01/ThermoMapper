#LSE: Estimate Change-Points in the Linear Spline Model
LSE<-function(X,Y,tau0,tol=1e-6,max.iter=1000)
{
  n_knots = length(tau0)
  erro = FALSE;count = 1;Test.0=TRUE
  while (Test.0)
  {
    ### update_step1: fix the change-points and estimate other coefficients ###
    X_tau = sapply(tau0,function(tau) {return(as.numeric(X>tau)*(X-tau))})
    lin.inital <- lm(Y~ cbind(X,X_tau)) #fit a linear regression
    theta0 = c(as.numeric(coef(lin.inital)),tau0)
    y_pred0 = predict(lin.inital)

    ### update_step2: update change-points via a modified NR procedure ###
    I_tau = sapply(tau0,function(tau) {return(as.numeric(X>tau))})
    sig_error = Y-y_pred0
    ### Calculate U and J ###
    Un = sapply(1:n_knots,function(i) {return(theta0[2+i]*mean(sig_error*I_tau[,i]))})
    Un = cbind(Un)

    Jn = matrix(nrow=n_knots,ncol=n_knots)
    for(i in c(1:n_knots)) {
      Jn[i,i] = -theta0[i+2]^2*mean(I_tau[,i])
      for(j in c(1:n_knots)[-i]) {
        Jn[i,j] = -theta0[i+2]*theta0[j+2]*mean(I_tau[,i]*I_tau[,j])
      }
    }
    ### Calculate the update_step ###
    tryCatch( {update_step = solve(Jn)%*%Un},error=function(cond) erro = TRUE)
    ### Update Tau ###
    if (erro == FALSE & !any(is.na(update_step)))
    {tau0 = tau0 + update_step}
    count=count+1

    Test.0=any(abs(update_step) > tol) & count<=max.iter & erro==FALSE #count: number of iterations
    if (is.na(Test.0)){Test.0=FALSE}#whether the algorithm is running properly
  }

  if(FALSE) {
  tau.out = rep(NA,n_knots)
  Test = (erro == FALSE & max(tau0,na.rm=T)<max(X) & min(tau0,na.rm=T)>min(X))
  if (Test & (count<=max.iter| (count==(max.iter+1) & all(abs(update_step) <= sqrt(tol)))))
  {tau.out = c(tau0)}
  }

  return(tau0)
}
