###### R function for determining the value of a in the 
###### Gamma(a, 1) prior for the range parameter phi in the correlation function


### M: the number of simulations in each homogeneous scenario
### N: the number of patients enrolled in each indication group
### alpha: the significance level of the confidence interval used in the first step of the algorithm
### I: the total number of indications in the trial, 
### varrho.lb: the lower bound of the correlation threshold varrho
### varrho.ub: the upper bound of the correlation threshold varrho

###### Recommended setting: M = 5000, alpha = 0.05, varrho.lb = 0.3, varrho.ub = 0.5.

a_select = function(M, N, q0, q1, alpha, I, varrho.lb, varrho.ub)
{
  ######## Step 1.1: generate sample distance under the two homogeneous scenarios
  dsample = NULL
  for (i in 1:I)
  {
    ##### Homogeneous scenario 1: p1 = p2 = q0[i]
    rr1 = rbinom(M,size=N[i],prob=q0[i])
    rr2 = rbinom(M,size=N[i],prob=q0[i])
    a1=1+rr1
    a2=1+rr2
    b1=1+N[i]-rr1
    b2=1+N[i]-rr2
    DB= -log(beta((a1+a2)/2,(b1+b2)/2)/sqrt(beta(a1,b1)*beta(a2,b2)))
    dsample = c(dsample, DB)
    ##### Homogeneous scenario 2: p1 = p2 = q1[i]
    rr1 = rbinom(M,size=N[i],prob=q1[i])
    rr2 = rbinom(M,size=N[i],prob=q1[i])
    a1=1+rr1
    a2=1+rr2
    b1=1+N[i]-rr1
    b2=1+N[i]-rr2
    DB= -log(beta((a1+a2)/2,(b1+b2)/2)/sqrt(beta(a1,b1)*beta(a2,b2)))
    dsample = c(dsample, DB)
  }
  ######## Step 1.2: calculate the d_t: the 1-alpha quantile of dsample
  d_t = quantile(x = dsample, probs = 1-alpha)
  ######## Step 2: Determine the lower and upper bound for a
  a.lb = -log(varrho.ub)/d_t
  a.ub = -log(varrho.lb)/d_t
  ######## Step 3: randomly select a value for a from Unif(a.lb,a.ub)
  a = runif(1,a.lb,a.ub)
  return(a)
}



########### example:
a = a_select(M = 5000, N = c(20,22,24,22,22), q0 = c(0.2, 0.15, 0.25, 0.2, 0.2), q1 = c(0.4, 0.35, 0.35, 0.45, 0.4), 
             alpha = 0.1, I = 5, varrho.lb = 0.3, varrho.ub = 0.5)
a

