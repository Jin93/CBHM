###################### an example code for CBHM using KL distance. 
##### the only differences from using B distance:
# 1. the specification of D matrix
# 2. the jags code for CBHM
# 3. Q.cbhm.b -> Q.cbhm.kl
library(entropy)
library(rjags)
library(mvtnorm)
library(boot) # for inverse logit function
library(nonpar) # for Cochran's Q test
library(distrEx) #for TotalVarD
library(nonpar) # for function cochrans.q

###### General Simulation Settings:
K=6 # number of indications
alpha=0.1 # significance level for the test
num.sim=5000 # the number of simulations per simulation setting
Ni=24 # the maximum of total sample size for each indication group
Ni1=14 # stage-one sample size for each indication group
nik=matrix(NA,2,K) # each row records the number of patients in indication k at stage i
rik=matrix(NA,2,K) # each row records the number of responders in indication k at stage i
nik[1,]=rep(Ni1,K) # number of patients enrolled at stage 1

q0=0.2 # standard of care (null) response rate
q1=0.4 # target response rate
##################################################
############### Calibration Stage: ###############
##################################################
p0=rep(q0,K) # true rr: set the true rr for all indications to null rr (null scenario)
############## calibrate the tuning parameters for each method so that 
############## the type I error rate is well controlled under the null scenario

########## Calibration for CBHM with B distance measure & expnential correlation function:
Qf=0.05 # probability cut-off for interim analysis
epsilon = 3*(q1-q0)/K # the small value added to the number of responsders for the indication groups that have equal sample rr

posterior.cbhm.b=matrix(NA,num.sim,K)
for (sim in 1:num.sim)
{
  ##### Stage 1:
  rik[1,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[1,x],prob=p0[x])})
  rikstar=rik
  ######### H distance matrix:
  D=matrix(NA,K,K)
  for (j in 1:(K-1))
  {
    for (k in (j+1):K)
    {
      if ((rik[1,j]/nik[1,j])==(rik[1,k]/nik[1,k]))
      {
        rik[1,k]=rik[1,j] + epsilon
      }
    }
  }
  for (j in 1:K)
  {
    for (k in j:K)
    {
      if (j == k)
      {
        D[j,k] = 0
      }
      if (j != k)
      {
        cdf1 = pbeta(seq(0,0.8,1/100),1+rik[1,j],1+nik[1,j]-rik[1,j])
        freqs1 = sapply(1:(length(cdf1)-1),FUN=function(x){cdf1[x+1]-cdf1[x]})
        cdf2 = pbeta(seq(0,0.8,1/100),1+rik[1,k],1+nik[1,k]-rik[1,k])
        freqs2 = sapply(1:(length(cdf2)-1),FUN=function(x){cdf2[x+1]-cdf2[x]})
        D[j,k] = D[k,j] = (KL.plugin(freqs1, freqs2)+KL.plugin(freqs2, freqs1))/2
      }
    }
  }
  rik=rikstar
  zero=rep(0,K)
  ## Jags model for CBHM:
  jags.data <- list("n"=nik[1,], "Y"=rik[1,], "D"=D, "K"=K, "zero"=zero,'mu0'=log(((q0+q1)/2)/(1-(q0+q1)/2)))
  jags.fit <- jags.model(file = "cbhm_kl.txt",data = jags.data,
                         n.adapt=1000,n.chains=1,quiet=T)
  update(jags.fit, 4000,quiet=T)
  cbhm.out <- coda.samples(jags.fit,variable.names = c("theta0","phi","tausq","tausq2","e","p"),n.iter=10000,quiet=T)
  ## Interim analysis:
  posterior=numeric()
  for (k in 1:K)
  {
    post.sample=cbhm.out[[1]][,paste("p[",k,"]",sep="")]
    posterior[k]=sum(post.sample>(q0+q1)/2)/length(post.sample)
  }
  ## Futility stop:
  stage2.stop=which(posterior<Qf)
  stage2.cont=which(posterior>=Qf)
  nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
  posterior.cbhm.b[sim,stage2.stop]=0
  ## Stage 2:
  if (length(stage2.cont)>0)
  {
    rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
    ri=colSums(as.matrix(rik[,stage2.cont]))
    ristar=ri
    ni=colSums(as.matrix(nik[,stage2.cont]))
    K1=length(stage2.cont)
    R=D=matrix(NA,K1,K1)
    if (K1>1)
    {
      for (j in 1:(K1-1))
      {
        for (k in (j+1):K1)
        {
          if (ri[j]/ni[j]==ri[k]/ni[k])
          {
            ri[k] = ri[j] + epsilon
          }
        }
      }
    }
    
    for (j in 1:K1)
    {
      for (k in j:K1)
      {
        if (j == k)
        {
          D[j,k] = 0
        }
        if (j != k)
        {
          cdf1 = pbeta(seq(0,0.8,1/100),1+ri[j],1+ni[j]-ri[j])
          freqs1 = sapply(1:(length(cdf1)-1),FUN=function(x){cdf1[x+1]-cdf1[x]})
          cdf2 = pbeta(seq(0,0.8,1/100),1+ri[k],1+ni[k]-ri[k])
          freqs2 = sapply(1:(length(cdf2)-1),FUN=function(x){cdf2[x+1]-cdf2[x]})
          D[j,k] = D[k,j] = (KL.plugin(freqs1, freqs2)+KL.plugin(freqs2, freqs1))/2
        }
      }
    }
    ri=ristar
    zero=rep(0,K1)
    jags.data <- list("n"=ni, "Y"=ri, "D"=D, "K"=K1, "zero"=zero,'mu0'=log(((q0+q1)/2)/(1-(q0+q1)/2)))
    jags.fit <- jags.model(file = "cbhm_kl.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000,quiet=T)
    cbhm.out <- coda.samples(jags.fit,variable.names = c("theta0","phi","tausq","tausq2","e","p"),n.iter=10000,quiet=T)
    ### Testing & estimation
    posterior=numeric()
    for (k in 1:K1)
    {
      if (K1>1)
      {
        post.sample=cbhm.out[[1]][,paste("p[",k,"]",sep="")]
      }
      if (K1==1)
      {
        post.sample=cbhm.out[[1]][,"p"]
      }
      posterior[k]=sum(post.sample>q0)/length(post.sample)
    }
    posterior.cbhm.b[sim,stage2.cont]=posterior
  }
  print(sim)
}
Q.cbhm.kl=quantile(posterior.cbhm.b,1-alpha) # probability cut-off for the final decision




################################ the sim-th simulation

########### Stage 1 data:
stage1resp=simdata[[scenario]][((sim-1)*Ni+1):((sim-1)*Ni+Ni1),]
rik[1,]=colSums(stage1resp)

rikstar=rik
######### H distance matrix:
D=matrix(NA,K,K)
for (j in 1:(K-1))
{
  for (k in (j+1):K)
  {
    if ((rik[1,j]/nik[1,j])==(rik[1,k]/nik[1,k]))
    {
      rik[1,k]=rik[1,j] + epsilon
    }
  }
}
for (j in 1:K)
{
  for (k in j:K)
  {
    if (j == k)
    {
      D[j,k] = 0
    }
    if (j != k)
    {
      cdf1 = pbeta(seq(0,0.8,1/100),1+rik[1,j],1+nik[1,j]-rik[1,j])
      freqs1 = sapply(1:(length(cdf1)-1),FUN=function(x){cdf1[x+1]-cdf1[x]})
      cdf2 = pbeta(seq(0,0.8,1/100),1+rik[1,k],1+nik[1,k]-rik[1,k])
      freqs2 = sapply(1:(length(cdf2)-1),FUN=function(x){cdf2[x+1]-cdf2[x]})
      D[j,k] = D[k,j] = (KL.plugin(freqs1, freqs2)+KL.plugin(freqs2, freqs1))/2
    }
  }
}
rik=rikstar
zero=rep(0,K)

############ Jags model for Spatial BHM:
jags.data <- list("n"=nik[1,], "Y"=rik[1,], "D"=D, "K"=K, "zero"=zero,'mu0'=log(((q0+q1)/2)/(1-(q0+q1)/2)))
jags.fit <- jags.model(file = "cbhm_kl.txt",data = jags.data,
                       n.adapt=1000,n.chains=1,quiet=T)
update(jags.fit, 4000,quiet=T)
sbhm.out <- coda.samples(jags.fit,variable.names = c("theta0","phi","tausq","tausq2","e","p"),n.iter=10000,quiet=T)
### Interim analysis:
posterior=numeric()
for (k in 1:K)
{
  post.sample=sbhm.out[[1]][,paste("p[",k,"]",sep="")]
  posterior[k]=sum(post.sample>(q0+q1)/2)/length(post.sample)
}
## Futility stop:
stage2.stop=which(posterior<Qf)
stage2.cont=which(posterior>=Qf)
nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
decision.sbhm[sim,stage2.stop]=-1
if(length(stage2.stop)>0)
{
  Bias.sbhm[sim,stage2.stop]=abs(summary(sbhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"Mean"]-p0[stage2.stop])
  pest.sbhm[sim,stage2.stop]=summary(sbhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"Mean"]
}

## Stage 2:
if (length(stage2.cont)>0)
{
  rik[2,]=colSums(simdata[[scenario]][((sim-1)*Ni+Ni1+1):((sim-1)*Ni+Ni),] %*% diag(ifelse(nik[2,]>0,1,0),K))
  ri=colSums(as.matrix(rik[,stage2.cont]))
  ni=colSums(as.matrix(nik[,stage2.cont]))
  K1=length(stage2.cont)
  ristar=ri
  ni=colSums(as.matrix(nik[,stage2.cont]))
  K1=length(stage2.cont)
  R=D=matrix(NA,K1,K1)
  if (K1>1)
  {
    for (j in 1:(K1-1))
    {
      for (k in (j+1):K1)
      {
        if (ri[j]/ni[j]==ri[k]/ni[k])
        {
          ri[k] = ri[j] + epsilon
        }
      }
    }
  }
  
  for (j in 1:K1)
  {
    for (k in j:K1)
    {
      if (j == k)
      {
        D[j,k] = 0
      }
      if (j != k)
      {
        cdf1 = pbeta(seq(0,0.8,1/100),1+ri[j],1+ni[j]-ri[j])
        freqs1 = sapply(1:(length(cdf1)-1),FUN=function(x){cdf1[x+1]-cdf1[x]})
        cdf2 = pbeta(seq(0,0.8,1/100),1+ri[k],1+ni[k]-ri[k])
        freqs2 = sapply(1:(length(cdf2)-1),FUN=function(x){cdf2[x+1]-cdf2[x]})
        D[j,k] = D[k,j] = (KL.plugin(freqs1, freqs2)+KL.plugin(freqs2, freqs1))/2
      }
    }
  }
  ri=ristar
  zero=rep(0,K1)
  jags.data <- list("n"=ni, "Y"=ri, "D"=D, "K"=K1, "zero"=zero,'mu0'=log(((q0+q1)/2)/(1-(q0+q1)/2)))
  jags.fit <- jags.model(file = "cbhm_kl.txt",data = jags.data,
                         n.adapt=1000,n.chains=1,quiet=T)
  update(jags.fit, 4000,quiet=T)
  sbhm.out <- coda.samples(jags.fit,variable.names = c("theta0","phi","tausq","tausq2","tausq3","e","p"),n.iter=10000,quiet=T)
  
  ### Final decision:
  posterior=numeric()
  if (K1==1)
  {
    post.sample=sbhm.out[[1]][,"p"]
    posterior=sum(post.sample>q0)/length(post.sample)
    Bias.sbhm[sim,stage2.cont]=abs(summary(sbhm.out[[1]])[[1]]["p","Mean"]-p0[stage2.cont])
    pest.sbhm[sim,stage2.cont]=summary(sbhm.out[[1]])[[1]]["p","Mean"]
  }
  if (K1>1)
  {
    for (k in 1:K1)
    {
      post.sample=sbhm.out[[1]][,paste("p[",k,"]",sep="")]
      posterior[k]=sum(post.sample>q0)/length(post.sample)
    }
    Bias.sbhm[sim,stage2.cont]=abs(summary(sbhm.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",x,"]",sep="")}),"Mean"]-p0[stage2.cont])
    pest.sbhm[sim,stage2.cont]=summary(sbhm.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",x,"]",sep="")}),"Mean"]
  }
  decision.sbhm[sim,stage2.cont]=ifelse(posterior>Q.sbhm.kl,1,0)
}
