load("spatial-model/simulation.RData")
library(entropy)
library(rjags)
library(mvtnorm)
library(boot) # for inverse logit function
library(nonpar) # for Cochran's Q test
library(distrEx) #for TotalVarD

###### Overall simulation settings:
K=6 # number of indications
alpha=0.1
num.sim=5000
Ni=24
Ni1=14
nik=matrix(NA,2,K) # number of patients in indication k at state i
rik=matrix(NA,2,K) # number of responders in indication k at state i
nik[1,]=rep(Ni1,K) # number of patients enrolled at stage 1

q0=0.2
q1=0.4
##################################################
############### Calibration Stage: ###############
##################################################
p_null=rep(q0,K)
p0=p_null
Qf=0.05
########## Calibration for CBHM with B distance measure 
# & expnential correlation function:
decision=matrix(NA,num.sim,K)
posterior.sbhm.b=matrix(NA,num.sim,K)
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
        rik[1,k]=rik[1,j] + 0.1
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
        #R[j,k] = 1
      }
      if (j != k)
      {
        a1=1+rik[1,j]
        a2=1+rik[1,k]
        b1=1+nik[1,j]-rik[1,j]
        b2=1+nik[1,k]-rik[1,k]
        D[j,k] = D[k,j] = -log(beta((a1+a2)/2,(b1+b2)/2)/sqrt(beta(a1,b1)*beta(a2,b2)))
      }
    }
  }
  
  rik=rikstar
  zero=rep(0,K)
  ############ Jags model for Spatial BHM:
  jags.data <- list("n"=nik[1,], "Y"=rik[1,], "D"=D, "K"=K, "zero"=zero)
  jags.fit <- jags.model(file = "~/Google Drive/Sanofi/Jin/R/spatial-model/s-bhm.txt",data = jags.data,
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
  decision[sim,stage2.stop]=0
  posterior.sbhm.b[sim,stage2.stop]=0
  
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
            ri[k] = ri[j] + 0.1
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
          #R[j,k] = 1
        }
        if (j != k)
        {
          a1=1+ri[j]
          a2=1+ri[k]
          b1=1+ni[j]-ri[j]
          b2=1+ni[k]-ri[k]
          D[j,k] = D[k,j] = -log(beta((a1+a2)/2,(b1+b2)/2)/sqrt(beta(a1,b1)*beta(a2,b2)))
        }
      }
    }
    ri=ristar
    zero=rep(0,K1)
    jags.data <- list("n"=ni, "Y"=ri, "D"=D, "K"=K1, "zero"=zero)
    jags.fit <- jags.model(file = "~/Google Drive/Sanofi/Jin/R/spatial-model/s-bhm.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000,quiet=T)
    sbhm.out <- coda.samples(jags.fit,variable.names = c("theta0","phi","tausq","tausq2","e","p"),n.iter=10000,quiet=T)
    
    ### Final decision:
    posterior=numeric()
    for (k in 1:K1)
    {
      if (K1>1)
      {
        post.sample=sbhm.out[[1]][,paste("p[",k,"]",sep="")]
      }
      if (K1==1)
      {
        post.sample=sbhm.out[[1]][,"p"]
      }
      posterior[k]=sum(post.sample>q0)/length(post.sample)
    }
    decision[sim,stage2.cont]=ifelse(posterior>Q.sbhm.b,1,0)
    posterior.sbhm.b[sim,stage2.cont]=posterior
  }
  print(sim)
}
Q.sbhm.b=quantile(posterior.sbhm.b,0.9)   ###0.8716 


########## Calibration for EXNEX:
decision=matrix(NA,num.sim,K)
posterior.exnex=matrix(NA,num.sim,K)
for (sim in 1:num.sim)
{
  ##### Stage 1:
  rik[1,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[1,x],prob=p0[x])})
  jags.data <- list(
    "Nexch"=1, "Nmix"=2,
    "dirichprior" = c(1,1),
    # alternative weights for EX and EXNEX-2
    "Nstrata"=K,
    "n" = nik[1,],
    # original data
    "r" = rik[1,],
    # prior means and precisions for EX parameter mu
    "mu.mean"=c(0,0), "mu.prec"=c(0.2,0.2),
    # scale parameter of Half-Normal prior for tau
    "tau.HN.scale"=c(1,1),
    # NEX priors; make them strata-specific if needed
    "nex.mean"=log(q0/(1-q0)), "nex.prec"=0.15,
    "p.cut" = 0.3
  )
  jags.fit <- jags.model(file = "~/Google Drive/Sanofi/Jin/R/exnex.txt",data = jags.data,
                         n.adapt=1000,n.chains=1,quiet=T)
  update(jags.fit, 4000)
  exnex.out <- coda.samples(jags.fit,variable.names = c("tau","mu","p","theta","pMix","exch.index"),n.iter=10000)
  
  ### Interim analysis:
  posterior=numeric()
  for (k in 1:K)
  {
    post.sample=exnex.out[[1]][,paste("p[",k,"]",sep="")]
    posterior[k]=sum(post.sample>(q0+q1)/2)/length(post.sample)
  }
  ## Futility stop:
  stage2.stop=which(posterior<Qf)
  stage2.cont=which(posterior>=Qf)
  nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
  decision[sim,stage2.stop]=0
  posterior.exnex[sim,stage2.stop]=0
  
  ## Stage 2:
  if (length(stage2.cont)>0)
  {
    rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
    ri=colSums(as.matrix(rik[,stage2.cont]))
    ni=colSums(as.matrix(nik[,stage2.cont]))
    K1=length(stage2.cont)
    
    ############ Jags model for BHM:
    jags.data <- list(
      "Nexch"=1, "Nmix"=2,
      "dirichprior" = c(1,1),
      # alternative weights for EX and EXNEX-2
      "Nstrata"=K1,
      "n" = ni,
      # original data
      "r" = ri,
      # prior means and precisions for EX parameter mu
      "mu.mean"=c(0,0), "mu.prec"=c(0.2,0.2),
      # scale parameter of Half-Normal prior for tau
      "tau.HN.scale"=c(1,1),
      # NEX priors; make them strata-specific if needed
      "nex.mean"=log(q0/(1-q0)), "nex.prec"=0.15,
      "p.cut" = 0.3
    )
    jags.fit <- jags.model(file = "~/Google Drive/Sanofi/Jin/R/exnex.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000)
    exnex.out <- coda.samples(jags.fit,variable.names = c("tau","mu","p","theta","pMix","exch.index"),n.iter=10000)
    
    ### Final decision:
    posterior=numeric()
    if (K1==1)
    {
      post.sample=exnex.out[[1]][,"p"]
      posterior=sum(post.sample>q0)/length(post.sample)
    }
    if (K1>1)
    {
      for (k in 1:K1)
      {
        post.sample=exnex.out[[1]][,paste("p[",k,"]",sep="")]
        posterior[k]=sum(post.sample>q0)/length(post.sample)
      }
    }
    ## Final decision:
    decision[sim,stage2.cont]=ifelse(posterior>Q.exnex,1,0)
    posterior.exnex[sim,stage2.cont]=posterior
  }
  print(sim)
}
Q.exnex=quantile(posterior.exnex,0.9) 

########## Calibration for Liu's Method:
##### FFR for homogeneity test: 0.2 from Liu's paper

############### Tuning homogeneity path:
C=0.5 # threshold for futility stopping
decision=matrix(NA,num.sim,K)
posterior.liu=matrix(0,num.sim,K)
for (sim in 1:num.sim)
{
  ##### Stage 1:
  stage1resp=sapply(1:K,FUN=function(x){rbinom(n=nik[1,x],size=1,prob=p0[x])})
  rik[1,]=colSums(stage1resp)
  ################ Homogeneity Path: ################
  ### Interim analysis:
  posterior=numeric()
  for (k in 1:K)
  {
    post.pk=rbeta(5000,(0.5+rik[1,k]),(0.5+nik[1,k]-rik[1,k]))
    post.r2k=sapply(1:5000,FUN=function(x){rbinom(1,Ni-nik[1,k],post.pk[x])})
    post.rr=(rik[1,k]+post.r2k)/Ni
    posterior[k]=sum(post.rr>q0)/length(post.rr)
  }
  ## Futility stop:
  stage2.stop=which(posterior<C)
  stage2.cont=which(posterior>=C)
  nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
  
  decision[sim,stage2.stop]=0
  posterior.liu[sim,stage2.stop]=0
  
  ## Stage 2:
  if (length(stage2.cont)>0)
  {
    rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
    ri=colSums(as.matrix(rik))
    ni=colSums(as.matrix(nik))
    K1=length(stage2.cont)
    
    jags.data <- list("n"=ni, "Y"=ri, "K"=K)
    jags.fit <- jags.model(file = "~/Google Drive/Sanofi/Jin/R/liu.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000)
    liu.out <- coda.samples(jags.fit,variable.names = c("tausq","p","theta0"),n.iter=10000)
    
    ### Final decision:
    posterior=numeric()
    for (k in 1:K1)
    {
      post.sample=liu.out[[1]][,paste("p[",stage2.cont[k],"]",sep="")]
      # posterior probability that rr_k > q0_k:
      posterior[k]=sum(post.sample>q0)/length(post.sample)
    }
    ## Final decision:
    decision[sim,stage2.cont]=ifelse(posterior>Q.liu,1,0)
    posterior.liu[sim,stage2.cont]=posterior
  }
  print(sim)
}
Q.liu=quantile(posterior.liu,0.9) # 0.0940 0.0965 0.1125 0.0900 0.1035 0.1030

########## Calibration for BHM:
decision=matrix(NA,num.sim,K)
posterior.bhm=matrix(0,num.sim,K)
for (sim in 1:num.sim)
{
  ##### Stage 1:
  rik[1,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[1,x],prob=p0[x])})
  
  ############ Jags model for BHM:
  jags.data <- list("n"=nik[1,], "Y"=rik[1,], "K"=K)
  jags.fit <- jags.model(file = "~/Google Drive/Sanofi/Jin/R/spatial-model/BHM.txt",data = jags.data,
                         n.adapt=1000,n.chains=1,quiet=T)
  update(jags.fit, 4000)
  bhm.out <- coda.samples(jags.fit,variable.names = c("theta0","tausq","p"),n.iter=10000)
  ### Interim analysis:
  posterior=numeric()
  for (k in 1:K)
  {
    post.sample=bhm.out[[1]][,paste("p[",k,"]",sep="")]
    posterior[k]=sum(post.sample>(q0+q1)/2)/length(post.sample)
  }
  ## Futility stop:
  stage2.stop=which(posterior<Qf)
  stage2.cont=which(posterior>=Qf)
  nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
  decision[sim,stage2.stop]=0
  posterior.bhm[sim,stage2.stop]= posterior[stage2.stop]#0
  ## Stage 2:
  if (length(stage2.cont)>0)
  {
    rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
    ri=colSums(as.matrix(rik[,stage2.cont]))
    ristar=ri
    ni=colSums(as.matrix(nik[,stage2.cont]))
    K1=length(stage2.cont)
    
    ############ Jags model for BHM:
    jags.data <- list("n"=ni, "Y"=ri, "K"=K1)
    jags.fit <- jags.model(file = "~/Google Drive/Sanofi/Jin/R/spatial-model/BHM.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000)
    bhm.out <- coda.samples(jags.fit,variable.names = c("theta0","tausq","p"),n.iter=10000)
    
    ### Final decision:
    posterior=numeric()
    if (K1==1)
    {
      post.sample=bhm.out[[1]][,"p"]
      posterior=sum(post.sample>q0)/length(post.sample)
    }
    if (K1>1)
    {
      for (k in 1:K1)
      {
        post.sample=bhm.out[[1]][,paste("p[",k,"]",sep="")]
        posterior[k]=sum(post.sample>q0)/length(post.sample)
      }
    }
    ## Final decision:
    decision[sim,stage2.cont]=ifelse(posterior>Q.bhm,1,0)
    posterior.bhm[sim,stage2.cont]=posterior
  }
  print(sim)
}
Q.bhm=quantile(posterior.bhm,0.9)
sapply(1:6,FUN=function(x){sum(decision[,x]>0)})/num.sim #0.0860 0.0925 0.0975 0.0990 0.0985 0.1015

########## Calibration for Independent Analysis:
Q=Qcandidate[Q.index]
decision=matrix(NA,num.sim,K)
posterior.ind=matrix(0,num.sim,K)
for (sim in 1:num.sim)
{
  ##### Stage 1:
  rik[1,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[1,x],prob=p0[x])})
  
  ############ Jags model for BHM:
  jags.data <- list("n"=nik[1,], "Y"=rik[1,], "K"=K)
  jags.fit <- jags.model(file = "~/Google Drive/Sanofi/Jin/R/spatial-model/independent.txt",data = jags.data,
                         n.adapt=1000,n.chains=1,quiet=T)
  update(jags.fit, 4000)
  independent.out <- coda.samples(jags.fit,variable.names = c("p"),n.iter=10000)
  ### Interim analysis:
  posterior=numeric()
  for (k in 1:K)
  {
    post.sample=independent.out[[1]][,paste("p[",k,"]",sep="")]
    posterior[k]=sum(post.sample>(q0+q1)/2)/length(post.sample)
  }
  ## Futility stop:
  stage2.stop=which(posterior<Qf)
  stage2.cont=which(posterior>=Qf)
  #nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),(Ni-Ni1)*K/length(stage2.cont),0)})
  nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
  decision[sim,stage2.stop]=0
  posterior.ind[sim,stage2.stop]=posterior[stage2.stop] #0
  ## Stage 2:
  if (length(stage2.cont)>0)
  {
    rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
    ri=colSums(as.matrix(rik[,stage2.cont]))
    ristar=ri
    ni=colSums(as.matrix(nik[,stage2.cont]))
    K1=length(stage2.cont)
    
    ############ Jags model for BHM:
    jags.data <- list("n"=ni, "Y"=ri, "K"=K1)
    jags.fit <- jags.model(file = "~/Google Drive/Sanofi/Jin/R/spatial-model/independent.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000)
    independent.out <- coda.samples(jags.fit,variable.names = c("p"),n.iter=10000)
    
    ### Final decision:
    posterior=numeric()
    if (K1==1)
    {
      post.sample=independent.out[[1]][,"p"]
      posterior=sum(post.sample>q0)/length(post.sample)
    }
    if (K1>1)
    {
      for (k in 1:K1)
      {
        post.sample=independent.out[[1]][,paste("p[",k,"]",sep="")]
        posterior[k]=sum(post.sample>q0)/length(post.sample)
      }
    }
    ## Futility stop:
    decision[sim,stage2.cont]=ifelse(posterior>Q.independent,1,0)
    posterior.ind[sim,stage2.cont]=posterior
  }
  print(sim)
}
Q.independent=quantile(posterior.ind,1-alpha) 


