library(entropy)
library(rjags)
library(mvtnorm)
library(boot) # for inverse logit function
library(nonpar) # for Cochran's Q test
library(distrEx) #for TotalVarD
library(nonpar) # for function cochrans.q
K=6 # the number of indications
num.sim=5000 # the number of simulations per setting
Ni=24 # the maximum of total sample size for each indication group
Ni1=14 # stage-one sample size for each indication group
nik=matrix(NA,2,K) # each row records the number of patients in indication k at stage i
rik=matrix(NA,2,K) # each row records the number of responders in indication k at stage i
nik[1,]=rep(Ni1,K) # number of patients enrolled at stage 1
q0=0.2 # standard of care (null) response rate
q1=0.4 # target response rate
###### p.scenario: each row stores the true rrs in one simulation setting
p.scenario = t(sapply(0:(K-1),function(x){c(rep(q1,x),rep(q0,K-x))}))
methods=c("CBHM","EXNEX","Liu","BHM","Independent")

###### parameters
Qf=0.05 # probability cut-off for interim analysis
epsilon = 3*(q1-q0)/K # the small value added to the number of responsders for the indication groups that have equal sample rr
C=0.5 # threshold for futility stopping for Liu's two-stage method, follows Liu et al. (2017)

#### the following Qs are calibrated to have the same level of type I error control under the null scenario. 
#### Please see calibration.R for details.
Q.cbhm.b # the calibrated prob cut-off for final decision for CBHM using B distance
Q.exnex # the calibrated prob cut-off for final decision for EXNEX
Q.liu # the calibrated prob cut-off for final decision for Liu's two-stage method
Q.bhm # the calibrated prob cut-off for final decision for BHM
Q.independent # the calibrated prob cut-off for final decision for independent analysis
#############################################
############## Full Simulations #############
#############################################

##################### Create simulated data: #####################
simdata=list()
for (scenario in 1:nrow(p.scenario))
{
  simdata[[scenario]]=matrix(NA,num.sim*Ni,K)
  p0=p.scenario[scenario,]
  for (sim in 1:num.sim)
  {
    simdata[[scenario]][((sim-1)*Ni+1):(sim*Ni),]=sapply(1:K,FUN=function(x){rbinom(n=Ni,size=1,prob=p0[x])})
  }
}

##################### Simulation: #####################
OC=list() # store operating charasteristics
Decision=list()
decisions=list() # save the testing results for all simulations, optional
Bias=list() # store absolute bias
MSE=list() # store MSE
samplesize=matrix(0,nrow(p.scenario),length(methods)) # average sample size considering interim stop
for (scenario in 1:nrow(p.scenario))
{
  p0=p.scenario[scenario,]
  # decision*: store the decision of each simulation. -1: interim stop, 0: do not reject, 1: reject
  decision.cbhm=decision.exnex=decision.liu=decision.bhm=decision.independent=matrix(NA,num.sim,K)
  Bias.cbhm=Bias.exnex=Bias.liu=Bias.bhm=Bias.independent=matrix(0,num.sim,K)
  pest.cbhm=pest.exnex=pest.liu=pest.bhm=pest.independent=matrix(0,num.sim,K)
  MSE.cbhm=MSE.exnex=MSE.liu=MSE.bhm=MSE.independent=matrix(0,num.sim,K)
  samplesize[scenario,]=rep(0,length(methods))
  Decision[[scenario]]=matrix(0,length(methods),5)
  colnames(Decision[[scenario]])=c("% Perfect","# TP","# TN","# FP","# FN")
  rownames(Decision[[scenario]])=methods
  tp=which(p0>=q1)
  tn=which(p0<q1)
  
  for (sim in 1:num.sim)
  {
    ########### Stage 1 data:
    stage1resp=simdata[[scenario]][((sim-1)*Ni+1):((sim-1)*Ni+Ni1),]
    #stage1resp=sapply(1:K,FUN=function(x){rbinom(n=nik[1,x],size=1,prob=p0[x])})
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
    jags.data <- list("n"=nik[1,], "Y"=rik[1,], "D"=D, "K"=K, "zero"=zero,'mu0'=log(((q0+q1)/2)/(1-(q0+q1)/2)))
    jags.fit <- jags.model(file = "cbhm.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000,quiet=T)
    cbhm.out <- coda.samples(jags.fit,variable.names = c("theta0","phi","tausq","tausq2","e","p"),n.iter=10000,quiet=T)
    ### Interim analysis:
    posterior=numeric()
    for (k in 1:K)
    {
      post.sample=cbhm.out[[1]][,paste0("p[",k,"]")]
      posterior[k]=sum(post.sample>(q0+q1)/2)/length(post.sample)
    }
    ## Futility stop:
    stage2.stop=which(posterior<Qf)
    stage2.cont=which(posterior>=Qf)
    nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
    decision.cbhm[sim,stage2.stop]=-1
    if(length(stage2.stop)>0)
    {
      Bias.cbhm[sim,stage2.stop]=abs(summary(cbhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste0("p[",stage2.stop[x],"]")}),"Mean"]-p0[stage2.stop])
      pest.cbhm[sim,stage2.stop]=summary(cbhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste0("p[",stage2.stop[x],"]")}),"Mean"]
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
      jags.data <- list("n"=ni, "Y"=ri, "D"=D, "K"=K1, "zero"=zero,'mu0'=log(((q0+q1)/2)/(1-(q0+q1)/2)))
      jags.fit <- jags.model(file = "cbhm.txt",data = jags.data,
                             n.adapt=1000,n.chains=1,quiet=T)
      update(jags.fit, 4000,quiet=T)
      cbhm.out <- coda.samples(jags.fit,variable.names = c("theta0","phi","tausq","tausq2","tausq3","e","p"),n.iter=10000,quiet=T)
      
      ### Final decision:
      posterior=numeric()
      if (K1==1)
      {
        post.sample=cbhm.out[[1]][,"p"]
        posterior=sum(post.sample>q0)/length(post.sample)
        Bias.cbhm[sim,stage2.cont]=abs(summary(cbhm.out[[1]])[[1]]["p","Mean"]-p0[stage2.cont])
        pest.cbhm[sim,stage2.cont]=summary(cbhm.out[[1]])[[1]]["p","Mean"]
      }
      if (K1>1)
      {
        for (k in 1:K1)
        {
          post.sample=cbhm.out[[1]][,paste0("p[",k,"]")]
          posterior[k]=sum(post.sample>q0)/length(post.sample)
        }
        Bias.cbhm[sim,stage2.cont]=abs(summary(cbhm.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste0("p[",x,"]")}),"Mean"]-p0[stage2.cont])
        pest.cbhm[sim,stage2.cont]=summary(cbhm.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste0("p[",x,"]")}),"Mean"]
      }
      decision.cbhm[sim,stage2.cont]=ifelse(posterior>Q.cbhm.b,1,0)
    }
    ##### total sample size:
    ni.cbhm=sum(nik)
    
    
    ########################### EXNEX
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
    jags.fit <- jags.model(file = "exnex.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000)
    exnex.out <- coda.samples(jags.fit,variable.names = c("tau","mu","p","theta","pMix","exch.index"),n.iter=10000)
    ### Interim analysis:
    posterior=numeric()
    for (k in 1:K)
    {
      post.sample=exnex.out[[1]][,paste0("p[",k,"]")]
      posterior[k]=sum(post.sample>(q0+q1)/2)/length(post.sample)
    }
    ## Futility stop:
    stage2.stop=which(posterior<Qf)
    stage2.cont=which(posterior>=Qf)
    nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
    decision.exnex[sim,stage2.stop]=-1
    
    if(length(stage2.stop)>0)
    {
      Bias.exnex[sim,stage2.stop]=abs(summary(exnex.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste0("p[",stage2.stop[x],"]")}),"Mean"]-p0[stage2.stop])
      pest.exnex[sim,stage2.stop]=summary(exnex.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste0("p[",stage2.stop[x],"]")}),"Mean"]
    }
    
    ## Stage 2:
    if (length(stage2.cont)>0)
    {
      rik[2,]=colSums(simdata[[scenario]][((sim-1)*Ni+Ni1+1):((sim-1)*Ni+Ni),] %*% diag(ifelse(nik[2,]>0,1,0),K))
      #rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
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
        Bias.exnex[sim,stage2.cont]=abs(summary(exnex.out[[1]])[[1]]["p","Mean"]-p0[stage2.cont])
        pest.exnex[sim,stage2.cont]=summary(exnex.out[[1]])[[1]]["p","Mean"]
      }
      if (K1>1)
      {
        for (k in 1:K1)
        {
          post.sample=exnex.out[[1]][,paste0("p[",k,"]")]
          posterior[k]=sum(post.sample>q0)/length(post.sample)
        }
        Bias.exnex[sim,stage2.cont]=abs(summary(exnex.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste0("p[",x,"]")}),"Mean"]-p0[stage2.cont])
        pest.exnex[sim,stage2.cont]=summary(exnex.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste0("p[",x,"]")}),"Mean"]
      }
      ## Final decision:
      decision.exnex[sim,stage2.cont]=ifelse(posterior>Q.exnex,1,0)
    }
    ni.exnex=sum(nik)
    
    #################### Liu's Design: ####################
    ########## Assess heterogeneity:
    testdata=stage1resp[1:7,] # 7 x K matrix
    if (sum(testdata)==0) {testdata=stage1resp[sample(1:Ni1,7),]}
    cochransq=cochrans.q(testdata)
    temp=strsplit(cochransq@PVal,split=" ")[[1]]
    pval=as.numeric(temp[length(temp)]) # extract the p-value from the output
    if (pval>=0.2)
    {
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
      
      decision.liu[sim,stage2.stop]=-1
      if(length(stage2.stop)>0)
      {
        for (k in stage2.stop)
        {
          post.pk=rbeta(5000,(0.5+rik[1,k]),(0.5+nik[1,k]-rik[1,k]))
          post.r2k=sapply(1:5000,FUN=function(x){rbinom(1,Ni-nik[1,k],post.pk[x])})
          post.rr=(rik[1,k]+post.r2k)/Ni
          Bias.liu[sim,k]=abs(mean(post.rr)-p0[k])
          pest.liu[sim,k]=mean(post.rr)
        }
      }
      
      ## Stage 2:
      if (length(stage2.cont)>0)
      {
        #rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
        #ri=colSums(as.matrix(rik[,stage2.cont]))
        #ni=colSums(as.matrix(nik[,stage2.cont]))
        #K1=length(stage2.cont)
        
        rik[2,]=colSums(simdata[[scenario]][((sim-1)*Ni+Ni1+1):((sim-1)*Ni+Ni),] %*% diag(ifelse(nik[2,]>0,1,0),K))
        #rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
        ri=colSums(as.matrix(rik))
        ni=colSums(as.matrix(nik))
        K1=length(stage2.cont)
        ############ Jags model for BHMM:
        jags.data <- list("n"=ni, "Y"=ri, "K"=K)
        jags.fit <- jags.model(file = "liu.txt",data = jags.data,
                               n.adapt=1000,n.chains=1,quiet=T)
        update(jags.fit, 4000)
        liu.out <- coda.samples(jags.fit,variable.names = c("tausq","p","theta0"),n.iter=10000)
        
        ### Final decision:
        posterior=numeric()
        for (k in 1:K1)
        {
          post.sample=liu.out[[1]][,paste0("p[",stage2.cont[k],"]")]
          posterior[k]=sum(post.sample>q0)/length(post.sample)
        }
        Bias.liu[sim,stage2.cont]=abs(summary(liu.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste0("p[",stage2.cont[x],"]")}),"Mean"]-p0[stage2.cont])
        pest.liu[sim,stage2.cont]=summary(liu.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste0("p[",stage2.cont[x],"]")}),"Mean"]
        decision.liu[sim,stage2.cont]=ifelse(posterior>Q.liu,1,0)
      }
    }
    
    if (pval<0.2)
    {
      ######### Heterogeneity Path: Simon's 2-stage Design
      ### Interim analysis:
      posterior=rik[1,]/nik[1,]
      ## Futility stop:
      stage2.stop=which(rik[1,]<2)
      stage2.cont=which(rik[1,]>=2)
      nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
      if(length(stage2.stop)>0)
      {
        Bias.liu[sim,stage2.stop]=abs(posterior[stage2.stop]-p0[stage2.stop])
        pest.liu[sim,stage2.stop]=posterior[stage2.stop]
      }
      decision.liu[sim,stage2.stop]=-1
      ## Stage 2:
      if (length(stage2.cont)>0)
      {
        rik[2,]=colSums(simdata[[scenario]][((sim-1)*Ni+Ni1+1):((sim-1)*Ni+Ni),] %*% diag(ifelse(nik[2,]>0,1,0),K))
        ri=colSums(as.matrix(rik[,stage2.cont]))
        ni=colSums(as.matrix(nik[,stage2.cont]))
        K1=length(stage2.cont)
        
        ############ Jags model for BHM:
        ### Final decision:
        posterior=ri/ni
        decision.liu[sim,stage2.cont]=ifelse(ri>7,1,0)
        Bias.liu[sim,stage2.cont]=abs(posterior-p0[stage2.cont])
        pest.liu[sim,stage2.cont]=posterior
      }
    }
    ni.liu=sum(nik)
    
    
    ########################## BHM: ##########################
    ############ Jags model for BHM:
    jags.data <- list("n"=nik[1,], "Y"=rik[1,], "K"=K)
    jags.fit <- jags.model(file = "BHM.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000)
    bhm.out <- coda.samples(jags.fit,variable.names = c("theta0","tausq","p"),n.iter=10000)
    ### Interim analysis:
    posterior=numeric()
    for (k in 1:K)
    {
      post.sample=bhm.out[[1]][,paste0("p[",k,"]")]
      posterior[k]=sum(post.sample>(q0+q1)/2)/length(post.sample)
    }
    ## Futility stop:
    stage2.stop=which(posterior<Qf)
    stage2.cont=which(posterior>=Qf)
    nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
    decision.bhm[sim,stage2.stop]=-1
    if(length(stage2.stop)>0)
    {
      Bias.bhm[sim,stage2.stop]=abs(summary(bhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste0("p[",stage2.stop[x],"]")}),"Mean"]-p0[stage2.stop])
      pest.bhm[sim,stage2.stop]=summary(bhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste0("p[",stage2.stop[x],"]")}),"Mean"]
    }
    ## Stage 2:
    if (length(stage2.cont)>0)
    {
      rik[2,]=colSums(simdata[[scenario]][((sim-1)*Ni+Ni1+1):((sim-1)*Ni+Ni),] %*% diag(ifelse(nik[2,]>0,1,0),K))
      ri=colSums(as.matrix(rik[,stage2.cont]))
      ristar=ri
      ni=colSums(as.matrix(nik[,stage2.cont]))
      K1=length(stage2.cont)
      
      ############ Jags model for BHM:
      jags.data <- list("n"=ni, "Y"=ri, "K"=K1)
      jags.fit <- jags.model(file = "BHM.txt",data = jags.data,
                             n.adapt=1000,n.chains=1,quiet=T)
      update(jags.fit, 4000)
      bhm.out <- coda.samples(jags.fit,variable.names = c("theta0","tausq","p"),n.iter=10000)
      
      ### Final decision:
      posterior=numeric()
      if (K1==1)
      {
        post.sample=bhm.out[[1]][,"p"]
        posterior=sum(post.sample>q0)/length(post.sample)
        Bias.bhm[sim,stage2.cont]=abs(summary(bhm.out[[1]])[[1]]["p","Mean"]-p0[stage2.cont])
        pest.bhm[sim,stage2.cont]=summary(bhm.out[[1]])[[1]]["p","Mean"]
      }
      if (K1>1)
      {
        for (k in 1:K1)
        {
          post.sample=bhm.out[[1]][,paste0("p[",k,"]")]
          posterior[k]=sum(post.sample>q0)/length(post.sample)
        }
        Bias.bhm[sim,stage2.cont]=abs(summary(bhm.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste0("p[",x,"]")}),"Mean"]-p0[stage2.cont])
        pest.bhm[sim,stage2.cont]=summary(bhm.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste0("p[",x,"]")}),"Mean"]
      }
      decision.bhm[sim,stage2.cont]=ifelse(posterior>Q.bhm,1,0)
    }
    ##### total sample size:
    ni.bhm=sum(nik)
    
    ############################## Independent analysis:
    ############ Jags model:
    jags.data <- list("n"=nik[1,], "Y"=rik[1,], "K"=K)
    jags.fit <- jags.model(file = "independent.txt",data = jags.data,
                           n.adapt=1000,n.chains=1,quiet=T)
    update(jags.fit, 4000)
    independent.out <- coda.samples(jags.fit,variable.names = c("p"),n.iter=10000)
    ### Interim analysis:
    posterior=numeric()
    for (k in 1:K)
    {
      post.sample=independent.out[[1]][,paste0("p[",k,"]")]
      posterior[k]=sum(post.sample>(q0+q1)/2)/length(post.sample)
    }
    ## Futility stop:
    stage2.stop=which(posterior<Qf)
    stage2.cont=which(posterior>=Qf)
    nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
    decision.independent[sim,stage2.stop]=-1
    if(length(stage2.stop)>0)
    {
      Bias.independent[sim,stage2.stop]=abs(summary(independent.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste0("p[",stage2.stop[x],"]")}),"Mean"]-p0[stage2.stop])
      pest.independent[sim,stage2.stop]=summary(independent.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste0("p[",stage2.stop[x],"]")}),"Mean"]
    }
    ## Stage 2:
    if (length(stage2.cont)>0)
    {
      rik[2,]=colSums(simdata[[scenario]][((sim-1)*Ni+Ni1+1):((sim-1)*Ni+Ni),] %*% diag(ifelse(nik[2,]>0,1,0),K))
      ri=colSums(as.matrix(rik[,stage2.cont]))
      ristar=ri
      ni=colSums(as.matrix(nik[,stage2.cont]))
      K1=length(stage2.cont)
      
      ############ Jags model for BHM:
      jags.data <- list("n"=ni, "Y"=ri, "K"=K1)
      jags.fit <- jags.model(file = "independent.txt",data = jags.data,
                             n.adapt=1000,n.chains=1,quiet=T)
      update(jags.fit, 4000)
      independent.out <- coda.samples(jags.fit,variable.names = c("p"),n.iter=10000)
      
      ### Final decision:
      posterior=numeric()
      if (K1==1)
      {
        post.sample=independent.out[[1]][,"p"]
        posterior=sum(post.sample>q0)/length(post.sample)
        ##### Bias:
        Bias.independent[sim,stage2.cont]=abs(summary(independent.out[[1]])[[1]]["Mean"]-p0[stage2.cont])
        pest.independent[sim,stage2.cont]=summary(independent.out[[1]])[[1]]["Mean"]
        #MSE.independent[sim,stage2.cont]=Bias.independent[sim,stage2.cont]^2 + (summary(independent.out[[1]])[[1]][,"SD"])^2
      }
      if (K1>1)
      {
        for (k in 1:K1)
        {
          post.sample=independent.out[[1]][,paste0("p[",k,"]")]
          posterior[k]=sum(post.sample>q0)/length(post.sample)
        }
        ##### Bias:
        Bias.independent[sim,stage2.cont]=abs(summary(independent.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste0("p[",x,"]")}),"Mean"]-p0[stage2.cont])
        pest.independent[sim,stage2.cont]=summary(independent.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste0("p[",x,"]")}),"Mean"]
      }
      decision.independent[sim,stage2.cont]=ifelse(posterior>Q.independent,1,0)############ Jags model for BHM:
    }
    ni.independent=sum(nik)
    
    # update average sample size
    samplesize[scenario,]=samplesize[scenario,]+c(ni.cbhm,ni.exnex,ni.liu,ni.bhm,ni.independent)
    
    # Summarize TP, TN, FP, FN:
    Decision[[scenario]]['CBHM',"% Perfect"]=Decision[[scenario]]['CBHM',"% Perfect"] + 1*(sum(c(decision.cbhm[sim,tp]==1,decision.cbhm[sim,tn]<=0))==K)
    Decision[[scenario]]['CBHM',"# TP"]=Decision[[scenario]]['CBHM',"# TP"] + sum(decision.cbhm[sim,tp]==1)
    Decision[[scenario]]['CBHM',"# TN"]=Decision[[scenario]]['CBHM',"# TN"] +  sum(decision.cbhm[sim,tn]<=0)
    Decision[[scenario]]['CBHM',"# FP"]=Decision[[scenario]]['CBHM',"# FP"] +  sum(decision.cbhm[sim,tn]==1)
    Decision[[scenario]]['CBHM',"# FN"]=Decision[[scenario]]['CBHM',"# FN"] +  sum(decision.cbhm[sim,tp]<=0)
    
    Decision[[scenario]]['EXNEX',"% Perfect"]=Decision[[scenario]]['EXNEX',"% Perfect"] + 1*(sum(c(decision.exnex[sim,tp]==1,decision.exnex[sim,tn]<=0))==K)
    Decision[[scenario]]['EXNEX',"# TP"]=Decision[[scenario]]['EXNEX',"# TP"] + sum(decision.exnex[sim,tp]==1)
    Decision[[scenario]]['EXNEX',"# TN"]=Decision[[scenario]]['EXNEX',"# TN"] +  sum(decision.exnex[sim,tn]<=0)
    Decision[[scenario]]['EXNEX',"# FP"]=Decision[[scenario]]['EXNEX',"# FP"] +  sum(decision.exnex[sim,tn]==1)
    Decision[[scenario]]['EXNEX',"# FN"]=Decision[[scenario]]['EXNEX',"# FN"] +  sum(decision.exnex[sim,tp]<=0)
    
    Decision[[scenario]]['Liu',"% Perfect"]=Decision[[scenario]]['Liu',"% Perfect"] + 1*(sum(c(decision.liu[sim,tp]==1,decision.liu[sim,tn]<=0))==K)
    Decision[[scenario]]['Liu',"# TP"]=Decision[[scenario]]['Liu',"# TP"] + sum(decision.liu[sim,tp]==1)
    Decision[[scenario]]['Liu',"# TN"]=Decision[[scenario]]['Liu',"# TN"] +  sum(decision.liu[sim,tn]<=0)
    Decision[[scenario]]['Liu',"# FP"]=Decision[[scenario]]['Liu',"# FP"] +  sum(decision.liu[sim,tn]==1)
    Decision[[scenario]]['Liu',"# FN"]=Decision[[scenario]]['Liu',"# FN"] +  sum(decision.liu[sim,tp]<=0)
    
    Decision[[scenario]]['BHM',"% Perfect"]=Decision[[scenario]]['BHM',"% Perfect"] + 1*(sum(c(decision.bhm[sim,tp]==1,decision.bhm[sim,tn]<=0))==K)
    Decision[[scenario]]['BHM',"# TP"]=Decision[[scenario]]['BHM',"# TP"] + sum(decision.bhm[sim,tp]==1)
    Decision[[scenario]]['BHM',"# TN"]=Decision[[scenario]]['BHM',"# TN"] +  sum(decision.bhm[sim,tn]<=0)
    Decision[[scenario]]['BHM',"# FP"]=Decision[[scenario]]['BHM',"# FP"] +  sum(decision.bhm[sim,tn]==1)
    Decision[[scenario]]['BHM',"# FN"]=Decision[[scenario]]['BHM',"# FN"] +  sum(decision.bhm[sim,tp]<=0)
    
    Decision[[scenario]]['Independent',"% Perfect"]=Decision[[scenario]]['Independent',"% Perfect"] + 1*(sum(c(decision.independent[sim,tp]==1,decision.independent[sim,tn]<=0))==K)
    Decision[[scenario]]['Independent',"# TP"]=Decision[[scenario]]['Independent',"# TP"] + sum(decision.independent[sim,tp]==1)
    Decision[[scenario]]['Independent',"# TN"]=Decision[[scenario]]['Independent',"# TN"] +  sum(decision.independent[sim,tn]<=0)
    Decision[[scenario]]['Independent',"# FP"]=Decision[[scenario]]['Independent',"# FP"] +  sum(decision.independent[sim,tn]==1)
    Decision[[scenario]]['Independent',"# FN"]=Decision[[scenario]]['Independent',"# FN"] +  sum(decision.independent[sim,tp]<=0)
    
    print(sim)
    # optional: save the results for all simulations
    #decisions[[scenario]]=list(decision.cbhm,decision.exnex,decision.liu,decision.bhm,decision.independent)
  }
  
  Decision[[scenario]]=Decision[[scenario]]/num.sim # average across all num.sim simulations
  
  OC[[scenario]]=matrix(NA,2*length(methods),K)
  colnames(OC[[scenario]])=sapply(1:ncol(OC[[scenario]]),FUN=function(x){paste0("cancer ",x)})
  rownames(OC[[scenario]])=sapply(1:nrow(OC[[scenario]]),FUN=function(x){paste0(methods[ceiling(x/2)]," - % ",ifelse(x/2!=floor(x/2),"reject","stop"))})
  OC[[scenario]][paste0(methods[1]," - % ","reject"),]=sapply(1:K,FUN=function(x){sum(decision.cbhm[,x]==1)/num.sim*100})
  OC[[scenario]][paste0(methods[1]," - % ","stop"),]=sapply(1:K,FUN=function(x){sum(decision.cbhm[,x]==-1)/num.sim*100})
  OC[[scenario]][paste0(methods[2]," - % ","reject"),]=sapply(1:K,FUN=function(x){sum(decision.exnex[,x]==1)/num.sim*100})
  OC[[scenario]][paste0(methods[2]," - % ","stop"),]=sapply(1:K,FUN=function(x){sum(decision.exnex[,x]==-1)/num.sim*100})
  OC[[scenario]][paste0(methods[3]," - % ","reject"),]=sapply(1:K,FUN=function(x){sum(decision.liu[,x]==1)/num.sim*100})
  OC[[scenario]][paste0(methods[3]," - % ","stop"),]=sapply(1:K,FUN=function(x){sum(decision.liu[,x]==-1)/num.sim*100})
  OC[[scenario]][paste0(methods[4]," - % ","reject"),]=sapply(1:K,FUN=function(x){sum(decision.bhm[,x]==1)/num.sim*100})
  OC[[scenario]][paste0(methods[4]," - % ","stop"),]=sapply(1:K,FUN=function(x){sum(decision.bhm[,x]==-1)/num.sim*100})
  OC[[scenario]][paste0(methods[5]," - % ","reject"),]=sapply(1:K,FUN=function(x){sum(decision.independent[,x]==1)/num.sim*100})
  OC[[scenario]][paste0(methods[5]," - % ","stop"),]=sapply(1:K,FUN=function(x){sum(decision.independent[,x]==-1)/num.sim*100})
  samplesize[scenario,]=samplesize[scenario,]/num.sim
  OC[[scenario]]=round(OC[[scenario]],1)
  Decision[[scenario]]=round(Decision[[scenario]],3)
  samplesize=round(samplesize,1)
  
  Bias[[scenario]]=matrix(NA,length(methods),K)
  colnames(Bias[[scenario]])=sapply(1:ncol(Bias[[scenario]]),FUN=function(x){paste0("cancer ",x)})
  rownames(Bias[[scenario]])=methods
  MSE[[scenario]]=matrix(NA,length(methods),K)
  colnames(MSE[[scenario]])=sapply(1:ncol(MSE[[scenario]]),FUN=function(x){paste0("cancer ",x)})
  rownames(MSE[[scenario]])=methods
  Bias[[scenario]][methods[1],]=colMeans(Bias.cbhm)
  Bias[[scenario]][methods[2],]=colMeans(Bias.exnex)# !!!!! 1001:2000
  Bias[[scenario]][methods[3],]=colMeans(Bias.liu)
  Bias[[scenario]][methods[4],]=colMeans(Bias.bhm)
  Bias[[scenario]][methods[5],]=colMeans(Bias.independent)
  Bias[[scenario]]=round(Bias[[scenario]],4)
  MSE[[scenario]][methods[1],]=Bias[[scenario]][methods[1],]^2 + apply(pest.cbhm,2,var)
  MSE[[scenario]][methods[2],]=Bias[[scenario]][methods[2],]^2 + apply(pest.exnex,2,var)
  MSE[[scenario]][methods[3],]=Bias[[scenario]][methods[3],]^2 + apply(pest.liu,2,var)
  MSE[[scenario]][methods[4],]=Bias[[scenario]][methods[4],]^2 + apply(pest.bhm,2,var)
  MSE[[scenario]][methods[5],]=Bias[[scenario]][methods[5],]^2 + apply(pest.independent,2,var)
  MSE[[scenario]]=round(MSE[[scenario]],4)
}
