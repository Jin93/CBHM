#############################################
############## Full Simulations #############
#############################################
##### each row: one setting of true RRs
p.scenario=matrix(c(rep(0.2,K),c(rep(0.4,4),rep(0.2,2)),c(rep(0.4,2),rep(0.2,4)),
                    c(rep(0.4,3),rep(0.2,3)),c(0.4,rep(0.2,5)),
                    c(0.5,rep(0.4,4),0.2),rep(0.4,6)),7,K,byrow=T)
methods=c("S-BHM","EXNEX","Liu","BHM","Independent")



##################### Create simulated data: #####################
K=6 # number of indications
num.sim=5000
Ni=24
Ni1=14
simdata=list()
for (scenario in 1:8)
{
  simdata[[scenario]]=matrix(NA,num.sim*Ni,K)
  p0=p.scenario[scenario,]
  for (sim in 1:num.sim)
  {
    simdata[[scenario]][((sim-1)*Ni+1):(sim*Ni),]=sapply(1:K,FUN=function(x){rbinom(n=Ni,size=1,prob=p0[x])})
  }
}

##################### Simulation: #####################
OC=list()
Decision=list()
decisons=list()
Bias=list()
MSE=list()
samplesize=matrix(0,nrow(p.scenario),length(methods))
for (scenario in 1:9) # needs to run by method: 1(bhm) 2(exnex,liu),3(exnex,liu,bhm,ind), 5(bhm)
{
  p0=p.scenario[scenario,]
  decision.sbhm=decision.exnex=decision.liu=decision.bhm=decision.independent=matrix(NA,num.sim,K)
  Bias.sbhm=Bias.exnex=Bias.liu=Bias.bhm=Bias.independent=matrix(0,num.sim,K)
  pest.sbhm=pest.exnex=pest.liu=pest.bhm=pest.independent=matrix(0,num.sim,K)
  MSE.sbhm=MSE.exnex=MSE.liu=MSE.bhm=MSE.independent=matrix(0,num.sim,K)
  samplesize[scenario,]=rep(0,length(methods))
  Decision[[scenario]]=matrix(0,length(methods),5)
  colnames(Decision[[scenario]])=c("% Perfect","# TP","# TN","# FP","# FN")
  rownames(Decision[[scenario]])=methods
  tp=which(p0>=0.4)
  tn=which(p0<0.4)
  
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
    decision.sbhm[sim,stage2.stop]=-1
    if(length(stage2.stop)>0)
    {
      Bias.sbhm[sim,stage2.stop]=abs(summary(sbhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"Mean"]-p0[stage2.stop])
      #MSE.sbhm[sim,stage2.stop]=Bias.sbhm[sim,stage2.stop]^2 + (summary(sbhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"SD"])^2
      pest.sbhm[sim,stage2.stop]=summary(sbhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"Mean"]
    }
    
    ## Stage 2:
    if (length(stage2.cont)>0)
    {
      rik[2,]=colSums(simdata[[scenario]][((sim-1)*Ni+Ni1+1):((sim-1)*Ni+Ni),] %*% diag(ifelse(nik[2,]>0,1,0),K))
      #rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
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
      decision.sbhm[sim,stage2.cont]=ifelse(posterior>Q.sbhm.b,1,0)
    }
    ##### total sample size:
    ni.sbhm=sum(nik)
    
    
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
    decision.exnex[sim,stage2.stop]=-1
    
    if(length(stage2.stop)>0)
    {
      Bias.exnex[sim,stage2.stop]=abs(summary(exnex.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"Mean"]-p0[stage2.stop])
      pest.exnex[sim,stage2.stop]=summary(exnex.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"Mean"]
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
          post.sample=exnex.out[[1]][,paste("p[",k,"]",sep="")]
          posterior[k]=sum(post.sample>q0)/length(post.sample)
        }
        Bias.exnex[sim,stage2.cont]=abs(summary(exnex.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",x,"]",sep="")}),"Mean"]-p0[stage2.cont])
        pest.exnex[sim,stage2.cont]=summary(exnex.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",x,"]",sep="")}),"Mean"]
      }
      ## Final decision:
      decision.exnex[sim,stage2.cont]=ifelse(posterior>Q.exnex,1,0)
    }
    ni.exnex=sum(nik)
    
    #################### Liu's Design: ####################
    ########## Assess heterogeneity:
    testdata=stage1resp[1:7,] # 7 x K matrix.
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
        jags.fit <- jags.model(file = "~/Google Drive/Sanofi/Jin/R/liu.txt",data = jags.data,
                               n.adapt=1000,n.chains=1,quiet=T)
        update(jags.fit, 4000)
        liu.out <- coda.samples(jags.fit,variable.names = c("tausq","p","theta0"),n.iter=10000)
        
        ### Final decision:
        posterior=numeric()
        for (k in 1:K1)
        {
          post.sample=liu.out[[1]][,paste("p[",stage2.cont[k],"]",sep="")]
          posterior[k]=sum(post.sample>q0)/length(post.sample)
        }
        Bias.liu[sim,stage2.cont]=abs(summary(liu.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",stage2.cont[x],"]",sep="")}),"Mean"]-p0[stage2.cont])
        pest.liu[sim,stage2.cont]=summary(liu.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",stage2.cont[x],"]",sep="")}),"Mean"]
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
        #rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
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
    decision.bhm[sim,stage2.stop]=-1
    if(length(stage2.stop)>0)
    {
      Bias.bhm[sim,stage2.stop]=abs(summary(bhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"Mean"]-p0[stage2.stop])
      pest.bhm[sim,stage2.stop]=summary(bhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"Mean"]
      #MSE.bhm[sim,stage2.stop]=Bias.bhm[sim,stage2.stop]^2 + (summary(bhm.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"SD"])^2
    }
    ## Stage 2:
    if (length(stage2.cont)>0)
    {
      rik[2,]=colSums(simdata[[scenario]][((sim-1)*Ni+Ni1+1):((sim-1)*Ni+Ni),] %*% diag(ifelse(nik[2,]>0,1,0),K))
      #rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
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
        Bias.bhm[sim,stage2.cont]=abs(summary(bhm.out[[1]])[[1]]["p","Mean"]-p0[stage2.cont])
        pest.bhm[sim,stage2.cont]=summary(bhm.out[[1]])[[1]]["p","Mean"]
      }
      if (K1>1)
      {
        for (k in 1:K1)
        {
          post.sample=bhm.out[[1]][,paste("p[",k,"]",sep="")]
          posterior[k]=sum(post.sample>q0)/length(post.sample)
        }
        Bias.bhm[sim,stage2.cont]=abs(summary(bhm.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",x,"]",sep="")}),"Mean"]-p0[stage2.cont])
        pest.bhm[sim,stage2.cont]=summary(bhm.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",x,"]",sep="")}),"Mean"]
      }
      decision.bhm[sim,stage2.cont]=ifelse(posterior>Q.bhm,1,0)
    }
    ##### total sample size:
    ni.bhm=sum(nik)
    
    ############################## Independent analysis:
    ############ Jags model:
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
    nik[2,]=sapply(1:K,FUN=function(x){ifelse(is.element(x,stage2.cont),Ni-Ni1,0)})
    decision.independent[sim,stage2.stop]=-1
    if(length(stage2.stop)>0)
    {
      Bias.independent[sim,stage2.stop]=abs(summary(independent.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"Mean"]-p0[stage2.stop])
      pest.independent[sim,stage2.stop]=summary(independent.out[[1]])[[1]][sapply(1:length(stage2.stop),FUN=function(x){paste("p[",stage2.stop[x],"]",sep="")}),"Mean"]
    }
    ## Stage 2:
    if (length(stage2.cont)>0)
    {
      rik[2,]=colSums(simdata[[scenario]][((sim-1)*Ni+Ni1+1):((sim-1)*Ni+Ni),] %*% diag(ifelse(nik[2,]>0,1,0),K))
      #rik[2,]=sapply(1:K,FUN=function(x){rbinom(n=1,size=nik[2,x],prob=p0[x])})
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
        ##### Bias:
        Bias.independent[sim,stage2.cont]=abs(summary(independent.out[[1]])[[1]]["Mean"]-p0[stage2.cont])
        pest.independent[sim,stage2.cont]=summary(independent.out[[1]])[[1]]["Mean"]
        #MSE.independent[sim,stage2.cont]=Bias.independent[sim,stage2.cont]^2 + (summary(independent.out[[1]])[[1]][,"SD"])^2
      }
      if (K1>1)
      {
        for (k in 1:K1)
        {
          post.sample=independent.out[[1]][,paste("p[",k,"]",sep="")]
          posterior[k]=sum(post.sample>q0)/length(post.sample)
        }
        ##### Bias:
        Bias.independent[sim,stage2.cont]=abs(summary(independent.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",x,"]",sep="")}),"Mean"]-p0[stage2.cont])
        pest.independent[sim,stage2.cont]=summary(independent.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",x,"]",sep="")}),"Mean"]
        #MSE.independent[sim,stage2.cont]=Bias.independent[sim,stage2.cont]^2 + (summary(independent.out[[1]])[[1]][sapply(1:K1,FUN=function(x){paste("p[",x,"]",sep="")}),"SD"])^2
      }
      decision.independent[sim,stage2.cont]=ifelse(posterior>Q.independent,1,0)############ Jags model for BHM:
    }
    ni.independent=sum(nik)
    
    print(sim)
    
    samplesize[scenario,1]=samplesize[scenario,1]+ni.sbhm
    samplesize[scenario,2]=samplesize[scenario,2]+ni.exnex
    samplesize[scenario,3]=samplesize[scenario,3]+ni.liu
    samplesize[scenario,4]=samplesize[scenario,4]+ni.bhm
    samplesize[scenario,5]=samplesize[scenario,5]+ni.independent
    
    ######## Summarize TP, TN, FP, FN:
    Decision[[scenario]][1,"% Perfect"]=Decision[[scenario]][1,"% Perfect"] + 1*(sum(c(decision.sbhm[sim,tp]==1,decision.sbhm[sim,tn]<=0))==K)
    Decision[[scenario]][1,"# TP"]=Decision[[scenario]][1,"# TP"] + sum(decision.sbhm[sim,tp]==1)
    Decision[[scenario]][1,"# TN"]=Decision[[scenario]][1,"# TN"] +  sum(decision.sbhm[sim,tn]<=0)
    Decision[[scenario]][1,"# FP"]=Decision[[scenario]][1,"# FP"] +  sum(decision.sbhm[sim,tn]==1)
    Decision[[scenario]][1,"# FN"]=Decision[[scenario]][1,"# FN"] +  sum(decision.sbhm[sim,tp]<=0)
    
    Decision[[scenario]][2,"% Perfect"]=Decision[[scenario]][2,"% Perfect"] + 1*(sum(c(decision.exnex[sim,tp]==1,decision.exnex[sim,tn]<=0))==K)
    Decision[[scenario]][2,"# TP"]=Decision[[scenario]][2,"# TP"] + sum(decision.exnex[sim,tp]==1)
    Decision[[scenario]][2,"# TN"]=Decision[[scenario]][2,"# TN"] +  sum(decision.exnex[sim,tn]<=0)
    Decision[[scenario]][2,"# FP"]=Decision[[scenario]][2,"# FP"] +  sum(decision.exnex[sim,tn]==1)
    Decision[[scenario]][2,"# FN"]=Decision[[scenario]][2,"# FN"] +  sum(decision.exnex[sim,tp]<=0)
    
    Decision[[scenario]][3,"% Perfect"]=Decision[[scenario]][3,"% Perfect"] + 1*(sum(c(decision.liu[sim,tp]==1,decision.liu[sim,tn]<=0))==K)
    Decision[[scenario]][3,"# TP"]=Decision[[scenario]][3,"# TP"] + sum(decision.liu[sim,tp]==1)
    Decision[[scenario]][3,"# TN"]=Decision[[scenario]][3,"# TN"] +  sum(decision.liu[sim,tn]<=0)
    Decision[[scenario]][3,"# FP"]=Decision[[scenario]][3,"# FP"] +  sum(decision.liu[sim,tn]==1)
    Decision[[scenario]][3,"# FN"]=Decision[[scenario]][3,"# FN"] +  sum(decision.liu[sim,tp]<=0)
    
    Decision[[scenario]][4,"% Perfect"]=Decision[[scenario]][4,"% Perfect"] + 1*(sum(c(decision.bhm[sim,tp]==1,decision.bhm[sim,tn]<=0))==K)
    Decision[[scenario]][4,"# TP"]=Decision[[scenario]][4,"# TP"] + sum(decision.bhm[sim,tp]==1)
    Decision[[scenario]][4,"# TN"]=Decision[[scenario]][4,"# TN"] +  sum(decision.bhm[sim,tn]<=0)
    Decision[[scenario]][4,"# FP"]=Decision[[scenario]][4,"# FP"] +  sum(decision.bhm[sim,tn]==1)
    Decision[[scenario]][4,"# FN"]=Decision[[scenario]][4,"# FN"] +  sum(decision.bhm[sim,tp]<=0)
    
    Decision[[scenario]][5,"% Perfect"]=Decision[[scenario]][5,"% Perfect"] + 1*(sum(c(decision.independent[sim,tp]==1,decision.independent[sim,tn]<=0))==K)
    Decision[[scenario]][5,"# TP"]=Decision[[scenario]][5,"# TP"] + sum(decision.independent[sim,tp]==1)
    Decision[[scenario]][5,"# TN"]=Decision[[scenario]][5,"# TN"] +  sum(decision.independent[sim,tn]<=0)
    Decision[[scenario]][5,"# FP"]=Decision[[scenario]][5,"# FP"] +  sum(decision.independent[sim,tn]==1)
    Decision[[scenario]][5,"# FN"]=Decision[[scenario]][5,"# FN"] +  sum(decision.independent[sim,tp]<=0)
    
    print(sim)
    decisions[[scenario]]=list(decision.sbhm,decision.exnex,decision.liu,decision.bhm,decision.independent)
  }
  #save(Decision,decisions,decision.sbhm,decision.exnex,decision.liu,decision.bhm,decision.independent,
  #     Bias.sbhm,Bias.exnex,Bias.liu,Bias.bhm,Bias.independent,
  #     pest.sbhm,pest.exnex,pest.liu,pest.bhm,pest.independent,samplesize,
  #     file=paste("~/Google Drive/Sanofi/Jin/R/simuresults/simuraw.RData",sep=""))
  
  Decision[[scenario]]=Decision[[scenario]]/num.sim
  
  OC[[scenario]]=matrix(NA,2*length(methods),K)
  colnames(OC[[scenario]])=sapply(1:ncol(OC[[scenario]]),FUN=function(x){paste("cancer ",x,sep="")})
  rownames(OC[[scenario]])=sapply(1:nrow(OC[[scenario]]),FUN=function(x){paste(methods[ceiling(x/2)]," - % ",ifelse(x/2!=floor(x/2),"reject","stop"),sep="")})
  OC[[scenario]][1,]=sapply(1:K,FUN=function(x){sum(decision.sbhm[,x]==1)/num.sim*100})
  OC[[scenario]][2,]=sapply(1:K,FUN=function(x){sum(decision.sbhm[,x]==-1)/num.sim*100})
  OC[[scenario]][3,]=sapply(1:K,FUN=function(x){sum(decision.exnex[,x]==1)/num.sim*100})
  OC[[scenario]][4,]=sapply(1:K,FUN=function(x){sum(decision.exnex[,x]==-1)/num.sim*100})
  OC[[scenario]][5,]=sapply(1:K,FUN=function(x){sum(decision.liu[,x]==1)/num.sim*100})
  OC[[scenario]][6,]=sapply(1:K,FUN=function(x){sum(decision.liu[,x]==-1)/num.sim*100})
  OC[[scenario]][7,]=sapply(1:K,FUN=function(x){sum(decision.bhm[,x]==1)/num.sim*100})
  OC[[scenario]][8,]=sapply(1:K,FUN=function(x){sum(decision.bhm[,x]==-1)/num.sim*100})
  OC[[scenario]][9,]=sapply(1:K,FUN=function(x){sum(decision.independent[,x]==1)/num.sim*100})
  OC[[scenario]][10,]=sapply(1:K,FUN=function(x){sum(decision.independent[,x]==-1)/num.sim*100})
  samplesize[scenario,]=samplesize[scenario,]/num.sim
  OC[[scenario]]=round(OC[[scenario]],1)
  Decision[[scenario]]=round(Decision[[scenario]],3)
  samplesize=round(samplesize,1)
  
  Bias[[scenario]]=matrix(NA,length(methods),K)
  colnames(Bias[[scenario]])=sapply(1:ncol(Bias[[scenario]]),FUN=function(x){paste("cancer ",x,sep="")})
  rownames(Bias[[scenario]])=methods
  MSE[[scenario]]=matrix(NA,length(methods),K)
  colnames(MSE[[scenario]])=sapply(1:ncol(MSE[[scenario]]),FUN=function(x){paste("cancer ",x,sep="")})
  rownames(MSE[[scenario]])=methods
  Bias[[scenario]][1,]=colMeans(Bias.sbhm)
  Bias[[scenario]][2,]=colMeans(Bias.exnex)# !!!!! 1001:2000
  Bias[[scenario]][3,]=colMeans(Bias.liu)
  Bias[[scenario]][4,]=colMeans(Bias.bhm)
  Bias[[scenario]][5,]=colMeans(Bias.independent)
  Bias[[scenario]]=round(Bias[[scenario]],4)
  MSE[[scenario]][1,]=Bias[[scenario]][1,]^2 + apply(pest.sbhm,2,var)
  MSE[[scenario]][2,]=Bias[[scenario]][2,]^2 + apply(pest.exnex,2,var)
  MSE[[scenario]][3,]=Bias[[scenario]][3,]^2 + apply(pest.liu,2,var)
  MSE[[scenario]][4,]=Bias[[scenario]][4,]^2 + apply(pest.bhm,2,var)
  MSE[[scenario]][5,]=Bias[[scenario]][5,]^2 + apply(pest.independent,2,var)
  MSE[[scenario]]=round(MSE[[scenario]],4)
}

