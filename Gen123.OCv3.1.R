
## V3.1: add the prob for go decision

library(rjags)
library(LaplacesDemon)
library(invgamma)
library(survival)  



Gen123.OC=function(n.rep,pT.true,pE.true,pT.true.c,pE.true.c,rho,beta0,betae,betat,gammad,aw,time.between,A1,A2,A3,ncohort1,cohortsize1,ncohort2,cohortsize2,n3,u00,u11,targetT,targetE,cutoff.eli.T,cutoff.eli.E,
                   targetHR,tauS,delta,cutoff.eli.TE){
  ## n.rep: number of replication
  ## pT.true: true toxicity rates for intervention doses 
  ## pE.true: true efficacy rates for intervention doses
  ## pT.true.c: true toxicity rate for control
  ## pE.true.c: true efficacy rate for control
  ## rho: latent multivariate-normal eff-tox correlation coefficient
  ## beta0: baseline log HR vector for pexp survival generating
  ## betae: efficacy log HR vector for pexp survival generating
  ## betat: toxicity log HR vector for pexp survival generating
  ## gammad: dose level log HR matrix for pexp survival generating, row for interval, column for dose level, last column for control (all 0) 
  ## aw: assess window for survival outcome
  ## time.between: the time between the phase II and III   
  ## A1: accrual rate per patient, 1 patient per A1 month in phase I  
  ## A2: accrual rate per patient, 1 patient per A2 month in phase II
  ## A3: accrual rate per patient, 1 patient per A3 month in phase III 
  ## ncohort1: number of cohorts in phase I
  ## cohortsize1: cohortsize in phase I
  ## ncohort2: number of cohorts in phase II
  ## cohortize2: cohortsize in phase II
  ## n3: sample size for phase III (include the number of patients treated at control and selected optimal dose in phase II) 
  ## u00: utility score for no toxicity, no efficacy
  ## u11: utility score for both toxicity and efficacy  
  ## targetT: highest acceptable toxicity rate for dose-finding
  ## targetE: lowest acceptable efficacy rate  for dose-finding
  ## cutoff.eli.T, cut-off probability for toxicity admissible
  ## cutoff.eli.E, cut-off probability for efficacy admissible 
  ## targetHR: highest acceptable HR for survial (compared to control)
  ## tauS: cut-off posterior probability for HR
  ## delta: cut-off posterior predictive probability for HR  
  ## cutoff.eli.TE: cut-off probability for utility, the ratio of the maximum value  
  
  
  ndose=length(pT.true) ## number of doses, without the control
  int=c(aw/2,aw/2) ## intervals for pexp survival generating
  u01=100 ## utility score for no toxicity, efficacy
  u10=0 ## utility score for toxicity, no efficacy
  utility=c(u11,u10,u01,u00)
  startdose=1
  n.earlystop=ncohort1*cohortsize1
  p.saf=0.6*targetT
  p.tox=1.4*targetT
  N1=6 
  N2=9
  phi=0.5 ## tunning parameter for RAR
  
  
  pTE=function(pT.true, pE.true, rho){
    f1<-function(x,bn.m1,bn.m2,rho){
      ff<-dnorm(x,bn.m1,1)*(1-pnorm(0,bn.m2+rho*(x-bn.m1),sqrt(1-rho^2)))
      return(ff)
    }
    ndose=length(pT.true)
    p10<-p01<-p00<-p11<-rep(0,ndose)
    for(d in 1:ndose){
      p11[d]<-integrate(f1,bn.m1=qnorm(pT.true[d]),bn.m2=qnorm(pE.true[d]),rho=rho,lower=0,upper=Inf)$value
      p10[d]<-pT.true[d]-p11[d]
      p01[d]<-pE.true[d]-p11[d]
      p00[d]<-1-p11[d]-p10[d]-p01[d]
    }
    return(list("pt1e1"=p11,"pt1e0"=p10,"pt0e1"=p01,"pt0e0"=p00))
  }
  
  uti=function(pT.true,pE.true,rho,u11,u00,targetT,targetE){
    u01=100
    u10=0
    utility=c(u11,u10,u01,u00)
    # Assume independence between toxicity and efficacy
    targetP<-c(targetE*targetT,targetT*(1-targetE),(1-targetT)*targetE,(1-targetT)*(1-targetE))
    
    # Calculate the benchmark utility
    uu = sum(targetP*utility) #highest unacceptable utility
    uu = uu+(100-uu)/2        # benchmark utility (i.e., desirable utility)
    p=pTE(pT.true, pE.true, rho)
    p11=p[[1]]
    p01=p[[3]]
    p00=p[[4]]
    u.true<-u11*p11+u00*p00+100*p01
    return(list( "benchmark"=uu, "true utility"=u.true   ))
    
  }
  
  get.boundary <- function(target, targetE, ncohort, cohortsize, p.saf=NA, p.tox=NA,  cutoff.eli,
                           cutoff.eli.E)
  {
    
    # if the user does not provide p.saf and p.tox, use the default values
    if(is.na(p.saf)) p.saf=0.6*target;
    if(is.na(p.tox)) p.tox=1.4*target;
    
    ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
    npts = ncohort*cohortsize;
    ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;elimE=NULL;
    for(n in (1:ncohort)*cohortsize)
    {
      error.min=3;
      for(m1 in 0:(n-1))
      {
        for(m2 in (m1+1):n)
        {
          
          error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
          error2 = 1-pbinom(m1, n, p.saf);
          error3 = pbinom(m2-1, n, p.tox);
          
          error=error1+error2+error3;
          if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
        }
      }
      ntrt = c(ntrt, n);
      b.e = c(b.e, cutoff1);
      b.d = c(b.d, cutoff2);
      
      elimineed=0; # indicating whether elimination is needed
      elimineedE=0
      if(n<3) { elim = c(elim, NA); elimE = c(elimE,NA)}  # require treating at least 3 patients before eliminating a dose
      else
      {
        for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
        {
          if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
        }
        if(elimineed==1) { elim = c(elim, ntox); }
        else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
        
        for(neff in n:0){
          if(pbeta(targetE,neff+1,n-neff+1)>cutoff.eli.E){elimineedE=1; break;}
        }
        if(elimineedE==1){elimE=c(elimE,neff)} else {elimE=c(elimE,NA)}
      }
    }
    for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
    boundaries = rbind(ntrt, elim, b.d, b.e,elimE);
    rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=", 
                             "Deescalate if # of DLT >=",  "Escalate if # of DLT <=", "Eliminate if # of Eff <=");
    colnames(boundaries) = rep("", ncohort);
    
    return(boundaries);
  }
  
  pi=pTE(pT.true, pE.true, rho)
  
  p11=pi[[1]]
  p10=pi[[2]]
  p01=pi[[3]]
  p00=pi[[4]]
  A1C=A1*cohortsize1 ## 1 cohort per A1C month in phase I
  u.t=uti(pT.true,pE.true,rho,u11,u00,targetT,targetE)
  uu=u.t[[1]]
  u.true=u.t[[2]]
  
  boundary.temp=get.boundary(targetT, targetE, ncohort1, cohortsize1,cutoff.eli=cutoff.eli.T, cutoff.eli.E = cutoff.eli.E) 
  b.e=boundary.temp[4,]   # escalation boundary
  b.d=boundary.temp[3,]   # deescalation boundary
  b.elim=boundary.temp[2,]  # elimination boundary
  b.elimE=boundary.temp[5,]
  
  phaseI=function(p11,p10,p01,p00){
    yT<-yE<-rep(0, ndose);    ## number of DLT/efficacy at each dose level
    y01<-y10<-y11<-y00<-rep(0,ndose); ## number of different outcomes at each dose level
    n<-rep(0, ndose);         ## number of patients treated at each dose level
    earlystop=0;              ## indicate whether the trial terminates early
    d=startdose;              ## starting dose level
    elimi = rep(0, ndose);    ## indicate whether doses are eliminated due to toxicity
    elimiE=  rep(0,ndose);    ## indicate whether doses are eliminated due to efficacy
    safe = 0
    posH<-rep(1-uu/100,ndose)
    duration=0
    for (i in 1: ncohort1){
      datate=rmultinom(1,cohortsize1,c(p11[d],p10[d],p01[d],p00[d])  )
      y11[d]=y11[d]+datate[1]
      y10[d]=y10[d]+datate[2]
      y01[d]=y01[d]+datate[3]
      y00[d]=y00[d]+datate[4]
      duration=duration+A1C
      n[d]=n[d]+cohortsize1
      yT<-y10+y11
      yE<-y01+y11
      nc<-n[d]/cohortsize1
      # determine whether current dose level is overly toxic
      if(!is.na(b.elim[nc]))
      {
        if(yT[d]>=b.elim[nc]) 
        {      
          elimi[d:ndose]=1;
          if(d==1) {earlystop=1; break;} 
        }
      }
      
      if(!is.na(b.elimE[nc]))
      {
        if(yE[d]<=b.elimE[nc]) 
        {      
          elimi[d]=1;
        }
      }
      
      if(sum(elimi==1)==ndose) {earlystop=1; break;} 
      
      u_curr<-(u01*y01[d]+u10*y10[d]+u11*y11[d]+u00*y00[d])/100
      
      posH[d] = 1-pbeta(uu/100,1+u_curr,n[d]-u_curr+1)
      
      posH <- posH*(1-elimi);
      if(n[d]>=N1){safe=1} else{safe=0}
      if(n[d]>=n.earlystop){break}
      
      if (yT[d]>=b.d[nc] && d!=1) {
        if(sum(elimi[1:(d-1)]==0)>0){d_opt=max(which(elimi[1:(d-1)]==0))} else {
          if(elimi[d]==1){earlystop=1;break} else{d_opt=d}}
      } else if (yT[d]>=b.d[nc] && d==1) {if(elimi[d]==0){d_opt=d} else{earlystop=1;break}
      } else{
        admi_set=d;
        if(d>1){
          if(sum(elimi[1:(d-1)]==0)>0){admi_set<-c(admi_set,max(which(elimi[1:(d-1)]==0)))} 
        }
        if(d<ndose){
          if(safe==0){
            if(sum(elimi[(d+1):ndose]==0)>0){admi_set<-c(admi_set,d+min(which(elimi[(d+1):ndose]==0)))}
          } else {
            if(yT[d]<=b.e[nc] & sum(elimi[(d+1):ndose]==0)>0){admi_set<-c(admi_set,d+min(which(elimi[(d+1):ndose]==0)))}
          }
        }
        #if(length(admi_set)>1 & eff_cut>0.5 & n[d]>=9){admi_set<-admi_set[-1]}
        temp.posH<-posH[admi_set]+runif(length(admi_set))*(10^-15)
        d_opt=admi_set[which.max(temp.posH)]
      }	
      
      if (elimi[d_opt]==1) {earlystop=1; break} 
      if (sum(elimi)==ndose) {earlystop=1; break}
      
      if (d<ndose){
        if(sum(elimi[(d+1):ndose]==0)>0){
          d_temp=d+min(which(elimi[(d+1):ndose]==0))
          if(n[d]>=N2 & n[min(d_temp,ndose)]==0 & yT[d]<b.d[n[d]/cohortsize1]){d_opt<-d_temp} 
        }  
        
      }
      d<-d_opt
    }
    return(list( "n"=n, "y"=rbind(y11,y10,y01,y00), "duration"=duration, "earlystop"=earlystop     )) 
    
  }
  A2C=A2*cohortsize2
  pi.c=pTE(pT.true.c, pE.true.c, rho)
  p11.c=pi.c[[1]]
  p10.c=pi.c[[2]]
  p01.c=pi.c[[3]]
  p00.c=pi.c[[4]]
  p11.all=c(p11,p11.c)
  p10.all=c(p10,p10.c)
  p01.all=c(p01,p01.c)
  p00.all=c(p00,p00.c)
  
  
  adm=function(y,targetT,targetE,cutoff.eli.T,cutoff.eli.E){
    y=as.matrix(y)
    n=apply(y,2,sum)
    n.A=ifelse(n>0,1,0)
    yT=y[1,]+y[2,]
    yE=y[1,]+y[3,]
    postT=1-pbeta(targetT,1+yT,1+n-yT)
    T.A=ifelse(postT<cutoff.eli.T,1,0)
    postE=pbeta(targetE,1+yE,1+n-yE)
    E.A=ifelse(postE<cutoff.eli.E,1,0)
    re=n.A*T.A*E.A
    re1=which(re==1)
    return(re1)
  }
  
  can=function(y,utility,cutoff.eli.TE){
    y=as.matrix(y)
    J=dim(y)[2]
    u.e=rep(0,J)
    for (j in 1:J){
      u.e[j]=sum(utility*y[,j])/100
    }
    n=apply(y,2,sum)
    temp=u.e/n
    mtemp=max(temp)
    re1=which( temp>=cutoff.eli.TE*mtemp   )
    re2=u.e[re1]/n[re1]
    return(list(re1,re2))
  }
  
  
  rpwexp <- function(n, rate, intervals=NULL, cumulative=FALSE){
    if(is.null(intervals)){
      if (cumulative){return(cumsum(rexp(n,rate[1])))}else
        return(rexp(n,rate[1]))}
    k <- length(rate)
    if (k==1){
      if(cumulative){return(cumsum(rexp(n,rate)))}else
        return(rexp(n,rate))
    }
    if (length(intervals) < k-1) stop("length(intervals) must be at least length(rate) - 1")
    tx <- 0
    j <- 1
    times <- array(0,n)
    timex <- cumsum(intervals)
    indx <- array(TRUE,n)
    for(i in 1:k){
      nindx <- sum(indx)
      if (nindx==0) break
      increment <- rexp(nindx,rate[i])
      if (cumulative) times[indx] <- tx + cumsum(increment)
      else times[indx] <- tx + increment
      if (i<k){
        tx <- timex[i]
        indx <- (times > timex[i])
      }
    }
    return(times)
  }
  
  
  
  pfs.gen=function(beta0,int,betae,betat,gammad,dose,ye,yt){
    ## gammad: row for interval; column for dose level
    
    ha=c(exp(beta0+betae*ye+betat*yt+gammad[,dose]) ,1)
    return(rpwexp(1,ha,int))
  }
  
  
  
  phaseII=function(data1){
    n.12=data1[[1]] 
    n.12=c(n.12,0) ## patient allocation combine phase I and II, the last cell is for control arm 
    y.12=data1[[2]]
    y.12=cbind(y.12, rep(0,4)) ## toxicity-efficacy distribution combine phase I and II, the last column is for control arm  
    dur.12=data1[[3]]
    earlystop=data1[[4]] 
    
    yT.2=NULL ## individual toxicity outcome in phase II
    yE.2=NULL ## individual efficacy outcome in phase II
    yS.2=NULL ## individual survival outcome in phase II
    enter=NULL ## individual enter time in phase II
    
    dose.2=NULL ## individual dose allocation in phase II
    
    if(earlystop==1){
      return(list( "n.12"=n.12, "y.12"=y.12, "dur.12"=dur.12, "earlystop"=earlystop, "yT.2"=NULL, "yE.2"=NULL,"yS.2"=NULL, "yS.2.obs"=NULL, "status.2"=NULL, "dose.2"=NULL, "can.set"=NULL, "follow"=NULL     )    )
      
    }
    else{
      for (i in 1:ncohort2){
        adset=adm(y.12[,1:ndose],targetT,targetE,cutoff.eli.T,cutoff.eli.E)
        if(length(adset)==0){can.set=NULL; earlystop=1; break}
        else{
          candidate=can(as.matrix(y.12[,adset]),utility,cutoff.eli.TE)
          can.set=adset[ candidate[[1]] ]
          can.uti=candidate[[2]]
          m=length(can.set)
          ran.ratio=c(  (m*(can.uti)^phi)/( (m+1)*sum((can.uti)^phi)  ), 1/(m+1)      )
          dose.temp=sample(c(can.set,(ndose+1)),cohortsize2,replace=T,prob=ran.ratio)
          dose.2=c(dose.2,dose.temp)
          enter.temp=rep(dur.12,cohortsize2)
          enter=c(enter,enter.temp)
          y.temp=matrix(0,nrow=4,ncol=cohortsize2)
          for (j in 1: cohortsize2){
            y.temp[,j]=rmultinom(1,1,c(p11.all[dose.temp[j]],p10.all[dose.temp[j]],p01.all[dose.temp[j]],p00.all[dose.temp[j]])  )
            y.12[,dose.temp[j]]=y.12[,dose.temp[j]]+y.temp[,j]
          }
          n.12=apply(y.12,2,sum)
          dur.12=dur.12+A2C
          yT.2.temp=y.temp[1,]+y.temp[2,]
          yT.2=c(yT.2,yT.2.temp)
          yE.2.temp=y.temp[1,]+y.temp[3,]
          yE.2=c(yE.2,yE.2.temp)
          yS.2.temp=NULL
          for (j in 1:cohortsize2){
            yS.2.temp[j]=pfs.gen(beta0,int,betae,betat,gammad,dose.temp[j],yE.2.temp[j],yT.2.temp[j])
          }
          yS.2=c(yS.2,yS.2.temp) 
        }
        
        
      }
      dur.12=dur.12+time.between-A2C
      follow=dur.12-enter
      follow[follow>aw]=aw ## actual follow-up time, cannot greater than assessment window
      status.2=rep(1,length(yS.2)) ## 1 for no censoring, 0 for censoring
      yS.2.obs=yS.2
      for(j in 1: length(yS.2)){
        if(yS.2.obs[j]>follow[j]){
          yS.2.obs[j]=follow[j]
          status.2[j]=0
        }
      }
      return(list( "n.12"=n.12, "y.12"=y.12, "dur.12"=dur.12, "earlystop"=earlystop, "yT.2"=yT.2, "yE.2"=yE.2,"yS.2"=yS.2, "yS.2.obs"=yS.2.obs, "status.2"=status.2, "dose.2"=dose.2, "can.set"=can.set, "follow"=follow   )    )
    }
  }
  
  
  
  model<-"model{
  
  
  for (i in 1:n.obs){
  t.obs[i]~dweib(alpha,lambda.obs[i])
  lambda.obs[i]<-(1/b)^alpha*exp(beta1*eff.obs[i]+beta2*tox.obs[i]+beta3[dose.obs[i]])    
  }
  
  for (i in 1:n.cen){
  status.cen[i]~dbern(S[i])
  S[i]<-pweib(t.cen[i],alpha,lambda.cen[i])
  lambda.cen[i]<-(1/b)^alpha*exp(beta1*eff.cen[i]+beta2*tox.cen[i]+beta3[dose.cen[i]]) 
  }
  
  beta1~dnorm(0.0,0.01)
  beta2~dnorm(0.0,0.01)
  alpha~dgamma(0.01,0.01)
  b~dgamma(0.01,0.01)
  for(i in 1: ndose){
  beta3[i]~dnorm(0.0,0.01)
  }
  beta3[ndose+1]~dnorm(0.0,100)
}"
  
  ## get mode
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  
  
  
  ## find on the optimal dose based on hazard ratio
  opt=function(data2){
    y.12=data2[[2]]
    earlystop=data2[[4]]
    yT.2=data2[[5]]
    yE.2=data2[[6]]
    yS.2.obs=data2[[8]]
    status.2=data2[[9]]
    dose.2=data2[[10]]
    can.set=data2[[11]]
    if (earlystop==1){
      opt.dose=0
      return(opt.dose) ## 0 is control arm
    }else {
      data.final=cbind(yS.2.obs,status.2,dose.2,yT.2,yE.2)
      data.obs<-subset(data.final,status.2==1)
      data.cen<-subset(data.final,status.2==0)
      n.obs<-length(data.obs[,1])
      n.cen<-length(data.cen[,1])
      data.mcmc<-list(n.obs=n.obs,n.cen=n.cen,t.obs=data.obs[,1],t.cen=data.cen[,1],
                      dose.obs=data.obs[,3], dose.cen=data.cen[,3],tox.obs=data.obs[,4], tox.cen=data.cen[,4],
                      eff.obs=data.obs[,5], eff.cen=data.cen[,5],status.cen=data.cen[,2],ndose=ndose)
      mcmc.final<- jags.model(textConnection(model), data=data.mcmc,  quiet=T)
      update(mcmc.final, 1000, progress.bar="none")
      post<- jags.samples(mcmc.final, n.iter=1000,variable.names=c("alpha","b","beta1","beta2","beta3"), progress.bar="none")
      post.alpha.sample=post$alpha[1:    1000   ]
      post.b.sample=post$b[1:    1000   ]
      post.beta1.sample=post$beta1[1:    1000   ]
      post.beta2.sample=post$beta2[1:    1000   ]
      post.beta3.sample=matrix(post$beta3,ncol=ndose+1,byrow=T)
      post.p11.sample=matrix(0, nrow=1000,ncol=ndose+1)
      post.p10.sample=matrix(0, nrow=1000,ncol=ndose+1)
      post.p01.sample=matrix(0, nrow=1000,ncol=ndose+1)
      post.p00.sample=matrix(0, nrow=1000,ncol=ndose+1)
      for (j in 1: (ndose+1)){
        p.temp=rdirichlet(1000, y.12[,j]+0.25)
        post.p11.sample[,j]=p.temp[,1]
        post.p10.sample[,j]=p.temp[,2]
        post.p01.sample[,j]=p.temp[,3]
        post.p00.sample[,j]=p.temp[,4]
      }
      post.s.sample=matrix(0, nrow=1000,ncol=ndose+1) ## calculate survival rate, row for MCMC samples, column for dose levels
      for (i in 1:1000){
        for (j in 1: (ndose+1)){
          post.s.sample[i,j]=post.p11.sample[i,j]*exp( -(aw/post.b.sample[i])^post.alpha.sample[i]*exp( post.beta1.sample[i]+post.beta2.sample[i]+post.beta3.sample[i,j]    )        )+
            post.p10.sample[i,j]*exp( -(aw/post.b.sample[i])^post.alpha.sample[i]*exp(post.beta2.sample[i]+post.beta3.sample[i,j]    )        )+
            post.p01.sample[i,j]*exp( -(aw/post.b.sample[i])^post.alpha.sample[i]*exp( post.beta1.sample[i]+post.beta3.sample[i,j]    )        )+
            post.p00.sample[i,j]*exp( -(aw/post.b.sample[i])^post.alpha.sample[i]*exp( post.beta3.sample[i,j]    )        )
        }
      }
      post.hr.sample=matrix(0, nrow=1000,ncol=length(can.set)) ## calculate hazard ratio compared with the control, row for MCMC samples, column for dose levels in candidate set
      for (i in 1:1000){
        for (j in 1: length(can.set)){
          post.hr.sample[i,j]=log( post.s.sample[i,can.set[j]]   )/log( post.s.sample[i,ndose+1]   )
        }
      }
      hr.min=NULL
      for (i in 1:1000){
        hr.min[i]=which( post.hr.sample[i,]==min(post.hr.sample[i,]     )      )[1]
      }
      opt.dose=can.set[getmode(hr.min)]
      return(opt.dose)
    }
    
    
  }
  
  
  ## no/no-go decision based on the predictive probability
  ## dose.sel: the selected optimal dose
  ## n3: the maximum sample size for phase 3, combined with the patients treated at optimal dose and control in phase II
  ## rinvgamma(n, shape, rate = 1, scale = 1/rate) for inverse gamma, shape=alpha, rate=beta
  
  pre.decision=function(data2,dose.sel,n3,targetHR,tauS,delta){
    yS.2.obs=data2[[8]]
    status.2=data2[[9]] ## 1 for no censoring, 0 for censoring
    dose.2=data2[[10]]
    if (dose.sel==0){
      decision=0
      return(decision)
    } else{
      yS.i=yS.2.obs[dose.2==dose.sel]
      yS.c=yS.2.obs[dose.2==(ndose+1)]
      status.i=status.2[dose.2==dose.sel]
      status.c=status.2[dose.2==(ndose+1)]
      alpha.i=sum(status.i)+1
      beta.i=sum(yS.i)+1
      alpha.c=sum(status.c)+1
      beta.c=sum(yS.c)+1
      alpha.i.cob=rep(alpha.i, 1000)
      beta.i.cob=rep(beta.i, 1000)
      alpha.c.cob=rep(alpha.c, 1000)
      beta.c.cob=rep(beta.c, 1000)
      for (i in 1:1000){
        mu.i=invgamma::rinvgamma(n3/2-length(yS.i), alpha.i, beta.i)
        pre.S.i=rexp( n3/2-length(yS.i), 1/mu.i  )
        pre.status.i=as.numeric(pre.S.i<aw)
        pre.S.i[pre.S.i>aw]=aw
        alpha.i.cob[i]=alpha.i.cob[i]+sum(pre.status.i)
        beta.i.cob[i]=beta.i.cob[i]+sum(pre.S.i)
        mu.c=rinvgamma(n3/2-length(yS.c), alpha.c, beta.c)
        pre.S.c=rexp( n3/2-length(yS.c), 1/mu.c  )
        pre.status.c=as.numeric(pre.S.c<aw)
        pre.S.c[pre.S.c>aw]=aw
        alpha.c.cob[i]=alpha.c.cob[i]+sum(pre.status.c)
        beta.c.cob[i]=beta.c.cob[i]+sum(pre.S.c)
        
      }
      
      temp1=NULL
      temp2=NULL
      for (i in 1:1000){
        mu.i.sample=rinvgamma(1000, alpha.i.cob[i],beta.i.cob[i])
        mu.c.sample=rinvgamma(1000, alpha.c.cob[i],beta.c.cob[i])
        hr.sample=mu.c.sample/mu.i.sample
        temp1[i]=mean(hr.sample<targetHR)
        temp2[i]=ifelse(temp1[i]>tauS,1,0)
      }
      if( mean(temp2)>=delta   ){
        decision=1
      } else{decision=0}
      return(decision)
      
    }
    
  }
  
  
  
  phaseIII=function(data2,dose.sel,decision){
    ##   return(list( "n.12"=n.12, "y.12"=y.12, "dur.12"=dur.12, "earlystop"=earlystop, "yT.2"=yT.2, "yE.2"=yE.2,"yS.2"=yS.2, "yS.2.obs"=yS.2.obs, "status.2"=status.2, "dose.2"=dose.2, "can.set"=can.set, "follow"=follow   )    )
    ##                 1            2            3                4                      5            6           7           8                     9                    10               11                 12           
    n.123=data2[[1]]
    y.123=data2[[2]]
    dur.123=data2[[3]]
    yS.23=data2[[7]]
    dose.23=data2[[10]]
    if (decision==0){
      final.decision=0 ## final decision; 0 for the control arm
      ##yS.23.obs=data2[[8]]
      ##status.23=data2[[9]] ## 1 for no censoring, 0 for censoring 
      ##fit=survfit(Surv(yS.23.obs, status.23)~1)
      ##survmedian.23=as.numeric(surv_median(fit))   ## median survival time combined phase 2 and 3, and all the doses and controls
      ss.123=sum(n.123) ## total sample size
      yT.123=y.123[1,]+y.123[2,]
      yE.123=y.123[1,]+y.123[3,]
      pT.123=sum(yT.123)/ss.123
      pE.123=sum(yE.123)/ss.123
      n.123=c( n.123[(ndose+1)], n.123[1:ndose]   )
      ##return(list("final decision"=final.decision, "patient allocation"=n.123, "median survival"=survmedian.23, "toxicity rate"=pT.123, "efficacy rate"=pE.123, "trial duration"=dur.123, "total sample size"=ss.123     ))
      return(list("optimal dose"=dose.sel,"final decision"=final.decision, "patient allocation"=n.123, "toxicity rate"=pT.123, "efficacy rate"=pE.123, "trial duration"=dur.123, "total sample size"=ss.123     ))
    } else{ ## interim analysis
      n2.i=length(dose.23[dose.23==dose.sel])
      n2.c=length(dose.23[dose.23==(ndose+1)])
      n2.i.stage1=n3/4-n2.i  
      n2.c.stage1=n3/4-n2.c 
      dur.123=dur.123+(n2.i.stage1+n2.c.stage1)*A3+aw
      dose.temp=rep( c(dose.sel, (ndose+1) ),c(  n2.i.stage1, n2.c.stage1   )   )
      dose.23=c(dose.23, dose.temp)
      y.temp=matrix(0,nrow=4,ncol=(n2.i.stage1+n2.c.stage1))
      for (j in 1: length(dose.temp) ){
        y.temp[,j]=rmultinom(1,1,c(p11.all[dose.temp[j]],p10.all[dose.temp[j]],p01.all[dose.temp[j]],p00.all[dose.temp[j]])  )
        y.123[,dose.temp[j]]=y.123[,dose.temp[j]]+y.temp[,j]
      }
      n.123=apply(y.123,2,sum)
      ss.123=sum(n.123) ## total sample size
      yT.123=y.123[1,]+y.123[2,]
      yE.123=y.123[1,]+y.123[3,]
      pT.123=sum(yT.123)/ss.123
      pE.123=sum(yE.123)/ss.123
      yT.temp=y.temp[1,]+y.temp[2,]
      yE.temp=y.temp[1,]+y.temp[3,]
      yS.temp=NULL
      for (j in 1:length(dose.temp)){
        yS.temp[j]=pfs.gen(beta0,int,betae,betat,gammad,dose.temp[j],yE.temp[j],yT.temp[j])
      }
      yS.23=c(yS.23, yS.temp)
      status.23=rep(1,length(yS.23)) ## 1 for no censoring, 0 for censoring
      yS.23.obs=yS.23
      for(j in 1: length(yS.23)){
        if(yS.23.obs[j]>aw){
          yS.23.obs[j]=aw
          status.23[j]=0
        }
      }
      yS.23.obs.i=yS.23.obs[dose.23==dose.sel]
      yS.23.obs.c=yS.23.obs[dose.23==(ndose+1)]
      status.23.i=status.23[dose.23==dose.sel]
      status.23.c=status.23[dose.23==(ndose+1)]
      
      yS.23.obs.comb=c(yS.23.obs.i,yS.23.obs.c)
      status.23.comb=c(status.23.i,status.23.c)
      group.23.comb=rep( c(1,0),c( n3/4, n3/4  )    )
      logrank.chisq=survdiff(Surv(yS.23.obs.comb, status.23.comb) ~ group.23.comb)$chisq
      
      if ( logrank.chisq<=0.54^2  ){  ## early stop and select the control
        final.decision=0
        n.123=c( n.123[(ndose+1)], n.123[1:ndose]   )
        return(list("optimal dose"=dose.sel,"final decision"=final.decision, "patient allocation"=n.123, "toxicity rate"=pT.123, "efficacy rate"=pE.123, "trial duration"=dur.123, "total sample size"=ss.123     ))
      } else if ( logrank.chisq>2.96^2  ) { ## early stop and select the intervention
        final.decision=dose.sel
        n.123=c( n.123[(ndose+1)], n.123[1:ndose]   )
        return(list("optimal dose"=dose.sel,"final decision"=final.decision, "patient allocation"=n.123, "toxicity rate"=pT.123, "efficacy rate"=pE.123, "trial duration"=dur.123, "total sample size"=ss.123     ))
      } else { ## enroll the second half 
        n2.i.stage2=n3/4  
        n2.c.stage2=n3/4 
        dur.123=dur.123+(n2.i.stage2+n2.c.stage2)*A3
        dose.temp2=rep( c(dose.sel, (ndose+1) ),c(  n2.i.stage2, n2.c.stage2   )   )
        dose.23=c(dose.23, dose.temp2)
        y.temp2=matrix(0,nrow=4,ncol=(n2.i.stage2+n2.c.stage2))
        for (j in 1: length(dose.temp2) ){
          y.temp2[,j]=rmultinom(1,1,c(p11.all[dose.temp2[j]],p10.all[dose.temp2[j]],p01.all[dose.temp2[j]],p00.all[dose.temp2[j]])  )
          y.123[,dose.temp2[j]]=y.123[,dose.temp2[j]]+y.temp2[,j]
        }
        n.123=apply(y.123,2,sum)
        ss.123=sum(n.123) ## total sample size
        yT.123=y.123[1,]+y.123[2,]
        yE.123=y.123[1,]+y.123[3,]
        pT.123=sum(yT.123)/ss.123
        pE.123=sum(yE.123)/ss.123
        yT.temp2=y.temp2[1,]+y.temp2[2,]
        yE.temp2=y.temp2[1,]+y.temp2[3,]
        yS.temp2=NULL
        for (j in 1:length(dose.temp2)){
          yS.temp2[j]=pfs.gen(beta0,int,betae,betat,gammad,dose.temp2[j],yE.temp2[j],yT.temp2[j])
        }
        yS.23=c(yS.23, yS.temp2)
        status.23=rep(1,length(yS.23)) ## 1 for no censoring, 0 for censoring
        yS.23.obs=yS.23
        for(j in 1: length(yS.23)){
          if(yS.23.obs[j]>aw){
            yS.23.obs[j]=aw
            status.23[j]=0
          }
        }
        
        yS.23.obs.i=yS.23.obs[dose.23==dose.sel]
        yS.23.obs.c=yS.23.obs[dose.23==(ndose+1)]
        status.23.i=status.23[dose.23==dose.sel]
        status.23.c=status.23[dose.23==(ndose+1)]
        
        yS.23.obs.comb=c(yS.23.obs.i,yS.23.obs.c)
        status.23.comb=c(status.23.i,status.23.c)
        group.23.comb=rep( c(1,0), c( n3/2, n3/2  )    )
        logrank.chisq=survdiff(Surv(yS.23.obs.comb, status.23.comb) ~ group.23.comb)$chisq
        if (logrank.chisq>1.94^2){
          final.decision=dose.sel
        } else { final.decision=0 }
        n.123=c( n.123[(ndose+1)], n.123[1:ndose]   )
        return(list("optimal dose"=dose.sel,"final decision"=final.decision, "patient allocation"=n.123, "toxicity rate"=pT.123, "efficacy rate"=pE.123, "trial duration"=dur.123, "total sample size"=ss.123     ))
      }
      
      
    }
    
  }
  
  upfs=function(p11.all,p10.all,p01.all,p00.all,beta0,int,betae,betat,gammad){
    re=NULL
    for (j in 1: (ndose+1)){
      re[j]=p11.all[j]*exp(-sum(exp(beta0+betat+betae+gammad[,j])*int))+p10.all[j]*exp(-sum(exp(beta0+betat+gammad[,j])*int))+p01.all[j]*exp(-sum(exp(beta0+betae+gammad[,j])*int))+p00.all[j]*exp(-sum(exp(beta0+gammad[,j])*int)) 
    }
    return(re)
  }
  
  pT.true.all=c(pT.true,pT.true.c)
  pE.true.all=c(pE.true,pE.true.c)
  u.true.all=c(u.true,u11*p11.c+u00*p00.c+100*p01.c)
  pS.true.all=upfs(p11.all,p10.all,p01.all,p00.all,beta0,int,betae,betat,gammad)
  
  dose.sel.m=rep(0,n.rep)
  decision.m=rep(0,n.rep)
  allocation.m=matrix(0, nrow=n.rep, ncol=(ndose+1))
  Trate.m=rep(0,n.rep)
  Erate.m=rep(0,n.rep)
  duration.m=rep(0,n.rep)
  ss.m=rep(0,n.rep)
  go.m=rep(0,n.rep)
  
  for (i in 1: n.rep){
    data1=phaseI(p11,p10,p01,p00)
    data2=phaseII(data1)
    dose.sel=opt(data2)
    decision=pre.decision(data2,dose.sel,n3,targetHR,tauS,delta)
    data3=phaseIII(data2,dose.sel,decision)
    if(dose.sel==0){
      go.m[i]=99
    } else{ go.m[i]=decision    }
    dose.sel.m[i]=data3[[1]]
    decision.m[i]=data3[[2]]
    allocation.m[i,]=data3[[3]]
    Trate.m[i]=data3[[4]]
    Erate.m[i]=data3[[5]]
    duration.m[i]=data3[[6]]
    ss.m[i]=data3[[7]]
  }
  return( list("true toxicity"=c(pT.true.all[ndose+1],pT.true.all[1:ndose]), "true efficacy"=c(pE.true.all[ndose+1],pE.true.all[1:ndose]),
               "true utility"=c(u.true.all[ndose+1],u.true.all[1:ndose]), "true suvival"=c(pS.true.all[ndose+1],pS.true.all[1:ndose]),optimal.dose=table(dose.sel.m)/n.rep,
               "final decision"=table(decision.m)/n.rep, "patient allocation"=apply( allocation.m,2,mean   ), "toxicity rate"=mean(Trate.m), "efficacy rate"=mean(Erate.m), "trial duration"=mean(duration.m), "total sample size"=mean(ss.m),
               "Prob of go decision"=  length(go.m[go.m==1])/length(go.m[go.m!=99])       )    )
  
}




# scenario 1.
 n.rep=5000
 pT.true=c(0.05,0.1,0.2,0.25,0.4)
 pE.true=c(0.1,0.35,0.4,0.5,0.5)
 pT.true.c=0.2
 pE.true.c=0.3
 rho=0.1
 beta0=c(1.1*log(-log(0.362)/6), 0.9*log(-log(0.362)/6))
 betae=rep(  -log(2),2)
 betat=rep(log(1.5),2)
 gammad=matrix(c(log( log(0.3)/log(0.373168)),log( log(0.3)/log(0.373168)),log( log(0.346)/log(0.4264148)),log( log(0.346)/log(0.4264148)),log( log(0.4)/log(0.4244721)),log( log(0.4)/log(0.4244721)),
                 log( log(0.398)/log(0.4417212)),log( log(0.398)/log(0.4417212)),log( log(0.4)/log(0.4207907)),log( log(0.4)/log(0.4207907)),0,0), nrow=2  )
 aw=6
 time.between=1
 A1=1/3
 A2=1/5
 A3=1/10
 ncohort1=10
 cohortsize1=3
 ncohort2=10
 cohortsize2=5
 n3=500
 u00=40
 u11=60
 targetT=0.35
 targetE=0.2
 cutoff.eli.T=0.95
 cutoff.eli.E=0.9
 targetHR=0.85
 tauS=0.8
 cutoff.eli.TE=0.5
 delta=0.5

 Gen123.OC(n.rep,pT.true,pE.true,pT.true.c,pE.true.c,rho,beta0,betae,betat,gammad,aw,time.between,A1,A2,A3,ncohort1,cohortsize1,ncohort2,cohortsize2,n3,u00,u11,targetT,targetE,cutoff.eli.T,cutoff.eli.E,
           targetHR,tauS,delta,cutoff.eli.TE)
#
#
#

# 
# ## scenario 2.
 n.rep=5000
 pT.true=c(0.02,0.05,0.1,0.15,0.2)
 pE.true=c(0.1,0.2,0.3,0.4,0.5)
 pT.true.c=0.1
 pE.true.c=0.3
 rho=0.1
 beta0=c(1.1*log(-log(0.242)/6), 0.9*log(-log(0.242)/6))   
 betae=rep(  -log(2),2)
 betat=rep(log(1.5),2)
 gammad=matrix(c(log( log(0.122)/log(0.3)),log( log(0.122)/log(0.3)),log( log(0.194)/log(0.28)),log( log(0.194)/log(0.28)),log( log(0.406)/log(0.3)),log( log(0.406)/log(0.3)),
                 log( log(0.598)/log(0.3)),log( log(0.598)/log(0.3)),log( log(0.262)/log(0.3)),log( log(0.262)/log(0.3)),0,0), nrow=2  )
 aw=6
 time.between=1
 A1=1/3
 A2=1/5
 A3=1/10
 ncohort1=10
 cohortsize1=3
 ncohort2=10
 cohortsize2=5
 n3=500
 u00=40
 u11=60
 targetT=0.35
 targetE=0.2
 cutoff.eli.T=0.95
 cutoff.eli.E=0.9
 targetHR=0.85
 tauS=0.8
 cutoff.eli.TE=0.5

 delta=0.5
 
 Gen123.OC(n.rep,pT.true,pE.true,pT.true.c,pE.true.c,rho,beta0,betae,betat,gammad,aw,time.between,A1,A2,A3,ncohort1,cohortsize1,ncohort2,cohortsize2,n3,u00,u11,targetT,targetE,cutoff.eli.T,cutoff.eli.E,
           targetHR,tauS,delta,cutoff.eli.TE)
 





# ## scenario 3.
 n.rep=5000
 pT.true=c(0.05,0.1,0.15,0.2,0.3)
 pE.true=c(0.15,0.3,0.4,0.35,0.3)
 pT.true.c=0.1
 pE.true.c=0.3
 rho=0.1
 beta0=c(1.1*log(-log(0.242)/6), 0.9*log(-log(0.242)/6))
 betae=rep(  -log(2),2)
 betat=rep(log(1.5),2)
 gammad=matrix(c(log( log(0.225)/log(0.3)),log( log(0.225)/log(0.3)),log( log(0.28)/log(0.28)),log( log(0.28)/log(0.28)),log( log(0.387)/log(0.3)),log( log(0.387)/log(0.3)),
                 log( log(0.714)/log(0.3)),log( log(0.714)/log(0.3)),log( log(0.33)/log(0.3)),log( log(0.33)/log(0.3)),0,0), nrow=2  )
 aw=6
 time.between=1
 A1=1/3
 A2=1/5
 A3=1/10
 ncohort1=10
 cohortsize1=3
 ncohort2=10
 cohortsize2=5
 n3=500
 u00=40
 u11=60
 targetT=0.35
 targetE=0.2
 cutoff.eli.T=0.95
 cutoff.eli.E=0.9
 targetHR=0.85
 tauS=0.8
 cutoff.eli.TE=0.5
 
 delta=0.5

 Gen123.OC(n.rep,pT.true,pE.true,pT.true.c,pE.true.c,rho,beta0,betae,betat,gammad,aw,time.between,A1,A2,A3,ncohort1,cohortsize1,ncohort2,cohortsize2,n3,u00,u11,targetT,targetE,cutoff.eli.T,cutoff.eli.E,
           targetHR,tauS,delta,cutoff.eli.TE)


# ## scenario 4.
 n.rep=5000
 pT.true=c(0.05,0.1,0.2,0.4,0.5)
 pE.true=c(0.15,0.4,0.45,0.5,0.5)
 pT.true.c=0.1
 pE.true.c=0.3
 rho=0.1
 beta0=c(1.1*log(-log(0.242)/6), 0.9*log(-log(0.242)/6))
 betae=rep(  -log(2),2)
 betat=rep(log(1.5),2)
 gammad=matrix(c(log( log(0.225)/log(0.3)),log( log(0.225)/log(0.3)),log( log(0.573)/log(0.28)),log( log(0.573)/log(0.28)),log( log(0.382)/log(0.3)),log( log(0.382)/log(0.3)),
                 log( log(0.29)/log(0.3)),log( log(0.29)/log(0.3)),log( log(0.305)/log(0.3)),log( log(0.305)/log(0.3)),0,0), nrow=2  )
 aw=6
 time.between=1
 A1=1/3
 A2=1/5
 A3=1/10
 ncohort1=10
 cohortsize1=3
 ncohort2=10
 cohortsize2=5
 n3=500
 u00=40
 u11=60
 targetT=0.35
 targetE=0.2
 cutoff.eli.T=0.95
 cutoff.eli.E=0.9
 targetHR=0.85
 tauS=0.8
 cutoff.eli.TE=0.5

 delta=0.5

 Gen123.OC(n.rep,pT.true,pE.true,pT.true.c,pE.true.c,rho,beta0,betae,betat,gammad,aw,time.between,A1,A2,A3,ncohort1,cohortsize1,ncohort2,cohortsize2,n3,u00,u11,targetT,targetE,cutoff.eli.T,cutoff.eli.E,
           targetHR,tauS,delta,cutoff.eli.TE)


## scenario 5.
 n.rep=5000
 pT.true=c(0.15,0.3,0.35,0.4,0.45)
 pE.true=c(0.4,0.6,0.4,0.2,0.1)
 pT.true.c=0.2
 pE.true.c=0.3
 rho=0.1
 beta0=c(1.1*log(-log(0.362)/6), 0.9*log(-log(0.362)/6))
 betae=rep(  -log(2),2)
 betat=rep(log(1.5),2)
 gammad=matrix(c(log( log(0.17)/log(0.3)),log( log(0.17)/log(0.3)),log( log(0.29)/log(0.35)),log( log(0.29)/log(0.35)),log( log(0.66)/log(0.7)),log( log(0.66)/log(0.7)),
                 log( log(0.775)/log(0.8)),log( log(0.775)/log(0.8)),log( log(0.73)/log(0.8)),log( log(0.73)/log(0.8)),0,0), nrow=2  )
 aw=6
 time.between=1
 A1=1/3
 A2=1/5
 A3=1/10
 ncohort1=10
 cohortsize1=3
 ncohort2=10
 cohortsize2=5
 n3=500
 u00=40
 u11=60
 targetT=0.35
 targetE=0.2
 cutoff.eli.T=0.95
 cutoff.eli.E=0.9
 targetHR=0.85
 tauS=0.8
 cutoff.eli.TE=0.5
 
 
 
 delta=0.5

 Gen123.OC(n.rep,pT.true,pE.true,pT.true.c,pE.true.c,rho,beta0,betae,betat,gammad,aw,time.between,A1,A2,A3,ncohort1,cohortsize1,ncohort2,cohortsize2,n3,u00,u11,targetT,targetE,cutoff.eli.T,cutoff.eli.E,
           targetHR,tauS,delta,cutoff.eli.TE)
 

## scenario 6.
 n.rep=5000
 pT.true=c(0.05,0.08,0.1,0.15,0.2)
 pE.true=c(0.3,0.5,0.3,0.2,0.1)
 pT.true.c=0.1
 pE.true.c=0.3
 rho=0.1
 beta0=c(1.1*log(-log(0.243)/6), 0.9*log(-log(0.243)/6))
 betae=rep(  -log(2),2)
 betat=rep(log(1.5),2)
 gammad=matrix(c(log( log(0.49)/log(0.5)),log( log(0.49)/log(0.5)),log( log(0.545)/log(0.5)),log( log(0.545)/log(0.5)),log( log(0.755)/log(0.5)),log( log(0.755)/log(0.5)),
                 log( log(0.62)/log(0.5)),log( log(0.62)/log(0.5)),log( log(0.56)/log(0.5)),log( log(0.56)/log(0.5)),0,0), nrow=2  )
 aw=6
 time.between=1
 A1=1/3
 A2=1/5
 A3=1/10
 ncohort1=10
 cohortsize1=3
 ncohort2=10
 cohortsize2=5
 n3=500
 u00=40
 u11=60
 targetT=0.35
 targetE=0.2
 cutoff.eli.T=0.95
 cutoff.eli.E=0.9
 targetHR=0.85
 tauS=0.8
 cutoff.eli.TE=0.5
 
 delta=0.5

 Gen123.OC(n.rep,pT.true,pE.true,pT.true.c,pE.true.c,rho,beta0,betae,betat,gammad,aw,time.between,A1,A2,A3,ncohort1,cohortsize1,ncohort2,cohortsize2,n3,u00,u11,targetT,targetE,cutoff.eli.T,cutoff.eli.E,
           targetHR,tauS,delta,cutoff.eli.TE)

## scenario 7.
 n.rep=5000
 pT.true=c(0.05,0.1,0.2,0.4,0.5)
 pE.true=c(0.3,0.4,0.3,0.2,0.1)
 pT.true.c=0.1
 pE.true.c=0.3
 rho=0.1
 beta0=c(1.1*log(-log(0.242)/6), 0.9*log(-log(0.242)/6))
 betae=rep(  -log(2),2)
 betat=rep(log(1.5),2)
 gammad=matrix(c(log( log(0.493)/log(0.5)),log( log(0.493)/log(0.5)),log( log(0.74)/log(0.5)),log( log(0.74)/log(0.5)),log( log(0.513)/log(0.5)),log( log(0.513)/log(0.5)),
                 log( log(0.456)/log(0.5)),log( log(0.456)/log(0.5)),log( log(0.36)/log(0.5)),log( log(0.36)/log(0.5)),0,0), nrow=2  )
 aw=6
 time.between=1
 A1=1/3
 A2=1/5
 A3=1/10
 ncohort1=10
 cohortsize1=3
 ncohort2=10
 cohortsize2=5
 n3=500
 u00=40
 u11=60
 targetT=0.35
 targetE=0.2
 cutoff.eli.T=0.95
 cutoff.eli.E=0.9
 targetHR=0.85
 tauS=0.8
 cutoff.eli.TE=0.5
 delta=0.5

 Gen123.OC(n.rep,pT.true,pE.true,pT.true.c,pE.true.c,rho,beta0,betae,betat,gammad,aw,time.between,A1,A2,A3,ncohort1,cohortsize1,ncohort2,cohortsize2,n3,u00,u11,targetT,targetE,cutoff.eli.T,cutoff.eli.E,
           targetHR,tauS,delta,cutoff.eli.TE)


## scenario 8.
 n.rep=5000
 pT.true=c(0.02,0.05,0.3,0.4,0.5)
 pE.true=c(0.2,0.4,0.4,0.4,0.4)
 pT.true.c=0.1
 pE.true.c=0.3
 rho=0.1
 beta0=c(1.1*log(-log(0.242)/6), 0.9*log(-log(0.242)/6))
 betae=rep(  -log(2),2)
 betat=rep(log(1.5),2)
 gammad=matrix(c(log( log(0.515)/log(0.5)),log( log(0.515)/log(0.5)),log( log(0.74)/log(0.5)),log( log(0.74)/log(0.5)),log( log(0.76)/log(0.5)),log( log(0.76)/log(0.5)),
                 log( log(0.765)/log(0.5)),log( log(0.765)/log(0.5)),log( log(0.775)/log(0.5)),log( log(0.775)/log(0.5)),0,0), nrow=2  )
 aw=6
 time.between=1
 A1=1/3
 A2=1/5
 A3=1/10
 ncohort1=10
 cohortsize1=3
 ncohort2=10
 cohortsize2=5
 n3=500
 u00=40
 u11=60
 targetT=0.35
 targetE=0.2
 cutoff.eli.T=0.95
 cutoff.eli.E=0.9
 targetHR=0.85
 tauS=0.8
 cutoff.eli.TE=0.5

 delta=0.5

 Gen123.OC(n.rep,pT.true,pE.true,pT.true.c,pE.true.c,rho,beta0,betae,betat,gammad,aw,time.between,A1,A2,A3,ncohort1,cohortsize1,ncohort2,cohortsize2,n3,u00,u11,targetT,targetE,cutoff.eli.T,cutoff.eli.E,
           targetHR,tauS,delta,cutoff.eli.TE)




