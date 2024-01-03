library ('nimble')
library('parallel')
load("/bsuscratch/brianrolek/gyps/data.RData")
#load("data\\data.RData")
set.seed(5757575)

run <- function(seed, datl){
  library('nimble')
  library('coda')

ifgreaterFun <- nimbleFunction(
    run = function(x = integer(0), 
                   cond1 = integer(0), cond2 = integer(0) # specify subadult stage here
    ){
      if(x < cond1){ ans <- 1}
      if(x >= cond1 & x < cond2){ 
          ans <- 2 
        } 
      if (x >= cond2){
          ans <- 3
        }
      return(ans)
      returnType(integer(0))
    })
  assign('ifgreaterFun', ifgreaterFun, envir = .GlobalEnv)
  
# impute median for unknown ages  
datl$first_age[is.na(datl$first_age)] <- median(1:6)  

code <- nimbleCode({ 
  # -------------------------------------------------
  # Parameters:
  # s: monthly survival probability intercept
  # tagfail: probability that tag will fail
  # p.dead: probability for carcass to be recovered
  # p.tagfail: probability of observing a tag failure
  # -------------------------------------------------
  # States (S):
  # 1 alive with functioning tag
  # 2 alive, tag failed or lost
  # 3 dead with functioning tag
  # 4 dead, tag failed or lost
  # 5 long dead
  
  # Observations (O):
  # 1 Observed alive, tag works 
  # 2 Observed alive, tag failed
  # 3 Observed dead, tag works
  # 4 Observed dead, tag failed- NO OCCURRENCES
  # 5 Not observed, uncertain tag status and fate
  # -------------------------------------------------
  # Priors and constraints
  for (xx in 1:3){
    mean.s[xx] ~ dbeta(1, 1)   # uninformative prior for all MONTHLY survival probabilities
    l.s[xx] <- logit(mean.s[xx])
  }# xx logit transformed survival intercept
for (xxx in 1:4){
  mean.tagfail[xxx] ~ dbeta(1, 1)   # uninformative prior for all MONTHLY survival probabilities
  l.tagfail[xxx] <- logit(mean.tagfail[xxx])
} # xxx
  mean.p.tagfail ~ dbeta(1, 1)   # uninformative prior for all MONTHLY survival probabilities
  l.p.tagfail <- logit(mean.p.tagfail)    # logit transformed survival intercept
  mean.p.dead ~ dbeta(1, 1)   # uninformative prior for all MONTHLY survival probabilities
  l.p.dead <- logit(mean.p.dead)    # logit transformed survival intercept

  for (x in 1:6){
    delta[x] ~ dnorm(0, sd=10) # covariates for survival        
  } # x
  for (xxxx in 1:2){
  beta[xxxx] ~ dnorm(0, sd=10) # covariates for tagloss  
} # xxxx
  gamma ~ dnorm(0, sd=10)
  eta ~ dnorm(0, sd=10)
  
  for (yr in 1:(nyears)){
    eps.s[yr] ~ dnorm(0, sd=sigma.s)
  }
  sigma.s ~ dexp(1)           
  # pi[1,] just simulates values but is unused in model
  # alpha[1:6] <- c(1,1,1,1,1,1)
  # pi[1,1:6] ~ ddirch(alpha[1:6])
  # # pi[2,] estimates age for unknown age subadults
  # pi[2,1] <- 0
  # pi[2,2:5] ~ ddirch(alpha[2:5])
  # pi[2,6] <- 0
  
  #### MONTHLY SURVIVAL PROBABILITY
  for (i in 1:nind){
   # first_age[i] ~ dcat( pi[ known[i], 1:6] )
    for (t in f[i]:ntime){
      age[i,t] <- ( first_age[i] + t/12 - f[i]/12 )
      age.class[i,t] <- ifgreaterFun( age[i,t], 1, 6 ) # translate to age classes: 0 yro first-year, 1-5 yro subadult, and >=6 yro adult
      
      logit(s[i,t]) <- l.s[ age.class[i,t] ] + 
                        #delta[1]*managed.cat[i,t] + 
                        delta[1]*year.cont[t] +
                        delta[2]*year.cont[t]^2 +
                        delta[3]*year.cont[t]^3 +
                        delta[4]*rehabbed[i] + 
                        delta[5]*sp[i] +
                        delta[6]*region[i] +
                        eps.s[ year.factor[t] ]

      logit(tagfailed[i,t]) <- l.tagfail[ study[i] ] + 
                                beta[1]*tag.age.sc[i,t] + 
                                beta[2]*tag.age.sc[i,t]^2 
      logit(p.tagfailed[i,t]) <- l.p.tagfail + gamma*region[i]  
      logit(p.dead[i,t]) <- l.p.dead + eta*region[i]  #### probability of dead recovery
    } #t
  } #i
  
  # -------------------------------------------------
  # Define state-transition and observation matrices 
  # -------------------------------------------------
  for (i in 1:nind){
    for (t in f[i]:(ntime-1)){
      # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
      ps[1,i,t,1]<-(1-tagfailed[i,t])*s[i,t]    ## dead birds stay dead
      ps[1,i,t,2]<-tagfailed[i,t]*s[i,t]
      ps[1,i,t,3]<-(1-tagfailed[i,t])*(1-s[i,t])
      ps[1,i,t,4]<-tagfailed[i,t]*(1-s[i,t])
      ps[1,i,t,5]<-0
      
      ps[2,i,t,1]<-0
      ps[2,i,t,2]<-s[i,t]
      ps[2,i,t,3]<-0
      ps[2,i,t,4]<-(1-s[i,t])
      ps[2,i,t,5]<-0
      
      ps[3,i,t,1]<-0
      ps[3,i,t,2]<-0
      ps[3,i,t,3]<-0
      ps[3,i,t,4]<-0
      ps[3,i,t,5]<-1
      
      ps[4,i,t,1]<-0
      ps[4,i,t,2]<-0
      ps[4,i,t,3]<-0
      ps[4,i,t,4]<-0
      ps[4,i,t,5]<-1
      
      ps[5,i,t,1]<-0
      ps[5,i,t,2]<-0
      ps[5,i,t,3]<-0
      ps[5,i,t,4]<-0
      ps[5,i,t,5]<-1
      
      # Define probabilities of O(t) [last dim] given S(t)  [first dim]
      po[1,i,t,1]<-1
      po[1,i,t,2]<-0
      po[1,i,t,3]<-0
      po[1,i,t,4]<-0
      po[1,i,t,5]<-0
      
      po[2,i,t,1]<-0
      po[2,i,t,2]<-p.tagfailed[i,t]
      po[2,i,t,3]<-0
      po[2,i,t,4]<-0
      po[2,i,t,5]<-(1-p.tagfailed[i,t])
      
      po[3,i,t,1]<-0
      po[3,i,t,2]<-0
      po[3,i,t,3]<-p.dead[i,t]
      po[3,i,t,4]<-0
      po[3,i,t,5]<-(1-p.dead[i,t])
      
      po[4,i,t,1]<-0
      po[4,i,t,2]<-0
      po[4,i,t,3]<-0
      po[4,i,t,4]<-p.tagfailed[i,t]*p.dead[i,t]
      po[4,i,t,5]<-(1-p.tagfailed[i,t])*(1-p.dead[i,t])
      
      po[5,i,t,1]<-0
      po[5,i,t,2]<-0
      po[5,i,t,3]<-0
      po[5,i,t,4]<-0
      po[5,i,t,5]<-1
      
    } #t
  } #i
  
  # Likelihood 
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- first_zs[i] ## alive with tag when first marked
    for (t in (f[i]+1):ntime){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:5])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,1:5])
    } #t
  } #i

} ) # nimbleCode

#fa.inits <- rep(NA, times=datl$nind)
#fa.inits[is.na(datl$first_age)] <- 3
inits <- function(){ list(z=datl$z.inits,
                         beta = rnorm(2,0,0.5),
                         delta = rnorm(6,0,0.5),
                         gamma = rnorm(1,0,0.5), 
                         eta = rnorm(1,0,0.5),
                         mean.s = runif(3,0,1), 
                         mean.tagfail = runif(4), 
                         mean.p.tagfail =  runif(1), 
                         mean.p.dead = runif(1),
                         sigma.s = runif(1) 
                         # pi = matrix(c(rdirch(1,alpha=c(1,1,1,1,1,1)), 
                         #               NA, rdirch(1,alpha=c(1,1,1,1)), NA), nrow=2, byrow=T ),
                         # first_age=fa.inits
                         
)}

pars <- c(  "beta", "delta", "gamma", "eta",
            "mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead",
            "sigma.s", "eps.s",
            "l.s", "l.tagfail", "l.p.tagfail", "l.p.dead") 

nimbleOptions(showCompilerOutput= TRUE)        
mod <- nimbleModel(code, calculate=T, constants = datl[c(4:6, 8:length(datl))],
                   data = datl[c(1,2,7)], inits = inits())
cmod <- compileNimble(mod )
conf <- configureMCMC(cmod, monitors=pars, print = TRUE)
mc <- buildMCMC(conf, project=cmod)
cmc <- compileNimble(mc, project=cmod, showCompilerOutput = TRUE)
printErrors()

nc <- 1; nt <- 50; ni <- 100000; nb <- 50000
#nc <- 1; nt <- 5; ni <- 200; nb <- 100

post <- runMCMC(cmc,
                niter = ni, 
                nburnin = nb,
                nchains = 1,
                thin = nt,
                samplesAsCodaMCMC = T)

return(post)
} # run model function end
# 
this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster,
                  X = 1:4,
                  fun = run,
                  datl=datl)
stopCluster(this_cluster)
# save(post, datl, file="outputs/gyps-nimble-study.RData")
# save(post, datl, file="outputs/gyps-nimble-age-managecat.RData")
save(post, datl, file="/bsuscratch/brianrolek/gyps/gyps-region.RData")
