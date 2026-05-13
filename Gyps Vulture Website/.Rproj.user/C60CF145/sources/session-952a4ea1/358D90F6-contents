# ---- sens analysis
library ('nimble')
library('parallel')
library('nimbleEcology')
library('MCMCvis')
load("data\\data.RData")
set.seed(5757575)

#************************
#* Create model and run function 
#* for sensitivity tests
#************************
# model code
code <- nimbleCode({ 
  # Priors and constraints
  mean.p.dead ~ dbeta(1, 1)   # uninformative prior for monthly probabilities
  l.p.dead <- logit(mean.p.dead)    # logit transformed intercept
  
  mean.tagfail ~ dbeta(1, 1)   # uninformative prior for monthly probabilities
  l.tagfail <- logit(mean.tagfail) # logit transformed intercept
  
  mean.p.tagfail ~ dbeta(1, 1)   # uninformative prior for monthly probabilities
  l.p.tagfail <- logit(mean.p.tagfail)    # logit transformed intercept
  
  for (x in 1:3){
    mean.s[x] ~ dbeta(1, 1)   # uninformative prior for monthly probabilities
    l.s[x] <- logit(mean.s[x])
  }# xx logit transformed intercept
  for (xxx in 1:2){
    beta[xxx] ~ dnorm(0, sd=10) # covariates for tagloss
  } # xxxx
  # ages 0,1,2,3,4,5,6
  alpha[1:7] <- c(1,1,1,1,1,1,1)
  pi[1,1:7] ~ ddirch(alpha[1:7])
  pi[2,1] <- 0
  pi[2,2:6] ~ ddirch(alpha[2:6])
  pi[2,7] <- 0
  
  #### MONTHLY SURVIVAL PROBABILITY
  for (i in 1:nind){
    first_age[i] ~ dcat(pi[known[i],1:7])
    for (t in f[i]:(ntime-1)){
      # translate to age classes: 
      # 0 year old = 1, 
      # 1-5 year old subadult = 2-6
      # and >=6 year old adult = 7
      age[i,t] <- ( (first_age[i]-1) + t/12 - f[i]/12 )
      age.class[i,t] <- ifgreaterFun( age[i,t], 1, 6 ) 
      
      logit(s[i,t]) <- l.s[ age.class[i,t] ] 
      
      logit(tagfailed[i,t]) <- l.tagfail + 
        beta[1]*tag.age.sc[i,t] +
        beta[2]*year.cont[t] 
      logit(p.tagfailed[i,t]) <- l.p.tagfail  # prob of detecting tag failure
      logit(p.dead[i,t]) <- l.p.dead  # probability of dead recovery
    } #t
  } #i
  
  # -------------------------------------------------
  # Define state-transition and observation matrices 
  # -------------------------------------------------
  for (i in 1:nind){
    for (t in f[i]:(ntime-1)){
      # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
      ps[i,1,1,t]<-(1-tagfailed[i,t])*s[i,t]  
      ps[i,1,2,t]<-tagfailed[i,t]*s[i,t]
      ps[i,1,3,t]<-(1-tagfailed[i,t])*(1-s[i,t])
      ps[i,1,4,t]<-tagfailed[i,t]*(1-s[i,t])
      
      ps[i,2,1,t]<-0
      ps[i,2,2,t]<-1
      ps[i,2,3,t]<-0
      ps[i,2,4,t]<-0
      
      ps[i,3,1,t]<-0
      ps[i,3,2,t]<-0
      ps[i,3,3,t]<-1
      ps[i,3,4,t]<-0
      
      ps[i,4,1,t]<-0
      ps[i,4,2,t]<-0
      ps[i,4,3,t]<-0
      ps[i,4,4,t]<-1
      
      # Define probabilities of O(t) [last dim] given S(t)  [first dim]
      po[i,1,1,t]<-1
      po[i,1,2,t]<-0
      po[i,1,3,t]<-0
      po[i,1,4,t]<-0
      po[i,1,5,t]<-0
      
      po[i,2,1,t]<-0
      po[i,2,2,t]<-p.tagfailed[i,t]
      po[i,2,3,t]<-0
      po[i,2,4,t]<-0
      po[i,2,5,t]<-(1-p.tagfailed[i,t])
      
      po[i,3,1,t]<-0
      po[i,3,2,t]<-0
      po[i,3,3,t]<-p.dead[i,t]
      po[i,3,4,t]<-0
      po[i,3,5,t]<-(1-p.dead[i,t])
      
      po[i,4,1,t]<-0
      po[i,4,2,t]<-p.tagfailed[i,t]*(1-p.dead[i,t])
      po[i,4,3,t]<-(1-p.tagfailed[i,t])*p.dead[i,t]
      po[i,4,4,t]<-p.tagfailed[i,t]*p.dead[i,t]
      po[i,4,5,t]<-(1-p.tagfailed[i,t])*(1-p.dead[i,t])
      
    } #t
  } #i
  
  # Likelihood 
  for (i in 1:nind){
    y[i, (f[i] + 1):last[i]] ~ 
      dDHMMo(init = ps[i, y.first[i], 1:4, f[i]], 
             probObs = po[i, 1:4, 1:5, f[i]:(last[i] - 1)], 
             probTrans = ps[i, 1:4, 1:4, (f[i] + 1):(last[i] - 1)], 
             len = last[i] - f[i], 
             checkRowSums = 0)
  } #i
  
} ) # nimbleCode

# function to run model
run <- function(seed, datl, constl, code){
  library('nimble')
  library('coda')
  library('nimbleEcology')
  source("R/functions.R")
  
  inits <- function(){ list(beta = rnorm(2,0,0.5),
                            mean.s = runif(3,0,1), 
                            mean.tagfail = runif(1), 
                            mean.p.tagfail = runif(1), 
                            mean.p.dead = runif(1),
                            alpha = rep(1/7, 7),
                            first_age = ifelse(is.na(datl$first_age), 3, NA)
  )}
  
  pars <- c(  "beta",  
              "mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead",
              "l.s", "l.tagfail", "l.p.tagfail", "l.p.dead") 
  
  nimbleOptions(showCompilerOutput= TRUE)        
  mod <- nimbleModel(code, calculate=T, 
                     constants = constl,
                     data = datl, 
                     inits = inits()) 
  
  cmod <- compileNimble( mod )
  conf <- configureMCMC(cmod, monitors=pars, print = TRUE)
  mc <- buildMCMC(conf, project=cmod)
  cmc <- compileNimble(mc, project=cmod, showCompilerOutput = TRUE)
  printErrors()
  
  nc <- 1; nt <- 20; ni <- 100000; nb <- 80000
  
  post <- runMCMC(cmc,
                  niter = ni, 
                  nburnin = nb,
                  nchains = 1,
                  thin = nt,
                  samplesAsCodaMCMC = T)
  
  return(post)
} # run model function end

#************************
#* Sensitivity test state 2 to state 3
#************************
load("data\\data.RData")
set.seed(5757575)
# which y's are 2s
# and which are in 2009-2011
which(datl$y==2, arr.ind=T)
# Sensitivity test
# reassign all 
# (2) alive without a functioning transmitter 
# to (3) recently dead with a functioning transmitter
datl$y[datl$y==2] <- 3

this_cluster <- makeCluster(4)
post1 <- parLapply(cl = this_cluster,
                   X = 1:4,
                   fun = run,
                   datl=datl,
                   constl=constl, 
                   code=code)

#************************
#* Sensitivity test state 2 to state 4
#************************
load("data\\data.RData")
#load("/bsuscratch/brianrolek/gyps/data.RData")
set.seed(5757575)
# which y's are 2s
# and which are in 2009-2011
which(datl$y==2, arr.ind=T)
# Sensitivity test
# reassign all 
# (2) alive without a functioning transmitter  
# to (4) recently dead without a functioning transmitter
datl$y[datl$y==2] <- 4
post2 <- parLapply(cl = this_cluster,
                  X = 1:4,
                  fun = run,
                  datl=datl, 
                  constl=constl, 
                  code=code)

#************************
#* Sensitivity test state 2 to state 5
#************************
load("data\\data.RData")
#load("/bsuscratch/brianrolek/gyps/data.RData")
set.seed(5757575)
# which y's are 2s
# and which are in 2009-2011
which(datl$y==2, arr.ind=T)
# Sensitivity test
# reassign all 
# (2) alive without a functioning transmitter  
# to (3) recently dead with a functioning transmitter 
datl$y[datl$y==2] <- 5

post3 <- parLapply(cl = this_cluster,
                  X = 1:4,
                  fun = run,
                  datl=datl, 
                  constl=constl,
                  code=code)
stopCluster(this_cluster)

save(post1, post2, post3, datl, file="outputs/gyps-28Apr2026-sensitivitytest.RData")

# ---- sens table
#********************
#* Sensitivity 
#********************
#* Compare results models with 
#* those having changed data
#* state two to three
library ('MCMCvis')
load("outputs/gyps-28Apr2026-sensitivitytest.RData")
pars <- c(  "beta",  
            "mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead",
            "l.s", "l.tagfail", "l.p.tagfail", "l.p.dead")
post.sens23 <- post1 
p23 <- MCMCpstr(post1, pars[-1], type="chains")

iters <- ncol(p23$beta)
sum95.sens23 <- MCMCsummary(post.sens23, pars[-c(1,7:10)], HPD=TRUE, digits=2, 
                            hpd_prob=0.95, pg0=TRUE, func=median, func_name="md")
coef.est.sens23 <- data.frame(Model="State 2 to 3",
                              Parameter= rownames(sum95.sens23),
                              Median=sum95.sens23$md, 
                              Mean=sum95.sens23$mean,
                              LHDI95=sum95.sens23$`95%_HPDL`, 
                              UHDI95=sum95.sens23$`95%_HPDU`,
                              p= sum95.sens23$`p>0`, 
                              Rhat=sum95.sens23$Rhat)

post.sens24 <- post2
p24 <- MCMCpstr(post2, pars[-1], type="chains")
iters <- ncol(p24$beta)
sum95.sens24 <- MCMCsummary(post.sens24, pars[-c(1,7:10)], HPD=TRUE, digits=2, 
                            hpd_prob=0.95, pg0=TRUE, func=median, func_name="md")
coef.est.sens24 <- data.frame( Model= "State 2 to 4",
                               Parameter= rownames(sum95.sens24),
                               Median=sum95.sens24$md, 
                               Mean=sum95.sens24$mean,
                               LHDI95=sum95.sens24$`95%_HPDL`, 
                               UHDI95=sum95.sens24$`95%_HPDU`,
                               p= sum95.sens24$`p>0`, 
                               Rhat=sum95.sens24$Rhat)

post.sens25 <- post3
p25 <- MCMCpstr(post3, pars[-1], type="chains")
iters <- ncol(p25$beta)
sum95.sens25 <- MCMCsummary(post.sens25, pars[-c(1,7:10)], HPD=TRUE, digits=2, 
                            hpd_prob=0.95, pg0=TRUE, func=median, func_name="md")
coef.est.sens25 <- data.frame(Model= "State 2 to 5",
                              Parameter= rownames(sum95.sens25),
                              Median=sum95.sens25$md, 
                              Mean=sum95.sens25$mean,
                              LHDI95=sum95.sens25$`95%_HPDL`, 
                              UHDI95=sum95.sens25$`95%_HPDU`,
                              p= sum95.sens25$`p>0`, 
                              Rhat=sum95.sens25$Rhat
)

load("outputs\\gyps-28Apr2026-marginalized-reduced.RData")
post.sens.red <- lapply(post, function(x){ x$samples })
pr <- MCMCpstr(post.sens.red, pars, type="chains")
iters <- ncol(pr$beta)
sum95.sens.r <- MCMCsummary(post.sens.red, pars[-c(1,7:10)], HPD=TRUE, digits=2, 
                            hpd_prob=0.95, pg0=TRUE, func=median, func_name="md")
coef.est.reduced <- data.frame(Model= "Reduced",
                              Parameter= rownames(sum95.sens.r),
                              Median=sum95.sens.r$md, 
                              Mean=sum95.sens.r$mean,
                              LHDI95=sum95.sens.r$`95%_HPDL`, 
                              UHDI95=sum95.sens.r$`95%_HPDU`,
                              p= sum95.sens.r$`p>0`, 
                              Rhat=sum95.sens.r$Rhat
)

# Compare outputs
df.sens.tab <- rbind(coef.est.sens23[1:3,1:6],
                    coef.est.sens24[1:3,1:6],
                    coef.est.sens25[1:3,1:6],
                    coef.est.reduced[1:3,1:6])
df.sens.tab[,3:6] <- df.sens.tab[,3:6] |> round(3)
df.sens.tab
# write.csv(file= "C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Gyps Vulture Survival in Africa\\docs\\sensitivity-table.csv", 
#           df.compare)
