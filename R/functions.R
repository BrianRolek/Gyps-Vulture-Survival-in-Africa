#*****************
#* Create a nimble function 
#* to reclassify ages
#* into 3 cohorts
#*****************
library('nimble')
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

#**************
#* Create a function to check 
#* that all rows sum to one
#* when generating random probabilities
#**************

sumtoone.func <- function(p.tagfailed= runif(1), tagfailed= runif(1), 
                          s=runif(1), p.dead=runif(1)){
  
  ps <- array(NA, dim=c(4,4))
  po <- array(NA, dim=c(4,5))
  # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
  ps[1,1]<-(1-tagfailed)*s    ## dead birds stay dead
  ps[1,2]<-tagfailed*s
  ps[1,3]<-(1-tagfailed)*(1-s)
  ps[1,4]<-tagfailed*(1-s)
  
  ps[2,1]<-0
  ps[2,2]<-1
  ps[2,3]<-0
  ps[2,4]<-0
  
  ps[3,1]<-0
  ps[3,2]<-0
  ps[3,3]<-1
  ps[3,4]<-0
  
  ps[4,1]<-0
  ps[4,2]<-0
  ps[4,3]<-0
  ps[4,4]<-1
  
  # Define probabilities of O(t) [last dim] given S(t)  [first dim]
  po[1,1]<-1
  po[1,2]<-0
  po[1,3]<-0
  po[1,4]<-0
  po[1,5]<-0
  
  po[2,1]<-0
  po[2,2]<-p.tagfailed
  po[2,3]<-0
  po[2,4]<-0
  po[2,5]<-(1-p.tagfailed)
  
  po[3,1]<-0
  po[3,2]<-0
  po[3,3]<-p.dead
  po[3,4]<-0
  po[3,5]<-(1-p.dead)
  
  po[4,1]<-0
  po[4,2]<-p.tagfailed*(1-p.dead)
  po[4,3]<-(1-p.tagfailed)*p.dead
  po[4,4]<-p.tagfailed*p.dead
  po[4,5]<-(1-p.tagfailed)*(1-p.dead)
  
  return(list(truestates.rowsums=rowSums(ps), 
              obsstates.rowsums=rowSums(po), 
              truestates.probs=ps, 
              obsstates.probs=po) )
}

#*********************
# Create a function to calculate
# WAIC
#*********************
waic.func <- function(flnm) {
  library ('nimble')
  library('nimbleEcology')
  load(flnm)
  mod <- nimbleModel(code, # code
                     constants = constl, 
                     data = datl, 
                     buildDerivs = FALSE,
                     calculate = T) 
  
  cmod <- compileNimble(mod, showCompilerOutput = TRUE)
  post1 <- list(post[[1]]$samples,
                post[[2]]$samples,
                post[[3]]$samples,
                post[[4]]$samples)
  
  samples <- do.call(rbind, post1)
  waic <- calculateWAIC(samples, cmod)
  return(waic) }