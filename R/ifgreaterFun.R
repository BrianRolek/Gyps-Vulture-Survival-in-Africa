library(nimble)
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