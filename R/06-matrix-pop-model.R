# ---- mpm
######################################################
# Vulture Trends Matrix Model
# Survival rates are for African white-backed and Ruppell's vultures
######################################################
# load required packages
library(popbio)
library (MCMCvis)
library(HDInterval)
# load output from reduced survival model
load("outputs\\gyps-28Apr2026-marginalized-reduced.RData")
postl.reduced <- lapply(post, function(x){ x$samples })
post.reduced <- do.call(rbind, postl.reduced)
p3 <- MCMCpstr(postl.reduced, "mean.s", type="chains")

#### Matrix values ----------------------------------------
# these are equal across the six matrices 
f1 <- 0 # n females fledged per female aged 1 to 4
f2 <- 0.6 / 2 # n females fledged per female aged >=5
# Convert monthly survival to yearly
s0 <- p3$mean.s[1,]^12 # survival prob. from fledging until the following year
s1 <- p3$mean.s[2,]^12 # survival prob. from age 1-5
s2 <- p3$mean.s[3,]^12# survival prob. for older birds

# Pre-breeding pulse matrix model
A.list <- list()
for (i in 1:length(s1)){
A <-matrix(c(
  s0[i] * f1, s0[i] * f1, s0[i] * f1, s0[i] * f1, s0[i] * f2,
  s1[i], 0, 0, 0, 0,
  0, s1[i], 0, 0, 0, 
  0, 0, s1[i], 0, 0,  
  0, 0, 0, s1[i], s2[i]),nrow=5,byrow=TRUE)

colnames(A) <- c("1","2","3","4","5+")
A.list[[i]] <- A
}
eigen_results <- lapply(A.list, eigen.analysis)

#' pull out the lambda values 
get_lambda <- function(x){
  lambdas <- eigen_results[[x]]$lambda1
}
runs <- c(1:length(s1)) # matrix indexed by number
lambda_values <- lapply(runs, get_lambda)

# pop.growth rate contains the posterior draws of lambda
# where one is a stable population
pop.growth.rate <- lambda_values |> unlist() 
paste0("Geometric mean population growth rate = ",
round(exp(mean(log(pop.growth.rate))), 3)) |> print() # Geometric mean

paste0("Median = ", round(median(pop.growth.rate),3)) |> print()
paste0("95% HDIs = ", round(hdi(pop.growth.rate),3)) |> print()
hist(pop.growth.rate, main="Population growth rate")
paste0("Probability of decline (lambda<1) = ",
round(mean(pop.growth.rate<1),3)) |> print() # probability of lambda <1
