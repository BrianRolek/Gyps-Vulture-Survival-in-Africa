# ---- mpm
######################################################
# Vulture Trends Matrix Model
# Survival rates are for African white-backed and Ruppell's vultures
######################################################
# load required packages
library(popbio)
library (MCMCvis)
library(HDInterval)
library(ggplot2)
library(ggpubr)
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
pg <- data.frame(pop.growth.rate= pop.growth.rate)
pg.hist <- ggplot(pg, aes(x=pop.growth.rate)) +
              geom_histogram() +
              geom_vline(xintercept=1.0, linewidth=1.5, linetype="dashed") +
              xlab("Population growth rate") +
              ylab("Count") +
              theme_classic() +
              annotate("text", x=0.75, y=275, size=7,
                       label= paste0("Prob<1 = ", round(mean(pg$pop.growth.rate<1),3))) +
              annotate("text", label="A", size=15, 
                       x = 0.7, y = 350) + #, hjust = 1, vjust = 1) +
              theme(axis.text=element_text(size=20), #change font size of axis text
                    axis.title=element_text(size=20)) #change font size of axis titles

paste0("Probability of decline (lambda<1) = ",
round(mean(pop.growth.rate<1),3)) |> print() # probability of lambda <1

# Calculate the expected percent decline from the 
# population growth rate. 
pg$perc.change.studyperiod <- (pg$pop.growth.rate^(2024-2009)-1)*100 # entire study period
pg$perc.change.annual <- (pg$pop.growth.rate-1)*100 # yearly
pd.hist <- ggplot(pg, aes(x=perc.change.studyperiod)) +
  geom_histogram() +
  xlab("Percent change over 15 years") +
  ylab("Count") +
  theme_classic() +
  annotate("text", label="B", size=15, 
           x = -Inf, y = Inf, hjust = -0.25, vjust = 1) +
  geom_vline(xintercept=0, linewidth=1.5, linetype="dashed") +
  theme(axis.text=element_text(size=20), #change font size of axis text
          axis.title=element_text(size=20)) #change font size of axis titles

# Create composite
pg_pd <- ggarrange(pg.hist, pd.hist, nrow=2)
pg_pd
# ggsave("figs\\pop-growth-both.tiff",
#        pg_pd, device="tiff",
#        width=6.5, height=6.5, units="in", dpi=300)
# Print annual percent change
print("median yearly percent change [95% HDIs] = ")
perc.ann <- paste0(median(pg$perc.change.annual) |> round(3), " [ ",
        HDInterval::hdi(pg$perc.change.annual)[1] |> round(3), ", ",
        HDInterval::hdi(pg$perc.change.annual)[2] |> round(3),
        " ]")
perc.ann
# Print study percent change
perc.period <- paste0( median(pg$perc.change.studyperiod) |> round(3), " [ ",
       HDInterval::hdi(pg$perc.change.studyperiod)[1] |> round(3), ", ",
       HDInterval::hdi(pg$perc.change.studyperiod)[2] |> round(3),
       " ]")
print("median percent change over full study [95% HDIs] = ")
perc.period

# Transform estimates from Virani et al. 2011
# their estimates are percent change 
# spanning 1976/1988 to 2003/2005
print("Conversions from Virani et al. 2011")
print("Mean percent change over longest duration of study (1976-2005) = ")
(( (1-0.52)^(1/(2005-1976))-1 )*100 ) |> round(digits=2) # longest span of years
print("Mean percent change over shortest duration of study (1988-2003) = ")
(( (1-0.52)^(1/(2003-1988))-1 )*100 ) |> round(digits=1) # briefest span of years
