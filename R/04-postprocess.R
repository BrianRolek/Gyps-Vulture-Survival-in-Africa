load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Gyps Vulture Survival in Africa\\outputs\\gyps-nimble-noage.RData")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Gyps Vulture Survival in Africa\\outputs\\gyps-nimble-age.RData")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Gyps Vulture Survival in Africa\\outputs\\gyps-nimble-age-managecont-yr2.RData")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Gyps Vulture Survival in Africa\\outputs\\gyps-study.Rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Gyps Vulture Survival in Africa\\outputs\\gyps-region.Rdata")
library(MCMCvis)
library (coda)
library (ggplot2)
library (reshape2)
library (tidybayes)

pars <- c(  "beta1", "beta2", "delta1", "delta2", "delta3", "delta4",
            "mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead",
            "sigma.s", "eps.s",
            "l.s", "l.tagfail", "l.p.tagfail", "l.p.dead")

pars <- c(  "beta", "delta", "gamma", "eta",
            "mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead",
            "sigma.s", "eps.s",
            "l.s", "l.tagfail", "l.p.tagfail", "l.p.dead")

MCMCsummary(post, pars, HPD=TRUE, digits=3, hpd_prob=0.95)
MCMCsummary(post, pars, HPD=TRUE, digits=3, hpd_prob=0.85)
MCMCtrace(post, pars, pdf=F)
p <- MCMCpstr(post, pars, type="chains")
p2 <- mcmc.list(post)
MCMCplot(p2, params=c("beta", "delta", 
                      "gamma", "eta"), 
         ISB=TRUE, 
         HPD=TRUE, ci=c(85, 95))

#***************
#* Plot effect from tag age on tag failure
#***************
ta <- seq(0,15, by=0.1)
ta.sc <- (ta-4.958314)/3.749984
ni <- 4000
pred.ta <- array(NA, dim=c(length(ta.sc), ni), dimnames=list(ta, 1:ni) )
for (i in 1:length(ta.sc)){
pred.ta[i,] <- p$l.tagfail[1,] + p$beta[1,]*ta.sc[i] + p$beta[2,]*ta.sc[i]^2
}
lp.ta <- melt(pred.ta)
colnames(lp.ta)[1:2] <- c("tagage", "iter" )
lp.ta$pred <- plogis(lp.ta$value)
mn <- plogis(apply(pred.ta, 1, mean, na.rm=T))
lhdi95 <- plogis(apply(pred.ta, 1, HDInterval::hdi, na.rm=T)[1,])
uhdi95 <- plogis(apply(pred.ta, 1, HDInterval::hdi, na.rm=T)[2,])
lhdi85 <- plogis(apply(pred.ta, 1, HDInterval::hdi, na.rm=T, credMass=0.85)[1,])
uhdi85 <- plogis(apply(pred.ta, 1, HDInterval::hdi, na.rm=T, credMass=0.85)[2,])
df <- data.frame(mn=mn, 
                 lhdi95=lhdi95, uhdi95=uhdi95, 
                 lhdi85=lhdi85, uhdi85=uhdi85,
                 ta=ta)

p1 <- ggplot() + theme_minimal() +
  geom_line(data=lp.ta, aes(x=tagage, y=pred, group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  geom_line(data=df, aes(x=ta, y=mn), linewidth=2) +
  geom_line(data=df, aes(x=ta, y=lhdi85), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=uhdi85), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=lhdi95), linewidth=1, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=uhdi95), linewidth=1, linetype="dashed") +
  ylab("Probability of tag failure (monthly)") + xlab("Tag age (years)")

#****************
#* plot survival by age class
#****************
ls <- melt(p$mean.s)

p2 <- ls |>
      ggplot(aes(x = value, y = Var1)) + theme_minimal() +
      scale_y_discrete(labels=c("First year", "Subadult", "Adult")) +
      stat_halfeye(.width=c(0.85, 0.95), point_interval="mode_hdi") +
      ylab("Age class") + xlab("Survival probability (monthly)")
p3 <- ls |>
      ggplot(aes(x = value^12, y = Var1)) + theme_minimal() +
      scale_y_discrete(labels=c("First year", "Subadult", "Adult")) +
      stat_halfeye(.width=c(0.85, 0.95), point_interval="mode_hdi") +
      ylab("Age class") + xlab("Survival probability (yearly)")

# Calculate yearly survival for results
ls$yr.s<- ls$value^12
tapply(ls$yr.s, ls$Var1, mean)
tapply(ls$yr.s, ls$Var1, HDInterval::hdi)

#****************
#* plot survival in response to management/time
#****************
yr <- seq(2009,2023, by=0.5)
yr.sc <- (yr-2016)/7
ni <- 4000
pred.yr <- array(NA, dim=c(3, length(yr.sc), ni), 
                 dimnames=list(c("First year", "Subadult", "Adult"), yr, 1:ni) )
for (a in 1:3){
for (i in 1:length(yr.sc)){
  pred.yr[a,i,] <- p$l.s[a,] + 
                    p$delta[1,]*yr.sc[i] + 
                    p$delta[2,]*yr.sc[i]^2 +
                    p$delta[3,]*yr.sc[i]^3
}}
lp.yr <- melt(pred.yr)
colnames(lp.yr)[1:3] <- c("Ageclass", "yr", "iter" )
lp.yr$pred <- plogis(lp.yr$value)

mn <- plogis(apply(pred.yr, c(1,2), mean, na.rm=T))
lhdi95 <- plogis(apply(pred.yr, c(1,2), HDInterval::hdi, na.rm=T)[1,,])
uhdi95 <- plogis(apply(pred.yr, c(1,2), HDInterval::hdi, na.rm=T)[2,,])
lhdi85 <- plogis(apply(pred.yr, c(1,2), HDInterval::hdi, na.rm=T, credMass=0.85)[1,,])
uhdi85 <- plogis(apply(pred.yr, c(1,2), HDInterval::hdi, na.rm=T, credMass=0.85)[2,,])

df <- data.frame(melt(mn), 
                 lhdi95=melt(lhdi95)[,3], uhdi95=melt(uhdi95)[,3],
                 lhdi85=melt(lhdi85)[,3], uhdi85=melt(uhdi85)[,3])
colnames(df)[1:3] <- c("Ageclass", "yr", "Mean") 

p4 <- ggplot() + theme_minimal() +
  geom_vline(xintercept=2012, col="blue", linewidth=2) +
  geom_line(data=lp.yr, aes(x=yr, y=pred^12, group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  facet_wrap(facets=vars(Ageclass)) +
  geom_line(data=df, aes(x=yr, y=Mean^12), linewidth=2) +
  geom_line(data=df, aes(x=yr, y=lhdi85^12), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=uhdi85^12), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=lhdi95^12), linewidth=1, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=uhdi95^12), linewidth=1, linetype="dashed") +
  ylab("Probability of survival (yearly)") + xlab("Year")


#*#****************
#* plot survival, historic v recent
#****************