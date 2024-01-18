
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\MunirVultures\\outputs\\gyps-cat2.RData")
library(MCMCvis)
library (coda)
library (ggplot2)
library (reshape2)
library (tidybayes)
options(scipen=999)

pars <- c(  "beta", "delta", #"gamma", "eta",
            "mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead",
            "sigma.s", 
            "l.s", "l.tagfail", "l.p.tagfail", "l.p.dead")
p <- MCMCpstr(post, pars, type="chains")
p2 <- mcmc.list(post)

MCMCsummary(post, pars, HPD=TRUE, digits=3, 
            hpd_prob=0.95, pg0=TRUE, func=median, func_name="md")
MCMCsummary(post, pars, HPD=TRUE, digits=3, 
            hpd_prob=0.85, pg0=TRUE, func=median, func_name="md")
MCMCtrace(post, pars, pdf=F)

par(mfrow=c(1,1))
MCMCplot(p2, params=c("beta", "delta", 
                      "gamma", "eta"), 
         ISB=TRUE, 
         HPD=TRUE, ci=c(85, 95))
MCMCplot(p2, params=c("mean.s", "mean.tagfail", 
                      "mean.p.tagfail", "mean.p.dead"), 
         ISB=TRUE, 
         HPD=TRUE, ci=c(85, 95))
MCMCplot(p2, params=c("mean.p.dead"), 
         ISB=TRUE, 
         HPD=TRUE, ci=c(85, 95))
#***************
#* Plot effect from tag age on tag failure
#***************
ta <- seq(0,5.3, by=0.1)
ta.sc <- (ta-5.436775)/3.997963
ni <- 4000
pred.ta <- array(NA, dim=c(length(ta.sc), ni), dimnames=list(ta, 1:ni) )
for (i in 1:length(ta.sc)){
pred.ta[i,] <- p$l.tagfail[1,] + p$beta[1,]*ta.sc[i] + p$beta[2,]*ta.sc[i]^2
}
lp.ta <- melt(pred.ta)
colnames(lp.ta)[1:2] <- c("tagage", "iter" )
lp.ta$pred <- plogis(lp.ta$value)
md <- plogis(apply(pred.ta, 1, median, na.rm=T))
lhdi95 <- plogis(apply(pred.ta, 1, HDInterval::hdi, na.rm=T)[1,])
uhdi95 <- plogis(apply(pred.ta, 1, HDInterval::hdi, na.rm=T)[2,])
lhdi85 <- plogis(apply(pred.ta, 1, HDInterval::hdi, na.rm=T, credMass=0.85)[1,])
uhdi85 <- plogis(apply(pred.ta, 1, HDInterval::hdi, na.rm=T, credMass=0.85)[2,])
df <- data.frame(md=md, 
                 lhdi95=lhdi95, uhdi95=uhdi95, 
                 lhdi85=lhdi85, uhdi85=uhdi85,
                 ta=ta)

pyta <- ggplot() + theme_minimal() +
  geom_line(data=lp.ta, aes(x=tagage, y=1-((1-pred)^12), group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  geom_line(data=df, aes(x=ta, y=1-((1-md)^12)), linewidth=2) +
  geom_line(data=df, aes(x=ta, y=1-((1-lhdi85)^12)), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=1-((1-uhdi85)^12)), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=1-((1-lhdi95)^12)), linewidth=1, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=1-((1-uhdi95)^12)), linewidth=1, linetype="dashed") +
  ylab("Transmitter failure (yearly probability)") + xlab("Transmitter age (years)")

pmta <- ggplot() + theme_minimal() +
  geom_line(data=lp.ta, aes(x=tagage, y=pred, group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  geom_line(data=df, aes(x=ta, y=md), linewidth=2) +
  geom_line(data=df, aes(x=ta, y=lhdi85), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=uhdi85), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=lhdi95), linewidth=1, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=uhdi95), linewidth=1, linetype="dashed") +
  ylab("Transmitter failure (monthly probability)") + xlab("Transmitter age (years)")

ggsave("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\MunirVultures\\docs\\figs\\transmitter age and failure.tiff",
       pyta, device="tiff", 
       width=6.5, height=4, units="in", dpi=300)

ggsave("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\MunirVultures\\docs\\figs\\transmitter age and failure-byMonth.tiff",
       pmta, device="tiff", 
       width=6.5, height=4, units="in", dpi=300)

df2 <- data.frame(md =1-((1-md)^12), 
                  lhdi95= 1-((1-lhdi95)^12),
                  uhdi95=1-((1-uhdi95)^12)    
)
#****************
#* plot survival by age class
#****************
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\MunirVultures\\outputs\\gyps-simplified.RData")
p3 <- MCMCpstr(post, pars, type="chains")
lss <- melt(p3$mean.s)

ps2 <- lss |>
      ggplot(aes(x = value, y = Var1)) + theme_minimal() +
      scale_y_discrete(labels=c("First year", "Subadult", "Adult")) +
      stat_halfeye(.width=c(0.85, 0.95), point_interval="median_hdi") +
      ylab("Age class") + xlab("Survival (monthly probability)") +
      coord_flip()
ps3 <- lss |>
      ggplot(aes(x = value^12, y = Var1)) + theme_minimal() +
      scale_y_discrete(labels=c("First year", "Subadult", "Adult")) +
      stat_halfeye(.width=c(0.85, 0.95), point_interval="median_hdi") +
      ylab("Age class") + xlab("Survival (yearly probability)") +
      coord_flip() 

# Calculate yearly survival for results
lss$yr.s<- lss$value^12
tapply(lss$yr.s, lss$Var1, median)
tapply(lss$yr.s, lss$Var1, HDInterval::hdi, credMass=0.95)
tapply(lss$yr.s, lss$Var1, HDInterval::hdi, credMass=0.85)

MCMCplot(p2, params=c("mean.s"), 
         ISB=TRUE, 
         HPD=TRUE, ci=c(85, 95))

#****************
#* plot survival in response to management/time
#****************
ni <- 4000
pred.man <- array(NA, dim=c(3, 2, ni), 
                 dimnames=list(c("First year", "Subadult", "Adult"), c("Early", "Late"), 1:ni) )
for (a in 1:3){
for (m in 1:2){
  pred.man[a,m,] <- p$l.s[a,] + 
                    p$delta[1,]*c(0,1)[m]
}}
lp.man <- melt(pred.man)
colnames(lp.man)[1:3] <- c("Ageclass", "man", "iter" )
lp.man$pred <- plogis(lp.man$value)

md <- plogis(apply(pred.man, c(1,2), median, na.rm=T))
lhdi95 <- plogis(apply(pred.man, c(1,2), HDInterval::hdi, na.rm=T)[1,,])
uhdi95 <- plogis(apply(pred.man, c(1,2), HDInterval::hdi, na.rm=T)[2,,])
lhdi85 <- plogis(apply(pred.man, c(1,2), HDInterval::hdi, na.rm=T, credMass=0.85)[1,,])
uhdi85 <- plogis(apply(pred.man, c(1,2), HDInterval::hdi, na.rm=T, credMass=0.85)[2,,])

df <- data.frame(melt(md), 
                 lhdi95=melt(lhdi95)[,3], uhdi95=melt(uhdi95)[,3],
                 lhdi85=melt(lhdi85)[,3], uhdi85=melt(uhdi85)[,3])
colnames(df)[1:3] <- c("Ageclass", "man", "median") 


p4 <- ggplot(data=lp.man, aes(x=pred^12, y=man)) + theme_minimal() +
  # geom_line(data=lp.man, aes(x=pred^12, y=man, group=iter),
  #           color="gray40", linewidth=0.5, alpha=0.05) +
  stat_halfeye(.width=c(0.85, 0.95), point_interval="median_hdi") +
  facet_wrap(facets=vars(Ageclass)) + 
  coord_flip() +
  ylab("Survival (yearly probability)") + xlab("Period")

ggsave("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\MunirVultures\\docs\\figs\\survival-ageclass-management.tiff",
       p4, device="tiff", 
       width=6.5, height=4, units="in", dpi=300)

#***************
#* Plot effect from year of study on tag failure
#***************
yr <- seq(2009,2023, by=0.1)
yr.sc <- (yr-2016)/7
ni <- 4000
pred.tf.yr <- array(NA, dim=c(length(yr.sc), ni), dimnames=list(yr, 1:ni) )
for (i in 1:length(yr.sc)){
  pred.tf.yr[i,] <- p$l.tagfail[1,] + p$beta[3,]*yr.sc[i] + p$beta[4,]*yr.sc[i]^2
}
lp.tf.yr <- melt(pred.tf.yr)
colnames(lp.tf.yr)[1:2] <- c("year", "iter" )
lp.tf.yr$pred <- plogis(lp.tf.yr$value)
md <- plogis(apply(pred.tf.yr, 1, median, na.rm=T))
lhdi95 <- plogis(apply(pred.tf.yr, 1, HDInterval::hdi, na.rm=T)[1,])
uhdi95 <- plogis(apply(pred.tf.yr, 1, HDInterval::hdi, na.rm=T)[2,])
lhdi85 <- plogis(apply(pred.tf.yr, 1, HDInterval::hdi, na.rm=T, credMass=0.85)[1,])
uhdi85 <- plogis(apply(pred.tf.yr, 1, HDInterval::hdi, na.rm=T, credMass=0.85)[2,])
df <- data.frame(md=md, 
                 lhdi95=lhdi95, uhdi95=uhdi95, 
                 lhdi85=lhdi85, uhdi85=uhdi85,
                 yr=yr, 
                 md12= 1-((1-md)^12) )

pyty <- ggplot() + theme_minimal() +
  geom_line(data=lp.tf.yr, aes(x=year, y=1-((1-pred)^12), group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  geom_line(data=df, aes(x=yr, y=1-((1-md)^12)), linewidth=2) +
  geom_line(data=df, aes(x=yr, y=1-((1-lhdi85)^12)), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=1-((1-uhdi85)^12)), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=1-((1-lhdi95)^12)), linewidth=1, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=1-((1-uhdi95)^12)), linewidth=1, linetype="dashed") +
  ylab("Transmitter failure (yearly probability)") + xlab("Year of study")

pmty <- ggplot() + theme_minimal() +
  geom_line(data=lp.tf.yr, aes(x=year, y=pred, group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  geom_line(data=df, aes(x=yr, y=md), linewidth=2) +
  geom_line(data=df, aes(x=yr, y=lhdi85), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=uhdi85), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=lhdi95), linewidth=1, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=uhdi95), linewidth=1, linetype="dashed") +
  ylab("Transmitter failure (monthly probability)") + xlab("Year of study")

ggsave("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\MunirVultures\\docs\\figs\\tagfailure-year.tiff",
       pyty, device="tiff", 
       width=6.5, height=4, units="in", dpi=300)

ggsave("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\MunirVultures\\docs\\figs\\tagfailure-month.tiff",
       pmty, device="tiff", 
       width=6.5, height=4, units="in", dpi=300)

df2 <- data.frame(md =1-((1-md)^12), 
                  lhdi95= 1-((1-lhdi95)^12),
                  uhdi95=1-((1-uhdi95)^12)    
)

#*#****************
#* survival by 
#****************