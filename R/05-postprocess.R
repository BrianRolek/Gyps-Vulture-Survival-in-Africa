# ---- setup2
library (MCMCvis)
library (coda)
library (ggplot2)
library (reshape2)
library (tidybayes)
library (bayestestR)
library (ggpubr)
library (scales)
options(scipen=999)
load("data//data.RData")
source("R//functions.R")
pars <- c(  "delta", "beta",
            "mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead",
            "l.s", "l.tagfail", "l.p.tagfail", "l.p.dead")
# load output from global model
load("outputs\\gyps-28Apr2026-marginalized-global.RData")
postl.global <- lapply(post, function(x){ x$samples })
post.global <- do.call(rbind, postl.global)
p1 <- MCMCpstr(postl.global, pars, type="chains")
p2 <- mcmc.list(postl.global)
# load output from reduced model
load("outputs\\gyps-28Apr2026-marginalized-reduced.RData")
postl.reduced <- lapply(post, function(x){ x$samples })
post.reduced <- do.call(rbind, postl.reduced)
p3 <- MCMCpstr(postl.reduced, pars[-1], type="chains")
p4 <- mcmc.list(postl.reduced)


# Model diagnostics
# ---- sumtoone
sumtoone.func()
# Check for convergence
# Priors are depicted in red
# ---- traceplots global
iters <- ncol(p1$delta)
MCMCtrace(postl.global, "delta", pdf=F, Rhat=T, 
          priors=rnorm(iters, 0, 10), post_zm = FALSE)
MCMCtrace(postl.global, "beta", pdf=F, Rhat=T, 
          priors=rnorm(iters, 0, 10), post_zm = FALSE)       
MCMCtrace(postl.global, c("mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead"), pdf=F, Rhat=T, 
          priors=rbeta(iters, 1, 1), post_zm = FALSE)  

# ---- traceplots reduced
iters <- ncol(p3$beta)
MCMCtrace(postl.reduced, "beta", pdf=F, Rhat=T, 
          priors=rnorm(iters, 0, 10), post_zm = FALSE)       
MCMCtrace(postl.reduced, c("mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead"), pdf=F, Rhat=T, 
          priors=rbeta(iters, 1, 1), post_zm = FALSE) 

# ---- table S1 
# Estimates
# Survival by management period
ps <- c(  "delta", "beta",
          "mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead")
pars <- list( ps[-c(1,2)], ps, ps[-1])
flnms <- c(null="outputs/gyps-28Apr2026-marginalized-null.RData", 
           global="outputs/gyps-28Apr2026-marginalized-global.RData",
           reduced="outputs/gyps-28Apr2026-marginalized-reduced.RData")
coef.ests <- list()
for (i in 1:length(flnms)){
  load(flnms[i])
  postl <- lapply(post, function(x){ x$samples })
  post.samps <- do.call(rbind, postl)
  sum95 <- MCMCsummary(postl, pars[[i]], HPD=TRUE, digits=2, 
                       hpd_prob=0.95, pg0=TRUE, 
                       func=median, func_name="median",
                       exact=FALSE, ISB=TRUE)
  
  coef.ests[[i]] <- cbind(Model=names(flnms)[i], sum95)
} # end loop i
names(coef.ests) <- c("null", "global", "reduced")
coef.ests2 <- do.call(rbind, coef.ests)
coef.ests2
# write.csv(file= "docs\\Coef_ests.csv",
#           coef.ests2)

# ---- cat plots
par(mfrow=c(1,1))
MCMCplot(p2, params=c("delta", 
                      "beta"), 
         ISB=TRUE, 
         HPD=TRUE, ci=c(85, 95))
MCMCplot(p2, params=c("mean.s", "mean.tagfail", 
                      "mean.p.tagfail", "mean.p.dead"), 
         ISB=TRUE, 
         HPD=TRUE, ci=c(85, 95))

# ---- tag failure age
#***************
#* Plot effect from tag age 
#* on tag failure
#* from reduced model
#***************
ta <- seq(0,5.3, by=0.1)
ta.sc <- (ta-5.436775)/3.997963
ni <- ncol(p3$beta)
yrs <- c(-0.73, 0.46) # set to 2011 for early period and 2020 for late
pred.ta <- array(NA, dim=c(length(ta.sc),2,ni), 
                 dimnames=list(tagage=ta, period=c("Early", "Late"), iter=1:ni) )
for (i in 1:length(ta.sc)){
  for (j in 1:2){
pred.ta[i,j,] <- p3$l.tagfail[1,] + p3$beta[1,]*ta.sc[i] + p3$beta[2,]*yrs[j]
}}
lp.ta <- as.data.frame.table(pred.ta, responseName = "value") 
lp.ta$tagage <- as.numeric(as.character(lp.ta$tagage))
lp.ta$pred <- plogis(lp.ta$value)
md <- plogis(apply(pred.ta, c(1,2), median, na.rm=T))
lhdi95 <- plogis(apply(pred.ta, c(1,2), HDInterval::hdi, na.rm=T)[1,,])
uhdi95 <- plogis(apply(pred.ta, c(1,2), HDInterval::hdi, na.rm=T)[2,,])
lhdi85 <- plogis(apply(pred.ta, c(1,2), HDInterval::hdi, na.rm=T, credMass=0.85)[1,,])
uhdi85 <- plogis(apply(pred.ta, c(1,2), HDInterval::hdi, na.rm=T, credMass=0.85)[2,,])
df <- data.frame(md=c(md[,1], md[,2]),
                 lhdi95=c(lhdi95[,1],lhdi95[,2]), 
                 uhdi95=c(uhdi95[,1],uhdi95[,2]),
                 lhdi85=c(lhdi85[,1],lhdi85[,2]), 
                 uhdi85=c(uhdi85[,1], uhdi85[,2]),
                 ta=c(ta, ta),
                 period=c(rep("Early", length(ta)), rep("Late", length(ta))) )

pyta <- ggplot() + theme_minimal() +
  geom_line(data=lp.ta, aes(x=tagage, y=1-((1-pred)^12), group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  scale_x_continuous(breaks=c(0,2,4,6)) +
  geom_line(data=df, aes(x=ta, y=1-((1-md)^12)), linewidth=1) +
  geom_line(data=df, aes(x=ta, y=1-((1-lhdi95)^12)), linewidth=0.5, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=1-((1-uhdi95)^12)), linewidth=0.5, linetype="dashed") +
  ylab("Transmitter failure (yearly probability)") + xlab("Transmitter age (years)") +
  facet_wrap("period")

pmta <- ggplot() + theme_minimal() +
  geom_line(data=lp.ta, aes(x=tagage, y=pred, group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  geom_line(data=df, aes(x=ta, y=md), linewidth=1) +
  geom_line(data=df, aes(x=ta, y=lhdi95), linewidth=0.5, linetype="dashed") +
  geom_line(data=df, aes(x=ta, y=uhdi95), linewidth=0.5, linetype="dashed") +
  ylab("Transmitter failure (monthly probability)") + xlab("Transmitter age (years)") +
  facet_wrap("period")


pyta
pmta
# ggsave("figs\\transmitter age and failure.tiff",
#        pyta, device="tiff",
#        width=6.5, height=4, units="in", dpi=300)
# 
# ggsave("figs\\transmitter age and failure-byMonth.tiff",
#        pmta, device="tiff",
#        width=6.5, height=4, units="in", dpi=300)

df2 <- data.frame(early.md = (1-((1-md)^12))[,1],
                  early.lhdi95= (1-((1-lhdi95)^12))[,1],
                  early.uhdi95= (1-((1-uhdi95)^12))[,1],
                  late.md = (1-((1-md)^12))[,2],
                  late.lhdi95= (1-((1-lhdi95)^12))[,2],
                  late.uhdi95= (1-((1-uhdi95)^12))[,2]
)
# df2 |> round(2)

# ---- tag failure year
#***************
#* Plot effect from year of study 
#* on tag failure
#***************
yr <- seq(2009,2025, by=0.1)
yr.sc <- (yr-2016)/7
ni <- ncol(p3$beta)
pred.tf.yr <- array(NA, dim=c(length(yr.sc), ni), 
                    dimnames=list(year=yr, iter=1:ni) )
for (i in 1:length(yr.sc)){
  pred.tf.yr[i,] <- p3$l.tagfail[1,] + p3$beta[1,]*-1.355 + p3$beta[2,]*yr.sc[i] 
}
lp.tf.yr <- as.data.frame.table(pred.tf.yr, responseName = "value") 
lp.tf.yr$year <- as.numeric(as.character(lp.tf.yr$year))
lp.tf.yr$pred <- plogis(lp.tf.yr$value)
lp.tf.yr <- lp.tf.yr[ (lp.tf.yr$year>=2009 & lp.tf.yr$year<=2011) |
                        (lp.tf.yr$year>=2017 & lp.tf.yr$year<=2024),  ]
lp.tf.yr$Period <- factor(ifelse(lp.tf.yr$year<=2011, "Early", "Late"), 
                          levels=c("Early", "Late"))
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
sub.df <- (df$yr >=2009 & df$yr<=2011) |
  (df$yr>=2017 & df$yr<=2024)
df <- df[sub.df, ]  
df$Period <- factor(ifelse(df$yr<=2011, "Early", "Late"), 
                    levels=c("Early", "Late"))

pyty <- ggplot() + theme_minimal() +
  geom_line(data=lp.tf.yr, aes(x=year, y=1-((1-pred)^12), group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  geom_line(data=df, aes(x=yr, y=1-((1-md)^12)), linewidth=1) +
  geom_line(data=df, aes(x=yr, y=1-((1-lhdi95)^12)), linewidth=0.5, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=1-((1-uhdi95)^12)), linewidth=0.5, linetype="dashed") +
  ylab("Transmitter failure (yearly probability)") + xlab("Year of study") +
  scale_x_continuous(breaks=c(2009:2011, 2018, 2020, 2022, 2024)) +
  facet_wrap("Period", scales="free_x")

pmty <- ggplot() + theme_minimal() +
  geom_line(data=lp.tf.yr, aes(x=year, y=pred, group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  geom_line(data=df, aes(x=yr, y=md), linewidth=1) +
  geom_line(data=df, aes(x=yr, y=lhdi95), linewidth=0.5, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=uhdi95), linewidth=0.5, linetype="dashed") +
  ylab("Transmitter failure (monthly probability)") + xlab("Year of study") +
  scale_x_continuous(breaks=c(2009:2011, 2018, 2020, 2022, 2024)) +
  facet_wrap("Period", scales="free_x")

pyty
pmty
# ggsave("figs\\tagfailure-year.tiff",
#        pyty, device="tiff", 
#        width=6.5, height=4, units="in", dpi=300)
# 
# ggsave("figs\\tagfailure-month.tiff",
#        pmty, device="tiff", 
#        width=6.5, height=4, units="in", dpi=300)

df2 <- data.frame(md =1-((1-md)^12), 
                  lhdi95= 1-((1-lhdi95)^12),
                  uhdi95=1-((1-uhdi95)^12) )

# ---- survival period
#****************
#* plot survival in response to period
#****************
lss <- as.data.frame.table(p3$mean.s, responseName = "value")

ps2 <- lss |>
  ggplot(aes(x = value, y = Var1)) + theme_minimal() +
  scale_y_discrete(labels=c("First year", "Subadult", "Adult")) +
  stat_halfeye(.width=c(0.85, 0.95), point_interval="median_hdci") +
  ylab("Age class") + xlab("Survival (monthly probability)") +
  coord_flip() 

ps3 <- lss |>
  ggplot(aes(x = value^12, y = Var1)) + 
  theme_minimal() +
  scale_y_discrete(labels=c("", "", "")) +
  stat_halfeye(.width=c(0.85, 0.95), point_interval="median_hdci") +
  ylab("") + xlab("Survival (yearly probability)") +
  xlim(0, 1) +
  coord_flip() +
  ggtitle("(B) Combined") 

ni <- ncol(p1$delta)
pred.man <- array(NA, dim=c(3, 2, ni), 
                 dimnames=list(c("First year", "Subadult", "Adult"), c("Early", "Late"), 1:ni) )
for (a in 1:3){
for (m in 1:2){
  pred.man[a,m,] <- p1$l.s[a,] + 
                    p1$delta[1,]*c(-1,1)[m]
}}
lp.man <- as.data.frame.table(pred.man, responseName = "value")
colnames(lp.man)[1:3] <- c("Ageclass", "Period", "iter" )
lp.man$pred <- plogis(lp.man$value)

p4 <- ggplot(data=lp.man, aes(x=pred^12, y=Period)) + theme_minimal() +
  geom_line(data=lp.man, aes(x=pred^12, y=Period, group=iter),
            color="gray40", linewidth=0.5, alpha=0.025) +
  stat_pointinterval(.width=c(0.85, 0.95), 
                     point_interval ="median_hdci") +
  facet_wrap(facets=vars(Ageclass)) + 
  xlim(0,1) +
  coord_flip() +
  ylab("Period") + xlab("Survival (yearly probability)") +
  ggtitle("(A) Period")

all_p <- ggarrange(p4, ps3, nrow=2)

all_p
# ggsave("figs\\survival-ageclass-period-combined.tiff",
#        all_p, device="jpeg", 
#        width=6, height=6, units="in", dpi=300)

# ---- survival estimates 
# survival by age class
# and management period
md <- plogis(apply(pred.man, c(1,2), median, na.rm=T))^12
mn <- plogis(apply(pred.man, c(1,2), mean, na.rm=T))^12
lhdi95 <- plogis(apply(pred.man, c(1,2), HDInterval::hdi, na.rm=T)[1,,])^12
uhdi95 <- plogis(apply(pred.man, c(1,2), HDInterval::hdi, na.rm=T)[2,,])^12

df <- data.frame(as.data.frame.table(md, responseName = "value"), 
                 Mean= as.data.frame.table(mn, responseName = "value")[,3],
                 lhdi95=as.data.frame.table(lhdi95, responseName = "value")[,3] |> round(2), 
                 uhdi95=as.data.frame.table(uhdi95, responseName = "value")[,3] |> round(2))
df <- df[order(df$Var1, df$Var2),]
colnames(df)[1:3] <- c("Age.class", "Period", "Median")
df$Median <- round(df$Median, 2)
df$Mean <- round(df$Mean, 2)

# get combined estimates from reduced model
df2 <- data.frame("Age.class"= c("First year", "Subadult", "Adult"),
                  "Period"= rep("Combined", 3),
                  "Median"= apply(p3$mean.s^12, 1, median) |> round(2),
                  "Mean"= apply(p3$mean.s^12, 1, mean) |> round(2),
                  "lhdi95"=apply(p3$mean.s^12, 1, HDInterval::hdi, credMass=0.95)[1,] |> round(2), 
                  "uhdi95"=apply(p3$mean.s^12, 1, HDInterval::hdi, credMass=0.95)[2,] |> round(2)
)
df3 <- rbind(df, df2)
df3 <- df3[order(df3$Age.class, df3$Period), ]
df3
# write.csv(file= "docs\\Survival_age_period.csv",
#           df3  )

# Generate probability of direction
# for comparing survival of age classes
# First we calculate the posterior survival differences
# Then we calculate the median and HDIs of survival differences
# Then we calculate PDs for age classes
s.diffs <- list()
s.diffs[[1]] <- p3$mean.s[3,] - p3$mean.s[2,]
s.diffs[[2]] <- p3$mean.s[3,] - p3$mean.s[1,]
s.diffs[[3]] <- p3$mean.s[2,] - p3$mean.s[1,]
surv.diffs <- data.frame( comparison = c("adults and subadults", "adults and first years", "subadults and first years"),
                          pd = lapply(s.diffs, function(x) { pd(x)$pd} ) |> unlist() |> round(2), 
                          median.diff = lapply(s.diffs, median ) |> unlist() |> round(3),
                          lapply(s.diffs, HDInterval::hdi, credMass=0.95) |> do.call(what=rbind) |> round(3) ) 
surv.diffs





