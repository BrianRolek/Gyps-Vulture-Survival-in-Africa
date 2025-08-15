library(MCMCvis)
library (coda)
library (ggplot2)
library (reshape2)
library (tidybayes)
library (bayestestR)
library (ggpubr)
options(scipen=999)
load("data//data.RData")
pars <- c(  "delta", "beta",#"gamma", "eta",
            "mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead",
            "l.s", "l.tagfail", "l.p.tagfail", "l.p.dead")
# Export list of individuals 
# included in survival data
# Allows Leah to make a map.
write.csv(file= "C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Gyps Vulture Survival in Africa\\docs\\individIDs_for_Map.csv",
          rownames(datl$y) )
# load output from global model
load("outputs\\gyps-11Aug2025-marginalized.RData")
post.global <- post
p <- MCMCpstr(post, pars, type="chains")
p2 <- mcmc.list(post)
# load output from reduced model
load("outputs\\gyps-12Aug2025-marginalized-reduced.RData")
post.reduced <- post
p3 <- MCMCpstr(post, pars[-1], type="chains")
# Calculate age sample sizes
first <- last <- c()
ages <- array(NA, dim=dim(datl$y))
for (i in 1:datl$nind){
  y.min <- which.min( datl$y[i,] )
  y.max <- which.max( datl$y[i,] )
  z.min <- which.min( datl$z[i,] )
  z.max <- which.max( datl$z[i,] )
  first[i] <- min(c(y.min, z.min))
  last[i] <- min(c(y.max, z.max))
  for (t in 1:datl$ntime){
    
ages[i, t] <- ( datl$first_age[i] + t/12 - datl$f[i]/12 )  
}}
ragged.melt <- function(x, first, last){
  val.list <- list()
  for (i in 1:nrow(x)){
    val.list[[i]] <- x[i, first[i]:last[i]]
  } # i
  all.vals <- do.call(c, val.list)
  return( all.vals )
} # function end
age.tab <- table(ragged.melt(floor(ages), first, last))
age.tab[1] # first year bird months
sum(age.tab[2:6]) # subadult bird months
sum(age.tab[7:length(age.tab)]) # adult bird months

# Model diagnostics
# Check for convergence
iters <- ncol(p$delta)
MCMCtrace(post.global, "delta", pdf=F, Rhat=T, 
          priors=rnorm(iters, 0, 10), post_zm = FALSE)
MCMCtrace(post.global, "beta", pdf=F, Rhat=T, 
          priors=rnorm(iters, 0, 10), post_zm = FALSE)       
MCMCtrace(post.global, c("mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead"), pdf=F, Rhat=T, 
          priors=rbeta(iters, 1, 1), post_zm = FALSE)  

iters <- ncol(p3$beta)
MCMCtrace(post.reduced, "beta", pdf=F, Rhat=T, 
          priors=rnorm(iters, 0, 10), post_zm = FALSE)       
MCMCtrace(post.reduced, c("mean.s", "mean.tagfail", "mean.p.tagfail", "mean.p.dead"), pdf=F, Rhat=T, 
          priors=rbeta(iters, 1, 1), post_zm = FALSE) 
# Table 1
# Summary stats
# Sample Sizes
ly <- melt(datl$y)
sp.df <- data.frame(Var1=rownames(datl$y), species=datl$sp)
ly2 <- merge(ly, sp.df, by="Var1")
ly2 <- ly2[!is.na(ly2$value),]
colSums(table(ly2$value, ly2$species))
# table 1.
t(table(ly2$value, ly2$species))
colSums(t(table(ly2$value, ly2$species)))
rowSums(t(table(ly2$value, ly2$species)))

# Table 2 
# Estimates
# Survival by management period
sum95.global <- MCMCsummary(post.global, pars[-c(7:10)], HPD=TRUE, digits=2, 
            hpd_prob=0.95, pg0=TRUE, func=median, func_name="md")

coef.est.global <- data.frame(Parameter= rownames(sum95.global),
                       Median=sum95.global$md, 
                       Mean=sum95.global$mean,
                       LHDI95=sum95.global$`95%_HPDL`, 
                       UHDI95=sum95.global$`95%_HPDU`,
                       p= sum95.global$`p>0`, 
                       Rhat=sum95.global$Rhat
                       )
coef.est.global <- cbind(Model="Global", coef.est.global)

sum95.reduced <- MCMCsummary(post.reduced, pars[-c(1,7:10)], HPD=TRUE, digits=2, 
                     hpd_prob=0.95, pg0=TRUE, func=median, func_name="md")
coef.est.reduced <- data.frame(Parameter= rownames(sum95.reduced),
                       Median=sum95.reduced$md, 
                       Mean=sum95.reduced$mean,
                       LHDI95=sum95.reduced$`95%_HPDL`, 
                       UHDI95=sum95.reduced$`95%_HPDU`,
                       p= sum95.reduced$`p>0`, 
                       Rhat=sum95.reduced$Rhat
)
coef.est.reduced <- cbind(Model="Reduced", coef.est.reduced)

coef.est <- rbind(coef.est.global, coef.est.reduced)
write.csv(file= "docs\\Coef_ests.csv",
          coef.est)

par(mfrow=c(1,1))
MCMCplot(p2, params=c("delta", 
                      "beta"), 
         ISB=TRUE, 
         HPD=TRUE, ci=c(85, 95))
MCMCplot(p2, params=c("mean.s", "mean.tagfail", 
                      "mean.p.tagfail", "mean.p.dead"), 
         ISB=TRUE, 
         HPD=TRUE, ci=c(85, 95))

#***************
#* Plot effect from tag age on tag failure
#* from reduced model
#***************
ta <- seq(0,5.3, by=0.1)
ta.sc <- (ta-5.436775)/3.997963
ni <- ncol(p3$beta)
pred.ta <- array(NA, dim=c(length(ta.sc), ni), dimnames=list(ta, 1:ni) )
for (i in 1:length(ta.sc)){
pred.ta[i,] <- p3$l.tagfail[1,] + p3$beta[1,]*ta.sc[i] #+ p$beta[2,]*ta.sc[i]^2
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

ggsave("figs\\transmitter age and failure.tiff",
       pyta, device="tiff",
       width=6.5, height=4, units="in", dpi=300)

ggsave("figs\\transmitter age and failure-byMonth.tiff",
       pmta, device="tiff",
       width=6.5, height=4, units="in", dpi=300)

df2 <- data.frame(md =1-((1-md)^12),
                  lhdi95= 1-((1-lhdi95)^12),
                  uhdi95=1-((1-uhdi95)^12)
)

#****************
#* plot survival by age class
#****************
# from reduced model, p3
lss <- melt(p3$mean.s)
#lss$Ageclass <- ifelse(lss$Var1 == "mean.s[1]", "First year",
#                       ifelse(lss$Var1 == "mean.s[2]"), "Subadult", "Adult")

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

# Calculate yearly survival for results
lss$yr.s<- lss$value^12
tapply(lss$yr.s, lss$Var1, median)
tapply(lss$yr.s, lss$Var1, mean)
tapply(lss$yr.s, lss$Var1, HDInterval::hdi, credMass=0.95)
tapply(lss$yr.s, lss$Var1, HDInterval::hdi, credMass=0.85)

ggsave("figs\\survival-ageclass.tiff",
       ps3, device="tiff", 
       width=6.5, height=4, units="in", dpi=300)

# Calculate PDs for age classes
s.diffs <- list()
s.diffs[[2]] <- p3$mean.s[3,] - p3$mean.s[2,]
s.diffs[[1]] <- p3$mean.s[3,] - p3$mean.s[1,]
s.diffs[[3]] <- p3$mean.s[2,] - p3$mean.s[1,]
lapply(s.diffs, pd)

#****************
#* plot survival in response to management/time
#****************
ni <- ncol(p$delta)
pred.man <- array(NA, dim=c(3, 2, ni), 
                 dimnames=list(c("First year", "Subadult", "Adult"), c("Early", "Late"), 1:ni) )
for (a in 1:3){
for (m in 1:2){
  pred.man[a,m,] <- p$l.s[a,] + 
                    p$delta[1,]*c(0,1)[m]
}}
lp.man <- melt(pred.man)
colnames(lp.man)[1:3] <- c("Ageclass", "Period", "iter" )
lp.man$pred <- plogis(lp.man$value)

p4 <- ggplot(data=lp.man, aes(x=pred^12, y=Period)) + theme_minimal() +
  geom_line(data=lp.man, aes(x=pred^12, y=Period, group=iter),
            color="gray40", linewidth=0.5, alpha=0.025) +
  stat_pointinterval(.width=c(0.85, 0.95), 
                     point_interval ="median_hdci") +
  # stat_halfeye(lss, aes(x = value^12, y = Var1),
  #              .width=c(0.85, 0.95), point_interval="median_hdci") +
  facet_wrap(facets=vars(Ageclass)) + 
  xlim(0,1) +
  coord_flip() +
  ylab("Period") + xlab("Survival (yearly probability)") +
  ggtitle("(A) Management")

ggsave("figs\\survival-ageclass-management.tiff",
       p4, device="tiff", 
       width=6.5, height=4, units="in", dpi=300)

all_p <- ggarrange(p4, ps3, nrow=2)
ggsave("figs\\survival-ageclass-management-combined.tiff",
       all_p, device="tiff", 
       width=6, height=6, units="in", dpi=300)

# Table 3 
# survival by age class
# and managment period
md <- plogis(apply(pred.man, c(1,2), median, na.rm=T))^12
mn <- plogis(apply(pred.man, c(1,2), mean, na.rm=T))^12
lhdi95 <- plogis(apply(pred.man, c(1,2), HDInterval::hdi, na.rm=T)[1,,])^12
uhdi95 <- plogis(apply(pred.man, c(1,2), HDInterval::hdi, na.rm=T)[2,,])^12
lhdi85 <- plogis(apply(pred.man, c(1,2), HDInterval::hdi, na.rm=T, credMass=0.85)[1,,])^12
uhdi85 <- plogis(apply(pred.man, c(1,2), HDInterval::hdi, na.rm=T, credMass=0.85)[2,,])^12

df <- data.frame(melt(md), 
                 Mean= melt(mn)[,3],
                 lhdi85=melt(lhdi85)[,3] |> round(2), 
                 uhdi85=melt(uhdi85)[,3] |> round(2),
                 lhdi95=melt(lhdi95)[,3] |> round(2), 
                 uhdi95=melt(uhdi95)[,3] |> round(2))
df <- df[order(df$Var1, df$Var2),]
colnames(df)[1:3] <- c("Age.class", "Period", "Median")
df$Median <- round(df$Median, 2)
df$Mean <- round(df$Mean, 2)

# get combined estimates from reduced model
df2 <- data.frame("Age.class"= c("First year", "Subadult", "Adult"),
                  "Period"= rep("Combined", 3),
                  "Median"= apply(p3$mean.s^12, 1, median) |> round(2),
                  "Mean"= apply(p3$mean.s^12, 1, mean) |> round(2),
                  "lhdi85"=apply(p3$mean.s^12, 1, HDInterval::hdi, credMass=0.85)[1,] |> round(2), 
                  "uhdi85"=apply(p3$mean.s^12, 1, HDInterval::hdi, credMass=0.85)[2,] |> round(2),
                  "lhdi95"=apply(p3$mean.s^12, 1, HDInterval::hdi, credMass=0.95)[1,] |> round(2), 
                  "uhdi95"=apply(p3$mean.s^12, 1, HDInterval::hdi, credMass=0.95)[2,] |> round(2)
)
df3 <- rbind(df, df2)
df3 <- df3[order(df3$Age.class, df3$Period), ]
write.csv(file= "docs\\Survival_age_period.csv",
          df3  )


#***************
#* Plot effect from year of study on tag failure
#***************
yr <- seq(2009,2025, by=0.1)
yr.sc <- (yr-2016)/7
ni <- 4000
pred.tf.yr <- array(NA, dim=c(length(yr.sc), ni), dimnames=list(yr, 1:ni) )
for (i in 1:length(yr.sc)){
  pred.tf.yr[i,] <- p3$l.tagfail[1,] + p3$beta[2,]*yr.sc[i] #+ p$beta[4,]*yr.sc[i]^2
}
lp.tf.yr <- melt(pred.tf.yr)
colnames(lp.tf.yr)[1:2] <- c("year", "iter" )
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
  geom_line(data=df, aes(x=yr, y=1-((1-md)^12)), linewidth=2) +
  geom_line(data=df, aes(x=yr, y=1-((1-lhdi85)^12)), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=1-((1-uhdi85)^12)), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=1-((1-lhdi95)^12)), linewidth=1, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=1-((1-uhdi95)^12)), linewidth=1, linetype="dashed") +
  ylab("Transmitter failure (yearly probability)") + xlab("Year of study") +
  facet_wrap("Period", scales="free_x")

pmty <- ggplot() + theme_minimal() +
  geom_line(data=lp.tf.yr, aes(x=year, y=pred, group=iter),
            color="gray40", linewidth=0.5, alpha=0.05) +
  geom_line(data=df, aes(x=yr, y=md), linewidth=2) +
  geom_line(data=df, aes(x=yr, y=lhdi85), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=uhdi85), linewidth=2, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=lhdi95), linewidth=1, linetype="dashed") +
  geom_line(data=df, aes(x=yr, y=uhdi95), linewidth=1, linetype="dashed") +
  ylab("Transmitter failure (monthly probability)") + xlab("Year of study") +
  facet_wrap("Period", scales="free_x")

ggsave("figs\\tagfailure-year.tiff",
       pyty, device="tiff", 
       width=6.5, height=4, units="in", dpi=300)

ggsave("figs\\tagfailure-month.tiff",
       pmty, device="tiff", 
       width=6.5, height=4, units="in", dpi=300)

df2 <- data.frame(md =1-((1-md)^12), 
                  lhdi95= 1-((1-lhdi95)^12),
                  uhdi95=1-((1-uhdi95)^12)    
)


