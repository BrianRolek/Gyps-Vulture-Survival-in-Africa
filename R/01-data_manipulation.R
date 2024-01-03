###############
# Data manipulation for discrete-time CJS 
# Survival models
################
library ("readxl")
library ("lubridate")
library ("data.table")
# data manip
dat1 <- read_xlsx("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\MunirVultures\\data\\VultureMortalityData RB LD.xlsx",
                  sheet="CombinedData",
                  col_types = c(rep("guess",11),
                                "date", "date", 
                                "guess", "date", "guess",
                                "logical", "logical", "logical", 
                                "logical","logical", 
                                rep("guess",4) 
                                ))
dat1 <- dat1[dat1$include==T, ]
dat1 <- dat1[dat1$Species %in% c("RUVU", "WBV") , ]
# Setup time bins for survival matrices
# seasonal dates based on equinoxes 
# "long rains" March-May, "short-rains" Oct-Dec
# Chamberlin and Wairoto 1997
# this pub cites Anyaba 1983 and Ogallo 1988
# Wet Dry season setup
# s1 <- sort(mdy(c(paste0(c("1/1/"), 2009:2023), 
#         paste0(c("3/1/"), 2009:2023),
#         paste0(c("6/1/"), 2009:2023),
#         paste0(c("10/1/"), 2009:2023))))
# e1 <- sort(mdy( c(ifelse(leap_year(2009:2023), 
#              paste0(c("2/29/"), 2009:2023),
#              paste0(c("2/28/"), 2009:2023)), # account for leap years
#              paste0(c("5/31/"), 2009:2023),
#              paste0(c("9/31/"), 2009:2023),
#              paste0(c("12/31/"), 2009:2023)) ))

# Monthly setup
smonth <- list()
for (i in 1:12){ 
  smonth[[i]] <- c(paste0(i,"/1/", 2009:2023)) 
}
s1 <- sort(mdy(do.call(c, smonth)))
mlength1 <- c(31,28,31,30,
              31,30,31,31,
              30,31,30,31)
emonth <- list()
dy <- ifelse(leap_year(2009:2023), 29, 28)
for (i in 1:12){
  if(i==2){
  emonth[[i]] <- c(paste0(i,"/",dy,"/", 2009:2023))
  } else{
    emonth[[i]] <- c(paste0(i,"/",mlength1[i],"/", 2009:2023))
  }
}
e1 <- sort(mdy(do.call(c, emonth)))

labels <- c(paste0("d", year(s1), "_", month(s1)))
dtbins <- data.frame(
  start = s1, # cutoffs between pre and post-hatching
  end = e1,
  labs=labels
  )
mndt <- round_date(min(ymd(dat1$DateAdded), na.rm=T), unit="months")
mxdt <- max(ymd(dat1$`DateOfLoss-last activity`), na.rm=T)
dtbins <- dtbins[dtbins$start >= mndt & dtbins$end <= mxdt,]
# dtbins$nms <-  paste0(c("D", "W"), year(dtbins$end))
dtbins$yr_int <- dtbins$start %--% dtbins$end 

int_overlaps_numeric <- function (int1, int2) {
  stopifnot(c(is.interval(int1), is.interval(int2)))
  x <- intersect(int1, int2)@.Data
  x[is.na(x)] <- 0
  as.duration(x)
}
# create survival intervals from data frame
# using banding through mortality or censor date
added <- ymd(dat1$DateAdded)
dead <- ymd(dat1$DateOfMortality)
dcens <- as_date(ymd(as.character(dat1$`DateOfLoss-last activity`)))
# sub in the max date of study if bird is not dead
end1 <- fifelse(is.na(dcens), dead, dcens )
end2 <- fifelse(is.na(end1), 
             max(dtbins$end, na.rm=T), 
             end1)
ntime <- nrow(dtbins)
nind <- nrow(dat1)
ed <- ad <- array(NA, dim=c(nind, ntime), dimnames=list(dat1$UnitID, dtbins$labs) )
# Create intervals for exposure
int <-  added %--% end2
for (t in 1:ntime){
  # total number of possible exposure days, accounts for hatch day in mid-year
  startd <- fifelse(dtbins$start[t]>added,
                   dtbins$start[t], added )
  ed[,t] <- int_overlaps_numeric(int, dtbins$yr_int[t]) / ddays(1) # year length in days
}
ad <- ifelse(ed>0, 1, 0)
colSums(ad)
rowSums(ad)

#************
#* Construct observation matrix
#************
ch <- array(NA, dim=dim(ad), dimnames=dimnames(ad))
ch[] <- ifelse(ad==1, 1, NA)
get.last <- function(x) max(which(x==1), na.rm=T)
last <- apply(ch, 1, get.last)
last.val <- with(dat1, ifelse(taglost==F & founddead==F & is.na(DateOfMortality), 1,
                  ifelse(taglost==T & founddead==F & is.na(DateOfMortality), 2,
                         ifelse(taglost==F & founddead==T & !is.na(DateOfMortality), 3,
                                ifelse(taglost==T & founddead==T & !is.na(DateOfMortality), 4, 
                                     5)))))
# replace last value with fate
last.val2 <- c()
for (i in 1:nrow(ch)){
  ch[i,last[i] ] <- last.val[i] 
  if(dat1$UncertainFate[i]==T){
    ch[i,last[i] ] <- 5
  }
  last.val2[i] <- ch[i,last[i] ]
}
get.last2 <- function(x) max(which(!is.na(x)), na.rm=T)
last2 <- apply(ch, 1, get.last2)

get.first <- function(x) min(which(x %in% c(1:4)), na.rm=T)
f <- apply(ch, 1, get.first)
first.val <- c(NA)
for (i in 1:nrow(ch)){
  first.val[i] <- ch[i,f[i] ]
}

live.seq <- list()
for (i in 1:nrow(ch)){
  live.seq[[i]] <- c(ch[i, f[i]:last2[i]])
}


#***************************
#* Calculate age of birds and 
#* Age of tags
#***************************
bird.age <- array(NA, dim=dim(ch), dimnames=dimnames(ch))
for (i in 1:nrow(ch)){
  bird.age[i, f[i] ] <- as.numeric(dat1$Age2[i]) 
}



df <- data.frame(time=dat1$DaysWithTransmitter, 
                 censored=ifelse(dat1$Censored=="Y", T, F), 
                 managed= ifelse(dat1$Manage=="Yes", T, F),
                 dead= ifelse(is.na(dat1$DateOfMortality), F, T),
                 species=as.factor(dat1$Species ), 
                 gyps=as.factor(ifelse(dat1$Species %in% c("RUVU", "WBV"), "gyps", "non-gyps") ) ) 
df$event <- ifelse(df$censored==T, 0, 1)
# too few LFV and RUVU during the managed period
# subset to WBV
table(df$species, df$managed)
table(df$species, df$dead, df$managed)
df <- df[df$species %in% c("WBV", "RUVU"),]
table(df$managed)

yr <- as.numeric(substr(colnames(ch), 2, 5))
yr.cont <- (yr-2016)/7 # for continuous covariate in survival
yr.factor <- as.numeric(factor(yr))
# assign a 1 for known age, and 2 for unknown age subadults.
# we're assuming survival is constant after 6 years so no need to track 
# age after 6
known1 <- as.numeric(!is.na(as.numeric(dat1$Age2)))
known <- ifelse(known1==1, 1, 2)
# z values for directly observed states
zmat <- ifelse(ch<=4, ch, NA)
first_zs <- c()
for (i in 1:nrow(zmat)){
  first_zs[i] <- zmat[i,f[i]] 
  zmat[i, f[i] ] <- NA # sub in NA for first capture
}

for (i in 1:nrow(zmat)){
if (last.val2[i]==3){
  zmat[i, (last2[i]+1) :ntime  ] <- 5
}}

lastz <- apply(zmat, 1, get.last2)

# z initial values for states that weren't observed
z.inits <- array(NA, dim(zmat), dimnames(zmat))
for (i in 1:nrow(z.inits)){ 
  indz1  <- ifelse(lastz[i]==ncol(z.inits), ncol(z.inits), lastz[i]+1)
  indz2  <- ifelse(lastz[i]==ncol(z.inits), ncol(z.inits), lastz[i]+2)
  ind1  <- ifelse(last2[i]==ncol(z.inits), ncol(z.inits), last2[i]+1)
  ind2 <- ifelse(last2[i]==ncol(z.inits), ncol(z.inits), last2[i]+2)
  if ( ch[i, last2[i] ] == 1 ) {
  z.inits[i, ind1] <- 3
  z.inits[i, (ind2):ncol(z.inits)] <- 5
  #next
  } # zlast=1
if ( ch[i, last2[i] ] == 2 ) {
  z.inits[i, (ind1)] <- 4
  z.inits[i, (ind2):ncol(z.inits)] <- 5
  #next
} # zlast=2
 if ( ch[i, last2[i] ] %in% c(3,4) ) {
   z.inits[i, (ind1):ncol(z.inits)] <- 5
   #next
 } # zlast=3 or 4
  if ( ch[i, last2[i] ] == 5 ) {
    z.inits[i, last2[i] ] <- 3
    z.inits[i, (ind1):ncol(z.inits)] <- 5
    #next
  } # zlast=5
}
z.inits[!is.na(zmat)] <- NA


tag.age <- array(NA, dim(ch), dimnames=dimnames(ch))
len <- sums <- c()
for (i in 1:nind){
  tag.age[i,f[i]] <- 0
for (t in f[i]:(ntime-1) ){    
  tag.age[i,t+1] <- tag.age[i,t] + 1/12
} }
tag.age.sc <- (tag.age-mean(tag.age, na.rm=T))/sd(tag.age, na.rm=T)
tag.age.sc[is.na(tag.age.sc) ] <- 0

man.cat <- array(0, dim(ch), dimnames=dimnames(ch))
man.cat[,31:ntime] <- 1

datl <- list(
  y=ch, # observation matrix
  z=zmat, # known state matrix
  z.inits = z.inits,
  f=as.numeric(f),# time interval of first capture
  k=last2, # last time observed prior to censured length nind. 
  known= known, # whether age was known or not, subadults are unknown
  first_age=as.numeric(dat1$Age2), # data matrix of known first ages, but NAs for unknown ages
  managed.cat=man.cat,
  year.cont=yr.cont,
  year.factor=yr.factor,
  rehabbed=ifelse(dat1$rehabbed==T, 1, 0),
  nind= nrow(ch),
  ntime= ncol(ch),
  nyears= length(unique(yr.factor)), 
  first_zs= first_zs, 
  tag.age.sc= tag.age.sc, 
  sp=ifelse(dat1$Species=="WBV", 0, 1),
  region= ifelse(dat1$Region=="N", 1, 0),
  study=as.numeric(as.factor(dat1$Dataset))
  )

save(datl=datl, get.last=get.last2, get.first=get.first, 
     file="data\\data.RData")
