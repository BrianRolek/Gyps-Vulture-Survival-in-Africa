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
                  col_types = c(rep("guess",10),
                                "date", "date", 
                                "guess", "date", "guess", "guess",
                                "guess", "date", "date",
                                rep("guess",4) 
                                ))
dat1 <- dat1[dat1$include=='T', ]

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
mndt <- min(ymd(dat1$DateAdded), na.rm=T)
mxdt <- max(ymd(dat1$`DateOfLoss-last activity`), na.rm=T)
# ignore that row 63 failed to parse because it is ymd_hms() format
# not the max value so irrelevant
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
#dcens[63] <- as_date(ymd_hms(as.character(dat1$`DateOfLoss-last activity`[63])))
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
  # calculate overlap length of intervals in days
  ad[, t] <- int_overlaps_numeric(int[t], dtbins$yr_int[t]) / ddays(1)
  # total number of possible exposure days, accounts for hatch day in mid-year
  startd <- fifelse(dtbins$start[t]>added,
                   dtbins$start[t], added )
  ed[,t] <- int_overlaps_numeric(int, dtbins$yr_int[t]) / ddays(1) # year length in days
}
ad <- ifelse(ed>15, 1, 0)
colSums(ad)

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

datal <- list(
  ch=, # survival matrix
  f=, # time interval of first capture
  k= , # last time observed prior to censured length nind. 
  known=, # whether age was known or not, subadults are unknown
  first_age=, # data matrix of known first ages, but NAs for unknown ages
  nind= nrow(ch)
  )