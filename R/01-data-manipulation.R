###############
# Data manipulation for discrete-time CJS 
# Survival models
################
library ("readxl")
library ("lubridate")
library ("data.table")
library ("tidyverse")
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
dat1 <- dat1[dat1$Dataset!="Ogada",]
dat1 <- dat1[dat1$Species %in% c("RUVU", "WBV") , ]
# print dataset for appendix
fields <- c("UnitID", "Dataset", "Species", "Stage", "WhereTrapped", 
             "DateAdded", "DateOfLoss-last activity", 
            "DaysWithTransmitter","include", "rehabbed",
             "taglost", "founddead", "CauseOfMortality",
            "UncertainFateTag")

toprint <- dat1[, fields]
write.csv(toprint,
 file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\MunirVultures\\docs\\Appendix1.csv")
dat1 <- dat1[dat1$include==T, ]


# Setup time bins for survival matrices
# Monthly setup
smonth <- list()
years <- 2009:2025
for (i in 1:12){ 
  smonth[[i]] <- c(paste0(i,"/1/", years)) 
}
s1 <- sort(mdy(do.call(c, smonth)))
mlength1 <- c(31,28,31,30,
              31,30,31,31,
              30,31,30,31)
emonth <- list()
dy <- ifelse(leap_year(years), 29, 28)
for (i in 1:12){
  if(i==2){
  emonth[[i]] <- c(paste0(i,"/",dy,"/", years))
  } else{
    emonth[[i]] <- c(paste0(i,"/",mlength1[i],"/", years))
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
  if(dat1$UncertainFateTag[i]==T){
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
# too few LFV during the later period
# subset to WBV
table(df$species, df$managed)
table(df$species, df$dead, df$managed)
df <- df[df$species %in% c("WBV", "RUVU"),]
table(df$managed)

yr <- as.numeric(substr(colnames(ch), 2, 5))
yr.cont <- (yr-median(2009:2024))/7.5 # for continuous covariate in survival
yr.factor <- as.numeric(factor(yr))
# assign a 1 for known age, and 2 for unknown age subadults.
# we're assuming survival is constant after 6 years so no need to track 
# age after 6
known1 <- as.numeric(!is.na(as.numeric(dat1$Age2)))
known <- ifelse(known1==1, 1, 2)

# Calculate tag age
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
# Remove individuals with <=2 occasions
# software can't run, 4 individs
tooshort <- c()
for (i in 1:nrow(ch)){
  tooshort[i] <- ((last2[i]-1)-(as.numeric(f[i])+1)) 
}
wtooshort <- which(tooshort<0)

ch[is.na(ch)] <- 5

datl <- list(
  y=ch[-wtooshort,], # observation matrix
  first_age=(as.numeric(dat1$Age2)+1)[-wtooshort] # data matrix of known first ages, but NAs for unknown ages
  )

constl <- list(
  f=as.numeric(f)[-wtooshort],# time interval of first capture
  last=last2[-wtooshort], # last time observed prior to censured length nind. 
  period.cat=ifelse(yr.cont<0,0,1),
  year.cont=yr.cont,
  rehabbed=ifelse(dat1$rehabbed==T, 1, 0)[-wtooshort],
  nind= nrow(datl$y),
  ntime= ncol(datl$y),
  nyears= length(unique(yr.factor)), 
  tag.age.sc= tag.age.sc[-wtooshort,], 
  sp=ifelse(dat1$Species=="WBV", 0, 1)[-wtooshort],
  known = ifelse(is.na(as.numeric(dat1$Age2)), 2, 1)[-wtooshort],
  y.first=first.val[-wtooshort]
)

save(datl=datl, constl, get.last=get.last2, get.first=get.first, 
     file="data\\data.RData")

# Export list of individuals 
# included in survival data
# for GIS map.
write.csv(file= "C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Gyps Vulture Survival in Africa\\docs\\individIDs_for_Map.csv",
          rownames(datl$y) )

# data summaries
# Setup data for calculating number of vultures
# and observation of each species during each period
y2 <- array(NA, dim=dim(datl$y), dimnames=dimnames(datl$y))
for(i in 1:constl$nind){
  for (t in constl$f[i]:constl$last[i]){
    y2[i,t] <- datl$y[i,t]
  }} 
ly <- y2 %>%
  as_tibble(rownames = "id") %>% # Convert to tibble, moving rownames to a column named "Var1"
  pivot_longer(
    cols = starts_with("d"),     # Specify columns to pivot (all starting with "col")
    names_to = "year_month",             # Name the new column for variable names "Var2"
    values_to = "state"            # Name the new column for values "value"
  )

species <- ifelse(constl$sp==0, "W", "R")
names(species) <- rownames(datl$y)  
period <- ifelse(yr.cont<0,"Early","Late")
names(period) <- colnames(datl$y) 
ly <- merge(ly, species, by.x="id", by.y=0)
ly <- merge(ly, period, by.x="year_month", by.y=0)
colnames(ly)[4:5] <- c("species","period")
# calculate ages
age <- array(999, dim=dim(datl$y), dimnames=dimnames(datl$y))
for(i in 1:constl$nind){
  age[i,constl$f[i]] <- datl$first_age[i]-1
  for (t in (constl$f[i]+1):constl$last[i]){
    age[i,t] <- floor(( (datl$first_age[i]-1) + t/12 - constl$f[i]/12 ))
  }} 
lage <- age %>%
  as_tibble(rownames = "id") %>% # Convert to tibble, moving rownames to a column named "Var1"
  pivot_longer(
    cols = starts_with("d"),     # Specify columns to pivot (all starting with "col")
    names_to = "year_month",             # Name the new column for variable names "Var2"
    values_to = "age"            # Name the new column for values "value"
  )
lage <- lage[lage$age<999,]

ly <- ly[ !is.na(ly$state), ]
ly <- merge(ly, lage, by=c("id", "year_month"), all.x=T)
ly$agec <- ifelse(is.na(ly$age), "SA-Unknown",
            ifelse(ly$age<1, "FY",
              ifelse(ly$age>=1 & ly$age<6, "SA-Known",
                ifelse(ly$age>=6, "A", NA  ))))
rh <- ifelse(constl$rehabbed==1, "rehabbed", "no rehab")
names(rh) <- rownames(datl$y)
ly <- merge(ly, rh, by.x="id", by.y=0, all.x=T)
colnames(ly)[8] <- "rehabbed"
# these sum to <1393 because unknown ages
# number of monthly observations

# number of individuals
li <- ly[!duplicated(ly$id),]

table(datl$y) ; sum(table(datl$y)[1:3])

table(ly$agec)
table(ly$period)
table(ly$species, ly$period)
table(ly$species, ly$period, ly$rehabbed)
table(ly$species, ly$period, ly$agec)
table(ly$species, ly$period, is.na(ly$age))


table(li$period)
table(li$species, li$period)
table(li$species, li$period, li$rehabbed)
table(li$species, li$period, li$agec)
table(li$species, li$period, is.na(li$age))

write.csv(rownames(datl$y), 
          file="docs/individIDs_for_Map.csv")

tdat <- read_xlsx("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\MunirVultures\\data\\VultureMortalityData_transmitter_harness _RBMV.xlsx")
tdat <- tdat[tdat$UnitID %in% rownames(datl$y),]
tdat$yr <- year(tdat$DateAdded)
tdat$period <-  ifelse(tdat$yr %in% c(2009, 2010, 2011),
                       "Early", "Late")

table(tdat$`Harness type1`)
table(tdat$`Transmitter type1`, tdat$period)
