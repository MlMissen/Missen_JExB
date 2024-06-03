# Missen_JExB

setwd("Desktop")

library(ggplot2)
library(dplyr)
library(tidyverse)
library(MASS)
library(car)
library(emmeans)
library(lme4)
library(rcompanion)
library(MuMIn)

#TasFACE2 SWC

SWC.freq <- read.csv("soil_water_content.csv", na.strings ="na")

SWC.freq$Month <- factor(SWC.freq$Month, levels= c("Jun-20","Jul-20", "Aug-20", "Sep-20", "Oct-20", "Nov-20", "Dec-20", "Jan-21", "Feb-21", "Mar-21", "Apr-21", "May-21", "Jun-21","Jul-21", "Aug-21"))
SWC.freq$CO2 <- factor(SWC.freq$CO2, levels=c("Elevated","Ambient"))
SWC.freq$Ring<- factor(SWC.freq$Ring)
SWC.freq$Plot_ID<- factor(SWC.freq$Plot_ID)
SWC.freq$Water <- factor(SWC.freq$Water, levels=c("3 days","5 days", "10 days"))
SWC.freq$Date <- as.Date(SWC.freq$Date , format = "%d/%m/%Y", tz = "Australia/Tasmania")


# Calculate SWC cycle means
SWC.ALL.CYCLE.MEAN <- SWC.freq %>%
  group_by(Plot_ID, Cycle) %>%
  dplyr::mutate(SWC.CYCLE= mean(SWC.MEAN)) %>%
  ungroup()

# delete duplicate rows
SWC.ALL.CYCLE.MEAN <- SWC.ALL.CYCLE.MEAN %>%
  distinct(Plot_ID, Cycle, .keep_all = TRUE)

# CALCULATE MONTHLY MEAN FROM CYCLE MEAN

SWC.ALL.CYCLE.MM <- SWC.ALL.CYCLE.MEAN %>%
  group_by(Plot_ID, Month, Ring, CO2, Water) %>%
  dplyr::summarise(SWC.MM= mean(SWC.CYCLE))


exclude_months <- c('Apr-21', 'May-21', 'Jun-21', 'Jul-21', 'Aug-21')
SWC.ALL.CYCLE.MM <- SWC.ALL.CYCLE.MM[!SWC.ALL.CYCLE.MM$Month %in% exclude_months, ]

# plot boxcox for SWC data and transform as necessary

SWC.NORMAL.MM.m <- aov(SWC.MM~CO2*Water*Month , data = SWC.ALL.CYCLE.MM)
boxcox(SWC.NORMAL.MM.m, lambda = seq(0,2, length=10)) #^4 transform

# run lmer on SWC data, testing for CO2, Water and Month effects

SWC.MM.lme1 <- lmer(log(SWC.MM)~CO2*Water*Month + (1|Ring/Plot_ID), SWC.ALL.CYCLE.MM)
Anova(SWC.MM.lme1, test.statistic = "F") 
emmeans(SWC.MM.lme1, "Water", "Month", pairwise = TRUE)

#####################################################################

### Relative SWC

# Get SWC average from June - Aug 2020, by plot_ID
SWC.ALL.CYCLE.MEAN.PRETREATMENT <- subset(SWC.ALL.CYCLE.MEAN, Month == "Jun-20" | Month == "Jul-20"|Month == "Aug-20")

SWC.NORMAL <<- SWC.ALL.CYCLE.MEAN.PRETREATMENT %>%
  group_by(Plot_ID) %>%
  dplyr::summarize(SWC.PRETREATMENT.AVERAGE = mean(SWC.CYCLE))%>%
  ungroup()

# SWC.AVERAGE OVER 3 MONTHS

# MARRY DATAFRAMES
SWC.NORMAL.ALL <-merge(x = SWC.ALL.CYCLE.MEAN, y = SWC.NORMAL, by = c("Plot_ID"))

SWC.NORMAL.ALL.1 <- SWC.NORMAL.ALL%>%
  mutate(NORMALISED.SWC=SWC.CYCLE/SWC.PRETREATMENT.AVERAGE)


SWC.NORMAL.ALL.MM <- SWC.NORMAL.ALL.1 %>%
  group_by(Plot_ID, Month, CO2, Water, Ring) %>%
  dplyr::summarize(NORMALISED.SWC.MM= mean(NORMALISED.SWC))

exclude_months <- c('Apr-21', 'May-21', 'Jun-21', 'Jul-21', 'Aug-21')
SWC.NORMAL.ALL.MM <- SWC.NORMAL.ALL.MM[!SWC.NORMAL.ALL.MM$Month %in% exclude_months, ]

# check distribution 
SWC.NORMAL.MM.m <- aov(NORMALISED.SWC.MM~CO2*Water*Month , data = SWC.NORMAL.ALL.MM)
boxcox(SWC.NORMAL.MM.m, lambda = seq(2,5, length=10)) #^4 transform

# run lmer on RSWC data, testing for CO2, Water and Month effects

SWC.NORMALISED.MM.lme1 <- lmer(NORMALISED.SWC.MM^4.5~CO2*Water*Month + (1|Ring/Plot_ID), SWC.NORMAL.ALL.MM)
Anova(SWC.NORMALISED.MM.lme1, test.statistic = "F") 



######################################################################

## Coefficient of variation of SWC

SWC.freq.SD.MEAN <- SWC.freq %>%
  group_by(Plot_ID, Month, Ring, CO2, Water) %>%
  dplyr::summarise(sd=sd(SWC.MEAN), SWC.MEAN.M=mean(SWC.MEAN))

SWC.freq.CV <- SWC.freq.SD.MEAN %>%
  group_by(Plot_ID, Month, Ring, CO2, Water) %>%
  dplyr::summarise(CV= sd/SWC.MEAN.M*100)

# plot boxcox for coefficient of variation of SWC and transform as necessary

CV.m <- aov(CV~CO2*Water*Month, data = SWC.freq.CV)
boxcox(CV.m, lambda = seq(0,1.5, length=10)) #log transform

# run lmer on coefficient of variation of SWC, testing for CO2, Water and Month effects

CV.lme <- lmer(log(CV)~CO2*Water*Month + (1|Ring/Plot_ID), SWC.freq.CV)
Anova(CV.lme, test.statistic = "F")


######################################################################

####  mean amplitude of change in SWC

SWC.freq.MAX.MIN <- SWC.freq %>%
  group_by(Plot_ID, Month, Ring, CO2,  Water, Cycle) %>%
  dplyr::summarise(MAX=max(SWC.MAX), MIN=min(SWC.MIN))

SWC.freq.MAX.MIN.CHANGE <- SWC.freq.MAX.MIN %>%
  group_by(Plot_ID, Month, Ring, CO2,  Water, Cycle) %>%
  dplyr::mutate(change=MAX - MIN)

SWC.freq.MAX.MIN.Percent <- SWC.freq.MAX.MIN.CHANGE %>%
  group_by(Plot_ID, Month, Ring, CO2,  Water, Cycle) %>%
  dplyr::mutate(percent=change/MAX*100)

# only include months in study period 

exclude_months <- c('Apr-21', 'May-21', 'Jun-21', 'Jul-21', 'Aug-21')
subsetted_df <- SWC.freq.MAX.MIN.Percent[!SWC.freq.MAX.MIN.Percent$Month %in% exclude_months, ]


# plot boxcox for amplitude of change of SWC and transform as necessary

PC.m <- aov(percent~CO2*Water*Month, data = subsetted_df)
boxcox(PC.m, lambda = seq(0,1.5, length=10)) #sqrt transform

# run lmer on amplitude of change of SWC, testing for CO2, Water and Month effects

PC.lme <- lmer(percent^0.5~CO2*Month*Water + (1|Ring/Plot_ID), subsetted_df)
Anova(PC.lme, test.statistic = "F")

PC.emm <- emmeans(PC.lme, ~ Month * Water)
pairs(PC.emm, by = "Month", adjust = "tukey")


#####################################################################
################### Gas exchange analysis data ######################

GE <- read.csv('gas_exchange.csv')

GE$Water <- factor(GE$Water, levels = c("3 days", "5 days", "10 days"))
GE$CO2 <- factor(GE$CO2, levels=c("elevated", "ambient"))
GE$Season <- factor(GE$Season, levels=c("Winter20", "Spring20", "Summer21", "Autumn21"))
GE$DSW <- factor(GE$DSW, levels=c("1", "2", "3", "4", "5", "6", "7", "9","10"))
GE$Plot_ID <- paste(GE$Ring,GE$Water, sep="_") #create Plot_ID

############## Photosynthesis #####################

# plot boxcox for Photosynthesis data and transform as necessary

GE.m <- aov(Anet~CO2*Water*Month, data = GE)
boxcox(GE.m, lambda = seq(0,1.5, length=10)) #sqrt transform

# run lmer on photosynthesis data, testing for CO2, Water and Month effects

GE.lme1 <- lmer(Anet~CO2*Month*Water + VpdL+(1|Ring/Plot_ID), GE)
Anova(GE.lme1, test.statistic = "F")
summary(GE.lme1)


# pairwise analysis - CO2 x Water
emm.photo <- emmeans(GE.lme1, ~ CO2|Month)
pairs(emm.photo, adjust = "tukey")

# pairwise analysis - Month x Water

emm.photo <- emmeans(GE.lme1, ~ Water|Month)
pairs(emm.photo, adjust = "tukey")


############ stomatal conductance ################

# plot boxcox for Gs data and transform as necessary
GE.Gs.m <- aov(Gs~CO2*Water*Month, data = GE)
boxcox(GE.Gs.m, lambda = seq(0,1.5, length=10)) #sqrt transform
plot(GE.Gs.m) 

# run lmer on Gs data, testing for CO2, Water and Month effects

GE.Gs.lme1 <- lmer(Gs^0.5~CO2*Month*Water + VpdL + (1|Ring/Plot_ID), GE)
Anova(GE.Gs.lme1, test.statistic = "F")

# pairwise analysis - CO2 x Month

emm.cond <- emmeans(GE.Gs.lme1, ~ CO2|Month)
pairs(emm.cond, adjust = "tukey")


###################################################################

# Days since water analysis

# plot boxcox for photosynthesis data and transform as necessary
photo.DSW.m <- aov(Anet~CO2*DSW, data = GE)
boxcox(photo.DSW.m, lambda = seq(-1,2, length=10)) #sqrt transform

# run lmer on photosynthesis data
photo.DSW.lme <- lmer(Anet~CO2*DSW +  (1|Ring/Plot_ID), GE)
Anova(photo.DSW.lme, test.statistic = "F")

# run pairwise means
emm.photo <- emmeans(photo.DSW.lme, ~ CO2 |DSW)
pairs(emm.photo, adjust = "tukey")

# plot boxcox for gs data and transform as necessary

cond.DSW.m <- aov(Gs~CO2*DSW, data = GE)
boxcox(cond.DSW.m, lambda = seq(-1,2, length=10)) #sqrt transform

cond.DSW.lme <- lmer(Gs^0.5~CO2*DSW +  (1|Ring/Plot_ID), GE)
Anova(cond.DSW.lme)

emm.rswc <- emmeans(cond.DSW.lme, ~ CO2|DSW)
pairs(emm.rswc, adjust = "tukey")


# Put days since watering into categories: 1-6 and 7-10 days
GE$DSW.CAT = ifelse(GE$DSW %in% c("1", "2","3", "4","5","6"), "1 - 6 days", "7 - 10 days")
GE$DSW.CAT <- factor(GE$DSW.CAT, levels=c("1 - 6 days", "7 - 10 days"))

# Run lmer for Gs category data

cond.cat.lme <- lmer(Gs^0.5~CO2*DSW.CAT +  (1|Ring/Plot_ID), GE)
Anova(cond.cat.lme)

emm.cond.cat <- emmeans(cond.cat.lme, ~ DSW.CAT)
pairs(emm.cond.cat, adjust = "tukey")

# Run lmer for photosynthesis category data

photo.cat.lme <- lmer(Anet~CO2*DSW.CAT + (1|Ring/Plot_ID), GE)
Anova(photo.cat.lme, test.statistic = "F")

# pairwise analysis
emm.photo.cat <- emmeans(photo.cat.lme, ~ DSW.CAT)
pairs(emm.photo.cat, adjust = "tukey")



######################################################################

#### Intrinsic water use efficiency code

GE.WUE <- GE %>%
  group_by(Plot_ID, CO2, Ring, Water, Month) %>%
  dplyr::mutate(WUE= Anet/Gs)

# plot boxcox for Anet data and transform as necessary

COND.PHOTO.FREQ.m <- aov(Anet~Gs*CO2*Water*Month , data = GE.WUE)
boxcox(COND.PHOTO.FREQ.m, lambda = seq(-1,1.5, length=10)) #log transform

# Run lmer for photosynthesis as a function of Gs category data

COND.PHOTO.FREQ.lme <- lmer(Anet^0.5~Gs + (1|Ring/Plot_ID), GE.WUE)
Anova(COND.PHOTO.FREQ.lme, test.statistic = "F")

# find R^2
r.squaredGLMM(COND.PHOTO.FREQ.lme)


# plot boxcox for iWUE data and transform as necessary

WUE.FREQ.m <- aov(WUE~CO2*Water*Month , data = GE.WUE)
boxcox(WUE.FREQ.m, lambda = seq(-1,1.5, length=10)) #log transform

# Run lmer for iWUE as a function of Gs

WUE.FREQ.lme <- lmer(log(WUE)~CO2*Water*Month + (1|Ring/Plot_ID), GE.WUE)
Anova(WUE.FREQ.lme, test.statistic = "F")


#####################################################################

Summer.WP <- read.csv("water_potential.csv")

Summer.WP$Water <- factor(Summer.WP$Water, levels = c("3 days", "5 days", "10 days"))
Summer.WP$CO2 <- factor(Summer.WP$CO2, levels=c("Elevated", "Ambient"))
Summer.WP$DSW <- factor(Summer.WP$DSW)
Summer.WP$Time.of.day <- factor(Summer.WP$Time.of.day, levels= c("Pre-dawn", "Midday"))
Summer.WP$Plot_ID <- paste(Summer.WP$Ring,Summer.WP$Water, sep="_")

# analyse pre-dawn and midday separately 

Summer.WP.midday <- subset(Summer.WP, Time.of.day =="Midday")
Summer.WP.predawn <- subset(Summer.WP, Time.of.day =="Pre-dawn")

# PREDAWN

# plot boxcox for Predawn WP  and transform as necessary
WP.m <- aov(Mpa.positive~CO2*Water, data = Summer.WP.predawn)
boxcox(WP.m, lambda = seq(-1.5,1.5, length=10)) #log transform

# Run lmer testing the effects of CO2 and WATER on predawn WP 

WP.lme1 <- lmer(Mpa.positive^0.5~CO2*Water + (1|Ring), Summer.WP.predawn)
Anova(WP.lme1, test.statistic = "F")

# Run lmer testing the effects of DSW on predawn WP 

WP.lme1 <- lmer(Mpa.positive^0.5~DSW + (1|Ring), Summer.WP.predawn)
Anova(WP.lme1, test.statistic = "F")

emm.wp <- emmeans(WP.lme1, ~ DSW)
pairs(emm.wp, adjust = "tukey")

# MIDDAY WP

# plot boxcox for midday WP  and transform as necessary
WP.m <- aov(Mpa.positive~CO2*Water, data = Summer.WP.midday)
boxcox(WP.m, lambda = seq(-1.5,1.5, length=10)) #log transform

# Run lmer testing the effects of CO2 and WATER on midday WP 

WP.lme1 <- lmer(Mpa.positive^0.5~CO2*Water + (1|Ring), Summer.WP.midday)
Anova(WP.lme1, test.statistic = "F")

# Run lmer testing the effects of DSW on midday WP 

WP.lme1 <- lmer(Mpa.positive^0.5~DSW + (1|Ring), Summer.WP.midday)
Anova(WP.lme1, test.statistic = "F")



########################################################################
############## Leaf Area Index analysis

# read file
LAI.WT <-read.csv("LAI.csv", na.strings = "na")

# make treatments factors

LAI.WT$CO2 <- factor(LAI.WT$CO2, levels=c("Elevated", "Ambient"))
LAI.WT$Ring <- factor(LAI.WT$Ring)
LAI.WT$Water <- factor(LAI.WT$Water, levels=c("3 days", "5 days", "10 days"))
LAI.WT$Month <- factor(LAI.WT$Month, levels= c("Jun '20","Aug '20", "Sep '20", "Nov '20", "Dec '20", "Jan '21", "Feb '21", "Mar '21", "May '21", "Jul '21", "Aug '21"))
LAI.WT$Plot_ID <- paste(LAI.WT$Ring,LAI.WT$Water, sep="_")

# plot boxcox for LAI  and transform as necessary

LAI.m <- aov(LAI~CO2*Water*Month, data = LAI.WT)
boxcox(LAI.m, lambda = seq(0,1, length=10)) #sqrt transform 

# Looking for effects on LAI
# run lmer on photosynthesis data, testing for CO2, Water and Month effects

LAI.me <- lmer(LAI^0.5~CO2*Month*Water + (1|Ring/Plot_ID), LAI.WT)
Anova(LAI.me, test.statistic = "F") ### significant interaction of CO2*Month > plot

emm.LAI <- emmeans(LAI.me, ~ CO2|Month)
pairs(emm.LAI, adjust = "tukey")
