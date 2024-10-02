# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# clear memory
rm(list = ls())

# 
library(dplyr)
library(lubridate)
library(hmmTMB)
library(HDInterval)


theme_set(theme_bw())
pal <- hmmTMB:::hmmTMB_cols
##### data preparation #####

sticklebacks <- read.csv("data/regularTS1min.csv")

# create time variables
sticklebacks$time <- strptime(sticklebacks$DateTime, 
                              format = "%Y-%m-%d %H:%M:%S")
sticklebacks$time <- as.POSIXct(sticklebacks$time)
colnames(sticklebacks)[which(colnames(sticklebacks) == "Time")] <- "Clock"
sticklebacks$Clock <- format(sticklebacks$time, 
                             format = "%H:%M:%S")
sticklebacks$tod <- as.numeric(difftime(
  as.POSIXct(sticklebacks$Clock, format = '%H:%M:%S'), 
  as.POSIXct('00:00', format = '%H:%M'), 
  units = 'min')) + 1

# change variable day
sticklebacks$Date <- as.Date(format(sticklebacks$time,
                                    format = "%Y-%m-%d"))
df_help <- sticklebacks %>% group_by(ID) %>% 
  mutate(days = as.numeric(difftime(Date, Date[1], 
                                    units = 'day')) + 1,
         index = 1:n()) %>% ungroup()
sticklebacks$days <- df_help$days
sticklebacks$index <- df_help$index
sticklebacks$index[which(sticklebacks$ID == "F47.1")] <- 
  2:(sum(sticklebacks$ID == "F47.1") + 1)
rm(df_help)

# covariates

sticklebacks$ID <- as.factor(sticklebacks$ID)
sticklebacks$Fish_id <- as.factor(sticklebacks$Fish_id)

sticklebacks$Early <- as.factor(sticklebacks$Early)
levels(sticklebacks$Early)
sticklebacks$Early <- relevel(sticklebacks$Early, ref = "1")

sticklebacks$Exposed <- as.factor(sticklebacks$Exposed)
levels(sticklebacks$Exposed)
sticklebacks$Infected <- as.factor(sticklebacks$Infected)
levels(sticklebacks$Infected)

sticklebacks <- sticklebacks %>% mutate(group = case_when(
  Infected == 0 & Exposed == 0 & Early == 1 ~ 1,
  Infected == 0 & Exposed == 1 & Early == 1 ~ 2,
  Infected == 1 & Exposed == 1 & Early == 1 ~ 3,
  Infected == 0 & Exposed == 0 & Early == 0 ~ 4,
  Infected == 0 & Exposed == 1 & Early == 0 ~ 5,
  Infected == 1 & Exposed == 1 & Early == 0 ~ 6
))
sticklebacks$group <- as.factor(sticklebacks$group)

sticklebacks$cosFull <- cos(2 * pi * sticklebacks$tod / 1440)
sticklebacks$sinFull <- sin(2 * pi * sticklebacks$tod / 1440)
sticklebacks$cosHalf <- cos(2 * pi * sticklebacks$tod / 720)
sticklebacks$sinHalf <- sin(2 * pi * sticklebacks$tod / 720)

# final HMM
hmm <- readRDS("outputs/3stateHMM_restrict_group&tod_int_reFishID.rds")


# Get the predictions for each treatment group


ControlEarly <- readRDS("outputs/probDraws_group1.rds")
ControlLate <- readRDS("outputs/probDraws_group4.rds")

ExposedEarly <- readRDS("outputs/probDraws_group2.rds")
ExposedLate <- readRDS("outputs/probDraws_group5.rds")

InfectedEarly <- readRDS("outputs/probDraws_group3.rds")
InfectedLate <- readRDS("outputs/probDraws_group6.rds")

# Get the predictions for each treatment group



LOS <- function(x){
   100* length(which(x >0))/length(x)
}

# estimate the proportion of overlaping states between the two groups
tod <- 1:1440 /60

Night_Time <- which(tod < 6 | tod > 21)


## Results for probability of sleep (state 1)


##########------------------------- Early infection stages (first 4 days)
EE_state1 <- (ExposedEarly[, , 1] - ControlEarly[, , 1]) *100  # as a matrix
IE_state1 <- (InfectedEarly[, , 1] - ControlEarly[, , 1]) *100  # as a matrix


#####-------------------------------- Differences during the night
# Differences in the probability of sleep during the night between infected and control individuals
mean(c(IE_state1[,Night_Time]))
hdi(c(IE_state1[,Night_Time]))
LOS(c(IE_state1[,Night_Time]))

# Differences between infected and exposed but not infected in early stages
mean(c(IE_state1[,Night_Time]- EE_state1[,Night_Time]))
hdi(c(IE_state1[,Night_Time]- EE_state1[,Night_Time]))

# Differences during the day

mean(c(IE_state1[,-Night_Time]))
hdi(c(IE_state1[,-Night_Time]))
LOS(c(IE_state1[,-Night_Time]))

mean(c(EE_state1[,-Night_Time]))
hdi(c(EE_state1[,-Night_Time]))
LOS(c(EE_state1[,-Night_Time]))


# Differences between infected and exposed but not infected in early stages
mean(c(IE_state1[,-Night_Time]- EE_state1[,-Night_Time]))
hdi(c(IE_state1[,-Night_Time]- EE_state1[,-Night_Time]))
LOS(c(IE_state1[,-Night_Time]- EE_state1[,-Night_Time]))





##--------------------------------- Late infection stages (day 29 to 32)
EL_state1 <- (ExposedLate[, , 1] - ControlLate[, , 1])  * 100  # as a matrix
IL_state1 <- (InfectedLate[, , 1] - ControlLate[, , 1]) * 100  # as a matrix


# Night time differences
mean(c(IL_state1[,Night_Time]))
hdi(c(IL_state1[,Night_Time]))
LOS(c(IL_state1[,Night_Time]))


mean(c(IL_state1[,Night_Time]- EL_state1[,Night_Time]))
hdi(c(IL_state1[,Night_Time]- EL_state1[,Night_Time]))



# Overall differences between infected and control from early to late stages

mean(c(InfectedEarly[, , 1][,Night_Time])*100)
hdi(c(InfectedEarly[, , 1][,Night_Time])*100)

mean(c(InfectedLate[, , 1][,Night_Time])*100)
hdi(c(InfectedLate[, , 1][,Night_Time])*100)


mean(c(InfectedEarly[, , 1][,-Night_Time])*100)
hdi(c(InfectedEarly[, , 1][,-Night_Time])*100)

mean(c(InfectedLate[, , 1][,-Night_Time])*100)
hdi(c(InfectedLate[, , 1][,-Night_Time])*100)



mean(c(IL_state1- IE_state1))
hdi(c(IL_state1- IE_state1))
LOS(c(IL_state1- IE_state1))

mean(c(IL_state1[,Night_Time]- IE_state1[,Night_Time]))
hdi(c(IL_state1[,Night_Time]- IE_state1[,Night_Time]))
LOS(c(IL_state1[,Night_Time] - IE_state1[,Night_Time]))



# run script to get the state probabilities
namePlot = "Prob_state1"
Y1axis_name = "Prob. of Sleep (% of state 1)"
Y2axis_name = "Ln(Treatment/Control)"

source("R/state1_plots.R")


# run script to get the state 2
namePlot = "Prob_state2"
Y1axis_name = "Prob. of state 2"
# Y2axis_name = "Ln(Treatment/Control)"
source("R/state2_plots.R")

# run script to get the state 3
namePlot = "Prob_state3"
Y1axis_name = "Prob. of state 3"
# Y2axis_name = "Ln(Treatment/Control)"
source("R/state3_plots.R")

###
## Results for Dwell time


ControlEarly <- readRDS("outputs/dwellDraws_group1.rds")
ControlLate <- readRDS("outputs/dwellDraws_group4.rds")

ExposedEarly <- readRDS("outputs/dwellDraws_group2.rds")
ExposedLate <- readRDS("outputs/dwellDraws_group5.rds")

InfectedEarly <- readRDS("outputs/dwellDraws_group3.rds")
InfectedLate <- readRDS("outputs/dwellDraws_group6.rds")



EE_state1 <- (ExposedEarly[, , 1] - ControlEarly[, , 1])   # as a matrix
IE_state1 <- (InfectedEarly[, , 1] - ControlEarly[, , 1])   # as a matrix
#####-------------------------------- Differences during the night
# Differences in the probability of sleep during the night between infected and control individuals
mean(c(IE_state1[,Night_Time]))
hdi(c(IE_state1[,Night_Time]))
LOS(c(IE_state1[,Night_Time]))

# Differences between infected and exposed but not infected in early stages
mean(c(IE_state1[,Night_Time]- EE_state1[,Night_Time]))
hdi(c(IE_state1[,Night_Time]- EE_state1[,Night_Time]))
LOS(c(IE_state1[,Night_Time]- EE_state1[,Night_Time]))

# Differences during the day

mean(c(IE_state1[,-Night_Time]))
hdi(c(IE_state1[,-Night_Time]))
LOS(c(IE_state1[,-Night_Time]))

# Differences between infected and exposed but not infected in early stages
mean(c(IE_state1[,-Night_Time]- EE_state1[,-Night_Time]))
hdi(c(IE_state1[,-Night_Time]- EE_state1[,-Night_Time]))
LOS(c(IE_state1[,-Night_Time]- EE_state1[,-Night_Time]))





########------------------------ Late infection stages (day 29 to 32)


EL_state1 <- (ExposedLate[, , 1] - ControlLate[, , 1])   # as a matrix
IL_state1 <- (InfectedLate[, , 1] - ControlLate[, , 1])   # as a matrix

# Night time differences
mean(c(IL_state1[,Night_Time]))
hdi(c(IL_state1[,Night_Time]))

# Overall differences between infected and control from early to late stages
mean(c(InfectedEarly[, , 1][,Night_Time]))
hdi(c(InfectedEarly[, , 1][,Night_Time]))

mean(c(InfectedLate[, , 1][,Night_Time]))
hdi(c(InfectedLate[, , 1][,Night_Time]))


# Overall differences between infected and control from early to late stages
mean(c(InfectedEarly[, , 1][,-Night_Time]))
hdi(c(InfectedEarly[, , 1][,-Night_Time]))

mean(c(InfectedLate[, , 1][,-Night_Time]))
hdi(c(InfectedLate[, , 1][,-Night_Time]))

# 


mean(c(IL_state1[,which(tod > 12 & tod < 21)]))
hdi(c(IL_state1[,which(tod > 12 & tod < 21)]))
LOS(c(IL_state1[,which(tod > 12 & tod < 21)]))

mean(c(IL_state1[,which(tod > 12 & tod < 21)]- EL_state1[,which(tod > 12 & tod < 21)]))
hdi(c(IL_state1[,which(tod > 12 & tod < 21)]- EL_state1[,which(tod > 12 & tod < 21)]))
LOS(c(IL_state1[,which(tod > 12 & tod < 21)]- EL_state1[,which(tod > 12 & tod < 21)]))


mean(c(IL_state1[,which(tod > 18 & tod < 21)]))
hdi(c(IL_state1[,which(tod > 18 & tod < 21)]))
LOS(c(IL_state1[,which(tod > 18 & tod < 21)]))


mean(c(IL_state1[,which(tod > 18 & tod < 21)]- EL_state1[,which(tod > 18 & tod < 21)]))
hdi(c(IL_state1[,which(tod > 18 & tod < 21)]- EL_state1[,which(tod > 18 & tod < 21)]))
LOS(c(IL_state1[,which(tod > 18 & tod < 21)]- EL_state1[,which(tod > 18 & tod < 21)]))



# run script to get the state probabilities
namePlot = "Dwell_Prob_state1"
Y1axis_name = "Dwell time (min)"
Y2axis_name = "Ln(Treatment/Control)"

source("R/state1_plots.R")


# run script to get the state 2
namePlot = "Dwell_Prob_state2"
# Y1axis_name = "Prob. of state 2"
# Y2axis_name = "Ln(Treatment/Control)"
source("R/state2_plots.R")

# run script to get the state 3
namePlot = "Dwell_Prob_state3"
# Y1axis_name = "Prob. of state 3"
# Y2axis_name = "Ln(Treatment/Control)"
source("R/state3_plots.R")

