# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(lubridate)
library(hmmTMB)


##### functions #####

getDelta <- function(Gamma, cycle = 1440, N = 3){
  
  deltas <- matrix(NA, nrow = cycle, ncol = N)
  index <- c(1:cycle,1:cycle) # index makes it easier to choose the correct Gamma in the matrix multiplications
  
  # get Gamma_star for periodic stationary distributions
  Gamma_star <- Gamma # Gamma_star will be the multiplication of the next 24 hourly(or 48 half-hourly) Gamma matrices
  for (t in 1:cycle){
    # loopi determines which gamma matrices will be multiplied for which t
    loopi <- index[seq(t+1,cycle+t-1, by = 1)]
    for (i in loopi){
      Gamma_star[, , t] <- Gamma_star[, , t]%*%Gamma[,,i]
    }
  }
  # get periodic stationary distributions
  for (t in 1:cycle) {
    deltas[t, ] <- solve(t(diag(N) - Gamma_star[, , t] + 1), rep(1, N))
  }
  return(deltas)
}

statProbs <- function(beta, tod = 1:1440, group = 1){
  group2 = ifelse(group == 2, 1, 0)
  group3 = ifelse(group == 3, 1, 0)
  group4 = ifelse(group == 4, 1, 0)
  group5 = ifelse(group == 5, 1, 0)
  group6 = ifelse(group == 6, 1, 0)
  n <- length(tod)
  N <- 3
  beta12 <- beta[1:18]
  beta13 <- beta[19:36]
  beta21 <- beta[37:54]
  beta23 <- beta[55:72]
  beta31 <- beta[73:90]
  beta32 <- beta[91:108]
  covsMat <- data.frame(Intercept = 1, 
                        group2 = group2,
                        group3 = group3,
                        group4 = group4,
                        group5 = group5,
                        group6 = group6,
                        cosFull = cos(2 * pi * tod / n),
                        sinFull = sin(2 * pi * tod / n),
                        int2_cosFull = cos(2 * pi * tod / n) * group2,
                        int3_cosFull = cos(2 * pi * tod / n) * group3,
                        int4_cosFull = cos(2 * pi * tod / n) * group4,
                        int5_cosFull = cos(2 * pi * tod / n) * group5,
                        int6_cosFull = cos(2 * pi * tod / n) * group6,
                        int2_sinFull = sin(2 * pi * tod / n) * group2,
                        int3_sinFull = sin(2 * pi * tod / n) * group3,
                        int4_sinFull = sin(2 * pi * tod / n) * group4,
                        int5_sinFull = sin(2 * pi * tod / n) * group5,
                        int6_sinFull = sin(2 * pi * tod / n) * group6)
  Gamma12 <- beta12 %*% t(as.matrix(covsMat))
  Gamma13 <- beta13 %*% t(as.matrix(covsMat))
  Gamma21 <- beta21 %*% t(as.matrix(covsMat))
  Gamma23 <- beta23 %*% t(as.matrix(covsMat))
  Gamma31 <- beta31 %*% t(as.matrix(covsMat))
  Gamma32 <- beta32 %*% t(as.matrix(covsMat))
  Gammas <- array(1, dim = c(N, N, n))
  Gammas[1, 2, ] <- exp(Gamma12)
  Gammas[1, 3, ] <- exp(Gamma13)
  Gammas[2, 1, ] <- exp(Gamma21)
  Gammas[2, 3, ] <- exp(Gamma23)
  Gammas[3, 1, ] <- exp(Gamma31)
  Gammas[3, 2, ] <- exp(Gamma32)
  for (i in 1:n) {
    Gammas[, , i] <- Gammas[, , i] / rowSums(Gammas[, , i])
  }
  delta <- getDelta(Gammas)
  # delta <- matrix(NA, n, N)
  # for(i in 1:n){
  #   delta[i, ] <- solve(t(diag(N) - Gammas[, , i] + 1), rep(1, N))
  # }
  return(delta)
}


ddwell_t_new = function(
    r, # vector of dwell-times to compute, must be of the form 1:R
    t, # time point where the state is entered
    state, # which state to compute
    Gamma # array of dim c(N,N,L)
){
  L = dim(Gamma)[3] # cycle length
  I = (t+r-1)%%L
  I[which(I==0)] = L
  gamma_ii = Gamma[state,state,I]
  pmf = c(1, cumprod(gamma_ii)[-length(I)])*(1-gamma_ii)
  return(pmf)
}

getGammas <- function(beta, tod = 1:1440, group = 1){
  group2 = ifelse(group == 2, 1, 0)
  group3 = ifelse(group == 3, 1, 0)
  group4 = ifelse(group == 4, 1, 0)
  group5 = ifelse(group == 5, 1, 0)
  group6 = ifelse(group == 6, 1, 0)
  n <- length(tod)
  N <- 3
  beta12 <- beta[1:18]
  beta13 <- beta[19:36]
  beta21 <- beta[37:54]
  beta23 <- beta[55:72]
  beta31 <- beta[73:90]
  beta32 <- beta[91:108]
  covsMat <- data.frame(Intercept = 1, 
                        group2 = group2,
                        group3 = group3,
                        group4 = group4,
                        group5 = group5,
                        group6 = group6,
                        cosFull = cos(2 * pi * tod / n),
                        sinFull = sin(2 * pi * tod / n),
                        int2_cosFull = cos(2 * pi * tod / n) * group2,
                        int3_cosFull = cos(2 * pi * tod / n) * group3,
                        int4_cosFull = cos(2 * pi * tod / n) * group4,
                        int5_cosFull = cos(2 * pi * tod / n) * group5,
                        int6_cosFull = cos(2 * pi * tod / n) * group6,
                        int2_sinFull = sin(2 * pi * tod / n) * group2,
                        int3_sinFull = sin(2 * pi * tod / n) * group3,
                        int4_sinFull = sin(2 * pi * tod / n) * group4,
                        int5_sinFull = sin(2 * pi * tod / n) * group5,
                        int6_sinFull = sin(2 * pi * tod / n) * group6)
  Gamma12 <- beta12 %*% t(as.matrix(covsMat))
  Gamma13 <- beta13 %*% t(as.matrix(covsMat))
  Gamma21 <- beta21 %*% t(as.matrix(covsMat))
  Gamma23 <- beta23 %*% t(as.matrix(covsMat))
  Gamma31 <- beta31 %*% t(as.matrix(covsMat))
  Gamma32 <- beta32 %*% t(as.matrix(covsMat))
  Gammas <- array(1, dim = c(N, N, n))
  Gammas[1, 2, ] <- exp(Gamma12)
  Gammas[1, 3, ] <- exp(Gamma13)
  Gammas[2, 1, ] <- exp(Gamma21)
  Gammas[2, 3, ] <- exp(Gamma23)
  Gammas[3, 1, ] <- exp(Gamma31)
  Gammas[3, 2, ] <- exp(Gamma32)
  for (i in 1:n) {
    Gammas[, , i] <- Gammas[, , i] / rowSums(Gammas[, , i])
  }
  return(Gammas)
}



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



##### simulate distribution of stationary state probabilities #####

# draws <- 1000
# betaDraws <- hmm$post_coeff(n_post = draws)
# betaDraws <- betaDraws[, which(colnames(betaDraws) == "coeff_fe_hid")]

# newdata <- data.frame(Late = rep(c(0, 1), each = 3), 
#                       Exposed = rep(c(0, 1, 1), 2),
#                       Infected = rep(c(0, 0, 1), 2) )

# for (j in 1:6) { 
#   deltaDraws <- array(dim = c(draws, 1440, 3)) 
#   for (i in 1:draws) {
#     beta <- as.vector(betaDraws[i, ])
#     deltaDraws[i, , ] <- statProbs(beta, group = j)
#   }
#   saveRDS(deltaDraws, file = paste0("probDraws_group", j, ".rds"))
# }

# # load data array (for example, group 6)
# test <- readRDS("probDraws_group6.rds")

# # note: which group number corresponds to which treatment can be checked in 
# # the data frame called "newdata" (see line 224)

# # get data frame for each state of group 6:
# test_state1 <- test[, , 1] # as a matrix or...

# dim(test_state1)

# test_state2 <- as.data.frame(test[, , 2]) # ...as a data frame
# test_state3 <- test[, , 3]



##### simulate distribution of state dwell times #####

draws <- 1000
betaDraws <- hmm$post_coeff(n_post = draws)
betaDraws <- betaDraws[, which(colnames(betaDraws) == "coeff_fe_hid")]

newdata <- data.frame(Late = rep(c(0, 1), each = 3), 
                      Exposed = rep(c(0, 1, 1), 2),
                      Infected = rep(c(0, 0, 1), 2) )

L <- 1440
r <- 1:L
for (j in 1:6) { 
  dwellDraws <- array(NA, dim = c(draws, 1440, 3)) 
  for (i in 1:draws) {
    beta <- as.vector(betaDraws[i, ])
    gammas <- getGammas(beta, group = j)
    for (l in 1:L) {
      tDwell <- ddwell_t_new(r, t = l, state = 1, Gamma = gammas)
      dwellDraws[i, l, 1] <- (L + sum(r * tDwell)) / sum(tDwell) - L
    }
    for (l in 1:L) {
      tDwell <- ddwell_t_new(r, t = l, state = 2, Gamma = gammas)
      dwellDraws[i, l, 2] <- (L + sum(r * tDwell)) / sum(tDwell) - L
    }
    for (l in 1:L) {
      tDwell <- ddwell_t_new(r, t = l, state = 3, Gamma = gammas)
      dwellDraws[i, l, 3] <- (L + sum(r * tDwell)) / sum(tDwell) - L
    }
  }
  saveRDS(dwellDraws, file = paste0("dwellDraws_group", j, ".rds"))
}

# load data array (for example, group 1)
example <- readRDS("outputs/dwellDraws_group1.rds")

# note: which group number corresponds to which treatment can be checked in 
# the data frame "newdata" (see line 256) --- same as for state probabilities

# get data frame for each state of group 1:
example_state1 <- example[, , 1] # as a matrix or...
example_state2 <- as.data.frame(example[, , 2]) # ...as a data frame
example_state3 <- example[, , 3]
