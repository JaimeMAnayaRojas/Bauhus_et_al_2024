# clear memory
rm(list = ls())

library(dplyr)
library(lubridate)
library(hmmTMB)
library(ggplot2)
library(gridExtra)
library(ggpubr)
theme_set(theme_bw())
pal <- hmmTMB:::hmmTMB_cols
cbPalette <- c("dimgrey", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")



##### data preparation #####

sticklebacks <- read.csv("data/regularTS1min.csv")

# create time variables
sticklebacks$time <- strptime(sticklebacks$DateTime, format = "%Y-%m-%d %H:%M:%S")
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
                                    units = 'day')) + 1) %>% ungroup()
sticklebacks$days <- df_help$days
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

sticklebacks$Late <- 0
sticklebacks$Late[which(sticklebacks$Early == 0)] <- 1
sticklebacks$Late <- as.factor(sticklebacks$Late)
levels(sticklebacks$Late)

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

# dummy variables

# sticklebacks <- sticklebacks %>% mutate(group2 = ifelse(group == 2, 1, 0),
#                                         group3 = ifelse(group == 3, 1, 0),
#                                         group4 = ifelse(group == 4, 1, 0),
#                                         group5 = ifelse(group == 5, 1, 0),
#                                         group6 = ifelse(group == 6, 1, 0))

head(sticklebacks)



##### 3-state basic HMM #####

tpm0 <- matrix(c(0.96, 0.04, 0,
                 0.02, 0.96, 0.02,
                 0, 0.06, 0.94),
               nrow = 3, byrow = TRUE)
hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                       initial_state = "stationary") 
obs <- Observation$new(data = sticklebacks, n_states = 3,
                       dists = list(Loco_sum = "gamma2"),
                       par = list(Loco_sum = list(mean = c(0.4, 1, 2.5),
                                                  sd = c(0.4, 0.8, 2))) )
fixpar <- list(hid = c("S1>S3.(Intercept)" = NA,
                       "S3>S1.(Intercept)" = NA))
hmm <- HMM$new(obs = obs, hid = hid, fixpar = fixpar)
hmm$fit(itnmax = 10000, control = list(eval.max = 10000))

saveRDS(hmm, "outputs/3stateHMM_restricted.rds")

### random starting values (not working?!)
# llks <- rep(NA, 10)
# mods <- list()
# tpm0 <- matrix(c(0.9, 0.1, 0,
#                  0.05, 0.9, 0.05,
#                  0, 0.1, 0.9),
#                nrow = 3, byrow = TRUE)
# for (k in 1:10){
#   hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
#                          initial_state = "stationary") 
#   par0 <- list(Loco_sum = list(mean = runif(3, 0.01, 5),
#                                sd = runif(3, 0.1, 2.5) ))
#   obs <- Observation$new(data = sticklebacks, n_states = 3,
#                          dists = list(Loco_sum = "gamma2"),
#                          par = par0)
#   hmm <- HMM$new(obs = obs, hid = hid, fixpar = fixpar)
#   hmm$fit(itnmax = 10000, control = list(eval.max = 10000))
#   mods[[k]] <- hmm
#   llks[k] <- hmm$llk()
# }



##### covariates: late, exposed, infected #####

hmm <- readRDS("outputs/3stateHMM_restricted.rds")

f <- ~ Late + Exposed + Infected
tpm0 <- matrix(c(0.96, 0.04, 0,
                 0.02, 0.96, 0.02,
                 0, 0.06, 0.94),
               nrow = 3, byrow = TRUE)
hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                       initial_state = "stationary", formula = f) 
obs <- Observation$new(data = sticklebacks, n_states = 3,
                       dists = list(Loco_sum = "gamma2"),
                       par = list(Loco_sum = list(mean = c(0.4, 1, 2.5),
                                                  sd = c(0.4, 0.8, 2))) )
fixpar <- list(hid = c("S1>S3.(Intercept)" = NA,
                       "S1>S3.Late1" = NA,
                       "S1>S3.Exposed1" = NA,
                       "S1>S3.Infected1" = NA,
                       "S3>S1.(Intercept)" = NA,
                       "S3>S1.Late1" = NA,
                       "S3>S1.Exposed1" = NA,
                       "S3>S1.Infected1" = NA
                       ))
hmm_new <- HMM$new(obs = obs, hid = hid, init = hmm, fixpar = fixpar)
# hmm_new <- update(hmm_new, type = "hid", i = 1, j = 3, change = ~ ., fit = FALSE)
# hmm_new <- update(hmm_new, type = "hid", i = 3, j = 1, change = ~ ., fit = FALSE)
hmm_new$fit(itnmax = 10000, control = list(eval.max = 10000))

saveRDS(hmm_new, "outputs/3stateHMM_restrict_covs.rds")
hmm_new <- readRDS("outputs/3stateHMM_restrict_covs.rds")

hmm_new$out()
hmm_new$AIC_marginal()
hmm_new$confint()

newdata <- data.frame(Late = rep(c(0, 1), each = 3), 
                      Exposed = rep(c(0, 1, 1), 2),
                      Infected = rep(c(0, 0, 1), 2) )
gammas <- hmm_new$predict("tpm", newdata = newdata) 


# random starting values 

runs <- 10
allAIC <- numeric(runs)
mods <- list()
for (k in 1:runs){
  tpm0 <- matrix(c(0.96, 0.04, 0,
                   0.02, 0.96, 0.02,
                   0, 0.06, 0.94),
                 nrow = 3, byrow = TRUE)
  hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                         initial_state = "stationary", formula = f)
  coeff_fe <- hid$coeff_fe()
  coeff_fe[which(coeff_fe != 0 & coeff_fe != -Inf)] <- runif(4, -5, -2.5)
  coeff_fe[c(2:4, 10:12, 14:16, 22:24)] <- runif(12, -3, 3)
  hid$update_coeff_fe(coeff_fe = coeff_fe)
  par0 <- list(Loco_sum = list(mean = c(rnorm(1, 0.4, 0.01),
                                        rnorm(1, 1, 0.05),
                                        rnorm(1, 2.5, 0.1)),
                               sd = c(rnorm(1, 0.4, 0.05),
                                      rnorm(1, 0.8, 0.1),
                                      rnorm(1, 2, 0.2)) ))
  obs <- Observation$new(data = sticklebacks, n_states = 3,
                         dists = list(Loco_sum = "gamma2"),
                         par = par0)
  hmm <- HMM$new(obs = obs, hid = hid, fixpar = fixpar)
  hmm$fit(silent = TRUE, itnmax = 10000, control = list(eval.max = 10000))
  allAIC[k] <- hmm$AIC_marginal()
  mods[[k]] <- hmm
  gc()
}
summary(allAIC)
plot(allAIC)



##### covariates: late, exposed, infected & tod #####

# previous HMM
hmm <- readRDS("outputs/3stateHMM_restrict_covs.rds")

f <- ~ Late + Exposed + Infected + cosFull + sinFull
tpm0 <- matrix(c(0.96, 0.04, 0,
                 0.02, 0.96, 0.02,
                 0, 0.06, 0.94),
               nrow = 3, byrow = TRUE)
hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                       initial_state = "stationary", formula = f)
obs <- Observation$new(data = sticklebacks, n_states = 3,
                       dists = list(Loco_sum = "gamma2"),
                       par = list(Loco_sum = list(mean = c(0.4, 1, 2.5),
                                                  sd = c(0.4, 0.8, 2))) )
fixpar <- list(hid = c("S1>S3.(Intercept)" = NA,
                       "S1>S3.Late1" = NA,
                       "S1>S3.Exposed1" = NA,
                       "S1>S3.Infected1" = NA,
                       "S1>S3.cosFull" = NA,
                       "S1>S3.sinFull" = NA,
                       "S3>S1.(Intercept)" = NA,
                       "S3>S1.Late1" = NA,
                       "S3>S1.Exposed1" = NA,
                       "S3>S1.Infected1" = NA,
                       "S3>S1.cosFull" = NA,
                       "S3>S1.sinFull" = NA
                       ))
hmm_new <- HMM$new(obs = obs, hid = hid, init = hmm, fixpar = fixpar)
hmm_new$fit(itnmax = 10000, control = list(eval.max = 10000))

saveRDS(hmm_new, "outputs/3stateHMM_restrict_covs&tod.rds")

hmm_new$out()
hmm_new$confint()

hmm_new$AIC_marginal()
hmm$AIC_marginal()



##### covariates: interactions of late, exposed, infected & tod #####

# previous HMM
hmm <- readRDS("outputs/3stateHMM_restrict_covs&tod.rds")

f <- ~ Late * cosFull + Exposed * cosFull + Infected * cosFull + 
  Late * sinFull + Exposed * sinFull + Infected * sinFull
tpm0 <- matrix(c(0.96, 0.04, 0,
                 0.02, 0.96, 0.02,
                 0, 0.06, 0.94),
               nrow = 3, byrow = TRUE)
hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                       initial_state = "stationary", formula = f)
obs <- Observation$new(data = sticklebacks, n_states = 3,
                       dists = list(Loco_sum = "gamma2"),
                       par = list(Loco_sum = list(mean = c(0.4, 1, 2.5),
                                                  sd = c(0.4, 0.8, 2))) )
fixpar <- list(hid = c("S1>S3.(Intercept)" = NA,
                       "S1>S3.Late1" = NA,
                       "S1>S3.cosFull" = NA,
                       "S1>S3.Exposed1" = NA,
                       "S1>S3.Infected1" = NA,
                       "S1>S3.sinFull" = NA,
                       "S1>S3.Late1:cosFull" = NA,   
                       "S1>S3.cosFull:Exposed1" = NA,
                       "S1>S3.cosFull:Infected1" = NA,
                       "S1>S3.Late1:sinFull" = NA,
                       "S1>S3.Exposed1:sinFull" = NA,
                       "S1>S3.Infected1:sinFull" = NA,
                       "S3>S1.(Intercept)" = NA,
                       "S3>S1.Late1" = NA,
                       "S3>S1.cosFull" = NA,
                       "S3>S1.Exposed1" = NA,
                       "S3>S1.Infected1" = NA,
                       "S3>S1.sinFull" = NA,
                       "S3>S1.Late1:cosFull" = NA,   
                       "S3>S1.cosFull:Exposed1" = NA,
                       "S3>S1.cosFull:Infected1" = NA,
                       "S3>S1.Late1:sinFull" = NA,
                       "S3>S1.Exposed1:sinFull" = NA,
                       "S3>S1.Infected1:sinFull" = NA ))
hmm_new <- HMM$new(obs = obs, hid = hid, init = hmm, fixpar = fixpar)
hmm_new$fit(itnmax = 10000, control = list(eval.max = 10000))

saveRDS(hmm_new, "outputs/3stateHMM_restrict_covs&tod_int.rds")

hmm_new$out()
hmm_new$confint()

hmm_new$AIC_marginal()
hmm$AIC_marginal()



##### covariates: random effects, interactions of covs & tod #####

# previous HMM
hmm <- readRDS("outputs/3stateHMM_restrict_covs&tod_int.rds")

f <- ~ Late * cosFull + Exposed * cosFull + Infected * cosFull + 
  Late * sinFull + Exposed * sinFull + Infected * sinFull + s(Fish_id, bs = "re")
tpm0 <- matrix(c(0.96, 0.04, 0,
                 0.02, 0.96, 0.02,
                 0, 0.06, 0.94),
               nrow = 3, byrow = TRUE)
hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                       initial_state = "stationary", formula = f)
obs <- Observation$new(data = sticklebacks, n_states = 3,
                       dists = list(Loco_sum = "gamma2"),
                       par = list(Loco_sum = list(mean = c(0.4, 1, 2.5),
                                                  sd = c(0.4, 0.8, 2))) )
fixpar <- list(hid = c("S1>S3.(Intercept)" = NA,
                       "S1>S3.Late1" = NA,
                       "S1>S3.cosFull" = NA,
                       "S1>S3.Exposed1" = NA,
                       "S1>S3.Infected1" = NA,
                       "S1>S3.sinFull" = NA,
                       "S1>S3.Late1:cosFull" = NA,   
                       "S1>S3.cosFull:Exposed1" = NA,
                       "S1>S3.cosFull:Infected1" = NA,
                       "S1>S3.Late1:sinFull" = NA,
                       "S1>S3.Exposed1:sinFull" = NA,
                       "S1>S3.Infected1:sinFull" = NA,
                       "S3>S1.(Intercept)" = NA,
                       "S3>S1.Late1" = NA,
                       "S3>S1.cosFull" = NA,
                       "S3>S1.Exposed1" = NA,
                       "S3>S1.Infected1" = NA,
                       "S3>S1.sinFull" = NA,
                       "S3>S1.Late1:cosFull" = NA,   
                       "S3>S1.cosFull:Exposed1" = NA,
                       "S3>S1.cosFull:Infected1" = NA,
                       "S3>S1.Late1:sinFull" = NA,
                       "S3>S1.Exposed1:sinFull" = NA,
                       "S3>S1.Infected1:sinFull" = NA ))
hmm_new <- HMM$new(obs = obs, hid = hid, init = hmm, fixpar = fixpar)
hmm_new$fit(itnmax = 10000, control = list(eval.max = 10000))

saveRDS(hmm_new, "outputs/3stateHMM_restrict_covs&tod_int_reFishID.rds")

hmm_new$out()
hmm_new$confint()

hmm_new$AIC_marginal()
hmm$AIC_marginal()



##### covariates: group #####

# previous HMM
hmm <- readRDS("outputs/3stateHMM_restricted.rds")

f <- ~ group 
tpm0 <- matrix(c(0.96, 0.04, 0,
                 0.02, 0.96, 0.02,
                 0, 0.06, 0.94),
               nrow = 3, byrow = TRUE)
hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                       initial_state = "stationary", formula = f)
obs <- Observation$new(data = sticklebacks, n_states = 3,
                       dists = list(Loco_sum = "gamma2"),
                       par = list(Loco_sum = list(mean = c(0.4, 1, 2.5),
                                                  sd = c(0.4, 0.8, 2))) )
fixpar <- list(hid = c("S1>S3.(Intercept)" = NA,
                       "S1>S3.group2" = NA,
                       "S1>S3.group3" = NA,
                       "S1>S3.group4" = NA,
                       "S1>S3.group5" = NA,
                       "S1>S3.group6" = NA,
                       "S3>S1.(Intercept)" = NA,
                       "S3>S1.group2" = NA,
                       "S3>S1.group3" = NA,
                       "S3>S1.group4" = NA,
                       "S3>S1.group5" = NA,
                       "S3>S1.group6" = NA ))
hmm_new <- HMM$new(obs = obs, hid = hid, init = hmm, fixpar = fixpar)
hmm_new$fit(itnmax = 10000, control = list(eval.max = 10000))

saveRDS(hmm_new, "outputs/3stateHMM_restrict_group.rds")

hmm_new$out()
hmm_new$confint()

hmm_new$AIC_marginal()
hmm$AIC_marginal()



##### covariates: group & tod #####

# previous HMM
hmm <- readRDS("outputs/3stateHMM_restrict_group.rds")

f <- ~ group + cosFull + sinFull
tpm0 <- matrix(c(0.96, 0.04, 0,
                 0.02, 0.96, 0.02,
                 0, 0.06, 0.94),
               nrow = 3, byrow = TRUE)
hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                       initial_state = "stationary", formula = f)
obs <- Observation$new(data = sticklebacks, n_states = 3,
                       dists = list(Loco_sum = "gamma2"),
                       par = list(Loco_sum = list(mean = c(0.4, 1, 2.5),
                                                  sd = c(0.4, 0.8, 2))) )
fixpar <- list(hid = c("S1>S3.(Intercept)" = NA,
                       "S1>S3.group2" = NA,
                       "S1>S3.group3" = NA,
                       "S1>S3.group4" = NA,
                       "S1>S3.group5" = NA,
                       "S1>S3.group6" = NA,
                       "S1>S3.cosFull" = NA,
                       "S1>S3.sinFull" = NA,
                       "S3>S1.(Intercept)" = NA,
                       "S3>S1.group2" = NA,
                       "S3>S1.group3" = NA,
                       "S3>S1.group4" = NA,
                       "S3>S1.group5" = NA,
                       "S3>S1.group6" = NA,
                       "S3>S1.cosFull" = NA,
                       "S3>S1.sinFull" = NA ))
hmm_new <- HMM$new(obs = obs, hid = hid, init = hmm, fixpar = fixpar)
hmm_new$fit(itnmax = 10000, control = list(eval.max = 10000))

saveRDS(hmm_new, "outputs/3stateHMM_restrict_group&tod.rds")

hmm_new$out()
hmm_new$confint()

hmm_new$AIC_marginal()
hmm$AIC_marginal()



##### covariates: interactions of group & tod #####

# previous HMM
hmm <- readRDS("outputs/3stateHMM_restrict_group&tod.rds")

f <- ~ group * cosFull + group * sinFull
tpm0 <- matrix(c(0.96, 0.04, 0,
                 0.02, 0.96, 0.02,
                 0, 0.06, 0.94),
               nrow = 3, byrow = TRUE)
hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                       initial_state = "stationary", formula = f)
obs <- Observation$new(data = sticklebacks, n_states = 3,
                       dists = list(Loco_sum = "gamma2"),
                       par = list(Loco_sum = list(mean = c(0.4, 1, 2.5),
                                                  sd = c(0.4, 0.8, 2))) )
fixpar <- list(hid = c("S1>S3.(Intercept)" = NA,
                       "S1>S3.group2" = NA,
                       "S1>S3.group3" = NA,
                       "S1>S3.group4" = NA,
                       "S1>S3.group5" = NA,
                       "S1>S3.group6" = NA,
                       "S1>S3.cosFull" = NA,
                       "S1>S3.sinFull" = NA,
                       "S1>S3.group2:cosFull" = NA,
                       "S1>S3.group3:cosFull" = NA,
                       "S1>S3.group4:cosFull" = NA,
                       "S1>S3.group5:cosFull" = NA,
                       "S1>S3.group6:cosFull" = NA,
                       "S1>S3.group2:sinFull" = NA,
                       "S1>S3.group3:sinFull" = NA,
                       "S1>S3.group4:sinFull" = NA,
                       "S1>S3.group5:sinFull" = NA,
                       "S1>S3.group6:sinFull" = NA,
                       "S3>S1.(Intercept)" = NA,
                       "S3>S1.group2" = NA,
                       "S3>S1.group3" = NA,
                       "S3>S1.group4" = NA,
                       "S3>S1.group5" = NA,
                       "S3>S1.group6" = NA,
                       "S3>S1.cosFull" = NA,
                       "S3>S1.sinFull" = NA,
                       "S3>S1.group2:cosFull" = NA,
                       "S3>S1.group3:cosFull" = NA,
                       "S3>S1.group4:cosFull" = NA,
                       "S3>S1.group5:cosFull" = NA,
                       "S3>S1.group6:cosFull" = NA,
                       "S3>S1.group2:sinFull" = NA,
                       "S3>S1.group3:sinFull" = NA,
                       "S3>S1.group4:sinFull" = NA,
                       "S3>S1.group5:sinFull" = NA,
                       "S3>S1.group6:sinFull" = NA ))
hmm_new <- HMM$new(obs = obs, hid = hid, init = hmm, fixpar = fixpar)
hmm_new$fit(itnmax = 10000, control = list(eval.max = 10000))

saveRDS(hmm_new, "outputs/3stateHMM_restrict_group&tod_int.rds")

hmm_new$out()
hmm_new$confint()

hmm_new$AIC_marginal()
hmm$AIC_marginal()



##### covariates: random effects, interactions of group & tod #####

# previous HMM
hmm <- readRDS("outputs/3stateHMM_restrict_group&tod_int.rds")

f <- ~ group * cosFull + group * sinFull  + s(Fish_id, bs = "re")
tpm0 <- matrix(c(0.96, 0.04, 0,
                 0.02, 0.96, 0.02,
                 0, 0.06, 0.94),
               nrow = 3, byrow = TRUE)
hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                       initial_state = "stationary", formula = f)
obs <- Observation$new(data = sticklebacks, n_states = 3,
                       dists = list(Loco_sum = "gamma2"),
                       par = list(Loco_sum = list(mean = c(0.4, 1, 2.5),
                                                  sd = c(0.4, 0.8, 2))) )
fixpar <- list(hid = c("S1>S3.(Intercept)" = NA,
                       "S1>S3.group2" = NA,
                       "S1>S3.group3" = NA,
                       "S1>S3.group4" = NA,
                       "S1>S3.group5" = NA,
                       "S1>S3.group6" = NA,
                       "S1>S3.cosFull" = NA,
                       "S1>S3.sinFull" = NA,
                       "S1>S3.group2:cosFull" = NA,
                       "S1>S3.group3:cosFull" = NA,
                       "S1>S3.group4:cosFull" = NA,
                       "S1>S3.group5:cosFull" = NA,
                       "S1>S3.group6:cosFull" = NA,
                       "S1>S3.group2:sinFull" = NA,
                       "S1>S3.group3:sinFull" = NA,
                       "S1>S3.group4:sinFull" = NA,
                       "S1>S3.group5:sinFull" = NA,
                       "S1>S3.group6:sinFull" = NA,
                       "S3>S1.(Intercept)" = NA,
                       "S3>S1.group2" = NA,
                       "S3>S1.group3" = NA,
                       "S3>S1.group4" = NA,
                       "S3>S1.group5" = NA,
                       "S3>S1.group6" = NA,
                       "S3>S1.cosFull" = NA,
                       "S3>S1.sinFull" = NA,
                       "S3>S1.group2:cosFull" = NA,
                       "S3>S1.group3:cosFull" = NA,
                       "S3>S1.group4:cosFull" = NA,
                       "S3>S1.group5:cosFull" = NA,
                       "S3>S1.group6:cosFull" = NA,
                       "S3>S1.group2:sinFull" = NA,
                       "S3>S1.group3:sinFull" = NA,
                       "S3>S1.group4:sinFull" = NA,
                       "S3>S1.group5:sinFull" = NA,
                       "S3>S1.group6:sinFull" = NA ))
hmm_new <- HMM$new(obs = obs, hid = hid, init = hmm, fixpar = fixpar)
hmm_new$fit(itnmax = 10000, control = list(eval.max = 10000))

saveRDS(hmm_new, "outputs/3stateHMM_restrict_group&tod_int_reFishID.rds")

hmm_new$out()
hmm_new$confint()

hmm_new$AIC_marginal()
hmm$AIC_marginal()



##### covariates: final model + additional trigo functions #####

# previous HMM
hmm <- readRDS("outputs/3stateHMM_restrict_group&tod_int_reFishID.rds")

f <- ~ group * cosFull + group * sinFull + 
  group * cosHalf + group * sinHalf + s(Fish_id, bs = "re")
tpm0 <- matrix(c(0.96, 0.04, 0,
                 0.02, 0.96, 0.02,
                 0, 0.06, 0.94),
               nrow = 3, byrow = TRUE)
hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                       initial_state = "stationary", formula = f)
obs <- Observation$new(data = sticklebacks, n_states = 3,
                       dists = list(Loco_sum = "gamma2"),
                       par = list(Loco_sum = list(mean = c(0.4, 1, 2.5),
                                                  sd = c(0.4, 0.8, 2))) )
fixpar <- list(hid = c("S1>S3.(Intercept)" = NA,
                       "S1>S3.group2" = NA,
                       "S1>S3.group3" = NA,
                       "S1>S3.group4" = NA,
                       "S1>S3.group5" = NA,
                       "S1>S3.group6" = NA,
                       "S1>S3.cosFull" = NA,
                       "S1>S3.sinFull" = NA,
                       "S1>S3.cosHalf" = NA,
                       "S1>S3.sinHalf" = NA,
                       "S1>S3.group2:cosFull" = NA,
                       "S1>S3.group3:cosFull" = NA,
                       "S1>S3.group4:cosFull" = NA,
                       "S1>S3.group5:cosFull" = NA,
                       "S1>S3.group6:cosFull" = NA,
                       "S1>S3.group2:sinFull" = NA,
                       "S1>S3.group3:sinFull" = NA,
                       "S1>S3.group4:sinFull" = NA,
                       "S1>S3.group5:sinFull" = NA,
                       "S1>S3.group6:sinFull" = NA,
                       "S1>S3.group2:cosHalf" = NA,
                       "S1>S3.group3:cosHalf" = NA,
                       "S1>S3.group4:cosHalf" = NA,
                       "S1>S3.group5:cosHalf" = NA,
                       "S1>S3.group6:cosHalf" = NA,
                       "S1>S3.group2:sinHalf" = NA,
                       "S1>S3.group3:sinHalf" = NA,
                       "S1>S3.group4:sinHalf" = NA,
                       "S1>S3.group5:sinHalf" = NA,
                       "S1>S3.group6:sinHalf" = NA,
                       "S3>S1.(Intercept)" = NA,
                       "S3>S1.group2" = NA,
                       "S3>S1.group3" = NA,
                       "S3>S1.group4" = NA,
                       "S3>S1.group5" = NA,
                       "S3>S1.group6" = NA,
                       "S3>S1.cosFull" = NA,
                       "S3>S1.sinFull" = NA,
                       "S3>S1.cosHalf" = NA,
                       "S3>S1.sinHalf" = NA,
                       "S3>S1.group2:cosFull" = NA,
                       "S3>S1.group3:cosFull" = NA,
                       "S3>S1.group4:cosFull" = NA,
                       "S3>S1.group5:cosFull" = NA,
                       "S3>S1.group6:cosFull" = NA,
                       "S3>S1.group2:sinFull" = NA,
                       "S3>S1.group3:sinFull" = NA,
                       "S3>S1.group4:sinFull" = NA,
                       "S3>S1.group5:sinFull" = NA,
                       "S3>S1.group6:sinFull" = NA,
                       "S3>S1.group2:cosHalf" = NA,
                       "S3>S1.group3:cosHalf" = NA,
                       "S3>S1.group4:cosHalf" = NA,
                       "S3>S1.group5:cosHalf" = NA,
                       "S3>S1.group6:cosHalf" = NA,
                       "S3>S1.group2:sinHalf" = NA,
                       "S3>S1.group3:sinHalf" = NA,
                       "S3>S1.group4:sinHalf" = NA,
                       "S3>S1.group5:sinHalf" = NA,
                       "S3>S1.group6:sinHalf" = NA ))
hmm_new <- HMM$new(obs = obs, hid = hid, init = hmm, fixpar = fixpar)
hmm_new$fit(itnmax = 10000, control = list(eval.max = 10000))

saveRDS(hmm_new, "outputs/3stateHMM_restrict_additionalTrigo.rds")

hmm_new$out()
hmm_new$confint()

hmm_new$AIC_marginal()
hmm$AIC_marginal()



##### random starting values #####

runs <- 20
aic_alt <- hmm$AIC_marginal()
allAIC <- numeric(runs)
for (k in 1:runs){
  tpm0 <- matrix(c(0.96, 0.04, 0,
                   0.02, 0.96, 0.02,
                   0, 0.06, 0.94),
                 nrow = 3, byrow = TRUE)
  hid <- MarkovChain$new(data = sticklebacks, n_states = 3, tpm = tpm0,
                         initial_state = "stationary", formula = f)
  coeff_fe <- hid$coeff_fe()
  coeff_fe[which(coeff_fe != 0 & coeff_fe != -Inf)] <- runif(4, -4, -1.5)
  coeff_fe[c(2:18, 38:54, 56:72, 92:108)] <- runif(68, -1.5, 1.5)
  hid$update_coeff_fe(coeff_fe = coeff_fe)
  par0 <- list(Loco_sum = list(mean = c(rnorm(1, 0.4, 0.01),
                                        rnorm(1, 1, 0.05),
                                        rnorm(1, 2.5, 0.1)),
                               sd = c(rnorm(1, 0.4, 0.05),
                                      rnorm(1, 0.8, 0.1),
                                      rnorm(1, 2, 0.2)) ))
  obs <- Observation$new(data = sticklebacks, n_states = 3,
                         dists = list(Loco_sum = "gamma2"),
                         par = par0)
  hmm <- HMM$new(obs = obs, hid = hid, fixpar = fixpar)
  hmm$fit(silent = TRUE, itnmax = 10000, control = list(eval.max = 10000))
  
  aic <- hmm$AIC_marginal()
  allAIC[k] <- aic
  # allAIC <- c(allAIC, aic)
  saveRDS(allAIC, "outputs/allAIC_group_reFishID.rds")
  saveRDS(hmm, paste0("outputs/3stateHMM_randomInitials_group_reFishID", k, ".rds"))
  # if(aic < aic_alt){saveRDS(hmm, paste0("outputs/3stateHMM_covsRI", k, ".rds"))
  #   # aic_alt <- aic
  # }
  gc()
}
summary(allAIC)
plot(allAIC)



##### tpm #####

hmm <- readRDS("outputs/3stateHMM_restrict_group&tod_int_reFishID.rds") # "3stateHMM_restrict_covs&tod_int_reFishID.rds"

tod <- 1:1440
newdata <- data.frame(tod = rep(tod, 6),
                      Late = rep(c(0, 1), each = 3 * length(tod)), 
                      Exposed = rep(rep(c(0, 1, 1), each = length(tod)), 2),
                      Infected = rep(rep(c(0, 0, 1), each = length(tod)), 2),
                      cosFull = rep(cos(2 * pi * tod / 1440), 6),
                      sinFull = rep(sin(2 * pi * tod / 1440), 6),
                      cosHalf = rep(cos(2 * pi * tod / 720), 6),
                      sinHalf = rep(sin(2 * pi * tod / 720), 6))

newdata <- newdata %>% mutate(group = case_when(
  Infected == 0 & Exposed == 0 & Late == 0 ~ 1,
  Infected == 0 & Exposed == 1 & Late == 0 ~ 2,
  Infected == 1 & Exposed == 1 & Late == 0 ~ 3,
  Infected == 0 & Exposed == 0 & Late == 1 ~ 4,
  Infected == 0 & Exposed == 1 & Late == 1 ~ 5,
  Infected == 1 & Exposed == 1 & Late == 1 ~ 6
))
newdata$group <- as.factor(newdata$group)
newdata$Fish_id <- "F"

gammas <- hmm$predict("tpm", newdata = newdata, n_post = 100) 

newdata$gamma11 <- gammas$mean[1, 1, ]
newdata$gamma22 <- gammas$mean[2, 2, ]
newdata$gamma33 <- gammas$mean[3, 3, ]
newdata$gamma12 <- gammas$mean[1, 2, ]
newdata$gamma13 <- gammas$mean[1, 3, ]
newdata$gamma21 <- gammas$mean[2, 1, ]
newdata$gamma23 <- gammas$mean[2, 3, ]
newdata$gamma31 <- gammas$mean[3, 1, ]
newdata$gamma32 <- gammas$mean[3, 2, ]
newdata$gamma11_lci <- gammas$lcl[1, 1, ]
newdata$gamma22_lci <- gammas$lcl[2, 2, ]
newdata$gamma33_lci <- gammas$lcl[3, 3, ]
newdata$gamma12_lci <- gammas$lcl[1, 2, ]
newdata$gamma13_lci <- gammas$lcl[1, 3, ]
newdata$gamma21_lci <- gammas$lcl[2, 1, ]
newdata$gamma23_lci <- gammas$lcl[2, 3, ]
newdata$gamma31_lci <- gammas$lcl[3, 1, ]
newdata$gamma32_lci <- gammas$lcl[3, 2, ]
newdata$gamma11_uci <- gammas$ucl[1, 1, ]
newdata$gamma22_uci <- gammas$ucl[2, 2, ]
newdata$gamma33_uci <- gammas$ucl[3, 3, ]
newdata$gamma12_uci <- gammas$ucl[1, 2, ]
newdata$gamma13_uci <- gammas$ucl[1, 3, ]
newdata$gamma21_uci <- gammas$ucl[2, 1, ]
newdata$gamma23_uci <- gammas$ucl[2, 3, ]
newdata$gamma31_uci <- gammas$ucl[3, 1, ]
newdata$gamma32_uci <- gammas$ucl[3, 2, ]

pBase <- ggplot(filter(newdata, group == 4 | group == 5 | group == 6),
                aes(x = tod, group = group, color = group, fill = group)) + 
  scale_x_continuous(name = "time of day (hours)",
                     breaks = seq(0, 1440, by = 120),
                     labels = seq(0, 24, by = 2)) + 
  # scale_y_continuous(limits = c(0, 1),
  #                    breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(name = "condition", values = cbPalette[c(1, 6, 7)],
                     labels = c("control", "exposed", "infected")) +
  scale_fill_manual(name = "condition", values = cbPalette[c(1, 6, 7)],
                    labels = c("control", "exposed", "infected"))

p11 <- pBase + geom_line(aes(y = gamma11)) +
  geom_ribbon(aes(ymin = gamma11_lci, ymax = gamma11_uci), alpha = 0.2)
p22 <- pBase + geom_line(aes(y = gamma22)) +
  geom_ribbon(aes(ymin = gamma22_lci, ymax = gamma22_uci), alpha = 0.2)
p33 <- pBase + geom_line(aes(y = gamma33)) +
  geom_ribbon(aes(ymin = gamma33_lci, ymax = gamma33_uci), alpha = 0.2)
p12 <- pBase + geom_line(aes(y = gamma12)) +
  geom_ribbon(aes(ymin = gamma12_lci, ymax = gamma12_uci), alpha = 0.2)
p13 <- pBase + geom_line(aes(y = gamma13)) +
  geom_ribbon(aes(ymin = gamma13_lci, ymax = gamma13_uci), alpha = 0.2)
p21 <- pBase + geom_line(aes(y = gamma21)) +
  geom_ribbon(aes(ymin = gamma21_lci, ymax = gamma21_uci), alpha = 0.2)
p23 <- pBase + geom_line(aes(y = gamma23)) +
  geom_ribbon(aes(ymin = gamma23_lci, ymax = gamma23_uci), alpha = 0.2)
p31 <- pBase + geom_line(aes(y = gamma31)) +
  geom_ribbon(aes(ymin = gamma31_lci, ymax = gamma31_uci), alpha = 0.2)
p32 <- pBase + geom_line(aes(y = gamma32)) +
  geom_ribbon(aes(ymin = gamma32_lci, ymax = gamma32_uci), alpha = 0.2)

p_all <- ggarrange(p11, p12, p13, p21, p22, p23, p31, p32, p33, nrow = 3, ncol = 3,
                   common.legend = TRUE, legend="bottom")
p_all

pdf("tpm_group456_group&tod_final.pdf", width = 12, height = 8)
p_all
dev.off()



##### data frame with proportions of decoded states per fish #####

sticklebacks$states <- hmm$viterbi()
df_ID <- sticklebacks %>% group_by(ID) %>% 
  summarise(p_state1 = sum(states == 1) / n(), 
            p_state2 = sum(states == 2) / n(), 
            p_state3 = sum(states == 3) / n())
rowSums(df_ID[, 2:4], na.rm = TRUE)

write.csv(df_ID, "outputs/df_stateProportions.csv", row.names = FALSE)



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

statProbs_covs <- function(beta, Late = 0, Exposed = 0, Infected = 0, 
                           tod = 1:1440){
  beta12 <- beta[1:12]
  beta13 <- beta[13:24]
  beta21 <- beta[25:36]
  beta23 <- beta[37:48]
  beta31 <- beta[49:60]
  beta32 <- beta[61:72]
  n <- length(tod)
  N <- 3
  covsMat <- data.frame(Intercept = 1, 
                        Late1 = Late,
                        cosFull = cos(2 * pi * tod / n),
                        Exposed1 = Exposed,
                        Infected1 = Infected, 
                        sinFull = sin(2 * pi * tod / n),
                        cosFull_Late = cos(2 * pi * tod / n) * Late,
                        cosFull_Exposed = cos(2 * pi * tod / n) * Exposed,
                        cosFull_Infected = cos(2 * pi * tod / n) * Infected,
                        sinFull_Late = sin(2 * pi * tod / n) * Late,
                        sinFull_Exposed = sin(2 * pi * tod / n) * Exposed,
                        sinFull_Infected = sin(2 * pi * tod / n) * Infected)
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
  # return(Gammas)
}


statProbs_group <- function(beta, tod = 1:1440, group = 1, simple = "yes"){
  group2 = ifelse(group == 2, 1, 0)
  group3 = ifelse(group == 3, 1, 0)
  group4 = ifelse(group == 4, 1, 0)
  group5 = ifelse(group == 5, 1, 0)
  group6 = ifelse(group == 6, 1, 0)
  n <- length(tod)
  N <- 3
  if(simple == "yes"){
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
  } else{
    beta12 <- beta[1:30]
    beta13 <- beta[31:60]
    beta21 <- beta[61:90]
    beta23 <- beta[91:120]
    beta31 <- beta[121:150]
    beta32 <- beta[151:180]
    covsMat <- data.frame(Intercept = 1, 
                          group2 = group2,
                          group3 = group3,
                          group4 = group4,
                          group5 = group5,
                          group6 = group6,
                          cosFull = cos(2 * pi * tod / n),
                          sinFull = sin(2 * pi * tod / n),
                          cosHalf = cos(2 * pi * tod / (n / 2)),
                          sinHalf = sin(2 * pi * tod / (n / 2)),
                          int2_cosFull = cos(2 * pi * tod / n) * group2,
                          int3_cosFull = cos(2 * pi * tod / n) * group3,
                          int4_cosFull = cos(2 * pi * tod / n) * group4,
                          int5_cosFull = cos(2 * pi * tod / n) * group5,
                          int6_cosFull = cos(2 * pi * tod / n) * group6,
                          int2_sinFull = sin(2 * pi * tod / n) * group2,
                          int3_sinFull = sin(2 * pi * tod / n) * group3,
                          int4_sinFull = sin(2 * pi * tod / n) * group4,
                          int5_sinFull = sin(2 * pi * tod / n) * group5,
                          int6_sinFull = sin(2 * pi * tod / n) * group6,
                          int2_cosHalf = cos(2 * pi * tod / (n / 2)) * group2,
                          int3_cosHalf = cos(2 * pi * tod / (n / 2)) * group3,
                          int4_cosHalf = cos(2 * pi * tod / (n / 2)) * group4,
                          int5_cosHalf = cos(2 * pi * tod / (n / 2)) * group5,
                          int6_cosHalf = cos(2 * pi * tod / (n / 2)) * group6,
                          int2_sinHalf = sin(2 * pi * tod / (n / 2)) * group2,
                          int3_sinHalf = sin(2 * pi * tod / (n / 2)) * group3,
                          int4_sinHalf = sin(2 * pi * tod / (n / 2)) * group4,
                          int5_sinHalf = sin(2 * pi * tod / (n / 2)) * group5,
                          int6_sinHalf = sin(2 * pi * tod / (n / 2)) * group6)
  }
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



tod_reFish <- function(beta, reInt, Fish_id, tod = 1:1440, group = 1){
  int <- reInt$int[which(reInt$Fish_id == Fish_id)]
  beta12 <- c(int[1], beta[1:18])
  beta13 <- c(int[2], beta[19:36])
  beta21 <- c(int[3], beta[37:54])
  beta23 <- c(int[4], beta[55:72])
  beta31 <- c(int[5], beta[73:90])
  beta32 <- c(int[6], beta[91:108])
  group2 = ifelse(group == 2, 1, 0)
  group3 = ifelse(group == 3, 1, 0)
  group4 = ifelse(group == 4, 1, 0)
  group5 = ifelse(group == 5, 1, 0)
  group6 = ifelse(group == 6, 1, 0)
  n <- length(tod)
  N <- 3
  covsMat <- data.frame(reIntercept = 1, 
                        Intercept = 1, 
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
  return(delta) # list(delta = delta, gammas = Gammas)
}



##### state probabilities with CIs as a function of tod (new method) #####

draws <- 100 # 200 
betaDraws <- hmm$post_coeff(n_post = draws)
betaDraws <- betaDraws[, which(colnames(betaDraws) == "coeff_fe_hid")]

newdata <- data.frame(Late = rep(c(0, 1), each = 3), 
                      Exposed = rep(c(0, 1, 1), 2),
                      Infected = rep(c(0, 0, 1), 2) )

df_all <- data.frame()
for (j in 1:6) { 
  print(j)
  df <- data.frame(tod = 1:1440, state1 = NA, state2 = NA, state3 = NA,
                   lci1 = NA, lci2 = NA, lci3 = NA, 
                   uci1 = NA, uci2 = NA, uci3 = NA, group = j)
  beta <- as.vector(hmm$coeff_fe()$hid)
  # df[, 2:4] <- statProbs_covs(beta, newdata$Late[j], newdata$Exposed[j],
  #                             newdata$Infected[j])
  df[, 2:4] <- statProbs_group(beta, group = j, simple = "no")
  deltaDraws <- array(NA, dim = c(1440, 3, draws)) 
  for (i in 1:draws) {
    beta <- as.vector(betaDraws[i, ])
    # deltaDraws[, , i] <- statProbs_covs(beta, newdata$Late[j], newdata$Exposed[j],
    #                                     newdata$Infected[j])
    deltaDraws[, , i] <- statProbs_group(beta, group = j, simple = "no")
  }
  for (i in 1:1440) {
    df$lci1[i] <- quantile(deltaDraws[i, 1, ], 0.025, na.rm = TRUE)
    df$uci1[i] <- quantile(deltaDraws[i, 1, ], 0.975, na.rm = TRUE)
    df$lci2[i] <- quantile(deltaDraws[i, 2, ], 0.025, na.rm = TRUE)
    df$uci2[i] <- quantile(deltaDraws[i, 2, ], 0.975, na.rm = TRUE)
    df$lci3[i] <- quantile(deltaDraws[i, 3, ], 0.025, na.rm = TRUE)
    df$uci3[i] <- quantile(deltaDraws[i, 3, ], 0.975, na.rm = TRUE)
  }
  df_all <- rbind(df_all, df)
}
head(df_all)
write.csv(df_all, "outputs/stateProbs_additionalTrigo.csv", row.names = FALSE)


df_all$Early <- ifelse(df_all$group == 1 | df_all$group == 2 | 
                          df_all$group == 3, 1, 0)
df_all$group <- as.factor(df_all$group)

pBase <- ggplot(mapping = aes(x = tod, group = group, 
                              color = group, fill = group)) + 
  scale_x_continuous(name = "time of day (hours)",
                     breaks = seq(0, 1440, by = 120),
                     labels = seq(0, 24, by = 2)) + 
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(name = "condition", values = cbPalette[c(1, 6, 7)],
                     labels = c("control", "exposed", "infected")) +
  scale_fill_manual(name = "condition", values = cbPalette[c(1, 6, 7)],
                    labels = c("control", "exposed", "infected"))

pDelta1_1 <- pBase + geom_line(aes(y = state1), linewidth = 0.8,
                               data = filter(df_all, Early == 1)) + 
  geom_ribbon(aes(ymin = lci1, ymax = uci1), alpha = 0.2,
              data = filter(df_all, Early == 1)) +
  ylab("probability of state 1")

pDelta2_1 <- pBase + geom_line(aes(y = state2), linewidth = 0.8,
                               data = filter(df_all, Early == 1)) + 
  geom_ribbon(aes(ymin = lci2, ymax = uci2), alpha = 0.2,
              data = filter(df_all, Early == 1)) +
  ylab("probability of state 2")

pDelta3_1 <- pBase + geom_line(aes(y = state3), linewidth = 0.8,
                               data = filter(df_all, Early == 1)) + 
  geom_ribbon(aes(ymin = lci3, ymax = uci3), alpha = 0.2,
              data = filter(df_all, Early == 1)) +
  ylab("probability of state 3")

pDelta1_0 <- pBase + geom_line(aes(y = state1), linewidth = 0.8,
                               data = filter(df_all, Early == 0)) + 
  geom_ribbon(aes(ymin = lci1, ymax = uci1), alpha = 0.2,
              data = filter(df_all, Early == 0)) +
  ylab("probability of state 1")

pDelta2_0 <- pBase + geom_line(aes(y = state2), linewidth = 0.8,
                               data = filter(df_all, Early == 0)) + 
  geom_ribbon(aes(ymin = lci2, ymax = uci2), alpha = 0.2,
              data = filter(df_all, Early == 0)) +
  ylab("probability of state 2")

pDelta3_0 <- pBase + geom_line(aes(y = state3), linewidth = 0.8,
                               data = filter(df_all, Early == 0)) + 
  geom_ribbon(aes(ymin = lci3, ymax = uci3), alpha = 0.2,
              data = filter(df_all, Early == 0)) +
  ylab("probability of state 3")

pDelta <- ggarrange(pDelta1_1, pDelta2_1, pDelta3_1, 
                    pDelta1_0, pDelta2_0, pDelta3_0, nrow = 2, ncol = 3,
                    common.legend = TRUE, legend="bottom")
pDelta

ggsave("delta_additionalTrigo.pdf", pDelta, width = 12, height = 8)



### for each fish (reID)

beta <- as.vector(hmm$coeff_fe()$hid)
reInt <- data.frame(Fish_id = rep(levels(sticklebacks$Fish_id), 6),
                    int = as.vector(hmm$coeff_re()$hid))

newdata <- sticklebacks %>% select(Fish_id, group) %>% distinct()
tod <- rep(1:1440, nrow(newdata))
df <- data.frame(tod = tod,
                 Fish_id = as.factor(rep(newdata$Fish_id, each = 1440)),
                 group = as.factor(rep(newdata$group, each = 1440)),
                 state1 = NA, state2 = NA, state3 = NA)

for (i in 1:nrow(newdata)) {
  print(i)
  delta <- tod_reFish(beta, reInt, Fish_id = newdata$Fish_id[i],
                      group = as.numeric(newdata$group[i]),)
  df[(1440 * (i - 1) + 1):(1440 * i), 4:6] <- delta
}

write.csv(df, "outputs/stateProbs_final_reID.csv", row.names = FALSE)



##### expected state dwell times #####

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


### expected dwell time

# tod <- 1:1440
# newdata <- data.frame(tod = rep(tod, 6),
#                       Late = rep(c(0, 1), each = 3 * length(tod)), 
#                       Exposed = rep(rep(c(0, 1, 1), each = length(tod)), 2),
#                       Infected = rep(rep(c(0, 0, 1), each = length(tod)), 2),
#                       cosFull = rep(cos(2 * pi * tod / 1440), 6),
#                       sinFull = rep(sin(2 * pi * tod / 1440), 6),
#                       cosHalf = rep(cos(2 * pi * tod / 720), 6),
#                       sinHalf = rep(sin(2 * pi * tod / 720), 6))
# 
# newdata <- newdata %>% mutate(group = case_when(
#   Infected == 0 & Exposed == 0 & Late == 0 ~ 1,
#   Infected == 0 & Exposed == 1 & Late == 0 ~ 2,
#   Infected == 1 & Exposed == 1 & Late == 0 ~ 3,
#   Infected == 0 & Exposed == 0 & Late == 1 ~ 4,
#   Infected == 0 & Exposed == 1 & Late == 1 ~ 5,
#   Infected == 1 & Exposed == 1 & Late == 1 ~ 6
# ))
# newdata$group <- as.factor(newdata$group)
# newdata$Fish_id <- "F"
# 
# gammas <- hmm$predict("tpm", newdata = newdata) 
# 
# L <- 1440
# r <- 1:L
# newdata$dwellTime1 <- 0
# newdata$dwellTime2 <- 0
# newdata$dwellTime3 <- 0
# for (i in 1:6) {
#   print(i)
#   gIndex <- which(newdata$group == i)
#   for (l in 1:L) {
#     tDwell <- ddwell_t_new(r, t = l, state = 1, Gamma = gammas[, , gIndex])
#     newdata$dwellTime1[gIndex[l]] <- (L + sum(r * tDwell)) / sum(tDwell) - L
#   }
#   for (l in 1:L) {
#     tDwell <- ddwell_t_new(r, t = l, state = 2, Gamma = gammas[, , gIndex])
#     newdata$dwellTime2[gIndex[l]] <- (L + sum(r * tDwell)) / sum(tDwell) - L
#   }
#   for (l in 1:L) {
#     tDwell <- ddwell_t_new(r, t = l, state = 3, Gamma = gammas[, , gIndex])
#     newdata$dwellTime3[gIndex[l]] <- (L + sum(r * tDwell)) / sum(tDwell) - L
#   }
# }
# write.csv(newdata, "expDwellTimes.csv", row.names = FALSE)


### expected dwell time & CIs

draws <- 100 # 200 
betaDraws <- hmm$post_coeff(n_post = draws)
betaDraws <- betaDraws[, which(colnames(betaDraws) == "coeff_fe_hid")]

newdata <- data.frame(Late = rep(c(0, 1), each = 3), 
                      Exposed = rep(c(0, 1, 1), 2),
                      Infected = rep(c(0, 0, 1), 2) )

df_all <- data.frame()
L <- 1440
r <- 1:L
for (j in 1:6) { 
  print(j)
  df <- data.frame(tod = 1:1440, dstate1 = NA, dstate2 = NA, dstate3 = NA,
                   lci1 = NA, lci2 = NA, lci3 = NA, 
                   uci1 = NA, uci2 = NA, uci3 = NA, group = j)
  beta <- as.vector(hmm$coeff_fe()$hid)
  gammas <- getGammas(beta, group = j)
  for (l in 1:L) {
    tDwell <- ddwell_t_new(r, t = l, state = 1, Gamma = gammas)
    df$dstate1[l] <- (L + sum(r * tDwell)) / sum(tDwell) - L
  }
  for (l in 1:L) {
    tDwell <- ddwell_t_new(r, t = l, state = 2, Gamma = gammas)
    df$dstate2[l] <- (L + sum(r * tDwell)) / sum(tDwell) - L
  }
  for (l in 1:L) {
    tDwell <- ddwell_t_new(r, t = l, state = 3, Gamma = gammas)
    df$dstate3[l] <- (L + sum(r * tDwell)) / sum(tDwell) - L
  }

  dwellDraws <- array(NA, dim = c(1440, 3, draws)) 
  for (i in 1:draws) {
    beta <- as.vector(betaDraws[i, ])
    gammas <- getGammas(beta, group = j)
    for (l in 1:L) {
      tDwell <- ddwell_t_new(r, t = l, state = 1, Gamma = gammas)
      dwellDraws[l, 1, i] <- (L + sum(r * tDwell)) / sum(tDwell) - L
    }
    for (l in 1:L) {
      tDwell <- ddwell_t_new(r, t = l, state = 2, Gamma = gammas)
      dwellDraws[l, 2, i] <- (L + sum(r * tDwell)) / sum(tDwell) - L
    }
    for (l in 1:L) {
      tDwell <- ddwell_t_new(r, t = l, state = 3, Gamma = gammas)
      dwellDraws[l, 3, i] <- (L + sum(r * tDwell)) / sum(tDwell) - L
    }
  }
  for (i in 1:1440) {
    df$lci1[i] <- quantile(dwellDraws[i, 1, ], 0.025, na.rm = TRUE)
    df$uci1[i] <- quantile(dwellDraws[i, 1, ], 0.975, na.rm = TRUE)
    df$lci2[i] <- quantile(dwellDraws[i, 2, ], 0.025, na.rm = TRUE)
    df$uci2[i] <- quantile(dwellDraws[i, 2, ], 0.975, na.rm = TRUE)
    df$lci3[i] <- quantile(dwellDraws[i, 3, ], 0.025, na.rm = TRUE)
    df$uci3[i] <- quantile(dwellDraws[i, 3, ], 0.975, na.rm = TRUE)
  }
  df_all <- rbind(df_all, df)
}
head(df_all)
write.csv(df_all, "outputs/expDwellTimes_CIs.csv", row.names = FALSE)


### overall dwell-time distributions (+ empirical dwell times)

ddwell_new = function(
    r, # vector of dwell-times to compute, must be of the form 1:R
    state, # which state to compute
    Gamma # array of dim c(N,N,L)
){
  L = dim(Gamma)[3]
  N = dim(Gamma)[1]
  # calculate p-stationary
  delta = matrix(NA, L, N)
  GammaT=Gamma[,,1]
  for (k in 2:L){ GammaT=GammaT%*%Gamma[,,k] }
  delta[1,] = solve(t(diag(N)-GammaT+1), rep(1,N))
  for(k in 2:L){ delta[k,] = delta[k-1,]%*%Gamma[,,k-1] }
  # calculate weights (only for state of interest)
  weights=numeric(L)
  weights[1] = sum(delta[L,-state] * Gamma[-state,state, L])
  for (k in 2:L){ weights[k]=sum(delta[k-1,-state] * Gamma[-state,state, k-1]) }
  weights = weights/sum(weights)
  # calculate all weighted d_i^t's
  pmfs_weighted = matrix(NA, L, length(r))
  tDwell_weighted = numeric(L)
  for(t in 1:L){ 
    tDwell = ddwell_t_new(r=r, t=t, state=state, Gamma=Gamma) 
    pmfs_weighted[t,] = weights[t] * tDwell
    tDwell_weighted[t] = weights[t] * ((L + sum(r * tDwell)) / sum(tDwell))
  }
  # overall dwell-time dist.
  pmf = apply(pmfs_weighted, 2, sum)
  # mean overall dwell time 
  meanOverall = (sum(tDwell_weighted) - L) 
  return(list(pmf = pmf, meanOverall = meanOverall))
}

L <- 1440
r <- 1:L
df_all <- data.frame()
beta <- as.vector(hmm$coeff_fe()$hid)
for (j in 1:6) {
  df <- data.frame(tod = 1:1440, dwellDist1 = NA, meanDwell1 = NA, group = j)
  gammas <- getGammas(beta, group = j)
  
  overall <- ddwell_new(r, state = 1, Gamma = gammas)
  df$dwellDist1 <- overall$pmf
  df$meanDwell1 <- rep(overall$meanOverall, L)
  df_all <- rbind(df_all, df)
}

df_all$period <- c(rep("Early", 3*L), rep("Late", 3*L))
df_all$treatment <- rep(c(rep("Control", L), rep("Exposed", L),
                          rep("Infected", L)), 2)
dwellMean <- df_all %>% select("meanDwell1", "period", "treatment") %>% 
  distinct()
ggplot(df_all, aes(tod, dwellDist1)) + 
  geom_bar(stat="identity", color = pal[1], fill = pal[1]) +
  geom_text(data = dwellMean,
            mapping = aes(x = 40, y = 0.125,
                          label = paste0("mean = ", round(meanDwell1, 1)) )) +
  labs(y = "overall dwell-time dist. in state 1",
       x = "time (in minutes)") +
  scale_x_continuous(limits = c(0, 60)) + 
  facet_grid(period ~ treatment) 

# empirical dwell times

ids <- unique(sticklebacks$ID)
df_all <- data.frame()
for (i in 1:length(ids)) {
  id <- unique(sticklebacks$ID)[i]
  rows <- which(sticklebacks$ID == id)
  x <- which(sticklebacks[rows, ]$states == 1)
  runs <- split(x, cumsum(c(0, diff(x) > 1))) # seq_along(x)
  dwellID <- as.vector(sapply(runs, length))
  df <- data.frame(ID = id, Fish_id = unique(sticklebacks$Fish_id[rows]),
                   Exposed = unique(sticklebacks$Exposed[rows]),
                   Infected = unique(sticklebacks$Infected[rows]),
                   Late = unique(sticklebacks$Late[rows]),
                   empDwell = dwellID, dwellMean = mean(dwellID))
  df_all <- rbind(df_all, df)
}
df_all <- df_all %>% mutate(period = ifelse(Late == 0, "Early", "Late"),
                            treatment = case_when(
                              Exposed == 0 & Infected == 0 ~ "Control",
                              Exposed == 1 & Infected == 0 ~ "Exposed",
                              Exposed == 1 & Infected == 1 ~ "Infected"
                            ))

df <- df_all[-which(df_all$empDwell == 0), ]
df_empDwell <- data.frame()
for (i in 1:60) {
  df_empDwell <- rbind(df_empDwell, df %>% group_by(period, treatment) %>%
                         summarise(y = sum(empDwell == i) / n(), x = i))
}
dwellMeans <- df_all %>% group_by(period, treatment) %>% 
  summarise(mean = mean(empDwell)) %>% ungroup() 

p_dwellDist <- ggplot(df_empDwell, aes(x, y)) + 
  geom_bar(stat="identity", color = pal[1], fill = pal[1]) +
  geom_text(aes(x = 40, y = 0.11,
                label = paste0("mean = ", round(mean, 1)) ),
            dwellMeans) +
  labs(y = "empirical dwell-time dist. in state 1",
       x = "time (in minutes)") +
  facet_grid(period ~ treatment) 
ggsave("fig_empirical-dwell-time-dist.pdf", p_dwellDist,
       width = 10, height = 6)

df_all <- df_all %>% group_by(Fish_id, period, treatment, dwellMean) %>% 
  summarise(n = n(), overallTime = sum(empDwell)) %>% ungroup()
x <- which(df_all$treatment == "Infected" & df_all$period == "Late")
sum(df_all$n[x] * df_all$dwellMean[x]) / sum(df_all$n[x])

# in diesem Plot sieht man nicht die Gewichtung der individual dwell times 
# für den overall mean!!!
p_empDwells <- ggplot(df_all, aes(period, dwellMean, color = Fish_id, group = Fish_id)) + 
  geom_point() + geom_line() + 
  ylab("average dwell time (in min) in state 1") + 
  facet_wrap(~ treatment) +
  theme(legend.position = "none") 
ggsave("fig_individual-dwell-times.pdf", p_empDwells,
       width = 8, height = 4)

ggplot(df_all, aes(period, n, color = Fish_id, group = Fish_id)) + 
  geom_point() + geom_line() + 
  ylab("number of times in state 1") + 
  facet_wrap(~ treatment) +
  theme(legend.position = "none") 



##### model checking #####

gof <- function(data) {
  s <- c(quantile(data$Loco_sum, seq(0, 1, by = 0.25), na.rm = TRUE),
         autocor = cor(data$Loco_sum[-1], data$Loco_sum[-nrow(data)])) # use = "complete.obs"
}
# Run posterior predictive checks
# checks <- hmm$check(check_fn = gof, silent = FALSE, nsims = 100)
# saveRDS(checks, "outputs/modelChecks.rds")
checks <- readRDS("outputs/modelChecks.rds")
checks$plot

head(checks$obs_stat)
cor(sticklebacks$Loco_sum[-1], sticklebacks$Loco_sum[-nrow(sticklebacks)], 
    use = "complete.obs")
head(checks$stats)


# pr <- hmm$pseudores()
# saveRDS(pr, "outputs/pseudoResiduals.rds")
pr <- readRDS("outputs/pseudoResiduals.rds")
qqnorm(pr, ylim = c(-4, 4))
abline(0, 1, col = 2, lwd = 2)
hist(pr, breaks = 40, prob = TRUE, main = "Histogram of pseudo-residuals",
     ylab = "density", xlab = "pseudo-residuals")
curve(dnorm(x), col = 2, lwd = 2, add = TRUE)
acf(pr, na.action = na.pass, lag.max = 10)
