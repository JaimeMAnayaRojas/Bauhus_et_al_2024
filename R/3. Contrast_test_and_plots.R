# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(lubridate)
library(hmmTMB)
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
# test_state1 <- test[, , 1] # as a matrix or...
# test_state2 <- as.data.frame(test[, , 2]) # ...as a data frame
# test_state3 <- test[, , 3]

# Do state 1 contrast first
EE_state1 <- log(ExposedEarly[, , 1] / ControlEarly[, , 1])   # as a matrix
IE_state1 <- log(InfectedEarly[, , 1] / ControlEarly[, , 1])   # as a matrix


library(HDInterval)


# 

EC_s1 <- data.frame(medians = apply(ControlEarly[, , 1], 2, median),
            lower = t(apply(ControlEarly[, , 1], 2, hdi))[,1],
            upper = t(apply(ControlEarly[, , 1], 2, hdi))[,2]
            )

EC_s1$Infection <- "Control"
EC_s1$Time <- "Early"
EC_s1$tod <- 1:1440

EE_s1 <- data.frame(medians = apply(ExposedEarly[, , 1], 2, median),
            lower = t(apply(ExposedEarly[, , 1], 2, hdi))[,1],
            upper = t(apply(ExposedEarly[, , 1], 2, hdi))[,2]
            )

EE_s1$Infection <- "Exposed"
EE_s1$Time <- "Early"
EE_s1$tod <- 1:1440




EI_s1 <- data.frame(medians = apply(InfectedEarly[, , 1], 2, median),
            lower = t(apply(InfectedEarly[, , 1], 2, hdi))[,1],
            upper = t(apply(InfectedEarly[, , 1], 2, hdi))[,2]
            )

EI_s1$Infection <- "Infected"
EI_s1$Time <- "Early"
EI_s1$tod <- 1:1440

#---------------


EE_state1 <- data.frame(medians = apply(EE_state1, 2, median),
            lower = t(apply(EE_state1, 2, hdi))[,1],
            upper = t(apply(EE_state1, 2, hdi))[,2]
            )

EE_state1$Infection <- "Exposed"
EE_state1$Time <- "Early"
EE_state1$tod <- 1:1440

IE_state1 <- data.frame(medians = apply(IE_state1, 2, median),
            lower = t(apply(IE_state1, 2, hdi))[,1],
            upper = t(apply(IE_state1, 2, hdi))[,2]
            )


IE_state1$Infection <- "Infected"
IE_state1$Time <- "Early"
IE_state1$tod <- 1:1440


## Late

LC_s1 <- data.frame(medians = apply(ControlLate[, , 1], 2, median),
            lower = t(apply(ControlLate[, , 1], 2, hdi))[,1],
            upper = t(apply(ControlLate[, , 1], 2, hdi))[,2]
            )

LC_s1$Infection <- "Control"
LC_s1$Time <- "Late"
LC_s1$tod <- 1:1440



LE_s1 <- data.frame(medians = apply(ExposedLate[, , 1], 2, median),
            lower = t(apply(ExposedLate[, , 1], 2, hdi))[,1],
            upper = t(apply(ExposedLate[, , 1], 2, hdi))[,2]
            )

LE_s1$Infection <- "Exposed"
LE_s1$Time <- "Late"
LE_s1$tod <- 1:1440


LI_s1 <- data.frame(medians = apply(InfectedLate[, , 1], 2, median),
            lower = t(apply(InfectedLate[, , 1], 2, hdi))[,1],
            upper = t(apply(InfectedLate[, , 1], 2, hdi))[,2]
            )

LI_s1$Infection <- "Infected"
LI_s1$Time <- "Late"
LI_s1$tod <- 1:1440


# put them together

dfs <- as.data.frame(rbind(EC_s1, EE_s1, EI_s1, LC_s1, LE_s1, LI_s1))

# Plot the results

dfs$medians <- dfs$medians * 100
dfs$lower <- dfs$lower * 100
dfs$upper <- dfs$upper * 100

library(ggplot2)
head(df)
p2c <- ggplot(dfs, aes(x = tod, y = medians, color = Infection)) +
  geom_line(linewidth = 1.25) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Infection), alpha = 0.3, colour =NA) +
  scale_fill_manual(values = c("black","#78baed", "#D55E00")) +
  scale_color_manual(values = c("black","#0072B2", "#D55E00")) + # use theme bw but remove grids 
    
    
    # scale x axis by dividing by 60
    scale_x_continuous(breaks = seq(0, 1440, 120), labels = seq(0, 24, 2)) + 
    facet_grid(Time ~ .) 

p2c <- p2c + ylab("Probability of state 1 (%)") + xlab("Time of day (Hours)") + 
    theme(legend.position = c(0.6,0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("plots/EE&IE_state1c.pdf", p2c, width = 15, height = 15, units = "cm")



#---


EL_state1 <- log(ExposedLate[, , 1] / ControlLate[, , 1])   # as a matrix
IL_state1 <- log(InfectedLate[, , 1] / ControlLate[, , 1])   # as a matrix


EL_state1 <- data.frame(medians = apply(EL_state1, 2, median),
            lower = t(apply(EL_state1, 2, hdi))[,1],
            upper = t(apply(EL_state1, 2, hdi))[,2]
            )

EL_state1$Infection <- "Exposed"
EL_state1$Time <- "Late"
EL_state1$tod <- 1:1440

IL_state1 <- data.frame(medians = apply(IL_state1, 2, median),
            lower = t(apply(IL_state1, 2, hdi))[,1],
            upper = t(apply(IL_state1, 2, hdi))[,2]
            )


IL_state1$Infection <- "Infected"
IL_state1$Time <- "Late"
IL_state1$tod <- 1:1440

df <- as.data.frame(rbind(EE_state1, IE_state1, EL_state1, IL_state1))
head(df)

# Plot the results

library(ggplot2)
head(df)
p2d <- ggplot(df, aes(x = tod, y = medians, color = Infection)) +
  geom_line(linewidth = 1.25) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Infection), alpha = 0.3, colour =NA) +
  scale_fill_manual(values = c("#78baed", "#D55E00")) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  geom_hline(yintercept = 0, linetype = "dashed") + # scale x axis by dividing by 60
  scale_x_continuous(breaks = seq(0, 1440, 120), labels = seq(0, 24, 2)) + 
  facet_grid(Time ~ .) 

p2d <- p2d + ylab("log-ratio [ln(Treatment/Control)]") + xlab("Time of day (Hours)") + 
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"))

p2d

ggsave("plots/EE&IE_state1.pdf",p2d, width = 15, height = 15, units = "cm")

# use ggarange to put them together

library(ggpubr)

theme_set(theme_bw())
figure <- ggarrange(p2c, p2d,
                    labels = c("D", "F"),
                    ncol = 2, nrow = 1)
figure

ggsave("plots/EE&IE_state1D.pdf",figure, width = 20, height = 15, units = "cm")

### Do the same for state 2 ### ### ### ### ### ### 





