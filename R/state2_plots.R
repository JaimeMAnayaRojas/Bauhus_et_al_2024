
# Do state 1 contrast first
EE_state2 <- log(ExposedEarly[, , 2] / ControlEarly[, , 2])   # as a matrix
IE_state2 <- log(InfectedEarly[, , 2] / ControlEarly[, , 2])   # as a matrix


# 

EC_s2 <- data.frame(medians = apply(ControlEarly[, , 2], 2, median),
            lower = t(apply(ControlEarly[, , 2], 2, hdi))[,1],
            upper = t(apply(ControlEarly[, , 2], 2, hdi))[,2]
            )

EC_s2$Infection <- "Control"
EC_s2$Time <- "Early"
EC_s2$tod <- 1:1440

EE_s2 <- data.frame(medians = apply(ExposedEarly[, , 2], 2, median),
            lower = t(apply(ExposedEarly[, , 2], 2, hdi))[,1],
            upper = t(apply(ExposedEarly[, , 2], 2, hdi))[,2]
            )

EE_s2$Infection <- "Exposed"
EE_s2$Time <- "Early"
EE_s2$tod <- 1:1440




EI_s2 <- data.frame(medians = apply(InfectedEarly[, , 2], 2, median),
            lower = t(apply(InfectedEarly[, , 2], 2, hdi))[,1],
            upper = t(apply(InfectedEarly[, , 2], 2, hdi))[,2]
            )

EI_s2$Infection <- "Infected"
EI_s2$Time <- "Early"
EI_s2$tod <- 1:1440

#---------------


EE_state2 <- data.frame(medians = apply(EE_state2, 2, median),
            lower = t(apply(EE_state2, 2, hdi))[,1],
            upper = t(apply(EE_state2, 2, hdi))[,2]
            )

EE_state2$Infection <- "Exposed"
EE_state2$Time <- "Early"
EE_state2$tod <- 1:1440

IE_state2 <- data.frame(medians = apply(IE_state2, 2, median),
            lower = t(apply(IE_state2, 2, hdi))[,1],
            upper = t(apply(IE_state2, 2, hdi))[,2]
            )


IE_state2$Infection <- "Infected"
IE_state2$Time <- "Early"
IE_state2$tod <- 1:1440


## Late

LC_s2 <- data.frame(medians = apply(ControlLate[, , 2], 2, median),
            lower = t(apply(ControlLate[, , 2], 2, hdi))[,1],
            upper = t(apply(ControlLate[, , 2], 2, hdi))[,2]
            )

LC_s2$Infection <- "Control"
LC_s2$Time <- "Late"
LC_s2$tod <- 1:1440



LE_s2 <- data.frame(medians = apply(ExposedLate[, , 2], 2, median),
            lower = t(apply(ExposedLate[, , 2], 2, hdi))[,1],
            upper = t(apply(ExposedLate[, , 2], 2, hdi))[,2]
            )

LE_s2$Infection <- "Exposed"
LE_s2$Time <- "Late"
LE_s2$tod <- 1:1440


LI_s2 <- data.frame(medians = apply(InfectedLate[, , 2], 2, median),
            lower = t(apply(InfectedLate[, , 2], 2, hdi))[,1],
            upper = t(apply(InfectedLate[, , 2], 2, hdi))[,2]
            )

LI_s2$Infection <- "Infected"
LI_s2$Time <- "Late"
LI_s2$tod <- 1:1440


# put them together

dfs <- as.data.frame(rbind(EC_s2, EE_s2, EI_s2, LC_s2, LE_s2, LI_s2))

# Plot the results

if(namePlot != "Dwell_Prob_state1"){
  dfs$medians <- dfs$medians * 100
  dfs$lower <- dfs$lower * 100
  dfs$upper <- dfs$upper * 100
}


p2c <- ggplot(dfs, aes(x = tod, y = medians, color = Infection)) +
  geom_line(linewidth = 1.25) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Infection), alpha = 0.3, colour =NA) +
  scale_fill_manual(values = c("black","#78baed", "#D55E00")) +
  scale_color_manual(values = c("black","#0072B2", "#D55E00")) + # use theme bw but remove grids 
    
    
    # scale x axis by dividing by 60
    scale_x_continuous(breaks = seq(0, 1440, 120), labels = seq(0, 24, 2)) + 
    facet_grid(Time ~ .) 
p2c <- p2c + ylab(Y1axis_name) + xlab("Time of day (Hours)") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
            ylim(0,max(dfs$upper)*1.1)
p2c
# ggsave("plots/EE&IE_state2c.pdf", p2c, width = 15, height = 15, units = "cm")



#---


EL_state2 <- log(ExposedLate[, , 2] / ControlLate[, , 2])   # as a matrix
IL_state2 <- log(InfectedLate[, , 2] / ControlLate[, , 2])   # as a matrix


EL_state2 <- data.frame(medians = apply(EL_state2, 2, median),
            lower = t(apply(EL_state2, 2, hdi))[,1],
            upper = t(apply(EL_state2, 2, hdi))[,2]
            )

EL_state2$Infection <- "Exposed"
EL_state2$Time <- "Late"
EL_state2$tod <- 1:1440

IL_state2 <- data.frame(medians = apply(IL_state2, 2, median),
            lower = t(apply(IL_state2, 2, hdi))[,1],
            upper = t(apply(IL_state2, 2, hdi))[,2]
            )


IL_state2$Infection <- "Infected"
IL_state2$Time <- "Late"
IL_state2$tod <- 1:1440

df <- as.data.frame(rbind(EE_state2, IE_state2, EL_state2, IL_state2))
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

p2d <- p2d + ylab(Y2axis_name) + xlab("Time of day (Hours)") + 
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"))

p2d

# ggsave("plots/EE&IE_state2.pdf",p2d, width = 15, height = 15, units = "cm")

# use ggarange to put them together

library(ggpubr)

theme_set(theme_bw())
figure <- ggarrange(p2c, p2d,
                    labels = c("A", "B"), common.legend = TRUE, legend = "bottom",
                    ncol = 2, nrow = 1)
figure
namePlot = "state2"
plotName = paste("plots/", namePlot, ".pdf", sep = "")

ggsave(plotName,figure, width = 20, height = 15, units = "cm")
