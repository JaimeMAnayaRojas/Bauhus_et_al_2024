
# use the value of a variable for a name to save the plot


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
# ggsave("plots/EE&IE_state1c.pdf", p2c, width = 15, height = 15, units = "cm")



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

p2d <- p2d + ylab(Y2axis_name) + xlab("Time of day (Hours)") + 
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"))

p2d

# ggsave("plots/EE&IE_state1.pdf",p2d, width = 15, height = 15, units = "cm")

# use ggarange to put them together

library(ggpubr)

theme_set(theme_bw())
figure <- ggarrange(p2c, p2d,
                    labels = c("A", "B"), common.legend = TRUE, legend = "bottom",
                    ncol = 2, nrow = 1)
figure

plotName = paste("plots/", namePlot, ".pdf", sep = "")
ggsave(plotName,figure, width = 20, height = 15, units = "cm")
