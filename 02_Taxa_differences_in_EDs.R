#Code for merging processed site level ED data

#load packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(rstudioapi)
library(RColorBrewer)
library(CBASSED50) #Voolstra R package
library(car) #Leveneʻs test
library(broom)
library(FSA)
library(lme4)
library(broom.mixed)
library(lmerTest)


rm(list = ls())
setwd("C:/github/CBASS/AS-CBASS.2026")
#####Load and Prep Data-----------------------------

S1 <- read_csv("Data/RAW_ESA-001_20260214_PAM_EDsdf.csv")
S3 <- read_csv("Data/RAW_ESA-003_20260219_PAM_EDsdf.csv")
S4 <- read_csv("Data/RAW_ESA-004_20260220_PAM_EDsdf.csv")
S5 <- read_csv("Data/RAW_ESA-005_20260223_PAM_EDsdf.csv")
S7 <- read_csv("Data/RAW_ESA-007_20260213_PAM_EDsdf.csv")
S9 <- read_csv("Data/RAW_ESA-009_20260216_PAM_EDsdf.csv")
S10 <- read_csv("Data/RAW_ESA-010_20260217_PAM_EDsdf.csv")
S11 <- read_csv("Data/RAW_ESA-011_20260210_PAM_EDsdf.csv")

EDdf <- rbind(S1,S3, S4, S5, S7, S9, S10, S11)
EDdf <- EDdf %>% mutate(across(where(is.character), as.factor))
colnames(EDdf)

write.csv(EDdf,
          "C:/github/CBASS/AS-CBASS.2026/Data/EDsdf-all.csv",
          row.names = FALSE)


#long format for facet plotting
ED_long <- EDdf %>%
  pivot_longer(
    cols = c(ED5, ED50, ED95),
    names_to = "ED_level",
    values_to = "ED_value"
  )

#summary stats; mean + SE
ED_summary <- ED_long %>%
  group_by(Species, ED_level) %>%
  summarise(
    Mean = mean(ED_value),
    SE = sd(ED_value) / sqrt(n()),
    .groups = "drop")


#set species order for plotting
ED_long <- ED_long %>% mutate(Species = factor(Species, levels = c("AABR", "AHYA", "AGLO", "ICRA")))

#boxplot of ED5, 50, 95-----------
g1 <- ggplot(ED_long, aes(x = Species, y = ED_value, fill = Species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(aes(color = Species),
              width = 0.15,
              alpha = 0.5,
              size = 1,
              show.legend = FALSE) +
  facet_wrap(~ ED_level, scales = "free_y") +
  theme_classic(base_size = 14) +
  labs(
    x = "Species",
    y = expression(paste("Temperature (", degree, "C)")),
    title = "CBASS-Derived Thermal Thresholds",
    subtitle = "American Samoa"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  theme(legend.position = "none")

g1

#Add Mean ± SE to Plot; drop boxplot

g2 <- ggplot(ED_long, aes(x = Species, y = ED_value)) +
  geom_jitter(aes(color = Species),width = 0.15, alpha = 0.4, size = 1, show.legend = FALSE) +
  geom_point(data = ED_summary,
             aes(x = Species, y = Mean, color = Species),
             size = 4,
             inherit.aes = FALSE) +
  
  geom_errorbar(data = ED_summary,
                aes(x = Species,
                    ymin = Mean - SE,
                    ymax = Mean + SE,
                    color = Species),
                width = 0.2,
                size = 0.8,
                inherit.aes = FALSE) +
  
  facet_wrap(~ ED_level, scales = "free_y") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Species",
    y = expression(paste("Temperature (", degree, "C)")),
    title = "CBASS-Derived Thermal Thresholds",
    subtitle = "American Samoa"
  )

g2


ggsave("C:/github/CBASS/AS-CBASS.2026/Plots/All_sites_boxplot.pdf", g1,  width = 16, height = 9, device = "pdf")
ggsave("C:/github/CBASS/AS-CBASS.2026/Plots/All_sites_points&means.pdf", g2,  width = 16, height = 9, device = "pdf")

#stats: does ED-value differ among species-------------

ED5_df  <- filter(ED_long, ED_level == "ED5")
ED50_df <- filter(ED_long, ED_level == "ED50")
ED95_df <- filter(ED_long, ED_level == "ED95")

#Shapiro-Wilk Normality
ED5_df %>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value)$p.value)

ED50_df %>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value)$p.value)

ED95_df %>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value_log)$p.value)
#not very normal --> log

ED95_df <- ED95_df %>%
  mutate(ED_value_log = log(ED_value))
ED95_df %>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value_log)$p.value)
#not very normal --> log transform didnʻt help, go non parametric

#Homogeneity of variance
leveneTest(ED_value ~ Species, data = ED5_df)
leveneTest(ED_value ~ Species, data = ED50_df)
leveneTest(ED_value ~ Species, data = ED95_df)
#all fine here


#quick table of ANOVA results for ED5, ED50, ED95
ED_long %>%
  group_by(ED_level) %>%
  do(tidy(aov(ED_value ~ Species, data = .)))

#now with site added as a random factor
ED_long %>%
  group_by(ED_level) %>%
  do(anova(lmer(ED_value ~ Species + (1 | Site), data = .)))

#non-parametric follow-up for ED95
kruskal.test(ED_value ~ Species, data = ED95_df)

