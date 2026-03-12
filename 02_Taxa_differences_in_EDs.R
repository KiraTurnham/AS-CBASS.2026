#Code for merging processed site level ED data

#load packages
library(tidyverse)
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
library(emmeans)
library(purrr)


rm(list = ls())
setwd("C:/github/CBASS/AS-CBASS.2026")
#####Load and Prep Data-----------------------------

EDdf <- read_csv("Data/EDsdf_all_ramps.csv") %>% 
  filter(Dataset =="All")%>%
  mutate(across(where(is.character), as.factor))

#long format 
ED_long <- EDdf %>%
  pivot_longer(
    cols = c(ED5, ED50, ED95),
    names_to = "ED_level",
    values_to = "ED_value"
  ) %>%
  mutate(Species = factor(Species, levels = c("AABR", "AHYA", "AGLO", "ICRA")))


#stats: does ED-value differ among species; mixed model for each ED level-------------

#subset by ED
ED5_df  <- filter(ED_long, ED_level == "ED5")
ED50_df <- filter(ED_long, ED_level == "ED50")
ED95_df <- filter(ED_long, ED_level == "ED95")


#evaluate parametric assumptions
#Shapiro-Wilk Normality
ED5_df %>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value)$p.value)

ED50_df %>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value)$p.value)

ED95_df %>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value)$p.value)
#not very normal --> log transform didnʻt help, go non parametric

#Homogeneity of variance
leveneTest(ED_value ~ Species, data = ED5_df)
leveneTest(ED_value ~ Species, data = ED50_df)
leveneTest(ED_value ~ Species, data = ED95_df)
#all fine here

#log-transformed ED95
ED_long <- ED_long %>%
  mutate(ED_value_log = ifelse(ED_level == "ED95", log(ED_value), ED_value))

models <- ED_long %>%
  split(.$ED_level) %>%
  map(~ lmer(ED_value_log ~ Species + (1|Site), data = .))

lapply(models, anova)
lapply(models, function(m) emmeans(m, pairwise ~ Species))

# extract ANOVA tables and add ED_level column
anova_table <- lapply(names(models), function(ed) {
  m <- models[[ed]]
  an <- anova(m)
  an_df <- as.data.frame(an) %>%
    rownames_to_column("term") %>%
    mutate(ED_level = ed)
  an_df
}) %>%
  bind_rows() %>%
  select(-term)

anova_table

# extract emmeans pairwise comparisons
emmeans_table <- lapply(names(models), function(ed) {
  m <- models[[ed]]
  em <- emmeans(m, pairwise ~ Species)
  tidy(em$contrasts) %>%
    mutate(ED_level = ed)
}) %>%
  bind_rows()

# view
emmeans_table

#check residual diagnostics for ED95
m_ED95 <- lmer(ED_value ~ Species + (1|Site), data = ED95_df)
# residual diagnostics
par(mfrow = c(1,2))
plot(resid(m_ED95) ~ fitted(m_ED95))
qqnorm(resid(m_ED95))
qqline(resid(m_ED95))

#alternatively, fit a combined model:
lmm_combined <- lmer(ED_value ~ Species * ED_level + (1 | Site), data = ED_long)
anova(lmm_combined)

emmeans_species <- emmeans(lmm_combined, pairwise ~ Species | ED_level)
emmeans_species

VarCorr(lmm_combined)

####plotting----------------------------------


#summary stats; mean + SE
ED_summary <- ED_long %>%
  group_by(Species, ED_level) %>%
  summarise(
    Mean = mean(ED_value),
    SE = sd(ED_value) / sqrt(n()),
    .groups = "drop")

site_summary <- ED_long %>%
  group_by(Site, Species, ED_level) %>%
  summarise(mean_ED = mean(ED_value), .groups="drop")


#boxplot of ED5, 50, 95-----------
g1 <- ggplot(ED_long, aes(x = Species, y = ED_value, fill = Species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(aes(color = Species), width = 0.15, alpha = 0.5, size = 1, show.legend = FALSE) +
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


