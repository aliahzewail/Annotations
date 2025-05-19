library(tidyverse)
library(ggpubr)
library(rstatix)
library(papaja)
library(patchwork)
library(brms)
library(lmerTest)

Annotation_data <- read_csv("Annotation/Data.csv")

names(Annotation_data)
class(Annotation_data$Group)
Annotation_data$Group <- as.factor(Annotation_data$Group)
Annotation_data$Time <- as.factor(Annotation_data$Time)
Annotation_data$ID <- as.factor(Annotation_data$ID)

Annotation_data %>%
  group_by(Time, Group) %>%
  identify_outliers(RWA)


Annotation_data <- subset(Annotation_data,Annotation_data$ID != 15)


Annotation_data$MH[is.na(Annotation_data$MH)] <- mean(Annotation_data$MH, na.rm = TRUE)

Annotation_data$PTSD[is.na(Annotation_data$PTSD)] <- mean(Annotation_data$PTSD, na.rm = TRUE)

#Linear Model

lmm_model <- lmer(MH~ Time * Group+Age+Cohort+ (1 | ID), data = Annotation_data)
summary(lmm_model)
lmm_model <- lmer(PTSD ~ Time * Group+Age+Cohort+ (1 | ID), data = Annotation_data)
summary(lmm_model)

  

p1 <- ggplot(Annotation_data, aes(x = Time, y = MH, fill = Group)) +
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) +
  geom_point(aes(x = as.numeric(Time)-.15, y = MH, colour = Group), position = position_jitter(width = .05), size = 1, shape = 20) +
  geom_boxplot(aes(x = Time, y = MH, fill = Group), outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +  theme_apa() +ylab("MH") + xlab("")+
  theme(legend.position = "none")
p1<-p1 +
  stat_summary(
    aes(group = Group, colour = Group),
    fun    = median,
    geom   = "line",
    size   = 1
  )


p2 <- ggplot(Annotation_data, 
             aes(x = factor(Time, levels = c("Pretest", "Posttest")),
                 y = PTSD, 
                 fill = Group)) +
  geom_flat_violin(
    aes(fill = Group), 
    position = position_nudge(x = .1, y = 0), 
    adjust = 1.5, trim = FALSE, alpha = .5, colour = NA
  ) +
  geom_point(
    aes(x = as.numeric(factor(Time, levels = c("Posttest", "Pretest"))) - .15, 
        y = PTSD, 
        colour = Group),
    position = position_jitter(width = .05), 
    size = 1, shape = 20
  ) +
  geom_boxplot(
    aes(x = factor(Time, levels = c("Pretest", "Posttest")), 
        y = PTSD, 
        fill = Group),
    outlier.shape = NA, alpha = .5, width = .1, colour = "black"
  ) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_apa() +
  ylab("PTSD") +
  xlab("") +
  theme(legend.position = "right")

p2<-p2 +
  stat_summary(
    aes(group = Group, colour = Group),
    fun    = median,
    geom   = "line",
    size   = 1
  )

p1 + p2+
  plot_layout(ncol = 2, widths = c(1, 1))

#Bayesian Model

bayesian_model <- brm(
MH~ Time * Group + Age + Cohort +(1 | ID),  
  data = Annotation_data,
  family = gaussian(), 
  prior = set_prior("normal(0, 10)", class = "b"),  
  chains = 4, iter = 2000, warmup = 1000, cores = 4,  
  save_pars = save_pars(all = TRUE)
)
summary(bayesian_model)

# a null model (without interaction)
null_model <- brm(
MH~ Time + Group + Age + Cohort + (1 | ID),
  data = Annotation_data,
  family = gaussian(),
  prior = set_prior("normal(0, 10)", class = "b"),
  chains = 4, iter = 2000, warmup = 1000, cores = 4,
  save_pars = save_pars(all = TRUE)
)
summary(null_model)
library(bridgesampling)

# Compute bridge sampling estimates
bridge_full <- bridge_sampler(bayesian_model)
bridge_null <- bridge_sampler(null_model)

bf(bridge_full, bridge_null)

library(ggplot2)
library(tidybayes)

# posterior samples for the interaction term
post_samples <- posterior_samples(bayesian_model, pars = "b_TimePretest:GroupTreatment")

# viz 
ggplot(post_samples, aes(x = `b_TimePretest:GroupTreatment`)) +
  geom_density(fill = "blue", alpha = 0.3) +  
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +  
  labs(title = "Posterior dist of the interaction Effect",
       x = "treatment Effect (Time × Group)",
       y = "") +
  theme_minimal()

#diagnostics 
ggstatsplot::ggwithinstats(data=Annotation_data%>%dplyr::filter(Group == "Treatment"),
                           y =MH,
                           x = Time,type = "paramteric")



bayesian_model <- brm(
  PTSD~ Time * Group + Age + Cohort +(1 | ID),  
  data = Annotation_data,
  family = gaussian(), 
  prior = set_prior("normal(0, 10)", class = "b"),  
  chains = 4, iter = 2000, warmup = 1000, cores = 4,  
  save_pars = save_pars(all = TRUE)
)
summary(bayesian_model)

# a null model (without interaction)
null_model <- brm(
  PTSD~ Time + Group + Age + Cohort + (1 | ID),
  data = Annotation_data,
  family = gaussian(),
  prior = set_prior("normal(0, 10)", class = "b"),
  chains = 4, iter = 2000, warmup = 1000, cores = 4,
  save_pars = save_pars(all = TRUE)
)
summary(null_model)
library(bridgesampling)

# Compute bridge sampling estimates
bridge_full <- bridge_sampler(bayesian_model)
bridge_null <- bridge_sampler(null_model)

bf(bridge_full, bridge_null)

library(ggplot2)
library(tidybayes)

# posterior samples for the interaction term
post_samples <- posterior_samples(bayesian_model, pars = "b_TimePretest:GroupTreatment")

# viz 
ggplot(post_samples, aes(x = `b_TimePretest:GroupTreatment`)) +
  geom_density(fill = "blue", alpha = 0.3) +  
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +  
  labs(title = "Posterior dist of the interaction Effect",
       x = "treatment Effect (Time × Group)",
       y = "") +
  theme_minimal()

#diagnostics 
ggstatsplot::ggwithinstats(data=Annotation_data%>%dplyr::filter(Group == "Treatment"),
                           y =MH,
                           x = Time,type = "paramteric")

