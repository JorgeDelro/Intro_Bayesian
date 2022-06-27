# Bayesian Data Analysis Using brms: A Tutorial for Sport Scientists

# libraries
library(tidyverse)
library(tidybayes)
library(brms)
library(bayestestR)
library(bayesplot)
library(loo)
library(ggpubr)
library(ggdist)

# Read Mergensen et al. database
db <- as_tibble(read.csv2(file = "db.csv",
                          header = T,
                          dec = "."))

# Hemoglobin database
dbHb <- db %>%
  select(ID, Group, ChangeWtr, HbmassPost, HbmassPre) %>%                  # Select variables for analysis
  mutate(HMabs = HbmassPost - HbmassPre,                                   # Compute Absolute values
         Group = factor(Group, levels = c("Placebo", "IHE", "LHTL"))) %>%  # Set Placebo level as reference
  na.omit()                                                                # Get complete cases

# Summary statistics
dbHb %>%
  group_by(Group) %>%                   # Results by Group
  summarise(mean_group = mean(HMabs),   # Summary by levels of Group variable
            sd_group = sd(HMabs))


# Pairwise differences between levels of Group variable
dbHb %>%
  summarize(mean_ChangeWtr = mean(ChangeWtr),
            sd_ChangeWtr = sd(ChangeWtr),
            mean_Placebo_IHE = mean(HMabs[Group == "Placebo"] - HMabs[Group == "IHE"]),
            sd_Placebo_IHE = sd(HMabs[Group == "Placebo"] - HMabs[Group == "IHE"]),
            mean_Placebo_LHTL = mean(HMabs[Group == "Placebo"] - HMabs[Group == "LHTL"]),
            sd_Placebo_LHTL = sd(HMabs[Group == "Placebo"] - HMabs[Group == "LHTL"]),
            mean_IHE_LHTL = mean(HMabs[Group == "IHE"] - HMabs[Group == "LHTL"]),
            sd_IHE_LHTL = sd(HMabs[Group == "IHE"] - HMabs[Group == "LHTL"]))

# Boxplot
dbHb %>% ggplot(aes(x = Group, y = HMabs, fill = Group)) +
  geom_boxplot() +
  labs(x = "Group", y = "Hemoglobin mass (g)") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14))

ggsave("Figure1.png", width = 7, height = 7, dpi= 600)


# Informative prior
# b_GroupLHTL
# Gore et al. (2013) - LHTL - increase for 240h of exposure 2.6% (2.3% - 2.9%)
dbHb %>%
  filter(Group == "LHTL") %>%
  mutate(HM_mean_change = (HbmassPre*2.6)/100,
         HM_low_CI = (HbmassPre*2.3)/100,
         HM_high_CI = (HbmassPre*2.9)/100) %>%
  summarize(Mean_exp_change = mean(HM_mean_change),
            LowCI_exp_change = mean(HM_low_CI),
            HighCI_exp_change = mean(HM_high_CI))


# Get priors
get_prior(formula = HMabs ~ 0 + Intercept + ChangeWtr + Group,
          data = dbHb,
          family = gaussian(link = "identity"))

# Set priors for bmod1 (Gaussian model)
bmod1Priors <- c(prior(normal(0, 2), class = "b", coef = "ChangeWtr"),    # prior for ChangeWtr effect
                 prior(normal(0, 2), class = "b", coef = "GroupIHE"),     # prior for IHE group effect
                 prior(normal(22.6, 1), class = "b", coef = "GroupLHTL"), # prior for LHTL group effect
                 prior(normal(0, 2), class = "b", coef = "Intercept"),  # prior for placebo level
                 prior(student_t(3, 0, 15), class = "sigma"))             # prior for LHTL effect


# Plot prior distributions
# Intercept
Intercept_prior <- ggplot(data.frame(x = c(-10, 10)), aes(x = x)) +
  stat_function(fun = dnorm,args = list(mean = 0, sd = 2), size = 1.5) +
  stat_function(fun = dnorm,geom="area", fill= "skyblue", args = list(mean = 0, sd = 2)) +
  labs(x = "Values", y = "Density") +
  ggtitle(expression(alpha~"~ Normal(0, 2)")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 12,
                                  hjust = 0.5))


# ChangeWtr
ChangeWtr_prior <- ggplot(data.frame(x = c(-10, 10)), aes(x = x)) +
  stat_function(fun = dnorm,args = list(mean = 0, sd = 2), size = 1.5) +
  stat_function(fun = dnorm,geom="area", fill= "skyblue", args = list(mean = 0, sd = 2)) +
  labs(x = "Values", y = "Density") +
  ggtitle(expression(beta[1]~"~ Normal(0, 2)")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 12,
                                  hjust = 0.5))


# IHE
IHE_prior <- ggplot(data.frame(x = c(-10, 10)), aes(x = x)) +
  stat_function(fun = dnorm,args = list(mean = 0, sd = 2), size = 1.5) +
  stat_function(fun = dnorm,geom="area", fill= "skyblue", args = list(mean = 0, sd = 2)) +
  labs(x = "Values", y = "Density") +
  ggtitle(expression(beta[IHE]~"~ Normal(0, 2)")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 12,
                                  hjust = 0.5))

# LHTL
LHTL_prior <- ggplot(data.frame(x = c(15, 30)), aes(x = x)) +
  stat_function(fun = dnorm,args = list(mean = 22.6, sd = 1), size = 1.5) +
  stat_function(fun = dnorm,geom="area", fill= "skyblue", args = list(mean = 22.6, sd = 1)) +
  labs(x = "Values", y = "Density") +
  ggtitle(expression(beta[LHTL]~"~ Normal(22.6, 1)")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 12,
                                  hjust = 0.5))

# sigma - residual variance
sigma_prior <- ggplot(data.frame(x = c(-100, 100)), aes(x = x)) +
  stat_function(fun = dstudent_t,args = list(df = 3, mu = 0, sigma = 15), size = 1.5) +
  stat_function(fun = dstudent_t,geom="area", fill= "skyblue", args = list(df = 3, mu = 0, sigma = 15)) +
  labs(x = "Values", y = "Density") +
  ggtitle(expression(sigma~"~ StudentT(0, 15, 3)")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 12,
                                  hjust = 0.5))

ggarrange(Intercept_prior,
          ChangeWtr_prior,
          IHE_prior,
          LHTL_prior,
          sigma_prior,
          nrow = 3,
          ncol = 2)

ggsave("Figure2.png", width = 7, height = 7, dpi= 600)

### Prior Predictive Checking
bmod1PPd <- brm(formula =  HMabs ~ 0 + Intercept + ChangeWtr + Group,
                data = dbHb,
                family = gaussian(link = "identity"),
                warmup = 1000,
                iter = 2000,
                chains = 4,
                seed = 1234,
                prior = bmod1Priors,                     # Prior distribution for population effects
                sample_prior = c("only"))                # sample only from priors

# Results bmod1Prior
summary(bmod1PPd)

# Check priors after fit the model
prior_summary(bmod1PPd)

# get parameters names
parnames(bmod1PPd)

# Set brightblue scheme
color_scheme_set("brightblue")

# Prior Parameter Distribution
# bmod1PPd %>%
#  plot(combo = c("dens", "trace"), widths = c(1, 1.5),
#       theme = theme_bw(base_size = 10))

# Prior predictive distribution
bmod1PPd %>%
  posterior_predict(draws = 50) %>%
  ppc_dens_overlay(y = dbHb$HMabs) +
  xlim(-750, 750)
# Prior predictive distribution by group condition
#bmod1Prior %>%
#  posterior_predict(draws = 50) %>%
#  ppc_stat_grouped(y = dbHb$HMabs,
#                   group = dbHb$Group,
#                   stat = "mean")


ggsave("Figure3.png", width = 5, height = 5, dpi= 600)

# Mergensen et al. (2015) model with weakly and informative priors on population effects
bmod1 <- brm(formula =  HMabs ~ 0 + Intercept + ChangeWtr + Group,
             data = dbHb,
             family = gaussian(link = "identity"),
             warmup = 1000,
             iter = 2000,
             chains = 4,
             seed = 1234,
             prior = bmod1Priors,                     # Prior distribution for population effects
             sample_prior = c("yes"))

summary(bmod1)

# Posterior parameter distribution
png('Figure4.png', units="in", width=6, height=5, res=600)
bmod1 %>%
  plot(combo = c("hist", "trace"), widths = c(1, 1.5),
       theme = theme_bw(base_size = 10))
dev.off()

# Posterior Predictive Distribution
ppc_dens_overlay(y = dbHb$HMabs, yrep = posterior_predict(bmod1, draws = 50)) +
  xlim(-150, 150)
ggsave("Figure5.png", width = 5, height = 5, dpi= 600)

## Model selection
# Change likelihood to Student-t
bmod2 <- brm(formula = HMabs ~ 0 + Intercept + ChangeWtr + Group,
             data = dbHb,
             family = student(link = "identity"),
             warmup = 1000,
             iter = 2000,
             chains = 4,
             seed = 1234,
             prior = bmod1Pdist,                     # Prior distribution for population effects
             sample_prior = "yes")

summary(bmod2)

# PSIS-LOO
loo1 <- loo(bmod1, save_psis = TRUE)
# print(loo1)
# plot(loo1, label_points = T)
loo2 <- loo(bmod2, save_psis = TRUE)
loo_compare(loo1, loo2)

# Post-hoc Contrasts
table_1 <- hypothesis(bmod1, c("Intercept = Intercept + GroupIHE",  # Placebo VS IHE
                    "Intercept = Intercept + GroupLHTL",             # Placebo VS LHTL
                    "Intercept + GroupIHE = Intercept + GroupLHTL")) # Placebo VS LHTL


######
###### Supplemental File Analysis
######

## Region of practical equivalence - bayestestR
equivalence_test(bmod1, range = c(-4, 4), ci = 0.95)

# ROPE
bmod1 %>%
  gather_draws(b_Intercept, b_GroupIHE, b_GroupLHTL) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < 4.0))) +
  stat_halfeye() +
  geom_vline(xintercept = c(-4.0, 4.0), linetype = "dashed") +
  scale_fill_manual(values = c("gray80", "skyblue")) +
  labs(x = "Hemoglobin mass (g)", y = "Condition") +
  scale_y_discrete(breaks = c("b_Intercept", "b_GroupIHE", "b_GroupLHTL"),
                   labels = c("Placebo", "IHE", "LHTL")) +
  annotate(geom = "text", x = 12, y = 3.5, label = 'bold("ROPE = (-4.0, 4.0)")',
           size = 5, parse = T) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=10),
    axis.title.x = element_text(size=12),
    axis.text.y = element_text(size=10),
    axis.title.y = element_text(size=12),
   legend.position = "none"
  )

# Save FigureS1
ggsave("FigureS1.png", width = 7, height = 5, dpi= 600)
