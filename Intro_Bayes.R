# An Introduction to Bayesian Data Analysis for Sport Scientists

library(tidyverse)
library(tidybayes)
library(brms)
# library(introbayes)
library(bayestestR)
library(effsize)
library(emmeans)
library(bayesplot)
library(loo)
library(posterior)
library(ggpubr)
library(ggdist)

# Set pseudo-random number generation
set.seed(7)

# ----------------------------------------------------------------------------
#  2.2.1. Metropolis Algorithm -----------------------------------------------
# ----------------------------------------------------------------------------



# Initialize the algorithm
theta_init <- 0
# Number of iterations
n_iterations <- 10^4
# Vector to store the results
theta <- rep(0,n_iterations)
# First value of the vector is the init
theta[1] <- theta_init

# Run the algorithm
for(i in 2:n_iterations){
  # Store the result at every iteration
  theta_current <- theta[i-1]
  # Propose a move
  theta_proposed <- theta_current + rnorm(1, mean = 0, sd = 1)
  # Probability of accept the move
  p_accept <- min(1, dnorm(x = theta_proposed,
                           mean = 10,
                           sd = 1) / dnorm(x = theta_current,
                                           mean = 10,
                                           sd = 1))

  # Generates a random value form Uniform(0,1)
  accept_value <- runif(1)
  # Accept or reject the proposed move
  if(accept_value < p_accept){
    theta[i] <- theta_proposed
  } else {
    theta[i] <- theta_current
  }
}

# Get the chain
Markov_chain <- data.frame(position = 1:10^4,
                           theta = theta)
Markov_chain_burn <- Markov_chain[1:100,]
Markov_chain_optim <- Markov_chain[9901:10^4,]

# Plot it
MC_1 <- Markov_chain_burn %>%
  ggplot(aes(x = position, y = theta)) +
  ylim(0, 15) +
  labs(y = expression(theta), x = "Number of chain steps (init)") +
  geom_line(size = 1.2, colour="skyblue") +
  theme_bw() +
  theme(panel.grid = element_blank())

MC_2 <-Markov_chain_optim %>%
  ggplot(aes(x = position, y = theta)) +
  ylim(0, 15) +
  labs(y = expression(theta), x = "Number of chain steps (final)") +
  geom_line(size = 1.2, colour="skyblue") +
  theme_bw() +
  theme(panel.grid = element_blank())

MC_3 <-Markov_chain %>%
  ggplot(aes(x = theta)) +
  geom_histogram(color = "white", fill="skyblue", binwidth = 0.3) +
  xlim(0, 15) +
  labs(x = expression(theta), y = expression("Number of accept move of " ~ theta)) +
  theme_bw() +
  theme(#axis.title.x = element_blank(),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.grid = element_blank()) +
  coord_flip()


ggarrange(MC_1, MC_2, MC_3,
          ncol = 3,
          nrow = 1)

# ----------------------------------------------------------------------------
#  2.2.2. Hamiltonian Monte Carlo Algorithm (Gelman, 2013)--------------------
# ----------------------------------------------------------------------------

# log posterior function
# X: the model matrix
# y: the target vector
# th: theta, the current parameter estimates
log_posterior <- function(X, y, theta) {


  beta <- th[-length(theta)]
  sigma <- th[length(theta)]
  sigma2 <- sigma^2
  mu <- X %*% beta

  # priors
  priorbvarinv <- diag(1/100, 4)
  prioralpha <- priorbeta = .001

  if (is.nan(sigma) | sigma<=0) {
    return(-Inf)
  }
  else {
    -.5*nrow(X)*log(sigma2) - (.5*(1/sigma2) * (crossprod(y-mu))) +
      -.5*ncol(X)*log(sigma2) - (.5*(1/sigma2) * (t(beta) %*% priorbvarinv %*% beta)) +
      -(prioralpha + 1)*log(sigma2) + log(sigma2) - priorbeta/sigma2
  }

} # end log posterior

# gradient
gradient_theta <- function(X, y, th) {
  d <- length(th)
  e <- .0001
  diffs <- numeric(d)

  for (k in 1:d) {
    th_hi <- th
    th_lo <- th
    th_hi[k] <- th[k] + e
    th_lo[k] <- th[k] - e
    diffs[k] <- (log_posterior(X, y, th_hi) - log_posterior(X, y, th_lo)) / (2 * e)
  }

  diffs
} # end gradient

# function for a single HMC iteration

hmc_iteration <- function(X, y, th, epsilon, L, M) {
  # Args
  # epsilon: the stepsize
  # L: the number of leapfrog steps
  # M: a diagonal mass matrix

  # initialization
  M_inv = 1/M
  d   = length(th)
  phi = rnorm(d, 0, sqrt(M))
  th_old = th

  log_p_old = log_posterior(X, y, th) - .5*sum(M_inv * phi^2)

  phi = phi + .5 * epsilon * gradient_theta(X, y, th)

  for (l in 1:L) {
    th  = th + epsilon*M_inv*phi
    phi = phi + ifelse(l == L, .5, 1) * epsilon * gradient_theta(X, y, th)
  }

  # here we get into standard MCMC stuff, jump or not based on a draw from a
  # proposal distribution
  phi = -phi
  log_p_star = log_posterior(X, y, th) - .5*sum(M_inv * phi^2)
  r = exp(log_p_star - log_p_old)

  if (is.nan(r)) r = 0

  p_jump = min(r, 1)

  if (runif(1) < p_jump) {
    th_new = th
  }
  else {
    th_new = th_old
  }

  # returns estimates and acceptance rate
  list(th = th_new, p_jump = p_jump)

}

# Main HMC function
hmc_run <- function(starts, iter, warmup, epsilon_0, L_0, M, X, y) {
  # # Args:
  # starts:  starting values
  # iter: total number of simulations for each chain (note chain is based on the dimension of starts)
  # warmup: determines which of the initial iterations will be ignored for inference purposes
  # epsilon0: the baseline stepsize
  # L0: the baseline number of leapfrog steps
  # M: is the mass vector
  chains = nrow(starts)
  d = ncol(starts)
  sims = array(NA,
               c(iter, chains, d),
               dimnames = list(NULL, NULL, colnames(starts)))
  p_jump = matrix(NA, iter, chains)

  for (j in 1:chains) {
    th = starts[j,]

    for (t in 1:iter) {
      epsilon = runif(1, 0, 2*epsilon_0)
      L    = ceiling(2*L_0*runif(1))

      temp = hmc_iteration(X, y, th, epsilon, L, M)

      p_jump[t,j] = temp$p_jump
      sims[t,j,]  = temp$th

      th = temp$th
    }
  }

  # acceptance rate
  acc = round(colMeans(p_jump[(warmup + 1):iter,]), 3)

  message('Avg acceptance probability for each chain: ',
          paste0(acc[1],', ',acc[2]), '\n')

  list(sims = sims, p_jump = p_jump)
}

# Simulate data


## Example
# Starting values and mcmc settings
parnames <- c(paste0('beta[', 1:4, ']'), 'sigma')
d <- length(parnames)

chains <- 2

theta_start <- t(replicate(chains, c(runif(d-1, -1, 1), 1)))
colnames(theta_start) <- parnames

nsim <- 1000
wu <- 500

stepsize <- .08
nLeap <- 10
vars <- rep(1, 5)
mass_vector <- 1 / vars

# Run the model
fit_hmc = hmc_run(
  starts    = theta_start,
  iter      = nsim,
  warmup    = wu,
  epsilon_0 = stepsize,
  L_0 = nLeap,
  M   = mass_vector,
  X   = X,
  y   = y
)

# ----------------------------------------------------------------------------
#  3. Applied Bayesian data analysis example   -------------------------------
# ----------------------------------------------------------------------------


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
#dbHb %>%
#  filter(Group == "LHTL") %>%
#  mutate(HM_mean_change = (HbmassPre*2.6)/100,
#         HM_low_CI = (HbmassPre*2.3)/100,
#         HM_high_CI = (HbmassPre*2.9)/100) %>%
#  summarize(Mean_exp_change = mean(HM_mean_change),
#            LowCI_exp_change = mean(HM_low_CI),
#            HighCI_exp_change = mean(HM_high_CI))


# Get priors
get_prior(formula = HMabs ~ 0 + Intercept + ChangeWtr + Group,
          data = dbHb,
          family = gaussian(link = "identity"))

# Set priors for bmod1 (Gaussian model)
bmod1Priors <- c(prior(normal(0, 2), class = "b", coef = "ChangeWtr"),    # prior for ChangeWtr effect
                 prior(normal(0, 2), class = "b", coef = "GroupIHE"),     # prior for IHE group effect
                 prior(normal(2.6, 0.5), class = "b", coef = "GroupLHTL"), # prior for LHTL group effect
                 prior(normal(0, 2), class = "b", coef = "Intercept"),  # prior for placebo level
                 prior(student_t(3, 0, 15), class = "sigma"))             # prior for LHTL effect


# Plot prior distributions ----------------------------------------------------
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
LHTL_prior <- ggplot(data.frame(x = c(-10, 10)), aes(x = x)) +
  stat_function(fun = dnorm,args = list(mean = 2.6, sd = 0.5), size = 1.5) +
  stat_function(fun = dnorm,geom="area", fill= "skyblue", args = list(mean = 2.6, sd = 0.5)) +
  labs(x = "Values", y = "Density") +
  ggtitle(expression(beta[LHTL]~"~ Normal(2.6, 0.5)")) +
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
  ggtitle(expression(sigma~"~ Half-StudentT(0, 15, 3)")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlim(c(0,NA)) +
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

### Prior Predictive Checking ----------------------------------------------------
bmod1_prior <- brm(formula = HMabs ~ 0 + Intercept + ChangeWtr + Group,
                data = dbHb,
                family = gaussian(link = "identity"),
                warmup = 1000,
                iter = 2000,
                chains = 4,
                seed = 1234,
                prior = bmod1Priors,          # Prior distribution for population effects
                sample_prior = c("only"))     # sample only from priors

# Results bmod1Prior
# summary(bmod1PPd)

# Check priors after fit the model
prior_summary(bmod1PPd)

# get parameters names
parnames(bmod1PPd)

# Set brightblue scheme
color_scheme_set("brightblue")

# Prior predictive distribution
bmod1_prior %>%
  posterior_predict(draws = 50) %>%
  ppc_dens_overlay(y = dbHb$HMabs) +
  xlim(-500, 500)

ggsave("Figure3.png", width = 5, height = 5, dpi= 600)

# Mergensen et al. (2015) model with weakly and informative priors on population effects
bmod1 <- brm(formula =  HMabs ~ 0 + Intercept + ChangeWtr + Group,
             data = dbHb,
             family = gaussian(link = "identity"),
             warmup = 1000,
             iter = 2000,
             chains = 4,
             seed = 1234,
             prior = bmod1Priors,   # Prior distribution for population effects
             sample_prior = c("yes"))


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
             prior = bmod1Priors,
             sample_prior = "yes")

# summary(bmod2)

# PSIS-LOO ------------------------------------------------------------------
loo1 <- loo(bmod1, save_psis = TRUE)
# print(loo1)
# plot(loo1, label_points = T)
loo2 <- loo(bmod2, save_psis = TRUE)
loo_compare(loo1, loo2)

# ROPE analysis - bayestestR --------------------------------------------------------
describe_posterior(bmod1, rope_range = c(-0.5, 0.5), ci_method = "HDI", ci = 1)
equivalence_test(bmod1, range = c(-0.5, 0.5), ci = 1)

# plot - ROPE - tidybayes
bmod1 %>%
  gather_draws(b_Intercept, b_GroupIHE, b_GroupLHTL) %>%
  ggplot(aes(y = .variable, x = .value, fill = stat(abs(x) < 0.5))) +
  stat_halfeye() +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  scale_fill_manual(values = c("skyblue", "gray80")) +
  labs(x = "Hemoglobin mass (g)", y = "Condition") +
  scale_y_discrete(breaks = c("b_Intercept", "b_GroupIHE", "b_GroupLHTL"),
                   labels = c("Placebo", "IHE", "LHTL")) +
  annotate(geom = "text", x = 4, y = 3.5, label = 'bold("ROPE = (-0.5, 0.5)")',
           size = 4.5, parse = T) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=10),
    axis.title.x = element_text(size=12),
    axis.text.y = element_text(size=10),
    axis.title.y = element_text(size=12),
    legend.position = "none"
  )

ggsave("Figure6.png", width = 5, height = 4.5, dpi= 600)


# Create table 2
# tbl_a <- summary(bmod1)
# tbl_b <- describe_posterior(bmod1, rope_range = c(-0.5, 0.5))
# tbl_c <- equivalence_test(bmod1, range = c(-0.5, 0.5))
# tbl2 <- data.frame(tbl_a$fixed, tbl_b$pd, tbl_b$ROPE_Percentage, tbl_c$ROPE_Equivalence)
# rownames(tbl2) <- c("alpha", "beta_1", "beta_2-IHE", "beta_2-LHTL")
# colnames(tbl2) <- c("mean","SE", "lower bound", "upper bound", "Rhat", "Bulk ESS", "Tail ESS", "pd",
#                    "% in ROPE", "H0")

# Post-hoc Contrasts - brms --------------------------------------------------------
hypothesis(bmod1, c("Intercept + GroupIHE = Intercept ",  # IHE VS Placebo
                    "Intercept + GroupLHTL = Intercept",  # LHTL VS Placebo
                    "Intercept + GroupIHE = Intercept + GroupLHTL")) # IHE VS LHTL

# Create table 3
tbl_a <- hypothesis(bmod1, c("Intercept + GroupIHE = Intercept ",  # IHE VS Placebo
                            "Intercept + GroupLHTL = Intercept",  # LHTL VS Placebo
                            "Intercept + GroupIHE = Intercept + GroupLHTL"), # IHE VS LHTL
                   alpha = 0.00) # Full posterior distribution

 bmod1_draws <- as_draws_df(bmod1)
bmod1_df <- tibble(Placebo = bmod1_draws$b_Intercept,
                   IHE = bmod1_draws$b_Intercept + bmod1_draws$b_GroupIHE,
                   LHTL = bmod1_draws$b_Intercept + bmod1_draws$b_GroupLHTL)
tbl_b <- data.frame(c(
  cohen.d(bmod1_df$IHE, bmod1_df$Placebo)$estimate,
  cohen.d(bmod1_df$LHTL, bmod1_df$Placebo)$estimate,
  cohen.d(bmod1_df$IHE, bmod1_df$LHTL)$estimate
))
colnames(tbl_b) <- "Effect size"

tbl3 <- data.frame(Hypothesis = c("1", "2", "3"),
                   tbl_a$hypothesis$Estimate,
                   tbl_a$hypothesis$Est.Error,
                   tbl_a$hypothesis$CI.Lower,
                   tbl_a$hypothesis$CI.Upper,
                   tbl_a$hypothesis$Evid.Ratio,
                   Evidence = c("Anecdotical", "Extreme", "Anecdotical"),
                   tbl_b$`Effect size`)

colnames(tbl3) <- c("hypothesis", "Mean diff", "EE",
                    "Lower bound", "Upper bound", "BF", "Evidence",
                    "Effect size")

tbl3[, 6] <- c("1.05", ">100", "0.91")


# Post-hoc Contrasts - bayestestR --------------------------------------------------------
# group_diff_prior <- emmeans(bmod1PPd, pairwise ~ Group)
# group_diff <- emmeans(bmod1, pairwise ~ Group)
# bayesfactor_parameters(group_diff, group_diff_prior)


# Sensitivity analysis -----------------------------------------------------
# Model with non-informative priors
# Set non-informative priors
bmod3Priors <- c(prior(normal(0, 10^5), class = "b", coef = "ChangeWtr"),
                 prior(normal(0, 10^5), class = "b", coef = "GroupIHE"),
                 prior(normal(0, 10^5), class = "b", coef = "GroupLHTL"),
                 prior(normal(0, 10^5), class = "b", coef = "Intercept"),
                 prior(normal(0, 10^5), class = "sigma"))


bmod3 <- brm(formula = HMabs ~ 0 + Intercept + ChangeWtr + Group,
             data = dbHb,
             family = gaussian(link = "identity"),
             warmup = 1000,
             iter = 2000,
             chains = 4,
             seed = 1234,
             prior = bmod3Priors,
             sample_prior = "yes")


# Sensitivity analysis ------------------------------------------------------
sensitivity_results <- sensitivity_analysis(bmodels = list(original_model = bmod1,
                                                  alternative_prior = bmod3),
                                   params = c("b_Intercept", "b_ChangeWtr",
                                              "b_GroupIHE", "b_GroupLHTL",
                                              "sigma"))
# Table
sensitivity_results$results_tables

# plot
ggarrange(
  # b_Intercept
  sensitivity_results$results_plots[[1]] +
    ggtitle(expression(alpha)) +
    scale_fill_manual(values=c("#FFCCCB", "skyblue")),
  # b_ChangeWtr
  sensitivity_results$results_plots[[2]] +
    ggtitle(expression(beta[1])) +
    scale_fill_manual(values=c("#FFCCCB", "skyblue")),
  # b_GroupIHE
  sensitivity_results$results_plots[[3]] +
    ggtitle(expression(beta[2-IHE])) +
    scale_fill_manual(values=c("#FFCCCB", "skyblue")),
  # b_GroupLHTL
  sensitivity_results$results_plots[[4]] +
    ggtitle(expression(beta[2-LHTL])) +
    scale_fill_manual(values=c("#FFCCCB", "skyblue")),
  # sigma
  sensitivity_results$results_plots[[5]] +
    ggtitle(expression(sigma)) +
    scale_fill_manual(values=c("#FFCCCB", "skyblue")),
  nrow = 3,
  ncol = 2)

ggsave("Figure6.png", width = 7, height = 7, dpi= 600)

