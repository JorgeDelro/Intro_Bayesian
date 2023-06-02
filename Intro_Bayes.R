# An Introduction to Bayesian Data Analysis for Sport Scientists

library(tidyverse)
library(tidybayes)
library(brms)
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
  geom_line(linewidth = 1.2, colour="skyblue") +
  theme_bw() +
  theme(panel.grid = element_blank())

MC_2 <-Markov_chain_optim %>%
  ggplot(aes(x = position, y = theta)) +
  ylim(0, 15) +
  labs(y = expression(theta), x = "Number of chain steps (final)") +
  geom_line(linewidth = 1.2, colour="skyblue") +
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

# Figure 1
ggarrange(MC_1, MC_2, MC_3,
          ncol = 3,
          nrow = 1)

# ----------------------------------------------------------------------------
#  2.2.2. Hamiltonian Monte Carlo Algorithm (McElreath, 2020) --------------------
# ----------------------------------------------------------------------------

# U - neg-log-probability
U <- function(q, a = 0, b = 1, k = 0, d = 1){
  muy <- q[1]
  mux <- q[2]
  U <- sum(dnorm(y, muy, 1, log = TRUE)) + sum(dnorm(x, mux, 1, log = TRUE)) +
    dnorm(muy, a, b, log = TRUE) + dnorm(mux, k, d, log = TRUE)
  return(-U)
}

# U_gradient - partial derivatives of U with respect to vector q
U_gradient <- function(q, a = 0, b = 1, k = 0, d = 1){
  muy <- q[1]
  mux <- q[2]
  G1 <- sum(y - muy) + (a - muy)/b^2 #dU/dmuy
  G2 <- sum(x - mux) + (a - mux)/d^2 #dU/dmux
  return( c(-G1, -G2)) # negative bc energy is neg-log-prob
}

# Hamiltonian Monte Carlo function
HMC <- function(U, grad_U, epsilon, L, current_q){
  q <- current_q
  p <- rnorm(length(q), 0, 1) # random flick - p is momentum
  current_p <- p
  # make half step for momentum at the beginning
  p <- p - epsilon * grad_U(q) / 2
  # initialize bookkeeping - saves trajectory
  qtraj <- matrix(NA, nrow =  L+1, ncol = length(q))
  ptraj <- qtraj
  qtraj[1,] <- current_q
  ptraj[1,] <- p

  # alternate full steps for position and momentum
  for (i in 1:L) {
    q <- q + epsilon * p # full step for the position
    # make a full step for the momentum, except at the end of the trajectory
    if(i!=L){
      p <- p - epsilon * grad_U(q)
      ptraj[i+1,] <- p
    }
    qtraj[i+1,] <- q
  }

  # make a half step for momentum at the end
  p <- p - epsilon * grad_U(q) / 2
  ptraj[L+1,] <- p
  # negative momentum at the end of trajectory to make the proposal symmetric
  p <- -p
  # evaluate potential and kinetic energies at start and end of trajectory
  current_U <- U(current_q)
  current_K <- sum(current_p^2) / 2
  proposed_U <- U(q)
  proposed_K <- sum(p^2) / 2
  # accept of reject the state at the end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  accept <- 0
  if(runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) {
    new_q <- q # accept
    accept <- 1
  } else new_q <- current_q # reject

  return(list(q = new_q,
              traj = qtraj,
              ptraj = ptraj,
              accept = accept))

}

# test data
set.seed(7)
y <- rnorm(50)
x <- rnorm(50)
y <- as.numeric(scale(y))
x <- as.numeric(scale(x))

## Plot a multivariate normal distribution
set.seed(7)

# Target parameters for univariate normal distributions
rho <- 0.0
mu1 <- 0; s1 <- 1
mu2 <- 0; s2 <- 1

# Parameters for bivariate normal distribution
mu <- c(mu1,mu2) # Mean
sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2) # Covariance matrix

# Store in a dataset
data.grid <- expand.grid(s.1 = seq(-0.4, 0.4, length.out=200), s.2 = seq(-0.4, 0.4, length.out=200))

# Add density
q1_samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu, sigma=sigma))

HMC_plot <- ggplot() +
  geom_contour(data=q1_samp,aes(x=s.1,y=s.2,z=prob),
               col="black",
               size = 0.5,
               binwidth = 0.004) +
  geom_point(aes(x = -0.1, y = 0.2), shape = 18, size = 7, col = "black") +
  ggtitle("2D Gaussian, L = 10") +
  ylab(expression(mu[2])) +
  xlab(expression(mu[1])) +
  ylim(-0.4,0.4) +
  xlim(-0.4,0.4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.position="none")

# Base plot HMC
HMC_plot

# Helper function to create a data frame
# for plotting segments in ggplot
# from HMC results
create_df_segments <- function(Q){

  # empty vector to store coord
  x <- c()
  y <- c()
  xend <- c()
  yend <- c()
  K0 <- c() # kinetic energy

  # loop over HMC results
  for (j in 1:L) {
    x[j] <- Q$traj[j,1]
    y[j] <- Q$traj[j,2]
    xend[j] <- Q$traj[j+1,1]
    yend[j] <- Q$traj[j+1,2]
    K0[j] <- sum(Q$ptraj[j,]^2)/2
  }

  # Create tibble with results
  df_seg <- tibble(
    x = x,
    y = y,
    xend = xend,
    yend = yend,
    seg = c(t(outer("segment", 1:11, paste, sep = '_'))),
    K0 = K0
  )

  return(df_seg)
}

# example
# Store results of HMC
Q <- list()
# Initial value
Q$q <- c(-0.1,0.2)
# Step size
step <- 0.03
# Number of steps
L <- 11
# Number of iterations
n_samples <- 4

# Add HMC samples to base plot
for (i in 1:n_samples) {

  Q <- HMC(U, U_gradient, step, L, Q$q)

  df_seg <- create_df_segments(Q)

  if(i == 1){
    HMC_plot_final <- HMC_plot +
      geom_segment(data = df_seg, aes(x = x, y = y, xend = xend, yend = yend,
                                      linewidth = 0.2 + K0),
                   col = "skyblue") +
      geom_point(data = df_seg[-1,], aes(x = x, y = y), col = "red") +
      geom_point(data = df_seg[-c(1:10),], aes(x = xend, y = yend), col = "black",
                 size = 5)
  }
  if(i > 1){
    HMC_plot_final <- HMC_plot_final +
      geom_segment(data = df_seg, aes(x = x, y = y, xend = xend, yend = yend,
                                      linewidth = 0.2 + K0),
                   col = "skyblue") +
      geom_point(data = df_seg[-1,], aes(x = x, y = y), col = "red") +
      geom_point(data = df_seg[-c(1:10),], aes(x = xend, y = yend), col = "black",
                 size = 5)
  }
}

# Final plot - HMC
HMC_plot_final

# FIGURE 2 - Re-drawn the points & Add number to sample points
# Re-drawn the points
#df_points <- tibble(
#  x = c(0.05676668, 0.26030178, -0.06749449, -0.1410822, 0.34581499),
#  y = c(0.10346451, 0.04858683, -0.23142964, -0.3517979, 0.01544388)
#)
#
# HMC_plot_final <- HMC_plot_final +
#  geom_point(aes(x = -0.1, y = 0.2), shape = 18, size = 7, col = "black") +
#  geom_point(data = df_points, aes(x = x, y = y), col = "black",
#             size = 5)

# Add number to sample points
# Figure 2
# HMC_plot_final <- HMC_plot_final +
#  annotate("text", x = -0.12, y = 0.24, label = "0", size = 5) +
#  annotate("text", x = 0.07, y = 0.14, label = "1", size = 5) +
#  annotate("text", x = 0.260, y = 0.1, label = "2", size = 5) +
#  annotate("text", x = -0.09, y = -0.20, label = "3", size = 5) +
#  annotate("text", x = -0.115, y = -0.38, label = "4", size = 5) +
#  annotate("text", x = 0.345, y = -0.035, label = "5", size = 5)

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
# Figure 1
dbHb %>% ggplot(aes(x = Group, y = HMabs, fill = Group)) +
  geom_boxplot() +
  labs(x = "Group", y = "Hemoglobin mass (g)") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14))

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

# Figure 4
ggarrange(Intercept_prior,
          ChangeWtr_prior,
          IHE_prior,
          LHTL_prior,
          sigma_prior,
          nrow = 3,
          ncol = 2)

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
prior_summary(bmod1_prior)

# get parameters names
variables(bmod1_prior)

# Set brightblue scheme
color_scheme_set("brightblue")

# Prior predictive distribution
# Figure 5
bmod1_prior %>%
  posterior_predict(draws = 50) %>%
  ppc_dens_overlay(y = dbHb$HMabs) +
  xlim(-500, 500) +
  theme_bw()

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
# Figure 6
png('Figure4.png', units="in", width=6, height=5, res=600)
bmod1 %>%
  plot(combo = c("hist", "trace"), widths = c(1, 1.5),
       theme = theme_bw(base_size = 10))
dev.off()

# Posterior Predictive Distribution
# Figure 7
ppc_dens_overlay(y = dbHb$HMabs, yrep = posterior_predict(bmod1, draws = 50)) +
  xlim(-150, 150) +
  theme_bw()

## Model selection
bmod2 <- brm(formula = HMabs ~ 0 + Intercept + ChangeWtr + Group,
             data = dbHb,
             family = student(link = "identity"),
             warmup = 1000,
             iter = 2000,
             chains = 4,
             seed = 1234,
             #prior = bmod1Priors,
             sample_prior = "yes")

# summary(bmod2)

# PSIS-LOO ------------------------------------------------------------------
loo1 <- loo(bmod1, save_psis = TRUE)
# print(loo1)
# plot(loo1, label_points = T)
loo2 <- loo(bmod2, save_psis = TRUE)
# Table 1
loo_compare(loo1, loo2)

# ROPE analysis - bayestestR --------------------------------------------------------
describe_posterior(bmod1, rope_range = c(-0.5, 0.5), ci_method = "HDI", ci = 1)
equivalence_test(bmod1, range = c(-0.5, 0.5), ci = 1)

# Figure 8
bmod1 %>%
  gather_draws(b_Intercept, b_GroupIHE, b_GroupLHTL) %>%
  ggplot(aes(y = .variable, x = .value, fill = after_stat(abs(x) < 0.5))) +
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

# Create table 2
tbl_a <- summary(bmod1)
tbl_b <- describe_posterior(bmod1, rope_range = c(-0.5, 0.5))
tbl_c <- equivalence_test(bmod1, range = c(-0.5, 0.5))
tbl2 <- data.frame(tbl_a$fixed, tbl_b$pd, tbl_b$ROPE_Percentage, tbl_c$ROPE_Equivalence)
rownames(tbl2) <- c("alpha", "beta_1", "beta_2-IHE", "beta_2-LHTL")
colnames(tbl2) <- c("mean","SE", "lower bound", "upper bound", "Rhat", "Bulk ESS", "Tail ESS", "pd",
                    "% in ROPE", "H0")

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

# Sensitivity analysis function
# https://github.com/JorgeDelro/sensitivity_analysis
source("Bayes_sensitivity_analysis.R")

sensitivity_results <- sensitivity_analysis(bmodels = list(original_model = bmod1,
                                                  alternative_prior = bmod3),
                                   params = c("b_Intercept", "b_ChangeWtr",
                                              "b_GroupIHE", "b_GroupLHTL",
                                              "sigma"))
# Table
sensitivity_results$results_tables

# Figure 9
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


