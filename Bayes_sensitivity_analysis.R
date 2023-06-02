
#
library(ggplot2) # ggplot2()
library(tibble) # tibble()
library(dplyr) # %>% / pull()
library(tidyr)
library(posterior) # as_draws_df() / summarise_draws()
library(RColorBrewer) # scale_fill_brewer()

# Helpers functions

# Check that the parameters were estimated in all the models
check_params <- function(bmodels,
                         params){


  for (i in seq_along(bmodels)) {

    # Get estimated parameter of each model
    model_parameters <- tidybayes::get_variables(bmodels[[i]])

    for (j in seq_along(params)) {

      # Check if the parameter is in the model
      if(!(params[j] %in% model_parameters)){
        # Stop if a parameter have not been estimated in one model
        stop(paste(params[j], "is not in model", names(bmodels[i]), sep = " "))
      }
    }
  }

} # end check_params

# Create tables to store results
create_tables <- function(bmodels,
                          params){

  list_tables <- list()
  # loop over models
  for (i in seq_along(bmodels)) {

    if(i == 1){

      # Get original draws
      original_draws <- as_draws_df(bmodels[[i]])

      original_means <- c()
      for (j in seq_along(params)) {
        # loop over parameters
        original_means[j] <- summarise_draws(original_draws, "mean") %>%
          select("variable", "mean") %>%
          filter(variable == params[j]) %>%
          pull()
      } # end loop params-original model

      next # go to the second model
    } # end get original model info


    # Initialize variables to obtain
    mod_means <- c()
    mod_sds <- c()
    quantile_5 <- c()
    quantile_50 <- c()
    quantile_95 <- c()
    p_deviation <- c()

    # Loop over the parameters
    for (j in seq_along(params)) {

      mod_means[j] <- summarise_draws(bmodels[[i]], "mean") %>%
        select("variable", "mean") %>%
        filter(variable == params[j]) %>%
        pull()

      mod_sds[j] <- summarise_draws(bmodels[[i]], "sd") %>%
        select("variable", "sd") %>%
        filter(variable == params[j]) %>%
        pull()

      # Percentage Deviation = (Mean - Original Mean)/Original Mean * 100
      p_deviation[j] <- (mod_means[j] - original_means[j]) / original_means[j]*100

      quantile_5[j] <- summarise_draws(bmodels[[i]], function(x) quantile(x, probs = c(0.05))) %>%
        select("variable", "5%") %>%
        filter(variable == params[j]) %>%
        pull()

      quantile_50[j] <- summarise_draws(bmodels[[i]], function(x) quantile(x, probs = c(0.5))) %>%
        select("variable", "50%") %>%
        filter(variable == params[j]) %>%
        pull()

      quantile_95[j] <- summarise_draws(bmodels[[i]], function(x) quantile(x, probs = c(0.95))) %>%
        select("variable", "95%") %>%
        filter(variable == params[j]) %>%
        pull()
    } # end loop params-alternative models


    # Built table with results
    sensitivity_result <- tibble(parameter = params,
                                 mean = mod_means,
                                 sd = mod_sds,
                                 '5%' = quantile_5,
                                 '50%' = quantile_50,
                                 '95%' = quantile_95,
                                 'original mean' = original_means,
                                 'percentage deviation' = p_deviation)

    list_tables[[i-1]] <- knitr::kable(sensitivity_result, caption = names(bmodels)[i])

  } # end loop over bmodels

  return(list_tables)

} # end create tables


# Create plots
create_plots <- function(bmodels,
                         params){

  list_plots <- list()

  # Get all models draws
  list_draws <- list()
  for (i in seq_along(bmodels)) {

    list_draws[[i]] <- as_draws_df(bmodels[[i]])

  } # end loop bmodels

  # Loop over the parameters
  for (j in seq_along(params)) {

    # Get the selected parameter of every model
    for (k in seq_along(bmodels)) {

      if(k == 1){
        db_plot <- data.frame(extract_variable(list_draws[[k]], variable = params[j]))
      }
      if(k > 1){
        db_plot_ex <- data.frame(extract_variable(list_draws[[k]], variable = params[j]))
        db_plot <- cbind(db_plot, db_plot_ex)
      }

    } # end loop bmodels

    # Add models names
    names(db_plot) <- names(bmodels)

    db_plot <- db_plot %>%
      pivot_longer(cols = everything(),
                   names_to = c("model"),
                   values_to = "value")

    # Create plot
    density_plot <- db_plot %>%
      ggplot(aes(x = value, fill = model)) +
      geom_density(alpha = 0.5) +
      ggtitle(params[j]) +
      #geom_vline(data=mu, aes(xintercept=grp.mean, color=model), linetype="dashed") +
      scale_fill_brewer(palette="Dark2") +
      theme_bw() +
      theme(legend.title=element_blank(),
            plot.title=element_text(hjust=0.5),
            legend.position="bottom")


    list_plots[[j]] <- density_plot

  } # end loop params

  return(list_plots)

} # end create_plots


# Sensitivity analysis function
sensitivity_analysis <- function(bmodels,
                                 params){

  # Check that models is a named list
  if(is.null(names(bmodels))) {
    stop("bmodels must be a named list")
  }

  if(is.null(params)) {
    stop("users must specify at least one parameter to perform sensitivity analysis")
  }

  # Check that the parameters were
  # estimated in all the models
  check_params(bmodels, params)

  # Get tables
  results_tables <- create_tables(bmodels,
                                  params)

  # Get plots
  results_plots <- create_plots(bmodels,
                                  params)

  return(list(
    results_tables = results_tables,
    results_plots = results_plots
  ))



} # end sensitivity analysis

