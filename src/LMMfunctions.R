# simple plot to visualize lmm interactions
simple_interaction <- function(xvar, yvar, tracevar, data) {
  ggplot(data, 
         aes(x = xvar,
             y = yvar,
             group = tracevar,
             linetype = tracevar)) +
    stat_summary(fun = "mean", # connect means
                 geom = "line")
}

# z-score a variable in dplyr pipe
# @vrbl (symbol, required): the variable
my_z_score <- function(vrbl) {
  vrbl_z = (vrbl - mean(vrbl, na.rm = TRUE)) / sd(vrbl, na.rm = TRUE)
  return(vrbl_z)
}

# Create summary tables (Mean & SD)
summary_tab <- function(df, dvs, factorstring) {
  
  factor = sym(factorstring)
  
  summary_group <- df %>%
    group_by(!!factor) %>%
    summarise(across(all_of(dvs), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)), .names = "{.col}-{.fn}"))
  
  
  summary_total = df %>%
    ungroup() %>%
    summarise(across(all_of(dvs), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)), .names = "{.col}-{.fn}")) %>%
    mutate(!!factorstring := "total") %>%
    rbind(summary_group) %>%
    pivot_longer(cols = ends_with(c("mean", "sd")), names_to = c("variable", ".value"), names_sep = "\\-") %>%
    pivot_wider(id_cols = variable, names_from = !!factor, values_from = c(mean, sd), names_vary = "slowest") %>%
    mutate_if(is.numeric, round, 2)
  
  # TODO: add independent t-tests for comparison between factor levels
  
  return(summary_total)
  
  # print as apa table
  # apa_table(
  # summary_total,
  # caption = paste(deparse(substitute(fun)), "of", dv),
  # escape = TRUE
  # )
}


# Create cross-correlation tables
cross_corrs <- function(df, vars, filtervar=c(), filterlevel=c()) {
  
  if(!is.null(filtervar)){
    df <- df %>% 
      filter(!!sym(filtervar) == filterlevel) # filter df
  } 
  
  corrtest_list <- corr.test(df[, vars], method = "spearman")
  corr_matrix <- matrix(str_c("r=", round(corrtest_list$r, 2), ", t=", round(corrtest_list$t, 2), ", p=", round(corrtest_list$p, 3), sep=""), 
                        nrow = length(vars), 
                        ncol = length(vars), 
                        byrow = TRUE, 
                        dimnames = list(vars, vars))
}

# Compute correlations between dv and vars
# @dv (character, required): the name of the dependent variable
# @vars (character list, required): a list of names of variables to correlate with dv
# @data (symbol, required): the data frame
# @method (character): force this correlation method
# @output (character): If "table", function will print a markdown table of correlations >=.30. Else defaults to returning a df
dv_corr <- function(dv, vars, data, method = c(), output = "df") {
  df <- data
  corr_results <- list()
  
  for (var in vars) {
    if (is.null(method)) {
      # if var is dichotomous, use point-biserial correlation; if var is continuous, use spearman
      unique_values <- unique(df[[var]])
      my_method <- ifelse(length(unique_values) <= 2, "pearson", "spearman")
    } else {
      my_method = method
    }
    
    # correlation test
    p <- corr.test(df[[var]], df[[dv]], method=my_method, use = "pairwise.complete.obs")
    
    # format result
    tmp <- data.frame(N = p$n, # sample size
                      r = p$r, # correlation coefficient 
                      p_value = p$p, # p-value 
                      stars = ifelse(p$p < 0.001, "***",
                                     ifelse(p$p < 0.01, "**",
                                            ifelse(p$p < 0.05, "*", ""))) )# significance stars
    
    # append result to list
    corr_results[[var]] <- tmp
  }
  
  # bind corr.test results for all control vars in one dataframe
  corr_df <- do.call(rbind, corr_results)
  
  if (output %in% "table") {
    # filter for r >= .30
    corr_df_filtered <- corr_df %>%
      filter(p_value < .05)
    
    if(nrow(corr_df_filtered) > 0) {
      print(corr_df_filtered)
    } else {
      paste("No correlations r >= .30 between", dv, "and control variables")
    }
  } else {
    return(corr_df)
  } 
}

# Plot dots and bars for mean & CI 
plot_dotbars <- function(df, dv, by_v, ylab, xlab, my_colour) {
  
  # remove missings
  data <- df %>%
    filter(!is.na(!!sym(dv)))
  
  # sample size
  sample_size = data %>% 
    group_by(!!sym(by_v)) %>% 
    summarise(num=n()) %>%
    mutate(myaxis = paste0(!!sym(by_v), "\n", "n=", num))
  
  plot <- ggplot(data = data,
   aes(x = !!sym(by_v),
       y = !!sym(dv))) +
    geom_hline(yintercept = 0, # add zero change line
               linewidth = 0.8) + 
    stat_summary(fun = "mean", # add mean bars
                 geom = "bar",
                 fill = my_colour,
                 alpha = .5,
                 width = .5) +
    geom_point(position = position_jitter(width = .2, height = 0), # add dots for individual participants
               shape = 16,
               colour = "darkgrey",
               alpha = 0.5) + 
    stat_summary(fun.data = mean_cl_boot, # add 95% bootstrapped CI with 1000 samples
                 geom = "errorbar", 
                 width = 0.18,
                 linewidth = 0.5,
                 color = "#696969",
                 show.legend = FALSE) +
    coord_fixed(ratio = 0.25, 
                ylim = c(-4, 6), 
                clip = "off") +
    scale_x_discrete(label = sample_size[["myaxis"]]) +
    theme_bw(base_size = 11) +
    labs(x = xlab,
         y = ylab) +
    theme(# axis.title = element_text(size=9.5),
          axis.text.x = element_text(angle=20, hjust=0.8))
  
  return(plot)
}

# Plot means and CI by order across timepoint function
plot_meanSE <- function(df, dv, facet_var=NULL, ylab, my_colours, show_labels=F) {
  
  if(endsWith(dv, "_z")) {
    my_yaxis <- c(-2, 2) # y axis limits
    my_ratio <- 1.05
    }
  else if(startsWith(dv, "log_")) {
    my_yaxis <- c(0, 2)
    my_ratio <- 2.1
    my_nudge <- data.frame("nudge" = c(0.1, 0.1, -0.1, -0.1))
  } else {
    my_yaxis <- c(1, 7)
    my_ratio <- 0.6
    my_nudge <- data.frame("nudge" = c(0.3, 0.3, -0.3, -0.3))
  }
  
  # remove missings
  data <- df %>%
    filter(!is.na(!!sym(dv))) # drop rows where no rating was made
  
  # if dv has less than 10 unique values, jitter the individual participants' lines to avoid overlapping
  if (nrow(unique(data[dv])) < 10) {
    my_linejitter <- 0.2
  } else {
    my_linejitter <- 0
  }
  
  dv = sym(dv)
  
  splot <- data %>%
    ggplot(aes(x = timepoint, 
               y = !!dv, 
               group = order,
               color = order,
               linetype = order,
               shape = medium)) + 
    geom_line(aes(group = p_id),
              alpha = 0.5,
              position = position_jitter(width = 0, height = my_linejitter)) +
    scale_colour_manual(values = c("darkgrey", "lightgrey"), guide = "none") +
    new_scale_color() +
    stat_summary(aes(color = order),
                 fun = "mean", # connect means
                 geom = "line",
                 linewidth = 0.9, 
                 position = position_dodge(0.5)) + 
    stat_summary(aes(color = order),
                 fun = "mean", # add mean
                 geom = "point",
                 size = 3, 
                 position = position_dodge(0.5)) + 
    stat_summary(aes(color = order),
                 fun.data = mean_cl_boot, # add 95% bootstrapped CI with 1000 samples
                 geom = "errorbar", 
                 width = 0.35,
                 linewidth = 0.9,
                 linetype = "solid", 
                 position = position_dodge(0.5),
                 show.legend = FALSE) + 
    scale_colour_manual(values = my_colours) +
    guides(
      colour = guide_legend("Order", override.aes = list(shape = NA)),
      linetype = guide_legend("Order"),
      shape = "none"
    ) +
    scale_y_continuous(limits = my_yaxis) +
    theme_bw(base_size = 11) +
    labs(x = "Time point",
         y = ylab) +
    theme(panel.grid.major.x = element_blank(), # remove vertical grid line,
          # legend.text=element_text(size=6),
          # legend.title=element_text(margin = margin(b = 0)),
          legend.background = element_rect(fill = "white", color = "black"),
          # axis.title = element_text(size=9.5)
          )
  
  if (!is.null(facet_var)) {
    facet_var = sym(facet_var)
    
    splot <- splot +
      facet_wrap(vars(!!facet_var), scales = "free_x")
  } else {
    # only add sample size if not faceting x-axis
    sample_size = data %>% 
      group_by(timepoint) %>% 
      summarise(num=n()) %>%
      mutate(myaxis = paste0(timepoint, "\n", "n=", num))
    
    splot <- splot +
      scale_x_discrete(label = sample_size[["myaxis"]]) +
      theme(axis.text.x = element_text(angle=20, hjust=0.8)) # tilt x-axis labels
  }
  
  if (show_labels==T) {
    label_position <- data %>%
      filter(!timepoint == "Baseline") %>%
      group_by(order, timepoint, medium) %>%
      summarise(mean = mean(!!dv, na.rm = TRUE), 
                .groups = "keep") %>% 
      bind_cols(my_nudge) %>%
      mutate(nudged_mean = mean + nudge) # calculate label positions as nudged group mean
    
    splot <- splot +
      geom_label(data = label_position, # add labels for film and performance
                 aes(x = timepoint, y = nudged_mean, label = medium),
                 position = position_dodge(1.7),
                 size = 5
      )
  }
  
  return(splot)
}

# Filter df according to exclusion criteria and create a lmer object
# @dv (character, required): the name of the dependent variable
# @formula (character, required): the main predictors to be tested in the model
# @data (symbol, required): the data frame 
# @control_vars (character list, optional): a list of names of the variables to include as controls (add control vars that are correlated with dv by r >= .30)
# @filter_sd (character, optional): the name of the variable that stores the variance of the items that contribute to the dv scale
# @extra_vars (character list, optional): a list of names of variables to include as additional predictors

lmm_analysis <- function(dv, formula=c(), data, control_vars=c(), extra_vars=c()) {
  df <- data %>%
    filter(!is.na(!!sym(dv))) # drop rows where dependent variable is missing
  
  # if using intermission survey data, drop row with time_spent outlier   
  # if("timepoint_3" %in% colnames(df)){
  #   df <- df %>%
  #     filter(is.na(check_tooslow1) | !(check_tooslow1 == TRUE & timepoint_3 == "Intermission")) 
  # } 
  
  # format formulas
  control_formula <- paste(control_vars, collapse = " + ")
  if(!is.null(control_vars)) {
    control_formula <- paste0(control_formula, " + ")
    control_vars <- trimws(strsplit(control_vars, "\\+")[[1]])
  }
  
  if(!is.null(formula)) {
    formula <- paste0(formula, " + ")
  }
  
  extra_formula <- paste(extra_vars, collapse = " + ")
  if(!is.null(extra_vars)) {
    extra_formula <- paste0(extra_formula, " + ")
    extra_vars <- trimws(strsplit(extra_vars, "\\+")[[1]])
  }
  
  # concatenate formulas
  formulas <- paste0(control_formula, formula, extra_formula, "(1|p_id)")
  
  # impute means in control and extra variables
  var_cols <- c(dv, control_vars, extra_vars) # get names of all measured variables (as opposed to design variables) in the model
  for(column in var_cols) {
    if(column %in% colnames(df)) {
      cat("\n Number of NAs imputed in ", column, ":", length(df[which(is.na(df[ , column])), ][[column]]), "\n")
      if(length(df[which(is.na(df[ , column])), ][[column]]) > 0) {
        df[which(is.na(df[ , column])) , column] <- mean(df[[column]], na.rm=TRUE)
      }
    }
  }
  
  lmm_formula <- as.formula(paste(dv, "~", formulas))
  
  print(lmm_formula)
  
  lmm_model <- lmer(lmm_formula, data=df)
  
  return(lmm_model)
}

# Compute likelihood ratio test of two nested models (lmer, lm, or polr), return result in apa format
# @model1 (symbol): The model that serves as comparison
# @model2 (symbol): The model that is compared against model1
# @greeks (BOOL): Whether to use latex greek characters
model_comp <- function(model1, model2, greeks = TRUE) {
  suppressMessages(tmp_chi2 <- anova(model1, model2))
  
  pval_column <- grep("^Pr", colnames(tmp_chi2), value = TRUE)
  df_column <- setdiff(grep("Df", colnames(tmp_chi2), value = TRUE), "Res.Df") # lm models have "Res.Df" and "Df"
  statistic_column <- ifelse("Chisq" %in% colnames(tmp_chi2), "Chisq", # lmerTest models
                             ifelse("LR stat." %in% colnames(tmp_chi2), "LR stat.", # polr models
                                    "F")) # lm models
  
  if (tmp_chi2[2, pval_column] < 0.001) {
    my_pval <- "$p < .001$"
  } else {
    my_pval <- paste0("$p = ", round(tmp_chi2[2, pval_column], 3), "$")
  }

  if (greeks) {
    # Create the LaTeX-compatible APA formatted result
    formatted_result <- paste0(
      "$\\chi^2 (", round(tmp_chi2[2, df_column], 0), ") = ", round(tmp_chi2[2, statistic_column], 2), 
      "$, ", my_pval
    )
  } else {
    # Create the APA formatted result without greek characters 
    formatted_result <- paste0(
      "$chi^2 (", round(tmp_chi2[2, df_column], 0), ") = ", round(tmp_chi2[2, statistic_column], 2), 
      "$, ", my_pval
    )
  } 
  
  
  # Return formatted result
  return(formatted_result)
}

# Compute Cohen's d for fixed effects in LMM
# @model (symbol): The lmerModLmerTest object
my_effectsize <- function(model) {
  my_effectnames <- rownames(coef(summary(model)))
  my_randvar <- VarCorr(model)$p_id[1]
  my_residvar <- (attr(VarCorr(model), "sc"))^2
  my_effectsizes <- list()
  
  for (effect in my_effectnames) {
    tmp_estimate <- coef(summary(model))[effect, "Estimate"]
    tmp_d <- (tmp_estimate*2) / sqrt(my_randvar + my_residvar) # estimate*2 = difference in means for sum-to-zero contrasts c(1, -1)
    
    my_effectsizes[[effect]] <- tmp_d
  }
  
  df_effectsizes <- do.call(rbind, my_effectsizes)
  colnames(df_effectsizes) <- c("CohensD")
  df_effectsizes <- abs(df_effectsizes)
  
  return(df_effectsizes)
}

# Print out the evaluation of a fixed effect using Type III Anova
# @model (character, required): The name of the lmerModLmerTest object
# @var (character, required): The name of the variable to be evaluated
eval_fixed <- function(model, var) {
  tmp_coeffs <- coef(summary(model))
  tmp_effects <- my_effectsize(model) 
  
  indx <- which(rownames(tmp_coeffs)==var) # get index of variable 
  
  if (length(indx) == 0) {
    stop("Variable not found in model.")
  } 
  
  tmp_est <- tmp_coeffs[indx, "Estimate"]
  tmp_tval <- tmp_coeffs[indx, "t value"]
  tmp_df <- tmp_coeffs[indx, "df"]
  tmp_pval <- tmp_coeffs[indx, "Pr(>|t|)"]
  tmp_cohen <- tmp_effects[indx]
  
  if (tmp_pval < 0.001) {
    my_pval <- "$p < .001$"
  } else {
    my_pval <- paste0("$p = ", round(tmp_pval, 3), "$")
  }
  
  # Create the LaTeX-compatible APA formatted result
  formatted_result <- paste0(
    "$\\beta = ", round(tmp_est, 2), 
    "$, $t (", round(tmp_df, 0), ") = ", round(tmp_tval, 2), 
    "$, ", my_pval,
    ", $d = ", round(tmp_cohen, 2), "$"
  )
  
  # Return formatted result
  return(formatted_result)
  
  # cat("$\\chi^2(`r printnum(tmp_anova$Df[indx_anova], digits = 0)`) = `r printnum(tmp_anova$Chisq[indx_anova], digits = 2)`,\\ p `r printp(tmp_anova$Pr[indx_anova])`,\\ \\omega^2 = `r printnum(tmp_effects$Omega2_partial[indx_effects], digits = 2)`$") # print the result in markdown maths
}


# Check model assumptions

# * Normality of residuals with qqnorm(residuals(problem_awareness_model)) and shapiro.test(residuals(problem_awareness_model))  
# * Normality of random effect: shapiro.test(random effect estimate)  
# * Homogeneity: plot(ios_anna_model) and bartlett.test(score~cells,df) or leveneTest(residuals(model) ~ df$cells), see https://stats.stackexchange.com/a/76228  
# * if there is an issue, try log(outcome) or reciprocal outcome: 1/outcome  
# cf. https://cran.r-project.org/doc/contrib/Faraway-PRA.pdf (chapter 16.2.3) and https://doi.org/10.1111/j.2041-210X.2009.00001.x
# * Multicollinearity: Myers (1990) suggests that a value of 10 is a good value at which to worry (Field, Discovering Statistics)

check_assumptions <- function(model, dv, df) {
  op <- par(mfrow = c(1, 3)) # set 2x2 pictures on one plot
  
  # normality of residuals
  qqnorm(residuals(model))
  qqline(residuals(model))
  title("Check normality of residuals:", line = 0.3)
  print(shapiro.test(residuals(model)))
  
  # normality of random effect
  if (class(model) == "lmerModLmerTest") {
    rand_intercepts <- ranef(model)$p_id$`(Intercept)`
    qqnorm(rand_intercepts)
    qqline(rand_intercepts)
    title("Check normality of random intercepts:", line = 0.3)
    print(shapiro.test(rand_intercepts))
  }
  
  # for homogeneity
  plot(fitted(model), residuals(model), xlab="Fitted", ylab="Residuals")
  abline(h = 0, col = "red")
  title("Check homogeneity of variances:")
  my_formula <- formula(paste0(dv, " ~ cell"))
  print(bartlett.test(my_formula, data = df))
  
  # for multicollinearity
  cat("\n       Variance Influence Factors: \n\n")
  print(vif(model, type="predictor"))
  
  par(op) # reset plotting window
}

# get OR and CIs from a logistic regression model
# @model (required): a polr model (from MASS:)
get_OR_CI <- function(model) {
  my_coefficients <- summary(model)$coefficients # extract all coefficients
  my_ORs <- exp(my_coefficients[, 1]) # compute odds ratios = exponents of coefficients
  CI_lower <- exp(my_coefficients[, 1] - 1.96 * my_coefficients[, 2]) # Compute 95% confidence intervals
  CI_upper <- exp(my_coefficients[, 1] + 1.96 * my_coefficients[, 2])
  
  # combine into dataframe
  my_table <- data.frame(
    OR = my_ORs,
    CI_lower = CI_lower,
    CI_upper = CI_upper
  )
  
  return(my_table)
}

calculate_scale_score <- function(items, unique_column, cutoff) {
  # @items (character vector, required): names of items to average over
  # @unique_column (character string, required): name of column that stores the number of unique values across these items
  # @cutoff (numeric): number of unique values that would be considered a lack of variance
  
  # calculate percentage of items that have been answered
  n_items <- length(items)
  valid_items <- sum(!is.na(items))
  threshold <- 0.8 * n_items
  
  if (valid_items >= threshold & unique_column > cutoff) {
    return(mean(items, na.rm = TRUE))
  } else {
    return(NA)
  }
}

# function for annotating faceted plots from https://stackoverflow.com/a/44897816
annotation_custom2 <- 
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){ 
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob, 
                                            xmin = xmin, xmax = xmax, 
                                            ymin = ymin, ymax = ymax))
  }