##############################################################################
# R code accompanying Nalborczyk, Bürkner, & Williams (2018)
# OSF projet: https://osf.io/nmp6x/
# Contact: ladislas.nalborczyk@gmail.com
# Last update: September 25, 2018
#############################################################

#####################################################
# regression example
#####################################

library(tidyverse)
library(brms)

set.seed(666) # setting the seed for reproducibility
nsims <- 1e4 # number of simulations

lower <- numeric()
upper <- numeric()

for(i in 1:nsims) {
    
    y <- data.frame(y = rnorm(n = 1e2, mean = 0, sd = 10) )
    
    if(i == 1) {
        
        m <- brm(y ~ 1, data = y)
        
    } else {
        
        m <- update(m, newdata = y)
        
    }
    
    intervals <- fixef(m)[, 3:4] %>% as.vector
    
    lower[i] <- intervals[1]
    upper[i] <- intervals[2]
    
    print(i)
    
}

sim <- 
    data.frame(lower = lower, upper = upper, id = 1:length(lower) ) %>%
    mutate(group = ifelse(0 > lower & 0 < upper, 1, 0) )

png(filename = "coverage1.png", width = 4500, height = 3000, res = 300)

sim %>%
    # selecting a subset of observations for readability
    # sample_frac(0.5) %>%
    ggplot(aes(x = id, colour = as.factor(group), size = as.factor(group) ) ) +
    # representing the population value of the intercept
    geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
    # plotting the intervals
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
    # ylim
    scale_y_continuous(limits = c(-8, 8) ) +
    # x-axis ticks
    scale_x_continuous(limits = c(0, 100), breaks = c(seq(0, 100, 10) ) ) +
    # scale_x_continuous(breaks = c(seq(0, nsims, nsims / 10) ) ) +
    # defining colors
    scale_color_manual(values = c("steelblue", "gray80"), guide = FALSE) +
    # and width of the intervals
    scale_size_manual(values = c(1, 0.75), guide = FALSE) +
    ylab("Estimate") +
    xlab("Simulation number") +
    theme_minimal(base_size = 20) +
    # representing the coverage proportion
    ggtitle(paste0("Coverage = ", round(mean(sim$group), 3) ) )

dev.off()

####################################################################
# meta-analysis example
################################################

library(metaBMA)
library(metafor)

set.seed(666) # setting the seed for reproducibility

nsims <- 1e4 # number of simulations
tau <- 0.1 # population value of tau

# initialising an empty dataframe
res <- data.frame(
    lower_freq = numeric(), upper_freq = numeric(),
    lower_bayes = numeric(), upper_bayes = numeric()
    )

for(i in 1:nsims) {
    
    print(i) # print current simulation run
    
    y <- data.frame(
        # generating data for a population d = 0.2 and a given tau
        effsize = rnorm(6, 0.2, tau),
        SE = runif(6, 0.05, 0.1),
        study = 1:6
        )
    
    # fitting a frequentist random-effects meta-analysis model
    fit_freq <- rma(effsize, sei = SE, data = y, method = "PM")
    
    lower_freq <- confint(fit_freq, level = 0.95)$random[2, 2]
    upper_freq <- confint(fit_freq, level = 0.95)$random[2, 3]
    
    # fitting a Bayesian random-effects meta-analysis model
    fit_bayes <- meta_random(
        y$effsize, y$SE, y$study,
        d = "norm", d.par = c(0, 10000),
        tau = "halfcauchy", tau.par = 10000,
        sample = 500, summarize = "jags"
        )
    
    lower_bayes <- fit_bayes$estimates[2, 4]
    upper_bayes <- fit_bayes$estimates[2, 5]
    
    res[i, ] <- cbind(lower_freq, upper_freq, lower_bayes, upper_bayes)
    
}

png(filename = "coverage2.png", width = 4500, height = 3000, res = 300)

res %>%
    gather(interval, value) %>%
    separate(col = interval, into = c("interval", "type"), sep = "_", remove = TRUE) %>%
    mutate(id = rep(1:(nrow(.) / 4), 4) ) %>%
    spread(interval, value) %>%
    # identifying misses and hits
    mutate(group = ifelse(tau > lower & tau < upper, 1, 0) ) %>%
    # computing coverage
    group_by(type) %>%
    mutate(coverage = mean(group) ) %>%
    # selecting a subset of observations for readability
    # sample_frac(0.1) %>%
    # plotting
    ggplot(aes(x = id, colour = as.factor(group), size = as.factor(group) ) ) +
    # plotting the population value of tau
    geom_hline(yintercept = tau, lty = 2, alpha = 0.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
    # y-axis
    scale_y_continuous(limits = c(0, 0.75) ) +
    # x-axis ticks
    scale_x_continuous(limits = c(0, 100), breaks = c(seq(0, 100, 10) ) ) +
    # scale_x_continuous(breaks = c(seq(0, nsims, nsims / 10) ) ) +
    # colors of intervals
    scale_color_manual(values = c("steelblue", "gray80"), guide = FALSE) +
    # width of intervals
    scale_size_manual(values = c(1, 0.75), guide = FALSE) +
    # facetting by type of interval
    facet_wrap(
        ~type, ncol = 1, scales = "free",
        labeller = as_labeller(c(
            "freq" = "Frequentist confidence interval",
            "bayes" = "Bayesian credible interval"
            ) )
        ) +
    # adding coverage proportion
    geom_label(
        aes(x = 90, y = 0.7, label = paste0("Coverage = ", round(coverage, 3) ) ),
        inherit.aes = FALSE, size = 6
        ) +
    ylab("Estimate") +
    xlab("Simulation number") +
    theme_minimal(base_size = 20)

dev.off()
