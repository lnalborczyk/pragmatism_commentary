library(tidyverse)
library(brms)

set.seed(666)
nsims <- 1e3

lower <- numeric()
upper <- numeric()

for(i in 1:nsims) {
    
    y <- data.frame(y = rnorm(1e2, mean = 100, sd = 15) )
    
    if(i == 1) {
        
        m <- brm(y ~ 1, data = y)
        
    } else {
        
        m <- update(m, newdata = y)
        
    }
    
    res <- fixef(m)[, 3:4] %>% as.vector
    lower[i] <- res[1]
    upper[i] <- res[2]
    
}

sim <- data.frame(lower = lower, upper = upper, id = 1:length(lower) ) %>%
    mutate(group = ifelse(100 > lower & 100 < upper, 1, 0) )

sim %>%
    ggplot(aes(x = id, colour = as.factor(group) ) ) +
    geom_hline(yintercept = 100, lty = 2, alpha = 0.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.25) +
    scale_y_continuous(limits = c(80, 120) ) +
    scale_x_continuous(breaks = c(seq(0, nsims, nsims / 10) ) ) +
    scale_color_discrete(guide = FALSE) +
    ylab("Estimate") +
    xlab("Simulation number") +
    theme_minimal(base_size = 12) +
    ggtitle(paste0("Coverage = ", mean(sim$group) ) )
