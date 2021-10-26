rm(list = ls())

# Purpose: ----
# Perform monte-carlo simulation to test if an indicator variable
# can help mitigate ovb when the bias is only on a few observations;
# for our purposes we know that the bias occurs likely among low levels
# of the covariate.

# Example: The CFS data does not differentiate between domestic and 
# international flows for some of its state dyads. If we know which 
# state dyads may be compromise, then we can simply use an indicator
# to control for unwanted variation. For our purposes, we know the 
# measurement error may occur in LA and WA where there is also low
# levels of gdp.

# Conclusion: Indicator can indeed improved our results but not at 
# the 95% confidence interval. No including an indicator, however, 
# is way worse.

n <- 10000
sims <- 10000

t.stat_boolean <- function(b, SE) {# Function to calculate confidence
  
  boolean <- 1.96 > abs((b - 1)/SE)
  
  return(boolean)
  
}

df_out <- data.frame()
for(s in 1:sims) {
  set.seed(s*123) # Moving seed across simulations
  
  z <- rnorm(n, mean = 1, sd = 1) # The Omitted variable
  x <- rnorm(n, mean = 1, sd = 1) + 0.25*z # x being correlated with z
  
  b1 <- 1 # The coefficient we are interested in
  y <- 1 + b1*x + rnorm(n, mean = 0, sd = 1)
  
  # Benchmark (1)
  out1 <- lm(y~x)
  b1   <- out1$coefficients[2]
  SE1  <- sqrt(diag(vcov(out1)))[2]
  val1 <- t.stat_boolean(b1,SE1)
  
  # Tweaks fo add the measurement error (OVB) for 25% lowest of x
  indices <- x %in% head(sort(x, decreasing = FALSE), 0.25*n)
  
  youts <- y
  youts[indices] <- y[indices] + .25*z[indices] # adding the measurement error
  
  ones <- rep(1, n)
  ones[!indices] <- 0 # Dummy variable to capture only the ones being affected
  
  # No indicator (2)
  out2 <- lm(youts~x)
  b2   <- out2$coefficients[2]
  SE2  <- sqrt(diag(vcov(out2)))[2]
  val2 <- t.stat_boolean(b2,SE2)
  
  # Indicator (3)
  out3 <- lm(youts~x + ones)
  b3   <- out3$coefficients[2]
  SE3  <- sqrt(diag(vcov(out3)))[2]
  val3 <- t.stat_boolean(b3,SE3)
  
  # Should we drop outliers (4)
  out4 <- lm(youts[!indices]~x[!indices])
  b4   <- out4$coefficients[2]
  SE4  <- sqrt(diag(vcov(out4)))[2]
  val4 <- t.stat_boolean(b4,SE4)
  
  row <- c(b1, val1, b2, val2, b3, val3, b4, val4)
  
  df_out <- rbind.data.frame(df_out, row)
  
}

names(df_out) <- c("b_unbiased", "Confidence1", 
                   "b_noindica", "Confidence2",
                   "b_indicato", "Confidence3",
                   "b_dropping", "Confidence4")
if(TRUE) {
cat("\n")
cat(":::: The level of confidence for each specification is ::::")
cat("\n")
cat("Unbiased:", sum(df_out$Confidence1)/sims, "\n")
cat("Biased (no indicator):", sum(df_out$Confidence2)/sims, "\n")
cat("Biased (indicator):", sum(df_out$Confidence3)/sims, "\n")
cat("Biased (dropping):", sum(df_out$Confidence4)/sims, "\n")
}
#end