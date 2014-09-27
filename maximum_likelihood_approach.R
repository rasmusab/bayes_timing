# This script documents an alternative method for estimating the mean and standard deviation
# of a asynchrony distribution when predictive responses above a certain limit (here in the 
# script set to 100 ms) are assumed to be contaminated by reactive responses.
# The this method fits a truncated normal distribution to the data and
# is described (in the context of reaction time distributions) by
# Ulrich, R., & Miller, J. (1994). 
# Effects of truncation on reaction time analysis. 
# Journal of Experimental Psychology: General, 123(1), 34.

# Librart that implements the truncated normal distribution...
library(truncnorm)

# First defining the function that estimates the parameters of a truncated normal
# using maximum likelihood
truncnorm_mle <- function(x, lower = -Inf, upper = Inf) {
  x_trunc <- as.numeric(x[x > lower & x < upper]) 
  opt <- optim(
    par = c(mean(x), sd(x)),
    fn = function(par) sum( log(dtruncnorm(x_trunc, lower, upper, mean = par[1], sd = par[2]))),
    control = list(fnscale = -1)) # to do maximization instead of minimization
  c(mean = opt$par[1], sd = opt$par[2])
}

# Estimating the mean and sd using the truncated normal maximum likelihood
truncnorm_mle(isi_3000, upper = cens)

# Estimating and plotting the sd for alla isi levels of participant 1
est_sd <- aggregate(d$asynch, list(isi = d$isi), function(x) {
  truncnorm_mle(x, upper = cens)["sd"]
})

plot(est_sd$isi, est_sd$x, ylim=c(0, max(est_sd$x)), xlim=c(0, max(est_sd$isi)), 
     type="b", xlab = "ISI in ms", ylab = "SD in ms")
