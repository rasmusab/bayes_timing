# This is a script implementing the model described in Baath (in preparation). The model
# below is the one descibed in the supplementary material which assumes that the 
# SD of the asynchrony is modeled as coming from a n-variate normal distribution where n
# is the number of ISI levels. This allows information regarding a participants performance 
# at one ISI level to inform the SD asynchrony at other ISI levels. The n-variate
# normal distribution is used to model the SDs on the log scale, that is
# log(sigma) ~ multivariate_norm(...).

# JAGS have a specific way of specifying censoring and a good guide is:
# http://doingbayesiandataanalysis.blogspot.se/2012/01/complete-example-of-right-censoring-in.html

# This script requres that you have JAGS and rjags installed, for information
# regarding JAGS see http://mcmc-jags.sourceforge.net/
library(rjags)
# We also use ddply for data manipulation
library(rjags)

# Reading in the data
d <- read.csv("tapping_data_baath_madison_2012.csv")

# Here the example data is from Baath and Madison (2012) 
# with the folowing variables/columns:

#   subject trial  isi asynch
# 1       1     1 3000   -383
# 2       1     1 3000   -302
# 3       1     1 3000    144
# 4       1     1 3000     -9
# ...     ...  ... ...    ...

# The right censoring limit in ms
cens = 100


# Calculating the mean and SD asynchrony per subject and isi using robust estimators (median and MAD).
# to be used to initialize the MCMC chains.
d_sum <- ddply(d, c("subject", "isi"), summarise, mean_asynch = median(asynch), sd_asynch = mad(asynch))
mean_asynch_per_subject <- daply(d_sum, "subject", function(d_sum_subj) {
  mean_asynch <- d_sum_subj$mean_asynch[order( d_sum_subj$isi )]
  names(mean_asynch) <- d_sum_subj$isi
  mean_asynch
}) 
sd_asynch_per_subject <- daply(d_sum, "subject", function(d_sum_subj) {
  sd_asynch <- d_sum_subj$sd_asynch[order( d_sum_subj$isi )]
  names(sd_asynch) <- d_sum_subj$isi
  sd_asynch
}) 


# Calculating which asynchronies that should be censored...
d$is_cens <- d$asynch > cens
# ... and creating a new variable where the censored asynchronies are set to NA.
d$asynch_cens <- ifelse(d$is_cens, NA, d$asynch)

# Generating the list of data to be sent to JAGS
isi_levels <- sort(unique(d$isi))
subjects <- sort(unique(d$subject))

data_list = list( 
  max_mean_sigma = isi_levels / 2,
  min_mean_mu = -isi_levels / 2,
  max_mean_mu = isi_levels / 2,
  max_sd_mu = isi_levels / 2,
  
  asynch = d$asynch_cens,
  subject = match(d$subject, subjects),
  n_subject = length(subjects),
  cens = cens,
  is_cens =d$is_cens,
  isi_level = match(d$isi, isi_levels),
  n_isi = length(isi_levels),
  identity_matrix = diag(1, length(isi_levels), length(isi_levels))
)

# Specifying the initial values for the MCMC chains using estimators

inits_list = list(
  mu = tapply(d$asynch, list(d$subject, d$isi), median),
  mean_mu = colMeans(mean_asynch_per_subject),
  sd_mu = apply(mean_asynch_per_subject, 2, sd),
  log_sigma = log(tapply(d$asynch, list(d$subject, d$isi), mad)),
  mean_log_sigma = colMeans(log(sd_asynch_per_subject)),
  inv_cov_mat_log_sigma = solve(cov(log(sd_asynch_per_subject)))
)


if(any(d$is_cens)) {
  # Initialize the censored asynchronies to be above the censoring threshold
  inits_list$asynch <- ifelse(d$is_cens, cens + 1, NA)
}

# The model specified in the JAGS language
model_string <- "model {
  for ( i in 1:length(asynch) ) {
    asynch[i] ~ dnorm(mu[subject[i], isi_level[i]], tau[subject[i], isi_level[i]])
    is_cens[i] ~ dinterval(asynch[i], cens)
  }

  for(subject_i in 1:n_subject) {
    log_sigma[subject_i, 1:n_isi] ~ dmnorm(mean_log_sigma[1:n_isi], inv_cov_mat_log_sigma[1:n_isi,1:n_isi])
    for(isi_i in 1:n_isi) {
      tau[subject_i, isi_i] <- 1 / pow(exp(log_sigma[subject_i, isi_i]), 2)
      mu[subject_i, isi_i] ~ dnorm(mean_mu[isi_i], 1 / pow(sd_mu[isi_i], 2) )
    }
  }

  # Specifying the hyper-priors
  for(isi_i in 1:n_isi) {
    mean_log_sigma[isi_i] ~ dunif(log(1), log(max_mean_sigma[isi_i]))
    mean_mu[isi_i]        ~ dunif(min_mean_mu[isi_i], max_mean_mu[isi_i])
    sd_mu[isi_i]          ~ dunif(0,max_sd_mu[isi_i])
  }

  inv_cov_mat_log_sigma ~ dwish(identity_matrix, n_isi)
  log_sigma_dist[1:n_isi] ~ dmnorm(mean_log_sigma, inv_cov_mat_log_sigma)
  for(isi_i in 1:n_isi) {
    sigma_dist[isi_i] <- exp(log_sigma_dist[isi_i])
  }
}"

# Fit the model using JAGS
# The parameters to be monitored
parameters = c("mu" , "log_sigma", 
                "mean_mu", "mean_log_sigma", "sigma_dist",
                "sd_mu", "inv_cov_mat_log_sigma")     
# Create, initialize, and adapt the model:
m = jags.model(textConnection(model_string), data=data_list, inits=inits_list, n.chains = 2 , n.adapt = 1000)
# Burn-in:
update(m, n.iter = 100 )
# Sampling from the model
s = coda.samples(m, parameters, n.iter = 1000)

# Looking at the posterior correlation matrix for the asynchrony SD.
smat <- as.matrix(s)
points_with_mean <- function(x, y, ...) {
  points(x, y, ...)
  points(mean(x), mean(y), col="red", pch=15)
}
pairs(smat[, c("sigma_dist[1]", "sigma_dist[2]", "sigma_dist[3]", "sigma_dist[4]")], 
      log="xy", panel=points_with_mean, col=rgb(0, 0, 0, 0.2),
      labels = paste("ISI", isi_levels[1:4]), xlim=c(15, 550), ylim=c(15, 550)) 

# Plotting point estimates for the group SD asynchrony
ssum <- summary(s)  
est <- ssum$quantiles

plot(isi_levels, est[grep(pattern = "^sigma_dist", rownames(est)), "50%"], type="b", ylim=c(0, 400), col="darkred", lwd=2)
for(i in seq_along(isi_levels)) {
  lines( rep(isi_levels[i], 2), est[paste0("sigma_dist[", i, "]"), c("25%", "75%")], lwd=4, col="darkred")
  lines( rep(isi_levels[i], 2), est[paste0("sigma_dist[", i, "]"), c("2.5%", "97.5%")], lwd=1, col="darkred")
}
