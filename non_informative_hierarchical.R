# This is a script implementing the model described in Baath (in preparation). The model
# below is the one refered to as the Non-informative hierarchical model in the manuscript.

# JAGS have a specific way of specifying censoring and a good guide is:
# http://doingbayesiandataanalysis.blogspot.se/2012/01/complete-example-of-right-censoring-in.html

# This script requres that you have JAGS and rjags installed, for information
# regarding JAGS see http://mcmc-jags.sourceforge.net/
require(rjags)

# The right censoring limit in ms
cens =  100

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


# Calculating which asynchronies that should be censored...
d$is_cens <- d$asynch > cens
# ... and creating a new variable where the censored asynchronies are set to NA.
d$asynch_cens <- ifelse(d$is_cens, NA, d$asynch)

# Generating the list of data to be sent to JAGS
isi_levels <- sort(unique(d$isi))
subjects <- sort(unique(d$subject))

data_list = list(
  max_mean_sigma = isi_levels / 2,
  max_sd_sigma = isi_levels / 4,
  min_mean_mu = -isi_levels / 2,
  max_mean_mu = isi_levels / 2,
  max_sd_mu = isi_levels / 2,
  
  asynch = d$asynch_cens,
  subject = match(d$subject, subjects),
  n_subject = length(subjects),
  cens = cens,
  is_cens =d$is_cens,
  isi = match(d$isi, isi_levels),
  n_isi = length(isi_levels)
)

# Specifying the initial values for the MCMC chains using estimators

inits_list = list(
  mu = tapply(d$asynch, list(d$subject, d$isi), mean),
  sigma = tapply(d$asynch, list(d$subject, d$isi), sd)
)

inits_list$log_mean_sigma <- log(apply(inits_list$sigma, 2, mean))
inits_list$sd_sigma <- apply(inits_list$sigma, 2, sd)
inits_list$mean_mu <- apply(inits_list$mu, 2, mean)
inits_list$sd_mu <- apply(inits_list$mu, 2, sd)

if(any(d$is_cens)) {
  # Initialize the censored asynchronies to be above the censoring threshold
  inits_list$asynch <- ifelse(d$is_cens, cens + 1, NA)
}

# The model specified in the JAGS language
model_string <- "model {
  for ( i in 1:length(asynch) ) {
    asynch[i] ~ dnorm(mu[subject[i], isi[i]], tau[subject[i], isi[i]])
    is_cens[i] ~ dinterval(asynch[i], cens)
  }
  
  # Specifying the hyper-priors
  for(isi_i in 1:n_isi) {
    for(subject_i in 1:n_subject) {
      mu[subject_i, isi_i] ~ dnorm(mean_mu[isi_i], 1/pow(sd_mu[isi_i], 2))
      tau[subject_i, isi_i] <- 1 / pow(sigma[subject_i, isi_i], 2)
      sigma[subject_i, isi_i] ~ dlnorm(logmu_sigma[isi_i], 1/pow(logsigma_sigma[isi_i], 2))
    }
    
    # Transformation for the log-normal on sigma
    logsigma_sigma[isi_i] <- sqrt(log( pow((sd_sigma[isi_i]/mean_sigma[isi_i]), 2) + 1 ))
    logmu_sigma[isi_i] <- log(mean_sigma[isi_i]) - pow(logsigma_sigma[isi_i], 2) / 2

    # Restricted Jeffreys prior on mean_sigma
    mean_sigma[isi_i] <- exp(log_mean_sigma[isi_i])
    log_mean_sigma[isi_i] ~ dunif(log(1), log(max_mean_sigma[isi_i]))

    sd_sigma[isi_i] ~ dunif(0, max_sd_sigma[isi_i])
    mean_mu[isi_i] ~ dunif(min_mean_mu[isi_i], max_mean_mu[isi_i])
    sd_mu[isi_i] ~ dunif(0,max_sd_mu[isi_i])
  }
}"

# Fit the model using JAGS
# The parameters to be monitored
parameters = c( "mu" , "sigma", "mean_mu", "sd_mu", "mean_sigma", "sd_sigma")     
# Create, initialize, and adapt the model:
m = jags.model(textConnection(model_string), data=data_list, inits=inits_list, n.chains = 3 , n.adapt = 1000)
# Burn-in:
update(m, n.iter = 1000 )
# Sampling from the model
s = coda.samples(m, parameters, n.iter = 5000)

# Inspecting the model fit. 
plot(s)
summary(s)  

# Here the output from JAGS can be slightly confusing in 
# that the subject and isi levels are not retained. For example mu[3, 4] refers to the 
# estimated mean of the third subject at the fourth isi level. The third subject might very 
# well not the subject with the id 3. To see what subjects and isi levels there are
# inspect the subjects and isi_levels variables.
isi_levels
subjects

# To do further calclulations on the MCMC samples it might be convenient to combine them into one matrix
s_mat <- as.matrix(s)

# To look at the posteriors for the group level mean sigma we could do the following:
# First extract the column with the samples of the group level mean sigma
mean_sigma <- s_mat[, grep("mean_sigma", colnames(s_mat)) ]
# Put back the isi levels as column names for easier comparison
colnames(mean_sigma) <- isi_levels
# Calculate the posterior medians and a 95% credible intervall
est_sd <- apply(mean_sigma, 2, quantile, probs = c(0.975, 0.5, 0.025))
# Plot it
plot(isi_levels, est_sd["50%",], xlim=c(0, max(isi_levels)), ylim=c(0, max(est_sd)), 
     type="b", lwd=2, xlab="ISI in ms", ylab="Mean SD in ms")
lines(isi_levels, est_sd["97.5%",], type = "b", col="blue")
lines(isi_levels, est_sd["2.5%",], type = "b", col="blue")

### The following replicates the plot that is found in section 3.3 in Baath (in preparation)
# ggplot2 is required for plotting
library(ggplot2)

# Create a data frame to contain the point estimates and credible 
# intervals for the asychrony SD for the three first participants
# plus the group level mean SD.
est_df <- expand.grid(subject = 1:3, isi_level = 1:5)
est_df$isi <- isi_levels[ est_df$isi_level]
for(row_i in seq_len(nrow(est_df))) {
  subject <- est_df$subject[row_i]  
  isi_level <- est_df$isi_level[row_i]  
  asynch <- s_mat[, paste0("sigma[", subject, ",", isi_level,"]")]
  est <- quantile(asynch, c(0.025, 0.5, 0.975))
  est_df$lower[row_i] <- est[1]
  est_df$median[row_i] <- est[2]
  est_df$upper[row_i] <- est[3]
}

est_df$name <- paste("Subject", est_df$subject, "SD")

mean_sigma <- s_mat[, grep("mean_sigma", colnames(s_mat)) ]
# Put back the isi levels as column names for easier comparison
colnames(mean_sigma) <- isi_levels
# Calculate the posterior medians and a 95% credible intervall
for(isi in colnames(mean_sigma)) {
  est <-quantile(mean_sigma[, isi], c(0.025, 0.5, 0.975))
  est_df <- rbind(est_df, data.frame(subject = NA, isi_level = NA, lower = est[1], median = est[2],
                                     upper = est[3], isi = as.numeric(isi), name = "Group mean SD"))
}

# Produce the actual plot

qplot(isi, median, ymin = lower, ymax = upper, data=est_df, geom=c("line", "linerange"),
      facets = . ~ name, xlab = "Interstimulus interval in ms",
      ylab = expression("Median posterior with 95% CI")) +
  theme_minimal() + scale_x_continuous(breaks = isi_levels) + 
  coord_cartesian(ylim=c(0, 390), xlim=c(200, 3100))

# The following would save the plot as a pdf
ggsave("hierarch_asynch_sd_plot.pdf", width=8, height=2.5)


