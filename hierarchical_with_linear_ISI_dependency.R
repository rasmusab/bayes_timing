# This is a script implementing the model described in Baath (in preparation). The model
# below is the one descibed in the supplementary material which assumes that the group mean
# and group SD of the asynchrony depends linearly on the ISI.

# JAGS have a specific way of specifying censoring and a good guide is:
# http://doingbayesiandataanalysis.blogspot.se/2012/01/complete-example-of-right-censoring-in.html

# This script requres that you have JAGS and rjags installed, for information
# regarding JAGS see http://mcmc-jags.sourceforge.net/
require(rjags)

# The right censoring limit in ms
cens =  100

# Reading in the data and picking out the 12 first participants
d <- read.csv("tapping_data_baath_madison_2012.csv")
d <- d[ d$subject %in% 1:12,] 

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
  max_sd_sigma = isi_levels / 4,
  max_sd_mu = isi_levels / 2,
  
  asynch = d$asynch_cens,
  subject = match(d$subject, subjects),
  n_subject = length(subjects),
  cens = cens,
  is_cens =d$is_cens,
  isi_level = match(d$isi, isi_levels),
  isi_in_ms = isi_levels,
  n_isi = length(isi_levels)
)

# Specifying the initial values for the MCMC chains using estimators

inits_list = list(
  mu = tapply(d$asynch, list(d$subject, d$isi), mean),
  sigma = tapply(d$asynch, list(d$subject, d$isi), sd)
)

inits_list$sd_sigma <- apply(inits_list$sigma, 2, sd)
inits_list$sd_mu <- apply(inits_list$mu, 2, sd)

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
  
  # Specifying the hyper-priors
  for(isi_i in 1:n_isi) {
    for(subject_i in 1:n_subject) {
      mu[subject_i, isi_i] ~ dnorm(mean_mu[isi_i], 1/pow(sd_mu[isi_i], 2))
      tau[subject_i, isi_i] <- 1 / pow(sigma[subject_i, isi_i], 2)
      sigma[subject_i, isi_i] ~ dlnorm(logmu_sigma[isi_i], 1/pow(logsigma_sigma[isi_i], 2))
    }
    
    
    # Restricted Jeffreys prior on mean_sigma
    mean_sigma[isi_i] <- intercept_mean_sigma + (isi_in_ms[isi_i] - min(isi_in_ms)) * isi_mean_sigma
    sd_sigma[isi_i] ~ dunif(0, max_sd_sigma[isi_i])
    # Transformation for the log-normal on sigma
    logsigma_sigma[isi_i] <- sqrt(log( pow((sd_sigma[isi_i]/mean_sigma[isi_i]), 2) + 1 ))
    logmu_sigma[isi_i] <- log(mean_sigma[isi_i]) - pow(logsigma_sigma[isi_i], 2) / 2
    sigma_dist[isi_i] ~ dlnorm(logmu_sigma[isi_i], 1/pow(logsigma_sigma[isi_i], 2))

    mean_mu[isi_i] <- intercept_mean_mu + (isi_in_ms[isi_i] - min(isi_in_ms)) * isi_mean_mu
    sd_mu[isi_i] ~ dunif(0,max_sd_mu[isi_i])
    mu_dist[isi_i] ~ dnorm(mean_mu[isi_i], 1/pow(sd_mu[isi_i], 2))
  }


  isi_mean_sigma ~ dnorm(0, 0.001) T(0,)
  intercept_mean_sigma ~ dnorm(0, 0.001) T(0,)
  isi_mean_mu ~ dnorm(0, 0.001)
  intercept_mean_mu ~ dnorm(0, 0.001)
}"

# Fit the model using JAGS
# The parameters to be monitored
parameters = c( "mu" , "sigma", 
                "mean_mu", "sd_mu",
                "mean_sigma", "sd_sigma",
                "mu_dist", "sigma_dist",
                "intercept_mean_mu", "isi_mean_mu", "isi_mean_sigma", "intercept_mean_sigma")     
# Create, initialize, and adapt the model:
m = jags.model(textConnection(model_string), data=data_list, inits=inits_list, n.chains = 2 , n.adapt = 1000)
# Burn-in:
update(m, n.iter = 1000 )
# Sampling from the model
s = coda.samples(m, parameters, n.iter = 4000)

# Inspecting the model fit. 
plot(s)


# Plotting point estimates for the mean asynchrony, both the group estimate (green)
# and the subject estimates (grey).
ssum <- summary(s)  
est <- ssum$quantiles

plot(isi_levels, est[grep(pattern = "^mu_dist", rownames(est)), "50%"], type="b", ylim=c(-200, 100), col="darkgreen",
     ylab="Mean asynchrony in ms", xlab="ISI in ms")
for(i in seq_along(isi_levels)) {
  lines( rep(isi_levels[i], 2), est[paste0("mu_dist[", i, "]"), c("25%", "75%")], lwd=4, col="darkgreen")
  lines( rep(isi_levels[i], 2), est[paste0("mu_dist[", i, "]"), c("2.5%", "97.5%")], lwd=1, col="darkgreen")
}
for(j in seq_along(subjects)) {
  points(isi_levels, est[grep(paste0("^mu\\[", j, ","), rownames(est)), "50%"], type="b", col="grey")
}

# Plotting point estimates for the SD asynchrony, both the group estimate (red)
# and the subject estimates (grey).
plot(isi_levels, est[grep(pattern = "^sigma_dist", rownames(est)), "50%"], type="b", ylim=c(0, 400), col="darkred", lwd=2,
     ylab="Asynchrony SD in ms", xlab="ISI in ms")
for(i in seq_along(isi_levels)) {
  lines( rep(isi_levels[i], 2), est[paste0("sigma_dist[", i, "]"), c("25%", "75%")], lwd=4, col="darkred")
  lines( rep(isi_levels[i], 2), est[paste0("sigma_dist[", i, "]"), c("2.5%", "97.5%")], lwd=1, col="darkred")
}
for(j in seq_along(subjects)) {
  points(isi_levels, est[grep(paste0("^sigma\\[", j, ","), rownames(est)), "50%"], type="b", col="grey")
}