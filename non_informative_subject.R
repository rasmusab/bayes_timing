# This is a script implementing the model described in Baath (in preparation). The model
# below is the one refered to as the Non-hierarchical model in the manuscript.

# JAGS have a specific way of specifying censoring and a good guide is:
# http://doingbayesiandataanalysis.blogspot.se/2012/01/complete-example-of-right-censoring-in.html

# This script requres that you have JAGS and rjags installed, for information
# regarding JAGS see http://mcmc-jags.sourceforge.net/
require(rjags)

# The right censoring limit in ms
cens =  100

# Reading in the data, and extract just the 3 first participants
d <- read.csv("tapping_data_baath_madison_2012.csv")
d <- d[ d$subject %in% 1:3,] 

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
  asynch = d$asynch_cens,
  subject = match(d$subject, subjects),
  n_subject = length(subjects),
  sigma_high = isi_levels / 2, # This is the assumed maximum possible SD
  mu_low = -isi_levels / 2,
  mu_high = isi_levels / 2,
  cens = cens,
  is_cens = d$is_cens,
  isi = match(d$isi, isi_levels),
  n_isi = length(isi_levels)
)

# Specifying the initial values for the MCMC chains using robust estimators

inits_list = list(
  mu = tapply(d$asynch, list(d$subject, d$isi), mean, trim=0.2),
  log.sigma = log(tapply(d$asynch, list(d$subject, d$isi), mad))
)

if(any(d$is_cens)) {
  # Initialize the censored asynchronies to be above the censoring threshold
  inits_list <- c(inits_list, list(asynch = ifelse(d$is_cens, cens + 1, NA) ))
}

# The model specified in the JAGS language
model_string <- "model {
  for ( i in 1:length(asynch) ) {
    asynch[i] ~ dnorm(mu[subject[i], isi[i]], tau[subject[i], isi[i]])
    is_cens[i] ~ dinterval(asynch[i], cens)
  }

  # Specifying the priors
  for(isi_i in 1:n_isi) {
    for(subject_i in 1:n_subject) {
      mu[subject_i, isi_i] ~ dunif(mu_low[isi_i], mu_high[isi_i])
      
      # A restricted Jeffrey's prior on sigma
      sigma[subject_i, isi_i] <- exp(log.sigma[subject_i, isi_i])
      log.sigma[subject_i, isi_i] ~ dunif(log(1), log(sigma_high[isi_i]))
      # transforming to precision in order to adhere to the parameterization 
      # of the Normal distribution in JAGS
      tau[subject_i, isi_i] <- 1/pow( sigma[subject_i, isi_i] , 2 )
    }
  }
}"

# Fit the model using JAGS
parameters = c( "mu" , "sigma")     # The parameters to be monitored
# Create, initialize, and adapt the model:
m = jags.model(textConnection(model_string), data=data_list , inits=inits_list,
                        n.chains = 3 , n.adapt = 1000)
# Burn-in:
update( m , n.iter = 1000  )
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
