data("ms_data")
mu <- ms_data$mu
basis <- ms_data$basis
coefs <- ms_data$coefs
time <- ms_data$time
lp <- matrix(runif(nrow(basis)),  nrow = nrow(basis))

sim_tte(pi = lp, mu = mu, basis = basis, coefs = coefs, time = time, type = "ms")
sim_tte(pi = lp, mu = -1, coefs = 1.1, time = time, type = "weibull", end = 100, obs.only = F, obs.aug = T, delta = 0.5, end = 100)
