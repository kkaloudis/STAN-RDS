// Stan model for logistic map reconstruction
data {
 int < lower = 1 > N; // Sample size
 vector[N] x; // map
}

parameters {
 real theta; // Control parameter
 real x0; // Initial condition
 real < lower = 0 > sigma; // Error SD
}

model {
  // priors
 theta ~ uniform(0,2);
 x0 ~ uniform(-1,1);
 sigma ~ gamma(1e-03,1e-03);
 // likelihood
 x[1] ~ normal(1 - theta * x0^2, sigma);
 for (n in 2:N) {
 x[n] ~ normal(1 - theta * x[n-1]^2, sigma);
 }
}

