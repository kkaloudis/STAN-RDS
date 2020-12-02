// Stan model for logistic map reconstruction & prediction
data {
  int < lower = 1 > N; // Sample size
  vector[N] x; // map
  int <lower=0> N_new; // prediction horizon
}

parameters {
  real theta; // Control parameter
  real < lower = -1, upper = 1 > x0; // Initial condition
  real < lower = 0 > sigma; // Error SD
  vector < lower = -1, upper = 1 >[N_new] x_pred; // predictions
}

model {
  // priors
  theta ~ uniform(0,2);
  x0 ~ uniform(-1,1);
  sigma ~ gamma(1e-03,1e-03);
  for (n in 1:N_new) {
    x_pred[n] ~ uniform(-1, 1);
  }  
  // likelihood
  x[1] ~ normal(1 - theta * x0^2, sigma);
  for (n in 2:N) {
    x[n] ~ normal(1 - theta * x[n-1]^2, sigma);
  }
  x_pred[1] ~ normal(1 - theta * x[N]^2, sigma);
  for (n in 2:N_new) {
    x_pred[n] ~ normal(1 - theta * x_pred[n-1]^2, sigma);
  }
}

