// Stan model for Henon map reconstruction
functions {
  vector funhenon(vector thetat, row_vector xt) {
    vector[2] y;
    y[1] = thetat[1] + thetat[2] * xt[1] + thetat[3] * xt[1]^2 + xt[2];
    y[2] = thetat[4] * xt[1];
    return y;
  }
}

data {
  int < lower = 1 > N; // Sample size
  matrix[N,2] x; // map
}

parameters {
  real theta[4]; // Control parameter
  row_vector[2] x0; // Initial condition
  corr_matrix[2] Omega;// prior correlation
  vector<lower=0>[2] tau; // prior scale
  // transformed parameters {
  //   matrix[2,2] Sigma;
  //   Sigma = diag_matrix(tau) * Sigma * diag_matrix(tau);
  // }
}

model {
  // priors
  for (n in 1:4){
    theta[n] ~ uniform(-2,2);
  }  
  for (n in 1:2){
    x0[n] ~ uniform(-2,2);
  }  
  tau ~ cauchy(0, 2.5);
  Omega ~ lkj_corr(1.5);
  // likelihood
  x[1,] ~ multi_normal_prec(funhenon(to_vector(theta),x0), quad_form_diag(Omega, tau));
  for (n in 2:N) {
    x[n,] ~ multi_normal_prec(funhenon(to_vector(theta),x[n-1,]), quad_form_diag(Omega, tau));
  }
}

// generated quantities {
//   matrix[2,2] Sigma;
//   Sigma = diag_matrix(tau) * Sigma * diag_matrix(tau);
// }

