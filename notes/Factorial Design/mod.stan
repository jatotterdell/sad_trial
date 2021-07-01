data {
  matrix[4, 4] X;
  vector[4] y;
  vector[4] mu;
  matrix[4, 4] Sigma;
  real<lower=0> sigma;
}
parameters {
  vector[4] beta;
}
model {
  beta ~ multi_normal(mu, Sigma);
  y ~ normal(X*beta, sigma);
}
