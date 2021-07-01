data {
  matrix[4, 4] X;
  int<lower=1> n[4];
  int<lower=0> y[4];
  vector[4] mu;
  matrix[4, 4] Sigma;
}
parameters {
  vector[4] beta;
}
model {
  beta ~ multi_normal(mu, Sigma);
  for(i in 1:4)
    y[i] ~ binomial_logit(n[i], X[i] * beta);
}
