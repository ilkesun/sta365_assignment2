data {
  int<lower=1> K;          
  int<lower=1> N;
  real wind[N];
  real temp[N];
  real y[N];
}

parameters {
  ordered[K] mu;             
  vector<lower=0>[K] sigma;
  real b0;
  real b1;
  real b2;
}

transformed parameters {
  real p;
  for (n in 1:N)
    p = inv_logit(b0 + b1 * wind[n] + b2 * temp[n]);
}

model {  
  mu ~ normal(0, 8);
  sigma ~ normal(0, 3);
  b0 ~ normal(0,5);
  b1 ~ normal(0,5);
  b2 ~ normal(0,5);
  
  for (n in 1:N) {
    target += log_sum_exp(log(p) + normal_lpdf(y[n]| mu[1], sigma[1]),
              log1m(p) + normal_lpdf(y[n]| mu[2], sigma[2]));
  }
}
