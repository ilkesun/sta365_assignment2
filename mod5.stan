data {
  int<lower=0> K;          
  int<lower=0> N;          
  vector [N]y;
  simplex [K]p;
}

parameters {
  ordered[K] mu;             
  vector<lower=0>[K] sigma;  
}

model { 
  mu ~ normal(0, 3);
  sigma ~ normal(0, 2);
  
  for (n in 1:N) {
    for (k in 1:K)
    target += log_sum_exp(log(p[k]) + normal_lpdf(y[n]| mu[1], sigma[1]),
                          log(p[2]) + normal_lpdf(y[n]| mu[2], sigma[2]),
                          log(p[3]) + normal_lpdf(y[n]| mu[3], sigma[3]));
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    for (k in 1:K)
      log_lik[n] = normal_lpdf(y[n]| mu[k], sigma[k]);
  }
}