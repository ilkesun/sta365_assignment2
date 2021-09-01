data {
  int<lower=0> K;          
  int<lower=0> N;          
  vector [N]y;
  real<lower=0, upper=1>p;
}

parameters {
  ordered[K] mu;             
  vector<lower=0>[K] sigma;  
}

model { 
  mu ~ normal(0, 8);
  sigma ~ normal(0, 3);
  
  for (n in 1:N) {
    target += log_mix(p,
                      normal_lpdf(y[n]| mu[1], sigma[1]),
                      normal_lpdf(y[n]| mu[2], sigma[2]));
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    for (k in 1:K)
      log_lik[n] = normal_lpdf(y[n]| mu[k], sigma[k]);
  }
}