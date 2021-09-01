data {
  int<lower=0> K;          
  int<lower=0> N;          
  real y[N];
}

parameters {
  simplex[K] p;          
  ordered[K] mu;             
  vector<lower=0>[K] sigma;  
}

model {
  vector[K] log_p = log(p);  
  sigma ~ normal(0, 3);
  mu ~ normal(0, 8);
  
  for (n in 1:N) {
    vector[K] lps = log_p;
    for (k in 1:K)
      lps[k] += normal_lpdf(y[n] | mu[k], sigma[k]);
    target += log_sum_exp(lps);
  }
}
