// Unanchored direct-effect surface model.
//
// Comparative estimates alone: y_n ~ N(theta_target - theta_comparator,
// se_n^2). Only differences are identified, so theta[1] is fixed at 0 as
// an arbitrary reference; the R adapter orders drugs so the chosen
// reference drug is index 1. The prior on the free effects is
// deliberately weak so that comparison against the frequentist engine is
// meaningful. (The design's pseudocode names the index arrays `target`
// and `comparator`; `target` is a reserved word in Stan, hence the
// `_idx` suffixes.)
data {
  int<lower=1> N;                             // number of comparisons
  int<lower=2> K;                             // number of drugs
  array[N] int<lower=1, upper=K> target_idx;
  array[N] int<lower=1, upper=K> comparator_idx;
  vector[N] y;                                // log-scale estimates
  vector<lower=0>[N] se;                      // their standard errors
}
parameters {
  vector[K - 1] theta_free;
}
transformed parameters {
  vector[K] theta;
  theta[1] = 0;
  theta[2:K] = theta_free;
}
model {
  theta_free ~ normal(0, 5);
  y ~ normal(theta[target_idx] - theta[comparator_idx], se);
}
