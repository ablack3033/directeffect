// Anchored direct-effect surface model.
//
// Comparative estimates constrain differences; absolute placebo anchors
// a_m ~ N(theta_drug, a_se_m^2) position the whole surface against
// placebo = 0. No arbitrary identification constraint is needed — the
// anchors, not a convention, determine the absolute location, and their
// uncertainty propagates into every drug's absolute effect. This is a
// separate program from the unanchored surface model by design: two
// simple deep modules rather than one model full of conditionals.
data {
  int<lower=1> N;                             // number of comparisons
  int<lower=2> K;                             // number of drugs
  array[N] int<lower=1, upper=K> target_idx;
  array[N] int<lower=1, upper=K> comparator_idx;
  vector[N] y;                                // log-scale estimates
  vector<lower=0>[N] se;                      // their standard errors
  int<lower=1> M;                             // number of anchors
  array[M] int<lower=1, upper=K> anchor_idx;
  vector[M] a;                                // anchor estimates vs placebo
  vector<lower=0>[M] a_se;                    // their standard errors
}
parameters {
  vector[K] theta;
}
model {
  theta ~ normal(0, 5);
  y ~ normal(theta[target_idx] - theta[comparator_idx], se);
  a ~ normal(theta[anchor_idx], a_se);
}
