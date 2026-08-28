# Fast demonstrations for the stats review of directeffect.
# A. anchor_surface() frequentist SEs vs exact mathematics + Monte Carlo.
# B. Bridge-edge residual is identically zero (vacuous diagnostic).
# C. validate_recovery() rank correlation includes the pinned reference.
# D. The draw-a-default-seed idiom consumes the caller's RNG stream.

# Run from the repository root.
for (f in list.files("R", full.names = TRUE)) source(f)

cat("=== Demo A: anchored frequentist SEs (mandatory example) ===\n")
comparisons <- data.frame(
  study_id = c("S1", "S2", "S3"),
  target = c("A", "A", "B"),
  comparator = c("B", "C", "C"),
  estimate = c(0.0, 0.4, 0.4),
  std_error = c(0.05, 0.05, 0.05))
anchors <- data.frame(
  study_id = "RCT1", drug = "C", reference = "placebo",
  estimate = 0.3, std_error = 0.04)
de <- direct_effect_network(comparisons, anchors = anchors,
                            effect_measure = "HR")

fit <- fit_surface(de, engine = "netmeta")
absolute <- anchor_surface(fit)
cat("Package (netmeta path) anchored effects:\n")
print(absolute$effects[, c("drug", "estimate", "std_error", "lower",
                           "upper")], digits = 4)

# Exact mathematics: the same Gaussian model, all three thetas free,
# anchor as a fourth observation row. GLS covariance is the true
# sampling covariance of the joint MLE; with a single anchor the
# package's two-stage estimator IS the joint MLE (the surface shifts by
# a - theta_hat_C, so C's absolute estimate is exactly the anchor), so
# these are the correct SEs for the package's own point estimates.
X <- rbind(
  c(1, -1, 0),   # y_AB = theta_A - theta_B
  c(1, 0, -1),   # y_AC = theta_A - theta_C
  c(0, 1, -1),   # y_BC = theta_B - theta_C
  c(0, 0, 1))    # a_C  = theta_C
se <- c(0.05, 0.05, 0.05, 0.04)
y <- c(0.0, 0.4, 0.4, 0.3)
XtW <- t(X) %*% diag(1 / se^2)
V <- solve(XtW %*% X)
cat("\nExact GLS (correct) absolute effects:\n")
print(data.frame(drug = c("A", "B", "C"),
                 estimate = drop(V %*% XtW %*% y),
                 std_error = sqrt(diag(V))), digits = 4)

# Monte Carlo under the truth the test suite itself uses:
# theta_A = 0.7, theta_B = 0.7, theta_C = 0.3. Redraw data, rerun the
# package's own two-stage pipeline (WLS surface identical to netmeta's
# common-effect estimates -- verified below -- then
# anchor_surface_offset), and measure the actual sampling SD of the
# reported absolute estimates and the actual coverage of the reported
# 95% intervals.
theta_true <- c(A = 0.7, B = 0.7, C = 0.3)

wls_surface_fit <- function(y_ab, y_ac, y_bc) {
  Xs <- rbind(c(-1, 0), c(0, -1), c(1, -1))       # free: B, C (ref A)
  Ws <- diag(1 / c(0.05, 0.05, 0.05)^2)
  Vs <- solve(t(Xs) %*% Ws %*% Xs)
  est <- drop(Vs %*% t(Xs) %*% Ws %*% c(y_ab, y_ac, y_bc))
  effects <- data.frame(
    drug = c("A", "B", "C"),
    estimate = c(0, est), std_error = c(0, sqrt(diag(Vs))),
    lower = NA_real_, upper = NA_real_, scale = "log",
    reference = "A", engine = "netmeta")
  new_directeffect_fit(effects, list(), "netmeta", NULL, de)
}

# Check the hand WLS reproduces the netmeta surface before relying on it.
hand <- wls_surface_fit(0.0, 0.4, 0.4)
stopifnot(max(abs(hand$effects$estimate - fit$effects$estimate)) < 1e-10,
          max(abs(hand$effects$std_error - fit$effects$std_error)) < 1e-10)

set.seed(1)
reps <- 20000
draws <- matrix(NA_real_, reps, 3, dimnames = list(NULL, c("A", "B", "C")))
covered <- matrix(NA, reps, 3, dimnames = list(NULL, c("A", "B", "C")))
for (r in seq_len(reps)) {
  y_ab <- rnorm(1, 0.0, 0.05)
  y_ac <- rnorm(1, 0.4, 0.05)
  y_bc <- rnorm(1, 0.4, 0.05)
  a <- rnorm(1, 0.3, 0.04)
  anch <- anchors; anch$estimate <- a
  res <- anchor_surface_offset(wls_surface_fit(y_ab, y_ac, y_bc), anch)
  est <- res$effects$estimate
  se_rep <- res$effects$std_error
  draws[r, ] <- est
  covered[r, ] <- abs(est - theta_true) <= 1.96 * se_rep
}
cat("\nMonte Carlo (20k reps), truth A=0.7 B=0.7 C=0.3:\n")
print(data.frame(
  drug = c("A", "B", "C"),
  reported_se = absolute$effects$std_error,
  true_sampling_sd = apply(draws, 2, sd),
  exact_gls_se = sqrt(diag(V)),
  empirical_coverage_of_nominal_95 = colMeans(covered)), digits = 4)

cat("\n=== Demo B: bridge-edge residual is identically zero ===\n")
bridge <- function(cd_estimate) {
  cmp <- data.frame(
    study_id = c("S1", "S2", "S3", "S4"),
    target = c("A", "A", "B", "C"),
    comparator = c("B", "C", "C", "D"),
    estimate = c(0.1, 0.3, 0.25, cd_estimate),
    std_error = c(0.05, 0.05, 0.05, 0.05))
  de_b <- direct_effect_network(cmp, effect_measure = "HR")
  edge_residuals(fit_surface(de_b, engine = "netmeta"))
}
r1 <- bridge(0.2)
r2 <- bridge(5.0)
cat("standardized residual of the C-D bridge, estimate = 0.2:",
    r1$standardized_residual[4], "\n")
cat("standardized residual of the C-D bridge, estimate = 5.0:",
    r2$standardized_residual[4], "\n")

cat("\n=== Demo C: rank correlation includes the pinned reference ===\n")
sim <- simulate_direct_effect_network(n_drugs = 6, n_comparisons = 10,
                                      n_anchors = 0, heterogeneity = 0,
                                      seed = 5)
fit_s <- fit_surface(sim$network, engine = "netmeta")
rec <- validate_recovery(fit_s, sim)
truth <- sim$truth
theta_t <- truth$theta[match(fit_s$effects$drug, truth$drug)]
theta_t <- theta_t - theta_t[fit_s$effects$drug == "drug_01"]
free <- fit_s$effects$drug != "drug_01"
cat("reported rank_correlation (includes pinned (0,0) row):",
    rec$rank_correlation, "\n")
cat("rank correlation over the", sum(free), "free drugs only:      ",
    cor(fit_s$effects$estimate[free], theta_t[free],
        method = "spearman"), "\n")

cat("\n=== Demo D: draw-a-default-seed idiom consumes caller RNG ===\n")
set.seed(7); expected <- rnorm(1)
set.seed(7)
invisible(sample.int(.Machine$integer.max, 1))  # fit_surface_stan default
after <- rnorm(1)
cat("rnorm(1) after set.seed(7):            ", expected, "\n")
cat("rnorm(1) after set.seed(7) + the idiom:", after, "\n")
