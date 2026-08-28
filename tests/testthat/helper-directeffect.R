# Shared across test files (testthat loads helper-*.R before all tests).

# The spec's 3-drug example: A, B, C compared pairwise, C anchored to
# placebo.
spec_comparisons <- function() {
  data.frame(
    study_id   = c("S1", "S2", "S3"),
    target     = c("A", "A", "B"),
    comparator = c("B", "C", "C"),
    estimate   = c(log(1.02), log(1.34), log(1.29)),
    std_error  = c(0.07, 0.09, 0.08)
  )
}

spec_anchors <- function() {
  data.frame(
    study_id  = "RCT1",
    drug      = "C",
    reference = "placebo",
    estimate  = log(1.20),
    std_error = 0.04
  )
}

# Hand oracle: the common-effect surface is weighted least squares of the
# comparison estimates on drug-difference contrasts with the reference
# drug fixed at 0. Computed directly from the normal equations so engine
# output can be checked against arithmetic, not against another engine.
# `covariance` is the full drug-by-drug covariance in `treatments`
# order, with the reference's row and column exactly 0.
wls_surface <- function(comparisons, treatments, reference) {
  free <- setdiff(treatments, reference)
  X <- matrix(0, nrow(comparisons), length(free),
              dimnames = list(NULL, free))
  for (k in seq_len(nrow(comparisons))) {
    target <- comparisons$target[k]
    comparator <- comparisons$comparator[k]
    if (target != reference) X[k, target] <- 1
    if (comparator != reference) X[k, comparator] <- -1
  }
  W <- diag(1 / comparisons$std_error^2, nrow = nrow(comparisons))
  V <- solve(t(X) %*% W %*% X)
  covariance <- matrix(0, length(treatments), length(treatments),
                       dimnames = list(treatments, treatments))
  covariance[free, free] <- V
  list(
    estimate = drop(V %*% t(X) %*% W %*% comparisons$estimate),
    std_error = sqrt(diag(V)),
    covariance = covariance,
    design = X,
    weights = W
  )
}

# Independent joint-GLS oracle for anchored fits: append the anchor rows
# to the weighted comparison design and solve the normal equations for
# the absolute effects directly. All drugs are free (the anchors, not a
# pinned reference, identify the absolute location), so this is the
# one-stage answer the two-stage surface-then-offset pipeline is
# checked against.
joint_gls <- function(comparisons, anchors, treatments) {
  n <- nrow(comparisons) + nrow(anchors)
  X <- matrix(0, n, length(treatments),
              dimnames = list(NULL, treatments))
  for (k in seq_len(nrow(comparisons))) {
    X[k, comparisons$target[k]] <- 1
    X[k, comparisons$comparator[k]] <- -1
  }
  for (m in seq_len(nrow(anchors))) {
    X[nrow(comparisons) + m, anchors$drug[m]] <- 1
  }
  y <- c(comparisons$estimate, anchors$estimate)
  W <- diag(1 / c(comparisons$std_error, anchors$std_error)^2)
  V <- solve(t(X) %*% W %*% X)
  dimnames(V) <- list(treatments, treatments)
  list(
    estimate = drop(V %*% t(X) %*% W %*% y),
    std_error = sqrt(diag(V)),
    covariance = V
  )
}

# An asymmetric four-drug network — unequal standard errors and a
# repeated pair — so covariance off-diagonals and per-row leverages are
# all distinct and informative for the matrix oracles.
asymmetric_comparisons <- function() {
  data.frame(
    study_id   = c("S1", "S2", "S3", "S4", "S5"),
    target     = c("A", "A", "B", "C", "B"),
    comparator = c("B", "C", "C", "D", "D"),
    estimate   = c(0.10, 0.32, 0.19, 0.05, 0.28),
    std_error  = c(0.05, 0.08, 0.06, 0.04, 0.09)
  )
}

# A single three-arm trial (arms A, B, C): three comparison rows sharing
# one study_id, with internally consistent contrasts (AB + BC = AC) and
# equal contrast standard errors (arm variance 0.005 each). The input
# class the simulator cannot generate.
three_arm_comparisons <- function() {
  data.frame(
    study_id   = c("T1", "T1", "T1"),
    target     = c("A", "A", "B"),
    comparator = c("B", "C", "C"),
    estimate   = c(0.1, 0.3, 0.2),
    std_error  = c(0.1, 0.1, 0.1)
  )
}

# Collect every warning an expression raises, muffling them, and return
# the messages; lets a test assert on one expected warning while an
# engine emits unrelated ones (e.g. rstan's own low-ESS warnings).
collect_warnings <- function(expr) {
  messages <- character()
  withCallingHandlers(
    expr,
    warning = function(cnd) {
      messages <<- c(messages, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  messages
}
