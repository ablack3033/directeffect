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
  list(
    estimate = drop(V %*% t(X) %*% W %*% comparisons$estimate),
    std_error = sqrt(diag(V))
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
