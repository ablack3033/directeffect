# Generates the packaged example network: a simulated evidence base of
# eight fictional drugs compared head-to-head, with two placebo anchors.
# Fully deterministic; re-run this script to regenerate data/*.rda.

devtools::load_all()

drug_names <- c("abrelor", "belvane", "cordexa", "delmiro",
                "envarin", "folzane", "gravix", "harlent")

simulation <- simulate_direct_effect_network(
  n_drugs = 8,
  n_comparisons = 16,
  n_anchors = 2,
  heterogeneity = 0.02,
  seed = 731,
  effect_sd = 0.35,
  se_range = c(0.05, 0.12),
  anchor_se_range = c(0.04, 0.08)
)

rename <- stats::setNames(drug_names,
                          sprintf("drug_%02d", seq_along(drug_names)))

example_comparisons <- simulation$comparisons
example_comparisons$study_id <- sprintf("study_%02d",
                                        seq_len(nrow(example_comparisons)))
example_comparisons$target <- unname(rename[example_comparisons$target])
example_comparisons$comparator <-
  unname(rename[example_comparisons$comparator])
example_comparisons$estimate <- round(example_comparisons$estimate, 3)
example_comparisons$std_error <- round(example_comparisons$std_error, 3)

example_anchors <- simulation$anchors
example_anchors$study_id <- sprintf("rct_%02d",
                                    seq_len(nrow(example_anchors)))
example_anchors$drug <- unname(rename[example_anchors$drug])
example_anchors$estimate <- round(example_anchors$estimate, 3)
example_anchors$std_error <- round(example_anchors$std_error, 3)

# The generating truth, for the record (log-HR scale, placebo = 0):
example_truth <- simulation$truth
example_truth$drug <- c("placebo", unname(rename[example_truth$drug[-1]]))
example_truth$theta <- round(example_truth$theta, 3)
print(example_truth)

save(example_comparisons, file = "data/example_comparisons.rda",
     compress = "bzip2")
save(example_anchors, file = "data/example_anchors.rda",
     compress = "bzip2")
