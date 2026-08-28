# Generates the packaged example network: a realistic but ENTIRELY
# SIMULATED statin evidence base. Real drug names, fake numbers: the
# true effects below are inventions calibrated to the magnitude of the
# statin trial literature (major cardiovascular events, hazard ratios
# roughly 0.7-0.85 versus placebo), and every estimate is drawn from
# them. Nothing here is a real trial result.
#
# The structure mirrors how statin evidence actually looks: head-to-head
# trials among active drugs (the comparisons), and a few landmark
# placebo-controlled trials (the anchors). Fully deterministic; re-run
# this script to regenerate data/*.rda. tests/testthat/
# test-example-data.R re-runs this generation and fails if the shipped
# data drifts from it.

# The invented truth: log hazard ratio versus placebo, placebo = 0.
truth <- c(
  atorvastatin = -0.30,
  fluvastatin  = -0.17,
  lovastatin   = -0.20,
  pravastatin  = -0.22,
  rosuvastatin = -0.32,
  simvastatin  = -0.26
)

# Head-to-head comparisons: newer / higher-intensity statins trialed
# against the older standards, with a few repeat pairs, forming a
# connected network with cycles.
design <- data.frame(
  study_id   = c("rct_ator_simv_1", "rct_ator_simv_2", "rct_ator_prav_1",
                 "rct_ator_prav_2", "rct_rosu_ator_1", "rct_rosu_ator_2",
                 "rct_rosu_simv_1", "rct_simv_prav_1", "rct_prav_lova_1",
                 "rct_lova_fluv_1", "rct_simv_fluv_1", "rct_prav_fluv_1"),
  target     = c("atorvastatin", "atorvastatin", "atorvastatin",
                 "atorvastatin", "rosuvastatin", "rosuvastatin",
                 "rosuvastatin", "simvastatin", "pravastatin",
                 "lovastatin", "simvastatin", "pravastatin"),
  comparator = c("simvastatin", "simvastatin", "pravastatin",
                 "pravastatin", "atorvastatin", "atorvastatin",
                 "simvastatin", "pravastatin", "lovastatin",
                 "fluvastatin", "fluvastatin", "fluvastatin"),
  std_error  = c(0.06, 0.09, 0.05, 0.10, 0.07, 0.11,
                 0.08, 0.07, 0.08, 0.10, 0.09, 0.12)
)

# Landmark placebo-controlled trials anchor the surface absolutely.
anchor_design <- data.frame(
  study_id  = c("rct_simv_plac", "rct_prav_plac", "rct_rosu_plac"),
  drug      = c("simvastatin", "pravastatin", "rosuvastatin"),
  reference = "placebo",
  std_error = c(0.06, 0.07, 0.08)
)

tau <- 0.02 # mild between-study heterogeneity
set.seed(12)

delta <- truth[design$target] - truth[design$comparator]
example_comparisons <- data.frame(
  study_id   = design$study_id,
  target     = design$target,
  comparator = design$comparator,
  estimate   = round(unname(stats::rnorm(
    nrow(design), delta, sqrt(design$std_error^2 + tau^2)
  )), 3),
  std_error  = design$std_error
)

example_anchors <- data.frame(
  study_id  = anchor_design$study_id,
  drug      = anchor_design$drug,
  reference = anchor_design$reference,
  estimate  = round(unname(stats::rnorm(
    nrow(anchor_design), truth[anchor_design$drug],
    anchor_design$std_error
  )), 3),
  std_error = anchor_design$std_error
)

# The generating truth ships with the package so recovery on the
# example network can be verified, not assumed.
example_truth <- data.frame(
  drug  = c("placebo", names(truth)),
  theta = c(0, unname(truth))
)

print(example_comparisons)
print(example_anchors)
print(example_truth)

save(example_comparisons, file = "data/example_comparisons.rda",
     compress = "bzip2")
save(example_anchors, file = "data/example_anchors.rda",
     compress = "bzip2")
save(example_truth, file = "data/example_truth.rda",
     compress = "bzip2")
