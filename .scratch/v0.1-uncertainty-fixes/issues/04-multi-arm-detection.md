# 04: Multi-arm evidence is detected loudly, never silently mishandled

**What to build:** A researcher whose evidence base contains a
multi-arm trial (several comparison rows sharing a `study_id`) is told
so at network construction, is refused by the Stan engine with an
explanation and a pointer to the netmeta engine (which handles the
within-study correlation correctly), and sees the same caveat from the
independence-based diagnostics. The README's and docs' model-identity
claim becomes exactly true: identical likelihood for
one-comparison-per-study networks, multi-arm supported by the
frequentist engine only in this version. Fail loudly, don't
approximate silently.

**Blocked by:** None (can start immediately).

**Status:** done

- [x] Constructing a direct-effect network where any `study_id`
      contributes more than one comparison raises a warning naming the
      studies and stating what differs between the engines for such
      input
- [x] Fitting a multi-arm network with the Stan engine fails with an
      informative error (its likelihood treats rows as independent)
      that points at the netmeta engine
- [x] `cycle_consistency()` and `edge_residuals()` warn on multi-arm
      networks that their variances assume independent comparisons
- [x] README, surface-fitting docs, and the spec language state the
      narrowed model-identity claim precisely
- [x] Fixture test: a single hand-built three-arm study is accepted by
      the netmeta engine and yields a contrast SE equal to the trial's
      own contrast SE (the √(2/3)-deflated independence answer is no
      longer reachable through any engine)
- [x] Networks with one comparison per study trigger none of the new
      warnings — the simulator-based tests and all existing examples
      stay silent

## Comments

Done. `direct_effect_network()` warns via `warn_on_multiarm_studies()`
(names the studies, states that netmeta models the within-study
correlation while the Stan engine refuses). `fit_surface_stan()` and
`anchor_surface_stan()` both call `refuse_multiarm_stan()` before
sampling: an error naming the studies, explaining the independence
likelihood, and pointing at `engine = "netmeta"`. `edge_residuals()`
and `cycle_consistency()` warn through `warn_on_multiarm_diagnostic()`.
README, `fit_surface()` docs, all three vignettes, and the v0.1-core
spec's "Two engines, one model" bullet now state the narrowed claim:
identical likelihood for one-comparison-per-study networks, multi-arm
via netmeta only in this version. Fixture (`three_arm_comparisons()`
helper, consistent contrasts, SE 0.1 each): netmeta accepts it and
returns the trial's own contrast SE 0.1 — the √(2/3)-deflated 0.0816 is
no longer reachable through any engine. Construction of every
one-comparison-per-study network in the suite (simulator, spec and
example data) stays warning-free, asserted explicitly.
