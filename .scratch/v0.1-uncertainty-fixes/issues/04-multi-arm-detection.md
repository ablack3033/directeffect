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

**Status:** ready-for-agent

- [ ] Constructing a direct-effect network where any `study_id`
      contributes more than one comparison raises a warning naming the
      studies and stating what differs between the engines for such
      input
- [ ] Fitting a multi-arm network with the Stan engine fails with an
      informative error (its likelihood treats rows as independent)
      that points at the netmeta engine
- [ ] `cycle_consistency()` and `edge_residuals()` warn on multi-arm
      networks that their variances assume independent comparisons
- [ ] README, surface-fitting docs, and the spec language state the
      narrowed model-identity claim precisely
- [ ] Fixture test: a single hand-built three-arm study is accepted by
      the netmeta engine and yields a contrast SE equal to the trial's
      own contrast SE (the √(2/3)-deflated independence answer is no
      longer reachable through any engine)
- [ ] Networks with one comparison per study trigger none of the new
      warnings — the simulator-based tests and all existing examples
      stay silent
