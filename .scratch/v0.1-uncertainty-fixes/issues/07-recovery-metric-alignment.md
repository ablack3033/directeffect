# 07: Recovery metrics computed over the same drugs

**What to build:** A methodologist reading `validate_recovery()` gets
all four summaries — bias, RMSE, coverage, rank correlation — over the
same set of drugs: the free ones. The pinned reference's exact (0, 0)
pair can no longer move Spearman's ρ by ±0.3 while the other metrics
exclude it.

**Blocked by:** None (can start immediately).

**Status:** done

- [x] Rank correlation is computed over the free drugs only for surface
      fits; anchored fits (no pinned drug) are unchanged
- [x] A test on a noisy simulated network asserts the free-drugs value
      and would fail under the old all-drugs computation (a case where
      the two demonstrably differ, as in the review)
- [x] Documentation states that the reference row is excluded from all
      four metrics for surface fits
- [x] The existing recovery tests (low-noise regime, anchored path)
      pass, with expectations updated only where the metric definition
      changed

## Comments

Done. `validate_recovery()` computes Spearman's ρ over `free` — the
same drugs as bias, RMSE, and coverage — so surface fits exclude the
pinned reference and anchored fits keep every drug (unchanged). Test
uses the 6-drug / 10-comparison / seed 8 simulation, where the
free-drugs ranking is perfect (ρ = 1) while the all-drugs computation
gives 0.943: the reported value must equal the free-only hand
computation and the test asserts the two definitions differ on this
case, so it fails under the old code. (The review's seed-5 case no
longer separates the definitions under the current fit path; seed 8
was measured to.) Docs state the reference row is excluded from all
four metrics for surface fits; existing low-noise and anchored
recovery tests pass unchanged.
