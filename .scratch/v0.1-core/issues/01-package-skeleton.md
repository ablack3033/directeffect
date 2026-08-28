# 01: Package skeleton with green CI

**What to build:** An installable `directeffect` R package that passes `R CMD check` cleanly, locally and in CI, with test infrastructure ready for every later ticket. From the user's perspective: `library(directeffect)` works, and every push to the repo shows a green check.

**Blocked by:** None (can start immediately)

**Status:** done

- [x] `R CMD check` passes locally with no errors or warnings
- [x] A GitHub Actions workflow runs `R CMD check` on every push and passes
- [x] testthat is wired up with at least one passing test
- [x] DESCRIPTION declares title, description, authors, and license; dependencies stay minimal — added only as later tickets need them

## Comments

Done in commit f18d64d. `R CMD check` Status: OK locally (no errors,
warnings, or notes; run under C.UTF-8 — the container lacks the
en_US.UTF-8 locale). CI run #1 (R-CMD-check on ubuntu-latest, R release)
completed green: https://github.com/ablack3033/directeffect/actions/runs/33147925739
Dependencies: none in Imports at this ticket; testthat (>= 3.0.0) in
Suggests. License: Apache 2.0.
