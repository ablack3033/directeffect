# 01: Package skeleton with green CI

**What to build:** An installable `directeffect` R package that passes `R CMD check` cleanly, locally and in CI, with test infrastructure ready for every later ticket. From the user's perspective: `library(directeffect)` works, and every push to the repo shows a green check.

**Blocked by:** None (can start immediately)

**Status:** claimed

- [ ] `R CMD check` passes locally with no errors or warnings
- [ ] A GitHub Actions workflow runs `R CMD check` on every push and passes
- [ ] testthat is wired up with at least one passing test
- [ ] DESCRIPTION declares title, description, authors, and license; dependencies stay minimal — added only as later tickets need them
