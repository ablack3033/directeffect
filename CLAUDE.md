# directeffect — working notes for Claude

## Git workflow (owner preference, recorded 2026-08-28)

- This is a single-user repository (ablack3033). Work directly on `main`:
  commit and push to `main`, do not create feature branches.
- Delete stale branches when you find them. (Note: `git push --delete` is
  blocked by the remote-session git proxy; ask the user to delete branches
  in the GitHub UI, or use the API if available.)

## Project layout

- v0.1 spec, design, and tickets live in `.scratch/v0.1-core/`.
- CI is `.github/workflows/R-CMD-check.yaml`; it must stay green on every
  push to `main`.
- Run `R CMD check` under the `C.UTF-8` locale in containers that lack
  `en_US.UTF-8`.
