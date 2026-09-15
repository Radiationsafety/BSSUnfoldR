## Summary

<!-- Brief description of what this PR does and why. -->

## Type of change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New unfolding method ported from Python `bssunfold`
- [ ] New response function / dose-coefficient dataset
- [ ] Enhancement of existing method or Detector class
- [ ] Documentation / vignette
- [ ] Refactor / internal cleanup
- [ ] Breaking change

## Related issues

<!-- e.g. Closes #12 -->

## Checklist

- [ ] Code follows the package style (Roxygen2 docs, 4-space indent, snake_case).
- [ ] New function has `@export` and Roxygen2 documentation with a runnable example.
- [ ] New unfolding method follows the `solve_<name>` + `unfold_<name>` + `Detector$unfold_<name>` pattern.
- [ ] `R CMD check --as-cran` passes locally with no ERRORs/WARNINGs.
- [ ] If porting from Python `bssunfold`, the corresponding module is referenced in the source file's header comment.
- [ ] Numeric data (constants, response functions) is generated via `scripts/extract_constants_to_r.py`, not hand-typed.

## Notes for reviewer

<!-- Anything you want the reviewer to focus on. -->
