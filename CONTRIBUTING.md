# Contributing to fastcramer

Bug reports, feature requests and questions are welcome as GitHub issues.
Please include a minimal reproducible example and `sessionInfo()`.

Pull requests: fork, create a branch, add or update tests under
`tests/testthat/`, make sure `R CMD check` passes with no errors or warnings,
and open the PR against `main`. New statistics or kernels must be checked
against an independent reference (for the existing kernels that reference is
`cramer::cramer.test(..., just.statistic = TRUE)`).

Please note that this project is released with a Contributor Code of Conduct
(`CODE_OF_CONDUCT.md`). By participating you agree to abide by its terms.
