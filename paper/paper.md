---
title: 'fastcramer: Fast exact computation of the Cramér multivariate two-sample statistic in R'
tags:
  - R
  - statistics
  - two-sample test
  - permutation test
  - embeddings
  - energy distance
authors:
  - name: Alexander C. Furnas
    orcid: 0000-0000-0000-0000
    affiliation: "1, 2"
affiliations:
  - name: Kellogg School of Management, Northwestern University, USA
    index: 1
  - name: Center for Science of Science and Innovation, Northwestern University, USA
    index: 2
date: 2 October 2026
bibliography: paper.bib
---

# Summary

The Cramér test [@baringhaus2004] is a nonparametric test of whether two
multivariate samples come from the same distribution. Its statistic is a
weighted contrast of mean pairwise distances between and within the two
samples, which makes it natural for high-dimensional data such as text or
image embeddings, where the question "do these two groups of documents occupy
different regions of the space?" has no parametric form. The `cramer` package
[@cramerpkg] has long provided the test in R, but its implementation builds
and loops over an $(m+n)\times(m+n)$ lookup matrix in interpreted R, so it
costs about four seconds for two samples of 500 observations in 768
dimensions and scales quadratically from there. In practice users subsample
heavily, which adds resampling noise to every inference.

`fastcramer` computes the identical statistic from its three mean pairwise
distances using BLAS matrix products and a single fused C++ pass that
evaluates the squared distance, the kernel and the running sum without
materialising large temporaries. Results agree with
`cramer::cramer.test(..., just.statistic = TRUE)` to floating-point rounding
for all five kernels the original package offers, at roughly a thousandth of
the cost: 0.3 seconds for two samples of 5,000 observations in 768
dimensions on a laptop, where the original would need several minutes. The
package also provides a label-permutation test with optional strata and a
bootstrap-versus-blocked-null design for samples too large to test in full.

# Statement of need

Dense vector representations of documents, molecules, images and genomes have
made distribution-level comparisons of two groups a routine question in the
social and natural sciences: do the documents produced by two groups, the
items funded by two agencies, or the products of two firms differ in content
once measured categories are held fixed? The Cramér statistic is well suited to
such questions because it is distribution-free, consistent against all
alternatives, and interpretable as an energy distance [@szekely2013], but
the available implementation forces a choice between tiny subsamples and
infeasible run times. In our own work on partisan differences in the content
of U.S. science [@furnas2026], the original implementation limited each
bootstrap replicate to 500 papers per group; the resampling noise of such
draws was large enough that intervals on the observed-minus-null difference
crossed zero in most fields even though the exceedance probabilities were
high. Re-running the same design with 5,000 papers per draw in `fastcramer`
removed that ambiguity in every field, and the full analysis of 26 fields,
five permutation nulls and 1,000 replicates runs in hours rather than months.

Three features address what users of the original package need:

1. **Exactness.** `cramer_statistic()` reproduces the original statistic for
   the `cramer`, `bahr`, `log`, `fraca` and `fracb` kernels; the test suite
   checks every kernel against `cramer.test()`.
2. **Permutation inference with strata.** `cramer_test()` reassigns group
   labels at random, optionally within user-supplied strata, so that nuisance
   structure such as time period or topic is preserved under the null. The
   within-stratum shuffle is vectorised and costs a fraction of a second for
   millions of observations.
3. **A design for very large samples.** `cramer_bootstrap()` implements
   repeated draws of $k$ observations per group, compared with draws made
   after permuting labels within strata, and summarises the comparison by the
   exceedance probability and the percentile interval of the pairwise
   differences, avoiding significance tests whose p-values scale with the
   number of replicates rather than with the data.

The implementation depends only on `Rcpp` and the BLAS shipped with R; a
tuned BLAS such as OpenBLAS or Apple's vecLib speeds the pure-R reference but
the fused kernel is memory-bound and performs similarly on either.

# Example

Two groups of 768-dimensional embeddings, with a nuisance category structure
(`block`), are compared on the full samples, then with labels permuted within
category so that category composition is held fixed under the null:

```r
library(fastcramer)
cramer_statistic(x, y)                                    # exact statistic, 0.3 s at 5,000 x 5,000
cramer_test(x, y, replicates = 999, block = block)        # permutation p-value, labels shuffled within category
boot <- cramer_bootstrap(x, y, replicates = 1000, draw = 5000, block = block)
boot$summary$exceedance                                   # share of observed > permuted draws
```

The package vignette (`vignette("fastcramer")`) walks through this design on
simulated embeddings, shows the agreement with `cramer::cramer.test()` for
every kernel, and discusses the interpretation of the exceedance probability
and the percentile interval.

# Acknowledgements

Development was supported by the Air Force Office of Scientific Research
under award FA9550-19-1-0354.

# References
