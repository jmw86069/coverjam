
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coverjam

<!-- badges: start -->
<!-- badges: end -->

The goal of coverjam is to provide flexible methods to create genome
sequence coverage heatmaps. These heatmaps typically represent sequence
coverage (depth):

- around a central genome feature, such as transcript start site (TSS)
  or ChIP-seq/ATAC-seq peaks, or
- across a scaled genome region, such as the gene locus “body” defined
  from the transcript start site (TSS) to transcript termination site
  (TTS).

This package is mainly a wrapper around fantastic Bioconductor packages
`EnrichedHeatmap` which extends `ComplexHeatmap` to render the heatmap.

The `coverjam` package provides additional methods:

- import coverage matrix data from multiple files
- define color gradients for each heatmap
- optionally handle coverage data in a group, which applies consistent
  numerical scales to the heatmap and corresponding meta-profile plot
- optionally apply row clustering by k-means or hierarchical clustering;
  offering correlation as an additional distance metric
- optionally split rows by aforementioned clustering, and/or by optional
  row annotations
- display row annotations alongside the heatmaps, in order to present
  additional context to the data displayed, for example chromosome,
  genome region type, additional metrics such as RNA-seq fold changes,
  etc.
- optionally apply numeric transformations (log2p1, sqrt, none)
- optionally sort rows by additional annotation, or by enriched signal
  using the `EnrichedHeatmap` metric.

Finally, `coverjam` intends to provide a re-usable command-line tool:

- apply an RMarkdown workflow to a set of coverage files, driven by
- configuration file used to customize the colors, numerical scales
  dimensions
- the method should be independent of user’s R environment, instead
  using the R environment of the installing user, setting the
  appropriate R_LIBS and RHOME variables to ensure the correct R is
  used.

## Installation

You can install the development version of coverjam with:

``` r
# install.packages("remotes")
remotes::install_github("jmw86069/coverjam")
```

## Notes

The methods are being migrated out of the `platjam` package, and into
this package for more focused maintenance.
