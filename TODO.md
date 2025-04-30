# TODO 23apr2024

* Complete the full migration so that `platjam` calls functions in this package

* For profile plots, the error bars sometimes exceed the y-axis range.

   * Either expand the y-axis range to accomodate the region required, or
   * Clip the error bar so it is not shown beyond the plot panel boundary.


# TODO 13nov2023

* Migrate fully functional `nmatlist2heatmaps()` and supporting functions.
* Debug Github pkgdown actions
* Plan shell script wrapper function similar to slicejam:

   * `setup_coverjam()` creates bash shell script to define RHOME,
   then wrappers a call to `Rscript` which calls a `.R` file, which
   then calls `.Rmd`
   * Create `run_coverjam.R` file that renders a `.Rmd` file.
   * Create `coverjam_analysis.Rmd` to create coverage heatmaps.
