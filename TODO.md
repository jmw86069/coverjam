# TODO 13nov2023

* Migrate fully functional `nmatlist2heatmaps()` and supporting functions.
* Debug Github pkgdown actions
* Plan shell script wrapper function similar to slicejam:

   * `setup_coverjam()` creates bash shell script to define RHOME,
   then wrappers a call to `Rscript` which calls a `.R` file, which
   then calls `.Rmd`
   * Create `run_coverjam.R` file that renders a `.Rmd` file.
   * Create `coverjam_analysis.Rmd` to create coverage heatmaps.
