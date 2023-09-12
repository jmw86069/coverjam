
## 07sep2023

Plan command-line interface.
Two options.

1. RMarkdown workflow, with shell script wrapper.

   * Shell script needs to be "installed" by the user who installed
   `coverjam`. It detects `RHOME` and `R_LIBS` so the script will
   call `Rscript` with the correct `R` executable, that user's R library
   files, and therefore is "portable" in that the script can be called
   by anybody else without having to install R and R package dependencies.
   * The shell script calls `Rscript` to render the `.Rmd` file.
   * User needs to supply config file, or command-line arguments in some way.

2. The tool [viash.io](https://viash.io) can create command-line executable.

   * It requires Java.
   * See funkyheatmap for example of deploying executables for re-use.
   * optionally creates Docker and Nextflow files.

