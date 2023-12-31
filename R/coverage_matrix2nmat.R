
#' Import genome coverage matrix files
#'
#' Import genome coverage matrix files
#'
#' This function imports genome coverage data matrix
#' and returns an object of class
#' `normalizedMatrix` compatible for use by the
#' package `"EnrichedHeatmap"`.
#'
#' There is a conversion function `EnrichedHeatmap::as.normalizedMatrix()`,
#' however this function does not call that function, in
#' favor of defining the attributes directly. In future, this
#' function may change to call that function.
#'
#' @family jam coverage heatmap functions
#' @family jam import functions
#'
#' @return `normalizedMatrix` numeric matrix, where additional
#'    metadata is stored in the object attributes. See
#'    `EnrichedHeatmap::as.normalizedMatrix()` for more
#'    details about the metadata. The `rownames` are defined
#'    by the first colname which does not match
#'    `mat_grep`, which by default is `"Gene ID"`,
#'    otherwise rownames are `NULL`.
#'
#' @param x `data.frame` or compatible object containing
#'    genome coverage data, or a character file path. When
#'    `x` is not supplied, `filename` is used to import
#'    data. When `x` is a filename, it is used to populate
#'    `filename`, then data is imported into `x`.
#' @param filename character path to a genome coverage file.
#'    When `x` is supplied, this argument is ignored. When
#'    `filename` is used, only the first file is imported.
#' @param signal_name The name of signal regions. It is only used
#'    for printing the object. When `signal_name` is `NULL`, the
#'    `signal_name` is derived from `names(filename)` if
#'    available, then `basename(filename)`, or `"signal"` then
#'    only `x` is supplied.
#' @param target_name The name of the target names. It is only used
#'    for printing the object.
#' @param background numeric value containing the background
#'    value in the matrix.
#' @param smooth logical whether to apply smoothing on rows.
#' @param target_is_single_point,signal_is_categorical logical
#'    indicating whether the target region is a single point,
#'    and whether signal matrix is categorical, respectively.
#' @param mat_grep character regular expression pattern used
#'    to identify colnames which contain coverage data. The
#'    default pattern expects the format `"-200:-100"`.
#' @param upstream_grep character regular expression pattern
#'    used to identify upstream colnames from values that
#'    match `mat_grep`. The default assumes any region
#'    beginning `"-"` is negative and upstream the central
#'    target region.
#' @param downstream_grep character regular expression pattern
#'    used to identify upstream colnames from values that
#'    match `mat_grep`. The default assumes all colnames which
#'    are not upstream are therefore downstream.
#' @param target_grep character regular expression pattern
#'    used to identify a colname referring to the `target`,
#'    which by default can only be `"0"`. Otherwise, no target
#'    region is defined.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @import data.table
#'
#' @examples
#' ## There are two small example files to use for testing
#' cov_file <- system.file("data", "tss_coverage.matrix",
#'    package="coverjam");
#' cov_file <- system.file("data", "h3k4me1_coverage.matrix",
#'    package="coverjam");
#' nmat <- coverage_matrix2nmat(cov_file);
#' jamba::printDebug("signal_name: ",
#'    attr(nmat, "signal_name"));
#'
#' #if (suppressPackageStartupMessages(require(EnrichedHeatmap))) {
#' color <- "red3";
#' colorramp <- colorjam::make_jam_divergent(color,
#'    color_complement(color));
#' signal_name <- attr(nmat, "signal_name");
#' k <- 6;
#' set.seed(123);
#' partition <- kmeans(log10(1+nmat), centers=k)$cluster;
#' EH <- EnrichedHeatmap::EnrichedHeatmap(
#'    log10(1+nmat),
#'    split=partition,
#'    pos_line=FALSE,
#'    use_raster=TRUE,
#'    col=color,
#'    #col=jamba::getColorRamp(color, n=10),
#'    top_annotation=ComplexHeatmap::HeatmapAnnotation(
#'       lines=EnrichedHeatmap::anno_enriched(
#'       gp=grid::gpar(col=colorjam::rainbowJam(k)))
#'    ),
#'    axis_name_gp=grid::gpar(fontsize=8),
#'    name=signal_name,
#'    column_title=signal_name);
#'
#'    PHM <- Heatmap(partition,
#'       use_raster=TRUE,
#'       col=structure(colorjam::rainbowJam(k),
#'          names=as.character(seq_len(k))),
#'       name="partition",
#'       show_row_names=FALSE,
#'       width=grid::unit(3, "mm"));
#'    draw(PHM + EH, main_heatmap=2);
#' }
#' }
#'
#' @export
coverage_matrix2nmat <- function
(x=NULL,
 filename=NULL,
 signal_name=NULL,
 target_name="target",
 background=0,
 smooth=FALSE,
 target_is_single_point=FALSE,
 signal_is_categorical=FALSE,
 mat_grep="[-0-9]+:[-0-9]+",
 upstream_grep="^[-]",
 downstream_grep="^[^-]",
 target_grep="^0$",
 verbose=FALSE,
 ...)
{
   ## read a matrix file, and convert to format EnrichedHeatmap understands

   ## Interpret vector of filenames as input
   if (length(filename) == 0 && length(x) > 1 && is.atomic(x)) {
      if (all(file.exists(x))) {
         xseq <- jamba::nameVector(seq_along(x), names(x));
         nmatlist <- lapply(xseq, function(i){
            coverage_matrix2nmat(x=x[i],
               filename=NULL,
               signal_name=names(x)[i],
               target_name=target_name,
               background=background,
               smooth=smooth,
               target_is_single_point=target_is_single_point,
               signal_is_categorical=signal_is_categorical,
               mat_grep=mat_grep,
               upstream_grep=upstream_grep,
               downstream_grep=downstream_grep,
               target_grep=target_grep,
               verbose=verbose,
               ...);
         });
         return(nmatlist);
      } else {
         stop("When x is supplied as a vector, it should contain files accessible with file.exists(x).");
      }
   }

   ## Other single-entry input types
   if (length(x) == 0) {
      if (length(filename) == 0) {
         stop("Must suppled x as data.frame, or filename.");
      } else if (!file.exists(head(x, 1))) {
         stop("No x supplied as data.frame, and filename is not found.");
      } else {
         if (verbose) {
            jamba::printDebug("matrix2nmat(): ",
               "Importing data from filename:",
               filename);
         }
         x <- data.table::fread(head(filename, 1),
            sep="\t",
            data.table=FALSE);
         if (length(signal_name) == 0) {
            if (length(names(filename)) > 0 && nchar(head(names(filename), 1)) > 0) {
               signal_name <- head(names(filename), 1);
            } else {
               signal_name <- basename(head(filename, 1));
            }
         }
      }
   } else if (jamba::igrepHas("character", class(x)) && file.exists(head(x, 1))) {
      if (verbose) {
         jamba::printDebug("matrix2nmat(): ",
            "Importing data from filename sent as x:",
            x);
      }
      filename <- x;
      x <- data.table::fread(head(filename, 1),
         data.table=FALSE,
         sep="\t");
      if (length(signal_name) == 0) {
         if (length(names(filename)) > 0 && nchar(head(names(filename), 1)) > 0) {
            signal_name <- head(names(filename), 1);
         } else {
            signal_name <- basename(head(filename, 1));
         }
      }
   }
   if (verbose) {
      jamba::printDebug("matrix2nmat(): ",
         "signal_name:",
         signal_name);
      jamba::printDebug("matrix2nmat(): ",
         "names(filename):",
         names(filename));
   }

   if (length(signal_name) == 0) {
      signal_name <- "signal";
   }

   if (!jamba::igrepHas("data.frame|dataframe|tibble|data.table", class(x))) {
      stop("Supplied x must be class data.frame, DataFrame, data.table, or tibble.");
   }
   if (jamba::igrepHas("data.table", class(x))) {
      x <- as.data.frame(x);
   }

   ## Extract numeric matrix
   mat_colnames <- jamba::vigrep(mat_grep, colnames(x));
   if (length(mat_colnames) == 0) {
      stop("No coordinate colnames found matching pattern.");
   }
   if (verbose) {
      jamba::printDebug("matrix2nmat(): ",
         "mat_colnames:",
         mat_colnames);
   }
   mat <- as.matrix(x[,mat_colnames,drop=FALSE]);
   name_colname <- head(setdiff(colnames(x), mat_colnames), 1);
   if (length(name_colname) > 0) {
      rownames(mat) <- x[,1];
   }

   ## Extract upstream and downstream
   colnames_up <- jamba::igrep(upstream_grep, mat_colnames);
   colnames_dn <- jamba::igrep(downstream_grep, mat_colnames);
   colnames_target <- jamba::igrep(target_grep, mat_colnames);
   colnames_dn <- setdiff(colnames_dn, colnames_target);
   attr(mat, "upstream_index") <- colnames_up;
   attr(mat, "target_index") <- colnames_target;
   attr(mat, "downstream_index") <- colnames_dn;

   ## Store filename if supplied
   if (length(filename) > 0) {
      attr(mat, "filename") <- filename;
   }

   ## Extract number of bases extended
   mat_extend <- abs(range(as.numeric(unlist(strsplit(mat_colnames, ":")))));
   attr(mat, "extend") <- mat_extend;

   ## Smooth
   attr(mat, "smooth") <- smooth;

   ## Signal and target name
   attr(mat, "signal_name") <- signal_name;
   attr(mat, "target_name") <- target_name;
   attr(mat, "target_is_single_point") <- target_is_single_point;
   attr(mat, "background") <- background;
   attr(mat, "signal_is_categorical") <- signal_is_categorical;
   class(mat) <- c("normalizedMatrix", "matrix");
   mat;
}
