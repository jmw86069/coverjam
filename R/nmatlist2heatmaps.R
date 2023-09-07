
#' Make multiple coverage heatmaps
#'
#' Make multiple coverage heatmaps
#'
#' This function takes a list of `normalizedMatrix` objects,
#' usually the output of `coverage_matrix2nmat()`, and
#' produces multiple heatmaps using
#' `EnrichedHeatmap`.
#'
#' This function is intended to be a convenient wrapper to
#' help keep each data matrix in order, to apply consistent
#' clustering and filtering across all data matrices,
#' and to enable optional multi-row heatmap layout.
#'
#' @param nmatlist `list` containing `normalizedMatrix` objects,
#'    usually the output from `coverage_matrix2nmat()`.
#' @param k_clusters integer number of k-means clusters to
#'    use to partition each heatmap. Use `0` or `NULL` for
#'    no clustering.
#' @param k_subset integer vector of k-means clusters to
#'    retain. Often one cluster contains mostly empty
#'    values, and can be removed using this mechanism.
#' @param k_colors vector of R colors, or `NULL` to use
#'    the output of `colorjam::rainbowJam(k_clusters)`.
#' @param k_width `unit` width of the k-means cluster color
#'    bar, used with `k_clusters`.
#' @param k_method `character` string indicating the distance
#'    used by k-means, where the common default is
#'    `"euclidean"`, however a useful alternative for
#'    sequence coverage data is `"correlation"` as implemented
#'    in `amap::Kmeans()`.
#' @param k_heatmap `integer` indicating which one or more
#'    `normalizedMatrix` objects in `nmatlist` will be used
#'    for k-means clustering, when `k_clusters` is defined
#'    more than 1.
#' @param partition `character` or `factor` vector used to split rows
#'    of each matrix in `nmatlist`, named by rownames. This
#'    value is ignored when `k_clusters` is supplied.
#' @param rows `character` vector of `rownames(nmatlist)` or
#'    `integer` vector with index of rows to keep from each
#'    matrix in `nmatlist`.
#' @param row_order integer vector used to order rows.
#'    When `TRUE` or `NULL` it uses
#'    the default for `EnrichedHeatmap::EnrichedHeatmap()`
#'    which is the `EnrichedHeatmap::enriched_score()`
#'    for the matrix `main_heatmap`. When `FALSE` the
#'    rows are ordered by the order they appear in `rows`,
#'    which is either the order they appear in `nmatlist`
#'    or the order after sorting `anno_df`. When
#'    `TRUE` the default
#' @param nmat_colors named character vector of R colors,
#'    to colorize each heatmap. When `NULL` then
#'    `colorjam::rainbowJam()` is used to create colors
#'    for each heatmap panel.
#' @param middle_color `character` R compatible color used
#'    when creating a divergent color gradient, this color
#'    is used as the middle color. Usually this color should
#'    be either `"white"` or `"black"`.
#' @param nmat_names `character` vector, or `NULL`, optional,
#'    used as custom names for each heatmap in `nmatlist`.
#'    When `nmat_names=NULL` the `signal_name` values are
#'    used from each `nmatlist` matrix.
#' @param main_heatmap integer index referring to the
#'    entry in `nmatlist` to use for clustering and row
#'    ordering.
#' @param anno_df `data.frame` or object that can be coerced,
#'    used to annotate rows of each matrix. It must have
#'    `rownames(anno_df)` that match `rownames(nmatlist)`.
#'    When supplied, data can be sorted using `byCols`.
#'    Note that only the `rownames(anno_df)`
#'    present in both `nmatlist` and `anno_df` are
#'    used to display the heatmaps. These rows
#'    may also be subsetted using argument `rows`.
#' @param byCols character vector of  values in
#'    `colnames(anno_df)` used to sort the data.frame
#'    via `jamba::mixedSortDF()`. Any colname with
#'    prefix `-` will be reverse-sorted.
#' @param color_sub `character vector` of R colors to be used
#'    as categorical colors, whose names match items to be
#'    colored. This argument is intended for `anno_df`,
#'    for any column in `anno_df` where all values in that
#'    column are also in `names(color_sub)` will be colorized
#'    using `color_sub` instead of generating new colors.
#'    Also colors for partition and kmeans clusters, usually
#'    defined with `k_colors` can be defined in color_sub,
#'    if `names(color_sub)` match the partition labels.
#' @param legend_max_ncol integer number indicating the maximum
#'    number of columns allowed for a categorical color legend.
#' @param legend_base_nrow integer number indicating the base
#'    number of rows used for a categorical color legend, before
#'    additional columns are added. Once the number of elements
#'    exceeds `(legend_max_ncol * legend_base_nrow)` then
#'    rows are added, but columns never exceed `legend_max_ncol`.
#' @param legend_max_labels `integer` to define the maximum labels
#'    to display as a color legend. When any `anno_df` column contains
#'    more than this number of categorical colors, the legend is
#'    not displayed (because it would prevent display of the
#'    heatmaps at all).
#' @param anno_row_marks character vector of `rownames`
#'    which will be labeled beside the heatmaps, using
#'    the `ComplexHeatmap::anno_mark()` method. It currently
#'    requires `anno_df` be defined, since it uses the
#'    first column in `anno_df` as a one-column heatmap,
#'    to anchor the labels.
#' @param anno_row_labels character vector of optional
#'    character labels to use instead of `rownames`.
#'    If `NULL` then `anno_row_marks` are used. Or
#'    `anno_row_labels` may contain a character vector
#'    of `colnames(anno_df)` which will create labels
#'    by concatenating each column value separated by
#'    space `" "`.
#' @param top_annotation `HeatmapAnnotation` or `logical` or `list`:
#'    * `TRUE` to use the default approach
#'    `EnrichedHeatmap::anno_enriched()`
#'    * `FALSE` to prevent the display of top annotation
#'    * `HeatmapAnnotation` which should be in the form
#'    `ComplexHeatmap::HeatmapAnnotation(EnrichedHeatmap::anno_enriched())`
#'    or equivalent. This form is required for the annotation
#'    function to be called on each coverage matrix heatmap.
#'    * `list` of objects suitable to be passed as a
#'    `top_annotation` argument for each coverage heatmap,
#'    in order of `nmatlist`.
#' @param top_anno_height `unit` object to define the default
#'    height of the `top_annotation`. When `top_annotation`
#'    is not defined, the default method uses
#'    `EnrichedHeatmap::anno_enriched()` with
#'    `height=top_anno_height`.
#' @param top_axis_side `character` value indicating which side
#'    of the top annotation to place the y-axis labels.
#'    When there is one value, it is repeated to `length(nmatlist)`,
#'    otherwise it is mainly used when `panel_groups` are
#'    provided, in which case only one top annotation is
#'    label per contiguous set of panels in the same panel group.
#'    In that case `"left"` will label the left side of the
#'    first panel in each group, `"right"` will label the
#'    right side of the last panel in each group.
#'    Values:
#'    "left", "right", "both", "none", "all".
#'
#' @param hm_nrow integer number of rows used to display
#'    the heatmap panels.
#' @param transform either `character` string referring to
#'    a numeric transformation, or a `function` that applies
#'    a numeric transformation. Valid `character` string values:
#'    `"log2signed"` applies `jamba::log2signed()` which applies
#'    `log2(1+x)` transform to the absolute value, then multiplies
#'    by the original `sign(x)`; `"sqrt"` applies square root;
#'    `"cubert"` applies cube root `x^(1/3)`; `"qrt"` applies
#'    fourth root `x^(1/4)`. When there are negative numeric
#'    values, the transformation is applied to absolute value,
#'    then multiplied by the original sign. Therefore, the
#'    transformation is applied to adjust the magnitude of
#'    the values. These values are passed to `get_numeric_transform()`
#'    which may have more information.
#' @param signal_ceiling `numeric` vector whose values are recycled
#'    to length `length(nmatlist)`. The signal_ceiling
#'    applies a maximum numeric value to the
#'    color ramp for each matrix in `nmatlist`. The value is
#'    passed to `get_nmat_ceiling()`, which recognizes three
#'    numeric forms:
#'    * `signal_ceiling > 1`: this specific numeric value
#'    is applied as the ceiling
#'    * `signal_ceiling > 0` and `signal_ceiling <= 1`: this numeric
#'    value is interpreted as a quantile threshold, for example
#'    `signal_ceiling=0.75` would calculate ceiling `quantile(x, probs=0.75)`.
#'    * `signal_ceiling` is `NULL`: the maximum absolute value of each
#'    matrix is used as the ceiling.
#'
#'    Note that the ceiling is only applied to color scale and
#'    not to the underlying data, which is useful to know because any
#'    clustering and row ordering steps will use the full data
#'    as needed.
#'
#'    If data needs to be strictly controlled to a
#'    numeric ceiling, that processing should take place
#'    on `nmatlist` before calling `nmatlist2heatmaps()`.
#' @param lens numeric value used to scale each heatmap
#'    color ramp, using `getColorRamp()`. Values above zero
#'    apply the color gradient more rapidly starting from the
#'    lowest value, making the color appear more intense for
#'    lower numeric values. Values below zero apply the color gradient
#'    less rapidly, which makes lower numeric values appear
#'    less intense. This adjustment is intended to help
#'    apply suitable color contrast depending upon the range
#'    of numeric values. The `lens` values are applied to
#'    each matrix in `nmatlist`, and so it is recycled to
#'    `length(nmatlist)` as needed. Note that `signal_ceiling`
#'    is also intended to help apply the color gradient to
#'    a suitable numeric range, and the `lens` argument is
#'    applied relative to the numeric range being used.
#' @param anno_lens numeric value used to scale the annotation
#'    heatmap color scales, see `lens` for details. Values
#'    higher than 1 make the color gradient more intense,
#'    values below -1 make the color gradient less intense.
#' @param axis_name_gp x-axis label graphic parameters,
#'    as output from `grid::gpar()`. For example to define
#'    the x-axis font size, use the form
#'    `grid::gpar(fontsize=8)`.
#' @param  axis_name_rot numeric value either `0` or `90` indicating
#'    whether to rotate the x-axis names, where `90` will rotate
#'    labels, and `0` will leave labels horizontal.
#' @param column_title_gp heatmap title graphic parameters,
#'    as output from `grid::gpar()`. For example to define
#'    the x-axis font size, use the form
#'    `grid::gpar(fontsize=8)`. This argument is passed
#'    directly to `ComplexHeatmap::Heatmap()`.
#' @param seed numeric value used with `set.seed()` to
#'    set the random seed. Set to `NULL` to avoid running
#'    `set.seed()`.
#' @param ht_gap `unit` size to specify the gap between multiple heatmaps.
#'    This argument is passed to `ComplexHeatmap::draw()`. An example
#'    is `grid::unit(8, "mm")` to specify 8 millimeters.
#' @param profile_value character string to define the type of numeric
#'    profile to display at the top of each heatmap. This argument is
#'    passed to `EnrichedHeatmap::anno_enriched()`. Values: `"mean"` the
#'    mean profile; `"sum"` the sum; `"abs_sum"` sum of absolute values;
#'    `"abs_mean"` the mean of absolute values.
#' @param ylims `vector` of maximum y-axis values for each heatmap profile;
#'    or `list`
#' @param border `logical` indicating whether to draw a border around the
#'    heatmap, which includes all heatmap panels in the event of
#'    splitting by clustering. The `border` can be supplied as a vector,
#'    so the `border` can be applied specifically to each heatmap
#'    if needed.
#' @param iter.max integer value indicating the maximum iterations
#'    performed by k-means clustering, only relevant when `k_clusters`
#'    is non-zero.
#' @param use_raster logical indicating whether to create heatmaps
#'    using raster resizing, almost always recommended `TRUE`.
#' @param raster_by_magick `logical` passed to `ComplexHeatmap::Heatmap()`,
#'    to enable ImageMagick use during rasterization. By default this
#'    option is `TRUE` and is only disabled when the R package
#'    `"magick"` is not installed, or not properly configured.
#'    If you see a warning "instalilng 'magick' will improve rasterization"
#'    then check the R package with `library(magick)` and see if
#'    there are error messages. When `"magick"` is not available,
#'    the rasterization is substantially slower, and may produce
#'    files much larger than normal.
#' @param raster_quality `logical` passed to `ComplexHeatmap::Heatmap()`,
#'    used when `use_raster=TRUE` and defines the level of detail retained,
#'    and is used only when `raster_by_magick=FALSE`. Using larger numbers
#'    decreases speed substantially.
#' @param do_plot logical indicating whether to draw the heatmaps,
#'    where `FALSE` will return the data used to create heatmaps
#'    without actually drawing the heatmaps.
#' @param return_type character string indicating the type of
#'    data to return: `"heatmaplist"` returns the list of heatmaps,
#'    which can separately be arranged together using
#'    `ComplexHeatmap::draw()` or `grid::grid.draw()`.
#' @param show_error logical indicating whether to add error
#'    bars to the profile plot at the top of each heatmap.
#'    These error bars are calculated by
#'    `EnrichedHeatmap::anno_enriched()` using
#'    `matrixStats::colSds(x)/nrow(x)`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to
#'    `EnrichedHeatmap::EnrichedHeatmap()` to allow greater
#'    customization of details. Note that many `...` arguments
#'    are also passed to `ComplexHeatmap::Heatmap()`.
#'
#' @family jam coverage heatmap functions
#'
#' @examples
#' ## There is a small example file to use for testing
#' library(jamba)
#' cov_file1 <- system.file("data", "tss_coverage.matrix", package="platjam");
#' cov_file2 <- system.file("data", "h3k4me1_coverage.matrix", package="platjam");
#' cov_files <- c(cov_file1, cov_file2);
#' names(cov_files) <- gsub("[.]matrix",
#'    "",
#'    basename(cov_files));
#' nmatlist <- lapply(cov_files, coverage_matrix2nmat);
#' nmatlist2heatmaps(nmatlist);
#'
#' # sometimes data transform can be helpful
#' nmatlist2heatmaps(nmatlist,
#'    transform=c("log2signed", "sqrt"));
#'
#' # k-means clusters, default uses euclidean distance
#' nmatlist2heatmaps(nmatlist, k_clusters=4,
#'    transform=c("log2signed", "sqrt"));
#'
#' # k-means clusters, "correlation" or "pearson" sometimes works better
#' nmatlist2heatmaps(nmatlist,
#'    k_clusters=4,
#'    k_method="pearson",
#'    transform=c("log2signed", "sqrt"));
#'
#' # example showing usage of top_axis_side
#' # and panel_groups
#' nmatlist2 <- nmatlist[c(1, 1, 1, 2, 2, 2)];
#' names(nmatlist2) <- jamba::makeNames(names(nmatlist2))
#' for (iname in names(nmatlist2)) {
#'    attr(nmatlist2[[iname]], "signal_name") <- gsub("coverage", "cov", iname);
#' }
#' # top_axis_side="left"
#' # assumes 12x7 figure size
#' nmatlist2heatmaps(nmatlist2,
#'    signal_ceiling=0.8,
#'    nmat_colors=rep(c("firebrick", "tomato"), each=3),
#'    panel_groups=rep(c("tss", "h3k4me1"), each=3),
#'    ht_gap=grid::unit(4, "mm"),
#'    top_axis_side="left",
#'    transform=rep(c("log2signed", "sqrt"), each=3));
#'
#' # top_axis_side="both"
#' nmatlist2heatmaps(nmatlist2,
#'    panel_groups=rep(c("tss", "h3k4me1"), each=3),
#'    ht_gap=grid::unit(6, "mm"),
#'    top_axis_side="both",
#'    transform=rep(c("log2signed", "sqrt"), each=3));
#'
#' # multiple heatmap rows
#' nmatlist2heatmaps(nmatlist2,
#'    k_clusters=4,
#'    k_method="pearson",
#'    hm_nrow=2,
#'    panel_groups=rep(c("tss", "h3k4me1"), each=3),
#'    ht_gap=grid::unit(6, "mm"),
#'    top_axis_side="both",
#'    top_anno_height=grid::unit(0.8, "cm"),
#'    transform=rep(c("log2signed", "sqrt"), each=3));
#'
#' # invent anno_df data.frame of additional annotations
#' anno_df <- data.frame(
#'    tss_score=EnrichedHeatmap::enriched_score(jamba::log2signed(nmatlist[[1]])),
#'    h3k4me1_score=EnrichedHeatmap::enriched_score(jamba::log2signed(nmatlist[[2]]))
#' );
#' rownames(anno_df) <- rownames(nmatlist[[1]]);
#' nmatlist2heatmaps(nmatlist,
#'    title="k-means clustering across both heatmaps",
#'    k_clusters=4,
#'    k_method="pearson",
#'    k_heatmap=c(1, 2),
#'    ht_gap=grid::unit(6, "mm"),
#'    top_axis_side="left",
#'    anno_df=anno_df,
#'    transform=rep(c("log2signed", "sqrt"), each=3));
#'
#' @export
nmatlist2heatmaps <- function
(nmatlist,
 panel_groups=NULL,
 title=NULL,
 caption=NULL,
 upstream_length=NULL,
 downstream_length=NULL,
 k_clusters=0,
 k_subset=NULL,
 k_colors=NULL,
 k_width=grid::unit(5, "mm"),
 k_method=c("euclidean", "pearson", "correlation"),
 k_heatmap=main_heatmap,
 partition=NULL,
 rows=NULL,
 row_order=NULL,
 nmat_colors=NULL,
 middle_color="white",
 nmat_names=NULL,
 main_heatmap=1,
 anno_df=NULL,
 byCols=NULL,
 color_sub=NULL,
 anno_row_marks=NULL,
 anno_row_labels=NULL,
 top_annotation=NULL,
 top_anno_height=grid::unit(3, "cm"),
 top_axis_side=c("right"),
 legend_max_ncol=2,
 legend_base_nrow=5,
 legend_max_labels=40,
 show_heatmap_legend=TRUE,
 hm_nrow=1,
 transform="none",
 #transform=jamba::log2signed,
 signal_ceiling=NULL,
 axis_name=NULL,
 axis_name_gp=grid::gpar(fontsize=8),
 axis_name_rot=90,
 column_title_gp=grid::gpar(fontsize=12),
 lens=-2,
 anno_lens=8,
 pos_line=FALSE,
 seed=123,
 ht_gap=grid::unit(3, "mm"),
 profile_value=c("mean", "sum", "abs_mean", "abs_sum"),
 ylims=NULL,
 border=TRUE,
 iter.max=20,
 use_raster=TRUE,
 raster_quality=1,
 raster_by_magick=TRUE,
 do_plot=TRUE,
 legend_width=grid::unit(3, "cm"),
 trim_legend_title=TRUE,
 heatmap_legend_param=NULL,
 annotation_legend_param=NULL,
 return_type=c("heatmaplist", "grid"),
 show_error=FALSE,
 verbose=FALSE,
 ...)
{
   #
   return_type <- match.arg(return_type);
   profile_value <- match.arg(profile_value);
   if (length(seed) > 0) {
      set.seed(seed);
   }
   if (length(main_heatmap) == 0 || main_heatmap > length(nmatlist)) {
      main_heatmap <- 1;
   }

   if (length(border) == 0) {
      border <- FALSE;
   }
   border <- rep(border, length.out=length(nmatlist));
   if (length(legend_width) == 0) {
      legend_width <- grid::unit(3, "cm");
   }

   ## optional coordinate range zoom for coverage data
   if (length(upstream_length) > 0 || length(downstream_length) > 0) {
      if (verbose > 1) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "call zoom_nmatlist()");
      }
      nmatlist <- zoom_nmatlist(nmatlist=nmatlist,
         upstream_length=upstream_length,
         downstream_length=downstream_length);
   }

   ## k_method
   kmeans <- stats::kmeans;
   k_method <- head(k_method, 1);
   if (jamba::igrepHas("pearson|correlation|spearman|maximum|manhattan", k_method)) {
      if (!suppressPackageStartupMessages(require(amap))) {
         k_method <- "euclidean";
         jamba::printDebug("nmatlist2heatmaps(): ",
            "k_method requires the ",
            "amap",
            " package, which is not installed. Setting k_method to ",
            "euclidean");
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Using amap::Kmeans()");
      }
      kmeans <- function(centers,...){amap::Kmeans(..., centers=centers,method=k_method)};
   }
   if (length(k_method) == 0 || nchar(k_method) == 0) {
      k_method <- "euclidean";
   }
   nmat_rows <- Reduce("intersect",
      lapply(nmatlist, rownames));
   if (length(rows) == 0) {
      ## Make sure rows are present in all nmatlist entries.
      rows <- nmat_rows;
   } else {
      ## Make sure rows are present in all rownames of nmatlist
      if (is.numeric(rows)) {
         rows <- rmNA(nmat_rows[rows]);
      } else {
         rows <- rows[rows %in% nmat_rows];
      }
   }
   ## Also optionally subset rows by rownames(anno_df)
   if (length(anno_df) > 0) {
      rows <- rows[rows %in%  rownames(anno_df)];
   }
   if (length(partition) > 0) {
      if (length(names(partition)) > 0) {
         rows <- rows[rows %in% names(partition)];
         partition <- partition[match(rows, names(partition))];
      } else {
         stop("names(partition) must match rownames in nmatlist.");
      }
   }

   ## row_order must be named by rows
   if (length(row_order) > 1) {
      if (length(names(row_order)) == 0) {
         names(row_order) <- rows;
      }
   }

   if (verbose) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "Recognized ",
         jamba::formatInt(length(rows)),
         " rows shared across all matrices.");
   }

   if (length(panel_groups) > 0) {
      panel_groups <- rep(panel_groups,
         length.out=length(nmatlist));
   }
   if (length(nmat_colors) == 0) {
      if (length(panel_groups) > 0) {
         nmat_colors <- colorjam::group2colors(panel_groups,
            ...);
      } else {
         nmat_colors <- colorjam::rainbowJam(length(nmatlist),
            ...);
      }
   }
   if (length(nmat_colors) < length(nmatlist)) {
      nmat_colors <- rep(nmat_colors,
         length.out=length(nmatlist));
   }

   ## Optional transformation of each matrix
   if (length(transform) == 0) {
      transform <- function(x){x}
   }
   transform <- get_numeric_transform(transform);
   if (!is.list(transform)) {
      transform <- list(transform);
   }
   if (length(transform) != length(nmatlist)) {
      transform <- rep(transform,
         length.out=length(nmatlist));
   }
   if (verbose > 1) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "str(transform):");
      print(str(transform));
   }

   ## pos_line
   if (length(pos_line) == 0) {
      pos_line <- FALSE;
   }
   pos_line <- rep(pos_line,
      length.out=length(nmatlist));

   ## optional signal_ceiling
   if (length(signal_ceiling) > 0) {
      signal_ceiling <- rep(signal_ceiling,
         length.out=length(nmatlist));
   }
   if (length(ylims) > 0) {
      if (is.list(ylims)) {
         ylims <- rep(ylims,
            length.out=length(nmatlist));
      } else {
         ylims <- lapply(rep(ylims, length.out=length(nmatlist)), function(ylim){
            range(c(0, 0.001, ylim))
         });
      }
   }
   if (length(ylims) > 0 && verbose) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "recognized ylims: ",
         paste0("(", jamba::cPaste(ylims), ")"),
         sep="; ");
   }

   ## Define some empty variables
   PHM <- NULL;
   AHM <- NULL;
   MHM <- NULL;

   ##################################
   ## Optional k-means clustering
   if (length(k_clusters) > 0 && k_clusters > 0) {
      if (length(k_colors) == 0) {
         k_colors <- colorjam::group2colors(
            seq_len(k_clusters),
            colorSub=color_sub);
      } else if (length(k_colors) < k_clusters) {
         ## Expand the given colors using color2gradient()
         k_multiplier <- ceiling(k_clusters / length(k_colors));
         k_colors <- jamba::nameVector(
            rev(head(
               jamba::color2gradient(k_colors,
                  n=k_multiplier,
                  gradientWtFactor=1/3),
               k_clusters)),
            seq_len(k_clusters));
      } else if (length(names(k_colors)) == 0) {
         names(k_colors) <- rev(seq_len(k_clusters));
      }
      if (length(k_heatmap) == 0) {
         k_heatmap <- main_heatmap;
      }
      if (length(k_heatmap) > 1) {
         imatrix <- do.call(cbind, lapply(k_heatmap, function(k_heatmap1){
            transform[[k_heatmap1]](nmatlist[[k_heatmap1]][rows,,drop=FALSE]);
         }));
      } else {
         imatrix <- transform[[k_heatmap]](nmatlist[[k_heatmap]][rows,,drop=FALSE]);
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Running kmeans, k_clusters:",
            k_clusters,
            ", k_method:",
            k_method);
      }

      ## 28aug2020 update to kmeans cluster within partitions
      if (length(partition) == 0) {
         kpartition <- kmeans(
            imatrix,
            iter.max=iter.max,
            centers=k_clusters)$cluster;
         ## Confirm that names(partition) match rows
         names(kpartition) <- rows;
         partition <- kpartition;
         if (verbose) {
            k_sizes <- table(kpartition);
            jamba::printDebug("nmatlist2heatmaps(): ",
               "k-means cluster sizes: ",
               paste0("cluster", names(k_sizes), "=", k_sizes), sep=", ");
         }
      } else {
         partition_rows_list <- split(rows, partition[rows]);
         itransform <- transform[[main_heatmap]];
         kpartitions_list <- lapply(partition_rows_list, function(prows){
            kpartition <- kmeans(
               itransform(nmatlist[[main_heatmap]][prows,]),
               iter.max=iter.max,
               centers=k_clusters)$cluster;
            names(kpartition) <- prows;
            kpartition;
         })
         kpartition <- unlist(unname(kpartitions_list))[rows];
         partition_df <- data.frame(partition=partition[rows],
            kpartition=kpartition[rows]);
         partition <- jamba::nameVector(
            jamba::pasteByRowOrdered(partition_df,
               sep=" - "),
            rows);
         k_colors <- colorjam::group2colors(levels(partition),
            colorSub=color_sub);
      }

      ## comment out the old method
      if (1 == 2) {
         ## Optionally combine with user-defined partition
         if (length(partition) == 0) {
            partition <- kpartition;
         } else {
            if (!all(rows %in% names(partition))) {
               jamba::printDebug("head(partition):");
               print(head(partition));
               jamba::printDebug("head(rows):");
               print(head(rows));
               print(table(all(rows) %in% names(partition)));
               stop("names(partition) must match rownames in nmatlist.");
            }
            partition_df <- data.frame(partition=partition[rows],
               kpartition=kpartition[rows]);
            partition <- jamba::nameVector(
               jamba::pasteByRowOrdered(partition_df,
                  sep=" - "),
               rows);
            k_colors <- colorjam::group2colors(levels(partition),
               colorSub=color_sub);
            if (verbose) {
               jamba::printDebug("k_colors:");
               print(k_colors);
               jamba::printDebugI(jamba::nameVector(k_colors));
            }
            if (verbose) {
               k_sizes <- table(partition);
               jamba::printDebug("nmatlist2heatmaps(): ",
                  "Combined partition and k-means cluster sizes: ",
                  paste0("partition_cluster", names(k_sizes), "=", k_sizes), sep=", ");
            }
         }
      }
   }
   # End k-means clustering
   ##################################

   ##################################
   ## Partition heatmap sidebar
   if (length(partition) > 0) {
      ## Make sure to use the partition values with the properly ordered rows
      if (!all(rows %in% names(partition))) {
         jamba::printDebug("head(partition):");
         print(head(partition));
         jamba::printDebug("head(rows):");
         print(head(rows));
         print(table(all(rows) %in% names(partition)));
         stop("names(partition) must match rownames in nmatlist.");
      }
      partition <- partition[match(rows, names(partition))];
      if (!is.factor(partition)) {
         partition <- factor(partition,
            levels=sort(unique(partition)));
      } else {
         partition <- factor(partition);
      }
      ## Define colors if not provided
      if (length(k_colors) == 0) {
         k_colors <- jamba::nameVector(
            colorjam::rainbowJam(length(levels(partition)),
               ...),
            levels(partition));
         k_colors <- k_colors[sort(names(k_colors))];
      }

      ## Optional subset of k-means clusters
      if (length(k_subset) > 0) {
         if (!any(as.character(partition) %in% as.character(k_subset))) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Warning: k_subset was supplied but does not match any partition values.");
            jamba::printDebug("nmatlist2heatmaps(): ",
               "head(levels(partition), 20):",
               paste0("'", head(levels(partition), 20), "'"));
         }
         partition_keep <- (as.character(partition) %in% as.character(k_subset));
         partition <- factor(partition[partition_keep],
            levels=intersect(as.character(k_subset),
               levels(partition)));
         rows <- names(partition);
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Subsetting partition for k_subset:",
               k_subset,
               " from ",
               jamba::formatInt(length(k_colors)),
               " rows down to ",
               jamba::formatInt(length(rows)));
         }
         k_colors <- k_colors[names(k_colors) %in% as.character(k_subset)];
         ## Subset AHM and MHM if defined
         if (length(AHM) > 0) {
            AHM <- AHM[rows,];
            if (verbose) {
               jamba::printDebug("nmatlist2heatmaps(): ",
                  "Taking k_subset rows for annotation heatmap.");
            }
         }
         if (length(MHM) > 0) {
            MHM <- MHM[rows,];
            if (verbose) {
               jamba::printDebug("nmatlist2heatmaps(): ",
                  "Taking k_subset rows for annotation mark heatmap.");
            }
         }
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "creating partition heatmap PHM.");
      }
      ##################################
      ## Partition Heatmap
      p_num <- length(levels(partition[rows]));
      p_ncol <- min(c(ceiling(p_num / legend_base_nrow), legend_max_ncol));
      p_nrow <- ceiling(p_num / p_ncol);
      p_at <- jamba::mixedSort(unique(partition[rows]));
      p_labels <- gsub("\n", " ", p_at);
      p_heatmap_legend_param <- list(
         title_position="topleft",
         border="black",
         nrow=p_nrow,
         at=p_at,
         labels=p_labels
      )
      PHM <- ComplexHeatmap::Heatmap(partition[rows],
         border=FALSE,
         heatmap_legend_param=p_heatmap_legend_param,
         use_raster=use_raster,
         raster_quality=raster_quality,
         raster_by_magick=raster_by_magick,
         col=k_colors,
         name="cluster",
         show_row_names=FALSE,
         width=k_width);
      PHM_rows <- rows;
   }


   ##########################################################
   ## Optional data.frame with additional annotations
   if (length(anno_df) > 0) {
      if (!jamba::igrepHas("data.frame|dataframe|data.table|tibble", class(anno_df))) {
         stop("anno_df must be data.frame, DataFrame, data.table, or tibble.");
      }
      if (!any(rownames(anno_df) %in% rows)) {
         stop("anno_df must contain rownames present in nmatlist.");
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Preparing anno_df.");
      }
      if (!all(rows %in% rownames(anno_df)) && verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Using subset of ",
            jamba::formatInt(sum(rows %in% rownames(anno_df))),
            " rows present in rownames(anno_df).");
      }
      if (length(byCols) > 0) {
         if (is.numeric(row_order)) {
            row_rank <- match(rows, rows[row_order]);
            anno_df <- jamba::mixedSortDF(
               data.frame(check.names=FALSE,
                  anno_df[jamba::rmNA(match(rows, rownames(anno_df))),, drop=FALSE],
                  row_rank_JAM=row_rank),
               byCols=c(byCols, "row_rank_JAM"));
            anno_df <- anno_df[,setdiff(
               colnames(anno_df),
               "row_rank_JAM"), drop=FALSE];
         } else {
            anno_df <- jamba::mixedSortDF(anno_df,
               byCols=byCols);
         }
         rows <- rownames(anno_df)[rownames(anno_df) %in% rows];
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Sorted rows by:",
               byCols);
         }
         row_order <- FALSE;
      } else {
         rows <- rows[rows %in% rownames(anno_df)];
      }
      anno_df <- anno_df[match(rows, rownames(anno_df)),, drop=FALSE];
      ## Determine a list of color functions, one for each column
      anno_colors_l <- lapply(jamba::nameVector(colnames(anno_df)), function(i){
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "anno_colors_l colname:", i);
         }
         i1 <- jamba::rmNA(anno_df[[i]]);
         if (any(c("integer", "numeric") %in% class(i1))) {
            if (min(i1, na.rm=TRUE) < 0 || max(abs(i1), na.rm=TRUE) < 50) {
               ## Bi-directional color scale
               #ibreaks1 <- max(abs(i1), na.rm=TRUE);
               ibreaks1 <- quantile(abs(i1), c(0.995));
               ibreaks <- unique(seq(from=-ibreaks1,
                  to=ibreaks1,
                  length.out=25));
               if (length(ibreaks) <= 1) {
                  ibreaks <- c(-1, 0, 1);
               }
               colBR <- jamba::getColorRamp("RdBu_r",
                  lens=anno_lens,
                  trimRamp=c(2, 2),
                  n=length(ibreaks));
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " bi-directional data");
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "ibreaks:");
                  jamba::printDebugI(
                     jamba::nameVector(colBR, round(digits=2, ibreaks)),
                     sep=", ");
               }
               cBR <- circlize::colorRamp2(breaks=ibreaks,
                  col=colBR);
            } else {
               iminmax <- quantile(i1,
                  c(0.005, 0.995));
               ibreaks <- unique(seq(from=iminmax[1],
                  to=iminmax[2],
                  length.out=15));
               if (max(ibreaks) == 0 || length(ibreaks) <= 1) {
                  ibreaks <- c(0, 1);
               }
               colBR <- jamba::getColorRamp("Purples",
                  n=length(ibreaks),
                  lens=anno_lens);
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " uni-directional data");
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "ibreaks:");
                  jamba::printDebugI(
                     jamba::nameVector(colBR, round(digits=2, ibreaks)),
                     sep=", ");
               }
               cBR <- circlize::colorRamp2(breaks=ibreaks,
                  col=colBR);
            }
         } else {
            i2 <- jamba::mixedSort(unique(jamba::rmNA(i1)));
            if (all(isColor(i2))) {
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " pre-defined colors");
               }
               cBR <- jamba::nameVector(i2);
            } else if (length(color_sub) > 0 && all(i2 %in% names(color_sub))) {
               color_match <- match(i2, names(color_sub));
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " categorical data using color_sub");
               }
               cBR <- color_sub[color_match];
            } else {
               if (!"factor" %in% class(i2)) {
                  i2 <- factor(i2,
                     levels=jamba::mixedSort(i2));
               }
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "anno_colors_l colname:", i,
                     " categorical data");
               }
               cBR <- colorjam::group2colors(
                  i2,
                  colorSub=color_sub);
               cBR <- jamba::rmNA(cBR);
            }
         }
         cBR;
      });
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "sdim(anno_colors_l):");
         print(sdim(anno_colors_l));
         str(anno_colors_l);
      }

      ## annotation_legend_param
      if (length(annotation_legend_param) == 0) {
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Defining annotation_legend_param.");
         }
         ## list show_legend
         annotation_show_legend <- sapply(anno_colors_l, function(ac1){
            if (!is.function(ac1) & length(ac1) > legend_max_labels) {
               FALSE;
            } else {
               TRUE;
            }
         });
         ## list of annotation_legend_param
         annotation_legend_param <- lapply(jamba::nameVector(colnames(anno_df)), function(i){
            i1 <- jamba::rmNA(anno_df[[i]]);
            a_num <- length(setdiff(unique(i1), c("", NA)));
            a_ncol <- min(c(ceiling(a_num / legend_base_nrow), legend_max_ncol));
            a_nrow <- ceiling(a_num / a_ncol);
            i_title <- jamba::cPaste(strwrap(i, width=15), sep="\n");
            if (a_num <= 10) {
               ## display distinct steps
               if (is.numeric(i1)) {
                  i1_at <- sort(unique(i1));
                  i1_labels <- jamba::formatInt(i1_at);
               } else {
                  i1_at <- jamba::mixedSort(unique(i1));
                  i1_labels <- gsub("\n", " ", i1_at);
               }
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "annotation_legend_param colname:", i,
                     " discrete numeric color legend, a_num:",
                     a_num,
                     ", a_nrow:",
                     a_nrow);
               }
               list(
                  title=i_title,
                  title_position="topleft",
                  at=i1_at,
                  labels=i1_labels,
                  color_bar="discrete",
                  border="black",
                  nrow=a_nrow
               );
            } else if (any(c("integer", "numeric") %in% class(i1))) {
               ## display continuous color gradient
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "annotation_legend_param colname:", i,
                     " continuous color legend");
               }
               list(
                  direction="horizontal",
                  title=i,
                  labels_rot=90,
                  legend_width=legend_width,
                  title_position="topleft",
                  border="black",
                  grid_width=grid::unit(1, "npc"));
            } else {
               if (verbose) {
                  jamba::printDebug("nmatlist2heatmaps(): ",
                     "annotation_legend_param colname:", i,
                     " discrete categorical color legend");
               }
               if (is.numeric(i1)) {
                  i1_at <- sort(unique(i1));
                  i1_labels <- jamba::formatInt(i1_at);
               } else {
                  i1_at <- jamba::mixedSort(unique(i1));
                  i1_labels <- gsub("\n", " ", i1_at);
               }
               list(
                  title=i,
                  title_position="topleft",
                  at=i1_at,
                  labels=i1_labels,
                  border="black",
                  nrow=a_nrow
               )
               #   grid_width=grid::unit(1, "npc")
            }
         });
      }
      for (jj in 1:ncol(anno_df)) {
         i1 <- anno_df[[jj]];
         if (any(c("character") %in% class(i1))) {
            i1 <- factor(i1,
               levels=jamba::mixedSort(unique(i1)));
            anno_df[[jj]] <- i1;
         }
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Creating anno_df heatmap AHM.");
      }
      ##################################
      ## Annotation heatmap
      AHM <- ComplexHeatmap::rowAnnotation(
         df=anno_df[rows,,drop=FALSE],
         annotation_legend_param=annotation_legend_param,
         show_legend=annotation_show_legend,
         name="Annotation",
         col=anno_colors_l);
      AHM_rows <- rows;

      ## Optional row marks
      anno_rows <- rows[rows %in% anno_row_marks];
      if (length(anno_rows) > 0) {
         anno_row_which <- match(anno_rows, rows);
         if (length(anno_row_labels) > 0 && all(anno_row_labels %in% colnames(anno_df))) {
            anno_row_labels <- pasteByRow(
               anno_df[anno_rows,anno_row_labels,drop=FALSE],
               sep=" ");
         } else if (length(anno_row_labels) >= length(anno_rows)) {
            anno_row_labels <- anno_row_labels[anno_rows];
         } else {
            anno_row_labels <- anno_rows;
         }
         ## Print optional verbose output
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Preparing row marks for ",
               jamba::formatInt(length(anno_rows)),
               " anno_rows found in anno_df, top 20 entries are shown:");
            print(head(
               data.frame(
                  anno_rows=anno_rows,
                  anno_row_which=anno_row_which,
                  anno_row_labels=anno_row_labels),
               20));
         }
         ##################################
         ## Mark Heatmap
         MHM <- ComplexHeatmap::Heatmap(jamba::nameVector(anno_df[rows,1], rows),
            col=anno_colors_l[[1]],
            name=colnames(anno_df)[1],
            show_row_names=FALSE,
            width=k_width,
            cluster_rows=FALSE,
            right_annotation=ComplexHeatmap::rowAnnotation(
               foo=anno_mark(at=anno_row_which,
                  labels=anno_row_labels)
            )
         );
         MHM_rows <- rows;
      } else {
         MHM <- NULL;
      }
   } else {
      AHM <- NULL;
      MHM <- NULL;
   }


   ##################################
   ## panel_groups
   if (length(panel_groups) > 0) {
      ## Make sure we have some duplicated panel_groups
      if (length(jamba::tcount(panel_groups, minCount=2)) > 0) {
         if (length(show_heatmap_legend) <= 1) {
            if (length(show_heatmap_legend) == 0) {
               show_heatmap_legend <- TRUE;
            }
            show_heatmap_legend <- ifelse(duplicated(panel_groups),
               FALSE,
               show_heatmap_legend);
         }
         panel_split <- split(seq_along(nmatlist), panel_groups);
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Defining ylims for panel_groups.");
            jamba::printDebug("nmatlist2heatmaps(): ",
               "panel_groups:",
               paste0("(", jamba::cPaste(panel_split), ")"),
               sep="; ")
         }
         ## Calculate shared ylims per panel_group
         panel_ylims <- lapply(panel_split, function(idxs){
            idx_ranges <- lapply(idxs, function(idx){
               nmat <- nmatlist[[idx]][rows,,drop=FALSE];
               itransform <- transform[[idx]];
               if (length(unique(partition)) > 1) {
                  if (verbose>1) {
                     jamba::printDebug("nmatlist2heatmaps(): ",
                        "Calculating colMeans() for panels (",
                        jamba::cPaste(idxs),
                        ") across row clusters.");
                  }
                  ## jamba::rmNULL() removes empty elements if factor levels
                  ## are not present in the rows of data being used
                  plist <- jamba::rmNULL(split(names(partition),
                     partition));
                  range(
                     unlist(
                        lapply(plist, function(prows){
                           range(colMeans(
                              itransform(nmat[prows,,drop=FALSE]),
                              na.rm=TRUE))
                        })
                     )
                  )
               } else {
                  if (verbose>1) {
                     jamba::printDebug("nmatlist2heatmaps(): ",
                        "Calculating colMeans() for panels (",
                        jamba::cPaste(idxs),
                        ") across all rows.");
                  }
                  range(
                     colMeans(
                        itransform(nmat),
                        na.rm=TRUE))
               }
            })
            idx_range <- range(pretty(unlist(idx_ranges)));
            idx_range;
         });
         if (length(ylims) == 0) {
            ylims <- panel_ylims[panel_groups];
         }
         if (length(ylims) > 0 && verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "panel_groups ylims: ",
               paste0("(", jamba::cPaste(ylims), ")"),
               sep="; ");
         }
         ## Calculate signal ceiling per panel_group
         panel_ceilings <- lapply(panel_split, function(idxs){
            idx_ceilings <- lapply(idxs, function(idx){
               nmat <- nmatlist[[idx]][rows,,drop=FALSE];
               itransform <- transform[[idx]];
               iceiling <- get_nmat_ceiling(itransform(nmat),
                  signal_ceiling[[idx]],
                  verbose=verbose>1);
            });
            idx_ceiling <- max(unlist(idx_ceilings), na.rm=TRUE);
            if (verbose) {
               jamba::printDebug("nmatlist2heatmaps(): ",
                  "Determined signal ceilings for panels (",
                  jamba::cPaste(idxs),
                  "), ceilings: (",
                  jamba::cPaste(round(digits=2, unlist(idx_ceilings))),
                  "), signal_ceiling=",
                  round(digits=2, idx_ceiling));
            }
            idx_ceiling;
         })
         signal_ceiling <- panel_ceilings[panel_groups];
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "panel_ceilings[panel_groups]: ",
               paste0("(", jamba::cPaste(panel_ceilings[panel_groups]), ")"),
               sep="; ");
         }
      }
   }
   if (length(show_heatmap_legend) == 0) {
      show_heatmap_legend <- TRUE;
   }
   if (length(show_heatmap_legend) != length(nmatlist)) {
      show_heatmap_legend <- rep(show_heatmap_legend,
         length.out=length(nmatlist));
   }

   #############################################
   ## Iterate each matrix to create heatmaps
   if (length(lens) == 0) {
      lens <- 0;
   }
   lens <- rep(lens,
      length.out=length(nmatlist));
   if ("gpar" %in% class(axis_name_gp)) {
      axis_name_gp <- rep(list(axis_name_gp),
         length.out=length(nmatlist))
   } else {
      axis_name_gp <- rep(axis_name_gp,
         length.out=length(nmatlist));
   }
   if ("gpar" %in% class(column_title_gp)) {
      column_title_gp <- rep(list(column_title_gp),
         length.out=length(nmatlist));
   } else {
      column_title_gp <- rep(column_title_gp,
         length.out=length(nmatlist));
   }
   if (is.list(axis_name)) {
      axis_name <- rep(axis_name,
         length.out=length(nmatlist));
   } else {
      axis_name <- rep(list(axis_name),
         length.out=length(nmatlist));
   }

   if (length(row_order) == 0) {
      row_order <- TRUE;
   }
   if (is.logical(row_order)) {
      if (any(row_order)) {
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ", sep="",
               c("Defining row_order with ",
                  "EnrichedHeatmap::enriched_score()"));
         }
         row_order <- order(
            EnrichedHeatmap::enriched_score(nmatlist[[main_heatmap]][rows,,drop=FALSE]),
            decreasing=TRUE);
         names(row_order) <- rows;
      } else {
         row_order <- jamba::nameVector(rows);
      }
   }
   if (length(row_order) > 1) {
      row_order <- row_order[rows];
   }
   if (any(is.na(row_order))) {
      jamba::printDebug("Fixed NA row_order by assigning rows.");
      row_order <- jamba::nameVector(rows);
   }
   #} else if (isFALSE(row_order)) {
   #   row_order <- seq_along(rows);
   #}

   ###############################################
   # expand heatmap_legend_param to each heatmap
   if (length(heatmap_legend_param) == 0) {
      # prototype default heatmap legend
      heatmap_legend_direction <- "horizontal";
      heatmap_legend_param_1 <- list(
         direction=heatmap_legend_direction,
         legend_width=legend_width,
         title_position="topleft",
         border="black",
         grid_width=grid::unit(1, "npc"));

      heatmap_legend_param <- lapply(seq_along(nmatlist), function(ipanel){
         hlp_1 <- heatmap_legend_param_1;
         if (length(panel_groups) > 0) {
            hlp_1$title <- panel_groups[[ipanel]];
         } else {
            hlp_1$title <- attr(nmatlist[[ipanel]], "signal_name");
         }
         if (trim_legend_title) {
            hlp_1$title <- gsub("\n.*",
               "",
               hlp_1$title)
         }
         hlp_1;
      })
      names(heatmap_legend_param) <- names(nmatlist);
   }
   if (any(c("legend_width", "border", "direction", "title_position") %in% names(heatmap_legend_param)) ||
         length(heatmap_legend_param) != length(nmatlist)) {
      if (any(c("legend_width", "border", "direction", "title_position") %in% names(heatmap_legend_param))) {
         heatmap_legend_param <- rep(list(heatmap_legend_param),
            length.out=length(nmatlist));
      } else {
         heatmap_legend_param <- rep(heatmap_legend_param,
            length.out=length(nmatlist));
      }
   }

   ## Handle top_annotation provided as a custom list
   top_annotation_list <- NULL;
   if (is.list(top_annotation)) {
      top_annotation_list <- top_annotation;
      if (!all(names(nmatlist) %in% names(top_annotation_list))) {
         if (length(top_annotation_list) != length(nmatlist)) {
            top_annotation_list <- rep(top_annotation_list,
               length.out=length(nmatlist));
            names(top_annotation_list) <- names(nmatlist);
         }
      }
      top_annotation_list <- top_annotation_list[names(nmatlist)];
   }
   top_axis <- rep(TRUE, length(nmatlist));
   if (length(panel_groups) > 0) {
      if ("left" %in% top_axis_side) {
         top_axis <- c("blahblah", head(panel_groups, -1)) != panel_groups;
         top_axis_side <- ifelse(top_axis, "left", "left");
      } else if ("right" %in% top_axis_side) {
         top_axis <- (panel_groups != tail(c(panel_groups, "blahblah"), -1))
         top_axis_side <- ifelse(top_axis, "right", "left");
      } else if ("both" %in% top_axis_side) {
         top_axis_left <- c("blahblah", head(panel_groups, -1)) != panel_groups;
         top_axis_right <- (panel_groups != tail(c(panel_groups, "blahblah"), -1));
         top_axis <- (top_axis_left | top_axis_right);
         top_axis_side <- ifelse(top_axis_right, "right",
            ifelse(top_axis_left, "left", "left"));
      } else if ("none" %in% top_axis_side) {
         top_axis <- rep(FALSE, length.out=length(panel_groups));
         top_axis_side <- rep("left", length.out=length(panel_groups));
      } else {
         top_axis <- rep(TRUE, length.out=length(panel_groups));
         top_axis_side <- rep("right", length.out=length(panel_groups));
      }
   } else {
      if (length(top_axis_side) == 0) {
         top_axis_side <- "left";
      }
      top_axis_side <- rep(top_axis_side,
         length.out=length(nmatlist));
   }

   #############################
   ## Iterate each heatmap
   if (verbose) {
      jamba::printDebug("nmatlist2heatmaps(): ",
         "Iterating each heatmap.");
   }

   EH_l <- lapply(seq_along(nmatlist), function(i){
      nmat <- nmatlist[[i]][rows,,drop=FALSE];
      signal_name <- attr(nmat, "signal_name");
      target_name <- attr(nmat, "target_name");
      s_name <- gsub("_at_", "\nat_", signal_name);
      if (length(nmat_names) > 0) {
         signal_name <- nmat_names[i];
      }

      color <- nmat_colors[[i]];
      if (length(color) == 0 || is.na(color)) {
         color <- "aquamarine4";
      }
      ## Define ylim
      if (length(ylims) > 0) {
         ylim <- sort(ylims[[i]]);
      } else {
         ylim <- NULL;
      }
      itransform <- transform[[i]];
      imat <- itransform(nmat);
      iceiling <- signal_ceiling[[i]];
      divergent <- FALSE;
      if (any(imat < 0)) {
         divergent <- TRUE;
         if (verbose) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               c("divergent=", TRUE), sep="")
         }
      }
      if (!is.function(color) && length(color) == 1 && divergent) {
         color2 <- color_complement(color, ...);
         color <- c(color2, middle_color, color);
         divergent <- TRUE;
      }
      if (verbose) {
         if (length(ylim) == 0) {
            ylim_txt <- "NULL";
         } else {
            ylim_txt <- format(ylim, digits=2);
         }
         if (is.function(color)) {
            if (divergent) {
               color_txt <- list(color(c(-10000, 0, 10000)));
            } else {
               color_txt <- list(color(0, 10000));
            }
         } else {
            color_txt <- color;
         }
         jamba::printDebug("nmatlist2heatmaps(): ",
            "signal_name:\n'",
            signal_name,
            "',\ntarget_name:",
            target_name,
            ",\ncolor:",
            color_txt,
            ",\nylim=(", jamba::cPaste(ylim_txt), ")",
            fgText=c(
               rep(list("darkorange",
                  "dodgerblue"),
                  length.out=6),
               NA),
            bgText=as.list(
               rep(list(NA),
                  length.out=6),
               color_txt)
         );
      }
      if (length(iceiling) > 0 && !is.na(iceiling)) {
         iceiling <- get_nmat_ceiling(imat,
            iceiling,
            verbose=verbose>1);
         if (divergent) {
            ibreaks <- seq(from=-iceiling,
               to=iceiling,
               length=21);
         } else {
            ibreaks <- seq(from=0,
               to=iceiling,
               length=21);
         }
         if (is.function(color)) {
            colramp <- color;
         } else {
            colramp <- circlize::colorRamp2(
               breaks=ibreaks,
               colors=jamba::getColorRamp(color,
                  defaultBaseColor=middle_color,
                  divergent=divergent,
                  n=21,
                  lens=lens[[i]]));
         }
      } else {
         if (is.function(color)) {
            colramp <- color;
         } else {
            colramp <- jamba::getColorRamp(color,
               n=21,
               defaultBaseColor=middle_color,
               divergent=divergent,
               lens=lens[[i]]);
         }
      }

      ## if partition is a factor, call factor() which forces
      ## it to drop any missing factor levels
      if (is.factor(partition[rows])) {
         partition <- factor(partition[rows]);
      }
      use_colors <- k_colors[levels(partition)];
      if (length(use_colors) == 0) {
         use_colors <- k_colors[unique(partition)];
         if (length(use_colors) == 0) {
            n <- max(c(1, length(unique(partition))));
            use_colors <- colorjam::rainbowJam(n);
         }
      }

      ######################################
      # process top_annotation
      if (length(top_annotation_list) > 0) {
         top_annotation <- top_annotation_list[[i]];
      }
      if (length(top_annotation) > 0 && is.logical(top_annotation) && !top_annotation) {
         top_annotation <- NULL;
      } else {
         top_annotation <- ComplexHeatmap::HeatmapAnnotation(
            lines=EnrichedHeatmap::anno_enriched(
               gp=grid::gpar(col=use_colors),
               value=profile_value,
               ylim=ylim,
               axis=top_axis[i],
               axis_param=list(side=top_axis_side[i]),
               height=top_anno_height,
               show_error=show_error)
         )
      }

      ######################################
      # create coverage heatmap
      EH <- EnrichedHeatmap::EnrichedHeatmap(imat[rows,],
         split=partition[rows],
         pos_line=pos_line[[i]],
         use_raster=use_raster,
         raster_quality=raster_quality,
         raster_by_magick=raster_by_magick,
         col=colramp,
         border=border[[i]],
         top_annotation=top_annotation,
         show_heatmap_legend=show_heatmap_legend[[i]],
         heatmap_legend_param=heatmap_legend_param[[i]],
         axis_name_gp=axis_name_gp[[i]],
         axis_name=axis_name[[i]],
         axis_name_rot=axis_name_rot,
         name=signal_name,
         column_title=signal_name,
         row_order=row_order[rows],
         column_title_gp=column_title_gp[[i]],
         ...);
      EH;
   });

   ######################################
   ## Layout heatmap panels
   ##
   HM_drawn <- NULL;
   if (hm_nrow > 1 && length(nmatlist) > 1) {
      ######################################
      ## Optional multi-row layout
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Applying multi-row layout.");
      }
      hm_split <- rep(
         rep(
            seq_len(hm_nrow),
            each=ceiling(length(nmatlist) / hm_nrow)),
         length.out=length(nmatlist));
      EH_l3 <- split(EH_l, hm_split);
      HM_temp <- NULL;
      main_heatmap_temp <- main_heatmap +
         (length(partition) > 0) +
         (length(AHM) > 0);
      ht_l <- lapply(EH_l3, function(EHs){
         HM_temp <- Reduce("+", EHs);
         main_heatmap_temp <- main_heatmap;
         if (length(partition) > 0) {
            HM_temp <- PHM[match(rows, PHM_rows),] + HM_temp;
            main_heatmap_temp <- main_heatmap_temp + 1;
         }
         if (length(AHM) > 0) {
            HM_temp <- AHM[match(rows, AHM_rows),] + HM_temp;
            main_heatmap_temp <- main_heatmap_temp + 1;
         }
         ht_1 <- grid::grid.grabExpr(
            ComplexHeatmap::draw(HM_temp,
               ht_gap=ht_gap,
               main_heatmap=main_heatmap_temp));
         ht_1;
      });
      if (do_plot) {
         l <- grid::grid.layout(hm_nrow, 1);
         vp <- grid::viewport(width=1, height=1, layout=l);
         grid::grid.newpage();
         grid::pushViewport(vp);
         for (i in seq_along(ht_l)) {
            grid::pushViewport(grid::viewport(layout.pos.row=i));
            grid::grid.draw(ht_l[[i]]);
            grid::popViewport();
         }
         grid::popViewport();
      }
      if ("grid" %in% return_type) {
         EH_l <- ht_l;
      }
   } else {
      ################################
      ## Single row layout
      HM_temp <- Reduce("+", EH_l);
      main_heatmap_temp <- main_heatmap;

      ht_gap <- rep(ht_gap,
         length.out=max(c(1, length(nmatlist)-1)));
      if (length(panel_groups) > 0) {
         ht_gap_adjust <- head(
            (panel_groups != tail(c(panel_groups, "blahblah"), -1)) * 2.5 + 0.5,
            -1);
         ht_gap <- ht_gap * ht_gap_adjust;
      }
      if (length(partition) > 0) {
         HM_temp <- PHM[match(rows, PHM_rows),] + HM_temp;
         main_heatmap_temp <- main_heatmap_temp + 1;
         ht_gap <- grid::unit.c(grid::unit(1, "mm"), ht_gap);
      }
      if (length(AHM) > 0) {
         HM_temp <- AHM[match(rows, AHM_rows),] + HM_temp;
         main_heatmap_temp <- main_heatmap_temp + 1;
         ht_gap <- grid::unit.c(grid::unit(1, "mm"), ht_gap);
      }
      if (length(MHM) > 0) {
         HM_temp <- HM_temp + MHM[match(rows, MHM_rows),];
         ht_gap <- grid::unit.c(ht_gap, grid::unit(1, "mm"));
      }
      if (verbose) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "ht_gap:");
         print(ht_gap);
      }
      if (do_plot &&
            (length(title) > 0 || length(caption) > 0)) {
         jamba::printDebug("nmatlist2heatmaps(): ",
            "Preparing ComplexHeatmap::draw(HeatmapList)");
         HM_drawn <- ComplexHeatmap::draw(HM_temp,
            column_title=title,
            ht_gap=ht_gap,
            adjust_annotation_extension=TRUE,
            main_heatmap=main_heatmap_temp)
         if (FALSE) {
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Preparing HeatmapList grob for grid_with_title()");
            HM_grob <- grid::grid.grabExpr(
               ComplexHeatmap::draw(HM_temp,
                  ht_gap=ht_gap,
                  main_heatmap=main_heatmap_temp)
            );
            jamba::printDebug("nmatlist2heatmaps(): ",
               "Calling grid_with_title()");
            multienrichjam::grid_with_title(HM_grob,
               title=title,
               caption=caption,
               verbose=verbose,
               ...);
         }
      } else {
         if (do_plot) {
            HM_drawn <- ComplexHeatmap::draw(HM_temp,
               ht_gap=ht_gap,
               main_heatmap=main_heatmap_temp);
         }
      }
   }
   ret_list <- list(AHM=AHM,
      PHM=PHM,
      EH_l=EH_l,
      MHM=MHM,
      draw=list(
         HM_temp=HM_temp,
         ht_gap=ht_gap,
         main_heatmap=main_heatmap_temp)
   );
   # return HM_drawn with the heatmap as drawn
   if (length(HM_drawn) > 0) {
      ret_list$HM_drawn <- HM_drawn;
   }
   invisible(ret_list);
}
