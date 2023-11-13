
#' Zoom the x-axis range for a list of normalizedMatrix coverage data
#'
#' Zoom the x-axis range for a list of normalizedMatrix coverage data
#'
#' This function filters the matrix columns by distance, and updates
#' important associated attributes:
#'
#' * `attr(nmat, "upstream_index")` - the column index positions upstream the target region
#' * `attr(nmat, "downstream_index")` - the column index positions downstream the target region
#' * `attr(nmat, "target_index")` - the column index positions representing the target region
#' * `attr(nmat, "extend")` - the genomic distance upstream and downstream the target region
#'
#' @family jam coverage heatmap functions
#'
#' @param nmatlist `list` of `normalizedMatrix` objects. Each
#'    `normalizedMatrix` is passed to `zoom_nmat()`.
#' @param upstream_length,downstream_length `numeric` vector whose
#'    values are recycled to length `length(nmatlist)`. Each value is
#'    passed to `zoom_nmat()` so each matrix can be zoomed to independent
#'    ranges.
#' @param ... additional arguments are passed to `zoom_nmat()`.
#'
#' @export
zoom_nmatlist <- function
(nmatlist,
   upstream_length=500,
   downstream_length=500,
   ...)
{
   #
   upstream_length <- rep(
      jamba::rmNULL(upstream_length, nullValue=NA),
      length.out=length(nmatlist));
   downstream_length <- rep(
      jamba::rmNULL(downstream_length, nullValue=NA),
      length.out=length(nmatlist));
   new_nmatlist <- lapply(seq_along(nmatlist), function(inmat){
      zoom_nmat(nmatlist[[inmat]],
         upstream_length=upstream_length[[inmat]],
         downstream_length=downstream_length[[inmat]],
         ...)
   })
   names(new_nmatlist) <- names(nmatlist);
   return(new_nmatlist);
}


#' Zoom the x-axis range for a normalizedMatrix coverage data
#'
#' Zoom the x-axis range for a normalizedMatrix coverage data
#'
#' This function is typically called by `zoom_nmatlist()` but can
#' be called on an individual `normalizedMatrix` object.
#'
#' This function filters the matrix columns by distance, and updates
#' important associated attributes:
#'
#' * `attr(nmat, "upstream_index")` - the column index positions upstream the target region
#' * `attr(nmat, "downstream_index")` - the column index positions downstream the target region
#' * `attr(nmat, "target_index")` - the column index positions representing the target region
#' * `attr(nmat, "extend")` - the genomic distance upstream and downstream the target region
#'
#' @family jam coverage heatmap functions
#'
#' @param nmat `normalizedMatrix` object, where the length extended from
#'    the target region is stored in `attr(nmat, "extend")` as a two-element
#'    integer vector representing upstream, and downstream length.
#'    Each column indicated in `attr(nmat, "upstream_index")` is expected
#'    to represent equal-sized bins spanning that range. Columns are
#'    retained if the farthest distance of the column is less
#'    than `upstream_length`.
#' @param upstream_length,downstream_length `numeric` coordinate maximum
#'    range from the target center region. When either value is `NULL`
#'    no threshold is applied, which is equivalent to `Inf`.
#'    The values are forced positive `abs(upstream_length)` as these
#'    are absolute magnitude length from the target region.
#' @param ... additional arguments are ignored.
#'
#' @export
zoom_nmat <- function
(nmat,
   upstream_length=500,
   downstream_length=500,
   ...)
{
   #
   if (length(upstream_length) == 0 || all(is.na(upstream_length))) {
      upstream_length <- Inf;
   }
   upstream_length <- abs(upstream_length);
   if (length(downstream_length) == 0 || all(is.na(downstream_length))) {
      downstream_length <- Inf;
   }
   downstream_length <- abs(downstream_length);

   # detect bin size
   bin_width <- NULL;
   if (length(attr(nmat, "upstream_index")) > 0) {
      bin_width <- round(attr(nmat, "extend")[1] /  length(attr(nmat, "upstream_index")))
      upstream_start <- rev(seq_along(attr(nmat, "upstream_index")) * bin_width);
   } else {
      upstream_start <- NULL;
   }
   upstream_keep <- (upstream_start <= upstream_length)
   new_extend1 <- max(upstream_start[upstream_keep], na.rm=TRUE);

   if (length(attr(nmat, "downstream_index")) > 0) {
      if (length(bin_width) == 0) {
         bin_width <- round(attr(nmat, "extend")[2] /  length(attr(nmat, "downstream_index")))
      }
      downstream_end <- seq_along(attr(nmat, "downstream_index")) * bin_width;
   } else {
      downstream_end <- NULL;
   }
   downstream_keep <- (downstream_end <= downstream_length)
   new_extend2 <- max(downstream_end[downstream_keep], na.rm=TRUE);
   new_extend <- c(new_extend1, new_extend2);

   target_keep <- rep(TRUE, length(attr(nmat, "target_index")));

   column_set <- c(attr(nmat, "upstream_index"),
      attr(nmat, "target_index"),
      attr(nmat, "downstream_index"))
   column_keep <- column_set[c(upstream_keep,
      target_keep,
      downstream_keep)];

   new_nmat <- nmat[, column_keep, drop=FALSE];
   attr_keep <- setdiff(names(attributes(nmat)),
      c("dim", "dimnames"));
   attributes(new_nmat)[attr_keep] <- attributes(nmat)[attr_keep];
   attr(new_nmat, "upstream_index") <- seq_len(sum(upstream_keep));
   attr(new_nmat, "downstream_index") <- tail(seq_len(ncol(new_nmat)), sum(downstream_keep));
   attr(new_nmat, "target_index") <- setdiff(seq_len(ncol(new_nmat)),
      c(attr(new_nmat, "upstream_index"),
         attr(new_nmat, "downstream_index")));
   attr(new_nmat, "extend") <- new_extend;
   return(new_nmat);
}
