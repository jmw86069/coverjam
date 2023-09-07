
#' Get appropriate numeric transformation function
#'
#' Get appropriate numeric transformation function
#'
#' This function recognizes numeric transformation functions by
#' name, or searches the attached R package environments for
#' a matching function by name.
#'
#' Recognized transform names:
#'
#' * `"none"` or `"linear"` returns the data without change
#' * `"log2signed"` applies `jamba::log2signed()`, defined as
#'      `log2(1+x)` transform to the absolute value, then multiplies
#'      by the original `sign(x)`
#' * `"exp2signed"` applies the inverse of `"log2signed"`, which
#'    exponentiates numeric values which were previously transformed
#'    with `log2(1+x)`. Note that value=0 when exponentiated becomes `+1`
#' * `"sqrt"` applies square root transform, equivalent to `sqrt()`
#' * `"cubert"` applies cube root `x^(1/3)`
#' * `"qrt"` applies fourth root `x^(1/4)`
#' * `"frt"` or `"fthrt"` applies fifth root `x^(1/5)`
#' * `"square"` applies `x^2` to absolute value, multiplied by the `sign(x)`
#' * `"cube"` applies `x^3`;
#'
#' Any other character name is used to find a function with the same
#' name, and if the function is found it is used. The function is
#' expected to take input `x` and return corresponding output with the
#' same length as the input. A function such as `max()` does not fit
#' this criteria, but a function such as `log2()` is acceptable.
#'
#' @returns `function` or `NULL` when no matching function is
#'    found, or `list` is returned when the input `transform`
#'    has multiple values.
#'
#' @family jam utility functions
#'
#' @param transform one of the following formats:
#'    * `character` string matching a recognized transformation:
#'
#'       * `"log2signed"`
#'       * `"exp2signed"`
#'       * `"sqrt"`
#'       * `"square"`
#'       * `"cubert"`
#'       * `"cube"`
#'       * `"qrt"`, `"quadrt"`
#'       * `"frt"`, `"fthrt"`
#'       * recognized mathematical function names (in scope), such as
#'       `log()`, `log10()`, etc. Any function retrieved by `get()`
#'       can be used, provided it returns data in the same dimensions
#'       as provided.
#'
#'    * `function`: this function is applied to the `numeric` matrix
#'    with no additional modifications. It is expected to return
#'    a `numeric` matrix with exact same dimensions.
#'    * `list` that contains a series of values with either
#'    `character` or `function` as described above.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' get_numeric_transform("log2signed")
#'
#' transform_list <- get_numeric_transform(
#'    c("none",
#'    "log2signed",
#'    log2,
#'    sqrt));
#' jamba::sdim(transform_list);
#'
#' x <- jamba::nameVector(0:16);
#' transform_list[[1]](x)
#'
#' transform_list[[2]](x)
#'
#' @export
get_numeric_transform <- function
(transform,
 ...)
{
   #
   if (length(transform) > 1) {
      if (length(names(transform)) == 0) {
         ## Determine names from the input functions
         sys <- sys.call();
         sys2 <- as.list(as.list(sys)[[2]])[-1];
         sys3 <- lapply(sys2, function(i){
            if ("call" %in% class(i)) {
               head(deparse(i), 1)
            } else {
               as.character(i)
            }
         })
         trans_names <- jamba::cPaste(sys3)
         names(transform) <- trans_names;
      }
      transform <- lapply(transform, function(it){
         get_numeric_transform(it, ...)
      });
      return(transform);
   }
   it <- transform;
   if (is.atomic(it)) {
      if ("none" %in% it) {
         it <- function(x)x;
      } else if ("log2signed" %in% it) {
         it <- jamba::log2signed;
      } else if ("exp2signed" %in% it) {
         it <- jamba::exp2signed;
      } else if ("sqrt" %in% it) {
         it <- function(x){sign(x)*sqrt(x)};
      } else if ("square" %in% it) {
         it <- function(x){sign(x)*(abs(x)^2)};
      } else if ("cubert" %in% it) {
         it <- function(x){sign(x)*abs(x)^(1/3)};
      } else if ("cube" %in% it) {
         it <- function(x){sign(x)*(abs(x)^3)};
      } else if (any(c("qrt", "quadrt") %in% it)) {
         it <- function(x){sign(x)*abs(x)^(1/4)};
      } else if (any(c("frt", "fthrt") %in% it)) {
         it <- function(x){sign(x)*abs(x)^(1/5)};
      } else {
         itf <- tryCatch({
            get(it);
         }, error=function(e){
            jamba::printDebug("get_numeric_transform(): ",
               "Error:",
               "transform name '",
               it,
               "' was not recognized, and not available on the search path:\n",
               search(),
               sep=", ")
            NULL;
         });
         if (is.function(itf)) {
            it <- function(x){sign(x)*itf(f)};
         } else {
            it <- NULL;
         }
      }
   }
   if (!is.function(it)) {
      return(NULL);
   }
   return(it);
}


#' Calculate signal ceiling of numeric matrix
#'
#' Calculate signal ceiling of numeric matrix,
#' called internally by `nmatlist2heatmaps()`.
#'
#' This function is called by `nmatlist2heatmaps()` and is not
#' intended to be called directly.
#'
#' It takes a `normalizedMatrix` or `numeric` matrix object, and
#' a ceiling value `iceiling` and determines an appropriate numeric
#' ceiling with the following rules:
#'
#' * if `iceiling=NULL` or `iceiling=NA` it returns the highest
#' absolute value in `imat`.
#' * if `iceiling > 0` and `iceiling <= 1`, it calculates the quantile
#' of the absolute values observed, using only non-zero values with
#' `quantile(abs(imat), probs=iceiling)`
#' * otherwise `iceiling` is interpreted as a fixed `numeric` ceiling
#'
#' In all cases, `iceiling` is rounded to 3 digits with
#' `round(iceiling, digits=3)`
#'
#' Also in all cases, `na.rm=TRUE` is used, to prevent returning `NA`.
#'
#' @family jam coverage heatmap functions
#'
#' @param imat `numeric` matrix or `normalizedMatrix` object.
#' @param iceiling `numeric` maximum value, interpreted as an absolute
#'    threshold for any value above `1`, and interpreted as a quantile
#'    threshold for any value above `0` and no greater than `1`.
#'    When `iceiling=NULL` or `iceiling=NA` the `numeric` maximum observed
#'    value is used as the ceiling.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
get_nmat_ceiling <- function
(imat,
 iceiling=NULL,
 verbose=TRUE,
 ...)
{
   if (length(iceiling) == 0 || any(is.na(iceiling))) {
      iceiling <- max(abs(imat),
         na.rm=TRUE);
      if (verbose) {
         jamba::printDebug("get_nmat_ceiling(): ",
            "   Applied max(nmat) ceiling=",
            round(digits=3, iceiling));
      }
   } else if (iceiling > 0 && iceiling <= 1) {
      # apply quantile
      imat_values <- setdiff(abs(imat), 0);
      iquantile <- quantile(imat_values,
         probs=iceiling,
         na.rm=TRUE);
      if (verbose) {
         jamba::printDebug("get_nmat_ceiling(): ",
            "Applied quantile=",
            iceiling,
            " and defined ceiling=",
            round(digits=3, iquantile));
      }
      iceiling <- iquantile;
   } else if (verbose) {
      jamba::printDebug("get_nmat_ceiling(): ",
         "   Applied ceiling=",
         round(digits=3, iceiling));
   }
   return(iceiling);
}
