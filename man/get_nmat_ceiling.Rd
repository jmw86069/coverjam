% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformation_ceiling.R
\name{get_nmat_ceiling}
\alias{get_nmat_ceiling}
\title{Calculate signal ceiling of numeric matrix}
\usage{
get_nmat_ceiling(imat, iceiling = NULL, verbose = TRUE, ...)
}
\arguments{
\item{imat}{\code{numeric} matrix or \code{normalizedMatrix} object.}

\item{iceiling}{\code{numeric} maximum value, interpreted as an absolute
threshold for any value above \code{1}, and interpreted as a quantile
threshold for any value above \code{0} and no greater than \code{1}.
When \code{iceiling=NULL} or \code{iceiling=NA} the \code{numeric} maximum observed
value is used as the ceiling.}

\item{verbose}{\code{logical} indicating whether to print verbose output.}

\item{...}{additional arguments are ignored.}
}
\description{
Calculate signal ceiling of numeric matrix,
called internally by \code{nmatlist2heatmaps()}.
}
\details{
This function is called by \code{nmatlist2heatmaps()} and is not
intended to be called directly.

It takes a \code{normalizedMatrix} or \code{numeric} matrix object, and
a ceiling value \code{iceiling} and determines an appropriate numeric
ceiling with the following rules:
\itemize{
\item if \code{iceiling=NULL} or \code{iceiling=NA} it returns the highest
absolute value in \code{imat}.
\item if \code{iceiling > 0} and \code{iceiling <= 1}, it calculates the quantile
of the absolute values observed, using only non-zero values with
\code{quantile(abs(imat), probs=iceiling)}
\item otherwise \code{iceiling} is interpreted as a fixed \code{numeric} ceiling
}

In all cases, \code{iceiling} is rounded to 3 digits with
\code{round(iceiling, digits=3)}

Also in all cases, \code{na.rm=TRUE} is used, to prevent returning \code{NA}.
}
\seealso{
Other jam coverage heatmap functions: 
\code{\link{coverage_matrix2nmat}()},
\code{\link{nmatlist2heatmaps}()},
\code{\link{zoom_nmatlist}()},
\code{\link{zoom_nmat}()}
}
\concept{jam coverage heatmap functions}
