#  File R/InitErgmReference.R in package ergm.rank, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

#' @templateVar name CompleteOrder
#' @title A uniform distribution over the possible complete orderings of the alters by each ego
#' @description The network must be directed or bipartite, and the ordering of the alters by each ego must be complete, i.e., there must not be any ties.
#'
#' @usage
#' # CompleteOrder
#'
#' @template ergmReference-general
#' @concept valued
#' @concept ordinal
InitErgmReference.CompleteOrder <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)

  if(!is.directed(nw) && !is.bipartite(nw))
    ergm_Init_stop("The network must be directed or bipartite.")

  nalters <- if (is.bipartite(nw)) network.size(nw) - nw %v% "bipartite"
             else network.size(nw) - 1L
  m <- as.matrix(nw, attrname = nw %ergmlhs% "response", matrix.type = "adjacency")
  diag(m) <- NA
  if (!all(apply(m, 1, function(x) length(unique(na.omit(x)))) == nalters))
    ergm_Init_stop("The rankings contain ties.")

  list(name="CompleteOrder", init_methods = c("CD","zeros"))
}



InitErgmReference.PartialOrder <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("classes"),
                      vartypes = c("numeric"),
                      defaultvalues = list(NULL),
                      required = c(FALSE)
                      )

  if(!is.directed(nw) && !is.bipartite(nw))
    ergm_Init_stop("The network must be directed or bipartite.")

  nalters <- if (is.bipartite(nw)) network.size(nw) - nw %v% "bipartite"
             else network.size(nw) - 1L

  if (!is.null(a$classes)) {
    classes <- a$classes
    if (classes > nalters) ergm_Init_stop("Number of distinct equivalence classes must not exceed the number of alters.")
  } else classes <- nalters

    # Extract the sociomatrix of edge weights
    v <- nw %e% (nw %ergmlhs% "response")

  if(!all(is.na(v))) {
    invalid <- abs(v - round(v)) > .Machine$double.eps^0.5 |
      v < 1 | v > nalters
    if(any(invalid, na.rm=TRUE)) {
      ergm_Init_stop("all edge weights must be integers between 1 and n-1 (inclusive) for partial order proposals.")
    }
  }

  ## TODO: Provide an API for locating and calling a reference.
  list(name = "DiscUnif", arguments = list(a = 1L, b = classes), init_methods = c("CD", "zeros"))
}

InitErgmConstraint.ranking <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(
    nw, arglist,
    varnames = c("M"),
    vartypes = c("numeric,matrix"),
    required = c(TRUE)
  )
  # M is the proposed sociomatrix
  M <- a$M

  if (!is.matrix(M))
    ergm_Init_stop("M must be a sociomatrix (numeric matrix).")

  n <- network.size(nw)
  if (!all(dim(M) == c(n, n)))
    ergm_Init_stop("M must be an n x n matrix, matching network size.")

  list(M = M)
}
