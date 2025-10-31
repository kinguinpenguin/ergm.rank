#' @templateVar name adjacent
#' @title Try to make jumps as small as possible
#' @description Statistic calculation is more efficient if the proposed value of an edge is as close as possible to its original value.
#'
#' @usage
#' # adjacent
#'
#' @template ergmHint-general
#'
#' @concept dyad-independent
InitErgmConstraint.adjacent <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist)
  list(priority = 10, constrain = "adjacent")
}

InitErgmConstraint.ranking <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(
    nw, arglist,
    varnames = c("M"),
    vartypes = c("numeric,matrix"),
    required = c(TRUE),
    defaultvalues = list(NULL)
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
