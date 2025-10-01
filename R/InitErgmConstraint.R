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
