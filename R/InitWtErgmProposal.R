#  File R/InitWtErgmProposal.R in package ergm.rank, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################


#' @templateVar name AlterSwap
#' @aliases InitWtErgmProposal.AlterSwap
#' @title A proposal that swaps values of two alters incident on an ego
#' @description This proposal randomly selects two dyads \eqn{(i,j)} and \eqn{(i,j')} with a common sender and proposes to swap their values if distinct.
#' @template ergmProposal-general
#' @concept valued
#' @concept ordinal
NULL

InitWtErgmProposal.AlterSwap <- function(arguments, nw) {
  MHproposal <- list(name = "AlterSwap", inputs=NULL)
  MHproposal
}

InitWtErgmProposal.PartialDisc  <- function(arguments, nw) {
  params <- with(arguments$reference$arguments,
                 switch(arguments$reference$name,
                        DiscUnif = list(3L, c(a, b))))

  list(name = "PartialDisc", iinputs = params[[1]], inputs=as.double(params[[2]]), dyadgen = ergm_dyadgen_select(arguments, nw))
}


#' @templateVar name AdjacentAlterSwap
#' @aliases InitWtErgmProposal.AdjacentAlterSwap
#' @title A proposal that swaps values of two adjacently-ranked alters incident on an ego
#' @description This proposal randomly selects a dyad \eqn{(i,j)}, then finds an alter \eqn{j'} ranked immediately above \eqn{j}, and proposes to swap their values.
#' @template ergmProposal-general
#' @concept valued
#' @concept ordinal
NULL

InitWtErgmProposal.AdjacentAlterSwap <- function(arguments, nw) {
  MHproposal <- list(name = "AdjacentAlterSwap", inputs = NULL, auxiliaries = ~ .sociomatrix("numeric") + .updown)
  MHproposal
}

InitErgmReference.PartialOrder <- function(nw, response=NULL, ...) {
  n <- network.size(nw)
  
  # The range of allowable values:
  min_val <- 1
  max_val <- n - 1
  
  # Return a reference list (like other InitErgmReference.* functions)
  list(
    name = "PartialOrder",
    parameters = c(min_val, max_val),
  )
}

InitWtErgmProposal.PartialDisc <- function(arguments, nw) {
  # Extract the sociomatrix of edge weights
  sm <- as.matrix(nw %n% "weights")
  if(is.null(sm)) sm <- as.matrix(as.sociomatrix(nw, attrname="weights", expand.bipartite=TRUE))

  # If the network has weights, check that they're all integers
  if(!all(is.na(sm))) {
    noninteger <- abs(sm - round(sm)) > .Machine$double.eps^0.5
    if(any(noninteger, na.rm=TRUE)) {
      stop("InitWtErgmProposal.PartialDisc: all edge weights must be integers for partial order proposals.")
    }
  }
  MHproposal <- list(
    name = "PartialDisc",
    inputs = NULL,
    auxiliaries = ~ .sociomatrix("numeric") + .updown
  )
  MHproposal
}

InitWtErgmProposal.AdjacentAlterSwapPartial <- function(arguments, nw) {
  MHproposal <- list(name = "MH_AdjacentAlterSwapPartial", inputs = NULL, auxiliaries = ~ .sociomatrix("numeric") + .updown)
  MHproposal
}
