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

#' @templateVar name AdjacentAlterSwap
#' @aliases InitWtErgmProposal.AdjacentAlterSwap
#' @title A proposal that swaps values of two adjacently-ranked alters incident on an ego
#' @description This proposal randomly selects a dyad \eqn{(i,j)}, then finds an alter \eqn{j'} ranked immediately above \eqn{j}, and proposes to swap their values.
#' @template ergmProposal-general
#' @concept valued
#' @concept ordinal
NULL

InitWtErgmProposal.AdjacentAlterSwap <- function(arguments, nw) {
  M <- arguments$constraints$ranking$M
  MHproposal <- list(name = "AdjacentAlterSwap", inputs = M, auxiliaries = ~ .sociomatrix("numeric") + .updown)
  MHproposal
}

InitWtErgmProposal.PartialDisc <- function(arguments, nw) {
  M <- arguments$constraints$ranking$M
  list(name = "PartialDisc", iinputs = unlist(arguments$reference$arguments),
       dyadgen = ergm_dyadgen_select(arguments, nw))
}
