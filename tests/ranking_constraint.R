#  File tests/termTests_rank.R in package ergm.rank, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################
library(ergm.rank)
n <- 4
S <- 2

M <- matrix(sample.int(n-1, n*n, replace=TRUE), n, n)
sm0 <- M + runif(n*n)
diag(sm0) <- 0
sm0 <- apply(sm0, 1, rank, ties.method = "random") |> t() |> (`-`)(1)

print(M)

nw0 <- as.network(sm0, directed = TRUE, ignore.eval = FALSE, names.eval = "r")

nws <- simulate(nw0~sum(),
                coef = 0, response = "r",
                reference = ~CompleteOrder, constraints = ~adjacent + ranking(M), nsim = S)

ms <- lapply(nws, as.matrix, attrname = "r")

consistent_ranking1 <- function(x, m) {
  sapply(seq_along(m), function(i)
    sapply(seq_along(m), function(j)
      i==j || (m[i]==m[j]) || ( (x[i]<x[j])==(m[i]<m[j])  && (x[i]>x[j])==(m[i]>m[j]) )
      )
    ) |> all()
}

consistent_ranking <- function(X, M) {
  sapply(seq_len(nrow(M)), function(i) consistent_ranking1(X[i,], M[i,])) |> all()
}

lapply(ms, consistent_ranking, M)
