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
n <- 7
S <- 20

set.seed(0)
M <- matrix(sample.int(n - 1, n * n, replace = TRUE), n, n)
diag(M) <- 0
sm0 <- M + runif(n * n)
diag(sm0) <- 0
sm0 <- apply(sm0, 1, rank, ties.method = "random") |> t() |> (`-`)(1)

print(M)
as.matrix(sm0, attrname = "r")

nw0 <- as.network(sm0, directed = TRUE, ignore.eval = FALSE, names.eval = "r")

nws <- simulate(nw0 ~ sum(),
                coef = 0, response = "r",
                reference = ~CompleteOrder,
                constraints = ~adjacent + ranking(M),
                nsim = S)

ms <- lapply(nws, as.matrix, attrname = "r")

consistent_ranking1 <- function(x, m) {
  outer(seq_along(m), seq_along(m),
        function(i, j) {
          i == j | m[i] == m[j] |
            ((x[i] < x[j]) == (m[i] < m[j]) &
             (x[i] > x[j]) == (m[i] > m[j]))
        }) |> all()
}

consistent_ranking <- function(X, M) {
  all(sapply(seq_len(nrow(M)), function(i) consistent_ranking1(X[i, ], M[i, ])))
}

stopifnot(!statnet.common::all_identical(ms))
stopifnot(all(sapply(ms, consistent_ranking, M)))
