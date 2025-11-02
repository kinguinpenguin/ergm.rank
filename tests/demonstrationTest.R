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
data("newcomb")

Y_full <- newcomb[[2]]
rank_matrix <- as.matrix(Y_full, attrname = "descrank")

n <- nrow(rank_matrix)

# Create top-5 observed matrix
newcomb2.top5 <- matrix(1, nrow = n, ncol = n)

# For each ego, keep only the top 5 highest descrank values
for (i in 1:nrow(rank_matrix)) {
  valid_indices <- which(!is.na(rank_matrix[i, ]))
  if (length(valid_indices) > 0) {
    ordered <- valid_indices[order(rank_matrix[i, valid_indices], decreasing = TRUE)]
    top5_indices <- head(ordered, 5)
    newcomb2.top5[i, top5_indices] <- rank_matrix[i, top5_indices]
  }
}

newcomb2.top5
nw2 <- simulate(newcomb[[2]] # Start with newcomb time point 2.
                ~ sum, coef = 0, # A term is required, but its coefficient is set to 0.
                response="descrank", # Use edge attribute "descrank" as response.
                reference=~CompleteOrder, # Complete ordering.
                constraints = ~ adjacent + ranking(newcomb2.top5)
)

fit2 <- ergm(nw2 ~ rank.nonconformity,
             response="descrank",
             reference=~CompleteOrder,
             constraints = ~ adjacent, # Sample space: make adjacent swap proposals.
             obs.constraints = ~ ranking(newcomb2.top5) # For the constrained sampler, *also* constrain ranking.
)
