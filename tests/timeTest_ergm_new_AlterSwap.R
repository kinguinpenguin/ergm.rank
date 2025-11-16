library(ergm.rank)
data("newcomb")
fit_true <- list()
fit_constrained <- list()
total_time <- 0
for (week in 1:15) {
  fit.time <- system.time(fit2.true <- ergm(newcomb[[week]] ~ rank.deference+rank.nonconformity("all")+
                      rank.nonconformity("localAND"),
              response="descrank",
              reference=~CompleteOrder, # Sample space: make adjacent swap proposals.
              constraints = ~ adjacent,
              control = snctrl()
  ))
  total_time <- total_time + fit.time
}
average_time <- total_time / 15
cat(sprintf("Total time: %.2f seconds\n", average_time))