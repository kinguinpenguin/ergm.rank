library(ergm)
library(ergm.rank)
library(network)

data(samplk)
samplk1
samplk2
samplk3
sampdlk1
sampdlk2
sampdlk3
adj1 <- as.matrix(samplk1, attrname = "score")
adj2 <- as.matrix(samplk2, attrname = "score")
adj3 <- as.matrix(samplk3, attrname = "score")
adj4 <- as.matrix(sampdlk1, attrname = "score")
adj5 <- as.matrix(sampdlk2, attrname = "score")
adj6 <- as.matrix(sampdlk3, attrname = "score")
write.table(adj1, file = "tests/sampson/samplk1.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(adj2, file = "tests/sampson/samplk2.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(adj3, file = "tests/sampson/samplk3.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(adj4, file = "tests/sampson/sampdlk1.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(adj5, file = "tests/sampson/sampdlk2.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(adj6, file = "tests/sampson/sampdlk3.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

fit <- ergm(samplk1 ~ rank.deference + rank.nonconformity("all")+rank.nonconformity("localAND"),
            response = "score",
            reference = ~CompleteOrder,
            constraints = ~ adjacent,
            obs.constraints = ~ ranking(adj1),
            control = snctrl(init.method = "zeros"))