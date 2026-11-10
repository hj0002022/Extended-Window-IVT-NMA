library(gemtc)
library(rjags)
library(coda)

comparisons <- data.frame(
  t1 = c("CONTROL", "CONTROL", "TNK"),
  t2 = c("ALT", "TNK", "ALT"),
  comparison = c("ALT vs CONTROL", "TNK vs CONTROL", "ALT vs TNK")
)

outcomes <- c("mRS01", "mRS02", "ENI", "SICH", "Death", "mRS56")
or_results <- vector("list", length(outcomes))
sucra_results <- vector("list", length(outcomes))

for (i in seq_along(outcomes)) {
  outcome <- outcomes[i]
  d <- read.csv(file.path("Data", paste0(outcome, ".csv")), stringsAsFactors = FALSE)
  d$treatment <- unname(c("ALT 0.9 mg/kg" = "ALT", "TNK 0.25 mg/kg" = "TNK", "CONTROL" = "CONTROL")[d$treatment])
  d <- d[, c("study", "treatment", "responders", "sampleSize")]
  network <- mtc.network(data.ab = d)
  model <- mtc.model(network, type = "consistency", n.chain = 4,
                     likelihood = "binom", link = "logit", linearModel = "random")
  set.seed(20261107 + i)
  fit <- mtc.run(model, n.adapt = 20000, n.iter = 50000, thin = 1)

  or_results[[i]] <- do.call(rbind, lapply(seq_len(nrow(comparisons)), function(j) {
    effect <- relative.effect(fit, t1 = comparisons$t1[j],
                              t2 = comparisons$t2[j], preserve.extra = FALSE)
    samples <- as.matrix(as.mcmc.list(effect))
    effect_column <- grep("^d\\.", colnames(samples), value = TRUE)
    stopifnot(length(effect_column) == 1L)
    q <- exp(quantile(samples[, effect_column], c(0.025, 0.5, 0.975), names = FALSE))
    data.frame(outcome = outcome, comparison = comparisons$comparison[j],
               OR = q[2], CrI_lower = q[1], CrI_upper = q[3])
  }))

  ranking <- rank.probability(fit, preferredDirection = 1)
  scores <- sucra(ranking)
  sucra_results[[i]] <- data.frame(outcome = outcome,
                                   treatment = names(scores),
                                   SUCRA = as.numeric(scores))
}

or_results <- do.call(rbind, or_results)
sucra_results <- do.call(rbind, sucra_results)

write.csv(or_results, "Main_Analysis_OR.csv", row.names = FALSE)
write.csv(sucra_results, "Main_Analysis_SUCRA.csv", row.names = FALSE)

print(or_results, row.names = FALSE)
print(sucra_results, row.names = FALSE)
