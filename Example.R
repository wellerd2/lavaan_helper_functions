source("/SEM_HelperFunctions.R")

res <- add_effects_for_prefix_outcomes(
  model_syntax = fib_manpract.base,
  prefixes = c("gw", "sw", "fecal"),
  max_path_len = Inf
)

fit <- sem(res$model, data=graz_cs, estimator="ML", missing="ml", fixed.x = FALSE, check.gradient = FALSE, em.h1.iter.max=10000, em.h1.tol = 1e-08, se = "standard", cluster="region")
summary(fit , standardize=TRUE, fit.measures=TRUE, rsquare=TRUE) # FIT

pe <- parameterEstimates(fit, standardize=TRUE)
pe <- parameterEstimates(fit, standardized = TRUE)

direct_paths <- pe[pe$op == "~", ]     # all direct regression paths
effects      <- pe[pe$op == ":=", ]    # all derived indirect / total effects
all_paths <- rbind(direct_paths, effects)
# initialize new columns WITHOUT touching lhs
all_paths$effect_class <- NA   # total / indirect / direct (derived only)
all_paths$X <- NA
all_paths$Y <- NA

idx <- all_paths$op == ":="

tmp <- do.call(rbind, strsplit(all_paths$lhs[idx], "__"))

all_paths$effect_class[idx] <- tmp[, 1]   # "total", "indirect", "direct"
all_paths$X[idx] <- tmp[, 2]
all_paths$Y[idx] <- tmp[, 3]
write.csv(all_paths, file = "C:/Users/welle/OneDrive/Documents/2025 Rotational Grazing/lavaan_effects_FIB_H2OModel.csv", row.names = FALSE)
