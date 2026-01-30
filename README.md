# lavaan_helper_functions

This repository provides helper functions to:
1) label all structural regression coefficients (`~`) in a lavaan model using Excel-style labels (`a`, `b`, …, `aa`, …),
2) build the corresponding directed regression graph,
3) define direct, indirect, and total effects for selected outcomes via derived parameters (`:=`) by enumerating all directed paths.

## Why the labels are correct

The labels are correct because they are mechanically tied to the fitted regression graph, not inferred heuristically. Specifically:

- Coefficient labels are assigned only to structural regression terms (`lhs ~ rhs`), not to loadings (`=~`), covariances (`~~`), intercepts (`~1`), or thresholds (`|`).
- The directed graph is constructed from the same labeled regression table, with edges `rhs → lhs`.  
  Therefore, each edge corresponds to exactly one labeled coefficient.
- Indirect effects are defined as formal products of labeled coefficients along directed paths.  
  Total effects are defined as the sum of direct and indirect effects.

## What is computed

For each selected outcome `Y`:
- Direct effect for each upstream `X` (if `Y ~ X` exists)
- Indirect effect for each upstream `X` (sum of products for all directed paths of length ≥ 2 from `X` to `Y`)
- Total effect = direct + indirect

## Assumptions and scope

- The structural model is a directed acyclic graph (no reciprocal paths).
- Indirect effects are computed only along directed regression paths (`~`).
- Correlated residuals (`~~`) are not treated as causal pathways.
- Categorical variables are modeled using the estimator specified by the user.

> Note: In dense DAGs, the number of indirect paths can be large. Use `max_path_len` to cap mediation chain length if needed.

## Installation

```r
install.packages(c("lavaan", "igraph", "stringr"))

