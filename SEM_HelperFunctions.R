# ============================================================
#' Helper functions for computing direct, indirect, and total effects in lavaan SEMs
#'
#' These functions provide a reproducible framework for labeling structural
#' regression coefficients, enumerating directed paths in acyclic SEMs, and
#' computing direct, indirect, and total effects using formal path products.
# ===========================================================
library(lavaan)
library(stringr)
library(igraph)

#' Generate Excel-style alphabetic coefficient labels
#'
#' Generates a sequence of alphabetic labels suitable for uniquely labeling
#' structural regression coefficients (e.g., a, b, ..., z, aa, ab, ...).
#' Labels are assigned deterministically in the order requested.
#'
#' @param n Integer giving the number of labels to generate.
#'
#' @return A character vector of length \code{n} containing alphabetic labels.
alpha_labels <- function(n) {
  stopifnot(n >= 1)
  lab_one <- function(i) {
    s <- ""
    while (i > 0) {
      i <- i - 1
      r <- i %% 26
      s <- paste0(letters[r + 1], s)
      i <- i %/% 26
    }
    s
  }
  vapply(seq_len(n), lab_one, character(1))
}

#' Compute direct, indirect, and total effects for selected outcomes
#'
#' Labels all structural regression coefficients in a lavaan model and computes
#' direct, indirect, and total effects for a specified set of outcome variables.
#' Indirect effects are defined as products of coefficients along directed paths.
#'
#' @param model_syntax lavaan model syntax (character string)
#' @param effects_to character vector of outcome variables
#' @param max_path_len maximum number of edges in indirect paths
#' @param keep_existing_labels logical; preserve existing coefficient labels
#' @param include_latents_in_paths logical; include latent variables as intermediates
#'
#' @return A list containing modified model syntax and path metadata
auto_label_and_effects_for_outcomes <- function(model_syntax,
                                                effects_to,
                                                max_path_len = Inf,
                                                keep_existing_labels = TRUE,
                                                include_latents_in_paths = TRUE) {
  
  pt <- lavaan::lavaanify(model_syntax, fixed.x = FALSE, auto = FALSE)
  
  # Identify latent variable names (LHS of "=~") (used only if include_latents_in_paths=FALSE)
  latents <- unique(pt$lhs[pt$op == "=~"])
  is_latent <- function(v) v %in% latents
  
  # Structural regressions only (exclude intercepts ~1)
  reg <- pt[pt$op == "~" & !is.na(pt$lhs) & !is.na(pt$rhs), , drop = FALSE]
  reg <- reg[reg$rhs != "1", , drop = FALSE]
  
  if (nrow(reg) == 0) {
    stop("No structural regressions (~) found in model_syntax after excluding intercepts (~1).")
  }
  
  # Ensure label column exists
  if (!("label" %in% names(reg))) reg$label <- NA_character_
  if (!keep_existing_labels) reg$label <- NA_character_
  
  # Assign alphabetic labels to unlabeled regression coefficients
  unlabeled_idx <- which(is.na(reg$label) | reg$label == "")
  if (length(unlabeled_idx) > 0) {
    reg$label[unlabeled_idx] <- alpha_labels(length(unlabeled_idx))
  }
  
  # Rebuild structural regressions with labels (lhs ~ a*x + b*z ...)
  reg_lines <- lapply(split(reg, reg$lhs), function(df) {
    terms <- paste0(df$label, "*", df$rhs)
    paste0(df$lhs[1], " ~ ", paste(terms, collapse = " + "))
  })
  reg_block <- paste(unlist(reg_lines), collapse = "\n")
  
  # Remove ONLY structural "~" lines from original syntax (keep =~, ~~ , ~1, |)
  orig_lines <- unlist(strsplit(model_syntax, "\n"))
  drop_struct <- stringr::str_detect(orig_lines, "\\s~\\s") &
    !stringr::str_detect(orig_lines, "=~") &
    !stringr::str_detect(orig_lines, "~~") &
    !stringr::str_detect(orig_lines, "~\\s*1") &
    !stringr::str_detect(orig_lines, "\\|")
  kept_lines <- orig_lines[!drop_struct]
  
  # Build directed graph: rhs -> lhs
  edges <- data.frame(from = reg$rhs, to = reg$lhs, label = reg$label, stringsAsFactors = FALSE)
  edges <- edges[!is.na(edges$from) & !is.na(edges$to), , drop = FALSE]
  
  g <- igraph::graph_from_data_frame(edges[, c("from", "to")], directed = TRUE)
  
  # Optionally remove latent nodes from allowed mediation chains
  allowed_nodes <- igraph::V(g)$name
  if (!include_latents_in_paths) {
    allowed_nodes <- allowed_nodes[!sapply(allowed_nodes, is_latent)]
  }
  g_allowed <- igraph::induced_subgraph(g, vids = allowed_nodes)
  
  # Helper: coefficient label for edge u -> v
  edge_label <- function(u, v) {
    idx <- which(edges$from == u & edges$to == v)
    if (length(idx) != 1) stop("Ambiguous/missing edge label for: ", u, " -> ", v)
    edges$label[idx]
  }
  
  safe <- function(s) gsub("[^A-Za-z0-9_]", "_", s)
  def_params <- character(0)
  
  for (y in effects_to) {
    if (!(y %in% igraph::V(g_allowed)$name)) next
    
    # IMPORTANT FIX:
    # subcomponent() returns a vertex sequence (igraph.vs), not a graph.
    anc_vs <- suppressWarnings(igraph::subcomponent(g_allowed, y, mode = "in"))
    anc_names <- setdiff(igraph::as_ids(anc_vs), y)
    
    for (x in anc_names) {
      paths <- igraph::all_simple_paths(g_allowed, from = x, to = y, mode = "out")
      if (length(paths) == 0) next
      
      path_list <- lapply(paths, igraph::as_ids)
      path_list <- Filter(function(pv) (length(pv) - 1) <= max_path_len, path_list)
      if (length(path_list) == 0) next
      
      # Direct effect if edge exists
      direct_label <- if (any(edges$from == x & edges$to == y)) edge_label(x, y) else NA_character_
      
      # Indirect paths have >= 2 edges
      indirect_paths <- Filter(function(pv) length(pv) >= 3, path_list)
      
      indirect_exprs <- character(0)
      if (length(indirect_paths) > 0) {
        indirect_exprs <- vapply(indirect_paths, function(pv) {
          labs <- vapply(seq_len(length(pv) - 1), function(i) edge_label(pv[i], pv[i + 1]), character(1))
          paste(labs, collapse = " * ")
        }, character(1))
      }
      
      xnm <- safe(x); ynm <- safe(y)
      direct_name   <- paste0("direct__",   xnm, "__", ynm)
      indirect_name <- paste0("indirect__", xnm, "__", ynm)
      total_name    <- paste0("total__",    xnm, "__", ynm)
      
      direct_def   <- if (!is.na(direct_label)) direct_label else "0"
      indirect_def <- if (length(indirect_exprs) > 0) paste(indirect_exprs, collapse = " + ") else "0"
      
      def_params <- c(def_params,
                      paste0(direct_name,   " := ", direct_def),
                      paste0(indirect_name, " := ", indirect_def),
                      paste0(total_name,    " := ", direct_name, " + ", indirect_name))
    }
  }
  
  new_model <- paste(
    paste(kept_lines, collapse = "\n"),
    "\n\n# --- Labeled structural regressions (auto-generated) ---\n",
    reg_block,
    "\n\n# --- Effects for selected outcomes (auto-generated) ---\n",
    paste(unique(def_params), collapse = "\n"),
    sep = ""
  )
  
  list(
    model = new_model,
    path_table = reg[, c("lhs", "op", "rhs", "label")],
    effects_defined = length(unique(def_params)),
    effects_to = effects_to
  )
}

#' Compute direct, indirect, and total effects for outcomes selected by name prefix
#'
#' Convenience wrapper around \code{auto_label_and_effects_for_outcomes()} that
#' selects outcome variables based on name prefixes (e.g., "gw", "sw", "fecal").
#' All upstream variables with directed paths into each selected outcome are
#' included when computing direct, indirect, and total effects.
#'
#' @param model_syntax lavaan model syntax (character string).
#' @param prefixes Character vector of prefixes used to select outcome variables.
#'   Outcomes are selected if their names begin with any of the supplied prefixes.
#' @param max_path_len Maximum number of edges allowed in indirect paths.
#'   Use to limit mediation chain length in dense DAGs.
#' @param ignore_case Logical; if \code{TRUE}, prefix matching is case-insensitive.
#' @param keep_existing_labels Logical; if \code{TRUE}, existing coefficient labels
#'   in the model syntax are preserved.
#' @param include_latents_in_paths Logical; if \code{FALSE}, latent variables are
#'   excluded as intermediates when enumerating indirect paths.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{model}: Modified lavaan model syntax with labeled regressions and
#'         derived effect definitions.
#'   \item \code{path_table}: Data frame mapping each structural regression to its
#'         assigned coefficient label.
#'   \item \code{effects_defined}: Number of derived direct, indirect, and total
#'         effects created.
#'   \item \code{effects_to}: Character vector of outcome variables used.
#' }
#'
#' @examples
#' res <- add_effects_for_prefix_outcomes(
#'   model_syntax = my_model,
#'   prefixes = c("gw", "sw"),
#'   max_path_len = 4
#' )
add_effects_for_prefix_outcomes <- function(model_syntax,
                                            prefixes = c("gw", "sw", "fecal"),
                                            max_path_len = Inf,
                                            ignore_case = FALSE,
                                            keep_existing_labels = TRUE,
                                            include_latents_in_paths = TRUE) {
  
  pt <- lavaan::lavaanify(model_syntax, fixed.x = FALSE, auto = FALSE)
  reg <- pt[pt$op == "~" & !is.na(pt$lhs) & !is.na(pt$rhs) & pt$rhs != "1", , drop = FALSE]
  
  outcomes <- sort(unique(reg$lhs))
  rx <- paste0("^(", paste(prefixes, collapse = "|"), ")")
  effects_to <- outcomes[grepl(rx, outcomes, ignore.case = ignore_case)]
  
  if (length(effects_to) == 0) {
    stop("No outcomes matched prefixes: ", paste(prefixes, collapse = ", "),
         ". Outcomes seen were: ", paste(head(outcomes, 50), collapse = ", "),
         if (length(outcomes) > 50) " ..." else "")
  }
  
  auto_label_and_effects_for_outcomes(
    model_syntax = model_syntax,
    effects_to = effects_to,
    max_path_len = max_path_len,
    keep_existing_labels = keep_existing_labels,
    include_latents_in_paths = include_latents_in_paths
  )
}