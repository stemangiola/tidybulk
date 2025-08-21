# Helper script to inspect the airway dataset and update roxygen examples

# 1) Inspect airway dataset structure (assays and colData)
cat("\n=== Inspecting 'airway' dataset ===\n")
if (!requireNamespace("airway", quietly = TRUE)) {
  stop("Package 'airway' is required to inspect and use the airway dataset. Please install it: BiocManager::install('airway')")
}

# Load the data without attaching the package
utils::data("airway", package = "airway", envir = environment())
airway_obj <- get("airway", envir = environment())

cat(sprintf("Class: %s\n", paste(class(airway_obj), collapse = ", ")))
if (methods::is(airway_obj, "SummarizedExperiment")) {
  cat(sprintf("Assays: %s\n", paste(SummarizedExperiment::assayNames(airway_obj), collapse = ", ")))
  cat(sprintf("Dimensions: %s features x %s samples\n",
              nrow(SummarizedExperiment::assay(airway_obj)),
              ncol(SummarizedExperiment::assay(airway_obj))))
  cd <- SummarizedExperiment::colData(airway_obj)
  cat("colData columns:\n")
  print(colnames(cd))
  cat("\ncolData snapshot (first 6 rows):\n")
  print(as.data.frame(cd)[seq_len(min(6, nrow(cd))), , drop = FALSE])
} else {
  warning("'airway' is not a SummarizedExperiment; skipping detailed structure output.")
}

# 2) Update roxygen examples: replace 'airway_mini' with 'airway' in roxygen lines only
#    and adjust 'tidybulk::airway_mini' to 'airway::airway'. Also, insert a minimal
#    preamble after each @examples tag if the example block references 'airway'.

root <- normalizePath(".", mustWork = TRUE)
r_dir <- file.path(root, "R")

stopifnot(dir.exists(r_dir))

r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)

insert_preamble_lines <- function(lines, start_idx) {
  # Insert after the @examples line (start_idx)
  preamble <- c(
    "#' ## Load airway dataset for examples",
    "#' library(airway)",
    "#' data('airway', package = 'airway')",
    "#' # Create/overwrite a 'condition' column used by examples",
    "#' SummarizedExperiment::colData(airway)$condition <- as.factor(SummarizedExperiment::colData(airway)$dex)"
  )
  # Only insert if not already present close to the start of the example block
  window_end <- min(length(lines), start_idx + 6L)
  window <- paste(lines[(start_idx + 1L):window_end], collapse = "\n")
  if (!grepl("library\\(airway\\)|data\\('airway'", window)) {
    lines <- append(lines, preamble, after = start_idx)
  }
  lines
}

replace_in_roxygen <- function(lines) {
  in_roxy <- FALSE
  in_examples <- FALSE
  example_starts <- integer(0)

  # First pass: identify roxygen and @examples blocks
  for (i in seq_along(lines)) {
    if (grepl("^#'", lines[i])) {
      if (!in_roxy) in_roxy <- TRUE
      if (grepl("^#'[[:space:]]*@examples\\b", lines[i])) {
        in_examples <- TRUE
        example_starts <- c(example_starts, i)
      }
    } else {
      in_roxy <- FALSE
      in_examples <- FALSE
    }
  }

  # Second pass: do replacements only on roxygen lines
  for (i in seq_along(lines)) {
    if (grepl("^#'", lines[i])) {
      # Normalize to using the object loaded by data('airway'), not namespaced object
      # Replace any previous mapping to namespaced airway
      lines[i] <- gsub("airway::airway", "airway", lines[i], fixed = TRUE)
      # Replace fully qualified old object to plain airway
      lines[i] <- gsub("tidybulk::airway_mini", "airway", lines[i], fixed = TRUE)
      # Replace standalone airway_mini tokens with airway (word boundary)
      lines[i] <- gsub("\\bse_mini\\b", "airway", lines[i])
    }
  }

  # Third pass: insert preamble in example blocks that reference airway
  if (length(example_starts) > 0) {
    # We iterate from the last to the first to keep indices valid while inserting
    for (idx in rev(example_starts)) {
      # If within the next ~30 lines of the examples block there is an 'airway' reference,
      # add the preamble lines.
      look_ahead_end <- min(length(lines), idx + 30L)
      block_text <- paste(lines[(idx + 1L):look_ahead_end], collapse = "\n")
      if (grepl("airway", block_text, fixed = TRUE)) {
        lines <- insert_preamble_lines(lines, idx)
      }
    }
  }

  # Fourth pass: remove any conditional guards in roxygen examples and make them unconditional
  for (i in seq_along(lines)) {
    if (grepl("^#'[[:space:]]*if[[:space:]]*\\(", lines[i])) {
      # Drop lines that open conditional guards
      lines[i] <- "#'"
    }
    if (grepl("^#'[[:space:]]*}[[:space:]]*$", lines[i])) {
      # Drop closing braces added by previous conditional blocks
      lines[i] <- "#'"
    }
  }

  # Fifth pass: within example blocks, normalize batch assignment to ensure two levels
  if (length(example_starts) > 0) {
    i <- 1
    while (i <= length(lines)) {
      if (grepl("^#'", lines[i])) {
        # Detect a line that zeros batch
        if (grepl("^#'\\s*cm\\$batch\\s*=\\s*0\\s*$", lines[i])) {
          # Replace with robust two-level assignment
          lines[i] <- "#' cm$batch <- rep(c('A','A','B','B'), length.out = ncol(cm))"
          # Remove any immediate subsequent line that sets batch by colnames subset
          if (i + 1 <= length(lines) && grepl("^#'\\s*cm\\$batch\\[", lines[i + 1])) {
            lines <- lines[-(i + 1)]
          }
        }
        # Upgrade method argument to combat_seq in examples
        if (grepl("method\\s*=\\s*\"combat\"", lines[i])) {
          lines[i] <- sub("method\\s*=\\s*\"combat\"", "method=\"combat_seq\"", lines[i])
        }
      }
      i <- i + 1
    }
  }

  # Sixth pass: Ensure AnnotationDbi/org.Hs.eg.db blocks are unconditional
  for (i in seq_along(lines)) {
    if (grepl("^#'[[:space:]]*if", lines[i]) && grepl("AnnotationDbi", lines[i])) {
      # Replace the guard with library calls
      lines[i] <- "#' library(AnnotationDbi)"
      # Insert org.Hs.eg.db library on next line
      lines <- append(lines, "#' library(org.Hs.eg.db)", after = i)
    }
  }

  lines
}

cat("\n=== Updating roxygen examples in R/ ===\n")
changed_files <- character(0)
for (f in r_files) {
  original <- readLines(f, warn = FALSE)
  modified <- replace_in_roxygen(original)
  if (!identical(original, modified)) {
    writeLines(modified, f)
    changed_files <- c(changed_files, f)
    cat(sprintf("Updated: %s\n", basename(f)))
  }
}

if (length(changed_files) == 0) {
  cat("No roxygen examples required changes.\n")
} else {
  cat(sprintf("\nTotal files updated: %d\n", length(changed_files)))
}

cat("\nDone. You can now regenerate Rd files with devtools::document() if needed.\n")


