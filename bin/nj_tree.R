#!/usr/bin/env Rscript

library(ape)

# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    result <- list(input = NULL, format = NULL, prefix = "nj_tree")
    i <- 1
    while (i <= length(args)) {
        if      (args[i] == "--input")  { result$input  <- args[i + 1]; i <- i + 2 }
        else if (args[i] == "--format") { result$format <- args[i + 1]; i <- i + 2 }
        else if (args[i] == "--prefix") { result$prefix <- args[i + 1]; i <- i + 2 }
        else i <- i + 1
    }
    result
}

params <- parse_args(args)

if (is.null(params$input))  stop("--input is required")
if (is.null(params$format)) stop("--format is required (fastani or ska2)")

# --------------------------------------------------------------------------
# Parsers: return a symmetric numeric distance matrix
# --------------------------------------------------------------------------

parse_fastani <- function(file) {
    df <- read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("query", "ref", "ani", "q_frags", "r_frags"))

    # Strip directory path and file extension to get clean sample names
    df$query <- tools::file_path_sans_ext(basename(df$query))
    df$ref   <- tools::file_path_sans_ext(basename(df$ref))

    # Remove self-comparisons
    df <- df[df$query != df$ref, ]

    # Convert ANI (%) to distance (0-1); average both directions for symmetry
    df$dist <- (100 - df$ani) / 100

    samples <- sort(unique(c(df$query, df$ref)))
    n       <- length(samples)
    mat     <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))

    for (i in seq_len(nrow(df))) {
        a <- df$query[i]
        b <- df$ref[i]
        # Both (A→B) and (B→A) are present; accumulate and later divide by 2
        mat[a, b] <- mat[a, b] + df$dist[i] / 2
        mat[b, a] <- mat[b, a] + df$dist[i] / 2
    }
    mat
}

parse_ska2 <- function(file) {
    df <- read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                     check.names = FALSE)

    # Use the normalised mismatch proportion as the distance
    dist_col <- "Mismatches (proportion)"
    if (!dist_col %in% colnames(df)) {
        stop(paste("Expected column not found:", dist_col,
                   "\nColumns present:", paste(colnames(df), collapse = ", ")))
    }

    samples <- sort(unique(c(df$Sample1, df$Sample2)))
    n       <- length(samples)
    mat     <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))

    for (i in seq_len(nrow(df))) {
        a <- df$Sample1[i]
        b <- df$Sample2[i]
        d <- df[[dist_col]][i]
        mat[a, b] <- d
        mat[b, a] <- d
    }
    mat
}

# --------------------------------------------------------------------------
# Build distance matrix
# --------------------------------------------------------------------------

mat <- if (params$format == "fastani") {
    parse_fastani(params$input)
} else if (params$format == "ska2") {
    parse_ska2(params$input)
} else {
    stop(paste("Unknown --format:", params$format, "(must be 'fastani' or 'ska2')"))
}

if (nrow(mat) < 3) {
    stop(paste("At least 3 samples are required for NJ tree construction; got", nrow(mat)))
}

# --------------------------------------------------------------------------
# NJ tree
# --------------------------------------------------------------------------

tree <- nj(as.dist(mat))

# Write Newick
write.tree(tree, file = paste0(params$prefix, ".nwk"))

# Write PDF plot
pdf(paste0(params$prefix, ".pdf"), width = 8, height = 6)
plot.phylo(tree,
           main  = paste("Neighbour-joining tree (", params$format, "distances)"),
           cex   = 0.8,
           type  = "phylogram")
add.scale.bar()
invisible(dev.off())

message("NJ tree written to: ", params$prefix, ".nwk")
message("Tree plot written to: ", params$prefix, ".pdf")
