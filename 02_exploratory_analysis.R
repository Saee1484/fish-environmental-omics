# ============================================================
# Project: fish-environmental-omics
# Script: 02_exploratory_analysis.R
# Dataset: GSE116687
# Purpose: Exploratory analysis of liver and testis expression
# ============================================================

# 1. Load packages --------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

# Confirm that the data objects already exist
stopifnot(exists("liver"), exists("testis"))

# Define output directory
project_dir <- file.path(
  path.expand("~"),
  "fish-environmental-omics"
)

figure_dir <- file.path(
  project_dir, "results", "figures"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Prepare expression data ----------------------------------

prepare_tissue <- function(dat, tissue_name) {

  stopifnot("GENE-ID" %in% names(dat))

  sample_cols <- setdiff(names(dat), "GENE-ID")

  if (!all(vapply(dat[sample_cols], is.numeric, logical(1)))) {
    stop("Expression columns must be numeric: ", tissue_name)
  }

  expr <- as.matrix(dat[sample_cols])
  rownames(expr) <- make.unique(as.character(dat[["GENE-ID"]]))

  if (any(!is.finite(expr)) || any(expr < 0)) {
    stop("Missing, non-finite or negative values found in ",
         tissue_name, ". Inspect the data before plotting.")
  }

  group <- ifelse(grepl("CTRL$", sample_cols), "Control",
                  ifelse(grepl("BaP$", sample_cols),
                         "BaP", NA_character_))

  if (anyNA(group)) {
    stop("Could not identify all sample groups in ", tissue_name)
  }

  metadata <- data.frame(
    sample = sample_cols,
    group = factor(group, levels = c("Control", "BaP")),
    tissue = tissue_name
  )

  list(expr = expr, metadata = metadata)
}

liver_data <- prepare_tissue(liver, "Liver")
testis_data <- prepare_tissue(testis, "Testis")

# 3. Save sample metadata -------------------------------------

project_dir <- "C:/Users/YASH CHAVAN/OneDrive/Documents/fish-environmental-omics"

table_dir <- file.path(project_dir, "results", "tables")

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Check whether the folder now exists
dir.exists(table_dir)
test_path <- file.path(table_dir, "write_test.txt")

writeLines("Testing project output.", test_path)
if (exists("metadata_all")) {
  write.csv(
    metadata_all,
    file.path(table_dir, "sample_metadata.csv"),
    row.names = FALSE
  )
  cat("Sample metadata saved successfully.\n")
} else {
  message("metadata_all does not exist yet. We need to recreate it.")
}
file.exists(test_path)

metadata_all <- bind_rows(
  liver_data$metadata,
  testis_data$metadata
)

write.csv(
  metadata_all,
  file.path(project_dir, "results", "tables", "sample_metadata.csv"),
  row.names = FALSE
)


# 4. Helper function for saving figures -----------------------

save_plot <- function(plot, folder, filename,
                      width = 8, height = 6) {

  out_dir <- file.path(figure_dir, folder)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ggsave(
    file.path(out_dir, paste0(filename, ".png")),
    plot, width = width, height = height,
    dpi = 300, bg = "white"
  )

  ggsave(
    file.path(out_dir, paste0(filename, ".pdf")),
    plot, width = width, height = height,
    bg = "white"
  )
}


# 5. Create exploratory figures -------------------------------

plot_tissue <- function(obj, tissue_name) {

  expr <- obj$expr
  metadata <- obj$metadata

  # Log transformation for visualisation only
  log_expr <- log2(expr + 1)

  # Convert to long format
  long_data <- as.data.frame(log_expr) |>
    tibble::rownames_to_column("gene_id") |>
    pivot_longer(
      cols = -gene_id,
      names_to = "sample",
      values_to = "log2_expression"
    ) |>
    left_join(metadata, by = "sample")

  prefix <- tolower(tissue_name)

  # Figure 1: Sample overview
  p1 <- ggplot(metadata, aes(x = sample, fill = group)) +
    geom_bar() +
    labs(
      title = paste(tissue_name, "- sample overview"),
      x = "Sample",
      y = "Number of samples",
      fill = "Group"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(p1, "01_data_overview",
            paste0(prefix, "_01_sample_overview"))

  # Figure 2: Expression boxplots
  p2 <- ggplot(long_data,
               aes(x = sample, y = log2_expression, fill = group)) +
    geom_boxplot(outlier.size = 0.25, linewidth = 0.35) +
    labs(
      title = paste(tissue_name, "- expression distributions"),
      subtitle = "log2(expression + 1); for visualisation",
      x = "Sample",
      y = "log2(expression + 1)",
      fill = "Group"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(p2, "02_expression_distributions",
            paste0(prefix, "_02_expression_boxplots"))

  # Figure 3: Density plots
  p3 <- ggplot(long_data,
               aes(x = logfurther2_expression, colour = group,
                   group = sample)) +
    geom_density(linewidth = 0.6, na.rm = TRUE) +
    facet_wrap(~sample, ncol = 2) +
    labs(
      title = paste(tissue_name, "- sample expression densities"),
      x = "log2(expression + 1)",
      y = "Density",
      colour = "Group"
    ) +
    theme_classic()

  save_plot(p3, "02_expression_distributions",
            paste0(prefix, "_03_expression_density"),
            width = 9, height = 8)

  # Figure 4: Percentage of zero-expression values
  zero_summary <- data.frame(
    sample = colnames(expr),
    zero_percent = colMeans(expr == 0) * 100
  ) |>
    left_join(metadata, by = "sample")

  write.csv(
    zero_summary,
    file.path(project_dir, "results", "tables",
              paste0(prefix, "_zero_expression_summary.csv")),
    row.names = FALSE
  )

  p4 <- ggplot(zero_summary,
               aes(x = sample, y = zero_percent, fill = group)) +
    geom_col() +
    labs(
      title = paste(tissue_name, "- zero-expression values"),
      x = "Sample",
      y = "Genes with value zero (%)",
      fill = "Group"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(p4, "01_data_overview",
            paste0(prefix, "_04_zero_expression"))

  # Figure 5: Sample correlation heatmap
  # Spearman correlation is used for exploratory comparison.
  cor_matrix <- cor(
    log_expr,
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  cor_data <- as.data.frame(as.table(cor_matrix))
  names(cor_data) <- c("sample1", "sample2", "correlation")

  p5 <- ggplot(cor_data,
               aes(x = sample1, y = sample2, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "navy", mid = "white", high = "firebrick",
      midpoint = 0.5, limits = c(-1, 1)
    ) +
    coord_equal() +
    labs(
      title = paste(tissue_name, "- sample correlation"),
      x = NULL, y = NULL, fill = "Spearman\ncorrelation"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  save_plot(p5, "03_sample_relationships",
            paste0(prefix, "_05_sample_correlation"))

  # Figure 6: PCA using the 500 most variable genes
  # PCA is exploratory; no differential-expression testing here.
  gene_variance <- apply(log_expr, 1, var)
  gene_variance[!is.finite(gene_variance)] <- -Inf

  valid_genes <- which(is.finite(gene_variance) &
                         gene_variance > 0)

  if (length(valid_genes) < 3) {
    warning("Too few variable genes for PCA: ", tissue_name)
    return(invisible(NULL))
  }

  top_genes <- valid_genes[
    order(gene_variance[valid_genes], decreasing = TRUE)
  ]

  top_genes <- head(top_genes, 500)

  pca <- prcomp(
    t(log_expr[top_genes, , drop = FALSE]),
    center = TRUE,
    scale. = FALSE
  )

  variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  pca_data <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  ) |>
    left_join(metadata, by = "sample")

  p6 <- ggplot(pca_data,
               aes(x = PC1, y = PC2, colour = group, label = sample)) +
    geom_point(size = 3.5) +
    geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    labs(
      title = paste(tissue_name, "- principal component analysis"),
      subtitle = "PCA of the 500 most variable genes",
      x = sprintf("PC1 (%.1f%%)", variance_percent[1]),
      y = sprintf("PC2 (%.1f%%)", variance_percent[2]),
      colour = "Group"
    ) +
    theme_classic()

  save_plot(p6, "03_sample_relationships",
            paste0(prefix, "_06_PCA"))

  # Save PCA coordinates and variance explained
  write.csv(
    pca_data,
    file.path(project_dir, "results", "tables",
              paste0(prefix, "_PCA_coordinates.csv")),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      component = paste0("PC", seq_along(variance_percent)),
      variance_percent = variance_percent
    ),
    file.path(project_dir, "results", "tables",
              paste0(prefix, "_PCA_variance.csv")),
    row.names = FALSE
  )

  invisible(list(
    expression = expr,
    metadata = metadata,
    correlation = cor_matrix,
    pca = pca
  ))
}

# 6. Run both tissue analyses ---------------------------------

liver_results <- plot_tissue(liver_data, "Liver")
testis_results <- plot_tissue(testis_data, "Testis")

cat("\nExploratory figures saved to:\n", figure_dir, "\n")


project_dir <- file.path(
  path.expand("~"),
  "fish-environmental-omics"
)

figure_dir <- file.path(
  project_dir, "results", "figures"
)

table_dir <- file.path(
  project_dir, "results", "tables"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(dir.exists(figure_dir), dir.exists(table_dir))

write.csv(
  metadata_all,
  file.path(table_dir, "sample_metadata.csv"),
  row.names = FALSE
)


exists("liver_results")
exists("testis_results")

project_dir <- file.path(
  path.expand("~"),
  "fish-environmental-omics"
)

list.files(
  file.path(project_dir, "results"),
  recursive = TRUE
)

head(liver$`GENE-ID`)
head(testis$`GENE-ID`)

# Check whether expression values are whole numbers
liver_values <- as.matrix(liver[, -1])
testis_values <- as.matrix(testis[, -1])

cat("Liver: all values are whole numbers:",
    all(liver_values == round(liver_values)), "\n")

cat("Testis: all values are whole numbers:",
    all(testis_values == round(testis_values)), "\n")
# Inspect the range of values
range(liver_values, na.rm = TRUE)
range(testis_values, na.rm = TRUE)



list.files(
  file.path(
    path.expand("~"),
    "fish-environmental-omics",
    "data"
  ),
  recursive = TRUE
)



library(GEOquery)

geo_files <- getGEOSuppFiles(
  "GSE116687",
  fetch_files = FALSE
)

print(geo_files)



# Find the downloaded files
geo_paths <- list.files(
  raw_dir,
  pattern = "\\.txt\\.gz$",
  full.names = TRUE
)

print(basename(geo_paths))

# Inspect the first six rows of each file
for (f in geo_paths) {
  cat("\n\nFILE:", basename(f), "\n")

  preview <- read.delim(
    gzfile(f),
    nrows = 6,
    check.names = FALSE
  )

  print(dim(preview))
  print(names(preview))
  print(preview)
}
project_dir <- file.path(
  path.expand("~"),
  "fish-environmental-omics"
)

raw_dir <- file.path(
  project_dir,
  "data",
  "raw",
  "geo_supplementary"
)

dir.exists(raw_dir)

list.files(
  raw_dir,
  recursive = TRUE,
  full.names = TRUE
)

geo_paths <- list.files(
  raw_dir,
  pattern = "\\.txt\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)
print(basename(geo_paths))
for (f in geo_paths) {
  cat("\n\nFILE:", basename(f), "\n")


  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )

  raw_dir <- file.path(
    project_dir,
    "data",
    "raw",
    "geo_supplementary"
  )

  geo_paths <- list.files(
    raw_dir,
    pattern = "\\.txt\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )

  print(geo_paths)

f <- geo_paths[1]
cat("Inspecting:", basename(f), "\n")
preview <- read.delim(
  gzfile(f),
  nrows = 6,
  check.names = FALSE
)
print(dim(preview))
print(names(preview))
print(preview)


for (i in seq_along(geo_paths)) {

  f <- geo_paths[i]

  cat("\n==============================\n")
  cat("FILE:", basename(f), "\n")

  dat <- read.delim(
    gzfile(f),
    check.names = FALSE,
    nrows = 6
  )

  cat("Preview dimensions:", nrow(dat), "rows x",
      ncol(dat), "columns\n")

  cat("\nColumn names:\n")
  print(names(dat))

  cat("\nFirst six rows:\n")
  print(dat)
}


for (f in geo_paths) {

  original <- read.delim(
    gzfile(f),
    check.names = FALSE
  )

  cat("\nFile:", basename(f), "\n")
  cat("Dimensions:", nrow(original), "x", ncol(original), "\n")

  if (all(names(original) %in% names(liver))) {
    cat("Columns match liver:",
        identical(names(original), names(liver)), "\n")
  }

  if (all(names(original) %in% names(testis))) {
    cat("Columns match testis:",
        identical(names(original), names(testis)), "\n")
  }
}

# Compare official downloaded files with existing R objects
for (f in geo_paths) {

  original <- read.delim(
    gzfile(f),
    check.names = FALSE
  )
# Compare official downloaded files with existing R objects
for (f in geo_paths) {

  original <- read.delim(
    gzfile(f),
    check.names = FALSE
  )
  cat("\nFile:", basename(f), "\n")

  if (grepl("_MTG", basename(f))) {
    cat("Comparing with testis object:\n")
    print(all.equal(original, testis))
  }
  if (grepl("_MTL", basename(f))) {
    cat("Comparing with liver object:\n")
    print(all.equal(original, liver))
  }
}

  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )
  table_dir <- file.path(
    project_dir,
    "results",
    "tables"
  )
  dir.create(
    table_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  provenance <- data.frame(
    GEO_accession = "GSE116687",
    organism = "Oreochromis niloticus",
    exposure = "Benzo(a)pyrene (BaP)",
    tissue = c("Testis", "Liver"),
    source_file = c(
      "GSE116687_MTG.txt.gz",
      "GSE116687_MTL.txt.gz"
    ),
    genes = c(nrow(testis), nrow(liver)),
    control_samples = c(3, 4),
    BaP_samples = c(3, 4),
    expression_format = "Not yet verified",
    stringsAsFactors = FALSE
  )

  write.csv(
    provenance,
    file.path(table_dir, "data_provenance.csv"),
    row.names = FALSE
  )
  print(provenance)
  library(GEOquery)

  geo_metadata <- getGEO(
    "GSE116687",
    GSEMatrix = FALSE
  )
  Meta(geo_metadata)$summary
  Meta(geo_metadata)$overall_design
  Meta(geo_metadata)$relation
  for (f in geo_paths) {

    dat <- read.delim(
      gzfile(f),
      check.names = FALSE
    )

    values <- as.matrix(dat[, -1])

    cat("\nFile:", basename(f), "\n")
    cat("Minimum:", min(values, na.rm = TRUE), "\n")
    cat("Median:", median(values, na.rm = TRUE), "\n")
    cat("Maximum:", max(values, na.rm = TRUE), "\n")
    cat("Fractional values:", sum(values != round(values)), "\n")
  }

  library(GEOquery)

  geo_metadata <- getGEO(
    "GSE116687",
    GSEMatrix = FALSE
  )

  Meta(geo_metadata)$summary
  Meta(geo_metadata)$overall_design
  Meta(geo_metadata)$supplementary_file


  library(GEOquery)

  sample_ids <- c(
    "GSM3258483",
    "GSM3258484",
    "GSM3258485",
    "GSM3258486",
    "GSM3258487",
    "GSM3258488",
    "GSM3258489",
    "GSM3258490",
    "GSM3258491",
    "GSM3258492",
    "GSM3258493",
    "GSM3258494"
  )

  for (id in sample_ids) {

    cat("\n========================\n")
    cat("Sample:", id, "\n")

    gsm <- getGEO(id, GSEMatrix = FALSE)
    meta <- Meta(gsm)

    cat("Title:\n")
    print(meta$title)

    cat("Data processing:\n")
    print(meta$data_processing)

    cat("Characteristics:\n")
    print(meta$characteristics_ch1)
  }
  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )
  raw_dir <- file.path(
    project_dir,
    "data",
    "raw",
    "geo_supplementary"
  )
  geo_paths <- list.files(
    raw_dir,
    pattern = "\\.txt\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  for (f in geo_paths) {
    cat("\n====================================\n")
    cat("File:", basename(f), "\n")
    cat("====================================\n")

    con <- gzfile(f, open = "rt")
    lines <- readLines(con, n = 8)
    close(con)
    cat(paste(lines, collapse = "\n"), "\n")
  }

  library(GEOquery)
  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )
  table_dir <- file.path(
    project_dir,
    "results",
    "tables"
  )
  dir.create(
    table_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  sample_ids <- c(
    "GSM3258483", "GSM3258484", "GSM3258485",
    "GSM3258486", "GSM3258487", "GSM3258488",
    "GSM3258489", "GSM3258490", "GSM3258491",
    "GSM3258492", "GSM3258493", "GSM3258494"
  )
  sample_metadata <- data.frame(
    GSM_ID = character(),
    title = character(),
    tissue = character(),
    data_processing = character(),
    stringsAsFactors = FALSE
  )
  for (id in sample_ids) {
    gsm <- getGEO(id, GSEMatrix = FALSE)
    meta <- Meta(gsm)

    tissue_text <- paste(
      meta$characteristics_ch1,
      collapse = "; "
    )

    processing_text <- paste(
      meta$data_processing,
      collapse = " "
    )

    sample_metadata <- rbind(
      sample_metadata,
      data.frame(
        GSM_ID = id,
        title = paste(meta$title, collapse = "; "),
        tissue = tissue_text,
        data_processing = processing_text,
        stringsAsFactors = FALSE
      )
    )
  }
  write.csv(
    sample_metadata,
    file.path(table_dir, "GSE116687_sample_metadata.csv"),
    row.names = FALSE
  )
  cat("Sample metadata saved successfully.\n")
  cat("Rows:", nrow(sample_metadata), "\n")
  cat("File:", file.path(
    table_dir,
    "GSE116687_sample_metadata.csv"
  ), "\n")

  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )
  metadata_file <- file.path(
    project_dir,
    "results",
    "tables",
    "GSE116687_sample_metadata.csv"
  )
  # Confirm the metadata file exists
  file.exists(metadata_file)

  # Check the sample titles and tissue information
  sample_metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
  )
  sample_metadata[, c("GSM_ID", "title", "tissue")]
  # Show the project output files created so far
  list.files(
    file.path(project_dir, "results"),
    recursive = TRUE
  )

  library(GEOquery)
  gsm <- getGEO(
    "GSM3258483",
    GSEMatrix = FALSE
  )

  meta <- Meta(gsm)
  cat("Sample title:\n")
  print(meta$title)
  cat("\nSupplementary files:\n")
  print(meta$supplementary_file)

  cat("\nData processing:\n")
  print(meta$data_processing)
  cat("\nOther metadata field names:\n")
  print(names(meta))
  library(GEOquery)
  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )
  metadata_file <- file.path(
    project_dir, "results", "tables",
    "GSE116687_sample_metadata.csv"
  )
  sample_metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
  )
  missing_ids <- c("GSM3258495", "GSM3258496")
  for (id in missing_ids) {
    gsm <- getGEO(id, GSEMatrix = FALSE)
    meta <- Meta(gsm)

    sample_metadata <- rbind(
      sample_metadata,
      data.frame(
        GSM_ID = id,
        title = paste(meta$title, collapse = "; "),
        tissue = paste(
          meta$characteristics_ch1,
          collapse = "; "
        ),
        data_processing = paste(
          meta$data_processing,
          collapse = " "
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  write.csv(
    sample_metadata,
    metadata_file,
    row.names = FALSE
  )
  cat("Total samples:", nrow(sample_metadata), "\n")
  print(sample_metadata[, c("GSM_ID", "title")])

  library(GEOquery)

  for (id in c("GSM3258495", "GSM3258496")) {
    cat("\n-----------------------------\n")
    cat("Checking:", id, "\n")

    tryCatch({
      gsm <- getGEO(id, GSEMatrix = FALSE)
      meta <- Meta(gsm)

      cat("Title:\n")
      print(meta$title)

      cat("Characteristics:\n")
      print(meta$characteristics_ch1)
    }, error = function(e) {
      cat("Retrieval error:", conditionMessage(e), "\n")
    })
  }

  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )
  metadata_file <- file.path(
    project_dir, "results", "tables",
    "GSE116687_sample_metadata.csv"
  )
  sample_metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
  )
  new_samples <- data.frame(
    GSM_ID = c("GSM3258495", "GSM3258496"),
    title = c("14L-BaP rep 3", "11L-BaP rep 4"),
    tissue = c(
      "tissue: Liver; age: One year old; Sex: male",
      "tissue: Liver; age: One year old; Sex: male"
    ),
    data_processing = rep(
      sample_metadata$data_processing[1],
      2
    ),
    stringsAsFactors = FALSE
  )


  sample_metadata <- rbind(
    sample_metadata,
    new_samples
  )
  # Load the existing metadata
  sample_metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
  )
  # Define the two missing samples
  new_samples <- data.frame(
    GSM_ID = c("GSM3258495", "GSM3258496"),
    title = c("14L-BaP rep 3", "11L-BaP rep 4"),
    tissue = c(
      "tissue: Liver; age: One year old; Sex: male",
      "tissue: Liver; age: One year old; Sex: male"
    ),
    data_processing = rep(
      sample_metadata$data_processing[1],
      2
    ),
    stringsAsFactors = FALSE
  )
  # Add only samples not already present
  new_samples <- new_samples[
    !new_samples$GSM_ID %in% sample_metadata$GSM_ID,
  ]

  if (nrow(new_samples) > 0) {
    sample_metadata <- rbind(
      sample_metadata,
      new_samples
    )
  }
  # Save the updated metadata
  write.csv(
    sample_metadata,
    metadata_file,
    row.names = FALSE
  )
  # Verify the result
  cat("Total samples:", nrow(sample_metadata), "\n")
  print(sample_metadata[, c("GSM_ID", "title")])


  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )

  metadata_file <- file.path(
    project_dir, "results", "tables",
    "GSE116687_sample_metadata.csv"
  )

  # Read the current file
  sample_metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
  )

  # Add missing sample records
  extra <- data.frame(
    GSM_ID = c("GSM3258495", "GSM3258496"),
    title = c("14L-BaP rep 3", "11L-BaP rep 4"),
    tissue = c(
      "tissue: Liver; age: One year old; Sex: male",
      "tissue: Liver; age: One year old; Sex: male"
    ),
    data_processing = c(
      sample_metadata$data_processing[1],
      sample_metadata$data_processing[1]
    ),
    stringsAsFactors = FALSE
  )

  # Keep existing rows and add only IDs that are missing
  sample_metadata <- sample_metadata[
    !sample_metadata$GSM_ID %in% extra$GSM_ID,
  ]

  sample_metadata <- rbind(sample_metadata, extra)

  # Save and read the file back to verify
  write.csv(sample_metadata, metadata_file, row.names = FALSE)

  verified <- read.csv(metadata_file, stringsAsFactors = FALSE)

  cat("Rows in saved file:", nrow(verified), "\n")
  cat("Missing sample IDs found:\n")
  print(verified[verified$GSM_ID %in% c(
    "GSM3258495", "GSM3258496"
  ), c("GSM_ID", "title")])

  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )

  metadata_file <- file.path(
    project_dir, "results", "tables",
    "GSE116687_sample_metadata.csv"
  )

  # Load the existing file
  sample_metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
  )
  # Define the two missing records
  sample_1 <- sample_metadata[1, ]
  sample_2 <- sample_metadata[1, ]
  sample_1$GSM_ID <- "GSM3258495"
  sample_1$title <- "14L-BaP rep 3"
  sample_1$tissue <- "tissue: Liver; age: One year old; Sex: male"
  sample_2$GSM_ID <- "GSM3258496"
  sample_2$title <- "11L-BaP rep 4"
  sample_2$tissue <- "tissue: Liver; age: One year old; Sex: male"

  # Add records only if they are missing
  if (!"GSM3258495" %in% sample_metadata$GSM_ID) {
    sample_metadata <- rbind(sample_metadata, sample_1)
  }
  if (!"GSM3258496" %in% sample_metadata$GSM_ID) {
    sample_metadata <- rbind(sample_metadata, sample_2)
  }
  # Save the updated file
  write.csv(
    sample_metadata,
    metadata_file,
    row.names = FALSE
  )

  # Verify
  cat("Total rows:", nrow(sample_metadata), "\n")
  print(tail(sample_metadata[, c("GSM_ID", "title")], 4))

  project_dir <- file.path(path.expand("~"), "fish-environmental-omics")

  metadata_file <- file.path(
    project_dir, "results", "tables",
    "GSE116687_sample_metadata.csv"
  )

  sample_metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
  )
  # Add the missing liver BaP sample
  if (!"GSM3258495" %in% sample_metadata$GSM_ID) {
    new_sample <- sample_metadata[1, ]
    new_sample$GSM_ID <- "GSM3258495"
    new_sample$title <- "14L-BaP rep 3"
    new_sample$tissue <- "tissue: Liver; age: One year old; Sex: male"

    sample_metadata <- rbind(sample_metadata, new_sample)
  }
  # Save the corrected metadata
  write.csv(sample_metadata, metadata_file, row.names = FALSE)
  # Verify all samples
  cat("Total samples:", nrow(sample_metadata), "\n")
  print(sample_metadata[, c("GSM_ID", "title")])

  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )
  metadata_file <- file.path(
    project_dir, "results", "tables",
    "GSE116687_sample_metadata.csv"
  )
  sample_metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
  )

  cat("Total samples:", nrow(sample_metadata), "\n")
  cat("Unique sample IDs:", length(unique(sample_metadata$GSM_ID)), "\n")
  cat("Missing values in sample IDs:", sum(is.na(sample_metadata$GSM_ID)), "\n")
  cat("Duplicate sample IDs:", sum(duplicated(sample_metadata$GSM_ID)), "\n")


  print(sample_metadata[, c("GSM_ID", "title")])

  # 1. Check liver data
  cat("LIVER DATA\n")
  print(dim(liver))
  print(colnames(liver))
  print(head(liver[, 1:min(5, ncol(liver))]))

  # 2. Check testis data
  cat("\nTESTIS DATA\n")
  print(dim(testis))
  print(colnames(testis))
  print(head(testis[, 1:min(5, ncol(testis))]))
  # 3. Check whether values are integers
  liver_values <- as.matrix(liver[, -1])
  testis_values <- as.matrix(testis[, -1])

  cat(
    "\nLiver values are all integers:",
    all(abs(liver_values - round(liver_values)) < 1e-8),
    "\n"
  )

  cat(
    "Testis values are all integers:",
    all(abs(testis_values - round(testis_values)) < 1e-8),
    "\n"
  )

  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"
  )
  liver_file <- file.path(
    project_dir, "data", "raw",
    "geo_supplementary", "GSE116687_MTL.txt.gz"
  )
  testis_file <- file.path(
    project_dir, "data", "raw",
    "geo_supplementary", "GSE116687_MTG.txt.gz"
  )
  # Display the first lines of each original file
  cat("ORIGINAL LIVER FILE\n")
  print(readLines(gzfile(liver_file), n = 5))
  cat("\nORIGINAL TESTIS FILE\n")
  print(readLines(gzfile(testis_file), n = 5))

  # Summarise the numerical values in each dataset

  summarise_expression <- function(data, tissue_name) {
    values <- as.data.frame(data[, -1])
    values[] <- lapply(values, as.numeric)

    summary_table <- data.frame(
      Sample = names(values),
      Sum = sapply(values, sum, na.rm = TRUE),
      Median = sapply(values, median, na.rm = TRUE),
      Maximum = sapply(values, max, na.rm = TRUE),
      Zero_values = sapply(
        values,
        function(x) sum(x == 0, na.rm = TRUE)
      )
    )

    cat("\n", tissue_name, "SUMMARY\n")
    print(summary_table, row.names = FALSE)

    return(summary_table)
  }

  liver_summary <- summarise_expression(liver, "LIVER")
  testis_summary <- summarise_expression(testis, "TESTIS")

  # Save diagnostic summaries
  project_dir <- file.path(
    path.expand("~"),
    "fish-environmental-omics"

    )
  write.csv(
    liver_summary,
    file.path(project_dir, "results", "tables",
              "liver_expression_diagnostics.csv"),
    row.names = FALSE
  )
  write.csv(
    testis_summary,
    file.path(project_dir, "results", "tables",
              "testis_expression_diagnostics.csv"),
    row.names = FALSE

    )
