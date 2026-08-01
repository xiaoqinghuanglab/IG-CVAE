#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 160
)

arguments <- commandArgs(trailingOnly = TRUE)

mode <- if (length(arguments)) {
  toupper(arguments[[1]])
} else {
  "ALL"
}

valid_modes <- c(
  "ALL",
  "ADNI",
  "ADDNEUROMED",
  "VALIDATE"
)

if (!mode %in% valid_modes) {
  stop(
    "Mode must be ALL, ADNI, AddNeuroMed, or VALIDATE."
  )
}

root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

local_library <- file.path(
  root,
  ".R",
  "library"
)

.libPaths(
  c(
    local_library,
    .libPaths()
  )
)

required_packages <- c(
  "limma",
  "data.table",
  "ggplot2",
  "ggrepel",
  "openxlsx"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages)) {
  stop(
    "Missing R packages: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages(
  {
    library(limma)
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(openxlsx)
  }
)

manifest_path <- file.path(
  root,
  "Output",
  "Intermediate",
  "Splits",
  "model_split_manifest.csv"
)

split_config_path <- file.path(
  root,
  "Output",
  "Intermediate",
  "Splits",
  "split_config.json"
)

cohort_specifications <- list(
  ADNI = list(
    expression = file.path(
      root,
      "data",
      "ADNI",
      "gene_expression_data.csv"
    ),
    metadata = file.path(
      root,
      "data",
      "ADNI",
      "Paired_data_metadata.csv"
    ),
    metadata_id = "subject_id",
    diagnosis = "diagnosis",
    age = "age",
    sex = "sex",
    apoe = "apoe4_allele_count",
    batch = NULL,
    expected_development = 242L,
    expected_locked_test = 60L,
    covariate_description = (
      "age, sex, and APOE e4 allele count"
    )
  ),
  AddNeuroMed = list(
    expression = file.path(
      root,
      "data",
      "AdNeuroMed",
      "gene_expression_data.csv"
    ),
    metadata = file.path(
      root,
      "data",
      "AdNeuroMed",
      "gene_expression_metadata.csv"
    ),
    metadata_id = "Subject_ID",
    diagnosis = "Diagnosis",
    age = "Age",
    sex = "Sex",
    apoe = "APOE",
    batch = "Batch",
    expected_development = 280L,
    expected_locked_test = 70L,
    covariate_description = (
      paste0("age, sex, APOE e4 allele count, ", "and microarray batch")
    )
  )
)

input_paths <- c(
  manifest_path,
  split_config_path,
  unlist(
    lapply(
      cohort_specifications,
      function(specification) {
        c(
          specification$expression,
          specification$metadata
        )
      }
    ),
    use.names = FALSE
  )
)

missing_inputs <- input_paths[
  !file.exists(input_paths)
]

if (length(missing_inputs)) {
  stop(
    "Missing inputs:\n  ",
    paste(
      missing_inputs,
      collapse = "\n  "
    )
  )
}

manifest <- data.table::fread(
  manifest_path,
  data.table = FALSE,
  check.names = FALSE
)

if (
  !identical(
    unique(manifest$split_version),
    "dual_cohort_independent_v1"
  )
) {
  stop(
    paste0("Expected split version ", "dual_cohort_independent_v1.")
  )
}

if (mode == "VALIDATE") {
  for (cohort in names(cohort_specifications)) {
    specification <- cohort_specifications[[cohort]]

    development <- manifest[
      manifest$cohort == cohort
      & manifest$model == "gene_branch"
      & manifest$role == "development",
      ,
      drop = FALSE
    ]

    locked_test <- manifest[
      manifest$cohort == cohort
      & manifest$model == "gene_branch"
      & manifest$role == "locked_test",
      ,
      drop = FALSE
    ]

    if (
      nrow(development)
      != specification$expected_development
    ) {
      stop(
        cohort,
        " development count mismatch."
      )
    }

    if (
      nrow(locked_test)
      != specification$expected_locked_test
    ) {
      stop(
        cohort,
        " locked-test count mismatch."
      )
    }

    if (
      length(
        intersect(
          development$subject_id,
          locked_test$subject_id
        )
      )
    ) {
      stop(
        cohort,
        " development/test overlap."
      )
    }
  }

  cat("Dual-cohort split validation: PASS\n")
  cat("DEG packages and inputs: PASS\n")
  cat(
    paste0("Locked-test labels are not used for DEG ", "fitting or feature ranking.\n")
  )
  cat("DUAL-COHORT DEG VALIDATION: PASS\n")

  quit(
    save = "no",
    status = 0
  )
}

table_directory <- file.path(
  root,
  "Output",
  "Table",
  "Supplementary"
)

figure_directory <- file.path(
  root,
  "Output",
  "Figure",
  "Supplementary"
)

model_input_directory <- file.path(
  root,
  "Output",
  "Model_input",
  "Gene"
)

qc_directory <- file.path(
  root,
  "Output",
  "QC",
  "DEG"
)

log_directory <- file.path(
  root,
  "Output",
  "Log",
  "DEG"
)

for (
  directory in c(
    table_directory,
    figure_directory,
    model_input_directory,
    qc_directory,
    log_directory
  )
) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

managed_patterns <- c(
  "^Supplementary_Table_[1-4]_",
  "^Supplementary_Figure_1[AB]_",
  "^(ADNI|AddNeuroMed)_top500_gene_",
  "^Primary_373_gene_order\\.txt$",
  "_373_gene_"
)

for (
  directory in c(
    table_directory,
    figure_directory,
    model_input_directory
  )
) {
  existing <- list.files(
    directory,
    full.names = TRUE
  )

  remove <- existing[
    vapply(
      basename(existing),
      function(filename) {
        any(
          vapply(
            managed_patterns,
            grepl,
            logical(1),
            x = filename
          )
        )
      },
      logical(1)
    )
  ]

  if (length(remove)) {
    unlink(
      remove,
      recursive = TRUE,
      force = TRUE
    )
  }
}

clean_text <- function(value) {
  trimws(as.character(value))
}

normalize_diagnosis <- function(value) {
  value <- toupper(clean_text(value))

  result <- ifelse(
    value %in% c(
      "AD",
      "ALZHEIMER'S DISEASE",
      "ALZHEIMERS DISEASE"
    ),
    "AD",
    ifelse(
      value %in% c(
        "CN",
        "CONTROL",
        "CTL",
        "COGNITIVELY NORMAL"
      ),
      "CN",
      NA_character_
    )
  )

  if (anyNA(result)) {
    stop(
      "Unrecognized diagnosis: ",
      paste(
        unique(value[is.na(result)]),
        collapse = ", "
      )
    )
  }

  factor(
    result,
    levels = c("CN", "AD")
  )
}

normalize_sex <- function(value) {
  value <- toupper(clean_text(value))

  result <- ifelse(
    value %in% c("F", "FEMALE"),
    "Female",
    ifelse(
      value %in% c("M", "MALE"),
      "Male",
      NA_character_
    )
  )

  factor(result)
}

parse_apoe4 <- function(value) {
  output <- rep(
    NA_real_,
    length(value)
  )

  for (index in seq_along(value)) {
    if (is.na(value[[index]])) {
      next
    }

    item <- toupper(
      trimws(
        as.character(value[[index]])
      )
    )

    if (
      !nzchar(item)
      || item %in% c(
        "NA",
        "N/A",
        "MISSING",
        "."
      )
    ) {
      next
    }

    if (item %in% c("0", "1", "2")) {
      output[[index]] <- as.numeric(item)
      next
    }

    cleaned <- gsub(
      "EPSILON|APOE|E",
      "",
      item
    )

    alleles <- regmatches(
      cleaned,
      gregexpr(
        "[234]",
        cleaned,
        perl = TRUE
      )
    )[[1]]

    if (length(alleles) == 2L) {
      output[[index]] <- sum(
        alleles == "4"
      )
    }
  }

  output
}

style_workbook <- function(
  workbook,
  sheet,
  rows,
  columns
) {
  header_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 10,
    textDecoration = "bold",
    fgFill = "#D9EAF7",
    border = "Bottom",
    borderColour = "#7F7F7F"
  )

  body_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9
  )

  openxlsx::addStyle(
    workbook,
    sheet,
    header_style,
    rows = 1,
    cols = seq_len(columns),
    gridExpand = TRUE
  )

  if (rows > 1L) {
    openxlsx::addStyle(
      workbook,
      sheet,
      body_style,
      rows = 2:rows,
      cols = seq_len(columns),
      gridExpand = TRUE
    )
  }

  openxlsx::freezePane(
    workbook,
    sheet,
    firstRow = TRUE
  )

  openxlsx::addFilter(
    workbook,
    sheet,
    rows = 1,
    cols = seq_len(columns)
  )

  openxlsx::setColWidths(
    workbook,
    sheet,
    cols = seq_len(columns),
    widths = "auto"
  )
}

write_supplementary_table <- function(
  path,
  title,
  description,
  data
) {
  workbook <- openxlsx::createWorkbook(
    creator = "Dual-cohort XIG-CVAE DEG pipeline"
  )

  openxlsx::addWorksheet(
    workbook,
    "README"
  )

  readme <- data.frame(
    Item = c(
      "Title",
      "Analysis population",
      "Ranking",
      "Multiple testing",
      "Locked-test policy"
    ),
    Value = c(
      title,
      description,
      (
        paste0("Genes ranked by moderated limma raw ", "P value for AD minus CN")
      ),
      (
        "Benjamini-Hochberg false-discovery rate"
      ),
      (
        paste0("Locked-test subjects were excluded ", "from model fitting and gene ranking")
      )
    )
  )

  openxlsx::writeData(
    workbook,
    "README",
    readme
  )

  style_workbook(
    workbook,
    "README",
    nrow(readme) + 1L,
    ncol(readme)
  )

  openxlsx::addWorksheet(
    workbook,
    "Data"
  )

  openxlsx::writeData(
    workbook,
    "Data",
    data
  )

  style_workbook(
    workbook,
    "Data",
    nrow(data) + 1L,
    ncol(data)
  )

  openxlsx::saveWorkbook(
    workbook,
    path,
    overwrite = TRUE
  )
}

make_volcano <- function(
  results,
  cohort,
  output_prefix,
  top500_cutoff
) {
  plot_data <- results

  plot_data$minus_log10_P <- -log10(
    pmax(
      plot_data$P.Value,
      .Machine$double.xmin
    )
  )

  plot_data$Category <- "Other genes"

  plot_data$Category[
    plot_data$adj.P.Val <= 0.05
  ] <- "BH FDR <= 0.05"

  plot_data$Category[
    plot_data$In_top500
  ] <- "Top 500"

  plot_data$Category <- factor(
    plot_data$Category,
    levels = c(
      "Other genes",
      "BH FDR <= 0.05",
      "Top 500"
    )
  )

  label_data <- head(
    plot_data[
      order(plot_data$P.Value),
    ],
    10L
  )

  fdr_raw_boundary <- if (
    any(plot_data$adj.P.Val <= 0.05)
  ) {
    max(
      plot_data$P.Value[
        plot_data$adj.P.Val <= 0.05
      ],
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  figure <- ggplot(
    plot_data,
    aes(
      x = logFC,
      y = minus_log10_P,
      colour = Category
    )
  ) +
    geom_point(
      alpha = 0.68,
      size = 1.25
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      colour = "#333333",
      linewidth = 0.55
    ) +
    geom_hline(
      yintercept = -log10(top500_cutoff),
      linetype = "dotdash",
      colour = "#009E73",
      linewidth = 0.65
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(label = Gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.2,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = c(
        "Other genes" = "#BDBDBD",
        "BH FDR <= 0.05" = "#0072B2",
        "Top 500" = "#D55E00"
      )
    ) +
    labs(
      title = paste0(
        cohort,
        " development-only differential expression"
      ),
      subtitle = (
        "AD minus CN; covariate-adjusted limma model"
      ),
      x = "log2 fold change",
      y = expression(-log[10](P)),
      colour = NULL
    ) +
    theme_classic(
      base_family = "sans",
      base_size = 11
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )

  if (is.finite(fdr_raw_boundary)) {
    figure <- figure +
      geom_hline(
        yintercept = -log10(fdr_raw_boundary),
        linetype = "longdash",
        colour = "#0072B2",
        linewidth = 0.6
      )
  }

  ggsave(
    paste0(output_prefix, ".png"),
    figure,
    width = 7.2,
    height = 5.7,
    dpi = 400,
    bg = "white"
  )

  ggsave(
    paste0(output_prefix, ".pdf"),
    figure,
    width = 7.2,
    height = 5.7,
    device = "pdf",
    bg = "white"
  )
}

run_cohort <- function(
  cohort,
  table_numbers,
  figure_letter
) {
  specification <- cohort_specifications[[cohort]]

  cat(
    "\n========== ",
    toupper(cohort),
    " DEVELOPMENT-ONLY LIMMA DEG ==========\n",
    sep = ""
  )

  expression_data <- data.table::fread(
    specification$expression,
    data.table = FALSE,
    check.names = FALSE
  )

  metadata <- data.table::fread(
    specification$metadata,
    data.table = FALSE,
    check.names = FALSE
  )

  expression_id <- names(expression_data)[[1]]

  expression_data[[expression_id]] <- clean_text(
    expression_data[[expression_id]]
  )

  metadata[[specification$metadata_id]] <- clean_text(
    metadata[[specification$metadata_id]]
  )

  development_split <- manifest[
    manifest$cohort == cohort
    & manifest$model == "gene_branch"
    & manifest$role == "development",
    ,
    drop = FALSE
  ]

  locked_split <- manifest[
    manifest$cohort == cohort
    & manifest$model == "gene_branch"
    & manifest$role == "locked_test",
    ,
    drop = FALSE
  ]

  if (
    nrow(development_split)
    != specification$expected_development
  ) {
    stop(
      cohort,
      " development size mismatch."
    )
  }

  if (
    nrow(locked_split)
    != specification$expected_locked_test
  ) {
    stop(
      cohort,
      " locked-test size mismatch."
    )
  }

  development_ids <- clean_text(
    development_split$subject_id
  )

  expression_index <- match(
    development_ids,
    expression_data[[expression_id]]
  )

  metadata_index <- match(
    development_ids,
    metadata[[specification$metadata_id]]
  )

  if (anyNA(expression_index)) {
    stop(
      cohort,
      " development IDs missing from expression matrix."
    )
  }

  if (anyNA(metadata_index)) {
    stop(
      cohort,
      " development IDs missing from metadata."
    )
  }

  development_expression <- expression_data[
    expression_index,
    ,
    drop = FALSE
  ]

  development_metadata <- metadata[
    metadata_index,
    ,
    drop = FALSE
  ]

  analysis <- data.frame(
    Subject_ID = development_ids,
    Diagnosis = normalize_diagnosis(
      development_metadata[[specification$diagnosis]]
    ),
    Age = suppressWarnings(
      as.numeric(
        development_metadata[[specification$age]]
      )
    ),
    Sex = normalize_sex(
      development_metadata[[specification$sex]]
    ),
    APOE4 = parse_apoe4(
      development_metadata[[specification$apoe]]
    )
  )

  if (!is.null(specification$batch)) {
    analysis$Batch <- factor(
      clean_text(
        development_metadata[[specification$batch]]
      )
    )
  }

  complete <- complete.cases(analysis)

  analysis_complete <- analysis[
    complete,
    ,
    drop = FALSE
  ]

  expression_complete <- development_expression[
    complete,
    ,
    drop = FALSE
  ]

  analysis_complete$Diagnosis <- droplevels(
    analysis_complete$Diagnosis
  )

  analysis_complete$Sex <- droplevels(
    analysis_complete$Sex
  )

  if ("Batch" %in% names(analysis_complete)) {
    analysis_complete$Batch <- droplevels(
      analysis_complete$Batch
    )
  }

  if (
    nlevels(analysis_complete$Diagnosis) != 2L
  ) {
    stop(
      cohort,
      " diagnosis has fewer than two levels."
    )
  }

  gene_names <- names(expression_complete)[-1]

  gene_matrix <- as.matrix(
    expression_complete[
      ,
      gene_names,
      drop = FALSE
    ]
  )

  storage.mode(gene_matrix) <- "double"

  if (any(!is.finite(gene_matrix))) {
    stop(
      cohort,
      " contains nonfinite expression values."
    )
  }

  design_formula <- if (cohort == "ADNI") {
    ~ Diagnosis + Age + Sex + APOE4
  } else {
    ~ Diagnosis + Age + Sex + APOE4 + Batch
  }

  design <- model.matrix(
    design_formula,
    data = analysis_complete
  )

  coefficient <- "DiagnosisAD"

  if (!coefficient %in% colnames(design)) {
    stop(
      "DiagnosisAD coefficient is absent for ",
      cohort
    )
  }

  fit <- limma::lmFit(
    t(gene_matrix),
    design
  )

  fit <- limma::eBayes(
    fit,
    trend = TRUE,
    robust = TRUE
  )

  results <- limma::topTable(
    fit,
    coef = coefficient,
    number = Inf,
    sort.by = "P"
  )

  results$Gene <- rownames(results)
  rownames(results) <- NULL

  results <- results[
    ,
    c(
      "Gene",
      setdiff(
        names(results),
        "Gene"
      )
    ),
    drop = FALSE
  ]

  results$rank_by_PValue <- seq_len(
    nrow(results)
  )

  results$Direction <- ifelse(
    results$logFC > 0,
    "Higher in AD",
    ifelse(
      results$logFC < 0,
      "Lower in AD",
      "No change"
    )
  )

  results$nominal_P_le_0.05 <- (
    results$P.Value <= 0.05
  )

  results$BH_FDR_le_0.05 <- (
    results$adj.P.Val <= 0.05
  )

  results$In_top500 <- (
    results$rank_by_PValue <= 500L
  )

  top500 <- results[
    seq_len(500L),
    ,
    drop = FALSE
  ]

  top500_cutoff <- max(
    top500$P.Value
  )

  full_path <- file.path(
    table_directory,
    sprintf(
      "Supplementary_Table_%d_%s_complete_DEG.xlsx",
      table_numbers[[1]],
      cohort
    )
  )

  top500_path <- file.path(
    table_directory,
    sprintf(
      "Supplementary_Table_%d_%s_top500_gene_panel.xlsx",
      table_numbers[[2]],
      cohort
    )
  )

  population_description <- sprintf(
    paste0(
      "%s gene-branch development subjects only; ",
      "%d eligible, %d complete-case analyzed, ",
      "%d locked-test subjects excluded"
    ),
    cohort,
    nrow(development_split),
    nrow(analysis_complete),
    nrow(locked_split)
  )

  write_supplementary_table(
    full_path,
    paste0(
      cohort,
      " complete differential-expression results"
    ),
    population_description,
    results
  )

  write_supplementary_table(
    top500_path,
    paste0(
      cohort,
      " independently selected top-500 gene panel"
    ),
    population_description,
    top500
  )

  full_expression_ids <- clean_text(
    expression_data[[expression_id]]
  )

  requested_ids <- c(
    clean_text(development_split$subject_id),
    clean_text(locked_split$subject_id)
  )

  model_index <- match(
    requested_ids,
    full_expression_ids
  )

  if (anyNA(model_index)) {
    stop(
      cohort,
      " model IDs missing from expression matrix."
    )
  }

  model_input <- expression_data[
    model_index,
    c(
      expression_id,
      top500$Gene
    ),
    drop = FALSE
  ]

  names(model_input)[[1]] <- "Subject_ID"

  model_input_path <- file.path(
    model_input_directory,
    paste0(
      cohort,
      "_top500_gene_model_input.csv"
    )
  )

  order_path <- file.path(
    model_input_directory,
    paste0(
      cohort,
      "_top500_gene_order.txt"
    )
  )

  data.table::fwrite(
    model_input,
    model_input_path
  )

  writeLines(
    top500$Gene,
    order_path
  )

  design_output <- data.frame(
    Subject_ID = analysis_complete$Subject_ID,
    design,
    check.names = FALSE
  )

  data.table::fwrite(
    design_output,
    file.path(
      qc_directory,
      paste0(
        cohort,
        "_DEG_design_matrix.csv"
      )
    )
  )

  figure_prefix <- file.path(
    figure_directory,
    paste0(
      "Supplementary_Figure_1",
      figure_letter,
      "_",
      cohort,
      "_volcano"
    )
  )

  make_volcano(
    results,
    cohort,
    figure_prefix,
    top500_cutoff
  )

  nominal_count <- sum(
    results$P.Value <= 0.05
  )

  fdr_count <- sum(
    results$adj.P.Val <= 0.05
  )

  ad_count <- sum(
    analysis_complete$Diagnosis == "AD"
  )

  cn_count <- sum(
    analysis_complete$Diagnosis == "CN"
  )

  qc_lines <- c(
    paste0("cohort: ", cohort),
    "contrast: AD - CN",
    paste0(
      "split_version: ",
      unique(manifest$split_version)
    ),
    paste0(
      "development_eligible: ",
      nrow(development_split)
    ),
    paste0(
      "development_analyzed: ",
      nrow(analysis_complete)
    ),
    paste0(
      "development_excluded_missing_covariates: ",
      sum(!complete)
    ),
    paste0(
      "development_AD: ",
      ad_count
    ),
    paste0(
      "development_CN: ",
      cn_count
    ),
    paste0(
      "locked_test_excluded_from_DEG: ",
      nrow(locked_split)
    ),
    paste0(
      "genes_tested: ",
      nrow(results)
    ),
    paste0(
      "covariates: ",
      specification$covariate_description
    ),
    paste0(
      "model_formula: ",
      paste(
        deparse(design_formula),
        collapse = ""
      )
    ),
    paste0(
      "nominal_P_le_0.05: ",
      nominal_count
    ),
    paste0(
      "BH_FDR_le_0.05: ",
      fdr_count
    ),
    paste0(
      "top500_raw_P_cutoff: ",
      format(
        top500_cutoff,
        scientific = TRUE,
        digits = 8
      )
    ),
    "locked_test_labels_used: false",
    paste0(
      "model_input: ",
      normalizePath(
        model_input_path,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )

  writeLines(
    qc_lines,
    file.path(
      qc_directory,
      paste0(
        cohort,
        "_DEG_QC.txt"
      )
    )
  )

  cat("Cohort: ", cohort, "\n", sep = "")

  cat(
    "Development eligible: ",
    nrow(development_split),
    "\n",
    sep = ""
  )

  cat(
    "Development analyzed: ",
    nrow(analysis_complete),
    " (AD=",
    ad_count,
    ", CN=",
    cn_count,
    ")\n",
    sep = ""
  )

  cat(
    "Excluded for missing covariates: ",
    sum(!complete),
    "\n",
    sep = ""
  )

  cat(
    "Locked test excluded from DEG: ",
    nrow(locked_split),
    "\n",
    sep = ""
  )

  cat(
    "Genes tested: ",
    nrow(results),
    "\n",
    sep = ""
  )

  cat(
    "Covariates: ",
    specification$covariate_description,
    "\n",
    sep = ""
  )

  cat(
    "Nominal P <= 0.05: ",
    nominal_count,
    "\n",
    sep = ""
  )

  cat(
    "BH FDR <= 0.05: ",
    fdr_count,
    "\n",
    sep = ""
  )

  cat(
    "Top-500 raw-P cutoff: ",
    format(
      top500_cutoff,
      scientific = TRUE,
      digits = 6
    ),
    "\n",
    sep = ""
  )

  cat("Top 10 genes:\n")

  print(
    top500[
      1:10,
      c(
        "rank_by_PValue",
        "Gene",
        "logFC",
        "P.Value",
        "adj.P.Val"
      )
    ],
    row.names = FALSE
  )

  created <- c(
    full_path,
    top500_path,
    paste0(figure_prefix, ".png"),
    paste0(figure_prefix, ".pdf"),
    model_input_path,
    order_path
  )

  cat(
    "Created:\n  ",
    paste(
      created,
      collapse = "\n  "
    ),
    "\n",
    sep = ""
  )

  list(
    cohort = cohort,
    results = results,
    top500 = top500,
    analyzed = nrow(analysis_complete),
    eligible = nrow(development_split),
    locked = nrow(locked_split),
    nominal = nominal_count,
    fdr = fdr_count,
    cutoff = top500_cutoff
  )
}

selected_cohorts <- switch(
  mode,
  ALL = c("ADNI", "AddNeuroMed"),
  ADNI = "ADNI",
  ADDNEUROMED = "AddNeuroMed"
)

results_by_cohort <- list()

for (cohort in selected_cohorts) {
  numbering <- if (cohort == "ADNI") {
    c(1L, 2L)
  } else {
    c(3L, 4L)
  }

  figure_letter <- if (cohort == "ADNI") {
    "A"
  } else {
    "B"
  }

  results_by_cohort[[cohort]] <- run_cohort(
    cohort,
    numbering,
    figure_letter
  )
}

if (
  all(
    c(
      "ADNI",
      "AddNeuroMed"
    ) %in% names(results_by_cohort)
  )
) {
  adni_top <- results_by_cohort$ADNI$top500
  anm_top <- results_by_cohort$AddNeuroMed$top500

  overlap <- intersect(
    adni_top$Gene,
    anm_top$Gene
  )

  overlap_table <- merge(
    adni_top[
      ,
      c(
        "Gene",
        "rank_by_PValue",
        "logFC",
        "P.Value",
        "adj.P.Val"
      )
    ],
    anm_top[
      ,
      c(
        "Gene",
        "rank_by_PValue",
        "logFC",
        "P.Value",
        "adj.P.Val"
      )
    ],
    by = "Gene",
    suffixes = c(
      "_ADNI",
      "_AddNeuroMed"
    )
  )

  overlap_table$Direction_concordant <- (
    sign(overlap_table$logFC_ADNI)
    == sign(overlap_table$logFC_AddNeuroMed)
  )

  data.table::fwrite(
    overlap_table,
    file.path(
      qc_directory,
      "Independent_top500_overlap_QC.csv"
    )
  )

  cat(
    paste0("\n========== INDEPENDENT TOP-500 ", "PANEL COMPARISON ==========\n")
  )

  cat("ADNI panel: 500 genes\n")
  cat("AddNeuroMed panel: 500 genes\n")

  cat(
    "Top-500 overlap: ",
    length(overlap),
    " genes\n",
    sep = ""
  )

  cat(
    "Direction-concordant overlap: ",
    sum(overlap_table$Direction_concordant),
    "/",
    nrow(overlap_table),
    "\n",
    sep = ""
  )

  cat(
    paste0("Overlap is descriptive only; panels remain ", "cohort-specific.\n")
  )
}

summary_lines <- c(
  "DUAL-COHORT DEVELOPMENT-ONLY DEG SUMMARY",
  paste0(
    "split_version: ",
    unique(manifest$split_version)
  ),
  paste0("mode: ", mode),
  (
    paste0("feature_panels: independent cohort-specific ", "top 500 by limma raw P value")
  ),
  (
    "locked_test_used_for_fitting_or_ranking: false"
  )
)

for (cohort in names(results_by_cohort)) {
  item <- results_by_cohort[[cohort]]

  summary_lines <- c(
    summary_lines,
    paste0(
      cohort,
      ": eligible=",
      item$eligible,
      "; analyzed=",
      item$analyzed,
      "; locked_test_excluded=",
      item$locked,
      "; nominal_P_le_0.05=",
      item$nominal,
      "; BH_FDR_le_0.05=",
      item$fdr,
      "; top500_cutoff=",
      format(
        item$cutoff,
        scientific = TRUE,
        digits = 8
      )
    )
  )
}

writeLines(
  summary_lines,
  file.path(
    qc_directory,
    "Dual_cohort_DEG_summary.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    log_directory,
    "DEG_sessionInfo.txt"
  )
)

cat(
  "\n========== FINAL NUMBERED DEG OUTPUTS ==========\n"
)

final_paths <- sort(
  c(
    list.files(
      table_directory,
      pattern = "^Supplementary_Table_[1-4]_",
      full.names = TRUE
    ),
    list.files(
      figure_directory,
      pattern = "^Supplementary_Figure_1[AB]_",
      full.names = TRUE
    ),
    list.files(
      model_input_directory,
      pattern = "top500_gene",
      full.names = TRUE
    )
  )
)

for (path in final_paths) {
  cat(
    path,
    "\t",
    file.info(path)$size,
    " bytes\n",
    sep = ""
  )
}

cat(
  paste0("Locked-test subjects were not used for ", "DEG fitting or ranking.\n")
)

cat("DUAL-COHORT DEG PIPELINE: PASS\n")
