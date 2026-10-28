#!/usr/bin/env Rscript
# Exploratory inference on continuous predicted KO abundances, using base R.
# This is a compositional log-ratio analysis, not a count/RNA-seq model.

ko_read_counts <- function(path) {
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                    colClasses = 'character', quote = '', comment.char = '', na.strings = NULL)
  if (ncol(tab) < 3L || nrow(tab) < 1L || !names(tab)[1] %in% c('FeatureID', 'feature_id', 'KO'))
    stop('Expected a KO-by-sample TSV with FeatureID, feature_id or KO as the first column')
  ids <- tab[[1]]; samples <- names(tab)[-1]
  if (any(!grepl('^K[0-9]{5}$', ids)) || anyDuplicated(ids))
    stop('KO identifiers must be unique canonical identifiers such as K00001')
  if (any(!nzchar(samples)) || anyDuplicated(samples) || any(grepl('[\t\r\n]', samples)))
    stop('Duplicate, empty or invalid count SampleID')
  counts <- matrix(suppressWarnings(as.numeric(as.matrix(tab[-1]))), nrow = nrow(tab),
                   dimnames = list(ids, samples))
  if (any(!is.finite(counts)) || any(counts < 0))
    stop('Predicted KO abundances must be finite nonnegative numbers; fractional values are allowed')
  if (any(!is.finite(colSums(counts))) || any(colSums(counts) <= 0))
    stop('Every sample must have a positive finite total KO abundance')
  counts
}

ko_read_metadata <- function(path, ids) {
  meta <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                   colClasses = 'character', na.strings = NULL)
  if (anyDuplicated(names(meta)) || !all(c('SampleID', 'Group') %in% names(meta)))
    stop('Metadata CSV must have unique columns including SampleID and Group')
  if (any(!nzchar(meta$SampleID)) || anyDuplicated(meta$SampleID) || !setequal(meta$SampleID, ids))
    stop('Metadata SampleID must match the exact count sample set without duplicates or omissions')
  meta <- meta[match(ids, meta$SampleID), , drop = FALSE]
  if (any(!nzchar(trimws(meta$Group)))) stop('Every sample requires a nonempty Group')
  if (any(grepl('[\t\r\n]', meta$Group))) stop('Group labels cannot contain tabs or line breaks')
  repeated_fields <- names(meta)[tolower(gsub('[^[:alnum:]]', '', names(meta))) %in% c('animalid', 'subjectid')]
  for (field in repeated_fields) {
    value <- meta[[field]]
    observed <- !toupper(trimws(value)) %in% c('', 'NA', 'N/A', 'NULL', 'NAN')
    if (anyDuplicated(value[observed]))
      stop('Independent samples are required: repeated ', field, ' values need a repeated-measures analysis')
  }
  rownames(meta) <- meta$SampleID
  meta
}

ko_covariates <- function(requested, meta) {
  if (is.null(requested) || !length(requested)) return(list())
  if (is.character(requested)) requested <- as.list(requested)
  result <- list()
  for (item in requested) {
    name <- if (is.character(item) && length(item) == 1L) item else item$name
    type <- if (is.list(item)) item$type else NULL
    if (is.null(name) || length(name) != 1L || !is.character(name) ||
        !name %in% names(meta) || name %in% c('SampleID', 'Group') || name %in% names(result))
      stop('Covariates must name distinct metadata columns other than SampleID and Group')
    raw <- meta[[name]]
    if (any(toupper(trimws(raw)) %in% c('', 'NA', 'N/A', 'NULL', 'NAN')))
      stop('Missing selected covariate values in ', name, '; supply a complete prespecified sample cohort')
    numeric <- suppressWarnings(as.numeric(raw))
    if (is.null(type)) type <- if (all(is.finite(numeric))) 'numeric' else 'categorical'
    if (length(type) != 1L || !type %in% c('numeric', 'categorical'))
      stop('Covariate type must be numeric or categorical')
    if (type == 'numeric' && any(!is.finite(numeric))) stop('Nonfinite or nonnumeric covariate: ', name)
    result[[name]] <- list(name = name, type = type, alias = sprintf('V%03d', length(result) + 1L),
                           value = if (type == 'numeric') numeric else factor(raw, levels = unique(raw)),
                           levels = if (type == 'categorical') as.list(unique(raw)) else NULL)
  }
  result
}

ko_write_table <- function(value, path) {
  write.table(value, path, sep = '\t', row.names = FALSE, quote = FALSE, na = 'NA')
}

ko_write_matrix <- function(value, path) {
  ko_write_table(data.frame(KO = rownames(value), value, check.names = FALSE), path)
}

ko_fit <- function(clr, retained, counts, model, formula, comparisons, contrasts, global_reason) {
  rows <- vector('list', nrow(counts) * nrow(comparisons)); position <- 0L
  for (i in seq_len(nrow(counts))) {
    reason <- ''; fit <- NULL; covariance <- NULL
    if (!retained[i]) {
      reason <- if (all(counts[i, ] == 0)) 'all_zero' else 'below_prevalence'
    } else if (nzchar(global_reason)) {
      reason <- global_reason
    } else {
      y <- as.numeric(clr[rownames(counts)[i], ])
      tolerance <- 1e-12 * max(1, max(abs(y)))
      if (diff(range(y)) <= tolerance) {
        reason <- 'constant_response'
      } else {
        model$response <- y
        fit <- stats::lm(formula, data = model, na.action = stats::na.fail,
                         contrasts = list(Group = stats::contr.treatment(nlevels(model$Group))))
        if (fit$rank < length(coef(fit))) reason <- 'rank_deficient_design'
        else if (df.residual(fit) <= 0) reason <- 'insufficient_residual_df'
        else if (sqrt(sum(residuals(fit)^2)) <= 1e-12 * max(1, sqrt(sum(y^2))))
          reason <- 'zero_residual_variance'
        else covariance <- stats::vcov(fit)
      }
    }
    for (j in seq_len(nrow(comparisons))) {
      effect <- se <- p <- NA_real_; contrast_reason <- reason
      if (!nzchar(contrast_reason)) {
        contrast <- contrasts[[j]]
        effect <- sum(contrast * coef(fit))
        variance <- as.numeric(t(contrast) %*% covariance %*% contrast)
        if (!is.finite(effect) || !is.finite(variance) || variance <= 0) {
          effect <- NA_real_; contrast_reason <- 'nonfinite_or_zero_contrast_variance'
        } else {
          se <- sqrt(variance)
          p <- 2 * stats::pt(-abs(effect / se), df = df.residual(fit))
          if (!is.finite(p)) { effect <- se <- p <- NA_real_; contrast_reason <- 'nonfinite_test' }
        }
      }
      position <- position + 1L
      rows[[position]] <- data.frame(KO = rownames(counts)[i], comparison = comparisons$id[j],
        numerator = comparisons$numerator[j], denominator = comparisons$denominator[j],
        effect = effect, SE = se, P = p, q = NA_real_, N = ncol(counts),
        Status = if (!retained[i]) 'excluded' else if (nzchar(contrast_reason)) 'not_estimable' else 'ok',
        Reason = contrast_reason, stringsAsFactors = FALSE)
    }
  }
  result <- do.call(rbind, rows)
  valid <- is.finite(result$P)
  result$q[valid] <- stats::p.adjust(result$P[valid], method = 'BH')
  result
}

run_ko_differential <- function(opts, run_dir) {
  if (!requireNamespace('jsonlite', quietly = TRUE)) stop('jsonlite is required')
  if (is.null(opts$counts_path) || is.null(opts$metadata_path)) stop('counts_path and metadata_path are required')
  if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE)
  destination <- file.path(run_dir, '02_statistics')
  lock <- file.path(run_dir, '.ko_statistics.lock')
  if (file.exists(destination)) stop('Saved 02_statistics is immutable; use a new analysis run')
  if (!dir.create(lock, showWarnings = FALSE)) stop('KO statistics are already running, or a previous lock remains')
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  staging <- tempfile('.ko_statistics-', tmpdir = run_dir)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)

  prevalence <- if (is.null(opts$feature_prevalence)) .1 else opts$feature_prevalence
  pseudo <- if (is.null(opts$pseudocount)) 1e-6 else opts$pseudocount
  if (!is.numeric(prevalence) || length(prevalence) != 1L || !is.finite(prevalence) || prevalence < 0 || prevalence > 1)
    stop('feature_prevalence must be a finite number between 0 and 1')
  if (!is.numeric(pseudo) || length(pseudo) != 1L || !is.finite(pseudo) || pseudo <= 0)
    stop('pseudocount must be a finite positive value in relative-proportion units')
  counts <- ko_read_counts(opts$counts_path)
  meta <- ko_read_metadata(opts$metadata_path, colnames(counts))
  groups <- unique(meta$Group)
  order <- if (is.null(opts$group_order)) groups else unlist(opts$group_order, use.names = FALSE)
  if (!is.character(order) || anyDuplicated(order) || !setequal(order, groups))
    stop('group_order must contain each actual Group exactly once')
  sizes <- table(factor(meta$Group, levels = order))
  if (length(order) < 2L || any(sizes < 2L)) stop('At least two independent samples in each of at least two groups are required')
  warnings <- character()
  if (any(sizes < 3L)) {
    warnings <- 'Some groups have only two samples; exploratory inference has limited precision.'
    warning(warnings, call. = FALSE)
  }
  covariates <- ko_covariates(opts$covariates, meta)
  aliases <- setNames(sprintf('G%03d', seq_along(order)), order)
  model <- data.frame(Group = factor(unname(aliases[meta$Group]), levels = unname(aliases)))
  constant_categorical <- FALSE
  for (covariate in covariates) {
    value <- covariate$value
    if (is.factor(value) && nlevels(value) < 2L) { constant_categorical <- TRUE; value <- rep(0, nrow(meta)) }
    model[[covariate$alias]] <- value
  }
  formula <- stats::reformulate(c('Group', vapply(covariates, function(x) x$alias, character(1))), response = 'response')
  design_formula <- stats::delete.response(stats::terms(formula))
  design <- stats::model.matrix(design_formula, data = model,
                                 contrasts.arg = list(Group = stats::contr.treatment(length(order))))
  reason <- if (constant_categorical || qr(design)$rank < ncol(design)) 'rank_deficient_design'
            else if (nrow(design) <= ncol(design)) 'insufficient_residual_df' else ''
  pairs <- utils::combn(seq_along(order), 2)
  comparisons <- data.frame(id = sprintf('C%03d', seq_len(ncol(pairs))),
                            numerator = order[pairs[2, ]], denominator = order[pairs[1, ]], stringsAsFactors = FALSE)
  comparisons$label <- paste(comparisons$numerator, 'versus', comparisons$denominator)
  newdata <- model[rep(1L, length(order)), , drop = FALSE]
  newdata$Group <- factor(unname(aliases), levels = unname(aliases))
  for (covariate in covariates) if (covariate$type == 'numeric')
    newdata[[covariate$alias]] <- mean(covariate$value)
  prediction_design <- stats::model.matrix(design_formula, data = newdata,
                                           contrasts.arg = list(Group = stats::contr.treatment(length(order))))
  contrasts <- lapply(seq_len(ncol(pairs)), function(j) prediction_design[pairs[2, j], ] - prediction_design[pairs[1, j], ])
  prevalence_actual <- rowMeans(counts > 0)
  retained <- prevalence_actual >= prevalence & rowSums(counts) > 0
  if (!any(retained)) stop('No positive KO features remain after prevalence filtering')
  filtered <- counts[retained, , drop = FALSE]
  totals <- colSums(filtered)
  if (any(!is.finite(totals)) || any(totals <= 0)) stop('Every sample must retain a positive total after feature filtering')
  relative <- sweep(filtered, 2, totals, '/')
  logs <- log2(relative + pseudo)
  clr <- sweep(logs, 2, colMeans(logs), '-')
  result <- ko_fit(clr, retained, counts, model, formula, comparisons, contrasts, reason)

  ko_write_table(result, file.path(staging, 'ko_differential.tsv'))
  ko_write_matrix(counts, file.path(staging, 'raw_abundances.tsv'))
  ko_write_matrix(sweep(counts, 2, colSums(counts), '/'), file.path(staging, 'relative_abundances.tsv'))
  ko_write_matrix(relative, file.path(staging, 'model_relative_abundances.tsv'))
  ko_write_matrix(clr, file.path(staging, 'clr_log2.tsv'))
  ko_write_table(data.frame(KO = rownames(counts), prevalence = prevalence_actual, retained = retained,
    total_abundance = rowSums(counts), reason = ifelse(retained, 'retained', ifelse(rowSums(counts) == 0, 'all_zero', 'below_prevalence'))),
    file.path(staging, 'feature_inclusion.tsv'))
  ko_write_table(data.frame(SampleID = meta$SampleID, Group = meta$Group,
    input_total = colSums(counts), model_total = totals, included = TRUE, reason = 'included_complete_metadata'),
    file.path(staging, 'sample_inclusion.tsv'))
  write.csv(meta, file.path(staging, 'metadata_aligned.csv'), row.names = FALSE, na = '')
  ko_write_table(data.frame(SampleID = meta$SampleID, design, check.names = FALSE), file.path(staging, 'model_matrix.tsv'))
  jsonlite::write_json(comparisons, file.path(staging, 'comparisons.json'), pretty = TRUE, auto_unbox = TRUE)
  covariate_report <- lapply(covariates, function(x) x[c('name', 'type', 'alias', 'levels')])
  group_report <- lapply(seq_along(order), function(i) list(label = order[i], alias = unname(aliases[i]), N = as.integer(sizes[i])))
  report <- list(method = 'exploratory_log2_CLR_linear_model',
    interpretation = 'Adjusted between-group difference in log2-CLR units; not an absolute abundance fold change. Predicted functional potential is not measured gene expression.',
    independent_samples_required = TRUE, sample_policy = 'Exact count/metadata sample-set match; no implicit exclusion or missing-value imputation',
    feature_prevalence = prevalence, prevalence_definition = 'fraction of all aligned samples with predicted abundance > 0',
    pseudocount = pseudo, pseudocount_units = 'relative proportion, added to every retained feature after TSS',
    transformation = 'Filter features by prevalence and positive total; TSS over retained KOs per sample; log2(proportion + pseudocount); subtract per-sample mean log2 over retained KOs',
    tss_denominator = 'sum of predicted abundances over retained KOs in each sample',
    clr_feature_universe = as.list(rownames(filtered)), input_feature_universe = as.list(rownames(counts)),
    sample_ids = as.list(colnames(counts)), included_sample_ids = as.list(colnames(counts)), excluded_sample_ids = list(),
    sample_mask = as.list(rep(TRUE, ncol(counts))), groups = group_report, covariates = covariate_report,
    effective_formula = paste(deparse(formula), collapse = ' '), group_aliases = as.list(aliases),
    design_columns = as.list(colnames(design)), design_rank = qr(design)$rank, residual_df = nrow(design) - qr(design)$rank,
    inference_status = if (nzchar(reason)) reason else 'per_feature_estimability_checks',
    contrast_definition = 'numerator minus denominator, covariates held equal; all unique unordered Group pairs',
    p_value = 'two-sided Student t from lm coefficient contrast and full coefficient covariance, residual degrees of freedom',
    multiplicity = list(method = 'BH', family = 'all valid KO-by-pairwise-comparison P values in this analysis, before any map or display filtering',
                        tested_hypotheses = sum(is.finite(result$P)), reference = 'Benjamini and Hochberg (1995), JRSS B 57:289-300, doi:10.1111/j.2517-6161.1995.tb02031.x'),
    excluded_feature_count = sum(!retained), retained_feature_count = sum(retained), warnings = as.list(warnings),
    sources = list(counts_path = normalizePath(opts$counts_path), metadata_path = normalizePath(opts$metadata_path),
                   counts_md5 = unname(tools::md5sum(opts$counts_path)), metadata_md5 = unname(tools::md5sum(opts$metadata_path))),
    software = list(R = R.version.string, jsonlite = as.character(utils::packageVersion('jsonlite'))))
  jsonlite::write_json(report, file.path(staging, 'transformation_report.json'), pretty = TRUE, auto_unbox = TRUE, null = 'null', digits = 16)
  jsonlite::write_json(opts, file.path(staging, 'input_options.json'), pretty = TRUE, auto_unbox = TRUE, null = 'null', digits = 16)
  writeLines(capture.output(utils::sessionInfo()), file.path(staging, 'sessionInfo.txt'))
  writeLines('complete', file.path(staging, 'COMPLETE'))
  if (file.exists(destination) || !file.rename(staging, destination)) stop('Cannot publish immutable KO statistics directory')
  invisible(report)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L || args[1] != '--run') stop('Usage: Rscript ko_differential.R --run OPTIONS.json RUN_DIR')
  run_ko_differential(jsonlite::read_json(args[2], simplifyVector = FALSE), args[3])
}
