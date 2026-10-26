#!/usr/bin/env Rscript
# Independent exploratory analyses. Statistical fits are cached once; --render
# reads only that cache and changes presentation, never samples or hypotheses.

`%||%` <- function(x, y) if (is.null(x)) y else x
adv_write <- function(x, file) {
  if (is.matrix(x)) x <- data.frame(FeatureID = rownames(x), x, check.names = FALSE)
  write.table(x, file, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA', fileEncoding = 'UTF-8')
}
adv_json <- function(x, file) jsonlite::write_json(x, file, pretty = TRUE, auto_unbox = TRUE,
                                                  null = 'null', na = 'null', digits = 12)
adv_label <- function(ids) {
  out <- vapply(strsplit(ids, ';', fixed = TRUE), function(x) {
    x <- trimws(x[nzchar(trimws(x))]); if (length(x)) tail(x, 1) else ''
  }, character(1))
  out[!nzchar(out)] <- ids[!nzchar(out)]
  # Full IDs remain the keys even when the last rank labels are identical.
  make.unique(out, sep = ' [lineage] ')
}
adv_defaults <- function(o) {
  defaults <- list(modules = c('core', 'ordination', 'association'), core_min_reads = 1,
    core_min_relative = .0001, core_prevalence = .8, feature_prevalence = .1,
    clr_pseudocount = .5, constrained_method = 'hellinger_rda', predictors = 'Group',
    covariates = character(), phenotypes = character(), association_method = 'partial_spearman', association_transform = 'clr',
    association_covariates = 'Group', permutations = 999L, seed = 42L, top_n = 20L,
    repeated_subjects = FALSE, font = 'DejaVu Serif')
  for (key in names(defaults)) if (is.null(o[[key]])) o[[key]] <- defaults[[key]]
  for (key in c('modules', 'predictors', 'covariates', 'phenotypes', 'association_covariates'))
    o[[key]] <- unlist(o[[key]], use.names = FALSE) %||% character()
  if (!all(o$modules %in% c('core', 'ordination', 'association')) || anyDuplicated(o$modules)) stop('Invalid advanced modules')
  for (key in c('core_min_relative', 'core_prevalence', 'feature_prevalence'))
    if (length(o[[key]]) != 1L || !is.finite(o[[key]]) || o[[key]] < 0 || o[[key]] > 1) stop('Invalid ', key)
  for (key in c('core_min_reads', 'permutations', 'seed', 'top_n'))
    if (length(o[[key]]) != 1L || !is.finite(o[[key]]) || o[[key]] != floor(o[[key]]) || o[[key]] < 0) stop('Invalid ', key)
  if (o$top_n < 1 || o$permutations < 1 || o$permutations > 99999 ||
      length(o$clr_pseudocount) != 1 || !is.finite(o$clr_pseudocount) || o$clr_pseudocount <= 0) stop('Invalid analysis sizing or pseudocount')
  if (!o$constrained_method %in% c('hellinger_rda', 'bray_dbrda', 'none') ||
      !o$association_method %in% c('spearman', 'partial_spearman') ||
      !o$association_transform %in% c('clr', 'relative')) stop('Invalid analysis method or association transform')
  o
}
adv_require <- function() {
  needed <- c('vegan', 'permute', 'ppcor', 'ggplot2', 'jsonlite', 'svglite', 'systemfonts')
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Missing advanced R packages: ', paste(missing, collapse = ', '))
  setNames(vapply(needed, function(p) as.character(utils::packageVersion(p)), character(1)), needed)
}
adv_counts <- function(path) {
  x <- read.delim(path, colClasses = 'character', check.names = FALSE, quote = '', comment.char = '', na.strings = NULL)
  if (ncol(x) < 3 || nrow(x) < 1 || names(x)[1] != 'FeatureID' || anyDuplicated(x[[1]]) ||
      anyDuplicated(names(x)[-1]) || any(!nzchar(x[[1]])) || any(!nzchar(names(x)[-1]))) stop('Expected unique FeatureID and sample count identifiers')
  m <- matrix(suppressWarnings(as.numeric(as.matrix(x[-1]))), nrow = nrow(x), dimnames = list(x[[1]], names(x)[-1]))
  if (any(!is.finite(m)) || any(m < 0) || any(m != floor(m))) stop('Advanced analysis requires nonnegative integer unrarefied counts')
  m
}
adv_metadata <- function(path, ids, o) {
  x <- read.csv(path, colClasses = 'character', check.names = FALSE, na.strings = NULL, fileEncoding = 'UTF-8')
  if (!all(c('SampleID', 'Group') %in% names(x)) || anyDuplicated(names(x)) || anyNA(x$SampleID) ||
      anyDuplicated(x$SampleID) || !setequal(x$SampleID, ids)) stop('Metadata SampleID must match count sample IDs exactly')
  x <- x[match(ids, x$SampleID), , drop = FALSE]; rownames(x) <- ids
  if (anyNA(x$Group) || any(!nzchar(trimws(x$Group)))) stop('Every sample requires Group metadata')
  for (key in setdiff(names(x), c('SampleID', 'Group'))) {
    x[[key]][x[[key]] %in% c('', 'NA')] <- NA_character_
    converted <- suppressWarnings(as.numeric(x[[key]]))
    if (all(is.na(x[[key]]) | is.finite(converted))) x[[key]] <- converted
  }
  requested <- unique(c(o$predictors, o$covariates, o$phenotypes, o$association_covariates))
  if (length(setdiff(requested, names(x)))) stop('Unknown metadata variables: ', paste(setdiff(requested, names(x)), collapse = ', '))
  if ('SampleID' %in% requested) stop('SampleID is an identifier and cannot be an explanatory variable')
  actual <- unique(x$Group)
  order <- c(intersect(unlist(o$group_order), actual), setdiff(actual, unlist(o$group_order)))
  x$Group <- factor(x$Group, levels = order)
  x
}
adv_design <- function(meta, variables) {
  variables <- unique(variables)
  used <- if (length(variables)) complete.cases(meta[, variables, drop = FALSE]) else rep(TRUE, nrow(meta))
  d <- data.frame(row.names = rownames(meta)[used])
  aliases <- if (length(variables)) paste0('v', seq_along(variables)) else character()
  for (i in seq_along(variables)) {
    z <- meta[[variables[i]]][used]
    if (is.character(z) || is.factor(z)) z <- droplevels(factor(z))
    if (length(unique(z)) < 2) stop('Constant metadata variable in complete cases: ', variables[i])
    d[[aliases[i]]] <- z
  }
  f <- if (length(aliases)) stats::reformulate(aliases) else stats::as.formula('~ 1')
  mm <- stats::model.matrix(f, data = d)
  rank <- qr(mm)$rank
  if (rank != ncol(mm)) stop('Rank deficient explanatory design; remove redundant variables')
  if (nrow(mm) - rank < 2L) stop('Insufficient residual degrees of freedom (minimum 2)')
  list(data = d, used = used, matrix = mm, rank = rank, residual_df = nrow(mm) - rank,
       aliases = setNames(aliases, variables), variables = variables)
}
adv_anova_table <- function(x, aliases = character(), family = '') {
  d <- as.data.frame(x); d$Term <- rownames(d); rownames(d) <- NULL
  index <- match(d$Term, unname(aliases)); hit <- !is.na(index)
  d$Term[hit] <- names(aliases)[index[hit]]
  # Canonical finite-only p-values prevent degenerate fits looking significant.
  p <- if ('Pr(>F)' %in% names(d)) d[['Pr(>F)']] else rep(NA_real_, nrow(d))
  valid <- is.finite(p) & p >= 0 & p <= 1
  if ('F' %in% names(d)) valid <- valid & is.finite(d$F)
  d$P <- ifelse(valid, p, NA_real_); d$Q <- NA_real_
  d$Q[valid] <- p.adjust(d$P[valid], method = 'BH')
  d$Status <- ifelse(valid, 'tested', 'not_tested')
  d$Family <- family
  if ('Pr(>F)' %in% names(d)) d[['Pr(>F)']] <- NULL
  d[, c('Term', setdiff(names(d), 'Term')), drop = FALSE]
}
adv_clr <- function(counts, pseudocount) {
  if (nrow(counts) < 2L) stop('CLR requires at least two prevalent features in the selected count subcomposition')
  logged <- log(counts + pseudocount)
  sweep(logged, 2, colMeans(logged), '-')
}
adv_association_definition <- function(transform, pseudocount) {
  if (transform == 'clr') paste0('CLR: log(count + ', pseudocount,
    ') minus the mean log(count + pseudocount) over all selected prevalent features per sample; ',
    'taxon log-ratio to the geometric mean of this count subcomposition; not absolute abundance') else
    'Relative: selected feature counts divided by full original chosen-rank sample totals; exploratory compositional association, not absolute abundance'
}
adv_fit <- function(o, run) {
  dir.create(file.path(run, '02_statistics'), recursive = TRUE, showWarnings = FALSE)
  folder <- file.path(run, '02_statistics'); versions <- adv_require()
  counts <- adv_counts(o$counts_path); total <- adv_counts(o$total_counts_path %||% o$counts_path)
  if (!identical(colnames(counts), colnames(total)) || !all(rownames(counts) %in% rownames(total)) ||
      !identical(unname(counts), unname(total[rownames(counts), , drop = FALSE]))) stop('Selected counts must be an unchanged subset of the full denominator count table')
  totals <- colSums(total)
  if (any(totals <= 0)) stop('All samples require a positive full-table library size')
  meta <- adv_metadata(o$metadata_path, colnames(counts), o)
  rel <- sweep(counts, 2, totals, '/')
  labels <- setNames(adv_label(rownames(counts)), rownames(counts))
  prevalence <- rowMeans(counts > 0)
  keep <- prevalence >= o$feature_prevalence & rowSums(counts) > 0
  diagnostics <- list(); data <- list(); models <- list()
  status <- function(module, state, reason = '', n = ncol(counts), rank = NA_integer_, residual_df = NA_integer_, variables = character(), samples = colnames(counts)) {
    diagnostics[[length(diagnostics) + 1L]] <<- data.frame(Module = module, Status = state, Reason = reason,
      N = n, Rank = rank, ResidualDF = residual_df, Variables = paste(variables, collapse = ';'), SampleIDs = paste(samples, collapse = ';'), stringsAsFactors = FALSE)
  }
  perform <- function(module, fn) tryCatch(fn(), error = function(e) { status(module, 'skipped', conditionMessage(e)); NULL })
  adv_write(data.frame(SampleID = colnames(counts), Group = as.character(meta$Group), FullLibrarySize = totals,
    SelectedLibrarySize = colSums(counts), stringsAsFactors = FALSE), file.path(folder, 'sample_audit.tsv'))
  adv_write(data.frame(FeatureID = rownames(counts), Label = unname(labels), Prevalence = prevalence,
    MeanRelativeAbundance = rowMeans(rel), RetainedForMultivariateAndAssociation = keep), file.path(folder, 'feature_filter.tsv'))
  adv_write(rel, file.path(folder, 'sample_relative_abundance.tsv'))
  set.seed(o$seed)
  if ('core' %in% o$modules) perform('core', function() {
    detected <- counts >= o$core_min_reads & rel >= o$core_min_relative & counts > 0
    chunks <- lapply(levels(meta$Group), function(g) {
      cols <- which(meta$Group == g)
      data.frame(FeatureID = rownames(counts), Label = unname(labels), Group = g, N = length(cols),
        DetectedSamples = rowSums(detected[, cols, drop = FALSE]), Prevalence = rowMeans(detected[, cols, drop = FALSE]),
        MeanRelativeAbundance = rowMeans(rel[, cols, drop = FALSE]), stringsAsFactors = FALSE)
    })
    d <- do.call(rbind, chunks); d$Core <- d$Prevalence >= o$core_prevalence
    adv_write(d, file.path(folder, 'core_by_group.tsv'))
    global <- data.frame(FeatureID = rownames(counts), Label = unname(labels), N = ncol(counts),
      DetectedSamples = rowSums(detected), Prevalence = rowMeans(detected), MeanRelativeAbundance = rowMeans(rel))
    global$Core <- global$Prevalence >= o$core_prevalence
    adv_write(global, file.path(folder, 'core_global.tsv'))
    data$core <<- d; data$occupancy <<- global
    status('core', 'complete', 'Descriptive detection requires both count and full-table relative-abundance thresholds')
  })
  if ('ordination' %in% o$modules) {
    perform('clr_pca', function() {
      if (sum(keep) < 2L || ncol(counts) < 3L) stop('CLR-PCA requires at least two retained features and three samples')
      clr <- t(adv_clr(counts[keep, , drop = FALSE], o$clr_pseudocount))
      if (sum(apply(clr, 2, var)) <= .Machine$double.eps) stop('CLR matrix has no between-sample variation')
      fit <- stats::prcomp(clr, center = TRUE, scale. = FALSE)
      coords <- data.frame(SampleID = rownames(clr), Group = as.character(meta$Group),
        Axis1 = fit$x[, 1], Axis2 = if (ncol(fit$x) >= 2) fit$x[, 2] else 0, stringsAsFactors = FALSE)
      variance <- data.frame(Axis = paste0('PC', seq_along(fit$sdev)), Eigenvalue = fit$sdev^2,
        Fraction = fit$sdev^2 / sum(fit$sdev^2))
      adv_write(t(clr), file.path(folder, 'clr_transformed.tsv')); adv_write(coords, file.path(folder, 'clr_pca_coordinates.tsv'))
      adv_write(variance, file.path(folder, 'clr_pca_variance.tsv'))
      data$pca <<- list(coordinates = coords, variance = variance); models$pca <<- fit
      status('clr_pca', 'complete', 'Pseudocount added to retained raw counts; CLR of selected subcomposition; PCA centered without scaling')
    })
    perform('permanova', function() {
      if (isTRUE(o$repeated_subjects)) stop('Independent samples only: repeated subjects disable inferential tests')
      if (!any(keep)) stop('No prevalent features remain')
      des <- adv_design(meta, unique(c('Group', o$covariates)))
      m <- meta[des$used, , drop = FALSE]
      if (nlevels(droplevels(m$Group)) < 2 || any(table(droplevels(m$Group)) < 2)) stop('PERMANOVA requires at least two samples per Group')
      distance <- vegan::vegdist(t(rel[keep, des$used, drop = FALSE]), method = 'bray')
      if (any(!is.finite(distance)) || max(distance) <= .Machine$double.eps) stop('Degenerate Bray-Curtis distances')
      formula <- stats::as.formula(paste('distance ~', paste(unname(des$aliases), collapse = ' + ')), env = environment())
      fit <- vegan::adonis2(formula, data = des$data, permutations = o$permutations, by = 'margin')
      tab <- adv_anova_table(fit, des$aliases, 'PERMANOVA marginal terms')
      if (!any(tab$Status == 'tested')) stop('PERMANOVA did not produce finite test statistics')
      adv_write(tab, file.path(folder, 'permanova.tsv'))
      distance_matrix <- as.matrix(distance)
      adv_write(data.frame(SampleID = rownames(distance_matrix), distance_matrix, check.names = FALSE),
        file.path(folder, 'bray_curtis_distance.tsv'))
      models$permanova <<- fit; models$bray_distance <<- distance
      data$permanova <<- tab
      status('permanova', 'complete', 'Bray-Curtis on full-table-normalized selected prevalent taxa; marginal Group and covariate tests',
        sum(des$used), des$rank, des$residual_df, des$variables, rownames(m))
    })
    perform('permdisp', function() {
      if (isTRUE(o$repeated_subjects)) stop('Independent samples only: repeated subjects disable inferential tests')
      if (!any(keep)) stop('No prevalent features remain')
      # Use the same complete-case cohort as the adjusted PERMANOVA.
      des <- adv_design(meta, unique(c('Group', o$covariates))); m <- meta[des$used, , drop = FALSE]
      groups <- droplevels(m$Group)
      if (nlevels(groups) < 2 || any(table(groups) < 3)) stop('PERMDISP requires at least three complete-case samples per Group')
      distance <- vegan::vegdist(t(rel[keep, des$used, drop = FALSE]), method = 'bray')
      if (any(!is.finite(distance)) || max(distance) <= .Machine$double.eps) stop('Degenerate Bray-Curtis distances')
      fit <- vegan::betadisper(distance, group = groups, type = 'median', bias.adjust = TRUE)
      if (any(!is.finite(fit$distances)) || var(fit$distances) <= .Machine$double.eps) stop('Degenerate distances to group medians')
      test <- vegan::permutest(fit, permutations = o$permutations)
      tab <- adv_anova_table(test$tab, family = 'PERMDISP global Group test')
      if (!any(tab$Status == 'tested')) stop('PERMDISP did not produce finite test statistics')
      adv_write(tab, file.path(folder, 'permdisp.tsv'))
      adv_write(data.frame(SampleID = rownames(m), Group = as.character(groups), DistanceToMedian = fit$distances), file.path(folder, 'permdisp_distances.tsv'))
      models$dispersion <<- fit; models$dispersion_test <<- test; data$permdisp <<- tab
      status('permdisp', 'complete', 'Unadjusted global Group dispersion diagnostic; spatial median; small-sample bias adjustment', sum(des$used), variables = 'Group', samples = rownames(m))
    })
    perform('constrained', function() {
      if (o$constrained_method == 'none') stop('Constrained ordination disabled by selection')
      if (isTRUE(o$repeated_subjects)) stop('Independent samples only: repeated subjects disable constrained inference')
      if (!length(o$predictors)) stop('No constrained predictors selected')
      if (length(intersect(o$predictors, o$covariates))) stop('Predictors and conditioning covariates must not overlap')
      if (sum(keep) < 2) stop('Constrained ordination requires at least two retained features')
      des <- adv_design(meta, c(o$predictors, o$covariates))
      response <- t(rel[keep, des$used, drop = FALSE])
      rhs <- paste(des$aliases[o$predictors], collapse = ' + ')
      if (length(o$covariates)) rhs <- paste(rhs, '+ Condition(', paste(des$aliases[o$covariates], collapse = ' + '), ')')
      formula <- stats::as.formula(paste('response ~', rhs), env = environment())
      if (o$constrained_method == 'hellinger_rda') {
        # Relative abundances retain the original full-table denominators.
        response <- sqrt(response)
        fit <- vegan::rda(formula, data = des$data)
      } else fit <- vegan::dbrda(formula, data = des$data, distance = 'bray')
      if (is.null(fit$CCA) || fit$CCA$rank < 1 || !is.finite(fit$CCA$tot.chi) || fit$CCA$tot.chi <= .Machine$double.eps) stop('Constrained fit has zero estimable constrained inertia')
      overall <- adv_anova_table(stats::anova(fit, permutations = o$permutations), family = 'Constrained overall test')
      marginal <- adv_anova_table(stats::anova(fit, permutations = o$permutations, by = 'margin'), des$aliases, 'Constrained marginal predictors')
      adv_write(overall, file.path(folder, 'constrained_overall.tsv')); adv_write(marginal, file.path(folder, 'constrained_marginal.tsv'))
      sites <- as.matrix(vegan::scores(fit, display = 'sites', choices = seq_len(min(2L, fit$CCA$rank)), scaling = 1))
      coords <- data.frame(SampleID = rownames(sites), Group = as.character(meta$Group[des$used]), Axis1 = sites[, 1],
        Axis2 = if (ncol(sites) >= 2) sites[, 2] else 0, stringsAsFactors = FALSE)
      eig <- fit$CCA$eig
      variance <- data.frame(Axis = names(eig), Eigenvalue = unname(eig), FractionOfTotal = unname(eig) / fit$tot.chi)
      adv_write(coords, file.path(folder, 'constrained_coordinates.tsv')); adv_write(variance, file.path(folder, 'constrained_variance.tsv'))
      data$constrained <<- list(coordinates = coords, variance = variance, method = o$constrained_method,
        overall = overall, marginal = marginal, single_axis = fit$CCA$rank == 1)
      models$constrained <<- fit
      status('constrained', 'complete', paste(o$constrained_method, '; covariates partialled out; overall and marginal permutation tests; separate BH families'),
        sum(des$used), des$rank, des$residual_df, des$variables, rownames(meta)[des$used])
    })
  }
  if ('association' %in% o$modules) perform('association', function() {
    if (!length(o$phenotypes)) stop('No quantitative phenotypes selected')
    if (!any(keep)) stop('No prevalent taxa remain')
    association_input <- if (o$association_transform == 'clr') adv_clr(counts[keep, , drop = FALSE], o$clr_pseudocount) else rel[keep, , drop = FALSE]
    definition <- adv_association_definition(o$association_transform, o$clr_pseudocount)
    adv_write(association_input, file.path(folder, 'association_transformed.tsv'))
    adv_json(list(transform = o$association_transform, definition = definition,
      pseudocount = if (o$association_transform == 'clr') o$clr_pseudocount else NULL,
      feature_ids = as.list(rownames(association_input)), sample_ids = as.list(colnames(association_input))),
      file.path(folder, 'association_transform.json'))
    data$association_inputs <<- list(transform = o$association_transform, definition = definition, values = association_input)
    taxa <- rownames(counts)[keep]; rows <- list(); samples <- list()
    covariates <- if (o$association_method == 'partial_spearman') o$association_covariates else character()
    for (trait in o$phenotypes) {
      reason <- ''; des <- NULL; used <- rep(FALSE, nrow(meta)); response <- meta[[trait]]
      if (!is.numeric(response)) reason <- 'Phenotype must be quantitative numeric metadata'
      if (trait %in% covariates) reason <- 'Phenotype cannot also be an adjustment covariate'
      if (isTRUE(o$repeated_subjects)) reason <- 'Independent samples only: repeated subjects disable correlation tests'
      if (!nzchar(reason)) {
        used <- is.finite(response) & if (length(covariates)) complete.cases(meta[, covariates, drop = FALSE]) else TRUE
        if (sum(used) < 5) reason <- 'Fewer than five complete-case samples'
      }
      if (!nzchar(reason)) {
        des <- tryCatch(adv_design(meta[used, , drop = FALSE], covariates), error = function(e) { reason <<- conditionMessage(e); NULL })
        if (!is.null(des)) {
          # ppcor ranks dummy columns too; inspect that exact adjustment design.
          z <- des$matrix[, -1, drop = FALSE]
          ranked <- if (ncol(z)) apply(z, 2, rank) else matrix(numeric(), nrow = sum(used), ncol = 0)
          if (is.null(dim(ranked))) ranked <- matrix(ranked, ncol = 1)
          if (qr(cbind(1, ranked))$rank != ncol(ranked) + 1) reason <- 'Rank deficient ranked adjustment design'
          if (sum(used) - ncol(z) - 2 < 3) reason <- 'Insufficient effective degrees of freedom (minimum 3)'
        }
      }
      samples[[length(samples) + 1L]] <- data.frame(Phenotype = trait, SampleID = rownames(meta), Included = used,
        Reason = ifelse(used, if (nzchar(reason)) reason else 'complete_case', if (nzchar(reason)) reason else 'missing_phenotype_or_covariate'))
      for (taxon in taxa) {
        r <- data.frame(FeatureID = taxon, Label = labels[[taxon]], Phenotype = trait,
          Method = o$association_method, Transform = o$association_transform,
          Covariates = paste(covariates, collapse = ';'), N = sum(used),
          CovariateDF = if (is.null(des)) NA_integer_ else ncol(des$matrix) - 1L,
          ResidualDF = if (is.null(des)) NA_integer_ else sum(used) - ncol(des$matrix) - 1L,
          Rho = NA_real_, P = NA_real_, Q = NA_real_, Status = 'skipped', Reason = reason,
          SampleIDs = paste(rownames(meta)[used], collapse = ';'), stringsAsFactors = FALSE)
        if (!nzchar(reason)) {
          x <- association_input[taxon, used]; y <- response[used]
          failure <- tryCatch({
            if (length(unique(x)) < 2 || length(unique(y)) < 2) stop('Constant feature or phenotype')
            if (o$association_method == 'partial_spearman' && ncol(z)) {
              ranked_all <- cbind(rank(x), rank(y), ranked)
              if (qr(cbind(1, ranked_all))$rank != ncol(ranked_all) + 1) stop('Singular ranked feature, phenotype and adjustment covariance')
              test <- ppcor::pcor.test(x, y, z, method = 'spearman')
              rho <- test$estimate; p <- test$p.value
            } else {
              test <- suppressWarnings(stats::cor.test(x, y, method = 'spearman', exact = FALSE))
              rho <- unname(test$estimate); p <- test$p.value
            }
            if (!is.finite(rho) || !is.finite(p) || p < 0 || p > 1 || abs(rho) > 1 + 1e-8) stop('Nonfinite or invalid correlation result')
            r$Rho <- rho; r$P <- p; r$Status <- 'tested'; r$Reason <- ''; NULL
          }, error = conditionMessage)
          if (!is.null(failure)) r$Reason <- failure
        }
        rows[[length(rows) + 1L]] <- r
      }
    }
    results <- do.call(rbind, rows); good <- results$Status == 'tested'
    results$Q[good] <- p.adjust(results$P[good], method = 'BH')
    adv_write(results, file.path(folder, 'phenotype_associations.tsv'))
    adv_write(do.call(rbind, samples), file.path(folder, 'association_sample_inclusion.tsv'))
    data$association <<- results
    message <- paste('Exploratory taxa-by-phenotype correlations;', definition,
      '; one BH family across all tested taxa and traits;', sum(good), 'tested of', nrow(results))
    if (!any(good)) message <- paste(message, ';', paste(unique(results$Reason), collapse = ' | '))
    status('association', if (any(good)) 'complete' else 'skipped', message,
      variables = c(o$phenotypes, covariates))
  })
  diag <- if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame(Module = character(), Status = character(), Reason = character())
  adv_write(diag, file.path(folder, 'diagnostics.tsv'))
  cache <- list(version = 1L, options = o, data = data, models = models, diagnostics = diag,
    sample_ids = colnames(counts), groups = levels(meta$Group), labels = labels, packages = versions,
    normalization = 'Selected feature counts divided by full original chosen-rank sample totals; never renormalized after filtering',
    independent_samples_only = TRUE)
  saveRDS(cache, file.path(folder, 'advanced_cache.rds'))
  adv_json(list(source = 'advanced_analysis', exploratory = TRUE, independent_samples_only = TRUE,
    repeated_subjects = o$repeated_subjects, sample_ids = as.list(colnames(counts)), package_versions = as.list(versions),
    normalization = cache$normalization, association_transform = o$association_transform,
    association_definition = adv_association_definition(o$association_transform, o$clr_pseudocount),
    diagnostics = diag, seed = o$seed, permutations = o$permutations,
    sources = c('https://vegandevs.github.io/vegan/reference/adonis.html',
      'https://vegandevs.github.io/vegan/reference/betadisper.html', 'https://vegandevs.github.io/vegan/reference/cca.html',
      'https://vegandevs.github.io/vegan/reference/anova.cca.html', 'https://cran.r-project.org/package=ppcor')),
    file.path(folder, 'analysis_summary.json'))
  cache
}

adv_config <- function(o, profile) {
  general <- list(dpi = 200, font_family = 'original', base_size = 11, title_size = NULL,
    axis_text_size = NULL, legend_text_size = NULL, legend_position = 'right', theme = 'bw',
    background = 'white', show_caption = TRUE, formats = 'png_svg_pdf')
  specific <- list(width = NULL, height = NULL, title = '', x_label = '', y_label = '')
  specific <- c(specific, switch(profile,
    beta = list(point_size = 3, point_alpha = .9, show_ellipse = TRUE, ellipse_level = .95,
      ellipse_alpha = .15, show_sample_labels = FALSE, show_stats = TRUE),
    heatmap = list(label_size = 8, row_order = 'abundance', show_values = FALSE, top_n = NULL),
    bubble = list(max_size = 8, point_alpha = .9, show_values = FALSE, top_n = NULL)))
  requested <- o$plot_settings %||% list()
  for (key in names(general)) if (!is.null(requested$general[[key]])) general[[key]] <- requested$general[[key]]
  for (key in names(specific)) if (!is.null(requested[[profile]][[key]])) specific[[key]] <- requested[[profile]][[key]]
  if (general$legend_position == 'original') general$legend_position <- 'right'
  if (general$theme == 'original') general$theme <- 'bw'
  font <- if (general$font_family == 'original') o$font else general$font_family
  if (!font %in% systemfonts::system_fonts()$family) stop('Advanced figure font is not installed: ', font)
  list(general = general, profile = specific, font = font)
}
adv_colors <- function(o) {
  colors <- list(heatmap_low = '#2166AC', heatmap_mid = '#F7F7F7', heatmap_high = '#B2182B',
    bubble_low = '#E8F1FA', bubble_high = '#2166AC')
  for (key in names(colors)) if (!is.null(o[[key]])) colors[[key]] <- o[[key]]
  if (any(!vapply(colors, function(x) length(x) == 1 && is.character(x) &&
    !is.na(x) && grepl('^#[[:xdigit:]]{6}$', x), logical(1)))) stop('Advanced colors must use #RRGGBB')
  colors
}
adv_group_style <- function(cache, o) {
  groups <- cache$groups
  order <- c(intersect(unlist(o$group_order), groups), setdiff(groups, unlist(o$group_order)))
  if (anyDuplicated(order)) stop('Duplicate group display order')
  colors <- setNames(grDevices::hcl.colors(length(order), 'Dark 3'), order)
  fixed <- c(CON = '#D55E00', LDC = '#0072B2', HDC = '#CC79A7')
  colors[intersect(names(fixed), order)] <- fixed[intersect(names(fixed), order)]
  custom <- unlist(o$group_colors)
  if (length(custom) && (is.null(names(custom)) || anyDuplicated(names(custom)) ||
      any(!grepl('^#[[:xdigit:]]{6}$', custom)))) stop('Invalid advanced group colors')
  colors[intersect(names(custom), order)] <- custom[intersect(names(custom), order)]
  list(group_order = order, group_colors = colors)
}
adv_theme <- function(p, config) {
  g <- config$general; s <- config$profile
  theme <- switch(g$theme, bw = ggplot2::theme_bw, classic = ggplot2::theme_classic,
    minimal = ggplot2::theme_minimal)
  if (is.null(theme)) stop('Invalid advanced plot theme')
  p <- p + theme(base_size = g$base_size, base_family = config$font) +
    ggplot2::theme(legend.position = g$legend_position)
  if (!is.null(g$title_size)) p <- p + ggplot2::theme(plot.title = ggplot2::element_text(size = g$title_size))
  if (!is.null(g$axis_text_size)) p <- p + ggplot2::theme(axis.text = ggplot2::element_text(size = g$axis_text_size))
  if (!is.null(g$legend_text_size)) p <- p + ggplot2::theme(legend.text = ggplot2::element_text(size = g$legend_text_size))
  if (!is.null(s$label_size)) p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(size = s$label_size))
  if (g$background == 'transparent') p <- p + ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = 'transparent', colour = NA),
    panel.background = ggplot2::element_rect(fill = 'transparent', colour = NA),
    legend.background = ggplot2::element_rect(fill = 'transparent', colour = NA),
    legend.key = ggplot2::element_rect(fill = 'transparent', colour = NA))
  for (key in c('title', 'x_label', 'y_label')) if (nzchar(s[[key]])) {
    target <- switch(key, x_label = 'x', y_label = 'y', title = 'title')
    p <- p + do.call(ggplot2::labs, setNames(list(s[[key]]), target))
  }
  if (!isTRUE(g$show_caption)) p <- p + ggplot2::labs(caption = NULL, subtitle = NULL)
  p
}
adv_target <- function(target) {
  profiles <- c(core_heatmap = 'heatmap', occupancy = 'bubble', CLR_PCA = 'beta', constrained = 'beta', association_heatmap = 'heatmap')
  if (is.null(target)) return(invisible(TRUE))
  if (!is.list(target) || length(target$stem) != 1 || !is.character(target$stem) ||
      !target$stem %in% paste0('03_figures/', names(profiles)) ||
      !identical(target$profile, unname(profiles[sub('^03_figures/', '', target$stem)]))) stop('Invalid advanced cached figure target')
  invisible(TRUE)
}
adv_render <- function(cache, o, run, cache_run) {
  adv_target(o$plot_target)
  folder <- file.path(run, '03_figures'); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  style <- adv_group_style(cache, o); colors <- adv_colors(o); figures <- list()
  save <- function(p, name, profile, width = 8, height = 5) {
    stem <- paste0('03_figures/', name)
    if (!is.null(o$plot_target) && o$plot_target$stem != stem) return(invisible(NULL))
    config <- adv_config(o, profile); p <- adv_theme(p, config)
    width <- config$profile$width %||% width; height <- config$profile$height %||% height
    formats <- if (config$general$formats == 'png_svg') c('png', 'svg') else c('png', 'svg', 'pdf')
    files <- paste0(file.path(run, stem), '.', formats)
    for (i in seq_along(files)) {
      device <- switch(formats[i], svg = svglite::svglite, pdf = grDevices::cairo_pdf, png = 'png')
      ggplot2::ggsave(files[i], plot = p, device = device, width = width, height = height,
        dpi = config$general$dpi, bg = config$general$background)
    }
    effective <- list(profile = profile, width = width, height = height, dpi = config$general$dpi,
      font = config$font, formats = as.list(formats), general = config$general,
      profile_settings = config$profile, colors = colors,
      labels = as.list(p$labels[c('title', 'x', 'y', 'caption', 'subtitle')]))
    adv_json(list(requested = o$plot_settings %||% list(), effective = effective),
      file.path(run, paste0(stem, '_plot_settings.json')))
    figures[[length(figures) + 1L]] <<- list(stem = stem, profile = profile,
      files = as.list(paste0(stem, '.', formats)), labels = effective$labels)
  }
  top <- function(d, n = cache$options$top_n) {
    global <- cache$data$occupancy
    sums <- if (!is.null(global)) setNames(global$MeanRelativeAbundance, global$FeatureID) else
      tapply(d$MeanRelativeAbundance, d$FeatureID, mean)
    names(head(sort(sums, decreasing = TRUE), n))
  }
  d <- cache$data$core
  if (!is.null(d)) {
    config <- adv_config(o, 'heatmap')
    keep <- top(d, config$profile$top_n %||% cache$options$top_n); d <- d[d$FeatureID %in% keep, , drop = FALSE]
    order <- if (config$profile$row_order == 'alphabetical') sort(keep) else keep
    d$Feature <- factor(d$FeatureID, levels = rev(order)); d$Group <- factor(d$Group, levels = style$group_order)
    config <- adv_config(o, 'heatmap')
    p <- ggplot2::ggplot(d, ggplot2::aes(x = Group, y = Feature, fill = Prevalence)) +
      ggplot2::geom_tile(colour = 'white') + ggplot2::scale_fill_gradient2(low = colors$heatmap_low,
        mid = colors$heatmap_mid, high = colors$heatmap_high, midpoint = .5, limits = c(0, 1)) +
      ggplot2::scale_y_discrete(labels = cache$labels) +
      ggplot2::labs(title = 'Core taxa prevalence by group', x = NULL, y = NULL, fill = 'Prevalence',
        caption = paste0('Descriptive. Detection: reads >= ', cache$options$core_min_reads,
          ' and full-table relative-abundance fraction >= ', cache$options$core_min_relative,
          '. Core prevalence fraction >= ', cache$options$core_prevalence,
          '. * marks core taxa; display limited to ', length(keep), ' taxa by mean abundance.'))
    d$Display <- if (isTRUE(config$profile$show_values)) paste0(round(d$Prevalence * 100), '%', ifelse(d$Core, '*', '')) else ifelse(d$Core, '*', '')
    p <- p + ggplot2::geom_text(data = d, ggplot2::aes(label = Display), size = 3, family = config$font)
    save(p, 'core_heatmap', 'heatmap', height = max(4, length(keep) * .25 + 2))
    d <- cache$data$core; d$Group <- factor(d$Group, levels = style$group_order)
    config <- adv_config(o, 'bubble')
    if (!is.null(config$profile$top_n)) d <- d[d$FeatureID %in% top(d, config$profile$top_n), , drop = FALSE]
    p <- ggplot2::ggplot(d, ggplot2::aes(x = MeanRelativeAbundance * 100, y = Prevalence * 100,
      colour = Group, fill = MeanRelativeAbundance * 100, size = MeanRelativeAbundance * 100)) +
      ggplot2::geom_hline(yintercept = cache$options$core_prevalence * 100, linetype = 2, colour = 'grey50') +
      ggplot2::geom_point(alpha = config$profile$point_alpha, shape = 21, stroke = .8) +
      ggplot2::scale_size_continuous(range = c(1, config$profile$max_size)) +
      ggplot2::scale_color_manual(values = style$group_colors, breaks = style$group_order) +
      ggplot2::scale_fill_gradient(low = colors$bubble_low, high = colors$bubble_high) +
      ggplot2::labs(title = 'Taxon prevalence and mean abundance', x = 'Mean full-table relative abundance (%)',
        y = 'Detected samples (%)', size = 'Mean abundance (%)', fill = 'Mean abundance (%)',
        caption = 'Descriptive; each point is one lineage within one Group. Dashed line is the selected core prevalence threshold.')
    if (isTRUE(config$profile$show_values)) {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = sprintf('%.2f%%', MeanRelativeAbundance * 100)),
        size = 2.5, family = config$font, vjust = -1, show.legend = FALSE)
    }
    save(p, 'occupancy', 'bubble')
  }
  ordination <- function(value, name, title, axes, note) {
    d <- value$coordinates; d$Group <- factor(d$Group, levels = style$group_order)
    config <- adv_config(o, 'beta'); s <- config$profile
    p <- ggplot2::ggplot(d, ggplot2::aes(x = Axis1, y = Axis2, colour = Group))
    if (isTRUE(s$show_ellipse)) {
      # Singular groups have no estimable covariance ellipse.
      good <- names(which(vapply(split(d, d$Group), function(z) nrow(z) >= 4 &&
        qr(cbind(z$Axis1 - mean(z$Axis1), z$Axis2 - mean(z$Axis2)))$rank == 2, logical(1))))
      ell <- d[as.character(d$Group) %in% good, , drop = FALSE]
      if (nrow(ell)) p <- p + ggplot2::stat_ellipse(data = ell, ggplot2::aes(fill = Group),
        geom = 'polygon', level = s$ellipse_level, alpha = s$ellipse_alpha, colour = NA, type = 'norm', show.legend = FALSE)
    }
    p <- p + ggplot2::geom_point(size = s$point_size, alpha = s$point_alpha) +
      ggplot2::scale_color_manual(values = style$group_colors, breaks = style$group_order) +
      ggplot2::scale_fill_manual(values = style$group_colors, breaks = style$group_order) +
      ggplot2::coord_equal() + ggplot2::labs(title = title, x = axes[1], y = axes[2], caption = note)
    if (isTRUE(s$show_sample_labels)) p <- p + ggplot2::geom_text(ggplot2::aes(label = SampleID),
      size = 2.5, vjust = -1, family = config$font, show.legend = FALSE)
    if (isTRUE(s$show_stats)) {
      tab <- if (name == 'CLR_PCA') cache$data$permanova else value$overall
      if (!is.null(tab)) {
        hit <- which(tab$Status == 'tested')
        if (name == 'CLR_PCA') hit <- hit[tab$Term[hit] == 'Group']
        if (length(hit)) p <- p + ggplot2::labs(subtitle = paste(if (name == 'CLR_PCA') 'Separate Bray-Curtis PERMANOVA Group' else 'Constrained overall permutation test',
          'P =', format.pval(tab$P[hit[1]], digits = 3)))
      }
    }
    save(p, name, 'beta')
  }
  if (!is.null(cache$data$pca)) {
    d <- cache$data$pca; v <- d$variance$Fraction
    ordination(d, 'CLR_PCA', 'CLR principal component analysis',
      sprintf('PC%d (%.1f%%)', 1:2, 100 * v[1:2]),
      paste0('Exploratory. Raw-count pseudocount = ', cache$options$clr_pseudocount,
        '; CLR of retained subcomposition; PCA centered, unscaled. Ellipses are descriptive; independent samples required for tests.'))
  }
  if (!is.null(cache$data$constrained)) {
    d <- cache$data$constrained; v <- d$variance$FractionOfTotal
    axes <- c(paste0(d$variance$Axis[1], ' (', round(100 * v[1], 1), '% total inertia)'),
      if (nrow(d$variance) >= 2) paste0(d$variance$Axis[2], ' (', round(100 * v[2], 1), '% total inertia)') else 'No second constrained axis')
    ordination(d, 'constrained', if (d$method == 'hellinger_rda') 'Hellinger redundancy analysis' else 'Bray-Curtis distance based RDA',
      axes, paste0('Exploratory; predictors: ', paste(cache$options$predictors, collapse = ', '),
        '; conditioning covariates: ', paste(cache$options$covariates, collapse = ', '),
        '. Axis percentages use total inertia. Permutation P-values and BH families are saved in tables.'))
  }
  d <- cache$data$association
  if (!is.null(d) && any(d$Status == 'tested')) {
    d <- d[d$Status == 'tested', , drop = FALSE]
    scores <- tapply(abs(d$Rho), d$FeatureID, max)
    config <- adv_config(o, 'heatmap')
    keep <- names(head(sort(scores, decreasing = TRUE), config$profile$top_n %||% cache$options$top_n))
    d <- d[d$FeatureID %in% keep, , drop = FALSE]
    if (config$profile$row_order == 'alphabetical') keep <- sort(keep)
    d$Feature <- factor(d$FeatureID, levels = rev(keep)); d$Phenotype <- factor(d$Phenotype, levels = cache$options$phenotypes)
    d$Display <- if (isTRUE(config$profile$show_values)) paste0(sprintf('%.2f', d$Rho), ifelse(d$Q < .05, '*', '')) else ifelse(d$Q < .05, '*', '')
    caption <- paste0(cache$options$association_method, '; ',
      adv_association_definition(cache$options$association_transform %||% 'relative', cache$options$clr_pseudocount), '; covariates: ',
      paste(if (cache$options$association_method == 'partial_spearman') cache$options$association_covariates else 'none', collapse = ', '),
      '. * BH Q < 0.05 across all tested taxa x traits. Missing or invalid pairs are absent; sample sizes are in tables.')
    p <- ggplot2::ggplot(d, ggplot2::aes(x = Phenotype, y = Feature, fill = Rho)) +
      ggplot2::geom_tile(colour = 'white') + ggplot2::scale_fill_gradient2(low = colors$heatmap_low,
        mid = colors$heatmap_mid, high = colors$heatmap_high, limits = c(-1, 1)) +
      ggplot2::scale_y_discrete(labels = cache$labels) +
      ggplot2::geom_text(ggplot2::aes(label = Display), size = 3, family = config$font) +
      ggplot2::labs(title = 'Exploratory phenotype associations', x = NULL, y = NULL, fill = 'Spearman rho',
        caption = paste(strwrap(caption, width = 110), collapse = '\n'))
    save(p, 'association_heatmap', 'heatmap', height = max(4, length(keep) * .25 + 2))
  }
  if (!is.null(o$plot_target) && length(figures) != 1) stop('Selected figure is absent from the saved advanced analysis: ', o$plot_target$stem)
  adv_json(list(source = 'advanced_analysis', advanced_run = cache_run, statistics_refit = FALSE,
    cache_file = file.path(cache_run, '02_statistics', 'advanced_cache.rds'), font = adv_config(o, 'beta')$font,
    top_n = o$top_n,
    plot_settings = o$plot_settings %||% list(), group_style = style, visual_options = colors,
    figures = figures), file.path(folder, 'render_summary.json'))
  invisible(figures)
}
adv_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3 || !args[1] %in% c('--check', '--run', '--render')) stop('Usage: advanced_analysis.R --check|--run|--render OPTIONS.json RUN')
  if (!requireNamespace('jsonlite', quietly = TRUE)) stop('jsonlite is required')
  request <- jsonlite::read_json(args[2], simplifyVector = FALSE)
  # JSON scalar lists are intentionally retained for nested plot settings.
  o <- adv_defaults(request); run <- normalizePath(args[3], mustWork = FALSE)
  if (args[1] == '--check') {
    versions <- adv_require(); adv_config(o, 'beta')
    dir.create(run, recursive = TRUE, showWarnings = FALSE)
    adv_json(list(status = 'ready', source = 'advanced_analysis', packages = as.list(versions),
      independent_samples_only = TRUE), file.path(run, 'advanced_check.json'))
    cat('Advanced R packages and font are ready\n'); return(invisible(TRUE))
  }
  if (args[1] == '--run') {
    cache <- adv_fit(o, run)
    adv_render(cache, o, run, run)
  } else {
    if (is.null(request$cache_run) || length(request$cache_run) != 1) stop('cache_run is required for cached rendering')
    cache_run <- normalizePath(request$cache_run, mustWork = TRUE)
    cache <- readRDS(file.path(cache_run, '02_statistics', 'advanced_cache.rds'))
    if (!identical(cache$version, 1L)) stop('Unsupported advanced statistics cache version')
    visual <- cache$options
    for (key in c('font', 'plot_settings', 'group_order', 'group_colors', 'plot_target', names(adv_colors(list()))))
      if (!is.null(request[[key]])) visual[[key]] <- request[[key]]
    adv_require(); adv_render(cache, visual, run, cache_run)
  }
  invisible(TRUE)
}
if (sys.nframe() == 0L) adv_main()
