#!/usr/bin/env Rscript
# Calls the official ANCOMBC::ancombc2; no reimplementation of its statistics.
read_count_input <- function(filename) {
  tab <- read.delim(filename, check.names = FALSE, stringsAsFactors = FALSE,
                    colClasses = 'character', quote = '', comment.char = '', na.strings = NULL)
  if (!identical(names(tab)[1], 'FeatureID') || ncol(tab) < 3L || nrow(tab) < 1L)
    stop('Canonical FeatureID count table is missing samples or taxa')
  if (anyDuplicated(tab[[1]]) || anyDuplicated(names(tab)[-1])) stop('Duplicate count identifiers')
  mat <- matrix(suppressWarnings(as.numeric(as.matrix(tab[-1]))), nrow = nrow(tab),
                dimnames = list(tab[[1]], names(tab)[-1]))
  if (any(!is.finite(mat)) || any(mat < 0) || any(mat != floor(mat)))
    stop('ANCOM-BC2 requires nonnegative integer unrarefied counts')
  mat
}

read_group_metadata <- function(filename, ids, reference, levels) {
  meta <- read.delim(filename, check.names = FALSE, colClasses = 'character',
                     stringsAsFactors = FALSE, quote = '', comment.char = '', na.strings = NULL)
  if (!identical(names(meta), c('sample-id', 'Group'))) stop('Expected sample-id and actual Group columns')
  meta <- meta[meta[['sample-id']] != '#q2:types', , drop = FALSE]
  if (anyDuplicated(meta[['sample-id']]) || !setequal(meta[['sample-id']], ids))
    stop('Metadata IDs do not match the actual count table')
  meta <- meta[match(ids, meta[['sample-id']]), , drop = FALSE]
  if (any(!nzchar(meta$Group)) || !reference %in% meta$Group) stop('Missing Group or reference group')
  meta$Group <- factor(meta$Group, levels = c(reference, setdiff(levels, reference)), ordered = FALSE)
  if (anyNA(meta$Group) || length(unique(meta$Group)) < 2L || any(table(meta$Group) < 2L))
    stop('At least two actual groups and two samples per group are required')
  rownames(meta) <- meta[['sample-id']]
  meta
}

write_result <- function(value, filename) {
  if (is.null(value)) return(FALSE)
  if (is.matrix(value)) value <- data.frame(taxon = rownames(value), value, check.names = FALSE)
  if (!is.data.frame(value)) return(FALSE)
  write.table(value, filename, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
  TRUE
}

fit_ancom_model <- function(counts, meta, opts, fit_function = ANCOMBC::ancombc2) {
  # Official matrix input avoids the optional microbiome dependency of phyloseq.
  if (!is.matrix(counts) || !is.numeric(counts) || !is.data.frame(meta) ||
      !identical(colnames(counts), rownames(meta)))
    stop('ANCOM-BC2 matrix columns and metadata rows must match exactly')
  fit_function(
    data = counts, taxa_are_rows = TRUE, meta_data = meta,
    tax_level = NULL, fix_formula = 'Group', rand_formula = NULL,
    p_adj_method = opts$p_adj_method, pseudo_sens = opts$pseudo_sens,
    prv_cut = opts$prv_cut, lib_cut = opts$lib_cut, group = 'Group',
    struc_zero = opts$struc_zero, neg_lb = opts$neg_lb, alpha = opts$alpha,
    n_cl = opts$threads, verbose = TRUE, global = opts$global, pairwise = opts$pairwise,
    dunnet = FALSE, trend = FALSE,
    iter_control = opts$iter_control, em_control = opts$em_control,
    mdfdr_control = opts$mdfdr_control)
}

effect_group_style <- function(actual, requested = NULL) {
  order <- if (is.null(requested$group_order)) character() else unlist(requested$group_order, use.names = FALSE)
  if (!is.character(order) || anyNA(order) || any(!nzchar(order)) || anyDuplicated(order))
    stop('Invalid group display order')
  order <- c(intersect(order, actual), setdiff(actual, order))
  colors <- c(CON = '#D55E00', LDC = '#0072B2', HDC = '#CC79A7')
  extra <- setdiff(actual, names(colors))
  if (length(extra)) colors <- c(colors, setNames(grDevices::hcl.colors(length(extra), 'Dark 3'), extra))
  custom <- if (is.null(requested$group_colors)) character() else unlist(requested$group_colors, use.names = TRUE)
  if (length(custom) && (is.null(names(custom)) || anyDuplicated(names(custom)) || anyNA(custom) ||
                        any(!grepl('^#[[:xdigit:]]{6}$', custom)))) stop('Invalid group display colors; use #RRGGBB')
  colors[intersect(names(custom), actual)] <- custom[intersect(names(custom), actual)]
  list(group_order = order, group_colors = colors[order])
}

plot_effects <- function(res, directory, reference, font, group_style = NULL) {
  fields <- if (is.null(res) || !'taxon' %in% names(res)) character() else grep('^lfc_Group', names(res), value = TRUE)
  targets <- sub('^lfc_Group', '', fields)
  style <- effect_group_style(unique(c(reference, targets)), group_style)
  fields <- fields[order(match(targets, style$group_order))]
  plots <- character()
  for (i in seq_along(fields)) {
    term <- sub('^lfc_', '', fields[i])
    qfield <- paste0('q_', term); sefield <- paste0('se_', term)
    if (!all(c(qfield, sefield) %in% names(res))) next
    dat <- data.frame(taxon = res$taxon, lfc = res[[fields[i]]], se = res[[sefield]], q = res[[qfield]])
    dat <- dat[is.finite(dat$lfc) & is.finite(dat$se) & is.finite(dat$q), , drop = FALSE]
    if (!nrow(dat)) next
    dat <- head(dat[order(dat$q, -abs(dat$lfc)), , drop = FALSE], 25L)
    target <- sub('^Group', '', term)
    # Only figure order/color changes. The fit, reference, effects and q-values remain intact.
    dat$direction_group <- ifelse(dat$lfc > 0, target, ifelse(dat$lfc < 0, reference, NA_character_))
    present <- style$group_order
    dat$direction_group <- factor(dat$direction_group, levels = present)
    dat$taxon <- factor(dat$taxon, levels = rev(dat$taxon))
    write_result(dat, file.path(directory, paste0('effect_', i, '_plot_input.tsv')))
    fig <- ggplot2::ggplot(dat, ggplot2::aes(x = lfc, y = taxon, color = direction_group)) +
      ggplot2::geom_vline(xintercept = 0, color = 'grey70') +
      ggplot2::geom_segment(ggplot2::aes(x = lfc - se, xend = lfc + se, yend = taxon)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = style$group_colors, breaks = present, na.value = '#777777') +
      ggplot2::labs(x = 'ANCOM-BC2 log fold change (± SE)', y = NULL, color = 'Effect direction',
                    title = paste(target, 'versus', reference),
                    caption = paste0('Positive: ', target, ' − ', reference,
                                     '; negative: reference is higher; zero effect is gray. Adjusted q-values remain in the exported tables.')) +
      ggplot2::theme_bw(base_family = font) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    for (ext in c('png', 'svg')) {
      name <- file.path(directory, paste0('effect_', i, '.', ext))
      ggplot2::ggsave(name, fig, width = 12, height = max(4, nrow(dat) * .25), dpi = 200)
      plots <- c(plots, name)
    }
  }
  if (!length(plots)) {
    fig <- ggplot2::ggplot() + ggplot2::annotate('text', x = 0, y = 0,
             label = 'No estimable finite Group effects.\nInspect model results, structural zeros and warnings.', family = font) +
      ggplot2::theme_void() + ggplot2::labs(title = 'ANCOM-BC2')
    for (ext in c('png', 'svg')) {
      filename <- file.path(directory, paste0('no_estimable_effects.', ext))
      ggplot2::ggsave(filename, fig, width = 10, height = 4, dpi = 200)
      plots <- c(plots, filename)
    }
  }
  plots
}

plot_qiime <- function(opts, run) {
  for (entry in opts$tables) {
    folder <- file.path(run, '02_differential', paste0('L', entry$level, '_', entry$name))
    result <- read.delim(file.path(folder, 'ancom_results.tsv'), check.names = FALSE,
                         stringsAsFactors = FALSE, colClasses = 'character')
    if (!all(c('W', 'Reject null hypothesis') %in% names(result))) stop('Official QIIME ANCOM W/decision columns missing')
    result$W <- suppressWarnings(as.numeric(result$W))
    if (any(!is.finite(result$W))) stop('QIIME ANCOM W must be finite in official export')
    result$taxon <- result[[1]]
    result$decision <- factor(result[['Reject null hypothesis']], levels = c('False', 'True', 'FALSE', 'TRUE'))
    dat <- head(result[order(-result$W), , drop = FALSE], 25L)
    dat$taxon <- factor(dat$taxon, levels = rev(dat$taxon))
    fig <- ggplot2::ggplot(dat, ggplot2::aes(x = W, y = taxon, fill = decision)) +
      ggplot2::geom_col() + ggplot2::labs(x = 'Official QIIME 2 ANCOM W', y = NULL,
                                        fill = 'Reject null', title = paste('QIIME 2 ANCOM:', entry$name)) +
      ggplot2::theme_bw(base_family = opts$font) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    for (ext in c('png', 'svg')) ggplot2::ggsave(file.path(folder, paste0('ancom_W.', ext)), fig,
                              width = 12, height = max(4, nrow(dat) * .25), dpi = 200)
  }
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L || !args[1] %in% c('--check', '--run', '--plot-qiime'))
    stop('Usage: differential_ancombc.R --check|--run|--plot-qiime options.json run_directory')
  required <- if (args[1] == '--plot-qiime') c('jsonlite', 'ggplot2', 'svglite') else
    c('ANCOMBC', 'jsonlite', 'ggplot2', 'svglite')
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (any(!available) && args[1] == '--plot-qiime') stop(paste('Missing isolated R plotting dependencies:',
                    paste(required[!available], collapse = ', '), '; run 安装依赖.sh for the selected plotting r_env.'))
  if (any(!available)) stop(paste('Missing isolated ANCOM-BC2 dependencies:',
                                  paste(required[!available], collapse = ', '),
                                  '; run install_differential.sh in Ubuntu/WSL, or select a prepared ancom_env.'))
  if (args[1] != '--plot-qiime' && !'ancombc2' %in% getNamespaceExports('ANCOMBC')) stop('Installed ANCOMBC does not provide ancombc2')
  opts <- jsonlite::fromJSON(args[2], simplifyVector = FALSE)
  run <- normalizePath(args[3], mustWork = TRUE)
  versions <- list(R = R.version.string, packages = as.list(vapply(required,
                    function(pkg) as.character(utils::packageVersion(pkg)), character(1))))
  jsonlite::write_json(versions, file.path(run, 'software_versions.json'), pretty = TRUE, auto_unbox = TRUE)
  if (args[1] == '--plot-qiime') { plot_qiime(opts, run); return(invisible(TRUE)) }
  if (args[1] == '--check') return(invisible(TRUE))
  destination <- file.path(run, '02_differential'); dir.create(destination, showWarnings = FALSE)
  summaries <- list()
  for (entry in opts$tables) {
    folder <- file.path(destination, paste0('L', entry$level, '_', entry$name))
    dir.create(folder, showWarnings = FALSE)
    counts <- read_count_input(entry$counts_tsv)
    meta <- read_group_metadata(opts$metadata_tsv, colnames(counts), opts$reference_group, unlist(opts$groups))
    if (any(colSums(counts) <= 0)) stop('Zero count library; no pseudocount replacement is permitted')
    warnings <- character()
    started <- format(Sys.time(), tz = 'UTC', usetz = TRUE)
    set.seed(opts$seed)
    model_error <- NULL
    fit <- tryCatch(withCallingHandlers(fit_ancom_model(counts, meta, opts),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)) }),
      error = function(e) { model_error <<- conditionMessage(e); NULL })
    if (!is.null(model_error)) {
      message('L', entry$level, ' (', entry$name, ') ANCOM-BC2 failed: ', model_error)
      failure <- list(level = entry$level, method = 'ANCOM-BC2', status = 'failed',
                       error = model_error, warnings = unique(warnings),
                       note = 'No p/q or effect estimates created for a failed model')
      jsonlite::write_json(failure, file.path(folder, 'diagnostics.json'), pretty = TRUE, auto_unbox = TRUE)
      write_result(data.frame(taxon = rownames(counts), status = 'model_failed',
                               lfc = NA_real_, se = NA_real_, p = NA_real_, q = NA_real_),
                   file.path(folder, 'failed_taxa.tsv'))
      summaries[[length(summaries) + 1L]] <- failure
      next
    }
    saveRDS(fit, file.path(folder, 'ancombc2_fit.rds'))
    for (key in c('res', 'res_global', 'res_pair', 'zero_ind', 'ss_tab', 'feature_table'))
      write_result(fit[[key]], file.path(folder, paste0(key, '.tsv')))
    fit_taxa <- if (!is.null(fit$res)) as.character(fit$res$taxon) else character()
    status <- data.frame(taxon = rownames(counts), prevalence = rowMeans(counts > 0),
                          total_count = rowSums(counts), in_model_results = rownames(counts) %in% fit_taxa,
                          check.names = FALSE)
    status$note <- ifelse(status$in_model_results, 'See official model result including NA',
                          'Not in model results; inspect prevalence, structural zero table and warnings')
    write_result(status, file.path(folder, 'taxon_status.tsv'))
    plots <- plot_effects(fit$res, folder, opts$reference_group, if (is.null(opts$font)) 'sans' else opts$font,
                          opts$group_style)
    diagnostics <- list(level = entry$level, method = 'ANCOM-BC2', status = 'executed', started_at = started,
                         finished_at = format(Sys.time(), tz = 'UTC', usetz = TRUE), warnings = unique(warnings),
                         convergence = 'Not inferred; use actual package warnings, saved fit and workflow.log',
                         missing_values = 'NA is preserved in all exported result tables',
                         primary_tests = 'Nonreference nominal Group levels versus explicit reference',
                         global = opts$global, pairwise = opts$pairwise,
                         pairwise_adjustment = 'Official ANCOM-BC2 mixed directional FDR with Holm FWER control and B=100',
                         structural_zeros = 'Reported separately by official zero_ind; no artificial finite p/q values',
                         feature_count = nrow(counts), modeled_features = length(fit_taxa), plots = plots)
    jsonlite::write_json(diagnostics, file.path(folder, 'diagnostics.json'), pretty = TRUE, auto_unbox = TRUE)
    summaries[[length(summaries) + 1L]] <- diagnostics
  }
  failed <- sum(vapply(summaries, function(x) identical(x$status, 'failed'), logical(1)))
  jsonlite::write_json(list(method = 'ANCOM-BC2', status = if (failed) 'partial_failure' else 'executed',
                            failed_levels = failed, levels = summaries), file.path(destination, 'summary.json'),
                        pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(run, 'R_session_info.txt'))
  if (failed) stop(paste(failed, 'taxonomy level(s) failed; actual outputs and failure diagnostics are preserved'))
}

if (sys.nframe() == 0L) main()
