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

plot_effects <- function(res, directory, reference, font, group_style = NULL, plot_settings = NULL) {
  opts <- list(font = font, plot_settings = plot_settings)
  font <- diff_plot_font(opts); opts$font <- font
  point <- diff_plot_config(opts, 'effect')$profile
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
      ggplot2::geom_vline(xintercept = 0, color = 'grey70')
    if (isTRUE(point$show_errorbars)) fig <- fig + ggplot2::geom_segment(ggplot2::aes(x = lfc - se, xend = lfc + se, yend = taxon))
    fig <- fig + ggplot2::geom_point(size = if (is.null(point$point_size)) 2 else point$point_size) +
      ggplot2::scale_color_manual(values = style$group_colors, breaks = present, na.value = '#777777') +
      ggplot2::labs(x = if (isTRUE(point$show_errorbars)) 'ANCOM-BC2 log fold change (± SE)' else 'ANCOM-BC2 log fold change', y = NULL, color = 'Effect direction',
                    title = paste(target, 'versus', reference),
                    caption = paste0('Positive: ', target, ' − ', reference,
                                     '; negative: reference is higher; zero effect is gray. Adjusted q-values remain in the exported tables.')) +
      ggplot2::theme_bw(base_family = font) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    plots <- c(plots, save_differential_figure(fig, file.path(directory, paste0('effect_', i)),
                                             height = max(4, nrow(dat) * .25), opts = opts, comparison = paste(target, 'versus', reference)))
  }
  if (!length(plots)) {
    fig <- ggplot2::ggplot() + ggplot2::annotate('text', x = 0, y = 0,
             label = 'No estimable finite Group effects.\nInspect model results, structural zeros and warnings.', family = font) +
      ggplot2::theme_void() + ggplot2::labs(title = 'ANCOM-BC2')
    attr(fig, 'differential_no_hit') <- TRUE
    plots <- c(plots, save_differential_figure(fig, file.path(directory, 'no_estimable_effects'), width = 10, height = 4, opts = opts))
  }
  plots
}

visual_options <- function(opts) {
  defaults <- list(diff_comparisons = 'all_pairs', diff_heatmap_scale = 'relative_abundance',
                   heatmap_low = '#2166AC', heatmap_mid = '#F7F7F7', heatmap_high = '#B2182B',
                   bubble_color_mode = 'abundance', bubble_low = '#E8F1FA', bubble_high = '#2166AC',
                   diff_plot_top_n = 30L)
  for (key in names(defaults)) if (!is.null(opts[[key]])) defaults[[key]] <- opts[[key]]
  if (!defaults$diff_comparisons %in% c('all_pairs', 'reference') ||
      !defaults$diff_heatmap_scale %in% c('row_zscore', 'relative_abundance') ||
      !defaults$bubble_color_mode %in% c('abundance', 'group')) stop('Invalid differential visualization option')
  for (key in c('heatmap_low', 'heatmap_mid', 'heatmap_high', 'bubble_low', 'bubble_high'))
    if (length(defaults[[key]]) != 1L || !grepl('^#[[:xdigit:]]{6}$', defaults[[key]])) stop('Use #RRGGBB plot colors: ', key)
  n <- defaults$diff_plot_top_n
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n != floor(n) || n < 1 || n > 200)
    stop('Differential plot top N must be an integer from 1 to 200')
  defaults$diff_plot_top_n <- as.integer(n)
  defaults
}

diff_plot_config <- function(opts = NULL, profile = 'effect') {
  requested <- if (is.null(opts$plot_settings)) list() else opts$plot_settings
  general <- list(dpi = NULL, font_family = 'original', base_size = NULL, title_size = NULL, axis_text_size = NULL,
                  legend_text_size = NULL, legend_position = 'original', theme = 'original',
                  background = 'white', show_caption = TRUE, formats = 'original')
  specific <- list(width = NULL, height = NULL, title = '', x_label = '', y_label = '')
  if (profile == 'heatmap') specific <- c(specific, list(label_size = NULL, row_order = 'abundance', show_values = FALSE))
  if (profile == 'bubble') specific <- c(specific, list(max_size = NULL, point_alpha = .9, show_values = FALSE))
  if (profile == 'effect') specific <- c(specific, list(point_size = NULL, show_errorbars = TRUE))
  for (key in names(general)) if (!is.null(requested$general[[key]])) general[[key]] <- requested$general[[key]]
  for (key in names(specific)) if (!is.null(requested[[profile]][[key]])) specific[[key]] <- requested[[profile]][[key]]
  list(general = general, profile = specific, requested = requested)
}

diff_plot_font <- function(opts) {
  selected <- diff_plot_config(opts)$general$font_family
  if (selected == 'original') return(if (is.null(opts$font)) 'sans' else opts$font)
  if (!requireNamespace('systemfonts', quietly = TRUE) || !selected %in% systemfonts::system_fonts()$family)
    stop('Selected plot font is not installed in this R environment: ', selected)
  selected
}

differential_figures <- list()
validate_diff_target <- function(target) {
  if (is.null(target)) return(invisible(TRUE))
  suffix <- c(heatmap='abundance/abundance_heatmap', bubble='abundance/abundance_bubble',
              effect='(comparisons/pair_[0-9]{3}/effect|ancom_W)')
  if (!is.list(target) || !is.character(target$profile) || length(target$profile)!=1L ||
      !target$profile %in% names(suffix) || !is.character(target$stem) || length(target$stem)!=1L ||
      !grepl(paste0('^02_differential/L[2-7]_[a-z]+/',suffix[[target$profile]],'$'),target$stem))
    stop('Invalid single differential figure target')
  if (!is.null(target$comparison_key) && (!identical(target$profile,'effect') ||
      !is.character(target$comparison_key) || length(target$comparison_key)!=1L || is.na(target$comparison_key) ||
      !nzchar(target$comparison_key))) stop('Invalid selected original comparison key')
  invisible(TRUE)
}

diff_relative_stem <- function(opts, stem) {
  root <- opts$plot_output_run
  if (is.null(root) || !startsWith(stem,paste0(root,'/'))) return(stem)
  substring(stem,nchar(root)+2L)
}

diff_target_wanted <- function(opts, stem, profile) {
  target <- opts$plot_target
  is.null(target) || (identical(target$profile,profile) && identical(target$stem,diff_relative_stem(opts,stem)))
}

diff_target_tables <- function(opts) {
  validate_diff_target(opts$plot_target)
  if (is.null(opts$plot_target)) return(opts)
  level <- strsplit(opts$plot_target$stem,'/',fixed=TRUE)[[1]][2]
  keep <- vapply(opts$tables,function(entry) identical(paste0('L',entry$level,'_',entry$name),level),logical(1))
  if (sum(keep)!=1L) stop('Requested figure level is absent from the original saved analysis')
  opts$tables <- opts$tables[keep]
  opts
}

write_diff_figure_manifest <- function(run, target=NULL) {
  if (!is.null(target) && length(differential_figures)!=1L) stop('The requested cached differential figure could not be produced: ',target$stem)
  jsonlite::write_json(list(figures=differential_figures),file.path(run,'figure_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
}

diff_label <- function(value, rank, comparison) {
  gsub('{comparison}', comparison, gsub('{rank}', rank, value, fixed = TRUE), fixed = TRUE)
}

style_differential_figure <- function(fig, config, font = 'sans', rank = '', comparison = '') {
  general <- config$general; profile <- config$profile
  no_hit <- isTRUE(attr(fig, 'differential_no_hit'))
  if (general$theme != 'original') {
    theme <- switch(general$theme, classic = ggplot2::theme_classic, minimal = ggplot2::theme_minimal, bw = ggplot2::theme_bw)
    if (is.null(theme)) stop('Invalid plot theme')
    # Preserve geometry layout safeguards while replacing the requested visual theme.
    original <- fig$theme
    fig <- fig + theme(base_size = if (is.null(general$base_size)) 11 else general$base_size, base_family = font)
    if (!is.null(original$axis.text.x)) fig <- fig + ggplot2::theme(axis.text.x = original$axis.text.x)
    if (!is.null(original$axis.text.y)) fig <- fig + ggplot2::theme(axis.text.y = original$axis.text.y)
    for (key in c('plot.margin', 'legend.key.height', 'legend.spacing.y', 'legend.position'))
      if (!is.null(original[[key]])) fig <- fig + do.call(ggplot2::theme, setNames(list(original[[key]]), key))
  }
  if (!is.null(general$base_size)) fig <- fig + ggplot2::theme(text = ggplot2::element_text(size = general$base_size))
  if (general$font_family != 'original') fig <- fig + ggplot2::theme(text = ggplot2::element_text(family = font))
  if (!is.null(general$title_size)) fig <- fig + ggplot2::theme(plot.title = ggplot2::element_text(size = general$title_size))
  if (!is.null(general$axis_text_size)) fig <- fig + ggplot2::theme(axis.text.x = ggplot2::element_text(size = general$axis_text_size),
                                                               axis.text.y = ggplot2::element_text(size = general$axis_text_size))
  if (!is.null(general$legend_text_size)) fig <- fig + ggplot2::theme(legend.text = ggplot2::element_text(size = general$legend_text_size))
  if (!is.null(profile$label_size)) fig <- fig + ggplot2::theme(axis.text.y = ggplot2::element_text(size = profile$label_size))
  if (general$legend_position != 'original') fig <- fig + ggplot2::theme(legend.position = general$legend_position)
  if (general$background == 'transparent') fig <- fig + ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = 'transparent', color = NA),
    panel.background = ggplot2::element_rect(fill = 'transparent', color = NA),
    legend.background = ggplot2::element_rect(fill = 'transparent', color = NA),
    legend.box.background = ggplot2::element_rect(fill = 'transparent', color = NA),
    legend.key = ggplot2::element_rect(fill = 'transparent', color = NA))
  if (!isTRUE(general$show_caption)) fig <- fig + ggplot2::labs(caption = NULL, subtitle = NULL)
  for (key in c('title', 'x_label', 'y_label')) if (nzchar(profile[[key]])) {
    label <- if (key == 'x_label') 'x' else if (key == 'y_label') 'y' else 'title'
    fig <- fig + do.call(ggplot2::labs, setNames(list(diff_label(profile[[key]], rank, comparison)), label))
  }
  if (no_hit) fig <- fig + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                                         axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  fig
}

save_differential_figure <- function(fig, stem, width = 12, height = 5, opts = NULL,
                                     profile = 'effect', rank = '', comparison = '') {
  if (!diff_target_wanted(opts,stem,profile)) return(character())
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  config <- diff_plot_config(opts, profile); general <- config$general
  font <- diff_plot_font(opts)
  fig <- style_differential_figure(fig, config, font, rank, comparison)
  if (!is.null(config$profile$width)) width <- config$profile$width
  if (!is.null(config$profile$height)) height <- config$profile$height
  dpi <- if (is.null(general$dpi)) 200 else general$dpi
  formats <- if (general$formats == 'png_svg_pdf') c('png', 'svg', 'pdf') else c('png', 'svg')
  files <- paste0(stem, '.', formats)
  for (file in files) {
    arguments <- list(filename = file, plot = fig, width = width, height = height, dpi = dpi, bg = general$background)
    if (endsWith(file, '.pdf')) {
      if (!capabilities('cairo')) stop('Cairo graphics are required for font-compatible PDF figures')
      arguments$device <- grDevices::cairo_pdf
    }
    do.call(ggplot2::ggsave, arguments)
  }
  effective <- list(profile = profile, width = width, height = height, dpi = dpi, formats = as.list(formats),
                    background = general$background, font = font,
                    general = general, profile_settings = config$profile,
                    labels = as.list(fig$labels[c('title', 'x', 'y', 'caption', 'subtitle')]),
                    rank = rank, comparison = comparison,
                    null_settings = 'Null sizing/font controls preserve the existing figure defaults')
  effective$legend_position <- fig$theme$legend.position
  if (profile == 'bubble' && is.null(effective$profile_settings$max_size)) effective$profile_settings$max_size <- 10
  if (profile == 'effect' && is.null(effective$profile_settings$point_size)) effective$profile_settings$point_size <- 2
  if (profile == 'heatmap' && is.null(effective$profile_settings$label_size)) effective$profile_settings$label_size <- 8
  jsonlite::write_json(list(requested = config$requested, effective = effective), paste0(stem, '_plot_settings.json'),
                       pretty = TRUE, auto_unbox = TRUE, null = 'null')
  attr(files, 'plot_settings') <- effective
  relative <- diff_relative_stem(opts,stem)
  differential_figures[[length(differential_figures)+1L]] <<- list(stem=relative,profile=profile,
                                                                  files=as.list(paste0(relative,'.',formats)))
  files
}

empty_differential_figure <- function(title, message, font) {
  fig <- ggplot2::ggplot() + ggplot2::annotate('text', x = 0, y = 0, label = paste(strwrap(message, 85), collapse = '\n'), family = font) +
    ggplot2::theme_void() + ggplot2::labs(title = title)
  attr(fig, 'differential_no_hit') <- TRUE
  fig
}

official_contrast_map <- function(meta, reference) {
  # Reproduce only the official contrast NAMES, never the fitted statistics.
  # RELEASE_3_20 .combn_fun names diff(c(beta_i,beta_j)) as coef_j_coef_i.
  groups <- levels(meta$Group)
  if (is.ordered(meta$Group) || !identical(groups[1], reference)) stop('Saved nominal model reference is inconsistent')
  model <- stats::model.matrix(~ Group, data = meta)
  coefficients <- colnames(model)[attr(model, 'assign') == 1L]
  nonreference <- groups[-1]
  if (length(coefficients) != length(nonreference)) stop('Cannot map the saved Group coefficients exactly')
  plan <- data.frame(key = coefficients, positive_group = nonreference,
                     negative_group = rep(reference, length(coefficients)), stringsAsFactors = FALSE)
  if (length(coefficients) > 1L) {
    pairs <- utils::combn(seq_along(coefficients), 2L)
    plan <- rbind(plan, data.frame(key = paste(coefficients[pairs[2, ]], coefficients[pairs[1, ]], sep = '_'),
                                  positive_group = nonreference[pairs[2, ]], negative_group = nonreference[pairs[1, ]],
                                  stringsAsFactors = FALSE))
  }
  if (anyDuplicated(plan$key)) stop('Original Group names produce ambiguous official contrast keys; no label or direction was guessed')
  plan
}

official_flag <- function(value) {
  if (is.logical(value)) return(!is.na(value) & value)
  if (is.numeric(value)) return(!is.na(value) & value == 1)
  !is.na(value) & value %in% c('TRUE', 'True', 'true', '1')
}

comparison_model_data <- function(res, contrast, opts, source) {
  if (!is.data.frame(res) || !'taxon' %in% names(res) || anyDuplicated(res$taxon)) stop('Saved official model results have invalid taxon identifiers')
  key <- contrast$key
  fields <- paste0(c('lfc_', 'se_', 'p_', 'q_', 'diff_'), key)
  if (!all(fields %in% names(res))) stop('Saved official comparison is unavailable: ', key)
  dat <- data.frame(taxon = as.character(res$taxon), lfc = res[[fields[1]]], se = res[[fields[2]]],
                    p = res[[fields[3]]], q = res[[fields[4]]], diff = res[[fields[5]]], stringsAsFactors = FALSE)
  if (!all(vapply(dat[c('lfc', 'se', 'p', 'q')], is.numeric, logical(1)))) stop('Official comparison estimates must be numeric')
  if (any(is.finite(dat$q) & (dat$q < 0 | dat$q > 1)) || any(is.finite(dat$se) & dat$se < 0))
    stop('Invalid saved official q-value or SE; no correction was substituted')
  ss_field <- paste0('passed_ss_', key); robust_field <- paste0('diff_robust_', key)
  dat$passed_ss <- if (ss_field %in% names(res)) res[[ss_field]] else rep(NA, nrow(dat))
  if (robust_field %in% names(res)) {
    robust <- official_flag(res[[robust_field]])
  } else {
    robust <- official_flag(dat$diff)
    if (ss_field %in% names(res)) robust <- robust & official_flag(dat$passed_ss)
    else if (isTRUE(opts$pseudo_sens)) stop('Saved sensitivity flag is missing for comparison: ', key)
  }
  alpha <- if (is.null(opts$alpha)) .05 else opts$alpha
  dat$robust_significant <- robust & is.finite(dat$q) & dat$q <= alpha
  dat$source <- rep(source, nrow(dat))
  dat$positive_group <- rep(contrast$positive_group, nrow(dat)); dat$negative_group <- rep(contrast$negative_group, nrow(dat))
  dat$structural_zero_group <- rep(NA_character_, nrow(dat))
  dat
}

structural_zero_data <- function(fit, contrast, counts, meta) {
  zero <- fit$zero_ind
  if (is.null(zero)) return(NULL)
  if (!is.data.frame(zero) || !'taxon' %in% names(zero) || anyDuplicated(zero$taxon) ||
      !all(as.character(zero$taxon) %in% rownames(counts))) stop('Saved structural-zero taxon identifiers do not match counts')
  # Official .get_struc_zero uses structural_zero (Group = ORIGINAL_LABEL).
  groups <- c(contrast$positive_group, contrast$negative_group)
  columns <- paste0('structural_zero (Group = ', groups, ')')
  if (!all(columns %in% names(zero))) stop('Saved structural-zero group columns do not match the original Groups')
  a <- zero[[columns[1]]]; b <- zero[[columns[2]]]
  explicit_false <- function(value) !is.na(value) & value %in% c(FALSE, 'FALSE', 'False', 'false', '0')
  hit <- (official_flag(a) & explicit_false(b)) | (explicit_false(a) & official_flag(b))
  ids <- as.character(zero$taxon[hit])
  if (!length(ids)) return(NULL)
  data.frame(taxon = ids, lfc = NA_real_, se = NA_real_, p = NA_real_, q = NA_real_, diff = NA,
             passed_ss = NA, robust_significant = FALSE, source = 'structural_zero_presence_absence',
             positive_group = contrast$positive_group, negative_group = contrast$negative_group,
             structural_zero_group = ifelse(official_flag(a[hit]), groups[1], groups[2]), stringsAsFactors = FALSE)
}

plot_actual_comparison <- function(dat, folder, contrast, style, opts) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  write_result(dat, file.path(folder, 'comparison.tsv'))
  model <- dat[dat$robust_significant & is.finite(dat$lfc) & is.finite(dat$se) & is.finite(dat$q), , drop = FALSE]
  model <- head(model[order(model$q, -abs(model$lfc), model$taxon), , drop = FALSE], opts$diff_plot_top_n)
  model$direction_group <- ifelse(model$lfc > 0, contrast$positive_group,
                                  ifelse(model$lfc < 0, contrast$negative_group, NA_character_))
  model$direction_group <- factor(model$direction_group, levels = style$group_order)
  write_result(model, file.path(folder, 'effect_plot_input.tsv'))
  title <- paste(contrast$positive_group, 'versus', contrast$negative_group)
  point <- diff_plot_config(opts, 'effect')$profile
  if (!nrow(model)) {
    fig <- empty_differential_figure(title, 'No robust significant taxa with finite model effects. Structural-zero findings, if present, have no model LFC/q and are shown in abundance plots and comparison.tsv.', opts$font)
  } else {
    model$taxon <- factor(model$taxon, levels = rev(model$taxon))
    fig <- ggplot2::ggplot(model, ggplot2::aes(lfc, taxon, color = direction_group)) +
      ggplot2::geom_vline(xintercept = 0, color = 'grey70')
    if (isTRUE(point$show_errorbars)) fig <- fig + ggplot2::geom_segment(ggplot2::aes(x = lfc - se, xend = lfc + se, yend = taxon))
    fig <- fig + ggplot2::geom_point(size = if (is.null(point$point_size)) 2 else point$point_size) +
      ggplot2::scale_color_manual(values = style$group_colors, breaks = style$group_order, na.value = '#777777') +
      ggplot2::labs(title = title, x = if (isTRUE(point$show_errorbars)) 'Official ANCOM-BC2 log fold change (± SE)' else 'Official ANCOM-BC2 log fold change', y = NULL, color = 'Effect direction',
                    caption = paste(strwrap(paste0('Positive = ', contrast$positive_group, ' − ', contrast$negative_group,
                      '. Robust significant model taxa only; ranked by the original comparison q. Structural-zero findings are separate; no LFC/q was invented.'), 105), collapse = '\n')) +
      ggplot2::theme_bw(base_family = opts$font) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
  }
  files <- save_differential_figure(fig, file.path(folder, 'effect'), height = max(4, nrow(model) * .25 + 1),
                                     opts = opts, profile = 'effect', rank = if (is.null(opts$plot_rank)) '' else opts$plot_rank, comparison = title)
  info <- list(contrast = contrast$key, positive_group = contrast$positive_group, negative_group = contrast$negative_group,
               sign = paste(contrast$positive_group, 'minus', contrast$negative_group),
               model_result_source = contrast$result_source, model_robust_taxa = sum(dat$robust_significant),
               structural_zero_findings = sum(dat$source == 'structural_zero_presence_absence'),
               effect_shown = nrow(model), q_recomputed = FALSE, files = files,
               plot_settings = attr(files, 'plot_settings'))
  jsonlite::write_json(info, file.path(folder, 'comparison.json'), pretty = TRUE, auto_unbox = TRUE)
  info
}

plot_abundance_views <- function(counts, meta, findings, folder, opts, title) {
  want_heat <- diff_target_wanted(opts,file.path(folder,'abundance_heatmap'),'heatmap')
  want_bubble <- diff_target_wanted(opts,file.path(folder,'abundance_bubble'),'bubble')
  if (!want_heat && !want_bubble) return(list(plots=character()))
  opts$font <- diff_plot_font(opts)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  if (!identical(colnames(counts), rownames(meta)) || any(colSums(counts) <= 0)) stop('Abundance counts and actual sample metadata must align with nonzero libraries')
  style <- effect_group_style(levels(meta$Group), opts$group_style)
  heat_config <- diff_plot_config(opts, 'heatmap'); bubble_config <- diff_plot_config(opts, 'bubble')
  relative <- sweep(counts, 2L, colSums(counts), '/')
  means <- vapply(style$group_order, function(group) rowMeans(relative[, as.character(meta$Group) == group, drop = FALSE]), numeric(nrow(relative)))
  if (!is.matrix(means)) means <- matrix(means, nrow = nrow(relative))
  dimnames(means) <- list(rownames(counts), style$group_order)
  write_result(relative, file.path(folder, 'sample_relative_abundance.tsv'))
  write_result(means, file.path(folder, 'group_mean_relative_abundance.tsv'))
  if (nrow(findings) && !all(findings$taxon %in% rownames(counts))) stop('Official differential findings do not match original count taxon IDs')
  model_ids <- unique(findings$taxon[findings$source %in% c('res', 'res_pair')])
  omnibus_ids <- unique(findings$taxon[findings$source == 'qiime_ancom'])
  zero_ids <- unique(findings$taxon[findings$source == 'structural_zero_presence_absence'])
  selection <- data.frame(taxon = rownames(counts), mean_relative_abundance = rowMeans(relative),
                           model_robust = rownames(counts) %in% model_ids,
                           omnibus_reject = rownames(counts) %in% omnibus_ids,
                           structural_zero = rownames(counts) %in% zero_ids, stringsAsFactors = FALSE)
  selection$selected_union <- selection$model_robust | selection$structural_zero | selection$omnibus_reject
  selected <- selection[selection$selected_union, , drop = FALSE]
  selected <- head(selected[order(-selected$mean_relative_abundance, selected$taxon), , drop = FALSE], opts$diff_plot_top_n)
  selection$shown <- selection$taxon %in% selected$taxon
  selection$source <- ifelse(selection$model_robust & selection$structural_zero, 'model_robust+structural_zero',
                      ifelse(selection$model_robust, 'model_robust', ifelse(selection$structural_zero, 'structural_zero_presence_absence',
                      ifelse(selection$omnibus_reject, 'qiime_ancom_omnibus_reject', 'not_selected'))))
  write_result(selection, file.path(folder, 'taxon_selection.tsv'))
  write_result(findings, file.path(folder, 'selection_findings.tsv'))
  note <- 'Each saved input sample is normalized by its own full unrarefied collapsed-count total; group abundance is the arithmetic mean of sample proportions. Abundance is descriptive and is not the bias-corrected model LFC.'
  zero_note <- if (any(selection$structural_zero)) '[SZ] marks an official structural-zero presence/absence finding in at least one requested pair. Its LFC/SE/p/q are not estimated; known TRUE vs FALSE is required, not TRUE vs TRUE or missing flags.' else ''
  if (!nrow(selected)) {
    empty <- data.frame(taxon = character(), Group = character(), mean_relative_abundance = numeric(), mean_percent = numeric(), row_zscore = numeric())
    if (want_heat) {
      write_result(empty,file.path(folder,'heatmap_plot_input.tsv'))
      heat <- empty_differential_figure(paste(title, 'abundance heatmap'), 'No official differential findings selected. No nonsignificant or synthetic taxa were substituted.', opts$font)
    }
    if (want_bubble) {
      write_result(empty,file.path(folder,'bubble_plot_input.tsv'))
      bubble <- empty_differential_figure(paste(title, 'abundance bubble plot'), 'No official differential findings selected. No nonsignificant or synthetic taxa were substituted.', opts$font)
    }
  } else {
    values <- means[selected$taxon, , drop = FALSE]; sd <- apply(values, 1L, stats::sd)
    z <- sweep(sweep(values, 1L, rowMeans(values), '-'), 1L, ifelse(is.finite(sd) & sd > 0, sd, 1), '/')
    labels <- paste0(selected$taxon, ifelse(selected$structural_zero, ' [SZ]', ''))
    labels <- vapply(labels, function(x) paste(strwrap(gsub(';[[:space:]]*', '; ', x), width = 80), collapse = '\n'), character(1))
    long <- data.frame(taxon = rep(selected$taxon, times = ncol(values)),
                       Group = rep(colnames(values), each = nrow(values)),
                       mean_relative_abundance = as.vector(values), mean_percent = as.vector(values) * 100,
                       row_zscore = as.vector(z), plot_y = rep(rev(seq_len(nrow(values))), times = ncol(values)),
                       structural_zero = rep(selected$structural_zero, times = ncol(values)), stringsAsFactors = FALSE)
    long$Group <- factor(long$Group, levels = style$group_order)
    long$heatmap_value <- if (opts$diff_heatmap_scale == 'row_zscore') long$row_zscore else long$mean_percent
    heat_order <- if (heat_config$profile$row_order == 'alphabetical') sort(selected$taxon) else selected$taxon
    heat_long <- long; heat_long$plot_y <- length(heat_order) - match(heat_long$taxon, heat_order) + 1L
    heat_labels <- labels[match(heat_order, selected$taxon)]
    if (want_heat) write_result(heat_long,file.path(folder,'heatmap_plot_input.tsv'))
    if (want_bubble) write_result(long,file.path(folder,'bubble_plot_input.tsv'))
    markers <- data.frame(Group = factor(style$group_order, levels = style$group_order), plot_y = .15)
    common <- ggplot2::theme_bw(base_family = opts$font) + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 8), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    y_scale <- function(axis_labels = labels) ggplot2::scale_y_continuous(breaks = seq_len(nrow(values)), labels = rev(axis_labels),
                                                    limits = c(0, nrow(values) + .5), expand = c(0, 0))
    value_size <- if (is.null(heat_config$general$base_size)) 3 else heat_config$general$base_size / (72.27 / 25.4)
    if (want_heat) {
    heat_label <- if (opts$diff_heatmap_scale == 'row_zscore') 'Row Z-score of group means' else 'Group mean relative abundance (%)'
    limit <- if (opts$diff_heatmap_scale == 'row_zscore') max(1, abs(long$row_zscore)) else max(long$mean_percent)
    midpoint <- if (opts$diff_heatmap_scale == 'row_zscore') 0 else limit / 2
    heat <- ggplot2::ggplot(heat_long, ggplot2::aes(Group, plot_y, fill = heatmap_value)) + ggplot2::geom_tile(width = .95, height = .95) +
      ggplot2::geom_point(data = markers, ggplot2::aes(Group, plot_y, color = Group), inherit.aes = FALSE, shape = 15, size = 3) +
      ggplot2::scale_fill_gradient2(low = opts$heatmap_low, mid = opts$heatmap_mid, high = opts$heatmap_high,
         midpoint = midpoint, limits = if (opts$diff_heatmap_scale == 'row_zscore') c(-limit, limit) else c(0, limit)) +
      ggplot2::scale_color_manual(values = style$group_colors, breaks = style$group_order) + y_scale(heat_labels) + common +
      ggplot2::labs(title = paste(title, 'abundance heatmap'), x = NULL, y = NULL, fill = heat_label, color = 'Group',
                    caption = paste(strwrap(paste(note, zero_note, if (opts$diff_heatmap_scale == 'row_zscore') 'Constant rows have Z=0; Z-scores do not represent LFC or significance.' else ''), 110), collapse = '\n'))
    if (isTRUE(heat_config$profile$show_values)) heat <- heat + ggplot2::geom_text(ggplot2::aes(label = sprintf('%.2f%%', mean_percent)),
                                                                               size = value_size, family = opts$font, color = 'black')
    }
    if (want_bubble) {
    bubble <- ggplot2::ggplot(long, ggplot2::aes(Group, plot_y, size = mean_percent))
    if (opts$bubble_color_mode == 'group') {
      bubble <- bubble + ggplot2::geom_point(ggplot2::aes(color = Group), alpha = bubble_config$profile$point_alpha, stroke = 0) +
        ggplot2::scale_color_manual(values = style$group_colors, breaks = style$group_order) +
        ggplot2::geom_point(data = markers, ggplot2::aes(Group, plot_y, color = Group), inherit.aes = FALSE, shape = 15, size = 3) +
        ggplot2::labs(color = 'Group')
    } else {
      bubble <- bubble + ggplot2::geom_point(ggplot2::aes(color = mean_percent), alpha = bubble_config$profile$point_alpha, stroke = 0) +
        ggplot2::scale_color_gradient(low = opts$bubble_low, high = opts$bubble_high) +
        ggplot2::geom_point(data = markers, ggplot2::aes(Group, plot_y, fill = Group), inherit.aes = FALSE, shape = 21, size = 3, color = 'grey30') +
        ggplot2::scale_fill_manual(values = style$group_colors, breaks = style$group_order) +
        ggplot2::labs(color = 'Group mean (%)', fill = 'Group')
    }
    size_breaks <- unique(signif(max(long$mean_percent) * c(.25, .5, 1), 2))
    short_bubble <- !is.null(bubble_config$profile$height) && bubble_config$profile$height < 6.5 &&
      bubble_config$general$legend_position == 'original'
    bottom_legend <- bubble_config$general$legend_position == 'bottom' || short_bubble
    bubble <- bubble + ggplot2::scale_size_area(max_size = if (is.null(bubble_config$profile$max_size)) 10 else bubble_config$profile$max_size,
                                               breaks = size_breaks) + y_scale() + common +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 8), legend.title = ggplot2::element_text(size = 9),
                    legend.key.height = grid::unit(.25, 'cm'), legend.spacing.y = grid::unit(.25, 'cm'),
                    plot.margin = ggplot2::margin(12, 12, 12, 12), plot.caption = ggplot2::element_text(size = 8)) +
      ggplot2::guides(size = ggplot2::guide_legend(order = 1, nrow = if (bottom_legend) 1 else NULL, title.position = 'top'),
                      fill = ggplot2::guide_legend(order = 2, nrow = if (bottom_legend) 1 else NULL, title.position = 'top'),
                      color = if (opts$bubble_color_mode == 'group') ggplot2::guide_legend(order = 2, nrow = if (bottom_legend) 1 else NULL, title.position = 'top') else
                        ggplot2::guide_colorbar(order = 3, direction = if (bottom_legend) 'horizontal' else 'vertical',
                          title.position = 'top',
                          barheight = grid::unit(if (bottom_legend) .4 else 2.8, 'cm'),
                          barwidth = grid::unit(if (bottom_legend) 2.8 else .4, 'cm'))) +
      ggplot2::labs(title = paste(title, 'abundance bubble plot'), x = NULL, y = NULL, size = 'Group mean (%)',
                    caption = paste(strwrap(paste(note, zero_note), 110), collapse = '\n'))
    if (short_bubble) bubble <- bubble + ggplot2::theme(legend.position = 'bottom')
    if (isTRUE(bubble_config$profile$show_values)) bubble <- bubble + ggplot2::geom_text(
      ggplot2::aes(label = sprintf('%.2f%%', mean_percent)), size = value_size, family = opts$font, color = 'black')
    }
  }
  lines <- if (nrow(selected)) max(lengths(strsplit(labels, '\n', fixed = TRUE))) else 1L
  height <- max(5, nrow(selected) * (.17 * lines + .12) + 2)
  rank <- if (is.null(opts$plot_rank)) '' else opts$plot_rank
  heat_files <- if (want_heat) save_differential_figure(heat, file.path(folder, 'abundance_heatmap'), height = height,
                                        opts = opts, profile = 'heatmap', rank = rank, comparison = 'all groups') else character()
  bubble_files <- if (want_bubble) save_differential_figure(bubble, file.path(folder, 'abundance_bubble'), width = 13, height = max(8, height),
                                          opts = opts, profile = 'bubble', rank = rank, comparison = 'all groups') else character()
  files <- c(heat_files, bubble_files)
  summary <- list(selected_total = sum(selection$selected_union), model_robust_total = sum(selection$model_robust),
                   omnibus_reject_total = sum(selection$omnibus_reject),
                   structural_zero_total = sum(selection$structural_zero), shown = nrow(selected),
                   selection = 'Union of actual robust model differential findings and pair-specific official structural-zero findings; top N by mean sample relative abundance',
                   table_units = 'relative abundance proportions (0 to 1); plot input mean_percent is percent',
                   heatmap_scale = opts$diff_heatmap_scale, bubble_color_mode = opts$bubble_color_mode,
                   group_style = list(group_order = as.list(style$group_order), group_colors = as.list(style$group_colors)),
                   visual_options = visual_options(opts), notes = c(note, zero_note), plots = files,
                   plot_settings = list(heatmap = attr(heat_files, 'plot_settings'), bubble = attr(bubble_files, 'plot_settings')))
  jsonlite::write_json(summary, file.path(folder, 'plot_summary.json'), pretty = TRUE, auto_unbox = TRUE)
  summary
}

plot_ancom_level <- function(fit, counts, meta, folder, opts, title, rank = '') {
  opts <- modifyList(opts, visual_options(opts)); if (is.null(opts$font)) opts$font <- 'sans'
  opts$font <- diff_plot_font(opts)
  opts$plot_rank <- rank
  plan <- official_contrast_map(meta, opts$reference_group)
  if (opts$diff_comparisons == 'reference') {
    plan <- plan[plan$negative_group == opts$reference_group, , drop = FALSE]
    result <- fit$res; result_source <- 'res'
  } else if (nlevels(meta$Group) == 2L) {
    result <- fit$res; result_source <- 'res'
  } else {
    if (is.null(fit$res_pair)) stop('Saved analysis has no official all-pairwise results. Select reference plots or rerun ANCOM-BC2 with pairwise tests; plots cannot create these statistics')
    result <- fit$res_pair; result_source <- 'res_pair'
  }
  style <- effect_group_style(levels(meta$Group), opts$group_style)
  plan <- plan[order(match(plan$positive_group, style$group_order), match(plan$negative_group, style$group_order)), , drop = FALSE]
  plan$result_source <- result_source
  comparisons <- list(); findings <- list()
  pair_indices <- seq_len(nrow(plan)); target_pair <- NULL
  if (!is.null(opts$plot_target) && identical(opts$plot_target$profile,'effect')) {
    comparison_key <- opts$plot_target$comparison_key
    if (is.null(comparison_key)) {
      selected_info <- file.path(opts$analysis_run,dirname(opts$plot_target$stem),'comparison.json')
      if (!file.exists(selected_info)) stop('The original comparison identity file is missing; refusing to substitute a different comparison')
      comparison_key <- jsonlite::fromJSON(selected_info,simplifyVector=FALSE)$contrast
    }
    pair_indices <- which(plan$key==comparison_key)
    if (length(pair_indices)!=1L) stop('Selected original comparison is not present in the saved official result plan')
    target_pair <- basename(dirname(opts$plot_target$stem))
  }
  for (i in pair_indices) {
    contrast <- as.list(plan[i, , drop = FALSE])
    dat <- comparison_model_data(result, contrast, opts, result_source)
    zero <- structural_zero_data(fit, contrast, counts, meta)
    if (!is.null(zero)) dat <- rbind(dat, zero)
    pair_folder <- file.path(folder,'comparisons',if (is.null(target_pair)) sprintf('pair_%03d',i) else target_pair)
    if (diff_target_wanted(opts,file.path(pair_folder,'effect'),'effect'))
      comparisons[[length(comparisons)+1L]] <- plot_actual_comparison(dat,pair_folder,contrast,style,opts)
    selected <- dat$robust_significant | dat$source == 'structural_zero_presence_absence'
    findings[[i]] <- data.frame(taxon = dat$taxon[selected], source = dat$source[selected],
                                comparison = rep(contrast$key, sum(selected)), stringsAsFactors = FALSE)
  }
  found <- if (length(findings)) do.call(rbind, findings) else data.frame(taxon = character(), source = character(), comparison = character())
  abundance <- plot_abundance_views(counts, meta, found, file.path(folder, 'abundance'), opts, title)
  list(status = 'rendered', comparisons_mode = opts$diff_comparisons, result_source = result_source,
       two_group_note = if (nlevels(meta$Group) == 2L) 'The official package disables multi-group pairwise for two groups; the primary result is the single actual comparison' else NULL,
       comparisons = comparisons, abundance = abundance,
       plots = c(unlist(lapply(comparisons, `[[`, 'files'), use.names = FALSE), abundance$plots),
       statistics_refit = FALSE, q_recomputed = FALSE)
}

plot_qiime <- function(opts, run, source_run = run) {
  opts$plot_output_run <- run
  opts <- diff_target_tables(opts)
  opts <- modifyList(opts, visual_options(opts)); if (is.null(opts$font)) opts$font <- 'sans'
  opts$font <- diff_plot_font(opts)
  summaries <- list()
  for (entry in opts$tables) {
    opts$plot_rank <- entry$name
    folder <- file.path(run, '02_differential', paste0('L', entry$level, '_', entry$name))
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    original <- file.path(source_run, '02_differential', paste0('L', entry$level, '_', entry$name), 'ancom_results.tsv')
    result <- read.delim(original, check.names = FALSE,
                         stringsAsFactors = FALSE, colClasses = 'character')
    if (!all(c('W', 'Reject null hypothesis') %in% names(result))) stop('Official QIIME ANCOM W/decision columns missing')
    result$W <- suppressWarnings(as.numeric(result$W))
    if (any(!is.finite(result$W))) stop('QIIME ANCOM W must be finite in official export')
    result$taxon <- result[[1]]
    result$decision <- factor(result[['Reject null hypothesis']], levels = c('False', 'True', 'FALSE', 'TRUE'))
    files <- character()
    if (diff_target_wanted(opts,file.path(folder,'ancom_W'),'effect')) {
    dat <- head(result[order(-result$W), , drop = FALSE], 25L)
    dat$taxon <- factor(dat$taxon, levels = rev(dat$taxon))
    fig <- ggplot2::ggplot(dat, ggplot2::aes(x = W, y = taxon, fill = decision)) +
      ggplot2::geom_col() + ggplot2::labs(x = 'Official QIIME 2 ANCOM W', y = NULL,
                                        fill = 'Reject null', title = paste('QIIME 2 ANCOM:', entry$name)) +
      ggplot2::theme_bw(base_family = opts$font) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    files <- save_differential_figure(fig, file.path(folder, 'ancom_W'), height = max(4, nrow(dat) * .25),
                                      opts = opts, profile = 'effect', rank = entry$name, comparison = 'omnibus ANCOM')
    }
    counts <- read_count_input(entry$counts_tsv)
    meta <- read_group_metadata(opts$metadata_tsv, colnames(counts), opts$reference_group, unlist(opts$groups))
    if (anyDuplicated(result$taxon) || !all(result$taxon %in% rownames(counts))) stop('Official QIIME ANCOM taxon IDs do not match original counts')
    hit <- official_flag(result[['Reject null hypothesis']])
    selected <- data.frame(taxon = result$taxon[hit], source = rep('qiime_ancom', sum(hit)),
                            comparison = rep('omnibus_ANCOM_reject_null', sum(hit)), stringsAsFactors = FALSE)
    abundance <- plot_abundance_views(counts, meta, selected, file.path(folder, 'abundance'), opts, paste('QIIME 2 ANCOM:', entry$name))
    abundance$selection <- 'Actual official omnibus ANCOM Reject null TRUE taxa; no pairwise inference or invented p/q values'
    if (length(abundance$plots)) jsonlite::write_json(abundance, file.path(folder, 'abundance', 'plot_summary.json'), pretty = TRUE, auto_unbox = TRUE)
    summaries[[length(summaries) + 1L]] <- list(level = entry$level, name = entry$name, status = 'rendered',
                                               comparisons = 'omnibus_ANCOM_only', abundance = abundance,
                                               plots = c(files, abundance$plots), statistics_refit = FALSE,
                                               plot_settings = attr(files, 'plot_settings'))
  }
  summary <- list(method = 'QIIME 2 ANCOM', status = 'rendered', levels = summaries,
                   statistics_refit = FALSE, original_analysis = normalizePath(source_run, mustWork = TRUE))
  jsonlite::write_json(summary, file.path(run, '02_differential', 'summary.json'), pretty = TRUE, auto_unbox = TRUE)
  summary
}

replot_saved_differential <- function(request, run) {
  if (is.null(request$analysis_run)) stop('The original differential analysis run must be provided')
  source <- normalizePath(request$analysis_run, mustWork = TRUE)
  if (identical(source, run)) stop('Cached differential figures must use a separate new output run')
  saved <- jsonlite::fromJSON(file.path(source, 'differential_options.json'), simplifyVector = FALSE)
  # Statistical controls, original Group/reference and actual count paths come only from the saved analysis.
  opts <- modifyList(saved, visual_options(request))
  if (!is.null(request$group_style)) opts$group_style <- request$group_style
  if (!is.null(request$font)) opts$font <- request$font
  if (!is.null(request$plot_settings)) opts$plot_settings <- modifyList(
    if (is.null(saved$plot_settings)) list() else saved$plot_settings, request$plot_settings, keep.null = TRUE)
  opts$plot_target <- request$plot_target
  opts$plot_output_run <- run
  opts$analysis_run <- source
  opts <- diff_target_tables(opts)
  if (identical(saved$method, 'qiime_ancom')) return(plot_qiime(opts, run, source))
  if (!is.null(saved$method) && !identical(saved$method, 'ancombc2')) stop('Unsupported saved differential method')
  summaries <- list()
  for (entry in opts$tables) {
    name <- paste0('L', entry$level, '_', entry$name)
    folder <- file.path(run, '02_differential', name); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    result <- tryCatch({
      fitfile <- file.path(source, '02_differential', name, 'ancombc2_fit.rds')
      if (!file.exists(fitfile)) stop('Saved model fit missing; original failed levels cannot be replotted: ', name)
      fit <- readRDS(fitfile)
      counts <- read_count_input(entry$counts_tsv)
      meta <- read_group_metadata(opts$metadata_tsv, colnames(counts), opts$reference_group, unlist(opts$groups))
      for (key in c('res', 'res_global', 'res_pair', 'zero_ind', 'ss_tab'))
        write_result(fit[[key]], file.path(folder, paste0(key, '.tsv')))
      report <- plot_ancom_level(fit, counts, meta, folder, opts, paste('ANCOM-BC2:', entry$name), rank = entry$name)
      report$level <- entry$level; report$name <- entry$name; report$original_fit <- fitfile
      report
    }, error = function(e) list(level = entry$level, name = entry$name, status = 'failed', error = conditionMessage(e), statistics_refit = FALSE))
    jsonlite::write_json(result, file.path(folder, 'diagnostics.json'), pretty = TRUE, auto_unbox = TRUE)
    summaries[[length(summaries) + 1L]] <- result
  }
  failed <- sum(vapply(summaries, function(x) identical(x$status, 'failed'), logical(1)))
  summary <- list(method = 'ANCOM-BC2', status = if (failed) 'partial_failure' else 'rendered', failed_levels = failed,
                   original_analysis = source, statistics_refit = FALSE, q_recomputed = FALSE, levels = summaries)
  jsonlite::write_json(summary, file.path(run, '02_differential', 'summary.json'), pretty = TRUE, auto_unbox = TRUE)
  if (failed) stop(failed, ' cached differential level(s) could not be plotted; actual results and diagnostics were preserved')
  summary
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L || !args[1] %in% c('--check', '--run', '--plot-qiime', '--plot-ancom'))
    stop('Usage: differential_ancombc.R --check|--run|--plot-qiime|--plot-ancom options.json run_directory')
  plotting_only <- args[1] %in% c('--plot-qiime', '--plot-ancom')
  required <- if (plotting_only) c('jsonlite', 'ggplot2', 'svglite') else
    c('ANCOMBC', 'jsonlite', 'ggplot2', 'svglite')
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (any(!available) && plotting_only) stop(paste('Missing isolated R plotting dependencies:',
                    paste(required[!available], collapse = ', '), '; run 安装依赖.sh for the selected plotting r_env.'))
  if (any(!available)) stop(paste('Missing isolated ANCOM-BC2 dependencies:',
                                  paste(required[!available], collapse = ', '),
                                  '; run install_differential.sh in Ubuntu/WSL, or select a prepared ancom_env.'))
  if (!plotting_only && !'ancombc2' %in% getNamespaceExports('ANCOMBC')) stop('Installed ANCOMBC does not provide ancombc2')
  opts <- jsonlite::fromJSON(args[2], simplifyVector = FALSE)
  run <- normalizePath(args[3], mustWork = TRUE)
  versions <- list(R = R.version.string, packages = as.list(vapply(required,
                    function(pkg) as.character(utils::packageVersion(pkg)), character(1))))
  jsonlite::write_json(versions, file.path(run, 'software_versions.json'), pretty = TRUE, auto_unbox = TRUE)
  if (args[1] == '--plot-qiime') { plot_qiime(opts, run); write_diff_figure_manifest(run,opts$plot_target); return(invisible(TRUE)) }
  if (args[1] == '--plot-ancom') { replot_saved_differential(opts, run); write_diff_figure_manifest(run,opts$plot_target); return(invisible(TRUE)) }
  if (args[1] == '--check') return(invisible(TRUE))
  opts$plot_output_run <- run
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
    plot_error <- NULL
    visual <- tryCatch(plot_ancom_level(fit, counts, meta, folder, opts, paste('ANCOM-BC2:', entry$name), rank = entry$name),
                       error = function(e) { plot_error <<- conditionMessage(e); NULL })
    if (!is.null(plot_error)) message('L', entry$level, ' model results preserved; plotting failed: ', plot_error)
    diagnostics <- list(level = entry$level, method = 'ANCOM-BC2', status = if (is.null(plot_error)) 'executed' else 'plot_failed', started_at = started,
                         finished_at = format(Sys.time(), tz = 'UTC', usetz = TRUE), warnings = unique(warnings),
                         convergence = 'Not inferred; use actual package warnings, saved fit and workflow.log',
                         missing_values = 'NA is preserved in all exported result tables',
                         primary_tests = 'Nonreference nominal Group levels versus explicit reference',
                         global = opts$global, pairwise = opts$pairwise,
                         pairwise_adjustment = 'Official ANCOM-BC2 mixed directional FDR with global screen and Holm within-taxon control',
                         structural_zeros = 'Reported separately by official zero_ind; no artificial finite p/q values',
                         feature_count = nrow(counts), modeled_features = length(fit_taxa),
                         visualizations = visual, plot_error = plot_error,
                         plots = if (is.null(visual)) character() else visual$plots)
    jsonlite::write_json(diagnostics, file.path(folder, 'diagnostics.json'), pretty = TRUE, auto_unbox = TRUE)
    summaries[[length(summaries) + 1L]] <- diagnostics
  }
  failed <- sum(vapply(summaries, function(x) x$status %in% c('failed', 'plot_failed'), logical(1)))
  jsonlite::write_json(list(method = 'ANCOM-BC2', status = if (failed) 'partial_failure' else 'executed',
                            failed_levels = failed, levels = summaries), file.path(destination, 'summary.json'),
                        pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(run, 'R_session_info.txt'))
  write_diff_figure_manifest(run)
  if (failed) stop(paste(failed, 'taxonomy level(s) failed; actual outputs and failure diagnostics are preserved'))
}

if (sys.nframe() == 0L) main()
