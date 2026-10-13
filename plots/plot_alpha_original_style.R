# source("qiime2_io.R") before sourcing this file.
# Original alpha plot aesthetics; QIIME alpha values are not recalculated.

alpha_kw <- function(value, group) {
  if (length(unique(value)) < 2L) return(c(H = NA_real_, P = NA_real_))
  fit <- stats::kruskal.test(value ~ droplevels(factor(group)))
  c(H = unname(fit$statistic), P = fit$p.value)
}

alpha_r_stats <- function(data) {
  global <- alpha_kw(data$Value, data$Group)
  groups <- utils::combn(levels(data$Group), 2)
  pair <- do.call(rbind, lapply(seq_len(ncol(groups)), function(j) {
    g <- groups[, j]
    selected <- data$Group %in% g
    test <- alpha_kw(data$Value[selected], data$Group[selected])
    data.frame(Group1 = g[1], Group2 = g[2], H = unname(test["H"]),
                P = unname(test["P"]), stringsAsFactors = FALSE)
  }))
  pair$q_within_metric <- NA_real_
  valid <- is.finite(pair$P)
  pair$q_within_metric[valid] <- stats::p.adjust(pair$P[valid], method = "BH")
  pair$Status <- ifelse(valid, "ok", "all_values_identical")
  list(overall = data.frame(KW_H = unname(global["H"]), KW_p = unname(global["P"]),
                             N = nrow(data)), pairwise = pair)
}

original_alpha_plot <- function(data, annotation,
                                 font_family = "Times New Roman",
                                 group_colors = c(CON = "#D55E00", LDC = "#0072B2",
                                                  HDC = "#CC79A7"),
                                 display_names = NULL, seed = 123L,
                                 point_size = 1.5, box_width = 0.6,
                                 jitter_width = 0.2, show_points = TRUE,
                                 show_letters = TRUE, show_outliers = FALSE,
                                 plot_type = "boxplot") {
  labels <- setNames(levels(data$Metric), levels(data$Metric))
  if (!is.null(display_names)) {
    known <- intersect(names(display_names), names(labels))
    labels[known] <- display_names[known]
  }
  p <- ggplot2::ggplot(data, ggplot2::aes(x = Group, y = Value, color = Group))
  p <- p + if (plot_type == "violin") {
    ggplot2::geom_violin(fill = NA, width = box_width)
  } else {
    ggplot2::geom_boxplot(fill = NA, outlier.shape = if (show_outliers) 19 else NA,
                          width = box_width)
  }
  if (show_points) p <- p + ggplot2::geom_jitter(
    position = ggplot2::position_jitter(width = jitter_width, seed = seed),
    size = point_size, alpha = 1)
  p <- p +
    ggplot2::facet_wrap(~ Metric, nrow = 1, scales = "free_y",
                        labeller = ggplot2::labeller(Metric = labels))
  if (show_letters) p <- p + ggplot2::geom_text(
      data = annotation, ggplot2::aes(x = Group, y = y, label = Letter),
      inherit.aes = FALSE, vjust = 0, family = font_family,
      fontface = "plain", size = 8
    )
  p +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "", x = "", y = "") +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 24, family = font_family),
      strip.background = ggplot2::element_rect(fill = "gray80", color = "black"),
      axis.title.y = ggplot2::element_text(size = 18, face = "bold", family = font_family),
      axis.text.x = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5,
                                         family = font_family),
      axis.text.y = ggplot2::element_text(size = 20, face = "bold", family = font_family),
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 1),
      legend.position = "none"
    ) +
    ggplot2::scale_color_manual(values = group_colors)
}

run_alpha_original_style <- function(
    project = "D:/qiime2/肥育猪试验-终版", input_csv = NULL,
    core_dir = "core-metrics-results", metadata_file = "metadata.txt",
    qiime_stats_dir = "alpha", stats_source = c("qiime2", "r"),
    output_dir = "R_original_style/alpha", metrics = NULL, metric_keys = NULL,
    vector_files = NULL, display_names = NULL,
    group_col = "Group", group_levels = c("CON", "LDC", "HDC"),
    group_colors = c(CON = "#D55E00", LDC = "#0072B2", HDC = "#CC79A7"),
    font_family = "Times New Roman", significance_level = 0.05, seed = 123L) {
  q2_packages(c("ggplot2", "svglite", "systemfonts", "multcompView",
                if (match.arg(stats_source) == "qiime2") "jsonlite"))
  stats_source <- match.arg(stats_source)
  font_family <- q2_font(font_family)
  if (!dir.exists(project)) stop("工作目录不存在：", project)
  resolve <- function(p) q2_resolve(project, p)
  if (!is.finite(significance_level) || significance_level <= 0 ||
      significance_level >= 1) stop("significance_level 必须位于 0 和 1 之间")
  if (!is.null(input_csv)) {
    df <- utils::read.csv(resolve(input_csv), check.names = FALSE,
                           colClasses = "character", fileEncoding = "UTF-8-BOM")
    if (!all(c("Sample", "Group") %in% names(df))) {
      stop("alpha CSV 需要 Sample、Group 两列，其余指定列为 alpha 指标")
    }
    if (is.null(metrics)) metrics <- setdiff(names(df), c("Sample", "Group"))
    if (!all(metrics %in% names(df))) stop("CSV 缺少指定指标")
    df <- df[c("Sample", "Group", metrics)]
  } else {
    if (is.null(metrics)) metrics <- c("observed_features", "shannon", "evenness", "faith_pd")
    metadata <- q2_read_metadata(resolve(metadata_file))
    if (!group_col %in% names(metadata)) stop("元数据没有分组列：", group_col)
    df <- NULL
    for (metric in metrics) {
      key <- if (!is.null(metric_keys) && metric %in% names(metric_keys))
        metric_keys[[metric]] else metric
      path <- if (!is.null(vector_files) && metric %in% names(vector_files))
        resolve(vector_files[[metric]]) else
        resolve(file.path(core_dir, paste0(key, "_vector.qza")))
      v <- q2_read_alpha(path)
      if (is.null(df)) {
        if (!all(v$Sample %in% metadata$SampleID)) stop("alpha 样本缺少元数据")
        df <- data.frame(Sample = v$Sample,
                          Group = metadata[[group_col]][match(v$Sample, metadata$SampleID)],
                          stringsAsFactors = FALSE)
      }
      if (!setequal(v$Sample, df$Sample)) stop("alpha 指标之间保留的样本不一致")
      df[[metric]] <- v$Value[match(df$Sample, v$Sample)]
    }
  }
  if (!length(metrics) || anyDuplicated(metrics)) stop("至少选择一个唯一的 alpha 指标")
  if (anyNA(df$Sample) || any(!nzchar(df$Sample)) || anyDuplicated(df$Sample)) {
    stop("Sample 不能为空或重复；样本 ID 必须唯一")
  }
  if (anyNA(df$Group) || any(!df$Group %in% group_levels)) stop("输入 Group 缺失或不在预设组中")
  df$Group <- factor(df$Group, levels = group_levels)
  if (any(table(df$Group) == 0L)) stop("某个预设组没有样本，请核对抽平后保留的样本")
  if (length(group_levels) < 2L || anyDuplicated(group_levels)) stop("需要至少两个不同的组")
  if (!all(group_levels %in% names(group_colors))) stop("请为每个组设置颜色")
  for (metric in metrics) {
    df[[metric]] <- suppressWarnings(as.numeric(df[[metric]]))
    if (any(!is.finite(df[[metric]]))) stop("alpha 指标存在缺失/非数值：", metric)
  }
  long <- data.frame(
    Sample = rep(df$Sample, each = length(metrics)),
    Group = rep(df$Group, each = length(metrics)),
    Metric = factor(rep(metrics, times = nrow(df)), levels = metrics),
    Value = as.vector(t(as.matrix(df[metrics])))
  )
  overall <- list(); pairwise <- list(); letters <- list()
  for (metric in metrics) {
    sub <- long[long$Metric == metric, , drop = FALSE]
    key <- if (!is.null(metric_keys) && metric %in% names(metric_keys))
      metric_keys[[metric]] else metric
    stat <- if (stats_source == "qiime2") {
      qf <- resolve(file.path(qiime_stats_dir, paste0(key, "-group-significance.qzv")))
      q2_read_alpha_stats(qf, key, sub, group_col)
    } else alpha_r_stats(sub)
    ov <- stat$overall; pr <- stat$pairwise
    ov$Metric <- metric; ov$Source <- stats_source
    pr$Metric <- metric; pr$Source <- stats_source
    overall[[metric]] <- ov; pairwise[[metric]] <- pr
    letter <- setNames(rep("", length(group_levels)), group_levels)
    if (all(is.finite(pr$q_within_metric))) {
      pm <- matrix(1, length(group_levels), length(group_levels),
                     dimnames = list(group_levels, group_levels))
      for (j in seq_len(nrow(pr))) {
        pm[pr$Group1[j], pr$Group2[j]] <- pm[pr$Group2[j], pr$Group1[j]] <-
          pr$q_within_metric[j]
      }
      letter <- multcompView::multcompLetters(pm, threshold = significance_level)$Letters
    }
    letters[[metric]] <- data.frame(Metric = metric, Group = names(letter),
                                    Letter = unname(letter), stringsAsFactors = FALSE)
  }
  annotation <- do.call(rbind, letters)
  annotation$Metric <- factor(annotation$Metric, levels = metrics)
  annotation$Group <- factor(annotation$Group, levels = group_levels)
  # Exactly the original per-group maximum * 1.02 placement.
  annotation$y <- vapply(seq_len(nrow(annotation)), function(i) {
    max(long$Value[long$Metric == annotation$Metric[i] &
                     long$Group == annotation$Group[i]]) * 1.02
  }, numeric(1))
  p <- original_alpha_plot(long, annotation, font_family, group_colors, display_names, seed)
  out <- resolve(output_dir)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(out, "α多样性_KW_pairwiseKW_BH.svg"), plot = p,
                  width = 19, height = 5, dpi = 300, device = svglite::svglite)
  q2_write_csv(df, file.path(out, "alpha_plot_input.csv"))
  q2_write_csv(do.call(rbind, overall), file.path(out, "alpha_KW.csv"))
  q2_write_csv(do.call(rbind, pairwise), file.path(out, "alpha_pairwise_KW_BH.csv"))
  q2_write_csv(annotation, file.path(out, "alpha_letters_positions.csv"))
  utils::capture.output(utils::sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
  message("alpha 绘图完成：", out)
  invisible(list(plot = p, data = long, annotation = annotation,
                 overall = overall, pairwise = pairwise))
}
