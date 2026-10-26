# source("qiime2_io.R") before sourcing this file.
# Reads QIIME 2 coordinates, explained proportions and saved PERMANOVA P.

original_pcoa_plot <- function(points, proportion, p_value = NULL,
                               title = "Weighted UniFrac",
                               font_family = "Times New Roman",
                               group_colors = c(CON = "#D55E00", LDC = "#0072B2",
                                                HDC = "#CC79A7"),
                               ellipse = TRUE, point_size = 2, point_alpha = NULL,
                               ellipse_level = 0.95, ellipse_alpha = 0.15,
                               show_sample_labels = FALSE) {
  p <- ggplot2::ggplot(points, ggplot2::aes(x = PC1, y = PC2,
                                          color = Group, fill = Group))
  if (!is.null(p_value)) {
    p <- p + ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label = paste0("PERMANOVA: P = ", signif(p_value, 3), "\n"),
      hjust = 1.05, vjust = 1.2, size = 5, family = font_family
    )
  }
  p <- p + if (is.null(point_alpha)) ggplot2::geom_point(size = point_size) else
    ggplot2::geom_point(size = point_size, alpha = point_alpha)
  if (ellipse) {
    # stat_ellipse requires >=4 points and a nonsingular 2D covariance matrix.
    good <- vapply(split(points, points$Group, drop = TRUE), function(g) {
      nrow(g) >= 4L && qr(scale(as.matrix(g[c("PC1", "PC2")]),
                                center = TRUE, scale = FALSE))$rank == 2L
    }, logical(1))
    if (any(!good)) message("以下组不满足椭圆绘制条件，仅显示散点：",
                            paste(names(good)[!good], collapse = ", "))
    ed <- points[as.character(points$Group) %in% names(good)[good], , drop = FALSE]
    if (nrow(ed)) {
      p <- p + ggplot2::stat_ellipse(data = ed, geom = "polygon", type = "t",
                                    color = NA, level = ellipse_level,
                                    alpha = ellipse_alpha)
    }
  }
  if (show_sample_labels) {
    if (!"SampleID" %in% names(points)) stop("SampleID is required for sample labels.")
    p <- p + ggplot2::geom_text(ggplot2::aes(label = SampleID), vjust = -0.8,
                                size = 3, family = font_family, show.legend = FALSE)
  }
  p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40",
                        linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40",
                        linewidth = 0.5) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0("PCoA (", title, ")"),
      x = paste0("PCoA1 (", round(proportion[1] * 100, 1), "%)"),
      y = paste0("PCoA2 (", round(proportion[2] * 100, 1), "%)"), color = "Group"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title.y = ggplot2::element_text(size = 20, face = "bold", family = font_family),
      axis.title.x = ggplot2::element_text(size = 20, face = "bold", family = font_family),
      axis.text.x = ggplot2::element_text(size = 16, face = "bold", family = font_family),
      axis.text.y = ggplot2::element_text(size = 16, face = "bold", family = font_family),
      legend.text = ggplot2::element_text(family = font_family, size = 14),
      legend.title = ggplot2::element_text(family = font_family, size = 14),
      plot.title = ggplot2::element_text(size = 20, face = "bold", family = font_family,
                                        hjust = 0.5),
      text = ggplot2::element_text(family = font_family)
    )
}

run_beta_original_style <- function(
    project = "D:/qiime2/肥育猪试验-终版",
    core_dir = "core-metrics-results", metadata_file = "metadata.txt",
    qiime_stats_dir = "beta", output_dir = "R_original_style/beta",
    group_col = "Group", group_levels = c("CON", "LDC", "HDC"),
    group_colors = c(CON = "#D55E00", LDC = "#0072B2", HDC = "#CC79A7"),
    font_family = "Times New Roman", annotate_p = TRUE, check_distance = TRUE,
    ellipse = TRUE,
    metric_titles = c(weighted_unifrac = "Weighted UniFrac",
                      unweighted_unifrac = "Unweighted UniFrac",
                      bray_curtis = "Bray-Curtis", jaccard = "Jaccard")) {
  q2_packages(c("ggplot2", "patchwork", "svglite", "xml2", "systemfonts"))
  font_family <- q2_font(font_family)
  if (!dir.exists(project)) stop("工作目录不存在：", project)
  resolve <- function(p) q2_resolve(project, p)
  metrics <- names(metric_titles)
  if (!length(metrics) || anyDuplicated(metrics)) stop("metric_titles 需要唯一的指标名称")
  if (anyDuplicated(group_levels)) stop("group_levels 不能重复")
  metadata <- q2_read_metadata(resolve(metadata_file))
  if (!group_col %in% names(metadata)) stop("元数据没有分组列：", group_col)
  metadata$Group <- trimws(metadata[[group_col]])
  out <- resolve(output_dir)
  dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  plots <- list(); details <- list(); overall <- list(); pairwise <- list()
  baseline_ids <- NULL
  for (metric in metrics) {
    pf <- resolve(file.path(core_dir, paste0(metric, "_pcoa_results.qza")))
    ord <- q2_read_pcoa(pf)
    # Retain the stored PCoA row order and coordinate signs.
    ids <- rownames(ord$coords)
    if (!is.null(baseline_ids) && !setequal(ids, baseline_ids)) {
      stop("四个指标使用的样本不一致：", metric)
    }
    baseline_ids <- ids
    if (!all(ids %in% metadata$SampleID)) stop("PCoA 中的样本缺少元数据：", metric)
    md <- metadata[match(ids, metadata$SampleID), , drop = FALSE]
    if (anyNA(md$Group) || any(!md$Group %in% group_levels)) {
      stop("PCoA 保留样本的 Group 缺失或不在 group_levels 中：", metric)
    }
    md$Group <- droplevels(factor(md$Group, levels = group_levels))
    if (nlevels(md$Group) < 2L) stop("需要至少两个保留的组：", metric)
    if (!all(levels(md$Group) %in% names(group_colors))) stop("请为每个组设置颜色")
    dm <- NULL
    if (check_distance) {
      df <- resolve(file.path(core_dir, paste0(metric, "_distance_matrix.qza")))
      dm <- q2_read_distance(df)
      if (!setequal(ids, rownames(dm))) stop("PCoA/距离矩阵样本不一致：", metric)
      dm <- dm[ids, ids, drop = FALSE]
    }
    saved_p <- NULL
    for (method in c("permanova", "permdisp")) {
      sf <- resolve(file.path(qiime_stats_dir, paste0(metric, "-", method, ".qzv")))
      if (!file.exists(sf)) {
        if (method == "permanova" && annotate_p) {
          stop("缺少已保存的 QIIME PERMANOVA：", sf,
               "。先运行 qiime_group_tests.sh；仅画坐标时可设 annotate_p=FALSE。")
        }
        next
      }
      stat <- q2_read_qiime_stats(sf, method, dm, md)
      ov <- stat$overall
      ov$Metric <- metric; ov$Method <- toupper(method); ov$Source <- "QIIME2_saved"
      overall[[paste(metric, method)]] <- ov
      if (nrow(stat$pairwise)) {
        pr <- stat$pairwise
        pr$Metric <- metric; pr$Method <- toupper(method); pr$Source <- "QIIME2_saved"
        pairwise[[paste(metric, method)]] <- pr
      }
      if (method == "permanova" && annotate_p) saved_p <- ov$P
    }
    points <- data.frame(SampleID = ids, ord$coords, Group = md$Group,
                          check.names = FALSE)
    p <- original_pcoa_plot(points, ord$proportion, saved_p, metric_titles[[metric]],
                             font_family, group_colors, ellipse)
    plots[[metric]] <- p
    details[[metric]] <- list(points = points, proportion = ord$proportion, P = saved_p)
    ggplot2::ggsave(file.path(out, "figures", paste0("PCoA_", metric, ".svg")),
                    plot = p, width = 5, height = 4, device = svglite::svglite)
    q2_write_csv(points, file.path(out, "tables", paste0(metric, "_PCoA_coordinates.csv")))
    q2_write_csv(data.frame(Axis = paste0("PC", seq_along(ord$proportion)),
                             ProportionExplained = ord$proportion,
                             PercentExplained = ord$proportion * 100),
                  file.path(out, "tables", paste0(metric, "_PCoA_axes.csv")))
  }
  if (length(overall)) q2_write_csv(do.call(rbind, overall),
                                    file.path(out, "tables", "beta_QIIME_overall.csv"))
  if (length(pairwise)) q2_write_csv(do.call(rbind, pairwise),
                                     file.path(out, "tables", "beta_QIIME_pairwise.csv"))
  combined <- patchwork::wrap_plots(plots, ncol = 2)
  ggplot2::ggsave(file.path(out, "figures", "β多样性_四指标.svg"), plot = combined,
                  width = 10, height = 4 * ceiling(length(plots) / 2),
                  device = svglite::svglite)
  if (all(c("weighted_unifrac", "unweighted_unifrac") %in% names(plots))) {
    unifrac_pair <- patchwork::wrap_plots(
      plots[c("weighted_unifrac", "unweighted_unifrac")], ncol = 2)
    ggplot2::ggsave(file.path(out, "figures", "β多样性_UniFrac.svg"), plot = unifrac_pair,
                    width = 10, height = 4, device = svglite::svglite)
  }
  utils::capture.output(utils::sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
  message("PCoA 绘图完成：", file.path(out, "figures"))
  invisible(list(plots = plots, combined = combined, data = details,
                 overall = overall, pairwise = pairwise))
}
