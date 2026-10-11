# QIIME 2 2025.7 standard QZA/QZV readers. No qiime2R installation needed.
# Read saved data only; no rarefaction, distance, PCoA, or permutation calculation.

q2_read_zip_text <- function(path, relative_path) {
  if (!file.exists(path)) stop("File not found: ", path)
  z <- utils::unzip(path, list = TRUE)$Name
  root_meta <- z[grepl("^[^/]+/metadata\\.yaml$", z)]
  if (length(root_meta) != 1L) stop("Not a standard QIIME 2 archive: ", path)
  root <- sub("/metadata\\.yaml$", "", root_meta)
  target <- paste0(root, "/", relative_path)
  if (!target %in% z) stop("Missing ", relative_path, " in ", path)
  td <- tempfile("q2-beta-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::unzip(path, files = target, exdir = td)
  readLines(file.path(td, target), warn = FALSE, encoding = "UTF-8")
}

q2_uuid <- function(path) {
  x <- q2_read_zip_text(path, "metadata.yaml")
  u <- trimws(sub("^uuid:\\s*", "", x[grepl("^uuid:", x)]))
  if (length(u) != 1L || !nzchar(u)) stop("Cannot read QZA UUID: ", path)
  u
}

q2_require_type <- function(path, expected) {
  x <- q2_read_zip_text(path, "metadata.yaml")
  ty <- trimws(sub("^type:\\s*", "", x[grepl("^type:", x)]))
  ty <- gsub("^['\"]|['\"]$", "", ty)
  if (length(ty) != 1L || ty != expected) {
    stop("Expected ", expected, " in ", path, "; found ", paste(ty, collapse = ","))
  }
}

q2_read_metadata <- function(path) {
  if (!file.exists(path)) stop("Metadata file not found: ", path)
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x <- x[nzchar(trimws(x))]
  if (!length(x)) stop("Empty metadata file")
  x[1] <- sub("^\ufeff", "", x[1])
  # Keep a #SampleID header, remove subsequent QIIME 2 directives/comments.
  x <- c(x[1], x[-1][!grepl("^\\s*#", x[-1])])
  m <- utils::read.delim(text = paste(x, collapse = "\n"),
                        check.names = FALSE, colClasses = "character",
                        quote = "", comment.char = "", na.strings = c("", "NA"))
  if (ncol(m) < 2L || anyDuplicated(names(m))) stop("Invalid metadata columns")
  names(m)[1] <- "SampleID"
  m$SampleID <- trimws(m$SampleID)
  if (anyNA(m$SampleID) || any(!nzchar(m$SampleID)) || anyDuplicated(m$SampleID)) {
    stop("Metadata SampleID must be nonempty and unique")
  }
  m
}

q2_read_distance <- function(path) {
  q2_require_type(path, "DistanceMatrix")
  x <- q2_read_zip_text(path, "data/distance-matrix.tsv")
  d <- utils::read.delim(text = paste(x, collapse = "\n"),
                        check.names = FALSE, colClasses = "character",
                        quote = "", comment.char = "")
  ids <- d[[1]]
  if (anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids))) stop("Invalid distance IDs")
  values <- as.matrix(d[-1])
  storage.mode(values) <- "double"
  rownames(values) <- ids
  if (nrow(values) != ncol(values) || !setequal(ids, colnames(values))) {
    stop("Distance row and column IDs do not match: ", path)
  }
  values <- values[ids, ids, drop = FALSE]
  if (any(!is.finite(values)) || any(values < 0) ||
      max(abs(diag(values))) > 1e-10 || max(abs(values - t(values))) > 1e-10) {
    stop("Distance must be finite, nonnegative, symmetric, with zero diagonal")
  }
  values
}

q2_read_pcoa <- function(path) {
  q2_require_type(path, "PCoAResults")
  x <- q2_read_zip_text(path, "data/ordination.txt")
  # Standard scikit-bio ordination format: use saved Site and Proportion explained.
  section <- function(name) {
    h <- which(startsWith(x, paste0(name, "\t")))
    if (length(h) != 1L) stop("Cannot locate ordination section: ", name)
    dims <- as.integer(strsplit(x[h], "\t", fixed = TRUE)[[1]][-1])
    nr <- if (length(dims) == 1L) as.integer(dims[1] > 0) else dims[1]
    if (anyNA(dims) || any(dims < 0)) stop("Invalid ordination dimensions")
    list(dims = dims, lines = if (nr > 0) x[h + seq_len(nr)] else character())
  }
  st <- section("Site")
  pe <- section("Proportion explained")
  ev <- section("Eigvals")
  if (length(st$dims) != 2L || st$dims[2] < 2L) stop("PCoA needs two stored axes")
  cells <- strsplit(st$lines, "\t", fixed = TRUE)
  if (any(lengths(cells) != st$dims[2] + 1L)) stop("Malformed PCoA Site section")
  ids <- vapply(cells, `[`, character(1), 1L)
  if (anyDuplicated(ids) || any(!nzchar(ids))) stop("Invalid PCoA sample IDs")
  coords <- t(vapply(cells, function(a) as.numeric(a[-1]), numeric(st$dims[2])))
  rownames(coords) <- ids
  colnames(coords) <- paste0("PC", seq_len(ncol(coords)))
  prop <- as.numeric(strsplit(pe$lines, "\t", fixed = TRUE)[[1]])
  eig <- as.numeric(strsplit(ev$lines, "\t", fixed = TRUE)[[1]])
  if (length(prop) != ncol(coords) || length(eig) != ncol(coords) ||
      any(!is.finite(coords)) || any(!is.finite(prop)) || any(!is.finite(eig))) {
    stop("Invalid PCoA coordinates or explained proportions")
  }
  list(coords = coords, proportion = prop, eigenvalues = eig)
}

q2_pair_key <- function(a, b) {
  mapply(function(x, y) paste(sort(c(x, y)), collapse = "\034"), a, b,
         USE.NAMES = FALSE)
}

q2_read_qiime_stats <- function(path, method, dm, meta, permutations = NULL) {
  q2_require_type(path, "Visualization")
  html <- q2_read_zip_text(path, "data/index.html")
  doc <- xml2::read_html(paste(html, collapse = "\n"))
  rows <- xml2::xml_find_all(doc, "//table//tr")
  kv <- list()
  for (tr in rows) {
    v <- trimws(xml2::xml_text(xml2::xml_find_all(tr, "./th|./td")))
    if (length(v) == 2L && nzchar(v[1])) kv[[v[1]]] <- v[2]
  }
  required <- c("method name", "sample size", "number of groups", "test statistic",
                "p-value", "number of permutations")
  if (!all(required %in% names(kv))) stop("Cannot parse QIIME 2 overview: ", path)
  if (tolower(kv[["method name"]]) != method) stop("Wrong test method in ", path)
  numeric_value <- function(key) {
    v <- as.numeric(kv[[key]])
    if (length(v) != 1L || !is.finite(v)) stop("Invalid QIIME 2 statistic: ", key)
    v
  }
  result <- data.frame(N = numeric_value("sample size"),
                       Groups = numeric_value("number of groups"),
                       F = numeric_value("test statistic"),
                       P = numeric_value("p-value"),
                       Permutations = numeric_value("number of permutations"))
  if (result$N != nrow(meta) || result$Groups != nlevels(meta$Group) ||
      (!is.null(permutations) && result$Permutations != permutations)) {
    stop("QZV sample/group/permutation counts disagree with this R analysis: ", path)
  }
  # raw_data.tsv records all unique sample pairs (some occur twice).
  # Validate the saved group assignments and distances, preventing stale QZV reuse.
  raw <- q2_read_zip_text(path, "data/raw_data.tsv")
  raw <- utils::read.delim(text = paste(raw, collapse = "\n"), check.names = FALSE,
                          colClasses = "character", quote = "", comment.char = "")
  req <- c("SubjectID1", "SubjectID2", "Group1", "Group2", "Distance")
  if (!all(req %in% names(raw))) stop("Unsupported QZV raw_data.tsv")
  gm <- unique(rbind(data.frame(ID = raw$SubjectID1, Group = raw$Group1),
                    data.frame(ID = raw$SubjectID2, Group = raw$Group2)))
  if (anyDuplicated(gm$ID) || !setequal(gm$ID, meta$SampleID)) {
    stop("QZV sample IDs do not match distance matrix")
  }
  if (any(gm$Group[match(meta$SampleID, gm$ID)] != as.character(meta$Group))) {
    stop("QZV saved groups differ from current metadata")
  }
  a <- match(raw$SubjectID1, meta$SampleID)
  b <- match(raw$SubjectID2, meta$SampleID)
  saved <- as.numeric(raw$Distance)
  if (anyNA(a) || anyNA(b) || any(!is.finite(saved)) ||
      (!is.null(dm) && max(abs(dm[cbind(a, b)] - saved)) > 1e-9) ||
      length(unique(q2_pair_key(raw$SubjectID1, raw$SubjectID2))) != choose(nrow(meta), 2)) {
    stop("QZV contains different or incomplete distances: ", path)
  }
  # CSV preserves the QIIME 2 pairwise P and q values: do not adjust q again.
  pair <- data.frame()
  z <- utils::unzip(path, list = TRUE)$Name
  if (any(endsWith(z, paste0("/data/", method, "-pairwise.csv")))) {
    px <- q2_read_zip_text(path, paste0("data/", method, "-pairwise.csv"))
    p <- utils::read.csv(text = paste(px, collapse = "\n"), check.names = FALSE,
                         colClasses = "character")
    need <- c("Group 1", "Group 2", "Sample size", "Permutations", "p-value", "q-value")
    if (!all(need %in% names(p))) stop("Unsupported QIIME 2 pairwise CSV")
    stat_col <- setdiff(names(p), need)
    if (length(stat_col) != 1L) stop("Cannot identify pairwise test statistic")
    pair <- data.frame(Group1 = p[["Group 1"]], Group2 = p[["Group 2"]],
                       N = as.numeric(p[["Sample size"]]),
                       F = as.numeric(p[[stat_col]]),
                       P = as.numeric(p[["p-value"]]),
                       Permutations = as.numeric(p[["Permutations"]]),
                       q_within_metric = as.numeric(p[["q-value"]]))
    expected <- utils::combn(levels(meta$Group), 2)
    keys <- q2_pair_key(pair$Group1, pair$Group2)
    target <- q2_pair_key(expected[1, ], expected[2, ])
    if (anyDuplicated(keys) || !setequal(keys, target) ||
        any(!is.finite(as.matrix(pair[c("N", "F", "P", "Permutations", "q_within_metric")])))) {
      stop("Invalid or incomplete QIIME 2 pairwise comparisons")
    }
    pair <- pair[match(target, keys), , drop = FALSE]
    if (any(pair$Permutations != result$Permutations)) stop("Pairwise permutation counts disagree")
  }
  list(overall = result, pairwise = pair)
}

q2_resolve <- function(project, path) {
  if (grepl("^(/|~|[A-Za-z]:[/\\\\])", path)) path.expand(path) else
    file.path(project, path)
}

q2_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("请先安装 R 包：", paste(missing, collapse = ", "))
}

q2_font <- function(family) {
  q2_packages("systemfonts")
  if (!family %in% c("serif", "sans", "mono") &&
      !family %in% systemfonts::system_fonts()$family) {
    stop("找不到字体：", family,
         "。请安装该字体并重启 R，或手动修改 font_family；代码不会自动替换字体。")
  }
  family
}

q2_write_csv <- function(x, path) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw("\ufeff"), con)
  utils::write.table(x, con, sep = ",", row.names = FALSE, col.names = TRUE,
                     quote = TRUE, na = "", qmethod = "double")
}

q2_read_alpha <- function(path) {
  q2_require_type(path, "SampleData[AlphaDiversity]")
  x <- q2_read_zip_text(path, "data/alpha-diversity.tsv")
  d <- utils::read.delim(text = paste(x, collapse = "\n"), check.names = FALSE,
                         colClasses = "character", quote = "", comment.char = "")
  if (ncol(d) != 2L || anyNA(d[[1]]) || anyDuplicated(d[[1]])) {
    stop("无效的 alpha vector：", path)
  }
  value <- as.numeric(d[[2]])
  if (any(!is.finite(value))) stop("alpha vector 含非数值：", path)
  data.frame(Sample = d[[1]], Value = value, stringsAsFactors = FALSE)
}

q2_read_alpha_stats <- function(path, metric_key, data, group_col = "Group") {
  q2_require_type(path, "Visualization")
  q2_packages("jsonlite")
  escaped <- utils::URLencode(group_col, reserved = TRUE)
  jp <- paste(q2_read_zip_text(path, paste0("data/column-", escaped, ".jsonp")),
              collapse = "\n")
  # Read the numeric JSON object only; never execute the JavaScript callback.
  objects <- regmatches(jp, gregexpr("\\{[^{}]*\\}", jp, perl = TRUE))[[1]]
  objects <- objects[grepl('"H"\\s*:', objects) & grepl('"p"\\s*:', objects)]
  if (length(objects) != 1L) stop("无法读取 QIIME alpha 全局 KW 结果：", path)
  global <- jsonlite::fromJSON(objects)
  if (length(global$H) != 1L || length(global$p) != 1L ||
      !is.finite(global$H) || !is.finite(global$p)) {
    stop("QIIME alpha KW 结果无效：", path)
  }
  # metadata.tsv contains the sample-level alpha values used in this QZV.
  saved_lines <- q2_read_zip_text(path, "data/metadata.tsv")
  saved_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(saved_file), add = TRUE)
  writeLines(saved_lines, saved_file, useBytes = TRUE)
  saved <- q2_read_metadata(saved_file)
  if (!setequal(saved$SampleID, data$Sample) || !group_col %in% names(saved)) {
    stop("alpha QZV 与绘图数据的样本/分组列不一致：", path)
  }
  saved <- saved[match(data$Sample, saved$SampleID), , drop = FALSE]
  if (anyNA(saved[[group_col]]) ||
      any(saved[[group_col]] != as.character(data$Group))) {
    stop("alpha QZV 保存的分组与当前输入不一致：", path)
  }
  candidates <- if (metric_key %in% names(saved)) metric_key else
    setdiff(names(saved), c("SampleID", group_col))
  matches <- vapply(candidates, function(nm) {
    v <- suppressWarnings(as.numeric(saved[[nm]]))
    all(is.finite(v)) &&
      max(abs(v - data$Value) / pmax(1, abs(data$Value))) < 1e-9
  }, logical(1))
  if (!any(matches)) {
    stop("alpha QZV 的指标值与当前输入不一致。请保留原始数值精度：", path)
  }
  px <- q2_read_zip_text(path,
                        paste0("data/kruskal-wallis-pairwise-", escaped, ".csv"))
  p <- utils::read.csv(text = paste(px, collapse = "\n"), check.names = FALSE,
                       colClasses = "character")
  required <- c("Group 1", "Group 2", "H", "p-value", "q-value")
  if (!all(required %in% names(p))) stop("无法读取 QIIME alpha 两两比较 CSV")
  pair <- data.frame(
    Group1 = sub(" \\(n=[0-9]+\\)$", "", p[["Group 1"]]),
    Group2 = sub(" \\(n=[0-9]+\\)$", "", p[["Group 2"]]),
    H = as.numeric(p[["H"]]), P = as.numeric(p[["p-value"]]),
    q_within_metric = as.numeric(p[["q-value"]]), stringsAsFactors = FALSE
  )
  pairs <- utils::combn(levels(droplevels(data$Group)), 2)
  target <- q2_pair_key(pairs[1, ], pairs[2, ])
  key <- q2_pair_key(pair$Group1, pair$Group2)
  if (anyDuplicated(key) || !all(key %in% target)) stop("alpha QZV 比较的组不一致")
  # QIIME can omit a comparison if all values in those two groups are identical.
  pos <- match(target, key)
  result <- data.frame(Group1 = pairs[1, ], Group2 = pairs[2, ],
                        H = pair$H[pos], P = pair$P[pos],
                        q_within_metric = pair$q_within_metric[pos])
  result$Status <- ifelse(is.na(pos), "not_in_qiime_output", "ok")
  list(overall = data.frame(KW_H = global$H, KW_p = global$p, N = nrow(data)),
       pairwise = result)
}

