# Presentation tables only: never replace count tables used in statistical models.
taxonomy_rank_label <- function(lineage, rank) {
  prefix <- c('k','p','c','o','f','g','s')[rank]
  vapply(strsplit(ifelse(is.na(lineage),'',lineage),';',fixed=TRUE), function(parts) {
    parts <- trimws(parts)
    hit <- grep(paste0('^(',prefix,'__|D_',rank-1,'__)'),parts,ignore.case=TRUE)
    if (length(hit)>1L) stop('Duplicate target rank in taxonomy lineage.')
    if (length(hit)!=1L) return('Others')
    label <- trimws(sub('^([A-Za-z]|D_[0-9]+)__','',parts[hit],ignore.case=TRUE))
    if (!nzchar(label) || tolower(label) %in% c('na','n/a','nan','none','null','missing','other','others','not assigned','not_assigned') ||
        grepl('incertae[[:space:]_-]*sedis|no[[:space:]_-]*rank|uncultured|unclassified|unassigned|unknown|unidentified|unresolved|ambiguous[[:space:]_-]*taxa',label,ignore.case=TRUE)) 'Others' else label
  },character(1),USE.NAMES=FALSE)
}

taxonomy_clean_counts <- function(counts, lineage, rank) {
  if (length(lineage)!=nrow(counts) || any(!is.finite(counts)) || any(counts<0) || any(counts!=floor(counts))) stop('Invalid taxonomy counts.')
  labels <- taxonomy_rank_label(lineage,rank)
  clean <- rowsum(counts,labels,reorder=FALSE)
  totals <- colSums(clean)
  relative <- sweep(clean,2,ifelse(totals>0,totals,1),'/')
  list(counts=clean,relative=relative,labels=labels)
}

taxonomy_top_table <- function(relative, rank, options) {
  limit <- switch(options$top_mode,by_rank=if (rank %in% 4:6) 20L else Inf,
                  custom=options$top_n,all=Inf)
  # Each sample has equal weight; library depth does not decide the ranking.
  score <- rowMeans(relative)
  names_ranked <- setdiff(rownames(relative)[order(-score,rownames(relative))],'Others')
  keep <- if (is.infinite(limit)) names_ranked else head(names_ranked,limit)
  rest <- setdiff(rownames(relative),keep)
  # A zero Others row is retained to make the meaning and color assignment stable.
  plotted <- rbind(relative[keep,,drop=FALSE],Others=colSums(relative[rest,,drop=FALSE]))
  list(relative=plotted,limit=if (is.infinite(limit)) NULL else limit,kept=keep,
       other_taxa=rest,ranking='mean per-sample relative abundance; ties by label')
}

taxonomy_fill_colors <- function(taxa, palette, registry) {
  definition <- registry[[palette]]
  if (is.null(definition)) stop('Unknown taxonomy palette.')
  known <- setdiff(taxa,'Others'); n <- length(known)
  if (identical(definition$mode,'hcl')) colors <- if(n) grDevices::hcl.colors(n,palette) else character()
  else if (identical(definition$mode,'hue')) colors <- if(n) scales::hue_pal()(n) else character()
  else {
    base <- unlist(definition$colors,use.names=FALSE)
    colors <- if (n<=length(base)) head(base,n) else c(base,grDevices::hcl.colors(n-length(base),'Dark 3'))
  }
  c(setNames(colors,known),Others='#857D83')[taxa]
}
