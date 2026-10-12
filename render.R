# Render saved QIIME outputs only: no rarefaction, ordination, or tests rerun in R.
args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop('Usage: Rscript render.R RUN_DIR FONT [PREPROCESS_DIR] [PLOTS_DIR]')
run <- normalizePath(args[1], mustWork=TRUE)
font <- if (length(args)>1) args[2] else 'DejaVu Serif'
parent <- if (length(args)>2) normalizePath(args[3], mustWork=TRUE) else dirname(run)
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
assets <- if (length(args)>3 && nzchar(args[4])) args[4] else file.path(dirname(normalizePath(script)),'plots')
needed <- c('ggplot2','svglite','systemfonts','patchwork','xml2','jsonlite','MASS')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop('Missing R packages: ',paste(missing,collapse=', '))
if (!font %in% c('serif','sans','mono') && !font %in% systemfonts::system_fonts()$family)
  stop('Font not installed: ',font,'. Select an installed font in the workstation.')
if (!capabilities('cairo')) stop('R needs Cairo graphics support for reliable PDF/PNG output.')
result_path <- function(base,...) {
  logical <- file.path(...); layout_file <- file.path(base,'layout.json')
  if (file.exists(layout_file)) {
    mapping <- jsonlite::read_json(layout_file)$paths
    keys <- names(mapping)[order(nchar(names(mapping)),decreasing=TRUE)]
    for (key in keys) if (logical==key || startsWith(logical,paste0(key,'/'))) {
      logical <- paste0(mapping[[key]],substring(logical,nchar(key)+1)); break
    }
  }
  file.path(base,logical)
}
out <- result_path(run,'R_original_style'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
notes <- character()
note <- function(x) { notes <<- c(notes,x); message(x) }
csv <- function(d,f) {
  dir.create(dirname(f),recursive=TRUE,showWarnings=FALSE)
  con <- file(f,'wb'); on.exit(close(con)); writeBin(charToRaw('\ufeff'),con)
  utils::write.table(d,con,sep=',',row.names=FALSE,quote=TRUE,na='')
}
tab <- function(f,metadata=FALSE) {
  x <- readLines(f,encoding='UTF-8',warn=FALSE)
  if (!length(x)) stop('Empty result: ',f)
  x[1] <- sub('^\ufeff','',x[1])
  if (metadata) x <- c(x[1],x[-1][!grepl('^\\s*#',x[-1])])
  utils::read.delim(text=paste(x,collapse='\n'),check.names=FALSE,colClasses='character',quote='',comment.char='',na.strings='')
}
md <- tab(result_path(run,'metadata.tsv'),TRUE); names(md)[1] <- 'SampleID'
if (!'Group' %in% names(md) || anyNA(md[c('SampleID','Group')]) || anyDuplicated(md$SampleID) || any(!nzchar(md$SampleID)) || any(!nzchar(md$Group))) stop('Invalid SampleID/Group metadata.')
groups <- unique(md$Group); colors <- c(CON='#D55E00',LDC='#0072B2',HDC='#CC79A7')
extra <- setdiff(groups,names(colors))
if (length(extra)) colors <- c(colors,setNames(grDevices::hcl.colors(length(extra),'Dark 3'),extra))
original <- all(file.exists(file.path(assets,c('qiime2_io.R','plot_alpha_original_style.R','plot_beta_original_style.R'))))
if (original) for (f in c('qiime2_io.R','plot_alpha_original_style.R','plot_beta_original_style.R')) source(file.path(assets,f),encoding='UTF-8')
if (!original) note('Original assets unavailable; embedded plots omit significance annotations.')
match_md <- function(ids) {
  if (anyDuplicated(ids) || anyNA(ids) || any(!nzchar(ids)) || !all(ids %in% md$SampleID)) stop('Saved result sample IDs invalid or absent from metadata.')
  d <- md[match(ids,md$SampleID),,drop=FALSE]; d$Group <- droplevels(factor(d$Group,levels=groups)); d
}
save_plot <- function(p,name,w,h) {
  f <- file.path(out,name); dir.create(dirname(f),recursive=TRUE,showWarnings=FALSE)
  ggplot2::ggsave(paste0(f,'.svg'),p,width=w,height=h,device=svglite::svglite,bg='white')
  ggplot2::ggsave(paste0(f,'.pdf'),p,width=w,height=h,device=grDevices::cairo_pdf,family=font,bg='white')
  ggplot2::ggsave(paste0(f,'.png'),p,width=w,height=h,dpi=300,type='cairo',bg='white')
}
theme_base <- function() ggplot2::theme_minimal(base_size=14,base_family=font)+
  ggplot2::theme(panel.grid=ggplot2::element_blank(),panel.border=ggplot2::element_rect(color='black',fill=NA),axis.text=ggplot2::element_text(face='bold'),axis.title=ggplot2::element_text(face='bold'),strip.background=ggplot2::element_rect(fill='gray80',color='black'))
alpha_titles <- c(observed_features='Observed features',shannon='Shannon diversity',evenness='Pielou evenness',faith_pd='Faith phylogenetic diversity',chao1='Chao1 richness',simpson='Simpson diversity',goods_coverage="Good's coverage")
alpha <- list(); annotations <- list(); overall <- list(); pairwise <- list(); alpha_ids <- NULL
for (metric in names(alpha_titles)) {
  f <- result_path(run,'exports','alpha',metric,'alpha-diversity.tsv')
  if (!file.exists(f)) stop('Required alpha export missing: ',metric,' (',f,'). Export the saved QIIME result before rendering.')
  v <- tab(f); if (ncol(v)!=2) stop('Invalid alpha vector: ',f)
  d <- match_md(v[[1]])
  if (!is.null(alpha_ids) && !setequal(alpha_ids,d$SampleID)) stop('Alpha metrics have different sample sets.')
  alpha_ids <- d$SampleID
  d <- data.frame(Sample=d$SampleID,Group=d$Group,Metric=metric,Value=as.numeric(v[[2]]))
  if (any(!is.finite(d$Value)) || any(d$Value<0)) stop('Invalid alpha values: ',metric)
  alpha[[metric]] <- d; an <- data.frame(Metric=metric,Group=levels(d$Group),Letter='')
  sf <- result_path(run,'alpha',paste0(metric,'-group-significance.qzv'))
  if (!file.exists(sf)) note(paste(metric,': no saved alpha group test; significance letters omitted.'))
  if (original && file.exists(sf) && nlevels(d$Group)>=2) {
    stat <- tryCatch(q2_read_alpha_stats(sf,metric,d),error=function(e) { note(paste(metric,'saved alpha tests not annotated:',e$message)); NULL })
    if (!is.null(stat)) {
      stat$overall$Metric <- metric; stat$pairwise$Metric <- metric
      overall[[metric]] <- stat$overall; pairwise[[metric]] <- stat$pairwise
      pr <- stat$pairwise
      valid <- all(is.finite(pr$q_within_metric)) && all(pr$q_within_metric>=0 & pr$q_within_metric<=1)
      if (valid && requireNamespace('multcompView',quietly=TRUE)) {
        pm <- matrix(1,nlevels(d$Group),nlevels(d$Group),dimnames=list(levels(d$Group),levels(d$Group)))
        for (j in seq_len(nrow(pr))) pm[pr$Group1[j],pr$Group2[j]] <- pm[pr$Group2[j],pr$Group1[j]] <- pr$q_within_metric[j]
        lt <- multcompView::multcompLetters(pm,threshold=.05)$Letters; an$Letter <- unname(lt[as.character(an$Group)])
      } else note(paste(metric,'letters omitted: incomplete saved tests or multcompView missing.'))
    }
  }
  an$y <- vapply(as.character(an$Group),function(g) max(d$Value[as.character(d$Group)==g])*1.02,numeric(1)); annotations[[metric]] <- an
}
if (length(alpha)) {
  a <- do.call(rbind,alpha); a$Metric <- factor(a$Metric,levels=names(alpha))
  an <- do.call(rbind,annotations); an$Metric <- factor(an$Metric,levels=names(alpha)); an$Group <- factor(an$Group,levels=groups)
  unlink(file.path(out,'alpha',paste0('alpha_diversity.',c('svg','pdf','png'))))
  for (metric in names(alpha)) {
    d <- droplevels(a[a$Metric==metric,,drop=FALSE]); annotation <- droplevels(an[an$Metric==metric,,drop=FALSE])
    if (original) p <- original_alpha_plot(d,annotation,font,colors,alpha_titles)+ggplot2::labs(y=alpha_titles[[metric]]) else p <-
      ggplot2::ggplot(d,ggplot2::aes(Group,Value,color=Group))+
      ggplot2::geom_boxplot(fill=NA,outlier.shape=NA,width=.6)+ggplot2::geom_jitter(position=ggplot2::position_jitter(width=.2,seed=123),size=1.5)+
      ggplot2::scale_color_manual(values=colors)+theme_base()+ggplot2::theme(legend.position='none')+ggplot2::labs(title=alpha_titles[[metric]],x='',y=alpha_titles[[metric]])
    save_plot(p,paste0('alpha/',metric),if (original) 4.75 else 6,5)
    csv(d,file.path(out,'alpha',paste0(metric,'_plot_input.csv')))
    csv(annotation,file.path(out,'alpha',paste0(metric,'_letters_positions.csv')))
  }
  csv(a,file.path(out,'alpha','alpha_plot_input.csv')); csv(an,file.path(out,'alpha','alpha_letters_positions.csv'))
  if (length(overall)) csv(do.call(rbind,overall),file.path(out,'alpha','alpha_KW_saved.csv'))
  if (length(pairwise)) csv(do.call(rbind,pairwise),file.path(out,'alpha','alpha_pairwise_KW_BH_saved.csv'))
}
read_ord <- function(f) {
  x <- readLines(f,warn=FALSE)
  section <- function(key) {
    h <- which(startsWith(x,paste0(key,'\t'))); if (length(h)!=1) stop('Missing/duplicate PCoA section: ',key)
    dims <- as.integer(strsplit(x[h],'\t',fixed=TRUE)[[1]][-1]); if (anyNA(dims) || any(dims<0)) stop('Invalid PCoA dimensions.')
    n <- if (length(dims)==1) as.integer(dims[1]>0) else dims[1]
    list(dims=dims,lines=if (n) x[h+seq_len(n)] else character())
  }
  s <- section('Site'); pe <- section('Proportion explained')
  if (length(s$dims)!=2 || s$dims[2]<2 || s$dims[1]<1) stop('PCoA needs two stored axes.')
  cells <- strsplit(s$lines,'\t',fixed=TRUE); if (any(lengths(cells)!=s$dims[2]+1)) stop('Malformed PCoA coordinates.')
  ids <- vapply(cells,`[`,character(1),1); coords <- t(vapply(cells,function(a) as.numeric(a[-1]),numeric(s$dims[2])))
  colnames(coords) <- paste0('PC',seq_len(ncol(coords))); rownames(coords) <- ids
  prop <- as.numeric(strsplit(pe$lines,'\t',fixed=TRUE)[[1]])
  if (length(prop)!=ncol(coords) || any(!is.finite(coords)) || any(!is.finite(prop)) || any(prop<0) || any(prop>1)) stop('Invalid stored PCoA coordinates/proportions.')
  list(ids=ids,coords=coords,proportion=prop)
}
read_dm <- function(f,ids) {
  d <- tab(f); m <- as.matrix(d[-1]); storage.mode(m) <- 'double'; rownames(m) <- d[[1]]
  if (anyDuplicated(d[[1]]) || nrow(m)!=ncol(m) || !setequal(ids,rownames(m)) || !setequal(ids,colnames(m))) stop('Distance sample IDs disagree.')
  m <- m[ids,ids,drop=FALSE]
  if (any(!is.finite(m)) || any(m<0) || max(abs(diag(m)))>1e-10 || max(abs(m-t(m)))>1e-10) stop('Invalid saved distance matrix.')
  m
}
titles <- c(weighted_unifrac='Weighted UniFrac',unweighted_unifrac='Unweighted UniFrac',bray_curtis='Bray-Curtis',jaccard='Jaccard')
plots <- list(); beta_stats <- list(); beta_ids <- NULL
for (metric in names(titles)) {
  f <- result_path(run,'exports','pcoa',metric,'ordination.txt')
  if (!file.exists(f)) { note(paste('PCoA export missing:',metric)); next }
  ord <- read_ord(f); d <- match_md(ord$ids)
  if (!is.null(beta_ids) && !setequal(beta_ids,ord$ids)) stop('PCoA metrics have different sample sets.')
  beta_ids <- ord$ids; dm <- read_dm(result_path(run,'exports','distance',metric,'distance-matrix.tsv'),ord$ids); saved_p <- NULL
  for (method in c('permanova','permdisp')) {
    sf <- result_path(run,'beta',paste0(metric,'-',method,'.qzv'))
    if (!file.exists(sf)) note(paste(metric,method,': no saved test; annotation omitted.'))
    if (!original || !file.exists(sf) || nlevels(d$Group)<2) next
    stat <- tryCatch(q2_read_qiime_stats(sf,method,dm,d),error=function(e) { note(paste(metric,method,'not annotated:',e$message)); NULL })
    if (!is.null(stat)) {
      ov <- stat$overall; ov$Metric <- metric; ov$Method <- method; beta_stats[[paste(metric,method)]] <- ov
      if (nrow(stat$pairwise)) csv(stat$pairwise,file.path(out,'beta',paste0(metric,'_',method,'_pairwise_saved.csv')))
      if (method=='permanova' && ov$P>=0 && ov$P<=1) saved_p <- ov$P
    }
  }
  pts <- data.frame(SampleID=ord$ids,ord$coords,Group=d$Group,check.names=FALSE)
  if (original) p <- original_pcoa_plot(pts,ord$proportion,saved_p,titles[[metric]],font,colors) else p <-
    ggplot2::ggplot(pts,ggplot2::aes(PC1,PC2,color=Group))+ggplot2::geom_point(size=2)+
    ggplot2::geom_hline(yintercept=0,linetype='dashed',color='grey40')+ggplot2::geom_vline(xintercept=0,linetype='dashed',color='grey40')+
    ggplot2::scale_color_manual(values=colors)+theme_base()+ggplot2::labs(title=paste0('PCoA (',titles[[metric]],')'),
      x=sprintf('PCoA1 (%.1f%%)',ord$proportion[1]*100),y=sprintf('PCoA2 (%.1f%%)',ord$proportion[2]*100))
  plots[[metric]] <- p; save_plot(p,paste0('beta/PCoA_',metric),6,5)
  csv(pts,file.path(out,'beta',paste0(metric,'_PCoA_coordinates.csv')))
  csv(data.frame(Axis=paste0('PC',seq_along(ord$proportion)),ProportionExplained=ord$proportion),file.path(out,'beta',paste0(metric,'_PCoA_axes.csv')))
}
if (length(plots)) save_plot(patchwork::wrap_plots(plots,ncol=2),'beta/PCoA_four_metrics',12,5*ceiling(length(plots)/2))
if (length(beta_stats)) csv(do.call(rbind,beta_stats),file.path(out,'beta','beta_QIIME_overall_saved.csv'))
# Composition uses filtered UNRAREFIED counts and each sample's relative abundance.
tf <- result_path(parent,'exported-table','feature-table.tsv'); tx <- result_path(parent,'exported-taxonomy','taxonomy.tsv')
if (file.exists(tf) && file.exists(tx)) {
  lines <- readLines(tf,warn=FALSE,encoding='UTF-8'); h <- which(grepl('^(#OTU ID|#?Feature ID|#?feature-id)\t',lines))
  if (length(h)!=1) stop('Cannot find BIOM TSV header.')
  tbl <- utils::read.delim(text=paste(lines[h:length(lines)],collapse='\n'),check.names=FALSE,colClasses='character',quote='',comment.char='')
  tax <- tab(tx); names(tax)[1] <- 'FeatureID'
  if (!'Taxon' %in% names(tax) || anyDuplicated(tax$FeatureID) || anyDuplicated(tbl[[1]])) stop('Invalid taxonomy/table feature IDs.')
  counts <- as.matrix(tbl[-1]); storage.mode(counts) <- 'double'
  if (any(!is.finite(counts)) || any(counts<0)) stop('Invalid feature counts.')
  pmd <- tab(result_path(parent,'metadata.tsv'),TRUE); names(pmd)[1] <- 'SampleID'
  if (!'Group' %in% names(pmd) || anyDuplicated(pmd$SampleID) || !all(colnames(counts) %in% pmd$SampleID)) stop('Taxonomy metadata mismatch.')
  pmd <- pmd[match(colnames(counts),pmd$SampleID),,drop=FALSE]
  if (anyNA(pmd$Group) || any(!nzchar(pmd$Group))) stop('Taxonomy groups missing.')
  valid <- colSums(counts)>0
  if (any(!valid)) note(paste('Zero-read samples omitted from composition:',paste(colnames(counts)[!valid],collapse=', ')))
  if (!any(valid)) stop('No nonzero samples for taxonomy composition.')
  counts <- counts[,valid,drop=FALSE]; pmd <- pmd[valid,,drop=FALSE]
  lineage <- tax$Taxon[match(tbl[[1]],tax$FeatureID)]; split_tax <- strsplit(ifelse(is.na(lineage),'',lineage),';',fixed=TRUE)
  for (rank_name in c('phylum','genus')) {
    rank <- if (rank_name=='phylum') 2 else 6
    labels <- vapply(split_tax,function(v) {
      hit <- grep(if (rank==2) '^\\s*(p__|D_1__)' else '^\\s*(g__|D_5__)',v)
      z <- if (length(hit)==1) trimws(v[hit]) else if (length(v)>=rank) trimws(v[rank]) else ''
      z <- sub('^([A-Za-z]|D_[0-9]+)__[[:space:]]*','',z)
      if (!nzchar(z) || tolower(z) %in% c('unassigned','uncultured','unknown')) 'Unclassified' else z
    },character(1))
    summed <- rowsum(counts,labels,reorder=FALSE); rel <- sweep(summed,2,colSums(summed),'/')
    csv(data.frame(Taxon=rownames(rel),rel,check.names=FALSE),file.path(out,'taxonomy',paste0(rank_name,'_sample_relative_abundance_unrarefied.csv')))
    all_long <- data.frame(Taxon=rep(rownames(rel),times=ncol(rel)),Abundance=as.vector(rel),Group=rep(pmd$Group,each=nrow(rel)))
    csv(stats::aggregate(Abundance~Taxon+Group,all_long,mean),file.path(out,'taxonomy',paste0(rank_name,'_group_mean_unrarefied.csv')))
    kept <- head(names(sort(rowMeans(rel),decreasing=TRUE)),15); plot_rel <- rel[kept,,drop=FALSE]
    rest <- setdiff(rownames(rel),kept); if (length(rest)) plot_rel <- rbind(plot_rel,Other=colSums(rel[rest,,drop=FALSE]))
    long <- data.frame(Taxon=rep(rownames(plot_rel),times=ncol(plot_rel)),SampleID=rep(colnames(plot_rel),each=nrow(plot_rel)),Abundance=as.vector(plot_rel),Group=rep(pmd$Group,each=nrow(plot_rel)))
    long$SampleID <- factor(long$SampleID,levels=colnames(plot_rel)); means <- stats::aggregate(Abundance~Taxon+Group,long,mean)
    csv(means,file.path(out,'taxonomy',paste0(rank_name,'_group_mean_top15_plot_input.csv')))
    for (kind in c('sample','group_mean')) {
      dat <- if (kind=='sample') long else means; dat$Label <- if (kind=='sample') dat$SampleID else factor(dat$Group,levels=unique(pmd$Group))
      p <- ggplot2::ggplot(dat,ggplot2::aes(Label,Abundance,fill=Taxon))+ggplot2::geom_col(width=.85)+theme_base()+
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))+ggplot2::labs(title=paste(rank_name,'composition (unrarefied filtered table)'),
          subtitle=if (kind=='sample') 'Top 15 taxa + Other; each sample normalized to 1' else 'Mean of per-sample relative abundances; top 15 + Other',x='',y='Relative abundance')
      save_plot(p,paste0('taxonomy/',rank_name,'_',kind,'_unrarefied'),max(8,min(20,nlevels(dat$Label)*.5+4)),6)
    }
  }
} else note('Taxonomy exports unavailable; composition plots omitted.')
writeLines(unique(notes),file.path(out,'render_notes.txt'),useBytes=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
jsonlite::write_json(list(original_assets=original,font=font,source='saved_QIIME2_outputs',notes=unique(notes)),file.path(out,'render_summary.json'),auto_unbox=TRUE,pretty=TRUE)
files <- list.files(out,recursive=TRUE,full.names=TRUE)
utils::write.table(data.frame(path=substring(files,nchar(run)+2),bytes=file.info(files)$size),file.path(out,'results.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
message('R figures completed: ',out)
