# Render saved QIIME outputs only: no rarefaction, ordination, or tests rerun in R.
args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop('Usage: Rscript render.R RUN_DIR FONT [PREPROCESS_DIR] [PLOTS_DIR] [GROUP_STYLE_JSON] [OUTPUT_DIR] [PLOT_SETTINGS_JSON] [PLOT_TARGET_JSON]')
run <- normalizePath(args[1], mustWork=TRUE)
font <- if (length(args)>1) args[2] else 'DejaVu Serif'
parent <- if (length(args)>2) normalizePath(args[3], mustWork=TRUE) else dirname(run)
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
assets <- if (length(args)>3 && nzchar(args[4])) args[4] else file.path(dirname(normalizePath(script)),'plots')
needed <- c('ggplot2','svglite','systemfonts','patchwork','xml2','jsonlite','MASS')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop('Missing R packages: ',paste(missing,collapse=', '))
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
out <- if (length(args)>5 && nzchar(args[6])) args[6] else result_path(run,'R_original_style')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,mustWork=TRUE)
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
style_file <- if (length(args)>4 && nzchar(args[5])) args[5] else file.path(run,'group_style.json')
if (length(args)>4 && nzchar(args[5]) && !file.exists(style_file)) stop('Group style JSON is missing: ',style_file)
requested_style <- if (file.exists(style_file)) jsonlite::read_json(style_file,simplifyVector=FALSE) else NULL
settings_file <- if (length(args)>6 && nzchar(args[7])) args[7] else file.path(run,'plot_settings.json')
if (length(args)>6 && nzchar(args[7]) && !file.exists(settings_file)) stop('Plot settings JSON is missing: ',settings_file)
plot_target <- NULL
if (length(args)>7 && nzchar(args[8])) {
  if (!file.exists(args[8])) stop('Plot target JSON is missing: ',args[8])
  plot_target <- jsonlite::read_json(args[8],simplifyVector=FALSE)
  patterns <- c(alpha='^alpha/(observed_features|shannon|evenness|faith_pd|chao1|simpson|goods_coverage)$',
                beta='^beta/PCoA_(weighted_unifrac|unweighted_unifrac|bray_curtis|jaccard|four_metrics)$',
                taxonomy='^taxonomy/(phylum|class|order|family|genus|species)_(sample|group_mean)_unrarefied$')
  if (!is.list(plot_target) || !is.character(plot_target$profile) || length(plot_target$profile)!=1 ||
      !plot_target$profile %in% names(patterns) || !is.character(plot_target$stem) || length(plot_target$stem)!=1 ||
      !grepl(patterns[[plot_target$profile]],plot_target$stem)) stop('Invalid single-figure plot target.')
}
want_profile <- function(profile) is.null(plot_target) || identical(plot_target$profile,profile)
want_figure <- function(stem,profile) is.null(plot_target) || (identical(plot_target$profile,profile) && identical(plot_target$stem,stem))
figures <- list()
plot_defaults <- list(
  general=list(dpi=NULL,base_size=NULL,title_size=NULL,axis_text_size=NULL,legend_text_size=NULL,
               legend_position='original',theme='original',background='white',show_caption=TRUE,formats='original',font_family='original'),
  alpha=list(width=NULL,height=NULL,title='',x_label='',y_label='',point_size=NULL,box_width=NULL,jitter_width=NULL,
             show_points=TRUE,show_letters=TRUE,show_outliers=FALSE,plot_type='boxplot'),
  beta=list(width=NULL,height=NULL,title='',x_label='',y_label='',point_size=NULL,point_alpha=NULL,
            show_ellipse=TRUE,ellipse_level=.95,ellipse_alpha=.15,show_sample_labels=FALSE,show_stats=TRUE,panel_columns=2),
  taxonomy=list(width=NULL,height=NULL,title='',x_label='',y_label='',top_n=20,top_mode='by_rank',levels=list(2,3,4,5,6),
                view='both',bar_width=.85,palette='user_R'))
plot_settings <- plot_defaults
requested_settings <- list()
if (file.exists(settings_file)) {
  requested <- jsonlite::read_json(settings_file,simplifyVector=FALSE)
  known_profiles <- c(names(plot_defaults),'heatmap','bubble','effect','lefse_bar','lefse_cladogram')
  if (!is.list(requested) || (length(requested) && is.null(names(requested))) || any(!names(requested) %in% known_profiles)) stop('Invalid plot settings sections.')
  requested_settings <- requested
  for (profile in intersect(names(requested),names(plot_defaults))) {
    section <- requested[[profile]]
    if (!is.list(section) || (length(section) && is.null(names(section))) || any(!names(section) %in% names(plot_defaults[[profile]]))) stop('Invalid plot settings: ',profile)
    plot_settings[[profile]][names(section)] <- section
  }
}
setting_number <- function(profile,key,lower,upper,nullable=FALSE,integer=FALSE) {
  x <- plot_settings[[profile]][[key]]
  if (is.null(x) && nullable) return(invisible(NULL))
  if (!is.numeric(x) || length(x)!=1 || !is.finite(x) || x<lower || x>upper || (integer && x!=floor(x))) stop('Invalid plot setting: ',profile,'.',key)
}
setting_choice <- function(profile,key,choices) {
  x <- plot_settings[[profile]][[key]]
  if (!is.character(x) || length(x)!=1 || is.na(x) || !x %in% choices) stop('Invalid plot setting: ',profile,'.',key)
}
for (profile in names(plot_settings)) {
  for (key in names(plot_settings[[profile]])) if (is.logical(plot_defaults[[profile]][[key]])) {
    x <- plot_settings[[profile]][[key]]
    if (!is.logical(x) || length(x)!=1 || is.na(x)) stop('Invalid plot setting: ',profile,'.',key)
  }
  if (profile!='general') {
    setting_number(profile,'width',2,30,nullable=TRUE)
    setting_number(profile,'height',2,40,nullable=TRUE)
    for (key in c('title','x_label','y_label')) {
      x <- plot_settings[[profile]][[key]]
      if (!is.character(x) || length(x)!=1 || is.na(x) || nchar(x)>300) stop('Invalid plot setting: ',profile,'.',key)
    }
  }
}
setting_number('general','dpi',72,600,nullable=TRUE,integer=TRUE)
setting_number('general','base_size',6,30,nullable=TRUE)
setting_number('general','title_size',6,40,nullable=TRUE)
for (key in c('axis_text_size','legend_text_size')) setting_number('general',key,5,30,nullable=TRUE)
setting_choice('general','legend_position',c('original','right','bottom','none'))
setting_choice('general','theme',c('original','classic','minimal','bw'))
setting_choice('general','background',c('white','transparent'))
setting_choice('general','formats',c('original','png_svg','png_svg_pdf'))
setting_choice('general','font_family',c('original','DejaVu Serif','DejaVu Sans','Times New Roman'))
if (plot_settings$general$font_family!='original') font <- plot_settings$general$font_family
if (!font %in% c('serif','sans','mono') && !font %in% systemfonts::system_fonts()$family)
  stop('Font not installed: ',font,'. Select an installed font in the workstation.')
setting_number('alpha','point_size',.1,12,nullable=TRUE)
setting_number('alpha','box_width',.1,1,nullable=TRUE)
setting_number('alpha','jitter_width',0,.5,nullable=TRUE)
setting_choice('alpha','plot_type',c('boxplot','violin'))
setting_number('beta','point_size',.1,12,nullable=TRUE)
setting_number('beta','point_alpha',.05,1,nullable=TRUE)
setting_number('beta','ellipse_level',.5,.99)
setting_number('beta','ellipse_alpha',0,1)
setting_number('beta','panel_columns',1,4,integer=TRUE)
setting_number('taxonomy','top_n',1,100,integer=TRUE)
setting_number('taxonomy','bar_width',.1,1)
setting_choice('taxonomy','view',c('both','sample','group'))
taxonomy_palettes <- jsonlite::read_json(file.path(dirname(normalizePath(script)),'taxonomy_palettes.json'))
setting_choice('taxonomy','palette',names(taxonomy_palettes))
setting_choice('taxonomy','top_mode',c('by_rank','custom','all'))
tax_levels <- unlist(plot_settings$taxonomy$levels,use.names=FALSE)
if (!is.numeric(tax_levels) || !length(tax_levels) || anyNA(tax_levels) || any(!tax_levels %in% 2:7) || anyDuplicated(tax_levels)) stop('Invalid plot setting: taxonomy.levels')
value_or <- function(x,default) if (is.null(x)) default else x
general <- plot_settings$general; alpha_options <- plot_settings$alpha; beta_options <- plot_settings$beta; taxonomy_options <- plot_settings$taxonomy
resolve_group_style <- function(actual,requested=NULL) {
  order <- if (is.null(requested$group_order)) character() else unlist(requested$group_order,use.names=FALSE)
  if (!is.character(order) || anyNA(order) || any(!nzchar(order)) || anyDuplicated(order)) stop('Invalid group display order.')
  order <- c(intersect(order,actual),setdiff(actual,order))
  palette <- c(CON='#D55E00',LDC='#0072B2',HDC='#CC79A7')
  extra <- setdiff(actual,names(palette))
  if (length(extra)) palette <- c(palette,setNames(grDevices::hcl.colors(length(extra),'Dark 3'),extra))
  custom <- if (is.null(requested$group_colors)) character() else unlist(requested$group_colors,use.names=TRUE)
  if (length(custom) && (is.null(names(custom)) || anyDuplicated(names(custom)) || anyNA(custom) || any(!grepl('^#[[:xdigit:]]{6}$',custom)))) stop('Invalid group display colors; use #RRGGBB.')
  palette[intersect(names(custom),actual)] <- custom[intersect(names(custom),actual)]
  list(group_order=order,group_colors=palette[order])
}
style <- resolve_group_style(unique(md$Group),requested_style)
groups <- style$group_order; colors <- style$group_colors
original <- all(file.exists(file.path(assets,c('qiime2_io.R','plot_alpha_original_style.R','plot_beta_original_style.R'))))
if (original) for (f in c('qiime2_io.R','plot_alpha_original_style.R','plot_beta_original_style.R')) source(file.path(assets,f),encoding='UTF-8')
if (!original) note('Original assets unavailable; embedded plots omit significance annotations.')
match_md <- function(ids) {
  if (anyDuplicated(ids) || anyNA(ids) || any(!nzchar(ids)) || !all(ids %in% md$SampleID)) stop('Saved result sample IDs invalid or absent from metadata.')
  d <- md[match(ids,md$SampleID),,drop=FALSE]; d$Group <- droplevels(factor(d$Group,levels=groups)); d
}
save_plot <- function(p,name,w,h,profile=NULL) {
  if (is.null(profile) && startsWith(name,'beta/')) profile <- 'beta'
  if (!want_figure(name,profile)) return(invisible(NULL))
  if (!is.null(profile) && name!='beta/PCoA_four_metrics') {
    w <- value_or(plot_settings[[profile]]$width,w)
    h <- value_or(plot_settings[[profile]]$height,h)
  }
  f <- file.path(out,name); dir.create(dirname(f),recursive=TRUE,showWarnings=FALSE)
  ggplot2::ggsave(paste0(f,'.svg'),p,width=w,height=h,device=svglite::svglite,bg=general$background)
  if (general$formats!='png_svg') ggplot2::ggsave(paste0(f,'.pdf'),p,width=w,height=h,device=grDevices::cairo_pdf,family=font,bg=general$background)
  else unlink(paste0(f,'.pdf'))
  ggplot2::ggsave(paste0(f,'.png'),p,width=w,height=h,dpi=value_or(general$dpi,300),type='cairo',bg=general$background)
  formats <- if (general$formats=='png_svg') c('png','svg') else c('png','svg','pdf')
  figures[[length(figures)+1L]] <<- list(stem=name,profile=profile,files=as.list(paste0(name,'.',formats)))
}
apply_plot_settings <- function(p,profile,metric='',rank='') {
  options <- plot_settings[[profile]]
  if (general$theme!='original') {
    make_theme <- switch(general$theme,classic=ggplot2::theme_classic,minimal=ggplot2::theme_minimal,bw=ggplot2::theme_bw)
    p <- p + make_theme(base_size=value_or(general$base_size,if (profile=='alpha') 11 else 14),base_family=font)
  }
  overrides <- list()
  if (!is.null(general$base_size)) for (key in c('text','axis.title','axis.title.x','axis.title.y','axis.text','axis.text.x','axis.text.y','legend.text','legend.title','strip.text','plot.title','plot.subtitle','plot.caption')) overrides[[key]] <- ggplot2::element_text(size=general$base_size)
  if (!is.null(general$title_size)) for (key in c('plot.title','strip.text')) overrides[[key]] <- ggplot2::element_text(size=general$title_size)
  if (!is.null(general$axis_text_size)) for (key in c('axis.text','axis.text.x','axis.text.y')) overrides[[key]] <- ggplot2::element_text(size=general$axis_text_size)
  if (!is.null(general$legend_text_size)) for (key in c('legend.text','legend.title')) overrides[[key]] <- ggplot2::element_text(size=general$legend_text_size)
  if (general$legend_position!='original') overrides$legend.position <- general$legend_position
  if (general$background=='transparent') for (key in c('plot.background','panel.background','legend.background','legend.box.background')) overrides[[key]] <- ggplot2::element_rect(fill='transparent',color=NA)
  if (length(overrides)) p <- p + do.call(ggplot2::theme,overrides)
  substitute_label <- function(x) gsub('{rank}',rank,gsub('{metric}',metric,x,fixed=TRUE),fixed=TRUE)
  labels <- list()
  for (key in c('title','x_label','y_label')) if (nzchar(options[[key]])) labels[[switch(key,x_label='x',y_label='y',title='title')]] <- substitute_label(options[[key]])
  if (!general$show_caption) { labels['caption'] <- list(NULL); labels['subtitle'] <- list(NULL) }
  if (length(labels)) p <- p + do.call(ggplot2::labs,labels)
  p
}
call_original_plot <- function(fun,arguments,options,profile) {
  supported <- names(formals(fun))
  unsupported <- if ('...' %in% supported) character() else setdiff(names(options),supported)
  for (key in unsupported) {
    control <- if (key=='ellipse') 'show_ellipse' else key
    if (!identical(plot_settings[[profile]][[control]],plot_defaults[[profile]][[control]]))
      note(paste('Custom plotting asset does not support',paste0(profile,'.',control),'; control omitted.'))
  }
  do.call(fun,c(arguments,options[setdiff(names(options),unsupported)]))
}
theme_base <- function() ggplot2::theme_minimal(base_size=14,base_family=font)+
  ggplot2::theme(panel.grid=ggplot2::element_blank(),panel.border=ggplot2::element_rect(color='black',fill=NA),axis.text=ggplot2::element_text(face='bold'),axis.title=ggplot2::element_text(face='bold'),strip.background=ggplot2::element_rect(fill='gray80',color='black'))
alpha_titles <- c(observed_features='Observed features',shannon='Shannon diversity',evenness='Pielou evenness',faith_pd='Faith phylogenetic diversity',chao1='Chao1 richness',simpson='Simpson diversity',goods_coverage="Good's coverage")
alpha <- list(); annotations <- list(); overall <- list(); pairwise <- list(); alpha_ids <- NULL
alpha_metrics <- if (!want_profile('alpha')) character() else if (is.null(plot_target)) names(alpha_titles) else sub('^alpha/','',plot_target$stem)
for (metric in alpha_metrics) {
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
    if (original) p <- call_original_plot(original_alpha_plot,list(d,annotation,font,colors,alpha_titles),
      list(point_size=value_or(alpha_options$point_size,1.5),box_width=value_or(alpha_options$box_width,.6),
           jitter_width=value_or(alpha_options$jitter_width,.2),show_points=alpha_options$show_points,
           show_letters=alpha_options$show_letters,show_outliers=alpha_options$show_outliers,plot_type=alpha_options$plot_type),'alpha')+
      ggplot2::labs(y=alpha_titles[[metric]]) else {
      p <- ggplot2::ggplot(d,ggplot2::aes(Group,Value,color=Group))
      p <- p + if (alpha_options$plot_type=='violin') ggplot2::geom_violin(fill=NA,width=value_or(alpha_options$box_width,.6)) else
        ggplot2::geom_boxplot(fill=NA,outlier.shape=if (alpha_options$show_outliers) 19 else NA,width=value_or(alpha_options$box_width,.6))
      if (alpha_options$show_points) p <- p + ggplot2::geom_jitter(position=ggplot2::position_jitter(width=value_or(alpha_options$jitter_width,.2),seed=123),size=value_or(alpha_options$point_size,1.5))
      p <- p + ggplot2::scale_color_manual(values=colors)+theme_base()+ggplot2::theme(legend.position='none')+ggplot2::labs(title=alpha_titles[[metric]],x='',y=alpha_titles[[metric]])
    }
    p <- apply_plot_settings(p,'alpha',metric=alpha_titles[[metric]])
    save_plot(p,paste0('alpha/',metric),if (original) 4.75 else 6,5,'alpha')
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
combined_target <- is.null(plot_target) || identical(plot_target$stem,'beta/PCoA_four_metrics')
beta_metrics <- if (!want_profile('beta')) character() else if (combined_target) names(titles) else sub('^beta/PCoA_','',plot_target$stem)
for (metric in beta_metrics) {
  # A zero distance matrix has no variance and therefore no defined PCoA
  # explained proportions. Preserve the other metrics and show the reason.
  dm <- read_dm(result_path(run,'exports','distance',metric,'distance-matrix.tsv'),md$SampleID)
  if (all(dm==0)) {
    reason <- 'PCoA undefined: all pairwise distances are zero.'
    note(paste(metric,reason))
    p <- ggplot2::ggplot()+ggplot2::annotate('text',x=0,y=0,label=reason,size=4,family=font)+
      ggplot2::theme_void(base_family=font)+ggplot2::labs(title=paste0('PCoA (',titles[[metric]],')'))
    p <- apply_plot_settings(p,'beta',metric=titles[[metric]])
    plots[[metric]] <- p; save_plot(p,paste0('beta/PCoA_',metric),6,5,'beta')
    csv(data.frame(Status='undefined_zero_distances',Reason=reason),file.path(out,'beta',paste0(metric,'_PCoA_status.csv')))
    next
  }
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
  if (original) p <- call_original_plot(original_pcoa_plot,list(pts,ord$proportion,if (beta_options$show_stats) saved_p else NULL,titles[[metric]],font,colors),
    list(ellipse=beta_options$show_ellipse,point_size=value_or(beta_options$point_size,2),point_alpha=beta_options$point_alpha,
         ellipse_level=beta_options$ellipse_level,ellipse_alpha=beta_options$ellipse_alpha,show_sample_labels=beta_options$show_sample_labels),'beta') else {
    point_layer <- if (is.null(beta_options$point_alpha)) ggplot2::geom_point(size=value_or(beta_options$point_size,2)) else
      ggplot2::geom_point(size=value_or(beta_options$point_size,2),alpha=beta_options$point_alpha)
    p <- ggplot2::ggplot(pts,ggplot2::aes(PC1,PC2,color=Group))+point_layer+
    ggplot2::geom_hline(yintercept=0,linetype='dashed',color='grey40')+ggplot2::geom_vline(xintercept=0,linetype='dashed',color='grey40')+
    ggplot2::scale_color_manual(values=colors)+theme_base()+ggplot2::labs(title=paste0('PCoA (',titles[[metric]],')'),
      x=sprintf('PCoA1 (%.1f%%)',ord$proportion[1]*100),y=sprintf('PCoA2 (%.1f%%)',ord$proportion[2]*100))+ggplot2::coord_equal()
    if (beta_options$show_sample_labels) p <- p + ggplot2::geom_text(ggplot2::aes(label=SampleID),vjust=-.8,size=3,family=font,show.legend=FALSE)
  }
  p <- apply_plot_settings(p,'beta',metric=titles[[metric]])
  plots[[metric]] <- p; save_plot(p,paste0('beta/PCoA_',metric),6,5,'beta')
  csv(pts,file.path(out,'beta',paste0(metric,'_PCoA_coordinates.csv')))
  csv(data.frame(Axis=paste0('PC',seq_along(ord$proportion)),ProportionExplained=ord$proportion),file.path(out,'beta',paste0(metric,'_PCoA_axes.csv')))
}
if (length(plots) && combined_target) {
  columns <- beta_options$panel_columns
  combined <- patchwork::wrap_plots(plots,ncol=columns)
  if (general$background=='transparent') combined <- combined + patchwork::plot_annotation(
    theme=ggplot2::theme(plot.background=ggplot2::element_rect(fill='transparent',color=NA)))
  combined_width <- value_or(beta_options$width,6)*columns
  combined_height <- value_or(beta_options$height,5)*ceiling(length(plots)/columns)
  if (!is.null(plot_target)) {
    combined_width <- value_or(beta_options$width,6*columns)
    combined_height <- value_or(beta_options$height,5*ceiling(length(plots)/columns))
  }
  save_plot(combined,'beta/PCoA_four_metrics',combined_width,combined_height,'beta')
}
if (length(beta_stats)) csv(do.call(rbind,beta_stats),file.path(out,'beta','beta_QIIME_overall_saved.csv'))
# Composition uses filtered UNRAREFIED counts and each sample's relative abundance.
tf <- result_path(parent,'exported-table','feature-table.tsv'); tx <- result_path(parent,'exported-taxonomy','taxonomy.tsv')
if (want_profile('taxonomy') && file.exists(tf) && file.exists(tx)) {
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
  composition_style <- resolve_group_style(unique(pmd$Group),requested_style)
  pmd$Group <- factor(pmd$Group,levels=composition_style$group_order)
  source(file.path(dirname(normalizePath(script)),'taxonomy_plot.R'),encoding='UTF-8')
  lineage <- tax$Taxon[match(tbl[[1]],tax$FeatureID)]
  rank_names <- c('kingdom','phylum','class','order','family','genus','species')
  if (!is.null(plot_target)) tax_levels <- match(sub('^taxonomy/([^_]+)_.*$','\\1',plot_target$stem),rank_names)
  for (rank in if (is.null(plot_target)) 2:7 else tax_levels) {
    rank_name <- rank_names[rank]
    cleaned <- taxonomy_clean_counts(counts,lineage,rank); rel <- cleaned$relative
    clean_dir <- file.path(out,'taxonomy','cleaned_tables',paste0('L',rank,'_',rank_name))
    dir.create(clean_dir,recursive=TRUE,showWarnings=FALSE)
    for (table_name in c('counts','relative')) {
      table <- data.frame(FeatureID=rownames(cleaned[[table_name]]),cleaned[[table_name]],check.names=FALSE)
      csv(table,file.path(clean_dir,paste0(table_name,'_cleaned.csv')))
      utils::write.table(table,file.path(clean_dir,paste0(table_name,'_cleaned.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
    }
    csv(data.frame(FeatureID=tbl[[1]],original_lineage=lineage,cleaned_taxon=cleaned$labels),file.path(clean_dir,'taxonomy_cleaning_mapping.csv'))
    jsonlite::write_json(list(level=rank,source='filtered_unrarefied',input_features=nrow(counts),
      cleaned_features=nrow(rel),per_sample_counts=as.list(colSums(counts)),
      per_sample_relative_sums=as.list(colSums(rel)),others_features=sum(cleaned$labels=='Others'),
      normalization='per-sample original filtered totals; group means give each sample equal weight',
      aggregation='requested rank name; missing or uncertain target rank becomes Others; originals retained'),
      file.path(clean_dir,'report.json'),pretty=TRUE,auto_unbox=TRUE)
    writeLines(c('清洗按当前层级的名称进行；未明确注释归为 Others，同名分类累加。',
      'counts_cleaned 保留未抽平整数计数；relative_cleaned 每列以该样本完整计数为分母重新计算。',
      'mapping 保存原 ASV 与原谱系，合并不会改变各样本的总 reads。',
      '完整清洗表中的 Others 仅为未识别分类；Top N 绘图表的 Others 还包括未入选的已识别分类。',
      '清洗表用于描述性丰度展示；ANCOM-BC2 继续使用完整谱系的未抽平计数，不把 Others 当作一个菌。'),
      file.path(clean_dir,'说明.md'),useBytes=TRUE)
    csv(data.frame(Taxon=rownames(rel),rel,check.names=FALSE),file.path(out,'taxonomy',paste0(rank_name,'_sample_relative_abundance_unrarefied.csv')))
    all_long <- data.frame(Taxon=rep(rownames(rel),times=ncol(rel)),Abundance=as.vector(rel),Group=rep(pmd$Group,each=nrow(rel)))
    csv(stats::aggregate(Abundance~Taxon+Group,all_long,mean),file.path(out,'taxonomy',paste0(rank_name,'_group_mean_unrarefied.csv')))
    if (!rank %in% tax_levels) next
    selection <- taxonomy_top_table(rel,rank,taxonomy_options); plot_rel <- selection$relative
    plot_taxa <- rownames(plot_rel)
    fill_colors <- taxonomy_fill_colors(plot_taxa,taxonomy_options$palette,taxonomy_palettes)
    suffix <- if (is.null(selection$limit)) 'all' else paste0('top',selection$limit)
    csv(data.frame(Taxon=plot_taxa,plot_rel,check.names=FALSE),file.path(out,'taxonomy',paste0(rank_name,'_sample_',suffix,'_plot_input.csv')))
    csv(data.frame(Taxon=plot_taxa,bottom_to_top=seq_along(plot_taxa),color=unname(fill_colors)),file.path(out,'taxonomy',paste0(rank_name,'_taxon_colors.csv')))
    long <- data.frame(Taxon=rep(rownames(plot_rel),times=ncol(plot_rel)),SampleID=rep(colnames(plot_rel),each=nrow(plot_rel)),Abundance=as.vector(plot_rel),Group=rep(pmd$Group,each=nrow(plot_rel)))
    sample_order <- pmd$SampleID[order(pmd$Group,seq_len(nrow(pmd)))]
    long$SampleID <- factor(long$SampleID,levels=sample_order); means <- stats::aggregate(Abundance~Taxon+Group,long,mean)
    means <- means[order(means$Group,match(means$Taxon,plot_taxa)),,drop=FALSE]
    csv(means,file.path(out,'taxonomy',paste0(rank_name,'_group_mean_',suffix,'_plot_input.csv')))
    jsonlite::write_json(list(level=rank,palette=taxonomy_options$palette,top_mode=taxonomy_options$top_mode,
      top_n=selection$limit,ranking=selection$ranking,bottom_to_top=as.list(plot_taxa),
      colors=as.list(fill_colors),others_includes=as.list(selection$other_taxa),
      value_unit='fraction (0 to 1); plot axis formatted as percent',grouping='saved metadata SampleID/Group',
      known_palette_extension='finite palettes retain first colors; extra taxa use HCL Dark 3; Others always #857D83'),
      file.path(out,'taxonomy',paste0(rank_name,'_plot_report.json')),pretty=TRUE,auto_unbox=TRUE,null='null')
    views <- switch(taxonomy_options$view,both=c('sample','group_mean'),sample='sample',group='group_mean')
    if (!is.null(plot_target)) views <- if (grepl('_group_mean_unrarefied$',plot_target$stem)) 'group_mean' else 'sample'
    for (kind in views) {
      dat <- if (kind=='sample') long else means; dat$Label <- if (kind=='sample') dat$SampleID else factor(dat$Group,levels=composition_style$group_order)
      dat$Taxon <- factor(dat$Taxon,levels=plot_taxa)
      markers <- unique(dat[c('Label','Group')])
      markers$Group <- factor(markers$Group,levels=composition_style$group_order)
      p <- ggplot2::ggplot(dat,ggplot2::aes(Label,Abundance,fill=Taxon))+ggplot2::geom_col(position=ggplot2::position_stack(reverse=TRUE),width=taxonomy_options$bar_width)+
        ggplot2::geom_point(data=markers,ggplot2::aes(x=Label,y=-.035,color=Group),inherit.aes=FALSE,shape=15,size=3)+
        ggplot2::scale_color_manual(values=composition_style$group_colors,breaks=composition_style$group_order)+
        ggplot2::scale_fill_manual(values=fill_colors,breaks=plot_taxa,drop=FALSE)+
        ggplot2::scale_y_continuous(breaks=seq(0,1,.2),labels=function(x) paste0(round(x*100),'%'))+theme_base()+
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))+ggplot2::labs(title=paste(rank_name,'composition (unrarefied filtered table)'),
          subtitle=paste0(if (is.null(selection$limit)) 'All identified taxa' else paste0('Top ',selection$limit,' identified taxa'),
            ' + Others; ',if (kind=='sample') 'per-sample relative abundance' else 'mean of per-sample relative abundances'),
          x='',y='Relative abundance (%)',color='Group',fill='Taxon')
      p <- apply_plot_settings(p,'taxonomy',rank=rank_name)
      save_plot(p,paste0('taxonomy/',rank_name,'_',kind,'_unrarefied'),max(8,min(20,nlevels(dat$Label)*.5+4)),6,'taxonomy')
    }
  }
} else if (want_profile('taxonomy')) note('Taxonomy exports unavailable; composition plots omitted.')
if (!is.null(plot_target) && length(figures)!=1L) stop('The requested cached figure could not be produced: ',plot_target$stem)
writeLines(unique(notes),file.path(out,'render_notes.txt'),useBytes=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
requested_settings[names(plot_settings)] <- plot_settings
jsonlite::write_json(list(original_assets=original,font=font,source='saved_QIIME2_outputs',input_run=run,
                         group_style=list(group_order=as.list(style$group_order),group_colors=as.list(style$group_colors)),
                         plot_settings=requested_settings,
                         plot_target=plot_target,figures=figures,
                         notes=unique(notes)),file.path(out,'render_summary.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
jsonlite::write_json(list(figures=figures),file.path(out,'figure_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
files <- list.files(out,recursive=TRUE,full.names=TRUE)
index_root <- if (length(args)>5 && nzchar(args[6])) out else run
utils::write.table(data.frame(path=substring(files,nchar(index_root)+2),bytes=file.info(files)$size),file.path(out,'results.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
message('R figures completed: ',out)
