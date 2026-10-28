# Render cached KO differential results at the source KGML coordinates only.
# CLI: Rscript --vanilla ko_map_render.R --render OPTIONS_JSON OUTPUT_DIR
# No statistical model, database query or network operation is performed here.
ko_map_render <- function(options_path, output_dir) {
  needed <- c('ggplot2', 'jsonlite', 'svglite', 'systemfonts')
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop('Missing R packages: ', paste(missing, collapse=', '))
  opts <- jsonlite::read_json(options_path, simplifyVector=FALSE)
  scalar <- function(x, default=NULL) if (is.null(x)) default else x
  scalar_text <- function(x, name, default='') {
    x <- scalar(x, default)
    if (!is.character(x) || length(x)!=1L || is.na(x) || nchar(x)>1000) stop('Invalid ', name)
    x
  }
  number <- function(x, name, default, lower=0, upper=Inf) {
    x <- scalar(x, default)
    if (!is.numeric(x) || length(x)!=1L || !is.finite(x) || x<lower || x>upper) stop('Invalid ', name)
    x
  }
  flag <- function(x, name, default=TRUE) {
    x <- scalar(x, default)
    if (!is.logical(x) || length(x)!=1L || is.na(x)) stop('Invalid ', name)
    x
  }
  geometry_path <- scalar_text(opts$geometry_path, 'geometry_path')
  hits_path <- scalar_text(opts$hits_path, 'hits_path')
  if (!nzchar(geometry_path) || !nzchar(hits_path)) stop('geometry_path and hits_path are required cached local files')
  geometry <- jsonlite::read_json(geometry_path, simplifyVector=FALSE)
  if (!identical(geometry$schema_version, 1L) && !identical(geometry$schema_version, 1)) stop('Unsupported KGML geometry schema')
  if (!identical(geometry$coordinate_system, 'KGML_original_y_down') ||
      !grepl('^ko[0-9]{5}$', scalar_text(geometry$map_id, 'map_id'))) stop('Invalid original KO KGML geometry')
  entries <- geometry$entries
  members <- function(entry) unlist(entry$ko_ids, use.names=FALSE)
  ko_entries <- Filter(function(entry) length(members(entry))>0L && length(entry$graphics)>0L, entries)
  if (!length(ko_entries)) stop('KGML contains no valid KO ortholog graphics')
  if (any(!grepl('^K[0-9]{5}$', unlist(lapply(ko_entries, members), use.names=FALSE)))) stop('Invalid KO geometry identifiers')
  map_kos <- unique(unlist(lapply(ko_entries, members), use.names=FALSE))

  table <- utils::read.delim(hits_path, check.names=FALSE, colClasses='character',
                            quote='', comment.char='', na.strings=c('', 'NA', 'NaN'), fileEncoding='UTF-8-BOM')
  if (!all(c('KO', 'effect', 'q', 'comparison') %in% names(table))) stop('Cached differential TSV requires KO, effect, q, comparison columns')
  table$KO <- toupper(sub('^ko:', '', trimws(table$KO), ignore.case=TRUE))
  if (anyNA(table$KO) || any(!grepl('^K[0-9]{5}$', table$KO)) || anyNA(table$comparison) || any(!nzchar(table$comparison))) stop('Invalid cached KO identifiers or comparison')
  for (key in c('effect', 'q')) {
    raw <- table[[key]]
    table[[key]] <- suppressWarnings(as.numeric(raw))
    if (any(!is.na(raw) & (is.na(table[[key]]) | !is.finite(table[[key]])))) stop('Invalid cached ', key)
  }
  if (any(!is.na(table$q) & (table$q<0 | table$q>1))) stop('Cached q values must be in [0,1]')
  comparison <- scalar_text(opts$comparison, 'comparison')
  if (!nzchar(comparison) && length(unique(table$comparison))==1L) comparison <- unique(table$comparison)
  if (!nzchar(comparison) || nrow(table)>0L && !comparison %in% table$comparison) stop('Select an available cached comparison')
  selected <- table[table$comparison==comparison,,drop=FALSE]
  # Repeated identical rows are retained in row reports but count once as a KO.
  # Conflicting copies would make the direction ambiguous and are rejected.
  repeated <- split(seq_len(nrow(selected)), selected$KO)
  for (indices in repeated) if (nrow(unique(selected[indices,c('effect','q'),drop=FALSE]))>1L)
    stop('Conflicting cached differential rows for KO ', selected$KO[indices[1]], ' in ', comparison)
  unique_selected <- selected[!duplicated(selected$KO),,drop=FALSE]

  requested <- scalar(opts$plot_settings, list())
  general <- list(dpi=NULL, font_family='original', base_size=NULL, title_size=NULL, axis_text_size=NULL,
                  legend_text_size=NULL, legend_position='original', theme='original', background='white',
                  show_caption=TRUE, formats='original')
  profile <- list(width=NULL, height=NULL, title='', up_color='#D73027', down_color='#4575B4',
                  mixed_color='#984EA3', background_color='#C8CDD0', line_width=1,
                  show_labels=TRUE, q_threshold=.05, effect_threshold=0)
  for (key in intersect(names(requested$general), names(general))) general[key] <- requested$general[key]
  for (key in intersect(names(requested$metabolic_map), names(profile))) profile[key] <- requested$metabolic_map[key]
  # Top-level parameters support initial rendering; saved profile settings own previews.
  for (key in names(profile)) if (!is.null(opts[[key]]) && !key %in% names(requested$metabolic_map)) profile[key] <- opts[key]
  for (key in c('up_color','down_color','mixed_color','background_color'))
    if (!is.character(profile[[key]]) || length(profile[[key]])!=1L || !grepl('^#[[:xdigit:]]{6}$', profile[[key]])) stop('Use #RRGGBB plot color: ', key)
  q_threshold <- number(profile$q_threshold, 'q_threshold', .05, 0, 1)
  effect_threshold <- number(profile$effect_threshold, 'effect_threshold', 0, 0)
  line_width <- number(profile$line_width, 'line_width', 1, .1, 10)
  show_labels <- flag(profile$show_labels, 'show_labels')
  show_title <- flag(opts$show_title, 'show_title')
  show_caption <- flag(general$show_caption, 'show_caption')
  font <- scalar_text(opts$font, 'font', 'sans')
  if (!identical(general$font_family, 'original')) font <- scalar_text(general$font_family, 'font_family')
  if (!font %in% c('serif','sans','mono') && !font %in% systemfonts::system_fonts()$family) stop('Font not installed: ',font)
  dpi <- number(general$dpi, 'dpi', scalar(opts$dpi, 200), 72, 600)
  if (dpi!=floor(dpi)) stop('DPI must be an integer')
  base_size <- number(general$base_size, 'base_size', 11, 6, 30)
  title_size <- number(general$title_size, 'title_size', base_size+2, 6, 40)
  legend_text_size <- number(general$legend_text_size, 'legend_text_size', base_size, 5, 30)
  if (!general$legend_position %in% c('original','right','bottom','none') ||
      !general$theme %in% c('original','classic','minimal','bw') ||
      !general$background %in% c('white','transparent') ||
      !general$formats %in% c('original','png_svg','png_svg_pdf')) stop('Invalid general metabolic map settings')
  legend_position <- if (general$legend_position=='original') 'bottom' else general$legend_position
  formats <- if (general$formats=='png_svg_pdf') c('png','svg','pdf') else c('png','svg')
  significant <- !is.na(unique_selected$q) & !is.na(unique_selected$effect) &
    unique_selected$q<=q_threshold & abs(unique_selected$effect)>=effect_threshold & unique_selected$effect!=0
  hit_kos <- unique_selected$KO[significant]
  up_kos <- unique_selected$KO[significant & unique_selected$effect>0]
  down_kos <- unique_selected$KO[significant & unique_selected$effect<0]
  join_ids <- function(x) paste(x, collapse=';')
  report_rows <- lapply(ko_entries, function(entry) {
    all <- members(entry); tested <- intersect(all, unique_selected$KO)
    up <- intersect(all, up_kos); down <- intersect(all, down_kos); hits <- intersect(all, hit_kos)
    status <- if (length(up) && length(down)) 'mixed' else if (length(up)) 'up' else if (length(down)) 'down' else 'background'
    data.frame(entry_id=entry$id, member_kos=join_ids(all), member_count=length(all),
               tested_kos=join_ids(tested), tested_count=length(tested), hit_kos=join_ids(hits), hit_count=length(hits),
               up_kos=join_ids(up), down_kos=join_ids(down), status=status,
               partial_hit=length(hits)>0L && length(hits)<length(all), stringsAsFactors=FALSE)
  })
  overlap <- do.call(rbind, report_rows)
  group_ids <- function(ko) join_ids(overlap$entry_id[vapply(ko_entries, function(entry) ko %in% members(entry), logical(1))])
  selected$map_entry_ids <- vapply(selected$KO, group_ids, character(1))
  selected$passes_display_threshold <- !is.na(selected$q) & !is.na(selected$effect) & selected$q<=q_threshold &
    abs(selected$effect)>=effect_threshold & selected$effect!=0
  mapped <- selected[nzchar(selected$map_entry_ids),,drop=FALSE]
  unmapped <- selected[!nzchar(selected$map_entry_ids),,drop=FALSE]

  # Coordinates below are copied from KGML; relation metadata supplies no paths.
  paths <- list(); boxes <- list(); circles <- list(); labels <- list()
  add <- function(collection, row) { collection[[length(collection)+1L]] <- row; collection }
  for (entry in entries) {
    position <- match(entry$id, overlap$entry_id)
    status <- if (is.na(position)) 'context' else overlap$status[position]
    for (graphic in entry$graphics) {
      x <- graphic$x; y <- graphic$y; width <- graphic$width; height <- graphic$height
      if (any(!is.finite(c(x,y,width,height)))) stop('Invalid cached KGML graphics coordinates')
      row <- data.frame(id=graphic$id, x=x, y=y, status=status, stringsAsFactors=FALSE)
      if (graphic$type=='line') {
        xy <- do.call(rbind, lapply(graphic$coords, unlist, use.names=FALSE))
        if (nrow(xy)<2L || ncol(xy)!=2L || any(!is.finite(xy))) stop('Invalid cached KGML line coordinates')
        paths <- add(paths, data.frame(id=graphic$id, x=xy[,1], y=xy[,2], status=status))
      } else if (graphic$type %in% c('rectangle','roundrectangle')) {
        row$xmin <- x-width/2; row$xmax <- x+width/2; row$ymin <- y-height/2; row$ymax <- y+height/2
        boxes <- add(boxes,row)
      } else if (graphic$type=='circle') {
        theta <- seq(0, 2*pi, length.out=49)
        circles <- add(circles, data.frame(id=graphic$id,x=x+cos(theta)*width/2,y=y+sin(theta)*height/2,status=status))
      } else stop('Unsupported cached KGML graphics type')
      label <- scalar_text(graphic$label, 'KGML label')
      if (show_labels && nzchar(label) && graphic$type!='line') {
        row$label <- label; labels <- add(labels,row[c('id','x','y','status','label')])
      }
    }
  }
  bind <- function(x) if (length(x)) do.call(rbind,x) else NULL
  paths <- bind(paths); boxes <- bind(boxes); circles <- bind(circles); labels <- bind(labels)
  palette <- c(up=profile$up_color, down=profile$down_color, mixed=profile$mixed_color,
               background=profile$background_color, context='#E1E4E6')
  legend_labels <- c(up='At least one up KO', down='At least one down KO', mixed='Both up and down KOs',
                     background='No KO passes thresholds')
  levels <- c('up','down','mixed','background')
  fig <- ggplot2::ggplot()
  # Separate layers keep threshold-passing groups above coincident backgrounds.
  for (highlight in c(FALSE,TRUE)) {
    layer_rows <- function(d) if (!is.null(d)) d[(d$status %in% c('up','down','mixed'))==highlight,,drop=FALSE] else NULL
    layer_paths <- layer_rows(paths); layer_boxes <- layer_rows(boxes); layer_circles <- layer_rows(circles)
    if (!is.null(layer_paths) && nrow(layer_paths)) fig <- fig + ggplot2::geom_path(data=layer_paths,
        ggplot2::aes(x=x,y=y,group=id,color=status),linewidth=line_width,lineend='round',show.legend=TRUE)
    if (!is.null(layer_boxes) && nrow(layer_boxes)) fig <- fig + ggplot2::geom_rect(data=layer_boxes,
        ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=status,color=status), linewidth=line_width*.45,show.legend=FALSE)
    if (!is.null(layer_circles) && nrow(layer_circles)) fig <- fig + ggplot2::geom_polygon(data=layer_circles,
        ggplot2::aes(x=x,y=y,group=id,fill=status,color=status), linewidth=line_width*.45,show.legend=FALSE)
  }
  if (!is.null(labels)) fig <- fig + ggplot2::geom_text(data=labels,
      ggplot2::aes(x=x,y=y,label=label), family=font,size=base_size/ggplot2::.pt*.72,color='#444444',show.legend=FALSE)
  # Train an explicit legend even for box-only maps or zero-hit comparisons.
  legend_data <- data.frame(x=rep(NA_real_,4),y=rep(NA_real_,4),status=factor(levels,levels=levels))
  fig <- fig + ggplot2::geom_point(data=legend_data,ggplot2::aes(x=x,y=y,color=status),na.rm=TRUE,show.legend=TRUE) +
    ggplot2::scale_color_manual(values=palette,breaks=levels,labels=legend_labels,limits=c(levels,'context'),drop=FALSE,name='KO feature-group overlap') +
    ggplot2::scale_fill_manual(values=palette,limits=c(levels,'context'),guide='none')
  fig <- fig + ggplot2::guides(color=ggplot2::guide_legend(title.position='top',byrow=TRUE,
      nrow=if(legend_position=='bottom') 2 else NULL,ncol=if(legend_position=='bottom') NULL else 1))
  bounds <- geometry$bounds
  extent <- unlist(bounds[c('xmin','xmax','ymin','ymax')],use.names=FALSE)
  if (length(extent)!=4L || any(!is.finite(extent)) || extent[2]<extent[1] || extent[4]<extent[3]) stop('Invalid cached KGML extent')
  span_x <- max(1, extent[2]-extent[1]); span_y <- max(1, extent[4]-extent[3])
  width <- number(profile$width, 'width', min(30,max(6,span_x/120)), 2, 30)
  height <- number(profile$height, 'height', min(40,max(4,width*span_y/span_x+1.1)), 2, 40)
  title <- scalar_text(profile$title, 'title')
  if (!nzchar(title)) title <- paste(geometry$map_id, '—', scalar_text(geometry$title, 'KGML title'))
  comparison_label <- scalar_text(opts$comparison_label, 'comparison_label', comparison)
  effect_label <- scalar_text(opts$effect_label, 'effect_label', 'cached effect')
  caption <- paste0('Cached comparison: ',comparison_label,'; q ≤ ',format(q_threshold),'; |',effect_label,'| ≥ ',format(effect_threshold),
                    '. ',length(intersect(hit_kos,map_kos)),' / ',length(map_kos),' map KOs highlighted.\n',
                    'A colored feature group has at least one matching KO; this does not establish activity of the complete reaction.\n',
                    'Mixed = qualifying KOs in both directions. Original KGML coordinates; relationships without paths are not drawn.')
  if (!length(hit_kos)) caption <- paste0('No KOs pass the display thresholds.\n',caption)
  if (length(hit_kos) && !length(intersect(hit_kos,map_kos))) caption <- paste0('No threshold-passing KOs overlap this map.\n',caption)
  # A fixed-aspect panel can occupy less than the full canvas width.  Wrap for
  # that narrower panel so pathway captions also fit a wide, short export.
  caption_width <- min(width,max(1.5,(height-1.8)*span_x/span_y))
  caption <- paste(vapply(strsplit(caption,'\n',fixed=TRUE)[[1]],function(line)
      paste(strwrap(line,width=max(24,floor(caption_width*105/max(6,base_size-3)))),collapse='\n'),character(1)),collapse='\n')
  theme <- switch(general$theme,classic=ggplot2::theme_classic,minimal=ggplot2::theme_minimal,bw=ggplot2::theme_bw,ggplot2::theme_void)
  fig <- fig + ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed(ratio=1,xlim=extent[1:2]+c(-1,1)*span_x*.015,ylim=rev(extent[3:4])+c(1,-1)*span_y*.015,clip='off') +
    theme(base_size=base_size,base_family=font) +
    ggplot2::theme(axis.text=ggplot2::element_blank(),axis.title=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),
      panel.grid=ggplot2::element_blank(),legend.position=legend_position,legend.text=ggplot2::element_text(size=legend_text_size),
      plot.title=ggplot2::element_text(size=title_size),plot.caption=ggplot2::element_text(hjust=0,size=max(6,base_size-3)),
      plot.title.position='plot',plot.caption.position='plot',plot.margin=ggplot2::margin(8,12,10,12),
      plot.background=ggplot2::element_rect(fill=general$background,color=NA),
      panel.background=ggplot2::element_rect(fill=general$background,color=NA),
      legend.background=ggplot2::element_rect(fill=general$background,color=NA),
      legend.key=ggplot2::element_rect(fill=general$background,color=NA)) +
    ggplot2::labs(title=if(show_title) title else NULL,caption=if(show_caption) caption else NULL)

  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  figure_dir <- file.path(output_dir,'03_figures'); dir.create(figure_dir,showWarnings=FALSE)
  stem <- file.path(figure_dir,'metabolic_map')
  for (format in formats) {
    arguments <- list(filename=paste0(stem,'.',format),plot=fig,width=width,height=height,dpi=dpi,bg=general$background,limitsize=FALSE)
    if (format=='pdf') {
      if (!capabilities('cairo')) stop('Cairo graphics are required for font-compatible PDF')
      arguments$device <- grDevices::cairo_pdf
    }
    do.call(ggplot2::ggsave,arguments)
  }
  mapping_dir <- file.path(output_dir,'02_mapping'); dir.create(mapping_dir,showWarnings=FALSE)
  write_tsv <- function(x,name) utils::write.table(x,file.path(mapping_dir,name),sep='\t',row.names=FALSE,quote=FALSE,na='',fileEncoding='UTF-8')
  write_tsv(mapped,'mapped_ko_rows.tsv'); write_tsv(unmapped,'unmapped_ko_rows.tsv'); write_tsv(overlap,'feature_group_overlap.tsv')
  estimable_kos <- unique(table$KO[is.finite(table$effect) & is.finite(table$q)])
  coverage <- list(input_rows_all_comparisons=nrow(table),all_input_kos=length(unique(table$KO)),all_tested_kos=length(estimable_kos),
                   mapped_all_tested_kos=length(intersect(estimable_kos,map_kos)),unmapped_all_tested_kos=length(setdiff(estimable_kos,map_kos)),
                   selected_rows=nrow(selected),selected_input_kos=nrow(unique_selected),
                   selected_tested_kos=sum(is.finite(unique_selected$effect) & is.finite(unique_selected$q)),
                   duplicate_selected_rows=nrow(selected)-nrow(unique_selected),
                   map_kos=length(map_kos),map_feature_groups=nrow(overlap),
                   mapped_selected_kos=length(intersect(unique_selected$KO,map_kos)),
                   unmapped_selected_kos=length(setdiff(unique_selected$KO,map_kos)),
                   selected_significant_kos=length(hit_kos),mapped_significant_kos=length(intersect(hit_kos,map_kos)),
                   unmapped_significant_kos=length(setdiff(hit_kos,map_kos)),
                   highlighted_feature_groups=sum(overlap$hit_count>0L),partial_hit_feature_groups=sum(overlap$partial_hit),
                   mixed_feature_groups=sum(overlap$status=='mixed'))
  effective <- list(profile='metabolic_map',width=width,height=height,dpi=dpi,font=font,formats=as.list(formats),
                    general=general,profile_settings=profile,comparison=comparison,
                    q_threshold=q_threshold,effect_threshold=effect_threshold)
  jsonlite::write_json(list(requested=requested,effective=effective),paste0(stem,'_plot_settings.json'),pretty=TRUE,auto_unbox=TRUE,null='null')
  summary <- list(source='ko_network',statistics_refit=FALSE,map_id=geometry$map_id,comparison=comparison,
                  coordinate_system=geometry$coordinate_system,geometry_sha256=geometry$source_sha256,
                  thresholds_display_only=TRUE,coverage=coverage,
                  feature_group_rule='At least one member KO passes q and absolute-effect thresholds with nonzero direction; mixed requires both up and down members. Partial groups do not establish complete reaction activity.',
                  comparison_label=comparison_label,effect_label=effect_label,effect_metadata=opts$effect_metadata,plot_settings=effective,
                  reports=list(mapped='02_mapping/mapped_ko_rows.tsv',unmapped='02_mapping/unmapped_ko_rows.tsv',
                               feature_groups='02_mapping/feature_group_overlap.tsv'))
  jsonlite::write_json(summary,file.path(output_dir,'render_summary.json'),pretty=TRUE,auto_unbox=TRUE,null='null')
  jsonlite::write_json(list(figures=list(list(stem='03_figures/metabolic_map',profile='metabolic_map',
    files=as.list(paste0('03_figures/metabolic_map.',formats))))),file.path(output_dir,'figure_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
  invisible(summary)
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args)!=3L || args[1]!='--render') stop('Usage: Rscript ko_map_render.R --render OPTIONS_JSON OUTPUT_DIR')
  ko_map_render(args[2],args[3])
}
