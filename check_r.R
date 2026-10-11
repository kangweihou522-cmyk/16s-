args <- commandArgs(trailingOnly=TRUE)
prefix <- Sys.getenv('CONDA_PREFIX')
if (!nzchar(prefix)) stop('绘图应在专用conda R环境中运行；请执行新版安装依赖.sh')
prefix <- normalizePath(prefix,mustWork=TRUE)
paths <- normalizePath(.libPaths(),mustWork=TRUE)
if (any(!startsWith(paste0(paths,'/'),paste0(prefix,'/')))) stop('R混入环境外包库: ',paste(paths[!startsWith(paste0(paths,'/'),paste0(prefix,'/'))],collapse=', '))
needed <- c('ggplot2','patchwork','svglite','xml2','jsonlite','systemfonts','multcompView','MASS')
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('专用R环境中缺少/无法加载包: ',paste(missing,collapse=', '),'; 请重跑新版安装依赖.sh')
if (!capabilities('cairo')) stop('专用R环境缺少Cairo图形支持')
font <- if (length(args)) args[1] else 'DejaVu Serif'
if (!tolower(font) %in% tolower(systemfonts::system_fonts()$family)) stop('字体未安装: ',font,'; 安装字体或在界面选择DejaVu Serif')
cat('R_HOME: ',R.home(),'\nR libraries: ',paste(paths,collapse='; '),'\n',sep='')
print(sessionInfo())
