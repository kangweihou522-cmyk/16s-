args <- commandArgs(trailingOnly=TRUE)
needed <- c('ggplot2','patchwork','svglite','xml2','jsonlite','systemfonts','multcompView')
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('缺少R包: ', paste(missing,collapse=', '), '; 请先运行安装依赖.sh')
if (!tolower(args[1]) %in% tolower(systemfonts::system_fonts()$family)) stop('字体未安装: ',args[1], '; 请安装字体或在界面明确选择 DejaVu Serif')
print(sessionInfo())
