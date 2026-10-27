# 16S 工作站方法审查与适用边界

本次审查针对现有源码及默认参数，对照 QIIME 2 2025.7 插件源码、DADA2/Cutadapt、ANCOM-BC2、LEfSe、PICRUSt2 与 vegan 官方资料。下表区分已经执行的保护、工作站经验规则和仍需实验设计判断的限制。它不代表对用户真实 FASTQ 或本机全部软件环境的认证。

## 审查结论

主线“按样本配对 → 有依据的去引物 → 处理后质量检查 → DADA2 降噪及去嵌合体 → 分类/过滤/建树 → 抽平多样性与未抽平丰度分析”符合这些工具的设计。不能把旧 OTU 流程的滑窗 Q20、10 bp 重叠、20% 错配、USEARCH+gold 去嵌合体参数直接套到 DADA2。DADA2 使用自身的质量过滤、误差学习、精确序列推断、配对合并和嵌合体识别；当前设置见每次任务参数与代码记录。

| 环节 | 本次核查与修正 | 解释和边界 |
|---|---|---|
| 样本/文件映射 | 四列显式映射与两列自动识别均拒绝同一真实 FASTQ 被复用；自动模式补上 inode/真实路径检查 | 防止链接造成伪重复；文件名相似不能证明实验样本独立 |
| 引物与附加序列 | 支持候选、自定义和 IUPAC；仅凭可靠匹配切除，未知引物后 barcode/spacer 不猜切；未检出不等于 clean | 现有清洗针对 5′ 引物及前方附加序列；短片段读穿的 3′ 对侧引物/接头未自动完整处理 |
| 快速 FASTQ 检查 | 保留前段记录、配对身份等输入检查 | clean 直通分支的快速检查不是全文件每条 ID/完整性证明；异常后段仍可能在下游报错。末尾 `.1/.2` 形式的部分 read ID 尚需转换或工具兼容确认 |
| 质量与截断建议 | 导出官方质量表、样本计数，结合长度/质量/已知扩增子重叠推荐；额外左端 trim 默认 0 | 全批抽样质量/长度覆盖不是每个样本的实际通过率；汇总图可能掩盖低深度样本问题。未知目标长度不凭质量图猜测 |
| DADA2 | 保留 maxEE、长度、误差学习、默认精确拼接及 consensus 去嵌合体；新增逐样本实际保留率 CSV/JSON 与界面表 | 相邻阶段保留率低于 50% 仅提示排查，不自动删除样本，也不是官方合格线。零分母留空 |
| 分类置信度 | 保留原浮点文本；允许仅边界舍入误差与官方禁用值 -1；拒绝实质越界/NaN | Confidence 不是丰度分母；`1.000000000000001` 不影响计数归一化 |
| 分类过滤 | 现有流程排除 Archaea、mitochondria、chloroplast；最低总丰度默认 1 | 是本项目目标选择，非官方统一必选项。并不等同于严格只保留 Bacteria；未分类/其他非目标需按实验核查。已排除古菌不能由后续预测恢复 |
| 系统发育树 | 保留 MAFFT → mask → FastTree → midpoint-root 的官方常见流程 | 短扩增子树具有分辨率限制；Faith PD/UniFrac 依赖树和完全匹配的特征集 |
| 抽平曲线 | 七项指标各自生成；请求 max-depth 超过实际最大计数时明确限至实际上限，并保存请求值、实际值和原因 | 曲线终点不是正式深度；极端最大计数可能使 10 个线性深度点过稀，应检查分布并手动调整 |
| 正式抽平 | 沿用用户指定“严格低于最低计数的最大整百数”，先预览保留样本；核心分析只抽平一次 | 这是保留全部样本的候选，不保证充分测序或曲线平台；不是 QIIME 官方最优公式。低质样本应按预先标准处理，不能按 P 值挑深度 |
| α 检验 | 修正某一常数配对导致整项指标被跳过的问题；全体常数明确跳过。无效检验保留 NA/诊断而非显著标签 | KW/Mann–Whitney 要求独立样本，小样本/大量并列值需谨慎；BH 按各指标的配对族执行，不是所有指标的共同 FDR |
| β 检验与图 | PERMANOVA 与 PERMDISP 分开报告；PERMDISP 任组少于 3 样本时跳过并说明；拒绝 NaN/Inf 统计量对应的假有效 P；PCoA 等比例坐标 | 分离不等于因果；组内离散度可能影响解释。上游对负特征值的处理需读 QIIME 日志，R 图不能恢复已丢失的负值；全零距离的上游核心流程未做真实 QIIME 验证 |
| 丰度/组成图 | 未抽平 ASV/L2–L7 计数与相对丰度；保留原表、清洗聚合表、Others、Top N 及作图表 | 图中 Others 合并/Top N 是展示操作；不能把合并后的模糊分类视为一个真实分类单元用于推断。相对量不是绝对菌量 |
| ANCOM-BC2 | 输入未抽平计数；按官方 prevalence 后 lib_cut 顺序估计样本，核对实际拟合样本；各层级保存 sample_inclusion.tsv/cohort.json，丰度图严格用实际模型样本 | 参考组影响系数方向；总体与所有两两比较另存。各层级单独分析可行，但门/纲特征少时偏差估计可能不稳；不把多个层级校正冒充共同 FDR |
| QIIME ANCOM | 保留 pseudocount → compositional ANCOM 的独立选项 | 与 ANCOM-BC2 不是同一算法；W 不是 LFC、LDA 或概率，也不产生相同的两两效应 |
| LEfSe | 调用已安装官方命令；Group 为 class，无虚构 subclass/subject；保存得分，预览不重拟合 | KW、官方一致性步骤、LDA 分阶段解释；官方小 subclass 采用中位数规则，不能保证每个 pair 都做了 Wilcoxon P<0.05。多组 one-against-all 与 all-against-all 不同；LDA 阈值不是多重校正 |
| PICRUSt2 | 未抽平 BIOM+匹配序列；保存版本、NSTI/比对过滤与用户 KO 映射来源 | 功能是参考基因组推断。KO→KEGG 求和是注释聚合，不能称为正式通路重建或活性；MetaCyc 官方推断另列；当前没有菌种贡献分层输出 |
| 新高级方法 | vegan：核心检出率、CLR-PCA、Bray PERMANOVA/PERMDISP、Hellinger-RDA/Bray-dbRDA；ppcor：表型偏 Spearman，另可选普通 Spearman | 原始计数与全表分母分别保留。关联默认 CLR，显式报告伪计数和筛选子组成；BH 覆盖本次全部 taxa×traits。相关/约束排序均不证明因果 |
| 设计保护 | 元数据按 SampleID 字面匹配；检查缺失、冲突、模型秩和自由度；检测重复 AnimalID/SubjectID 后高级推断跳过，仅保留描述分析 | 没有 ID 列不证明独立；配对、笼位、纵向/重复测量和复杂混杂尚需相应模型，不能强行按独立样本检验 |
| 重绘与复现 | 图后独立改尺寸/配色/顺序，读取缓存并另存；高级缓存 SHA 校验；保存真实外部命令、源码和参数 | 修改图形不改变统计；修改样本/模型/过滤标准需新分析。随机抽平和置换重跑可能不同，源码相同不等于逐数值复现 |
| 独立输入与恢复 | 外部计数/QZA按模块核验类型、样本/特征覆盖，保留原件与对齐报告；新旧任务明确区分；失败只复用已核验检查点 | 整数值不能证明未抽平，需用户确认来源；外部输入完成不表示做过DADA2或此前质控。旧版本无指纹时不猜测跳过阶段，失败源目录不覆盖 |

## 本次验证范围

验证分为样本与文件保护、命令编排、真实 R 合成数据计算/图形、HTTP/浏览器与 Windows 部署脚本检查。没有在云端重跑用户的完整 QIIME 2、ANCOM-BC2、LEfSe 生物分析，也没有 Windows/WSL2 桌面实机。具体执行数量及结果见 [VALIDATION.md](../VALIDATION.md)。新手解释、输入示例和决策步骤见 [USER_GUIDE.md](USER_GUIDE.md)。

## 优先参考的官方资料与原始论文

1. QIIME 2 2025.7：[q2-dada2 包装参数](https://github.com/qiime2/q2-dada2/blob/2025.7.0/q2_dada2/_denoise.py)、[实际 R 调用](https://github.com/qiime2/q2-dada2/blob/2025.7.0/q2_dada2/assets/run_dada.R)、[分类过滤](https://github.com/qiime2/q2-taxa/blob/2025.7.0/q2_taxa/_method.py)、[抽平曲线与 α 统计](https://github.com/qiime2/q2-diversity/blob/2025.7.0/q2_diversity/_alpha/_visualizer.py)、[β 统计](https://github.com/qiime2/q2-diversity/blob/2025.7.0/q2_diversity/_beta/_visualizer.py)、[core-metrics](https://github.com/qiime2/q2-diversity/blob/2025.7.0/q2_diversity/_core_metrics.py)。
2. [Cutadapt 官方指南](https://cutadapt.readthedocs.io/en/stable/guide.html)：配对处理、简并碱基、锚定/位置与 read-through；[DADA2 官方教程](https://benjjneb.github.io/dada2/tutorial.html)；Callahan 等，2016，DADA2，[doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)。
3. [QIIME Moving Pictures 官方教程](https://docs.qiime2.org/2024.10/tutorials/moving-pictures/)与[持续更新源码](https://github.com/qiime2/moving-pictures-tutorial/blob/main/book/tutorial.md)：抽平深度与样本保留的权衡。旧版本文档用于概念解释，不替代 2025.7 已安装 CLI 参数。
4. ANCOM-BC2：[官方 vignette](https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html)、[实际预处理顺序](https://github.com/FrederickHuangLin/ANCOMBC/blob/RELEASE_3_20/R/ancombc_prep.R)、Lin 与 Peddada，2024，[doi:10.1038/s41592-023-02092-7](https://doi.org/10.1038/s41592-023-02092-7)。
5. LEfSe：[官方仓库](https://github.com/biobakery/lefse)、[1.1.2 实现](https://github.com/biobakery/lefse/blob/1.1.2/lefse/lefse.py)；Segata 等，2011，[doi:10.1186/gb-2011-12-6-r60](https://doi.org/10.1186/gb-2011-12-6-r60)。
6. PICRUSt2：[官方流程](https://github.com/picrust/picrust2/wiki/Full-pipeline-script)、[FAQ](https://github.com/picrust/picrust2/wiki/Frequently-Asked-Questions)；Douglas 等，2020，[doi:10.1038/s41587-020-0548-6](https://doi.org/10.1038/s41587-020-0548-6)。
7. vegan 官方：[adonis2](https://vegandevs.github.io/vegan/reference/adonis.html)、[betadisper](https://vegandevs.github.io/vegan/reference/betadisper.html)、[rda/dbrda](https://vegandevs.github.io/vegan/reference/cca.html)、[约束排序检验](https://vegandevs.github.io/vegan/reference/anova.cca.html)。Legendre 与 Gallagher，2001，[doi:10.1007/s004420100716](https://doi.org/10.1007/s004420100716)。
8. [ppcor 官方包](https://cran.r-project.org/package=ppcor)；Kim，2015，[doi:10.5351/CSAM.2015.22.6.665](https://doi.org/10.5351/CSAM.2015.22.6.665)。Gloor 等，2017，组成数据解释，[doi:10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224)。

引用用于对应方法和限制，不表示某篇论文证明工作站所有默认阈值普适。核心出现率、质量截断启发式和整百推荐深度均须结合实验判断。
