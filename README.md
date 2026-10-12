# Windows / WSL 16S 本地分析工作站 v0.4

浏览器操作，WSL Ubuntu 中使用你自己安装的 **QIIME 2 和 R** 计算。原始测序数据只读，每次生成独立结果目录；分类器由用户选择本地 `.qza` 文件。程序不调用任何 AI API，不需保持 ChatGPT/Codex 会话运行。

## 获取、启动和更新

在 Ubuntu 中执行：

```bash
git clone https://github.com/kangweihou522-cmyk/16s-.git ~/16S-workstation-v02
cd ~/16S-workstation-v02
bash 启动.sh
```

Windows 浏览器打开 `http://localhost:8765`。保持 Ubuntu 终端打开；Ctrl+C 停止。若原基础版仍运行，先停止它，避免端口占用。本版本安装在新目录，保留 `~/16S-workstation-basic` 旧版。

若端口被其他程序占用：`bash 启动.sh --port 8766`，对应打开 `http://localhost:8766`。

更新源码后必须停止并重新启动后台，已经运行的Python不会自动载入新代码。网页显示服务的启动版本、代码目录和启动时间；这些标识在启动时记录，磁盘上执行`git pull`不会让旧进程的版本标识变新。若仍显示已修复的旧报错，查看服务标识并使用新端口启动更新后的目录，避免连到旧进程或旧任务。

以后只需停止旧进程，再执行 `bash 启动.sh`：启动器会自动检查并快进更新 main，然后启动最新代码。离线、认证失败或存在本地改动时使用现有版本；不强制覆盖、不自动合并冲突。跳过网络检查可用 `bash 启动.sh --no-update`。

首次克隆的权限与更新相同：公开仓库无需登录；私有仓库需要在本机先配置一次 GitHub 认证。ChatGPT 云端的平台授权不会自动传到WSL，浏览器GitHub登录也不是Git命令的认证。`git clone`用于首次获取；后续更新由启动器运行 `git pull --ff-only` 完成。安装依赖只在首次或缺少新依赖时执行。

Windows 双击入口 `启动工作站.bat` 针对你已安装的 `Ubuntu-24.04`。首次推荐从 Ubuntu 启动，便于看到报错。

## 软件与分类器

- 界面和文件检查：Python 3 标准库，无需 Node.js、Shiny 或第三方 Python 包。
- 分析：以用户已安装的 `qiime2-amplicon-2025.7` 为参数依据，WSL 内需有 dada2、cutadapt、feature-classifier、phylogeny、diversity、taxa、feature-table 等插件。其他版本先核对是否支持 `--p-max-merge-mismatch`。
- 分类器：在网页选择你自己的兼容 `.qza`。程序不下载、不替换分类器。`tools peek` 可以查看类型，但完整 sklearn/区域兼容性需实际分类时核验。
- 绘图：独立 conda 环境 `16s-workstation-r`，固定 R 4.4 并一次安装 ggplot2、patchwork、svglite、xml2、jsonlite、systemfonts、multcompView、MASS 的兼容二进制包。与QIIME环境分开，也不调用Windows R/RStudio。网页中的“绘图专用R环境名”默认为该环境。
- 选择你已安装的 QIIME 环境名。`conda env list` 可查看；环境检查会报告名称、插件、R依赖和分类器情况。没有选择环境也能返回发现结果。不同 QIIME 版本可能改变命令或文件布局，不能仅凭名称认定兼容。

如果缺R依赖，在项目目录执行 `bash 安装依赖.sh`。需要网络和 apt 权限；apt只准备Python/字体，脚本不再更新系统R，不修改QIIME环境。它从单一 conda-forge 渠道创建或修复 `16s-workstation-r`，保留TLS及包校验。后台通过 `conda run` 调用专用R，并在激活前后隔离系统/用户旧R包路径和启动配置。原先失败安装的 `~/R/16S-workstation-library` 保留，工作站不再加载。首次创建可能下载数百MB，具体随依赖变化。

绘图默认 Times New Roman；依赖脚本尝试复用已有Windows字体，并把系统DejaVu字体复制到工作站本地字体目录。没有Times时，请明确选择DejaVu Serif，不静默替换。

若早期安装日志出现 `SETLENGTH` / `SET_FRAME`、旧 `cpp11` 编译错误或旧 `gtable` 版本错误，先 `git pull --ff-only`，再重跑新版 `安装依赖.sh`。这是系统R与旧编译包的混用问题，不应删除原始数据或修改全局R配置来绕过。专用环境安装失败时不会退回混用系统R。

## 数据命名与分组

目前支持 Phred33 编码的 `.fastq`、`.fq` 和相应 gzip 压缩格式，单端或未拼接双端。推荐四列 CSV `SampleID,Group,forward,reverse`：`SampleID` 是本次分析的样本名，`Group` 是分组，`forward` 和 `reverse` 分别指定该样本的 R1 和 R2 文件。样本名直接按 CSV 中的文本读取，无需与文件名相同，前导零保留。

例如，选择的原始目录中有 `raw.A1.1.fq.gz` 和 `raw.A1.2.fq.gz`，可以填写：

```csv
SampleID,Group,forward,reverse
A1,CON,raw.A1.1.fq.gz,raw.A1.2.fq.gz
A2,CON,raw.A2.1.fq.gz,raw.A2.2.fq.gz
CAP1,LDC,CAP1.raw_1.fastq,CAP1.raw_2.fastq
CAP2,LDC,CAP2.raw_1.fastq,CAP2.raw_2.fastq
```

以上是格式示例，不代表实际测序数据。仓库中的 `分组模板.csv` 可用作四列模板，替换为你的样本、分组与真实文件名。`forward` / `reverse` 可填写文件名、相对于所选原始目录的子目录路径（如 `batch1/001_R1.fastq.gz`），或该目录内的 Windows / WSL 绝对路径。路径含逗号时，按标准 CSV 格式用双引号包围单元格。

双端必须填写两端；单端选择“真正的单端 FASTQ”后，保留四列表头，仅填写 `forward`，`reverse` 留空。方向以列为准，不依赖文件名末尾的 R1/R2 标记。每行对应一个样本，每个方向对应一个文件，多 lane 不会自动合并。

四列模式只分析 CSV 中列出的文件，目录中未列出的 FASTQ 不纳入本次分析。程序检查重复样本 ID、空分组、文件缺失、同一文件被重复使用（包括不同路径别名），以及文件是否位于所选原始目录内；目录外的路径会被拒绝。点击“检查配对与引物”可核对每个样本的分组、方向和实际文件名。程序将按此映射生成 QIIME 2 manifest，以及只含样本和分组的 QIIME metadata，后续分析和绘图使用相同的样本关系。

Excel 编辑时将 `SampleID` 列设为文本，并另存为 CSV UTF-8；如果 Excel 已去掉前导零，工作站无法恢复。CSV 路径按原始文件目录解析，而不是按 CSV 所在目录解析。

旧版两列 `SampleID,Group` 仍兼容：程序根据分组表中的样本名自动匹配以下常见文件命名。

| 双端文件示例 | CSV SampleID |
|---|---|
| `raw.A1.1.fq.gz` / `raw.A1.2.fq.gz` | `A1` |
| `CAP1.raw_1.fastq` / `CAP1.raw_2.fastq` | `CAP1` |
| `001_R1.fastq.gz` / `001_R2.fastq.gz` | `001` |
| `001_S7_L001_R1_001.fastq.gz` / 同样格式的 R2 | `001` |

两列模式使用分组表中的样本名验证匹配；识别 raw 标记和末尾方向号。真实样本名中的数字保留。如果同名可对应多个CSV样本，则报告歧义；重复方向、缺少另一端、多lane不静默合并。单端也支持 `raw.A1.fq.gz` 或 `A1.raw.fastq` 对应 `A1`。对于不符合命名规则的文件，使用四列映射直接指定即可。

旧版两列示例（已保留在 `分组模板_两列.csv`）：

```csv
SampleID,Group
A1,CON
A2,CON
CAP1,LDC
CAP2,LDC
```

目录选择器跳过无权读取的 Windows 系统项目，如 `DumpStack.log.tmp`，不修改系统权限。FASTQ本身若不可读会报告错误。

## 分阶段运行

1. **输入和环境检查**：选择原始目录、分组CSV、输出根目录（允许原始目录下新建 results 子目录），确认配对和软件依赖。可以输入 Windows 路径或 `/mnt/d/...`，支持中文目录。
2. **引物与质量**：前1000条 reads 的开头匹配常见338F/341F/806R（考虑简并碱基、偏移和1错配）。338F/341F有共有序列，未匹配不能证明已去除。用户确认已去引物、选择候选组合或输入真实序列后，生成 `demux.qzv` / `trimmed.qzv`，并导出各阶段的 `per-sample-fastq-counts.tsv`、`forward-seven-number-summaries.tsv` 和双端 reverse 表。cutadapt 使用0.1错误率、至少较短引物长度的匹配，丢弃未匹配引物的 reads；报告保留前后各样本数量。
3. **预处理**：默认从进入DADA2的数据质量表自动推荐截断、trim、线程、最低丰度；支持切换手动，手动修改会被保留。预计扩增子长度改为可选：有可靠的去引物后长度时提前核对至少12bp重叠，未知时不编造长度，DADA2仍按实际序列检查拼接。它不是250/300bp测序read长度。继续预处理复用已经检查的导入数据；改变输入文件或元数据会阻止旧结果复用。旧版质量目录可点击“从已有质量结果导出表格并推荐参数”，不重新导入或运行DADA2。
4. **抽平和多样性**：曲线最大深度可留空，按过滤后最大样本总量生成并记录实际值。预处理完成后自动推荐正式抽平深度，显示曲线依据及每组将保留/丢弃的样本，支持手动修改。已有预处理结果可点击“读取本批结果并推荐抽平深度”，无需重跑DADA2。每次核心分析只抽平一次，之后指标、统计及R图使用本次保存数据。
5. **查看结果**：界面按阶段折叠显示，支持筛选、下载和SVG/PDF/PNG预览；各步生成 `说明.md` / `report.json`，可预览说明和TSV/CSV/JSON。输入检查与环境检查在已选择的输出根目录保存独立check目录。QZV可在QIIME 2 View中查看；该查看器需联网，本工作站不内置完整QZV查看器。

## 质控依据与自动参数

用户提供的Q20/10bp滑窗、10bp overlap、20%拼接错配、USEARCH/gold是另一套处理流程；本工作站采用DADA2 ASV流程。按照QIIME 2 2025.7官方实现，明确使用正反向预期错误数上限2、`trunc_q=2`、最小重叠12bp、最大拼接错配0、`chimera_method=consensus`。DADA2在降噪后拼接、检测并去除嵌合体，原流程的默认值已有去嵌合体；新版将设置和实际统计明确保存，不额外叠加USEARCH或gold数据库。输入须是已拆分到样本的FASTQ，CSV关联样本与文件，不在这里重新按barcode拆样；需要原始barcode拆样的数据应先完成相应步骤。

自动截断是透明的保守启发式，**不是官方固定最优值**：每端先取至少99%抽样reads能覆盖的最长位置，再在第一个连续10个位置的Q25中位数低于Q30的窗口前截断。这个窗口只用于评估汇总质量图，不对每条read执行滑窗修剪。推荐长度不足50bp时要求人工判断，50bp只是自动推荐的适用范围，不是新增的DADA2通用过滤标准。记录抽样数、长度覆盖上限、质量窗口、选定截断和重叠检查；抽样长度覆盖比例不能当成实际质控保留率。

引物处理确认后，额外trim默认0；线程推荐不超过可用CPU的4线程；ASV全表最低丰度默认1，避免额外删除低丰度ASV影响丰富度与Chao1，手动可提高。曲线最大深度留空时仍由过滤后最大样本计数确定；正式抽平深度须结合曲线和各组样本保留情况选择。质量汇总无法推断真实扩增子长度、引物是否已去除或每条read的预期错误数。参数效果以各样本DADA2的filtered、merged、non-chimeric统计核验。

官方依据：[DADA2 2025.7参数及默认值](https://github.com/qiime2/q2-dada2/blob/2025.7.0/q2_dada2/_denoise.py)、[参数含义](https://github.com/qiime2/q2-dada2/blob/2025.7.0/q2_dada2/plugin_setup.py)、[质量统计与TSV导出](https://github.com/qiime2/q2-demux/blob/2025.7.0/q2_demux/_summarize/_visualizer.py)。

正式分组检验至少要求每组保留2个样本，且不允许某组因抽平全部丢失。该运行检查不替代研究设计。已拼接 reads 不能自动还原成双端。

## 正式抽平深度推荐与独立α图

正式深度是每个样本随机保留的有效序列数，依据DADA2及后续过滤后的计数；不是read长度，也不是曲线最大深度。低于所选深度的样本不参与本次核心分析。

推荐先要求每组至少保留90%的原样本且不少于2个，然后读取QIIME抽平曲线CSV中每个样本在各深度、各重复的指标均值，寻找最早的三点窗口：末点深度至少为首点的1.2倍；观测ASV与Faith PD变化≤5%、Shannon变化≤2%；每组至少90%原样本支持该窗口。上述比例和变化阈值是工作站启发式，不是官方统一最佳深度，也不是平台的统计显著性检验。

若没有兼顾曲线和留样的窗口，返回满足组内留样的最大深度，并明确标为“仅按样本保留给出的初步推荐”；缺少曲线时不宣称已确认平台。无法达到留样条件时不编造值。新run在 `06_rarefaction/curve_data/` 保存原始曲线CSV，推荐另存 `depth_recommendation.json` 并写入阶段报告；旧QZV也能直接读取，无需重新抽平。推荐不读取分组检验P值。

R把七项已计算α指标分别保存为 `alpha/observed_features`、`shannon`、`evenness`、`faith_pd`、`chao1`、`simpson`、`goods_coverage` 的SVG/PDF/PNG，每项有对应输入及注释CSV。全部指标不再合在一个α图里；缺少指标会明确报告，已有统计不可读取时仍不生成虚构显著性字母。

v0.4修复路径转换把 `qiime diversity alpha` 的动作名误改成 `03_statistics/alpha` 的错误：只转换输入/输出和元数据文件参数，动作名及普通参数值不转换；输出文件父目录自动创建，QIIME要求不存在的 `--output-dir` 目标由动作自行创建。已完成的预处理可继续用于新的多样性任务，失败的旧任务保留。空数值会提示具体参数名称，如“请填写本次抽平深度”，不再只显示Python的int转换错误。

## 保存的分析结果

本版执行以下流程，保存**每一步实际生成的所有 QZA/QZV 和中间结果**，并非运行QIIME全部可选方法：

- 导入、质量概览、可选cutadapt、DADA2表/代表序列/降噪统计。
- 用户分类器sklearn分类和分类结果可视化。
- 排除 Archaea、mitochondria、chloroplast，ASV全表最低总丰度默认1；过滤序列、分类堆积图、MAFFT/FastTree树。最低丰度可在网页修改；当前排除项固定，若研究包含古菌需另行调整。
- 抽平曲线（Shannon、observed features、Faith PD）；核心抽平表。
- 七项α：observed features、Shannon、evenness、Faith PD、Chao1、Simpson、Good’s coverage；QIIME总体/两两KW和保存的BH q值。
- 四项β：weighted/unweighted UniFrac、Bray-Curtis、Jaccard，保存距离、PCoA、Emperor及总体/两两PERMANOVA与PERMDISP。
- BIOM/TSV、分类TSV、序列FASTA、树、alpha向量、PCoA、距离及统计QZV的标准导出文件。
- R的α箱线图、四项PCoA和组合图、门/属的每样本组成与组均值图，输出SVG/PDF/PNG及绘图输入CSV。组成图使用**未抽平的过滤表**，先按每样本归一化，再求组均值；Top15+Other仅用于图，全表也保存。

已内置用户提供的QIIME结果读取/原始风格函数；无需本机额外保存原来的3个脚本。若已有自定义版本，可选择脚本目录覆盖内置函数。缺失或无法验证的保存统计不补造P值或显著性字母。R不重新抽平、计算PCoA、运行置换检验或再次BH校正已有q值。

`workflow.log`、`commands.jsonl`、`parameters.json`、输入清单、样本保留情况、软件版本和 `outputs.tsv` 供复核。QIIME本身在QZA/QZV中保存provenance。

如果某指标全部相同、两组值完全相同，或全部距离为0，相应检验会记录明确跳过原因；其他未经诊断的命令失败会停下并保留产物。统计状态保存在 `statistics_status.json`。R失败不会删除已完成的QIIME结果。随机抽平/置换重跑可能不同，因此保留每次产物，不能宣称固定种子复现。

## 目录和继续工作

输出根目录可以是原始数据目录内的新子目录，例如 `D:\16S\raw\results`。工作站每次在其下新建唯一的 `run_<时间>_<编号>`，不修改原始FASTQ或覆盖旧结果。仍要求输出不是原始目录本身或其父目录；结果目录中的文件不参与两列模式的原始文件自动识别，四列模式仅读取CSV指定文件。

```text
输出根目录/run_<时间>_<编号>/
  00_inputs/              manifest、metadata、输入说明
  01_import/              demux QZA/QZV、原始质量统计表
  02_primers/             去引物 QZA/QZV、处理后质量统计表
  03_dada2/               ASV、代表序列、去嵌合体统计与完整ASV表
  04_taxonomy_filter/     分类、过滤后表、代表序列、分类堆积图
  05_phylogeny/           比对、掩蔽、树及树导出
  06_rarefaction/         多深度抽平曲线
  08_parameter_advice/    自动推荐依据
  07_diversity/diversity_<深度>_<编号>/
    00_inputs/            本次保留样本分组与retention
    01_core_metrics/      本次抽平表、指标、距离、PCoA、Emperor
    02_exports/           抽平TSV、指标、坐标及统计导出
    03_statistics/        分组检验QZV及执行/跳过状态
    04_R_figures/         按alpha/beta/taxonomy分组的SVG/PDF/PNG/CSV
  处理说明.md、参数、布局映射、运行状态、日志
```

每个阶段目录都有 `说明.md` 和 `report.json`，记录参数、实际数量、保存文件和完成/跳过/失败状态。旧版已有结果保留原位置，新程序通过有版本号的布局映射兼容读取；不搬动或覆盖旧分析。

本版本不是任意步骤断点恢复工具：质量阶段可继续到预处理，完成预处理可创建多次不同深度分析；失败的部分结果不直接重用。重启网页后台后原文件保留，内存任务状态丢失；使用历史结果目录查看文件或开始下一阶段。

## 验证

详见 `VALIDATION.md`。后端/浏览器功能和真实R绘图用合成数据验证，QIIME命令编排测试采用模拟执行器。**未在你的电脑或真实FASTQ上完成端到端分析。** 首次使用先环境检查，再逐阶段运行；实际环境失败应依据日志处理，不把源码交付视为真实数据分析成功。
