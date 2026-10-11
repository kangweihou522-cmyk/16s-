# Windows / WSL 16S 本地分析工作站 v0.2

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

以后只需停止旧进程，再执行 `bash 启动.sh`：启动器会自动检查并快进更新 main，然后启动最新代码。离线、认证失败或存在本地改动时使用现有版本；不强制覆盖、不自动合并冲突。跳过网络检查可用 `bash 启动.sh --no-update`。

首次克隆的权限与更新相同：公开仓库无需登录；私有仓库需要在本机先配置一次 GitHub 认证。ChatGPT 云端的平台授权不会自动传到WSL，浏览器GitHub登录也不是Git命令的认证。`git clone`用于首次获取；后续更新由启动器运行 `git pull --ff-only` 完成。安装依赖只在首次或缺少新依赖时执行。

Windows 双击入口 `启动工作站.bat` 针对你已安装的 `Ubuntu-24.04`。首次推荐从 Ubuntu 启动，便于看到报错。

## 软件与分类器

- 界面和文件检查：Python 3 标准库，无需 Node.js、Shiny 或第三方 Python 包。
- 分析：WSL 内已有 QIIME 2 amplicon conda 环境，含 dada2、cutadapt、feature-classifier、phylogeny、diversity、taxa、feature-table 等插件。
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

1. **输入和环境检查**：选择原始目录、分组CSV、与原始目录分开的输出根目录，确认配对和软件依赖。可以输入 Windows 路径或 `/mnt/d/...`，支持中文目录。
2. **引物与质量**：前1000条 reads 的开头匹配常见338F/341F/806R（考虑简并碱基、偏移和1错配）。338F/341F有共有序列，未匹配不能证明已去除。用户确认已去引物、选择候选组合或输入真实序列后，生成 `demux.qzv` / `trimmed.qzv`。cutadapt 模式丢弃未匹配引物的 reads，数量保留在日志。
3. **预处理**：依据质量图填写DADA2截断/trim参数；双端还需填写进入DADA2的预计扩增子长度，用于核对至少12bp重叠。去引物时填去除两条引物后的扩增子长度，不是250/300bp测序read长度。继续预处理会复用已经检查的导入数据；改变输入文件或元数据会阻止旧结果复用。
4. **抽平和多样性**：曲线最大深度可留空，按过滤后最大样本总量生成并记录实际值。结合曲线和样本深度，选择正式抽平深度；界面显示每组将保留/丢弃的样本。每次核心分析只抽平一次，之后指标、统计及R图使用本次保存数据。
5. **查看结果**：界面列出结果文件，支持筛选、下载和SVG/PDF/PNG预览。QZV可在QIIME 2 View中查看；该查看器需联网，本工作站不内置完整QZV查看器。

正式分组检验至少要求每组保留2个样本，且不允许某组因抽平全部丢失。该运行检查不替代研究设计。已拼接 reads 不能自动还原成双端。

## 保存的分析结果

本版执行以下流程，保存**每一步实际生成的所有 QZA/QZV 和中间结果**，并非运行QIIME全部可选方法：

- 导入、质量概览、可选cutadapt、DADA2表/代表序列/降噪统计。
- 用户分类器sklearn分类和分类结果可视化。
- 排除 Archaea、mitochondria、chloroplast，ASV全表最低总丰度默认10；过滤序列、分类堆积图、MAFFT/FastTree树。阈值和排除项依研究目标还可调整。
- 抽平曲线（Shannon、observed features、Faith PD）；核心抽平表。
- 七项α：observed features、Shannon、evenness、Faith PD、Chao1、Simpson、Good’s coverage；QIIME总体/两两KW和保存的BH q值。
- 四项β：weighted/unweighted UniFrac、Bray-Curtis、Jaccard，保存距离、PCoA、Emperor及总体/两两PERMANOVA与PERMDISP。
- BIOM/TSV、分类TSV、序列FASTA、树、alpha向量、PCoA、距离及统计QZV的标准导出文件。
- R的α箱线图、四项PCoA和组合图、门/属的每样本组成与组均值图，输出SVG/PDF/PNG及绘图输入CSV。组成图使用**未抽平的过滤表**，先按每样本归一化，再求组均值；Top15+Other仅用于图，全表也保存。

已内置用户提供的QIIME结果读取/原始风格函数；无需本机额外保存原来的3个脚本。若已有自定义版本，可选择脚本目录覆盖内置函数。缺失或无法验证的保存统计不补造P值或显著性字母。R不重新抽平、计算PCoA、运行置换检验或再次BH校正已有q值。

`workflow.log`、`commands.jsonl`、`parameters.json`、输入清单、样本保留情况、软件版本和 `outputs.tsv` 供复核。QIIME本身在QZA/QZV中保存provenance。

如果某指标全部相同、两组值完全相同，或全部距离为0，相应检验会记录明确跳过原因；其他未经诊断的命令失败会停下并保留产物。统计状态保存在 `statistics_status.json`。R失败不会删除已完成的QIIME结果。随机抽平/置换重跑可能不同，因此保留每次产物，不能宣称固定种子复现。

## 目录和继续工作

```text
输出根目录/run_<时间>_<编号>/
  QIIME预处理QZA/QZV、导出表、manifest、metadata、参数、日志
  diversity_<深度>_<编号>/
    core-metrics-results/   抽平表、指标、距离、PCoA、Emperor
    alpha/ beta/            分组检验QZV
    exports/                标准导出数据和统计可视化
    R_original_style/       SVG、PDF、PNG、CSV、R版本
    statistics_status.json  已执行/跳过/未完成状态
```

本版本不是任意步骤断点恢复工具：质量阶段可继续到预处理，完成预处理可创建多次不同深度分析；失败的部分结果不直接重用。重启网页后台后原文件保留，内存任务状态丢失；使用历史结果目录查看文件或开始下一阶段。

## 验证

详见 `VALIDATION.md`。后端/浏览器功能和真实R绘图用合成数据验证，QIIME命令编排测试采用模拟执行器。**未在你的电脑或真实FASTQ上完成端到端分析。** 首次使用先环境检查，再逐阶段运行；实际环境失败应依据日志处理，不把源码交付视为真实数据分析成功。
