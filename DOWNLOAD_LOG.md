# PISA 2025 微观数据下载日志

> 论文包路径 B：下载 OECD 官方学生问卷公用文件，用于 SSCI 论文实证分析前置准备。
> 本日志如实记录实际下载与核验结果，未下载成功的项目不做任何成功标注。

## 执行信息

| 项目 | 内容 |
|---|---|
| 执行日期 | 2026-11-17（GMT+8） |
| 主数据下载开始 | 14:29:21 +0800 |
| 主数据下载结束 | 14:52:15 +0800 |
| 主数据下载耗时 | 22 分 53 秒 |
| 目标目录 | `C:/Users/chioerik/WorkBuddy/2026-11-17-12-50-12/ai_frequency_paradox/data` |
| 下载工具 | curl 8.22.0（断点续传 `-L -C -`，失败重试 3 次），unzip 6.00 |
| 网络中断 | 无。首次下载即成功（CURL_EXIT=0），未触发重试 |

## 一、下载成功的文件

### 1. 学生问卷公用文件（主数据）

| 项目 | 内容 |
|---|---|
| 文件名 | `CY09_MS_STU_PUF.zip` |
| 绝对路径 | `C:/Users/chioerik/WorkBuddy/2026-11-17-12-50-12/ai_frequency_paradox/data/CY09_MS_STU_PUF.zip` |
| 来源 URL | `https://webfs.oecd.org/pisa2022/2025/CY09_MS_STU_PUF.zip` |
| 大小（ls -lh） | 947M |
| 精确字节数 | 992,910,196 bytes |
| 服务器声明大小 | 992,910,196 bytes（HEAD 探测 Content-Length，完全一致） |
| 服务器元信息 | HTTP 200；Content-Type: application/x-zip-compressed；Last-Modified: Mon, 24 Aug 2026 13:04:09 GMT；Accept-Ranges: bytes；ETag: "814a1cdc933dd1:0" |
| 下载命令 | `curl -L -C - -f --retry 3 --retry-delay 5 --retry-all-errors --progress-bar -o CY09_MS_STU_PUF.zip <URL>` |
| curl 退出码 | 0 |
| zip 完整性 | `unzip -t` → `testing: CY09_MS_STU_PUF.sav OK`；`No errors detected in compressed data`；UNZIP_EXIT=0 |
| 压缩包内容 | 含 `CY09_MS_STU_PUF.sav`（SPSS 格式学生问卷数据，内含水问卷/家长问卷数据、学生成绩估计值；含 AI 使用频率题、AI 素养机会题、ESCS 与 5 个科学 plausible values） |
| 校验状态 | **通过**（字节数精确匹配 + zip 结构完整） |

### 2. 官方 Codebook

| 项目 | 内容 |
|---|---|
| 文件名 | `PISA2025_Codebook.xlsx` |
| 绝对路径 | `C:/Users/chioerik/WorkBuddy/2026-11-17-12-50-12/ai_frequency_paradox/data/PISA2025_Codebook.xlsx` |
| 来源 URL | `https://webfs.oecd.org/pisa2022/2025/PISA2025_Codebook.xlsx`（OECD 官方文件服务器，HEAD 探测 200 / 985,161 bytes） |
| 大小（ls -lh） | 963K |
| 魔数校验 | `PK\x03\x04`（有效 xlsx/zip） |
| 工作表 | STU、COG、LDW、SCH、TCH、STU_TT、COG_PROCESS |
| 关键变量核验 | AI 使用频率 `ST438Q01DA`–`ST438Q04DA`（4 项）；AI 素养机会 `ST437Q01DA`–`ST437Q02DA`（2 项）；AI 工具使用 `IC171Q01JA` 等（10 项）；AI 兴趣/素养 `IC400Q01DA` 等（6 项）；`ESCS`；5 个科学 plausible values `PV1SCIE`–`PV5SCIE`（全部齐全） |
| 校验状态 | **通过**（文件可读 + 关键变量全部命中） |

### 3. 学生问卷文档（官方 PDF，降级方案来源）

OECD 数据库页面为 AEM 动态渲染且主站有 Cloudflare 挑战，无法提取问卷直链；改用 NCES（美国国家教育统计中心）官方发布的 PISA 2025 主调查英文问卷作为文档来源。

| 文件名 | 绝对路径 | 来源 URL | 大小 | 校验 |
|---|---|---|---|---|
| `PISA2025_Student_Questionnaire_NCES.pdf` | `C:/Users/chioerik/WorkBuddy/2026-11-17-12-50-12/ai_frequency_paradox/data/PISA2025_Student_Questionnaire_NCES.pdf` | `https://nces.ed.gov/surveys/pisa/pisa2025/questionnaires/P25_StuQ-USA.pdf` | 700K | PDF 1.6 头；首页文本确认为 "Student QUESTIONNAIRE / National Version / PISA 2025 / Main Survey / March, 2026" |
| `PISA2025_ICT_Questionnaire_NCES.pdf` | `C:/Users/chioerik/WorkBuddy/2026-11-17-12-50-12/ai_frequency_paradox/data/PISA2025_ICT_Questionnaire_NCES.pdf` | `https://nces.ed.gov/surveys/pisa/pisa2025/questionnaires/P25_ICTQ-USA.pdf` | 531K | PDF 1.7 头，文件有效 |
| `PISA2025_Student_UH_Questionnaire_NCES.pdf` | `C:/Users/chioerik/WorkBuddy/2026-11-17-12-50-12/ai_frequency_paradox/data/PISA2025_Student_UH_Questionnaire_NCES.pdf` | `https://nces.ed.gov/surveys/pisa/pisa2025/questionnaires/P25_StuQ-UH-USA.pdf` | 579K | PDF 1.7 头，文件有效 |

## 二、未下载成功 / 未获取的项目（如实记录）

| 项目 | 状态 | 原因 |
|---|---|---|
| PISA 2025 Technical Report 草稿章节 | **未获取** | OECD 数据库页面（https://www.oecd.org/en/data/datasets/pisa-2025-database.html）虽有 "Download the draft technical report" 入口，但页面为动态渲染，curl 直取返回 Cloudflare "Just a moment..." 挑战页（5512 字节），WebFetch 与只读文本代理均未返回该链接的 href；文件服务器 `/2025/` 目录列表返回 HTTP 403；按官方命名惯例探测 9 个候选文件名（PISA2025_Technical_Report.zip 等）全部 HTTP 404。已按任务降级方案下载学生问卷 PDF 三份。 |
| 认知过程数据 COG_PROCESS | 未下载（任务明确排除） | — |
| LDW 数据 | 未下载（任务明确排除） | — |
| 学校问卷 SCH | 未下载（任务明确排除） | — |
| 教师问卷 TCH | 未下载（任务明确排除） | — |

> 非致命性探测错误：若干 codebook / 技术报告候选 URL 返回 404/403，属于定位直链过程中的探测响应，不是主数据下载失败。主数据下载本身无错误、无重试。

## 三、核验汇总

- `ls -lh`：5 个文件均已落盘，主数据 947M，总占用约 950M。
- 字节级比对：主数据 992,910,196 bytes 与服务器 Content-Length 完全一致。
- `unzip -t CY09_MS_STU_PUF.zip`：`No errors detected in compressed data of CY09_MS_STU_PUF.zip.`
- Codebook 关键变量（AI 使用频率、AI 素养机会、ESCS、5 个科学 PV）已逐一核对存在。

## 四、下一步

1. 本机需安装 R（及 RStudio），并安装 `haven` / `readxl` / `dplyr` 等包。
2. 按 `PISA2025_Codebook.xlsx` 的 STU 工作表建立变量映射：
   - AI 使用频率：`ST438Q01DA`–`ST438Q04DA`
   - AI 素养机会：`ST437Q01DA`–`ST437Q02DA`；AI 工具：`IC171Q*`、`IC400Q*`
   - 控制变量：`ESCS`
   - 因变量：5 个科学 plausible values `PV1SCIE`–`PV5SCIE`（按 OECD 推荐用重复估计合并，不能只取单个 PV）
3. 使用 `haven::read_sav()` 读取 `CY09_MS_STU_PUF.sav`，按官方 codebook 处理反向题、缺失值与权重（学生权重 W_FSTUWT + 重复权重），再运行后续分析脚本。
4. 本轮仅完成下载与核验，未修改论文正文，未运行 R。

## 五、复核轮记录（2026-11-17 15:54 +0800，自动化再次触发）

> 本轮为自动化任务再次触发。因上一轮已成功下载完整文件，为避免对 993MB 完整副本重复拉取（`curl -C -` 对完整文件只会返回 HTTP 416「Range Not Satisfiable」，无实际意义），本轮执行落盘复核而非重新下载。所有结论均来自对本机真实文件的检测，未伪造。

| 复核项 | 结果 |
|---|---|
| 主数据本地字节数 | 992,910,196 bytes（`stat` 实测） |
| 服务器 HEAD Content-Length | 992,910,196 bytes（HTTP 200，Last-Modified: Mon, 24 Aug 2026 13:04:09 GMT，Accept-Ranges: bytes） |
| 字节级比对 | **一致** |
| 主数据 zip 完整性 | `unzip -t` → `testing: CY09_MS_STU_PUF.sav OK`；`No errors detected in compressed data`；退出码 0 → **通过** |
| Codebook 有效性 | 魔数 `50 4b 03 04`（有效 xlsx）；7 个工作表 STU/COG/LDW/SCH/TCH/STU_TT/COG_PROCESS |
| Codebook 关键变量 | `ST438Q01DA`/`ST438Q04DA`/`ST437Q01DA`/`ST437Q02DA`/`IC171Q01JA`/`IC400Q01DA`/`ESCS`/`PV1SCIE`/`PV5SCIE` 全部 **FOUND** |
| 3 份问卷 PDF | `%PDF-1.6`、`%PDF-1.7`、`%PDF-1.7` 文件头有效 |
| 本轮下载动作 | 无（文件已完整且校验通过，跳过重复下载） |
| 本轮错误 | 无 |

**复核结论**：5 个文件全部完整有效，主数据与服务器字节级一致，可直接进入后续 R 分析。下一步仍需安装 R 与 `haven`/`readxl` 并按官方 codebook 映射变量（见第四节）。
