# Working preferences

The user prioritizes low AI usage and small, purpose-driven changes.

- Modify the existing implementation. Inspect only the relevant symbols and files; do not rediscover the whole project for each request.
- For small changes, use one agent. Do not spawn subagents by default. Delegate only when the user asks or a substantial task benefits enough to justify the added AI work.
- Keep the change within the requested scope. Do not add optional features, diagnostics, refactors, or documentation rewrites to a small adjustment.
- Run focused verification for the affected behavior. Reuse previous validation evidence; do not repeat successful checks without a new change, failure, or unresolved concern.
- Run the complete backend suite, browser tests, real R rendering, or dependency installation only when the affected workflow requires it. Avoid recreating environments or downloading packages that already work.
- Do not add tests for simple text, styling, or documentation edits. For data handling, retain checks that protect sample mapping, raw files, and saved analysis results.
- Keep progress updates and final answers concise: what changed, what was checked, and the local update command when needed.
- Never save time by weakening TLS/package verification, overwriting user files, inventing statistics, or claiming unexecuted analysis passed.

## Relevant files

- `app.py`: CSV/file mapping, QIIME commands, isolated R calls, local HTTP service.
- `index.html`: browser interface.
- `render.R` and `plots/`: plots from saved QIIME outputs.
- `launch.py` and `启动.sh`: startup and conservative source updates.
- `安装依赖.sh`: separate R environment; preserve the QIIME environment.
- `README.md` / `VALIDATION.md`: workflow and previous validation evidence.

- Extended modules: `abundance_tables.py` (ASV/L2–L7 exports), `differential_analysis.py` + `differential_ancombc.R` (explicit ANCOM-BC2/QIIME ANCOM), `lefse_analysis.py` (cached threshold plots), `functional_prediction.py` (PICRUSt2/KO maps). `app.extension` owns isolation and completion markers. Do not refit LEfSe for LDA previews or use relative/rarefied counts for ANCOM-BC2.
