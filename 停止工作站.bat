@echo off
chcp 65001 >nul
setlocal
if exist "%~dp0Station.ps1" (
    powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0Station.ps1" -Action Stop %*
) else (
    powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0windows\Workstation.ps1" -Action Stop %*
)
set "station_exit=%errorlevel%"
if not "%station_exit%"=="0" echo 工作站操作失败，请查看上方错误和工作站日志。
pause
exit /b %station_exit%
